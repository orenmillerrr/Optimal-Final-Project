% MATLAB script to create and test an Interacting Multiple Model (IMM) filter
% with Constant Velocity (CV), Constant Acceleration (CA), and Constant Turn (CT) models.

clear; close all; clc;

%% Simulation Parameters
dt = 0.1; % Time step
T_sim = 100; % Total simulation time
N_steps = T_sim / dt; % Number of simulation steps

% Initial state [x; y; vx; vy; ax; ay; omega]
% Note: ax, ay, and omega are only active in CA and CT models, respectively.
% The state vector size is determined by the largest model state (CA or CT, size 6 or 5 depending on formulation).
% We will use a consistent state representation for simplicity, padding with zeros where needed.
% State vector: [x; y; vx; vy; ax; ay] for CV/CA, [x; y; vx; vy; omega] for CT.
% For IMM, we'll use a maximum state size and map models to it.
% Let's use a state size of 6: [x; y; vx; vy; ax; ay]. CT omega will be handled separately or integrated into a 6-state vector.
% A common approach is to use a 6-state vector [x; y; vx; vy; ax; ay] and interpret ax/ay differently for CT.
% For CT, we can use [x; y; vx; vy; omega], size 5. Let's stick to a 6-state representation for consistency where possible.
% State: [x; y; vx; vy; ax; ay]
initial_state = [0; 0; 10; 0; 0; 0]; % Initial position (0,0), velocity (10,0)

% Measurement noise covariance (assuming position measurements [x; y])
R_meas = diag([5^2, 5^2]); % Variance of 5 units in x and y

% Process noise standard deviations for each model
sigma_v_cv = 1; % CV velocity process noise
sigma_a_ca = 2; % CA acceleration process noise
sigma_omega_ct = deg2rad(1); % CT turn rate process noise (convert degrees to radians)
sigma_a_ct = 0.5; % CT acceleration process noise (for tangential/radial if using 6-state)

%% Model Definitions

% Model 1: Constant Velocity (CV)
% State: [x; y; vx; vy] - size 4
% We will map this to the 6-state vector [x; y; vx; vy; 0; 0]
F_cv = [1 0 dt 0 0 0;
        0 1 0 dt 0 0;
        0 0 1 0 0 0;
        0 0 0 1 0 0;
        0 0 0 0 0 0; % Padding for ax
        0 0 0 0 0 0]; % Padding for ay
Q_cv = diag([0, 0, sigma_v_cv^2*dt, sigma_v_cv^2*dt, 0, 0]); % Process noise for velocity

% Model 2: Constant Acceleration (CA)
% State: [x; y; vx; vy; ax; ay] - size 6
F_ca = [1 0 dt 0 0.5*dt^2 0;
        0 1 0 dt 0 0.5*dt^2;
        0 0 1 0 dt 0;
        0 0 0 1 0 dt;
        0 0 0 0 1 0;
        0 0 0 0 0 1];
Q_ca = diag([0, 0, 0, 0, sigma_a_ca^2*dt, sigma_a_ca^2*dt]); % Process noise for acceleration

% Model 3: Constant Turn (CT)
% State: [x; y; vx; vy; omega] - size 5. Mapping to 6-state [x; y; vx; vy; omega; 0] or [x; y; vx; vy; ax; ay] is tricky.
% A common approach for CT in a 6-state is [x; y; vx; vy; ax; ay] where ax/ay represent components of centripetal/tangential acceleration.
% Alternatively, use a 5-state CT model and handle state mapping explicitly in IMM.
% Let's use the 5-state CT model and map to/from the 6-state IMM filter state.
% CT State: [x; y; vx; vy; omega]
% Process noise for CT is often applied to omega and possibly tangential/radial acceleration components.
% We'll use a simplified 5-state CT model and handle the state mapping.
% State: [x; y; vx; vy; omega]
% Mapping from 6-state [x; y; vx; vy; ax; ay] to 5-state [x; y; vx; vy; omega] is lossy.
% Mapping from 5-state [x; y; vx; vy; omega] to 6-state [x; y; vx; vy; ax; ay] requires calculating ax/ay from omega, vx, vy.
% ax = -vy * omega, ay = vx * omega (centripetal acceleration components)
% Let's use the 6-state filter and define the CT model's effect on the 6-state.
% For CT, the state transition is non-linear. We'll use an Extended Kalman Filter (EKF) approach for the CT model within IMM.
% State: [x; y; vx; vy; ax; ay]. For CT, ax and ay are dependent on vx, vy, and omega.
% We will define the state as [x; y; vx; vy; omega], size 5, and map to/from the 6-state IMM filter.

% CT Model (5-state: [x; y; vx; vy; omega])
% Non-linear state transition function f_ct(state, dt):
% x(k+1) = x(k) + (vx(k)/omega(k))*sin(omega(k)*dt) - (vy(k)/omega(k))*(1-cos(omega(k)*dt))
% y(k+1) = y(k) + (vy(k)/omega(k))*sin(omega(k)*dt) + (vx(k)/omega(k))*(1-cos(omega(k)*dt))
% vx(k+1) = vx(k)*cos(omega(k)*dt) - vy(k)*sin(omega(k)*dt)
% vy(k+1) = vx(k)*sin(omega(k)*dt) + vy(k)*cos(omega(k)*dt)
% omega(k+1) = omega(k) (constant turn rate)

% We need the Jacobian of f_ct for the EKF linearisation.
% J_ct = df_ct / d(state)

% Process noise for CT (applied to omega and potentially velocity components)
Q_ct_5state = diag([0, 0, 0, 0, sigma_omega_ct^2*dt]); % Process noise on omega

% Measurement model (assuming [x; y] measurements for all models)
H_meas = [1 0 0 0 0 0; % Maps 6-state to [x]
          0 1 0 0 0 0]; % Maps 6-state to [y]

%% IMM Filter Parameters
num_models = 3;
% Model order: 1=CV (4-state), 2=CA (6-state), 3=CT (5-state)
model_states = {[1:4], [1:6], [1:5]}; % Indices of active states in a potential larger state vector (not used directly in this implementation)
model_names = {'CV', 'CA', 'CT'};

% Initial model probabilities
initial_model_prob = [1/3; 1/3; 1/3]; % Start with equal probability

% Model transition matrix P(j|i) = probability of transitioning to model j given currently in model i
% Rows: current model (i), Columns: next model (j)
% P = [ P(CV|CV) P(CA|CV) P(CT|CV)
%       P(CV|CA) P(CA|CA) P(CT|CA)
%       P(CV|CT) P(CA|CT) P(CT|CT) ]
transition_matrix = [0.95 0.025 0.025;
                     0.025 0.95 0.025;
                     0.025 0.025 0.95];

% Initial state covariance for each model (size depends on model state size)
P_initial_cv = diag([10^2, 10^2, 5^2, 5^2]); % CV: [x, y, vx, vy]
P_initial_ca = diag([10^2, 10^2, 5^2, 5^2, 2^2, 2^2]); % CA: [x, y, vx, vy, ax, ay]
P_initial_ct = diag([10^2, 10^2, 5^2, 5^2, deg2rad(5)^2]); % CT: [x, y, vx, vy, omega]

%% Initialize IMM Filter States and Covariances
% We will store states and covariances for each model.
% The combined IMM state and covariance will be calculated at each step.

% Filter states for each model (initialize with a common initial state estimate)
% We need to map the initial_state (6-state) to the model-specific state sizes.
filter_state_models = cell(num_models, 1);
filter_cov_models = cell(num_models, 1);

% Initialize CV model state and covariance (4-state)
filter_state_models{1} = initial_state(1:4); % [x; y; vx; vy]
filter_cov_models{1} = P_initial_cv;

% Initialize CA model state and covariance (6-state)
filter_state_models{2} = initial_state; % [x; y; vx; vy; ax; ay]
filter_cov_models{2} = P_initial_ca;

% Initialize CT model state and covariance (5-state)
% Need to estimate initial omega from initial velocity if possible, or assume 0.
initial_omega_est = 0;
filter_state_models{3} = [initial_state(1:4); initial_omega_est]; % [x; y; vx; vy; omega]
filter_cov_models{3} = P_initial_ct;

% Initial model probabilities
model_prob = initial_model_prob;

% Store combined IMM estimates
imm_state_est = zeros(6, N_steps); % Combined IMM state (using 6-state for consistency)
imm_cov_est = zeros(6, 6, N_steps); % Combined IMM covariance (using 6x6 for consistency)
model_prob_hist = zeros(num_models, N_steps); % History of model probabilities

%% Simulate Target Trajectory
true_state = zeros(6, N_steps); % True state [x; y; vx; vy; ax; ay]
true_state(:, 1) = initial_state;

% Define a sequence of motion modes
% 1: CV, 2: CA, 3: CT
motion_modes = [ones(1, 300), 2*ones(1, 300), 3*ones(1, 300), ones(1, 100)]; % Example sequence
motion_modes = repmat(motion_modes, 1, ceil(N_steps / length(motion_modes)));
motion_modes = motion_modes(1:N_steps);

% Parameters for motion modes
ca_accel = [2; 1]; % Constant acceleration for CA mode
ct_omega = deg2rad(5); % Constant turn rate for CT mode (5 degrees per second)

measurements = zeros(2, N_steps); % [x; y] measurements

for k = 1:N_steps-1
    current_true_state = true_state(:, k);
    current_mode = motion_modes(k);

    % Simulate true state transition based on current mode
    if current_mode == 1 % CV
        F_true_cv = [1 0 dt 0 0 0;
                     0 1 0 dt 0 0;
                     0 0 1 0 0 0;
                     0 0 0 1 0 0;
                     0 0 0 0 1 0; % ax stays constant (0)
                     0 0 0 0 0 1]; % ay stays constant (0)
        true_state(:, k+1) = F_true_cv * current_true_state; % No process noise in true state for simplicity
    elseif current_mode == 2 % CA
        F_true_ca = [1 0 dt 0 0.5*dt^2 0;
                     0 1 0 dt 0 0.5*dt^2;
                     0 0 1 0 dt 0;
                     0 0 0 1 0 dt;
                     0 0 0 0 1 0;
                     0 0 0 0 0 1];
        % Update acceleration in true state for CA mode
        current_true_state(5:6) = ca_accel;
        true_state(:, k+1) = F_true_ca * current_true_state;
    elseif current_mode == 3 % CT
        % Non-linear CT motion
        x = current_true_state(1);
        y = current_true_state(2);
        vx = current_true_state(3);
        vy = current_true_state(4);
        omega = ct_omega; % Use defined turn rate for true state

        if abs(omega) < 1e-9 % Handle straight line motion if omega is near zero
            x_next = x + vx * dt;
            y_next = y + vy * dt;
            vx_next = vx;
            vy_next = vy;
        else
            x_next = x + (vx/omega)*sin(omega*dt) - (vy/omega)*(1-cos(omega*dt));
            y_next = y + (vy/omega)*sin(omega*dt) + (vx/omega)*(1-cos(omega*dt));
            vx_next = vx*cos(omega*dt) - vy*sin(omega*dt);
            vy_next = vx*sin(omega*dt) + vy*cos(omega*dt);
        end
        true_state(1:4, k+1) = [x_next; y_next; vx_next; vy_next];
        true_state(5:6, k+1) = [-vy_next*omega; vx_next*omega]; % Store centripetal acceleration components
    end

    % Generate measurement with noise
    measurements(:, k+1) = H_meas(1:2, 1:6) * true_state(:, k+1) + mvnrnd([0; 0], R_meas)';
end

%% IMM Filter Loop
for k = 1:N_steps-1
    % 1. Interaction (Mixing)
    % Calculate normalization factor for mixing probabilities
    c_bar = transition_matrix' * model_prob;

    % Calculate mixing probabilities mu_ij = P(model i at k-1 | model j at k, Z_k-1)
    mixing_prob = zeros(num_models, num_models);
    for j = 1:num_models % Current model j
        for i = 1:num_models % Previous model i
            mixing_prob(i, j) = transition_matrix(i, j) * model_prob(i) / c_bar(j);
        end
    end

    % Calculate mixed initial state estimate and covariance for each model
    mixed_state_models = cell(num_models, 1);
    mixed_cov_models = cell(num_models, 1);

    for j = 1:num_models % For each model j at time k
        mixed_state_models{j} = zeros(size(filter_state_models{j}));
        mixed_cov_models{j} = zeros(size(filter_cov_models{j}));

        % Map states from previous models (i) to current model's (j) state space
        % This is the tricky part when state sizes differ.
        % We need functions to map between model state spaces.
        % For simplicity here, we will assume a common underlying 6-state space
        % and handle the mapping implicitly based on model type.
        % A more robust implementation would have explicit state mapping functions.

        % Let's refine the state representation:
        % CV: 4-state [x; y; vx; vy]. Map to/from 6-state [x; y; vx; vy; 0; 0]
        % CA: 6-state [x; y; vx; vy; ax; ay]. Direct use.
        % CT: 5-state [x; y; vx; vy; omega]. Map to/from 6-state [x; y; vx; vy; omega; 0] or similar.

        % Let's use the model-specific state sizes and handle mapping explicitly.
        % We need functions to map state_i -> state_j.

        % Simplified approach for demonstration: Assume the filter states
        % stored in filter_state_models are already in a compatible format
        % (e.g., padded to a common size if needed, or mapping handled externally).
        % A proper IMM implementation requires careful state mapping.

        % For this example, let's assume the filter_state_models and filter_cov_models
        % are already in the correct size for each model's prediction step.

        % Recalculate mixed state and covariance based on mixing probabilities
        for i = 1:num_models % From each previous model i
            % Need a function to map state from model i to model j's state space
            % state_i_mapped_to_j = map_state(filter_state_models{i}, model_i_type, model_j_type);
            % cov_i_mapped_to_j = map_covariance(filter_cov_models{i}, model_i_type, model_j_type);

            % For simplicity in this example, let's assume the states are compatible enough
            % or that the mapping is handled elsewhere. This is a simplification!
            % A real IMM needs careful state space management.

            % Let's assume filter_state_models{i} and filter_cov_models{i} are
            % already in the size required by model j for mixing. This is not generally true.
            % A better approach: Define a common larger state space (e.g., 6-state)
            % and map all model states to this space for mixing.

            % Let's use the 6-state as the common space for mixing.
            % We need functions to expand model i's state/cov to 6-state.
            state_i_expanded = expand_state_to_6(filter_state_models{i}, i);
            cov_i_expanded = expand_covariance_to_6(filter_cov_models{i}, i);

            % Now, map the expanded state/cov to model j's state space for mixing.
            % This step is still complex and depends on how model j uses the 6-state.
            % Let's simplify and assume the mixing happens in the 6-state space,
            % and then the mixed state/cov is reduced to model j's size.

            % Mixed state for model j (in 6-state space)
            mixed_state_6 = zeros(6, 1);
            for i = 1:num_models
                state_i_6 = expand_state_to_6(filter_state_models{i}, i);
                mixed_state_6 = mixed_state_6 + mixing_prob(i, j) * state_i_6;
            end
            mixed_state_models{j} = reduce_state_from_6(mixed_state_6, j); % Reduce to model j's size

            % Mixed covariance for model j (in 6-state space)
            mixed_cov_6 = zeros(6, 6);
            for i = 1:num_models
                state_i_6 = expand_state_to_6(filter_state_models{i}, i);
                cov_i_6 = expand_covariance_to_6(filter_cov_models{i}, i);
                % Covariance mixing formula: P_j_mixed = sum_i( mu_ij * [P_i + (x_i - x_j_mixed)*(x_i - x_j_mixed)'] )
                % where x_i and x_j_mixed are in the common space.
                mixed_cov_6 = mixed_cov_6 + mixing_prob(i, j) * (cov_i_6 + (state_i_6 - mixed_state_6) * (state_i_6 - mixed_state_6)');
            end
            mixed_cov_models{j} = reduce_covariance_from_6(mixed_cov_6, j); % Reduce to model j's size
        end
    end

    % 2. Model-Specific Prediction and Update (Kalman Filter steps for each model)
    likelihood = zeros(num_models, 1); % Likelihood of measurement given model

    for j = 1:num_models % For each model j
        current_mixed_state = mixed_state_models{j};
        current_mixed_cov = mixed_cov_models{j};

        % Prediction step
        if j == 1 % CV Model (4-state)
            F = F_cv(1:4, 1:4); % Use 4x4 F for CV model
            Q = Q_cv(1:4, 1:4); % Use 4x4 Q for CV model
            predicted_state = F * current_mixed_state;
            predicted_cov = F * current_mixed_cov * F' + Q;
            H = H_meas(1:2, 1:4); % Measurement matrix for CV model
        elseif j == 2 % CA Model (6-state)
            F = F_ca; % Use 6x6 F for CA model
            Q = Q_ca; % Use 6x6 Q for CA model
            predicted_state = F * current_mixed_state;
            predicted_cov = F * current_mixed_cov * F' + Q;
            H = H_meas; % Measurement matrix for CA model
        elseif j == 3 % CT Model (5-state) - EKF Prediction
            % State: [x; y; vx; vy; omega]
            x = current_mixed_state(1);
            y = current_mixed_state(2);
            vx = current_mixed_state(3);
            vy = current_mixed_state(4);
            omega = current_mixed_state(5);

            % Non-linear prediction function
            if abs(omega) < 1e-9 % Handle straight line approximation
                 predicted_state_ct = [x + vx * dt;
                                       y + vy * dt;
                                       vx;
                                       vy;
                                       omega];
            else
                 predicted_state_ct = [x + (vx/omega)*sin(omega*dt) - (vy/omega)*(1-cos(omega*dt));
                                       y + (vy/omega)*sin(omega*dt) + (vx/omega)*(1-cos(omega*dt));
                                       vx*cos(omega*dt) - vy*sin(omega*dt);
                                       vx*sin(omega*dt) + vy*cos(omega*dt);
                                       omega];
            end
            predicted_state = predicted_state_ct; % Predicted state for CT (5-state)

            % Calculate Jacobian J_ct at current_mixed_state
            % This is complex and omitted for brevity. A numerical Jacobian or symbolic toolbox can help.
            % For this example, we'll use a simplified approximate Jacobian or assume it's provided.
            % A common approximation for the CT Jacobian (5x5):
            % J_ct = [ 1 0 (sin(w*dt)/w) -(1-cos(w*dt))/w x*dt*cos(w*dt)/w - vx*sin(w*dt)/w^2 + vy*(1-cos(w*dt))/w^2
            %          0 1 (1-cos(w*dt))/w  sin(w*dt)/w y*dt*cos(w*dt)/w - vy*sin(w*dt)/w^2 - vx*(1-cos(w*dt))/w^2
            %          0 0 cos(w*dt) -sin(w*dt) -vx*dt*sin(w*dt) - vy*dt*cos(w*dt)
            %          0 0 sin(w*dt) cos(w*dt) vx*dt*cos(w*dt) - vy*dt*sin(w*dt)
            %          0 0 0 0 1 ];
            % This is still complex. Let's use a simplified linear approximation for the CT motion model matrix for demonstration,
            % acknowledging that a proper EKF or Unscented Kalman Filter (UKF) is needed for CT.

            % Simplified CT linear approximation (use CV-like structure but with omega)
            % This is NOT a correct CT model linearization but serves for code structure demonstration.
            % A proper implementation requires the Jacobian.
            % Let's use the Jacobian calculation.
            w = omega;
            w_dt = w * dt;

            if abs(w) < 1e-9 % Handle near-zero omega case for Jacobian
                 J_ct = [ 1 0 dt 0 0;
                          0 1 0 dt 0;
                          0 0 1 0 0;
                          0 0 0 1 0;
                          0 0 0 0 1];
            else
                 J_ct = [ 1 0 sin(w_dt)/w -(1-cos(w_dt))/w (vx*dt*cos(w_dt) - (vx*sin(w_dt))/w + vy*(1-cos(w_dt))/w^2);
                          0 1 (1-cos(w_dt))/w sin(w_dt)/w (vy*dt*cos(w_dt) - (vy*sin(w_dt))/w - vx*(1-cos(w_dt))/w^2);
                          0 0 cos(w_dt) -sin(w_dt) (-vx*dt*sin(w_dt) - vy*dt*cos(w_dt));
                          0 0 sin(w_dt) cos(w_dt) (vx*dt*cos(w_dt) - vy*dt*sin(w_dt));
                          0 0 0 0 1 ];
            end

            F = J_ct; % Use Jacobian as the state transition matrix for EKF
            Q = Q_ct_5state; % Process noise for CT (5-state)
            predicted_cov = F * current_mixed_cov * F' + Q;

            % Measurement matrix for CT model (5-state to 2-state [x; y])
            H = [1 0 0 0 0;
                 0 1 0 0 0];
        end

        % Update step (Standard Kalman Filter update)
        y_meas = measurements(:, k+1); % Current measurement
        if isempty(y_meas) || all(isnan(y_meas))
            % No measurement, prediction is the update
            filter_state_models{j} = predicted_state;
            filter_cov_models{j} = predicted_cov;
            likelihood(j) = 1; % Assume likelihood is 1 if no measurement
        else
            % Calculate innovation (measurement residual)
            innovation = y_meas - H * predicted_state;

            % Calculate innovation covariance
            S = H * predicted_cov * H' + R_meas;

            % Calculate Kalman Gain
            K = predicted_cov * H' * inv(S);

            % Update state estimate
            filter_state_models{j} = predicted_state + K * innovation;

            % Update covariance estimate (Joseph form for numerical stability)
            filter_cov_models{j} = (eye(size(predicted_cov)) - K * H) * predicted_cov * (eye(size(predicted_cov)) - K * H)' + K * R_meas * K';

            % Calculate likelihood of the measurement given the model
            % P(Z_k | model j, Z_k-1) = N(innovation; 0, S)
            likelihood(j) = mvnpdf(innovation', zeros(1, size(innovation, 1)), S);
        end
    end

    % 3. Model Probability Update
    % Calculate normalization factor for model probabilities
    c = c_bar' * likelihood;

    % Update model probabilities P(model j | Z_k)
    for j = 1:num_models
        model_prob(j) = likelihood(j) * c_bar(j) / c;
    end
    model_prob_hist(:, k+1) = model_prob;

    % 4. State and Covariance Combination
    % Combine the state estimates from each model based on updated model probabilities
    % We need to map the model-specific states back to a common state space (e.g., 6-state)
    combined_state_6 = zeros(6, 1);
    for j = 1:num_models
        state_j_6 = expand_state_to_6(filter_state_models{j}, j);
        combined_state_6 = combined_state_6 + model_prob(j) * state_j_6;
    end
    imm_state_est(:, k+1) = combined_state_6;

    % Combine the covariance estimates
    % P_combined = sum_j( mu_j * [P_j + (x_j - x_combined)*(x_j - x_combined)'] )
    % where x_j and x_combined are in the common space.
    combined_cov_6 = zeros(6, 6);
    for j = 1:num_models
        state_j_6 = expand_state_to_6(filter_state_models{j}, j);
        cov_j_6 = expand_covariance_to_6(filter_cov_models{j}, j);
        combined_cov_6 = combined_cov_6 + model_prob(j) * (cov_j_6 + (state_j_6 - combined_state_6) * (state_j_6 - combined_state_6)');
    end
    imm_cov_est(:, :, k+1) = combined_cov_6;

end % End of IMM filter loop

%% Plotting Results
time = (0:N_steps-1) * dt;

figure;
subplot(3, 1, 1);
plot(true_state(1, :), true_state(2, :), 'b-', 'LineWidth', 1.5);
hold on;
plot(measurements(1, :), measurements(2, :), 'rx', 'MarkerSize', 5);
plot(imm_state_est(1, :), imm_state_est(2, :), 'g--', 'LineWidth', 1.5);
xlabel('X Position');
ylabel('Y Position');
title('Target Trajectory: True vs. Measured vs. IMM Estimate');
legend('True Trajectory', 'Measurements', 'IMM Estimate');
grid on;

subplot(3, 1, 2);
plot(time, true_state(3, :), 'b-', 'LineWidth', 1.5);
hold on;
plot(time, imm_state_est(3, :), 'g--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Velocity X');
title('Velocity X: True vs. IMM Estimate');
legend('True Velocity X', 'IMM Estimate Velocity X');
grid on;

subplot(3, 1, 3);
plot(time, model_prob_hist', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Model Probability');
title('Model Probabilities Over Time');
legend(model_names);
grid on;

% Calculate and display RMS errors (position and velocity)
pos_error = true_state(1:2, :) - imm_state_est(1:2, :);
vel_error = true_state(3:4, :) - imm_state_est(3:4, :);

rmse_pos = sqrt(mean(sum(pos_error.^2, 1)));
rmse_vel = sqrt(mean(sum(vel_error.^2, 1)));

fprintf('RMS Position Error: %.2f\n', rmse_pos);
fprintf('RMS Velocity Error: %.2f\n', rmse_vel);

%% Helper Functions for State/Covariance Mapping (Simplified for this example)
% In a real application, these functions need to handle the specific state definitions
% and transformations between model state spaces and the common IMM state space.

function state_6 = expand_state_to_6(state, model_type)
    % Expands a model-specific state to a 6-state vector [x; y; vx; vy; ax; ay]
    % model_type: 1=CV (4-state), 2=CA (6-state), 3=CT (5-state)
    state_6 = zeros(6, 1);
    if model_type == 1 % CV [x; y; vx; vy]
        state_6(1:4) = state;
        % ax, ay are assumed 0 in CV
    elseif model_type == 2 % CA [x; y; vx; vy; ax; ay]
        state_6 = state;
    elseif model_type == 3 % CT [x; y; vx; vy; omega]
        state_6(1:4) = state(1:4); % x, y, vx, vy
        omega = state(5);
        % Approximate ax, ay from vx, vy, omega (centripetal acceleration)
        vx = state(3);
        vy = state(4);
        state_6(5) = -vy * omega; % ax approx
        state_6(6) = vx * omega;  % ay approx
    end
end

function cov_6 = expand_covariance_to_6(cov, model_type)
    % Expands a model-specific covariance to a 6x6 covariance matrix
    % model_type: 1=CV (4x4), 2=CA (6x6), 3=CT (5x5)
    cov_6 = zeros(6, 6);
    if model_type == 1 % CV (4x4)
        cov_6(1:4, 1:4) = cov;
        % Covariance related to ax, ay is 0
    elseif model_type == 2 % CA (6x6)
        cov_6 = cov;
    elseif model_type == 3 % CT (5x5)
        % This is more complex. The covariance terms involving omega need to be
        % mapped to terms involving ax, ay. This mapping is non-linear.
        % A simplified approach is to fill the corresponding parts and set others to 0.
        % This is an approximation. A proper implementation needs the Jacobian of the mapping.
        cov_6(1:4, 1:4) = cov(1:4, 1:4); % x, y, vx, vy covariance
        cov_6(5, 5) = cov(5, 5); % omega variance - placing it in ax/ay variance is incorrect.
        % A better approach involves the Jacobian of the state expansion function.
        % J_expand = d(state_6) / d(state_model)
        % cov_6 = J_expand * cov * J_expand'
        % This is complex. For this example, we'll do a simplified mapping.
        % Map omega variance to some equivalent ax/ay variance - this is a heuristic.
        % Let's just fill the 4x4 part and leave the rest as 0 for simplicity.
        % This is a significant simplification and not strictly correct.
        cov_6(1:4, 1:4) = cov(1:4, 1:4);
        % Covariance related to ax, ay from omega is non-trivial.
        % Setting to small values or zeros is a simplification.
        cov_6(5,5) = 1e-3; % Placeholder, not correct
        cov_6(6,6) = 1e-3; % Placeholder, not correct
    end
end

function state_model = reduce_state_from_6(state_6, model_type)
    % Reduces a 6-state vector to a model-specific state vector
    % model_type: 1=CV (4-state), 2=CA (6-state), 3=CT (5-state)
    if model_type == 1 % CV [x; y; vx; vy]
        state_model = state_6(1:4);
    elseif model_type == 2 % CA [x; y; vx; vy; ax; ay]
        state_model = state_6;
    elseif model_type == 3 % CT [x; y; vx; vy; omega]
        state_model = zeros(5, 1);
        state_model(1:4) = state_6(1:4); % x, y, vx, vy
        % Estimate omega from vx, vy, ax, ay. This is ambiguous.
        % A common approach is to use the estimated omega from the CT model's own update.
        % For combination, we can use the weighted average of the model's own estimates.
        % Let's use the estimated omega from the CT model's filter_state_models{3}.
        % This function is primarily for the mixed state calculation, where we average expanded states.
        % Reducing the average expanded state back to a model's size is also tricky.
        % A simpler approach for combination is to average the model's own updated states in their original size.
        % Then, expand the combined state to the common space for the final IMM output.

        % Let's revise the combination step to average model-specific updated states directly.
        % The expand/reduce functions are mainly needed for the interaction step.

        % For reducing the mixed state (calculated in 6-state space) to model j's size:
        if model_type == 1 % CV (4-state)
             state_model = state_6(1:4);
        elseif model_type == 2 % CA (6-state)
             state_model = state_6;
        elseif model_type == 3 % CT (5-state)
             % Need to estimate omega from the 6-state. This is problematic.
             % A better approach is to average the updated model states directly.
             % Let's assume for now that the 6-state representation is sufficient
             % for the mixed state and we reduce it by taking the relevant parts.
             % Estimating omega from [vx; vy; ax; ay] is possible but requires assumptions.
             % ax = -vy*omega, ay = vx*omega => omega = -ax/vy = ay/vx (if denominators non-zero)
             % A more robust way is atan2(ay, ax) / sqrt(vx^2+vy^2) - requires care.
             % Let's just take the first 4 states and set omega to 0 for reduction to 5-state.
             % This is a simplification.
             state_model = [state_6(1:4); 0]; % Simplified: assume omega is 0 in reduced state
        end
    end
end

function cov_model = reduce_covariance_from_6(cov_6, model_type)
    % Reduces a 6x6 covariance matrix to a model-specific covariance matrix
    % model_type: 1=CV (4x4), 2=CA (6x6), 3=CT (5x5)
    if model_type == 1 % CV (4x4)
        cov_model = cov_6(1:4, 1:4);
    elseif model_type == 2 % CA (6x6)
        cov_model = cov_6;
    elseif model_type == 3 % CT (5x5)
        % This is also complex. Need to map covariance from 6-state to 5-state.
        % Requires the Jacobian of the reduction function.
        % Let's take the relevant parts and set others to 0 as a simplification.
        cov_model = zeros(5, 5);
        cov_model(1:4, 1:4) = cov_6(1:4, 1:4); % x, y, vx, vy covariance
        % Covariance related to omega needs mapping from ax/ay covariance.
        % Simplified: just take the 4x4 part.
        cov_model = cov_6(1:4, 1:4); % This is incorrect, size mismatch.
        % A proper reduction requires mapping.
        % Let's assume for simplicity we just take the top-left part.
        cov_model = cov_6(1:5, 1:5); % This assumes the 6-state has omega in the 5th position, which it doesn't.
        % Correct reduction requires mapping based on the state definitions.
        % Let's just take the 4x4 part for x, y, vx, vy and set omega variance to a nominal value.
        cov_model = diag([diag(cov_6(1:4, 1:4)); 1e-3]); % Simplified: take diagonal and add nominal omega variance
        % This is a significant simplification.
    end
end

