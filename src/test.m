% IMM_FILTER Implements an Interacting Multiple Model (IMM) filter.
%
% This is a basic structure for an IMM filter. You need to define
% your specific models (state transition, process noise, measurement,
% measurement noise matrices) and the model transition matrix.
%
% Inputs:
%   models: A cell array of structures, where each structure defines a model:
%           models{i}.A - State transition matrix (size state_dim_i x state_dim_i)
%           models{i}.H - Measurement matrix (size meas_dim x state_dim_i)
%           models{i}.Q - Process noise covariance (size state_dim_i x state_dim_i)
%           models{i}.R - Measurement noise covariance (size meas_dim x meas_dim)
%           models{i}.B - Control input matrix (optional, size state_dim_i x control_dim, set to [] if not used)
%   P_trans: Model transition probability matrix (size num_models x num_models, P_trans(i,j) = P(model j | model i at previous step))
%   x_init: Initial state estimate (a cell array of size num_models, x_init{i} is size state_dim_i x 1)
%   P_init: Initial state covariance (a cell array of size num_models, P_init{i} is size state_dim_i x state_dim_i)
%   mu_init: Initial model probabilities (a row vector of size 1 x num_models, sums to 1)
%   measurements: A matrix of measurements, where each row is a measurement vector (size num_measurements x meas_dim)
%   dt: Time step (scalar, used if models are time-varying, or for calculating A/Q/R)
%   u: Control input (optional, a matrix where each row is the control vector, size num_measurements x control_dim)
%
% Outputs:
%   x_est: Combined state estimate history (size num_measurements x state_dim_max)
%   P_est: Combined state covariance history (size state_dim_max x state_dim_max x num_measurements)
%   mu_est: Model probability history (size num_measurements x num_models)
%   x_model_est: State estimate history for each individual model (size num_measurements x state_dim_max x num_models)
%   P_model_est: State covariance history for each individual model (size state_dim_max x state_dim_max x num_measurements x num_models)
%
% Note: This implementation assumes all models have the same measurement dimension (meas_dim).
% It also assumes a common state dimension for simplicity in output storage (state_dim).
% If models have different state dimensions, you would need to adjust the storage
% and combination steps accordingly, potentially padding smaller state vectors/covariances
% or using cell arrays for outputs. The current code uses the state dimension of the first model.

function [x_est, P_est, mu_est, x_model_est, P_model_est] = imm_filter(models, P_trans, x_init, P_init, mu_init, measurements, dt, u)

    % --- Initialization ---
    num_models = length(models);
    num_measurements = size(measurements, 1);
    % Assuming all models have the same measurement dimension
    meas_dim = size(measurements, 2);
    % Assuming a common state dimension for output storage - using the first model's dimension
    state_dim = size(x_init{1}, 1);

    % Initialize storage for results
    x_est = zeros(num_measurements, state_dim);
    P_est = zeros(state_dim, state_dim, num_measurements);
    mu_est = zeros(num_measurements, num_models);
    % Note: If models have different state dimensions, x_model_est and P_model_est
    % storage would need to be adjusted (e.g., using cell arrays).
    x_model_est = zeros(num_measurements, state_dim, num_models);
    P_model_est = zeros(state_dim, state_dim, num_measurements, num_models);


    % Current state and covariance for each model (using cell arrays to handle potentially different state dimensions internally)
    x_current = x_init;
    P_current = P_init;
    mu_current = mu_init;

    % --- IMM Loop ---
    for k = 1:num_measurements
        y_k = measurements(k, :)'; % Current measurement (column vector)

        % --- 1. Mixing (Interaction) ---
        % Calculate predicted model probabilities
        mu_predicted = mu_current * P_trans;

        % Calculate mixing probabilities (probability of model i at k-1 given model j at k)
        % mu_ij = P(model_i(k-1) | model_j(k), y_k)
        % This is calculated using the predicted model probabilities and transition matrix
        % mu_ij = P(model_j(k) | model_i(k-1)) * P(model_i(k-1)) / P(model_j(k))
        % mu_ij = P_trans(i,j) * mu_current(i) / mu_predicted(j)
        mixing_probs = zeros(num_models, num_models);
        for j = 1:num_models % For each filter j (model j at current step)
            if mu_predicted(j) > eps % Use eps for small number comparison to avoid division by near-zero
                for k = 1:num_models % From each filter i (model i at previous step)
                    mixing_probs(k, j) = P_trans(k, j) * mu_current(k) / mu_predicted(j);
                end
            else
                 % If predicted probability is negligible, assume equal mixing from previous models
                 mixing_probs(:, j) = 1/num_models;
            end
            % Ensure mixing probabilities sum to 1 for numerical stability
            mixing_probs(:, j) = mixing_probs(:, j) / sum(mixing_probs(:, j));
        end

        % Calculate mixed initial state and covariance for each filter j
        x_mixed = cell(num_models, 1);
        P_mixed = cell(num_models, 1);
        for j = 1:num_models % For each filter j
            state_dim_j = size(models{j}.A, 1); % Get state dimension for model j
            x_mixed{j} = zeros(state_dim_j, 1);
            % Initialize P_mixed{j} with zeros, it will be built up
            P_mixed{j} = zeros(state_dim_j, state_dim_j);

            for k = 1:num_models % Sum over all previous models i
                 % Need to handle different state dimensions if they exist.
                 % This basic implementation assumes compatibility or requires careful
                 % definition of models such that the mixing makes sense.
                 % A more robust implementation would require state-space conversions
                 % between models if dimensions differ significantly.
                 % For models with different state dimensions, this mixing step
                 % as written might not be directly applicable without state augmentation
                 % or other techniques. Assuming compatible state spaces for now.
                 state_dim_i = size(models{k}.A, 1);
                 if state_dim_i == state_dim_j % Simple case: dimensions match
                    x_mixed{j} = x_mixed{j} + x_current{k} * mixing_probs(k, j);
                    % Calculate covariance contribution from model i
                    cov_i = P_current{k} + (x_current{k} - x_mixed{j}) * (x_current{k} - x_mixed{j})';
                    P_mixed{j} = P_mixed{j} + mixing_probs(k, j) * cov_i;
                 else
                     % --- Handle different state dimensions ---
                     % This is a simplification. A proper IMM with different state
                     % dimensions requires careful consideration of how states
                     % are related and potentially state augmentation or mapping.
                     % For this basic example, we will issue a warning and skip
                     % mixing from incompatible models. In a real application,
                     % you'd implement a specific state conversion or augmentation.
                     warning('IMM:StateDimMismatch', 'State dimensions of models %d and %d do not match for mixing.', k, j);
                     % A more advanced approach would involve mapping x_current{i}
                     % and P_current{i} into the state space of model j before mixing.
                     % For now, we'll just skip this contribution.
                 end
            end
        end

        % --- 2. Model-Conditioned Filtering (Kalman Filter) ---
        x_filtered = cell(num_models, 1);
        P_filtered = cell(num_models, 1);
        likelihood = zeros(1, num_models);

        for j = 1:num_models % For each filter j
            A = models{j}.A;
            H = models{j}.H;
            Q = models{j}.Q;
            R = models{j}.R;
            B = models{j}.B;
            state_dim_j = size(A, 1); % State dimension for model j
            meas_dim_j = size(H, 1); % Measurement dimension for model j (should match meas_dim)

            % Check measurement dimension consistency
            if meas_dim_j ~= meas_dim
                 error('IMM:MeasurementDimMismatch', 'Measurement dimension of model %d (%d) does not match overall measurement dimension (%d).', j, meas_dim_j, meas_dim);
            end

            % Prediction Step
            if ~isempty(B) && size(u,1) >= k
                % Ensure control input dimension matches B
                if size(u, 2) ~= size(B, 2)
                     error('IMM:ControlDimMismatch', 'Control input dimension (%d) does not match control matrix B dimension (%d) for model %d.', size(u, 2), size(B, 2), j);
                end
                x_pred = A * x_mixed{j} + B * u(k,:)';
            else
                 x_pred = A * x_mixed{j};
            end
            P_pred = A * P_mixed{j} * A' + Q;

            % Update Step
            y_pred = H * x_pred;
            innovation = y_k - y_pred;
            S = H * P_pred * H' + R; % Innovation covariance

            % Add a check for singular S matrix before inversion
            if det(S) < eps * max(abs(diag(S))) % Check for near-singularity
                 warning('IMM:SingularInnovationCovariance', 'Innovation covariance matrix S is singular or near-singular for model %d at step %d.', j, k);
                 % In case of singularity, the likelihood calculation and Kalman gain
                 % will fail. A common approach is to skip the update for this model
                 % or add a small regularization term to S. For this basic code,
                 % we'll set likelihood to a very small value and skip the update
                 % for this model in this step.
                 likelihood(j) = eps; % Assign a very small likelihood
                 x_filtered{j} = x_pred; % Keep predicted state
                 P_filtered{j} = P_pred; % Keep predicted covariance
                 continue; % Skip the rest of the update step for this model
            end

            % Calculate likelihood (assuming Gaussian distribution)
            % Likelihood = N(innovation; 0, S)
            likelihood(j) = (1 / sqrt((2 * pi)^meas_dim * det(S))) * exp(-0.5 * innovation' * inv(S) * innovation);

            K = P_pred * H' * inv(S); % Kalman gain
            x_filtered{j} = x_pred + K * innovation;
            % Joseph form for numerical stability: (I - KH)P(I - KH)' + KRK'
            P_filtered{j} = (eye(state_dim_j) - K * H) * P_pred * (eye(state_dim_j) - K * H)' + K * R * K';
            % A simpler form (less numerically stable): P_filtered{j} = (eye(state_dim_j) - K * H) * P_pred;
        end

        % --- 3. Model Probability Update ---
        c_bar = sum(likelihood .* mu_predicted); % Normalization constant

        mu_updated = zeros(1, num_models);
        if c_bar > eps % Use eps for small number comparison to avoid division by near-zero
            for j = 1:num_models
                mu_updated(j) = (likelihood(j) * mu_predicted(j)) / c_bar;
            end
        else
            % If normalization constant is negligible, assign equal probability
            % This can happen if all models have very low likelihoods for the measurement
            warning('IMM:NormalizationConstantZero', 'Normalization constant c_bar is zero or near-zero at step %d. Assigning equal model probabilities.', k);
            mu_updated = ones(1, num_models) / num_models;
        end
         % Ensure updated probabilities sum to 1 due to potential floating point issues
        mu_updated = mu_updated / sum(mu_updated);


        % --- 4. Combination ---
        % Note: This combination step assumes all models have the same state dimension
        % as 'state_dim'. If not, this needs significant modification.
        x_combined = zeros(state_dim, 1);
        P_combined = zeros(state_dim, state_dim);

        for j = 1:num_models
             state_dim_j = size(models{j}.A, 1);
             if state_dim_j == state_dim % Only combine if dimensions match the assumed output dimension
                x_combined = x_combined + mu_updated(j) * x_filtered{j};
             else
                 % If state dimensions differ, this combination is not straightforward
                 % and requires mapping or augmentation.
                 warning('IMM:CombinationDimMismatch', 'State dimension of model %d (%d) does not match the assumed output state dimension (%d) for combination.', j, state_dim_j, state_dim);
                 % Skipping contribution from this model to the combined state/covariance
             end
        end

        for j = 1:num_models
             state_dim_j = size(models{j}.A, 1);
             if state_dim_j == state_dim % Only combine if dimensions match the assumed output dimension
                % Calculate covariance contribution from model j to the combined state
                cov_j = P_filtered{j} + (x_filtered{j} - x_combined) * (x_filtered{j} - x_combined)';
                P_combined = P_combined + mu_updated(j) * cov_j;
             end
        end

        % --- Store Results and Update for Next Iteration ---
        x_est(k, :) = x_combined';
        P_est(:, :, k) = P_combined;
        mu_est(k, :) = mu_updated;
        for j = 1:num_models
            % Store individual model estimates - requires careful handling if state dimensions differ
            state_dim_j = size(models{j}.A, 1);
             if state_dim_j == state_dim % Only store if dimensions match the assumed output dimension
                x_model_est(k, :, j) = x_filtered{j}';
                P_model_est(:, :, k, j) = P_filtered{j};
             else
                 % If state dimensions differ, storing in a fixed-size array is not possible.
                 % You might need to store these in a cell array instead of a fixed matrix.
                 warning('IMM:StoringDimMismatch', 'State dimension of model %d (%d) does not match the assumed output state dimension (%d) for storing individual model estimates.', j, state_dim_j, state_dim);
                 % You would need alternative storage here, e.g.,
                 % x_model_est_cell{k, j} = x_filtered{j};
                 % P_model_est_cell{k, j} = P_filtered{j};
             end
        end


        x_current = x_filtered; % Update current state estimates for next iteration
        P_current = P_filtered; % Update current covariances for next iteration
        mu_current = mu_updated; % Update current model probabilities for next iteration
    end

end

% --- Example Usage (You would define your own models and data) ---
% Example: Two models - Constant Velocity (CV) and Constant Acceleration (CA)
% Note: This example has models with different state dimensions (CV=2, CA=3).
% The provided IMM function has limitations in handling different state dimensions
% in the mixing, combination, and output storage steps as noted in the comments.
% For this example to run without errors related to dimension mismatches in the
% current function implementation, you would need to either:
% 1. Use models with the same state dimension.
% 2. Modify the IMM function significantly to handle state-space conversions
%    or augmentation for mixing and combination.
% The example below will likely trigger the dimension mismatch warnings/errors
% in the current function if run directly.

dt = 0.1; % Time step

% CV Model (State: [position; velocity])
A_cv = [1 dt; 0 1];
H_cv = [1 0]; % Measure position
Q_cv = [0.1*dt^3/3 0.1*dt^2/2; 0.1*dt^2/2 0.1*dt]; % Example process noise
R_cv = 0.5; % Example measurement noise
model_cv = struct('A', A_cv, 'H', H_cv, 'Q', Q_cv, 'R', R_cv, 'B', []);

% CA Model (State: [position; velocity; acceleration])
A_ca = [1 dt 0.5*dt^2; 0 1 dt; 0 0 1];
H_ca = [1 0 0]; % Measure position
Q_ca = [0.5*dt^5/20 0.5*dt^4/8 0.5*dt^3/6; ...
        0.5*dt^4/8 0.5*dt^3/3 0.5*dt^2/2; ...
        0.5*dt^3/6 0.5*dt^2/2 0.5*dt]; % Example process noise
R_ca = 0.5; % Example measurement noise (same as CV for simplicity)
model_ca = struct('A', A_ca, 'H', H_ca, 'Q', Q_ca, 'R', R_ca, 'B', []);

models = {model_cv, model_ca}; % Cell array of models

% Model Transition Matrix (Example: 95% chance of staying in the same model)
P_trans = [0.95 0.05;
           0.05 0.95];

% Initial Conditions
x_init_cv = [0; 0]; % Initial state for CV [position; velocity]
P_init_cv = eye(2); % Initial covariance for CV

x_init_ca = [0; 0; 0]; % Initial state for CA [position; velocity; acceleration]
P_init_ca = eye(3); % Initial covariance for CA

% Initial states and covariances for each model in a cell array
x_init = {x_init_cv, x_init_ca};
P_init = {P_init_cv, P_init_ca};

mu_init = [0.5 0.5]; % Initial model probabilities (sum to 1)

% Generate some example measurements (e.g., noisy position data)
t = 0:dt:10;
true_pos = 0.5 * 1 * t.^2; % Example: Constant acceleration of 1
measurements = true_pos' + 0.5 * randn(length(t), 1); % Add noise

% Run the IMM filter
% Note: Running this with the current function and the CV/CA models
% will likely trigger dimension mismatch warnings/errors.
[x_est, P_est, mu_est, x_model_est, P_model_est] = imm_filter(models, P_trans, x_init, P_init, mu_init, measurements, dt, []);

% Plot results (example: position estimate and model probabilities)
figure;
subplot(2,1,1);
plot(t, measurements, 'o', t, x_est(:, 1), '-');
legend('Measurements', 'IMM Position Estimate');
xlabel('Time');
ylabel('Position');
title('IMM Filter Position Estimate');

subplot(2,1,2);
plot(t, mu_est);
legend('Model 1 (CV)', 'Model 2 (CA)');
xlabel('Time');
ylabel('Model Probability');
title('IMM Model Probabilities');
