% createCompositeTrajectory: Generates a composite trajectory from a sequence of motion segments and saves data.

% Define the time step
dt = 0.01; % seconds

% Define the sequence of motion segments
% Each segment is a struct with 'type', 'duration', and specific parameters.
% The first segment MUST define 'start_pos'.
% Subsequent segments' initial conditions are derived from the end of the previous.
segments = {
    struct('type', 'CV', 'duration', 10, 'velocity', [10; 0], 'start_pos', [0; 0]); % Start at (0,0), move right at speed 1 (defined by velocity vector)
    struct('type', 'CTR', 'duration', 10, 'turnRate', deg2rad(30)); % Turn 30 deg/s for 3s (total 90 deg turn), maintain speed 1
    %struct('type', 'CV', 'duration', 30, 'speed', 10)
    struct('type', 'CA', 'duration', 4, 'acceleration', [5; 5]); 
    %struct('type', 'CA', 'duration', 10, 'acceleration', [3; 3]);
    %struct('type', 'CTR', 'duration', 10, 'turnRate', deg2rad(-18));
    % struct('type', 'SINUSOIDAL', 'duration', 6, 'amplitude_vx', 0.2, 'amplitude_vy', 0.5, 'frequency', 0.5); % Sinusoidal velocity variation
    % struct('type', 'CONSTANT_JERK', 'duration', 3, 'jerk', [0.1; -0.2]); % Constant jerk segment
    %struct('type', 'CV', 'duration', 30, 'speed', 30); % Move at new constant speed 1.5, using current heading
    %struct('type', 'CTR', 'duration', 10, 'turnRate', deg2rad(-2), 'speed', []); % Turn -20 deg/s for 4s (total -80 deg turn), maintain speed (derived from previous segment)
};

% Initialize arrays to store data
all_positions = [];
all_velocities = [];
all_accelerations = [];
all_times = [];

% Initialize current state variables for segment transitions
current_pos = [];
current_velocity = [];
current_acceleration = []; % Added for jerk segment
current_heading = []; % Radians
current_speed = [];
current_time = 0;

% Loop through each segment
for i = 1:length(segments)
    segment = segments{i};
    segment_type = segment.type;
    duration = segment.duration;

    % Set initial conditions for the current segment based on the end of the previous
    if i == 1
        if ~isfield(segment, 'start_pos')
            error('First segment must specify start_pos.');
        end
        current_pos = segment.start_pos;

        % Initial velocity/heading/speed/acceleration for the first segment based on its type
        if strcmp(segment_type, 'CV')
            % For the first CV segment, velocity vector OR speed + heading must be provided
            if isfield(segment, 'velocity')
                 current_velocity = segment.velocity;
                 current_speed = norm(current_velocity);
                 current_heading = atan2(current_velocity(2), current_velocity(1)); % Calculate initial heading
            elseif isfield(segment, 'speed') && isfield(segment, 'start_heading')
                 current_speed = segment.speed;
                 current_heading = segment.start_heading;
                 current_velocity = [current_speed * cos(current_heading);
                                     current_speed * sin(current_heading)];
            else
                 error('First CV segment must specify either ''velocity'' or ''speed'' and ''start_heading''.');
            end
             current_acceleration = [0; 0]; % CV starts with zero acceleration

        elseif strcmp(segment_type, 'CA')
             if ~isfield(segment, 'start_velocity')
                 error('First CA segment must specify start_velocity.');
             end
             current_velocity = segment.start_velocity;
             current_speed = norm(current_velocity);
             current_heading = atan2(current_velocity(2), current_velocity(1)); % Calculate initial heading
             if ~isfield(segment, 'acceleration')
                  error('First CA segment must specify acceleration.');
             end
             current_acceleration = segment.acceleration; % CA starts with its defined acceleration


        elseif strcmp(segment_type, 'CTR')
             if ~isfield(segment, 'start_speed') || ~isfield(segment, 'start_heading')
                 error('First CTR segment must specify start_speed and start_heading.');
             end
             current_speed = segment.start_speed;
             current_heading = segment.start_heading;
             current_velocity = [current_speed * cos(current_heading);
                                 current_speed * sin(current_heading)]; % Calculate initial velocity
             % Initial acceleration for CTR is centripetal
             omega = segment.turnRate;
             v = current_speed;
             accel_magnitude = v * abs(omega);
             accel_heading = current_heading + sign(omega) * (pi/2);
             current_acceleration = [accel_magnitude * cos(accel_heading);
                                     accel_magnitude * sin(accel_heading)];

        elseif strcmp(segment_type, 'SINUSOIDAL')
             % For the first sinusoidal segment, start_velocity must be provided
             if ~isfield(segment, 'start_velocity')
                 error('First SINUSOIDAL segment must specify start_velocity.');
             end
             current_velocity = segment.start_velocity;
             current_speed = norm(current_velocity);
             current_heading = atan2(current_velocity(2), current_velocity(1)); % Calculate initial heading
             % Initial acceleration for Sinusoidal (at t=0)
             amplitude_vx = segment.amplitude_vx;
             amplitude_vy = segment.amplitude_vy;
             frequency = segment.frequency;
             omega = 2 * pi * frequency;
             current_acceleration = [amplitude_vx * omega * cos(0); % cos(omega*0) = 1
                                     amplitude_vy * omega * cos(0)];

        elseif strcmp(segment_type, 'CONSTANT_JERK')
            % For the first constant jerk segment, initial velocity and acceleration must be provided
            if ~isfield(segment, 'start_velocity') || ~isfield(segment, 'start_acceleration')
                error('First CONSTANT_JERK segment must specify start_velocity and start_acceleration.');
            end
            current_velocity = segment.start_velocity;
            current_acceleration = segment.start_acceleration;
            current_speed = norm(current_velocity);
            current_heading = atan2(current_velocity(2), current_velocity(1)); % Calculate initial heading

        else
            error('Unknown segment type in the first segment.');
        end
    else
        % For subsequent segments, initial conditions are the end conditions of the previous segment
        % current_pos, current_velocity, current_acceleration, current_speed, current_heading are already set from the end of the previous segment's loop
    end

    % Generate points and full state data for the current segment using the updated functions
    if strcmp(segment_type, 'CV')
        % Check if speed is specified for this CV segment, otherwise use velocity vector
        if isfield(segment, 'speed') && ~isempty(segment.speed)
            segment_speed = segment.speed;
            % Calculate velocity vector using the specified speed and current heading
            segment_velocity = [segment_speed * cos(current_heading);
                                segment_speed * sin(current_heading)];
        elseif isfield(segment, 'velocity')
             segment_velocity = segment.velocity; % Use velocity vector if provided
             segment_speed = norm(segment_velocity); % Update speed based on this velocity
        else
             error('CV segment must specify either ''velocity'' or ''speed''.');
        end

        [segment_points, segment_velocities, segment_accelerations] = generateConstantVelocitySegment(current_pos, segment_velocity, duration, dt);

        % Update state for the next segment (end of this segment)
        current_pos = segment_points(end, :)';
        current_velocity = segment_velocities(end, :)'; % End velocity is the segment velocity
        current_speed = norm(current_velocity);
        current_heading = atan2(current_velocity(2), current_velocity(1)); % Update heading
        current_acceleration = segment_accelerations(end, :)'; % Acceleration is zero at the end of CV


    elseif strcmp(segment_type, 'CA')
        segment_acceleration = segment.acceleration; % Use acceleration for this segment
        % Use current_velocity as the start_velocity for this CA segment
        [segment_points, segment_velocities, segment_accelerations, end_velocity_ca] = generateConstantAccelerationSegment(current_pos, current_velocity, segment_acceleration, duration, dt);

        % Update state for the next segment (end of this segment)
        current_pos = segment_points(end, :)';
        current_velocity = end_velocity_ca; % End velocity is calculated by the function
        current_speed = norm(current_velocity);
        current_heading = atan2(current_velocity(2), current_velocity(1)); % Update heading
        current_acceleration = segment_accelerations(end, :)'; % Acceleration is constant at the end of CA


    elseif strcmp(segment_type, 'CTR')
         % If speed is specified for this CTR segment, use it, otherwise use current_speed
        if isfield(segment, 'speed') && ~isempty(segment.speed)
            segment_speed = segment.speed;
        else
            segment_speed = current_speed; % Maintain speed from previous segment
        end
        segment_turn_rate = segment.turnRate; % Use turn rate for this segment

        % Use current_pos, current_speed, current_heading as start conditions
        [segment_points, segment_velocities, segment_accelerations, end_pos_ctr, end_heading_ctr] = generateConstantTurnRateSegment(current_pos, segment_speed, current_heading, segment_turn_rate, duration, dt);

        % Update state for the next segment (end of this segment)
        current_pos = end_pos_ctr; % End position from the function
        current_heading = end_heading_ctr; % End heading from the function
        current_speed = segment_speed; % Speed is constant during CTR
        current_velocity = [current_speed * cos(current_heading);
                            current_speed * sin(current_heading)]; % Update velocity vector
        current_acceleration = segment_accelerations(end, :)'; % Acceleration is centripetal at the end of CTR


    elseif strcmp(segment_type, 'SINUSOIDAL')
        % Use current_pos and current_velocity as start conditions
        amplitude_vx = segment.amplitude_vx;
        amplitude_vy = segment.amplitude_vy;
        frequency = segment.frequency;

        [segment_points, segment_velocities, segment_accelerations, end_velocity_sin] = generateSinusoidalSegment(current_pos, current_velocity, amplitude_vx, amplitude_vy, frequency, duration, dt);

        % Update state for the next segment (end of this segment)
        current_pos = segment_points(end, :)';
        current_velocity = end_velocity_sin; % End velocity is calculated by the function
        current_speed = norm(current_velocity);
        current_heading = atan2(current_velocity(2), current_velocity(1)); % Update heading
        current_acceleration = segment_accelerations(end, :)'; % Acceleration at the end of Sinusoidal


    elseif strcmp(segment_type, 'CONSTANT_JERK')
        segment_jerk = segment.jerk; % Use jerk for this segment
        % Use current_pos, current_velocity, current_acceleration as start conditions
        [segment_points, segment_velocities, segment_accelerations, end_velocity_jerk, end_acceleration_jerk] = generateConstantJerkSegment(current_pos, current_velocity, current_acceleration, segment_jerk, duration, dt);

        % Update state for the next segment (end of this segment)
        current_pos = segment_points(end, :)';
        current_velocity = end_velocity_jerk; % End velocity from the function
        current_acceleration = end_acceleration_jerk; % End acceleration from the function
        current_speed = norm(current_velocity);
        current_heading = atan2(current_velocity(2), current_velocity(1)); % Update heading


    else
        error('Unknown segment type.');
    end

    % Append the generated segment data to the full trajectory arrays
    % For the first segment, include the initial point (index 1).
    % For subsequent segments, exclude the first point (index 1) to avoid duplicates,
    % as it's the same as the last point of the previous segment.
    if i == 1
        all_positions = [all_positions; segment_points];
        all_velocities = [all_velocities; segment_velocities];
        all_accelerations = [all_accelerations; segment_accelerations];
        % Generate time vector for this segment and append
        segment_times = current_time + (0:size(segment_points, 1)-1)' * dt;
        all_times = [all_times; segment_times];
    else
        all_positions = [all_positions; segment_points(2:end, :)];
        all_velocities = [all_velocities; segment_velocities(2:end, :)];
        all_accelerations = [all_accelerations; segment_accelerations(2:end, :)];
         % Generate time vector for this segment (excluding the first point's time) and append
        segment_times = current_time + (1:size(segment_points, 1)-1)' * dt;
        all_times = [all_times; segment_times];
    end

    % Update current_time for the next segment
    current_time = all_times(end);


end % End of segments loop

filename = 'C:\Users\acr0093\Documents\UCAH\Optimal-Final-Project\src\trajectory\trajectoryDataCVCTCA';
traj = [all_positions,all_velocities,all_accelerations];
save(filename, 'traj');

% Plot the trajectory
figure;
plot(all_positions(:, 1), all_positions(:, 2), '-o');
xlabel('X Position');
ylabel('Y Position');
grid on;
axis equal; % Keep x and y scales the same

% You can now load the data in another script using:
% load('trajectory_data.mat');
% This will load variables all_positions, all_velocities, all_accelerations, all_times
