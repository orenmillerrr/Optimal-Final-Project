% createCompositeTrajectory: Generates a composite trajectory from a sequence of motion segments and saves data.
clear;clc;close all
% Define the time step
dt = 0.01; % seconds

% Define the sequence of motion segments
% Each segment is a struct with 'type', 'duration', and specific parameters.
% The first segment MUST define 'start_pos'.
% Subsequent segments' initial conditions are derived from the end of the previous.
segments = {
    struct('type', 'CV', 'duration', 10, 'velocity', [100; 0], 'start_pos', [0; 0]); % Start at (0,0), move right at speed 1
    struct('type', 'CTR', 'duration', 5, 'turnRate', deg2rad(90/5)); % Turn 30 deg/s for 3s (total 90 deg turn), maintain speed 1
    struct('type', 'CV', 'duration', 10, 'velocity', [0; 100]); % Move at new constant velocity
    struct('type', 'CTR', 'duration', 5, 'turnRate', deg2rad(90/5)); % Turn 30 deg/s for 3s (total 90 deg turn), maintain speed 1
    struct('type', 'CV', 'duration', 10, 'velocity', [-100; 0]); % Move at new constant velocity
    % struct('type', 'CA', 'duration', 10, 'acceleration', [0; -9.8]); % Accelerate upwards
    % struct('type', 'CV', 'duration', 5, 'velocity', [0.5; 9.8]); % Move at new constant velocity
    struct('type', 'CTR', 'duration',25, 'turnRate', deg2rad(90/5), 'speed', []); % Turn -20 deg/s for 4s (total -80 deg turn), maintain speed (derived from previous segment)
};



% Initialize arrays to store data
all_positions = [];
all_velocities = [];
all_accelerations = [];
all_times = [];

% Initialize current state variables for segment transitions
current_pos = [];
current_velocity = [];
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

        % Initial velocity/heading/speed for the first segment based on its type
        if strcmp(segment_type, 'CV')
            if ~isfield(segment, 'velocity')
                 error('First CV segment must specify velocity.');
            end
            current_velocity = segment.velocity;
            current_speed = norm(current_velocity);
            current_heading = atan2(current_velocity(2), current_velocity(1)); % Calculate initial heading

        elseif strcmp(segment_type, 'CA')
             if ~isfield(segment, 'start_velocity')
                 error('First CA segment must specify start_velocity.');
             end
             current_velocity = segment.start_velocity;
             current_speed = norm(current_velocity);
             current_heading = atan2(current_velocity(2), current_velocity(1)); % Calculate initial heading

        elseif strcmp(segment_type, 'CTR')
             if ~isfield(segment, 'start_speed') || ~isfield(segment, 'start_heading')
                 error('First CTR segment must specify start_speed and start_heading.');
             end
             current_speed = segment.start_speed;
             current_heading = segment.start_heading;
             current_velocity = [current_speed * cos(current_heading);
                                 current_speed * sin(current_heading)]; % Calculate initial velocity
        else
            error('Unknown segment type in the first segment.');
        end
    else
        % For subsequent segments, initial conditions are the end conditions of the previous segment
        % current_pos, current_velocity, current_speed, current_heading are already set from the end of the previous segment's loop
    end

    % Generate points and full state data for the current segment using the updated functions
    if strcmp(segment_type, 'CV')
        segment_velocity = segment.velocity; % Use velocity for this segment
        [segment_points, segment_velocities, segment_accelerations] = generateConstantVelocitySegment(current_pos, segment_velocity, duration, dt);

        % Update state for the next segment (end of this segment)
        current_pos = segment_points(end, :)';
        current_velocity = segment_velocities(end, :)'; % End velocity is the segment velocity
        current_speed = norm(current_velocity);
        current_heading = atan2(current_velocity(2), current_velocity(1)); % Update heading


    elseif strcmp(segment_type, 'CA')
        segment_acceleration = segment.acceleration; % Use acceleration for this segment
        % Use current_velocity as the start_velocity for this CA segment
        [segment_points, segment_velocities, segment_accelerations, end_velocity_ca] = generateConstantAccelerationSegment(current_pos, current_velocity, segment_acceleration, duration, dt);

        % Update state for the next segment (end of this segment)
        current_pos = segment_points(end, :)';
        current_velocity = end_velocity_ca; % End velocity is calculated by the function
        current_speed = norm(current_velocity);
        current_heading = atan2(current_velocity(2), current_velocity(1)); % Update heading


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

% Save the data to a .mat file
filename = 'D:\Classes\MECH 7710 - Optimal Control\Optimal-Final-Project';
traj = [all_positions(:,1),all_velocities(:,1),all_accelerations(:,1),all_positions(:,2),all_velocities(:,2),all_accelerations(:,2)];
save(filename, 'traj');

fprintf('Trajectory data saved to %s\n', filename);

% Plot the trajectory
figure;
plot(all_positions(:, 1), all_positions(:, 2), '-o');
title('Composite Trajectory');
xlabel('X Position');
ylabel('Y Position');
grid on;
axis equal; % Keep x and y scales the same

% You can now load the data in another script using:
% load('trajectory_data.mat');
% This will load variables all_positions, all_velocities, all_accelerations, all_times
