function [trajectory_segment, velocities, accelerations] = generateConstantVelocitySegment(start_pos, velocity, duration, dt)
% generateConstantVelocitySegment: Generates points, velocities, and accelerations for a constant velocity trajectory segment.
%
%   [trajectory_segment, velocities, accelerations] = generateConstantVelocitySegment(start_pos, velocity, duration, dt)
%
%   Inputs:
%       start_pos (2x1 vector): Initial position [x0; y0].
%       velocity (2x1 vector): Constant velocity [vx; vy].
%       duration (scalar): Duration of the segment.
%       dt (scalar): Time step.
%
%   Outputs:
%       trajectory_segment (Nx2 matrix): A matrix where each row is a point [x, y].
%       velocities (Nx2 matrix): A matrix where each row is a velocity vector [vx, vy].
%       accelerations (Nx2 matrix): A matrix where each row is an acceleration vector [ax, ay].

    num_steps = floor(duration / dt);
    trajectory_segment = zeros(num_steps + 1, 2);
    velocities = zeros(num_steps + 1, 2);
    accelerations = zeros(num_steps + 1, 2);

    x0 = start_pos(1);
    y0 = start_pos(2);
    vx = velocity(1);
    vy = velocity(2);

    % Velocity and acceleration are constant for this segment
    constant_velocity = [vx, vy];
    constant_acceleration = [0, 0];

    for i = 0:num_steps
        t = i * dt;
        current_x = x0 + vx * t;
        current_y = y0 + vy * t;

        trajectory_segment(i + 1, :) = [current_x, current_y];
        velocities(i + 1, :) = constant_velocity;
        accelerations(i + 1, :) = constant_acceleration;
    end

end
