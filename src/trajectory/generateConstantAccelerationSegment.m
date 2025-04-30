function [trajectory_segment, velocities, accelerations, end_velocity] = generateConstantAccelerationSegment(start_pos, start_velocity, acceleration, duration, dt)
% generateConstantAccelerationSegment: Generates points, velocities, and accelerations for a constant acceleration trajectory segment.
%
%   [trajectory_segment, velocities, accelerations, end_velocity] = generateConstantAccelerationSegment(start_pos, start_velocity, acceleration, duration, dt)
%
%   Inputs:
%       start_pos (2x1 vector): Initial position [x0; y0].
%       start_velocity (2x1 vector): Initial velocity [vx0; vy0].
%       acceleration (2x1 vector): Constant acceleration [ax; ay].
%       duration (scalar): Duration of the segment.
%       dt (scalar): Time step.
%
%   Outputs:
%       trajectory_segment (Nx2 matrix): A matrix where each row is a point [x, y].
%       velocities (Nx2 matrix): A matrix where each row is a velocity vector [vx, vy].
%       accelerations (Nx2 matrix): A matrix where each row is an acceleration vector [ax, ay].
%       end_velocity (2x1 vector): The velocity at the end of the segment [vxf; vyf].

    num_steps = floor(duration / dt);
    trajectory_segment = zeros(num_steps + 1, 2);
    velocities = zeros(num_steps + 1, 2);
    accelerations = zeros(num_steps + 1, 2);

    x0 = start_pos(1);
    y0 = start_pos(2);
    vx0 = start_velocity(1);
    vy0 = start_velocity(2);
    ax = acceleration(1);
    ay = acceleration(2);

    % Acceleration is constant for this segment
    constant_acceleration = [ax, ay];

    for i = 0:num_steps
        t = i * dt;
        current_x = x0 + vx0 * t + 0.5 * ax * t^2;
        current_y = y0 + vy0 * t + 0.5 * ay * t^2;
        current_vx = vx0 + ax * t;
        current_vy = vy0 + ay * t;

        trajectory_segment(i + 1, :) = [current_x, current_y];
        velocities(i + 1, :) = [current_vx, current_vy];
        accelerations(i + 1, :) = constant_acceleration;
    end

    % Calculate end velocity explicitly
    end_velocity = start_velocity + acceleration * duration;

end
