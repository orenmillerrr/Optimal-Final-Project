function [trajectory_segment, velocities, accelerations, end_pos, end_heading] = generateConstantTurnRateSegment(start_pos, start_speed, start_heading, turn_rate, duration, dt)
% generateConstantTurnRateSegment: Generates points, velocities, and accelerations for a constant turn rate (constant speed) trajectory segment.
%
%   [trajectory_segment, velocities, accelerations, end_pos, end_heading] = generateConstantTurnRateSegment(start_pos, start_speed, start_heading, turn_rate, duration, dt)
%
%   Inputs:
%       start_pos (2x1 vector): Initial position [x0; y0].
%       start_speed (scalar): Constant speed magnitude.
%       start_heading (scalar): Initial heading in radians (angle from positive x-axis).
%       turn_rate (scalar): Constant turn rate in radians per second.
%       duration (scalar): Duration of the segment.
%       dt (scalar): Time step.
%
%   Outputs:
%       trajectory_segment (Nx2 matrix): A matrix where each row is a point [x, y].
%       velocities (Nx2 matrix): A matrix where each row is a velocity vector [vx, vy].
%       accelerations (Nx2 matrix): A matrix where each row is an acceleration vector [ax, ay].
%       end_pos (2x1 vector): The position at the end of the segment [xf; yf].
%       end_heading (scalar): The heading at the end of the segment in radians.

    num_steps = floor(duration / dt);
    trajectory_segment = zeros(num_steps + 1, 2);
    velocities = zeros(num_steps + 1, 2);
    accelerations = zeros(num_steps + 1, 2);

    x0 = start_pos(1);
    y0 = start_pos(2);
    v = start_speed;
    psi0 = start_heading; % Initial heading
    omega = turn_rate;

    % Handle the special case of zero turn rate (straight line motion)
    if abs(omega) < 1e-9 % Use a small tolerance for floating point comparison
        constant_velocity = [v * cos(psi0), v * sin(psi0)];
        constant_acceleration = [0, 0];
        for i = 0:num_steps
            t = i * dt;
            current_x = x0 + v * cos(psi0) * t;
            current_y = y0 + v * sin(psi0) * t;
            trajectory_segment(i + 1, :) = [current_x, current_y];
            velocities(i + 1, :) = constant_velocity;
            accelerations(i + 1, :) = constant_acceleration;
        end
        end_pos = trajectory_segment(end, :)';
        end_heading = psi0; % Heading remains constant
    else
        % Calculate points using the analytical solution for CTR
        for i = 0:num_steps
            t = i * dt;
            % Analytical equations for position
            current_x = x0 - (v/omega) * sin(psi0) + (v/omega) * sin(psi0 + omega * t);
            current_y = y0 + (v/omega) * cos(psi0) - (v/omega) * cos(psi0 + omega * t);

            % Calculate instantaneous velocity
            current_heading = psi0 + omega * t;
            current_vx = v * cos(current_heading);
            current_vy = v * sin(current_heading);

            % Calculate instantaneous acceleration (centripetal)
            accel_magnitude = v * abs(omega);
            % Acceleration vector points towards the center of the turn, perpendicular to velocity
            accel_heading = current_heading + sign(omega) * (pi/2); % Add 90 deg for left turn, subtract for right
            current_ax = accel_magnitude * cos(accel_heading);
            current_ay = accel_magnitude * sin(accel_heading);


            trajectory_segment(i + 1, :) = [current_x, current_y];
            velocities(i + 1, :) = [current_vx, current_vy];
            accelerations(i + 1, :) = [current_ax, current_ay];
        end
        end_pos = trajectory_segment(end, :)';
        end_heading = psi0 + omega * duration; % Calculate final heading
    end

end
