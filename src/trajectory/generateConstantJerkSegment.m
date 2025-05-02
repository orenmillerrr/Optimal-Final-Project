function [trajectory_segment, velocities, accelerations, end_velocity, end_acceleration] = generateConstantJerkSegment(start_pos, start_velocity, start_acceleration, jerk, duration, dt)
% generateConstantJerkSegment: Generates points, velocities, and accelerations for a constant jerk trajectory segment.
%
%   [trajectory_segment, velocities, accelerations, end_velocity, end_acceleration] = generateConstantJerkSegment(start_pos, start_velocity, start_acceleration, jerk, duration, dt)
%
%   Inputs:
%       start_pos (2x1 vector): Initial position [x0; y0].
%       start_velocity (2x1 vector): Initial velocity [vx0; vy0].
%       start_acceleration (2x1 vector): Initial acceleration [ax0; ay0].
%       jerk (2x1 vector): Constant jerk [jx; jy].
%       duration (scalar): Duration of the segment.
%       dt (scalar): Time step.
%
%   Outputs:
%       trajectory_segment (Nx2 matrix): A matrix where each row is a point [x, y].
%       velocities (Nx2 matrix): A matrix where each row is a velocity vector [vx, vy].
%       accelerations (Nx2 matrix): A matrix where each row is an acceleration vector [ax, ay].
%       end_velocity (2x1 vector): The velocity at the end of the segment [vxf; vyf].
%       end_acceleration (2x1 vector): The acceleration at the end of the segment [axf; ayf].

    num_steps = floor(duration / dt);
    trajectory_segment = zeros(num_steps + 1, 2);
    velocities = zeros(num_steps + 1, 2);
    accelerations = zeros(num_steps + 1, 2);

    x0 = start_pos(1);
    y0 = start_pos(2);
    vx0 = start_velocity(1);
    vy0 = start_velocity(2);
    ax0 = start_acceleration(1);
    ay0 = start_acceleration(2);
    jx = jerk(1);
    jy = jerk(2);

    for i = 0:num_steps
        t = i * dt;

        % Kinematic equations for constant jerk
        current_ax = ax0 + jx * t;
        current_ay = ay0 + jy * t;

        current_vx = vx0 + ax0 * t + 0.5 * jx * t^2;
        current_vy = vy0 + ay0 * t + 0.5 * jy * t^2;

        current_x = x0 + vx0 * t + 0.5 * ax0 * t^2 + (1/6) * jx * t^3;
        current_y = y0 + vy0 * t + 0.5 * ay0 * t^2 + (1/6) * jy * t^3;

        trajectory_segment(i + 1, :) = [current_x, current_y];
        velocities(i + 1, :) = [current_vx, current_vy];
        accelerations(i + 1, :) = [current_ax, current_ay];
    end

    % Calculate end velocity and acceleration explicitly
    t_end = duration;
    end_ax = ax0 + jx * t_end;
    end_ay = ay0 + jy * t_end;
    end_acceleration = [end_ax; end_ay];

    end_vx = vx0 + ax0 * t_end + 0.5 * jx * t_end^2;
    end_vy = vy0 + ay0 * t_end + 0.5 * jy * t_end^2;
    end_velocity = [end_vx; end_vy];

end
