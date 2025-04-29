% kalmanFilter Class Definition
%
% This class implements a standard linear Kalman filter.
% It encapsulates the state, covariance, system matrices, and noise parameters.
% It provides methods for the prediction and update steps.

classdef kalmanFilter < handle
    
    properties

        A       % State Transition Matrix
        B       % Input Matrix
        H       % Observation Matrix
        Q       % Process Noise Covariance Matrix
        R       % Measurement Noise Covariance Matrix
        x_hat   % State Estimate Matrix
        P       % State Estimate Covariance Matrix

        Z       % State Innovation
        S       % Inovation Covariance

        m       % Number of States
        n       % Number of Inputs
        p       % Number of Measurements


    end
    
    methods

        function obj = kalmanFilter(A,B,H,Q,R,x0,P0)
            % Construct an instance of this class
            %
            % Inputs:
            %   A  - State transition matrix
            %   B  - Control input matrix (can be empty [])
            %   H  - Observation matrix
            %   Q  - Process noise covariance matrix
            %   R  - Measurement noise covariance matrix
            %   x0 - Initial state estimate vector
            %   P0 - Initial error covariance matrix
            
            % ------- Input Validatation ------
            % Inputs
            if nargin < 7
                error('KalmanFilter:NotEnoughInputs', 'Requires at least 7 input arguments (A, B, H, Q, R, x0, P0).')
            end

            % State Transition Matrix
            [obj.n, A_cols] = size(A);
            if obj.n ~= A_cols
                error('KalmanFilter:InvalidDim', 'Matrix A must be square (n x n).');
            end

            % Input Matrix
            if ~isempty(B)
                [n_B, obj.m] = size(B);
                if n_B ~= obj.n
                    error('KalmanFilter:InvalidDim', 'Matrix B must have n rows.');
                end
            else
                obj.m = 0; % No control input
            end

            % Measurement Matrix
            [obj.p, n_H] = size(H);
            if n_H ~= obj.n
                error('KalmanFilter:InvalidDim', 'Matrix H must have n columns.');
            end

            % Process Noise
            [n_Q, n_Q_cols] = size(Q);
            if n_Q ~= obj.n || n_Q_cols ~= obj.n
                error('KalmanFilter:InvalidDim', 'Matrix Q must be square (n x n).');
            end

            % Measurement Noise
            [p_R, p_R_cols] = size(R);
            if p_R ~= obj.p || p_R_cols ~= obj.p
                error('KalmanFilter:InvalidDim', 'Matrix R must be square (p x p).');
            end

            % Initial States
            if length(x0) ~= obj.n || ~iscolumn(x0)
                 error('KalmanFilter:InvalidDim', 'Initial state x0 must be a column vector of size n.');
            end

            [n_P, n_P_cols] = size(P0);
             if n_P ~= obj.n || n_P_cols ~= obj.n
                 error('KalmanFilter:InvalidDim', 'Initial covariance P0 must be square (n x n).');
             end

            % ------ Assign properties to object ------
            obj.A = A;
            obj.B = B;
            obj.H = H;
            obj.Q = Q;
            obj.R = R;
            obj.x_hat = x0; % Initial state estimate
            obj.P = P0;     % Initial error covariance

        end % obj = kalmanFilter(A,B,H,Q,R,x0,P0)
        
        function propagate(obj,u)
            % Perform the Kalman filter prediction step
            %
            % Inputs:
            %   obj - kalmanFilter object instance
            %   u   - Control input vector (m x 1), use [] or omit if no 
            %         control input
            % ** NOTE: This assumes that control input vector (u) is     **
            %          deterministic

            % Check control input dimension
            if obj.m > 0 % If there is a control input
                if nargin < 2 || isempty(u)
                    error('KalmanFilter:MissingInput', 'Control input u is required because B is defined.');
                end
                if length(u) ~= obj.m || ~iscolumn(u)
                    error('KalmanFilter:InvalidDim', 'Control input u must be a column vector of size m.');
                end
                % Predict state: x_hat_minus = A * x_hat + B * u
                obj.x_hat = obj.A * obj.x_hat + obj.B * u;
            else % No control input
                 if nargin > 1 && ~isempty(u)
                    warning('KalmanFilter:UnusedInput', 'Control input u provided but B is empty. Ignoring u.');
                 end
                % Predict state: x_hat_minus = A * x_hat
                obj.x_hat = obj.A * obj.x_hat;
            end

            % Predict error covariance: P_minus = A * P * A' + Q
            obj.P = obj.A * obj.P * obj.A' + obj.Q;

        end % predict(obj,u)

        function update(obj,y)
            % Perform the Kalman filter update step
            %
            % Inputs:
            %   obj - kalmanFilter object instance
            %   y   - Measurement vector (p x 1)

            % Check measurement dimension
             if length(y) ~= obj.p || ~iscolumn(y)
                 error('KalmanFilter:InvalidDim', 'Measurement z must be a column vector of size p.');
             end

            % Calculate Kalman Gain
            obj.S = obj.H * obj.P * obj.H' + obj.R; % Innovation (or residual) covariance
            K = (obj.P * obj.H') / obj.S; % More numerically stable than inv(S)

            % Update state estimate
            obj.Z = y - obj.H * obj.x_hat; % Measurement residual (innovation)
            obj.x_hat = obj.x_hat + K * obj.Z;

            % Update error covariance: P = (I - K * H) * P_minus
            I = eye(obj.n); % Identity matrix
            obj.P = (I - K * obj.H) * obj.P;

        end % update(obj,y)

        function immUpdate(obj,x0,P0)
            % Get data pertaining to kalmanFilter state estimate
            %
            % Inputs:
            %   obj - kalmanFilter object instance
            %   x0 - New x_hat for kalmanFilter (IMM)
            %   P0 - New P for kalmanFilter (IMM)

            obj.x_hat = x0;
            obj.P = P0;
        end
        

        function x_hat = getState(obj)
            % Get data pertaining to kalmanFilter state estimate
            %
            % Inputs:
            %   obj - kalmanFilter object instance
            % Outputs:
            %   x_hat - State Estimate

            x_hat = obj.x_hat;
        end


        function P = getCovariance(obj)
            % Get data pertaining to kalmanFilter covariance
            %
            % Inputs:
            %   obj - kalmanFilter object instance
            % Outputs:
            %   P - Covariance

            P = obj.P;
        end

        function [Z,S] = getInnovation(obj)
            % Get data pertaining to kalmanFilter innovation
            %
            % Inputs:
            %   obj - kalmanFilter object instance
            % Outputs:
            %   Z - Innovation
            %   S - Innovation Covariance

            Z = obj.Z;
            S = obj.S;

        end

        function [L] = computeLikelihood(obj)
            % Get data pertaining to kalmanFilter innovation
            %
            % Inputs:
            %   obj - kalmanFilter object instance
            % Outputs:
            %   L = Likelihood

            L = 1/sqrt(2*pi)*obj.S^(-1/2) * exp(-.5*obj.Z'*obj.S^-1*obj.Z);
        end
    end
end

