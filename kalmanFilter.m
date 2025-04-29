classdef kalmanFilter < handle
    %KF Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

        % System Matrices
        A
        B
        H

        % Noise Matrices
        Q 
        R

        % State Estimate 
        x_hat
        P

        % Dimensions
        m 
        n
        p

    end
    
    methods
        function obj = kalmanFilter(A,B,H,Q,R,x0,P0)

            % ------ Input Validatation -----
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



        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

