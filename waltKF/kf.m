classdef kf < handle
    %KF A generic Kalman Filter class
    % This class implements all the basic Kalman Filter equations, but uses
    % State Block and Measurement Block objects to define all the matrices
    % and state propagation functions.
    %
    % Author: Walter Livingston

    properties
        state_block         (1,1)           % State Block Object
        measurement_block   (1,1)           % Measurement Block Object
        covariance_block    (1,1)           % Covariance Block Object
        x                   (:,1) double    % State Vector
        X                   (:,1) double    % Nominal State Vector
        P                   (:,:) double    % State Covariance Matrix
        z                   (:,1) double    % Innovation Vector
        mode                (1,1) string    % Kalman Filter Mode
    end
    
    methods
        function obj = kf(state_block, measurement_block, covariance_block, options)
            arguments
                state_block         (1,1) stateblock
                measurement_block   (1,1) measurementblock
                covariance_block    (1,1) covarianceblock
                options.xi          (:,1) double = zeros(state_block.num_states,1);
                options.X_i         (:,1) double = zeros(state_block.num_states,1);
                options.dt          (1,1) double = 1
                options.mode        (1,1) string = 'kf'
            end

            obj.state_block = state_block;
            obj.measurement_block = measurement_block;
            obj.covariance_block = covariance_block;
            obj.mode = options.mode;
            obj.x = options.xi;

            if strcmp(obj.mode,"error")
                obj.X = options.X_i;

                sAux = obj.state_block.aux;
                sAux.X = obj.X;

                mAux = obj.measurement_block.aux;
                mAux.X = obj.X;
                
                obj.state_block.applyError(obj.x, obj.X, options.dt);
                obj.measurement_block.applyError(obj.x, obj.X);

                obj.state_block = obj.state_block.processAuxData(sAux);
                obj.measurement_block = obj.measurement_block.processAuxData(mAux);
            end
            
            obj.P = covariance_block.calcProcessCovarianceMatrix(obj, options.dt);

            switch obj.mode
                case 'kf'
                    obj.state_block.F = obj.state_block.updateStateTransitionMatrix(options.xi);
                    obj.measurement_block.H = ...
                        obj.measurement_block.updateObservationMatrix(options.xi, ...
                            0);
                case 'extended'

            end
        end
        
        function [obj] = process(obj, dt, customStateBlock, customInnBlock)
        %PROCESS Implements the Time Update of a Kalman Filter
        % This function runs the time update of a Kalman filter. Using the
        % stateblock passed into the constructor, this function creates the
        % state transition and process covariance matrices and propagates
        % the state estimates and state covariance matrix.
            if exist('customStateBlock', 'var')
                sBlock = customStateBlock;
                obj.state_block = customStateBlock;
            else
                sBlock = obj.state_block;
            end

            if exist('customInnBlock', 'var')
                cBlock = customInnBlock;
            else
                cBlock = obj.covariance_block;
            end

            [obj.x, Phi] = sBlock.propagate(obj.x, dt, ...
                (strcmp(obj.mode, 'extended')) || (strcmp(obj.mode, 'error')));
            Qd = cBlock.calcProcessCovarianceMatrix(obj, dt, sBlock);
            obj.P = Phi*obj.P*Phi' + Qd;
            
            if strcmp(obj.mode,'error')
                obj.X = sBlock.applyError(obj.x, obj.X, dt);
                aux = sBlock.aux;
                aux.X = obj.X;
                obj.state_block = sBlock.processAuxData(aux);
            end
        end

        function [obj] = update(obj, y, customMeasBlock, customInnBlock)
            if exist('customMeasBlock', 'var')
                mBlock = customMeasBlock;
                aux = mBlock.aux;
                aux.X = obj.X;
                aux.y = y;
                obj.measurement_block = customMeasBlock.processAuxData(aux);
            else
                mBlock = obj.measurement_block;
            end

            if exist('customInnBlock', 'var')
                cBlock = customInnBlock;
            else
                cBlock = obj.covariance_block;
            end

            [yhat, H] = mBlock.update(obj.x, y, ...
                (strcmp(obj.mode, 'extended')) || (strcmp(obj.mode, 'error')));
            if any(isnan(obj.P)); error(''); end
            R = cBlock.calcMeasurementCovarianceMatrix(obj, mBlock);
            oldz = obj.z;
            obj.z = (y - yhat);
            [S,reject] = cBlock.calcInnovationCovarianceMatrix(obj);
            if ~reject
                L = obj.P*H'/S;
    
                obj.x = obj.x + L*obj.z;
                obj.P = (eye(obj.state_block.num_states) - L*H)*obj.P*(eye(obj.state_block.num_states) - L*H)' + L*R*L';
    
                if strcmp(obj.mode, 'error')
                    [obj.x, obj.X] = mBlock.applyError(obj.x, obj.X);
                    aux = mBlock.aux;
                    aux.X = obj.X;
                    aux.y = y;
                    obj.measurement_block = mBlock.processAuxData(aux);
                end
            else
                obj.z = oldz;
            end
        end
    end
end

