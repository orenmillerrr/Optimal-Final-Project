classdef (Abstract) measurementblock < handle
    %MEASUREMENTBLOCK A generic Kalman Filter Measurement Block class
    % This class is an abstract implementation of the measurement block of 
    % a Kalman Filter. This class is used to form the observation,
    % innovation covariance, and measurement covariance matrices. It also
    % stands to correct an estimate with incoming measurements using the
    % measurement model.
    %
    % Author: Walter Livingston
    
    properties
        meas_sigmas     (1,:) double = 0;       % Vector of State Sigmas
        H               (:,:) double = NaN;     % Observation Matrix
        aux             (1,1) struct            % Struct of Aux Data
    end
    
    methods
        function obj = measurementblock(meas_sigmas)
            obj.meas_sigmas = meas_sigmas;
        end
    end

    methods(Abstract)
        [H]         = updateObservationMatrix(obj, x_, y_);
        [R]         = calcMeasurementCovarianceMatrix(obj);
        [S]         = calcInnovationCovarianceMatrix(obj, x_, P_);
        [yhat, H]   = update(obj, x_, y_, relinearize);
    end

    methods
        function [obj] = processAuxData(obj, aux)
            obj.aux = aux;       
        end
        function [x, X] = applyError(obj, x_, X_)
            error(['process_aux_data() for the state block has not ' ...
                'been implemented!']);        
        end
    end
end

