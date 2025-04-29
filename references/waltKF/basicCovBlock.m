classdef basicCovBlock < covarianceblock
    % 
    %
    % Author: Walter Livingston
    
    properties
    end

    methods
        function [S, reject] = calcInnovationCovarianceMatrix(obj, filter)
            H = filter.measurement_block.H;
            P = filter.P;
            R = filter.measurement_block.calcMeasurementCovarianceMatrix();
            z = filter.z;
            S = H*P*H' + R;
            obj.S = S;
            obj.r = sqrt(z'*(S^(-1))*z);

            reject = false;
        end
    end
end

