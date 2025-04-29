classdef (Abstract) covarianceblock < handle
    % 
    %
    % Author: Walter Livingston
    
    properties
        S (:,:) double = NaN;
        r (1,1) double = 0;
    end

    methods (Abstract)
        [S,reject]  = calcInnovationCovarianceMatrix(obj, filter);
    end
    
    methods
        function [Q] = calcProcessCovarianceMatrix(obj, filter, dt, customStateBlock)
            if exist('customStateBlock', 'var')
                sBlock = customStateBlock;
            else
                sBlock = filter.state_block;
            end

            Q = sBlock.calcProcessCovarianceMatrix(dt);
        end
        function [R] = calcMeasurementCovarianceMatrix(obj, filter, customMeasBlock)
            if exist('customMeasBlock', 'var')
                mBlock = customMeasBlock;
            else
                mBlock = filter.state_block;
            end
            
            R = mBlock.calcMeasurementCovarianceMatrix();
        end
    end
end

