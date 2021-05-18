classdef DAC < matlab.System
    % DAC class with scaling.
    
    properties	
        % number of bits, 0 indicates no quantization
        nbits = 6;    

        % Output type:
        % "int":  signed int of the form,
        %     q = 2*k+1, k=-M/2 to M/2-1, M=2^nbits
        % "float":  floating point output of the form:
        %     q = stepSize*(k+0.5)
        % If isComplex, then this representation is used in I and Q
        outputType = "int";

        step_bit = [1.5958,0.9957,0.586,0.3352,0.1881,0.1041];
        isComplex = true;    % complex input
        stepSize = 1.0;          % step size
        dither = false;      % enable dithering
                
		% DAC scaling parameters
		nscal = 10000;       % number of samples used for calibration
        
        % Parameters for linear model:
        %    Q(x) = linGain*x + N(0,quantVar),   x~N(0,inputVar)
		inputVar = 1;       % input variance        
        linGain = 1;    
        quantVar = 1;
        quantVar1 = 1;
        mseOpt = 0;         % Optimal MSE in dB
	end
    
    % Quantizer methods.  All methods are static
    methods 
        
        function obj = DAC(varargin)
            % Constructor
            % Set parameters from constructor arguments
            if nargin >= 1
                obj.set(varargin{:});
            end      
        end
        
        % Performs the quantization providing integer and floating
        % point values
        function [qfloat,qint] = qsat(obj, x0)            
            if (obj.nbits == 0)
                % No quantization
                qint = x0;
                qfloat = x0;
                return
            end
            obj.stepSize = obj.step_bit(obj.nbits);
            
            
            % Scale the input
            x = x0 / obj.stepSize;            

            % Dithering
            if obj.dither
                if obj.isComplex
                    d = (rand(size(x))) + 1i*(2*rand(size(x)));
                else
                    d = rand(size(x))-0.5;
                end
                x = x + d;
            end

            M2 = 2^(obj.nbits-1);
            if obj.isComplex
                % Perform quantization for complex signals
                xr = floor(real(x));
                xi = floor(imag(x));
                qr = max(min(xr,M2-1), -M2)+0.5;
                qi = max(min(xi,M2-1), -M2)+0.5;
                qint = qr + 1i*qi;                
            else                          
                % Perform quantization for complex signals
                x = floor(real(x));
                qint = max(min(x,M2-1), -M2)+0.5;
            end
            
            % Scale output to integer or float
            qint = 2*qint;
            qfloat = 0.5*obj.stepSize*qint;
            
            % Remove dithering
            if obj.dither
               qfloat = qfloat - d*obj.stepSize;
            end
                
        end
        
        % Sets scaling values from a different DAC
        function copyScale(obj, dac1)
            obj.stepSize = dac1.stepSize;
            obj.inputVar = dac1.inputVar;
            obj.quantVar =dac1.quantVar;
            obj.linGain = dac1.linGain;
            obj.mseOpt = dac1.mseOpt;
        end
               
    end
    
     methods (Access = protected)
       		
        function q = stepImpl(obj, x)
            % Step function:  Performs the quantization
            % Scale each stream by its variance
            % note that we need a factor of $2$ for when 
            % the signal is complex 
            sqrt_covx = sqrt(diag(x'*x/size(x,1)))';
%             size(sqrt_covx)
            if obj.isComplex
                sqrt_covx = sqrt_covx/sqrt(2);
            end
            x = x./sqrt_covx;
            
            [qf, qi] = obj.qsat(x);
            if strcmp(obj.outputType, 'int')
                q = qi;
            else
                q = qf;
            end
            q = q./sqrt(diag(q'*q/size(q,1)))'.*sqrt_covx*sqrt(2);
        end
     end
end