classdef IIR_filter < matlab.System
    % IIR_filter. The low pass filter after the D/A
    
    properties
        NRB;			% num of resource blocks
        SCS;			% subcarrier spacing
        isFixPoint;		% use low resolution IIR_filter-layer processing
        
        % CC Filter design specifications
        Fp;				% pass-band frequency
        Fst;			% stop-band frequency
        Ap = 1;			% pass-band ripple (dB)
        Ast = 50;		% stop-band attenuation (dB)
        nbcoeff = 12;	% number of bits for the filter coefficient
        nfilt;			% filter order
        bfilt;			% num coeffs
        afilt;			% denom coeffs
        hd;				% filter design
        spec;
        hd2;
        
        % Carrier aggregation
        enableCA;
        ncc;		% number of component carriers
        fcc;		% component carrier center frequency
        
        % DAC
        fsamp;
        nbdac;
    end
    
    methods
        function obj = IIR_filter(varargin)
            % Constructor
            
            % Set key-value pair arguments
            if nargin >= 1
                obj.set(varargin{:});
            end
        end
    end
    
    methods (Access = protected)
        function setupImpl(obj)
            
            % Create a fixed-point IIR low pass filter for low-resolution
            % baseband processing.
            
            % Find the effective signal bandwidth:
            % ResourceBlocks * 12 * SubCarrierSpacing
            fsig = obj.NRB * 12 * obj.SCS * 1e3;
            
            obj.Fp = fsig/obj.fsamp;	% pass-band frequency
            obj.Fst = 1/obj.ncc;		% stop-band frequency
            
            
            obj.spec = fdesign.lowpass('Fp,Fst,Ap,Ast', ...
                obj.Fp, obj.Fst, obj.Ap, obj.Ast);
            f = design(obj.spec, 'ellip');

            % Design a Fixed-Point Filter
            obj.hd = f;
                
            if obj.isFixPoint
                
                f = design(obj.spec,'ellip','SystemObject',true);
                
                sosMatrix = f.SOSMatrix;
                sclValues = f.ScaleValues;

                b = repmat(sclValues(1:(end-1)),1,3) .* sosMatrix(:,(1:3));
                a = sosMatrix(:,(5:6));
                num = b'; % matrix of scaled numerator sections
                den = a'; % matrix of denominator sections
                
                obj.bfilt = fi(num, true, obj.nbcoeff, ...
                    'RoundingMethod','Nearest', 'OverflowAction', 'Saturate');
                obj.afilt = fi(den, true, obj.nbcoeff, ...
                    'RoundingMethod','Nearest', 'OverflowAction', 'Saturate');
                obj.hd = dsp.BiquadFilter('SOSMatrixSource','Input port', ...
                        'ScaleValuesInputPort',false);

            end
        end
        
        function y = stepImpl(obj, x)
            
            if obj.isFixPoint
                % Integer real input from DAC with nbits resolution
                xi = fi(x, 1, max(obj.nbcoeff,obj.nbdac));
                y = double(obj.hd(xi, obj.bfilt, obj.afilt));
                % Power scale
                y = y./sqrt(diag(y'*y))'.*sqrt(diag(x'*x))';
            else
                y=  filter(obj.hd,x);
            end
        end
    end
end
