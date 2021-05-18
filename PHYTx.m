classdef PHYTx < matlab.System
    % PhyTx. The physical layer processing befor zoo the D/A
    
    properties
        NRB;			% num of resource blocks
        SCS;			% subcarrier spacing
        isFixPoint;		% use low resolution PHY-layer processing
        
        % CC Filter design specifications
        Fp;				% pass-band frequency
        Fst;			% stop-band frequency
        Ap = 1;			% pass-band ripple (dB)
        Ast = 50;		% stop-band attenuation (dB)
        nbcoeff = 6;	% number of bits for the filter coefficient
        nfilt;			% filter order
        bfilt;			% num coeffs
        afilt;			% denom coeffs
        hd;				% filter design
        
        % Carrier aggregation
        enableCA;
        ncc;		% number of component carriers
        fcc;		% component carrier center frequency
        
        % DAC
        fsamp;
        nbdac;
    end
    
    methods
        function obj = PHYTx(varargin)
            % Constructor
            
            % Set key-value pair arguments
            if nargin >= 1
                obj.set(varargin{:});
            end
        end
    end
    
    methods (Access = protected)
        function setupImpl(obj)
            
            % Create a fixed-point FIR low pass filter for low-resolution
            % baseband processing.
            if obj.isFixPoint
                % Find the effective signal bandwidth:
                % ResourceBlocks * 12 * SubCarrierSpacing
                fsig = obj.NRB * 12 * obj.SCS * 1e3;
                
                obj.Fp = fsig/obj.fsamp;	% pass-band frequency
                obj.Fst = 1/obj.ncc;		% stop-band frequency
                               
				% Design a Fixed-Point Filter
				spec = fdesign.lowpass('Fp,Fst,Ap,Ast', ...
					obj.Fp, obj.Fst, obj.Ap, obj.Ast);
				f = design(spec, 'minphase', false, 'SystemObject', true);

				coef = coeffs(f);

				bq = fi(coef.Numerator, true, obj.nbcoeff, ...
					'RoundingMethod','Nearest', 'OverflowAction', 'Saturate');
				L = bq.FractionLength;
				bsc = coef.Numerator*2^L;

				obj.hd = dfilt.dffir(bsc);
				obj.hd.Arithmetic = 'fixed';
				obj.hd.CoeffWordLength = obj.nbcoeff;

				% Integer real input from DAC with nbits resolution
                obj.hd.FilterInternals = 'SpecifyPrecision';
    
				obj.hd.OutputWordLength = obj.nbdac;
				obj.hd.OutputFracLength = 0;

				obj.bfilt = obj.hd.Numerator;
				obj.afilt = 1;
                
            end
        end
        
        function y = stepImpl(obj, x)
            
            if obj.isFixPoint
                % upsample and then filter the input signal with
                % a filter that has low-resolution coefficients
                xup = sqrt(obj.ncc)*upsample(x, obj.ncc);
                y = filter(obj.bfilt, obj.afilt, xup);
                % Power scale
                y = y./sqrt(diag(y'*y))'.*sqrt(diag(x'*x))';

            else
                % Use an High order low-pass filter with upsample
                y=  sqrt(1/obj.ncc)*resample(x,obj.ncc,1,50, 20); %
            end
        end
    end
end
