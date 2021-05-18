classdef NRgNBTx < matlab.System
    % 5G NR gNB transmitter class
    properties
        carrierConfig;	% carrier configuration
        pdschConfig;	% default PDSCH configuration
        waveformConfig;	% waveform configuration
        ofdmGridLayer;	% before pre-coding nsc x nsym x nlayers
        
        % Transmitted data
        bits;		% tx data bits
        pdschSym;	% tx data symbols
        dmrsSym;	% demodulation reference signals (DM-RS)
        ptrsSym;	% phase tracking reference signals (PT-RS)
        
        Nslot = 0;	% slot number
        fsamp;		% sample frequency of the DAC
        ncc;		% number of component carriers
        fcc;		% component carrier center frequency
        wtx;        % tx beamforming vector
        ratio;      % adjacent band to desired band power ratio 
        
        % DAC
        dac
        nbits_dac = 0;
		adcFOM = 65e-12;	% FOM of the ADC 65 fJ/conv

        % PHY
		phy;			% Physical layer object
		isFixPoint;		% use fixed-point arithmetic
		enableCA;		% enable carrier aggregation
        
        % IIR_filter
        iirf;

    end
    
    methods
        function obj = NRgNBTx(simParam, varargin)
            % Constructor
            
            % Get parameters from simulation parameters
            % Many 5G Toolbox routines do not take classes, the
            % objects need to be converted to older structures.
            obj.carrierConfig = mmwsim.nr.objtostruct(simParam.carrierConfig);
            obj.pdschConfig = mmwsim.nr.objtostruct(simParam.pdschConfig);
            obj.waveformConfig = mmwsim.nr.objtostruct(simParam.waveformConfig);
            
            % Set parameters from constructor arguments
            if nargin >= 1
                obj.set(varargin{:});
            end
            
            % set the sample frequency.
            obj.fsamp = simParam.waveformConfig.SamplingRate;
            
            if obj.ncc == 4
                obj.fcc = [-150e6, -50e6, 50e6, 150e6];
            elseif obj.ncc == 8
                obj.fcc = [-700e6, -500e6, -300e6, -100e6, ...
                    100e6, 300e6, 500e6, 700e6];
            elseif obj.ncc == 2
                obj.fcc = [-100e6, 100e6];
            else
                obj.fcc = 0;
            end
        end
    end
    
    methods (Access = protected)
        
        function setupImpl(obj)
            
            % Create the TX DAC
            obj.dac = DAC('nbits', obj.nbits_dac);
			
            % Create the TX PHY
            obj.phy = PHYTx(...
                'ncc', obj.ncc, ...
                'isFixPoint', obj.isFixPoint, ...
                'enableCA', obj.enableCA, ...
                'NRB', obj.carrierConfig.NRB, ...
                'SCS', obj.carrierConfig.SubcarrierSpacing, ...
                'fsamp', obj.waveformConfig.SamplingRate, ....
                'nbdac', obj.nbits_dac*(obj.nbits_dac~=0)+ 12*(obj.nbits_dac==0));
            % Create the TX IIR_filter
            
            obj.iirf = IIR_filter(...
                'ncc', obj.ncc, ...
                'isFixPoint', obj.isFixPoint, ...
                'enableCA', obj.enableCA, ...
                'NRB', obj.carrierConfig.NRB, ...
                'SCS', obj.carrierConfig.SubcarrierSpacing, ...
                'fsamp', obj.waveformConfig.SamplingRate, ....
                'nbdac', obj.nbits_dac*(obj.nbits_dac~=0)+ 12*(obj.nbits_dac==0));
        end
        
        function x = stepImpl(obj)
            % Set the slot number in the PDSCH config
            obj.pdschConfig.Nslot = obj.Nslot;
            
            % Get information for PDSCH and DM-RS allocations
            [pdschIndices, dmrsIndices, obj.dmrsSym, ...
                ptrsIndices, obj.ptrsSym, ~] = ...
                mmwsim.nr.hPDSCHResources(obj.carrierConfig, obj.pdschConfig);
            
            cplength = sum(obj.waveformConfig.CyclicPrefixLengths(1:obj.waveformConfig.SymbolsPerSlot));
            
            x = zeros(obj.ncc*(cplength + ...
                obj.waveformConfig.Nfft*obj.waveformConfig.SymbolsPerSlot), ...
                size(obj.wtx,1), obj.ncc);
            
            % Loop over the component carriers
            for icc = 1:obj.ncc
                % Create the PDSCH grid before pre-coding
                obj.ofdmGridLayer = zeros(...
                    obj.waveformConfig.NSubcarriers,...
                    obj.waveformConfig.SymbolsPerSlot, ...
                    obj.pdschConfig.NLayers);
                
                % Generate random bits
                bitsPerSym = mmwsim.nr.NRConst.bitsPerSym(...
                    obj.pdschConfig.Modulation);
                nsym = length(pdschIndices);
                nbits = bitsPerSym * nsym;
                obj.bits = randi([0 1], nbits, 1);
                
                % Modulate the bits to symbols
                M = mmwsim.nr.NRConst.modOrder(obj.pdschConfig.Modulation);
                obj.pdschSym(:,icc) = qammod(obj.bits,M,'InputType','bit',...
                    'UnitAveragePower',true);
                
                % Map symbols to OFDM grid
                obj.ofdmGridLayer(pdschIndices) = obj.pdschSym(:,icc);
                obj.ofdmGridLayer(dmrsIndices) = obj.dmrsSym;
                obj.ofdmGridLayer(ptrsIndices) = obj.ptrsSym;
                
                % Perform the OFDM modulation
                xlayer = mmwsim.nr.hOFDMModulate(obj.carrierConfig, ...
                    obj.ofdmGridLayer);
                
                % Pass the data from the PHY                
                xup = double(obj.phy.step(xlayer));

                % Pass the data from the DAC
                xdac = obj.dac.step(double(xup));

                xiir = double(obj.iirf.step(xdac));

                xnco = mmwsim.nr.hCarrierAggregationModulate(xiir, ...
                    obj.fsamp, obj.fcc(icc));

                % Perform the TX beamforming.
                if icc == 1
                    x(:,:,icc) = xnco*obj.wtx';
                else
                    x(:,:,icc) = obj.ratio*xnco*obj.wtx';
                end
            end
            x = sum(x, 3);
            % Increment the slot number
            obj.Nslot = obj.Nslot + 1;
        end
    end
end