classdef NumericMintun1984 < handle & mloxygen.Mintun1984
	%% NUMERICMINTUN1984  

	%  $Revision$
 	%  was created 07-Dec-2020 21:36:02 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.9.0.1524771 (R2020b) Update 2 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    methods (Static)
        function [this,oo,aif] = createFromDeviceKit(devkit, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param fs_Raichle_Martin ~ [fs PS lambda Dt v1]
            %  @param oo is numeric, default from devkit.
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}, default := 'simulanneal'.
            %  @param roi ...
            %  @param blurOo := {[], 0, 4.3, ...}          
            %  @return this.
            %  @return ho, blurred by ipr.blurHo.
            %  @return aif.
            
            import mloxygen.NumericMintun1984.reshapeArterial
            import mloxygen.NumericMintun1984.reshapeScanner
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'fs_Raichle_Martin', [], @(x) ~isempty(x))
            addParameter(ip, 'oo', [], @isnumeric)
            addParameter(ip, 'roi', 'brain_222.4dfp.hdr')
            addParameter(ip, 'blurOo', 4.3, @isnumeric)
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            
            % prepare atlas data
            
            sesd = ipr.devkit.sessionData;
            sesd.jitOnAtlas(sesd.wmparc1OnAtlas())
            sesd.jitOnAtlas(sesd.ooOnAtlas())
            
            % scanner provides calibrations, ancillary data
            
            ipr.roi = mlfourd.ImagingContext2(ipr.roi);
            roibin = ipr.roi.binarized();            
            scanner = ipr.devkit.buildScannerDevice();
            scanner = scanner.blurred(ipr.blurOo);
            scanner = scanner.volumeAveraged(roibin);
            [oo,timesMid1] = reshapeScanner(scanner.timesMid, scanner.activityDensity()); % calibrated, decaying
                        
            % AIF          
            
            arterial = ipr.devkit.buildArterialSamplingDevice(scanner);
            Ddatetime = seconds(scanner.datetime0 - arterial.datetime0); % Ddatetime ~ 62
            Dt = ipr.fs_Raichle_Martin(4);
            arterialTimes = arterial.times(arterial.index0:arterial.indexF) - arterial.time0 + Dt - Ddatetime; % ~ [-52 140]
            aif = reshapeArterial(arterialTimes, arterial.activityDensity(), scanner.timesMid);            
            fp = sprintf('mloygen_NumericMintun1984_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS'));  
            this = mloxygen.NumericMintun1984( ...
                'devkit', devkit, ...
                'fs_Raichle_Martin', ipr.fs_Raichle_Martin, ...
                'oo', oo, ...
                'solver', 'simulanneal', ...
                'times_sampled', timesMid1, ...
                'artery_interpolated', aif, ...
                'Dt', Dt, ...
                'fileprefix', fp, ...
                'roi', roibin, ...
                varargin{:});
        end
        function aif = reshapeArterial(times, activityDensity, scannerTimesMid)
            aif = mloxygen.NumericRaichle1983.reshapeArterial(times, activityDensity, scannerTimesMid);
        end
        function [ho,timesMid1] = reshapeScanner(timesMid, activityDensity)
            [ho,timesMid1] = mloxygen.NumericRaichle1983.reshapeScanner(timesMid, activityDensity);
        end
    end

	methods (Access = protected)		  
 		function this = NumericMintun1984(varargin)
 			%% NUMERICMINTUN1984
            %  @param oo is numeric.
            %  @param solver is in {'simulanneal'}.
            %  @param fs_Raichle_Martin ~ [fs PS lambda v1].

 			this = this@mloxygen.Mintun1984(varargin{:});
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'oo', [], @(x) isnumeric(x))
            addParameter(ip, 'solver', 'simulanneal', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            this.measurement = ipr.oo;
            this.model = mloxygen.Mintun1984Model(varargin{:});            
            switch lower(ipr.solver)
                case 'simulanneal'
                    this.strategy_ = mloxygen.Mintun1984SimulAnneal( ...
                        'context', this, varargin{:});
                otherwise
                    error('mloxygen:NotImplementedError', 'NumericMintun1984.ipr.solver->%s', ipr.solver)
            end
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

