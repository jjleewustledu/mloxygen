classdef NumericRaichle1983 < handle & mloxygen.Raichle1983
	%% NUMERICRAICHLE1983  
    %  bulds Raichle models for regions such as FreeSurfer regions.

	%  $Revision$
 	%  was created 10-Sep-2020 16:21:21 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.
    
    methods (Static)
        function [this,ho,aif] = createFromDeviceKit(devkit, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param ho is numeric, default from devkit.
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}, default := 'simulanneal'.
            %  @param roi ...
            %  @param blurHo := {[], 0, 4.3, ...}          
            %  @return this.
            %  @return ho, blurred by ipr.blurHo.
            %  @return aif.
            
            import mloxygen.NumericRaichle1983.DTimeToShift
            import mloxygen.NumericRaichle1983.reshapeArterial
            import mloxygen.NumericRaichle1983.reshapeScanner
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'ho', [], @isnumeric)
            addParameter(ip, 'roi', 'brain_222.4dfp.hdr')
            addParameter(ip, 'blurHo', 4.3, @isnumeric)
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            
            % prepare atlas data
            
            sesd = ipr.devkit.sessionData;
            sesd.jitOnAtlas(sesd.wmparc1OnAtlas())
            sesd.jitOnAtlas(sesd.hoOnAtlas())
            
            % scanner provides calibrations, ancillary data
            
            ipr.roi = mlfourd.ImagingContext2(ipr.roi);
            roibin = ipr.roi.binarized();            
            scanner = ipr.devkit.buildScannerDevice();
            scanner = scanner.blurred(ipr.blurHo);
            scanner = scanner.volumeAveraged(roibin);
            [ho,timesMid1] = reshapeScanner(scanner.timesMid, scanner.activityDensity()); % calibrated, decaying
                        
            % AIF            
            % Dt shifts the AIF in time:  Dt < 0 shifts left; Dt > 0 shifts right.            
            
            arterial = ipr.devkit.buildArterialSamplingDevice(scanner);
            Ddatetime = seconds(scanner.datetime0 - arterial.datetime0); % Ddatetime ~ 62
            Dt = DTimeToShift(arterial, scanner); % Dt ~ 5
            arterialTimes = arterial.times(arterial.index0:arterial.indexF) - arterial.time0 + Dt - Ddatetime; % ~ [-57 145]
            aif = reshapeArterial(arterialTimes, arterial.activityDensity(), scanner.timesMid);            
            fp = sprintf('mloygen_NumericRaichle1983_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS'));  
            this = mloxygen.NumericRaichle1983( ...
                'devkit', devkit, ...
                'ho', ho, ...
                'solver', 'simulanneal', ...
                'times_sampled', timesMid1, ...
                'artery_interpolated', aif, ...
                'Dt', Dt, ...
                'fileprefix', fp, ...
                'roi', roibin, ...
                varargin{:});
        end
        function Dt = DTimeToShift(varargin)
            %% to apply to arterial to match diff(scanner), in sec.
            
            ip = inputParser;
            addRequired(ip, 'arterial')
            addRequired(ip, 'scanner')
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            Ddatetime = seconds(ipr.scanner.datetime0 - ipr.arterial.datetime0); % Ddatetime ~ 62
            
            a          = ipr.arterial;
            t_a        = asrow(ipr.arterial.times(a.index0:a.indexF) - a.time0);
            activity_a = asrow(ipr.arterial.activityDensity());
            t_s        = asrow(ipr.scanner.timesMid);
            activity_s = asrow(ipr.scanner.activityDensity());
            
            unif_t = 0:max([t_a t_s]);
            unif_activity_s = pchip(t_s, activity_s, unif_t);
            d_activity_s = diff(unif_activity_s); % uniformly sampled time-derivative
            
            % shift dcv in time to match inflow with dtac    
            % use 0.1 of max since counting SNR >> 10 and idx_scanner ~ 1
            [~,idx_a] = max(activity_a > 0.5*max(activity_a));
            [~,idx_s] = max(d_activity_s > 0.5*max(d_activity_s));
            Dt = unif_t(idx_s) - t_a(idx_a); % Dt ~ -57
            if Dt < -t_a(idx_a)
                warning('mloxygen:ValueError', ...
                    'NumericRaichle1983.DTimeToShift.Dt -> %g; forcing -> %g', Dt, -t_a(idx_a))
                Dt = -t_a(idx_a);
            end
            if Dt > t_a(end)
                warning('mloxygen:ValueError', ...
                    'NumericRaichle1983.DTimeToShift.Dt -> %g; forcing -> 0', Dt)
                Dt = 0;
            end
            
            Dt = Dt + Ddatetime; % Dt ~ 5
        end
        function aif = reshapeArterial(times, activityDensity, scannerTimesMid)
            len07 = floor(0.7*length(times));
            len08 = floor(0.8*length(times));
            activityDensityLast = mean(activityDensity(len07:len08));
            times1 = times(1:len08);
            activityDensity1 = activityDensity(1:len08);
            for t = times(len08)+1:scannerTimesMid(end)
                times1 = [times1 t]; %#ok<AGROW>
                activityDensity1 = [activityDensity1 activityDensityLast*2^(-(t - times(len08))/122.2416)]; %#ok<AGROW>
            end
            aif = makima(times1, activityDensity1, ...
                         -scannerTimesMid(2):scannerTimesMid(end));
        end
        function [ho,timesMid1] = reshapeScanner(timesMid, activityDensity)
            timesMid1 = [-flip(timesMid(1:2)) timesMid];
            ho = makima([-timesMid(2) timesMid], [0 activityDensity], timesMid1);
        end
    end

    %% PROTECTED
    
	methods (Access = protected)
 		function this = NumericRaichle1983(varargin)
 			%% NUMERICRAICHLE1983
            %  @param ho is numeric.
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}.

 			this = this@mloxygen.Raichle1983(varargin{:});	
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'ho', [], @(x) isnumeric(x))
            addParameter(ip, 'solver', 'simulanneal', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            this.measurement = ipr.ho;
            this.model = mloxygen.Raichle1983Model(varargin{:});
            switch lower(ipr.solver)
                case 'nest'
                    this.strategy_ = mloxygen.Raichle1983Nest( ...
                        'context', this, varargin{:});
                case 'simulanneal'
                    this.strategy_ = mloxygen.Raichle1983SimulAnneal( ...
                        'context', this, varargin{:});
                case 'hmc'
                    this.strategy_ = mloxygen.Raichle1983HMC( ...
                        'context', this, varargin{:});
                case 'lm'
                    this.strategy_ = mloxygen.Raichle1983LM( ...
                        'context', this, varargin{:});
                case 'bfgs'
                    this.strategy_ = mloxygen.Raichle1983BFGS( ...
                        'context', this, varargin{:});
                otherwise
                    error('mloxygen:NotImplementedError', 'NumericRaichle1983.ipr.solver->%s', ipr.solver)
            end
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

