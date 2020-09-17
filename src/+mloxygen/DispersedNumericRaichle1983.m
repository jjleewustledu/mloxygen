classdef DispersedNumericRaichle1983 < handle & mloxygen.Raichle1983
	%% DISPERSEDNUMERICRAICHLE1983S  

	%  $Revision$
 	%  was created 10-Sep-2020 22:24:50 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    properties (Constant)
        LENK = 4
    end
    
    properties (Dependent)
        artery_interpolated
        times_sampled
    end
    
    methods (Static)
        function [this,ho] = createFromDeviceKit(devkit, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param ho is numeric, default from devkit.
            %  @param roi ...
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}, default := 'simulanneal'.            
            %  @param blurHo := {[], 0, 4.3, ...}            
            %  @param map, default := mloxygen.Raichle1983Model.preferredMap().
            %  @param times_sampled non-uniformly scheduled by the time-resolved PET reconstruction.
            %  @param artery_interpolated, default from devkit.
            %  @param sigma0, default from mloptimization.SimulatedAnnealing.
            %  @param fileprefix, default from devkit.    
            %  @param histology is char.
            %  @return this.
            %  @return ho, blurred by ipr.blurHo.
            
            import mloxygen.DispersedNumericRaichle1983.reshapeArterial
            import mloxygen.DispersedNumericRaichle1983.reshapeScanner
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'ho', [], @isnumeric)
            addParameter(ip, 'roi', 'brain.4dfp.hdr', @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'blurHo', 4.3, @isnumeric)
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            
            % prepare atlas data
            
            sesd = ipr.devkit.sessionData;
            sesd.jitOn222(sesd.wmparc1OnAtlas())
            sesd.jitOn222(sesd.hoOnAtlas())
            
            % scanner provides calibrations, decay, ancillary data
            
            ipr.roi = mlfourd.ImagingContext2(ipr.roi);
            roibin = ipr.roi.binarized();            
            scanner = ipr.devkit.buildScannerDevice();
            scanner = scanner.blurred(ipr.blurHo);
            scanner = scanner.volumeAveraged(roibin);
            [hoTimesMid,ho] = reshapeScanner(scanner);
            
            % arterial reshaping makes it quantitatively comparable to scanner
            
            arterial = ipr.devkit.buildArterialSamplingDevice(scanner);
            [~,aif,Dt] = reshapeArterial(arterial, scanner);
            
            % assemble this
            
            fp = sprintf('mloygen_DispersedRaichle1983_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS'));  
            this = mloxygen.DispersedNumericRaichle1983( ...
                devkit, ...
                'ho', ho, ...
                'solver', 'simulanneal', ...
                'times_sampled', hoTimesMid, ...
                'artery_interpolated', aif, ...
                'Dt', Dt, ...
                'fileprefix', fp, ...
                varargin{:});
        end
        function Dt = DTimeToShift(varargin)
            %% Dt by which to shift arterial to match diff(scanner):  Dt < 0 will shift left; Dt > 0 will shift right.
            
            ip = inputParser;
            addRequired(ip, 't_a')
            addRequired(ip, 'activity_a')
            addRequired(ip, 't_s')
            addRequired(ip, 'activity_s')
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            t_a        = ipr.t_a;
            activity_a = ipr.activity_a;
            t_s        = ipr.t_s;
            activity_s = ipr.activity_s;
            
            unif_t = 0:max([t_a t_s]);
            unif_activity_s = pchip(t_s, activity_s, unif_t);
            d_activity_s = diff(unif_activity_s); % uniformly sampled time-derivative
            
            % shift dcv in time to match inflow with dtac
            [~,idx_a] = max(activity_a > 0.95*max(activity_a)); % idx_a ~ 35
            [~,idx_s] = max(d_activity_s > 0.95*max(d_activity_s)); % idx_s ~ 35
            Dt = unif_t(idx_s) - t_a(idx_a); % Dt ~ 5
            if Dt < -t_a(idx_a)
                warning('mloxygen:ValueError', ...
                    'DispersedNumericRaichle1983.DTimeToShift.Dt -> %g; forcing -> %g', Dt, -t_a(idx_a))
                Dt = -t_a(idx_a);
            end
            if Dt > t_a(end)
                warning('mloxygen:ValueError', ...
                    'DispersedNumericRaichle1983.DTimeToShift.Dt -> %g; forcing -> 0', Dt)
                Dt = 0;
            end
        end        
        function aif = extrapolateEarlyLateAif(aif__)
            [~,idx0] = max(aif__ > 0.1*max(aif__));
            idx0 = idx0 - 2;
            baseline = mean(aif__(1:idx0));
            [~,idxF] = max(flip(aif__) > baseline);
            idxF = length(aif__) - idxF + 1;
            
            aif__ = aif__ - baseline;
            aif__(aif__ < 0) = 0;
            
            selection = zeros(size(aif__));
            selection(idx0:idxF) = 1;
            aif = selection .* aif__;
            lenLate = length(aif) - idxF;
            halflife = 122.2416;
            aif(idxF+1:end) = aif(idxF)*2.^(-(1:lenLate)/halflife);
        end
        function [aifTimes,aif,Dt] = reshapeArterial(arterial, scanner)
            %% 1.  form uniform time coordinates consistent with scanner
            %  2.  sample aif on uniform time coordinates
            %  3.  infer & apply shift of worldline, Dt, for aif
            
            import mloxygen.DispersedNumericRaichle1983.DTimeToShift 
            import mloxygen.DispersedNumericRaichle1983.extrapolateEarlyLateAif  
            import mloxygen.DispersedNumericRaichle1983.reshapeScanner   
            import mloxygen.DispersedNumericRaichle1983.shiftWorldlines            
                        
            [hoTimesMid,ho] = reshapeScanner(scanner); % hoTimesMid(1) ~ -5
            aifTimes = hoTimesMid(1):hoTimesMid(end); % aifTimes ~ [-5 -4 -3 ... 569]       
            hoDatetime0 = scanner.datetime0 + seconds(hoTimesMid(1)); % hoDatetime0 ~ 23-May-2019 12:04:20
            aifDatetime0 = arterial.datetime0; % aifDatetime0 ~ 23-May-2019 12:03:59
            Ddatetime = seconds(aifDatetime0 - hoDatetime0); % Ddatetime ~ -21
            aifTimes__ = (arterial.time0:arterial.timeF) - arterial.time0 + Ddatetime + hoTimesMid(1); 
            % aifTimes__ ~ [-26 -25 -24 ... 181]; satisfies 1.
            
            aif__ = arterial.activityDensity();
            aif__ = extrapolateEarlyLateAif(aif__);
            aif = makima([aifTimes__ aifTimes(end)], [aif__ 0], aifTimes); % satisfies 2.             
            
            Dt = DTimeToShift(aifTimes, aif, hoTimesMid, ho);
            aif = shiftWorldlines(aif, Dt);
        end
        function [timesMid,ho] = reshapeScanner(scanner)
            %% prepends frames to scanner.activityDensity, the resamples
            
            timesMid_ = scanner.timesMid;
            ho = scanner.activityDensity();
            timesMid = [-flip(timesMid_(1:2)) timesMid_];
            ho = makima([-timesMid_(2) timesMid_], [0 ho], timesMid);
        end
        function aif = shiftWorldlines(aif__, Dt)
            if Dt == 0
                aif = aif__;
                return
            end
            
            halflife = 122.2416;
            if Dt < 0
                aif = aif__(end)*ones(size(aif__));
                aif(1:(length(aif__)+Dt)) = aif__((1-Dt):end);
                selection = aif > 0.01*max(aif);
                aif(selection) = 2^(-Dt/halflife)*aif(selection);
                return
            end
            aif = aif__(1)*ones(size(aif__));
            aif((1+Dt):end) = aif__(1:(end-Dt));
            selection = aif > 0.01*max(aif);
            aif(selection) = 2^(-Dt/halflife)*aif(selection);
        end
    end

	methods 
        
        %% GET
        
        function g = get.artery_interpolated(this)
            g = this.strategy_.artery_interpolated;
        end
        function g = get.times_sampled(this)
            g = this.strategy_.times_sampled;
        end
        
        %%
		  
        function [k,sk] = k4(this, varargin)
            [k,sk] = k4(this.strategy_, varargin{:});
        end
        function [k,sk] = ks(this, varargin)
            k = zeros(1,this.LENK);
            sk = zeros(1,this.LENK);
            [k(1),sk(1)] = k1(this.strategy_, varargin{:});
            [k(2),sk(2)] = k2(this.strategy_, varargin{:});
            [k(3),sk(3)] = k3(this.strategy_, varargin{:});
            [k(4),sk(4)] = k4(this.strategy_, varargin{:});
        end
        function ho  = checkSimulated(this, varargin)
            %% CHECKSIMULATED simulates tissue activity with passed and internal parameters without changing state.
            %  @param required ks is [k1 k2 k3 k4 Dt].
            %  @param aif is numeric; default is this.artery_interpolated for model state.
            %  @return ho simulation is numeric.
        
            ip = inputParser;
            addOptional(ip, 'ks', this.ks(), @isnumeric)
            addParameter(ip, 'aif', this.artery_interpolated, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;  
            
            ho = this.model.simulated(ipr.ks, 'aif', ipr.aif, 'Dt', this.Dt);
        end
        
 		function this = DispersedNumericRaichle1983(devkit, varargin)
 			%% DISPERSEDNUMERICRAICHLE1983
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param ho is numeric.
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}.
            %  @param blurHo := {[], 0, 4.3, ...}            
            %  @param map, default := mloxygen.Raichle1983Model.preferredMap().
            %  @param times_sampled non-uniformly scheduled by the time-resolved PET reconstruction.
            %  @param artery_interpolated, default from devkit.
            %  @param sigma0, default from mloptimization.SimulatedAnnealing.
            %  @param fileprefix, default from devkit. 
            %  @param histology is char.

 			this = this@mloxygen.Raichle1983(devkit, varargin{:});            
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'ho', [], @(x) isnumeric(x))
            addParameter(ip, 'solver', 'simulanneal', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            this.measurement = ipr.ho;
            this.model = mloxygen.DispersedRaichle1983Model(varargin{:});
            switch lower(ipr.solver)
                case 'simulanneal'
                    this.strategy_ = mloxygen.DispersedRaichle1983SimulAnneal( ...
                        'context', this, varargin{:});
                otherwise
                    error('mloxygen:NotImplementedError', 'Raichle1983.ipr.solver->%s', ipr.solver)
            end
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

