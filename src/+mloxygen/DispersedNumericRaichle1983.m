classdef DispersedNumericRaichle1983 < handle & mloxygen.Raichle1983
	%% DISPERSEDNUMERICRAICHLE1983S  

	%  $Revision$
 	%  was created 10-Sep-2020 22:24:50 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    properties (Constant)
        LENK = 4
    end
    
    methods (Static)
        function [this,ho,aif] = createFromDeviceKit(devkit, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param ho is numeric, default from devkit.
            %  @param model is an mloxygen.DispersedRaichle1983Model.
            %  @param roi ...
            %  @param solver is in {'simulanneal'}, default := 'simulanneal'.            
            %  @param blurHo := {[], 0, 4.3, ...}            
            %  @param sigma0, default from mloptimization.SimulatedAnnealing.
            %  @return this.
            %  @return ho, blurred by ipr.blurHo.
            %  @return aif.
            
            reshapeArterial = @mloxygen.DispersedNumericRaichle1983.reshapeArterial;
            reshapeScanner = @mloxygen.DispersedNumericRaichle1983.reshapeScanner;
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'ho', [], @isnumeric)
            addParameter(ip, 'model', [], @(x) isa(x, 'mloxygen.DispersedRaichle1983Model'))
            addParameter(ip, 'roi', 'brain_222.4dfp.hdr')
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
            
            ipr.model = set_times_sampled(ipr.model, hoTimesMid);
            ipr.model = set_artery_interpolated(ipr.model, aif);
            this = mloxygen.DispersedNumericRaichle1983( ...
                'devkit', devkit, ...
                'ho', ho, ...
                'solver', 'simulanneal', ...
                'Dt', Dt, ...
                'model', ipr.model, ...
                'times_sampled', hoTimesMid, ...
                'artery_interpolated', aif, ...
                'roi', ipr.roi, ...
                varargin{:});
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
            aifTransition = mean(aif(idxF-10:idxF));
            aif(idxF+1:end) = aifTransition*2.^(-(1:lenLate)/halflife);
        end
        function [aifTimes,aif,Dt] = reshapeArterial(arterial, scanner)
            %% 1.  form uniform time coordinates consistent with scanner
            %  2.  sample aif on uniform time coordinates
            %  3.  infer & apply shift of worldline, Dt, for aif
            
            extrapolateEarlyLateAif = @mloxygen.DispersedNumericRaichle1983.extrapolateEarlyLateAif;
            reshapeScanner = @mloxygen.DispersedNumericRaichle1983.reshapeScanner;
             
            [hoTimesMid,~,hoTimeMin] = reshapeScanner(scanner); % hoTimesMid(1) ~ 0; hoTimeMin ~ -5
            %[hoTimesMid,ho,hoTimeMin] = reshapeScanner(scanner); % hoTimesMid(1) ~ 0; hoTimeMin ~ -5
            aifTimes = hoTimesMid(1):hoTimesMid(end); % aifTimes ~ [0 ... 574]
            hoDatetime0 = scanner.datetime0 + seconds(hoTimeMin); % hoDatetime0 ~ 23-May-2019 12:04:20
            aifDatetime0 = arterial.datetime0; % aifDatetime0 ~ 23-May-2019 12:03:59
            Ddatetime = seconds(aifDatetime0 - hoDatetime0); % Ddatetime ~ -21
            aifTimes__ = (arterial.time0:arterial.timeF) - arterial.time0 + Ddatetime + hoTimesMid(1); 
            % aifTimes__ ~ [-26 -25 -24 ... 181]; satisfies 1.
            
            aif__ = arterial.activityDensity();
            aif__ = extrapolateEarlyLateAif(aif__);
            aif = makima([aifTimes__ aifTimes(end)], [aif__ 0], aifTimes); % satisfies 2.             
            
            Dt = arterial.Dt;
        end
        function [timesMid,ho,timeMin] = reshapeScanner(scanner)
            %% prepends frames to scanner.activityDensity, then resamples
            
            timesMid_ = scanner.timesMid;
            ho = scanner.activityDensity();

            % create timesMid < 0 & interpolated ho
            timesMid = [-flip(timesMid_(1:2)) timesMid_];
            ho = makima([-timesMid_(2) timesMid_], [0 ho], timesMid);
            
            % reset timesMid(1) := 0 while preserving sampling structure
            timeMin = min(timesMid);
            timesMid = timesMid - timeMin;
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
        function ks = buildKs(this, varargin)
            this = solve(this, varargin{:});
            ks = [k1(this) k2(this) k3(this) k4(this) k5(this)];
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
        function fs_ = fs(this, varargin)
            %% fs == [f PS lambda Delta Dt]
            %  @param 'typ' is char, understood by imagingType.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'typ', 'mlfourd.ImagingContext2', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            k(1) = k1(this.strategy_, varargin{:});
            k(2) = k2(this.strategy_, varargin{:});
            k(3) = k3(this.strategy_, varargin{:});
            k(4) = k4(this.strategy_, varargin{:});
            k(5) = this.Dt;
             
            roibin = logical(this.roi);
            fs_ = copy(this.roi.fourdfp);
            fs_.img = zeros([size(this.roi) length(k)]);
            for t = 1:length(k)
                img = zeros(size(this.roi), 'single');
                img(roibin) = k(t);
                fs_.img(:,:,:,t) = img;
            end
            fs_.fileprefix = this.sessionData.fsOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
            fs_ = imagingType(ipr.typ, fs_);
        end
        function [k,sk] = k4(this, varargin)
            [k,sk] = k4(this.strategy_, varargin{:});
        end
        function [k,sk] = k5(this, varargin)
            k = this.Dt;
            sk = nan;
        end
        function ks_ = ks(this, varargin)
            ks_ = this.fs(varargin{:});
        end
    end
    
    %% PROTECTED
    
    methods (Access = protected)        
 		function this = DispersedNumericRaichle1983(varargin)
 			%% DISPERSEDNUMERICRAICHLE1983
            %  @param devkit is mlpet.IDeviceKit.
            %  @param ho is numeric.
            %  @param solver is in {'simulanneal'}.
            %  @param blurHo := {[], 0, 4.3, ...}            
            %  @param sigma0, default from mloptimization.SimulatedAnnealing.

 			this = this@mloxygen.Raichle1983(varargin{:});            
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'ho', [], @(x) isnumeric(x))
            addParameter(ip, 'solver', 'simulanneal', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            this.measurement = ipr.ho;
            switch lower(ipr.solver)
                case 'simulanneal'
                    this.strategy_ = mloxygen.DispersedRaichle1983SimulAnneal( ...
                        'context', this, varargin{:});
                otherwise
                    error('mloxygen:NotImplementedError', 'DispersedNumericRaichle1983.ipr.solver->%s', ipr.solver)
            end
        end
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

