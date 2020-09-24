classdef DispersedImagingRaichle1983 < handle & matlab.mixin.Copyable
	%% DispersedImagingRaichle1983  
    %  builds Raichle models for imaging voxels.
    %  For performance considerations, see also https://blogs.mathworks.com/loren/2012/03/26/considering-performance-in-object-oriented-matlab-code/

	%  $Revision$
 	%  was created 11-Sep-2020 00:49:25 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.
 	
	properties (Constant)
        HALFLIFE = 122.2416; % s
        JITTER = 0.0 % > 0 to aid deep learning
        LENF = 4
        MAX_NORMAL_BACKGROUND = 20 % Bq/mL 		
    end
    
	properties
        artery_plasma_interpolated
        devkit
        ho % scanner.activityDensity(), decaying & calibrated
        metric % flow 1/s
        meanAbsError
        measurement % expose for performance when used by mloxygen.DispersedRaichle1983Strategy
        model       %
        normMeanAbsError
        Nroi
        prediction
        residual
        roi
        roibin % is logical
        taus
 		times_sampled
    end
    
    properties (Dependent) 
        petPointSpread
        metricTag
        regionTag
        sessionData
    end
    
    methods (Static)
        function this = createFromDeviceKit(devkit, varargin)
            %% makes no adjustments of AIF timing
            %  @param devkit is mlpet.IDeviceKit.
            %  @param ho is numeric, default from devkit.
            %  @param roi is understood by mlfourd.ImagingContext2.
            %  @param blurHo := {[], 0, 4.3, ...}
            
            import mloxygen.DispersedImagingRaichle1983.reshapeArterial
            import mloxygen.DispersedImagingRaichle1983.reshapeScanner
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'ho', [], @isnumeric)
            addParameter(ip, 'roi', 'brain.4dfp.hdr', @(x) ~isempty(x))
            addParameter(ip, 'blurHo', 4.3, @isnumeric)
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            
            % scanner provides calibrations, ancillary data
             
            scanner = ipr.devkit.buildScannerDevice();
            scanner = scanner.blurred(ipr.blurHo);
            [hoTimesMid,ho] = reshapeScanner(scanner);
            
            % AIF  
            
            arterial = ipr.devkit.buildArterialSamplingDevice(scanner);
            [~,aif] = reshapeArterial(arterial, scanner);
            
            hoTaus = diff(hoTimesMid);
            hoTaus = [hoTaus hoTaus(end)];
            this = mloxygen.DispersedImagingRaichle1983( ...
                'devkit', ipr.devkit, ...
                'ho', ho, ...
                'taus', hoTaus, ...
                'times_sampled', hoTimesMid, ...
                'artery_sampled', aif, ...
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
                    'DispersedImagingRaichle1983.DTimeToShift.Dt -> %g; forcing -> %g', Dt, -t_a(idx_a))
                Dt = -t_a(idx_a);
            end
            if Dt > t_a(end)
                warning('mloxygen:ValueError', ...
                    'DispersedImagingRaichle1983.DTimeToShift.Dt -> %g; forcing -> 0', Dt)
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
            
            import mloxygen.DispersedImagingRaichle1983.DTimeToShift 
            import mloxygen.DispersedImagingRaichle1983.extrapolateEarlyLateAif  
            import mloxygen.DispersedImagingRaichle1983.reshapeScanner   
            import mloxygen.DispersedImagingRaichle1983.shiftWorldlines            
                        
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
            
            ho_avgt = asrow(mean(mean(mean(ho, 1), 2), 3));
            Dt = DTimeToShift(aifTimes, aif, hoTimesMid, ho_avgt);
            aif = shiftWorldlines(aif, Dt);
        end
        function [timesMid,ho] = reshapeScanner(scanner)
            %% prepends frames to scanner.activityDensity, the resamples
            
            timesMid_ = scanner.timesMid;
            ho = scanner.activityDensity();
            timesMid = [-flip(timesMid_(1:2)) timesMid_];
            sz = size(ho);
            sz_ = sz; 
            sz_(end) = sz_(end) + 1;
            ho_ = zeros(sz_);
            ho_(:,:,:,2:sz_(end)) = ho;
            ho = makima([-timesMid_(2) timesMid_], ho_, timesMid);
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
        
        function g = get.metricTag(this)
            g = this.sessionData.metricTag;
        end
        function g = get.petPointSpread(this)
            g = this.sessionData.petPointSpread;
        end
        function g = get.regionTag(this)
            g = this.sessionData.regionTag;
        end
        function g = get.sessionData(this)
            g = this.devkit.sessionData;
        end
        
        %%		
        
        function [ic,nic] = buildMeanAbsError(this)
            % @return \Sigma_i activity_i \frac{tau_i}/{T}
            
            if isempty(this.residual)
                this.buildResidual()
            end
            assert(~isempty(this.taus))
            sesd = this.devkit.sessionData;
            hoDCorr = sesd.hoOnAtlas('typ', 'mlfourd.ImagingContext2');
            
            % MAE
            ic = abs(copy(this.residual));
            ic = ic.timeAveraged('taus', this.taus);
            ic.fileprefix = [hoDCorr.fileprefix this.metricTag this.regionTag '_MAE'];
            ic = ic.blurred(this.petPointSpread);
            this.meanAbsError = ic;
            
            % NMAE
            hoAvgt = copy(hoDCorr);
            hoAvgt = hoAvgt.timeAveraged('taus', this.taus);         
            nic = copy(this.meanAbsError) ./ hoAvgt;
            nic.fileprefix = [hoDCorr.fileprefix this.metricTag this.regionTag '_NMAE'];
            nic = nic.blurred(this.petPointSpread);
            this.normMeanAbsError = nic;
        end
        function ic = buildPrediction(this, varargin)
            %% @param reuseExisting is logical; default is false.
            
            ip = inputParser;
            addParameter(ip, 'reuseExisting', false, @islogical)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            assert(~isempty(this.metric))
            sesd = this.devkit.sessionData;
            scanner = this.devkit.buildScannerDevice();
            
            % init working_ifc from hoOnAtlas
            % return existing if reuseExisting
            working_ifc = sesd.hoOnAtlas('typ', 'mlfourd.ImagingFormatContext');
            fileprefix0 = [working_ifc.fileprefix this.metricTag this.regionTag '_predicted'];
            working_ifc.fileprefix = fileprefix0;
            sz = size(working_ifc);
            if ipr.reuseExisting && isfile(working_ifc.fqfilename)
                ic = mlfourd.ImagingContext2(working_ifc.fqfilename);
                this.prediction = ic;
                return
            end            
            
            % represent prediction in R^2:  voxels^1 (x) times            
            fs_ = zeros(this.Nroi,this.LENF+1);
            for i = 1:this.LENF+1     
                rate = this.metric.fourdfp.img(:,:,:,i);                
                fs_(:, i) = rate(this.roibin);
            end         
            this.ensureModel() % without voxelwise adjustments of AIF timings  
            img2d = zeros(this.Nroi, sz(4)); 
            idx0 = length(this.model.times_sampled) - sz(4); % N.B. includes extra early times lab frame
            for vxl = 1:this.Nroi                
                sim = this.model.simulated(fs_(vxl,:)); % adjust AIF timings
                img2d(vxl,:) = sim(1+idx0:sz(4)+idx0);
            end
            
            % embed prediction in R^4:  voxels^3 (x) times
            working_ifc.img = zeros(sz);
            for t = 1:sz(4)
                frame = zeros(sz(1:3));
                frame(this.roibin) = img2d(:,t);
                working_ifc.img(:,:,:,t) = frame;
            end
            
            % decay-correct & remove calibrations
            ic = scanner.decayCorrectLike(working_ifc);
            ic = ic ./ scanner.invEfficiencyf(sesd);
            ic.fileprefix = fileprefix0;
            this.prediction = ic;
        end
        function ic = buildResidual(this)
            if isempty(this.prediction)
                this.buildPrediction()
            end            
            sesd = this.devkit.sessionData;
            hoDCorr = sesd.hoOnAtlas('typ', 'mlfourd.ImagingContext2');
            ic = this.prediction - hoDCorr;
            ic.fileprefix = [hoDCorr.fileprefix this.metricTag this.regionTag  '_residual'];
            ic = ic.blurred(6); % inspired by BOLD spatial scales
            ic = ic .* this.prediction.binarized();
            %this.residual = ic;
        end
        function ic = buildResidualWmparc1(this)
            if isempty(this.prediction)
                this.buildPrediction()
            end            
            sesd = this.devkit.sessionData;
            hoDCorr = sesd.hoOnAtlas('typ', 'mlfourd.ImagingContext2');
            ic = this.prediction - this.parcellateWithWmparc1(hoDCorr);
            ic.fileprefix = [hoDCorr.fileprefix this.metricTag this.regionTag  '_residual'];
            ic = ic .* this.prediction.binarized();
            this.residual = ic;
        end
        function this = ensureModel(this)
            %% without voxelwise adjustments of AIF timings
            
            if isempty(this.model)                               
                this.model = mloxygen.Raichle1983Model.createFromMetric( ...
                    this.sessionData.metric, ...
                    'times_sampled', this.times_sampled, ...
                    'artery_interpolated', this.artery_plasma_interpolated);
            end            
        end
        function ic = parcellateWithWmparc1(this, ic)
            
            sesd = this.devkit.sessionData;
            sesd.region = 'wmparc1';
            wmparc1 = sesd.wmparc1OnAtlas('typ', 'mlfourd.ImagingContext2');
            wmparc1_ = wmparc1.fourdfp;
            
            ic = ic.blurred(4.3);
            ic = ic.masked(wmparc1.binarized);
            ic_ = copy(ic.fourdfp);
            ic_.fileprefix = [ic_.fileprefix '_parcWithWmparc1'];
            img_ = reshape(ic_.img, [128*128*75 size(ic,4)]); % 2D
            img__ = zeros(size(img_));
            
            indices = [1000:1035 2000:2035 3000:3035 4000:4035 5001:5002 6000 1:85 192:255];
            for idx = indices
                roibin_ = wmparc1_.img == idx;
                roivec_ = reshape(roibin_, [128*128*75 1]);
                if 0 == dipsum(roivec_)
                    continue
                end
                
                for t = 1:size(ic,4)
                    img__(roivec_,t) = mean(img_(roivec_,t));
                end
            end
            ic_.img = reshape(img__, [128 128 75 size(ic,4)]);
            ic = mlfourd.ImagingContext2(ic_);
        end
  	end 

    %% PROTECTED
    
	methods (Access = protected)
        function this = DispersedImagingRaichle1983(varargin)
            %% DispersedImagingRaichle1983
            %  @param devkit is mlpet.IDeviceKit.
            %  @param ho is understood by mlfourd.ImagingContext2.
            %  @param times_sampled from scanner is numeric.
            %  @param artery_sampled from counter is numeric.
            %  @param roi is understood by mlfourd.ImagingContext2.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'ho', [])
            addParameter(ip, 'taus', [], @isnumeric)
            addParameter(ip, 'times_sampled', [], @isnumeric)
            addParameter(ip, 'artery_sampled', [], @isnumeric)
            addParameter(ip, 'roi', [], @(x) ~isempty(x))
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.devkit = ipr.devkit;
            this.ho = mlfourd.ImagingContext2(ipr.ho);
            assert(4 == ndims(this.ho))
            this.taus = ipr.taus;
            this.times_sampled = ipr.times_sampled;  
            
            % artery management            
            t = this.times_sampled(1):this.times_sampled(end);
            this.artery_plasma_interpolated = pchip(0:length(ipr.artery_sampled)-1, ipr.artery_sampled, t);
            
            this.roi = mlfourd.ImagingContext2(ipr.roi);
            this.roibin = this.roi.fourdfp.img > 0;
            this.ho = this.masked(this.ho);
            this.Nroi = dipsum(this.roibin);
        end
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
        end
        function ic = masked(this, ic)
            %% retains original fileprefix.
            
            fp = ic.fileprefix;
            ic = ic.masked(this.roi);
            ic.fileprefix = fp;
        end
  	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

