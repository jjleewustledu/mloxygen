classdef DispersedMintun1984Model < mloxygen.Mintun1984Model
	%% DISPERSEDMINTUN1984MODEL  
    %   uses variational ks := [metabf oef] and deterministic fs := [f PS lambda Delta Dt v1].
    %   metabf \in [0 1] is the fraction of the oxygen AIF metabolized to water at 90 s after bolus arrival.
    %   lambda is mL/g.
    %   v1 is mL/mL <= 1.

	%  $Revision$
 	%  was created 05-Dec-2020 20:46:46 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.9.0.1524771 (R2020b) Update 2 for MACI64.  Copyright 2020 John Joowon Lee.

    
	methods (Static)
        function this = createFromSession(sesd, varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'sesd', @(x) isa(x, 'mlpipeline.ISessionData'))
            addParameter(ip, 'roi', [])
            parse(ip, sesd, varargin{:})
            ipr = ip.Results;
            roi = mlfourd.ImagingContext2(ipr.roi);
            roi = roi.binarized();
            
            fs_R_M_ = mloxygen.DispersedMintun1984Model.fs_R_M(sesd, roi);
            this = mloxygen.DispersedMintun1984Model('fs_Raichle_Martin', fs_R_M_, varargin{:});
        end
        function        ensureModelPrereq(sesd, ic)
            if isfile(ic.fqfilename)
                return
            end
            
            ds8 = datestr(sesd.datetime, 'yyyymmdd');
            re = regexp(ic.fileprefix, ['\S+dt' ds8 '(?<time>\d{6})\S+'], 'names');
            ds14 = [ds8 re.time];
            
            pwd0 = pushd(ic.filepath);
            globbed = globT(strrep(ic.filename, ds14, [ds8 '*']));
            ifc_avgens = mlfourd.ImagingFormatContext(globbed{1}); % average ensemble
            if length(globbed) > 1
                for ig = 2:length(globbed)
                    ifc = mlfourd.ImagingFormatContext(globbed{ig});
                    ifc_avgens.img = ifc_avgens.img + ifc.img;
                end
                ifc_avgens.img = ifc_avgens.img / length(globbed);
            end
            ifc_avgens.fileprefix = ic.fileprefix;
            ifc_avgens.save()
            popd(pwd0)
        end
        function vecT = fs_R_M(sesd, roi)
            import mloxygen.DispersedMintun1984Model.ensureModelPrereq
            
            blurTag = mlraichle.StudyRegistry.instance.blurTag;
            fs = sesd.fsOnAtlas('typ', 'mlfourd.ImagingContext2', ...
                'dateonly', true, 'tags', [blurTag sesd.regionTag]);
            ensureModelPrereq(sesd, fs)
            fs_avgxyz = fs.volumeAveraged(roi);  
            cbv = sesd.cbvOnAtlas('typ', 'mlfourd.ImagingContext2', ...
                'dateonly', true, 'tags', [blurTag sesd.regionTag]);
            ensureModelPrereq(sesd, cbv)
            v1_avgxyz = cbv.volumeAveraged(roi) * 0.0105;
            vecT = [fs_avgxyz.fourdfp.img v1_avgxyz.fourdfp.img];
        end
        function loss = loss_function(ks, fs, artery_interpolated, times_sampled, measurement, timeCliff)
            import mloxygen.DispersedMintun1984Model.sampled  
            estimation  = sampled(ks, fs, artery_interpolated, times_sampled);
            measurement = measurement(1:length(estimation));
            positive    = measurement > 0.05*max(measurement); % & times_sampled < timeCliff;
            eoverm      = estimation(positive)./measurement(positive);
            Q           = mean(abs(1 - eoverm));
            %Q           = sum((1 - eoverm).^2);
            loss        = Q; % 0.5*Q/sigma0^2 + sum(log(sigma0*measurement)); % sigma ~ sigma0*measurement
        end  
        function m    = preferredMap()
            %% init from Mintun J Nucl Med 25:177-187, 198.
            %  metabf described in Fig. 7.
            
            m = containers.Map;
            m('k1') = struct('min', 0.26, 'max', 0.62, 'init', 0.44,  'sigma', 0.01); % oef +/- 3 sd
            m('k2') = struct('min', eps,  'max', 1.5,  'init', 0.5,   'sigma', 0.01); % activity(HO(end))/activity(HO(90))
            m('k3') = struct('min', 0.2,  'max', 0.8,  'init', 0.5,   'sigma', 0.1); % activity(HO)/(activity(HO) + activity(OO)) at 90 sec
            m('k4') = struct('min', 0.25, 'max', 3,    'init', 0.835, 'sigma', 0.1); % v_post + 0.5 v_cap
        end
        function qs   = sampled(ks, fs, artery_interpolated, times_sampled)
            %  @param artery_interpolated is uniformly sampled at high sampling freq.
            %  @param times_sampled are samples scheduled by the time-resolved PET reconstruction
            
            import mloxygen.DispersedMintun1984Model.solution
            import mlpet.TracerKineticsModel.solutionOnScannerFrames 
            qs = solution(ks, fs, artery_interpolated);
            qs = solutionOnScannerFrames(qs, times_sampled);
        end
        function y    = sigmoid(x)
            %if x >= 0
            %    z = exp(-x);
            %    y = 1./(1 + z);
            %else
                z = exp(x);
                y = z./(1 + z);
            %end            
        end
        function loss = simulanneal_objective(ks, fs, artery_interpolated, times_sampled, qs0, sigma0)
            import mloxygen.DispersedMintun1984Model.sampled          
            qs = sampled(ks, fs, artery_interpolated, times_sampled);            
            loss = 0.5 * sum((1 - qs ./ qs0).^2) / sigma0^2; % + sum(log(sigma0*qs0)); % sigma ~ sigma0 * qs0
        end 
        function rho   = solution(ks, fs, artery_interpolated)
            %  @param artery_interpolated is uniform with high sampling freq. starting at time = 0.

            import mlpet.AerobicGlycolysisKit
            import mloxygen.DispersedMintun1984Model.sigmoid
            import mlpet.TracerKinetics
            
            RR = mlraichle.RaichleRegistry.instance();
            tBuffer = RR.tBuffer;
            ALPHA = 0.005670305; % log(2)/halflife in 1/s
            [~,idx0] = max(artery_interpolated > 0.05*max(artery_interpolated));
            idxU = idx0 + 90; % cf. Mintun1984
            
            oef = ks(1);
            metabTail = ks(2); 
            metabFrac = ks(3); 
            v_post_cap = ks(4);
            f = fs(1);
            PS = fs(2);
            lambda = fs(3); 
            v1 = fs(6);
            m = 1 - exp(-PS/f);
            n = length(artery_interpolated);
            times = 0:1:n-1;
            timesb = times; % - tBuffer;

            %% estimate shape of water of metabolism
            shape = zeros(1, n);
            shape(idx0:idxU) = linspace(0, 1, 91); % max(shape) == 1 
            shape(idxU+1:end) = linspace(1, metabTail, n-idxU); % min(shape) == metab2
            ductimes = zeros(1,n);
            ductimes(idx0+1:n) = 0:(n-idx0-1);
            ducshape = shape .* 2.^(-ductimes/122.2416); % decay-uncorrected
            
            %% set scale of artery_h2o
            metabScale = metabFrac*artery_interpolated(idxU); % activity water of metab \approx activity of oxygen after 90 sec
            metabScale = metabScale*TracerKinetics.DENSITY_PLASMA/TracerKinetics.DENSITY_BLOOD;
            artery_h2o = metabScale*ducshape;                     
            
            %% compartment 2, using m, f, lambda
            artery_o2 = artery_interpolated - artery_h2o;
            artery_o2(artery_o2 < 0) = 0;   
            kernel = exp(-m*f*timesb/lambda - ALPHA*timesb);
            rho2 =  m*f*conv(kernel, artery_h2o) + ...
                oef*m*f*conv(kernel, artery_o2);
            
            %% compartment 1
            % v_post = 0.83*v1;
            % v_cap = 0.01*v1;
            R = 0.85; % ratio of small-vessel to large-vessel Hct
            rho1 = v1*R*(1 - oef*v_post_cap)*artery_o2;
            
            rho = rho1(1:n) + rho2(1:n);        
            rho = rho(tBuffer+1:n);
        end        
    end

	methods		  
 		function this = DispersedMintun1984Model(varargin)
 			this = this@mloxygen.Mintun1984Model(varargin{:});
        end        
        function ho   = simulated(this, varargin)
            %% SIMULATED simulates tissue activity with passed and internal parameters.
            %  @param required ks is variational [metabf oef].
            %  @param required fs is deterministic [f PS lambda Delta Dt].
            %  @param aif is numeric; default is this.artery_interpolated for model state.
        
            ip = inputParser;
            addRequired(ip, 'ks', @isnumeric)
            addRequired(ip, 'fs', @isnumeric)
            addParameter(ip, 'aif', this.artery_interpolated, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;            
            
            Dt = ipr.fs(5);
            if Dt ~= 0
                times = 0:length(ipr.aif)-1;
                aif = makima(times - Dt, ipr.aif, times); % remove the delay Dt found by model
            else
                aif = ipr.aif;
            end
            ho = mloxygen.DispersedMintun1984Model.sampled(ipr.ks, ipr.fs, aif, this.times_sampled);
        end
    end 
    
    %% HIDDEN
    
    methods (Hidden)
        function rho   = solution_backup(ks, fs, artery_interpolated)
            %  @param artery_interpolated is uniform with high sampling freq. starting at time = 0.

            import mlpet.AerobicGlycolysisKit
            import mloxygen.DispersedMintun1984Model.sigmoid
            
            RR = mlraichle.RaichleRegistry.instance();
            tBuffer = RR.tBuffer;
            ALPHA = 0.005670305; % log(2)/halflife in 1/s
            [~,idx0] = max(artery_interpolated > 0.05*max(artery_interpolated));
            idxU = idx0 + 90; % cf. Mintun1984
            
            oef = ks(1);
            metab2 = ks(2); 
            metab3 = ks(3); 
            Delta = ks(4);
            f = fs(1);
            PS = fs(2);
            lambda = fs(3); 
            v1 = fs(6);
            m = 1 - exp(-PS/f);
            n = length(artery_interpolated);
            times = 0:1:n-1;
            timesb = times; % - tBuffer;
             
            % use Delta, metabf
            auc0 = trapz(artery_interpolated);
            artery_interpolated1 = conv(artery_interpolated, exp(-Delta*times));
            artery_interpolated1 = artery_interpolated1(1:n);
            artery_interpolated1 = artery_interpolated1*auc0/trapz(artery_interpolated1);

            shape = zeros(1, n);
            %stimes = times0(1:n-idx0) + 1e-6;
            %shape(idx0+1:n) = stimes.^(metab2-1) .* exp(-metab3*stimes);
            stimes1 = times(1:idxU-idx0) - times(idxU-idx0); % stimes1 := -90:0
            x = metab2*stimes1; % x := -4.5:0
            shape(idx0+1:idxU) = (sigmoid(x)      - sigmoid(x(1)))/ ...
                                 (sigmoid(x(end)) - sigmoid(x(1)));
            
            stimes2 = times(1:n-idxU);
            shape(idxU+1:n) = exp(-metab3*stimes2);
            
            ductimes = zeros(1,n);
            ductimes(idx0+1:n) = 0:(n-idx0-1);
            ducshape = shape .* 2.^(-ductimes/122.2416); % decay-uncorrected
            
            %metabScale = 0.1*artery_interpolated1(idxU)/ducshape(idxU);
            metabScale = 0.06*trapz(artery_interpolated1(1:idxU)) / ...
                              trapz(ducshape(1:idxU)); % water of metab is 6% of cumulative activity
            artery_h2o = metabScale*ducshape;
            artery_o2 = artery_interpolated1 - artery_h2o;
            artery_o2(artery_o2 < 0) = 0;            
            
            % compartment 2
            kernel = exp(-m*f*timesb/lambda - ALPHA*timesb);
            rho2 =  m*f*conv(kernel, artery_h2o) + ...
                oef*m*f*conv(kernel, artery_o2);
            
            % compartment 1
            % v_post = 0.83*v1;
            % v_cap = 0.01*v1;
            R = 0.85; % ratio of small-vessel to large-vessel Hct
            rho1 = v1*R*(1 - oef*0.835)*artery_o2;
            
            % use E, f, lambda
            rho = rho1(1:n) + rho2(1:n);
        
            rho = rho(tBuffer+1:n);
        end  
    end
    

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

