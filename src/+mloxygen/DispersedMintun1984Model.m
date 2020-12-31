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
            
            fs = sesd.fsOnAtlas('typ', 'mlfourd.ImagingContext2', ...
                'dateonly', true, 'tags', [sesd.petPointSpreadTag sesd.regionTag]);
            ensureModelPrereq(sesd, fs)
            fs_avgxyz = fs.volumeAveraged(roi);  
            cbv = sesd.cbvOnAtlas('typ', 'mlfourd.ImagingContext2', ...
                'dateonly', true, 'tags', [sesd.petPointSpreadTag sesd.regionTag]);
            ensureModelPrereq(sesd, cbv)
            v1_avgxyz = cbv.volumeAveraged(roi) * 0.0105;
            vecT = [fs_avgxyz.fourdfp.img v1_avgxyz.fourdfp.img];
        end
        function loss = loss_function(ks, fs, artery_interpolated, times_sampled, measurement, sigma0)
            import mloxygen.DispersedMintun1984Model.sampled  
            estimation  = sampled(ks, fs, artery_interpolated, times_sampled);
            measurement = measurement(1:length(estimation));
            taus        = diff(times_sampled);
            taus        = [taus taus(end)];
            taus        = taus(1:length(estimation));
            positive    = measurement > 0;
            e           = estimation .* taus;
            m           = measurement .* taus;
            eoverm      = e(positive)./m(positive);
            Q           = sum((1 - eoverm).^2);
            loss        = 0.5*Q/sigma0^2; % + sum(log(sigma0*measurement)); % sigma ~ sigma0*measurement
        end  
        function m    = preferredMap()
            %% init from Mintun J Nucl Med 25:177-187, 198.
            %  metabf described in Fig. 7.
            
            m = containers.Map;
            m('k1') = struct('min', 0, 'max', 1, 'init', 0.5,  'sigma', 0.01); % unit-less metabf
            m('k2') = struct('min', 0, 'max', 1, 'init', 0.44, 'sigma', 0.01); % unit-less oef
        end
        function qs   = sampled(ks, fs, artery_interpolated, times_sampled)
            %  @param artery_interpolated is uniformly sampled at high sampling freq.
            %  @param times_sampled are samples scheduled by the time-resolved PET reconstruction
            
            import mloxygen.DispersedMintun1984Model.solution
            qs = solution(ks, fs, artery_interpolated);
            n = length(artery_interpolated);
            qs = makima(0:n-1, qs, times_sampled);
        end
        function loss = simulanneal_objective(ks, fs, artery_interpolated, times_sampled, qs0, sigma0)
            import mloxygen.DispersedMintun1984Model.sampled          
            qs = sampled(ks, fs, artery_interpolated, times_sampled);            
            loss = 0.5 * sum((1 - qs ./ qs0).^2) / sigma0^2; % + sum(log(sigma0*qs0)); % sigma ~ sigma0 * qs0
        end 
        function rho   = solution(ks, fs, artery_interpolated)
            %  @param artery_interpolated is uniform with high sampling freq. starting at time = 0.

            import mlpet.AerobicGlycolysisKit
            
            ALPHA = 0.005670305; % log(2)/halflife in 1/s
            
            metabf = ks(1);
            oef = ks(2);
            f = fs(1);
            PS = fs(2);
            lambda = fs(3); 
            Delta = fs(4);
            Dt = round(fs(5));
            v1 = fs(6);
            m = max(1 - exp(-PS/f), AerobicGlycolysisKit.E_MIN);
            m = min(m, AerobicGlycolysisKit.E_MAX);
            n = length(artery_interpolated);
            times = 0:1:n-1;
             
            % use Delta, metabf
            auc0 = trapz(artery_interpolated);
            artery_interpolated1 = conv(artery_interpolated, exp(-Delta*times));
            artery_interpolated1 = artery_interpolated1(1:n);
            artery_interpolated1 = artery_interpolated1*auc0/trapz(artery_interpolated1);
            
            shape = ones(size(artery_interpolated1));
            if Dt < 0
                shape(1:91) = linspace(0, metabf, 91);
                shape(92:end) = linspace(metabf, 1, length(shape)-91);
                artery_h2o_dc = metabf*artery_interpolated1(91)*2^(90/122.2416)*shape;
            else
                shape(1:1+Dt) = 0;
                shape(1+Dt:91+Dt) = linspace(0, metabf, 91);
                shape(92+Dt:end) = linspace(metabf, 1, length(shape)-91-Dt);
                artery_h2o_dc = metabf*artery_interpolated1(91+Dt)*2^((90+Dt)/122.2416)*shape;
            end
            artery_h2o = artery_h2o_dc .* 2.^(-times/122.2416);
            artery_o2 = artery_interpolated1 - artery_h2o;
            artery_o2(artery_o2 < 0) = 0;
            
            % compartment 2
            rho2 =  m*f*conv(exp(-m*f*times/lambda - ALPHA*times), artery_h2o) + ...
                oef*m*f*conv(exp(-m*f*times/lambda - ALPHA*times), artery_o2);
            
            % compartment 1
            % v_post = 0.83*v1;
            % v_cap = 0.01*v1;
            R = 0.85; % ratio of small-vessel to large-vessel Hct
            rho1 = v1*R*(1 - oef*0.835)*artery_o2;
            
            % use E, f, lambda
            rho = rho1(1:n) + rho2(1:n);
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
                aif = pchip(times - Dt, ipr.aif, times); % remove the delay Dt found by model
            else
                aif = ipr.aif;
            end
            ho = mloxygen.DispersedMintun1984Model.sampled(ipr.ks, ipr.fs, aif, this.times_sampled);
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

