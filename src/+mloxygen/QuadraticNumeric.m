classdef (Abstract) QuadraticNumeric < handle & mlpet.TracerKinetics
	%% QUADRATICNUMERIC  
    %  1. Videen, T. O., Perlmutter, J. S., Herscovitch, P. & Raichle, M. E.
    %     Brain Blood Volume, Flow, and Oxygen Utilization Measured with 15 O Radiotracers and Positron Emission Tomography: Revised Metabolic Computations. 
    %     J Cereb Blood Flow Metab 7, 513–516 (1987).
    %  2. Herscovitch, P., Mintun, M. A. & Raichle, M. E. 
    %     Brain Oxygen Utilization Measured with Oxygen-15 Radiotracers and Positron Emission Tomography: Generation of Metabolic Images. 
    %     J Nucl Med 26, 416–417 (1985).

	%  $Revision$
 	%  was created 12-Jun-2021 21:31:16 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.10.0.1669831 (R2021a) Update 2 for MACI64.  Copyright 2021 John Joowon Lee.
 	
    methods (Abstract)
        this = solve(this)
    end
    
	properties (Dependent)
        a
        a1
        a2
        b
        b1
        b2
        b3
        b4
 		alpha_decay
        blurTag
        canonical_f
        regionTag
        roi % mlfourd.ImagingContext2
        tF
    end
    
    properties
        artery_interpolated
        canonical_cbf = 8:2:110 % mL/hg/min
                                 % not to exceed 110 per Cook's distance from buildQuadraticModel()
        measurement
        modelA
        modelB12
        modelB34
        t0
        taus
        timesMid
        tObs
    end
    
    methods (Static) 
        function [tac__,timesMid__,t0__,aif__,Dt,datetimePeak] = mixTacAif(devkit, varargin)
            %  @return taus__ are the durations of emission frames
            %  @return t0__ is the start of the first emission frame.
            %  @return aif__ is aligned in time to emissions.
            %  @return Dt is the time shift needed to align aif to emissions.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'scanner', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'arterial', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            if strcmp(class(ipr.scanner), class(ipr.arterial))
                [tac__,timesMid__,t0__,aif__,Dt,datetimePeak] = ...
                    mloxygen.QuadraticNumeric.mixTacIdif(devkit, varargin{:});
                return
            end
            ipr.roi = ipr.roi.binarized();
            ad = mlaif.AifData.instance();
            
            % scannerDevs provide calibrations & ROI-masking    
            %blur = ipr.devkit.sessionData.petPointSpread;
            s = ipr.scanner; %.blurred(blur);
            s = s.masked(ipr.roi);
            tac = s.activityDensity();
            tac(tac < 0) = 0;    
            tac = ad.normalizationFactor*tac; % empirical normalization                   
            tac__ = tac;
            taus__ = s.taus;
            timesMid__ = s.timesMid;
            Nt = ceil(timesMid__(end));
            
            % estimate t0__
            tac_avgxyz = squeeze(mean(mean(mean(tac__, 1), 2), 3));
            dtac_avgxyz = diff(tac_avgxyz);
            [~,idx] = max(dtac_avgxyz > 0.05*max(dtac_avgxyz));
            t0__ = timesMid__(idx) - taus__(idx)/2;
            
            % arterialDevs calibrate & align arterial times-series to localized scanner time-series            
            a0 = ipr.arterial;
            [a, datetimePeak] = devkit.alignArterialToScanner( ...
                a0, s, ipr.aifdata, 'sameWorldline', false);
            aif = a.activityDensity('Nt', Nt);
            switch class(a)
                case 'mlswisstrace.TwiliteDevice'
                    t = (0:Nt-1) - seconds(s.datetime0 - a.datetime0);
                case 'mlcapintec.CapracDevice'
                    t = a.times - seconds(s.datetime0 - a.datetime0);
                otherwise
                    error('mloxygen:ValueError', ...
                        'class(AugmentedData.mixTacAif.a) = %s', class(a))
            end
            
            % adjust aif__, get Dt
            if min(t) > 0
                aif = interp1([0 t], [0 aif], 0:s.timesMid(end), 'linear', 0);
            else                
                aif = interp1(t, aif, 0:s.timesMid(end), 'linear', 0);
            end
            aif(aif < 0) = 0;            
            aif__ = aif;            
            Dt = a.Dt;
        end
        function [tac__,timesMid__,t0__,idif__,Dt,datetimePeak] = mixTacIdif(devkit, varargin)
            %  @return taus__ are the durations of emission frames
            %  @return t0__ is the start of the first emission frame.
            %  @return idif__ is aligned in time to emissions.
            %  @return Dt is the time shift needed to align aif to emissions.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'scanner', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'arterial', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            ipr.roi = ipr.roi.binarized();
            
            % scannerDevs provide calibrations & ROI-masking    
            %blur = ipr.devkit.sessionData.petPointSpread;
            s = ipr.scanner; %.blurred(blur);
            s = s.masked(ipr.roi);
            tac = s.activityDensity();
            tac(tac < 0) = 0;                       
            tac__ = tac;
            taus__ = s.taus;
            timesMid__ = s.timesMid;
            
            % estimate t0__
            tac_avgxyz = squeeze(mean(mean(mean(tac__, 1), 2), 3));
            dtac_avgxyz = diff(tac_avgxyz);
            [~,idx] = max(dtac_avgxyz > 0.05*max(dtac_avgxyz));
            t0__ = timesMid__(idx) - taus__(idx)/2;
            
            % idif
            idif = a.activityDensity();            
            t = a.timesMid;
            
            % adjust aif__, get Dt
            idif = interp1([0 t], [0 idif], 0:s.timesMid(end), 'linear', 0);
            idif(idif < 0) = 0;            
            idif__ = idif; 

            % trivial values
            Dt = 0;
            datetimePeak = NaT;
        end
    end

	methods
        
        %% GET
        
        function g = get.a(this)
            g = [this.a1 this.a2];
        end
        function g = get.a1(this)            
            g = this.modelA.Coefficients{1, 'Estimate'};
        end
        function g = get.a2(this)            
            g = this.modelA.Coefficients{2, 'Estimate'};
        end
        function g = get.b(this)
            g = [this.b1 this.b2 this.b3 this.b4];
        end
        function g = get.b1(this)            
            g = this.modelB12.Coefficients{1, 'Estimate'};
        end
        function g = get.b2(this)            
            g = this.modelB12.Coefficients{2, 'Estimate'};
        end
        function g = get.b3(this)            
            g = this.modelB34.Coefficients{1, 'Estimate'};
        end
        function g = get.b4(this)            
            g = this.modelB34.Coefficients{2, 'Estimate'};
        end
        function g = get.alpha_decay(~)
            g = mlpet.Radionuclides.decayConstantOf('15O');
        end
        function g = get.blurTag(~)
            g = mlraichle.StudyRegistry.instance.blurTag;
        end
        function g = get.canonical_f(this)
            g = this.cbfToF1(this.canonical_cbf);
        end
        function g = get.regionTag(this)
            g = this.sessionData.regionTag;
        end
        function g = get.roi(this)
            g = this.roi_;
        end    
        function     set.roi(this, s)
            if ~isempty(s)
                this.roi_ = mlfourd.ImagingContext2(s);
                this.roi_ = this.roi_.binarized();
            end
        end
        function g = get.tF(this)
            g = this.t0 + this.tObs;
        end
        
        %%
        
        function obs = obsFromAif(this, varargin)
            %% Cf. Videen 1987, Eq. 2.
            %  @param required aif is a vector containing arterial input, sampled at 1 Hz.
            %  @param required cbf is the vector of flows to model, expressed as Hz.
            %  @return obs is the RHS of Eq. 2, a scalar \int_0^T dt \text{tracer density}.
            
            ip = inputParser;
            addRequired(ip, 'aif', @isnumeric)
            addRequired(ip, 'f', @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;            
            ipr.aif = asrow(ipr.aif);
            ipr.f = ascolumn(ipr.f);
            
            lam = this.LAMBDA;
            alph = this.alpha_decay;  
            N = length(ipr.aif);
            M = length(ipr.f);
            times = 0:(N-1);
            rho = zeros(M, this.tF-this.t0+1);
            for r = 1:M
                rho_ = ipr.f(r)*conv(ipr.aif, exp(-(ipr.f(r)/lam + alph)*times));
                rho(r,:) = rho_(this.t0+1:this.tF+1);
            end

            obs = trapz(rho, 2); % \int dt
        end
        function obs = obsFromTac(this, varargin)
            %% Cf. Videen 1987, Eq. 2.
            %  @param required tac is numeric, containing scanner TAC in R^(3+1) with native time-frames.
            %  @param optional taus is the vector of sampling durations for native time-frames.
            %  @param optional t0 is the start time of observation in seconds.
            %  @param optional tF is the finish time of observation in seconds.
            %  @return obs is the RHS of Eq. 2, in R^3 := \int_0^T dt \text{tracer density}.            
            
            ip = inputParser;
            addRequired(ip, 'tac', @isnumeric)
            addParameter(ip, 'timesMid', this.timesMid, @isnumeric)
            addParameter(ip, 't0', this.t0, @isscalar)
            addParameter(ip, 'tF', this.tF, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;            
            assert(size(ipr.tac,4) == length(ipr.timesMid))
            
            window = ipr.t0 <= ipr.timesMid & ipr.timesMid <= ipr.tF;            
            obs = trapz(ipr.timesMid(window), ipr.tac(:,:,:,window), 4);
        end
    end
    
    %% PROTECTED
    
    properties (Access = protected)
        img_
        roi_
    end
    
    methods (Static, Access = protected)        
        function cbf = cbfQuadraticModel(As, obs)
            %% @return cbf in mL/hg/min.
            
            cbf = obs.^2*As(1) + obs*As(2);
        end
    end
    
    methods (Access = protected)	
        function mdl = buildQuadraticModel(~, obs, cbf)
            %% BUILDQUADRATICMODEL 
            %  @param obs are numeric PET_{obs} := \int_{t \in \text{obs}} dt' \varrho(t').
            %  @param cbf are numeric CBF or similar flows in mL/min/hg.
            %  @returns mdl.  A1, A2 are in mdl.Coefficients{:,'Estimate'}.
            %  See also:  https://www.mathworks.com/help/releases/R2016b/stats/nonlinear-regression-workflow.html
            
            fprintf('QuadraticNumeric.buildQuadraticModel ..........\n');
            mdl = fitnlm( ...
                ascolumn(obs), ...
                ascolumn(cbf), ...
                @mloxygen.QuadraticNumeric.cbfQuadraticModel, ...
                [1 1]);
            disp(mdl)
            fprintf('mdl.RMSE -> %g, min(rho) -> %g, max(rho) -> %g\n', mdl.RMSE, min(obs), max(obs));
            if isempty(getenv('NOPLOT')) || ~isempty(getenv('DEBUG'))
                plotResiduals(mdl);
                plotDiagnostics(mdl, 'cookd');
                plotSlice(mdl);
            end
        end
        function savefig(this, varargin)
            ip = inputParser;
            addRequired(ip, 'handle', @ishandle) % fig handle
            addParameter(ip, 'tags', '', @istext) % for filenames
            parse(ip, varargin{:})
            ipr = ip.Results;            
            
            if ~isempty(ipr.tags)
                ipr.tags = strcat("_", strip(ipr.tags, "_"));
            end          
            dbs = dbstack;
            client = dbs(2).name;
            client_ = strrep(dbs(2).name, '.', '_');
            dtStr = datestr(this.sessionData.datetime);
            title(sprintf('%s\n%s %s', client, ipr.tags, dtStr))
            try
                dtTag = lower(this.sessionData.doseAdminDatetimeTag);
                savefig(ipr.handle, ...
                    fullfile(this.dataPath, ...
                    sprintf('%s%s_%s.fig', client_, tags_, dtTag)))
                figs = get(0, 'children');
                saveas(figs(1), ...
                    fullfile(this.dataPath, ...
                    sprintf('%si%s_%s.png', client_, tags_, dtTag)))
                close(figs(1))
            catch ME
                handwarning(ME)
            end
        end
        
 		function this = QuadraticNumeric(varargin)
 			%% QUADRATICNUMERIC
            %  @param devkit is mlpet.IDeviceKit.
            %  @param sessionData is mlpipeline.ISessionData, but defers to devkit.sessionData.
            %  @param roi is understood by mlfourd.ImagingContext2; will be binarized.

 			this = this@mlpet.TracerKinetics(varargin{:});
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            ip.PartialMatching = false;
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'artery_interpolated', [], @isvector)
            addParameter(ip, 'timesMid', [], @isvector)
            addParameter(ip, 't0', 0, @isscalar)
            addParameter(ip, 'tObs', 40, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            if ~isempty(ipr.roi)
                this.roi_ = mlfourd.ImagingContext2(ipr.roi);
                this.roi_ = this.roi_.binarized();
            end
            this.artery_interpolated = ipr.artery_interpolated;
            this.timesMid = ipr.timesMid;
            this.t0 = ipr.t0;
            this.tObs = ipr.tObs;
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

