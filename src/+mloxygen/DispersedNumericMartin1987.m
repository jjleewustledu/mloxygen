classdef DispersedNumericMartin1987 < handle & mlpet.AugmentedData & mloxygen.Martin1987_20210210
	%% DISPERSEDNUMERICMARTIN1987

	%  $Revision$
 	%  was created 29-Apr-2020 23:31:38 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1319299 (R2019b) Update 5 for MACI64.  Copyright 2020 John Joowon Lee.

    
    properties (Constant)   
        LENK = 1
    end
    
    properties (Dependent)
        artery_sampled
        Delta
        Dt % time-shift for AIF; Dt < 0 shifts backwards in time.
    end
    
    methods (Static)
        function [this,tac_,aif_] = createFromDeviceKit(devkit, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param required scanner is an mlpet.AbstractDevice.
            %  @param required arterial is an mlpet.AbstractDevice.
            %  @param roi is mlfourd.ImagingContext2.
            %  @param T0 isscalar.
            %  @param Tf isscalar.
            %  @return this.
            %  @return tac_, blurred by ipr.blurFdg.
            %  @return aif_.
            
            import mloxygen.DispersedNumericMartin1987
            import mloxygen.DispersedNumericMartin1987.mixTacAif
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'scanner', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'arterial', [], @(x) isa(x, 'mlpet.AbstractDevice'))          
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'T0', 120, @isscalar)
            addParameter(ip, 'Tf', 240, @isscalar)
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
                         
            [tac_,timesMid_,aif_,Dt] = mixTacAif(devkit, varargin{:});
            fp = sprintf('mlglucose_DispersedNumericMartin1987_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS')); 
            
            this = DispersedNumericMartin1987( ...
                'oc', tac_, ...
                'solver', '', ...
                'devkit', devkit, ...
                'Dt', Dt, ...
                'times_sampled', timesMid_, ...
                'artery_interpolated', aif_, ...
                'fileprefix', fp, ...
                'T0', ipr.T0, ...
                'Tf', min([ipr.Tf max(timesMid_) length(aif_)-1]), ...
                'timeCliff', ipr.arterial.timeCliff, ...
                varargin{:});
        end
    end
    
    methods 
        
        %% GET
        
        function g = get.artery_sampled(this)
            g = this.model.solutionOnScannerFrames( ...
                this.artery_interpolated(this.tBuffer+1:end), this.times_sampled);
        end
        function g = get.Delta(~)
            g = 0;
        end
        function g = get.Dt(this)
            g = this.strategy_.Dt;
        end
        
        %%
        
        function a = artery_global(this, varargin)
            %% ARTERY_GLOBAL
            %  @required roi is understood by mlfourd.ImagingContext2.
            %  @param typ is understood by imagingType.
            %  @return a is an imagingType.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'roi')
            addParameter(ip, 'typ', 'mlfourd.ImagingContext2', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            ipr.roi = mlfourd.ImagingContext2(ipr.roi);
            ipr.roi = ipr.roi.binarized();
            
            artery_interpolated1 = this.artery_interpolated(this.tBuffer+1:end);
            avec = this.model.solutionOnScannerFrames(artery_interpolated1, this.times_sampled);
            
            roibin = logical(ipr.roi);
            a = copy(ipr.roi.fourdfp);
            a.img = zeros([size(ipr.roi) length(avec)]);
            for t = 1:length(avec)
                img = zeros(size(ipr.roi), 'single');
                img(roibin) = avec(t);
                a.img(:,:,:,t) = img;
            end
            a.fileprefix = this.sessionData.aifsOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
            a = imagingType(ipr.typ, a);
        end
		function a = artery_local(this, varargin)
            %% ARTERY_LOCAL
            %  @param typ is understood by imagingType.
            %  @return a is an imagingType.
            %  See also ml*.Dispersed*Model.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'typ', 'mlfourd.ImagingContext2', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            n = length(this.artery_interpolated);
            if this.Dt < 0 % shift back to right
                artery_interpolated1 = zeros(1, n);
                artery_interpolated1(-this.Dt+1:end) = this.artery_interpolated(1:n+this.Dt);
            elseif this.Dt > 0 % shift back to left
                artery_interpolated1 = this.artery_interpolated(1)*ones(1, n);
                artery_interpolated1(1:end-this.Dt) = this.artery_interpolated(1+this.Dt:n);
            else
                artery_interpolated1 = this.artery_interpolated;
            end 
            artery_interpolated1 = artery_interpolated1(this.tBuffer+1:end);
            avec = this.model.solutionOnScannerFrames(artery_interpolated1, this.times_sampled);
            
            roibin = logical(this.roi);
            a = copy(this.roi.fourdfp);
            a.img = zeros([size(this.roi) length(avec)]);
            for t = 1:length(avec)
                img = zeros(size(this.roi), 'single');
                img(roibin) = avec(t);
                a.img(:,:,:,t) = img;
            end
            a.fileprefix = this.sessionData.aifsOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
            a = imagingType(ipr.typ, a);
        end 
        function ks = buildKs(this, varargin)
            this = solve(this, varargin{:});
            ks = [k1(this) k2(this)];
        end 
        function [k,sk] = k2(this, varargin)
            k = this.Dt;
            sk = nan;
        end
        function ks_ = ks(this, varargin)
            %% ks == [v1 Dt]
            %  @param 'typ' is char, understood by imagingType.            
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'typ', 'mlfourd.ImagingContext2', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            k(1) = k1(this.strategy_, varargin{:});
            k(2) = k2(this);
             
            roibin = logical(this.roi);
            ks_ = copy(this.roi.fourdfp);
            ks_.img = zeros([size(this.roi) this.LENK]);
            for t = 1:length(k)
                img = zeros(size(this.roi), 'single');
                img(roibin) = k(t);
                ks_.img(:,:,:,t) = img;
            end
            ks_.fileprefix = this.sessionData.vsOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
            ks_ = imagingType(ipr.typ, ks_);
        end 
    end
    
    %% PROTECTED

	methods (Access = protected)		  
 		function this = DispersedNumericMartin1987(varargin)
 			%% DISPERSEDNUMERICMARTIN1987
            %  @param ho is numeric.
            %  @param devkit is mlpet.IDeviceKit.   
            %  
            %  for mloxygen.DispersedMartin1987Model: 
            %  @param times_sampled for scanner is typically not uniform.
            %  @param artery_interpolated must be uniformly interpolated.
 			%  @param T0 is numeric.
            %  @param Tf is numeric.
            %
            %  for mloxygen.DispersedMartin1987Solver:
            %  @param context is mloxygen.Raichle1983.
            %  @param sigma0.
            %  @param fileprefix.
            
            this = this@mloxygen.Martin1987_20210210( ...
                'model', mloxygen.DispersedMartin1987Model(varargin{:}), ...
                varargin{:});

            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'oc', [], @isnumeric)
            addParameter(ip, 'solver', '', @ischar)
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            this.measurement = ipr.oc;
            this.strategy_ = mloxygen.DispersedMartin1987Solver( ...
                'context', this, varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

