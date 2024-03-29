classdef QuadraticNumericMartin1987 < handle & mloxygen.QuadraticNumeric
	%% QUADRATICNUMERICMARTIN1987  

	%  $Revision$
 	%  was created 12-Jun-2021 21:12:51 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.10.0.1669831 (R2021a) Update 2 for MACI64.  Copyright 2021 John Joowon Lee.
 	    
    properties (Constant)
        T0_DELAY = 60
    end

    methods (Static)
        function this = createFromDeviceKit(devkit, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param required scanner is an mlpet.AbstractDevice.
            %  @param required arterial is an mlpet.AbstractDevice.
            %  @param roi is mlfourd.ImagingContext2.
            %  @return this.
            
            import mloxygen.QuadraticNumericMartin1987
            import mlkinetics.ScannerKit.mixTacAif       
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'scanner', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'arterial', [], @(x) isa(x, 'mlpet.AbstractDevice')) 
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            parse(ip, devkit, varargin{:})
            
            % mix components for augmentation       
            [tac_,timesMid_,t0_,aif_] = mixTacAif(devkit, varargin{:});
            
            %
            
            fp = sprintf('mlglucose_QuadraticNumericMartin1987_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS')); 
            
            this = QuadraticNumericMartin1987( ...
                'oc', tac_, ...
                'devkit', devkit, ...
                'timesMid', timesMid_, ...
                't0', t0_ + QuadraticNumericMartin1987.T0_DELAY, ...
                'tObs', 240, ...
                'artery_interpolated', aif_, ...
                'fileprefix', fp, ...
                varargin{:});      
        end
    end

	methods 
 		function this = QuadraticNumericMartin1987(varargin)
 			%% QUADRATICNUMERICMARTIN1987
 			%  @param .

 			this = this@mloxygen.QuadraticNumeric(varargin{:});
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            ip.PartialMatching = false;
            addParameter(ip, 'oc', [], @(x) isnumeric(x))
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            this.measurement = ipr.oc; 			
        end
        function vs_ = vs(this, varargin)
            %% @param typ is forwarded to imagingType(), e.g., 'single', 'ImagingContext2', ...  
            %  @returns estimated v map in [0 1].
            
            ip = inputParser;
            addParameter(ip, 'typ', 'mlfourd.ImagingContext2', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            vs_ = copy(this.roi.imagingFormat);
            vs_.img = single(this.img_);
            vs_.fileprefix = this.vsOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
            vs_ = imagingType(ipr.typ, vs_);
        end
        function fp = vsOnAtlas(this, varargin)
            if isa (this.sessionData, 'mlnipet.SessionData')
                fp = this.sessionData.vsOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
                return
            end
            if isa (this.sessionData, 'mlpipeline.ImagingMediator')

                tags = strip([this.blurTag this.regionTag], '_');
                ic = this.sessionData.metricOnAtlas('vs', tags);
                fp = ic.fileprefix;
                return
            end
            error('mloxygen:RuntimeError', stackstr())
        end
        function this = solve(this)
            tF_brink = min(this.tF, this.timeCliff);
            obsPet = this.obsFromTac(this.measurement, t0=this.t0, tF=tF_brink);
            integralAif = trapz(this.artery_interpolated(this.t0+1:tF_brink+1));
            this.img_ = obsPet/(this.RATIO_SMALL_LARGE_HCT*this.DENSITY_BRAIN*integralAif);
        end
    end		  
    
    %% PROTECTED
    
	methods (Access = protected)
        function tcliff = timeCliff(this)
            artery = this.artery_interpolated;
            [~,tmax] = max(artery);
            artery1 = artery(tmax:end);
            [~,tcliff] = min(artery1);
            tcliff = tmax + tcliff - 1;
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

