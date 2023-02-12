classdef QuadraticNumericRaichle1983 < handle & mloxygen.QuadraticNumeric
	%% QUADRATICNUMERICRAICHLE1983  

	%  $Revision$
 	%  was created 12-Jun-2021 21:12:36 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.10.0.1669831 (R2021a) Update 2 for MACI64.  Copyright 2021 John Joowon Lee.
    
    methods (Static)
        function this = createFromDeviceKit(devkit, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param required scanner is an mlpet.AbstractDevice.
            %  @param required arterial is an mlpet.AbstractDevice.
            %  @param roi is mlfourd.ImagingContext2.
            %  @return this.
            
            import mloxygen.QuadraticNumericRaichle1983
            import mloxygen.QuadraticNumeric.mixTacAif       
            
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
            
            fp = sprintf('mlglucose_QuadraticNumericRaichle1983_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS')); 
            
            this = QuadraticNumericRaichle1983( ...
                'ho', tac_, ...
                'devkit', devkit, ...
                'timesMid', timesMid_, ...
                't0', t0_, ...
                'tObs', 60, ...
                'artery_interpolated', aif_, ...
                'fileprefix', fp, ...
                varargin{:});      
        end
    end
    
	methods  
        function fs_ = fs(this, varargin)
            %% @param typ is forwarded to imagingType(), e.g., 'single', 'ImagingContext2', ...  
            %  @returns estimated f map in Hz.
            
            ip = inputParser;
            addParameter(ip, 'typ', 'mlfourd.ImagingContext2', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            fs_ = copy(this.roi.imagingFormat);
            fs_.img = single(this.img_);
            fs_.fileprefix = this.fsOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
            fs_ = imagingType(ipr.typ, fs_);
        end 
        function fp = fsOnAtlas(this, varargin)
            if isa (this.sessionData, 'mlnipet.SessionData')
                fp = this.sessionData.fsOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
                return
            end
            if isa (this.sessionData, 'mlpipeline.ImagingMediator')
                tags = strip([this.blurTag this.regionTag], '_');
                ic = this.sessionData.metricOnAtlas('fs', tags);
                fp = ic.fileprefix;
                return
            end
            error('mloxygen:RuntimeError', stackstr())
        end
        function this = solve(this, varargin)
            obsAif = this.obsFromAif(this.artery_interpolated, this.canonical_f);
            this.modelA = this.buildQuadraticModel(obsAif, this.canonical_cbf);
            obsPet = this.obsFromTac(this.measurement);
            this.img_ = this.a1*obsPet.^2 + this.a2*obsPet;
            this.img_ = this.cbfToF1(this.img_);
        end
    end
    
    %% PROTECTED
    
	methods (Access = protected)
 		function this = QuadraticNumericRaichle1983(varargin)
 			%% QUADRATICNUMERICRAICHLE1983
 			%  @param .

 			this = this@mloxygen.QuadraticNumeric(varargin{:});
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            ip.PartialMatching = false;
            addParameter(ip, 'ho', [], @(x) isnumeric(x))
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            this.measurement = ipr.ho;
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

