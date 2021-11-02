classdef FieldNumericRaichle1983 < handle & mloxygen.FieldNumeric
	%% FIELDNUMERICRAICHLE1983  

	%  $Revision$
 	%  was created 23-Aug-2021 20:06:16 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.10.0.1710957 (R2021a) Update 4 for MACI64.  Copyright 2021 John Joowon Lee.
 	
    methods (Static)
        function this = createFromDeviceKit(devkit, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param required scanner is an mlpet.AbstractDevice.
            %  @param required arterial is an mlpet.AbstractDevice.
            %  @param roi is mlfourd.ImagingContext2.
            %  @return this.
            
            import mloxygen.FieldNumericRaichle1983
            import mloxygen.FieldNumeric.mixTacAif       
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'scanner', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'arterial', [], @(x) isa(x, 'mlpet.AbstractDevice')) 
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            parse(ip, devkit, varargin{:})
            
            % mix components for augmentation       
            [tacs_,timesMid_,t0_] = prepareTacs(devkit, varargin{:});
            
            %
            
            fp = sprintf('mlglucose_FieldNumericRaichle1983_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS')); 
            
            this = FieldNumericRaichle1983( ...
                'ho', tacs_, ...
                'devkit', devkit, ...
                'timesMid', timesMid_, ...
                't0', t0_, ...
                'tObs', 60, ...
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
            
            fs_ = copy(this.roi.fourdfp);
            fs_.img = single(this.img_);
            fs_.fileprefix = this.sessionData.fsOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
            fs_ = imagingType(ipr.typ, fs_);
        end 
        function this = solve(this, varargin)
        end
    end
    
    %% PROTECTED

	methods (Access = protected)		  
 		function this = FieldNumericRaichle1983(varargin)
 			%% FIELDNUMERICRAICHLE1983
 			%  @param .

 			this = this@mloxygen.FieldNumeric(varargin{:});
            
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

