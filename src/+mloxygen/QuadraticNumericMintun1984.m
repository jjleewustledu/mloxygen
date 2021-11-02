classdef QuadraticNumericMintun1984 < handle & mloxygen.QuadraticNumeric
	%% QUADRATICNUMERICMINTUN1984  

	%  $Revision$
 	%  was created 12-Jun-2021 21:13:11 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.10.0.1669831 (R2021a) Update 2 for MACI64.  Copyright 2021 John Joowon Lee.
 	
    properties (Constant)
        NOMINAL_O2_CONTENT = 0.17 % mL O2.  Estimated from Ito, Eur J Nucl Med Mol Imaging (2004) 31:635-643.
                                  % DOI 10.1007/s00259-003-1430-8
    end
    
    properties (Dependent)
        artery_oxygen
        artery_water_metab
        cbf
        cbv
        integral_artery_oxygen
    end
    
    methods (Static)
        function this = createFromDeviceKit(devkit, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param required scanner is an mlpet.AbstractDevice.
            %  @param required arterial is an mlpet.AbstractDevice.
            %  @param roi is mlfourd.ImagingContext2.
            %  @return this.
            
            import mloxygen.QuadraticNumericMintun1984
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
            
            fp = sprintf('mlglucose_QuadraticNumericMintun1984_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS')); 
            
            this = QuadraticNumericMintun1984( ...
                'oo', tac_, ...
                'devkit', devkit, ...
                'timesMid', timesMid_, ...
                't0', t0_, ...
                'tObs', 90, ...
                'artery_interpolated', aif_, ...
                'fileprefix', fp, ...
                varargin{:}); 
            sd = this.sessionData;
            this.cbf_ = sd.cbfOnAtlas('typ', 'ImagingContext2', 'dateonly', true, 'tags', sd.regionTag);
            this.cbv_ = sd.cbvOnAtlas('typ', 'ImagingContext2', 'dateonly', true, 'tags', sd.regionTag);
            this = metabolize(this);
        end
    end

	methods 
        
        %% GET
        
        function g = get.artery_oxygen(this)
            g = this.artery_oxygen_;
            assert(~isempty(g));
        end
        function g = get.artery_water_metab(this)
            g = this.artery_water_metab_;
            assert(~isempty(g));
        end
        function g = get.cbf(this)
            g = this.cbf_;
            assert(~isempty(g));            
        end
        function g = get.cbv(this)
            g = this.cbv_;
            assert(~isempty(g));            
        end
        function g = get.integral_artery_oxygen(this)
            g = this.integral_artery_oxygen_;
            assert(~isempty(g));
        end
        
        %%
        
        function this = metabolize(this)
            [~,idx0] = max(this.artery_interpolated > 0.05*max(this.artery_interpolated));
            idxU = idx0 + 90;
            metabTail = 0.5; % activity(HO(end))/activity(HO(90))
            metabFrac = 0.5; % activity(HO)/(activity(HO) + activity(OO)) at 90 sec
            n = length(this.artery_interpolated);
            
            %% estimate shape of water of metabolism
            shape = zeros(1, n);
            shape(idx0:idxU) = linspace(0, 1, 91); % max(shape) == 1 
            shape(idxU+1:end) = linspace(1, metabTail, n-idxU); % min(shape) == metab2
            ductimes = zeros(1,n);
            ductimes(idx0+1:n) = 0:(n-idx0-1);
            ducshape = shape .* 2.^(-ductimes/122.2416); % decay-uncorrected
            
            %% set scale of artery_h2o
            metabScale = metabFrac*this.artery_interpolated(idxU); % activity water of metab \approx activity of oxygen after 90 sec
            metabScale = metabScale*this.DENSITY_PLASMA/this.DENSITY_BLOOD;
            
            %% set internal params
            this.artery_water_metab_ = metabScale*ducshape;            
            this.artery_oxygen_ = this.artery_interpolated - this.artery_water_metab_;
            this.integral_artery_oxygen_ = ...
                0.01*this.RATIO_SMALL_LARGE_HCT*this.DENSITY_BRAIN* ...
                trapz(this.artery_oxygen_(this.t0+1:this.tF+1));            
        end
        function os_ = os(this, varargin)
            %% @param typ is forwarded to imagingType(), e.g., 'single', 'ImagingContext2', ...  
            %  @returns estimated f map in Hz.
            
            ip = inputParser;
            addParameter(ip, 'typ', 'mlfourd.ImagingContext2', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            os_ = copy(this.roi.fourdfp);
            os_.img = single(this.img_);
            os_.fileprefix = this.sessionData.osOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
            os_ = imagingType(ipr.typ, os_);
        end
        function this = solve(this)
            obsWaterMetab = this.obsFromAif(this.artery_water_metab, this.canonical_f);
            this.modelB12 = this.buildQuadraticModel(this.canonical_cbf, obsWaterMetab); % N.B. mapping CBF -> obs
            
            obsOxygen = this.obsFromAif(this.artery_oxygen, this.canonical_f);
            this.modelB34 = this.buildQuadraticModel(this.canonical_cbf, obsOxygen); % N.B. mapping CBF -> obs
            
            obsPet = this.obsFromTac(this.measurement);
            cbfImg = this.cbf.fourdfp.img;
            cbvImg = this.cbv.fourdfp.img;
            poly12 = this.b1*cbfImg.^2 + this.b2*cbfImg;
            poly34 = this.b3*cbfImg.^2 + this.b4*cbfImg;
            numerator = obsPet - poly12 - this.integral_artery_oxygen*cbvImg;
            denominator = poly34 - 0.835*this.integral_artery_oxygen*cbvImg;
            this.img_ = numerator ./ denominator;
            this.img_(isnan(this.img_)) = 0;
            this.img_(this.img_ < 0) = 0;
            this.img_(this.img_ > 1) = 1;
        end
    end
    
    %% PROTECTED
    
    properties (Access = protected)
        artery_oxygen_
        artery_water_metab_
        cbf_
        cbv_
        fs_
        integral_artery_oxygen_
        vs_
    end
    
	methods (Access = protected)		  
 		function this = QuadraticNumericMintun1984(varargin)
 			%% QUADRATICNUMERICMINTUN1984
 			%  @param .

 			this = this@mloxygen.QuadraticNumeric(varargin{:});
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            ip.PartialMatching = false;
            addParameter(ip, 'oo', [], @(x) isnumeric(x))
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            this.measurement = ipr.oo; 			
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

