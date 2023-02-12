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
                'tObs', 40, ...
                'artery_interpolated', aif_, ...
                'fileprefix', fp, ...
                varargin{:}); 
            this.cbf_ = this.cbfOnAtlas();
            this.cbv_ = this.cbvOnAtlas();
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
        
        function obj = cbfOnAtlas(this, varargin)
            sd = this.sessionData;
            if isa (this.sessionData, 'mlnipet.SessionData')
                obj = this.sessionData.cbfOnAtlas('typ', 'ImagingContext2', 'dateonly', true, 'tags', sd.regionTag);
                return
            end
            if isa (this.sessionData, 'mlpipeline.ImagingMediator')
                g = glob(fullfile(this.sessionData.sessionsPath, 'ses-*', 'pet', '*_cbf_*_voxels.nii.gz'));
                obj = mlfourd.ImagingContext2(g{1});
                return
            end
            error('mloxygen:RuntimeError', stackstr())
        end        
        function obj = cbvOnAtlas(this, varargin)
            sd = this.sessionData;
            if isa (this.sessionData, 'mlnipet.SessionData')
                obj = this.sessionData.cbvOnAtlas('typ', 'ImagingContext2', 'dateonly', true, 'tags', sd.regionTag);
                return
            end
            if isa (this.sessionData, 'mlpipeline.ImagingMediator')
                g = glob(fullfile(this.sessionData.sessionsPath, 'ses-*', 'pet', '*_cbv_*_voxels.nii.gz'));
                obj = mlfourd.ImagingContext2(g{1});
                return
            end
            error('mloxygen:RuntimeError', stackstr())
        end
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
            
            os_ = copy(this.roi.imagingFormat);
            os_.img = single(this.img_);
            os_.fileprefix = this.osOnAtlas.fileprefix;
            os_ = imagingType(ipr.typ, os_);
        end        
        function obj = osOnAtlas(this, varargin)
            if isa (this.sessionData, 'mlnipet.SessionData')
                obj = this.sessionData.osOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
                return
            end
            if isa (this.sessionData, 'mlpipeline.ImagingMediator')
                tags = strip([this.blurTag this.regionTag], '_');
                obj = this.sessionData.metricOnAtlas('os', tags);
                return
            end
            error('mloxygen:RuntimeError', stackstr())
        end
        function this = solve(this)
            this = this.solve_voxels();
        end
        function this = solve_voxels(this)
            obsWaterMetab = this.obsFromAif(this.artery_water_metab, this.canonical_f); % time series
            this.modelB12 = this.buildQuadraticModel(this.canonical_cbf, obsWaterMetab); % N.B. nonlin model mapping CBF -> obs
            
            obsOxygen = this.obsFromAif(this.artery_oxygen, this.canonical_f); % time series
            this.modelB34 = this.buildQuadraticModel(this.canonical_cbf, obsOxygen); % N.B. nonlin model mapping CBF -> obs

            ps = sqrt(0.5)*this.sessionData.petPointSpread(); % for applying twice

            obsPet_ifc = copy(this.roi.imagingFormat);
            obsPet_ifc.img = this.obsFromTac(this.measurement);
            obsPet_ic = mlfourd.ImagingContext2(obsPet_ifc);
            obsPet_ic = obsPet_ic.blurred(ps);

            cbf_ic = this.cbf.blurred(ps); % pre-blur polynomial
            cbv_ic = this.cbv.blurred(ps); % pre-blur polynomial
            poly12_ic = cbf_ic.^2.*this.b1 + cbf_ic.*this.b2;
            poly34_ic = cbf_ic.^2.*this.b3 + cbf_ic.*this.b4;
            numerator_ic = obsPet_ic - poly12_ic - cbv_ic*this.integral_artery_oxygen;
            numerator_ic = numerator_ic.blurred(ps); % post-blur
            denominator_ic = poly34_ic - cbv_ic*0.835*this.integral_artery_oxygen;
            denominator_ic = denominator_ic.blurred(ps); % post-blur
            ratio_ic = numerator_ic ./ denominator_ic;            
            ratio_ic = ratio_ic.thresh(0);
            ratio_ic = ratio_ic.uthresh(1);
            ratio_ic = ratio_ic.scrubNanInf();

            this.img_ = ratio_ic.imagingFormat.img;
        end
        function this = solve_on_wmparc1(this)
            obsWaterMetab = this.obsFromAif(this.artery_water_metab, this.canonical_f);
            this.modelB12 = this.buildQuadraticModel(this.canonical_cbf, obsWaterMetab); % N.B. mapping CBF -> obs
            
            obsOxygen = this.obsFromAif(this.artery_oxygen, this.canonical_f);
            this.modelB34 = this.buildQuadraticModel(this.canonical_cbf, obsOxygen); % N.B. mapping CBF -> obs
            
            obsPet = mlfourd.ImagingContext2(this.obsFromTac(this.measurement));
            obsPet.fileprefix = clientname(true, 2);
            poly12 = this.cbf.^2*this.b1 + this.cbf*this.b2;
            poly34 = this.cbf.^2*this.b3 + this.cbf*this.b4;


            
            wmparc1 = this.sessionData.wmparc1OnAtlas('typ', 'mlfourd.ImagingContext2');
            this.img_ = zeros(size(wmparc1), 'single');

            numerator = obsPet - poly12 - this.cbv*this.integral_artery_oxygen;
            denominator = poly34 - this.cbv*this.integral_artery_oxygen*0.835;
            fraction = numerator ./ denominator;
            this.img_ = this.img_ + fraction.selectFourdfp.img;
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

