classdef DispersedNumericRaichle1983 < handle & mlpet.AugmentedData & mloxygen.Raichle1983
	%% DISPERSEDNUMERICRAICHLE1983  

	%  $Revision$
 	%  was created 10-Sep-2020 22:24:50 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    
    properties (Constant)
        LENK = 4
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
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}, default := 'simulanneal'.
            %  @param roi is mlfourd.ImagingContext2.
            %  @param histology in {'g' 'w' 's'}
            %  @param sigma0, default from mloptimization.SimulatedAnnealing.
            %  @param fileprefix, default from devkit.
            %  @return this.
            %  @return tac_, blurred by ipr.blurFdg.
            %  @return aif_.
            
            import mloxygen.DispersedNumericRaichle1983
            import mlpet.AugmentedData.mixTacAif       
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'scanner', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'arterial', [], @(x) isa(x, 'mlpet.AbstractDevice'))           
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            parse(ip, devkit, varargin{:})
            
            % mix components for augmentation       
            [tac_,timesMid_,aif_,Dt] = mixTacAif(devkit, varargin{:});
            
            %
            
            fp = sprintf('mlglucose_DispersedNumericRaichle1983_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS')); 
            
            this = DispersedNumericRaichle1983( ...
                'ho', tac_, ...
                'solver', 'simulanneal', ...
                'devkit', devkit, ...
                'Dt', Dt, ...
                'times_sampled', timesMid_, ...
                'artery_interpolated', aif_, ...
                'fileprefix', fp, ...
                varargin{:});      
        end
    end

	methods 
        
        %% GET
        
        function g = get.artery_sampled(this)
            g = this.model.solutionOnScannerFrames( ...
                this.artery_interpolated(this.tBuffer+1:end), this.times_sampled);
        end
        function g = get.Delta(this)
            g = this.k4();
        end
        function g = get.Dt(this)
            g = this.strategy_.Dt;
        end
        
        %%
        
        function ks = buildKs(this, varargin)
            this = solve(this, varargin{:});
            ks = [k1(this) k2(this) k3(this) k4(this) k5(this)];
        end
        function ho = checkSimulated(this, varargin)
            %% CHECKSIMULATED simulates tissue activity with passed and internal parameters without changing state.
            %  @param required ks is [k1 k2 k3 k4 Dt].
            %  @param aif is numeric; default is this.artery_interpolated for model state.
            %  @return ho simulation is numeric.
        
            ip = inputParser;
            addOptional(ip, 'ks', this.ks(), @isnumeric)
            addParameter(ip, 'aif', this.artery_interpolated, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;  
            
            ho = this.model.simulated(ipr.ks, 'aif', ipr.aif, 'Dt', this.Dt);
        end
        function [k,sk] = k4(this, varargin)
            [k,sk] = k4(this.strategy_, varargin{:});
        end
        function [k,sk] = k5(this, varargin)
            k = this.Dt;
            sk = nan;
        end
        function ks_ = ks(this, varargin)
            %% ks == [f PS lambda Delta Dt]
            %  @param 'typ' is char, understood by imagingType.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'typ', 'mlfourd.ImagingContext2', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            k(1) = k1(this.strategy_, varargin{:});
            k(2) = k2(this.strategy_, varargin{:});
            k(3) = k3(this.strategy_, varargin{:});
            k(4) = k4(this.strategy_, varargin{:});
            k(5) = k5(this);
             
            roibin = logical(this.roi);
            ks_ = copy(this.roi.fourdfp);
            ks_.img = zeros([size(this.roi) length(k)]);
            for t = 1:length(k)
                img = zeros(size(this.roi), 'single');
                img(roibin) = k(t);
                ks_.img(:,:,:,t) = img;
            end
            ks_.fileprefix = this.sessionData.fsOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
            ks_ = imagingType(ipr.typ, ks_);
        end
    end
    
    %% PROTECTED
    
    methods (Access = protected)        
 		function this = DispersedNumericRaichle1983(varargin)
 			%% DISPERSEDNUMERICRAICHLE1983
            %  @param ho is numeric.
            %  @param solver is in {'simulanneal'}.
            %  @param devkit is mlpet.IDeviceKit.   
            %  
            %  for mloxygen.DispersedRaichle1983Model: 
            %  @param histology in {'g' 'w' 's'}
            %  @param map is a containers.Map.  Default := DispersedHuang1980Model.preferredMap.
            %  @param times_sampled for scanner is typically not uniform.
            %  @param artery_interpolated must be uniformly interpolated.
            %
            %  for mloxygen.DispersedRaichle1983SimulAnneal:
            %  @param context is mloxygen.Raichle1983.
            %  @param sigma0.
            %  @param fileprefix.

 			this = this@mloxygen.Raichle1983( ...
                'model', mloxygen.DispersedRaichle1983Model(varargin{:}), ...
                varargin{:});            
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'ho', [], @(x) isnumeric(x))
            addParameter(ip, 'solver', 'simulanneal', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            this.measurement = ipr.ho;
            switch lower(ipr.solver)
                case 'simulanneal'
                    this.strategy_ = mloxygen.DispersedRaichle1983SimulAnneal( ...
                        'context', this, varargin{:});
                otherwise
                    error('mloxygen:NotImplementedError', 'DispersedNumericRaichle1983.ipr.solver->%s', ipr.solver)
            end
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

