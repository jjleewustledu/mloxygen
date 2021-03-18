classdef AugmentedImagingMartin1987 < handle & mlpet.AugmentedData & mlpet.TracerKinetics
	%% AUGMENTEDIMAGINGMARTIN1987

	%  $Revision$
 	%  was created 29-Apr-2020 23:31:38 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1319299 (R2019b) Update 5 for MACI64.  Copyright 2020 John Joowon Lee.

    
    properties (Constant)        
        T = mlpet.TracerKineticsModel.T
    end
    
    methods (Static)        
        function this = createFromDeviceKit(devkit, devkit2, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param required devkit2 is mlpet.IDeviceKit.
            %  @param required scanner is an mlpet.AbstractDevice.
            %  @param required scanner2 is an mlpet.AbstractDevice.
            %  @param required arterial is an mlpet.AbstractDevice.
            %  @param required arterial2 is an mlpet.AbstractDevice.
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}, default := 'simulanneal'.
            %  @param roibin is mlfourd.ImagingContext2.
            %  @param roibin2 is mlfourd.ImagingContext2.
            %  @param DtMixing isscalar.
            %  @param fracMixing in [0 1] for mixing tacs and aifs.
            %  @param T0 isscalar.
            %  @param Tf isscalar.
            %  @return this.
            
            import mloxygen.AugmentedImagingMartin1987
            import mlpet.AugmentedData.mixScannersAifs
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addRequired(ip, 'devkit2', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'scanner', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'scanner2', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'arterial', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'arterial2', [], @(x) isa(x, 'mlpet.AbstractDevice'))            
            addParameter(ip, 'roibin', [], @islogical)
            addParameter(ip, 'roibin2', [], @islogical)
            addParameter(ip, 'T', AugmentedImagingMartin1987.T, @isscalar0)
            addParameter(ip, 'DtMixing', 0, @isscalar)
            addParameter(ip, 'fracMixing', 0.5, @isscalar)
            addParameter(ip, 'T0', 120, @isscalar)
            addParameter(ip, 'Tf', 240, @isscalar)
            parse(ip, devkit, devkit2, varargin{:})
            ipr = ip.Results;
            
            [scan_,timesMid_,aif_] = mixScannersAifs(varargin{:});
            sz = size(ipr.scanner.imagingContext);
            scanmat_ = reshape(scan_, [sz(1)*sz(2)*sz(3) sz(4)]);
            
            this = AugmentedImagingMartin1987( ...
                'devkit', ipr.devkit, ...
                'ocmat', scanmat_, ...
                'times_sampled', timesMid_, ...
                'artery_interpolated', aif_, ...
                'T0', 120, ...
                'Tf', min([ipr.Tf max(timesMid_) length(aif_)-1]), ...
                varargin{:});
            this.DtMixing = ipr.DtMixing;
        end
    end
    
    methods 
        function this = solve(this, roivecbin)
            scale = 1/this.RATIO_SMALL_LARGE_HCT;
            tac = this.ocmat_(roivecbin, :);            
            this.product_ = ...
                scale * ...
                trapz(this.tacTimeRange_, tac) / ...
                trapz(this.aifTimeRange_, this.artery_interpolated_);
        end
        function val = v1(this, varargin)
            val = this.product_;
        end
        
    end
    
    %% PROTECTED

	methods (Access = protected)		  
 		function this = AugmentedImagingMartin1987(varargin)
 			%% AUGMENTEDIMAGINGMARTIN1987
            %  @param devkit is mlpet.IDeviceKit.
            %  @param oc is numeric ~ [Nvoxels Ntimes].
 			%  @param T0 is numeric.
            %  @param Tf is numeric.
            
            this = this@mlpet.TracerKinetics(varargin{:});

            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'ocmat', [], @isnumeric)
            addParameter(ip, 'times_sampled', [], @isnumeric)
            addParameter(ip, 'artery_interpolated', [], @isnumeric)
            addParameter(ip, 'T0', nan, @isfinite);
            addParameter(ip, 'Tf', nan, @isfinite);
            parse(ip, varargin{:});
            ipr = ip.Results;
            this.ocmat_ = ipr.ocmat;
            this.times_sampled_ = ipr.times_sampled;
            this.artery_interpolated_ = ipr.artery_interpolated;
            this.T0_ = ipr.T0;
            this.Tf_ = ipr.Tf;
            
            ts = this.times_sampled_;
            this.tacTimeRange_ = ts(this.T0_ <= ts & ts <= this.Tf_);
            this.ocmat_ = this.ocmat_(:, this.T0_ <= ts & ts <= this.Tf_);
            
            this.aifTimeRange_ = this.T0_:this.Tf_;
            this.artery_interpolated_ = this.artery_interpolated_(this.T0_+1:this.Tf_+1);
 		end
 	end 
    
    %% PRIVATE
    
    properties (Access = private)
        aifTimeRange_
        artery_interpolated_
        ocmat_
        product_  % := v1 as scalar
        tacTimeRange_
        times_sampled_
        T0_
        Tf_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

