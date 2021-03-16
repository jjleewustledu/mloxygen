classdef AugmentedNumericRaichle1983 < handle & mloxygen.DispersedNumericRaichle1983
	%% AugmentedNUMERICRAICHLE1983S  

	%  $Revision$
 	%  was created 10-Sep-2020 22:24:50 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.

    
    methods (Static)
        function [this,tac_,aif_] = createFromDeviceKit(devkit, devkit2, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param required devkit2 is mlpet.IDeviceKit.
            %  @param required scanner is an mlpet.AbstractDevice.
            %  @param required scanner2 is an mlpet.AbstractDevice.
            %  @param required arterial is an mlpet.AbstractDevice.
            %  @param required arterial2 is an mlpet.AbstractDevice.
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}, default := 'simulanneal'.
            %  @param roi is mlfourd.ImagingContext2.
            %  @param roi2 is mlfourd.ImagingContext2.
            %  @param histology in {'g' 'w' 's'}
            %  @param Dt_aif isscalar.
            %  @param sigma0, default from mloptimization.SimulatedAnnealing.
            %  @param fileprefix, default from devkit.
            %  @param fracMixing in [0 1] for mixing tacs and aifs.
            %  @return this.
            %  @return tac_, blurred by ipr.blurFdg.
            %  @return aif_.
            
            import mloxygen.AugmentedNumericRaichle1983
            import mlpet.AugmentedData.mixTacsAifs          
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addRequired(ip, 'devkit2', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'scanner', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'scanner2', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'arterial', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'arterial2', [], @(x) isa(x, 'mlpet.AbstractDevice'))            
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'roi2', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'Dt_aif', 0, @isscalar)
            addParameter(ip, 'fracMixing', 0.9, @isscalar)
            parse(ip, devkit, devkit2, varargin{:})
            ipr = ip.Results;
            
            % mix components for augmentation
            [tac_,timesMid_,aif_,Dt] = mixTacsAifs(devkit, devkit2, varargin{:});
            fp = sprintf('mlglucose_AugmentedNumericRaichle1983_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS')); 
            
            this = AugmentedNumericRaichle1983( ...
                'ho', tac_, ...
                'solver', 'simulanneal', ...
                'devkit', devkit, ...
                'Dt', Dt, ...
                'times_sampled', timesMid_, ...
                'artery_interpolated', aif_, ...
                'fileprefix', fp, ...
                varargin{:});      
            this.Dt_aif = ipr.Dt_aif;
        end
    end
    
    %% PROTECTED
    
    methods (Access = protected)        
 		function this = AugmentedNumericRaichle1983(varargin)
 			%% DISPERSEDNUMERICRAICHLE1983
            %  @param ho is numeric.
            %  @param solver is in {'simulanneal'}.
            %  @param devkit is mlpet.IDeviceKit.   
            %  @param Dt is numeric, s of time-shifting for AIF.
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

 			this = this@mloxygen.DispersedNumericRaichle1983(varargin{:});
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

