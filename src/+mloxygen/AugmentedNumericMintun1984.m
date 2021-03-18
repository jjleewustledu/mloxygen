classdef AugmentedNumericMintun1984 < handle & mloxygen.DispersedNumericMintun1984
	%% AUGMENTEDNUMERICMINTUN1984  

	%  $Revision$
 	%  was created 09-Feb-2021 23:27:40 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.9.0.1570001 (R2020b) Update 4 for MACI64.  Copyright 2021 John Joowon Lee.
 	   
    
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
            %  @param DtMixing isscalar.
            %  @param sigma0, default from mloptimization.SimulatedAnnealing.
            %  @param fileprefix, default from devkit.
            %  @param fracMixing in [0 1] for mixing tacs and aifs.
            %  @return this.
            %  @return tac_, blurred by ipr.blurFdg.
            %  @return aif_.
            
            import mloxygen.AugmentedNumericMintun1984
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
            addParameter(ip, 'DtMixing', 0, @isscalar)
            addParameter(ip, 'fracMixing', 0.9, @isscalar)
            parse(ip, devkit, devkit2, varargin{:})
            ipr = ip.Results;
            
            % mix components for augmentation
            [tac_,timesMid_,aif_,Dt] = mixTacsAifs(devkit, devkit2, varargin{:});
            fp = sprintf('mlglucose_AugmentedNumericMintun1984_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS'));             
            fs_R_M_ = mloxygen.DispersedMintun1984Model.fs_R_M(devkit.sessionData, ipr.roi);
            
            this = AugmentedNumericMintun1984( ...
                'oo', tac_, ...
                'solver', 'simulanneal', ...
                'devkit', devkit, ...
                'Dt', Dt, ...
                'times_sampled', timesMid_, ...
                'artery_interpolated', aif_, ...
                'fileprefix', fp, ...
                'fs_Raichle_Martin', fs_R_M_, ...
                varargin{:});      
            this.DtMixing = ipr.DtMixing;
        end
    end

    %% PROTECTED
    
	methods (Access = protected)		  
 		function this = AugmentedNumericMintun1984(varargin)
 			%% AUGMENTEDNUMERICMINTUN1984
            %  @param oo is numeric.
            %  @param solver is in {'simulanneal'}.
            %  @param devkit is mlpet.IDeviceKit.   
            %  @param Dt is numeric, s of time-shifting for AIF.
            %  
            %  for mloxygen.DispersedMintun1984Model: 
            %  @param fs_Raichle_Martin is numeric, containing model parameters from Raichle1983 & Martin1987.
            %  @param map is a containers.Map.  Default := DispersedHuang1980Model.preferredMap.
            %  @param times_sampled for scanner is typically not uniform.
            %  @param artery_interpolated must be uniformly interpolated.
            %
            %  for mloxygen.DispersedMintun1984SimulAnneal:
            %  @param context is mloxygen.Mintun1984.
            %  @param sigma0.
            %  @param fileprefix.

 			this = this@mloxygen.DispersedNumericMintun1984(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

