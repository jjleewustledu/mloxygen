classdef AugmentedNumericMartin1987 < handle & mloxygen.DispersedNumericMartin1987
	%% AUGMENTEDNUMERICMARTIN1987

	%  $Revision$
 	%  was created 29-Apr-2020 23:31:38 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1319299 (R2019b) Update 5 for MACI64.  Copyright 2020 John Joowon Lee.

    
    methods (Static)        
        function [this,tac_,aif_] = createFromDeviceKit(devkit, devkit2, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param required devkit2 is mlpet.IDeviceKit.
            %  @param required scanner is an mlpet.AbstractDevice.
            %  @param required scanner2 is an mlpet.AbstractDevice.
            %  @param required arterial is an mlpet.AbstractDevice.
            %  @param required arterial2 is an mlpet.AbstractDevice.
            %  @param roi is mlfourd.ImagingContext2.
            %  @param roi2 is mlfourd.ImagingContext2.
            %  @param DtMixing isscalar.
            %  @param fracMixing in [0 1] for mixing tacs and aifs.
            %  @param T0 isscalar.
            %  @param Tf isscalar.
            %  @return this.
            %  @return tac_, blurred by ipr.blurFdg.
            %  @return aif_.
            
            import mloxygen.AugmentedNumericMartin1987
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
            addParameter(ip, 'T0', 120, @isscalar)
            addParameter(ip, 'Tf', 240, @isscalar)
            addParameter(ip, 'DtMixing', 0, @isscalar)
            addParameter(ip, 'fracMixing', 0.5, @isscalar)
            parse(ip, devkit, devkit2, varargin{:})
            ipr = ip.Results;
            
            % mix components for augmentation
            [tac_,timesMid_,aif_,Dt] = mixTacsAifs(devkit, devkit2, varargin{:});
            fp = sprintf('mlglucose_AugmentedNumericMartin1987_createFromDeviceKit_dt%s', datestr(now, 'yyyymmddHHMMSS')); 
            
            this = AugmentedNumericMartin1987( ...
                'oc', tac_, ...
                'solver', '', ...
                'devkit', devkit, ...
                'Dt', Dt, ...
                'times_sampled', timesMid_, ...
                'artery_interpolated', aif_, ...
                'fileprefix', fp, ...
                'T0', ipr.T0, ...
                'Tf', min([ipr.Tf max(timesMid_) length(aif_)-1]), ...
                varargin{:});
            this.DtMixing = ipr.DtMixing;
        end
    end
    
    %% PROTECTED

	methods (Access = protected)		  
 		function this = AugmentedNumericMartin1987(varargin)
 			%% AUGMENTEDNUMERICMARTIN1987
            %  @param ho is numeric.
            %  @param devkit is mlpet.IDeviceKit.   
            %  @param Dt is numeric, s of time-shifting for AIF.
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
            
            this = this@mloxygen.DispersedNumericMartin1987(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

