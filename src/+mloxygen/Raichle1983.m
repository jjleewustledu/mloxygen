classdef Raichle1983 < handle & matlab.mixin.Copyable 
	%% RAICHLE1983 is the context to a strategy design patterns which implements:
    %  mloxygen.{Raichle1983Nest, Raichle1983SimulAnneal, Raichle1983HMC, Raichle1983LM, Raichle1983BFGS}.
    %  For performance considerations, see also https://blogs.mathworks.com/loren/2012/03/26/considering-performance-in-object-oriented-matlab-code/

	%  $Revision$
 	%  was created 30-May-2018 01:54:09 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlpet/src/+mlpet.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
        devkit
        Dt          % time-shift for AIF; Dt < 0 shifts backwards in time.
        measurement % expose for performance when used by mlglucose.Huang1980Strategy
        model       %
        regionTag 		
    end    
    
    methods (Static)
        function this = createFromDeviceKit(devkit, varargin)
            %% 
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param Dt is numeric, s of time-shifting for AIF.
            %  @param ho is numeric.
            %  @param solver is in {'nest' 'simulanneal' 'hmc' 'lm' 'bfgs'}.
            %  @param blurHo := {[], 0, 4.3, ...}            
            %  @param map, default := mloxygen.Raichle1983Model.preferredMap().
            %  @param times_sampled non-uniformly scheduled by the time-resolved PET reconstruction.
            %  @param artery_interpolated, default from devkit.
            %  @param sigma0, default from mloptimization.SimulatedAnnealing.
            %  @param fileprefix, default from devkit. 
            %  @return this.
            
            this = mloxygen.NumericRaichle1983.createFromDeviceKit(devkit, varargin{:});
        end
    end
    
    methods
        function ks = buildKs(this, varargin)
            this = solve(this, varargin{:});
            ks = [k1(this) k2(this) k3(this)];
        end
        function [k,sk] = k1(this, varargin)
            [k,sk] = k1(this.strategy_, varargin{:});
        end
        function [k,sk] = k2(this, varargin)
            [k,sk] = k2(this.strategy_, varargin{:});
        end
        function [k,sk] = k3(this, varargin)
            [k,sk] = k3(this.strategy_, varargin{:});
        end
        function [k,sk] = ks(this, varargin)
            k = zeros(1,3);
            sk = zeros(1,3);
            [k(1),sk(1)] = k1(this.strategy_, varargin{:});
            [k(2),sk(2)] = k2(this.strategy_, varargin{:});
            [k(3),sk(3)] = k3(this.strategy_, varargin{:});
        end
        function h = plot(this, varargin)
            h = this.strategy_.plot(varargin{:});
        end
        function this = simulated(this, varargin)
            this.measurement = this.model.simulated(varargin{:});
            this.strategy_.Measurement = this.measurement; % strategy_ needs value copies for performance
        end
        function this = solve(this, varargin)
            this.strategy_ = solve(this.strategy_, varargin{:});
        end
    end
    
    %% PROTECTED
    
    properties (Access = protected)
        strategy_
    end

	methods (Access = protected)		  
 		function this = Raichle1983(devkit, varargin)
            %% RAICHLE1983
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param Dt is numeric, s of time-shifting for AIF.
            
            ip = inputParser;            
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'Dt', 0, @isscalar)
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            
            this.devkit = ipr.devkit;
            this.Dt = ipr.Dt;
            this.regionTag = this.devkit.sessionData.regionTag;
        end
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

