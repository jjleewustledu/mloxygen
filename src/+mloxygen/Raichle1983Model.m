classdef Raichle1983Model 
	%% RAICHLE1983MODEL provides model data and methods to the strategy design pattern comprising
    %  mloxygen.{Raichle1983, Raichle1983Nest, Raichle1983SimulAnneal, Raichle1983HMC, Raichle1983 LM, Raichle1983BFGS}.
    %  It operates on single voxels or regions.

	%  $Revision$
 	%  was created 10-Sep-2020 19:42:35 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.
 	
	properties
        artery_interpolated
 		map
        times_sampled 		
    end
    
    methods (Static)
        function qs       = solution(ks, artery_interpolated)
            %  @param artery_interpolated is uniformly with at high sampling freq. starting at time = 0.
            
            ALPHA = 0.005670305; % log(2)/halflife in 1/s
            
            f = ks(1);
            PS = ks(2);
            lambda = ks(3); 
            E = max(1 - exp(-PS/f), eps);
            [~,idx] = max(artery_interpolated > 0.1*max(artery_interpolated));
            n = min(length(artery_interpolated), idx+59); % limit duration of scan sampling
            times = 0:1:n-1;
             
            qs = E*f*conv(exp(-E*f*times/lambda - ALPHA*times), artery_interpolated);
            qs = qs(1:n);
        end
        function qs       = sampled(ks, artery_interpolated, times_sampled)
            %  @param artery_interpolated is uniformly sampled at high sampling freq.
            %  @param times_sampled are samples scheduled by the time-resolved PET reconstruction
            
            import mloxygen.Raichle1983Model.solution  
            qs = solution(ks, artery_interpolated);
            times_sampled = floor(times_sampled - times_sampled(1)) + 1; % times_sampled(1) == 1
            times_sampled = times_sampled(times_sampled <= length(qs));
            qs = qs(times_sampled);
        end
         
        function logp = log_likelihood(Z, Sigma)
            %% for Raichle1983HMC
            
            logp = sum(-log(Sigma) - 0.5*log(2*pi) - 0.5*Z.^2); % scalar
        end
        function loss = simulanneal_objective(ks, artery_interpolated, times_sampled, qs0, sigma0)
            import mloxygen.Raichle1983Model.sampled          
            qs = sampled(ks, artery_interpolated, times_sampled);            
            loss = 0.5 * sum((1 - qs ./ qs0).^2) / sigma0^2; % + sum(log(sigma0*qs0)); % sigma ~ sigma0 * qs0
        end   
        function m = preferredMap()
            %% init from Raichle J Nucl Med 24:  790-798, 1983.
            %  PS in [0.014 0.0245 0.0588] Hz for white, brain, grey;
            %  PS min := PS(1) - (PS(2) - PS(1))
            %  PS max := PS(3) + (PS(3) - PS(2))
            %  PS init := PS(2)
            %  PS sigma := 0.08*(PS init)
            
            m = containers.Map;
            m('k1')  = struct('min', 0.00175, 'max', 0.028,  'init', 0.00777, 'sigma', 6e-4); % f / s
            m('k2')  = struct('min', 0.0035,  'max', 0.0931, 'init', 0.0245,  'sigma', 0.00196); % PS / s
            m('k3')  = struct('min', 0.8,     'max', 1,      'init', 0.95,    'sigma', 0.05); % lambda in [0 1]
        end
    end

	methods		  
 		function this = Raichle1983Model(varargin)
 			
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'map', this.preferredMap(), @(x) isa(x, 'containers.Map'))
            addParameter(ip, 'times_sampled', [], @isnumeric)
            addParameter(ip, 'artery_interpolated', [], @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.map = ipr.map;
            this.times_sampled = ipr.times_sampled;
            if this.times_sampled(end)+1 ~= length(ipr.artery_interpolated)
                this.artery_interpolated = ...
                    pchip(0:length(ipr.artery_interpolated)-1, ipr.artery_interpolated, 0:this.times_sampled(end));
            else
                this.artery_interpolated = ipr.artery_interpolated;
            end
        end
        
        function ho  = simulated(this, varargin)
            %% SIMULATED simulates tissue activity with passed and internal parameters.
            %  @param required ks is [k1 k2 k3 Dt].
            %  @param aif is numeric; default is this.artery_interpolated for model state.
        
            ip = inputParser;
            addRequired(ip, 'ks', @isnumeric)
            addParameter(ip, 'aif', this.artery_interpolated, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;            
            
            ks = ipr.ks(1:3);
            Dt = ipr.ks(4);
            if Dt ~= 0
                times = 0:length(ipr.aif)-1;
                aif = pchip(times+Dt, ipr.aif, times);
            else
                aif = ipr.aif;
            end
            ho = mloxygen.Raichle1980Model.sampled(ks, aif, this.times_sampled);
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

