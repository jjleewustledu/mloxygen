classdef DispersedRaichle1983Model 
	%% DISPERSEDRAICHLE1983MODEL provides model data and methods to the strategy design pattern comprising
    %  mloxygen.{Raichle1983, DispersedRaichle1983SimulAnneal}.
    %  It operates on single voxels or regions.  It includes a dispersion parameter $\Delta$:
    %  aif_\text{disp}(t) = aif(t) \otimes e^{-\Delta t}. 

	%  $Revision$
 	%  was created 10-Sep-2020 22:23:42 by jjlee,
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
            Delta = ks(4);
            E = max(1 - exp(-PS/f), 0.7);
            E = min(E, 0.93);
            %[~,idx] = max(artery_interpolated > 0.1*max(artery_interpolated));            
            %n = min(length(artery_interpolated), idx+119); % limit duration of scan sampling
            n = length(artery_interpolated);
            times = 0:1:n-1;
             
            % use Delta
            auc0 = trapz(artery_interpolated);
            artery_interpolated1 = conv(artery_interpolated, exp(-Delta*times));
            artery_interpolated1 = artery_interpolated1(1:n);
            artery_interpolated1 = artery_interpolated1*auc0/trapz(artery_interpolated1);
            
            % use E, f, lambda
            qs = E*f*conv(exp(-E*f*times/lambda - ALPHA*times), artery_interpolated1);
            qs = qs(1:n);
        end
        function qs       = sampled(ks, artery_interpolated, times_sampled)
            %  @param artery_interpolated is uniformly sampled at high sampling freq.
            %  @param times_sampled are samples scheduled by the time-resolved PET reconstruction
            
            import mloxygen.DispersedRaichle1983Model.solution
            qs = solution(ks, artery_interpolated);
            idx_sampled = floor(times_sampled - times_sampled(1)) + 1; % times_sampled(1) == 1
            idx_sampled = idx_sampled(idx_sampled <= length(qs));
            qs = qs(idx_sampled);
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
            %% init from Raichle J Nucl Med 24:790-798, 1983; Herscovitch JCBFM 5:65-69 1985; Herscovitch JCBFM 7:527s-542 1987
            %  PS in [0.0140 0.0245 0.0588] Hz for white, brain, grey;
            %  PS min := PS(1) - (PS(2) - PS(1))
            %  PS max := PS(3) + (PS(3) - PS(2))
            %  PS init := PS(2)
            %  PS sigma := 0.08*(PS init)
            %  lambda described in Table 2
            
            m = containers.Map;
            m('k1') = struct('min', 0.000833, 'max', 0.028,  'init', 0.00777, 'sigma', 3.89e-4); % f / s
            m('k2') = struct('min', 0.00928,  'max', 0.0368, 'init', 0.0194,  'sigma', 0.002); % PS / s
            %m('k2') = struct('min', 0.7,      'max', 0.93,   'init', 0.825,   'sigma', 0.05); % E_w
            m('k3') = struct('min', 0.797,    'max', 1.09,   'init', 0.945,   'sigma', 0.05); % lambda in mL/mL
            m('k4') = struct('min', 0.333,    'max', 2,      'init', 1,       'sigma', 0.1); % Delta for cerebral dispersion
        end
    end

	methods		  
 		function this = DispersedRaichle1983Model(varargin)
 			
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'map', this.preferredMap(), @(x) isa(x, 'containers.Map'))
            addParameter(ip, 'times_sampled', [], @isnumeric)
            addParameter(ip, 'artery_interpolated', [], @isnumeric)
            addParameter(ip, 'histology', '', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.map = ipr.map;
            this = this.adjustMapForHistology(ipr.histology);
            this.times_sampled = ipr.times_sampled;
            if this.times_sampled(end)-this.times_sampled(1)+1 ~= length(ipr.artery_interpolated)
                this.artery_interpolated = ...
                    pchip(0:length(ipr.artery_interpolated)-1, ipr.artery_interpolated, this.times_sampled(1):this.times_sampled(end));
            else
                this.artery_interpolated = ipr.artery_interpolated;
            end
        end
        
        function this = adjustMapForHistology(this, histology)
            %% use PS ranges from Herscovitch et al 1987 Table 2
            %  use 2*sigma(lambda_gray) ~ 2*std(Herscovitch & Raichle 1987 Table 2)
            
            switch histology
                case 'g'
                    this.map('k2') = struct('min', 0.00788,  'max', 0.0368,  'init', 0.0209,  'sigma', 0.002); % PS / s
                    this.map('k3') = struct('min', 0.955,    'max', 1.09,    'init', 1.02,    'sigma', 0.05); % lambda in mL/mL                    
                case 'w'
                    this.map('k2') = struct('min', 0.0068,  'max', 0.0215,  'init', 0.0139,  'sigma', 0.002); % PS / s
                    this.map('k3') = struct('min', 0.797,    'max', 0.905,   'init', 0.851,   'sigma', 0.05); % lambda in mL/mL
                case 's' % subcortical
                    this.map('k2') = struct('min', 0.0109,   'max', 0.0282,  'init', 0.0193,  'sigma', 0.002); % PS / s
                otherwise
                    % noninformative
            end
        end
        function ho   = simulated(this, varargin)
            %% SIMULATED simulates tissue activity with passed and internal parameters.
            %  @param required ks is [k1 k2 k3 k4 Dt].
            %  @param aif is numeric; default is this.artery_interpolated for model state.
            %  @param Dt is numeric, in sec.
        
            ip = inputParser;
            addRequired(ip, 'ks', @isnumeric)
            addParameter(ip, 'aif', this.artery_interpolated, @isnumeric)
            addParameter(ip, 'Dt', 0, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;            
            
            ks = ipr.ks(1:4);
            if length(ks) > 4
                ipr.Dt = ks(5);
            end
            if ipr.Dt ~= 0
                times = 0:length(ipr.aif)-1;
                aif = pchip(times+ipr.Dt, ipr.aif, times);
            else
                aif = ipr.aif;
            end
            ho = mloxygen.DispersedRaichle1983Model.sampled(ks, aif, this.times_sampled);
        end
 	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

