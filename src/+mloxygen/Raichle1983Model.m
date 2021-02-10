classdef (Abstract) Raichle1983Model < mlpet.TracerKineticsModel
	%% RAICHLE1983MODEL provides model data and methods to the strategy design pattern comprising
    %  mloxygen.{Raichle1983, Raichle1983Nest, Raichle1983SimulAnneal, Raichle1983HMC, Raichle1983 LM, Raichle1983BFGS}.
    %  It operates on single voxels or regions.

	%  $Revision$
 	%  was created 10-Sep-2020 19:42:35 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.
    
    methods (Abstract, Static)
        solution
    end
    
    methods (Static)
        function this = createFromMetric(m, varargin)
            assert(ischar(m))
            switch m
                case {'' 'fs'}
                    this = mloxygen.DispersedRaichle1983Model(varargin{:});
                case 'gs'
                    this = mloxygen.ConstantERaichle1983Model(varargin{:});
                otherwise
                    error('mloxygen:ValueError', 'Raichle1983Model.createFromMetric.m -> %s', m)
            end
        end
        function logp = log_likelihood(Z, Sigma)
            %% for Raichle1983HMC
            
            logp = sum(-log(Sigma) - 0.5*log(2*pi) - 0.5*Z.^2); % scalar
        end
        function loss = loss_function(ks, artery_interpolated, times_sampled, measurement, sigma0)
            import mloxygen.Raichle1983Model.sampled  
            estimation  = sampled(ks, artery_interpolated, times_sampled);
            measurement = measurement(1:length(estimation));
            taus        = diff(times_sampled);
            taus        = [taus taus(end)];
            taus        = taus(1:length(estimation));
            positive    = measurement > 0;
            e           = estimation .* taus;
            m           = measurement .* taus;
            eoverm      = e(positive)./m(positive);
            Q           = sum((1 - eoverm).^2);
            loss        = 0.5*Q/sigma0^2; % + sum(log(sigma0*measurement)); % sigma ~ sigma0*measurement
        end  
        function m    = preferredMap()
            %% init from Raichle J Nucl Med 24:790-798, 1983; Herscovitch JCBFM 5:65-69 1985; Herscovitch JCBFM 7:527s-542 1987
            %  PS in [0.0140 0.0245 0.0588] Hz for white, brain, grey;
            %  PS min := PS(1) - (PS(2) - PS(1))
            %  PS max := PS(3) + (PS(3) - PS(2))
            %  PS init := PS(2)
            %  PS sigma := 0.08*(PS init)
            %  lambda described in Table 2
            
            m = containers.Map;
            m('k1') = struct('min', 0.0043, 'max', 0.0155, 'init', 0.00777, 'sigma', 3.89e-4); % f / s, max ~ 0.0155
            m('k2') = struct('min', 0.0137, 'max', 0.0266, 'init', 0.0228,  'sigma', 0.002); % PS / s
            m('k3') = struct('min', 0.608,  'max', 1.06,   'init', 0.945,   'sigma', 0.05); % lambda in mL/mL
            m('k4') = struct('min', 0.01,   'max', 2,      'init', 1,       'sigma', 0.1); % Delta for cerebral dispersion
        end
        function qs   = sampled(ks, artery_interpolated, times_sampled)
            %  @param artery_interpolated is uniformly sampled at high sampling freq.
            %  @param times_sampled are samples scheduled by the time-resolved PET reconstruction
            
            import mloxygen.Raichle1983Model.solution 
            import mlpet.TracerKineticsModel.solutionOnScannerFrames  
            qs = solution(ks, artery_interpolated);
            qs = solutionOnScannerFrames(qs, times_sampled);
        end
        function loss = simulanneal_objective(ks, artery_interpolated, times_sampled, qs0, sigma0)
            import mloxygen.Raichle1983Model.sampled          
            qs = sampled(ks, artery_interpolated, times_sampled);            
            loss = 0.5 * sum((1 - qs ./ qs0).^2) / sigma0^2; % + sum(log(sigma0*qs0)); % sigma ~ sigma0 * qs0
        end
    end

	methods		  
 		function this = Raichle1983Model(varargin)
            this = this@mlpet.TracerKineticsModel(varargin{:});            
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'histology', '', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this = this.adjustMapForHistology(ipr.histology);
        end
        function this = adjustMapForHistology(this, histology)
            %% use PS ranges from Herscovitch et al 1987 Table 2
            
            switch histology
                case 'g'
                    this.map('k2') = struct('min', 0.017,  'max', 0.0266, 'init', 0.0218,  'sigma', 0.002); % PS / s
                    this.map('k3') = struct('min', 0.738,  'max', 1.06,   'init', 1.02,    'sigma', 0.05); % lambda in mL/mL  
                    %this.map('k3') = struct('min', 0.987, 'max', 1.06,  'init', 1.02,    'sigma', 0.05); % lambda in mL/mL                    
                case 'w'
                    this.map('k2') = struct('min', 0.0137, 'max', 0.0142, 'init', 0.014,   'sigma', 0.002); % PS / s
                    this.map('k3') = struct('min', 0.608,  'max', 0.882,  'init', 0.851,   'sigma', 0.05); % lambda in mL/mL
                    %this.map('k3') = struct('min', 0.819,  'max', 0.882, 'init', 0.851,   'sigma', 0.05); % lambda in mL/mL
                case 's' % subcortical
                    this.map('k2') = struct('min', 0.0159, 'max', 0.0215, 'init', 0.0187,  'sigma', 0.002); % PS / s
                    this.map('k3') = struct('min', 0.738,  'max', 0.97,   'init', 0.924,   'sigma', 0.05); % lambda in mL/mL
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
            if length(ipr.ks) > 4
                ipr.Dt = ipr.ks(5);
            end
            if ipr.Dt ~= 0
                times = 0:length(ipr.aif)-1;
                aif = pchip(times - ipr.Dt, ipr.aif, times); % remove the delay Dt found by model
            else
                aif = ipr.aif;
            end
            ho = mloxygen.Raichle1983Model.sampled(ks, aif, this.times_sampled);
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

