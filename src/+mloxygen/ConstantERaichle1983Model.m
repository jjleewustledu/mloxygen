classdef ConstantERaichle1983Model < mloxygen.Raichle1983Model
	%% CONSTANTERAICHLE1983MODEL provides model data and methods to the strategy design pattern comprising
    %  mloxygen.{Raichle1983, DispersedRaichle1983SimulAnneal}.
    %  It operates on single voxels or regions.  It includes a dispersion parameter $\Delta$:
    %  aif_\text{disp}(t) = aif(t) \otimes e^{-\Delta t}. 

	%  $Revision$
 	%  was created 22-Sep-2020 23:58:16 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    methods (Static)
        function qs       = solution(ks, artery_interpolated)
            %  @param artery_interpolated is uniformly with at high sampling freq. starting at time = 0.

            import mlpet.AerobicGlycolysisKit
            
            ALPHA = 0.005670305; % log(2)/halflife in 1/s
            
            f = ks(1);
            lambda = ks(3); 
            Delta = ks(4);
            E = 1;
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
            
            import mloxygen.ConstantERaichle1983Model.solution
            qs = solution(ks, artery_interpolated);
            idx_sampled = floor(times_sampled - times_sampled(1)) + 1;
            idx_sampled = idx_sampled(idx_sampled <= length(qs));
            qs = qs(idx_sampled);
        end
        function loss = simulanneal_objective(ks, artery_interpolated, times_sampled, qs0, sigma0)
            import mloxygen.ConstantERaichle1983Model.sampled          
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
            m('k1') = struct('min', 0.0043, 'max', 0.0155, 'init', 0.00777, 'sigma', 3.89e-4); % f / s
            m('k2') = struct('min', 0.0137, 'max', 0.0266, 'init', 0.0228,  'sigma', 0.002); % PS / s
            m('k3') = struct('min', 0.608,  'max', 1.06,   'init', 0.945,   'sigma', 0.05); % lambda in mL/mL
            m('k4') = struct('min', 0.2,    'max', 3,      'init', 2,       'sigma', 0.1); % Delta for cerebral dispersion
        end
        function loss = loss_function(ks, artery_interpolated, times_sampled, measurement, sigma0)
            import mloxygen.ConstantERaichle1983Model.sampled            
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
    end
    
	methods 		  
 		function this = ConstantERaichle1983Model(varargin)
            
            this = this@mloxygen.Raichle1983Model(varargin{:});
            
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
                    this.map('k2') = struct('min', 0.017,    'max', 0.0266,  'init', 0.0218,  'sigma', 0.002); % PS / s
                    this.map('k3') = struct('min', 0.738,    'max', 1.06,    'init', 1.02,    'sigma', 0.05); % lambda in mL/mL  
                    %this.map('k3') = struct('min', 0.987,    'max', 1.06,    'init', 1.02,    'sigma', 0.05); % lambda in mL/mL                    
                case 'w'
                    this.map('k2') = struct('min', 0.0137,   'max', 0.0142,  'init', 0.014,   'sigma', 0.002); % PS / s
                    this.map('k3') = struct('min', 0.608,    'max', 0.882,   'init', 0.851,   'sigma', 0.05); % lambda in mL/mL
                    %this.map('k3') = struct('min', 0.819,    'max', 0.882,   'init', 0.851,   'sigma', 0.05); % lambda in mL/mL
                case 's' % subcortical
                    this.map('k2') = struct('min', 0.0159,   'max', 0.0215,  'init', 0.0187,  'sigma', 0.002); % PS / s
                    this.map('k3') = struct('min', 0.738,    'max', 0.97,    'init', 0.924,   'sigma', 0.05); % lambda in mL/mL
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
            ho = mloxygen.DispersedRaichle1983Model.sampled(ks, aif, this.times_sampled);
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

