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
        function loss = loss_function(ks, artery_interpolated, times_sampled, measurement, sigma0)
            import mloxygen.ConstantERaichle1983Model.sampled  
            estimation  = sampled(ks, artery_interpolated, times_sampled);
            measurement = measurement(1:length(estimation));
            positive    = measurement > 0.05*max(measurement);            
            eoverm      = estimation(positive)./measurement(positive);
            Q           = mean(abs(1 - eoverm));
            %Q           = sum((1 - eoverm).^2);
            loss        = 0.5*Q; % /sigma0^2; % + sum(log(sigma0*measurement)); % sigma ~ sigma0*measurement
        end  
        function qs   = sampled(ks, artery_interpolated, times_sampled)
            %  @param artery_interpolated is uniformly sampled at high sampling freq.
            %  @param times_sampled are samples scheduled by the time-resolved PET reconstruction
            
            import mloxygen.ConstantERaichle1983Model.solution 
            import mlpet.TracerKineticsModel.solutionOnScannerFrames  
            qs = solution(ks, artery_interpolated);
            qs = solutionOnScannerFrames(qs, times_sampled);
        end
        function loss = simulanneal_objective(ks, artery_interpolated, times_sampled, qs0, sigma0)
            import mloxygen.ConstantERaichle1983Model.sampled          
            qs = sampled(ks, artery_interpolated, times_sampled);            
            loss = 0.5 * sum((1 - qs ./ qs0).^2) / sigma0^2; % + sum(log(sigma0*qs0)); % sigma ~ sigma0 * qs0
        end
        function qs   = solution(ks, artery_interpolated)
            %  @param artery_interpolated is uniformly with at high sampling freq. starting at time = 0.
            
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
    end
    
	methods 		  
 		function this = ConstantERaichle1983Model(varargin)
            %  @param histology is:  'g', 'w', 's', else histology information is not used.
            
            this = this@mloxygen.Raichle1983Model(varargin{:});
        end        
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

