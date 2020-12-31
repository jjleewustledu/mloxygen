classdef (Abstract) Raichle1983Model < mlpet.TracerKineticsModel
	%% RAICHLE1983MODEL provides model data and methods to the strategy design pattern comprising
    %  mloxygen.{Raichle1983, Raichle1983Nest, Raichle1983SimulAnneal, Raichle1983HMC, Raichle1983 LM, Raichle1983BFGS}.
    %  It operates on single voxels or regions.

	%  $Revision$
 	%  was created 10-Sep-2020 19:42:35 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.
    
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
    end

	methods		  
 		function this = Raichle1983Model(varargin) 			
            this = this@mlpet.TracerKineticsModel(varargin{:});
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

