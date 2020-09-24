classdef (Abstract) Raichle1983Model 
	%% RAICHLE1983MODEL provides model data and methods to the strategy design pattern comprising
    %  mloxygen.{Raichle1983, Raichle1983Nest, Raichle1983SimulAnneal, Raichle1983HMC, Raichle1983 LM, Raichle1983BFGS}.
    %  It operates on single voxels or regions.

	%  $Revision$
 	%  was created 10-Sep-2020 19:42:35 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    methods (Abstract, Static)
        preferredMap
        sampled
        simulanneal_objective
        solution
    end
    
    methods (Abstract)
        simulated(this)
    end
    
	properties
        artery_interpolated
 		map
        times_sampled 		
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
    end

	methods		  
 		function this = Raichle1983Model(varargin)
            %  @param times_sampled for scanner is typically not uniform
            %  @param artery_interpolated must be uniformly interpolated
 			
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'map', this.preferredMap(), @(x) isa(x, 'containers.Map'))
            addParameter(ip, 'times_sampled', [], @isnumeric)
            addParameter(ip, 'artery_interpolated', [], @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.map = ipr.map;
            this = set_times_sampled(this, ipr.times_sampled);
            this = set_artery_interpolated(this, ipr.artery_interpolated);
        end
        function this = set_times_sampled(this, s)
            if isempty(s)
                return
            end
            this.times_sampled = s;
        end
        function this = set_artery_interpolated(this, s)
            if isempty(s)
                return
            end
            % artery_interpolated may be shorter than scanner times_sampled
            assert(~isempty(this.times_sampled))
            if length(s)-1 ~= this.times_sampled(end)
                this.artery_interpolated = ...
                    makima(0:length(s)-1, s, this.times_sampled(1):this.times_sampled(end));
            else
                this.artery_interpolated = s;
            end
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

