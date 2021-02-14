classdef DispersedMartin1987Model < mlpet.TracerKineticsModel
	%% DISPERSEDMARTIN1987MODEL  

	%  $Revision$
 	%  was created 10-Feb-2021 17:27:21 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.9.0.1570001 (R2020b) Update 4 for MACI64.  Copyright 2021 John Joowon Lee.
 	
    properties 
        T0
        Tf
    end
    
    methods (Static)
        function loss = loss_function(varargin)
            loss = [];
        end
        function m = preferredMap()
            m = [];
        end
        function qs = sampled(varargin)
            qs = [];
        end
        function loss = simulanneal_objective(varargin)
            loss = [];
        end
        function qs = solution(varargin)
            qs = [];
        end
    end

	methods 		  
 		function this = DispersedMartin1987Model(varargin)
 			this = this@mlpet.TracerKineticsModel(varargin{:});
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            ip.PartialMatching = false;
            addParameter(ip, 'T0', 120, @isfinite);
            addParameter(ip, 'Tf', 240, @isfinite);
            parse(ip, varargin{:});
            ipr = ip.Results;
            this.T0 = ipr.T0;
            this.Tf = ipr.Tf; 
        end
        function oc = simulated(~, varargin)
            oc = [];
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

