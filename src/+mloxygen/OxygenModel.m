classdef OxygenModel < mlkinetics.KineticsModel
	%% OXYGENMODEL  

	%  $Revision$
 	%  was created 13-Dec-2017 17:43:25 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
    properties (Dependent)
        cbf
        cbv
    end
    
	properties 		
        oef
        cmro2
        waterMetab
 	end

	methods 
        
        %% GET
        
        function g = get.cbf(this)
            g = this.hemodynamicsModel_.cbf;
        end
        function g = get.cbv(this)
            g = this.hemodynamicsModel_.cbv;
        end
        
        %%
		  
        function prms = solverParameters(this)
            import mloxygen.*;
            switch (class(this.solver_))
                case 'mlanalysis.LevenbergMarquardt'
                    prms = OxygenLMParameters;
                case 'mlbayesian.BretthorstMcmc'
                    prms = OxygenMcmcParameters;
                case 'mlnest.NestedSamplingMain'
                    prms = OxygenNestParameters;
                case 'mlstan.FlatHMC'
                    prms = OxygenFlatHMCParameters;
                case 'mlstan.HierarchicalHMC'
                    prms = OxygenHierarchicalHMCParameters;
                otherwise
                    error('mlraichle:unsupportedSwitchStrategy', ...
                        'OxygenModel.solverParameters:class(solver)->%s', this.solver_);
            end
        end 
		  
 		function this = OxygenModel(varargin)
 			%% OXYGENMODEL

 			this = this@mlkinetics.KineticsModel(varargin{:});
            this.hemodynamicsModel_ = mlhemodynamics.HemodynamicsModel(varargin{:});
 		end
    end 
    
    %% PRIVATE
    
    properties
        hemodynamicsModel_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

