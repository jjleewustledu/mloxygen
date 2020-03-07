classdef MetabDirector < handle & mlkinetics.MetabDirector
	%% METABDIRECTOR  

	%  $Revision$
 	%  was created 16-Dec-2018 23:43:17 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
	properties
 		
    end
    
    methods (Static)
        function this = constructRaichle(varargin)
            this = mloxygen.MetabDirector( ...
                mloxygen.RaichleBuilder(varargin{:}));
        end
        function this = constructMartin(varargin)
            this = mloxygen.MetabDirector( ...
                mloxygen.MartinBuilder(varargin{:}));
        end
        function this = constructMintun(varargin)
            this = mloxygen.MetabDirector( ...
                mloxygen.MintunBuilder(varargin{:}));
        end
        function this = constructFokkerPlanck(varargin)
            this = mloxygen.MetabDirector( ...
                mloxygen.FokkerPlanckBuilder(varargin{:}));
        end
        function this = constructKPZ(varargin)
            this = mloxygen.MetabDirector( ...
                mloxygen.KPZBuilder(varargin{:}));
        end
    end

	methods		  
 		function this = MetabDirector(varargin)
 			%% METABDIRECTOR
 			%  @param required builder is mloxygen.MetabBuilder.

 			this = this@mlkinetics.MetabDirector(varargin{:});	
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

