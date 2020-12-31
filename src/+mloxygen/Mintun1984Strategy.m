classdef (Abstract) Mintun1984Strategy 
	%% MINTUN1984STRATEGY is the interface for concrete strategies used with mloxygen.Mintun1984.

	%  $Revision$
 	%  was created 07-Dec-2020 23:28:15 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.9.0.1524771 (R2020b) Update 2 for MACI64.  Copyright 2020 John Joowon Lee.
 	
	properties (Abstract)
        artery_interpolated
        context
        fileprefix
        map                  % containers.Map containing model params as structs with fields:  min, max, init
        model                %
        Measurement          % external data
        results
        sigma0               % fraction of Measurement < 1
        times_sampled        % numeric times for Measurement; midpoints of frames for PET
        visualize 		
 	end

	methods (Abstract)
        k1(this)
        k2(this)
        solve(this)
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

