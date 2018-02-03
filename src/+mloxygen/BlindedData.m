classdef BlindedData 
	%% BLINDEDDATA  

	%  $Revision$
 	%  was created 14-Jan-2018 18:10:06 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2018 John Joowon Lee.
 	
    properties (Dependent)
        meanArterialPressure
    end
    
	properties
 		heartrate
        bloodPressure % [systolic diastolic]
 	end

	methods 
        
        %% GET
        
        function g = get.meanArterialPressure(this)
            g = this.bloodPressure(1)/3 + this.bloodPressure(2)*2/3;
        end
		  
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

