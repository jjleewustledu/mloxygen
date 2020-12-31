classdef (Abstract) Mintun1984Model < mlpet.TracerKineticsModel
	%% MINTUN1984MODEL  

	%  $Revision$
 	%  was created 05-Dec-2020 20:47:11 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.9.0.1524771 (R2020b) Update 2 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    properties (Constant)
        NOMINAL_O2_CONTENT = 0.17 % mL O2.  Estimated from Ito, Eur J Nucl Med Mol Imaging (2004) 31:635-643.
                                  % DOI 10.1007/s00259-003-1430-8
    end
    
    properties
        fs_Raichle_Martin % [fs PS lambda Dt v1] | [f PS lambda Delta Dt v1]
    end
    
	methods 		  
 		function this = Mintun1984Model(varargin)
            %  @param fs_Raichle_Martin ~ [fs PS lambda Dt v1] | [f PS lambda Delta Dt v1].
            
 			this = this@mlpet.TracerKineticsModel(varargin{:});
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'fs_Raichle_Martin', [], @isnumeric)
            parse(ip, varargin{:})
            
            this.fs_Raichle_Martin = ip.Results.fs_Raichle_Martin;
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end
