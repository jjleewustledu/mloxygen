classdef Chaojie2020 < mloxygen.Herscovitch1985
	%% CHAOJIE2020  

	%  $Revision$
 	%  was created 24-Jul-2020 10:32:57 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.8.0.1417392 (R2020a) Update 4 for MACI64.  Copyright 2020 John Joowon Lee.
 	
	properties
 		
 	end

	methods 
		  
 		function this = Chaojie2020(varargin)
 			%% CHAOJIE2020
 			%  @param .

 			this = this@mloxygen.Herscovitch1985(varargin{:});
 		end
    end 
    
    methods (Access = protected)        
        function this = readdata(this, tra)
            %% assumes that tac acquisition begins with tac.Start(1) == 0
            
            cache = readtable(fullfile(this.sessionPath, sprintf('%s%s1%s.dcv', this.pnumber, tra, this.studyTag)), 'FileType', 'text');
            msk = ~isnan(cache.Var2);
            dcv = table(cache.Var1(msk), cache.Var2(msk), cache.Var2(msk), 'VariableNames', {'time', 'dcv', 'dcv_well'});
            
            try
                tac = readtable(fullfile(this.sessionPath, sprintf('%s_%s1_frame_tac.csv', this.pnumber, tra)));
            catch ME
                handwarning(ME)
                tac = readtable(fullfile(this.sessionPath, sprintf('%s%s1_frame_tac.csv', this.pnumber, tra)));                
            end
            this.startTime = tac.Start(1);
            tac.Start = tac.Start - this.startTime;
            tac.End = tac.End - this.startTime;
            prop = [tra 'data'];
            this.(prop) = struct('dcv', dcv, 'tac', tac);
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

