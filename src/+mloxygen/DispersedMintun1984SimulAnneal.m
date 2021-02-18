classdef DispersedMintun1984SimulAnneal < mloxygen.Mintun1984SimulAnneal
	%% DISPERSEDMINTUN1984SIMULANNEAL  

	%  $Revision$
 	%  was created 07-Dec-2020 23:20:53 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.9.0.1524771 (R2020b) Update 2 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    properties
        Dt
    end
    
	methods 		  
 		function this = DispersedMintun1984SimulAnneal(varargin)
 			this = this@mloxygen.Mintun1984SimulAnneal(varargin{:});
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            ip.PartialMatching = false;
            addParameter(ip, 'Dt', [], @isnumeric)
            parse(ip, varargin{:})
            this.Dt = ip.Results.Dt;
        end
        
        function fprintfModel(this)
            fprintf('Simulated Annealing:\n');            
            for ky = 1:length(this.ks)
                fprintf('\tk%i = %f\n', ky, this.ks(ky));
            end     
            fprintf('\tDt = %f\n', this.Dt);       
            
            fs = this.model.fs_Raichle_Martin;
            for ky = 1:length(fs)-1
                fprintf('\tfs(%i) = %f\n', ky, fs(ky));
            end
            fprintf('\tv1 = %f\n', fs(end));
            fprintf('\tE = 1 - exp(-PS/f) = %f\n', 1 - exp(-fs(2)/fs(1)));
            
            fprintf('\tsigma0 = %f\n', this.sigma0);
            for ky = this.map.keys
                fprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})));
            end
        end
        function s = sprintfModel(this)
            s = sprintf('Simulated Annealing:\n');
            for ky = 1:length(this.ks)
                s = [s sprintf('\tk%i = %f\n', ky, this.ks(ky))]; %#ok<AGROW>
            end
            s = [s sprintf('\tDt = %f\n', this.Dt)];
            
            fs = this.model.fs_Raichle_Martin;
            for ky = 1:length(fs)-1
                s = [s sprintf('\tfs(%i) = %f\n', ky, fs(ky))]; %#ok<AGROW>
            end            
            s = [s sprintf('\tv1 = %f\n', fs(end))];            
            s = [s sprintf('\tE = 1 - exp(-PS/f) = %f\n', 1 - exp(-fs(2)/fs(1)))];
            
            s = [s sprintf('\tsigma0 = %f\n', this.sigma0)];
            for ky = this.map.keys
                s = [s sprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})))]; %#ok<AGROW>
            end
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

