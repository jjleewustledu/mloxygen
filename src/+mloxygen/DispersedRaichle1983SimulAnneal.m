classdef DispersedRaichle1983SimulAnneal < mloxygen.Raichle1983SimulAnneal
	%% DISPERSEDRAICHLE1983SIMULANNEAL  

	%  $Revision$
 	%  was created 10-Sep-2020 22:24:15 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.

    properties
        Dt
    end
    
	methods		  
 		function this = DispersedRaichle1983SimulAnneal(varargin)
 			this = this@mloxygen.Raichle1983SimulAnneal(varargin{:});
            
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
            fprintf('\tE = 1 - exp(-PS/f) = %f\n', 1 - exp(-this.ks(2)/this.ks(1)))
            fprintf('\tsigma0 = %f\n', this.sigma0);
            fprintf('\tDt = %f\n', this.Dt);
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
            s = [s sprintf('\tE = 1 - exp(-PS/f) = %f\n', 1 - exp(-this.ks(2)/this.ks(1)))];
            s = [s sprintf('\tsigma0 = %f\n', this.sigma0)];
            for ky = this.map.keys
                s = [s sprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})))]; %#ok<AGROW>
            end
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

