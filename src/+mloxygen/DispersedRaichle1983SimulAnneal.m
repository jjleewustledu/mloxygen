classdef DispersedRaichle1983SimulAnneal < mloxygen.Raichle1983SimulAnneal
	%% DISPERSEDRAICHLE1983SIMULANNEAL  

	%  $Revision$
 	%  was created 10-Sep-2020 22:24:15 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.

    properties
        Dt
        registry
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
            this.registry = mlraichle.RaichleRegistry.instance();
        end    
        
        function fprintfModel(this)
            fprintf('Simulated Annealing:\n');            
            for ky = 1:length(this.ks)
                fprintf('\tk%i = %g\n', ky, this.ks(ky));
            end 
            fprintf('\tE = 1 - exp(-PS/f) = %g\n', 1 - exp(-this.ks(2)/this.ks(1)))          
            fprintf('\tDt = %g\n', this.Dt);
            fprintf('\ttBuffer = %g\n', this.registry.tBuffer)
            fprintf('\ttimeCliff = %g\n', this.timeCliff)
            fprintf('\tloss = %g\n', this.loss())
            fprintf('\tsigma0 = %g\n', this.sigma0);
            for ky = this.map.keys
                fprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})));
            end
        end
        function s = sprintfModel(this)
            s = sprintf('Simulated Annealing:\n');
            for ky = 1:length(this.ks)
                s = [s sprintf('\tk%i = %g\n', ky, this.ks(ky))]; %#ok<AGROW>
            end
            s = [s sprintf('\tE = 1 - exp(-PS/f) = %g\n', 1 - exp(-this.ks(2)/this.ks(1)))];
            s = [s sprintf('\tDt = %g\n', this.Dt)];
            s = [s sprintf('\ttBuffer = %g\n', this.registry.tBuffer)];
            s = [s sprintf('\ttimeCliff = %g\n', this.timeCliff)];
            s = [s sprintf('\tloss = %g\n', this.loss())];
            s = [s sprintf('\tsigma0 = %g\n', this.sigma0)];
            for ky = this.map.keys
                s = [s sprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})))]; %#ok<AGROW>
            end
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

