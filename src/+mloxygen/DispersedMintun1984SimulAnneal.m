classdef DispersedMintun1984SimulAnneal < mloxygen.Mintun1984SimulAnneal
	%% DISPERSEDMINTUN1984SIMULANNEAL  

	%  $Revision$
 	%  was created 07-Dec-2020 23:20:53 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.9.0.1524771 (R2020b) Update 2 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    properties
        Dt
        DtMixing
        fracMixing
        aifdata
    end
    
	methods 		  
 		function this = DispersedMintun1984SimulAnneal(varargin)
 			this = this@mloxygen.Mintun1984SimulAnneal(varargin{:});
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            ip.PartialMatching = false;
            addParameter(ip, 'Dt', [], @isnumeric)
            addParameter(ip, 'DtMixing', [], @isnumeric)
            addParameter(ip, 'fracMixing', [], @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.Dt = ipr.Dt;
            this.DtMixing = ipr.DtMixing;
            this.fracMixing = ipr.fracMixing;
            this.aifdata = mlaif.AifData.instance();
        end
        
        function fprintfModel(this)
            fprintf('Simulated Annealing:\n');            
            for ky = 1:length(this.ks)
                fprintf('\tk%i = %g\n', ky, this.ks(ky));
            end     
            fprintf('\tDt = %f\n', this.Dt); 
            fprintf('\ttBuffer = %g\n', this.aifdata.tBuffer)  
            fprintf('\tDtMixing = %g\n', this.DtMixing)
            fprintf('\tfracMixing = %g\n', this.fracMixing)
            %fprintf('\ttimeCliff = %g\n', this.timeCliff)
            fprintf('\tloss = %g\n', this.loss())      
            
            fs = this.model.fs_Raichle_Martin;
            for ky = 1:length(fs)-1
                fprintf('\tfs(%i) = %g\n', ky, fs(ky));
            end
            fprintf('\tE = 1 - exp(-PS/f) = %g\n', 1 - exp(-fs(2)/fs(1)));
            fprintf('\tv1 = %g\n', fs(end));
            
            %fprintf('\tsigma0 = %g\n', this.sigma0);
            for ky = this.map.keys
                fprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})));
            end
        end
        function s = sprintfModel(this)
            s = sprintf('Simulated Annealing:\n');
            for ky = 1:length(this.ks)
                s = [s sprintf('\tk%i = %g\n', ky, this.ks(ky))]; %#ok<AGROW>
            end
            s = [s sprintf('\tDt = %g\n', this.Dt)];
            s = [s sprintf('\ttBuffer = %g\n', this.aifdata.tBuffer)];            
            s = [s sprintf('\tDtMixing = %g\n', this.DtMixing)];
            s = [s sprintf('\tfracMixing = %g\n', this.fracMixing)];
            %s = [s sprintf('\ttimeCliff = %g\n', this.timeCliff)];
            s = [s sprintf('\tloss = %g\n', this.loss())];
            
            fs = this.model.fs_Raichle_Martin;
            for ky = 1:length(fs)-1
                s = [s sprintf('\tfs(%i) = %g\n', ky, fs(ky))]; %#ok<AGROW>
            end                      
            s = [s sprintf('\tE = 1 - exp(-PS/f) = %g\n', 1 - exp(-fs(2)/fs(1)))];
            s = [s sprintf('\tv1 = %g\n', fs(end))];  
            
            %s = [s sprintf('\tsigma0 = %g\n', this.sigma0)];
            for ky = this.map.keys
                s = [s sprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})))]; %#ok<AGROW>
            end
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

