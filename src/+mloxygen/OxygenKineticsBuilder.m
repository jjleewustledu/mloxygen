classdef OxygenKineticsBuilder < mlkinetics.AbstractKineticsBuilder 
	%% OxygenKineticsBuilder  

	%  $Revision$
 	%  was created 04-Dec-2017 13:53:52 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlhemodynamics/src/+mlhemodynamics.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
	properties 		
 	end

	methods 
		  
        %%
        
        function cbf = buildCbf(this)
            this.solver_ = this.solver_.estimateParameters('cbf');
            assert(this.solver_.isfinished);
            cbf = this.solver_.model.cbf;
        end
        function cbv = buildCbv(this)
            this.solver_ = this.solver_.estimateParameters('cbv');
            assert(this.solver_.isfinished);
            cbv = this.solver_.model.cbv;
        end
        function oef = buildOef(this, cbf, cbv)
            this.solver_.model.cbf = cbf;
            this.solver_.model.cbv = cbv;
            this.solver_ = this.solver_.estimateParameters('oef');
            assert(this.solver_.isfinished);
            oef = this.solver_.model.oef;
        end
        function cmro2 = buildCmro2(this, cbf, oef)
            this.solver_.model.cbf = cbf;
            this.solver_.model.oef = oef;
            this.solver_ = this.solver_.estimateParameters('cmro2');
            assert(this.solver_.isfinished);
            cmro2 = this.solver_.model.cmro2;
        end
        function waterMetab = buildWaterMetab(this)
            this.solver_ = this.solver_.estimateParameters('waterMetab');
            assert(this.solver_.isfinished);
            waterMetab = this.solver_.model.waterMetab;
        end
        
 		function this = OxygenKineticsBuilder(varargin)
 			%% OxygenKineticsBuilder
 			%  @param named solver is an mlanalysis.ISolver
            
            this = this@mlkinetics.AbstractKineticsBuilder(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

