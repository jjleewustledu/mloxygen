classdef OxygenKineticsDirector 
	%% OXYGENKINETICSDIRECTOR  

	%  $Revision$
 	%  was created 04-Dec-2017 13:53:35 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/Local/src/mlcvl/mlhemodynamics/src/+mlhemodynamics.
 	%% It was developed on Matlab 9.3.0.713579 (R2017b) for MACI64.  Copyright 2017 John Joowon Lee.
 	
    properties (Dependent)
        product
    end
    
	methods 
        
        %% GET/SET
        
        function g = get.product(this)
            g = mloxygen.OxygenView( ...
                'cbf', this.cbf_, ...
                'cbv', this.cbv_, ...
                'oef', this.oef_, ...
                'cmro2', this.cmro2_, ...
                'waterMetab', this.waterMetab_);
        end
        
        %%
        		  
        function this = constructRates(this)
            this = this.constructCbf;
            this = this.constructCbv;
            this = this.constructOef;
        end
        function this = constructPhysiological(this)
            this = this.constructCmro2;
            this = this.constructWaterMetab;
        end
        
 		function this = OxygenKineticsDirector(varargin)
 			%% OXYGENKINETICSDIRECTOR
 			%  @param oxygenBldr is an mlkinetics.IKineticsBuilder.
            %  @param named model is an mlhemodynamics.HemodynamicsModel.

 			ip = inputParser;
            addRequired( ip, 'oxygenBldr', @(x) isa(x, 'mloxygen.OxygenKineticsBuilder'));
            addParameter(ip, 'roisBldr', @(x) isa(x, 'mlrois.IRoisBuilder'));
            parse(ip, varargin{:});
            
            this.kineticsBuilder_ = ip.Results.oxygenBldr;	
            this.roisBuilder_ = ip.Results.roisBldr;
 		end
 	end 

    %% PRIVATE
    
    properties (Access = private)
        cbf_
        cbv_
        cmro2_
        oef_
        waterMetab_
    end
    
    methods (Access = private)        
        function this = constructCbf(this)
            this.cbf_ = this.kineticsBuilder_.buildCbf;
        end
        function this = constructCbv(this)
            this.cbv_ = this.kineticsBuilder_.buildCbv;
        end
        function this = constructOef(this)
            if (isempty(this.cbf_))
                this.cbf_ = this.constructCbf;
            end
            if (isempty(this.cbv_))
                this.cbv_ = this.constructCbv;
            end
            this.oef_ = this.kineticsBuilder_.buildOef(this.cbf_, this.cbv_);
        end
        function this = constructCmro2(this)
            if (isempty(this.cbf_))
                this.cbf_ = this.constructCbf;
            end
            if (isempty(this.oef_))
                this.oef_ = this.constructOef;
            end
            this.cmro2_ = this.kineticsBuilder_.buildCmro2(this.cbf_, this.oef_);
        end
        function this = constructWaterMetab(this)
            this.waterMetab_ = this.kineticsBuilder_.buildWaterMetab;
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

