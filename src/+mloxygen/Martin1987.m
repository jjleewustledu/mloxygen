classdef Martin1987 < handle
	%% Martin1987  

	%  $Revision$
 	%  was created 31-Oct-2018 15:20:13 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
    properties (Constant)
        DENSITY_BRAIN = 1.05 % assumed mean brain density, g/mL
        RATIO_SMALL_LARGE_HCT = 0.85 % Grubb, et al., 1978
    end
    
	properties (Dependent)
        T0 % integration limits
        Tf % integration limits
 		aif
        ispercent
        tac
        W % := 1; aif and tac manage their own calibrations
        
        subject
        session
        scan
 	end

	methods 
        
        %% GET
        
        function g = get.T0(this)
            g = this.T0_;
        end
        function g = get.Tf(this)
            g = this.Tf_;
        end
        function g = get.aif(this)
            g = this.aif_;
        end
        function g = get.ispercent(this)
            g = this.ispercent_;
        end
        function g = get.tac(this)
            g = this.tac_;
        end
        function g = get.W(~)
            g = 1;
        end
        
        function g = get.subject(this)
            g = this.xnatinfo_.subject;
        end
        function g = get.session(this)
            g = this.xnatinfo_.session;
        end
        function g = get.scan(this)
            g = this.xnatinfo_.scan;
        end
        
        %%
        
        function dt = datetime(this)
            dt = datetime(this.tac);
            diff = seconds(dt - datetime(this.aif));
            assert(abs(diff) < seconds(10), ...
                'mloxygen:ValueError', 'Martin1987.datetime:  difference of datetimes TAC - AIF => %g sec', diff);
        end
        function e = estimate(this)
            assert(~this.tac.isDecayCorrected);
            assert(~this.aif.isDecayCorrected);

            scale = this.W / (this.RATIO_SMALL_LARGE_HCT * this.DENSITY_BRAIN);
            if (this.ispercent)
                scale = 100*scale;
            end
            
            e = scale * trapz(this.tac, this.T0, this.Tf) / trapz(this.aif, this.T0, this.Tf);
        end
        function plotWorldlines(this, varargin)
            a = this.aif;
            t = this.tac;
            plot(a.times, a.activities, t.times, t.activities, varargin{:});
            xlabel('t / sec');
            ylabel('activity / (Bq/mL)');
            legend('AIF', 'TAC');
            title(this);
        end
        function str = title(this)
            str = sprintf('Martin1987 CBV %s %s %s %s', ...
                char(this.subject), char(this.session), char(this.scan), datetime(this));
        end
		  
 		function this = Martin1987(varargin)
 			%% Martin1987
 			%  @param T0 is numeric.
            %  @param Tf is numeric.
            %  @param aif is mlpet.IAifData.
            %  @param ispercent is logical.
            %  @param tac is mlpet.IScannerData.
            %  @param xnatinfo is mlxnat.<>.

            ip = inputParser;
            addParameter(ip, 'T0', nan, @isnumeric);
            addParameter(ip, 'Tf', nan, @isnumeric);
            addParameter(ip, 'ispercent', true, @islogical);
            addParameter(ip, 'aif', [], @(x) isa(x, 'mlpet.IAifData'));
            addParameter(ip, 'tac', [], @(x) isa(x, 'mlpet.IScannerData'));
            addParameter(ip, 'xnatinfo', [], @(x) isa(x, 'mlxnat.<>'));
            parse(ip, varargin{:});
 			
            this.T0_        = ip.Results.T0;
            this.Tf_        = ip.Results.Tf;
            this.aif_       = ip.Results.aif;
            this.ispercent_ = ip.Results.ispercent;
            this.tac_       = ip.Results.tac;
            this.xnatinfo_  = ip.Results.xnatinfo;
 		end
    end 
    
    %% PRIVATE
    
    properties (Access = private)
        T0_
        Tf_
        aif_
        ispercent_
        tac_
        xnatinfo_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

