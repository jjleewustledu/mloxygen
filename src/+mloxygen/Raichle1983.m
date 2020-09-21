classdef Raichle1983 < handle & mlpet.TracerKinetics
	%% RAICHLE1983 is the context to a strategy design patterns which implements:
    %  mloxygen.{Raichle1983Nest, Raichle1983SimulAnneal, Raichle1983HMC, Raichle1983LM, Raichle1983BFGS}.
    %  For performance considerations, see also https://blogs.mathworks.com/loren/2012/03/26/considering-performance-in-object-oriented-matlab-code/

	%  $Revision$
 	%  was created 30-May-2018 01:54:09 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlpet/src/+mlpet.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	    
    properties        
        Dt          % time-shift for AIF; Dt < 0 shifts backwards in time.
        measurement % expose for performance when used by mloxygen.Raichle1983Strategy
        model       %
        regionTag 		
    end   
    
	properties (Dependent)
        averageVoxels
    end 
    
    methods (Static)
        function this = createFromDeviceKit(devkit, varargin)
            sesd = devkit.sessionData;
            sesd.jitOn222(sesd.ocOnAtlas())
            this = mloxygen.Raichle1983('devkit', devkit, varargin{:});
        end
    end
    
    methods
        
        %% GET
        
        function g = get.averageVoxels(this)
            g = this.averageVoxels_;
        end
        function     set.averageVoxels(this, s)
            assert(islogical(s))
            this.averageVoxels_ = s;
        end
        
        %%
        
        function ks = buildKs(this, varargin)
            this = solve(this, varargin{:});
            ks = [k1(this) k2(this) k3(this)];
        end
        function [k,sk] = k1(this, varargin)
            [k,sk] = k1(this.strategy_, varargin{:});
        end
        function [k,sk] = k2(this, varargin)
            [k,sk] = k2(this.strategy_, varargin{:});
        end
        function [k,sk] = k3(this, varargin)
            [k,sk] = k3(this.strategy_, varargin{:});
        end
        function [k,sk] = ks(this, varargin)
            k = zeros(1,3);
            sk = zeros(1,3);
            [k(1),sk(1)] = k1(this.strategy_, varargin{:});
            [k(2),sk(2)] = k2(this.strategy_, varargin{:});
            [k(3),sk(3)] = k3(this.strategy_, varargin{:});
        end
        function h = plot(this, varargin)
            h = this.strategy_.plot(varargin{:});
        end
        function this = simulated(this, varargin)
            this.measurement = this.model.simulated(varargin{:});
            this.strategy_.Measurement = this.measurement; % strategy_ needs value copies for performance
        end
        function this = solve(this, varargin)
            this.strategy_ = solve(this.strategy_, varargin{:});
        end
    end
    
    %% PROTECTED
    
    properties (Access = protected)
        averageVoxels_
        strategy_
    end

	methods (Access = protected)		  
 		function this = Raichle1983(varargin)
            %% RAICHLE1983
            %  @param devkit is mlpet.IDeviceKit.
            %  @param Dt is numeric, s of time-shifting for AIF.
            %  @param averageVoxels is logical, choosing creation of scalar results.
            
            this = this.mlpet.TracerKinetics(varargin{:});
            
            ip = inputParser;            
            ip.KeepUnmatched = true;
            addParameter(ip, 'Dt', 0, @isscalar)
            addParameter(ip, 'averageVoxels', false, @islogical);
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            
            this.Dt = ipr.Dt;
            this.averageVoxels_ = ipr.averageVoxels;
            this.regionTag = this.devkit_.sessionData.regionTag;
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

