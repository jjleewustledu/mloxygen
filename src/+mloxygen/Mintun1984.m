classdef Mintun1984 < handle & mlpet.TracerKineticsStrategy
	%% MINTUN1984 is the context to a strategy design patterns which implements:
    %  mloxygen.{Mintun1984SimulAnneal,Mintun1984LM}.
    %  For performance considerations, see also https://blogs.mathworks.com/loren/2012/03/26/considering-performance-in-object-oriented-matlab-code/

	%  $Revision$
 	%  was created 30-May-2018 01:53:47 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlpet/src/+mlpet.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.

    methods (Static)
        function this = createFromDeviceKit(devkit, varargin)
            sesd = devkit.sessionData;
            sesd.jitOn222(sesd.ooOnAtlas())
            this = mloxygen.Mintun1984('devkit', devkit, varargin{:});
        end
    end
    
    methods
        function ks = buildKs(this, varargin)
            this = solve(this, varargin{:});
            ks = [k1(this) k2(this) k3(this)];
        end
        function [k,sk] = k1(this, varargin)
            %% k1 == oef \in [0 1]
            
            [k,sk] = k1(this.strategy_, varargin{:});
        end
        function [k,sk] = k2(this, varargin)
            %% k2 == metab rate, 1/s
            
            [k,sk] = k2(this.strategy_, varargin{:});
        end
        function [k,sk] = k3(this, varargin)
            %% k3 == metab amplitude, unit-less
            
            [k,sk] = k3(this.strategy_, varargin{:});
        end
        function [k,sk] = k4(this, varargin)
            %% k4 == Delta, for dispersion of aif
            
            [k,sk] = k4(this.strategy_, varargin{:});
        end
        function ks_ = ks(this, varargin)
            %% ks == [oef meta-rate meta-amplitude]
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'typ', 'mlfourd.ImagingContext2', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            k(1) = k1(this.strategy_, varargin{:});
            k(2) = k2(this.strategy_, varargin{:});
            k(3) = k3(this.strategy_, varargin{:});
            k(4) = k4(this.strategy_, varargin{:});
             
            roibin = logical(this.roi);
            ks_ = copy(this.roi.fourdfp);
            ks_.img = zeros([size(this.roi) length(k)]);
            for t = 1:length(k)
                img = zeros(size(this.roi), 'single');
                img(roibin) = k(t);
                ks_.img(:,:,:,t) = img;
            end
            ks_.fileprefix = this.sessionData.osOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
            ks_ = imagingType(ipr.typ, ks_);
        end 
        function os_ = os(this, varargin)
            os_ = this.ks(varargin{:});
        end   
    end
    
    %% PROTECTED

	methods (Access = protected)		  
 		function this = Mintun1984(varargin)
            this = this@mlpet.TracerKineticsStrategy(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

