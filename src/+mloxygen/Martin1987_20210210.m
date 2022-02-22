classdef Martin1987_20210210 < handle & mlpet.TracerKineticsStrategy
	%% MARTIN1987_20210210  

	%  $Revision$
 	%  was created 10-Feb-2021 18:39:26 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.9.0.1570001 (R2020b) Update 4 for MACI64.  Copyright 2021 John Joowon Lee. 	
	
    methods (Static)
        function this = createFromDeviceKit(devkit, varargin)
            sesd = devkit.sessionData;
            sesd.jitOnAtlas(sesd.ocOnAtlas())
            this = mloxygen.Martin1987_20210210('devkit', devkit, varargin{:});
        end
    end

	methods         
        function ks = buildKs(this, varargin)
            this = solve(this, varargin{:});
            ks = k1(this);
        end 
        function [k,sk] = k1(this, varargin)
            [k,sk] = k1(this.strategy_, varargin{:});
        end
        function ks_ = ks(this, varargin)
            %% ks == v1
            %  @param 'typ' is char, understood by imagingType.            
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'typ', 'mlfourd.ImagingContext2', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            k(1) = k1(this.strategy_, varargin{:});
             
            roibin = logical(this.roi);
            ks_ = copy(this.roi.fourdfp);
            ks_.img = zeros([size(this.roi) this.LENK]);
            for t = 1:length(k)
                img = zeros(size(this.roi), 'single');
                img(roibin) = k(t);
                ks_.img(:,:,:,t) = img;
            end
            ks_.fileprefix = this.sessionData.vsOnAtlas('typ', 'fp', 'tags', [this.blurTag this.regionTag]);
            ks_ = imagingType(ipr.typ, ks_);
        end 
        function vs_ = vs(this, varargin)
            vs_ = this.ks(varargin{:});
        end 
    end
    
    %% PROTECTED
    
    methods (Access = protected)
		  
 		function this = Martin1987_20210210(varargin)
            this = this@mlpet.TracerKineticsStrategy(varargin{:});
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

