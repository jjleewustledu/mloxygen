classdef Martin1987 < handle & mlpet.TracerKinetics
	%% Martin1987 

	%  $Revision$
 	%  was created 31-Oct-2018 15:20:13 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
    
	properties (Dependent)
        aif
        averageVoxels % logical
        roi % mlfourd.ImagingContext2
        T0 % integration limits
        Tf % integration limits
    end
    
    methods (Static)
        function [this,oo,aif] = createFromDeviceKit(devkit, varargin)
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'roi', 'brain_222.4dfp.hdr')
            addParameter(ip, 'averageVoxels', true, @islogical)
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            
            sesd = ipr.devkit.sessionData;
            sesd.jitOn111(sesd.ocOnAtlas())
            this = mloxygen.Martin1987( ...
                'devkit', ipr.devkit, 'averageVoxels', ipr.averageVoxels, 'roi', ipr.roi, 'T0', 120, 'TF', 240);
            oo = this.scanner_.activityDensity();
            aif = this.aif;
        end
    end

	methods 
        
        %% GET, SET
        
        function g = get.aif(this)
            g = this.arterial_.activityDensity();
        end
        function g = get.averageVoxels(this)
            g = this.averageVoxels_;
        end
        function     set.averageVoxels(this, s)
            assert(islogical(s))
            this.averageVoxels_ = s;
        end
        function g = get.roi(this)
            g = this.roi_;
        end    
        function     set.roi(this, s)
            if ~isempty(s)
                this.roi_ = mlfourd.ImagingContext2(s);
                this.roi_ = this.roi_.binarized();
            end
        end
        function g = get.T0(this)
            g = this.T0_;
        end
        function g = get.Tf(this)
            g = this.Tf_;
        end
        
        %%
        
        function cbv = buildCbv(this, varargin)
            %  @param roi is ImagingContext2.
            %  @param averageVoxels is logical.
            %  @return cbv := numeric if averageVoxels;
            %              := ImagingContext2 if not averageVoxels.  (DEPRECATED)
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'roi', this.roi, @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'averageVoxels', this.averageVoxels, @islogical)
            parse(ip, varargin{:})
            ipr = ip.Results;

            scale = 100/(this.RATIO_SMALL_LARGE_HCT * this.BRAIN_DENSITY);
            aif_ = this.aif;
            this.boundTfByAif(aif_);
            aif_ = aif_(this.T0+1:this.Tf+1);
            
            if ipr.averageVoxels
                s = this.scanner_.volumeAveraged(ipr.roi);
                times = s.times;
                times = times(this.T0 <= s.times & s.times <= this.Tf);
                tac = s.activityDensity();
                tac = tac(this.T0 <= s.times & s.times <= this.Tf);
                cbv = scale * trapz(times, tac)/ ...
                    trapz(this.T0:this.Tf, aif_);
            else
                
                %% DEPRECATED
                
                s = this.scanner_.masked(ipr.roi);
                times = s.times;
                times = times(this.T0 <= s.times & s.times <= this.Tf);
                tac = s.activityDensity();
                tac = tac(:,:,:,this.T0 <= s.times & s.times <= this.Tf);
                tac = permute(tac, [4 1 2 3]);
                img = scale * trapz(times, tac)/ ...
                    trapz(this.T0:this.Tf, aif_);
                img = squeeze(img);
                ifc = this.scanner_.imagingContext.fourdfp;
                ifc.img = img;
                ifc.fileprefix = this.sessionData.cbvOnAtlas('typ', 'fp');
                cbv = mlfourd.ImagingContext2(ifc);
            end
        end
        function buildQC(~, varargin)
        end
        function h = plot(this, varargin)
            ip = inputParser;
            addParameter(ip, 'index', [], @isnumeric)
            addParameter(ip, 'roi', this.roi)
            addParameter(ip, 'zoom', 10, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;            
            ipr.roi = mlfourd.ImagingContext2(ipr.roi);
            
            a = this.arterial_;
            s = this.scanner_;
            tac = s.imagingContext;
            if rank(tac) == 4
                tac = tac.volumeAveraged(ipr.roi);
            end
            if rank(tac) == 3
                error('mloxygen:ValueError', 'Martin1987.plot:  rank(tac) == 3 is not supported')
            end
            h = figure; 
            plot(a.datetimes, a.activityDensity(), ':o', ...
                 s.datetimes, ipr.zoom*tac.fourdfp.img, ':o')
            ylabel('activity / (Bq/mL)')
            legend('aif', [num2str(ipr.zoom) 'x oc'])            
            dbs = dbstack;
            title(sprintf('%s: indices->%s\n%s', dbs(1).name, mat2str(ipr.index), ipr.roi.fileprefix))
        end
        function this = solve(this, varargin)
            %  @param tags is char, e.g., 'wmparc1', 'wholebrain', 'voxel'
            
            ip = inputParser;
            addParameter(ip, 'tags', ['_' this.sessionData.region], @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            cbv = this.roi.zeros;
            cbv = cbv.fourdfp;
            cbv.fileprefix = this.sessionData.cbvOnAtlas('typ', 'fp', 'tags', ipr.tags);            
            cbv.img(logical(this.roi)) = this.buildCbv(varargin{:});
            this.product_ = mlfourd.ImagingContext2(cbv);            
        end
        
        %% products
        
        function obj = cbv(this, varargin)
            ip = inputParser;
            addParameter(ip, 'typ', 'mlfourd.ImagingContext2', @ischar)
            parse(ip, varargin{:})
            obj = imagingType(ip.Results.typ, this.product_);
        end
        function obj = v1(this, varargin)
            ip = inputParser;
            addParameter(ip, 'typ', 'mlfourd.ImagingContext2', @ischar)
            parse(ip, varargin{:})
            
            v1 = copy(this.product_.fourdfp);
            v1.img = this.cbvToV1(v1.img);
            v1.fileprefix = strrep(this.product_.fileprefix, 'cbv', 'v1');
            obj = imagingType(ip.Results.typ, v1);
        end
    end
    
    %% PROTECTED
    
    methods (Access = protected)		  
 		function this = Martin1987(varargin)
 			%% Martin1987
            %  @param devkit is mlpet.IDeviceKit.
            %  @param averageVoxels is logical.
            %  @param roi is understood by ImagingContext2.
 			%  @param T0 is numeric.
            %  @param Tf is numeric.
            
            this = this@mlpet.TracerKinetics(varargin{:});

            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'averageVoxels', true, @islogical);
            addParameter(ip, 'roi', [])
            addParameter(ip, 'T0', nan, @isnumeric);
            addParameter(ip, 'Tf', nan, @isnumeric);
            parse(ip, varargin{:});
            ipr = ip.Results;
 			
            this.averageVoxels_ = ipr.averageVoxels;
            if ~isempty(ipr.roi)
                this.roi_ = mlfourd.ImagingContext2(ipr.roi);
                this.roi_ = this.roi_.binarized();
            end
            this.T0_        = ipr.T0;
            this.Tf_        = ipr.Tf;
            
            this.scanner_  = this.devkit_.buildScannerDevice();
            this.arterial_ = this.devkit_.buildArterialSamplingDevice(this.scanner_);
            this.counting_ = this.devkit_.buildCountingDevice();
        end
        function boundTfByAif(this, aif)
            if (this.Tf+1)/length(aif) > 1.25
                error('mloxygen:ValueError', ...
                    'Martin1987.boundTfByAif:  Tf->%g exceeds length(aif)->%g', this.Tf, length(aif))
            end
            if (this.Tf+1) > length(aif)
                warning('mloxygen:ValueWarning', ...
                    'Martin1987.boundTfByAif:  Tf->%g exceeds length(aif)->%g', this.Tf, length(aif))
                this.Tf_ = length(aif)-1;
            end
        end
    end 
    
    %% PRIVATE
    
    properties (Access = private)
        averageVoxels_
        arterial_ % := devkit
        counting_ % := devkit
        product_  % := cbv as mlfourd.ImagingContext2
        roi_
        scanner_ % := devkit
        T0_
        Tf_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

