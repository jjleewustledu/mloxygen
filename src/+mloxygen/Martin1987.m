classdef Martin1987 < handle & mlpet.TracerKinetics
	%% Martin1987  

	%  $Revision$
 	%  was created 31-Oct-2018 15:20:13 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
    
	properties (Dependent)
        averageVoxels
        ispercent
        T0 % integration limits
        Tf % integration limits
    end
    
    methods (Static)
        function this = createFromDeviceKit(devkit)
            this = mloxygen.Martin1987('devkit', devkit, 'T0', 120, 'TF', 240);
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
        function g = get.ispercent(this)
            g = this.ispercent_;
        end
        function g = get.T0(this)
            g = this.T0_;
        end
        function g = get.Tf(this)
            g = this.Tf_;
        end
        
        %%
        
        function v = buildCbv(this, varargin)
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'averageVoxels', this.averageVoxels, @islogical)
            parse(ip, varargin{:})
            ipr = ip.Results;

            scale = 1/(this.RATIO_SMALL_LARGE_HCT * this.BRAIN_DENSITY * this.BLOOD_DENSITY );
            if this.ispercent
                scale = 100*scale;
            end
            aif = this.arterial_.activityDensity();
            this.boundTfByAif(aif);
            aif = aif(this.T0+1:this.Tf+1);
            
            if ipr.averageVoxels
                s = this.scanner_.volumeAveraged(ipr.roi);
                times = s.times;
                times = times(this.T0 <= s.times & s.times <= this.Tf);
                tac = s.activityDensity();
                tac = tac(this.T0 <= s.times & s.times <= this.Tf);
                v = scale * trapz(times, tac)/ ...
                    trapz(this.T0:this.Tf, aif);
            else
                s = this.scanner_.masked(ipr.roi);
                % s = s.blurred(4.3); % for displaying AIF
                times = s.times;
                times = times(this.T0 <= s.times & s.times <= this.Tf);
                tac = s.activityDensity();
                tac = tac(:,:,:,this.T0 <= s.times & s.times <= this.Tf);
                tac = permute(tac, [4 1 2 3]);
                img = scale * trapz(times, tac)/ ...
                    trapz(this.T0:this.Tf, aif);
                img = squeeze(img);
                ifc = this.scanner_.imagingContext.fourdfp;
                ifc.img = img;
                if this.ispercent
                    ifc.fileprefix = this.sessionData.cbvOnAtlas('typ', 'fp');
                else
                    ifc.fileprefix = this.sessionData.v1OnAtlas('typ', 'fp');
                end
                v = mlfourd.ImagingContext2(ifc);
            end
        end
        function buildQC(this, varargin)

            return
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'cbv', [], @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.scanner_.masked(ipr.roi)            
            disp(this.arterial_)
            this.arterial_.plot()
            title('Martin1987.buildQC().this.arterial_.plot()')
            disp(this.scanner_)
            this.scanner_.plot()
            title('Martin1987.buildQC().this.scanner_.plot()')
            this.scanner_.imagingContext.fsleyes()
        end
        function h = plot(this, varargin)
            ip = inputParser;
            addParameter(ip, 'index', [], @isnumeric)
            addParameter(ip, 'roi', [])
            parse(ip, varargin{:})
            ipr = ip.Results;
            ipr.roi = mlfourd.ImagingContext2(ipr.roi);
            
            a = this.arterial_;
            s = this.scanner_;
            
            h = figure; 
            plot(a.datetimes, a.activityDensity(), ':o', ...
                s.datetimes, 10*s.activityDensity(), ':o')
            ylabel('activity / (Bq/mL)')
            legend('aif', '10x oc')
            title(sprintf('mloxygen.Margin1987.plot():  index->%g', ipr.index))
        end
    end
    
    %% PROTECTED
    
    methods (Access = protected)		  
 		function this = Martin1987(varargin)
 			%% Martin1987
            %  @param devkit is mlpet.IDeviceKit.
 			%  @param T0 is numeric.
            %  @param Tf is numeric.
            %  @param ispercent is logical.
            %  @param averageVoxels is logical.
            
            this = this@mlpet.TracerKinetics(varargin{:});

            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'T0', nan, @isnumeric);
            addParameter(ip, 'Tf', nan, @isnumeric);
            addParameter(ip, 'ispercent', true, @islogical);
            addParameter(ip, 'averageVoxels', true, @islogical);
            parse(ip, varargin{:});
            ipr = ip.Results;
 			
            this.T0_        = ipr.T0;
            this.Tf_        = ipr.Tf;
            this.ispercent_ = ipr.ispercent;
            this.averageVoxels_ = ipr.averageVoxels;
            
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
        arterial_
        counting_
        ispercent_
        scanner_
        T0_
        Tf_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

