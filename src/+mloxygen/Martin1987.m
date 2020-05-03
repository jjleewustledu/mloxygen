classdef Martin1987 < handle & matlab.mixin.Copyable 
	%% Martin1987  

	%  $Revision$
 	%  was created 31-Oct-2018 15:20:13 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.4.0.813654 (R2018a) for MACI64.  Copyright 2018 John Joowon Lee.
 	
    properties (Constant)
        BLOOD_DENSITY = 1.06 % https://hypertextbook.com/facts/2004/MichaelShmukler.shtml; human whole blood 37 C
        BRAIN_DENSITY = 1.05 % Torack et al., 1976
        RATIO_SMALL_LARGE_HCT = 0.85 % Grubb, et al., 1978
        PLASMA_DENSITY = 1.03
    end
    
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
            aif = aif(this.T0:this.Tf);
            
            if ipr.averageVoxels
                s = this.scanner_.volumeAveraged(ipr.roi);
                times = s.times;
                times = times(this.T0 < s.times & s.times < this.Tf);
                tac = s.activityDensity();
                tac = tac(this.T0 < s.times & s.times < this.Tf);
                v = scale * trapz(times, tac)/ ...
                    trapz(this.T0:this.Tf, aif);
            else
                s = this.scanner_.masked(ipr.roi);
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
                    ifc.fileprefix = strrep(ifc.fileprefix, 'oc', 'cbv');
                else
                    ifc.fileprefix = strrep(ifc.fileprefix, 'oc', 'v1');
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
    end
    
    %% PROTECTED
    
    methods (Access = protected)		  
 		function this = Martin1987(varargin)
 			%% Martin1987
            %  @param devkit is mlpet.IDeviceKit.
 			%  @param T0 is numeric.
            %  @param Tf is numeric.
            %  @param ispercent is logical.

            ip = inputParser;
            addParameter(ip, 'devkit', [], @(x) isa(x, 'mlpet.IDeviceKit'));
            addParameter(ip, 'T0', nan, @isnumeric);
            addParameter(ip, 'Tf', nan, @isnumeric);
            addParameter(ip, 'ispercent', true, @islogical);
            addParameter(ip, 'averageVoxels', true, @islogical);
            parse(ip, varargin{:});
            ipr = ip.Results;
 			
            this.devkit_    = ipr.devkit;
            this.T0_        = ipr.T0;
            this.Tf_        = ipr.Tf;
            this.ispercent_ = ipr.ispercent;
            this.averageVoxels_ = ipr.averageVoxels;
            
            this.scanner_  = this.devkit_.buildScannerDevice();
            this.arterial_ = this.devkit_.buildArterialSamplingDevice(this.scanner_);
            this.counting_ = this.devkit_.buildCountingDevice();
 		end
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
        end
    end 
    
    %% PRIVATE
    
    properties (Access = private)
        averageVoxels_
        arterial_
        counting_
        devkit_
        ispercent_
        scanner_
        T0_
        Tf_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

