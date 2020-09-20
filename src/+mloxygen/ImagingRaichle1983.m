classdef ImagingRaichle1983 < handle & matlab.mixin.Copyable
	%% IMAGINGRAICHLE1983  
    %  builds Raichle models for imaging voxels.
    %  For performance considerations, see also https://blogs.mathworks.com/loren/2012/03/26/considering-performance-in-object-oriented-matlab-code/

	%  $Revision$
 	%  was created 11-Sep-2020 00:49:25 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.
 	
	properties (Constant)
        HALFLIFE = 122.2416; % s
        JITTER = 0.0 % > 0 to aid deep learning
        LENF = 3
        MAX_NORMAL_BACKGROUND = 20 % Bq/mL 		
    end
    
	properties
        artery_plasma_interpolated
        devkit
        ho % scanner.activityDensity(), decaying & calibrated
        fs % flow 1/s
        meanAbsError
        measurement % expose for performance when used by mloxygen.Raichle1983Strategy
        model       %
        normMeanAbsError
        Nroi
        prediction
        residual
        roi
        roibin % is logical
        taus
 		times_sampled
    end
    
    properties (Dependent)        
        regionTag
        sessionData
    end
    
    methods (Static)
        function this = createFromDeviceKit(devkit, varargin)
            %% makes no adjustments of AIF timing
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param ho is numeric, default from devkit.
            %  @param roi is understood by mlfourd.ImagingContext2.
            %  @param blurHo := {[], 0, 4.3, ...}
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'ho', [], @isnumeric)
            addParameter(ip, 'roi', [], @(x) ~isempty(x))
            addParameter(ip, 'blurHo', 4.3, @isnumeric)
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            
            % scanner provides calibrations, ancillary data
            
            scanner = ipr.devkit.buildScannerDevice();
            scanner = scanner.blurred(ipr.blurHo);
            ho = scanner.activityDensity(); % calibrated, decaying
            
            % AIF  
            
            arterial = ipr.devkit.buildArterialSamplingDevice(scanner);
            aif = pchip(arterial.times, arterial.activityDensity(), 0:scanner.times(end));
            
            this = mloxygen.ImagingRaichle1983( ...
                'devkit', devkit, ...
                'ho', ho, ...
                'taus', scanner.taus, ...
                'times_sampled', scanner.timesMid, ...
                'artery_sampled', aif, ...
                varargin{:});
        end
    end

	methods 
        
        %% GET
        
        function g = get.regionTag(this)
            g = this.sessionData.regionTag;
        end
        function g = get.sessionData(this)
            g = this.devkit.sessionData;
        end
        
        %%		
        
        function [ic,nic] = buildMeanAbsError(this)
            % @return \Sigma_i activity_i \frac{tau_i}/{T}
            
            if isempty(this.residual)
                this.buildResidual()
            end
            assert(~isempty(this.taus))
            sesd = this.devkit.sessionData;
            hoDCorr = sesd.hoOnAtlas('typ', 'mlfourd.ImagingContext2');
            
            % MAE
            ic = abs(copy(this.residual));
            ic = ic.timeAveraged('taus', this.taus);
            ic.fileprefix = [hoDCorr.fileprefix this.regionTag '_MAE'];
            this.meanAbsError = ic;
            
            % NMAE
            hoAvgt = copy(hoDCorr);
            hoAvgt = hoAvgt.timeAveraged('taus', this.taus);         
            nic = copy(this.meanAbsError) ./ hoAvgt;
            nic.fileprefix = [hoDCorr.fileprefix this.regionTag '_NMAE'];
            this.normMeanAbsError = nic;
        end
        function ic = buildPrediction(this, varargin)
            %% @param reuseExisting is logical; default is false.
            
            ip = inputParser;
            addParameter(ip, 'reuseExisting', false, @islogical)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            assert(~isempty(this.fs))
            sesd = this.devkit.sessionData;
            scanner = this.devkit.buildScannerDevice();
            
            % init working_ifc from hoOnAtlas
            % return existing if reuseExisting
            working_ifc = sesd.hoOnAtlas('typ', 'mlfourd.ImagingFormatContext');
            fileprefix0 = [working_ifc.fileprefix this.regionTag '_predicted'];
            working_ifc.fileprefix = fileprefix0;
            sz = size(working_ifc);
            if ipr.reuseExisting && isfile(working_ifc.fqfilename)
                ic = mlfourd.ImagingContext2(working_ifc.fqfilename);
                this.prediction = ic;
                return
            end            
            
            % represent prediction in R^2:  voxels^1 (x) times            
            fs_ = zeros(this.Nroi,this.LENF+1);
            for i = 1:this.LENF+1     
                rate = this.fs.fourdfp.img(:,:,:,i);                
                fs_(:, i) = rate(this.roibin);
            end         
            this.ensureModel() % without voxelwise adjustments of AIF timings  
            img2d = zeros(this.Nroi, sz(4));          
            for vxl = 1:this.Nroi                
                img2d(vxl,:) = this.model.simulated(fs_(vxl,:)); % adjust AIF timings
            end
            
            % embed prediction in R^4:  voxels^3 (x) times
            working_ifc.img = zeros(sz);
            for t = 1:sz(4)
                frame = zeros(sz(1:3));
                frame(this.roibin) = img2d(:,t);
                working_ifc.img(:,:,:,t) = frame;
            end
            
            % decay-correct & remove calibrations
            ic = scanner.decayCorrectLike(working_ifc);
            ic = ic ./ scanner.invEfficiencyf(sesd);
            ic.fileprefix = fileprefix0;
            this.prediction = ic;
        end
        function ic = buildResidual(this)
            if isempty(this.prediction)
                this.buildPrediction()
            end            
            sesd = this.devkit.sessionData;
            hoDCorr = sesd.hoOnAtlas('typ', 'mlfourd.ImagingContext2');
            ic = this.prediction - hoDCorr;
            ic.fileprefix = [hoDCorr.fileprefix this.regionTag  '_residual'];
            this.residual = ic;
        end
        function this = ensureModel(this)
            if isempty(this.model)                               
                map = mloxygen.Raichle1983Model.preferredMap();
                this.model = mloxygen.Raichle1983Model( ...
                    'map', map, ...
                    'times_sampled', this.times_sampled, ...
                    'artery_interpolated', this.artery_plasma_interpolated);
            end            
        end
        function this = solve(this)
            ho_img_2d = this.projectedHoArray();
            fs_img_2d = zeros(dipsum(this.roibin), 4);
            map = mloxygen.Raichle1983Model.preferredMap();
            
            for vxl = 1:this.Nroi
                try
                    tic
                    fprintf('mloxygen.ImagingRaichle1983.solve():  vxl->%i this.Nroi->%i\n', vxl, this.Nroi)
                    this.measurement = ho_img_2d(vxl, :);
                    if this.sufficientData(this.measurement)
                        this.model = mloxygen.Raichle1983Model( ...
                            'map', map, ...
                            'times_sampled', this.times_sampled, ...
                            'artery_interpolated', this.artery_plasma_interpolated);
                        strategy = mloxygen.Raichle1983SimulAnneal('context', this);
                        strategy = solve(strategy);
                        
                        % store latest solutions
                        fs_img_2d(vxl, 1) = k1(strategy);
                        fs_img_2d(vxl, 2) = k2(strategy);
                        fs_img_2d(vxl, 3) = k3(strategy);
                        
                        % use latest solutions for initial conditions for solving neighboring voxels
                        map_k1 = map('k1'); % cache mapped struct
                        map_k2 = map('k2');
                        map_k3 = map('k3');                      
                        map_k1.init = fs_img_2d(vxl, 1); % cached struct.init := latest solutions with jitter
                        map_k2.init = fs_img_2d(vxl, 2);
                        map_k3.init = fs_img_2d(vxl, 3);
                        map('k1') = map_k1; % update mapped struct with adjusted cache
                        map('k2') = map_k2;
                        map('k3') = map_k3;
                    
                        % else fs_img_2d retains zeros                        
                    end                    
                    toc
                catch ME
                    %handexcept(ME)
                    handwarning(ME)
                end
            end
            
            this.fs = this.invProjectedFsArray(fs_img_2d);
        end
        function tf = sufficientData(this, measurement)
            tf = mean(measurement) > this.MAX_NORMAL_BACKGROUND;
        end
  	end 

    %% PROTECTED
    
	methods (Access = protected)
        function this = ImagingRaichle1983(varargin)
            %% IMAGINGRaichle1983
            %  @param devkit is mlpet.IDeviceKit.
            %  @param ho is understood by mlfourd.ImagingContext2.
            %  @param times_sampled from scanner is numeric.
            %  @param artery_sampled from counter is numeric.
            %  @param roi is understood by mlfourd.ImagingContext2.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'ho', [])
            addParameter(ip, 'taus', [], @isnumeric)
            addParameter(ip, 'times_sampled', [], @isnumeric)
            addParameter(ip, 'artery_sampled', [], @isnumeric)
            addParameter(ip, 'roi', [], @(x) ~isempty(x))
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.devkit = ipr.devkit;
            this.ho = mlfourd.ImagingContext2(ipr.ho);
            assert(4 == ndims(this.ho))
            this.taus = ipr.taus;
            this.times_sampled = ipr.times_sampled;  
            this.hct = ipr.hct;
            
            % artery management            
            t = 0:this.times_sampled(end);
            this.artery_plasma_interpolated = pchip(0:length(ipr.artery_sampled)-1, ipr.artery_sampled, t);
            
            this.roi = mlfourd.ImagingContext2(ipr.roi);
            this.roibin = this.roi.fourdfp.img > 0;
            this.ho = this.masked(this.ho);
            this.Nroi = dipsum(this.roibin);
        end
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
        end
        function ic = invProjectedFsArray(this, arr)
            img = zeros([size(this.roibin) this.LENF]);
            for i = 1:this.LENF
                cube = img(:,:,:,i);
                cube(this.roibin) = arr(:,i); 
                img(:,:,:,i) = cube;
            end
            fp = strrep(this.ho.fileprefix, 'ho', 'fs'); % sprintf('fs_%s_', this.roi.fileprefix));
            ic = mlfourd.ImagingContext2(img, 'filename', [fp '_' this.roi.fileprefix '.4dfp.hdr'], 'mmppix', [2 2 2]);
        end
        function ic = masked(this, ic)
            %% retains original fileprefix.
            
            fp = ic.fileprefix;
            ic = ic.masked(this.roi);
            ic.fileprefix = fp;
        end
        function arr = projectedHoArray(this)
            Nt = size(this.ho, 4);            
            arr = zeros(dipsum(this.roibin), Nt);
            for it = 1:Nt
                cube = this.ho.fourdfp.img(:,:,:,it);
                arr(:,it) = cube(this.roibin);
            end
        end
  	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

