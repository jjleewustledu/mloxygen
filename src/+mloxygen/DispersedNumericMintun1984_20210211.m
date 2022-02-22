classdef DispersedNumericMintun1984_20210211 < handle & mloxygen.Mintun1984
	%% DISPERSEDNUMERICMINTUN1984_20210211  

	%  $Revision$
 	%  was created 05-Dec-2020 20:45:51 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.9.0.1524771 (R2020b) Update 2 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    methods (Static)
        function [this,oo,aif] = createFromDeviceKit(devkit, varargin)
            %% adjusts AIF timings for coincidence of inflow with tissue activity from scanner
            %  @param required devkit is mlpet.IDeviceKit.
            %  @param oo is numeric, default from devkit.
            %  @param model is an mloxygen.DispersedMintun1984Model.
            %  @param roi ...
            %  @param solver is in {'simulanneal'}, default := 'simulanneal'.            
            %  @param blurOo := {[], 0, 4.3, ...}            
            %  @param sigma0, default from mloptimization.SimulatedAnnealing.
            %  @return this.
            %  @return oo, blurred by ipr.blurOo.
            %  @return aif.
            
            reshapeArterial = @mloxygen.DispersedNumericMintun1984_20210211.reshapeArterial;
            reshapeScanner = @mloxygen.DispersedNumericMintun1984_20210211.reshapeScanner;
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'oo', [], @isnumeric)
            addParameter(ip, 'model', [], @(x) isa(x, 'mloxygen.DispersedMintun1984Model'))
            addParameter(ip, 'roi', 'brain_222.4dfp.hdr')
            addParameter(ip, 'blurOo', 4.3, @isnumeric)
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            
            % prepare atlas data
            
            sesd = ipr.devkit.sessionData;
            sesd.jitOnAtlas(sesd.wmparc1OnAtlas())
            sesd.jitOnAtlas(sesd.ooOnAtlas())
            
            % scanner provides calibrations, decay, ancillary data
            
            ipr.roi = mlfourd.ImagingContext2(ipr.roi);
            roibin = ipr.roi.binarized();            
            scanner = ipr.devkit.buildScannerDevice();
            scanner = scanner.blurred(ipr.blurOo);
            scanner = scanner.volumeAveraged(roibin);
            [ooTimesMid,oo] = reshapeScanner(scanner);
            
            % arterial reshaping makes it quantitatively comparable to scanner
            
            arterial = ipr.devkit.buildArterialSamplingDevice(scanner);
            [~,aif,Dt] = reshapeArterial(arterial, scanner, ipr.model);
            
            % assemble this
            
            ipr.model = set_times_sampled(ipr.model, ooTimesMid);
            ipr.model = set_artery_interpolated(ipr.model, aif);
            this = mloxygen.DispersedNumericMintun1984_20210211( ...
                'devkit', devkit, ...
                'oo', oo, ...
                'solver', 'simulanneal', ...
                'Dt', Dt, ...
                'model', ipr.model, ...
                'times_sampled', ooTimesMid, ...
                'artery_interpolated', aif, ...
                'roi', ipr.roi, ...
                varargin{:});
        end    
        function aif = extrapolateEarlyLateAif(aif__)
            aif = mloxygen.DispersedNumericRaichle1983.extrapolateEarlyLateAif(aif__);
        end
        function [aifTimes,aif,Dt] = reshapeArterial(arterial, scanner, model)
            %% 1.  form uniform time coordinates consistent with scanner
            %  2.  sample aif on uniform time coordinates
            %  3.  infer & apply shift of worldline, Dt, for aif
            
            extrapolateEarlyLateAif = @mloxygen.DispersedNumericMintun1984_20210211.extrapolateEarlyLateAif;
            reshapeScanner = @mloxygen.DispersedNumericMintun1984_20210211.reshapeScanner;
            shiftWorldlines = @mloxygen.DispersedNumericMintun1984_20210211.shiftWorldlines;
                        
            [ooTimesMid,~,ooTimeMin] = reshapeScanner(scanner); % ooTimesMid(1) ~ 0; ooTimeMin ~ -5
            aifTimes = ooTimesMid(1):ooTimesMid(end); % aifTimes ~ [0 ... 574]
            ooDatetime0 = scanner.datetime0 + seconds(ooTimeMin); % ooDatetime0 ~ 23-May-2019 12:04:20
            aifDatetime0 = arterial.datetime0; % aifDatetime0 ~ 23-May-2019 12:03:59
            Ddatetime = seconds(aifDatetime0 - ooDatetime0); % Ddatetime ~ -21
            aifTimes__ = (arterial.time0:arterial.timeF) - arterial.time0 + Ddatetime + ooTimesMid(1); 
            % aifTimes__ ~ [-26 -25 -24 ... 181]; satisfies 1.
            
            aif__ = arterial.activityDensity();
            aif__ = extrapolateEarlyLateAif(aif__);
            aif = makima([aifTimes__ aifTimes(end)], [aif__ 0], aifTimes); % satisfies 2.             
            
            Dt = model.fs_Raichle_Martin(5);
            aif = shiftWorldlines(aif, Dt);
        end
        function [timesMid,ho,timeMin] = reshapeScanner(scanner)
            [timesMid,ho,timeMin] = mloxygen.DispersedNumericRaichle1983.reshapeScanner(scanner);
        end
        function aif = shiftWorldlines(aif__, Dt)
            aif = mloxygen.DispersedNumericRaichle1983.shiftWorldlines(aif__, Dt);
        end
    end
    
    methods
        function ho  = checkSimulated(this, varargin)
            %% CHECKSIMULATED simulates tissue activity with passed and internal parameters without changing state.
            %  @param required ks is variational [k1 k2].
            %  @param aif is numeric; default is this.artery_interpolated for model state.
            %  @return ho simulation is numeric.
        
            ip = inputParser;
            addOptional(ip, 'ks', this.ks(), @isnumeric)
            addParameter(ip, 'aif', this.artery_interpolated, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;  
            
            ho = this.model.simulated(ipr.ks, this.model.fs_Raichle_Martin, 'aif', ipr.aif, 'Dt', this.Dt);
        end
    end
    
    %% PROTECTED
    
    methods (Access = protected) 
 		function this = DispersedNumericMintun1984_20210211(varargin)
 			%% DISPERSEDNUMERICMINTUN1984_20210211
            %  @param devkit is mlpet.IDeviceKit.
            %  @param oo is numeric.
            %  @param solver is in {'simulanneal'}.
            %  @param blurOo := {[], 0, 4.3, ...}            
            %  @param sigma0, default from mloptimization.SimulatedAnnealing.

 			this = this@mloxygen.Mintun1984(varargin{:});
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'oo', [], @(x) isnumeric(x))
            addParameter(ip, 'solver', 'simulanneal', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
                        
            this.measurement = ipr.oo;
            switch lower(ipr.solver)
                case 'simulanneal'
                    this.strategy_ = mloxygen.DispersedMintun1984SimulAnneal( ...
                        'context', this, varargin{:});
                otherwise
                    error('mloxygen:NotImplementedError', 'DispersedNumericMintun1984_20210211.ipr.solver->%s', ipr.solver)
            end
 		end
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

