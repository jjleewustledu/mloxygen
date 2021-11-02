classdef (Abstract) FieldNumeric < handle & mlpet.TracerKinetics
	%% FIELDNUMERIC  

	%  $Revision$
 	%  was created 23-Aug-2021 20:13:57 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.10.0.1710957 (R2021a) Update 4 for MACI64.  Copyright 2021 John Joowon Lee.
 	
    methods (Abstract)
        this = solve(this)
    end
    
	properties (Dependent)
 		alpha_decay
        blurTag
        regionTag
        roi % mlfourd.ImagingContext2
        tF
    end
    
    properties
        artery_interpolated
        measurement
        t0
        taus
        timesMid
        tObs
    end
    
    methods (Static)
        function [tacs__,timesMid__,t0__] = prepareTacs(devkit, varargin)
            %  @return taus__ are the durations of emission frames
            %  @return t0__ is the start of the first emission frame.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addRequired(ip, 'devkit', @(x) isa(x, 'mlpet.IDeviceKit'))
            addParameter(ip, 'scanner', [], @(x) isa(x, 'mlpet.AbstractDevice'))
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            parse(ip, devkit, varargin{:})
            ipr = ip.Results;
            ipr.roi = ipr.roi.binarized();
            
            % scannerDevs provide calibrations & ROI-masking    
            blur = ipr.devkit.sessionData.petPointSpread;
            s = ipr.scanner.blurred(blur);
            s = s.masked(ipr.roi);
            tacs = s.activityDensity();
            tacs(tacs < 0) = 0;                       
            tacs__ = tacs;
            taus__ = s.taus;
            timesMid__ = s.timesMid;
            Nt = ceil(timesMid__(end));
            
            % estimate t0__
            tac_avgxyz = squeeze(mean(mean(mean(tacs__, 1), 2), 3));
            dtac_avgxyz = diff(tac_avgxyz);
            [~,idx] = max(dtac_avgxyz > 0.05*max(dtac_avgxyz));
            t0__ = timesMid__(idx) - taus__(idx)/2;
        end
    end
    
    methods
        
        %% GET
        
        function g = get.alpha_decay(~)
            g = mlpet.Radionuclides.decayConstantOf('15O');
        end
        function g = get.blurTag(~)
            g = mlraichle.StudyRegistry.instance.blurTag;
        end
        function g = get.regionTag(this)
            g = this.sessionData.regionTag;
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
        function g = get.tF(this)
            g = this.t0 + this.tObs;
        end        
    end
    
    %% PROTECTED
    
    properties (Access = protected)
        roi_
    end

	methods (Access = protected)        
        function savefig(this, varargin)
            ip = inputParser;
            addRequired(ip, 'handle', @ishandle) % fig handle
            addParameter(ip, 'tags', '', @ischar) % for filenames
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            tags = ipr.tags;
            if ~isempty(tags)
                tags_ = ['_' strrep(tags, ' ', '_')];
            else
                tags_ = '';
            end            
            dbs = dbstack;
            client = dbs(2).name;
            client_ = strrep(dbs(2).name, '.', '_');
            dtStr = datestr(this.sessionData.datetime);
            title(sprintf('%s\n%s %s', client, tags, dtStr))
            try
                dtTag = lower(this.sessionData.doseAdminDatetimeTag);
                savefig(ipr.handle, ...
                    fullfile(this.dataPath, ...
                    sprintf('%s%s_%s.fig', client_, tags_, dtTag)))
                figs = get(0, 'children');
                saveas(figs(1), ...
                    fullfile(this.dataPath, ...
                    sprintf('%si%s_%s.png', client_, tags_, dtTag)))
                close(figs(1))
            catch ME
                handwarning(ME)
            end
        end
        
 		function this = FieldNumeric(varargin)
 			%% FIELDNUMERIC
            %  @param devkit is mlpet.IDeviceKit.
            %  @param sessionData is mlpipeline.ISessionData, but defers to devkit.sessionData.
            %  @param roi is understood by mlfourd.ImagingContext2; will be binarized.

 			this = this@mlpet.TracerKinetics(varargin{:});
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            ip.PartialMatching = false;
            addParameter(ip, 'roi', [], @(x) isa(x, 'mlfourd.ImagingContext2'))
            addParameter(ip, 'timesMid', [], @isvector)
            addParameter(ip, 't0', 0, @isscalar)
            addParameter(ip, 'tObs', 40, @isscalar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            if ~isempty(ipr.roi)
                this.roi_ = mlfourd.ImagingContext2(ipr.roi);
                this.roi_ = this.roi_.binarized();
            end
            this.timesMid = ipr.timesMid;
            this.t0 = ipr.t0;
            this.tObs = ipr.tObs;
 		end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

