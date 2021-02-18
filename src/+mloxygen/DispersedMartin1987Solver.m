classdef DispersedMartin1987Solver 
	%% DISPERSEDMARTIN1987SOLVER  

	%  $Revision$
 	%  was created 10-Feb-2021 18:48:53 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.9.0.1570001 (R2020b) Update 4 for MACI64.  Copyright 2021 John Joowon Lee.
 	
	properties        
        artery_interpolated
        quiet = false
        visualize = false
        
        context
        fileprefix
        Measurement          % external data
        map                  % containers.Map containing model params as structs with fields:  min, max, init
        model                %
        sigma0 = 0.05        % fraction of Measurement < 1
        times_sampled        % numeric times for Measurement; midpoints of frames for PET 
        zoom = 1
        
        Dt
 	end
    
	properties (Dependent)   
        ks
        results
    end

	methods 
        
        %% GET
        
        function g = get.ks(this)
            g = this.results_.ks;
        end
        function g = get.results(this)
            g = this.results_;
        end
        
        %%
		  
 		function this = DispersedMartin1987Solver(varargin)
            this.fileprefix = ...
                sprintf('%s_ctor_%s', strrep(class(this), '.', '_'), datestr(now, 'yyyymmddHHMMSS'));
            
            ip = inputParser;            
            ip.KeepUnmatched = true;
            addParameter(ip, 'context', [], @(x) ~isempty(x))
            addParameter(ip, 'sigma0', this.sigma0, @(x) isscalar(x) && x <= 1)
            addParameter(ip, 'fileprefix', this.fileprefix, @ischar)
            addParameter(ip, 'Dt', [], @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;            
            
            this.context = ipr.context;
            this.model = this.context.model;                           % copy objects for speed    
            this.artery_interpolated = this.model.artery_interpolated; %
            this.map = this.model.map;                                 %
            this.Measurement = this.context.measurement;               %
            this.times_sampled = this.model.times_sampled;             %  
            this.sigma0 = ipr.sigma0;            
            this.fileprefix = ipr.fileprefix;
            this.Dt = ipr.Dt;
        end
        
        function fprintfModel(this)
            fprintf('DispersedMartin1987Solver:\n');            
            for ky = 1:length(this.ks)
                fprintf('\tk%i = %f\n', ky, this.ks(ky));
            end 
            fprintf('\tT0 => %g\n', this.model.T0); 
            fprintf('\tTf => %g\n', this.model.Tf); 
            fprintf('\tDt = %f\n', this.Dt);
        end
        function [k,sk] = k1(this, varargin)
            [k,sk] = find_result(this, 'k1');
        end
        function h = plot(this, varargin)
            ip = inputParser;
            addParameter(ip, 'showAif', true, @islogical)
            addParameter(ip, 'xlim', [-10 500], @isnumeric)            
            addParameter(ip, 'ylim', [], @isnumeric)
            addParameter(ip, 'zoom', 10, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.zoom = ipr.zoom;
            
            RR = mlraichle.RaichleRegistry.instance();
            tBuffer = RR.tBuffer;
            aif = this.artery_interpolated;
            h = figure;
            times = this.times_sampled;
            
            if ipr.showAif
                plot(times, ipr.zoom*this.Measurement, ':o', ...
                    0:length(aif)-tBuffer-1, aif(tBuffer+1:end), '--') 
                if ipr.zoom > 1
                    leg_aif = sprintf('aif x%i', ipr.zoom);
                else
                    leg_aif = 'aif';
                end
                legend('measurement', leg_aif)
            else
                plot(times, ipr.zoom*this.Measurement, 'o')                
                legend('measurement')
            end
            if ~isempty(ipr.xlim); xlim(ipr.xlim); end
            if ~isempty(ipr.ylim); ylim(ipr.ylim); end
            xlabel('times / s')
            ylabel('activity / (Bq/mL)')
            annotation('textbox', [.5 .3 .3 .3], 'String', sprintfModel(this), 'FitBoxToText', 'on', 'FontSize', 7, 'LineStyle', 'none')
            dbs = dbstack;
            title(dbs(1).name)
        end
        function save(this)
            save([this.fileprefix '.mat'], this);
        end
        function saveas(this, fn)
            save(fn, this);
        end   
        function this = solve(this, varargin)            
            [~,idx0] = max(this.Measurement > 0.1*max(this.Measurement));
            ts = this.times_sampled;
            T0 = this.model.T0 + ts(idx0);
            Tf = this.model.Tf + ts(idx0);
            RR = mlraichle.RaichleRegistry.instance();
            tBuffer = RR.tBuffer;
            
            tsWindow = T0 <= ts & ts <= Tf;
            tac = this.Measurement(tsWindow);
            tacTimeRange = ts(tsWindow);  
            rTacTimeRange = ceil(tacTimeRange);
            aif_0_f = this.artery_interpolated(tBuffer+rTacTimeRange(1)+1:tBuffer+rTacTimeRange(end)+1);
            
            scale = 1 / mlpet.TracerKinetics.RATIO_SMALL_LARGE_HCT;
            ks_ = scale * ...
                trapz(tacTimeRange, tac) / ...
                trapz(aif_0_f);
            
            this.results_ = struct('ks', ks_); 
            if ~this.quiet
                fprintfModel(this)
            end
            if this.visualize
                plot(this)
            end
        end 
        function this = solve_with_makima(this, varargin)            
            T0 = this.model.T0;
            Tf = this.model.Tf;
            
            dMeasurement = diff(this.Measurement);
            [~,idxdM] = max(dMeasurement > 0.1*max(dMeasurement));
            Measurement_ = makima( ...
                this.times_sampled(idxdM:end), ...
                this.Measurement(idxdM:end), ...
                this.times_sampled(idxdM):this.times_sampled(end));
            [~,idxai] = max(this.artery_interpolated > 0.1*max(this.artery_interpolated));
            artery_interpolated_ = this.artery_interpolated(idxai:end);
            
            scale = 1 / mlpet.TracerKinetics.RATIO_SMALL_LARGE_HCT;
            ks_ = scale * ...
                trapz(Measurement_(T0+1:Tf+1)) / ...
                trapz(artery_interpolated_(T0+1:Tf+1));
            
            this.results_ = struct('ks', ks_); 
            if ~this.quiet
                fprintfModel(this)
            end
            if this.visualize
                plot(this)
            end
        end 
        function s = sprintfModel(this)
            s = sprintf('DispersedMartin1987Solver:\n');
            for ky = 1:length(this.ks)
                s = [s sprintf('\tk%i = %f\n', ky, this.ks(ky))]; %#ok<AGROW>
            end
            s = [s sprintf('\tT0 => %g\n', this.model.T0)]; 
            s = [s sprintf('\tTf => %g\n', this.model.Tf)]; 
            s = [s sprintf('\tDt = %f\n', this.Dt)];
        end 
 	end 
    
    %% PROTECTED
    
    properties (Access = protected)        
        results_
    end
    
    methods (Access = protected)
        function [m,sd] = find_result(this, lbl)
            ks_ = this.ks;
            assert(strcmp(lbl(1), 'k'))
            ik = str2double(lbl(2));
            m = ks_(ik);
            sd = 0;
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

