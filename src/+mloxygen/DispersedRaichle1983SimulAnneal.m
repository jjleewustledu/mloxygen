classdef DispersedRaichle1983SimulAnneal < mloxygen.Raichle1983SimulAnneal
	%% DISPERSEDRAICHLE1983SIMULANNEAL  

	%  $Revision$
 	%  was created 10-Sep-2020 22:24:15 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.

    properties
        aifdata
        Dt
        DtMixing
        fracMixing
        v1
    end
    
	methods		  
 		function this = DispersedRaichle1983SimulAnneal(varargin)
 			this = this@mloxygen.Raichle1983SimulAnneal(varargin{:});
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            ip.PartialMatching = false;
            addParameter(ip, 'Dt', [], @isnumeric)
            addParameter(ip, 'DtMixing', [], @isnumeric)
            addParameter(ip, 'fracMixing', [], @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.Dt = ipr.Dt;
            this.DtMixing = ipr.DtMixing;
            this.fracMixing = ipr.fracMixing;
            this.v1 = this.model.v1;
            this.aifdata = mlaif.AifData.instance();
        end    
        
        function fprintfModel(this)
            fprintf('Simulated Annealing:\n');            
            for ky = 1:length(this.ks)
                fprintf('\tk%i = %g\n', ky, this.ks(ky));
            end 
            fprintf('\tE = 1 - exp(-PS/f) = %g\n', 1 - exp(-this.ks(2)/this.ks(1)))          
            fprintf('\tDt = %g\n', this.Dt);
            fprintf('\tv1 = %f\n', this.v1);
            fprintf('\ttBuffer = %g\n', this.aifdata.tBuffer)
            fprintf('\tDtMixing = %g\n', this.DtMixing)
            fprintf('\tfracMixing = %g\n', this.fracMixing)
            %fprintf('\ttimeCliff = %g\n', this.timeCliff)
            fprintf('\tloss = %g\n', this.loss())
            %fprintf('\tsigma0 = %g\n', this.sigma0);
            for ky = this.map.keys
                fprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})));
            end
        end
        function h = plot(this, varargin)
            ip = inputParser;
            addParameter(ip, 'showAif', true, @islogical)
            addParameter(ip, 'xlim', [-10 500], @isnumeric)
            addParameter(ip, 'ylim', [], @isnumeric)
            addParameter(ip, 'zoom', 3, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.zoom = ipr.zoom;
            
            tBuffer = this.aifdata.tBuffer;
            aif = this.dispersedAif(this.artery_interpolated, this.ks(end));
            h = figure;
            times = this.times_sampled;
            sampled = this.model.sampled(this.ks, this.v1, aif, times);
            
            if ipr.zoom > 1
                leg_meas = sprintf('measurement x%i', ipr.zoom);
            else
                leg_meas = 'measurement';
            end
            if ipr.zoom > 1
                leg_est = sprintf('estimation x%i', ipr.zoom);
            else
                leg_est = 'estimation';
            end
            if ipr.showAif
                plot(times, ipr.zoom*this.Measurement, ':o', ...
                    times(1:length(sampled)), ipr.zoom*sampled, '-', ...
                    -tBuffer:length(aif)-tBuffer-1, aif, '--')
                legend(leg_meas, leg_est, 'aif')
            else
                plot(times, ipr.zoom*this.Measurement, 'o', ...
                    times(1:length(sampled)), ipr.zoom*sampled, '-')                
                legend(leg_meas, leg_est)
            end
            if ~isempty(ipr.xlim); xlim(ipr.xlim); end
            if ~isempty(ipr.ylim); ylim(ipr.ylim); end
            xlabel('times / s')
            ylabel('activity / (Bq/mL)')
            annotation('textbox', [.25 .5 .3 .3], 'String', sprintfModel(this), 'FitBoxToText', 'on', 'FontSize', 8, 'LineStyle', 'none')
            dbs = dbstack;
            title(dbs(1).name)
        end 
        function this = solve(this, varargin)
            ip = inputParser;
            addRequired(ip, 'loss_function', @(x) isa(x, 'function_handle'))
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            options_fmincon = optimoptions('fmincon', ...
                'FunctionTolerance', 1e-12, ...
                'OptimalityTolerance', 1e-12, ...
                'TolCon', 1e-14, ...
                'TolX', 1e-14);
            if this.visualize_anneal
                options = optimoptions('simulannealbnd', ...
                    'AnnealingFcn', 'annealingboltz', ...
                    'FunctionTolerance', eps, ...
                    'HybridFcn', {@fmincon, options_fmincon}, ...
                    'InitialTemperature', 20, ...
                    'MaxFunEvals', 50000, ...
                    'ReannealInterval', 200, ...
                    'TemperatureFcn', 'temperatureexp', ...
                    'Display', 'diagnose', ...
                    'PlotFcns', {@saplotbestx,@saplotbestf,@saplotx,@saplotf,@saplotstopping,@saplottemperature});
            else
                options = optimoptions('simulannealbnd', ...
                    'AnnealingFcn', 'annealingboltz', ...
                    'FunctionTolerance', eps, ...
                    'HybridFcn', {@fmincon, options_fmincon}, ...
                    'InitialTemperature', 20, ...
                    'MaxFunEvals', 50000, ...
                    'ReannealInterval', 200, ...
                    'TemperatureFcn', 'temperatureexp');
            end
 			[ks_,sse,exitflag,output] = simulannealbnd( ...
                @(ks__) ipr.loss_function( ...
                       ks__, double(this.v1), this.artery_interpolated, this.times_sampled, double(this.Measurement), this.timeCliff), ...
                this.ks0, this.ks_lower, this.ks_upper, options); 
            
            this.product_ = struct('ks0', this.ks0, 'ks', ks_, 'sse', sse, 'exitflag', exitflag, 'output', output); 
            if ~this.quiet
                fprintfModel(this)
            end
            if this.visualize
                plot(this)
            end
        end 
        function s = sprintfModel(this)
            s = sprintf('Simulated Annealing:\n');
            for ky = 1:length(this.ks)
                s = [s sprintf('\tk%i = %g\n', ky, this.ks(ky))]; %#ok<AGROW>
            end
            s = [s sprintf('\tE = 1 - exp(-PS/f) = %g\n', 1 - exp(-this.ks(2)/this.ks(1)))];
            s = [s sprintf('\tDt = %g\n', this.Dt)];
            s = [s sprintf('\tv1 = %f\n', this.v1)];
            s = [s sprintf('\ttBuffer = %g\n', this.aifdata.tBuffer)];            
            s = [s sprintf('\tDtMixing = %g\n', this.DtMixing)];
            s = [s sprintf('\tfracMixing = %g\n', this.fracMixing)];
            %s = [s sprintf('\ttimeCliff = %g\n', this.timeCliff)];
            s = [s sprintf('\tloss = %g\n', this.loss())];
            %s = [s sprintf('\tsigma0 = %g\n', this.sigma0)];
            for ky = this.map.keys
                s = [s sprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})))]; %#ok<AGROW>
            end
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

