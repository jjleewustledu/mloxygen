classdef Mintun1984SimulAnneal < mlpet.TracerSimulAnneal & mloxygen.Mintun1984Strategy
	%% MINTUN1984SIMULANNEAL operates on single voxels/regions. 

	%  $Revision$
 	%  was created 07-Dec-2020 23:22:23 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.9.0.1524771 (R2020b) Update 2 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    properties
        timeCliff
    end
    
	methods		  
 		function this = Mintun1984SimulAnneal(varargin)
 			%% MINTUN1984SIMULANNEAL
            %  @param context is mloxygen.Mintun1984.
            %  @param sigma0.
            %  @param fileprefix.
 			
 			this = this@mlpet.TracerSimulAnneal(varargin{:});
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'timeCliff', Inf, @isscalar)
            parse(ip, varargin{:})
            this.timeCliff = ip.Results.timeCliff;
            
            [this.ks_lower,this.ks_upper,this.ks0] = remapper(this);
            this.artery_interpolated = this.model.artery_interpolated;
        end                
        
        function fprintfModel(this)
            fprintf('Simulated Annealing:\n');            
            for ky = 1:length(this.ks)
                fprintf('\tk%i = %f\n', ky, this.ks(ky));
            end            
            
            fs = this.model.fs_Raichle_Martin;
            for ky = 1:length(fs)-1
                fprintf('\tfs(%i) = %f\n', ky, fs(ky));
            end
            fprintf('\tv1 = %f\n', fs(end));
            fprintf('\tE = 1 - exp(-PS/f) = %f\n', 1 - exp(-fs(2)/fs(1)));
            
            fprintf('\tsigma0 = %f\n', this.sigma0);
            for ky = this.map.keys
                fprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})));
            end
        end
        function [k,sk] = k1(this, varargin)
            [k,sk] = find_result(this, 'k1');
        end
        function [k,sk] = k2(this, varargin)
            [k,sk] = find_result(this, 'k2');
        end
        function [k,sk] = k3(this, varargin)
            [k,sk] = find_result(this, 'k3');
        end 
        function [k,sk] = k4(this, varargin)
            [k,sk] = find_result(this, 'k4');
        end        
        function h = plot(this, varargin)
            ip = inputParser;
            addParameter(ip, 'showAif', true, @islogical)
            addParameter(ip, 'xlim', [-10 500], @isnumeric)            
            addParameter(ip, 'ylim', [], @isnumeric)
            addParameter(ip, 'zoom', 4, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.zoom = ipr.zoom;            
            
            RR = mlraichle.RaichleRegistry.instance();
            tBuffer = RR.tBuffer;
            aif = this.dispersedAif(this.artery_interpolated);
            h = figure;
            times = this.times_sampled;
            sampled = this.model.sampled(this.ks, this.model.fs_Raichle_Martin, aif, times);            
            
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
                       ks__, double(this.model.fs_Raichle_Martin), this.artery_interpolated, this.times_sampled, double(this.Measurement), this.timeCliff), ...
                this.ks0, this.ks_lower, this.ks_upper, options); 
            
            this.results_ = struct('ks0', this.ks0, 'ks', ks_, 'sse', sse, 'exitflag', exitflag, 'output', output); 
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
                s = [s sprintf('\tk%i = %f\n', ky, this.ks(ky))]; %#ok<AGROW>
            end
            
            fs = this.model.fs_Raichle_Martin;
            for ky = 1:length(fs)-1
                s = [s sprintf('\tfs(%i) = %f\n', ky, fs(ky))]; %#ok<AGROW>
            end            
            s = [s sprintf('\tv1 = %f\n', fs(end))];            
            s = [s sprintf('\tE = 1 - exp(-PS/f) = %f\n', 1 - exp(-fs(2)/fs(1)))];
            
            s = [s sprintf('\tsigma0 = %f\n', this.sigma0)];
            for ky = this.map.keys
                s = [s sprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})))]; %#ok<AGROW>
            end
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

