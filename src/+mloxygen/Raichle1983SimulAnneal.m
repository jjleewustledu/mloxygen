classdef Raichle1983SimulAnneal < mlpet.TracerSimulatedAnneal & mloxygen.Raichle1983Strategy
	%% RAICHLE1983SIMULANNEAL operates on single voxels/regions.

	%  $Revision$
 	%  was created 10-Sep-2020 19:43:31 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.

	methods
 		function this = Raichle1983SimulAnneal(varargin)
 			%% RAICHLE1983SIMULANNEAL
            %  @param context is mloxygen.Raichle1983.
            %  @param sigma0.
            %  @param fileprefix.

 			this = this@mlpet.TracerSimulatedAnneal(varargin{:});
            
            [this.ks_lower,this.ks_upper,this.ks0] = remapper(this);
            this.artery_interpolated = this.model.artery_interpolated;
        end        
        
        function disp(this)
            fprintf('\n')
            fprintf(class(this))
            if isempty(this.results_)
                return
            end
            fprintf('initial ks0: '); disp(this.results_.ks0)
            fprintf('est.     ks: '); disp(this.results_.ks)
            fprintf('        sse: '); disp(this.results_.sse)
            fprintf('   exitflag: '); disp(this.results_.exitflag)
            disp(this.results_.output)
            disp(this.results_.output.rngstate)
            disp(this.results_.output.temperature)
        end
        function fprintfModel(this)
            fprintf('Simulated Annealing:\n');
            %fprintf('\tE = 1 - exp(-PS/f) = %f\n', 1 - exp(-this.ks(2)/this.ks(1)))
            for ky = 1:length(this.ks)
                fprintf('\tk%i = %f\n', ky, this.ks(ky));
            end
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
            addParameter(ip, 'xlim', [-5 500], @isnumeric)            
            addParameter(ip, 'ylim', [], @isnumeric)
            addParameter(ip, 'zoom', 1, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            aif = this.artery_interpolated;
            h = figure;
            times = this.times_sampled;
            sampled = this.model.sampled(this.ks, aif, times);
            if ipr.showAif
                plot(times, ipr.zoom*this.Measurement, ':o', ...
                    times(1:length(sampled)), ipr.zoom*sampled, '-', ...
                    0:length(aif)-1, aif, '--')                
                legend('measurement', 'estimation', 'aif')
            else
                plot(times, ipr.zoom*this.Measurement, 'o', ...
                    times(1:length(sampled)), ipr.zoom*sampled, '-')                
                legend('measurement', 'estimation')
            end
            if ~isempty(ipr.xlim); xlim(ipr.xlim); end
            if ~isempty(ipr.ylim); ylim(ipr.ylim); end
            xlabel('times / s')
            ylabel('activity / (Bq/mL)')
            annotation('textbox', [.175 .25 .3 .3], 'String', sprintfModel(this), 'FitBoxToText', 'on', 'FontSize', 7, 'LineStyle', 'none')
            title('Raichle1983SimulAnneal.plot()')
        end 
        function        save(this)
            save([this.fileprefix '.mat'], this);
        end
        function        saveas(this, fn)
            save(fn, this);
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
                       ks__, this.artery_interpolated, this.times_sampled, double(this.Measurement), this.sigma0), ...
                this.ks0, this.ks_lower, this.ks_upper, options); 
            
            this.results_ = struct('ks0', this.ks0, 'ks', ks_, 'sse', sse, 'exitflag', exitflag, 'output', output); 
            if ~this.quiet
                fprintfModel(this)
            end
            if this.visualize
                plot(this)
            end
        end       
        function s    = sprintfModel(this)
            s = sprintf('Simulated Annealing:\n');
            %s = [s sprintf('\tE = 1 - exp(-PS/f) = %f\n', 1 - exp(-this.ks(2)/this.ks(1)))];
            for ky = 1:length(this.ks)
                s = [s sprintf('\tk%i = %f\n', ky, this.ks(ky))]; %#ok<AGROW>
            end
            s = [s sprintf('\tsigma0 = %f\n', this.sigma0)];
            for ky = this.map.keys
                s = [s sprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})))]; %#ok<AGROW>
            end
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

