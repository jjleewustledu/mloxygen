classdef Raichle1983SimulAnneal < mlpet.TCSimulAnneal
	%% RAICHLE1983SIMULANNEAL operates on single voxels/regions.
    %  It overloads fprintfModel(), solve() for higher accuracies, sprintfModel().

	%  $Revision$
 	%  was created 10-Sep-2020 19:43:31 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.

    properties
    end
    
	methods
 		function this = Raichle1983SimulAnneal(varargin)
 			%% RAICHLE1983SIMULANNEAL
            %  @param context is mloxygen.Raichle1983.
            %  @param sigma0.
            %  @param fileprefix.

 			this = this@mlpet.TCSimulAnneal(varargin{:});
        end        
        
        function fprintfModel(this)
            fprintf('Simulated Annealing:\n');            
            for ky = 1:length(this.ks)
                fprintf('\tk%i = %f\n', ky, this.ks(ky));
            end           
            PS = this.ks(3);
            f = this.ks(1);
            fprintf('\tE = 1 - exp(-PS/f) = %f\n', this.context.E(PS, f))
            fprintf('\tsigma0 = %f\n', this.sigma0);
            for ky = this.map.keys
                fprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})));
            end
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
                        ks__, this.Data, this.ArteryInterpolated, this.TimesSampled, double(this.Measurement)), ...
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
                s = [s sprintf('\tk%i = %f\n', ky, this.ks(ky))]; %#ok<AGROW>
            end
            PS = this.ks(3);
            f = this.ks(1);
            s = [s sprintf('\tE = 1 - exp(-PS/f) = %f\n', this.context.E(PS, f))];
            s = [s sprintf('\tsigma0 = %f\n', this.sigma0)];
            for ky = this.map.keys
                s = [s sprintf('\tmap(''%s'') => %s\n', ky{1}, struct2str(this.map(ky{1})))]; %#ok<AGROW>
            end
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

