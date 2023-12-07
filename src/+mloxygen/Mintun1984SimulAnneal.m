classdef Mintun1984SimulAnneal < mlpet.TCSimulAnneal
	%% MINTUN1984SIMULANNEAL operates on single voxels/regions. 
    %  It overloads fprintfModel(), solve() for higher accuracies, sprintfModel().

	%  $Revision$
 	%  was created 07-Dec-2020 23:22:23 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.9.0.1524771 (R2020b) Update 2 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    properties
    end
    
	methods		  
 		function this = Mintun1984SimulAnneal(varargin)
 			%% MINTUN1984SIMULANNEAL
            %  @param context is mloxygen.Mintun1984.
            %  @param sigma0.
            %  @param fileprefix.
 			
 			this = this@mlpet.TCSimulAnneal(varargin{:});
        end                
        
        function fprintfModel(this)
            fprintf('%s:\n', stackstr());
            for ky = 1:length(this.ks)
                fprintf('\t%s = %g\n', this.ks_names{ky}, this.ks(ky));
            end 

            raichleks = this.Data.raichleks;
            for ky = 1:length(raichleks)
                fprintf('\t%s = %g\n', mlkinetics.Raichle1983Model.ks_names{ky}, raichleks(ky));
            end
            fprintf('\tE = 1 - exp(-PS/f) = %f\n', 1 - exp(-raichleks(2)/raichleks(1)));

            fprintf('\tloss = %g\n', this.loss())
            keys = natsort(this.map.keys);
            for ky = 1:length(this.ks)
                fprintf('\tmap(''%s'') => %s\n', this.ks_names{ky}, ...
                    join(struct2str(this.map(keys{ky}), orientation='horz')));
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
 			[ks_,loss,exitflag,output] = simulannealbnd( ...
                @(ks__) ipr.loss_function( ...
                        ks__, this.Data, this.ArteryInterpolated, this.TimesSampled, double(this.Measurement)), ...
                        this.ks0, this.ks_lower, this.ks_upper, options);
            
            this.product_ = struct('ks0', this.ks0, 'ks', ks_, 'loss', loss, 'exitflag', exitflag, 'output', output); 
            if ~this.quiet
                fprintfModel(this)
            end
            if this.visualize
                plot(this)
            end
        end   
        function s = sprintfModel(this)
            s = sprintf('%s:\n', stackstr());
            for ky = 1:length(this.ks)
                s = [s sprintf('\t%s = %g\n', this.ks_names{ky}, this.ks(ky))]; %#ok<AGROW>
            end

            raichleks = this.Data.raichleks;
            for ky = 1:length(raichleks)
                s = [s sprintf('\t%s = %g\n', mlkinetics.Raichle1983Model.ks_names{ky}, raichleks(ky))]; %#ok<AGROW>
            end
            s = [s sprintf('\tE = 1 - exp(-PS/f) = %f\n', 1 - exp(-raichleks(2)/raichleks(1)))];

            s = [s sprintf('\tloss = %g\n', this.loss())];
            keys = natsort(this.map.keys);
            for ky = 1:length(this.ks)
                s = [s sprintf('\tmap(''%s'') => %s\n', this.ks_names{ky}, ...
                    join(struct2str(this.map(keys{ky}), orientation='horz')))]; %#ok<AGROW>
            end
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

