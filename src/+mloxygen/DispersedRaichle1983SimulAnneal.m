classdef DispersedRaichle1983SimulAnneal < mloxygen.Raichle1983SimulAnneal
	%% DISPERSEDRAICHLE1983SIMULANNEAL  

	%  $Revision$
 	%  was created 10-Sep-2020 22:24:15 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.7.0.1434023 (R2019b) Update 6 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    methods  (Static)
        function loss = loss_function(ks, artery_interpolated, times_sampled, measurement, sigma0)
            import mloxygen.DispersedRaichle1983Model.sampled            
            estimation  = sampled(ks, artery_interpolated, times_sampled);
            measurement = measurement(1:length(estimation));
            taus        = diff(times_sampled);
            taus        = [taus taus(end)];
            taus        = taus(1:length(estimation));
            positive    = measurement > 0;
            e           = estimation .* taus;
            m           = measurement .* taus;
            eoverm      = e(positive)./m(positive);
            Q           = sum((1 - eoverm).^2);
            loss        = 0.5*Q/sigma0^2; % + sum(log(sigma0*measurement)); % sigma ~ sigma0*measurement
        end 
    end

	methods		  
 		function this = DispersedRaichle1983SimulAnneal(varargin)
 			this = this@mloxygen.Raichle1983SimulAnneal(varargin{:});
 		end
        
        function [k,sk] = k4(this, varargin)
            [k,sk] = find_result(this, 'k4');
        end 
        function this = solve(this, varargin)
            ip = inputParser;
            parse(ip, varargin{:})
            ipr = ip.Results; %#ok<NASGU>
            
            import mloxygen.DispersedRaichle1983SimulAnneal.loss_function   
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
                @(ks__) loss_function( ...
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
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

