classdef Herscovitch1985 < handle & matlab.mixin.Copyable 
	%% HERSCOVITCH1985 
    %  Uses:  asrow()

	%  $Revision$
 	%  was created 15-Jun-2020 19:29:51 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.8.0.1396136 (R2020a) Update 3 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    properties (Constant)
        LAMBDA = 0.95                % brain-blood equilibrium partition coefficient, mL/mL, Herscovitch, Raichle, JCBFM (1985) 5:65
        DENSITY_BRAIN = 1.05         % assumed mean brain density, g/mL
        RBC_FACTOR = 0.766           % per Tom Videen, metproc.inc, line 193  
        RATIO_SMALL_LARGE_HCT = 0.85 % Grubb, et al., 1978               
        
        CBF_UTHRESH = 200
        CBV_UTHRESH = 100
    end
    
	properties
        a1 % for H2[15O] model
        a2
        b1 % for O[15O] model
        b2
        b3
        b4
        canonFlows = 10:100 % mL/100 g/min
        fracHOMetab = 153/263
        halflife = 122.2416 % sec
        hodata
        imagingContextCbf
        imagingContextCbv
        imagingContextCmro2
        imagingContextOef
        ocdata
        oodata
        ooFracTime  = 120
        petPointSpread = [7.311094309335641 7.311094309335641 5.330000000000000]
        pie % see also:  man pie       
 		sessionPath % contains all data specific to scan session.
        startTime
        verbose
    end
    
    properties (Dependent)
        aifHOMetab
        aifOO
        aifOOIntegral
        decayRate
        imagingContextHO
        imagingContextOC
        imagingContextOO
        imagingContextWmparc
        ooPeakTime
        pnumber
        videenBlur
        W % well factor
    end
    
    methods (Static) 
        function metric = polynomialMetric(As, petobs)
            metric = petobs.^2*As(1) + petobs*As(2);
        end  
        function cbf = invsToCbf(f)
            cbf = 6000*f/mloxygen.Herscovitch1985.DENSITY_BRAIN;
        end
        function f   = cbfToInvs(cbf)
            f = cbf*mloxygen.Herscovitch1985.DENSITY_BRAIN/6000;
        end
    end
    
	methods 
        
        %% GET
        
        function g = get.aifHOMetab(this)
            if isempty(this.aifHOMetab_)
                this.aifHOMetab_ = this.estimateAifHOMetab();
            end
            g = this.aifHOMetab_;
        end
        function g = get.aifOO(this)
            if isempty(this.aifOO_)
                this.aifOO_ = this.estimateAifOO();
            end
            g = this.aifOO_;
        end
        function g = get.aifOOIntegral(this)
            if isempty(this.aifOOIntegral_)
                this.aifOOIntegral_ = this.estimateAifOOIntegral();
            end
            g = this.aifOOIntegral_;
        end
        function g = get.decayRate(this)
            g = log(2)/this.halflife;
        end
        function g = get.imagingContextHO(this)
            g = mlfourd.ImagingContext2([this.pnumber 'ho1_reg_2fdg.nii']);
        end
        function g = get.imagingContextOC(this)
            g = mlfourd.ImagingContext2([this.pnumber 'oc1_reg_2fdg.nii']);
        end
        function g = get.imagingContextOO(this)
            g = mlfourd.ImagingContext2([this.pnumber 'oo1_reg_2fdg.nii']);
        end
        function g = get.imagingContextWmparc(this)
            g = mlfourd.ImagingContext2([this.pnumber 'fdg1_wmparc.nii']);
        end
        function g = get.ooPeakTime(this)
            [~,g] = max(this.oodata.dcv.dcv_well);
        end
        function g = get.pnumber(this)
            [~,g] = fileparts(this.sessionPath);
        end
        function g = get.videenBlur(this)
            fhalf = 0.3; % half-wave number in cm^{-1}; cf. gauss_4dfp
            fwhh  = 10*2*log(2)/(pi*fhalf);
            g     = fwhh*[1 1 1];
        end
        function g = get.W(this)
            if ~isempty(this.W_)
                g = this.W_;
                return
            end
            g = mean(this.hodata.dcv.dcv_well ./ this.hodata.dcv.dcv);
        end
        
        %%
		  
 		function this = Herscovitch1985(varargin)
 			%% HERSCOVITCH1985
 			%  @param sessionPath is folder.
            %  @param pie is numeric from pie-phantom calibration.
            %  @param startTime is seconds, obtained from p1234tra1.4dfp.img.rec from 
            %         $(962to4dfp_framecheck p1234tra1.v)
            %  @param verbose is logical and provides quality controls.

            ip = inputParser;
            addParameter(ip, 'sessionPath', pwd, @isfolder)
            addParameter(ip, 'pie', 5.2705, @isnumeric)
            addParameter(ip, 'verbose', false, @islogical)
            parse(ip, varargin{:})
 			ipr = ip.Results;
            this.sessionPath = ipr.sessionPath;
            this.pie = ipr.pie;
            this.verbose = ipr.verbose;
            
            for tra = {'ho' 'oo' 'oc'}
                this = this.readdata(tra{1});
            end
            this.W_ = mean(this.hodata.dcv.dcv_well ./ this.hodata.dcv.dcv);
        end
        
        function this = buildA1A2(this)
            this = decayUncorrect(this, 'ho');
            this = shiftAif(this, 'ho');
            petobs = this.estimatedPetobs(this.hodata.dcv, this.canonFlows);
            model = this.modelOfCbf(petobs, this.canonFlows);
            this.a1 = model.Coefficients{1, 'Estimate'};
            this.a2 = model.Coefficients{2, 'Estimate'};
        end
        function this = buildB1B2(this)
            this = decayUncorrect(this, 'oo');
            this = shiftAif(this, 'oo');
            flowHOMetab_ = 60*this.pie*this.estimatedPetobs(this.aifHOMetab, this.canonFlows);
            model = this.modelOfOOFlow(flowHOMetab_, this.canonFlows);
            this.b1 = model.Coefficients{1, 'Estimate'};
            this.b2 = model.Coefficients{2, 'Estimate'};
        end
        function this = buildB3B4(this)
            flowOO_ = 60*this.pie*this.estimatedPetobs(this.aifOO, this.canonFlows);
            model = this.modelOfOOFlow(flowOO_, this.canonFlows);
            this.b3 = model.Coefficients{1, 'Estimate'};
            this.b4 = model.Coefficients{2, 'Estimate'};
        end        
        function cbf  = buildCbf(this)
            if isempty(this.a1) || isempty(this.a2)
                this.buildA1A2();
            end            
            petobs = this.tracer2petobs('ho');
            cbf = this.polynomialMetric([this.a1 this.a2], petobs);
            cbf = mlfourd.ImagingContext2(cbf, 'filename', [this.pnumber 'cbf_reg_2fdg.nii']);
            this.imagingContextCbf = cbf;
        end
        function cbv  = buildCbv(this)
            dcv = this.ocdata.dcv;
            tac = this.ocdata.tac;
            petobs = this.tracer2petobs('oc');
            int = trapz(dcv.dcv_well(1:tac.Dur_sec-1));
            cbv = 100*petobs*60*this.pie/(this.RBC_FACTOR*this.DENSITY_BRAIN*int);
            cbv = mlfourd.ImagingContext2(cbv, 'filename', [this.pnumber 'cbv_reg_2fdg.nii']);
            this.imagingContextCbv = cbv;
        end  
        function cmro2 = buildCmro2(this, o2content)
            cmro2 = this.imagingContextCbf .* this.imagingContextOef .* o2content;
            this.imagingContextCmro2 = cmro2;
        end
        function oef  = buildOef(this)
            if isempty(this.b1) || isempty(this.b2)
                this.buildB1B2();
                this.buildB3B4();
            end
            if isempty(this.imagingContextCbf)
                this.buildCbf()
            end
            if isempty(this.imagingContextCbv)
                this.buildCbv()
            end            
            petobs = this.tracer2petobs('ho');
            nimg = this.buildOefNumer(petobs);
            dimg = this.buildOefDenom();
            oef = this.ensure0to1(nimg./dimg);
            oef = mlfourd.ImagingContext2(oef, 'filename', [this.pnumber 'oef_reg_2fdg.nii']);
            this.imagingContextOef = oef;
        end  
        function img  = buildOefNumer(this, petobs) 
            %% @return numeric
            
            cbv = this.imagingContextCbv.nifti.img;
            int = this.estimateAifOOIntegral();
            img = 60*this.pie*petobs - this.flowHOMetab() - int*cbv;
        end
        function img  = buildOefDenom(this)
            %% @return numeric
            
            cbv = this.imagingContextCbv.nifti.img;
            int = this.estimateAifOOIntegral();
            img = this.flowOO() - 0.835*int*cbv;
        end
        function this = decayUncorrect(this, varargin)
            ip = inputParser;
            addRequired(ip, 'tracer', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;            
            
            dat = [ipr.tracer 'data'];
            dcv = this.(dat).dcv;
            tac = this.(dat).tac;
            midt = (tac.Start+tac.End)/2 + this.startTime;
            
            this.(dat).dcv.dcv_well = dcv.dcv_well.*2.^(-dcv.time/this.halflife);
            this.(dat).tac.whole_brain = tac.whole_brain.*2.^(-midt/this.halflife);            
        end
        function dcv  = estimateAifOO(this)
            dcv = this.oodata.dcv;
            dcv.dcv_well = dcv.dcv_well - this.aifHOMetab.dcv_well;
        end
        function dcv  = estimateAifHOMetab(this)
            dcv      = this.oodata.dcv;
            [~,idxP] = max(dcv.time > this.ooPeakTime);
            dfrac_dt = this.fracHOMetab/(this.ooFracTime - dcv.time(idxP));
            fracVec  = zeros(size(dcv.time));
            fracVec(idxP:end) = dfrac_dt*(dcv.time(idxP:end) - dcv.time(idxP));            
            dcv.dcv_well = dcv.dcv_well.*fracVec;
        end
        function int  = estimateAifOOIntegral(this)
            tac = this.oodata.tac;
            midt = (tac.Start+tac.End)/2;
            range = 0 <= midt & midt <= this.oodata.dcv.time(end-2);
            int = trapz(this.aifOO.dcv_well(range));
            int = 0.01*this.RATIO_SMALL_LARGE_HCT*this.DENSITY_BRAIN*int;
        end
        function f    = flowHOMetab(this)
            %% @return numeric flow
            
            if isempty(this.imagingContextCbf)
                this.buildCbf()
            end
            img = this.imagingContextCbf.nifti.img;
            f = this.b1*img.^2 + this.b2*img;
        end
        function f    = flowOO(this)
            %% @return numeric flow
            
            if isempty(this.imagingContextCbf)
                this.buildCbf()
            end
            img = this.imagingContextCbf.nifti.img;
            f = this.b3*img.^2 + this.b4*img;
        end 
        function h    = plot(this, varargin)
            ip = inputParser;
            addRequired(ip, 'tracer', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            dat = [ipr.tracer 'data'];
            dcv = this.(dat).dcv;
            tac = this.(dat).tac;
            if 1 == length(tac.whole_brain)
                time = ((tac.Start+tac.End)/2:dcv.time(end))';
                whole_brain = tac.whole_brain*ones(size(time));
            else
                time = (tac.Start+tac.End)/2;
                whole_brain = tac.whole_brain;
            end
            h = plot(dcv.time, dcv.dcv_well, time, whole_brain);
            xlabel('time / sec')
            ylabel('activity / well-counter cps')
            legend([ipr.tracer ' dcv'], [ipr.tracer ' tac'])
        end      
        function f    = plotall(this)
            f = figure;
            subplot(3, 1, 1)
            plot(this, 'ho')
            subplot(3, 1, 2)
            plot(this, 'oc')
            subplot(3, 1, 3)
            plot(this, 'oo')
        end   
        function this = shiftAif(this, varargin)
            ip = inputParser;
            addRequired(ip, 'tracer', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            dat = [ipr.tracer 'data'];
            dcv = this.(dat).dcv;
            tac = this.(dat).tac;
            midt = (tac.Start+tac.End)/2;
            unif_midt = midt(1):midt(end);
            unif_whole_brain = pchip(midt, tac.whole_brain, unif_midt);
            dtac = diff(unif_whole_brain); % uniformly sampled time-derivative
            
            % shift dcv in time to match inflow with dtac            
            [~,idx_dcv_inflow] = max(dcv.dcv_well > max(dcv.dcv_well)/2);
            [~,idx_dtac_inflow] = max(dtac > max(dtac)/2);
            Dt = dcv.time(idx_dcv_inflow) - unif_midt(idx_dtac_inflow); % Dt > 0 typically for NNICU ECAT EXACT HR+
            this.(dat).dcv.time = dcv.time - Dt;
            
            % rescale dcv to adjust for alterations of activity worldlines
            this.(dat).dcv.dcv_well = dcv.dcv_well*2^(Dt/this.halflife);
        end    
        function petobs = tracer2petobs(this, tra)
            assert(ischar(tra))
            tra = lower(tra);
            dat = [tra 'data'];
            
            tac = this.(dat).tac;
            midt = (tac.Start+tac.End)/2;
            switch tra
                case 'ho'
                    range = 0 <= midt & midt <= this.(dat).dcv.time(end-2);
                    midt = midt(range);
                case 'oc'
                case 'oo'
                    range = 0 <= midt & midt <= this.(dat).dcv.time(end-2);
                    midt = midt(range);
                otherwise
                    error('mloxygen:ValueError', 'tracer %s is not supported', tra)
            end
            
            ict = ['imagingContext' upper(tra)];
            ic2 = this.(ict);
            mask = this.imagingContextWmparc;
            mask = mask.binarized;
            ic2 = ic2.blurred(this.videenBlur);
            ic2 = ic2.masked(mask);
            switch tra
                case 'ho' 
                    petCounts = ic2.nifti.img(:,:,:,range)/1000;
                    for it = 1:length(midt)
                        petCounts(:,:,:,it) = petCounts(:,:,:,it)*2^(-midt(it)/this.halflife);
                    end                   
                    petobs = trapz(midt, petCounts, 4);
                case 'oc'
                    petobs = ic2.nifti.img/7;
                case 'oo' 
                    petCounts = ic2.nifti.img(:,:,:,range)/365;
                    for it = 1:length(midt)
                        petCounts(:,:,:,it) = petCounts(:,:,:,it)*2^(-midt(it)/this.halflife);
                    end                   
                    petobs = trapz(midt, petCounts, 4);
                otherwise
                    error('mloxygen:ValueError', 'tracer %s is not supported', tra)
            end 
        end 
    end
    
    %% PRIVATE   
    
    properties (Access = private)
        aifHOMetab_
        aifOO_
        aifOOIntegral_
        W_
    end
    
    methods (Access = private)
        function img  = ensure0to1(~, img)
            %% @return numeric
                      
            img(~isfinite(img)) = 0;
            img(isnan(img)) = 0;
            img(img < 0) = 0;
            img(img > 1) = 1;
        end
        function petdyn = estimatedPetdyn(this, dcv, cbf)
            assert(istable(dcv))
            assert(isnumeric(cbf));
            
            import mloxygen.*;
            f = Herscovitch1985.cbfToInvs(cbf);
            [~,idx0] = max(dcv.dcv_well' > 2*dcv.dcv_well(1));
            time = dcv.time(idx0:end)'; % row
            activity = dcv.dcv_well(idx0:end)'; % row
            petdyn = zeros(length(f), length(time));
            for r = 1:size(petdyn,1)
                petdyn_ = (1/(60*this.pie))*f(r)*conv(activity, exp(-(f(r)/this.LAMBDA + this.decayRate)*time)); 
                petdyn(r,:) = petdyn_(1:length(time));
            end
        end        
        function petobs = estimatedPetobs(this, dcv, cbf)
            assert(istable(dcv))
            assert(isnumeric(cbf))
            
            rho = this.estimatedPetdyn(dcv, cbf);
            dt = dcv.time(2) - dcv.time(1);
            petobs = dt*trapz(rho, 2);
        end
        function mdl = modelOfCbf(this, petobs, cbf)
            %% modelOfCbf 
            %  @param petobs are numeric PETobs := \int_{t \in \text{obs}} dt' \varrho(t').
            %  @param model from fitnlm()
            %  @returns this with this.product := mdl.  A1, A2 are in mdl.Coefficients{:,'Estimate'}.
            %  See also:  https://www.mathworks.com/help/releases/R2016b/stats/nonlinear-regression-workflow.html
            
            fprintf('Herscovitch1985.modelOfCbf ..........\n');
            mdl = fitnlm( ...
                petobs, ...
                cbf', ...
                @mloxygen.Herscovitch1985.polynomialMetric, ...
                [1 1]);
            disp(mdl)
            fprintf('mdl.RMSE -> %g, min(rho) -> %g, max(rho) -> %g\n', mdl.RMSE, min(petobs), max(petobs));
            if (this.verbose)
                figure; plotResiduals(mdl);
                figure; plotDiagnostics(mdl, 'cookd');
                plotSlice(mdl);
            end
        end
        function mdl = modelOfOOFlow(this, flows, cbf)
            %% modelOfOOFlow 
            %  @param petobs are numeric PETobs := \int_{t \in \text{obs}} dt' \varrho(t').
            %  @param cbf are numeric CBF.
            %  @returns this with this.product := mdl.  A1, A2 are in mdl.Coefficients{:,'Estimate'}.
            %  See also:  https://www.mathworks.com/help/releases/R2016b/stats/nonlinear-regression-workflow.html
            
            fprintf('Herscovitch1985.modelOfOOFlow ..........\n');
            mdl    = fitnlm( ...
                cbf', ...
                flows, ...
                @mloxygen.Herscovitch1985.polynomialMetric, ...
                [1 1]);            
            disp(mdl)
            fprintf('mdl.RMSE -> %g, min(rho) -> %g, max(rho) -> %g\n', mdl.RMSE, min(flows), max(flows));
            if (this.verbose)
                figure; plotResiduals(mdl);
                figure; plotDiagnostics(mdl, 'cookd');
                plotSlice(mdl);
            end
        end
        function this = readdata(this, tra)
            %% assumes that tac acquisition begins with tac.Start(1) == 0
            
            dcv = readtable(fullfile(this.sessionPath, sprintf('p7667%s1.dcv', tra)), 'FileType', 'text');
            tac = readtable(fullfile(this.sessionPath, sprintf('p7667%s1_frame_tac.csv', tra)));
            this.startTime = tac.Start(1);
            tac.Start = tac.Start - this.startTime;
            tac.End = tac.End - this.startTime;
            prop = [tra 'data'];
            this.(prop) = struct('dcv', dcv, 'tac', tac);
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

