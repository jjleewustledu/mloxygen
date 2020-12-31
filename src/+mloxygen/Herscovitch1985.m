classdef Herscovitch1985 < handle& mlpet.TracerKinetics
	%% HERSCOVITCH1985 
    %  Uses:  asrow()

	%  $Revision$
 	%  was created 15-Jun-2020 19:29:51 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/src/+mloxygen.
 	%% It was developed on Matlab 9.8.0.1396136 (R2020a) Update 3 for MACI64.  Copyright 2020 John Joowon Lee.
 	
    properties (Constant)
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
        o2content
        ocdata
        oodata
        ooFracTime  = 120
        pie % see also:  man pie  
        scaleHO = 1000
        scaleOC = 1000
        scaleOO = 1000
 		sessionPath % contains all data specific to scan session.
        shiftWorldline = false
        startTime
        studyTag
        verbose
    end
    
    properties (Dependent)
        aifHOMetab
        aifOO
        aifOOIntegral
        decayRate
        imagingContextBrainmask
        imagingContextHO
        imagingContextOC
        imagingContextOO
        imagingContextWmparc
        ooPeakTime
        petPointSpread
        pnumber
        videenBlur
        W % well factor
    end
    
    methods (Static)
        function metric = polynomialMetric(As, petobs)
            metric = petobs.^2*As(1) + petobs*As(2);
        end  
        function replaceHeader(in, hdr)
            hdr = mlfourd.ImagingContext2(hdr);
            in  = mlfourd.ImagingContext2(in);
            hdr = hdr.nifti;
            hdr.img = in.nifti.img;
            hdr.fqfilename = in.fqfilename;
            hdr.save()
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
        function g = get.imagingContextBrainmask(this)
            g = mlfourd.ImagingContext2([this.pnumber '_brainmask_reg_2fdg.nii']);
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
        function g = get.petPointSpread(this)
            g = [7.311094309335641 7.311094309335641 5.330000000000000];
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
            %  @param o2content in [0 1] is numeric.
            %  @param verbose is logical and provides quality controls.
            %  @param studyTag is char

            ip = inputParser;
            addParameter(ip, 'sessionPath', pwd, @isfolder)
            addParameter(ip, 'pie', 5.2705, @isnumeric)
            addParameter(ip, 'o2content', 0.2, @isnumeric)
            addParameter(ip, 'verbose', false, @islogical)
            addParameter(ip, 'studyTag', '', @ischar)
            parse(ip, varargin{:})
 			ipr = ip.Results;
            this.sessionPath = ipr.sessionPath;
            this.pie = ipr.pie;
            this.o2content = ipr.o2content;
            this.verbose = ipr.verbose;
            this.studyTag = ipr.studyTag;
            
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
            if isempty(this.a1) && isempty(this.a2)
                this.buildA1A2();
            end            
            petobs = this.tracer2petobs('ho');
            cbf = this.polynomialMetric([this.a1 this.a2], petobs);
            cbf = this.buildImagingContext2(cbf, 'filename', [this.pnumber 'cbf_reg_2fdg.nii']);
            this.imagingContextCbf = cbf;
        end
        function cbv  = buildCbv(this)
            dcv = this.ocdata.dcv;
            tac = this.ocdata.tac;
            petobs = this.tracer2petobs('oc');
            int = trapz(dcv.dcv_well(1:tac.Dur_sec-1));
            cbv = 100*petobs*60*this.pie*tac.Dur_sec/(this.RBC_FACTOR*this.DENSITY_BRAIN*int);
            cbv = this.buildImagingContext2(cbv, 'filename', [this.pnumber 'cbv_reg_2fdg.nii']);
            this.imagingContextCbv = cbv;
        end  
        function cmro2 = buildCmro2(this, varargin)
            ip = inputParser;
            addParameter(ip, 'o2content', this.o2content, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            if isempty(this.imagingContextCbf)
                this.buildCbf()
            end
            if isempty(this.imagingContextCbv)
                this.buildCbv()
            end             
            if isempty(this.imagingContextOef)
                this.buildOef()
            end            
            cmro2 = this.imagingContextCbf .* this.imagingContextOef .* ipr.o2content;
            cmro2.fileprefix = [this.pnumber 'cmro2_reg_2fdg'];
            cmro2.filesuffix = '.nii';
            this.imagingContextCmro2 = cmro2;
        end
        function ic   = buildImagingContext2(this, varargin)
            ip = inputParser;
            addRequired(ip, 'img', @isnumeric)
            addParameter(ip, 'filename', @ischar)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            ifc = this.imagingContextBrainmask.nifti;
            ifc.img = ipr.img;
            ifc.filename = ipr.filename;
            ic = mlfourd.ImagingContext2(ifc);
        end
        function oef  = buildOef(this)
            if isempty(this.imagingContextCbf)
                this.buildCbf()
            end
            if isempty(this.imagingContextCbv)
                this.buildCbv()
            end 
            if isempty(this.b1) && isempty(this.b2)
                this.buildB1B2();
                this.buildB3B4();
            end           
            petobs = this.tracer2petobs('oo');
            nimg = this.buildOefNumer(petobs);
            dimg = this.buildOefDenom();
            oef = this.ensure0to1(nimg./dimg);
            oef = this.buildImagingContext2(oef, 'filename', [this.pnumber 'oef_reg_2fdg.nii']);
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
        function tbl  = buildWholebrain(this, varargin)
            import mlfourd.ImagingContext2
            
            ip = inputParser;
            addParameter(ip, 'checkmask', false, @islogical)
            parse(ip, varargin{:})
            
            %mask = ImagingContext2([this.pnumber 'fdg1_wmparc.nii']);
            mask = ImagingContext2([this.pnumber '_brainmask_reg_2fdg.nii']);
            mask = mask.thresh(40); % CSF intensity
            mask = mask.binarized();
            if ip.Results.checkmask
                mask.fsleyes(this.imagingContextBrainmask.fqfilename)
            end
            CBF1 = ImagingContext2([this.pnumber 'cbf_reg_2fdg.nii']);
            CBF1 = CBF1.volumeAveraged(mask);
            CBF1 = CBF1.nifti.img;
            CBV1 = ImagingContext2([this.pnumber 'cbv_reg_2fdg.nii']);
            CBV1 = CBV1.volumeAveraged(mask);
            CBV1 = CBV1.nifti.img;
            CMRO2_1 = ImagingContext2([this.pnumber 'cmro2_reg_2fdg.nii']);
            CMRO2_1 = CMRO2_1.volumeAveraged(mask);
            CMRO2_1 = CMRO2_1.nifti.img;
            OEF1 = ImagingContext2([this.pnumber 'oef_reg_2fdg.nii']);
            OEF1 = OEF1.volumeAveraged(mask);
            OEF1 = OEF1.nifti.img;
            
            PID = str2double(this.pnumber(2:end));
            HO_shift = this.hodata.Dt;
            OO_shift = this.oodata.Dt;
            tbl = table(PID, CBF1, CBV1, CMRO2_1, OEF1, HO_shift, OO_shift);
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
        function        flirt2atl(this)
            
            % align fdg
            if ~isfile([this.pnumber 'fdg1_avgt.nii'])
                fdg = mlfourd.ImagingContext2([this.pnumber 'fdg1.nii']);
                fdg = fdg.timeAveraged();
                fdg.save
            end            
            std = '/usr/local/fsl/data/standard/MNI152_T1_2mm';
            [~,stdbase] = fileparts(std);
            if ~isfile([this.pnumber 'fdg1_avgt_on_' stdbase '.nii'])
                mlbash(sprintf( ...
                    'flirt -in %sfdg1_avgt.nii -ref %s -out %sfdg1_avgt_on_%s.nii -omat %sfdg1_avgt_on_%s.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear', ...
                    this.pnumber, std, this.pnumber, stdbase, this.pnumber, stdbase))
            end
            
            % align physiological
            hiresT1 = fullfile(getenv('REFDIR'), 'mni_icbm152_t1_tal_nlin_asym_09b_hires.nii');
            refprefix = [this.pnumber 'fdg1_avgt_on_MNI152_T1_2mm'];
            for phy = {'cbf' 'cbv' 'oef' 'cmro2'}
                in  = [this.pnumber phy{1} '_reg_2fdg.nii'];
                out = [this.pnumber phy{1} '_on_MNI152_T1_2mm.nii.gz'];
                this.replaceHeader(in, [this.pnumber 'fdg1_avgt.nii'])
                mlbash(sprintf( ...
                    'flirt -in %s -applyxfm -init %s.mat -out %s -paddingsize 0.0 -interp trilinear -ref %s.nii.gz', ...
                    in, refprefix, out, std))
                mlbash(sprintf('fsleyes %s -cm brain_colours_nih_new %s -a 60', out, hiresT1))
            end
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
                time = ((tac.Start+tac.End)/2-1:(tac.Start+tac.End)/2+1)'; % dither for visibility
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
            this.(dat).Dt = dcv.time(idx_dcv_inflow) - unif_midt(idx_dtac_inflow); % Dt > 0 typically for NNICU ECAT EXACT HR+
            this.(dat).dcv.time = dcv.time - this.(dat).Dt;
            
            % rescale dcv to adjust for alterations of activity worldlines
            if this.shiftWorldline
                this.(dat).dcv.dcv_well = dcv.dcv_well*2^(this.(dat).Dt/this.halflife);
            else
                this.(dat).dcv.dcv_well = dcv.dcv_well;
            end
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
            mask = this.imagingContextBrainmask;
            mask = mask.binarized;
            ic2 = ic2.blurred(this.petPointSpread);
            ic2 = ic2.masked(mask);
            switch tra
                case 'ho' 
                    petCounts = ic2.nifti.img(:,:,:,range)/this.scaleHO;
                    for it = 1:length(midt)
                        petCounts(:,:,:,it) = petCounts(:,:,:,it)*2^(-midt(it)/this.halflife);
                    end                   
                    petobs = trapz(midt, petCounts, 4);
                case 'oc'
                    petobs = ic2.nifti.img/this.scaleOC;
                case 'oo' 
                    petCounts = ic2.nifti.img(:,:,:,range)/this.scaleOO;
                    for it = 1:length(midt)
                        petCounts(:,:,:,it) = petCounts(:,:,:,it)*2^(-midt(it)/this.halflife);
                    end                   
                    petobs = trapz(midt, petCounts, 4);
                otherwise
                    error('mloxygen:ValueError', 'tracer %s is not supported', tra)
            end 
        end 
    end
    
    %% PROTECTED
    
    properties (Access = protected)
        aifHOMetab_
        aifOO_
        aifOOIntegral_
        W_
    end
    
    methods (Access = protected)
        function img  = ensure0to1(~, img)
            %% @return numeric
                      
            img(~isfinite(img)) = 0;
            img(isnan(img)) = 0;
            img(img < 0) = 0;
            img(img > 1) = 0;
        end
        function petdyn = estimatedPetdyn(this, dcv, cbf)
            assert(istable(dcv))
            assert(isnumeric(cbf));
            
            import mloxygen.*;
            f = Herscovitch1985.cbfToF1(cbf);
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
            
            dcv = readtable(fullfile(this.sessionPath, sprintf('%s%s1%s.dcv', this.pnumber, tra, this.studyTag)), 'FileType', 'text');
            try
                tac = readtable(fullfile(this.sessionPath, sprintf('%s_%s1_frame_tac.csv', this.pnumber, tra)));
            catch ME
                handwarning(ME)
                tac = readtable(fullfile(this.sessionPath, sprintf('%s%s1_frame_tac.csv', this.pnumber, tra)));                
            end
            this.startTime = tac.Start(1);
            tac.Start = tac.Start - this.startTime;
            tac.End = tac.End - this.startTime;
            prop = [tra 'data'];
            this.(prop) = struct('dcv', dcv, 'tac', tac);
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

