classdef Test_Herscovitch1985 < matlab.unittest.TestCase
	%% TEST_HERSCOVITCH1985 

	%  Usage:  >> results = run(mloxygen_unittest.Test_Herscovitch1985)
 	%          >> result  = run(mloxygen_unittest.Test_Herscovitch1985, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 15-Jun-2020 19:29:51 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mloxygen/test/+mloxygen_unittest.
 	%% It was developed on Matlab 9.8.0.1396136 (R2020a) Update 3 for MACI64.  Copyright 2020 John Joowon Lee.
 	
	properties
        pwd0
 		testObj        
        xlsx
        
        a1 = 1.44003892234215e-06 
        a2 = 0.0039731579544141
        b1 = -15.802457226999099
        b2 = 3.007954577587040e+03
        b3 = -2.866745949342757e+02
        b4 = 6.999099137376797e+04
 	end

	methods (Test)
		function test_ctor(this)
 			disp(this.testObj)
 			disp(this.testObj.hodata)
            disp(this.testObj.hodata.tac)
 			disp(this.testObj.ocdata)
 			disp(this.testObj.oodata)
        end
        function test_W(this)
            this.verifyEqual(this.testObj.W, 16.27000, 'RelTol', 1e-6)
        end
        function test_plotall(this)
            plotall(this.testObj)
        end
        function test_decayUncorrect(this)
            subplot(2,1,1)
            plot(this.testObj, 'ho')
            this.testObj = decayUncorrect(this.testObj, 'ho');
            subplot(2,1,2)
            plot(this.testObj, 'ho')
        end
        function test_shiftAif(this)
            this.testObj = decayUncorrect(this.testObj, 'ho');
            subplot(2,1,1)
            plot(this.testObj, 'ho')
            this.testObj = shiftAif(this.testObj, 'ho');
            subplot(2,1,2)
            plot(this.testObj, 'ho')
        end
        function test_buildA1A2(this)
            this.testObj.verbose = true;
            obj = this.testObj.buildA1A2;
            this.verifyEqual(obj.a1, this.a1, 'RelTol', 0.01);
            this.verifyEqual(obj.a2, this.a2, 'RelTol', 0.01);
        end
        function test_wholeBrainCbf(this)
            obj = this.testObj.buildA1A2;
            tac = obj.hodata.tac;
            midt = (tac.Start+tac.End)/2;
            range = 0 <= midt & midt <= obj.hodata.dcv.time(end-2);
            
            petCounts = tac.whole_brain/1000;
            petCounts = petCounts.*2.^(obj.startTime/obj.halflife);
            %figure; plot(midt, wellCounts); title('test_wholeBrainCbf')
            
            petobs = trapz(midt(range), petCounts(range));
            cbf = obj.polynomialMetric([obj.a1 obj.a2], petobs);
            this.verifyEqual(cbf, 44.066954720038126, 'RelTol', 0.01)
            fprintf('test_wholeBrainCbf.obj.a1->%g\n', obj.a1)
            fprintf('test_wholeBrainCbf.obj.a2->%g\n', obj.a2)
            fprintf('test_wholeBrainCbf.petobs->%g\n', petobs)
        end
        function test_wholeBrainCbf2(this)
            obj = this.testObj.buildA1A2;
            tac = obj.hodata.tac;
            midt = (tac.Start+tac.End)/2;
            range = 0 <= midt & midt <= obj.hodata.dcv.time(end-2);
            
            ho = mlfourd.ImagingContext2([obj.pnumber 'ho1_reg_2fdg.nii']);
            wmparc = mlfourd.ImagingContext2([obj.pnumber 'fdg1_wmparc.nii']);
            mask = wmparc.binarized;
            ho = ho.volumeAveraged(mask);
            petCounts = ho.nifti.img'/1000;
            petCounts = petCounts.*2.^(-midt/obj.halflife);
            %figure; plot(midt, wellCounts); title('test_wholeBrainCbf2')
            
            petobs = trapz(midt(range), petCounts(range));
            cbf = obj.polynomialMetric([obj.a1 obj.a2], petobs);
            this.verifyEqual(double(cbf), 44.066954720038126, 'RelTol', 0.01)
            fprintf('test_wholeBrainCbf.obj.a1->%g\n', obj.a1)
            fprintf('test_wholeBrainCbf.obj.a2->%g\n', obj.a2)
            fprintf('test_wholeBrainCbf2.petobs->%g\n', petobs)
        end
        function test_buildB1B2(this)
            this.testObj.verbose = true;
            obj = this.testObj.buildB1B2;
            this.verifyEqual(obj.b1, this.b1, 'RelTol', 0.01);
            this.verifyEqual(obj.b2, this.b2, 'RelTol', 0.01);
        end
        function test_buildB3B4(this)
            this.testObj.verbose = true;
            obj = this.testObj.buildB3B4;
            this.verifyEqual(obj.b3, this.b3, 'RelTol', 0.01);
            this.verifyEqual(obj.b4, this.b4, 'RelTol', 0.01);
        end
        function test_plotAifHOMetab(this)
            figure;
            obj = this.testObj;
            plot(obj.aifHOMetab.time, obj.aifHOMetab.dcv_well);
        end
        function test_plotAifOO(this)
            figure;
            obj = this.testObj;
            plot(obj.aifOO.time, obj.aifOO.dcv_well);
        end  
        function test_buildCbf(this)
            cbf = this.testObj.buildCbf();
            cbf.fsleyes()
        end
        function test_buildCbv(this)
            cbv = this.testObj.buildCbv();
            cbv.fsleyes()
        end
        function test_buildOef(this)
            oef = this.testObj.buildOef();
            oef.fsleyes()
        end
        function test_buildCmro2(this)
            cmro2 = this.testObj.buildCmro2();
            cmro2.fsleyes()
        end
        function test_buildAll(this)
            this.testObj.buildCmro2();
            
            this.testObj.imagingContextCbf.fsleyes()
            this.testObj.imagingContextCbf.save()
            this.testObj.imagingContextCbv.fsleyes()
            this.testObj.imagingContextCbv.save()
            this.testObj.imagingContextCmro2.fsleyes()
            this.testObj.imagingContextCmro2.save()
            this.testObj.imagingContextOef.fsleyes()
            this.testObj.imagingContextOef.save()
        end
        function test_buildWholebrain(this)
            import mlfourd.ImagingContext2
            o = this.testObj;
            tbl = o.buildWholebrain('checkmask', true);
            disp(tbl)
            
            %     PID       CBF1      CBV1     CMRO2_1       OEF1
            %     ____    ________    _____    ________    _________
            % 
            %     7667    42.80746    5.358    2.485797    0.3601757
            
        end
        function test_buildCohort(this)
            cd('/Users/jjlee/Tmp/population_data_Herscovitch')
            tbl = [];
            for g = globFoldersT('p*')
                pwd1 = pushd(g{1});
                gnum = str2double(g{1}(2:end));
                o2c = this.xlsx.O2Cont(this.xlsx.PID == gnum);
                pie = this.xlsx.x2DPie(this.xlsx.PID == gnum);
                o = mloxygen.Herscovitch1985('sessionPath', pwd, 'pie', pie, 'o2content', o2c);    
                if abs(o.W - this.xlsx.wellcal(this.xlsx.PID == gnum)) < 0.2
                    warning('mloxygen:ValueWarning', ...
                        'discrepancy this.W -> %g', o.W - this.xlsx.wellcal(this.xlsx.PID == gnum))
                end
                o.buildCmro2();
                o.imagingContextCbf.save()
                o.imagingContextCbv.save()
                o.imagingContextCmro2.save()
                o.imagingContextOef.save()                
                tbl = [tbl; o.buildWholebrain]; %#ok<AGROW>    
                disp(tbl)
                popd(pwd1);
            end
            writetable(tbl, 'test_buildCohort_Herscovitch.xlsx')
        end
        function test_buildCohort_Chaojie(this)
            cd('/Users/jjlee/Tmp/population_data_Chaojie')
            tbl = [];
            for g = globFoldersT('p*')
                pwd1 = pushd(g{1});
                gnum = str2double(g{1}(2:end));
                o2c = this.xlsx.O2Cont(this.xlsx.PID == gnum);
                pie = this.xlsx.x2DPie(this.xlsx.PID == gnum);
                o = mloxygen.Chaojie2020('sessionPath', pwd, 'pie', pie, 'o2content', o2c, 'studyTag', '_pop');    
                o.buildCmro2();
                o.imagingContextCbf.save()
                o.imagingContextCbv.save()
                o.imagingContextCmro2.save()
                o.imagingContextOef.save()                
                tbl = [tbl; o.buildWholebrain]; %#ok<AGROW>    
                disp(tbl)
                popd(pwd1);
            end
            writetable(tbl, 'test_buildCohort_Chaojie.xlsx')
        end
        function test_buildSingle(this)
            cd('/Users/jjlee/Tmp/population_data_Herscovitch')
            for g = {'p7667'}
                pwd1 = pushd(g{1});
                gnum = str2double(g{1}(2:end));
                o2c = this.xlsx.O2Cont(this.xlsx.PID == gnum);
                pie = this.xlsx.x2DPie(this.xlsx.PID == gnum);
                o = mloxygen.Herscovitch1985('sessionPath', pwd, 'pie', pie, 'o2content', o2c);
                o.buildCmro2();
                o.imagingContextCbf.save()
                o.imagingContextCbv.save()
                o.imagingContextCmro2.save()
                o.imagingContextOef.save()
                tbl = o.buildWholebrain('checkmask', true);
                disp(tbl)
                
                %     PID       CBF1      CBV1     CMRO2_1       OEF1   
                %     ____    ________    _____    ________    _________
                % 
                %     7667    42.80746    5.358    2.485797    0.3601757
    
                popd(pwd1);
            end            
        end
        function test_buildSingle_Chaojie(this)
            cd('/Users/jjlee/Tmp/population_data_Chaojie')
            for g = {'p7667'}
                pwd1 = pushd(g{1});
                gnum = str2double(g{1}(2:end));
                o2c = this.xlsx.O2Cont(this.xlsx.PID == gnum);
                pie = this.xlsx.x2DPie(this.xlsx.PID == gnum);
                o = mloxygen.Chaojie2020('sessionPath', pwd, 'pie', pie, 'o2content', o2c, 'studyTag', '_pop');
                o.buildCmro2();
                o.imagingContextCbf.save()
                o.imagingContextCbv.save()
                o.imagingContextCmro2.save()
                o.imagingContextOef.save()
                tbl = o.buildWholebrain('checkmask', true);
                disp(tbl)
                
                %     PID       CBF1        CBV1      CMRO2_1       OEF1       HO_shift    OO_shift
                %     ____    ________    ________    ________    _________    ________    ________
                % 
                %     7667    29.31678    5.546002    1.415667    0.2920172       7           -2   
    
                popd(pwd1);
            end            
        end
        function test_tracerScalingInNativeSpaces(this)
            tra = 'fdg';
            
            % assemble smooth dynamic tracer962 and make time-summed
            tracerfp = [this.testObj.pnumber tra '1'];
            tracer962 = mlfourd.ImagingContext2([tracerfp '.4dfp.hdr']);
            tracer962.fsleyes()
            tracer962.fileprefix = [tracer962.fileprefix '_962'];
            tracer962.filesuffix = '.nii';
            tracer962 = tracer962.timeSummed();
            tracer962.fsleyes()
            tracer962.save()
            
            % make dynamic tracernii time-summed            
            tracernii = mlfourd.ImagingContext2([tracerfp '.nii']);
            tracernii = tracernii.timeSummed();
            tracernii.fsleyes()
            tracernii.save()
            
            %  plot vectors of time-summed for inspecting correlations
            tracer962_vec = reshape(tracer962.nifti.img, [128*128*63 1]);
            tracernii_vec = reshape(tracernii.nifti.img, [128*128*63 1]);
            plot(tracer962_vec, tracernii_vec, '.')
            ylabel('PETobs NIfTI')
            xlabel('PETobs 962to4dfp_framecheck')
            title([upper(tra) ' time-summed'])
            saveFigures(pwd, 'prefix', [this.testObj.pnumber '_' upper(tra) '_time_summed'], 'closeFigure', false)
        end
        function test_tacScaling(this)
            tra = 'oo';
            dat = [tra 'data'];
            tac = this.testObj.(dat).tac;
            
            % assemble smooth dynamic tracer962 and make time-summed
            tracerfp = [this.testObj.pnumber tra '1'];
            tracer962 = mlfourd.ImagingContext2([tracerfp '.4dfp.hdr']);
            %tracer962.fsleyes()
            if ~strcmp(tra, 'oc')
                peek = tracer962.volumeAveraged();
                plot(tac.Start, peek.fourdfp.img)
                tmp = tracer962.fourdfp;
                try
                    assert(length(tac.Dur_sec) == size(tracer962, 4))
                    for t = 1:length(tac.Dur_sec)
                        tmp.img(:,:,:,t) = tmp.img(:,:,:,t)/tac.Dur_sec(t);
                    end
                catch ME
                    handwarning(ME)
                    warning('mloxygen:RuntimeWarning', 'trying tac.Dur since tac.Dur_sec not found')
                    assert(length(tac.Dur) == size(tracer962, 4))
                    for t = 1:length(tac.Dur)
                        tmp.img(:,:,:,t) = tmp.img(:,:,:,t)/tac.Dur(t);
                    end
                end
                tracer962 = mlfourd.ImagingContext2(tmp);
                peek = tracer962.volumeAveraged();
                plot(tac.Start, peek.fourdfp.img)
                tracer962 = tracer962.timeSummed();
            end
            tracer962.fsleyes()
            
            % make dynamic tracernii time-summed            
            tracernii = mlfourd.ImagingContext2([tracerfp '_reg_2fdg.nii']);
            %tracernii.fsleyes()
            tracernii = tracernii.timeSummed();
            tracernii.save()
            tracernii.fsleyes()
            
            % register tracer962 to tracernii, all time-summed
            tracer962.fileprefix = [tracerfp '_962_sumt']; tracer962.filesuffix = '.nii'; tracer962.save();
            mlbash(sprintf( ...
                'flirt -in %s -ref %s_reg_2fdg_sumt.nii -out %s_reg_2fdg.nii.gz -omat %s_reg_2fdg.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6  -interp trilinear', ...
                tracer962.filename, tracerfp, tracer962.fileprefix, tracer962.fileprefix))
            tracer962 = mlfourd.ImagingContext2([tracer962.fileprefix '_reg_2fdg.nii.gz']);
            %tracer962.fsleyes()
            
            %  plot vectors of time-summed for inspecting correlations
            tracer962_vec = reshape(tracer962.nifti.img, [128*128*63 1]);
            tracernii_vec = reshape(tracernii.nifti.img, [128*128*63 1]);
            plot(tracer962_vec, tracernii_vec, '.')
            ylabel('PETobs NIfTI')
            xlabel('PETobs 962to4dfp_framecheck')
            title([upper(tra) ' time-summed'])
            saveFigures(pwd, 'prefix', [this.testObj.pnumber '_' upper(tra) '_time_summed'], 'closeFigure', false)
        end
        function test_flirt2atl(this)
            this.testObj.flirt2atl()
        end
	end

 	methods (TestClassSetup)
		function setupHerscovitch1985(this)
            sessp = fullfile(getenv('HOME'), 'Tmp', 'population_data_Chaojie', 'p7667');
            this.testObj_ = mloxygen.Chaojie2020('sessionPath', sessp, 'pie', 5.2705, 'o2content', 0.167, 'studyTag', '_pop');
            %this.testObj_ = mloxygen.Herscovitch1985('sessionPath', sessp, 'pie', 5.2705, 'o2content', 0.167, 'studyTag', '');
            %sessp = fullfile(getenv('HOME'), 'Tmp', 'p7757');
            %this.testObj_ = mloxygen.Herscovitch1985('sessionPath', sessp, 'pie', 5.35, 'o2content', 0.161);
            %this.testObj_ = [];
            
            this.xlsx = readtable('/Users/jjlee/Box/Chaojie/np674_de-id_jjlee.xlsx');
 		end
	end

 	methods (TestMethodSetup)
		function setupHerscovitch1985Test(this)
            this.pwd0 = pwd;
            if ~isempty(this.testObj_)
                cd(this.testObj_.sessionPath)
                this.testObj = copy(this.testObj_);
            end
 			this.addTeardown(@this.cleanTestMethod);
 		end
	end

	properties (Access = private)
 		testObj_
 	end

	methods (Access = private)
		function cleanTestMethod(this)
            cd(this.pwd0)
 		end
	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

