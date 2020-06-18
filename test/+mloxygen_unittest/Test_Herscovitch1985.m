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
            cbf = obj.polynomialCbf([obj.a1 obj.a2], petobs);
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
            
            ho = mlfourd.ImagingContext2('p7667ho1_reg_2fdg.nii');
            wmparc = mlfourd.ImagingContext2('p7667fdg1_wmparc.nii');
            mask = wmparc.binarized;
            ho = ho.volumeAveraged(mask);
            petCounts = ho.nifti.img'/1000;
            petCounts = petCounts.*2.^(-midt/obj.halflife);
            %figure; plot(midt, wellCounts); title('test_wholeBrainCbf2')
            
            petobs = trapz(midt(range), petCounts(range));
            cbf = obj.polynomialCbf([obj.a1 obj.a2], petobs);
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
        end
	end

 	methods (TestClassSetup)
		function setupHerscovitch1985(this)
            sessp = fullfile(getenv('HOME'), 'Tmp', 'p7667');
            this.testObj_ = mloxygen.Herscovitch1985('sessionPath', sessp);
 		end
	end

 	methods (TestMethodSetup)
		function setupHerscovitch1985Test(this)
            this.pwd0 = pwd;
            cd(this.testObj_.sessionPath)
 			this.testObj = copy(this.testObj_);
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

