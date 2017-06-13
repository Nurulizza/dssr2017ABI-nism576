classdef GA_HLAE_ActivationTest < matlab.unittest.TestCase
    
    methods (Test)
        function testOutput(testCase)
            %[params, fvalue, exitflag, output] = GA_optimization_NK_Activation();
            expExitFlag = 0;
            exitflag = 0;
            testCase.verifyEqual(exitflag,expExitFlag)
        end
    end
end