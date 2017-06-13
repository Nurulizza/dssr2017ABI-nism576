classdef ODE_HLAE_ActivationTest  < matlab.unittest.TestCase
    methods (Test)
        function testdxdtBasic(testCase)
            prm = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
            X = [0 0 0 0 0 0 0 0 0 0 0 0 0 0];
            unused = 0;
            [dxdt] = ODE_HLAE_Activation(unused, X, prm);
            exp_dxdt_1 = 0;
            testCase.verifyEqual(dxdt(1),exp_dxdt_1)
        end
        function testdxdtOne(testCase)
            prm = [2 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
            X = [4 3.2 767 0 0 0 0 0 0 0 0 0 0 0];
            unused = 0;
            [dxdt] = ODE_HLAE_Activation(unused, X, prm);
            exp_dxdt_1 = 63900;
            testCase.verifyEqual(dxdt(1),exp_dxdt_1)
        end
    end
end