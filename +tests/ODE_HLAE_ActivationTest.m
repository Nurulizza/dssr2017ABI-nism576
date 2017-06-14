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
        function testdxdtTwo(testCase)
            prm = [2 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
            X = [4 3.2 767 0 0 0 0 0 0 0 0 0 0 0];
            unused = 0;
            [dxdt] = ODE_HLAE_Activation(unused, X, prm);
            exp_dxdt_2 = 63900;
            testCase.verifyEqual(dxdt(2),exp_dxdt_2)
        end
        function testdxdtThree(testCase)
            prm = [2 3 2 3 2 0 0 0 0 0 0 0 0 0 0 0 ];
            X = [4 3.2 767 0 0 0 5.6 8.2 0 0 0 14.7 0];
            unused = 0;
            [dxdt] = ODE_HLAE_Activation(unused, X, prm);
            exp_dxdt_3 = -1184970;
            testCase.verifyEqual(dxdt(3),exp_dxdt_3)
        end
        function testdxdtFour(testCase)
            prm = [2 3 2 3 2 1 1 2 1 3 0 0 0 0 0 0 ];
            X = [4 3.2 767 56.7 2 0 5.6 8.2 0 0 10 14.7 22];
            unused = 0;
            [dxdt] = ODE_HLAE_Activation(unused, X, prm);
            exp_dxdt_4 = -15268;

            testCase.verifyEqual(dxdt(4),exp_dxdt_4)
        end
    end
end