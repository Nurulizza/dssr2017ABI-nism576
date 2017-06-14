classdef nvarTest < matlab.unittest.TestCase
    
    methods (Test)
        function testNVarOutput(testCase)
            explength = 25;
            output = Sensitivity_analysis_MICA_Activation;
            testCase.verifyEqual(length(output),explength);
        end
    end
end