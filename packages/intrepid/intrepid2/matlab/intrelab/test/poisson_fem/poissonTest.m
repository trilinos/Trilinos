function tests = poissonTest
  tests = functiontests(localfunctions);
end

function testConvergence(testCase)
  act_rate = conv_test;
  exp_rate = 1.9;
  verifyGreaterThanOrEqual(testCase,act_rate,exp_rate);
end
