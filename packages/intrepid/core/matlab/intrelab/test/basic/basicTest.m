function tests = basicTest
  tests = functiontests(localfunctions);
end

function testIntrepid(testCase)
  act_fdiff = m2i_test;
  exp_fdiff = 0;
  verifyEqual(testCase,act_fdiff,exp_fdiff,'AbsTol',1e-8);
end

function testIntrepidAndML(testCase)
  act_solnorm = m2ml_test;
  exp_solnorm = 2.5888736e+03;
  verifyEqual(testCase,act_solnorm,exp_solnorm,'RelTol',1e-6);
end
