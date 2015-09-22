import sys

sys.path.insert(0,"../build/build.dir/packages/PyTrilinos/src/stk/PyPercept")

from math import *
from numpy import *
import unittest
import time
import print_table 

from PerceptMesh import *

class StringFunctionUnitTests(unittest.TestCase):

    def setUp(self):
        self.testpoints = [ [0.1234,  -0.5678, 0.9, 0.812],
                          [0.1234e-3,  -0.5678e-5, 0.9e+8, 0.812e-4],
                          [.101, 102., 10201.0, 0.0122],
                          [0.003, -100001.1, 44.1, 3.0]
                        ]

        self.testpoints_fd = [ [0.1234,  -0.5678, 0.9, 0.812],
                            [0.1234e-3,  -0.5678e-5, 0.9e-3, 0.812e-4],
                            [101.0, 102.0, 10.2, 0.0122],
                            [0.003, .002, -0.0011, 0.0]
                          ]
   
    def test_stringFunction_xy_basic(self):
      x=1.234
      y=2.345
      z=0.0

      sf = StringFunction(" x - y ")
      
      input_array = array([x, y, z])
      
      time = 0.0
      output_array = sf.value(input_array, time)
      print output_array

      eval_print(x, y, z, time, sf)

    def test_stringFunction_xy_basic_1(self):
       sfx = StringFunction("x")
       sfy = StringFunction("y")
       sfxy = StringFunction("x-y")

       x = 1.234
       y = 5.678
       z = 0.0
       t = 0.0
       xy = x-y

       eval_print(1,2,3,0, sfxy)
       vx = eval_func(x, y, z, t, sfx)
       print "x = ", x, "vx = ", vx
       vy = eval_func(x, y, z, t, sfy)
       vxy = eval_func(x, y, z, t, sfxy)
       
       print "y = ", y, "vy = ", vy
       print "xy = ", xy, "vxy = ", vxy

       self.assertEqual(y, vy)
       self.assertEqual(xy, vxy)

    def test_stringFunction_xy_basic_2(self):
       sftestNA = StringFunction("x", "sftestNA", Dimensions(3), Dimensions(2, 3))
       sftestNA.setDomainDimensions(Dimensions(3))
       sftest = StringFunction("x", "sftestNA", Dimensions(3), Dimensions(2, 3))
       sftest_domain = sftest.getNewDomain()
       sftest_codomain = sftest.getNewCodomain()

       sfx = StringFunction("x", "sfx")

       sfy = StringFunction("y")
       sfxy = StringFunction("x-y")

       x = 1.234
       y = 5.678
       z = 0.0
       t = 0.0
       xy = x-y

       eval_print(1,2,3,0, sfxy)
       vx = eval_func(x,y,z,t, sfx)
       print "x = ", x, "vx = ", vx
       vy = eval_func(x, y, z, t, sfy)
       print "y = ", y, "vy = ", vy
       vxy = eval_func(x, y, z, t, sfxy)
       print "xy = ", xy, "vxy = ", vxy
          
       self.assertEqual(x, vx)
       self.assertEqual(y, vy)
       self.assertEqual(xy, vxy)

    def test_stringFunction_test_alias(self):
       sfx = StringFunction("x", "sfx", Dimensions(3), Dimensions(1))
       sfy = StringFunction("y", "sfy", Dimensions(3), Dimensions(1))
       sfxy = StringFunction("x-y", "sfxy", Dimensions(3), Dimensions(1))
       sfembedded = StringFunction("sfxy", "sfembedded", Dimensions(3), Dimensions(1))

       x = 1.234
       y = 5.678
       z = 0.0
       t = 0.0
       xy = x-y

       eval_print(1,2,3,0, sfxy)
       vx = eval_func(x,y,z,t, sfx)
       print "x = ", x, "vx = ", vx
       vy = eval_func(x, y, z, t, sfy)
       print "y = ", y, "vy = ", vy
       vxy = eval_func(x, y, z, t, sfxy)
       print "xy = ", xy, "vxy = ", vxy

       self.assertEqual(x, vx)
       self.assertEqual(y, vy)
       self.assertEqual(xy, vxy)

       print "sfembedded = ...", sfembedded
       eval_print(1,2,3,0,sfembedded)
       print "sfembedded = ", eval_func(x,y,z,t,sfembedded)

       vxy1 = eval_func(x,y,z,t,sfembedded)
       sfembedded.add_alias("sfalias")
       sftestalias = StringFunction("sfalias", "sftestalias")
       vxy2 = eval_func(x,y,z,t,sftestalias)
       print "sftestalias = ", vxy2

    def test_stringFunction_vector_valued(self):
       x = 1.234
       y = 5.678
       z = 3.456
       t = 0.0
      
       didCatch = 0
       try:
          sfv0 = StringFunction("v[0]=x; v[1]=y; v[2]=z; x", "sfv", Dimensions(1,4), Dimensions(1,3))
          eval_vec3_print(1,2,3,0, sfv0)
       except:
          didCatch = 1
          print "TEST::function::stringFunctionVector: expected to catch this since dom/codomain dimensions should be rank-1"
       sfv = StringFunction("v[0]=x*y*z; v[1]=y; v[2]=z; x", "sfv", Dimensions(3), Dimensions(3))
       eval_vec3_print(1.234, 2.345e-3, 3.456e+5, 0.0, sfv)
       vec = eval_vec3(x,y,z,t,sfv)
       print "x = ", x
       print "y = ", y
       print "z = ", z
       print "val = ", (vec[0]*vec[1]*vec[2]) 
          
       self.assertEqual(vec[0], (x*y*z))
       self.assertEqual(vec[1], (y))
       self.assertEqual(vec[2], (z))

    def test_stringFunction_constants(self):
       x = 1.234
       y = 5.678
       z = 3.456
       t = 0.0
       myC = 4.5678

       # dummy return value sf_myC, but could be used in table printing, or other pythonic uses
       sf_myC = StringFunction(str(myC), "myC", Dimensions(3), Dimensions(1)); 

       # alternative
       # sf_myC = StringFunction("4.5678", "myC", Dimensions(3), Dimensions(1));

       # this string function refers to the other through "myC"
       sfv = StringFunction("x+myC", "sfv", Dimensions(3), Dimensions(1))
      
       #eval_print(x,y,z, 0.0, sfv)
       vec = eval_func(x,y,z,t,sfv)
       print "x = ", x
       print "y = ", y
       print "z = ", z
       print "constants test val = ", vec, " expected = ", (myC + x)
         
       self.assertEqual(vec, (myC + x))

       # more...
       myConstants = {"C":1.234,"rho":1.e-5}
       sf_myC1 = []
       for cname, cvalue in myConstants.items():   # note: this could become a python function
           sf_myC1.append( StringFunction(str(cvalue),cname,Dimensions(3),Dimensions(1)) )
       
       sfv1 = StringFunction("x + C*rho", "sfv1", Dimensions(3), Dimensions(1))
      
       #eval_print(x,y,z, 0.0, sfv1)
       vec = eval_func(x,y,z,t,sfv1)
       expected = (x + myConstants["C"]*myConstants["rho"])
       print "constants test val1 = ", vec, " expected = ", expected
         
       self.assertEqual(vec, expected)

    def test_stringFunction_arithmetic_ops(self):
       for xyzt in self.testpoints:
         x = xyzt[0]
         y = xyzt[1]
         z = xyzt[2]
         t = xyzt[3]

         sfx = StringFunction("x")
         sfy = StringFunction("y")
         sfxy = StringFunction("x-y")

         sfxy2 = sfx - sfy
         xy = x - y
         vxy = eval_func(x,y,z,t,sfxy)
         vxy2 = eval_func(x,y,z,t,sfxy2)

         sfx1 = StringFunction("x")
         sfy2 = StringFunction("y")
         vx = eval_func(x,y,z,t, sfx1)
         vy = eval_func(x,y,z,t, sfy2)

         print "vxy2 = ", vxy2, " == vxy = ", vxy
         print "xy = ", xy, " == vxy = ", vxy
         print "x = ", x, " == vx = ", vx
         print "y = ", y, " == y = ", vy

         self.assertEqual(x, vx)
         self.assertEqual(y, vy)
         self.assertEqual(xy, vxy)
         self.assertEqual(vxy2, vxy)

         sfxy_minus = sfx - sfy
         xy_minus = x - y
         vxy_minus = eval_func(x,y,z,t,sfxy_minus)
         vxy1_minus = eval_func(x,y,z,t,sfxy_minus)
         print "xy_minus = ", xy_minus, " == vxy_minus = ", vxy_minus
         print "xy_minus = ", xy_minus, " == vxy1_minus = ", vxy1_minus

         self.assertEqual(xy_minus, vxy_minus)
         self.assertEqual(vxy_minus, vxy1_minus)

         sfxy_plus = sfx + sfy
         xy_plus = x + y
         vxy_plus = eval_func(x,y,z,t,sfxy_plus)
         vxy1_plus = eval_func(x,y,z,t,sfxy_plus)
         print "xy_plus = ", xy_plus, " == vxy_plus = ", vxy_plus
         print "xy_plus = ", xy_plus, " == vxy1_plus = ", vxy1_plus

         self.assertEqual(xy_plus, vxy_plus)
         self.assertEqual(vxy_plus, vxy1_plus)

         sfxy_mult = sfx * sfy
         xy_mult = x * y
         vxy_mult = eval_func(x,y,z,t,sfxy_mult)
         vxy1_mult = eval_func(x,y,z,t,sfxy_mult)
         print "xy_mult = ", xy_mult, " == vxy_mult = ", vxy_mult
         print "xy_mult = ", xy_mult, " == vxy1_mult = ", vxy1_mult

         self.assertEqual(xy_mult, vxy_mult)
         self.assertEqual(vxy_mult, vxy1_mult)

         sfxy_div = sfx / sfy
         xy_div = x / y
         vxy_div = eval_func(x,y,z,t,sfxy_div)
         vxy1_div = eval_func(x,y,z,t,sfxy_div)
         print "xy_div = ", xy_div, " == vxy_div = ", vxy_div
         print "xy_div = ", xy_div, " == vxy1_div = ", vxy1_div

         self.assertEqual(xy_div, vxy_div)
         self.assertEqual(vxy_div, vxy1_div)

    def test_stringFunction_derivative(self):
       for xyzt in self.testpoints:
         x = xyzt[0]
         y = xyzt[1]
         z = xyzt[2]
         t = xyzt[3]

         sfxy = StringFunction("x-y")
         dsfxy_y = StringFunction("-1")
         dy = array([["y"]])
         #input_array = array([x, y, z])
         print "dy= " , dy , " dy.ndim= " , dy.ndim, " dy.dtype= " , dy.dtype, " dy.itemsize= ", dy.itemsize , " dy.size= " , dy.size
         #sys.exit(1)
         dsfxy_y_1 = sfxy.derivative_test(dy)
         
         dvxy = eval_func(x,y,z,t,dsfxy_y_1)
         dvxy1 = eval_func(x,y,z,t,dsfxy_y)
         print "dvxy = ", dvxy, " == dvxy1 = ", dvxy1
         print "-1.0 = -1 == dvxy = ", dvxy 

         self.assertEqual(dvxy, dvxy1)
         self.assertEqual(-1, dvxy)
         print dsfxy_y_1 

    def test_stringFunction_derivative_1(self):
       for xyzt in self.testpoints:
         x = xyzt[0]
         y = xyzt[1]
         z = xyzt[2]
         t = xyzt[3]

         print "here 1"
         eps = 1.e-6
         eps_loc = eps*(fabs(x)+fabs(y)+fabs(z)+fabs(t))/4.0

         sfxy = StringFunction("x-y")
         dsfxy_grad = StringFunction("v[0]=1; v[1]= -1; v[2]=0", "test", Dimensions(3), Dimensions(3))
         dxyz = array([["x"],["y"],["z"]])    #new simpler user-interface 
         #dxyz = array([["x","y","z"]])    #new simpler user-interface 
         print "dxyz.shape= " , dxyz.shape

         grad = array(["1","-1","0"])
         sfxy.set_gradient_strings(grad)
         dsfxy_grad_1 = sfxy.derivative_test(dxyz)
         dsfxy_grad_fd = sfxy.derivative_test_fd(dxyz, eps_loc)
         dsfxy_grad_2 = sfxy.derivative(dxyz)

         dvxy1 = eval_vec3(x,y,z,t,dsfxy_grad_1)
         dvxy_fd = eval_vec3(x,y,z,t,dsfxy_grad_fd)
         dvxy2 = eval_vec3(x,y,z,t,dsfxy_grad_2)
         dvxy = eval_vec3(x,y,z,y,dsfxy_grad)

         i = 0
         while i < 3:
           self.assertEqual(dvxy[i], dvxy1[i])
           self.assertEqual(dvxy[i], dvxy2[i])
           self.assertAlmostEqual(dvxy[i], dvxy_fd[i])
           i = i + 1            

         self.assertEqual(dvxy[0], 1.0)
         self.assertEqual(dvxy[1], -1.0)

    def test_stringFunction_derivative_2(self):    
       for xyzt in self.testpoints_fd:
         x = xyzt[0]
         y = xyzt[1]
         z = xyzt[2]
         t = xyzt[3]
         eps = 1.e-10
         eps_loc = eps*(fabs(x)+fabs(y)+fabs(z)+fabs(t))/4.0
         sf = StringFunction(" sin(x*y*z*z) " )
         grad = array(["y*z*z*cos(x*y*z*z)", "x*z*z*cos(x*y*z*z)", "2*x*y*z*cos(x*y*z*z)"])
         gradv = "v[0]="+grad[0]+"; v[1]="+grad[1]+" ; v[2]="+grad[2]+";"
         dsf_grad = StringFunction(gradv, "test", Dimensions(3), Dimensions(3))
         #dxyz = array([["x","y","z"]])
         dxyz = array([["x"],["y"],["z"]])    #new simpler user-interface 
         sf.set_gradient_strings(grad)
         dsf_grad_fd = sf.derivative_test_fd(dxyz, eps_loc)
         dsf_grad_2 = sf.derivative(dxyz)

         dv_fd = eval_vec3(x,y,z,t,dsf_grad_fd)
         dv2 = eval_vec3(x,y,z,t,dsf_grad_2)
         dv = eval_vec3(x,y,z,t,dsf_grad)

         i = 0
         while i < 3:
           print "dv2[i] = ", dv2[i], " == dv[i] = ", dv[i]
           self.assertEqual(dv[i], dv2[i])
           if fabs(dv[i]-dv_fd[i]) > 0.5*(fabs(dv_fd[i])+fabs(dv[i]))*1.e-6:
             print "\n i = ", i, "x= ", x, "y= ", y, "z= ", z, "expected= ", dv[i], "actual = ", dv_fd[i]
           self.assertAlmostEqual(dv[i], dv_fd[i], delta = 1.e-1)
           i = i + 1  

    def test_stringFunction_multiplePoints(self):
       points = zeros(shape=(4,3))
       output = zeros(shape=(4,1))
       output_expect = zeros(shape=(4,1))

       sf1 = StringFunction("x+y*z")
       i = 0
       for xyzt in self.testpoints:
         x = xyzt[0]
         y = xyzt[1]
         z = xyzt[2]
         t = xyzt[3]
         points[i][0] = x
         points[i][1] = y
         points[i][2] = z
         vx = eval_func(x,y,z,t,sf1)
         output_expect[i][0] = vx 
         print "x+y*z = ", x+y*z, " == vx = ", vx
         self.assertEqual((x+y*z), vx)
         i = i + 1

       sf2 = StringFunction(str(sf1.getFunctionString()), "sf2", Dimensions(3), Dimensions(1))
       output = sf2.value(points, output, 0.0)
       
       i = 0 
       while i < 4:
         print "output_expect(i, 0)  = ", output_expect[i][0] , " == output(i, 0) = ", output[i][0]
         self.assertEqual(output_expect[i][0], output[i][0])
         i = i + 1
         
    def test_stringFunction_expressions(self):
       x = 0.1234
       y = -0.5678
       z = 0.9
       t = 0.812
       global PI
       global E
       PI = pi
       E = e

       sf1 = StringFunction("x+y")
       ve = x + y
       v = eval_func(x,y,z,t,sf1)
       print "x = ", x, "y = ", y, "v = ", v, "ve = ", ve

       EXPR_TO_TEST1 = "(exp(x)+log(x)+log10(x)+pow(x,y)+sqrt(x)+erfc(x)+erf(x)+acos(x)+asin(x)+atan(x)+atan2(x,z)+cos(x)+cosh(x)+sin(x)+sinh(x)+tan(x)+tanh(x)+abs(y)+fabs(y))"
       EXPR_TO_TEST2 = "(x/y*z-t+(4*x)-(1.23e-3/z))"
       EXPR_TO_TEST3 = "(4 % 2)"
       EXPR_TO_TEST4 = "(-z)"
       EXPR_TO_TEST5 = "(exp(E))"
       EXPR_TO_TEST6 = "(PI)"

       def DO_SF_STKUNIT_UNIT_TEST(expr,x,y,z,t):
         
         sf = StringFunction(expr)
         v_loc = eval_func(x,y,z,t,sf)
         if isinf(v_loc):       #this is kind of wierd but Python doesn't handle infinite values like C++ and otherwise generates OverflowError
           ve_loc = v_loc
         else:   
           ve_loc = eval(expr)
         print "ve_loc = ", ve_loc, " == v_loc = ", v_loc
         self.assertEqual(ve_loc, v_loc)

       DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST1,x,y,z,t)
       DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST2,x,y,z,t)
       DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST3,x,y,z,t)
       DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST4,x,y,z,t)
       DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST5,x,y,z,t)
       DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST6,x,y,z,t)

       for xyzt in self.testpoints:
         x = xyzt[0]
         y = xyzt[1]
         z = xyzt[2]
         t = xyzt[3]
     
         DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST1,x,y,z,t)
         DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST2,x,y,z,t)
         DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST3,x,y,z,t)
         DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST4,x,y,z,t)
         DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST5,x,y,z,t)
         DO_SF_STKUNIT_UNIT_TEST(EXPR_TO_TEST6,x,y,z,t)

    def test_stringFunction_timing(self):
       numIt = 1024
       EXPR_TO_TEST1A = "(exp(x)+log(x)+log10(x)+pow(x,y)+sqrt(x)+erfc(x)+erf(x)+acos(x)+asin(x))"
       EXPR_TO_TEST1B = "(atan(x)+atan2(x,z)+cos(x)+cosh(x)+sin(x)+sinh(x)+tan(x)+tanh(x)+abs(y)+fabs(y))"
       EXPR_TO_TEST2 = "(x/y*z-t+(4*x)-(1.23e-3/z))"
       EXPR_TO_TEST8 = "(sin(x+y))"
       
       def DO_SF_TIMING_TEST_CPP(expr, numIt, x, y, z, t):
         val = 0.0
         for it in range(numIt):
           try:
             val = val + eval(expr)
           except:
             pass

       def DO_SF_TIMING_TEST_STRING(expr, numIt, x, y, z, t):
         sf = StringFunction(expr)
         val = 0.0
         for it in range(numIt):
           try:
             val = val + eval_func(x,y,z,t,sf)
           except:
             pass
  

       def TIME_IT1(expr, numIt, x, y, z, t):
         t_cpp = time.time()
         DO_SF_TIMING_TEST_CPP(expr, numIt, x, y, z, t)
         t_cpp = time.time() - t_cpp
         t_string = time.time()
         DO_SF_TIMING_TEST_STRING(expr, numIt, x, y, z, t)
         t_string = time.time() - t_string
         ratio = t_string/t_cpp
         time_it = [expr, t_cpp, t_string, ratio]
         return time_it

       print ""
       i = 0
       for xyzt in self.testpoints:
         x = xyzt[0]
         y = xyzt[1]
         z = xyzt[2]
         t = xyzt[3]
         
         print "\n Timings for ", numIt, " iterations for point #", i, " x,y,z,t= ", x, y, z, t, " \n"       
         headers = ["expression", "cpp time", "string time", "ratio"]
          
         t1 = TIME_IT1(EXPR_TO_TEST1A, numIt, x, y, z, t)    
         t2 = TIME_IT1(EXPR_TO_TEST1B, numIt, x, y, z, t)    
         t3 = TIME_IT1(EXPR_TO_TEST2, numIt, x, y, z, t)    
         t4 = TIME_IT1(EXPR_TO_TEST8, numIt, x, y, z, t)    

         table = [headers, t1, t2, t3, t4]
         out = sys.stdout
         print_table.print_table(out, table)       
         i = i + 1

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(StringFunctionUnitTests)
    unittest.TextTestRunner(verbosity=2).run(suite)
