import sys
sys.path.insert(0,"/scratch/srkenno/Trilinos-BUILDS/build-PyPercept-scratch-srkenno-code/packages/PyTrilinos/src/stk/PyPercept")
#sys.path.append("/scratch/srkenno/Trilinos-BUILDS/build11-090711/packages/PyTrilinos/src/stk/PyPercept")

from mpi4py import MPI
from PerceptMesh import *
from math import *
import random
import time
import unittest
from numpy import *

class CheckCoordMag(GenericFunction):

   def __init__(self, name=""):
     self.name = name
     self.error = False
   
   def evaluate(self, domain, codomain, time_value_optional=0.0):
     x = domain(0)
     y = domain(1)
     z = domain(2)
     v = sqrt(x*x + y*y + z*z)
     if fabs(v-cmag_field_node > 1.e-6):
        print "CheckCoordMag:: ", self.name, "v= ", v, "x= ", x, "y= ", y, "z= ", z, "cmag_field_node= ", cmag_field_node
        self.assertAlmostEqual(v, cmag_field_node)
        self.error = True 

class FieldFunctionUnitTests(unittest.TestCase):

    def setUp(self):
        self.testpoints = [[0.1234,     0.5678,    0.9,    0.812],
                           [0.1234e-3,  0.5678e-5, 0.97,   0.01],
                           [0.101,      0.02,      0.1020, 0.0122],
                           [0.003,      0.89,      0.01,   0.5] ]
   
    def test_fieldFunction_multiplePoints(self):
       print "start..."
       num_x = 3
       num_y = 3
       num_z = 3
       config_mesh = str(num_x) + "x" + str(num_y) + "x" + str(num_z) + "|bbox:0,0,0,1,1,1"

       eMesh = PerceptMesh()
       eMesh.new_mesh(GMeshSpec(config_mesh))
       vectorDimension = 0
       eMesh.add_field("coords_mag_field", FEMMetaData.NODE_RANK, vectorDimension)
       eMesh.commit()
       f_coords = eMesh.get_field("coordinates")
       ff_coords = FieldFunction("ff_coords", f_coords, eMesh, Dimensions(3), Dimensions(3), FieldFunction.SIMPLE_SEARCH)
       val1 = eval_vec3(0.2, 0.3, 0.4, 0.0, ff_coords)
       print "val1= ", val1
     
       points = zeros(shape=(4,3))
       output_expect = zeros(shape=(4,3))    
       output = zeros(shape=(4,3))
       
       print "here 1"
       i = 0 
       for xyzt in self.testpoints:
          x = xyzt[0]
          y = xyzt[1]
          z = xyzt[2]
          t = xyzt[3]
          points[i][0] = x
          points[i][1] = y
          points[i][2] = z

          vec = eval_vec3(x,y,z,t,ff_coords)

          tol0 = fabs(1.e-5*x)
          tol1 = fabs(1.e-5*y)
          tol2 = fabs(1.e-5*z)

          print "vec(0) = ", vec[0], " == x = ", x
          print "vec(1) = ", vec[1], " == y = ", y
          print "vec(2) = ", vec[2], " == z = ", z
          
          self.assertAlmostEqual(x, vec[0], delta=tol0)
          self.assertAlmostEqual(y, vec[1], delta=tol1)
          self.assertAlmostEqual(z, vec[2], delta=tol2)         

          output_expect[i][0] = x
          output_expect[i][1] = y
          output_expect[i][2] = z
          i = i + 1

       print "field_op: NPTS= 4"
       ff_coords.setDomainDimensions(Dimensions(3))
       ff_coords.setCodomainDimensions(Dimensions(3))
       #output = ff_coords.evaluate(points)
       # pass in the output array to ensure result is properly dimensioned
       output = ff_coords.value(points, output)
       print "here 2, output= ", output
 
       for j in range(4): #NLM
           output_expect_j = output_expect[j][0]
           output_j = output[j][0]

           tol = 1.e-5*(fabs(output_expect_j))
           print "output[j] = ", output_j, " == output_expect[j] = ", output_expect_j  , " points[j] = ", points[j]
           self.assertAlmostEqual(output_j, output_expect_j, delta = tol) 
       print "start...done"

    def test_fieldFunction_demo_1_0_0(self):
       eMesh = PerceptMesh(3)
       eMesh.new_mesh(GMeshSpec("3x3x3|bbox:0,0,0,1,1,1"))
       eMesh.commit()
       eMesh.print_info("fieldFunction_demo_1_0_0", 2)

       f_coords = eMesh.get_field("coordinates")

       ff_coords = FieldFunction("ff_coords", f_coords, eMesh, 3, 3)
       x = 0.123
       y = 0.234
       z = 0.345
       time = 0.0

       eval_vec3_print(x, y, z, time, ff_coords)

    def test_fieldFunction_read_print(self):
       print_info = 0

       x = 3
       y = 3
       z = 3
       config_mesh = str(x) + "x" + str(y) + "x" + str(z) + "|bbox:0,0,0,1,1,1"

       eMesh = PerceptMesh()
       eMesh.new_mesh_read_only(GMeshSpec(config_mesh))

       metaData = eMesh.get_fem_meta_data()
       
       parts = metaData.get_parts()       
       nparts = len(parts)
       
       if print_info == 1:
          print "Number of parts = ", nparts
       fields = metaData.get_fields()
       nfields = len(fields)
       if print_info == 1:
          print "Number of fields = ", nfields      
          for i in range(nfields):
             field = fields[i] 
             # here we have a FieldBase* which we cannot access from Python
             # is this a problem and will we want this functionality?
                 
    def test_fieldFunction_demo_1(self):
       gms = Gmesh_STKmesh_Fixture(MPI.COMM_WORLD, "3x3x3|bbox:0,0,0,1,1,1")
       print "gms = ", gms
       print "gms= end"

       eMesh = PerceptMesh()
       eMesh.new_mesh(GMeshSpec("3x3x3|bbox:0,0,0,1,1,1"))
       eMesh.commit()

       f_coords = eMesh.get_field("coordinates")
       ff_coords = FieldFunction("ff_coords", f_coords, eMesh, 3, 3)

       x = 0.123
       y = 0.234
       z = 0.345
       time = 0.0
       eval_vec3_print(x,y,z,time,ff_coords)

    def test_fieldFunction_demo_2(self):   
       eMesh = PerceptMesh()
       eMesh.new_mesh(GMeshSpec("3x3x3|bbox:0,0,0,1,1,1"))

       vectorDimension = 0
       eMesh.add_field("coords_mag_field", FEMMetaData.NODE_RANK, vectorDimension)
       eMesh.commit()

       f_coords = eMesh.get_field("coordinates")
       coords_mag_field = eMesh.get_field("coords_mag_field")

       ff_coords = FieldFunction("ff_coords", f_coords, eMesh, 3, 3)
       eval_vec3_print(0.1,0.1,0.1,0.0,ff_coords)
       
       coords_mag_sf = StringFunction("sqrt(x*x + y*y + z*z)" , "coords_mag_sf", 3, 1)
       x = 0.123
       y = 0.234
       z = 0.345
       vv = sqrt(x*x + y*y + z*z)
       v1 = eval_func(x,y,z,0,coords_mag_sf)
       print "vv = ", vv, "== v1 = ", v1
       self.assertEqual(vv, v1)
       
       coords_mag_field_function = FieldFunction("coords_mag_field_function", coords_mag_field, eMesh, 3, 1)

       coords_mag_field_function.interpolateFrom(coords_mag_sf)

       eMesh.save_as("./cubehex8_withCoordMag_out.e")

       ff_coords.add_alias("mc")
       
       sfcm = StringFunction("sqrt(mc[0]*mc[0]+mc[1]*mc[1]+mc[2]*mc[2])", "sfcm", 3, 1)

    def test_fieldFunction_readMesh_createField_interpolateFrom(self):
       num_x = 3
       num_y = 3
       num_z = 3
       config_mesh = str(num_x) + "x" + str(num_y) + "x" + str(num_z) + "|bbox:0,0,0,1,1,1"

       eMesh = PerceptMesh()
       eMesh.new_mesh(GMeshSpec(config_mesh))
       vectorDimension = 0
       eMesh.add_field("coords_mag_field", FEMMetaData.NODE_RANK, vectorDimension)
       eMesh.commit()

       #p_rank = eMesh.get_bulk_data().parallel_rank()
       #setRank(p_rank)        
          #from Util 
       f_coords = eMesh.get_field("coordinates")

       coords_mag_field = eMesh.get_field("coords_mag_field")
       #VERIFY_OP_ON      Here the unit test does something
       
       ff_coords = FieldFunction("ff_coords", f_coords, eMesh, 3, 3, FieldFunction.SIMPLE_SEARCH)

       #here we could evaluate the function
       #eval_vec3_print(0.1,0.2,0.3,0.0,ff_coords)
       
       coords_mag_sf = StringFunction("sqrt(x*x + y*y + z*z)", "coords_mag_sf", 3, 1)
       coords_mag_field_function = FieldFunction("coords_mag_field_function", coords_mag_field, eMesh, 3, 3, FieldFunction.SIMPLE_SEARCH)
       coords_mag_field_function.interpolateFrom(coords_mag_sf)

       #The following is not doable from Python
      
       checkCoordMag = CheckCoordMag()
       #eMesh.nodalOpLoop(checkCoordMag, coords_mag_field)
       print checkCoordMag.error   

       ff_coords.add_alias("mc")
       sfcm = StringFunction("sqrt(mc[0]*mc[0]+mc[1]*mc[1]+mc[2]*mc[2])", "sfcm", Dimensions(3), Dimensions(1))
      
       tol1 = 1.e-12
      
       vv = eval_vec3(0.1, 0.2, 0.3, 0.0, ff_coords)
       print 
       print "0.1 == vv[0] = ", vv[0], "passed"
       print "0.2 == vv[1] = ", vv[1], "passed"
       print "0.3 == vv[2] = ", vv[2], "passed"

       self.assertAlmostEqual(.1, vv[0], delta=tol1)
       self.assertAlmostEqual(.2, vv[1], delta=tol1)
       self.assertAlmostEqual(.3, vv[2], delta=tol1)

       vv = eval_func(0.1, 0.2, 0.3, 0.0, sfcm)
       v_expect = sqrt(0.1*0.1+0.2*0.2+0.3*0.3)

       if ((vv-v_expect) < tol1):
          print "vv = ", vv, " == v_expect = ", v_expect, "passed"

       coords_mag_field_function.interpolateFrom(sfcm)

             
 
    def test_fieldFunction_point_eval_verify(self):
       num_x = 3
       num_y = 3
       num_z = 3
       config_mesh = str(num_x) + "x" + str(num_y) + "x" + str(num_z) + "|bbox:0,0,0,1,1,1"

       eMesh = PerceptMesh()
       eMesh.new_mesh(GMeshSpec(config_mesh))
       eMesh.commit()

       f_coords = eMesh.get_field("coordinates") 
       
       ff_coords = FieldFunction("ff_coords", f_coords, eMesh, Dimensions(3), Dimensions(3), FieldFunction.SIMPLE_SEARCH)

       val1 = eval_vec3_print(0.2,0.3,0.4,0.0,ff_coords)

       bulkData = eMesh.get_bulk_data()

       try:
          val10 = eval_print_vec3(1.2, 1.3, 1.4, 0.0, ff_coords)
       except:
          print "expected to catch this exception: "

       pts = array([0.2, 0.3, 0.4])
       output_pts = array([0.0, 0.0, 0.0])
       output_pts = ff_coords.value(pts, output_pts)
       
       tol = 1.e-9
         
       print "output(0) = ", pts[0], " == output_pts(0) = ", output_pts[0]
       print "output(1) = ", pts[1], " == output_pts(1) = ", output_pts[1]
       print "output(2) = ", pts[2], " == output_pts(2) = ", output_pts[2]

       self.assertAlmostEqual(pts[0], output_pts[0], delta = tol)
       self.assertAlmostEqual(pts[1], output_pts[1], delta = tol)
       self.assertAlmostEqual(pts[2], output_pts[2], delta = tol)

    def test_fieldFunction_point_eval_timing(self):
       num_x = 3
       num_y = 3
       num_z = 3
       config_mesh = str(num_x) + "x" + str(num_y) + "x" + str(num_z) + "|bbox:0,0,0,1,1,1"

       eMesh = PerceptMesh()
       eMesh.new_mesh(GMeshSpec(config_mesh))
       eMesh.commit()

       #FIXME
       #p_size = eMesh.get_bulk_data->parallel_size()

       f_coords = eMesh.get_field("coordinates")

       for iSearchType in range(2):
          if iSearchType == 0:
             search_type = FieldFunction.SIMPLE_SEARCH
             search_type_name = "SIMPLE_SEARCH"
          else:
             search_type = FieldFunction.STK_SEARCH
             search_type_name = "STK_SEARCH"
          ff_coords = FieldFunction("ff_coords", f_coords, eMesh, Dimensions(3), Dimensions(3), search_type)
         
          t1st = time.time()
          val1 = eval_vec3(0.2,0.3,0.4,0.0,ff_coords)
          val1 = eval_vec3(0.2,0.3,0.4,0.0,ff_coords) #evaluated twice???
          t1st = time.time() - t1st

          numIter = 10000
          random.seed(12345)      
          total_time = time.time()
          max_rand = 32767

          for iter in range(numIter):
             num0 = random.randint(1, max_rand)*1.0
             num1 = random.randint(1, max_rand)*1.0
             num2 = random.randint(1, max_rand)*1.0
             pts = array([(num0/max_rand), (num1/max_rand), (num2/max_rand)])
             output_pts = array([0.0,0.0,0.0])
             output_pts = ff_coords.value(pts, output_pts, 0.0)

          total_time = time.time() - total_time

          print "TEST::function::fieldFunction_point_eval_timing: "
          print " for search_type= ", search_type_name
          print "    time for 1st eval= ", t1st
          print "    for ", numIter, "iterations, evaluating field(x,y,z) time = ", total_time
          print "    average per point lookup and eval time = ", (total_time/numIter)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(FieldFunctionUnitTests)
    unittest.TextTestRunner(verbosity=2).run(suite)

