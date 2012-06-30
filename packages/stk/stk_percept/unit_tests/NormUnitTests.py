import sys

sys.path.insert(0,"../build/build.dir/packages/PyTrilinos/src/stk/PyPercept")

from mpi4py import MPI
import unittest
from math import *
from random import *
from numpy import *
from PerceptMesh import *


class LocalFixture:
  
  def __init__(self, num_xyz, num_y, num_z):
    
    self.eMesh = PerceptMesh()
    self.num_x = num_xyz
    self.num_y = num_y
    self.num_z = num_z
    config_mesh = str(self.num_x) + "x" + str(self.num_y) + "x" + str(self.num_z) + "|bbox:-0.5,-0.5,-0.5,0.5,0.5,0.5"

    self.eMesh.new_mesh(GMeshSpec(config_mesh))
    self.eMesh.commit()
    self.metaData = self.eMesh.get_fem_meta_data()
    self.bulkData = self.eMesh.get_bulk_data()
    self.coords_field = self.metaData.get_field("coordinates")
    self.sfx = StringFunction("x", "sfx", Dimensions(3), Dimensions(1))
    self.sfx_res = ConstantFunction(0.0, "sfx_res")

def rotationMatrix(axis, angle):
     """
     Create a rotation matrix corresponding to the rotation around a general
     axis by a specified angle.

     R = dd^T + cos(a) (I - dd^T) + sin(a) skew(d)

     Parameters:

         angle : float a
         direction : array d
     """
     if axis == 0:
       direction = array([1,0,0])
     elif axis == 1:
       direction = array([0,1,0])
     elif axis == 2:
       direction = array([0,0,1])
     else:
       print "ERROR: could not determine axis"
     d = array(direction, dtype=float64)
     d /= linalg.norm(d)

     iden = eye(3, dtype=float64)
     ddt = outer(d, d)
     skew = array([[    0,  d[2],  -d[1]],
                      [-d[2],     0,  d[0]],
                      [d[1], -d[0],    0]], dtype=float64)

     mtx = ddt + cos(angle) * (iden - ddt) + sin(angle) * skew
     return mtx 

def scalingMatrix(axis, scale):
    
    sm = eye(3, dtype=float64)
    sm[axis][axis] = scale
    return sm

class MeshTransformer(GenericFunction):
     
    def __init__(self, m):
       self.rotMat = m
 
    def transform(self, domain, time_value_optional):
       x = domain[0]
       y = domain[1]
       z = domain[2]
       v = array([x, y, z])

       v = self.rotMat * v
       codomain[0] = v[0]
       codomain[1] = v[1]
       codomain[2] = v[2]

def random01():
   rnd = random.randint(0, 32678)
   rnd = (rnd+1.0)/2.0
   return rnd

class NormUnitTests(unittest.TestCase):
    
    def test_norm_volume(self):
       fix = LocalFixture(3,3,12)
       metaData = fix.metaData
       bulkData = fix.bulkData
       eMesh = fix.eMesh
       coords_field = fix.coords_field

       ff_coords = FieldFunction("ff_coords", coords_field, bulkData, Dimensions(3), Dimensions(3), FieldFunction.SIMPLE_SEARCH)
       
       identity = ConstantFunction(1.0, "identity") 
       sqrt_volume = ConstantFunction(0.0, "sqrt_volume")

       l2Norm = L2Norm(bulkData) 

       result = l2Norm.evaluate(identity) 
       #result = eval_norm(bulkData, identity, 2)
       
       self.assertAlmostEqual(1.0, result)

       print "TEST.norm.volume: writing gmesh_hex8_original_out.e .."
       eMesh.save_as("./gmesh_hex8_original_out.e")
       print "TEST.norm.volume: writing gmesh_hex8_original_out.e done"

       rmx = rotationMatrix(0, 30)
       rmy = rotationMatrix(1, -45)
       rmz = rotationMatrix(2, 30)   
       rm = rmy*rmz
       rm = rmx*rm

       print "TEST.norm.volume: writing gmesh_hex8_rotated_out.e .."
       eMesh.save_as("./gmesh_hex8_rotated_out.e")
       print "TEST.norm.volume: writing gmesh_hex8_rotated_out.e done"

       scx = pi 
       scy = e
       scz = sqrt(3.0)      
       sc = scx*scy*scz
       smx = scalingMatrix(0, scx)
       smy = scalingMatrix(1, scy)
       smz = scalingMatrix(2, scz)
       sm = smy*smz
       sm = smx*sm
       meshScale = MeshTransformer(sm)
       #eMesh.nodalOpLoop(meshScale, coords_field)
       
    def test_norm_string_function(self):
      fix = LocalFixture(4,1,1)
      metaData = fix.metaData
      bulkData = fix.bulkData
      eMesh = fix.eMesh
      coords_field = fix.coords_field
      sfx = fix.sfx
      sfx_res = fix.sfx_res

      l2Norm = L2Norm(bulkData)
      result = l2Norm.evaluate(sfx)

      sfx_expect = sqrt((0.25/3.0))
      self.assertAlmostEqual(sfx_expect, result)     

      l1Norm = L1Norm(bulkData)
      result = l1Norm.evaluate(sfx)
      sfx_expect = 0.25
      self.assertAlmostEqual(sfx_expect, result)

      sfxyz = StringFunction("x*y*z", "sfxyz", Dimensions(3), Dimensions(1))
      result = l2Norm.evaluate(sfxyz)
      sfx_expect = 0.0240562612162344
      self.assertAlmostEqual(sfx_expect, result)

      rmz = rotationMatrix(2, 30)
      rm = rmz
      meshRotate = MeshTransformer(rm)
      #eMesh.nodalOpLoop(meshRotate, coords_field) #FIXME

      result = l2Norm.evaluate(sfxyz)
      sfx_expect = 0.0178406008037016
      #if fabs(result-sfx_expect) > 0.01*sfx_expect:
         #self.assertAlmostEqual(sfx_expect, result)
         #self.assertTrue(False)

    def string_function_turbo_verify_correctness(self):
       fix = LocalFixture(4,1,1)
       metaData = fix.metaData
       bulkData = fix.bulkData
       eMesh = fix.eMesh
       coords_field = fix.coords_field
       sfx = fix.sfx
       sfx_res = fix.sfx_res
 
       ff_coords = FieldFunction("ff_coords", coords_field, bulkData, Dimensions(3), Dimensions(3), FieldFunction.SIMPLE_SEARCH)
       ff_coords.add_alias("mc")

       sfx_mc = StringFunction("mc[0]", "sfx_mc", Dimensions(3), Dimensions(1))
       sfx_mc1 = StringFunction("mc[0]", "sfx_mc1", Dimensions(3), Dimensions(1))
     
       x = -0.49+.98*random01()
       y = -0.49+.98*random01()
       z = -0.49+.98*random01()

       self.assertAlmostEqual(eval_func(x,y,z,0.0,sfx), eval_func(x,y,z,0.0,sfx_mc))
       self.assertAlmostEqual(eval_func(0.34,0,0,0.0,sfx), eval_func(0.34,0,0,0.0,sfx_mc))

       sfx_res_turbo(0.0, "sfx_res_turbo")

    





    def test_norm_field_function(self):
      fix = LocalFixture(4,1,1)
      metaData = fix.metaData
      bulkData = fix.bulkData
      coords_field = fix.coords_field

      l2Norm = L2Norm(bulkData, metaData.universal_part())

      ff_coords = FieldFunction("ff_coords", coords_field, bulkData, Dimensions(3), Dimensions(3), FieldFunction.SIMPLE_SEARCH)
      
      vals = [3,0.0]
      #sfx_res_vec = ConstantFunctionVec(vals, "sfx_res_vec")
      result = l2Norm.evaluate(ff_coords)
      sfx_expect = sqrt(.25/3.0)
      self.assertAlmostEqual(sfx_expect, result)

      l1Norm = L1Norm(bulkData, metaData.universal_part())
      result_1 = l1Norm.evaluate(ff_coords)
      sfx_expect_1 = .25      
      self.assertAlmostEqual(sfx_expect_1, result_1)




if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(NormUnitTests)
    unittest.TextTestRunner(verbosity=2).run(suite)
