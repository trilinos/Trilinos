import sys

sys.path.insert(0,"../build/build.dir/packages/PyTrilinos/src/stk/PyPercept")

from mpi4py import MPI
from PerceptMesh import *
import unittest
from numpy import *
import subprocess

class UseCases(unittest.TestCase):

   def test_use_case_1(self):

    pMesh = PerceptMesh()
    pMesh.new_mesh(GMeshSpec("3x3x3|bbox:0,0,0,2,2,2"))
    field = pMesh.add_field("coordinates", 1)
    pMesh.commit()

    input_array = array([1.0, 0.5, 0.5])
    input_array_2 = array([1.0, 1.5, 1.5])

    ff = FieldFunction("ff", field, pMesh, 3, 3)
    ff.add_alias("myalias")
    ff_output = ff.evaluate(input_array)

    f2 = FieldFunction("f2", field, pMesh, 3, 3)
    f2_output = f2.evaluate(input_array_2)    

    sf = StringFunction("x+y+z", "myname", 3, 1)
    sf_output = sf.evaluate(input_array)

    sf_diff = StringFunction("ff-f2", "myname") 

    norm = L1Norm(pMesh.get_bulk_data())
    value = norm.evaluate(ff)
    diffnorm = norm.evaluate(sf_diff) 


    #Now use a helper function to evaluate the norm
    #eval_norm(bulkData, Function, power)

    #value1 = eval_norm(pMesh.get_bulk_data(), ff, 1)
    #diffnorm1 = eval_norm(pMesh.get_bulk_data(), sf_diff, 1)

    #self.assertEqual(value, value1)
    #self.assertEqual(diffnorm, diffnorm1)

   def test_use_case_2(self):
     
     pMesh = PerceptMesh()
     pMesh.open("exodus_files/tet-mesh.e")

     uniform_refiner = Refiner(pMesh, TET4_TET4_8)
     pMesh.commit()

     i = 0
     while i < 3:
       uniform_refiner.doBreak()
       i = i + 1

     pMesh.save_as("tet-mesh-refined-3-times.e")

   def test_use_case_3(self):
    try:
     subprocess.call("sierra aria -i tet_mesh.e -o result_0.e") 
 
     pMesh = PerceptMesh()
     pMesh.open("tet-mesh.e")     
     uniform_refiner = Refiner(pMesh, TET4_TET4_8)
     pMesh.commit()

     uniform_refiner.doBreak()
     pMesh.save_as("tet-mesh_refined.e")

     subprocess.call("sierra aria -i tet_mesh_refined.e -o result_1.e")
 
     pMesh_0 = PerceptMesh()     
     pMesh_1 = PerceptMesh()     
     pMesh_0.open_read_only("result_0.e")
     pMesh_1.open_read_only("result_1.e")

     ff_0 = Field_Function(pMesh_0)
     ff_1 = Field_Function(pMesh_1)
     diff = StringFunction("ff_0 - ff_1")

     #diffnorm = eval_norm(pMesh.get_bulk_data, diff, 2)
     #print "diffnorm = ", diffnorm
    except:
     print "Sierra not found."

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(UseCases)
    unittest.TextTestRunner(verbosity=2).run(suite)
