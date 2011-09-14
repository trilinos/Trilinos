import sys
sys.path.append("/scratch/srkenno/Trilinos-BUILDS/build11-090711/packages/PyTrilinos/src/stk/PyPercept")

from mpi4py import MPI
from PerceptMesh import *
import unittest
from numpy import *


class PerceptMeshUnitTests(unittest.TestCase):

    is_setup = 0
    pm = MPI.COMM_WORLD
    input_files = "./input_files/"

    def fixture_setup_0(self):
       eMesh = PerceptMesh()

       p_size = eMesh.getParallelSize()
      
       gmesh_spec = "4x4x" + str((4*p_size)) + "|bbox:0,0,0,1,1,1"
       eMesh.newMesh(GMeshSpec(gmesh_spec))
       eMesh.commit()
       eMesh.saveAs("hex_fixture.e")

    def fixture_setup_1(self):

       p_size = self.pm.size 
       if p_size <= 2:
         n = 12
         nx = n
         ny = n
         fixture = QuadFixture_4(self.pm, nx, ny, 1) 
         fixture.meta_data.commit()
         fixture.generate_mesh()
         
         eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data)
         eMesh.printInfo("quad fixture", 2)
         eMesh.saveAs("quad_fixture.e")
       
       if p_size <= 2:
         n = 12
         nx = n
         ny = n
       
         fixture = QuadFixture_4(self.pm, nx, ny, 0)
         fixture.meta_data.commit()
         fixture.generate_mesh()

         eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data)
         eMesh.printInfo("quad fixture no sidesets", 2)
         eMesh.saveAs("quad_fixture_no_sidesets.e")

    def fixture_setup(self): 
       if self.is_setup == 1:
         pass
       else:
         self.fixture_setup_0()
         self.fixture_setup_1()
         self.is_setup = 1
          
    def test_unit_perceptMesh_wedge6_1(self):     
       eMesh = PerceptMesh()
       p_size = eMesh.getParallelSize()
       if p_size == 1:
         wedgeFixture = WedgeFixture()
         wedgeFixture.createMesh(self.pm, 4, 3, 2, 0, 1, 0, 1, 0, 1, "swept-wedge_0.e")

    def test_perceptMesh_walk_nodes(self):
       self.fixture_setup()
       p_size = self.pm.size
       p_rank = self.pm.rank
       if p_size <= 2:
         n = 12
         nx = n
         ny = n

         sidesets_on = 1
         fixture = QuadFixture_4(self.pm, nx, ny, sidesets_on)
         fixture.meta_data.commit()
         fixture.generate_mesh()

         eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data)
         eMesh.printInfo("quad fixture", 2)

         metaData = eMesh.getFEM_meta_data()
         parts = metaData.get_parts()

         nparts = len(parts)
         print "Number of parts = ", nparts

         surface_id = 2
         surface_name = "surface_" + str(surface_id)
         part = eMesh.getNonConstPart(surface_name)
         in_surface_selector = Selector(part)
         bulkData = eMesh.getBulkData()
         coordField = eMesh.getCoordinatesField()

         if eMesh.getSpatialDim() == 2:
           buckets_arg = eMesh.edge_rank()
         else:
           buckets_arg = eMesh.face_rank
         buckets = bulkData.buckets(buckets_arg)
         sum = 0.0
       
       #  for bucket in buckets: #FIXME
       #    if in_surface_selector(bucket) == 1:
       #       cell_topo_data = PerceptMesh.get_cell_topology(bucket)
       #       cell_topo(cell_topo_data)  #FIXME 
       #       num_elements_in_bucket = bucket.size()

       #       for iElement in range(num_elements_in_bucket):
       #          element = bucket[iElement] #FIXME
       #          elem_nodes = element.relations(FEMMetaData.NODE_RANK) #FIXME  
       #          num_node = elem_nodes.size() #FIXME
       #  
       #          for inode in range(num_node):
       #             node = elem_nodes[inode].entity()
   
    def test_mesh_diff(self):
      self.fixture_setup()
      p_size = parallel_machine_size(self.pm)

      if p_size <= 2:
        eMesh_0 = PerceptMesh(2)
        eMesh_0.openReadOnly("quad_fixture.e")
        eMesh_0.saveAs("quad_fixture_readwrite.e")

        eMesh_1 = PerceptMesh(2)
        eMesh_2 = PerceptMesh(2)
        eMesh_1.openReadOnly("quad_fixture_readwrite.e")
        eMesh_2.openReadOnly("quad_fixture.e")

      if p_size == 1:
        add_newlines = False
        eMesh_1.printInfo("quad fixture", 2, add_newlines)
        eMesh_2.printInfo("quad fixture", 2, add_newlines)
        #Here the unit test compares an expected output string with the output of the printInfos

      diff_msg = "diff report: "
      diff = PerceptMesh.mesh_difference(eMesh_1, eMesh_2, diff_msg, True )
      self.assertFalse(diff)      
      
      metaData_1 = eMesh_1.getFEM_meta_data()
      metaData_2 = eMesh_2.getFEM_meta_data()
      bulkData_1 = eMesh_1.getBulkData()
      bulkData_2 = eMesh_2.getBulkData()
      coordField_1 = eMesh_1.getCoordinatesField()
      coordField_2 = eMesh_2.getCoordinatesField()

      diff = PerceptMesh.mesh_difference(metaData_1, metaData_2, bulkData_1, bulkData_2, diff_msg, True)
      self.assertFalse(diff)

      buckets = bulkData_1.buckets(FEMMetaData.NODE_RANK)
    #  for bucket in buckets: #FIXME
    #    num_elements_in_bucket = bucket.size()
    #    for iEntity in range(num_elements_in_bucket):
    #      entity = bucket[iEntity]
    #      coord = field_data(coordField_1, entity)
    #      coord[0] = coord[0] + 0.01
    #  diff_msg = "diff report after mod: "
    #  diff = PerceptMesh.mesh_difference(eMesh_1, eMesh_2, diff_msg, True)
    #  self.assertTrue(diff)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(PerceptMeshUnitTests)
    unittest.TextTestRunner(verbosity=2).run(suite)
