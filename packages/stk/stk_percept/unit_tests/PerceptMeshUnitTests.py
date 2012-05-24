import sys

sys.path.insert(0,"../build/build.dir/packages/PyTrilinos/src/stk/PyPercept")

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

       p_size = eMesh.get_parallel_size()
      
       gmesh_spec = "4x4x" + str((4*p_size)) + "|bbox:0,0,0,1,1,1"
       eMesh.new_mesh(GMeshSpec(gmesh_spec))
       eMesh.commit()
       eMesh.save_as("./exodus_files/hex_fixture.e")

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
         eMesh.print_info("quad fixture", 2)
         eMesh.save_as("./exodus_files/quad_fixture.e")
       
       if p_size <= 2:
         n = 12
         nx = n
         ny = n
       
         fixture = QuadFixture_4(self.pm, nx, ny, 0)
         fixture.meta_data.commit()
         fixture.generate_mesh()

         eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data)
         eMesh.print_info("quad fixture no sidesets", 2)
         eMesh.save_as("./exodus_files/quad_fixture_no_sidesets.e")

    def fixture_setup(self): 
       if self.is_setup == 1:
         pass
       else:
         self.fixture_setup_0()
         self.fixture_setup_1()
         self.is_setup = 1
          
    def test_unit_perceptMesh_wedge6_1(self):     
       eMesh = PerceptMesh()
       p_size = eMesh.get_parallel_size()
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
         eMesh.print_info("quad fixture", 2)

         metaData = eMesh.get_fem_meta_data()
         parts = metaData.get_parts()

         nparts = len(parts)
         print "Number of parts = ", nparts

         surface_id = 2
         surface_name = "surface_" + str(surface_id)
         part = eMesh.get_part(surface_name)
         in_surface_selector = Selector(part)
         bulkData = eMesh.get_bulk_data()
         coordField = eMesh.get_coordinates_field()

         if eMesh.get_spatial_dim() == 2:
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
        eMesh_0.open_read_only("./exodus_files/quad_fixture.e")
        eMesh_0.save_as("./exodus_files/quad_fixture_readwrite.e")

        eMesh_1 = PerceptMesh(2)
        eMesh_2 = PerceptMesh(2)
        eMesh_1.open_read_only("./exodus_files/quad_fixture_readwrite.e")
        eMesh_2.open_read_only("./exodus_files/quad_fixture.e")

      if p_size == 1:
        add_newlines = False
        eMesh_1.print_info("quad fixture", 2, add_newlines)
        eMesh_2.print_info("quad fixture", 2, add_newlines)
        #Here the unit test compares an expected output string with the output of the print_infos

      diff_msg = "diff report: "
      diff = PerceptMesh.mesh_difference(eMesh_1, eMesh_2, diff_msg, True )
      self.assertFalse(diff)      
      
      #metaData_1 = eMesh_1.get_fem_meta_data()
      #metaData_2 = eMesh_2.get_fem_meta_data()
      bulkData_1 = eMesh_1.get_bulk_data()
      #bulkData_2 = eMesh_2.get_bulk_data()
      #coordField_1 = eMesh_1.get_coordinates_field()
      #coordField_2 = eMesh_2.get_coordinates_field()

      #diff = PerceptMesh.mesh_difference(metaData_1, metaData_2, bulkData_1, bulkData_2, diff_msg, True)
      #self.assertFalse(diff)

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

    # new tests of simplified high-level PerceptMesh interface
    def test_high_level_interface(self):
      self.fixture_setup()
      p_size = parallel_machine_size(self.pm)

      if p_size <= 2:
        eMesh = PerceptMesh(2)
        eMesh.open("./exodus_files/quad_fixture.e")

        vectorDimension = 0
        eMesh.add_field("coords_mag_field", FEMMetaData.NODE_RANK, vectorDimension)
        eMesh.commit()

        f_coords = eMesh.get_field("coordinates")
        coords_mag_field = eMesh.get_field("coords_mag_field")

        ff_coords = FieldFunction("ff_coords", f_coords, eMesh, 2, 2)
        #eval_vec3_print(0.1,0.1,0.1,0.0,ff_coords)

        coords_mag_sf = StringFunction("sqrt(x*x + y*y )" , "coords_mag_sf", 2, 1)
        x = 0.123
        y = 0.234
        vv = sqrt(x*x + y*y )
        v1 = eval_func2(x,y,0,coords_mag_sf)
        print "vv = ", vv, "== v1 = ", v1
        self.assertEqual(vv, v1)

        coords_mag_field_function = FieldFunction("coords_mag_field_function", coords_mag_field, eMesh, 2, 1)

        coords_mag_field_function.interpolateFrom(coords_mag_sf)

        eMesh.save_as("./exodus_files/quad_fixture_with_coords_mag.e")

        ff_coords.add_alias("mc")

        sfcm = StringFunction("sqrt(mc[0]*mc[0]+mc[1]*mc[1]+mc[2]*mc[2])", "sfcm", 3, 1)

        add_newlines = True
        eMesh.print_info("quad fixture", 2, add_newlines)

        self.assertTrue(eMesh.get_spatial_dim() == 2)
        self.assertTrue(eMesh.get_number_elements() == 12*12)
        self.assertTrue(eMesh.get_number_nodes() == 13*13)

        self.assertTrue(eMesh.get_parallel_size() == p_size)

        self.assertTrue(eMesh.get_bulk_data() != 0)
        self.assertTrue(eMesh.get_fem_meta_data() != 0)

        # // entity data setter/getters
        node = eMesh.get_node(1)
        self.assertTrue(node != 0)
        cm1 = eMesh.get_field_data(coords_mag_field, node)
        co1 = [0,0]
        co1[0] = eMesh.get_field_data(f_coords, node, 0)
        co1[1] = eMesh.get_field_data(f_coords, node, 1)
        print "cm1= ", cm1, " co1= ", co1
        eMesh.set_field_data(123.0, f_coords, node, 0)
        co1[0] = eMesh.get_field_data(f_coords, node, 0)
        print " co1= ", co1
        
        element = eMesh.get_element(1)
        self.assertTrue(element != 0)

        element1 = eMesh.get_entity(eMesh.element_rank(), 1)
        self.assertTrue(element == element1)
        
        #/// find node closest to given point
        node = eMesh.get_node(0,0)
        self.assertTrue(node != 0)

        #/// find element that contains given point
        element = eMesh.get_element(0.01, 0.01)
        self.assertTrue(element != 0)
        


        #self.assertFalse(diff)

    # 
    def test_time_dep_interface(self):
      p_size = parallel_machine_size(self.pm)

      if p_size <= 2:
        eMesh = PerceptMesh(2)
        eMesh.open("./exodus_files/time-dep.e")
        eMesh.commit()

        Tnd_field = eMesh.get_field("Tnd")

        # // entity data setter/getters
        eMesh.read_database_at_time(0.0)
        node = eMesh.get_node(2,2)
        self.assertTrue(node != 0)
        t0 = eMesh.get_field_data(Tnd_field, node)
        self.assertTrue(t0 == 0.0)
        eMesh.read_database_at_time(1.0)
        t1 = eMesh.get_field_data(Tnd_field, node)
        self.assertTrue(t1 == 11.0)
        
        print "t0= " , t0, " t1= " , t1


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(PerceptMeshUnitTests)
    unittest.TextTestRunner(verbosity=2).run(suite)
