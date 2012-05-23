import sys
#sys.path.append("/scratch/srkenno/Trilinos-BUILDS/build11-090711/packages/PyTrilinos/src/stk/PyPercept")
sys.path.insert(0,"/scratch/srkenno/Trilinos-BUILDS/build-PyPercept-scratch-srkenno-code/packages/PyTrilinos/src/stk/PyPercept")

from mpi4py import MPI
from PerceptMesh import *
import unittest
from numpy import *

input_files_loc = "./input_files_"

def fixture_setup_0():
  eMesh = PerceptMesh()
  p_size = eMesh.get_parallel_size()
  gmesh_spec = "4x4x"+str(4*p_size)+"|bbox:0,0,0,1,1,1"
  eMesh.new_mesh(GMeshSpec(gmesh_spec))
  eMesh.commit()
  eMesh.save_as("hex_fixture.e")

  eMesh = PerceptMesh()
  eMesh.open("exodus_files/"+input_files_loc+"hex_fixture.e")
  scalarDimension = 0
  proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
  breaker = Refiner(eMesh, HEX8_TET4_24, proc_rank_field)
  eMesh.commit()

  breaker.doBreak()
  eMesh.save_as("tet_fixture.e")
   
def fixture_setup_1():
  pm = MPI.COMM_WORLD
  p_size = parallel_machine_size(pm)
  if p_size <= 3:
    n = 12
    nx = n
    ny = n
    fixture = QuadFixture_4(pm, nx, ny, True)
    fixture.meta_data.commit()
    fixture.generate_mesh()
    eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data)
    eMesh.print_info("quad fixture", 2)
    eMesh.save_as("quad_fixture.e")

    fixture = QuadFixture_4(pm, nx, ny, False)
    fixture.meta_data.commit()
    fixture.generate_mesh()
    eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data)
    eMesh.print_info("quad fixture no sidesets", 2)
    eMesh.save_as("quad_fixture_no_sidesets.e")

is_setup = False
def fixture_setup():
  if is_setup == True:
    pass
  else:
    fixture_setup_0()
    fixture_setup_1()
    is_setup == True
    


class UniformRefinerUnitTests(unittest.TestCase):

  def test_break_quad_to_quad_sierra_1(self):
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size <= 2:
      n = 2
      nx = 1
      ny = n
      createEdgeSets = False
      fixture = QuadFixture_4(pm,nx,ny,createEdgeSets)
      isCommited = False
      eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data, isCommited)
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      eMesh.add_field("proc_rank_edge", eMesh.edge_rank(), scalarDimension)
      breaker = Refiner(eMesh, QUAD4_QUAD4_4, proc_rank_field)
      eMesh.commit()
      fixture.generate_mesh()
      breaker.setRemoveOldElements(False)
      breaker.doBreak()
      eMesh.dump_elements_compact()

  def test_break_tri_to_tri_sierra_1_test(self):
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size <= 2:
      n = 2
      nx = n
      ny = n
      createEdgeSets = False
      fixture = QuadFixture_4(pm,nx,ny,createEdgeSets)
      isCommited = False
      eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data, isCommited)
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      eMesh.add_field("proc_rank_edge", eMesh.edge_rank(), scalarDimension)
      breaker = Refiner(eMesh, TRI3_TRI3_4, proc_rank_field)
      eMesh.commit()
      fixture.generate_mesh()
      breaker.setRemoveOldElements(False)
      breaker.doBreak()
      eMesh.dump_elements_compact()

  def test_break_quad_to_quad_sierra(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size <= 3:
      n = 12
      nx = n
      ny = n
      fixture = QuadFixture_4(pm, nx, ny, True)
      eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data, False)
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, QUAD4_QUAD4_4_SIERRA, proc_rank_field)
      eMesh.commit()
      fixture.generate_mesh()
      eMesh.print_info("quad mesh")
      breaker.doBreak()

  def test_hex8_hex8_8_1(self):
    fixture_setup()
    eMesh = PerceptMesh()
    p_size = eMesh.get_parallel_size()
    gmesh_spec = "4x4x"+str(4*p_size)+"|bbox:0,0,0,1,1,1"
    eMesh.new_mesh(GMeshSpec(gmesh_spec))
    scalarDimension = 0
    proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
    breaker = Refiner(eMesh, HEX8_HEX8_8, proc_rank_field)
    eMesh.commit()
    breaker.doBreak

  def test_wedge6_1(self):
    fixture_setup()
    eMesh = PerceptMesh()
    p_size = eMesh.get_parallel_size()
    if p_size == 1:
      wedgeFixture = WedgeFixture()
      wedgeFixture.createMesh(MPI.COMM_WORLD, 4,3,2,0,1,0,1,0,1, "swept_wedge_0.e")

  def test_beam_enrich(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size <= 1:
 
      eMesh = PerceptMesh()
      eMesh.open("exodus_files/beam.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      eMesh.save_as("beam_enrich_0.e")
      breaker = Refiner(eMesh, BEAM2_BEAM3_1, proc_rank_field)
      eMesh.commit()
      breaker.setIgnoreSideSets(True)
      breaker.doBreak()

      eMesh.save_as("beam_enrich.e")

  def test_beam_refine(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size <= 1:
      mesh = BeamFixture(pm, False)
      mesh.m_metaData.commit()
      mesh.populate()

      isCommited = True
      em1 = PerceptMesh(mesh.m_metaData, mesh.m_bulkData, isCommited)
      em1.save_as("beam_0.e")

      eMesh = PerceptMesh()
      eMesh.open("beam_0.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, BEAM2_BEAM2_2, proc_rank_field)
      eMesh.commit()
      breaker.setIgnoreSideSets(True)
      breaker.doBreak()

  def test_break_quad_to_tri_6(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size == 1 or p_size == 3:
      eMesh = PerceptMesh(2)
      eMesh.open("quad_fixture.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, QUAD4_TRI3_6, proc_rank_field)
      eMesh.commit()
      eMesh.print_info("quad mesh")
      breaker.setIgnoreSideSets(True)
      breaker.setRemoveOldElements(False)                
      breaker.doBreak()
      eMesh.save_as("square_quad4_out.e")

  def test_break_quad_to_tri_4(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size == 1 or p_size == 3:
      eMesh = PerceptMesh(2)
      eMesh.open("quad_fixture.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, QUAD4_TRI3_4, proc_rank_field)
      eMesh.commit()
      eMesh.print_info("quad mesh")
      breaker.setIgnoreSideSets(True)
      breaker.setRemoveOldElements(False)
      breaker.doBreak()
      eMesh.save_as("square_quad4_tri3_4_out.e")

  def test_break_quad_to_quad(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size == 1 or p_size == 3:
      eMesh = PerceptMesh(2)
      eMesh.open("quad_fixture_no_sidesets.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, QUAD4_QUAD4_4, proc_rank_field)
      eMesh.commit()
      eMesh.print_info("quad mesh")
      breaker.setIgnoreSideSets(True)
      breaker.doBreak()
      eMesh.save_as("square_quad4_ref_out.e")

  def test_break_quad_to_quad_sierra_unit1(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size == 1 or p_size == 3:
      eMesh = PerceptMesh(2)
      eMesh.open("quad_fixture.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, QUAD4_QUAD4_4_SIERRA, proc_rank_field)
      eMesh.commit()
      eMesh.print_info("quad mesh")
      breaker.doBreak()
      eMesh.save_as("square_quad4_sierra_ref_out.e")

  def test_break_quad_to_quad_sierra_sidesets(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size == 1 or p_size == 2:
      eMesh = PerceptMesh(2)
      eMesh.open("quad_fixture.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, QUAD4_QUAD4_4_SIERRA, proc_rank_field)
      eMesh.commit()
      eMesh.print_info("after refinement break_quad_to_quad_sierra_sidesets")
      breaker.doBreak()
      eMesh.save_as("quad_sidesets_sierra_out.e")

  def test_break_hex8_tet4_24_1(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    eMesh = PerceptMesh(3)
    gmesh_spec = "1x1x" + str(p_size) + "|bbox:0,0,0,1,1," + str(p_size)
    eMesh.new_mesh(GMeshSpec(gmesh_spec))
    scalarDimension = 0
    proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
    breaker = Refiner(eMesh, HEX8_TET4_24, proc_rank_field)
    eMesh.commit()
    eMesh.print_info()
    breaker.setRemoveOldElements(True)
    breaker.doBreak()
    eMesh.save_as("hex_tet_24_cube1x1x1.e")

  def test_hex8_tet4_6_12_1(self):
    fixture_setup()
    eMesh = PerceptMesh(3)
    p_size = eMesh.get_parallel_size() 
    gmesh_spec = "1x1x" + str(p_size) + "|bbox:0,0,0,1,1," + str(p_size)
    eMesh.new_mesh(GMeshSpec(gmesh_spec))
    scalarDimension = 0
    proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
    breaker = Refiner(eMesh, HEX8_TET4_6_12, proc_rank_field)
    eMesh.commit()
    eMesh.print_info()
    breaker.setRemoveOldElements(True)
    breaker.doBreak()
    eMesh.save_as("hex_tet_6_12_cube1x1x1.e")

  def test_hex8_tet4_6_12_2(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    eMesh = PerceptMesh(3)
    eMesh.open("hex_fixture.e")
    scalarDimension = 0
    proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
    breaker = Refiner(eMesh, HEX8_TET4_6_12, proc_rank_field)
    eMesh.commit()
    eMesh.print_info()
    breaker.doBreak()
    eMesh.save_as("hex_tet_6_12_1.e")

  def test_quad4_quad4_test_1(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size <= 3:
      n = 12
      nx = n
      ny = n
      fixture = QuadFixture_4(pm,nx,ny,True)
      fixture.meta_data.commit()
      fixture.generate_mesh()
      eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data)
      eMesh.print_info("quad fixture")
      eMesh.save_as("quad_fixture_test_1.e")

  def test_break_quad_to_quad_sierra_1_unit1(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    doGenSideSets = True
    if p_size <= 3:
       n = 12
       nx = n
       ny = n
       fixture = QuadFixture_4(pm,nx,ny,doGenSideSets)
       isCommited = False
       eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data, isCommited)
       eMesh.commit()
       fixture.generate_mesh()
       eMesh.save_as("quad_fixture_0.e")
       eMesh.close()
       
       i = 0
       while i < 2:
         print "\n\n\n ================ tmp Refine Pass = ", i
         eMesh1 = PerceptMesh(2)
         eMesh1.open("quad_fixture_" + str(i) + ".e")
         scalarDimension = 0
         proc_rank_field = eMesh1.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
         breaker = Refiner(eMesh1, QUAD4_QUAD4_4_SIERRA, proc_rank_field)
         eMesh1.commit()
         breaker.doBreak()
         eMesh1.save_as("quad_fixture_" + str(i+1) + ".e")
         eMesh1.close()
         i = i + 1

  def test_break_quad_to_quad_sierra_2(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    doGenSideSets = True
    if p_size <= 3:
       n = 12
       nx = n
       ny = n
       fixture = QuadFixture_4(pm,nx,ny,doGenSideSets)
       isCommited = False
       eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data, isCommited)
       eMesh.commit()
       fixture.generate_mesh()
       eMesh.save_as("quad_fixture_mbreak_0.e")
       eMesh.close()

       eMesh1 = PerceptMesh(2)
       eMesh1.open("quad_fixture_mbreak_0.e")
       scalarDimension = 0
       proc_rank_field = eMesh1.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
       breaker = Refiner(eMesh1, QUAD4_QUAD4_4_SIERRA, proc_rank_field)
       eMesh1.commit()
       
       i = 0
       while i < 2:
         print "\n\n\n ================ tmp Refine Pass = ", i
         breaker.doBreak()
         eMesh1.save_as("quad_fixture_mbreak_" + str(i) + ".e")
         i = i + 1

  def test_break_quad4_to_quad9(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size <= 3:
       n = 12
       nx = n
       ny = n
       doGenSideSets = True
       fixture = QuadFixture_4(pm,nx,ny,doGenSideSets)
       isCommited = False
       eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data, isCommited)
       scalarDimension = 0
       proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
       breaker = Refiner(eMesh, QUAD4_QUAD9_1, proc_rank_field)
       eMesh.commit()
       fixture.generate_mesh()
       eMesh.save_as("quad_fixture_quad9_0.e")
       breaker.doBreak()
       eMesh.save_as("quad_fixture_quad9_1.e")

  def test_break_quad4_to_quad9(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size <= 3:
       n = 12
       nx = n
       ny = n
       doGenSideSets = True
       fixture = QuadFixture_4(pm,nx,ny,doGenSideSets)
       isCommited = False
       eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data, isCommited)
       scalarDimension = 0
       proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
       breaker = Refiner(eMesh, QUAD4_QUAD8_1, proc_rank_field)
       eMesh.commit()
       fixture.generate_mesh()
       eMesh.save_as("quad_fixture_quad8_0.e")
       breaker.doBreak()
       eMesh.save_as("quad_fixture_quad8_1.e")
       eMesh.save_as("quad_fixture_quad8_quad8_0.e")

  def test_quad8_to_quad8(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size <= 3:
      eMesh = PerceptMesh(2)
      eMesh.open("quad_fixture_quad8_quad8_0.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, QUAD8_QUAD8_4, proc_rank_field)
      eMesh.commit()
      breaker.setIgnoreSideSets(False)
      breaker.doBreak()
      eMesh.save_as("quad_fixture_quad8_quad8_1.e")

  def test_break_quad4_to_quad9_to_quad9_0(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    doGenSideSets = False
    if p_size <= 1:
      n = 12
      nx = n
      ny = n
      fixture = QuadFixture_4(pm,nx,ny,doGenSideSets)
      isCommited = False
      eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data, isCommited)
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, QUAD4_QUAD9_1, proc_rank_field)
      eMesh.commit()
      fixture.generate_mesh()
      breaker.doBreak()
      eMesh.save_as("quad_1x1x_quad9_quad9_0.e")
     
      em1 = PerceptMesh(2)
      em1.open("quad_1x1x_quad9_quad9_0.e")
      proc_rank_field = em1.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(em1, QUAD9_QUAD9_4, proc_rank_field)
      em1.commit()
      breaker.setIgnoreSideSets(True)
      breaker.doBreak()
      em1.save_as("quad_1x1x_quad9_quad9_1.e")

  def test_break_quad4_to_quad9_to_quad9(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    doGenSideSets = True
    if p_size <= 3:
      n = 12
      nx = n
      ny = n
      fixture = QuadFixture_4(pm, nx, ny, doGenSideSets)
      isCommited = False
      eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data, isCommited)
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, QUAD4_QUAD9_1, proc_rank_field)
      eMesh.commit()
      fixture.generate_mesh()
      breaker.doBreak()
      eMesh.save_as("quad_fixture_quad9_quad9_0.e")

      em1 = PerceptMesh(2)
      em1.open("quad_fixture_quad9_quad9_0.e")
      scalarDimension = 0
      proc_rank_field = em1.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(em1, QUAD9_QUAD9_4, proc_rank_field)
      em1.commit()
      breaker.doBreak()
      em1.save_as("quad_fixture_quad9_quad9_1.e")

  def test_break_tri_to_tri_sierra_0(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size <= 3:
      n = 12
      nx = n
      ny = n
      fixture = QuadFixture_3(pm, nx, ny, True)
      eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data, False)
      eMesh.commit()
      fixture.generate_mesh()
      eMesh.print_info("tri mesh")
      eMesh.save_as("quad_fixture_tri3.e")      

  def test_break_tri_to_tri_sierra_1(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size <= 3:
      n = 12
      nx = n
      ny = n
      createEdgeSets = True
      fixture = QuadFixture_3(pm, nx, ny, createEdgeSets)
      isCommited = False
      eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data, isCommited)
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      eMesh.add_field("proc_rank_ede", eMesh.edge_rank(), scalarDimension)
      breaker = Refiner(eMesh, TRI3_TRI3_4, proc_rank_field)
      eMesh.commit()
      fixture.generate_mesh()
      eMesh.print_info("tri mesh")
      eMesh.save_as("quad_fixture_tri3_0.e")
      breaker.doBreak()
      eMesh.print_info("tri mesh refined")
      eMesh.save_as("quad_fixture_tri3_1.e")

  def test_break_tri3_to_tri6_sierra(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)

    print p_size, "++++++++++++++++++++"
    for i in range(100):
      print "-------------------"
    if p_size <= 3:
      n = 12
      nx = n
      ny = n
      createEdgeSets = True
      fixture = QuadFixture_3(pm, nx, ny, createEdgeSets)
      isCommited = False
      eMesh = PerceptMesh(fixture.meta_data, fixture.bulk_data, isCommited)
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      eMesh.add_field("proc_rank_ede", eMesh.edge_rank(), scalarDimension)
      breaker = Refiner(eMesh, TRI3_TRI6_1, proc_rank_field)
      eMesh.commit()
      fixture.generate_mesh()
      eMesh.print_info("tri mesh tri6")
      eMesh.save_as("quad_fixture_tri3_tri6_0.e")
      breaker.doBreak()
      eMesh.print_info("tri mesh enriched")
      eMesh.save_as("quad_fixture_tri6_tri6_0.e")
      eMesh.save_as("quad_fixture_tri3_tri6_1.e")
 
  def test_break_tri3_to_tri6_to_tri6_sierra(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size <= 3:
      eMesh = PerceptMesh(2)
      eMesh.open("quad_fixture_tri6_tri6_0.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      eMesh.add_field("proc_rank_ede", eMesh.edge_rank(), scalarDimension)
      breaker = Refiner(eMesh, TRI6_TRI6_4, proc_rank_field)
      eMesh.commit()  
      eMesh.print_info("tri mesh tri6")
      eMesh.save_as("quad_fixture_tri6_tri6_0.e")
      breaker.doBreak()
      eMesh.print_info("tri mesh refined")
      eMesh.save_as("quad_fixture_tri6_tri6_1.e")

  def test_break_tet4_tet4_1(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size == 1 or p_size == 3:
      eMesh = PerceptMesh(3)
      eMesh.open_read_only("tet_fixture.e")
      eMesh.save_as("tet_from_hex_fixture_0.e")
    if p_size == 1 or p_size == 3:
      eMesh = PerceptMesh(3)
      eMesh.open("tet_from_hex_fixture_0.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, TET4_TET4_8, proc_rank_field)
      eMesh.commit()
      eMesh.print_info("tet mesh")
      breaker.doBreak()
      eMesh.save_as("tet4_refined_1.e")
      breaker.doBreak()
      eMesh.save_as("tet4_refined_2.e")
      
  def test_break_tet4_tet10_1(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size == 1 or p_size == 3:
      eMesh = PerceptMesh(3)
      eMesh.open_read_only("tet_fixture.e")
      eMesh.save_as("tet_from_hex_fixture_0.e")
    if p_size == 1 or p_size == 3:
      eMesh = PerceptMesh(3)
      eMesh.open("tet_from_hex_fixture_0.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, TET4_TET10_1, proc_rank_field)
      eMesh.commit()
      eMesh.print_info("tet mesh")
      breaker.doBreak()
      eMesh.save_as("tet10_1.e")

  def test_break_tet4_tet10_tet10_1(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size == 1 or p_size == 3:
      eMesh = PerceptMesh(3)
      eMesh.open("tet_from_hex_fixture_0.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, TET4_TET10_1, proc_rank_field)
      eMesh.commit()
      eMesh.print_info("tet mesh")
      breaker.doBreak()
      eMesh.save_as("tet10_1.e")
      eMesh.print_info("tet10_1")
    if p_size == 1 or p_size == 3:
      eMesh = PerceptMesh(3)
      eMesh.open("tet10_1.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, TET10_TET10_8, proc_rank_field)
      eMesh.commit()
      breaker.doBreak()
      eMesh.save_as("tet10_tet10_1.e")
   
  def test_hex8_hex8_8_1_unit1(self):
    fixture_setup()
    eMesh = PerceptMesh(3)
    p_size = eMesh.get_parallel_size()
    gmesh_spec = "1x1x" + str(p_size) + "|bbox:0,0,0,1,1," + str(p_size)
    eMesh.new_mesh(GMeshSpec(gmesh_spec))
    scalarDimension = 0
    proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
    breaker = Refiner(eMesh, HEX8_HEX8_8, proc_rank_field) 
    eMesh.commit()
    eMesh.print_info()
    eMesh.save_as("hex_hex_cube1x1x" + str(p_size) + "-orig.e")
    breaker.setRemoveOldElements(True)
    breaker.doBreak()
    eMesh.save_as( "hex_hex_cube1x1x" + str(p_size)+".e")

  def test_hex8_hex8_8_2_unit1(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size == 1 or p_size == 3:
      eMesh = PerceptMesh(3)
      eMesh.open("hex_fixture.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, HEX8_HEX8_8, proc_rank_field)
      eMesh.commit()
      eMesh.print_info()
      eMesh.save_as("hex8_0.e")
      breaker.doBreak()
      eMesh.save_as("hex8_1.e")
      breaker.doBreak()
      eMesh.save_as("hex8_2.e")
      
  def test_hex8_hex27_1_2(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size == 1 or p_size == 3:
      eMesh = PerceptMesh(3)
      eMesh.open("hex_fixture.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, HEX8_HEX27_1, proc_rank_field)
      eMesh.commit()
      eMesh.print_info()
      eMesh.save_as("hex27_0.e")
      breaker.doBreak()
      eMesh.save_as("hex27_1.e")

  def test_hex08_hex20_1_1(self):
    fixture_setup()
    eMesh = PerceptMesh(3)
    p_size = eMesh.get_parallel_size()
    gmesh_spec = "1x1x" + str(p_size) + "|bbox:0,0,0,1,1," + str(p_size)
    eMesh.new_mesh(GMeshSpec(gmesh_spec))
    scalarDimension = 0
    proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
    breaker = Refiner(eMesh, HEX8_HEX20_1, proc_rank_field)
    eMesh.commit()
    eMesh.print_info()
    eMesh.save_as("hex8_hex_20_cube1x1x"+str(p_size) + "-orig.e")
    breaker.setRemoveOldElements(True)
    breaker.doBreak()
    eMesh.save_as("hex8_hex20_cube1x1x"+str(p_size) + ".e")
    eMesh.save_as("hex20_hex20_cube1x1x"+str(p_size) + "_0.e")
    
  def test_hex08_hex20_1_2(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size == 1 or p_size == 3:
      eMesh = PerceptMesh(3)
      eMesh.open("hex_fixture.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, HEX8_HEX20_1, proc_rank_field)
      eMesh.commit()
      eMesh.print_info()
      eMesh.save_as("hex20_0.e")
      breaker.doBreak()
      eMesh.save_as("hex20_1.e")
      eMesh.save_as("hex20_hex20_0.e")

  def test_hex20_hex20_1(self):
    fixture_setup()
    eMesh = PerceptMesh(3)
    p_size = eMesh.get_parallel_size()
    if p_size <= 3:
      eMesh.open("hex20_hex20_cube1x1x"+str(p_size) + "_0.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, HEX20_HEX20_8, proc_rank_field)
      eMesh.commit()
      breaker.setRemoveOldElements(True)
      breaker.doBreak()
      eMesh.save_as("hex20_hex20_cube1x1x" + str(p_size) + "_1.e")

  def test_hex20_hex20_1_2(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size == 1 or p_size == 3:
      eMesh = PerceptMesh(3)
      eMesh.open("hex20_hex20_0.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, HEX20_HEX20_8, proc_rank_field)
      eMesh.commit()
      eMesh.save_as("hex20_hex20_0.e")
      breaker.doBreak()
      eMesh.save_as("hex20_hex20_1.e")

  def test_hex27_hex27_0(self):
    fixture_setup()
    scalarDimension = 0
    eMesh = PerceptMesh()
    p_size = eMesh.get_parallel_size()
    gmesh_spec = "1x1x" + str(p_size) + "|bbox:0,0,0,1,1," + str(p_size)
    eMesh.new_mesh(GMeshSpec(gmesh_spec))
    
    proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
    breaker = Refiner(eMesh, HEX8_HEX27_1, proc_rank_field)
    eMesh.commit()
    eMesh.print_info()
    eMesh.save_as("hex27_hex27_cube1x1x" + str(p_size) + "-orig.e")
    breaker.setRemoveOldElements(True)
    breaker.doBreak()
    eMesh.save_as("hex27_hex27_cube1x1x" + str(p_size) + "_0.e")
    
    em1 = PerceptMesh(3)
    p_size = em1.get_parallel_size()
    em1.open("hex27_hex27_cube1x1x" + str(p_size) + "_0.e")
    proc_rank_field = em1.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
    breaker = Refiner(em1, HEX27_HEX27_8, proc_rank_field)
    em1.commit()
    breaker.setIgnoreSideSets(True)
    breaker.setRemoveOldElements(True)
    breaker.doBreak()
    em1.save_as("hex27_hex27_cube1x1x" + str(p_size) + "_1.e")

  def test_hex27_hex27_1(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size == 1 or p_size == 3:
      eMesh = PerceptMesh(3)
      eMesh.open("hex_fixture.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, HEX8_HEX27_1, proc_rank_field)
      eMesh.commit()
      eMesh.print_info()
      eMesh.save_as("hex8_hex27_0.e")
      breaker.doBreak()
      eMesh.save_as("hex8_27_1.e")
    if p_size == 1 or p_size == 3:
      eMesh = PerceptMesh(3)
      eMesh.open("hex8_27_1.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, HEX27_HEX27_8, proc_rank_field)
      eMesh.commit()
      breaker.setRemoveOldElements(True)
      breaker.doBreak()
      eMesh.save_as("hex8_hex27_hex27_1.e")

  def test_wedge6_2(self):
    fixture_setup()
    eMesh = PerceptMesh(3)
    p_size = eMesh.get_parallel_size()
    if p_size == 1:
      wedgeFixture = WedgeFixture()
      wedgeFixture.createMesh(MPI.COMM_WORLD, 4,3,2,0,1,0,1,0,1, "swept_wedge_0.e")
      eMesh.open("swept_wedge_0.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, WEDGE6_WEDGE6_8, proc_rank_field)
      eMesh.commit()
      breaker.doBreak()
      eMesh.save_as("swept-wedge_1.e")
     
  def test_wedge6_enrich_1(self):
    fixture_setup()
    eMesh = PerceptMesh(3)
    p_size = eMesh.get_parallel_size()
    if p_size == 1:
      wedgeFixture = WedgeFixture()
      wedgeFixture.createMesh(MPI.COMM_WORLD, 4,3,2,0,1,0,1,0,1, "swept_wedge_enrich_0.e")
      eMesh.open("swept_wedge_enrich_0.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, WEDGE6_WEDGE15_1, proc_rank_field)
      eMesh.commit()
      breaker.doBreak()
      eMesh.save_as("swept-wedge_enrich_1.e")
      eMesh.save_as("swept-wedge_enrich_refine_0.e")

  def test_wedge6_enrich_refine(self):
    fixture_setup()
    p_size = parallel_machine_size(MPI.COMM_WORLD)
    if p_size == 1:
      eMesh = PerceptMesh(3)
      wedgeFixture = WedgeFixture()
      wedgeFixture.createMesh(MPI.COMM_WORLD, 4,2,2,0,1,0,1,0,1, "tmp-swept-wedge_enrich_0.e")
      eMesh.open("tmp-swept-wedge_enrich_0.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, WEDGE6_WEDGE15_1, proc_rank_field)
      eMesh.commit()
      breaker.doBreak()
      eMesh.save_as("swept-wedge_2_enrich_refine_0.e")
    if p_size == 1:
      eMesh = PerceptMesh(3)
      eMesh.open("swept-wedge_2_enrich_refine_0.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, WEDGE15_WEDGE15_8, proc_rank_field)
      eMesh.commit()
      breaker.setIgnoreSideSets(True)
      breaker.doBreak()
      eMesh.save_as("swept-wedge_2_enrich_refine_1.e")

  def test_heterogeneous_mesh(self):
    fixture_setup()
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(MPI.COMM_WORLD)
    if p_size <= 1:
      mesh = HeterogeneousFixture(MPI.COMM_WORLD, False)
      #put_io_part_attribute(mesh.m_block_hex)
      #put_io_part_attribute(mesh.m_block_wedge)
      #put_io_part_attribute(mesh.m_block_tet)
      mesh.m_metaData.commit()
      mesh.populate()
      isCommited = True
      em1 = PerceptMesh(mesh.m_metaData, mesh.m_bulkData, isCommited)
      em1.save_as("heterogeneous_0.e")
      em1.close()
      
      eMesh = PerceptMesh(3)
      eMesh.open("heterogeneous_0.e")
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      #INCOMPLETE

  def test_wedge6_wedge18_enrich(self):
    pm = MPI.COMM_WORLD
    p_size = parallel_machine_size(pm)
    if p_size == 1:
      wedgeFixture = WedgeFixture()
      bulk = wedgeFixture.createMesh(MPI.COMM_WORLD, 4,3,2,0,1,0,1,0,1,"")
      eMesh = PerceptMesh(wedgeFixture.getMetaData(), bulk, False)
      scalarDimension = 0
      proc_rank_field = eMesh.add_field("proc_rank", eMesh.element_rank(), scalarDimension)
      breaker = Refiner(eMesh, WEDGE6_WEDGE18_1, proc_rank_field)
      eMesh.commit()
      wedgeFixture.createBulkAfterMetaCommit(MPI.COMM_WORLD)
      breaker.doBreak()
      

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(UniformRefinerUnitTests)
    unittest.TextTestRunner(verbosity=2).run(suite)

