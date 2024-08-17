#include "FaceCreatorFixture.hpp"
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/QuadFixture.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>

namespace
{

void convert_quad_fixture_to_my_bulk_data_flavor(unsigned numX, unsigned numY, stk::mesh::BulkData& bulkData)
{
  stk::mesh::fixtures::QuadFixture fixture(bulkData.parallel(), numX, numY, false);

  stk::mesh::Field<double> &coordField = fixture.m_meta.declare_field<double>(stk::topology::NODE_RANK, "model_coordinates");
  stk::mesh::put_field_on_mesh(coordField, fixture.m_meta.universal_part(), fixture.m_meta.spatial_dimension(), nullptr);
  stk::mesh::Part& block_1 = fixture.m_meta.declare_part_with_topology("block_1", stk::topology::QUADRILATERAL_4_2D);
  stk::io::set_field_output_type(coordField, stk::io::FieldOutputType::VECTOR_3D);
  stk::io::put_io_part_attribute(block_1);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  std::vector<double> x;
  std::vector<double> y;

  for(unsigned j=0;j<=numY;++j)
  {
    for(unsigned i=0;i<=numX;i++)
    {
      x.push_back(i); // 0 1 2, 0 1 2, 0 1 2, ...
      y.push_back(j); // 0 0 0, 1 1 1
    }
  }

  stk::mesh::EntityVector nodes;
  stk::mesh::get_selected_entities(fixture.m_meta.universal_part(), fixture.m_bulk_data.buckets(stk::topology::NODE_RANK), nodes);
  for(stk::mesh::Entity node : nodes )
  {
    double* coords = stk::mesh::field_data(coordField, node);
    unsigned id = fixture.m_bulk_data.identifier(node);
    coords[0] = x[id-1];
    coords[1] = y[id-1];
  }

  fixture.m_bulk_data.modification_begin();
  stk::mesh::EntityVector elems;
  stk::mesh::get_selected_entities(fixture.m_meta.locally_owned_part(), fixture.m_bulk_data.buckets(stk::topology::ELEM_RANK), elems);
  for(stk::mesh::Entity element : elems)
  {
    fixture.m_bulk_data.change_entity_parts(element, stk::mesh::ConstPartVector{&block_1});
  }
  fixture.m_bulk_data.modification_end();

  std::ostringstream os;
  const std::string file_temp("testadfasdasdfas.exo");
  stk::io::write_mesh(file_temp, fixture.m_bulk_data);
  stk::io::fill_mesh(file_temp, bulkData);

  STK_ThrowRequireMsg(fixture.m_bulk_data.parallel_size()<10, "Testing assumption violated.");
  os << file_temp << "." << fixture.m_bulk_data.parallel_size() << "." << fixture.m_bulk_data.parallel_rank();
  unlink(os.str().c_str());
}

class FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester : public FaceCreatorFixture
{
protected:

  FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester() : FaceCreatorFixture(2) {}

  void setup_2x1_2d_mesh(stk::mesh::BulkData::AutomaticAuraOption aura_option)
  {
    setup_empty_mesh(aura_option);
    unsigned numX = 2, numY = 1;
    convert_quad_fixture_to_my_bulk_data_flavor(numX, numY, get_bulk());
  }

  virtual stk::mesh::EntityVector get_nodes_of_face_for_this_proc()
  {
    std::vector<unsigned> face_node_ids = { 2, 5 };
    return get_nodes_for_proc(face_node_ids);
  }

  virtual unsigned get_permuted_index(unsigned i)
  {
    std::vector<std::vector<unsigned> > index_for_proc = {
      {0, 1},
      {1, 0}
    };
    return index_for_proc[get_bulk().parallel_rank()][i];
  }
};


TEST_F(FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester, twoHexesTwoProcsCreateTwoFacesWithAura)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    setup_2x1_2d_mesh(stk::mesh::BulkData::AUTO_AURA);
    test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
  }
}

TEST_F(FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester, twoHexesTwoProcsCreateTwoFacesWithoutAura)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    setup_2x1_2d_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    test_that_one_face_exists_after_both_procs_create_face_on_proc_boundary();
  }
}

TEST_F(FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester, twoHexesTwoProcsCreateOneFaceWithAura)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    setup_2x1_2d_mesh(stk::mesh::BulkData::AUTO_AURA);
    test_that_one_face_exists_on_both_procs_after_only_one_proc_makes_face();
  }
}

TEST_F(FaceCreator2DElemGraphUsingBDElemGraphFaceSharingTester, twoHexesTwoProcsCreateOneFaceWithoutAura)
{
  if(stk::parallel_machine_size(get_comm())==2)
  {
    setup_2x1_2d_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    test_that_one_face_exists_on_both_procs_after_only_one_proc_makes_face();
  }
}

}
