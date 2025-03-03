#include <iostream>

#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/getOption.h>

#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/baseImpl/EquivalentEntityBlocks.hpp>

namespace
{

class StkToolsA : public stk::unit_test_util::MeshFixture
{

};

TEST_F(StkToolsA, WriteTextMeshDescFromExodusFile)
{
  const std::string unNamed = "mesh not specified";
  const std::string meshName = stk::unit_test_util::get_option("-i", unNamed);
  STK_ThrowRequireMsg(meshName!=unNamed, "Please specify mesh with -i option.");
  setup_mesh(meshName, stk::mesh::BulkData::NO_AUTO_AURA);

  std::ostringstream os;
  os.precision(16);

  const stk::mesh::BucketVector &buckets = get_bulk().get_buckets(stk::topology::ELEM_RANK, get_meta().locally_owned_part());

  os << "const std::string meshDesc = \"";

  std::map<unsigned, int> localToGlobalId;
  std::vector<stk::mesh::PartOrdinal> partOrdinals;

  unsigned counter = 0;
  for(size_t i=0;i<buckets.size();++i)
  {
    const stk::mesh::Bucket& bucket = *buckets[i];
    for(size_t j=0;j<bucket.size();j++)
    {
      stk::mesh::Entity element = bucket[j];
      stk::mesh::impl::get_element_block_part_ordinals(element, get_bulk(), partOrdinals);
      unsigned partOrdinal = partOrdinals[0];
      if(localToGlobalId.find(partOrdinal)==localToGlobalId.end())
      {
        localToGlobalId[partOrdinal]=counter;
        counter++;
      }
    }
  }

  for(size_t i=0;i<buckets.size();++i)
  {
    const stk::mesh::Bucket& bucket = *buckets[i];
    for(size_t j=0;j<bucket.size();j++)
    {
      stk::mesh::Entity element = bucket[j];
      os << stk::parallel_machine_rank(get_comm()) << "," << get_bulk().identifier(element) << ",HEX_8";
      unsigned num_nodes = get_bulk().num_nodes(element);
      const stk::mesh::Entity *nodes = get_bulk().begin_nodes(element);
      for(unsigned k=0;k<num_nodes;k++)
      {
        os << "," << get_bulk().identifier(nodes[k]);
      }

      stk::mesh::impl::get_element_block_part_ordinals(element, get_bulk(), partOrdinals);
      unsigned part_id = partOrdinals[0];
      unsigned local_id = localToGlobalId.at(part_id);
      os << ",block_" << local_id+1 << "\\n\\\n";
    }
  }
  os << "\";\n";

  stk::mesh::Field<double>* coordField = get_meta().get_field<double>(stk::topology::NODE_RANK, "coordinates");
  ASSERT_TRUE(coordField!=nullptr);

  stk::mesh::Selector nodeSelector = get_meta().locally_owned_part() | get_meta().globally_shared_part();
  stk::mesh::EntityVector nodes;
  stk::mesh::get_entities(get_bulk(), stk::topology::NODE_RANK, nodes);

  os << "std::vector<double> coordinates{\n";
  for(stk::mesh::Entity node : nodes)
  {
    double *coord = stk::mesh::field_data(*coordField, node);
    os << coord[0] << ", " << coord[1] << ", " << coord[2] << "," << std::endl;
  }
  os << "};\n";

  std::ofstream outfile("test.out.txt");
  outfile << os.str();
}
}
