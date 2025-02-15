#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/getOption.h>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/DestroyElements.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_util/util/string_utils.hpp>

namespace
{

class StkToolsC : public stk::unit_test_util::MeshFixture
{};

TEST_F(StkToolsC, DeleteMeshExceptSpecifiedElems)
{
  const std::string unNamed = "mesh not specified";
  const std::string meshName = stk::unit_test_util::get_option("-i", unNamed);
  STK_ThrowRequireMsg(meshName!=unNamed, "Please specify mesh with -i option.");
  setup_mesh(meshName, stk::mesh::BulkData::NO_AUTO_AURA);

  std::string invalidElemId = "-1";
  std::string inputElemIds = stk::unit_test_util::get_command_line_option("-e", invalidElemId);
  STK_ThrowRequireMsg(inputElemIds != invalidElemId, "Please specify element list with -e.");

  std::set<stk::mesh::EntityId> elemIdsToKeep;
  std::vector<std::string> elemIdSegments = stk::split_csv_string(inputElemIds);
  for (const std::string & elemIdSegment : elemIdSegments) {
    const std::string trimmedElemId = stk::trim_string(elemIdSegment);
    const int elemId = std::stoi(trimmedElemId);
    elemIdsToKeep.insert(elemId);
  }

  const stk::mesh::BucketVector &buckets = get_bulk().get_buckets(stk::topology::ELEM_RANK, get_meta().locally_owned_part());

  std::vector<size_t> entityCounts;
  stk::mesh::comm_mesh_counts(get_bulk(), entityCounts);

  std::cout << "num entities = " << entityCounts[stk::topology::ELEMENT_RANK] << std::endl;

  stk::mesh::EntityVector elementsToDestroy;
  for(size_t i=0;i<buckets.size();++i)
  {
    const stk::mesh::Bucket& bucket = *buckets[i];
    for(size_t j=0;j<bucket.size();j++)
    {
      stk::mesh::Entity element = bucket[j];
      stk::mesh::EntityId elemId = get_bulk().identifier(element);
      if (elemIdsToKeep.find(elemId) == elemIdsToKeep.end())
      {
        elementsToDestroy.push_back(element);
      }
      else
      {
        std::cout << "keeping element ID: " << elemId << std::endl;
      }
    }
  }
  destroy_elements(get_bulk(), elementsToDestroy);
  stk::io::write_mesh("modified.e", get_bulk());
}

void stk_determine_centroid(const unsigned spatial_dim, stk::mesh::Entity element,
                            const stk::mesh::FieldBase& nodal_coord, std::vector<double> & centroid)
{
  centroid.resize(spatial_dim);
  std::fill(centroid.begin(), centroid.end(), 0.0);

  const stk::mesh::BulkData& bulk_data = nodal_coord.get_mesh();
  const stk::mesh::Entity* const node_vec = bulk_data.begin_nodes(element);
  const unsigned num_nodes = bulk_data.num_nodes(element);

  for (unsigned iNode = 0; iNode < num_nodes; ++iNode) {
    stk::mesh::Entity node = node_vec[iNode];
    double* coor = static_cast<double*>(stk::mesh::field_data(nodal_coord, node));
    STK_ThrowRequire(coor);
    for (unsigned i = 0; i < spatial_dim; ++i) {
      centroid[i] += coor[i];
    }
  }

  for(unsigned i = 0; i < spatial_dim; ++i) {
    centroid[i] /= num_nodes;
  }
}

TEST_F(StkToolsC, DeleteMeshExceptWithinBoundingBox)
{
  const std::string unNamed = "mesh not specified";
  const std::string inputMeshName = stk::unit_test_util::get_option("-i", unNamed);
  STK_ThrowRequireMsg(inputMeshName!=unNamed, "Please specify mesh with -i option.");
  setup_mesh(inputMeshName, stk::mesh::BulkData::NO_AUTO_AURA);

  const std::string outputMeshName = stk::unit_test_util::get_option("-o", "modified.g");

  double xLo = stk::unit_test_util::get_command_line_option("-x", std::numeric_limits<double>::lowest());
  double xHi = stk::unit_test_util::get_command_line_option("-X", std::numeric_limits<double>::max());
  double yLo = stk::unit_test_util::get_command_line_option("-y", std::numeric_limits<double>::lowest());
  double yHi = stk::unit_test_util::get_command_line_option("-Y", std::numeric_limits<double>::max());
  double zLo = stk::unit_test_util::get_command_line_option("-z", std::numeric_limits<double>::lowest());
  double zHi = stk::unit_test_util::get_command_line_option("-Z", std::numeric_limits<double>::max());

  const stk::mesh::BucketVector &buckets = get_bulk().get_buckets(stk::topology::ELEM_RANK, get_meta().locally_owned_part());
  const stk::mesh::FieldBase * coordinates = get_meta().coordinate_field();

  stk::mesh::EntityVector elementsToDestroy;
  std::vector<double> centroid(3);

  for (size_t i = 0; i < buckets.size(); ++i) {
    const stk::mesh::Bucket& bucket = *buckets[i];
    for (size_t j = 0; j < bucket.size(); ++j) {
      stk::mesh::Entity element = bucket[j];
      stk_determine_centroid(3, element, *coordinates, centroid);

      if ((centroid[0] < xLo) || (centroid[0] > xHi) ||
          (centroid[1] < yLo) || (centroid[1] > yHi) ||
          (centroid[2] < zLo) || (centroid[2] > zHi)) {
        elementsToDestroy.push_back(element);
      }
    }
  }

  destroy_elements(get_bulk(), elementsToDestroy);
  stk::io::write_mesh(outputMeshName, get_bulk());
}

TEST_F(StkToolsC, FlipElementConnectivity)
{
  const std::string unNamed = "mesh not specified";
  const std::string meshName = stk::unit_test_util::get_option("-i", unNamed);
  STK_ThrowRequireMsg(meshName!=unNamed, "Please specify mesh with -i option.");
  setup_mesh(meshName, stk::mesh::BulkData::NO_AUTO_AURA);

  int invalidBlockId = -1;
  int inputBlockId = stk::unit_test_util::get_command_line_option("-b", invalidBlockId);
  STK_ThrowRequireMsg(inputBlockId!=invalidBlockId, "Please specify block with -b.");

  std::ostringstream os;
  os << "block_" << inputBlockId;
  stk::mesh::Part* elemBlock = get_meta().get_part(os.str());

  stk::mesh::EntityVector elems;
  stk::mesh::BucketVector buckets = get_bulk().buckets(stk::topology::ELEM_RANK);
  stk::mesh::get_selected_entities(stk::mesh::Selector(*elemBlock), buckets, elems);

  get_bulk().modification_begin();
  for (auto elem : elems)
  {
    stk::topology topology = get_bulk().bucket(elem).topology();
    STK_ThrowRequireMsg(topology == stk::topology::HEX_8, "Input block must have HEX_8 topology but found topology " << topology);

    stk::mesh::EntityVector storedNodes;
    const stk::mesh::Entity* node = get_bulk().begin_nodes(elem);
    unsigned numNodes = get_bulk().num_nodes(elem);
    for (unsigned i = 0; i < numNodes; ++i)
    {
      storedNodes.push_back(node[i]);
    }
    for (unsigned i = 0; i < numNodes; ++i)
    {
      get_bulk().destroy_relation(elem, storedNodes[i], i);
    }
    for (unsigned i = 0; i < numNodes; ++i)
    {
      unsigned correctNodeId = (i+4)%8;
      get_bulk().declare_relation(elem, storedNodes[correctNodeId], i);
    }
  }
  get_bulk().modification_end();
  stk::io::write_mesh("flipped.e", get_bulk());
}
}
