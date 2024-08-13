#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/DestroyRelations.hpp>
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities
#include <stk_mesh/base/SkinBoundary.hpp> 
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraphUpdater.hpp>
#include <tokenize.h>                   // for tokenize
#include <cmath>                        // for atan2, cos, sin
#include <cstdlib>                      // for strtod, nullptr, strtol, exit, etc
#include <cstring>                      // for memcpy
#include <iomanip>                      // for operator<<, setw
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <string>                       // for string, operator==, etc
#include <vector>
#include <Ioss_EntityType.h>            // for EntityType, etc
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stddef.h>                     // for size_t
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker

namespace
{

std::string get_file_name(const std::string &parameters)
{
  std::string filename = parameters;
  size_t colon = filename.find(':');
  if (colon != std::string::npos && colon > 0) {
    filename = filename.substr(colon+1);
  }

  return filename;
}

void get_mesh_dimensions(const std::string &parameters, int &numX, int &numY, int &numZ)
{
  std::string filename = get_file_name(parameters);

  auto params = Ioss::tokenize(filename, "/");
  auto groups = Ioss::tokenize(params[params.size()-1], "|+");

  // First 'group' is the interval specification -- IxJxK
  auto tokens = Ioss::tokenize(groups[0], "x");
  assert(tokens.size() == 3);

  numX = std::strtol(tokens[0].c_str(), nullptr, 10);
  numY = std::strtol(tokens[1].c_str(), nullptr, 10);
  numZ = std::strtol(tokens[2].c_str(), nullptr, 10);
}

struct DeletedMeshInfo
{
  std::vector<stk::mesh::EntityId> elemIds;
  std::vector< std::vector<stk::mesh::Entity> > elemNodes;
};

struct NodeConnectivityInfo
{
  std::map<stk::mesh::Entity, int> nodeConnMap;
  std::vector<int> nodeConn;
};

typedef stk::mesh::Field<int> ScalarIntField;

class AnimateMeshDeletion
{
public:
  AnimateMeshDeletion(stk::mesh::BulkData &bulk, ScalarIntField &statusField, const std::string &animationFileName) :
    m_bulkData(bulk),
    m_stkIo(bulk.parallel()),
    m_fileHandler(0),
    m_statusField(statusField),
    m_animationFileName(animationFileName)
  {
    initialize_status_field();
    create_output_mesh();
  }

  ~AnimateMeshDeletion() {}

  void initialize_status_field()
  {
    stk::mesh::put_field_on_mesh(m_statusField, m_bulkData.mesh_meta_data().universal_part(), nullptr);
  }

  void create_output_mesh()
  {
    m_stkIo.set_bulk_data(m_bulkData);
    m_fileHandler = m_stkIo.create_output_mesh(m_animationFileName, stk::io::WRITE_RESULTS);
    m_stkIo.add_field(m_fileHandler, m_statusField);
    m_stkIo.write_output_mesh(m_fileHandler);
  }

  void deactivate_element(stk::mesh::EntityId elemId)
  {
    stk::mesh::Entity elem = m_bulkData.get_entity(stk::topology::ELEM_RANK, elemId);
    *stk::mesh::field_data(m_statusField, elem) = 1;
  }

  void animate(const std::vector<stk::mesh::EntityId> elemIds, int numSteps)
  {
    initialize_status_field();

    int numElem = elemIds.size();

    for (int m_time = 0; m_time <= numSteps; m_time++)
    {
      if((0 != m_time) && (m_time <= numElem))
      {
        stk::mesh::EntityId elemId = elemIds[m_time-1];
        deactivate_element(elemId);
      }

      m_stkIo.begin_output_step(m_fileHandler, m_time);
      m_stkIo.write_defined_output_fields(m_fileHandler);
      m_stkIo.end_output_step(m_fileHandler);
    }
  }

private:
  stk::mesh::BulkData &m_bulkData;
  stk::io::StkMeshIoBroker m_stkIo;
  size_t m_fileHandler;
  ScalarIntField &m_statusField;
  std::string m_animationFileName;
};

void declare_animation_field(stk::mesh::MetaData &metaData, std::string &animationFile)
{
  animationFile = stk::unit_test_util::get_option("-animationFile", "");

  if("" != animationFile)
  {
    ScalarIntField &deathField = metaData.declare_field<int>(stk::topology::ELEMENT_RANK, "death");
    stk::mesh::put_field_on_mesh(deathField, metaData.universal_part(), nullptr);
  }
}

void animate_death(stk::mesh::BulkData &bulkData, const std::vector<stk::mesh::EntityId> elemIds, const std::string &animationFile = "")
{
  if("" != animationFile)
  {
    int numSteps, numElem = elemIds.size();

    stk::all_reduce_max(bulkData.parallel(), &numElem, &numSteps, 1);

    ScalarIntField &deathField = static_cast<ScalarIntField &> (*bulkData.mesh_meta_data().get_field(stk::topology::ELEM_RANK, "death"));
    AnimateMeshDeletion deathAnimation(bulkData, deathField, animationFile);

    deathAnimation.animate(elemIds, numSteps);
  }
}

void print_delete_info(const stk::mesh::BulkData &bulkData, DeletedMeshInfo &meshInfo)
{
  int data[2];
  int data_sum[2], data_min, data_max;

  data[0] = meshInfo.elemIds.size();
  data[1] = stk::mesh::count_selected_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData.buckets(stk::topology::ELEM_RANK)); // localElems.size();

  stk::all_reduce_sum(bulkData.parallel(), data, data_sum, 2);
  stk::all_reduce_min(bulkData.parallel(), data, &data_min, 1);
  stk::all_reduce_max(bulkData.parallel(), data, &data_max, 1);

  if(0 == bulkData.parallel_rank())
  {
    double del_sum = (100.0*data_sum[0])/data_sum[1];
    std::cout << "[Delete Elements] " << data_sum[0] << " elements (" << del_sum << " % of mesh) Min : " << data_min << " Max: " << data_max << std::endl;
  }
}

void print_skin_info(const stk::mesh::BulkData &bulkData, const stk::mesh::Part& surface)
{
  stk::mesh::Selector selector(surface & bulkData.mesh_meta_data().locally_owned_part());
  std::vector<size_t> skin_count;
  stk::mesh::count_entities(selector, bulkData, skin_count);
  int minFaces, maxFaces, sumFaces, localFaces = skin_count[2];
  stk::all_reduce_min(bulkData.parallel(), &localFaces, &minFaces, 1);
  stk::all_reduce_max(bulkData.parallel(), &localFaces, &maxFaces, 1);
  stk::all_reduce_sum(bulkData.parallel(), &localFaces, &sumFaces, 1);

  if(bulkData.parallel_rank() == 0)
  {
    std::cout << "[Skinning] Min faces: " << minFaces << " Max faces: " << maxFaces << " Total faces: " << sumFaces << std::endl;
  }
}

void print_memory(stk::ParallelMachine communicator, const std::string &preamble)
{
  double maxHwmInMB = stk::get_max_hwm_across_procs(communicator) / (1024.0 * 1024.0);

  size_t curr_max, curr_min, curr_avg;
  stk::get_current_memory_usage_across_processors(communicator, curr_max, curr_min, curr_avg);
  double max = curr_max / (1024.0 * 1024.0);
  curr_min /= (1024.0 * 1024.0);
  curr_avg /= (1024.0 * 1024.0);

  if(stk::parallel_machine_rank(communicator) == 0)
  {
    std::ostringstream os;
    os << preamble << " high water mark = " << maxHwmInMB << " Mb, current max = " << max << " Mb" << std::endl;
    std::cerr << os.str();
  }
}


void get_ids_along_boundary(const stk::mesh::BulkData &bulkData, const std::string &meshSpec, DeletedMeshInfo &meshInfo)
{
  int nx, ny, nz;

  get_mesh_dimensions(meshSpec, nx, ny, nz);

  meshInfo.elemIds.resize(ny*nz);

  int count = 0;

  for(int z = 0; z < nz; ++z)
  {
    for(int y = 0; y < ny; ++y)
    {
      if((y%2 == 0) || ((ny-1) == y)) continue;

      int id = 1 + nx*(y + ny*z);

      stk::mesh::EntityId elemId = static_cast<stk::mesh::EntityId> (id);
      stk::mesh::Entity elem = bulkData.get_entity(stk::topology::ELEM_RANK, elemId);

      if(bulkData.is_valid(elem) && bulkData.bucket(elem).owned())
      {
        meshInfo.elemIds[count] = static_cast<stk::mesh::EntityId> (id);
        ++count;
      }
    }
  }

  meshInfo.elemIds.resize(count);

  print_delete_info(bulkData, meshInfo);
}

bool is_element_on_processor_boundary(const stk::mesh::ElemElemGraph& eeGraph, stk::mesh::Entity localElement)
{
  size_t numConn = eeGraph.get_num_connected_elems(localElement);

  bool onProcBoundary = false;
  for(size_t i=0; i<numConn; ++i)
  {
    if(!eeGraph.is_connected_elem_locally_owned(localElement, i))
      onProcBoundary = true;
  }
  return onProcBoundary;
}

bool will_deleting_element_create_orphaned_nodes(const stk::mesh::BulkData &bulkData,
                                                 stk::mesh::Entity localElement,
                                                 const NodeConnectivityInfo &nodeConnInfo)
{
  const stk::mesh::ConnectedEntities nodes = bulkData.get_connected_entities(localElement, stk::topology::NODE_RANK);
  const unsigned numNodes = nodes.size();

  for(unsigned i=0; i<numNodes; ++i)
  {
    stk::mesh::Entity node = nodes[i];
    int index = -1;
    std::map<stk::mesh::Entity, int>::const_iterator it = nodeConnInfo.nodeConnMap.find(node);

    if(nodeConnInfo.nodeConnMap.end() != it)
      index = it->second;

    if((-1 != index) && (nodeConnInfo.nodeConn[index] <= 1))
      return true;
  }

  return false;
}

void decrement_node_connectivity(const stk::mesh::BulkData &bulkData,
                                 stk::mesh::Entity localElement,
                                 NodeConnectivityInfo &nodeConnInfo)
{
  const stk::mesh::ConnectedEntities nodes = bulkData.get_connected_entities(localElement, stk::topology::NODE_RANK);
  const unsigned numNodes = nodes.size();

  for(unsigned i=0; i<numNodes; ++i)
  {
    stk::mesh::Entity node = nodes[i];

    int index = nodeConnInfo.nodeConnMap[node];
    nodeConnInfo.nodeConn[index]--;

    EXPECT_TRUE( nodeConnInfo.nodeConn[index] > 0);
  }
}

int num_local_elements(const stk::mesh::BulkData &bulkData, stk::mesh::Entity node)
{
  const stk::mesh::ConnectedEntities elements = bulkData.get_connected_entities(node, stk::topology::ELEM_RANK);
  int numElems = elements.size();
  int numLocalElems = numElems;

  for(int i=0; i<numElems; ++i)
  {
    stk::mesh::Entity elem = elements[i];
    if(!bulkData.bucket(elem).owned())
      numLocalElems--;
  }

  return numLocalElems;
}

void initialize_node_connectivity_map(const stk::mesh::BulkData &bulkData,
                                      NodeConnectivityInfo &nodeConnInfo)
{
  stk::mesh::EntityVector allNodes;
  stk::mesh::get_selected_entities(bulkData.mesh_meta_data().universal_part(), bulkData.buckets(stk::topology::NODE_RANK), allNodes);

  nodeConnInfo.nodeConn.resize(allNodes.size());

  for(unsigned i=0; i<allNodes.size(); ++i)
  {
    stk::mesh::Entity node = allNodes[i];
    nodeConnInfo.nodeConnMap[node] = i;
    nodeConnInfo.nodeConn[i] = num_local_elements( bulkData, node );
  }
}

void populate_interior_perforated_mesh_ids(const stk::mesh::BulkData &bulkData,
                                           const stk::mesh::ElemElemGraph& eeGraph,
                                           const stk::mesh::EntityVector &localElems,
                                           NodeConnectivityInfo &nodeConnInfo,
                                           DeletedMeshInfo &meshInfo)
{
  for(unsigned lastIndex = 0; lastIndex < localElems.size(); ++lastIndex)
  {
    stk::mesh::Entity localElement = localElems[lastIndex];
    stk::mesh::EntityId elemId = bulkData.identifier( localElement );

    if(!is_element_on_processor_boundary(eeGraph, localElement) &&
       !will_deleting_element_create_orphaned_nodes(bulkData, localElement, nodeConnInfo))
    {
      meshInfo.elemIds.push_back(elemId);
      decrement_node_connectivity(bulkData, localElement, nodeConnInfo);
    }
  }
}

void populate_processor_boundary_perforated_mesh_ids(const stk::mesh::BulkData &bulkData,
                                                     const stk::mesh::ElemElemGraph& eeGraph,
                                                     const stk::mesh::EntityVector &localElems,
                                                     NodeConnectivityInfo &nodeConnInfo,
                                                     DeletedMeshInfo &meshInfo)
{
  for(unsigned lastIndex = 0; lastIndex < localElems.size(); ++lastIndex)
  {
    stk::mesh::Entity localElement = localElems[lastIndex];
    stk::mesh::EntityId elemId = bulkData.identifier( localElement );

    if(is_element_on_processor_boundary(eeGraph, localElement) &&
       !will_deleting_element_create_orphaned_nodes(bulkData, localElement, nodeConnInfo))
    {
      meshInfo.elemIds.push_back(elemId);
      decrement_node_connectivity(bulkData, localElement, nodeConnInfo);
    }
  }
}

void get_perforated_mesh_ids(const stk::mesh::BulkData &bulkData, const stk::mesh::ElemElemGraph& eeGraph, DeletedMeshInfo &meshInfo)
{
  stk::mesh::EntityVector localElems;
  stk::mesh::get_selected_entities(bulkData.mesh_meta_data().locally_owned_part(), bulkData.buckets(stk::topology::ELEM_RANK), localElems);

  meshInfo.elemIds.reserve(localElems.size());

  NodeConnectivityInfo nodeConnInfo;

  initialize_node_connectivity_map(bulkData, nodeConnInfo);

  populate_interior_perforated_mesh_ids(bulkData, eeGraph, localElems, nodeConnInfo, meshInfo);

  populate_processor_boundary_perforated_mesh_ids(bulkData, eeGraph, localElems, nodeConnInfo, meshInfo);

  print_delete_info(bulkData, meshInfo);
}

class ElementGraphPerformance : public stk::unit_test_util::PerformanceTester
{
public:
  ElementGraphPerformance(stk::mesh::BulkData &bulk, const std::string &fileSpec)
    : stk::unit_test_util::PerformanceTester(bulk.parallel()),
      bulkData(bulk),
      meshSpec(fileSpec)
  {
  }

  virtual ~ElementGraphPerformance()
  {
  }

protected:
  void delete_elements(DeletedMeshInfo &meshInfo)
  {
    meshInfo.elemNodes.resize(meshInfo.elemIds.size());

    bulkData.modification_begin();
    for(size_t i = 0; i < meshInfo.elemIds.size(); ++i)
    {
      stk::mesh::EntityId elemId = meshInfo.elemIds[i];
      stk::mesh::Entity elem = bulkData.get_entity(stk::topology::ELEM_RANK, elemId);

      meshInfo.elemNodes[i].assign(bulkData.begin_nodes(elem), bulkData.end_nodes(elem));

      bulkData.destroy_entity(elem);
    }
    bulkData.modification_end();
  }

  void add_elements(DeletedMeshInfo &meshInfo)
  {
    bulkData.modification_begin();
    stk::mesh::Part &hexPart = bulkData.mesh_meta_data().get_topology_root_part(stk::topology::HEX_8);
    stk::mesh::Part *block_1 = bulkData.mesh_meta_data().get_part("block_1");
    stk::mesh::PartVector parts{&hexPart, block_1};
    for(size_t i = 0; i < meshInfo.elemIds.size(); ++i)
    {
      stk::mesh::EntityId elemId = meshInfo.elemIds[i];
      stk::mesh::Entity elem = bulkData.declare_element(elemId, parts);

      for(size_t j = 0; j<meshInfo.elemNodes[i].size(); ++j)
      {
        stk::mesh::Entity node = meshInfo.elemNodes[i][j];
        bulkData.declare_relation(elem, node, j);
      }
    }
    bulkData.modification_end();
  }

  virtual void run_algorithm_to_time()
  {
    // Add graph creation (and, consequently, deletion) to the timing
    bulkData.delete_face_adjacent_element_graph();
    bulkData.initialize_face_adjacent_element_graph();

    stk::mesh::ElemElemGraph & elemGraph = bulkData.get_face_adjacent_element_graph();
    DeletedMeshInfo meshInfo;

    unsigned oldGraphEdgeCount = elemGraph.num_edges();
    get_ids_along_boundary(bulkData, meshSpec, meshInfo);

    delete_elements(meshInfo);
    //stk::io::write_mesh("afterDelete.g", bulkData, bulkData.parallel());
    unsigned newGraphEdgeCount = elemGraph.num_edges();
    ASSERT_TRUE(oldGraphEdgeCount > newGraphEdgeCount);

    add_elements(meshInfo);

    unsigned finalGraphEdgeCount = elemGraph.num_edges();
    ASSERT_EQ(oldGraphEdgeCount, finalGraphEdgeCount);
  }
  virtual size_t get_value_to_output_as_iteration_count()
  {
    return 1;
  }

  stk::mesh::BulkData &bulkData;
  const std::string meshSpec;
};



class PerforatedElementGraphPerformance : public ElementGraphPerformance
{
public:
  PerforatedElementGraphPerformance(stk::mesh::BulkData &bulk, const std::string &fileSpec, const std::string &animateFile = "") :
    ElementGraphPerformance(bulk, fileSpec),
    animationFile(animateFile)
  {

  }

  ~PerforatedElementGraphPerformance()
  {

  }

  void delete_graph()
  {
    bulkData.delete_face_adjacent_element_graph();
  }

  void create_graph()
  {
    bulkData.initialize_face_adjacent_element_graph();
  }

  void delete_sides()
  {
    const stk::mesh::MetaData& meta = bulkData.mesh_meta_data();
    stk::mesh::EntityVector elems;
    stk::mesh::get_selected_entities(meta.locally_owned_part(),
                                     bulkData.buckets(stk::topology::ELEM_RANK), elems);
    bulkData.modification_begin();

    for(stk::mesh::Entity elem : elems)
    {
      stk::mesh::destroy_relations(bulkData, elem, meta.side_rank());
    }

    stk::mesh::EntityVector sides;
    stk::mesh::get_selected_entities(meta.locally_owned_part(),
                                     bulkData.buckets(meta.side_rank()), sides);
    for(stk::mesh::Entity side : sides)
    {
      bulkData.destroy_entity(side);
    }

    bulkData.modification_end();
  }

protected:
  virtual void run_algorithm_to_time()
  {
    stk::mesh::ElemElemGraph & elemGraph = bulkData.get_face_adjacent_element_graph();
    DeletedMeshInfo meshInfo;

    stk::mesh::Part& surface = bulkData.mesh_meta_data().declare_part("Skinned Surface");
    create_exposed_boundary(surface);

    unsigned oldGraphEdgeCount = elemGraph.num_edges();
    get_perforated_mesh_ids(bulkData, elemGraph, meshInfo);

    delete_elements(meshInfo);
    //stk::io::write_mesh("afterDelete.g", bulkData, bulkData.parallel());
    unsigned newGraphEdgeCount = elemGraph.num_edges();
    ASSERT_TRUE(oldGraphEdgeCount > newGraphEdgeCount);

    add_elements(meshInfo);
    //stk::io::write_mesh("afterAdd.g", bulkData, bulkData.parallel());
    unsigned finalGraphEdgeCount = elemGraph.num_edges();
    ASSERT_EQ(oldGraphEdgeCount, finalGraphEdgeCount);

    animate_death(bulkData, meshInfo.elemIds, animationFile);
  }

  void create_exposed_boundary(stk::mesh::Part& surface)
  {
    print_memory(bulkData.parallel(), "Before create exposed boundary sides");
    stk::mesh::create_exposed_block_boundary_sides(bulkData, bulkData.mesh_meta_data().universal_part(), {&surface});
    print_skin_info(bulkData, surface);
  }

  std::string animationFile;
};

class ElementGraphPerformanceTest : public stk::unit_test_util::MeshFixture
{
protected:
  void run_element_graph_perf_test()
  {
    ElementGraphPerformance perfTester(get_bulk(), get_mesh_spec());
    const int numTimes = 30;
    for (int i=0; i < numTimes; ++i) {
      perfTester.run_performance_test();
    }
  }
  std::string get_mesh_spec()
  {
    return stk::unit_test_util::get_option("-file", "NO_FILE_SPECIFIED");
  }
};

class PerforatedElementGraphPerformanceTest : public stk::unit_test_util::MeshFixture
{
protected:
  void run_element_graph_perf_test()
  {
    PerforatedElementGraphPerformance perfTester(get_bulk(), get_mesh_spec(), animationFile);
    const int numTimes = 12;
    for (int i=0; i < numTimes; ++i) {
      perfTester.delete_sides();
      perfTester.delete_graph();
      perfTester.create_graph();
      perfTester.run_performance_test();
    }
  }
  std::string get_mesh_spec()
  {
    return stk::unit_test_util::get_option("-file", "NO_FILE_SPECIFIED");
  }

  std::string animationFile = "";
};

TEST_F(ElementGraphPerformanceTest, read_mesh)
{
  print_memory(MPI_COMM_WORLD, "Before mesh creation");
  setup_mesh(get_mesh_spec(), stk::mesh::BulkData::AUTO_AURA);

  run_element_graph_perf_test();
}

TEST_F(ElementGraphPerformanceTest, read_mesh_no_aura)
{
  print_memory(MPI_COMM_WORLD, "Before mesh creation");
  setup_mesh(get_mesh_spec(), stk::mesh::BulkData::NO_AUTO_AURA);

  run_element_graph_perf_test();
}

TEST_F(PerforatedElementGraphPerformanceTest, read_mesh_with_auto_decomp)
{
  print_memory(MPI_COMM_WORLD, "Before mesh creation");
  allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
  declare_animation_field(get_meta(), animationFile);
  stk::unit_test_util::read_from_serial_file_and_decompose(get_mesh_spec(), get_bulk(), "rcb");

  run_element_graph_perf_test();
}

TEST_F(PerforatedElementGraphPerformanceTest, read_mesh_with_auto_decomp_no_aura)
{
  print_memory(MPI_COMM_WORLD, "Before mesh creation");
  allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
  declare_animation_field(get_meta(), animationFile);
  stk::unit_test_util::read_from_serial_file_and_decompose(get_mesh_spec(), get_bulk(), "rcb");

  run_element_graph_perf_test();
}

}

