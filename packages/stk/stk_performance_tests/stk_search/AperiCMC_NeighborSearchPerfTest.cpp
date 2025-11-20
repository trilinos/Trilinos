#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/OutputStreams.hpp>

#include <Kokkos_Core.hpp>
#include <memory>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_mesh/base/NgpForEachEntity.hpp>
#include <stk_mesh/base/NgpMesh.hpp>
#include <stk_mesh/base/NgpFieldBLAS.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_search/BoxIdent.hpp>
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/LocalCoarseSearch.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/Point.hpp>
#include <stk_search/SearchMethod.hpp>
#include <stk_search/Sphere.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_unit_test_utils/timer.hpp>

static constexpr size_t MAX_NEIGHBORS = 40;

class NeighborSearch
{
  using ExecSpace = stk::ngp::ExecSpace;
  using NodeIdent = stk::mesh::Entity::entity_value_type;
  using NodeIdentProc = stk::search::IdentProc<stk::mesh::EntityId, int>;
  using SphereIdent = stk::search::BoxIdent<stk::search::Sphere<double>, NodeIdent>;
  using SphereIdentProc = stk::search::BoxIdentProc<stk::search::Sphere<double>, NodeIdentProc>;
  using PointIdent = stk::search::BoxIdent<stk::search::Point<double>, NodeIdent>;
  using PointIdentProc = stk::search::BoxIdentProc<stk::search::Point<double>, NodeIdentProc>;
  using LocalIntersection = stk::search::IdentIntersection<NodeIdent, NodeIdent>;
  using Intersection = stk::search::IdentProcIntersection<NodeIdentProc, NodeIdentProc>;

  using LocalRangeViewType = Kokkos::View<SphereIdent *, ExecSpace>;
  using LocalDomainViewType = Kokkos::View<PointIdent *, ExecSpace>;
  using LocalResultViewType = Kokkos::View<LocalIntersection *, ExecSpace>;

  using RangeViewType = Kokkos::View<SphereIdentProc *, ExecSpace>;
  using DomainViewType = Kokkos::View<PointIdentProc *, ExecSpace>;
  using ResultViewType = Kokkos::View<Intersection *, ExecSpace>;

public:
  using DoubleField = stk::mesh::Field<double>;
  using NgpDoubleField = stk::mesh::NgpField<double>;

  NeighborSearch(std::shared_ptr<stk::mesh::BulkData> bulkData,
                 const std::vector<std::string> &blockNames = {})
  : m_bulkData(bulkData)
  {
    stk::mesh::MetaData& meta = m_bulkData->mesh_meta_data();

    if (blockNames.size() == 0) {
      m_selector = stk::mesh::Selector(meta.universal_part());
    } else {
      stk::mesh::PartVector parts;
      for (const auto &blockName : blockNames) {
        stk::mesh::Part *part = meta.get_part(blockName);
        STK_ThrowRequireMsg(part != nullptr, "Block " << blockName << " not found.");
        parts.push_back(part);
      }
      m_selector = stk::mesh::selectUnion(parts);
    }

    m_owned_selector = m_selector & m_bulkData->mesh_meta_data().locally_owned_part();

    // Get the node number of neighbors field
    m_numNeighborsField = &meta.get_field<double>(stk::topology::NODE_RANK, "num_neighbors")->field_of_state(stk::mesh::StateNone);

    // Get the node neighbors field
    m_neighborsField = &meta.get_field<double>(stk::topology::NODE_RANK, "neighbors")->field_of_state(stk::mesh::StateNone);

    // Get the coordinates field
    m_coordField = static_cast<const DoubleField*>(meta.coordinate_field());

    // Get the kernel radius field
    m_kernelRadiusField = &meta.get_field<double>(stk::topology::NODE_RANK, "kernel_radius")->field_of_state(stk::mesh::StateNone);

    // Get the function values field
    m_functionValuesField = &meta.get_field<double>(stk::topology::NODE_RANK, "function_values")->field_of_state(stk::mesh::StateNone);
  }

  void compute_kernel_radius()
  {
    const stk::mesh::NgpMesh &ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulkData);
    // Get the ngp fields
    auto ngpCoordField = stk::mesh::get_updated_ngp_field<double>(*m_coordField);
    ngpCoordField.sync_to_device();
    auto ngpKernelRadiusField = stk::mesh::get_updated_ngp_field<double>(*m_kernelRadiusField);
    ngpKernelRadiusField.clear_sync_state();
    const double tolerance = std::numeric_limits<double>::epsilon();

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, m_selector,
      KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex &node_index) {
        // Get the node's coordinates
        stk::mesh::EntityFieldData<double> coordinates = ngpCoordField(node_index);

        // Get the kernel radius
        double kernel_radius_squared = 0.0;
        stk::mesh::NgpMesh::ConnectedEntities connected_entities = ngpMesh.get_connected_entities(stk::topology::NODE_RANK, node_index, stk::topology::ELEMENT_RANK);
        for (size_t i = 0; i < connected_entities.size(); ++i) {
          stk::mesh::FastMeshIndex elem_index = ngpMesh.fast_mesh_index(connected_entities[i]);
          stk::mesh::NgpMesh::ConnectedNodes connected_nodes = ngpMesh.get_nodes(stk::topology::ELEM_RANK, elem_index);
          for (size_t j = 0; j < connected_nodes.size(); ++j) {
            stk::mesh::FastMeshIndex neighbor_index = ngpMesh.fast_mesh_index(connected_nodes[j]);
            stk::mesh::EntityFieldData<double> neighbor_coordinates = ngpCoordField(neighbor_index);
            double length_squared = 0;
            for (size_t k = 0; k < 3; ++k) {
              const double value = coordinates[k] - neighbor_coordinates[k];
              length_squared += value * value;
            }
            kernel_radius_squared = Kokkos::max(kernel_radius_squared, length_squared);
          }
        }
        const double kernel_radius = Kokkos::sqrt(kernel_radius_squared);
        ngpKernelRadiusField(node_index, 0) = kernel_radius + tolerance;
      });
    ngpKernelRadiusField.modify_on_device();
  }

  // The domain will be the nodes. Search will find the nodes within the spheres from above.
  // The identifiers will be the node entity-offsets
  LocalDomainViewType create_local_node_points() {
    const stk::mesh::MetaData &meta = m_bulkData->mesh_meta_data();
    stk::mesh::Selector ownedOrShared = m_owned_selector | meta.globally_shared_part();
    const unsigned num_local_nodes = stk::mesh::count_entities(*m_bulkData, stk::topology::NODE_RANK, ownedOrShared);
    LocalDomainViewType node_points("node_points", num_local_nodes);

    auto& ngpCoordField = stk::mesh::get_updated_ngp_field<double>(*m_coordField);
    ngpCoordField.sync_to_device();
    const stk::mesh::NgpMesh &ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulkData);

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, ownedOrShared,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& node_index) {
      stk::mesh::EntityFieldData<double> coords = ngpCoordField(node_index);
      stk::mesh::Entity node = ngpMesh.get_entity(stk::topology::NODE_RANK, node_index);
      const unsigned nodeLocalId = ngpMesh.local_id(node);
      node_points(nodeLocalId) = PointIdent{stk::search::Point<double>(coords[0], coords[1], coords[2]), NodeIdent(node.local_offset())};
    });

    return node_points;
  }

  // The domain will be the nodes. Search will find the nodes within the spheres from above.
  // The identifiers will be the global node ids.
  DomainViewType create_node_points() {
    const stk::mesh::MetaData &meta = m_bulkData->mesh_meta_data();
    stk::mesh::Selector ownedOrShared = m_owned_selector | meta.globally_shared_part();
    const unsigned num_local_nodes = stk::mesh::count_entities(*m_bulkData, stk::topology::NODE_RANK, ownedOrShared);
    DomainViewType node_points("node_points", num_local_nodes);

    auto& ngpCoordField = stk::mesh::get_updated_ngp_field<double>(*m_coordField);
    ngpCoordField.sync_to_device();
    const stk::mesh::NgpMesh &ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulkData);
    const int myRank = m_bulkData->parallel_rank();

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, ownedOrShared,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& node_index) {
      stk::mesh::EntityFieldData<double> coords = ngpCoordField(node_index);
      stk::mesh::Entity node = ngpMesh.get_entity(stk::topology::NODE_RANK, node_index);
      const unsigned nodeLocalId = ngpMesh.local_id(node);
      node_points(nodeLocalId) = PointIdentProc{stk::search::Point<double>(coords[0], coords[1], coords[2]), NodeIdentProc(node.local_offset(), myRank)};
    });

    return node_points;
  }

  // Sphere range. Will be used to find the nodes within a ball defined by the sphere.
  // The identifiers will be the entity local_offsets
  LocalRangeViewType create_local_node_spheres(const stk::mesh::Selector& selector) {
    const unsigned numNodes = stk::mesh::count_entities(*m_bulkData, stk::topology::NODE_RANK, selector);
    LocalRangeViewType node_spheres("node_spheres", numNodes);

    auto& ngpCoordField = stk::mesh::get_updated_ngp_field<double>(*m_coordField);
    ngpCoordField.sync_to_device();
    auto& ngpKernelRadiusField = stk::mesh::get_updated_ngp_field<double>(*m_kernelRadiusField);
    ngpKernelRadiusField.sync_to_device();

    const stk::mesh::NgpMesh &ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulkData);

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, selector,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& node_index) {
      stk::mesh::EntityFieldData<double> coords = ngpCoordField(node_index);
      stk::search::Point<double> center(coords[0], coords[1], coords[2]);
      stk::mesh::Entity node = ngpMesh.get_entity(stk::topology::NODE_RANK, node_index);
      const unsigned nodeLocalId = ngpMesh.local_id(node);
      double radius = ngpKernelRadiusField(node_index, 0);
      node_spheres(nodeLocalId) = SphereIdent{stk::search::Sphere<double>(center, radius), NodeIdent(node.local_offset())};
    });

    return node_spheres;
  }

  // Sphere range. Will be used to find the nodes within a ball defined by the sphere.
  // The identifiers will be the entity local_offsets
  RangeViewType create_node_spheres(const stk::mesh::Selector& selector) {
    const unsigned numNodes = stk::mesh::count_entities(*m_bulkData, stk::topology::NODE_RANK, selector);
    RangeViewType node_spheres("node_spheres", numNodes);

    auto& ngpCoordField = stk::mesh::get_updated_ngp_field<double>(*m_coordField);
    ngpCoordField.sync_to_device();
    auto& ngpKernelRadiusField = stk::mesh::get_updated_ngp_field<double>(*m_kernelRadiusField);
    ngpKernelRadiusField.sync_to_device();

    const stk::mesh::NgpMesh &ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulkData);
    const int myRank = m_bulkData->parallel_rank();

    stk::mesh::for_each_entity_run(ngpMesh, stk::topology::NODE_RANK, selector,
    KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& node_index) {
      stk::mesh::EntityFieldData<double> coords = ngpCoordField(node_index);
      stk::search::Point<double> center(coords[0], coords[1], coords[2]);
      stk::mesh::Entity node = ngpMesh.get_entity(stk::topology::NODE_RANK, node_index);
      const unsigned nodeLocalId = ngpMesh.local_id(node);
      double radius = ngpKernelRadiusField(node_index, 0);
      node_spheres(nodeLocalId) = SphereIdentProc{stk::search::Sphere<double>(center, radius), NodeIdentProc(node.local_offset(), myRank)};
    });

    return node_spheres;
  }

  // Ghost the neighbors to the nodes processor
  void ghost_node_neighbors(const ResultViewType::host_mirror_type &host_search_results)
  {
    m_bulkData->modification_begin();
    stk::mesh::Ghosting &neighbor_ghosting = m_bulkData->create_ghosting("neighbors");
    std::vector<stk::mesh::EntityProc> nodes_to_ghost;

    const int myRank = m_bulkData->parallel_rank();

    for (size_t i = 0; i < host_search_results.size(); ++i) {
      auto result = host_search_results(i);
      if (result.domainIdentProc.proc() != myRank && result.rangeIdentProc.proc() == myRank) {
        stk::mesh::Entity node(result.rangeIdentProc.id());
        nodes_to_ghost.emplace_back(node, result.domainIdentProc.proc());
      }
    }

    m_bulkData->change_ghosting(neighbor_ghosting, nodes_to_ghost);
    m_bulkData->modification_end();
  }

  // Put the search results into the neighbors field. The neighbors field is a field of entity offsets. The neighbors are sorted by distance. Near to far.
  void unpack_search_results_into_field(const LocalResultViewType& search_results)
  {
    auto ngpCoordField = stk::mesh::get_updated_ngp_field<double>(*m_coordField);
    auto ngpNeighborsField = stk::mesh::get_updated_ngp_field<double>(*m_neighborsField);
    auto ngpNumNeighborsField = stk::mesh::get_updated_ngp_field<double>(*m_numNeighborsField);
    auto ngpFunctionValuesField = stk::mesh::get_updated_ngp_field<double>(*m_functionValuesField);

    ngpCoordField.sync_to_device();
    ngpNeighborsField.sync_to_device();
    ngpNumNeighborsField.sync_to_device();
    ngpFunctionValuesField.sync_to_device();

    const stk::mesh::NgpMesh &ngpMesh = stk::mesh::get_updated_ngp_mesh(*m_bulkData);

    const size_t numResults = search_results.size();
    Kokkos::parallel_for("unpack_search_results", Kokkos::RangePolicy<stk::ngp::ExecSpace>(0, numResults),
      KOKKOS_LAMBDA(const int& i)
      {
        auto result = search_results(i);
        unsigned domainId = result.domainIdent;
        const bool shouldUnpack = i==0 || search_results(i-1).domainIdent != domainId;
        if (shouldUnpack) {
          unsigned ii = i;
          while(ii<search_results.size() && result.domainIdent == domainId) {
            stk::mesh::Entity node(result.domainIdent);
            stk::mesh::Entity neighbor(result.rangeIdent);
            stk::mesh::FastMeshIndex nodeIndex = ngpMesh.device_mesh_index(node);
            stk::mesh::FastMeshIndex neighborIndex = ngpMesh.device_mesh_index(neighbor);

            const stk::mesh::EntityFieldData<double> neighborCoords = ngpCoordField(neighborIndex);
            const stk::mesh::EntityFieldData<double> nodeCoords = ngpCoordField(nodeIndex);
            stk::mesh::EntityFieldData<double> p_neighbor_data = ngpNeighborsField(nodeIndex);
            double& num_neighbors = ngpNumNeighborsField(nodeIndex,0);
            stk::mesh::EntityFieldData<double> p_function_values = ngpFunctionValuesField(nodeIndex);  // Using the function values field as a temporary storage for the squared distances

            // Calculate the squared distance between the node and the neighbor
            double distance_squared = 0.0;
            for (size_t j = 0; j < 3; ++j) {
              const double value = neighborCoords[j] - nodeCoords[j];
              distance_squared += value * value;
            }
            // Find where to insert the neighbor, based on the distance
            unsigned insert_index = (unsigned)num_neighbors;  // Default to the end of the list
            for (size_t j = 0; j < insert_index; ++j) {
              if (distance_squared < p_function_values[j]) {
                insert_index = j;
                break;
              }
            }

            // Shift the function values and neighbors to make room for the new neighbor
            unsigned reverse_start_index = (unsigned)num_neighbors;
            if (reverse_start_index == MAX_NEIGHBORS) {
              Kokkos::printf("Node %ld has too many neighbors. The furthest neighbor will be removed.\n", ngpMesh.identifier(node));
              --reverse_start_index;
            } else {
              num_neighbors += 1;
            }
            for (size_t j = reverse_start_index; j > insert_index; --j) {
              p_function_values[j] = p_function_values[j - 1];
              p_neighbor_data[j] = p_neighbor_data[j - 1];
            }

            // Insert the new neighbor
            p_function_values[insert_index] = distance_squared;
            p_neighbor_data[insert_index] = (double)neighbor.local_offset();
            ++ii;
            result = search_results(ii);
          }
        }
      });
    m_neighborsField->modify_on_device();
    m_numNeighborsField->modify_on_device();
  }

  void do_ball_search()
  {
    DomainViewType node_points = create_node_points();
    RangeViewType node_spheres = create_node_spheres(m_owned_selector);
    m_kernelRadiusField->clear_sync_state();

    ResultViewType search_results;
    stk::search::SearchMethod search_method = stk::search::MORTON_LBVH;

    stk::ngp::ExecSpace exec_space = Kokkos::DefaultExecutionSpace{};
    const bool results_parallel_symmetry = true;

    stk::search::coarse_search(node_points, node_spheres, search_method, m_bulkData->parallel(), search_results, exec_space, results_parallel_symmetry);

    ResultViewType::host_mirror_type host_search_results = Kokkos::create_mirror_view(search_results);
    Kokkos::deep_copy(host_search_results, search_results);

    ghost_node_neighbors(host_search_results);

    stk::mesh::Selector all = m_bulkData->mesh_meta_data().universal_part();
    LocalDomainViewType local_node_points = create_local_node_points();
    LocalRangeViewType local_node_spheres = create_local_node_spheres(all);

    LocalResultViewType local_search_results;
    const bool sortResults = true;
    stk::search::local_coarse_search(local_node_points, local_node_spheres, search_method, local_search_results, exec_space, sortResults);

    unpack_search_results_into_field(local_search_results);
  }

  void add_nodes_neighbors_within_variable_ball() {
    compute_kernel_radius();
    do_ball_search();
  }

  void check_neighbor_stats(int nElemX, int nElemY, int nElemZ,
                            int expectedMinNeighbors, int expectedMaxNeighbors)
  {
    // Initialize the min and max values
    double max_num_neighbors = 0;
    double min_num_neighbors = std::numeric_limits<double>::max();
    double total_num_neighbors = 0;
    double num_entities = stk::mesh::count_entities(*m_bulkData, stk::topology::NODE_RANK, m_owned_selector);
    NgpDoubleField& ngpNumNeighborsField = stk::mesh::get_updated_ngp_field<double>(*m_numNeighborsField);
    NgpDoubleField& ngpNeighborsField = stk::mesh::get_updated_ngp_field<double>(*m_neighborsField);
    const bool newestDataIsOnDevice = !ngpNumNeighborsField.need_sync_to_device()
                                   && !ngpNeighborsField.need_sync_to_device();

    if (newestDataIsOnDevice) {
      ngpNumNeighborsField.sync_to_host();
      ngpNeighborsField.sync_to_host();
    }

    const stk::mesh::BucketVector& nodeBuckets = m_bulkData->get_buckets(stk::topology::NODE_RANK, m_owned_selector);
    auto numNeighborsData = m_numNeighborsField->data();
    for(const stk::mesh::Bucket* bptr : nodeBuckets) {
      for(stk::mesh::Entity node : *bptr) {
        auto numNeighbors = numNeighborsData.entity_values(node);
        max_num_neighbors = std::max(max_num_neighbors, numNeighbors(0_comp));
        min_num_neighbors = std::min(min_num_neighbors, numNeighbors(0_comp));
        total_num_neighbors += numNeighbors(0_comp);
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &max_num_neighbors, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &min_num_neighbors, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &total_num_neighbors, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &num_entities, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    EXPECT_EQ(min_num_neighbors, expectedMinNeighbors);//neighbors include 'self'
    EXPECT_EQ(max_num_neighbors, expectedMaxNeighbors);
  //  const double avg_num_neighbors = total_num_neighbors / num_entities;
  //  EXPECT_NEAR(avg_num_neighbors, 14.518518, 0.001);
    const size_t expectedTotalNodes = (nElemX + 1) * (nElemY + 1) * (nElemZ + 1);
    EXPECT_EQ(num_entities, expectedTotalNodes);
  }

  void zero_fields()
  {
    stk::ngp::ExecSpace execSpace = Kokkos::DefaultExecutionSpace{};

    constexpr double zero = 0.0;
    stk::mesh::field_fill(zero, *m_numNeighborsField, execSpace);
    stk::mesh::field_fill(zero, *m_neighborsField, execSpace);
    stk::mesh::field_fill(zero, *m_kernelRadiusField, execSpace);
    stk::mesh::field_fill(zero, *m_functionValuesField, execSpace);
  }

private:
  std::shared_ptr<stk::mesh::BulkData> m_bulkData;
  stk::mesh::Selector m_selector;
  stk::mesh::Selector m_owned_selector;
  stk::mesh::NgpMesh m_ngpMesh;
  const DoubleField *m_coordField;
  DoubleField *m_numNeighborsField;
  DoubleField *m_neighborsField;
  DoubleField *m_kernelRadiusField;
  DoubleField *m_functionValuesField;
};

class NeighborSearchTestFixture : public ::testing::Test
{
protected:
  void SetUp() override { }

  void create_mesh_and_search_manager(MPI_Comm comm, const std::string &mesh_spec)
  {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr =
        stk::mesh::MeshBuilder(comm)
             .set_spatial_dimension(3)
             .set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA)
             .set_maintain_local_ids(true)
             .create();
    stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();

    NeighborSearch::DoubleField& nodeNumNeighborsField = meta.declare_field<double>(stk::topology::NODE_RANK, "num_neighbors", 1);
    stk::mesh::put_field_on_entire_mesh(nodeNumNeighborsField, 1);

    NeighborSearch::DoubleField& nodeNeighborsField = meta.declare_field<double>(stk::topology::NODE_RANK, "neighbors", 1);
    stk::mesh::put_field_on_entire_mesh(nodeNeighborsField, MAX_NEIGHBORS);

    NeighborSearch::DoubleField& functionValuesField = meta.declare_field<double>(stk::topology::NODE_RANK, "function_values", 1);
    stk::mesh::put_field_on_entire_mesh(functionValuesField, MAX_NEIGHBORS);

    NeighborSearch::DoubleField& kernelRadiusField = meta.declare_field<double>(stk::topology::NODE_RANK, "kernel_radius", 1);
    stk::mesh::put_field_on_entire_mesh(kernelRadiusField, 1);

    stk::io::fill_mesh(mesh_spec, *bulkPtr);

    m_searchManager = std::make_shared<NeighborSearch>(bulkPtr, std::vector<std::string>{"block_1"});
  }

  std::shared_ptr<NeighborSearch> m_searchManager;
};

std::string get_mesh_spec(int nX, int nY, int nZ, double bboxScaleFactor = 1.0)
{
  std::string meshSpec = "generated:" + std::to_string(nX) + "x" + std::to_string(nY) + "x" + std::to_string(nZ);

  if (bboxScaleFactor != 1.0) {
    meshSpec += "|bbox:";
    meshSpec += std::to_string(-bboxScaleFactor * nX) + ","
              + std::to_string(-bboxScaleFactor * nY) + ","
              + std::to_string(-bboxScaleFactor * nZ);
    meshSpec += ","
              + std::to_string(bboxScaleFactor * nX) +  ","
              + std::to_string(bboxScaleFactor * nY) +  ","
              + std::to_string(bboxScaleFactor * nZ);
  }

  return meshSpec;
}

TEST_F(NeighborSearchTestFixture, VariableBallSearchUnitCubes) {
  MPI_Comm comm = stk::parallel_machine_world();
  const int numProcs = stk::parallel_machine_size(comm);
  if (numProcs > 5) {
      GTEST_SKIP_("Test only runs with 5 or fewer processes.");
  }

  int numElemX = 50;
  int numElemY = 50;
  int numElemZ = 10*numProcs;
  std::string meshSpec = get_mesh_spec(numElemX, numElemY, numElemZ);
  stk::outputP0()<<"meshSpec: "<<meshSpec<<std::endl;

  create_mesh_and_search_manager(comm, meshSpec);

  stk::unit_test_util::BatchTimer batchTimer(comm);
  batchTimer.initialize_batch_timer();
  const unsigned NUM_RUNS = 50;

  for(unsigned i=0; i<NUM_RUNS; ++i) {
    batchTimer.start_batch_timer();

    m_searchManager->zero_fields();
    m_searchManager->add_nodes_neighbors_within_variable_ball();

    batchTimer.stop_batch_timer();
  }

  int minNeighbors = 8;
  int maxNeighbors = 27;
  m_searchManager->check_neighbor_stats(numElemX, numElemY, numElemZ, minNeighbors, maxNeighbors);

  batchTimer.print_batch_timing(NUM_RUNS);
}

TEST_F(NeighborSearchTestFixture, VariableBallSearchScaledCubes) {
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 5) {
      GTEST_SKIP_("Test only runs with 5 or fewer processes.");
  }

  int numElemX = 3;
  int numElemY = 3;
  int numElemZ = 4;
  double bboxScaleFactor = 1.0 / 7.0;
  std::string meshSpec = get_mesh_spec(numElemX, numElemY, numElemZ, bboxScaleFactor);
  stk::outputP0()<<"meshSpec: "<<meshSpec<<std::endl;

  create_mesh_and_search_manager(comm, meshSpec);

  m_searchManager->add_nodes_neighbors_within_variable_ball();
  int minNeighbors = 8;
  int maxNeighbors = 27;
  m_searchManager->check_neighbor_stats(numElemX, numElemY, numElemZ, minNeighbors, maxNeighbors);
}
