#include <gtest/gtest.h>
#include <mpi.h>

#include <Kokkos_Core.hpp>
#include <array>
#include <memory>
#include <stk_io/StkMeshIoBroker.hpp>
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
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_search/BoxIdent.hpp>
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/Point.hpp>
#include <stk_search/SearchMethod.hpp>
#include <stk_search/Sphere.hpp>
#include <stk_topology/topology.hpp>
#include <stk_util/ngp/NgpSpaces.hpp>

using DoubleField = stk::mesh::Field<double>;
using NgpDoubleField = stk::mesh::NgpField<double>;
static constexpr size_t MAX_NEIGHBORS = 40;

class NeighborSearch {
    using ExecSpace = stk::ngp::ExecSpace;
    using NodeIdentProc = stk::search::IdentProc<stk::mesh::EntityId, int>;
    using SphereIdentProc = stk::search::BoxIdentProc<stk::search::Sphere<double>, NodeIdentProc>;
    using PointIdentProc = stk::search::BoxIdentProc<stk::search::Point<double>, NodeIdentProc>;
    using Intersection = stk::search::IdentProcIntersection<NodeIdentProc, NodeIdentProc>;

    using RangeViewType = Kokkos::View<SphereIdentProc *, ExecSpace>;
    using DomainViewType = Kokkos::View<PointIdentProc *, ExecSpace>;
    using ResultViewType = Kokkos::View<Intersection *, ExecSpace>;

    using FastMeshIndicesViewType = Kokkos::View<stk::mesh::FastMeshIndex *, ExecSpace>;

   public:
    NeighborSearch(stk::mesh::BulkData *bulk_data, const std::vector<std::string> &sets = {}) : m_bulk_data(bulk_data), m_sets(sets) {
        m_ngp_mesh = stk::mesh::get_updated_ngp_mesh(*m_bulk_data);
        stk::mesh::MetaData *meta_data = &m_bulk_data->mesh_meta_data();

        if (sets.size() == 0) {
            m_selector = stk::mesh::Selector(meta_data->universal_part());
        } else {
            stk::mesh::PartVector parts;
            for (const auto &set : sets) {
                stk::mesh::Part *part = meta_data->get_part(set);
                if (part == nullptr) {
                    throw std::runtime_error("Set " + set + " not found.");
                }
                parts.push_back(part);
            }
            m_selector = stk::mesh::selectUnion(parts);
        }
        // Warn if the selector is empty.
        if (m_selector.is_empty(stk::topology::ELEMENT_RANK)) {
            std::cout << "Warning: NeighborSearch selector is empty." << std::endl;
        }

        stk::mesh::Selector full_owned_selector = m_bulk_data->mesh_meta_data().locally_owned_part();
        m_owned_selector = m_selector & full_owned_selector;

        // Get the node number of neighbors field
        m_node_num_neighbors_field = &meta_data->get_field<double>(stk::topology::NODE_RANK, "num_neighbors")->field_of_state(stk::mesh::StateNone);
        m_ngp_node_num_neighbors_field = &stk::mesh::get_updated_ngp_field<double>(*m_node_num_neighbors_field);

        // Get the node neighbors field
        m_node_neighbors_field = &meta_data->get_field<double>(stk::topology::NODE_RANK, "neighbors")->field_of_state(stk::mesh::StateNone);
        m_ngp_node_neighbors_field = &stk::mesh::get_updated_ngp_field<double>(*m_node_neighbors_field);

        // Get the coordinates field
        m_coordinates_field = &meta_data->get_field<double>(stk::topology::NODE_RANK, m_bulk_data->mesh_meta_data().coordinate_field_name())->field_of_state(stk::mesh::StateNone);
        m_ngp_coordinates_field = &stk::mesh::get_updated_ngp_field<double>(*m_coordinates_field);

        // Get the kernel radius field
        m_kernel_radius_field = &meta_data->get_field<double>(stk::topology::NODE_RANK, "kernel_radius")->field_of_state(stk::mesh::StateNone);
        m_ngp_kernel_radius_field = &stk::mesh::get_updated_ngp_field<double>(*m_kernel_radius_field);

        // Get the function values field
        m_function_values_field = &meta_data->get_field<double>(stk::topology::NODE_RANK, "function_values")->field_of_state(stk::mesh::StateNone);
    }

    void ComputeKernelRadius(double scale_factor) {
        auto ngp_mesh = m_ngp_mesh;
        // Get the ngp fields
        auto ngp_coordinates_field = *m_ngp_coordinates_field;
        auto ngp_kernel_radius_field = *m_ngp_kernel_radius_field;
        const double tolerance = std::numeric_limits<double>::epsilon();

        stk::mesh::for_each_entity_run(
            ngp_mesh, stk::topology::NODE_RANK, m_selector,
            KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex &node_index) {
                // Get the node's coordinates
                stk::mesh::EntityFieldData<double> coordinates = ngp_coordinates_field(node_index);

                // Get the kernel radius
                double kernel_radius_squared = 0.0;
                stk::mesh::NgpMesh::ConnectedEntities connected_entities = ngp_mesh.get_connected_entities(stk::topology::NODE_RANK, node_index, stk::topology::ELEMENT_RANK);
                for (size_t i = 0; i < connected_entities.size(); ++i) {
                    stk::mesh::FastMeshIndex elem_index = ngp_mesh.fast_mesh_index(connected_entities[i]);
                    stk::mesh::NgpMesh::ConnectedNodes connected_nodes = ngp_mesh.get_nodes(stk::topology::ELEM_RANK, elem_index);
                    for (size_t j = 0; j < connected_nodes.size(); ++j) {
                        stk::mesh::FastMeshIndex neighbor_index = ngp_mesh.fast_mesh_index(connected_nodes[j]);
                        stk::mesh::EntityFieldData<double> neighbor_coordinates = ngp_coordinates_field(neighbor_index);
                        double length_squared = 0;
                        for (size_t k = 0; k < 3; ++k) {
                            const double value = coordinates[k] - neighbor_coordinates[k];
                            length_squared += value * value;
                        }
                        kernel_radius_squared = Kokkos::max(kernel_radius_squared, length_squared);
                    }
                }
                const double kernel_radius = Kokkos::sqrt(kernel_radius_squared);
                ngp_kernel_radius_field(node_index, 0) = kernel_radius * scale_factor + tolerance;
            });
        ngp_kernel_radius_field.clear_sync_state();
        ngp_kernel_radius_field.modify_on_device();
        ngp_kernel_radius_field.sync_to_host();
    }

    // Create local entities on host and copy to device
    FastMeshIndicesViewType GetLocalEntityIndices(stk::mesh::EntityRank rank, stk::mesh::Selector selector) {
        std::vector<stk::mesh::Entity> local_entities;
        stk::mesh::get_entities(*m_bulk_data, rank, selector, local_entities);

        FastMeshIndicesViewType mesh_indices("mesh_indices", local_entities.size());
        FastMeshIndicesViewType::HostMirror host_mesh_indices = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, mesh_indices);

        for (size_t i = 0; i < local_entities.size(); ++i) {
            const stk::mesh::MeshIndex &mesh_index = m_bulk_data->mesh_index(local_entities[i]);
            host_mesh_indices(i) = stk::mesh::FastMeshIndex{mesh_index.bucket->bucket_id(), mesh_index.bucket_ordinal};
        }

        Kokkos::deep_copy(mesh_indices, host_mesh_indices);
        return mesh_indices;
    }

    // The domain will be the nodes. Search will find the nodes within the spheres from above.
    // The identifiers will be the global node ids.
    DomainViewType CreateNodePoints() {
        const stk::mesh::MetaData &meta = m_bulk_data->mesh_meta_data();
        const unsigned num_local_nodes = stk::mesh::count_entities(*m_bulk_data, stk::topology::NODE_RANK, m_owned_selector | meta.globally_shared_part());
        DomainViewType node_points("node_points", num_local_nodes);

        auto ngp_coordinates_field = *m_ngp_coordinates_field;
        const stk::mesh::NgpMesh &ngp_mesh = m_ngp_mesh;

        // Slow host operation that is needed to get an index. There is plans to add this to the stk::mesh::NgpMesh.
        FastMeshIndicesViewType node_indices = GetLocalEntityIndices(stk::topology::NODE_RANK, m_owned_selector | meta.globally_shared_part());
        const int my_rank = m_bulk_data->parallel_rank();

        Kokkos::parallel_for(
            stk::ngp::DeviceRangePolicy(0, num_local_nodes), KOKKOS_LAMBDA(const unsigned &i) {
                stk::mesh::EntityFieldData<double> coords = ngp_coordinates_field(node_indices(i));
                stk::mesh::Entity node = ngp_mesh.get_entity(stk::topology::NODE_RANK, node_indices(i));
                node_points(i) = PointIdentProc{stk::search::Point<double>(coords[0], coords[1], coords[2]), NodeIdentProc(ngp_mesh.identifier(node), my_rank)};
            });

        return node_points;
    }

    // Sphere range. Will be used to find the nodes within a ball defined by the sphere.
    // The identifiers will be the global node ids.
    RangeViewType CreateNodeSpheres() {
        const unsigned num_local_nodes = stk::mesh::count_entities(*m_bulk_data, stk::topology::NODE_RANK, m_owned_selector);
        RangeViewType node_spheres("node_spheres", num_local_nodes);

        auto ngp_coordinates_field = *m_ngp_coordinates_field;
        auto ngp_kernel_radius_field = *m_ngp_kernel_radius_field;
        const stk::mesh::NgpMesh &ngp_mesh = m_ngp_mesh;

        // Slow host operation that is needed to get an index. There is plans to add this to the stk::mesh::NgpMesh.
        FastMeshIndicesViewType node_indices = GetLocalEntityIndices(stk::topology::NODE_RANK, m_owned_selector);
        const int my_rank = m_bulk_data->parallel_rank();

        Kokkos::parallel_for(
            stk::ngp::DeviceRangePolicy(0, num_local_nodes), KOKKOS_LAMBDA(const unsigned &i) {
                stk::mesh::EntityFieldData<double> coords = ngp_coordinates_field(node_indices(i));
                stk::search::Point<double> center(coords[0], coords[1], coords[2]);
                stk::mesh::Entity node = ngp_mesh.get_entity(stk::topology::NODE_RANK, node_indices(i));
                double radius = ngp_kernel_radius_field(node_indices(i), 0);
                node_spheres(i) = SphereIdentProc{stk::search::Sphere<double>(center, radius), NodeIdentProc(ngp_mesh.identifier(node), my_rank)};
            });

        return node_spheres;
    }

    // Ghost the neighbors to the nodes processor
    void GhostNodeNeighbors(const ResultViewType::HostMirror &host_search_results) {
        m_bulk_data->modification_begin();
        stk::mesh::Ghosting &neighbor_ghosting = m_bulk_data->create_ghosting("neighbors");
        std::vector<stk::mesh::EntityProc> nodes_to_ghost;

        const int my_rank = m_bulk_data->parallel_rank();

        for (size_t i = 0; i < host_search_results.size(); ++i) {
            auto result = host_search_results(i);
            if (result.domainIdentProc.proc() != my_rank && result.rangeIdentProc.proc() == my_rank) {
                stk::mesh::Entity node = m_bulk_data->get_entity(stk::topology::NODE_RANK, result.rangeIdentProc.id());
                nodes_to_ghost.emplace_back(node, result.domainIdentProc.proc());
            }
        }

        m_bulk_data->change_ghosting(neighbor_ghosting, nodes_to_ghost);
        m_bulk_data->modification_end();
    }

    // Put the search results into the neighbors field. The neighbors field is a field of global node ids. The neighbors are sorted by distance. Near to far.
    void UnpackSearchResultsIntoField(const ResultViewType::HostMirror &host_search_results) {
        const int my_rank = m_bulk_data->parallel_rank();

        for (size_t i = 0; i < host_search_results.size(); ++i) {
            auto result = host_search_results(i);
            if (result.domainIdentProc.proc() == my_rank) {
                stk::mesh::Entity node = m_bulk_data->get_entity(stk::topology::NODE_RANK, result.domainIdentProc.id());
                stk::mesh::Entity neighbor = m_bulk_data->get_entity(stk::topology::NODE_RANK, result.rangeIdentProc.id());
                const double *p_neighbor_coordinates = stk::mesh::field_data(*m_coordinates_field, neighbor);
                const double *p_node_coordinates = stk::mesh::field_data(*m_coordinates_field, node);
                double *p_neighbor_data = stk::mesh::field_data(*m_node_neighbors_field, node);
                double &num_neighbors = *stk::mesh::field_data(*m_node_num_neighbors_field, node);
                double *p_function_values = stk::mesh::field_data(*m_function_values_field, node);  // Using the function values field as a temporary storage for the squared distances

                // Calculate the squared distance between the node and the neighbor
                double distance_squared = 0.0;
                for (size_t j = 0; j < 3; ++j) {
                    const double value = p_neighbor_coordinates[j] - p_node_coordinates[j];
                    distance_squared += value * value;
                }

                // Find where to insert the neighbor, based on the distance
                size_t insert_index = (size_t)num_neighbors;  // Default to the end of the list
                for (size_t j = 0; j < insert_index; ++j) {
                    if (distance_squared < p_function_values[j]) {
                        insert_index = j;
                        break;
                    }
                }

                // Shift the function values and neighbors to make room for the new neighbor
                size_t reverse_start_index = (size_t)num_neighbors;
                if (reverse_start_index == MAX_NEIGHBORS) {
                    Kokkos::printf("Node %ld has too many neighbors. The furthest neighbor will be removed.\n", m_bulk_data->identifier(node));
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
            }
        }
        // Never communicate the neighbors field. The shared nodes need to have a processor local value and not the value of the owning processor.
        m_node_neighbors_field->modify_on_host();
        m_node_num_neighbors_field->modify_on_host();
        m_node_neighbors_field->sync_to_device();
        m_node_num_neighbors_field->sync_to_device();
    }

    void DoBallSearch() {
        DomainViewType node_points = CreateNodePoints();
        RangeViewType node_spheres = CreateNodeSpheres();

        ResultViewType search_results;
        stk::search::SearchMethod search_method = stk::search::MORTON_LBVH;

        stk::ngp::ExecSpace exec_space = Kokkos::DefaultExecutionSpace{};
        const bool results_parallel_symmetry = true;

        stk::search::coarse_search(node_points, node_spheres, search_method, m_bulk_data->parallel(), search_results, exec_space, results_parallel_symmetry);

        ResultViewType::HostMirror host_search_results = Kokkos::create_mirror_view(search_results);
        Kokkos::deep_copy(host_search_results, search_results);

        // Print sizes
        std::cout << "Neighborhood Search Information:" << std::endl;
        std::cout << "\n  Search Point-Sphere Pair Results Size: " << host_search_results.size()
                  << "\n  Evaluation Points Size: " << node_points.size()
                  << "\n  Neighbor Spheres Size: " << node_spheres.size() << std::endl;

        GhostNodeNeighbors(host_search_results);

        UnpackSearchResultsIntoField(host_search_results);
    }

    void add_nodes_neighbors_within_variable_ball(double scale_factor) {
        ComputeKernelRadius(scale_factor);
        DoBallSearch();
    }

    std::map<std::string, double> GetNumNeighborStats() {
        // Initialize the min and max values
        double max_num_neighbors = 0;
        double min_num_neighbors = std::numeric_limits<double>::max();
        double total_num_neighbors = 0;
        double num_entities = 0;
        NgpDoubleField ngp_num_neighbors_field;

        num_entities = stk::mesh::count_entities(*m_bulk_data, stk::topology::NODE_RANK, m_owned_selector);
        ngp_num_neighbors_field = *m_ngp_node_num_neighbors_field;
        ngp_num_neighbors_field.sync_to_host();

        const stk::mesh::FieldBase* hostField = ngp_num_neighbors_field.get_field_base();
        const stk::mesh::BucketVector& nodeBuckets = m_bulk_data->get_buckets(stk::topology::NODE_RANK, m_owned_selector);
        for(const stk::mesh::Bucket* bptr : nodeBuckets) {
          for(stk::mesh::Entity node : *bptr) {
            const double* numNeighbors = reinterpret_cast<const double*>(stk::mesh::field_data(*hostField, node));
            max_num_neighbors = std::max(max_num_neighbors, numNeighbors[0]);
            min_num_neighbors = std::min(min_num_neighbors, numNeighbors[0]);
            total_num_neighbors += numNeighbors[0];
          }
        }

        // Use MPI_Allreduce to calculate the min, max, and sum across all MPI ranks
        MPI_Allreduce(MPI_IN_PLACE, &max_num_neighbors, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &min_num_neighbors, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &total_num_neighbors, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &num_entities, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        std::map<std::string, double> stats;
        stats["max_num_neighbors"] = max_num_neighbors;
        stats["min_num_neighbors"] = min_num_neighbors;
        stats["avg_num_neighbors"] = total_num_neighbors / num_entities;
        stats["num_entities"] = num_entities;
        return stats;
    }

    void PrintNumNeighborsStats() {
        // Node
        std::map<std::string, double> node_stats = GetNumNeighborStats();

        std::cout << "Node Stats: " << std::endl;
        std::cout << "    Total Num Nodes: " << node_stats["num_entities"] << std::endl;
        std::cout << "  Max Num Neighbors: " << node_stats["max_num_neighbors"] << std::endl;
        std::cout << "  Min Num Neighbors: " << node_stats["min_num_neighbors"] << std::endl;
        std::cout << "  Avg Num Neighbors: " << node_stats["avg_num_neighbors"] << std::endl
                  << std::endl;  // Add a new line for readability
    }

    void SyncFieldsToHost() {
        m_ngp_node_num_neighbors_field->sync_to_host();
        m_ngp_kernel_radius_field->clear_sync_state();
        m_ngp_kernel_radius_field->sync_to_host();
    }

   private:
    stk::mesh::BulkData *m_bulk_data;                // The bulk data object.
    std::vector<std::string> m_sets;                 // The sets to process.
    stk::mesh::Selector m_selector;                  // The selector
    stk::mesh::Selector m_owned_selector;            // The local selector
    stk::mesh::NgpMesh m_ngp_mesh;                   // The ngp mesh object.
    DoubleField *m_coordinates_field;                // The coordinates field
    DoubleField *m_node_num_neighbors_field;         // The number of neighbors field
    DoubleField *m_node_neighbors_field;             // The neighbors field
    DoubleField *m_kernel_radius_field;              // The kernel radius field
    DoubleField *m_function_values_field;            // The function values field
    NgpDoubleField *m_ngp_coordinates_field;         // The ngp coordinates field
    NgpDoubleField *m_ngp_node_num_neighbors_field;  // The ngp number of neighbors field
    NgpDoubleField *m_ngp_node_neighbors_field;      // The ngp neighbors field
    NgpDoubleField *m_ngp_kernel_radius_field;       // The ngp kernel radius field
};

class NeighborSearchTestFixture : public ::testing::Test {
   protected:
    void SetUp() override {
    }

    void CreateMeshAndProcessors(const std::string &mesh_spec) {
        MPI_Comm p_communicator = MPI_COMM_WORLD;
        m_bulk_data = stk::mesh::MeshBuilder(p_communicator).create();
        stk::mesh::MetaData *p_meta_data = &m_bulk_data->mesh_meta_data();

        stk::io::StkMeshIoBroker mesh_reader;
        mesh_reader.set_bulk_data(*m_bulk_data);
        mesh_reader.add_mesh_database(mesh_spec, stk::io::READ_MESH);
        mesh_reader.create_input_mesh();
        mesh_reader.add_all_mesh_fields_as_input_fields();

        // Create the fields, start with nodes
        m_node_num_neighbors_field = &p_meta_data->declare_field<double>(stk::topology::NODE_RANK, "num_neighbors", 1);
        stk::mesh::put_field_on_entire_mesh(*m_node_num_neighbors_field, 1);

        m_node_neighbors_field = &p_meta_data->declare_field<double>(stk::topology::NODE_RANK, "neighbors", 1);
        stk::mesh::put_field_on_entire_mesh(*m_node_neighbors_field, MAX_NEIGHBORS);

        m_node_neighbors_function_values_field = &p_meta_data->declare_field<double>(stk::topology::NODE_RANK, "function_values", 1);
        stk::mesh::put_field_on_entire_mesh(*m_node_neighbors_function_values_field, MAX_NEIGHBORS);

        m_kernel_radius_field = &p_meta_data->declare_field<double>(stk::topology::NODE_RANK, "kernel_radius", 1);
        stk::mesh::put_field_on_entire_mesh(*m_kernel_radius_field, 1);

        mesh_reader.populate_bulk_data();

        // Create the NeighborSearch
        m_search_processor = std::make_shared<NeighborSearch>(m_bulk_data.get(), std::vector<std::string>{"block_1"});
    }

    std::shared_ptr<stk::mesh::BulkData> m_bulk_data;
    DoubleField *m_node_num_neighbors_field;
    DoubleField *m_node_neighbors_field;
    DoubleField *m_node_neighbors_function_values_field;
    DoubleField *m_kernel_radius_field;
    std::shared_ptr<NeighborSearch> m_search_processor;
};

TEST_F(NeighborSearchTestFixture, VariableBallSearchUnitCubes) {
    // Unit cube elements. Should give same answer on CPU and GPU.
    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    if (num_procs > 5) {
        GTEST_SKIP_("Test only runs with 5 or fewer processes.");
    }
    int num_elements_x = 1;
    int num_elements_y = 1;
    int num_elements_z = 5;
    std::string mesh_spec = "generated:" + std::to_string(num_elements_x) + "x" + std::to_string(num_elements_y) + "x" + std::to_string(num_elements_z);
    std::cout<<"mesh_spec: "<<mesh_spec<<std::endl;
    CreateMeshAndProcessors(mesh_spec);
    double ball_scale_factor = 1.0;
    m_search_processor->add_nodes_neighbors_within_variable_ball(ball_scale_factor);
    m_search_processor->SyncFieldsToHost();

    // Check the neighbor stats
    std::map<std::string, double> node_neighbor_stats = m_search_processor->GetNumNeighborStats();
    // Expected results are hard-coded to CPU results
    EXPECT_EQ(node_neighbor_stats["min_num_neighbors"], 7);
    EXPECT_EQ(node_neighbor_stats["max_num_neighbors"], 10);
    EXPECT_NEAR(node_neighbor_stats["avg_num_neighbors"], 9.0, 0.001);
    size_t expected_num_nodes = (num_elements_x + 1) * (num_elements_y + 1) * (num_elements_z + 1);
    EXPECT_EQ(node_neighbor_stats["num_entities"], expected_num_nodes);
}

TEST_F(NeighborSearchTestFixture, VariableBallSearchScaledCubes) {
    // Scaled cube elements. Noticing different answers on CPU and GPU.
    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    if (num_procs > 5) {
        GTEST_SKIP_("Test only runs with 5 or fewer processes.");
    }
    int num_elements_x = 1;
    int num_elements_y = 1;
    int num_elements_z = 5;
    std::string mesh_spec = "generated:" + std::to_string(num_elements_x) + "x" + std::to_string(num_elements_y) + "x" + std::to_string(num_elements_z);
    double bbox_scale_factor = 1.0 / 7.0; // Rational factor to exacerbate the CPU vs GPU differences. Adding a tolerance in the radius calculation fixes the differences.
    mesh_spec += "|bbox:";
    mesh_spec += "-" + std::to_string(bbox_scale_factor * num_elements_x) + ",-" + std::to_string(bbox_scale_factor * num_elements_y) + ",-" + std::to_string(bbox_scale_factor * num_elements_z);
    mesh_spec += "," + std::to_string(bbox_scale_factor * num_elements_x) +  "," + std::to_string(bbox_scale_factor * num_elements_y) +  "," + std::to_string(bbox_scale_factor * num_elements_z);
    std::cout<<"mesh_spec: "<<mesh_spec<<std::endl;
    CreateMeshAndProcessors(mesh_spec);
    double ball_scale_factor = 1.0;
    m_search_processor->add_nodes_neighbors_within_variable_ball(ball_scale_factor);
    m_search_processor->SyncFieldsToHost();

    // Check the neighbor stats
    std::map<std::string, double> node_neighbor_stats = m_search_processor->GetNumNeighborStats();
    // Expected results are hard-coded to CPU results
    EXPECT_EQ(node_neighbor_stats["min_num_neighbors"], 7);
    EXPECT_EQ(node_neighbor_stats["max_num_neighbors"], 11);
    EXPECT_NEAR(node_neighbor_stats["avg_num_neighbors"], 9.166667, 0.001);
    size_t expected_num_nodes = (num_elements_x + 1) * (num_elements_y + 1) * (num_elements_z + 1);
    EXPECT_EQ(node_neighbor_stats["num_entities"], expected_num_nodes);
}
