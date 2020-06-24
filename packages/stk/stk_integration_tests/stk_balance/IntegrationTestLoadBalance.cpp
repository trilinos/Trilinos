#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include <cmath>
#include <gtest/gtest.h>
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_unit_test_utils/StkMeshFromGeneratedMesh.h>
#include <stk_util/environment/WallTime.hpp>
#include <test_utils/OptionsForTesting.hpp>
#include <integrationtest/MeshUtilsForBoundingVolumes.hpp>
#include <unit_tests/UnitTestUtils.hpp>

#include <Teuchos_ParameterList.hpp>
#include <stk_balance/internal/StkMeshAdapterForZoltan2.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_ColoringProblem.hpp>

#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/StkBalanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>

#include <zoltan.h>

#include <algorithm>

#include <exodusII.h>

#include <test_utils/BalanceTestUtilities.hpp>

#include <stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp>

#include <stk_tools/mesh_clone/MeshClone.hpp>

namespace
{

struct LoadBalanceDiagnostics
{
    std::vector<double> sumOfVertexWeightsPerProc;
    std::vector<double> sumOfCutEdgeWeightsPerProc;
    std::vector<int> numElementsPerProc;
};

void writeParFiles(stk::io::StkMeshIoBroker &ioBroker, const std::string &output_file_name);
void fillIoBroker(MPI_Comm communicator, const std::string &generatedMeshSpec, stk::io::StkMeshIoBroker &ioBroker);
void setUpDefaultColoring(const stk::mesh::BulkData &stkMeshBulkData, std::vector<int>& coloring);
void createMockElementDecompositon(const int procId, stk::mesh::EntityProcVec &mockDecomposition, const stk::mesh::EntityVector& entities);
void verifyMeshAfterRebalance(stk::mesh::BulkData &stkMeshBulkData);
void verifyMeshPriorToRebalance(stk::mesh::BulkData &stkMeshBulkData);

template<typename GlobalId, typename LocalNumber>
void writeDotFile(const std::string &fileName, const std::vector<GlobalId>& globalIds, const std::vector<LocalNumber> &offsets, const std::vector<GlobalId>& adjacency);

std::string getSubdomainPartName(int subdomainId);

template <typename GlobalIds>
void gatherLoadBalanceDiagnostics(const std::vector<double> &vertexWeights, const std::vector<double> &edgeWeights, const std::vector<GlobalIds> &adjacency, stk::mesh::BulkData &stkMeshBulkData, MPI_Comm communicator, struct LoadBalanceDiagnostics &diagnostics);
void printLoadBalanceDiagnostics(const struct LoadBalanceDiagnostics &loadBalanceDiagnostics);
void checkMeshIsLoadBalanced(const stk::balance::BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData);

//====================


class ColoringSettings : public stk::balance::BalanceSettings
{
public:
    ColoringSettings() {}

    virtual size_t getNumNodesRequiredForConnection(stk::topology element1Topology, stk::topology element2Topology) const
    {
        return 1;
    }
    virtual double getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const
    {
        return 1.0;
    }
    virtual int getGraphVertexWeight(stk::topology type) const
    {
        return 1;
    }
    virtual GraphOption getGraphOption() const
    {
        return BalanceSettings::LOAD_BALANCE;
    }

    virtual double getGraphVertexWeight(stk::mesh::Entity entity, int criteria_index = 0) const
    {
        return 1.0;
    }
};

class OurColoringSettings : public ColoringSettings
{
public:
    OurColoringSettings() {}

    virtual GraphOption getGraphOption() const
    {
        return BalanceSettings::COLOR_MESH;
    }
};

TEST(LoadBalance, writeMesh)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);
    int myProcId = -1;
    MPI_Comm_rank(communicator, &myProcId);
    Options options = getOptionsForTest("generated:1x1x10|sideset:xXyYzZ");

    if(numProcs == 2)
    {
        stk::io::StkMeshIoBroker ioBroker(communicator);
        fillIoBroker(communicator, options.getMeshFileName(), ioBroker);
        stk::mesh::BulkData &outputBulkData = ioBroker.bulk_data();

        const std::string output_file_name = "output.exo";
        if(myProcId == 0)
        {
            balance_utils::clearFiles(output_file_name, numProcs);
        }

        writeParFiles(ioBroker, output_file_name);
        MPI_Barrier(communicator);

        stk::io::StkMeshIoBroker inputBroker(communicator);
        fillIoBroker(communicator, output_file_name, inputBroker);
        stk::mesh::BulkData &inputBulkData = inputBroker.bulk_data();

        std::vector<size_t> outputCounts;
        stk::mesh::count_entities(outputBulkData.mesh_meta_data().universal_part(), outputBulkData, outputCounts);
        std::vector<size_t> inputCounts;
        stk::mesh::count_entities(inputBulkData.mesh_meta_data().universal_part(), inputBulkData, inputCounts);

        EXPECT_TRUE(inputCounts == outputCounts);

        if(myProcId == 0)
        {
            balance_utils::clearFiles(output_file_name, numProcs);
        }
    }
}

TEST(LoadBalance, DISABLED_moveElementToAnotherProcessor)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);
    int myProcId = -1;
    MPI_Comm_rank(communicator, &myProcId);
    Options options = getOptionsForTest("generated:1x1x10|sideset:xXyYzZ");

    if(numProcs == 2)
    {
        stk::io::StkMeshIoBroker ioBroker(communicator);
        fillIoBroker(communicator, options.getMeshFileName(), ioBroker);

        stk::mesh::BulkData &stkMeshBulkData = ioBroker.bulk_data();

        // Move element 9 from proc 2 to proc 0
        std::vector<std::pair<stk::mesh::Entity, int> > entityProcPairs;
        stk::mesh::EntityKey entityKey(stk::topology::ELEMENT_RANK, 9);
        stk::mesh::Entity elementToMove = stkMeshBulkData.get_entity(entityKey);

        if(myProcId == 0)
        {
            EXPECT_FALSE(stkMeshBulkData.is_valid(elementToMove));
        }
        else
        {
            EXPECT_TRUE(stkMeshBulkData.is_valid(elementToMove));
            EXPECT_TRUE(stkMeshBulkData.bucket(elementToMove).owned());
        }

        if(myProcId == 1)
        {
            entityProcPairs.push_back(std::make_pair(elementToMove, 0));
            unsigned numFaces = stkMeshBulkData.num_faces(elementToMove);
            const stk::mesh::Entity *faces = stkMeshBulkData.begin_faces(elementToMove);
            for(unsigned int i = 0; i < numFaces; i++)
            {
                if(stkMeshBulkData.bucket(faces[i]).owned())
                {
                    entityProcPairs.push_back(std::make_pair(faces[i], 0));
                }
            }
            unsigned numEdges = stkMeshBulkData.num_edges(elementToMove);
            const stk::mesh::Entity *edges = stkMeshBulkData.begin_edges(elementToMove);
            for(unsigned int i = 0; i < numEdges; i++)
            {
                if(stkMeshBulkData.bucket(edges[i]).owned())
                {
                    entityProcPairs.push_back(std::make_pair(edges[i], 0));
                }
            }
        }

        stkMeshBulkData.change_entity_owner(entityProcPairs);

        stk::mesh::Entity elementAfterMove = stkMeshBulkData.get_entity(entityKey);

        if(myProcId == 0)
        {
            EXPECT_TRUE(stkMeshBulkData.is_valid(elementAfterMove));
            EXPECT_TRUE(stkMeshBulkData.bucket(elementAfterMove).owned());
        }
        else
        {
            EXPECT_TRUE(stkMeshBulkData.is_valid(elementAfterMove));
            EXPECT_TRUE(stkMeshBulkData.bucket(elementAfterMove).in_aura());
        }
    }
}

TEST(LoadBalance, DISABLED_specifyWhichProcessorYouWantEachElementToBeOnAndWriteOutParFiles)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);
    int myProcId = -1;
    MPI_Comm_rank(communicator, &myProcId);
    const std::string generatedMeshSpec = "generated:1x1x4|sideset:xXyYzZ";

    if(numProcs == 2)
    {
        stk::io::StkMeshIoBroker ioBroker(communicator);
        fillIoBroker(communicator, generatedMeshSpec, ioBroker);

        stk::mesh::BulkData &stkMeshBulkData = ioBroker.bulk_data();
        verifyMeshPriorToRebalance(stkMeshBulkData);

        stk::mesh::EntityVector elements;
        stk::mesh::get_selected_entities(stkMeshBulkData.mesh_meta_data().locally_owned_part(), stkMeshBulkData.buckets(stk::topology::ELEM_RANK), elements);

        stk::mesh::EntityProcVec mockDecomposition;
        createMockElementDecompositon(stkMeshBulkData.parallel_rank(), mockDecomposition, elements);

        stk::balance::internal::rebalance(stkMeshBulkData, mockDecomposition);

        verifyMeshAfterRebalance(stkMeshBulkData);
    }
}



TEST(LoadBalance, DISABLED_Zoltan2Parmetis)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);
    int procId;
    MPI_Comm_rank(communicator, &procId);

    Options options = getOptionsForTest("generated:1x1x4|sideset:xXyYzZ");

    if(numProcs == 2 || options.overRideTest())
    {
        stk::balance::internal::logMessage(communicator, "Creating mesh");

        stk::io::StkMeshIoBroker ioBroker(communicator);
        fillIoBroker(communicator, options.getMeshFileName(), ioBroker);

        stk::balance::internal::logMessage(communicator, "Finished creating mesh");

        stk::mesh::BulkData &stkMeshBulkData = ioBroker.bulk_data();
        if(!options.overRideTest())
        {
            verifyMeshPriorToRebalance(stkMeshBulkData);
        }

        stk::balance::internal::logMessage(communicator, "Starting to balance mesh");

        //stk::balance::BasicZoltan2Settings defaultSettings;
        stk::balance::GraphCreationSettings graphOptions;
        graphOptions.setToleranceForFaceSearch( options.getToleranceForFaceSearch() );
        graphOptions.setToleranceForParticleSearch( options.getToleranceForParticleSearch() );
        graphOptions.setDecompMethod( options.getPartmetisMethod() );

        stk::balance::balanceStkMesh(graphOptions, stkMeshBulkData);

        if ( !options.overRideTest() ) checkMeshIsLoadBalanced(graphOptions, stkMeshBulkData);

        stk::balance::internal::logMessage(communicator, "Finished balancing mesh");


        stk::balance::internal::logMessage(communicator, "Writing files");

        double time_start  = stk::wall_time();

        std::string nonConstString = options.getOutputFilename();
        if( nonConstString.empty() ) nonConstString = "output.exo";

        const std::string output_file_name = nonConstString;
        writeParFiles(ioBroker, output_file_name);

        if(procId == 0 && options.deleteFiles()) balance_utils::clearFiles(output_file_name, numProcs);


        double time_end = stk::wall_time();
        double writetime = time_end-time_start;
        std::ostringstream os;
        os << "IO Time: " << writetime;
        stk::balance::internal::logMessage(communicator, os.str());
    }
}

// This test balances then colors a graph.
// If balancing is optional, then it could be deleted,
// simplifying the test.

TEST(LoadBalance, zoltan2Coloring)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);
    int procId;
    MPI_Comm_rank(communicator, &procId);

    Options options = getOptionsForTest("generated:1x1x4|sideset:xXyYzZ");

    if(options.overRideTest())
    {
        stk::balance::internal::logMessage(communicator, "Creating mesh");

        stk::io::StkMeshIoBroker ioBroker(communicator);
        fillIoBroker(communicator, options.getMeshFileName(), ioBroker);

        stk::balance::internal::logMessage(communicator, "Finished creating mesh");

        stk::mesh::BulkData &stkMeshBulkData = ioBroker.bulk_data();
        if(!options.overRideTest())
        {
            verifyMeshPriorToRebalance(stkMeshBulkData);
        }

        stk::balance::internal::logMessage(communicator, "Starting to balance mesh");

        stk::balance::BasicZoltan2Settings graphSettings;
        stk::balance::balanceStkMesh(graphSettings, stkMeshBulkData);

        if ( !options.overRideTest() ) checkMeshIsLoadBalanced(graphSettings, stkMeshBulkData);

        stk::balance::internal::logMessage(communicator, "Finished balancing mesh");

        std::vector<int> coloring;
        setUpDefaultColoring(stkMeshBulkData, coloring);
        balance_utils::putFieldDataOnMesh(stkMeshBulkData, coloring);


        stk::balance::internal::logMessage(communicator, "Writing files");

        double time_start  = stk::wall_time();

        std::string nonConstString = options.getOutputFilename();
        if( nonConstString.empty() ) nonConstString = "output.exo";

        const std::string output_file_name = nonConstString;
        writeParFiles(ioBroker, output_file_name);

        if(procId == 0 && options.deleteFiles()) balance_utils::clearFiles(output_file_name, numProcs);


        double time_end = stk::wall_time();
        double writetime = time_end-time_start;
        std::ostringstream os;
        os << "IO Time: " << writetime;
        stk::balance::internal::logMessage(communicator, os.str());
    }

}

TEST(LoadBalance, zoltan2Adapter)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);

    if(numProcs == 1)
    {
        const std::string generatedMeshSpec = "generated:1x1x1";
        stk::io::StkMeshIoBroker ioBroker(communicator);
        fillIoBroker(communicator, generatedMeshSpec, ioBroker);

        stk::mesh::BulkData &stkMeshBulkData = ioBroker.bulk_data();

        stk::mesh::impl::LocalIdMapper localIds(stkMeshBulkData, stk::topology::ELEM_RANK);

        stk::balance::GraphCreationSettings graphSettings;
        Zoltan2ParallelGraph zoltan2Graph;
        std::vector<int> adjacencyProcs;
        zoltan2Graph.fillZoltan2AdapterDataFromStkMesh(stkMeshBulkData,
                                                       graphSettings,
                                                       adjacencyProcs,
                                                       stkMeshBulkData.mesh_meta_data().locally_owned_part(),
                                                       localIds);

        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts(stkMeshBulkData, counts);
        zoltan2Graph.set_num_global_elements(counts[stk::topology::ELEM_RANK]);

        zoltan2Graph.set_spatial_dim(stkMeshBulkData.mesh_meta_data().spatial_dimension());
        StkMeshZoltanAdapter stkMeshAdapter(zoltan2Graph);

        const int numEntityRanks = 4;
        const size_t numEntitiesOfRank[numEntityRanks] = {
                0, //nodes
                0, //edges
                0, //faces
                1}; //elements

        const enum Zoltan2::MeshEntityType zoltan2MeshEntityTypes[numEntityRanks] = {
                Zoltan2::MESH_VERTEX, //nodes
                Zoltan2::MESH_EDGE, //edges
                Zoltan2::MESH_FACE, //faces
                Zoltan2::MESH_REGION}; //elements

        for(int i = 0; i < numEntityRanks; i++)
        {
            ASSERT_EQ(numEntitiesOfRank[i], stkMeshAdapter.getLocalNumOf(zoltan2MeshEntityTypes[i]));
            const BalanceGlobalNumber *ids;
            stkMeshAdapter.getIDsViewOf(zoltan2MeshEntityTypes[i], ids);
            std::vector<int> sortedIds(ids, ids + numEntitiesOfRank[i]);
            std::sort(sortedIds.begin(), sortedIds.end());
            for(int j = 0; j < static_cast<int>(numEntitiesOfRank[i]); j++)
            {
                EXPECT_EQ(j+1, sortedIds[j])<< "failed for zoltan2 MeshEntityType " << zoltan2MeshEntityTypes[i];
            }
        }

        int fromTo[4][4] = {
                { 0, 0, 0, 0},
                { 0, 0, 0, 0},
                { 0, 0, 0, 0},
                { 0, 0, 1, 0}
        };

        for(int i = 0; i < numEntityRanks; i++)
        {
            for (int j=0;j<numEntityRanks;j++)
            {
                int numWeightsPer2ndAdj = stkMeshAdapter.getNumWeightsPer2ndAdj(zoltan2MeshEntityTypes[i], zoltan2MeshEntityTypes[j]);
                EXPECT_EQ(fromTo[i][j], numWeightsPer2ndAdj);
            }
        }

        const double *coordinates;
        int stride;

        for(int coordinateDirection = 0; coordinateDirection < stkMeshAdapter.getDimension(); coordinateDirection++)
        {
            stkMeshAdapter.getCoordinatesViewOf(Zoltan2::MESH_REGION, coordinates, stride, coordinateDirection);
            EXPECT_EQ(0.5, coordinates[0]);
            EXPECT_EQ(3, stride);
        }
    }
}

TEST(LoadBalance, DISABLED_createGraphEdgesUsingNodeConnectivity)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);
    int me;
    MPI_Comm_rank(communicator, &me);

    Options options = getOptionsForTest("generated:1x1x4");

    if(numProcs == 2)
    {
        stk::io::StkMeshIoBroker ioBroker(communicator);
        fillIoBroker(communicator, options.getMeshFileName(), ioBroker);

        stk::mesh::BulkData &stkMeshBulkData = ioBroker.bulk_data();
        verifyMeshPriorToRebalance(stkMeshBulkData);

        stk::mesh::impl::LocalIdMapper localIds(stkMeshBulkData, stk::topology::ELEM_RANK);

        Zoltan2ParallelGraph myGraph;

        size_t numElements = 0;
        std::vector<stk::balance::GraphEdge> graphEdges;
        stk::balance::GraphCreationSettings graphSettings;
        std::vector<int> adjacencyProcs;

        stk::mesh::Selector mySelector = stkMeshBulkData.mesh_meta_data().locally_owned_part();
        myGraph.fillZoltan2AdapterDataFromStkMesh(stkMeshBulkData,
                                                  graphSettings,
                                                  adjacencyProcs,
                                                  mySelector,
                                                  localIds);
        ASSERT_EQ(2u, numElements);
        std::vector<int> goldOffsets(numElements + 1);
        if(me == 0)
        {
            goldOffsets[0] = 0;
            goldOffsets[1] = 1;
            goldOffsets[2] = 3;
        }
        else
        {
            goldOffsets[0] = 0;
            goldOffsets[1] = 2;
            goldOffsets[2] = 3;
        }

        ASSERT_EQ(3u, graphEdges.size());
        std::vector<BalanceGlobalNumber> goldAdjacency(graphEdges.size());
        if(me == 0)
        {                         //    *-------*-------*-------*-------*
            goldAdjacency[0] = 2; //    |       |       |       |       |
            goldAdjacency[1] = 1; //    |   1   |   2   |   3   |   4   |
            goldAdjacency[2] = 3; //    |       |       |       |       |
                                  //    *-------*-------*-------*-------*
        }                         //
        else                      //        on proc 0:                     on proc 1:
        {                         //         *-------*-------* - - - *           * - - - *-------*-------*
            goldAdjacency[0] = 4; //         |       |       |       !           !       |       |       |
            goldAdjacency[1] = 2; //         |   1   |   2   |   3   !           !   2   |   3   |   4   |
            goldAdjacency[2] = 3; //         |       |       |       !           !       |       |       |
                                  // bucket  *-------*-------* - - - *           * - - - *-------*-------*
        }                         //   order:    1       2       3                   3       1       2

        std::vector<double> goldEdgeWeights(graphEdges.size(), 1.0);

        EXPECT_TRUE(goldOffsets == myGraph.get_offsets());
        EXPECT_TRUE(goldAdjacency == myGraph.get_adjacency());
        EXPECT_TRUE(goldEdgeWeights == myGraph.get_edge_weights());
    }
}

TEST(LoadBalance, zoltan2coloring)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);
    int me;
    MPI_Comm_rank(communicator, &me);

    Options options = getOptionsForTest("generated:3x3x3");

    if(numProcs == 2 || options.overRideTest())
    {
        stk::io::StkMeshIoBroker ioBroker(communicator);
        fillIoBroker(communicator, options.getMeshFileName(), ioBroker);

        stk::mesh::BulkData &stkMeshBulkData = ioBroker.bulk_data();

        Teuchos::ParameterList params("test params");

        if(options.debugZoltan())
        {
            params.set("timer_output_stream", "std::cout");
            params.set("debug_level", "verbose_detailed_status");
            params.set("debug_output_file", "kdd");
            params.set("debug_procs", "all");
        }

        stk::mesh::impl::LocalIdMapper localIds(stkMeshBulkData, stk::topology::ELEM_RANK);

        ColoringSettings coloringSettings;
        Zoltan2ParallelGraph zoltan2Graph;
        std::vector<int> adjacencyProcs;
        zoltan2Graph.fillZoltan2AdapterDataFromStkMesh(stkMeshBulkData,
                                                       coloringSettings,
                                                       adjacencyProcs,
                                                       stkMeshBulkData.mesh_meta_data().locally_owned_part(),
                                                       localIds);

        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts(stkMeshBulkData, counts);
        zoltan2Graph.set_num_global_elements(counts[stk::topology::ELEM_RANK]);

        zoltan2Graph.set_spatial_dim(stkMeshBulkData.mesh_meta_data().spatial_dimension());
        StkMeshZoltanAdapter stkMeshAdapter(zoltan2Graph);

        Zoltan2::ColoringProblem<StkMeshZoltanAdapter> problem(&stkMeshAdapter, &params);
        std::srand(stkMeshBulkData.parallel_rank()); // KHP: Temporary until an API is added to Zoltan2 for random seeds.
        problem.solve();

        if(options.debugZoltan())
        {
            problem.printTimers();
        }

        Zoltan2::ColoringSolution<StkMeshZoltanAdapter> *soln = problem.getSolution();
        size_t checkLength = soln->getColorsSize();
        int* checkColoring = soln->getColors();

        EXPECT_EQ(zoltan2Graph.get_vertex_ids().size(), checkLength);

        std::vector<int> coloring(checkColoring, checkColoring + checkLength);
        const std::string output_file_name = "output.exo";
        balance_utils::putFieldDataOnMesh(stkMeshBulkData, coloring);
        writeParFiles(ioBroker, output_file_name);

        if(me == 0)
        {
            balance_utils::clearFiles(output_file_name, numProcs);
        }

        if(!options.overRideTest())
        {
            std::vector<int>::iterator maxValueIter = std::max_element(coloring.begin(), coloring.end());
            int numberOfColors = *maxValueIter;
            int goldNumColors[] = {8, 4};
            EXPECT_EQ(goldNumColors[me], numberOfColors)<< " failed for processor " << me;
        }
    }
}

TEST(LoadBalance, ourColoring)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);
    int me;
    MPI_Comm_rank(communicator, &me);

    Options options = getOptionsForTest("generated:3x3x3");

    if(numProcs == 2)
    {
        stk::io::StkMeshIoBroker ioBroker(communicator);
        fillIoBroker(communicator, options.getMeshFileName(), ioBroker);

        stk::mesh::BulkData &stkMeshBulkData = ioBroker.bulk_data();

        stk::mesh::impl::LocalIdMapper localIds(stkMeshBulkData, stk::topology::ELEM_RANK);

        std::vector<int> adjacencyProcs;
        OurColoringSettings coloringSettings;

        Zoltan2ParallelGraph myGraph;
        myGraph.fillZoltan2AdapterDataFromStkMesh(stkMeshBulkData,
                                                  coloringSettings,
                                                  adjacencyProcs,
                                                  stkMeshBulkData.mesh_meta_data().locally_owned_part(),
                                                  localIds);

        std::vector<int> coloring(myGraph.get_vertex_ids().size(), std::numeric_limits<int>::max());

        for(size_t elementIndex = 0; elementIndex < coloring.size(); elementIndex++)
        {
            int numAdjElements = myGraph.get_offsets()[elementIndex + 1] - myGraph.get_offsets()[elementIndex];
            std::vector<int> adjColors(numAdjElements, -1);
            for(int i = 0; i < numAdjElements; i++)
            {
                int localId = myGraph.get_adjacency()[myGraph.get_offsets()[elementIndex] + i];

                adjColors[i] = coloring[localId];
            }
            std::sort(adjColors.begin(), adjColors.end());
            std::vector<int>::iterator endIter = std::unique(adjColors.begin(), adjColors.end());
            int numAdjColors = endIter - adjColors.begin();
            int findColor = numAdjColors;
            for(int i = 0; i < numAdjColors; i++)
            {
                if(adjColors[i] != i)
                {
                    findColor = i;
                    break;
                }
            }
            coloring[elementIndex] = findColor;
        }

        balance_utils::putFieldDataOnMesh(stkMeshBulkData, coloring);

        const std::string output_file_name = "output.exo";
        writeParFiles(ioBroker, output_file_name);

        if(me == 0 && options.deleteFiles())
        {
            balance_utils::clearFiles(output_file_name, numProcs);
        }

        std::vector<int>::iterator maxValueIter = std::max_element(coloring.begin(), coloring.end());
        int numberOfColors = *maxValueIter + 1;
        int goldNumColors[] = {8, 4};
        EXPECT_EQ(goldNumColors[me], numberOfColors);
    }
}

void checkMeshIsLoadBalanced(const stk::balance::BalanceSettings& balanceSettings, stk::mesh::BulkData& stkMeshBulkData)
{
    stk::mesh::impl::LocalIdMapper localIds(stkMeshBulkData, stk::topology::ELEM_RANK);

    Zoltan2ParallelGraph zoltan2Graph;
    std::vector<int> adjacencyProcs;
    zoltan2Graph.fillZoltan2AdapterDataFromStkMesh(stkMeshBulkData,
                                                   balanceSettings,
                                                   adjacencyProcs,
                                                   stkMeshBulkData.mesh_meta_data().locally_owned_part(),
                                                   localIds);

    struct LoadBalanceDiagnostics diagnostics;
    gatherLoadBalanceDiagnostics(zoltan2Graph.get_vertex_weights(),
                                 zoltan2Graph.get_edge_weights(),
                                 zoltan2Graph.get_adjacency(),
                                 stkMeshBulkData,
                                 stkMeshBulkData.parallel(), diagnostics);

    if(stkMeshBulkData.parallel_rank() == 0)
    {
      printLoadBalanceDiagnostics(diagnostics);
    }
}

TEST(LoadBalance, DISABLED_zoltan1decomposition)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);
    int procId;
    MPI_Comm_rank(communicator, &procId);

    Options options = getOptionsForTest("generated:3x3x3");

    if(numProcs == 2 || options.overRideTest())
    {
        stk::io::StkMeshIoBroker ioBroker(communicator);
        fillIoBroker(communicator, options.getMeshFileName(), ioBroker);

        stk::mesh::BulkData &stkMeshBulkData = ioBroker.bulk_data();

        stk::balance::GraphCreationSettingsWithCustomTolerances loadBalanceSettings;

        loadBalanceSettings.setToleranceForFaceSearch(options.getToleranceForFaceSearch());
        loadBalanceSettings.setToleranceForParticleSearch(options.getToleranceForParticleSearch());

        stk::balance::balanceStkMesh(loadBalanceSettings, stkMeshBulkData);

        checkMeshIsLoadBalanced(loadBalanceSettings, stkMeshBulkData);

        std::string output_file_name = "output.exo";
        std::vector<int> coloring(27, 0);
        balance_utils::putFieldDataOnMesh(stkMeshBulkData, coloring);
        writeParFiles(ioBroker, output_file_name);
    }
}

TEST(LoadBalance, MxN_decomposition)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);
    int procId;
    MPI_Comm_rank(communicator, &procId);

    Options options = getOptionsForTest("generated:3x3x3");

    if(numProcs == 2 || options.overRideTest())
    {
        stk::io::StkMeshIoBroker ioBroker(communicator);
        fillIoBroker(communicator, options.getMeshFileName(), ioBroker);

        stk::mesh::BulkData &stkMeshBulkData = ioBroker.bulk_data();

        stk::balance::GraphCreationSettingsWithCustomTolerances loadBalanceSettings;

        loadBalanceSettings.setToleranceForFaceSearch(options.getToleranceForFaceSearch());
        loadBalanceSettings.setToleranceForParticleSearch(options.getToleranceForParticleSearch());

        unsigned num_procs_decomp = static_cast<unsigned>(options.numSubdomains());
        stk::mesh::EntityProcVec decomp;
        std::vector<stk::mesh::Selector> selectors = {stkMeshBulkData.mesh_meta_data().locally_owned_part()};
        stk::balance::internal::calculateGeometricOrGraphBasedDecomp(loadBalanceSettings, num_procs_decomp, decomp, stkMeshBulkData, selectors);

        std::string output_file_name = "output.exo";
        balance_utils::putEntityProcOnMeshField(stkMeshBulkData, decomp);
        writeParFiles(ioBroker, output_file_name);

        std::vector<unsigned> mappings(num_procs_decomp, 0);
        int procCounter = 0;
        for(unsigned i = 0; i < num_procs_decomp; i++)
        {
            mappings[i] = procCounter;
            procCounter++;
            if(procCounter >= numProcs)
            {
                procCounter = 0;
            }
        }

        stk::balance::internal::rebalance(stkMeshBulkData, mappings, decomp);

        const std::string fieldName2 = "Coloring";
        stk::mesh::FieldBase *subdomainField = stkMeshBulkData.mesh_meta_data().get_field(stk::topology::ELEMENT_RANK, fieldName2);
        for(size_t i = 0; i < num_procs_decomp; i++)
        {
            std::vector<stk::mesh::Entity> entities;
            stkMeshBulkData.modification_begin();

            if(mappings[i] == static_cast<unsigned>(procId))
            {
                const stk::mesh::BucketVector &buckets = stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK);
                for(size_t j = 0; j < buckets.size(); j++)
                {
                    stk::mesh::Bucket &bucket = *buckets[j];
                    if(bucket.owned())
                    {
                        double *bucketSubdomainData = static_cast<double*>(stk::mesh::field_data(*subdomainField, bucket));
                        for(size_t k = 0; k < bucket.size(); k++)
                        {
                            if(bucketSubdomainData[k] == static_cast<double>(i))
                            {
                                entities.push_back(bucket[k]);
                            }
                        }
                    }
                }
            }

            stk::mesh::PartVector partVector;
            std::string partNameForSubdomain = getSubdomainPartName(i);
            stk::mesh::Part& subdomain = stkMeshBulkData.mesh_meta_data().declare_part(partNameForSubdomain, stk::topology::ELEMENT_RANK);
            partVector.push_back(&subdomain);

            for(size_t j = 0; j < entities.size(); j++)
            {
                stkMeshBulkData.change_entity_parts(entities[j], partVector);
            }

            stkMeshBulkData.modification_end();
        }

        //////////////// write out the files

        std::string filename = "subdomain.exo";

        for(size_t i = 0; i < num_procs_decomp; i++)
        {
            if(mappings[i] == static_cast<unsigned>(procId))
            {
                std::string partNameForSubdomain = getSubdomainPartName(i);
                stk::mesh::MetaData &stkMeshMetaData = stkMeshBulkData.mesh_meta_data();
                stk::mesh::Part& subdomain = *stkMeshMetaData.get_part(partNameForSubdomain);

                stk::mesh::MetaData newMeta;
                stk::mesh::BulkData newBulkData(newMeta, MPI_COMM_SELF);
                stk::tools::copy_mesh(stkMeshBulkData, subdomain, newBulkData);

                if(!options.overRideTest())
                {
                    EXPECT_EQ(9u, stk::mesh::count_selected_entities(newMeta.universal_part(), newBulkData.buckets(stk::topology::ELEMENT_RANK)));
                }

                std::string localFilename = balance_utils::getFilename(filename, num_procs_decomp, i);
                stk::io::StkMeshIoBroker meshIO(MPI_COMM_SELF);
                meshIO.set_bulk_data(newBulkData);
                size_t index = meshIO.create_output_mesh(localFilename, stk::io::WRITE_RESULTS);
                meshIO.write_output_mesh(index);
            }
        }

        if(procId == 0)
        {
            if(options.deleteFiles())
            {
                balance_utils::clearFiles(output_file_name, numProcs);
                balance_utils::clearFiles(filename, num_procs_decomp);
            }
        }
    }
}

TEST(LoadBalance, findBoundaryNodesAndFaces)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);
    int me;
    MPI_Comm_rank(communicator, &me);

    Options options = getOptionsForTest("generated:3x3x3");

    if(numProcs == 1 || options.overRideTest())
    {
        stk::io::StkMeshIoBroker ioBroker(communicator);
        fillIoBroker(communicator, options.getMeshFileName(), ioBroker);
        stk::mesh::BulkData &stkMeshBulkData = ioBroker.bulk_data();

        stk::mesh::Part & skin_part = stkMeshBulkData.mesh_meta_data().declare_part("SkinPart", stk::topology::FACE_RANK);
        stk::mesh::PartVector add_parts(1, &skin_part);
        stk::mesh::skin_mesh(stkMeshBulkData, stkMeshBulkData.mesh_meta_data().locally_owned_part(), add_parts);

        const stk::mesh::BucketVector &nodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, skin_part);
        const stk::mesh::BucketVector &faceBuckets = stkMeshBulkData.get_buckets(stk::topology::FACE_RANK, skin_part);

        unsigned numNodes = 0;
        for(size_t i = 0; i < nodeBuckets.size(); i++)
        {
            numNodes += nodeBuckets[i]->size();
        }

        unsigned numFaces = 0;
        for(size_t i = 0; i < faceBuckets.size(); i++)
        {
            numFaces += faceBuckets[i]->size();
        }

        std::vector<size_t> counts;
        stk::mesh::count_entities(stkMeshBulkData.mesh_meta_data().locally_owned_part(), stkMeshBulkData, counts);

        if(!options.overRideTest())
        {
            EXPECT_EQ(56u, numNodes);
            EXPECT_EQ(54u, numFaces);
            EXPECT_EQ(27u, counts[stk::topology::ELEMENT_RANK]);
        }
    }
}

TEST(LoadBalance, checkBBOnFace)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);

    Options options = getOptionsForTest("generated:1x1x1|sideset:x");

    if(numProcs == 1)
    {
        stk::io::StkMeshIoBroker ioBroker(communicator);
        fillIoBroker(communicator, options.getMeshFileName(), ioBroker);
        stk::mesh::BulkData &stkMeshBulkData = ioBroker.bulk_data();

        const stk::mesh::BucketVector &faceBuckets = stkMeshBulkData.buckets(stk::topology::FACE_RANK);
        ASSERT_EQ(1u, faceBuckets.size());
        ASSERT_EQ(1u, faceBuckets[0]->size());

        stk::balance::internal::StkBox goldFaceBoundingBox(-0.1, -0.1, -0.1, 0.1, 1.1, 1.1);
        stk::balance::internal::SearchBoxIdentProcs faceBoundingBoxesWithIdents;
        stk::mesh::Entity face = (*faceBuckets[0])[0];
        const double eps = 0.1;
        const stk::mesh::FieldBase * coord = stkMeshBulkData.mesh_meta_data().get_field(stk::topology::NODE_RANK,"coordinates");
        stk::balance::internal::addBoxForFace(stkMeshBulkData, face, eps, faceBoundingBoxesWithIdents, coord);

        ASSERT_EQ(1u, faceBoundingBoxesWithIdents.size());
        stk::balance::internal::StkBox resultBox = faceBoundingBoxesWithIdents[0].first;
        EXPECT_EQ(goldFaceBoundingBox, resultBox);

        unsigned numElements = stkMeshBulkData.num_elements(face);
        ASSERT_EQ(1u, numElements);

        const stk::mesh::Entity *element = stkMeshBulkData.begin_elements(face);
        stk::balance::internal::SearchIdentProc resultIdent = faceBoundingBoxesWithIdents[0].second;
        EXPECT_EQ(stkMeshBulkData.identifier(*element), resultIdent.id());

        EXPECT_EQ(stkMeshBulkData.parallel_rank(), resultIdent.proc());
    }
}

TEST(LoadBalance, doOneElementSearch)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);

    Options options = getOptionsForTest("generated:1x1x1");

    if(numProcs == 1)
    {
        stk::io::StkMeshIoBroker ioBroker(communicator);
        fillIoBroker(communicator, options.getMeshFileName(), ioBroker);
        stk::mesh::BulkData &stkMeshBulkData = ioBroker.bulk_data();

        stk::balance::internal::SearchBoxIdentProcs faceBoxes;
        const stk::mesh::FieldBase* coord = stkMeshBulkData.mesh_meta_data().get_field(stk::topology::NODE_RANK, "coordinates");
        stk::balance::GraphCreationSettings settings;
        settings.setToleranceForFaceSearch( 0.1 );
        stk::balance::internal::fillFaceBoxesWithIds(stkMeshBulkData, settings, coord, faceBoxes, stkMeshBulkData.mesh_meta_data().locally_owned_part());

        size_t goldNumFaceBoxes = 6u;
        EXPECT_EQ(goldNumFaceBoxes, faceBoxes.size());

        stk::balance::internal::SearchElemPairs searchResults;
        stk::search::coarse_search(faceBoxes, faceBoxes, stk::search::KDTREE, communicator, searchResults);

        size_t numNonSelfFaceInteractions = 24u;
        size_t numSelfInteractions = goldNumFaceBoxes;
        size_t totalNumInterations = numNonSelfFaceInteractions + numSelfInteractions;
        EXPECT_EQ(totalNumInterations, searchResults.size());

        stk::balance::internal::SearchElemPairs::iterator iter = std::unique(searchResults.begin(), searchResults.end());
        searchResults.resize(iter - searchResults.begin());
        size_t numUniqueInteractions = 1;
        EXPECT_EQ(numUniqueInteractions, searchResults.size());

        std::vector<stk::balance::GraphEdge> graphEdges;
        for(size_t i = 0; i < searchResults.size(); i++)
        {
            stk::mesh::Entity element1;
            element1 = static_cast<stk::mesh::Entity::Entity_t>(searchResults[i].first.id());
                    stk::mesh::Entity element2;
            element2 = static_cast<stk::mesh::Entity::Entity_t>(searchResults[i].second.id());
                    if(element1 != element2)
            {
                double edgeWeight = 1.0;
                const int owningProcElement2 = 0;
                graphEdges.push_back(stk::balance::GraphEdge(element1, stkMeshBulkData.identifier(element2), owningProcElement2, edgeWeight));
            }
        }
        EXPECT_TRUE(graphEdges.empty());
    }
}

template <typename GlobalIds>
void gatherLoadBalanceDiagnostics(const std::vector<double> &vertexWeights, const std::vector<double> &edgeWeights, const std::vector<GlobalIds> &adjacency, stk::mesh::BulkData &stkMeshBulkData, MPI_Comm communicator, struct LoadBalanceDiagnostics &diagnostics)
{

    double totalVertexWeight = 0;
    for(size_t i = 0; i < vertexWeights.size(); i++)
    {
        totalVertexWeight += vertexWeights[i];
    }

    double totalEdgeWeight = 0;
    for(size_t i = 0; i < adjacency.size(); i++)
    {
        stk::mesh::EntityKey entityKeyElement(stk::topology::ELEMENT_RANK, adjacency[i]);
        stk::mesh::Entity element = stkMeshBulkData.get_entity(entityKeyElement);
        if(!stkMeshBulkData.is_valid(element) || !stkMeshBulkData.bucket(element).owned())
        {
            totalEdgeWeight += edgeWeights[i];
        }
    }

    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);
    diagnostics.sumOfVertexWeightsPerProc.clear();
    diagnostics.sumOfVertexWeightsPerProc.resize(numProcs, 0);
    diagnostics.sumOfCutEdgeWeightsPerProc.clear();
    diagnostics.sumOfCutEdgeWeightsPerProc.resize(numProcs, 0);
    std::vector<int> counts(numProcs, 1);
    std::vector<int> displs(numProcs, 0);
    for(size_t i = 0; i < displs.size(); i++)
    {
        displs[i] = i;
    }

    int numElementsThisProc = vertexWeights.size();
    diagnostics.numElementsPerProc.clear();
    diagnostics.numElementsPerProc.resize(numProcs, 0);
    MPI_Allgatherv(&totalVertexWeight, 1, MPI_DOUBLE, &diagnostics.sumOfVertexWeightsPerProc[0], &counts[0], &displs[0], MPI_DOUBLE, communicator);
    MPI_Allgatherv(&totalEdgeWeight, 1, MPI_DOUBLE, &diagnostics.sumOfCutEdgeWeightsPerProc[0], &counts[0], &displs[0], MPI_DOUBLE, communicator);
    MPI_Allgatherv(&numElementsThisProc, 1, MPI_INT, &diagnostics.numElementsPerProc[0], &counts[0], &displs[0], MPI_INT, communicator);
}

void printLoadBalanceDiagnostics(const struct LoadBalanceDiagnostics &loadBalanceDiagnostics)
{
    for(size_t i = 0; i < loadBalanceDiagnostics.sumOfVertexWeightsPerProc.size(); i++)
    {
        std::cerr << "sumOfVertexWeightsPerProc(" << i + 1 << ") = " << loadBalanceDiagnostics.sumOfVertexWeightsPerProc[i] << ";" << std::endl;
    }
    for(size_t i = 0; i < loadBalanceDiagnostics.sumOfCutEdgeWeightsPerProc.size(); i++)
    {
        std::cerr << "sumOfCutEdgeWeightsPerProc(" << i + 1 << ") = " << loadBalanceDiagnostics.sumOfCutEdgeWeightsPerProc[i] << ";" << std::endl;
    }
    for(size_t i = 0; i < loadBalanceDiagnostics.numElementsPerProc.size(); i++)
    {
        std::cerr << "numElementsPerProc(" << i + 1 << ") = " << loadBalanceDiagnostics.numElementsPerProc[i] << ";" << std::endl;
    }
}

TEST(LoadBalance, doSearch)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);
    int me;
    MPI_Comm_rank(communicator, &me);

    Options options = getOptionsForTest("generated:3x3x3");

    if(numProcs == 2 || options.overRideTest())
    {
        stk::io::StkMeshIoBroker ioBroker(communicator);
        fillIoBroker(communicator, options.getMeshFileName(), ioBroker);
        stk::mesh::BulkData &stkMeshBulkData = ioBroker.bulk_data();

        stk::balance::internal::SearchBoxIdentProcs faceBoxes;
        const stk::mesh::FieldBase* coord = stkMeshBulkData.mesh_meta_data().get_field(stk::topology::NODE_RANK, "coordinates");
        stk::balance::GraphCreationSettings settings;
        settings.setToleranceForFaceSearch( 0.1 );
        stk::balance::internal::fillFaceBoxesWithIds(stkMeshBulkData, settings, coord, faceBoxes, stkMeshBulkData.mesh_meta_data().locally_owned_part());

        std::vector<stk::balance::internal::StkBox> faceItems(faceBoxes.size());
        for(size_t i = 0; i < faceBoxes.size(); i++)
        {
            faceItems[i] = faceBoxes[i].first;
        }
        std::string file = "faceBoxes.exo";
        std::string filename = balance_utils::getFilename(file, numProcs, me);
        writeExodusFileUsingBoxes(faceItems, filename);

        if(me == 0 && options.deleteFiles())
        {
            balance_utils::clearFiles(file, numProcs);
        }

        stk::balance::internal::SearchElemPairs searchResults;
        stk::search::coarse_search(faceBoxes, faceBoxes, stk::search::KDTREE, communicator, searchResults);

        stk::balance::internal::SearchElemPairs::iterator iter = std::unique(searchResults.begin(), searchResults.end());
        searchResults.resize(iter - searchResults.begin());

        for(size_t i = 0; i < searchResults.size(); i++)
        {
            stk::balance::internal::SearchIdentProc item1 = searchResults[i].first;
            stk::balance::internal::SearchIdentProc item2 = searchResults[i].second;
            std::ostringstream os;
            os << "For interaction " << i << " on processor " << me << " element " << item1.id() << " on proc " << item1.proc() <<
                    " iteracts with element " << item2.id() << " from proc " << item2.proc() << std::endl;
            if(me == 0 && options.overRideTest())
            {
                std::cerr << os.str();
            }

        }
    }
}

void create2DisconnectedHex8sStkMesh(stk::mesh::MetaData &meta, stk::mesh::BulkData &mesh);

TEST(LoadBalance, testGraphCreationUsingSearchForContact)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);

    Options options = getOptionsForTest("twoDis.exo");

    if(numProcs == 1)
    {
        const size_t spatialDim = 3;
        stk::mesh::MetaData meta(spatialDim);
        stk::mesh::BulkData stkMeshBulkData(meta, MPI_COMM_WORLD);
        create2DisconnectedHex8sStkMesh(meta, stkMeshBulkData);

        std::vector<size_t> counts;
        stk::mesh::count_entities(meta.universal_part(), stkMeshBulkData, counts);
        EXPECT_EQ(2u, counts[stk::topology::ELEMENT_RANK]);
        EXPECT_EQ(16u, counts[stk::topology::NODE_RANK]);

        if(options.overRideTest())
        {
            stk::io::StkMeshIoBroker meshIO(MPI_COMM_WORLD);
            meshIO.set_bulk_data(stkMeshBulkData);
            unsigned index = meshIO.create_output_mesh("twoDis.exo", stk::io::WRITE_RESULTS);
            meshIO.write_output_mesh(index);
        }

        std::vector<stk::balance::GraphEdge> graphEdges;
        stk::balance::GraphCreationSettingsWithCustomTolerances loadBalanceSettings;

        loadBalanceSettings.setToleranceForFaceSearch(options.getToleranceForFaceSearch());
        loadBalanceSettings.setToleranceForParticleSearch(options.getToleranceForParticleSearch());

        size_t numElements = stk::mesh::count_selected_entities(stkMeshBulkData.mesh_meta_data().locally_owned_part(),
                                                                stkMeshBulkData.buckets(stk::topology::ELEM_RANK));

        std::vector<double> vertexWeights(numElements, 1);
        stk::balance::internal::addGraphEdgesUsingBBSearch(stkMeshBulkData, loadBalanceSettings, graphEdges, meta.locally_owned_part());

        unsigned numEdgesCreated = 2;
        EXPECT_EQ(numEdgesCreated, graphEdges.size());

        if(options.overRideTest())
        {
            for(size_t i = 0; i < graphEdges.size(); i++)
            {
                std::cerr << "Edge " << i << " is(" << stkMeshBulkData.identifier(graphEdges[i].vertex1()) << ", " << graphEdges[i].vertex2() << ")" << std::endl;
            }

        }
    }
}

TEST(LoadBalance, createNewMeshFromPart_2DisconnectedHex8s)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);

    Options options = getOptionsForTest("twoDis.exo");

    if(numProcs == 1)
    {
        const size_t spatialDim = 3;
        stk::mesh::MetaData meta(spatialDim);
        stk::mesh::BulkData stkMeshBulkData(meta, MPI_COMM_WORLD);
        create2DisconnectedHex8sStkMesh(meta, stkMeshBulkData);
        const std::string blockName = "block_1";

        stk::mesh::Part &outputPart = *meta.get_part(blockName);

        std::vector<size_t> counts;
        stk::mesh::count_entities(meta.universal_part(), stkMeshBulkData, counts);
        EXPECT_EQ(2u, counts[stk::topology::ELEMENT_RANK]);
        EXPECT_EQ(16u, counts[stk::topology::NODE_RANK]);

        stk::mesh::MetaData newMeta;
        stk::mesh::BulkData newBulkData(newMeta, MPI_COMM_SELF);
        stk::tools::copy_mesh(stkMeshBulkData, outputPart, newBulkData);

        stk::mesh::count_entities(newMeta.universal_part(), newBulkData, counts);
        EXPECT_EQ(1u, counts[stk::topology::ELEMENT_RANK]);
        EXPECT_EQ(8u, counts[stk::topology::NODE_RANK]);
    }
}

TEST(LoadBalance, createNewMeshFromPart_writeFilesTillWeDontNeedThisTestAsADriverCodeAnymore)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);

    Options options = getOptionsForTest("generated:1x1x1");

    if(numProcs == 1 && options.overRideTest())
    {
        stk::io::StkMeshIoBroker ioBroker(communicator);
        fillIoBroker(communicator, options.getMeshFileName(), ioBroker);
        stk::mesh::BulkData &stkMeshBulkData = ioBroker.bulk_data();
        stk::mesh::MetaData &meta = stkMeshBulkData.mesh_meta_data();
        const std::string blockName = stk::unit_test_util::get_option("-b", "block_1");

        stk::mesh::Part &outputPart = *meta.get_part(blockName);

        stk::mesh::MetaData newMeta;
        stk::mesh::BulkData newBulkData(newMeta, MPI_COMM_SELF);
        stk::tools::copy_mesh(stkMeshBulkData, outputPart, newBulkData);

        std::string filename("block1.exo");
        stk::io::StkMeshIoBroker meshIO(MPI_COMM_SELF);
        meshIO.set_bulk_data(newBulkData);
        size_t index = meshIO.create_output_mesh(filename, stk::io::WRITE_RESULTS);
        meshIO.write_output_mesh(index);
    }
}

void create2ParticlesStkMesh(stk::mesh::MetaData &meta, stk::mesh::BulkData &mesh);

TEST(LoadBalance, testGraphCreationUsingSearchWithParticles)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);

    Options options = getOptionsForTest("twoDis.exo");

    if(numProcs == 1)
    {
        const size_t spatialDim = 3;
        stk::mesh::MetaData meta(spatialDim);
        stk::mesh::BulkData stkMeshBulkData(meta, MPI_COMM_WORLD);
        create2ParticlesStkMesh(meta, stkMeshBulkData);

        std::vector<size_t> counts;
        stk::mesh::count_entities(meta.universal_part(), stkMeshBulkData, counts);
        EXPECT_EQ(2u, counts[stk::topology::ELEMENT_RANK]);
        EXPECT_EQ(2u, counts[stk::topology::NODE_RANK]);

        if(options.overRideTest())
        {
            stk::io::StkMeshIoBroker meshIO(MPI_COMM_WORLD);
            meshIO.set_bulk_data(stkMeshBulkData);
            unsigned index = meshIO.create_output_mesh("twoDis.exo", stk::io::WRITE_RESULTS);
            meshIO.write_output_mesh(index);
        }

        std::vector<stk::balance::GraphEdge> graphEdges;
        stk::balance::GraphCreationSettingsWithCustomTolerances loadBalanceSettings;

        loadBalanceSettings.setToleranceForParticleSearch(2.01);

        size_t numElements = stk::mesh::count_selected_entities(stkMeshBulkData.mesh_meta_data().locally_owned_part(),
                                                                        stkMeshBulkData.buckets(stk::topology::ELEM_RANK));
        std::vector<double> vertexWeights(numElements, 1);
        stk::balance::internal::addGraphEdgesUsingBBSearch(stkMeshBulkData, loadBalanceSettings, graphEdges, meta.locally_owned_part());

        unsigned numEdgesCreated = 2;
        EXPECT_EQ(numEdgesCreated, graphEdges.size());

        if(options.overRideTest())
        {
            for(size_t i = 0; i < graphEdges.size(); i++)
            {
                std::cerr << "Edge " << i << " is(" << stkMeshBulkData.identifier(graphEdges[i].vertex1()) << ", " << graphEdges[i].vertex2() << ")" << std::endl;
            }

        }
    }
}

void create2Particles2HexStkMesh(stk::mesh::MetaData &meta, stk::mesh::BulkData &mesh);

TEST(LoadBalance, testGraphCreationUsingSearchWithParticlesAndSkin)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);

    Options options = getOptionsForTest("twoDis.exo");

    if(numProcs == 1)
    {
        const size_t spatialDim = 3;
        stk::mesh::MetaData meta(spatialDim);
        stk::mesh::BulkData stkMeshBulkData(meta, MPI_COMM_WORLD);
        create2Particles2HexStkMesh(meta, stkMeshBulkData);

        std::vector<size_t> counts;
        stk::mesh::count_entities(meta.universal_part(), stkMeshBulkData, counts);
        EXPECT_EQ(4u, counts[stk::topology::ELEMENT_RANK]);
        EXPECT_EQ(18u, counts[stk::topology::NODE_RANK]);

        if(options.overRideTest())
        {
            stk::io::StkMeshIoBroker meshIO(MPI_COMM_WORLD);
            meshIO.set_bulk_data(stkMeshBulkData);
            unsigned index = meshIO.create_output_mesh("twoDis.exo", stk::io::WRITE_RESULTS);
            meshIO.write_output_mesh(index);
        }

        std::vector<stk::balance::GraphEdge> graphEdges;
        stk::balance::GraphCreationSettingsWithCustomTolerances loadBalanceSettings;

        loadBalanceSettings.setToleranceForFaceSearch(0.21);
        loadBalanceSettings.setToleranceForParticleSearch(2.01);

        size_t numElements = stk::mesh::count_selected_entities(stkMeshBulkData.mesh_meta_data().locally_owned_part(),
                                                                        stkMeshBulkData.buckets(stk::topology::ELEM_RANK));

        std::vector<double> vertexWeights(numElements, 1);
        stk::balance::internal::addGraphEdgesUsingBBSearch(stkMeshBulkData, loadBalanceSettings, graphEdges, meta.locally_owned_part());

        unsigned numEdgesCreated = 12;
        EXPECT_EQ(numEdgesCreated, graphEdges.size());

        if(options.overRideTest())
        {
            for(size_t i = 0; i < graphEdges.size(); i++)
            {
                std::cerr << "Edge " << i << " is(" << stkMeshBulkData.identifier(graphEdges[i].vertex1()) << ", " << graphEdges[i].vertex2() << ")" << std::endl;
            }
        }
    }
}

TEST(LoadBalance, testSortingOfEdges)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);
    if(numProcs == 1)
    {
        std::vector<stk::balance::GraphEdge> graphEdges;
        stk::mesh::Entity vertex1Entity;
        unsigned globalIdOfEntity1 = 1;
        vertex1Entity = static_cast<stk::mesh::Entity::Entity_t>(globalIdOfEntity1);

        double largeWeight = 1000.0;
        double smallWeight = 1.0;
        int owningProcOfVertex2 = 0;
        unsigned globalIdOfVertex2 = 2;

        graphEdges.push_back(stk::balance::GraphEdge(vertex1Entity, globalIdOfVertex2, owningProcOfVertex2, largeWeight));
        graphEdges.push_back(stk::balance::GraphEdge(vertex1Entity, globalIdOfVertex2, owningProcOfVertex2, smallWeight));
        graphEdges.push_back(stk::balance::GraphEdge(vertex1Entity, globalIdOfVertex2, owningProcOfVertex2, largeWeight));

        std::sort(graphEdges.begin(), graphEdges.end());
        std::vector<stk::balance::GraphEdge>::iterator iter = std::unique(graphEdges.begin(), graphEdges.end());
        graphEdges.erase(iter, graphEdges.end());

        size_t goldNumberOfUniqueEdges = 1;
        EXPECT_EQ(goldNumberOfUniqueEdges, graphEdges.size());
        EXPECT_EQ(smallWeight, graphEdges[0].weight());
    }
}

namespace OurCartesianField
{
typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorField;
}

void putCoordinatesOnElement(stk::mesh::BulkData& bulk, stk::mesh::Entity element, double* nodalCoords, OurCartesianField::VectorField& coordField)
{
    unsigned spatialDim = bulk.mesh_meta_data().spatial_dimension();
    const stk::mesh::Entity* nodes = bulk.begin_nodes(element);
    unsigned numNodes = bulk.num_nodes(element);
    for(size_t j = 0; j < numNodes; j++)
    {
        double* coordsXyzForNode = stk::mesh::field_data(coordField, nodes[j]);

        coordsXyzForNode[0] = nodalCoords[j * spatialDim + 0];
        coordsXyzForNode[1] = nodalCoords[j * spatialDim + 1];
        coordsXyzForNode[2] = nodalCoords[j * spatialDim + 2];
    }
}

void create2DisconnectedHex8sStkMesh(stk::mesh::MetaData &meta, stk::mesh::BulkData &mesh)
{
    stk::mesh::Part& block1 = meta.declare_part_with_topology("block_1", stk::topology::HEX_8);
    stk::io::put_io_part_attribute(block1);
    stk::mesh::Part& block2 = meta.declare_part_with_topology("block_2", stk::topology::HEX_8);
    stk::io::put_io_part_attribute(block2);

    OurCartesianField::VectorField &coordField = meta.declare_field<OurCartesianField::VectorField>(stk::topology::NODE_RANK, "coordinates");

    stk::mesh::put_field_on_mesh(coordField, meta.universal_part(), meta.spatial_dimension(), nullptr);

    meta.commit();

    mesh.modification_begin();

    size_t offset = 8e9;
    stk::mesh::EntityId elem1GlobalId = 1 + offset;
    stk::mesh::EntityIdVector nodeIdsElement1 {1, 2, 3, 4, 5, 6, 7, 8};
    stk::mesh::Entity element1 = stk::mesh::declare_element(mesh, block1, elem1GlobalId, nodeIdsElement1);

    stk::mesh::EntityId elem2GlobalId = 2 + offset;
    stk::mesh::EntityIdVector nodeIdsElement2 {9, 10, 11, 12, 13, 14, 15, 16};
    stk::mesh::Entity element2 = stk::mesh::declare_element(mesh, block2, elem2GlobalId, nodeIdsElement2);

    mesh.modification_end();

    const size_t numNodesPerElement = 8;
    const size_t spatialDim = 3;
    double coordsElement1[numNodesPerElement * spatialDim] =
    {
        0,0,0,
        1,0,0,
        1,1,0,
        0,1,0,
        0,0,1,
        1,0,1,
        1,1,1,
        0,1,1
    };
    putCoordinatesOnElement(mesh, element1, coordsElement1, coordField);

    double coordsElement2[numNodesPerElement * spatialDim] =
    {
        1,0,0,
        2,0,0,
        2,1,0,
        1,1,0,
        1,0,1,
        2,0,1,
        2,1,1,
        1,1,1
    };
    putCoordinatesOnElement(mesh, element2, coordsElement2, coordField);
}

void create2ParticlesStkMesh(stk::mesh::MetaData &meta, stk::mesh::BulkData &mesh)
{
    stk::mesh::Part& block1 = meta.declare_part_with_topology("block_1", stk::topology::PARTICLE);
    stk::io::put_io_part_attribute(block1);

    OurCartesianField::VectorField &coordField = meta.declare_field<OurCartesianField::VectorField>(stk::topology::NODE_RANK, "coordinates");

    stk::mesh::put_field_on_mesh(coordField, meta.universal_part(), meta.spatial_dimension(), nullptr);

    meta.commit();

    mesh.modification_begin();

    stk::mesh::EntityId elem1GlobalId = 1;
    stk::mesh::EntityIdVector nodeIdsElement1 {1};
    stk::mesh::Entity element1 = stk::mesh::declare_element(mesh, block1, elem1GlobalId, nodeIdsElement1);

    stk::mesh::EntityId elem2GlobalId = 2;
    stk::mesh::EntityIdVector nodeIdsElement2 {2};
    stk::mesh::Entity element2 = stk::mesh::declare_element(mesh, block1, elem2GlobalId, nodeIdsElement2);

    mesh.modification_end();

    double coordsElement1[] = {0, 0, 0};
    putCoordinatesOnElement(mesh, element1, coordsElement1, coordField);

    double coordsElement2[] = {1, 0, 0};
    putCoordinatesOnElement(mesh, element2, coordsElement2, coordField);
}

void create2Particles2HexStkMesh(stk::mesh::MetaData &meta, stk::mesh::BulkData &mesh)
{
    stk::mesh::Part& block1 = meta.declare_part_with_topology("block_1", stk::topology::HEX_8);
    stk::io::put_io_part_attribute(block1);

    stk::mesh::Part& block2 = meta.declare_part_with_topology("block_2", stk::topology::PARTICLE);
    stk::io::put_io_part_attribute(block2);

    OurCartesianField::VectorField &coordField = meta.declare_field<OurCartesianField::VectorField>(stk::topology::NODE_RANK, "coordinates");

    stk::mesh::put_field_on_mesh(coordField, meta.universal_part(), meta.spatial_dimension(), nullptr);

    meta.commit();

    mesh.modification_begin();

    stk::mesh::EntityId elem1GlobalId = 1;
    stk::mesh::EntityIdVector nodeIdsElement1 {1, 2, 3, 4, 5, 6, 7, 8};
    stk::mesh::Entity element1 = stk::mesh::declare_element(mesh, block1, elem1GlobalId, nodeIdsElement1);

    stk::mesh::EntityId elem2GlobalId = 2;
    stk::mesh::EntityIdVector nodeIdsElement2 {9, 10, 11, 12, 13, 14, 15, 16};
    stk::mesh::Entity element2 = stk::mesh::declare_element(mesh, block1, elem2GlobalId, nodeIdsElement2);

    stk::mesh::EntityId elem3GlobalId = 3;
    stk::mesh::EntityIdVector nodeIdsElement3 {17};
    stk::mesh::Entity element3 = stk::mesh::declare_element(mesh, block2, elem3GlobalId, nodeIdsElement3);

    stk::mesh::EntityId elem4GlobalId = 4;
    stk::mesh::EntityIdVector nodeIdsElement4 {18};
    stk::mesh::Entity element4 = stk::mesh::declare_element(mesh, block2, elem4GlobalId, nodeIdsElement4);

    mesh.modification_end();

    const size_t numNodesPerElement = 8;
    const size_t spatialDim = 3;
    double coordsElement1[numNodesPerElement * spatialDim] =
    {
        0,0,0,
        1,0,0,
        1,1,0,
        0,1,0,
        0,0,1,
        1,0,1,
        1,1,1,
        0,1,1
    };
    putCoordinatesOnElement(mesh, element1, coordsElement1, coordField);

    const double translation = 0.001;
    double coordsElement2[numNodesPerElement * spatialDim] =
    {
        1+translation,0,0,
        2+translation,0,0,
        2+translation,1,0,
        1+translation,1,0,
        1+translation,0,1,
        2+translation,0,1,
        2+translation,1,1,
        1+translation,1,1
    };
    putCoordinatesOnElement(mesh, element2, coordsElement2, coordField);

    double coordsElement3[] = {0.95, 1.5, 0};
    putCoordinatesOnElement(mesh, element3, coordsElement3, coordField);

    double coordsElement4[] = {1.05, 1.5, 0};
    putCoordinatesOnElement(mesh, element4, coordsElement4, coordField);
}

void writeParFiles(stk::io::StkMeshIoBroker &ioBroker, const std::string &output_file_name)
{
    size_t index = ioBroker.create_output_mesh(output_file_name, stk::io::WRITE_RESULTS);

    ioBroker.set_active_mesh(index);
    ioBroker.write_output_mesh(index);

    double time = 0.5;

    const std::string fieldName1 = "DomainProc";
    stk::mesh::Field<double> &field1 = *ioBroker.meta_data().get_field<stk::mesh::Field<double> >(stk::topology::ELEMENT_RANK, fieldName1);
    ioBroker.add_field(index, field1);

    const std::string fieldName2 = "Coloring";
    stk::mesh::Field<double> &field2 = *ioBroker.meta_data().get_field<stk::mesh::Field<double> >(stk::topology::ELEMENT_RANK, fieldName2);
    ioBroker.add_field(index, field2);

    ioBroker.begin_output_step(index, time);
    ioBroker.write_defined_output_fields(index);
    ioBroker.end_output_step(index);
}

void fillIoBroker(MPI_Comm communicator, const std::string &generatedMeshSpec, stk::io::StkMeshIoBroker &ioBroker)
{
    std::string doDecomp = stk::unit_test_util::get_option("-decomp", "no");

    if ( doDecomp != "no")
    {
        //ioBroker.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
        //ioBroker.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RANDOM"));
        ioBroker.property_add(Ioss::Property("DECOMPOSITION_METHOD", "LINEAR"));
    }

    std::string useLargeInt = stk::unit_test_util::get_option("-lint", "no");
    if(useLargeInt != "no")
    {
        ioBroker.property_add(Ioss::Property("INTEGER_SIZE_API", 8));
    }

    ioBroker.add_mesh_database(generatedMeshSpec, stk::io::READ_MESH);
    ioBroker.create_input_mesh();

    {
        const std::string fieldName = "DomainProc";
        stk::mesh::Field<double> &field =
                ioBroker.meta_data().declare_field<stk::mesh::Field<double> >(stk::topology::ELEMENT_RANK,
                        fieldName, 1);
        stk::mesh::put_field_on_mesh(field, ioBroker.meta_data().universal_part(), nullptr);
    }

    {
        const std::string fieldName = "Coloring";
        stk::mesh::Field<double> &field =
                ioBroker.meta_data().declare_field<stk::mesh::Field<double> >(stk::topology::ELEMENT_RANK,
                        fieldName, 1);
        stk::mesh::put_field_on_mesh(field, ioBroker.meta_data().universal_part(), nullptr);
    }

    ioBroker.populate_bulk_data();
}

void setUpDefaultColoring(const stk::mesh::BulkData &stkMeshBulkData, std::vector<int>& coloring)
{
    coloring.clear();
    const stk::mesh::BucketVector &buckets = stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK);
    int elementCounter = 0;
    for(size_t i = 0; i < buckets.size(); i++)
    {
        const stk::mesh::Bucket &bucket = *buckets[i];
        if(bucket.owned())
        {
            for(size_t j = 0; j < bucket.size(); j++)
            {
                coloring.push_back(elementCounter);
                elementCounter++;
            }
        }
    }
}

void verifyMeshPriorToRebalance(stk::mesh::BulkData &stkMeshBulkData)
{
    stk::mesh::Selector locallyOwnedSelector = stkMeshBulkData.mesh_meta_data().locally_owned_part();
    const stk::mesh::BucketVector &buckets = stkMeshBulkData.get_buckets(stk::topology::ELEMENT_RANK, locallyOwnedSelector);
    size_t goldNumBuckets = 1u;
    EXPECT_EQ(goldNumBuckets, buckets.size());

    size_t goldNumLocalElementsInBucket = 2u;
    stk::mesh::Bucket &bucket1 = *buckets[0];
    EXPECT_EQ(goldNumLocalElementsInBucket, bucket1.size());

    stk::mesh::EntityId goldElementId = 1 + 2 * stkMeshBulkData.parallel_rank();
    EXPECT_EQ(goldElementId, stkMeshBulkData.identifier(bucket1[0]));
    EXPECT_EQ(goldElementId+1, stkMeshBulkData.identifier(bucket1[1]));
}

void verifyMeshAfterRebalance(stk::mesh::BulkData &stkMeshBulkData)
{
    stk::mesh::Selector locallyOwnedSelector = stkMeshBulkData.mesh_meta_data().locally_owned_part();
    const stk::mesh::BucketVector &buckets = stkMeshBulkData.get_buckets(stk::topology::ELEMENT_RANK, locallyOwnedSelector);
    size_t goldNumBuckets = 1u;
    EXPECT_EQ(goldNumBuckets, buckets.size());

    size_t goldNumLocalElementsInOnlyBucket = 2u;
    stk::mesh::Bucket &onlyBucket = *buckets[0];
    EXPECT_EQ(goldNumLocalElementsInOnlyBucket, onlyBucket.size());

    stk::mesh::EntityId goldElementId[2][2] = { {2, 3}, {1, 4}};
    EXPECT_EQ(goldElementId[stkMeshBulkData.parallel_rank()][0], stkMeshBulkData.identifier(onlyBucket[0]));
    EXPECT_EQ(goldElementId[stkMeshBulkData.parallel_rank()][1], stkMeshBulkData.identifier(onlyBucket[1]));
}

void createMockElementDecompositon(const int procId, stk::mesh::EntityProcVec &mockDecomposition, const stk::mesh::EntityVector& entities)
{
    int destProc = -1;
    if(procId == 0)
    {
        destProc = 1;
    }
    else
    {
        destProc = 0;
    }

    mockDecomposition.clear();
    mockDecomposition.resize(entities.size());
    for(size_t i=0;entities.size();++i)
    {
        mockDecomposition[i] = std::make_pair(entities[i], procId);
    }
    mockDecomposition[0].second = destProc;
}

template<typename GlobalId, typename LocalNumber>
void writeDotFile(const std::string &fileName, const std::vector<GlobalId>& globalIds, const std::vector<LocalNumber> &offsets, const std::vector<GlobalId>& adjacency)
{
    size_t numElem = offsets.size() - 1;

    std::ofstream out(fileName.c_str());
    out << "graph a {\n";
    for(size_t i = 0; i < numElem; i++)
    {
        out << globalIds[i] << "  ";
        out << std::endl;
        for(int j = offsets[i]; j < offsets[i + 1]; j++)
        {
            out << globalIds[i] << " -- " << adjacency[j] << "  ";
            out << "\n";
        }
    }
    out << "}\n";
    out.close();
}

std::string getSubdomainPartName(int subdomainId)
{
    std::ostringstream out;
    out << "subdomain_" << subdomainId;
    return out.str();
}

// 10. 0+1 = 1
// 100. 1+1 = 2

TEST(LoadBalance, checkWidth)
{
    EXPECT_EQ(1, balance_utils::getWidth(9));
    EXPECT_EQ(1, balance_utils::getWidth(10));
    EXPECT_EQ(2, balance_utils::getWidth(11));
    EXPECT_EQ(2, balance_utils::getWidth(99));
    EXPECT_EQ(2, balance_utils::getWidth(100));
    EXPECT_EQ(3, balance_utils::getWidth(101));
    EXPECT_EQ(3, balance_utils::getWidth(1000));
}

}

