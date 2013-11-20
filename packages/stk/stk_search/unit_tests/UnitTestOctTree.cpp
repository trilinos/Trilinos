#include <gtest/gtest.h>
#include <mpi.h>
#include <stk_search/OctTreeOps.hpp>
#include <stk_search/OctTree.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

namespace
{

STKUNIT_UNIT_TEST(stk_search_not_boost, testCalculationOfKeyUsingOffset)
{
    int procId=-1;
    MPI_Comm_rank(MPI_COMM_WORLD, &procId);
    if ( procId == 0 )
    {
        unsigned depth = 4;
        unsigned offset = 4;
        stk::OctTreeKey key;
        stk::search::calculate_key_using_offset(depth, offset, key);
        std::cerr << "key = " << key << std::endl;
    }
}

void getCutsForProcessorCount(unsigned numProcsLocal, const float * const weights, stk::OctTreeKey *cuts)
{
        unsigned depth = 4;
        stk::search::partition_oct_tree(numProcsLocal, depth, weights, numProcsLocal, cuts);
}


STKUNIT_UNIT_TEST(stk_search_not_boost, testPartitioningOfPhysicalTreeForVaryingNumberOfProcsAndWeights)
{
    int procId=-1;
    MPI_Comm_rank(MPI_COMM_WORLD, &procId);

    if ( procId == 0 )
    {
        std::vector<unsigned> ordinals;
        unsigned ordinalOnBottomRightOfTree = 1;
        unsigned ordinalOnBottomLeftOfTree = 4680;
        ordinals.push_back(ordinalOnBottomRightOfTree);
        ordinals.push_back(ordinalOnBottomLeftOfTree);
        ordinals.push_back(2000);

        unsigned depth = 4;
        unsigned tree_size = stk::oct_tree_size(depth);
        float * weights = new float[tree_size*2];
        for (size_t ord=0;ord<ordinals.size();ord++)
        {
            for (unsigned i=0;i<tree_size*2;i++)
            {
                weights[i] = 0.0;
            }
            unsigned goldOrdinal=ordinals[ord];
            weights[2*goldOrdinal]=1;
            std::vector<unsigned> numProcs;
            numProcs.push_back(2);
            numProcs.push_back(100);
            numProcs.push_back(1000000);
            for (size_t i=0;i<numProcs.size();i++)
            {
                unsigned numProcsLocal = numProcs[i];

                stk::OctTreeKey *cuts = new stk::OctTreeKey[numProcsLocal];
                getCutsForProcessorCount(numProcsLocal, weights, cuts);
                unsigned procIdToTest=0;
                unsigned ordinalToTest = stk::oct_tree_offset(depth, cuts[procIdToTest+1]);
                unsigned addOffset = 1;
                if ( goldOrdinal == tree_size-1 )
                {
                    addOffset = 0;
                }
                EXPECT_EQ(goldOrdinal, ordinalToTest-addOffset) << "Failed for proc count = " << numProcsLocal << " with ordinal " << goldOrdinal << std::endl;
                delete [] cuts;
            }
        }
        delete [] weights;
    }
}

STKUNIT_UNIT_TEST(stk_search_not_boost, stressTestPartitioningUpToOneMillionProcessors)
{
    unsigned depth = 4;
    unsigned tree_size = stk::oct_tree_size(depth);
    float * weights = new float[tree_size*2];

    int procId=-1;
    MPI_Comm_rank(MPI_COMM_WORLD, &procId);
    std::vector<unsigned> numProcs;
    numProcs.push_back(2);
    numProcs.push_back(8);
    numProcs.push_back(29);
    numProcs.push_back(4096);
    numProcs.push_back(10111);
    numProcs.push_back(28023);
    numProcs.push_back(49023);
    numProcs.push_back(102321);
    numProcs.push_back(480321);
    numProcs.push_back(1023845);

    if ( procId == 0 )
    {
        for (unsigned i=0;i<tree_size*2;i++)
        {
            weights[i] = 1.0;
        }
        float weightPerNode = 2.0;

        float totalWeight = float(2*tree_size);
        for (size_t i=0;i<numProcs.size();i++)
        {
            unsigned numProcsLocal = numProcs[i];
            float targetWeight = totalWeight/numProcsLocal;
            unsigned firstTargetWeightOrdinal = targetWeight/weightPerNode + 1;

            unsigned cuts_length = numProcsLocal;
            stk::OctTreeKey *cuts = new stk::OctTreeKey[cuts_length];

            stk::search::partition_oct_tree(numProcsLocal, depth, weights, cuts_length, cuts);

            unsigned procIdToTest=0;
            unsigned ordinal = stk::oct_tree_offset(depth, cuts[procIdToTest+1]);
            EXPECT_EQ(firstTargetWeightOrdinal, ordinal) << "Failed for proc count = " << numProcsLocal << std::endl;

            delete [] cuts;
        }
    }

    delete [] weights;
}

} // end namespace

