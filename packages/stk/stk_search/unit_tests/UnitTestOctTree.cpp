// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <gtest/gtest.h>
#include <stk_search/OctTreeOps.hpp>
#include <stk_search/OctTree.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc

#include <gtest/gtest.h>

namespace
{

void checkKeyAgainstOffset(const int key[], const int offset_key_1);
void getCutsForProcessorCount(unsigned numProcsLocal, const float * const weights, stk::OctTreeKey *cuts);

typedef stk::search::IdentProc<uint64_t, unsigned> Ident;
typedef stk::search::Point<double> Point;
typedef stk::search::Box<double> Box;

TEST(stk_search_oct_tree, checkCuts)
{
    typedef std::pair<Box, Ident> BoxIdent;
    typedef std::vector<BoxIdent> BoxVector;

    MPI_Comm comm = MPI_COMM_WORLD;
    int proc_id = stk::parallel_machine_rank(comm);
    int num_procs = stk::parallel_machine_size(comm);

    double offsetFromEdgeOfProcessorBoundary=0.1;
    double sizeOfDomainPerProcessor=1.0;
    double boxSize=0.8;
    ASSERT_TRUE(offsetFromEdgeOfProcessorBoundary<=sizeOfDomainPerProcessor-offsetFromEdgeOfProcessorBoundary);
    double min=offsetFromEdgeOfProcessorBoundary;
    double max=boxSize+offsetFromEdgeOfProcessorBoundary+1;

    if(num_procs >= 4)
    {
        Point min_corner, max_corner;

        min_corner[0] = offsetFromEdgeOfProcessorBoundary;
        min_corner[1] = offsetFromEdgeOfProcessorBoundary;
        min_corner[2] = offsetFromEdgeOfProcessorBoundary;
        max_corner[0] = boxSize+offsetFromEdgeOfProcessorBoundary;
        max_corner[1] = boxSize+offsetFromEdgeOfProcessorBoundary;
        max_corner[2] = boxSize+offsetFromEdgeOfProcessorBoundary;
        if (proc_id == 1)
        {
            min_corner[0] += sizeOfDomainPerProcessor;
            max_corner[0] += sizeOfDomainPerProcessor;
        }
        else if (proc_id == 3)
        {
            min_corner[1] += sizeOfDomainPerProcessor;
            max_corner[1] += sizeOfDomainPerProcessor;
        }
        else if (proc_id == 2)
        {
            min_corner[0] += sizeOfDomainPerProcessor;
            min_corner[1] += sizeOfDomainPerProcessor;
            max_corner[0] += sizeOfDomainPerProcessor;
            max_corner[1] += sizeOfDomainPerProcessor;
        }

        BoxVector local_domain;
        Ident domainBox(proc_id, proc_id);
        if ( proc_id < 4 )
        {
          local_domain.push_back(std::make_pair(Box(min_corner, max_corner), domainBox));
        }

        BoxVector local_range;
        local_range = local_domain;

        std::vector<float> global_box(6);
        stk::search::box_global_bounds(comm,
            local_domain.size(),
            &local_domain[0] ,
            local_range.size(),
            &local_range[0],
            &global_box[0]);

        float tolerance = 2*std::numeric_limits<float>::epsilon();
        EXPECT_NEAR(float(min), global_box[0], tolerance);
        EXPECT_NEAR(float(min), global_box[1], tolerance);
        EXPECT_NEAR(float(min), global_box[2], tolerance);
        EXPECT_NEAR(float(max), global_box[3], tolerance);
        EXPECT_NEAR(float(max), global_box[4], tolerance);
        EXPECT_NEAR(float(boxSize+offsetFromEdgeOfProcessorBoundary), global_box[5], tolerance);

        typedef std::map< stk::OctTreeKey, std::pair< std::list< BoxIdent >, std::list< BoxIdent > > > SearchTree ;

        SearchTree searchTree;

        bool local_violations = true;
        unsigned Dim = 3;
        stk::search::createSearchTree(&global_box[0], local_domain.size(), &local_domain[0],
                local_range.size(), &local_range[0], Dim, proc_id, local_violations,
                searchTree);

        const double tol = 0.001 ;
        std::vector< stk::OctTreeKey > cuts ;

        stk::search::oct_tree_partition( comm , searchTree , tol , cuts );

        const int tree_depth = 4;
        EXPECT_EQ(0u, stk::oct_tree_offset(tree_depth, cuts[0]));
        EXPECT_EQ(2u, stk::oct_tree_offset(tree_depth, cuts[1]));
        EXPECT_EQ(587u, stk::oct_tree_offset(tree_depth, cuts[2]));
        EXPECT_EQ(1172u, stk::oct_tree_offset(tree_depth, cuts[3]));

    }
    else
    {
        std::cerr << "WARNING: Test not setup for anything other than 4 processors; ran with " << num_procs << "." << std::endl;
    }
}

TEST(stk_search_oct_tree, testCalculationOfKeyUsingOffset)
{
    int procId=stk::parallel_machine_rank(MPI_COMM_WORLD);
    if ( procId == 0 )
    {
        int key_1[4] = { 1, 1, 1, 1 };
        int offset_key_1 = 4;
        checkKeyAgainstOffset(key_1, offset_key_1);
        int key_2[4] = { 8, 8, 8, 8 };
        int offset_key_2 = 4680;
        checkKeyAgainstOffset(key_2, offset_key_2);
        int key_3[4] = {1, 0, 0, 0 };
        int offset_key_3 = 1;
        checkKeyAgainstOffset(key_3, offset_key_3);
        int key_4[4] = { 1, 1, 1, 2 };
        int offset_key_4 = 5;
        checkKeyAgainstOffset(key_4, offset_key_4);
    }
}

TEST(stk_search_oct_tree, testPartitioningOfPhysicalTreeForVaryingNumberOfProcsAndWeights)
{
    int procId=stk::parallel_machine_rank(MPI_COMM_WORLD);
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

TEST(stk_search_oct_tree, stressTestPartitioningUpToOneMillionProcessors)
{
    unsigned depth = 4;
    unsigned tree_size = stk::oct_tree_size(depth);
    float * weights = new float[tree_size*2];

    int procId=stk::parallel_machine_rank(MPI_COMM_WORLD);
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

        float totalWeight = static_cast<float>(2*tree_size);
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

void getCutsForProcessorCount(unsigned numProcsLocal, const float * const weights, stk::OctTreeKey *cuts)
{
        unsigned depth = 4;
        stk::search::partition_oct_tree(numProcsLocal, depth, weights, numProcsLocal, cuts);
}

void setKeyUsingIntArray(const unsigned depth, const int key[], stk::OctTreeKey &goldKey)
{
    for (unsigned int i=0;i<depth;i++)
    {
        bool isLeadingBitIndexZero = ( key[i] == 0 );
        if ( isLeadingBitIndexZero ) break;
        goldKey.set_index(i+1, key[i]);
    }
}

void checkKeyAgainstOffset(const int key[], const int offset_key_1)
{
    unsigned depth = 4;
    stk::OctTreeKey octTreeKey;
    stk::search::calculate_key_using_offset(depth, offset_key_1, octTreeKey);

    stk::OctTreeKey goldKey;
    setKeyUsingIntArray(depth, key, goldKey);

    EXPECT_EQ( goldKey, octTreeKey);
}

} // end namespace

