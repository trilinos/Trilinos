#include <gtest/gtest.h>
#include <stk_topology/topology.hpp>
#include <Akri_DiagWriter.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <Akri_MOAB_TetRefiner.hpp>
#include <Akri_TriRefiner.hpp>

namespace krino {

static stk::topology get_refinement_topology(const stk::topology baseTopology)
{
  switch(baseTopology())
  {
    case stk::topology::TRIANGLE_3_2D:
      return stk::topology::TRIANGLE_6_2D;
    case stk::topology::TETRAHEDRON_4:
      return stk::topology::TETRAHEDRON_10;
    default:
        ThrowRuntimeError("Element topology not found in get_refinement_topology: " << baseTopology.name());
  }
}

void fill_permutations_and_permuted_case_ids(const stk::topology baseTopo, std::vector<unsigned> & permutations, std::vector<unsigned> & permutedCaseIds)
{
  const stk::topology topo = get_refinement_topology(baseTopo);
  const unsigned numEdges = topo.num_edges();
  const unsigned numNodes = topo.num_nodes();

  std::cout << "Building permutation tables for " << topo.name() << " with " << numEdges << " edges." << std::endl;

  const unsigned numCases = 1<<numEdges;
  permutations.resize(numCases);
  permutedCaseIds.resize(numCases);

  const unsigned numBaseNodes = numNodes-numEdges;
  std::vector<unsigned> permutedNodeIds(numNodes);

  for (unsigned i=0; i<numCases; ++i)
  {
    unsigned lowestCaseId = numCases;
    unsigned lowestPermutation = topo.num_positive_permutations();
    for (unsigned iPerm=0; iPerm<topo.num_positive_permutations(); ++iPerm)
    {
      topo.permutation_node_ordinals(iPerm, permutedNodeIds.data());

      unsigned permutedCaseId = 0;
      for (unsigned iEdge=0; iEdge<numEdges; ++iEdge)
      {
        const int permutedEdge = permutedNodeIds[iEdge+numBaseNodes]-numBaseNodes;
        if (i & (1<<permutedEdge))
          permutedCaseId += 1<<iEdge;
      }

      if (permutedCaseId < lowestCaseId)
      {
        lowestCaseId = permutedCaseId;
        lowestPermutation = iPerm;
      }
    }

    permutations[i] = lowestPermutation;
    permutedCaseIds[i] = lowestCaseId;
  }

  std::cout << "  Permutations = ";
  for (auto perm : permutations) std::cout << perm << ",";
  std::cout << std::endl;

  std::cout << "  Permuted case ids = ";
  for (auto permutedCaseId : permutedCaseIds) std::cout << permutedCaseId << ",";
  std::cout << std::endl;
}

TEST(Refinement,permutationTables)
{
  stk::ParallelMachine comm{MPI_COMM_WORLD};
  if(stk::parallel_machine_size(comm) == 1)
  {
    std::vector<unsigned> permutations;
    std::vector<unsigned> permutedCaseIds;
    fill_permutations_and_permuted_case_ids(stk::topology::TRIANGLE_3_2D, permutations, permutedCaseIds);
    for (unsigned iCase=0; iCase<permutations.size(); ++iCase)
    {
      EXPECT_EQ(permutations[iCase], TriRefiner::determine_permutation_tri3(iCase));
      EXPECT_EQ(permutedCaseIds[iCase], TriRefiner::determine_permuted_case_id_tri3(iCase));
    }
    fill_permutations_and_permuted_case_ids(stk::topology::TETRAHEDRON_4, permutations, permutedCaseIds);
    for (unsigned iCase=0; iCase<permutations.size(); ++iCase)
    {
      EXPECT_EQ(permutations[iCase], moab::SimplexTemplateRefiner::determine_permutation_tet4(iCase));
      EXPECT_EQ(permutedCaseIds[iCase], moab::SimplexTemplateRefiner::determine_permuted_case_id_tet4(iCase));
    }
  }
}

}


