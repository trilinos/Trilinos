// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   TensorBasisTagTests.cpp
    \brief  Tests of the DoF tag queries of tensor-product bases, which are answered from the component bases or from lazily
            built tables.
    \author Nathan V. Roberts
*/

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_DerivedBasisFamily.hpp"
#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"
#include "Intrepid2_TestUtils.hpp"

#include <map>
#include <thread>
#include <tuple>
#include <vector>

namespace
{
  using namespace Intrepid2;

  //! Checks that getDofCount(), getDofOrdinal(), getDofTag(), getAllDofTags(), and getAllDofOrdinal() agree with each other.
  template<class BasisPtrType>
  void testTagConsistency(const BasisPtrType &basis, Teuchos::FancyOStream &out, bool &success)
  {
    const int cardinality = basis->getCardinality();
    out << "Testing tags for " << basis->getName() << " (cardinality " << cardinality << ")\n";

    // tags from getDofTag(), copied, since each returned row is a view
    std::vector<std::array<ordinal_type,4>> tags(cardinality);
    for (int i=0; i<cardinality; i++)
    {
      auto tag = basis->getDofTag(i);
      for (int k=0; k<4; k++) tags[i][k] = tag(k);
    }

    // rows returned by getDofTag() are distinct: an earlier row is unchanged by a later query
    if (cardinality > 1)
    {
      auto first  = basis->getDofTag(0);
      auto second = basis->getDofTag(cardinality-1);
      for (int k=0; k<4; k++)
      {
        TEST_EQUALITY(first(k),  tags[0][k]);
        TEST_EQUALITY(second(k), tags[cardinality-1][k]);
      }
    }

    // getAllDofTags() agrees with getDofTag()
    auto allTags = basis->getAllDofTags();
    TEST_EQUALITY(int(allTags.extent(0)), cardinality);
    for (int i=0; i<cardinality; i++)
      for (int k=0; k<4; k++)
        TEST_EQUALITY(allTags(i,k), tags[i][k]);

    // subcell dof counts implied by the tags
    std::map<std::pair<int,int>,int> dofCounts;
    for (int i=0; i<cardinality; i++)
    {
      const auto key = std::make_pair(int(tags[i][0]), int(tags[i][1]));
      dofCounts[key]++;
      // the tag's dof count field agrees with the number of tags for the subcell (checked after the loop)
    }
    for (int i=0; i<cardinality; i++)
    {
      const auto key = std::make_pair(int(tags[i][0]), int(tags[i][1]));
      TEST_EQUALITY(int(tags[i][3]), dofCounts[key]);
      // getDofOrdinal() inverts getDofTag()
      TEST_EQUALITY(basis->getDofOrdinal(tags[i][0], tags[i][1], tags[i][2]), i);
    }
    // getDofCount() agrees with the tags, including subcells without dofs
    const int spaceDim = basis->getDomainDimension();
    auto allDofOrdinals = basis->getAllDofOrdinal();
    for (int d=0; d<=spaceDim; d++)
    {
      for (int s=0; s<int(allDofOrdinals.extent(1)); s++)
      {
        const auto entry = dofCounts.find({d,s});
        const int expected = (entry == dofCounts.end()) ? 0 : entry->second;
        TEST_EQUALITY(basis->getDofCount(d,s), expected);
      }
    }
  }

  //! Checks that concurrent tag queries from several host threads, on a basis whose tag tables have not been built, agree.
  template<class BasisPtrType>
  void testConcurrentTagQueries(const BasisPtrType &basis, Teuchos::FancyOStream &out, bool &success)
  {
    const int cardinality = basis->getCardinality();
    const int numThreads = 4;
    std::vector<std::vector<int>> results(numThreads, std::vector<int>(4 * cardinality + 1));
    std::vector<std::thread> threads;
    for (int t=0; t<numThreads; t++)
    {
      threads.emplace_back([&basis, &results, t, cardinality]()
      {
        // alternate between the scalar queries and the tag tables, so that both lazy paths race
        if (t % 2 == 0) results[t][4 * cardinality] = basis->getDofCount(0,0);
        auto allTags = basis->getAllDofTags();
        for (int i=0; i<cardinality; i++)
          for (int k=0; k<4; k++)
            results[t][4*i+k] = allTags(i,k);
        if (t % 2 == 1) results[t][4 * cardinality] = basis->getDofCount(0,0);
      });
    }
    for (auto &thread : threads) thread.join();
    for (int t=1; t<numThreads; t++)
    {
      TEST_COMPARE_ARRAYS(results[t], results[0]);
    }
  }

  using HierarchicalFamily = HierarchicalBasisFamily<Kokkos::DefaultHostExecutionSpace>;
  using NodalFamily        = DerivedNodalBasisFamily<Kokkos::DefaultHostExecutionSpace>;

  TEUCHOS_UNIT_TEST( TensorBasisTags, Consistency_Quad_Hex )
  {
    const std::vector<EFunctionSpace> spaces {FUNCTION_SPACE_HGRAD, FUNCTION_SPACE_HCURL, FUNCTION_SPACE_HDIV, FUNCTION_SPACE_HVOL};
    for (const auto &topo : {shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<> >()),
                             shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<> >())})
    {
      for (const auto fs : spaces)
      {
        for (int p=1; p<=3; p++)
        {
          testTagConsistency(getBasis<HierarchicalFamily>(topo, fs, p), out, success);
          testTagConsistency(getBasis<NodalFamily>       (topo, fs, p), out, success);
        }
      }
    }
  }

  TEUCHOS_UNIT_TEST( TensorBasisTags, Consistency_Hypercube )
  {
    for (int spaceDim=1; spaceDim<=5; spaceDim++)
    {
      for (int p=1; p<=2; p++)
      {
        testTagConsistency(getHypercubeBasis_HGRAD<HierarchicalFamily>(p, spaceDim), out, success);
        testTagConsistency(getHypercubeBasis_HVOL <HierarchicalFamily>(p, spaceDim), out, success);
      }
    }
  }

  TEUCHOS_UNIT_TEST( TensorBasisTags, ConcurrentQueries )
  {
#ifdef KOKKOS_ENABLE_ATOMICS_BYPASS
    // Kokkos atomics (including View reference counts) are not atomic in this build (Trilinos enables
    // Kokkos_ENABLE_ATOMICS_BYPASS by default when Serial is the only Kokkos backend), so Views, including those returned by
    // the tag queries, cannot be shared between host threads.
    out << "Skipping: Kokkos is built with KOKKOS_ENABLE_ATOMICS_BYPASS, so Views cannot be shared between host threads.\n";
#else
    const auto hex = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<> >());
    testConcurrentTagQueries(getBasis<HierarchicalFamily>(hex, FUNCTION_SPACE_HCURL, 3), out, success);
    testConcurrentTagQueries(getHypercubeBasis_HGRAD<HierarchicalFamily>(2, 4), out, success);
#endif
  }
} // namespace
