// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#include "gtest/gtest.h"                // for EXPECT_TRUE, TEST, etc
#include "stk_util/parallel/Parallel.hpp"  // for parallel_machine_size, etc
#include <stk_mesh/base/Types.hpp>      // for PartVector, BucketVector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include <stk_mesh/base/Part.hpp>       // for Part
#include <stk_mesh/base/Selector.hpp>   // for Selector, operator|, etc
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names
#include "stk_mesh/base/SkinBoundary.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/stk_mesh_fixtures/SelectorFixture.hpp"  // for SelectorFixture
#include "stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp"
#include <algorithm>                    // for find
#include <array>                        // for array
#include <iostream>                     // for ostringstream, ostream, etc
#include <stddef.h>                     // for size_t
#include <stdexcept>                    // for logic_error
#include <string>                       // for basic_string, operator==, etc
#include <utility>                      // for pair
#include <vector>                       // for vector
namespace stk { namespace mesh { class Bucket; } }

namespace {

using stk::mesh::fixtures::SelectorFixture;

void test_selector_with_buckets(const SelectorFixture &selectorFixture, const stk::mesh::Selector &selector, bool gold_shouldEntityBeInSelector[]);
void test_selector_with_partitions(const SelectorFixture &selectorFixture, const stk::mesh::Selector &selector, bool gold_shouldEntityBeInSelector[]);

void initialize(SelectorFixture& fixture)
{
  fixture.m_meta_data.commit();
  fixture.m_bulk_data.modification_begin();
  fixture.generate_mesh();
  ASSERT_TRUE(fixture.m_bulk_data.modification_end());
}

/** \defgroup stk_mesh_selector_unit "stk::mesh::Selector Unit Testing"
  * \addtogroup stk_mesh_selector_unit
  * \{
  *
  * Selector unit testing environment. <br>
  * A special set of mesh parts and entities are set up in the
  * following configuration for the Selector unit testing.<br>
  * Parts:  PartA, PartB, PartC, PartD, PartU <br>
  * PartU = MetaData.universal_part() <br>
  * Entities:  Entity1, Entity2, Entity3, Entity4, Entity5 <br>
  *
  * PartA contains Entity1, Entity2 <br>
  * PartB contains Entity2, Entity3 <br>
  * PartC contains Entity3, Entity4 <br>
  * PartD contains no entities <br>
  * Entity5 is not contained in any Part <br>
  *
  * <PRE>
  * |----------|--|-------|--|----------|    |-------------|
  * |<--PartA---->|       |<--PartC---->|    |   PartD     |
  * |          |<---PartB--->|          |    |             |
  * |  1       |2 |       |3 |       4  | 5  |             |
  * |          |  |       |  |          |    |             |
  * |          |  |       |  |          |    |             |
  * |----------|--|-------|--|----------|    |-------------|
  * </PRE>
  *
  * Note:  The unit test names use the convention of "i" for
  * intersection, "u" for union, and "c" for complement.
  *
  * */

TEST(Verify, selectorFixtureDoesNotSegFault)
{
  {
    SelectorFixture fix;
    initialize(fix);
  }
  EXPECT_TRUE(true);
}

TEST(Verify, twoSelectorFixturesCreatedSequentiallyDoNotSegFault)
{
  {
    SelectorFixture fix;
    initialize(fix);
  }
  {
    SelectorFixture fix;
    initialize(fix);
  }
  EXPECT_TRUE(true);
}


TEST(Verify, partASelector)
{
  SelectorFixture fix;
  initialize(fix);

  stk::mesh::Selector partASelector(fix.m_partA);

  const size_t numExpectedBuckets = 2;
  EXPECT_EQ(numExpectedBuckets, partASelector.get_buckets(stk::topology::NODE_RANK).size());

  EXPECT_FALSE(partASelector.is_empty(stk::topology::NODE_RANK));

  const int numEntities = 5;
  bool gold_shouldEntityBeInSelector[numEntities] = {true, true, false, false, false};

  test_selector_with_buckets(fix, partASelector, gold_shouldEntityBeInSelector);
  test_selector_with_partitions(fix, partASelector, gold_shouldEntityBeInSelector);
}

TEST(Verify, notPartASelector)
{
  SelectorFixture fix;
  initialize(fix);

  stk::mesh::Selector notPartASelector = ! fix.m_partA;

  const int numEntities = 5;
  bool gold_shouldEntityBeInSelector[numEntities] = {false, false, true, true, true};

  test_selector_with_buckets(fix, notPartASelector, gold_shouldEntityBeInSelector);
  test_selector_with_partitions(fix, notPartASelector, gold_shouldEntityBeInSelector);
}

TEST(Verify, emptyPartSelector)
{
  SelectorFixture fix;
  initialize(fix);

  stk::mesh::Selector selector( fix.m_partD );

  EXPECT_TRUE(selector.is_empty(stk::topology::NODE_RANK));

  const int numEntities = 5;
  bool gold_shouldEntityBeInSelector[numEntities] = {false, false, false, false, false};

  test_selector_with_buckets(fix, selector, gold_shouldEntityBeInSelector);
  test_selector_with_partitions(fix, selector, gold_shouldEntityBeInSelector);
}

TEST(Verify, selectorEmptyDuringMeshMod)
{
  const unsigned spatialDim=3;
  stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
  builder.set_spatial_dimension(spatialDim);
  builder.set_entity_rank_names(stk::mesh::entity_rank_names());
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = builder.create();
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  stk::mesh::Part& block1 = meta.declare_part_with_topology("block_1", stk::topology::HEX_8);
  meta.commit();

  stk::mesh::Selector block1Selector = block1;
  EXPECT_TRUE(block1Selector.is_empty(stk::topology::NODE_RANK));
  EXPECT_TRUE(block1Selector.is_empty(stk::topology::ELEM_RANK));

  bulk.modification_begin();

  if (bulk.parallel_rank()==0) {

    stk::mesh::EntityId elem1Id = 1;
    stk::mesh::Entity elem1 = bulk.declare_element(elem1Id, stk::mesh::ConstPartVector{&block1});

    EXPECT_FALSE(block1Selector.is_empty(stk::topology::ELEM_RANK));

    stk::mesh::PartVector addParts;
    stk::mesh::PartVector removeParts(1, &block1);
    bulk.change_entity_parts(elem1, addParts, removeParts);

    EXPECT_TRUE(block1Selector.is_empty(stk::topology::ELEM_RANK));

    addParts.push_back(&block1);
    removeParts.clear();

    bulk.change_entity_parts(elem1, addParts, removeParts);

    EXPECT_FALSE(block1Selector.is_empty(stk::topology::ELEM_RANK));
  }

  //cannot call modification_end which requires all elements, faces, edges to have connected nodes - which we have not defined
  //bulk.modification_end();

  if (bulk.parallel_rank()==0) {
    EXPECT_FALSE(block1Selector.is_empty(stk::topology::ELEM_RANK));
  }
  else {
    EXPECT_TRUE(block1Selector.is_empty(stk::topology::ELEM_RANK));
  }
}

TEST(Verify, selector_declarePartSubset)
{
  SelectorFixture fix;
  stk::mesh::Part& partAsub = fix.m_meta_data.declare_part_with_topology("subPart", stk::topology::NODE);
  fix.m_meta_data.declare_part_subset(fix.m_partA, partAsub);
  EXPECT_EQ(1u, fix.m_fieldABC->restrictions().size());
  const stk::mesh::Selector& fieldABCselector = fix.m_fieldABC->restrictions()[0].selector();
  EXPECT_TRUE(fieldABCselector(partAsub));
}

TEST(Verify, complementOfPartASelector)
{
  SelectorFixture fix;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA;
  stk::mesh::Selector partASelector(partA);

  size_t expected_num_buckets = 2;
  EXPECT_EQ(expected_num_buckets, partASelector.get_buckets(stk::topology::NODE_RANK).size());
  EXPECT_FALSE(partASelector.is_empty(stk::topology::NODE_RANK));
  EXPECT_TRUE(partASelector.is_empty(stk::topology::FACE_RANK));

  stk::mesh::Selector partAComplementSelector = partASelector.complement();

  expected_num_buckets = 3;
  EXPECT_EQ(expected_num_buckets, partAComplementSelector.get_buckets(stk::topology::NODE_RANK).size());

  const int numEntities = 5;
  bool gold_shouldEntityBeInSelector[numEntities] = {false, false, true, true, true};

  test_selector_with_buckets(fix, partAComplementSelector, gold_shouldEntityBeInSelector);
  test_selector_with_partitions(fix, partAComplementSelector, gold_shouldEntityBeInSelector);
}

TEST(Verify, andedSelector)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA;
  stk::mesh::Part & partB = fix.m_partB;

  stk::mesh::Selector andedSelector = partA & partB;
  const int numEntities = 5;
  bool gold_shouldEntityBeInPartASelector[numEntities] = {true , true, false, false, false};
  bool gold_shouldEntityBeInPartBSelector[numEntities] = {false, true, true , false, false};
  bool gold_shouldEntityBeInAndedSelector[numEntities] = {false, true, false, false, false};
  test_selector_with_buckets(fix, partA, gold_shouldEntityBeInPartASelector);
  test_selector_with_buckets(fix, partB, gold_shouldEntityBeInPartBSelector);
  test_selector_with_buckets(fix, andedSelector, gold_shouldEntityBeInAndedSelector);
}

TEST(Verify, oredSelector)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA;
  stk::mesh::Part & partB = fix.m_partB;

  stk::mesh::Selector oredSelector = partA | partB;
  const int numEntities = 5;
  bool gold_shouldEntityBeInPartASelector[numEntities] = {true , true, false, false, false};
  bool gold_shouldEntityBeInPartBSelector[numEntities] = {false, true, true , false, false};
  bool gold_shouldEntityBeInOredSelector[numEntities] =  {true , true, true , false, false};
  test_selector_with_buckets(fix, partA, gold_shouldEntityBeInPartASelector);
  test_selector_with_buckets(fix, partB, gold_shouldEntityBeInPartBSelector);
  test_selector_with_buckets(fix, oredSelector, gold_shouldEntityBeInOredSelector);
}

TEST(Verify, notAndedSelector)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA;
  stk::mesh::Part & partB = fix.m_partB;

  stk::mesh::Selector notAndedSelector = !(partA & partB);

  const int numEntities = 5;
  bool gold_shouldEntityBeInPartASelector[numEntities] = {true , true , false, false, false};
  bool gold_shouldEntityBeInPartBSelector[numEntities] = {false, true , true , false, false};
  bool gold_shouldEntityBeInAndedSelector[numEntities] = {true , false, true , true , true };
  test_selector_with_buckets(fix, partA, gold_shouldEntityBeInPartASelector);
  test_selector_with_buckets(fix, partB, gold_shouldEntityBeInPartBSelector);
  test_selector_with_buckets(fix, notAndedSelector, gold_shouldEntityBeInAndedSelector);
}

TEST(Verify, noredSelector)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA;
  stk::mesh::Part & partB = fix.m_partB;

  stk::mesh::Selector norSelector = !(partA | partB);
  const int numEntities = 5;
  bool gold_shouldEntityBeInPartASelector[numEntities] = {true , true , false, false, false};
  bool gold_shouldEntityBeInPartBSelector[numEntities] = {false, true , true , false, false};
  bool gold_shouldEntityBeInOredSelector[numEntities]  = {false, false, false, true , true };
  test_selector_with_buckets(fix, partA, gold_shouldEntityBeInPartASelector);
  test_selector_with_buckets(fix, partB, gold_shouldEntityBeInPartBSelector);
  test_selector_with_buckets(fix, norSelector, gold_shouldEntityBeInOredSelector);
}

TEST(Verify, differenceSelector)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA;
  stk::mesh::Part & partB = fix.m_partB;

  stk::mesh::Selector differenceSelector = partA - partB;
  stk::mesh::Selector equivSelector = partA & !partB;
  const int numEntities = 5;
  bool gold_shouldEntityBeInPartASelector[numEntities]      = {true , true , false, false, false};
  bool gold_shouldEntityBeInPartBSelector[numEntities]      = {false, true , true , false, false};
  bool gold_shouldEntityBeInDifferenceSelector[numEntities] = {true , false, false, false, false};
  test_selector_with_buckets(fix, partA, gold_shouldEntityBeInPartASelector);
  test_selector_with_buckets(fix, partB, gold_shouldEntityBeInPartBSelector);
  test_selector_with_buckets(fix, differenceSelector, gold_shouldEntityBeInDifferenceSelector);
  test_selector_with_buckets(fix, equivSelector, gold_shouldEntityBeInDifferenceSelector);
}

TEST(Verify, variousSelectorCombinations)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA;
  stk::mesh::Part & partB = fix.m_partB;
  stk::mesh::Part & partC = fix.m_partC;

  stk::mesh::Selector complexSelector = partA & !(partB | partC);

  const int numEntities = 5;
  bool gold_shouldEntityBeInPartBSelector[numEntities]           = {false, true , true , false, false};
  bool gold_shouldEntityBeInPartCSelector[numEntities]           = {false, false, true , true , false};
  bool gold_shouldEntityBeInPartBOrPartCSelector[numEntities]    = {false, true , true , true , false};
  bool gold_shouldEntityNotBeInPartBOrPartCSelector[numEntities] = {true , false, false, false, true };
  bool gold_shouldEntityBeInPartASelector[numEntities]           = {true , true , false, false, false};
  bool gold_shouldEntityBeInComplexSelector[numEntities]         = {true , false, false, false, false};
  test_selector_with_buckets(fix, partB, gold_shouldEntityBeInPartBSelector);
  test_selector_with_buckets(fix, partC, gold_shouldEntityBeInPartCSelector);
  test_selector_with_buckets(fix, partB|partC, gold_shouldEntityBeInPartBOrPartCSelector);
  test_selector_with_buckets(fix, !(partB|partC), gold_shouldEntityNotBeInPartBOrPartCSelector);
  test_selector_with_buckets(fix, partA, gold_shouldEntityBeInPartASelector);
  test_selector_with_buckets(fix, complexSelector, gold_shouldEntityBeInComplexSelector);
}

TEST(Verify, complementOfSelectorComplement)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA ;
  stk::mesh::Part & partB = fix.m_partB ;
  stk::mesh::Selector norSelector = (partA | partB).complement();

  const int numEntities = 5;
  bool gold_shouldEntityBeInPartASelector[numEntities]    = {true , true , false, false, false};
  bool gold_shouldEntityBeInPartBSelector[numEntities]    = {false, true , true , false, false};
  bool gold_shouldEntityBeInNoredSelector[numEntities]    = {false, false, false, true,  true };
  bool gold_shouldEntityBeInNotNoredSelector[numEntities] = {true , true , true , false, false};

  test_selector_with_buckets(fix, partA, gold_shouldEntityBeInPartASelector);
  test_selector_with_buckets(fix, partB, gold_shouldEntityBeInPartBSelector);
  test_selector_with_buckets(fix, norSelector, gold_shouldEntityBeInNoredSelector);

  stk::mesh::Selector selector = norSelector.complement();
  test_selector_with_buckets(fix, selector, gold_shouldEntityBeInNotNoredSelector);
  test_selector_with_buckets(fix, (partA | partB), gold_shouldEntityBeInNotNoredSelector);
}

TEST(Verify, complementOfDefaultConstructedSelector)
{
  SelectorFixture fix ;
  initialize(fix);

  const int numEntities = 5;
  bool goldEntityInDefaultCtor[numEntities]           = {false, false, false, false, false};
  bool goldEntityInDefaultCtorComplement[numEntities] = {true , true , true , true , true };

  stk::mesh::Selector defaultConstructedSelector;
  test_selector_with_buckets(fix, defaultConstructedSelector, goldEntityInDefaultCtor);

  stk::mesh::Selector complementOfDefault = defaultConstructedSelector.complement();
  test_selector_with_buckets(fix, complementOfDefault, goldEntityInDefaultCtorComplement);
}

TEST(Verify, usingPartVectorToSelectIntersection)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::PartVector parts ;
  parts.push_back( & fix.m_partA );
  parts.push_back( & fix.m_partB );
  stk::mesh::Selector selector = selectIntersection(parts);

  const int numEntities = 5;
  bool gold_shouldEntityBeInPartASelector[numEntities] = {true , true, false, false, false};
  bool gold_shouldEntityBeInPartBSelector[numEntities] = {false, true, true , false, false};
  bool gold_shouldEntityBeInAndedSelector[numEntities] = {false, true, false, false, false};
  test_selector_with_buckets(fix, fix.m_partA, gold_shouldEntityBeInPartASelector);
  test_selector_with_buckets(fix, fix.m_partB, gold_shouldEntityBeInPartBSelector);
  test_selector_with_buckets(fix, selector, gold_shouldEntityBeInAndedSelector);

  std::ostringstream msg;
  msg << selector;
  EXPECT_EQ("(PartA & PartB)", msg.str());
}

TEST(Verify, selectUnionIgnoresNullParts)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::PartVector parts ;
  parts.push_back( & fix.m_partA );
  parts.push_back( & fix.m_partB );
  parts.push_back( & fix.m_partC );
  stk::mesh::Selector selector = selectUnion(parts);

  parts.insert(parts.begin(), nullptr);
  parts.push_back(nullptr);
  stk::mesh::Selector selector2 = selectUnion(parts);
  EXPECT_TRUE(selector == selector2);
}

TEST(Verify, usingPartVectorToSelectUnion)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::PartVector parts ;
  parts.push_back( & fix.m_partA );
  parts.push_back( & fix.m_partB );
  parts.push_back( & fix.m_partC );
  stk::mesh::Selector selector = selectUnion(parts);

  const int numEntities = 5;
  bool gold_shouldEntityBeInPartASelector[numEntities]        = {true , true , false, false, false};
  bool gold_shouldEntityBeInPartBSelector[numEntities]        = {false, true , true , false, false};
  bool gold_shouldEntityBeInPartAOrPartBSelector[numEntities] = {true , true , true , false, false};
  bool gold_shouldEntityBeInPartCSelector[numEntities]        = {false, false, true , true , false};
  bool gold_shouldEntityBeInPartOredSelector[numEntities]     = {true , true , true , true , false};
  test_selector_with_buckets(fix, fix.m_partA, gold_shouldEntityBeInPartASelector);
  test_selector_with_buckets(fix, fix.m_partB, gold_shouldEntityBeInPartBSelector);
  test_selector_with_buckets(fix, (fix.m_partA | fix.m_partB), gold_shouldEntityBeInPartAOrPartBSelector);
  test_selector_with_buckets(fix, fix.m_partC, gold_shouldEntityBeInPartCSelector);
  test_selector_with_buckets(fix, selector, gold_shouldEntityBeInPartOredSelector);

  std::ostringstream msg;
  msg << selector;
  std::cout << "msg.str() = " << msg.str() << std::endl;
  EXPECT_EQ( "(PartA | PartB | PartC)", msg.str() );
}

TEST(Verify, defaultConstructorForSelector)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Selector defaultConstructedSelector;

  EXPECT_TRUE(defaultConstructedSelector.is_empty(stk::topology::NODE_RANK));

  const stk::mesh::BucketVector& nodeBuckets = defaultConstructedSelector.get_buckets(stk::topology::NODE_RANK);
  EXPECT_TRUE(nodeBuckets.empty());

  const int numEntities = 5;
  bool goldEntityInDefaultCtor[numEntities] = {false, false, false, false, false};

  test_selector_with_buckets(fix, defaultConstructedSelector, goldEntityInDefaultCtor);
}

TEST(Verify, usingEqualityOperator)
{
  SelectorFixture fix;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA;
  stk::mesh::Selector partASelector(partA);
  stk::mesh::Selector anotherPartASelector(partA);
  EXPECT_TRUE(partASelector == anotherPartASelector);

  stk::mesh::Part & partB = fix.m_partB;
  stk::mesh::Selector unionPartAPartB(partA | partB);
  stk::mesh::Selector anotherUnionPartAPartB(partA | partB);
  EXPECT_TRUE(unionPartAPartB == anotherUnionPartAPartB);

  stk::mesh::Selector intersectionPartAPartB(partA & partB);
  stk::mesh::Selector anotherIntersectionPartAPartB(partA & partB);
  EXPECT_TRUE(intersectionPartAPartB == anotherIntersectionPartAPartB);

  EXPECT_FALSE(unionPartAPartB == intersectionPartAPartB);

  stk::mesh::Selector complementPartA(!partA);
  stk::mesh::Selector anotherComplementPartA(!partA);
  EXPECT_TRUE(complementPartA == anotherComplementPartA);

  EXPECT_FALSE(partASelector == complementPartA);
  EXPECT_FALSE(complementPartA == partASelector);

  stk::mesh::Selector complementPartB(!partB);
  EXPECT_FALSE(complementPartA == complementPartB);

  stk::mesh::Selector notNotPartA(!!partA);
  EXPECT_FALSE(partASelector == notNotPartA);
}

TEST(Verify, usingCopyConstructor)
{
  SelectorFixture fix;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA;
  stk::mesh::Part & partB = fix.m_partB;
  stk::mesh::Part & partC = fix.m_partC;
  stk::mesh::Selector selector = (partA & partB) | partC;
  stk::mesh::Selector anotherSelector(selector);
  EXPECT_TRUE(selector == anotherSelector);
}


TEST(Verify, usingSelectField)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Selector selectFieldA = stk::mesh::selectField(*fix.m_fieldA);
  stk::mesh::Selector selectFieldABC = stk::mesh::selectField(*fix.m_fieldABC);

  stk::mesh::Part & partA = fix.m_partA;
  stk::mesh::Part & partB = fix.m_partB;
  stk::mesh::Part & partC = fix.m_partC;
  stk::mesh::Part & partD = fix.m_partD;

  //
  //  Test selection of buckets
  //
  const int numEntities = 5;
  bool gold_shouldEntityBeInPartASelector[numEntities]         = {true , true , false, false, false};
  bool gold_shouldEntityBeInPartBSelector[numEntities]         = {false, true , true , false, false};
  bool gold_shouldEntityBeInPartCSelector[numEntities]         = {false, false, true , true , false};
  bool gold_shouldEntityBeInPartsABCUnionSelector[numEntities] = {true , true , true , true , false};

  test_selector_with_buckets(fix, selectFieldA, gold_shouldEntityBeInPartASelector);
  test_selector_with_buckets(fix, fix.m_partB, gold_shouldEntityBeInPartBSelector);
  test_selector_with_buckets(fix, fix.m_partC, gold_shouldEntityBeInPartCSelector);
  test_selector_with_buckets(fix, selectFieldABC, gold_shouldEntityBeInPartsABCUnionSelector);

  //
  //  Test selection of parts.  Part selection returns true if there is overlap between the test part set
  //  and any part the field was registered on.
  //
  EXPECT_TRUE (selectFieldA(partA));
  EXPECT_FALSE(selectFieldA(partB));
  EXPECT_FALSE(selectFieldA(partC));
  EXPECT_FALSE(selectFieldA(partD));

  EXPECT_TRUE (selectFieldABC(partA));
  EXPECT_TRUE (selectFieldABC(partB));
  EXPECT_TRUE (selectFieldABC(partC));
  EXPECT_FALSE(selectFieldABC(partD));


  stk::mesh::PartVector partsAB;
  partsAB.push_back(&partA);
  partsAB.push_back(&partB);

  stk::mesh::PartVector partsCD;
  partsAB.push_back(&partC);
  partsAB.push_back(&partD);

  EXPECT_TRUE (selectFieldA(partsAB));
  EXPECT_FALSE(selectFieldA(partsCD));

  //
  //  Check selection of buckets also works in a mesh modification cycle
  //
  fix.get_NonconstBulkData().modification_begin("TEST MODIFICATION");

  test_selector_with_buckets(fix, selectFieldA, gold_shouldEntityBeInPartASelector);
  test_selector_with_buckets(fix, fix.m_partB, gold_shouldEntityBeInPartBSelector);
  test_selector_with_buckets(fix, fix.m_partC, gold_shouldEntityBeInPartCSelector);
  test_selector_with_buckets(fix, selectFieldABC, gold_shouldEntityBeInPartsABCUnionSelector);


  fix.get_NonconstBulkData().modification_end();


}

TEST(Verify, selectorContainsPart)
{
  SelectorFixture fix;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA;
  stk::mesh::Part & partB = fix.m_partB;
  stk::mesh::Part & partC = fix.m_partC;
  stk::mesh::Part & partD = fix.m_partD;
  stk::mesh::Selector selector = partA | partB | (!partC) | partD;
  std::cout << "select_part selector = " << selector << std::endl;
  EXPECT_TRUE(selector(partA));
  EXPECT_TRUE(selector(partB));
  EXPECT_FALSE(selector(partC));
  EXPECT_TRUE(selector(partD));

  selector = partA | ((!((partA & partB) | partC)) & ((!partD) | partB));
  EXPECT_TRUE(selector(partA));
  EXPECT_TRUE(selector(partB));
  EXPECT_FALSE(selector(partC));
  EXPECT_FALSE(selector(partD));

  selector = partC & (!partD);
  EXPECT_FALSE(selector(partA));
  EXPECT_FALSE(selector(partB));
  EXPECT_TRUE(selector(partC));
  EXPECT_FALSE(selector(partD));
}

TEST(Verify, printingOfSelectorUnion)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA ;
  stk::mesh::Part & partB = fix.m_partB ;
  stk::mesh::Part & partC = fix.m_partC ;
  stk::mesh::Part & partD = fix.m_partD ;
  stk::mesh::Selector selector = partA | partB | partC | partD;
  std::cout << "A|B|C|D = " << selector << std::endl;
  std::ostringstream msg;
  msg << selector;
  EXPECT_EQ( "(PartA | PartB | PartC | PartD)" , msg.str() );
}

TEST(Verify, printingOfSelectorIntersection)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA ;
  stk::mesh::Part & partB = fix.m_partB ;
  stk::mesh::Part & partC = fix.m_partC ;
  stk::mesh::Selector selector = partA & partB & partC;
  std::cout << "A&B&C = " << selector << std::endl;
  std::ostringstream msg;
  msg << selector;
  EXPECT_TRUE( msg.str() == "(PartA & PartB & PartC)" );
}

TEST(Verify, printingOfGeneralSelector)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA ;
  stk::mesh::Part & partB = fix.m_partB ;
  stk::mesh::Part & partC = fix.m_partC ;
  stk::mesh::Part & partD = fix.m_partD ;
  stk::mesh::Selector selector =  partA | ( (!((partA & partB) | partC))  & ((!partD) | partB));
  std::cout << "complicated selector = " << selector << std::endl;
  std::ostringstream msg;
  msg << selector;
  EXPECT_EQ( "(PartA | (!(((PartA & PartB) | PartC)) & (!(PartD) | PartB)))" , msg.str() );
}

TEST(Verify, printingOfNothingForComplementOfDefaultSelector)
{
  stk::mesh::Selector selectAll;
  selectAll.complement();

  stk::mesh::Selector anotherSelectAll;
  anotherSelectAll.complement();

  {
    stk::mesh::Selector selectAllANDAll = selectAll & anotherSelectAll;
    std::ostringstream description;
    description << selectAllANDAll;
    EXPECT_EQ( "(!(NOTHING) & !(NOTHING))", description.str());

    //will throw because the selector doesn't have access to a mesh
    EXPECT_THROW(selectAllANDAll.get_buckets(stk::topology::NODE_RANK), std::logic_error);
  }
  {
    stk::mesh::Selector selectAllORAll = selectAll | anotherSelectAll;
    std::ostringstream description;
    description << selectAllORAll;
    EXPECT_EQ( "(!(NOTHING) | !(NOTHING))", description.str());
  }
}

TEST(Verify, usingLessThanOperator)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA ;
  stk::mesh::Part & partB = fix.m_partB ;
  stk::mesh::Part & partC = fix.m_partC ;
  stk::mesh::Part & partD = fix.m_partD ;

  const unsigned ordA = partA.mesh_meta_data_ordinal();
  const unsigned ordB = partB.mesh_meta_data_ordinal();
  const unsigned ordC = partC.mesh_meta_data_ordinal();
  const unsigned ordD = partD.mesh_meta_data_ordinal();

  ASSERT_TRUE(ordA < ordB && ordB < ordC && ordC < ordD);

  stk::mesh::Selector partASelector = partA;
  EXPECT_FALSE(partASelector < partASelector);

  stk::mesh::Selector partBSelector = partB;
  EXPECT_TRUE(partASelector < partBSelector);

  stk::mesh::Selector partAIntersectPartDSelector = partA & partD;
  stk::mesh::Selector partBIntersectPartCSelector = partB & partC;
  EXPECT_TRUE(partAIntersectPartDSelector < partBIntersectPartCSelector);

  stk::mesh::Selector partCIntersectPartBSelectorOrderMatters = partC & partB;
  stk::mesh::Selector partDIntersectPartASelectorOrderMatters = partD & partA;
  EXPECT_TRUE(partDIntersectPartASelectorOrderMatters < partCIntersectPartBSelectorOrderMatters);

  stk::mesh::Selector partAUnionPartBSelector     = partA | partB;
  stk::mesh::Selector partAIntersectPartBSelector = partA & partB;
  EXPECT_TRUE(partASelector < partAUnionPartBSelector);
  EXPECT_TRUE(partAUnionPartBSelector < partAIntersectPartBSelector);

  stk::mesh::Selector partAUnionPartBIntersectPartCSelector = partAUnionPartBSelector & partC;
  EXPECT_TRUE(partAUnionPartBSelector < partAUnionPartBIntersectPartCSelector);
}

void test_selector_with_buckets(const SelectorFixture &selectorFixture, const stk::mesh::Selector &selector, bool gold_shouldEntityBeInSelector[])
{
  const stk::mesh::BulkData& stkMeshBulkData = selectorFixture.get_BulkData();
  {
    const stk::mesh::Bucket & bucket = stkMeshBulkData.bucket(selectorFixture.m_entity1);
    bool result = selector(bucket);
    EXPECT_EQ(gold_shouldEntityBeInSelector[0], result);
  }
  {
    const stk::mesh::Bucket & bucket = stkMeshBulkData.bucket(selectorFixture.m_entity2);
    bool result = selector(bucket);
    EXPECT_EQ(gold_shouldEntityBeInSelector[1], result);
  }
  {
    const stk::mesh::Bucket & bucket = stkMeshBulkData.bucket(selectorFixture.m_entity3);
    bool result = selector(bucket);
    EXPECT_EQ(gold_shouldEntityBeInSelector[2], result);
  }
  {
    const stk::mesh::Bucket & bucket = stkMeshBulkData.bucket(selectorFixture.m_entity4);
    bool result = selector(bucket);
    EXPECT_EQ(gold_shouldEntityBeInSelector[3], result);
  }
  {
    const stk::mesh::Bucket & bucket = stkMeshBulkData.bucket(selectorFixture.m_entity5);
    bool result = selector(bucket);
    EXPECT_EQ(gold_shouldEntityBeInSelector[4], result);
  }
}

void test_selector_with_partitions(const SelectorFixture &selectorFixture, const stk::mesh::Selector &selector, bool gold_shouldEntityBeInSelector[])
{
  const stk::mesh::BulkData& stkMeshBulkData = selectorFixture.get_BulkData();
  {
    const stk::mesh::impl::Partition* partition = stkMeshBulkData.bucket(selectorFixture.m_entity1).getPartition();
    bool result = selector(*partition);
    EXPECT_EQ(gold_shouldEntityBeInSelector[0], result);
  }
  {
    const stk::mesh::impl::Partition* partition = stkMeshBulkData.bucket(selectorFixture.m_entity2).getPartition();
    bool result = selector(*partition);
    EXPECT_EQ(gold_shouldEntityBeInSelector[1], result);
  }
  {
    const stk::mesh::impl::Partition* partition = stkMeshBulkData.bucket(selectorFixture.m_entity3).getPartition();
    bool result = selector(*partition);
    EXPECT_EQ(gold_shouldEntityBeInSelector[2], result);
  }
  {
    const stk::mesh::impl::Partition* partition = stkMeshBulkData.bucket(selectorFixture.m_entity4).getPartition();
    bool result = selector(*partition);
    EXPECT_EQ(gold_shouldEntityBeInSelector[3], result);
  }
  {
    const stk::mesh::impl::Partition* partition = stkMeshBulkData.bucket(selectorFixture.m_entity5).getPartition();
    bool result = selector(*partition);
    EXPECT_EQ(gold_shouldEntityBeInSelector[4], result);
  }
}

std::shared_ptr<stk::mesh::BulkData> create_mesh(stk::ParallelMachine comm,
                                                 unsigned spatialDim)
{
  stk::mesh::MeshBuilder builder(comm);
  builder.set_spatial_dimension(spatialDim);
  std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();
  return bulk;
}

void check_selector_does_not_return_root_topology_parts(stk::ParallelMachine pm, const std::string &part_name, stk::topology topo )
{
  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> meshPtr = create_mesh(pm, spatialDim);
  stk::mesh::BulkData& mesh = *meshPtr;
  stk::mesh::MetaData& meta = mesh.mesh_meta_data();

  stk::mesh::Part * part = &meta.declare_part_with_topology(part_name, topo);
  meta.commit();

  stk::mesh::Selector selector(*part);

  stk::mesh::PartVector selectorParts;

  selector.get_parts(selectorParts);

  stk::mesh::Part &rootPart= meta.get_topology_root_part(topo);

  auto iterator = std::find(selectorParts.begin(), selectorParts.end(), &rootPart );

  EXPECT_EQ(1u, selectorParts.size());

  EXPECT_TRUE(selectorParts.end() == iterator);

  iterator = std::find(selectorParts.begin(), selectorParts.end(), part );

  EXPECT_TRUE(selectorParts.end() != iterator);
}

TEST( UnitTestRootTopology, getPartsDoesNotFindAutoCreatedRootParts )
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);

  if(p_size > 1)
  {
    return;
  }

  std::vector < stk::topology > test_topologies {
    stk::topology::NODE
        //EDGE_RANK
        , stk::topology::LINE_2
        , stk::topology::LINE_3
        //FACE_RANK
        , stk::topology::TRI_3
        , stk::topology::TRI_4
        , stk::topology::TRI_6
        , stk::topology::QUAD_4
        , stk::topology::QUAD_6
        , stk::topology::QUAD_8
        , stk::topology::QUAD_9
        //ELEMENT_RANK
        , stk::topology::PARTICLE
        , stk::topology::BEAM_2
        , stk::topology::BEAM_3
        , stk::topology::SHELL_TRI_3
        // Cannot create SHELL_TRI_4 entities because Shards does not support!
        // , stk::topology::SHELL_TRI_4
        , stk::topology::SHELL_TRI_6
        , stk::topology::SHELL_QUAD_4
        , stk::topology::SHELL_QUAD_8
        , stk::topology::SHELL_QUAD_9
        , stk::topology::TET_4
        , stk::topology::TET_8
        , stk::topology::TET_10
        , stk::topology::TET_11
        , stk::topology::PYRAMID_5
        , stk::topology::PYRAMID_13
        , stk::topology::PYRAMID_14
        , stk::topology::WEDGE_6
        , stk::topology::WEDGE_15
        , stk::topology::WEDGE_18
        , stk::topology::HEX_8
        , stk::topology::HEX_20
        , stk::topology::HEX_27};

  for (unsigned i = 0; i < test_topologies.size(); ++i)
  {
    check_selector_does_not_return_root_topology_parts(pm, "topo_part" , test_topologies[i]);
  }
}

TEST(Selector, get_parts_intersection_ranked)
{
  const unsigned spatialDim = 3;
  stk::mesh::MetaData meta(spatialDim);

  stk::mesh::Part& rankedPart = meta.declare_part("ranked", stk::topology::ELEM_RANK);
  stk::mesh::Part& unrankedPart = meta.declare_part("unranked");
  meta.commit();

  EXPECT_TRUE(unrankedPart.primary_entity_rank() == stk::topology::INVALID_RANK);

  //It doesn't seem right that an intersection selector should only return the ranked
  //part in get_parts, but applications have come to depend on this. So this test
  //codifies the behavior to make sure that stk testing will catch any unintentional
  //breakages in stk changes.
  stk::mesh::Selector selector1 = rankedPart & unrankedPart;
  stk::mesh::PartVector selector1Parts;
  selector1.get_parts(selector1Parts);
  EXPECT_EQ(1u, selector1Parts.size());
  EXPECT_EQ(rankedPart.mesh_meta_data_ordinal(), selector1Parts[0]->mesh_meta_data_ordinal());

  stk::mesh::PartVector selector2Parts;
  selector1.get_parts(selector2Parts);
  EXPECT_EQ(1u, selector2Parts.size());
  EXPECT_EQ(rankedPart.mesh_meta_data_ordinal(), selector2Parts[0]->mesh_meta_data_ordinal());
}

TEST(Selector, partIntersection)
{
  stk::ParallelMachine comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) != 1) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = stk::mesh::MeshBuilder(comm).set_spatial_dimension(3).create();
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::Part& sidePart = meta.declare_part_with_topology("sidePart", stk::topology::QUAD_4);
  constexpr int NX=2, NY=2, NZ=2;
  stk::mesh::fixtures::HexFixture::fill_mesh(NX, NY, NZ, *bulkPtr);
  stk::mesh::create_exposed_block_boundary_sides(*bulkPtr, meta.universal_part(), {&sidePart});
  EXPECT_EQ(24u, stk::mesh::count_entities(*bulkPtr, stk::topology::FACE_RANK, sidePart));

  stk::mesh::Entity node1 = bulkPtr->get_entity(stk::topology::NODE_RANK, 1);
  EXPECT_TRUE(bulkPtr->is_valid(node1));
  const stk::mesh::Bucket& node1Bucket = bulkPtr->bucket(node1);

  stk::mesh::Selector validIntersection = sidePart & meta.locally_owned_part();
  EXPECT_TRUE(validIntersection(node1Bucket));
  EXPECT_TRUE(validIntersection(node1Bucket.supersets()));
  stk::mesh::impl::Partition& node1Partition = *node1Bucket.getPartition();
  EXPECT_TRUE(validIntersection(node1Partition));

  stk::mesh::Selector invalidIntersection = stk::mesh::Selector() & meta.locally_owned_part();
  EXPECT_FALSE(invalidIntersection(node1Bucket));
  EXPECT_FALSE(invalidIntersection(node1Bucket.supersets()));
  EXPECT_FALSE(invalidIntersection(node1Partition));
}

TEST(Selector, selectUnion_clone_for_different_mesh)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }

  constexpr unsigned spatialDim = 3;
  stk::mesh::MetaData meta1(spatialDim);
  std::string part1Name("myPart1"), part2Name("myPart2");
  stk::mesh::PartVector parts1 = { &meta1.declare_part(part1Name), &meta1.declare_part(part2Name)};
  stk::mesh::Selector select1 = stk::mesh::selectUnion(parts1);

  stk::mesh::MetaData meta2;
  meta2.declare_part(part1Name);
  meta2.declare_part(part2Name);
  meta2.initialize(spatialDim);

  stk::mesh::Selector select2 = select1.clone_for_different_mesh(meta2);
  stk::mesh::PartVector parts2;
  select2.get_parts(parts2);

  EXPECT_EQ(2u, parts2.size());
  EXPECT_NE(parts2[0]->mesh_meta_data_ordinal(), parts1[0]->mesh_meta_data_ordinal());
  EXPECT_NE(parts2[1]->mesh_meta_data_ordinal(), parts1[1]->mesh_meta_data_ordinal());
  EXPECT_EQ(parts2[0]->name(), part1Name);
  EXPECT_EQ(parts2[1]->name(), part2Name);
}

TEST( UnitTestRootTopology, bucketAlsoHasAutoCreatedRootParts )
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(pm) > 1) { GTEST_SKIP(); }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> meshPtr = create_mesh(pm, spatialDim);
  stk::mesh::BulkData& mesh = *meshPtr;
  stk::mesh::MetaData& meta = mesh.mesh_meta_data();

  stk::mesh::Part * triPart = &meta.declare_part_with_topology("tri_part", stk::topology::TRI_3);
  stk::mesh::Part * shellPart = &meta.declare_part_with_topology("shell_part", stk::topology::SHELL_TRI_3);
  meta.commit();

  mesh.modification_begin();

  std::array<int, 3> node_ids = {{1,2,3}};
  stk::mesh::ConstPartVector empty;
  stk::mesh::Entity shell3 = mesh.declare_element(1u, stk::mesh::ConstPartVector{shellPart});
  for(unsigned i = 0; i<node_ids.size(); ++i) {
    stk::mesh::Entity node = mesh.declare_node(node_ids[i], empty);
    mesh.declare_relation(shell3, node, i);
  }
  mesh.declare_element_side(shell3, 0u, stk::mesh::ConstPartVector{triPart} );

  mesh.modification_end();

  stk::mesh::Selector triSelector(*triPart);

  const stk::mesh::BucketVector & buckets = mesh.get_buckets(stk::topology::FACE_RANK, triSelector);

  EXPECT_EQ(1u, buckets.size());

  {
    const stk::mesh::PartVector &triBucketParts = buckets[0]->supersets();

    stk::mesh::Part &triRootPart= meta.get_topology_root_part(stk::topology::TRI_3);

    auto triIterator = std::find(triBucketParts.begin(), triBucketParts.end(), &triRootPart );

    EXPECT_TRUE(triBucketParts.end() != triIterator);

    EXPECT_TRUE(triRootPart.contains(*triPart));
  }
  {
    stk::mesh::OrdinalVector ordinals;
    buckets[0]->supersets(ordinals);

    stk::mesh::Part &triRootPart= meta.get_topology_root_part(stk::topology::TRI_3);

    auto triIterator = std::find(ordinals.begin(), ordinals.end(), triRootPart.mesh_meta_data_ordinal());

    EXPECT_TRUE(ordinals.end() != triIterator);
  }
  {
    std::pair<const unsigned *, const unsigned *> ords_range = buckets[0]->superset_part_ordinals();

    stk::mesh::Part &triRootPart= meta.get_topology_root_part(stk::topology::TRI_3);

    auto triIterator = std::find(ords_range.first, ords_range.second, triRootPart.mesh_meta_data_ordinal());

    EXPECT_TRUE(ords_range.second != triIterator);
  }
}

} // namespace
