/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>                     // for ostringstream, ostream, etc
#include <stk_mesh/base/Part.hpp>       // for Part
#include <stk_mesh/base/Selector.hpp>   // for Selector, operator|, etc
#include <stk_mesh/base/Types.hpp>      // for PartVector
#include <stk_mesh/fixtures/SelectorFixture.hpp>  // for SelectorFixture
#include <gtest/gtest.h>
#include <string>                       // for basic_string, operator==, etc
#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
namespace stk { namespace mesh { class Bucket; } }




// Unit test the Selector in isolation

namespace {

using stk::mesh::fixtures::SelectorFixture ;

void testSelectorWithBuckets(const SelectorFixture &selectorFixture, const stk::mesh::Selector &selector, bool gold_shouldEntityBeInSelector[]);

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

TEST(Verify, twoSelectorFixturesCreatedSequetiallyDoNotSegFault)
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

  testSelectorWithBuckets(fix, partASelector, gold_shouldEntityBeInSelector);
}

TEST(Verify, notPartASelector)
{
  SelectorFixture fix;
  initialize(fix);

  stk::mesh::Selector notPartASelector = ! fix.m_partA;

  const int numEntities = 5;
  bool gold_shouldEntityBeInSelector[numEntities] = {false, false, true, true, true};

  testSelectorWithBuckets(fix, notPartASelector, gold_shouldEntityBeInSelector);
}

TEST(Verify, emptyPartSelector)
{
  SelectorFixture fix;
  initialize(fix);

  stk::mesh::Selector selector( fix.m_partD );

  EXPECT_TRUE(selector.is_empty(stk::topology::NODE_RANK));

  const int numEntities = 5;
  bool gold_shouldEntityBeInSelector[numEntities] = {false, false, false, false, false};

  testSelectorWithBuckets(fix, selector, gold_shouldEntityBeInSelector);
}

TEST(Verify, selectorEmptyDuringMeshMod)
{
    const unsigned spatialDim=3;
    stk::mesh::MetaData meta(spatialDim, stk::mesh::entity_rank_names());
    stk::mesh::Part& block1 = meta.declare_part_with_topology("block_1", stk::topology::HEX_8);
    meta.commit();
    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);

    stk::mesh::Selector block1Selector = block1;
    EXPECT_TRUE(block1Selector.is_empty(stk::topology::NODE_RANK));
    EXPECT_TRUE(block1Selector.is_empty(stk::topology::ELEM_RANK));

    bulk.modification_begin();

    if (bulk.parallel_rank()==0) {

        stk::mesh::EntityId elem1Id = 1;
        stk::mesh::Entity elem1 = bulk.declare_entity(stk::topology::ELEM_RANK, elem1Id, block1);

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

  testSelectorWithBuckets(fix, partAComplementSelector, gold_shouldEntityBeInSelector);
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
  testSelectorWithBuckets(fix, partA, gold_shouldEntityBeInPartASelector);
  testSelectorWithBuckets(fix, partB, gold_shouldEntityBeInPartBSelector);
  testSelectorWithBuckets(fix, andedSelector, gold_shouldEntityBeInAndedSelector);
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
  testSelectorWithBuckets(fix, partA, gold_shouldEntityBeInPartASelector);
  testSelectorWithBuckets(fix, partB, gold_shouldEntityBeInPartBSelector);
  testSelectorWithBuckets(fix, oredSelector, gold_shouldEntityBeInOredSelector);
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
  testSelectorWithBuckets(fix, partA, gold_shouldEntityBeInPartASelector);
  testSelectorWithBuckets(fix, partB, gold_shouldEntityBeInPartBSelector);
  testSelectorWithBuckets(fix, notAndedSelector, gold_shouldEntityBeInAndedSelector);
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
  testSelectorWithBuckets(fix, partA, gold_shouldEntityBeInPartASelector);
  testSelectorWithBuckets(fix, partB, gold_shouldEntityBeInPartBSelector);
  testSelectorWithBuckets(fix, norSelector, gold_shouldEntityBeInOredSelector);
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
  testSelectorWithBuckets(fix, partA, gold_shouldEntityBeInPartASelector);
  testSelectorWithBuckets(fix, partB, gold_shouldEntityBeInPartBSelector);
  testSelectorWithBuckets(fix, differenceSelector, gold_shouldEntityBeInDifferenceSelector);
  testSelectorWithBuckets(fix, equivSelector, gold_shouldEntityBeInDifferenceSelector);
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
    testSelectorWithBuckets(fix, partB, gold_shouldEntityBeInPartBSelector);
    testSelectorWithBuckets(fix, partC, gold_shouldEntityBeInPartCSelector);
    testSelectorWithBuckets(fix, partB|partC, gold_shouldEntityBeInPartBOrPartCSelector);
    testSelectorWithBuckets(fix, !(partB|partC), gold_shouldEntityNotBeInPartBOrPartCSelector);
    testSelectorWithBuckets(fix, partA, gold_shouldEntityBeInPartASelector);
    testSelectorWithBuckets(fix, complexSelector, gold_shouldEntityBeInComplexSelector);
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

  testSelectorWithBuckets(fix, partA, gold_shouldEntityBeInPartASelector);
  testSelectorWithBuckets(fix, partB, gold_shouldEntityBeInPartBSelector);
  testSelectorWithBuckets(fix, norSelector, gold_shouldEntityBeInNoredSelector);

  stk::mesh::Selector selector = norSelector.complement();
  testSelectorWithBuckets(fix, selector, gold_shouldEntityBeInNotNoredSelector);
  testSelectorWithBuckets(fix, (partA | partB), gold_shouldEntityBeInNotNoredSelector);
}

TEST(Verify, complementOfDefaultConstructedSelector)
{
    SelectorFixture fix ;
    initialize(fix);

    const int numEntities = 5;
    bool goldEntityInDefaultCtor[numEntities]           = {false, false, false, false, false};
    bool goldEntityInDefaultCtorComplement[numEntities] = {true , true , true , true , true };

    stk::mesh::Selector defaultConstructedSelector;
    testSelectorWithBuckets(fix, defaultConstructedSelector, goldEntityInDefaultCtor);

    stk::mesh::Selector complementOfDefault = defaultConstructedSelector.complement();
    testSelectorWithBuckets(fix, complementOfDefault, goldEntityInDefaultCtorComplement);
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
  testSelectorWithBuckets(fix, fix.m_partA, gold_shouldEntityBeInPartASelector);
  testSelectorWithBuckets(fix, fix.m_partB, gold_shouldEntityBeInPartBSelector);
  testSelectorWithBuckets(fix, selector, gold_shouldEntityBeInAndedSelector);

  std::ostringstream msg;
  msg << selector;
  EXPECT_EQ("(PartA & PartB)", msg.str());
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
  testSelectorWithBuckets(fix, fix.m_partA, gold_shouldEntityBeInPartASelector);
  testSelectorWithBuckets(fix, fix.m_partB, gold_shouldEntityBeInPartBSelector);
  testSelectorWithBuckets(fix, (fix.m_partA | fix.m_partB), gold_shouldEntityBeInPartAOrPartBSelector);
  testSelectorWithBuckets(fix, fix.m_partC, gold_shouldEntityBeInPartCSelector);
  testSelectorWithBuckets(fix, selector, gold_shouldEntityBeInPartOredSelector);

  std::ostringstream msg;
  msg << selector;
  std::cout << "msg.str() = " << msg.str() << std::endl;
  EXPECT_EQ( "((PartA | PartB) | PartC)", msg.str() );
}

TEST(Verify, defaultConstructorForSelector)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Selector defaultConstructedSelector;

  const int numEntities = 5;
  bool goldEntityInDefaultCtor[numEntities] = {false, false, false, false, false};

  testSelectorWithBuckets(fix, defaultConstructedSelector, goldEntityInDefaultCtor);
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

  stk::mesh::Selector selectFieldA = stk::mesh::selectField(fix.m_fieldA);
  stk::mesh::Selector selectFieldABC = stk::mesh::selectField(fix.m_fieldABC);

  const int numEntities = 5;
  bool gold_shouldEntityBeInPartASelector[numEntities]         = {true , true , false, false, false};
  bool gold_shouldEntityBeInPartBSelector[numEntities]         = {false, true , true , false, false};
  bool gold_shouldEntityBeInPartCSelector[numEntities]         = {false, false, true , true , false};
  bool gold_shouldEntityBeInPartsABCUnionSelector[numEntities] = {true , true , true , true , false};

  testSelectorWithBuckets(fix, selectFieldA, gold_shouldEntityBeInPartASelector);
  testSelectorWithBuckets(fix, fix.m_partB, gold_shouldEntityBeInPartBSelector);
  testSelectorWithBuckets(fix, fix.m_partC, gold_shouldEntityBeInPartCSelector);
  testSelectorWithBuckets(fix, selectFieldABC, gold_shouldEntityBeInPartsABCUnionSelector);
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

    selector = partA | (!((partA & partB) | partC) & (!partD | partB));
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
  EXPECT_EQ( "(((PartA | PartB) | PartC) | PartD)" , msg.str() );
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
  EXPECT_TRUE( msg.str() == "((PartA & PartB) & PartC)" );
}

TEST(Verify, printingOfGeneralSelector)
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA ;
  stk::mesh::Part & partB = fix.m_partB ;
  stk::mesh::Part & partC = fix.m_partC ;
  stk::mesh::Part & partD = fix.m_partD ;
  stk::mesh::Selector selector =  partA | ( !((partA & partB) | partC)  & (!partD | partB));
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
  EXPECT_TRUE(partCIntersectPartBSelectorOrderMatters < partDIntersectPartASelectorOrderMatters);

  stk::mesh::Selector partAUnionPartBSelector     = partA | partB;
  stk::mesh::Selector partAIntersectPartBSelector = partA & partB;
  EXPECT_TRUE(partASelector < partAUnionPartBSelector);
  EXPECT_TRUE(partAUnionPartBSelector < partAIntersectPartBSelector);

  stk::mesh::Selector partAUnionPartBIntersectPartCSelector = partAUnionPartBSelector & partC;
  EXPECT_TRUE(partAUnionPartBSelector < partAUnionPartBIntersectPartCSelector);
}

void testSelectorWithBuckets(const SelectorFixture &selectorFixture, const stk::mesh::Selector &selector, bool gold_shouldEntityBeInSelector[])
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

} // namespace
