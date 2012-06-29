/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_util/environment/WallTime.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/fixtures/SelectorFixture.hpp>

// Unit test the Selector in isolation

namespace {

using stk::mesh::fixtures::SelectorFixture ;

void initialize(SelectorFixture& fixture)
{
  fixture.m_meta_data.commit();
  fixture.m_bulk_data.modification_begin();
  fixture.generate_mesh();
  STKUNIT_ASSERT(fixture.m_bulk_data.modification_end());
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

/** \brief Verify we can construct the selector unit testing fixture.
 *
 * */
STKUNIT_UNIT_TEST( UnitTestSelector, one_SelectorFixture )
{
  {
    SelectorFixture fix;
    initialize(fix);
  }
  STKUNIT_EXPECT_TRUE(true);
}


/** \brief Verify we can construct two selector unit testing fixtures one after another.
 *
 * */
STKUNIT_UNIT_TEST( UnitTestSelector, two_SelectorFixture )
{
  {
    SelectorFixture fix;
    initialize(fix);
  }
  {
    SelectorFixture fix;
    initialize(fix);
  }
  STKUNIT_EXPECT_TRUE(true);
}



/** \brief Test containment directly.
 *
 * Verify PartA contains Entity1 and Entity2, and does not contain
 * Entity3, Entity4, or Entity5.
 * */
STKUNIT_UNIT_TEST( UnitTestSelector, A_12345 )
{
  SelectorFixture fix;
  initialize(fix);

  stk::mesh::Selector selector( fix.m_partA );
  //std::cout << "Selector = " << selector << std::endl;

  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket() ;
    bool result = selector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity2->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity3->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity4->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity5->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
}

/** \brief Test containment with the complement.
 *
 * Verify !PartA does not contain Entity1 and Entity2, and does
 * contain Entity3, Entity4, and Entity5.
 ** */
STKUNIT_UNIT_TEST( UnitTestSelector, Ac_12345 )
{
  SelectorFixture fix;
  initialize(fix);

  stk::mesh::Selector selector = ! fix.m_partA ;
  //std::cout << "Selector = " << selector << std::endl;

  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity2->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity3->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity4->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity5->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }

}

/** \brief Verify PartD does not contain Entity5.
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, D_5 )
{
  SelectorFixture fix;
  initialize(fix);

  stk::mesh::Selector selector( fix.m_partD );

  const stk::mesh::Bucket & bucket = fix.m_entity5->bucket();

  bool result = selector(bucket);
  STKUNIT_EXPECT_FALSE(result);
}

/** \brief Verify PartA.complement contains Entity1 and Entity5.
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, Ac_15 )
{
  SelectorFixture fix;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA ;

  stk::mesh::Selector selector(partA);
  selector.complement();
  //std::cout << "Selector = " << selector << std::endl;

  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity5->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
}

/** \brief Verify (PartA AND PartB) does not contain Entity1 and does
 * contain Entity2.
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, AiB_12 )
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA ;
  stk::mesh::Part & partB = fix.m_partB ;

  stk::mesh::Selector selector = partA & partB;
  //std::cout << "Selector = " << selector << std::endl;

  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity2->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
}


/** \brief Verify (PartA OR PartB) contains Entity1 but not Entity4
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, AuB_14 )
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA ;
  stk::mesh::Part & partB = fix.m_partB ;

  stk::mesh::Selector selector = partA | partB;
  //std::cout << "Selector = " << selector << std::endl;

  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity4->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
}


/** \brief Verify !(PartA AND PartB) contains Entity1 but not Entity2.
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, AiBc_12 )
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA ;
  stk::mesh::Part & partB = fix.m_partB ;

  stk::mesh::Selector selector = partA & !partB;
  //std::cout << "Selector = " << selector << std::endl;

  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity2->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
}


/** \brief Verify !(PartA OR PartB) contains Entity1 but not Entity3.
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, AuBc_13 )
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA ;
  stk::mesh::Part & partB = fix.m_partB ;

  stk::mesh::Selector selector = partA | !partB;
  //std::cout << "Selector = " << selector << std::endl;

  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity3->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
}


// plus copy constructor
/** \brief Verify (PartA AND !(PartB OR PartC)) contains Entity1 but
 * not Entity2.
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, Ai_BuC_c_12 )
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA ;
  stk::mesh::Part & partB = fix.m_partB ;
  stk::mesh::Part & partC = fix.m_partC ;

  stk::mesh::Selector selector = partA & !(partB | partC);
  //std::cout << "Selector = " << selector << std::endl;

  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity2->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }

  stk::mesh::Selector newSelector(selector);
  // Should be the same:
  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = newSelector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity2->bucket();
    bool result = newSelector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }

}


/** \brief test on Selector operator for Entity
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, entityTest )
{
  {
    SelectorFixture fix ;
    initialize(fix);

    stk::mesh::Part & partA = fix.m_partA ;
    stk::mesh::Part & partB = fix.m_partB ;
    stk::mesh::Selector selector = partA & !partB;

    const stk::mesh::Entity & pEntity = *fix.m_entity5;
    bool result = selector(pEntity);
    STKUNIT_EXPECT_FALSE(result);
  }

}

/** \brief Verify the default constructor does not contain Entity1.
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, defaultConstructor )
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Selector selector;
  //std::cout << "Selector = " << selector << std::endl;

  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
}


/** \brief Verify flipping the complement bit on an OR expression
 * works correctly.
 * Verify !(PartA OR PartB) does not contain Entity1, Entity2,
 * or Entity3, and does contain Entity4.  Then check that !!(PartA OR
 * PartB) does contain Entity1, Entity2, and Entity3, and not Entity4.
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, flipComplement_AuB_c )
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA ;
  stk::mesh::Part & partB = fix.m_partB ;
  stk::mesh::Selector notOrSelector = partA | partB;
  //std::cout << "Or Selector = " << notOrSelector << std::endl;
  notOrSelector.complement();
  //std::cout << "Not Or Selector = " << notOrSelector << std::endl;

  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = notOrSelector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity2->bucket();
    bool result = notOrSelector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity3->bucket();
    bool result = notOrSelector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity4->bucket();
    bool result = notOrSelector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }

  stk::mesh::Selector notNotOrSelector = !notOrSelector;
  //std::cout << "Not Not Or Selector = " << notNotOrSelector << std::endl;
  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = notNotOrSelector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity2->bucket();
    bool result = notNotOrSelector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity3->bucket();
    bool result = notNotOrSelector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity4->bucket();
    bool result = notNotOrSelector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
}

/** \brief Verify flipping the complement bit on an AND expression
 * works correctly.
 * Verify !(PartA AND PartB) does not contain Entity2, and does
 * contain Entity1, Entity3, and Entity4.
 * Then check that !!(PartA AND PartB) does contain Entity2, but not Entity1, Entity3, or Entity4.
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, flipComplement_AiB_c )
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA ;
  stk::mesh::Part & partB = fix.m_partB ;
  stk::mesh::Selector notAndSelector = partA & partB;
  //std::cout << "And Selector = " << notAndSelector << std::endl;
  notAndSelector.complement();
  //std::cout << "Not And Selector = " << notAndSelector << std::endl;

  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = notAndSelector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity2->bucket();
    bool result = notAndSelector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity3->bucket();
    bool result = notAndSelector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity4->bucket();
    bool result = notAndSelector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }

  stk::mesh::Selector notNotAndSelector = !notAndSelector;
  //std::cout << "Not Not And Selector = " << notNotAndSelector << std::endl;
  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = notNotAndSelector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity2->bucket();
    bool result = notNotAndSelector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity3->bucket();
    bool result = notNotAndSelector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity4->bucket();
    bool result = notNotAndSelector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
}

/** \brief Verify flipping the complement bit on an empty Selector
 * works.
 * () does not contain Entity1.
 * !() contains Entity1.
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, complementEmpty ) {
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Selector selector;
  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  selector.complement();
  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
}


/** \brief Verify the fancy output for (PartA OR PartB OR PartC OR
 * PartD).
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, AuBuCuD )
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
  STKUNIT_EXPECT_EQUAL( "!(!PartA AND !PartB AND !PartC AND !PartD)" , msg.str() );
}


/** \brief Verify the fancy output for (PartA AND PartB AND PartC).
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, AiBiC )
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
  STKUNIT_EXPECT_TRUE( msg.str() == "PartA AND PartB AND PartC" );
}


/** \brief Verify the fancy output for a complex expression.
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, complicated )
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
  STKUNIT_EXPECT_EQUAL( "!(!PartA AND !((!(PartA AND PartB) AND !PartC) AND !(PartD AND !PartB)))" , msg.str() );
}


/** \brief Verify \ref stk::mesh::selectIntersection
 * "selectIntersection" works correctly.
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, selectIntersection )
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::PartVector parts ;
  parts.push_back( & fix.m_partA );
  parts.push_back( & fix.m_partB );
  stk::mesh::Selector selector = selectIntersection(parts);

  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = selector(bucket);
    STKUNIT_ASSERT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity2->bucket();
    bool result = selector(bucket);
    STKUNIT_ASSERT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity3->bucket();
    bool result = selector(bucket);
    STKUNIT_ASSERT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity4->bucket();
    bool result = selector(bucket);
    STKUNIT_ASSERT_FALSE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity5->bucket();
    bool result = selector(bucket);
    STKUNIT_ASSERT_FALSE(result);
  }

  std::ostringstream msg;
  msg << selector;
  STKUNIT_EXPECT_TRUE( msg.str() == "PartA AND PartB");
}


/** \brief Verify \ref stk::mesh::selectUnion "selectUnion" works
 * correctly.
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, selectUnion )
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::PartVector parts ;
  parts.push_back( & fix.m_partA );
  parts.push_back( & fix.m_partB );
  parts.push_back( & fix.m_partC );
  stk::mesh::Selector selector = selectUnion(parts);

  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = selector(bucket);
    STKUNIT_ASSERT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity2->bucket();
    bool result = selector(bucket);
    STKUNIT_ASSERT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity3->bucket();
    bool result = selector(bucket);
    STKUNIT_ASSERT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity4->bucket();
    bool result = selector(bucket);
    STKUNIT_ASSERT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity5->bucket();
    bool result = selector(bucket);
    STKUNIT_ASSERT_FALSE(result);
  }

  std::ostringstream msg;
  msg << selector;
  std::cout << "msg.str() = " << msg.str() << std::endl;
  STKUNIT_EXPECT_EQUAL( "!(!PartA AND !PartB AND !PartC)", msg.str() );
}

// Intersection first then union
// & before |
/** \brief Verify order of operations works correctly.
 * (PartA OR PartB AND PartC) = (PartA OR (PartB AND PartC)).
 * (PartB AND PartC OR PartA) = ((PartB AND PartC) OR PartA).
 *
 */
//STKUNIT_UNIT_TEST( UnitTestSelector, orderOfOperations )
//{
//  SelectorFixture fix ;
//  stk::mesh::Part & partA = fix.m_partA ;
//  stk::mesh::Part & partB = fix.m_partB ;
//  stk::mesh::Part & partC = fix.m_partC ;
//  {
//    stk::mesh::Selector selector = partA | partB & partC;
//    //std::cout << "A|B&C selector = " << selector << std::endl;
//    std::ostringstream msg;
//    msg << selector;
//    STKUNIT_EXPECT_EQUAL( "!(!PartA AND !(PartB AND PartC))", msg.str() );
//    {
//      const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
//      bool result = selector(bucket);
//      STKUNIT_EXPECT_TRUE(result);
//    }
//    {
//      const stk::mesh::Bucket & bucket = fix.m_entity2->bucket();
//      bool result = selector(bucket);
//      STKUNIT_EXPECT_TRUE(result);
//    }
//    {
//      const stk::mesh::Bucket & bucket = fix.m_entity3->bucket();
//      bool result = selector(bucket);
//      STKUNIT_EXPECT_TRUE(result);
//    }
//    {
//      const stk::mesh::Bucket & bucket = fix.m_entity4->bucket();
//      bool result = selector(bucket);
//      STKUNIT_EXPECT_FALSE(result);
//    }
//  }
//  {
//    stk::mesh::Selector selector = partB & partC | partA;
//    //std::cout << "B&C|A selector = " << selector << std::endl;
//    std::ostringstream msg;
//    msg << selector;
//    STKUNIT_EXPECT_EQUAL( "!(!(PartB AND PartC) AND !PartA)", msg.str() );
//    {
//      const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
//      bool result = selector(bucket);
//      STKUNIT_EXPECT_TRUE(result);
//    }
//    {
//      const stk::mesh::Bucket & bucket = fix.m_entity2->bucket();
//      bool result = selector(bucket);
//      STKUNIT_EXPECT_TRUE(result);
//    }
//    {
//      const stk::mesh::Bucket & bucket = fix.m_entity3->bucket();
//      bool result = selector(bucket);
//      STKUNIT_EXPECT_TRUE(result);
//    }
//    {
//      const stk::mesh::Bucket & bucket = fix.m_entity4->bucket();
//      bool result = selector(bucket);
//      STKUNIT_EXPECT_FALSE(result);
//    }
//  }
//}


/** \brief Verify unions and intersections of default constructors and
 * their complements.
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, ZeroiuZero ) {
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Selector selectNone;
  stk::mesh::Selector selectAll;
  selectAll.complement();
  {
    stk::mesh::Selector selector = selectNone & selectAll;
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    stk::mesh::Selector selector = selectNone | selectAll;
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = selector(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }

}

/** \brief Verify default \ref stk::mesh::Selector "selectors" work
 * well with \ref stk::mesh::Part "mesh part" instantiated \ref
 * stk::mesh::Selector "selectors".
 *
 * In particular, check that (() OR
 * PartA) contains Entity1 and (!() AND PartA) contains Entity1.
 *
 */

STKUNIT_UNIT_TEST( UnitTestSelector, ZeroiuA )
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA ;
  stk::mesh::Selector selectNone;
  stk::mesh::Selector selectAll;
  selectAll.complement();
  stk::mesh::Selector selectA = partA;
  stk::mesh::Selector selectNoneOrA = selectNone | selectA;
  stk::mesh::Selector selectAllAndA = selectAll & selectA;
  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = selectNoneOrA(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = selectAllAndA(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
}


/** \brief Verify copy constructed \ref stk::mesh::Selector "selector" creates same pretty print output.
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, copyConstructor )
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Part & partA = fix.m_partA ;
  stk::mesh::Part & partB = fix.m_partB ;
  stk::mesh::Part & partC = fix.m_partC ;
  stk::mesh::Selector selectA = (partA & partB) | partC;
  stk::mesh::Selector anotherSelectA(selectA);
  std::ostringstream descriptionA;
  descriptionA << selectA;
  std::ostringstream descriptionAnotherA;
  descriptionAnotherA << anotherSelectA;
  STKUNIT_EXPECT_EQUAL( descriptionA.str() == descriptionAnotherA.str(), true );
}


/** \brief Verify pretty printing of default constructors and their
 * complements works well.
 * In particular !() AND !() and (!() OR !()) == !(() AND ()).
 *
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, AlliuAll )
{
  stk::mesh::Selector selectAll;
  selectAll.complement();

  stk::mesh::Selector anotherSelectAll;
  anotherSelectAll.complement();

  {
    stk::mesh::Selector selectAllANDAll = selectAll & anotherSelectAll;
    std::ostringstream description;
    description << selectAllANDAll;
    STKUNIT_EXPECT_EQUAL( "!() AND !()", description.str() );
  }
  {
    stk::mesh::Selector selectAllORAll = selectAll | anotherSelectAll;
    std::ostringstream description;
    description << selectAllORAll;
    STKUNIT_EXPECT_EQUAL( "!(())", description.str() );
  }
}

/** \brief Verify that the 'selectField' selector works correctly.
 *
 */
STKUNIT_UNIT_TEST( UnitTestSelector, selectField )
{
  SelectorFixture fix ;
  initialize(fix);

  stk::mesh::Selector selectA = stk::mesh::selectField(fix.m_fieldA);
  stk::mesh::Selector selectABC = stk::mesh::selectField(fix.m_fieldABC);
  {
    //entity1 is in partA, so entity1's bucket should be selected by selectA:
    const stk::mesh::Bucket & bucket = fix.m_entity1->bucket();
    bool result = selectA(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
  {
    //entity3 is not in partA, so entity3's bucket should not be selected by selectA:
    const stk::mesh::Bucket & bucket = fix.m_entity3->bucket();
    bool result = selectA(bucket);
    STKUNIT_EXPECT_FALSE(result);
  }
  {
    //entity3 is in partB, so entity3's bucket should be selected by selectABC:
    const stk::mesh::Bucket & bucket = fix.m_entity3->bucket();
    bool result = selectABC(bucket);
    STKUNIT_EXPECT_TRUE(result);
  }
}

/** \} */


} // namespace
