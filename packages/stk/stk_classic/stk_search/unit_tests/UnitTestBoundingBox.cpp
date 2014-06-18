/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <mpi.h>

#include <sstream>

#include <stk_search/BoundingBoxCompare.hpp>
#include <stk_search/BoundingBox.hpp>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>
#include <stk_search_util/stk_mesh/PrintBoundingBox.hpp>

namespace stk_classic {
namespace search {
namespace box {

void tests() {
  float center[3] = {0.0f, 0.0f, 0.0f};
  float one[3] = {1.0f, 1.0f, 1.0f};
  float radius = 1.0f;
  float box[6] = {-1.0f,-1.0f,-1.0f,1.0f,1.0f,1.0f};

  //test constructors
  PointBoundingBox<uint64_t,float,3>        pa,  pb(center,1);
  SphereBoundingBox<uint64_t,float,3>       sa,  sb(center,radius,1);
  AxisAlignedBoundingBox<uint64_t,float,3>  aa,  ab(box,1);

  std::vector<PointBoundingBox<uint64_t,float,3> >       paa(3);
  std::vector<SphereBoundingBox<uint64_t,float,3> >       saa(3);
  std::vector<AxisAlignedBoundingBox<uint64_t,float,3> >  aaa(3);

  //test assignment
  pa = pb;
  sa = sb;
  aa = ab;

  //test copy constructors
  PointBoundingBox<uint64_t,float,3> pc(pa);
  SphereBoundingBox<uint64_t,float,3> sc(sa);
  AxisAlignedBoundingBox<uint64_t,float,3> ac(aa);

  SphereBoundingBox<uint64_t,float,3> sd(pa);
  SphereBoundingBox<uint64_t,float,3> se(aa);

  AxisAlignedBoundingBox<uint64_t,float,3> ad(pa);
  AxisAlignedBoundingBox<uint64_t,float,3> ae(sa);

  //test member functions
  for (int i = 0; i<3; ++i) {
    STKUNIT_ASSERT_EQUAL( pa.lower(i),  center[i]);
    STKUNIT_ASSERT_EQUAL( pa.middle(i), center[i]);
    STKUNIT_ASSERT_EQUAL( pa.upper(i),  center[i]);
    STKUNIT_ASSERT_EQUAL( pa.length(i),  0.0f);

    STKUNIT_ASSERT_EQUAL( sa.lower(i),  center[i] - radius);
    STKUNIT_ASSERT_EQUAL( sa.middle(i), center[i] );
    STKUNIT_ASSERT_EQUAL( sa.upper(i),  center[i] + radius );
    STKUNIT_ASSERT_EQUAL( sa.length(i),  2.0f * radius);

    STKUNIT_ASSERT_EQUAL( sa.lower(i),  box[i]);
    STKUNIT_ASSERT_EQUAL( sa.middle(i), box[i] + (box[i+3]-box[i])/2.0f);
    STKUNIT_ASSERT_EQUAL( sa.upper(i),  box[i+3] );
    STKUNIT_ASSERT_EQUAL( sa.length(i),  box[i+3]-box[i]);
  }

  STKUNIT_ASSERT_EQUAL( pa.intersect(pb), true);
  STKUNIT_ASSERT_EQUAL( pa.intersect(sb), true);
  STKUNIT_ASSERT_EQUAL( pa.intersect(ab), true);

  STKUNIT_ASSERT_EQUAL( sa.intersect(pb), true);
  STKUNIT_ASSERT_EQUAL( sa.intersect(sb), true);
  STKUNIT_ASSERT_EQUAL( sa.intersect(ab), true);

  STKUNIT_ASSERT_EQUAL( aa.intersect(pb), true);
  STKUNIT_ASSERT_EQUAL( aa.intersect(sb), true);
  STKUNIT_ASSERT_EQUAL( aa.intersect(ab), true);


  //test bounding box compare
  PointBoundingBox<uint64_t,float,3> a(center,0), b(one,1);

  compare::Compare<PointBoundingBox<uint64_t,float,3>,compare::LOWER> compare_lower(0);
  STKUNIT_ASSERT_EQUAL(compare_lower(a,b),true);

  compare::Compare<PointBoundingBox<uint64_t,float,3>,compare::MIDDLE> compare_middle(0);
  STKUNIT_ASSERT_EQUAL(compare_middle(a,b),true);

  compare::Compare<PointBoundingBox<uint64_t,float,3>,compare::UPPER> compare_upper(0);
  STKUNIT_ASSERT_EQUAL(compare_upper(a,b),true);

  compare::Compare<PointBoundingBox<uint64_t,float,3>,compare::LENGTH> compare_length(0);
  STKUNIT_ASSERT_EQUAL(compare_length(a,b),false);
  STKUNIT_ASSERT_EQUAL(compare_length(b,a),false);

  compare::Compare<PointBoundingBox<uint64_t,float,3>,compare::KEY> compare_key;
  STKUNIT_ASSERT_EQUAL(compare_key(a,b),true);

  compare::Partition<PointBoundingBox<uint64_t,float,3>,compare::LOWER> partition_lower(0,2.0f);
  STKUNIT_ASSERT_EQUAL(partition_lower(a),true);

  compare::Partition<PointBoundingBox<uint64_t,float,3>,compare::MIDDLE> partition_middle(0,2.0f);
  STKUNIT_ASSERT_EQUAL(partition_middle(a),true);

  compare::Partition<PointBoundingBox<uint64_t,float,3>,compare::UPPER> partition_upper(0,2.0f);
  STKUNIT_ASSERT_EQUAL(partition_upper(a),true);

  compare::Partition<PointBoundingBox<uint64_t,float,3>,compare::LENGTH> partition_length(0,2.0f);
  STKUNIT_ASSERT_EQUAL(partition_length(a),true);

  compare::Partition<PointBoundingBox<uint64_t,float,3>,compare::KEY> partition_key(2);
  STKUNIT_ASSERT_EQUAL(partition_key(a),true);

 // use_case::dw() << "Test diag writer for BoundingBox" << std::endl;
 // use_case::dw() << pa << std::endl;
 // use_case::dw() << sa << std::endl;
 // use_case::dw() << aa << std::endl;


  std::stringstream output;

  output << "Test writer for BoundingBox" << std::endl;
  output  << pa << std::endl;
  output << sa << std::endl;
  output << aa << std::endl;

  output << paa << std::endl;
  output << saa << std::endl;
  output << aaa << std::endl;

}

} // namespace box
} // namespace search
} // namespace stk_classic

STKUNIT_UNIT_TEST(UnitTestingOfSearchBox, testUnit)
{
  MPI_Barrier( MPI_COMM_WORLD );
  stk_classic::search::box::tests();
}
