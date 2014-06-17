/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <mpi.h>

#include <stk_search/BihTree.hpp>
#include <stk_search/BoundingBox.hpp>

typedef stk_classic::search::box::AxisAlignedBoundingBox<int,float,3> BoundingVolume;

namespace stk_classic {
namespace search {
namespace bih {

class UnitTestBihTree {
public:
  UnitTestBihTree() {}
  ~UnitTestBihTree() {}
  void testNode();
  void testTree();
};


void UnitTestBihTree::testNode()
{
  BihTree<BoundingVolume>::Node test;

  STKUNIT_ASSERT_EQUAL( test.is_leaf() , true);
  STKUNIT_ASSERT_EQUAL( test.size() , 0);

  //make the node an interior node
  float clip_planes[4] = {-1,1};
  test.make_interior_node( 0 , clip_planes);

  //the num of items in an interior node is -1
  STKUNIT_ASSERT_EQUAL( test.size() , -1);

  STKUNIT_ASSERT_EQUAL( test.is_leaf() , false);
  STKUNIT_ASSERT( test.axis() == 0 );
  STKUNIT_ASSERT_EQUAL( test.left_max(), -1.0f);
  STKUNIT_ASSERT_EQUAL( test.right_min(), 1.0f);

  //check that the axis is being encoded in the pointer correctly
  int axis = static_cast<int>(test.m_children_ptr - reinterpret_cast<uint64_t>(test.children()));
  STKUNIT_ASSERT_EQUAL( axis , test.axis());


  //test leaf node
  stk_classic::search::bih::BihTree<BoundingVolume>::Node & leftChild  = test.left();
  stk_classic::search::bih::BihTree<BoundingVolume>::Node & rightChild  = test.right();


  leftChild.make_leaf_node(0,0);

  STKUNIT_ASSERT_EQUAL( leftChild.is_leaf() , true);
  STKUNIT_ASSERT_EQUAL( leftChild.size() , 0);

  std::vector<BoundingVolume> m_test_keys(4);


  for (unsigned i=0; i<4; ++i) {
   m_test_keys[i].key = i;
  }

  //make empty left child
  rightChild.make_leaf_node(&m_test_keys[0],4);

  STKUNIT_ASSERT_EQUAL( rightChild.is_leaf() , true);
  STKUNIT_ASSERT_EQUAL( rightChild.size() , 4);

}

void test_tree(
  std::vector<BoundingVolume> & range,
  const std::vector<BoundingVolume> & domain,
  unsigned max_depth,
  unsigned items_per_leaf)
{
  stk_classic::search::bih::BihTree<BoundingVolume> tree(range.begin(),range.end(),max_depth,items_per_leaf);

  int extra_intersections = 0;

  std::vector< std::pair< BoundingVolume::Key, BoundingVolume::Key> > relations_vector;
  std::set< std::pair< BoundingVolume::Key, BoundingVolume::Key> > relations_set;
  std::vector<BoundingVolume> range_overlap;

  for (std::vector<BoundingVolume>::const_iterator i = domain.begin(); i != domain.end(); ++i) {
    relations_vector.clear();
    relations_set.clear();
    range_overlap.clear();

    tree.intersect(*i,relations_vector);
    tree.intersect(*i,relations_set);
    tree.intersect(*i,range_overlap);

    extra_intersections += range_overlap.size() - 64;

    STKUNIT_ASSERT(relations_vector.size() >= 64);
    STKUNIT_ASSERT(relations_set.size() >= 64);
    STKUNIT_ASSERT(range_overlap.size() >= 64);
  }

  std::set< std::pair< BoundingVolume::Key, BoundingVolume::Key> > self_relations_set;
  //test for self intersection
  tree.self_intersect(self_relations_set);
}

void UnitTestBihTree::testTree( )
{
  std::vector<BoundingVolume> range;

  const int END = 25;

  for (int i = 0; i<END; ++i) {
    for (int j = 0; j<END; ++j) {
      for (int k = 0; k<END; ++k) {
        BoundingVolume box;
        box.key = 16*i+4*j+k;
        box.box[0] = i-1; box.box[1] = j-1; box.box[2] = k-1;
        box.box[3] = i+1; box.box[4] = j+1; box.box[5] = k+1;
        range.push_back(box);
      }
    }
  }

  std::vector<BoundingVolume> domain;
  for (int i = 2; i<END-2; i += 2) {
    for (int j = 2; j<END-2; j += 2) {
      for (int k = 2; k<END-2; k += 2) {
        BoundingVolume box;
        box.key = -16*i+4*j+k;
        box.box[0] = i-1.5; box.box[1] = j-1.5; box.box[2] = k-1.5;
        box.box[3] = i+0.5; box.box[4] = j+0.5; box.box[5] = k+0.5;
        domain.push_back(box);
      }
    }
  }

  for (int i = 1<<6; i >= 1; i >>= 1) {
    test_tree(range,domain,64,i);
  }

  {
    //test bih tree with single box
    stk_classic::search::bih::BihTree<BoundingVolume> tree(range.begin(),range.begin()+1);
  }

  {
    //test a very unbalanced tree
    std::vector<BoundingVolume> unbalanced;

    for (int i=0; i<END; ++i) {
      BoundingVolume box;
      box.key = i;
      box.box[0] = -1; box.box[1] = -1; box.box[2] = -1;
      box.box[3] = 1; box.box[4] = 1; box.box[4] = 1;
      unbalanced.push_back(box);
    }

    {
      BoundingVolume box;
      box.key = 25;
      box.box[0] = 99; box.box[1] = 99; box.box[2] = 99;
      box.box[3] = 100; box.box[4] = 100; box.box[4] = 100;
      unbalanced.push_back(box);
    }
    {
      BoundingVolume box;
      box.key = 26;
      box.box[0] = 0; box.box[1] = 0; box.box[2] = 0;
      box.box[3] = 100; box.box[4] = 100; box.box[4] = 100;
      unbalanced.push_back(box);
    }

    stk_classic::search::bih::BihTree<BoundingVolume> tree(unbalanced.begin(),unbalanced.end());

    std::set< std::pair< BoundingVolume::Key, BoundingVolume::Key> > self_relations_set;
    //test for self intersection
    tree.self_intersect(self_relations_set);
  }
}

} // namespace bih
} // namespace search
} // namespace stk_classic

STKUNIT_UNIT_TEST(UnitTestingOfSearchBih, testUnit)
{
  MPI_Barrier( MPI_COMM_WORLD );
  stk_classic::search::bih::UnitTestBihTree unitTest;
  unitTest.testNode();
  unitTest.testTree();
}

