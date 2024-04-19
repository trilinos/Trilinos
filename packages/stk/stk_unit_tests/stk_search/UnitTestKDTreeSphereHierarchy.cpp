/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <gtest/gtest.h>
#include <cmath>
#include <stk_search/kdtree/KDTree_BoundingSpheres.hpp>

namespace {

template <typename NumT>
void test_stk_search_sphere_hier_UnionBoxes()
{
  using namespace stk::search;

  Point<NumT> ptA( 2, 100, 1000);
  Point<NumT> ptB(11, 100, 1000);

  const NumT radiusA = 1;
  const NumT radiusB = 5;

  Sphere<NumT> sphereA(ptA, radiusA);
  Sphere<NumT> sphereB(ptB, radiusB);

  Sphere<NumT> resultSph;
  UnionBoxes(sphereA, sphereB, resultSph);

  Point<NumT> resultCtr = resultSph.center();

  EXPECT_FLOAT_EQ(8.5, resultCtr[0]);
  EXPECT_FLOAT_EQ(100,   resultCtr[1]);
  EXPECT_FLOAT_EQ(1000,   resultCtr[2]);
  EXPECT_FLOAT_EQ(7.5, resultSph.radius());
}

TEST( stk_search_sphere_hier, UnionBoxes_float )
{
  test_stk_search_sphere_hier_UnionBoxes<float>();
}

TEST( stk_search_sphere_hier, UnionBoxes_double )
{
  test_stk_search_sphere_hier_UnionBoxes<double>();
}

TEST( stk_search_sphere_hier, UnionBoxesASwallowsB )
{
  using namespace stk::search;

  Point<float> ptA(2, 100, 1000);
  Point<float> ptB(4, 100, 1000);

  const float radiusA = 5;
  const float radiusB = 1;

  Sphere<float> sphereA(ptA, radiusA);
  Sphere<float> sphereB(ptB, radiusB);

  Sphere<float> resultSph;
  UnionBoxes(sphereA, sphereB, resultSph);

  Point<float> resultCtr = resultSph.center();

  EXPECT_FLOAT_EQ(2, resultCtr[0]);
  EXPECT_FLOAT_EQ(100,   resultCtr[1]);
  EXPECT_FLOAT_EQ(1000,   resultCtr[2]);
  EXPECT_FLOAT_EQ(5, resultSph.radius());
}

TEST( stk_search_sphere_hier, UnionBoxesBSwallowsA )
{
  using namespace stk::search;

  Point<float> ptA(2, 100, 1000);
  Point<float> ptB(4, 100, 1000);

  const float radiusA = 1;
  const float radiusB = 5;

  Sphere<float> sphereA(ptA, radiusA);
  Sphere<float> sphereB(ptB, radiusB);

  Sphere<float> resultSph;
  UnionBoxes(sphereA, sphereB, resultSph);

  Point<float> resultCtr = resultSph.center();

  EXPECT_FLOAT_EQ(4, resultCtr[0]);
  EXPECT_FLOAT_EQ(100,   resultCtr[1]);
  EXPECT_FLOAT_EQ(1000,   resultCtr[2]);
  EXPECT_FLOAT_EQ(5, resultSph.radius());
}

std::vector<stk::search::Sphere<float> > makeTestSphereArray(const int numSpheres_X, const int numSpheres_Y,
                                                             const int numSpheres_Z) {
  using namespace stk::search;
  typedef float                             num_type;
  typedef Sphere<num_type>                  sphere_type;
  typedef typename sphere_type::point_type  point_type;

  const num_type startX = 0;
  const num_type startY = 0;
  const num_type startZ = 0;

  const num_type spacing      = 2.0;
  const num_type sphereRadius = 0.5;

  std::vector<sphere_type> spheres;
  num_type curr_z = startZ;
  for (int i_z = 0; i_z < numSpheres_Z; ++i_z, curr_z += spacing) {
    num_type curr_y = startY;
    for (int i_y = 0; i_y < numSpheres_Y; ++i_y, curr_y += spacing) {
      num_type curr_x = startX;
      for (int i_x = 0; i_x < numSpheres_X; ++i_x, curr_x += spacing) {
        point_type ctr(curr_x, curr_y, curr_z);
        spheres.push_back(sphere_type(ctr, sphereRadius));
      }
    }
  }

  return spheres;
}



TEST( stk_search_sphere_hier, KDTreeSpheresFromSpheres) {
  using namespace stk::search;
  typedef float                             num_type;
  typedef Sphere<num_type>                  sphere_type;

  const int numSpheres_X = 4;
  const int numSpheres_Y = 3;
  const int numSpheres_Z = 2;

  std::vector<sphere_type> spheres =  makeTestSphereArray(numSpheres_X, numSpheres_Y, numSpheres_Z);
  ProximitySearchTree_T<sphere_type> sphereTree;
  makeSphereNodalKDTree(spheres, sphereTree);

  int numSpheres = spheres.size();
  for (int i = 0; i < numSpheres; ++i) {
    std::vector<int> returnIndexList;
    std::vector<sphere_type> returnBoxList;
    sphereTree.SearchForOverlap(spheres[i], returnIndexList, returnBoxList);

    EXPECT_EQ(1u, returnIndexList.size());
    EXPECT_EQ(1u, returnBoxList.size());
    EXPECT_EQ(i, returnIndexList[0]);
    EXPECT_TRUE(spheres[i] == returnBoxList[0]);
  }
}


TEST( stk_search_sphere_hier, TightenedKDTreeSpheresFromSpheres) {
  using namespace stk::search;
  typedef float                             num_type;
  typedef Sphere<num_type>                  sphere_type;
  typedef typename sphere_type::point_type  point_type;

  const int numSpheres_X = 4;
  const int numSpheres_Y = 3;
  const int numSpheres_Z = 2;

  std::vector<sphere_type> spheres =  makeTestSphereArray(numSpheres_X, numSpheres_Y, numSpheres_Z);
  ProximitySearchTree_T<sphere_type> sphereTree;
  makeAABBTightenedSphereKDTree(spheres, sphereTree);

  point_type firstCtr = spheres.front().center();
  point_type lastCtr = spheres.back().center();
  num_type cx = 0.5*(firstCtr[0] + lastCtr[0]);
  num_type cy = 0.5*(firstCtr[1] + lastCtr[1]);
  num_type cz = 0.5*(firstCtr[2] + lastCtr[2]);

  sphere_type rootSphere = sphereTree.BoundingBox();
  point_type rootCtr = rootSphere.center();
  EXPECT_FLOAT_EQ(rootCtr[0], cx);
  EXPECT_FLOAT_EQ(rootCtr[1], cy);
  EXPECT_FLOAT_EQ(rootCtr[2], cz);

  num_type leafRadius = spheres.front().radius();
  num_type rx = 0.5*(lastCtr[0] - firstCtr[0]) + leafRadius;
  num_type ry = 0.5*(lastCtr[1] - firstCtr[1]) + leafRadius;
  num_type rz = 0.5*(lastCtr[2] - firstCtr[2]) + leafRadius;

  num_type cr = sqrt(rx*rx + ry*ry + rz*rz);
  EXPECT_FLOAT_EQ(rootSphere.radius(), cr);

  int numSpheres = spheres.size();
  for (int i = 0; i < numSpheres; ++i) {
    std::vector<int> returnIndexList;
    std::vector<sphere_type> returnBoxList;
    sphereTree.SearchForOverlap(spheres[i], returnIndexList, returnBoxList);

    EXPECT_EQ(1u, returnIndexList.size());
    EXPECT_EQ(1u, returnBoxList.size());
    EXPECT_EQ(i, returnIndexList[0]);
    EXPECT_TRUE(spheres[i] == returnBoxList[0]);
  }
}

}
