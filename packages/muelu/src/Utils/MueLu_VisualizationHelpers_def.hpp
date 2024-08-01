// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_VISUALIZATIONHELPERS_DEF_HPP_
#define MUELU_VISUALIZATIONHELPERS_DEF_HPP_

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include "MueLu_VisualizationHelpers_decl.hpp"
#include "MueLu_Level.hpp"

#include <vector>
#include <list>
#include <algorithm>
#include <string>

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<ParameterList> VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<std::string>("visualization: output filename", "viz%LEVELID", "filename for VTK-style visualization output");
  validParamList->set<int>("visualization: output file: time step", 0, "time step variable for output file name");       // Remove me?
  validParamList->set<int>("visualization: output file: iter", 0, "nonlinear iteration variable for output file name");  // Remove me?
  validParamList->set<std::string>("visualization: style", "Point Cloud", "style of aggregate visualization for VTK output. Can be 'Point Cloud', 'Jacks', 'Convex Hulls'");
  validParamList->set<bool>("visualization: build colormap", false, "Whether to build a random color map in a separate xml file.");
  validParamList->set<bool>("visualization: fine graph edges", false, "Whether to draw all fine node connections along with the aggregates.");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doPointCloud(std::vector<int>& vertices, std::vector<int>& geomSizes, LO /* numLocalAggs */, LO numFineNodes) {
  vertices.reserve(numFineNodes);
  geomSizes.reserve(numFineNodes);
  for (LocalOrdinal i = 0; i < numFineNodes; i++) {
    vertices.push_back(i);
    geomSizes.push_back(1);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doJacks(std::vector<int>& vertices, std::vector<int>& geomSizes, LO numLocalAggs, LO numFineNodes, const std::vector<bool>& isRoot, const ArrayRCP<LO>& vertex2AggId) {
  // For each aggregate, find the root node then connect it with the other nodes in the aggregate
  // Doesn't matter the order, as long as all edges are found.
  vertices.reserve(vertices.size() + 3 * (numFineNodes - numLocalAggs));
  geomSizes.reserve(vertices.size() + 2 * (numFineNodes - numLocalAggs));
  int root = 0;
  for (int i = 0; i < numLocalAggs; i++)  // TODO: Replace this O(n^2) with a better way
  {
    while (!isRoot[root])
      root++;
    int numInAggFound = 0;
    for (int j = 0; j < numFineNodes; j++) {
      if (j == root)  // don't make a connection from the root to itself
      {
        numInAggFound++;
        continue;
      }
      if (vertex2AggId[root] == vertex2AggId[j]) {
        vertices.push_back(root);
        vertices.push_back(j);
        geomSizes.push_back(2);
        // Also draw the free endpoint explicitly for the current line
        vertices.push_back(j);
        geomSizes.push_back(1);
        numInAggFound++;
        // if(numInAggFound == aggSizes_[vertex2AggId_[root]]) //don't spend more time looking if done with that root
        //   break;
      }
    }
    root++;  // get set up to look for the next root
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doConvexHulls2D(std::vector<int>& vertices, std::vector<int>& geomSizes, LO numLocalAggs, LO numFineNodes, const std::vector<bool>& /* isRoot */, const ArrayRCP<LO>& vertex2AggId, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& xCoords, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& yCoords) {
  // This algorithm is based on Andrew's Monotone Chain variant of the Graham Scan for Convex Hulls.  It adds
  // a colinearity check which we'll need for our viz.
  for (int agg = 0; agg < numLocalAggs; agg++) {
    std::vector<int> aggNodes;
    for (int i = 0; i < numFineNodes; i++) {
      if (vertex2AggId[i] == agg)
        aggNodes.push_back(i);
    }

    // have a list of nodes in the aggregate
    TEUCHOS_TEST_FOR_EXCEPTION(aggNodes.size() == 0, Exceptions::RuntimeError,
                               "CoarseningVisualization::doConvexHulls2D: aggregate contains zero nodes!");
    if (aggNodes.size() == 1) {
      vertices.push_back(aggNodes.front());
      geomSizes.push_back(1);
      continue;
    }
    if (aggNodes.size() == 2) {
      vertices.push_back(aggNodes.front());
      vertices.push_back(aggNodes.back());
      geomSizes.push_back(2);
      continue;
    } else {
      int N        = (int)aggNodes.size();
      using MyPair = std::pair<myVec2, int>;
      std::vector<MyPair> pointsAndIndex(N);
      for (int i = 0; i < N; i++) {
        pointsAndIndex[i] = std::make_pair(myVec2(xCoords[aggNodes[i]], yCoords[aggNodes[i]]), i);
      }

      // Sort by x coordinate
      std::sort(pointsAndIndex.begin(), pointsAndIndex.end(), [](const MyPair& a, const MyPair& b) {
        return a.first.x < b.first.x || (a.first.x == b.first.x && a.first.y < b.first.y);
      });

      // Colinearity check
      bool colinear = true;
      for (int i = 0; i < N; i++) {
        if (ccw(pointsAndIndex[i].first, pointsAndIndex[(i + 1) % N].first, pointsAndIndex[(i + 2) % N].first) != 0) {
          colinear = false;
          break;
        }
      }

      if (colinear) {
        vertices.push_back(aggNodes[pointsAndIndex.front().second]);
        vertices.push_back(aggNodes[pointsAndIndex.back().second]);
        geomSizes.push_back(2);
      } else {
        std::vector<int> hull(2 * N);
        int count = 0;

        // Build lower hull
        for (int i = 0; i < N; i++) {
          while (count >= 2 && ccw(pointsAndIndex[hull[count - 2]].first, pointsAndIndex[hull[count - 1]].first, pointsAndIndex[i].first) <= 0) {
            count--;
          }
          hull[count++] = i;
        }

        // Build the upper hull
        int t = count + 1;
        for (int i = N - 1; i > 0; i--) {
          while (count >= t && ccw(pointsAndIndex[hull[count - 2]].first, pointsAndIndex[hull[count - 1]].first, pointsAndIndex[i - 1].first) <= 0) {
            count--;
          }
          hull[count++] = i - 1;
        }

        // Remove the duplicated point
        hull.resize(count - 1);

        // Verify: Check that hull retains CCW order
        for (int i = 0; i < (int)hull.size(); i++) {
          TEUCHOS_TEST_FOR_EXCEPTION(ccw(pointsAndIndex[hull[i]].first, pointsAndIndex[hull[(i + 1) % hull.size()]].first, pointsAndIndex[hull[(i + 1) % hull.size()]].first) == 1, Exceptions::RuntimeError, "CoarseningVisualization::doConvexHulls2D: counter-clockwise verification fails");
        }

        // We now need to undo the indices from the sorting
        for (int i = 0; i < (int)hull.size(); i++) {
          hull[i] = aggNodes[pointsAndIndex[hull[i]].second];
        }

        vertices.reserve(vertices.size() + hull.size());
        for (size_t i = 0; i < hull.size(); i++) {
          vertices.push_back(hull[i]);
        }
        geomSizes.push_back(hull.size());
      }  // else colinear
    }    // else 3 + nodes
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doConvexHulls3D(std::vector<int>& vertices, std::vector<int>& geomSizes, LO numLocalAggs, LO numFineNodes, const std::vector<bool>& /* isRoot */, const ArrayRCP<LO>& vertex2AggId, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& xCoords, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& yCoords, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& zCoords) {
  // Use 3D quickhull algo.
  // Vector of node indices representing triangle vertices
  // Note: Calculate the hulls first since will only include point data for points in the hulls
  // Effectively the size() of vertIndices after each hull is added to it
  typedef std::list<int>::iterator Iter;
  for (int agg = 0; agg < numLocalAggs; agg++) {
    std::list<int> aggNodes;  // At first, list of all nodes in the aggregate. As nodes are enclosed or included by/in hull, remove them
    for (int i = 0; i < numFineNodes; i++) {
      if (vertex2AggId[i] == agg)
        aggNodes.push_back(i);
    }
    // First, check anomalous cases
    TEUCHOS_TEST_FOR_EXCEPTION(aggNodes.size() == 0, Exceptions::RuntimeError,
                               "CoarseningVisualization::doConvexHulls3D: aggregate contains zero nodes!");
    if (aggNodes.size() == 1) {
      vertices.push_back(aggNodes.front());
      geomSizes.push_back(1);
      continue;
    } else if (aggNodes.size() == 2) {
      vertices.push_back(aggNodes.front());
      vertices.push_back(aggNodes.back());
      geomSizes.push_back(2);
      continue;
    }
    // check for collinearity
    bool areCollinear = true;
    {
      Iter it = aggNodes.begin();
      myVec3 firstVec(xCoords[*it], yCoords[*it], zCoords[*it]);
      myVec3 comp;
      {
        it++;
        myVec3 secondVec(xCoords[*it], yCoords[*it], zCoords[*it]);  // cross this with other vectors to compare
        comp = vecSubtract(secondVec, firstVec);
        it++;
      }
      for (; it != aggNodes.end(); it++) {
        myVec3 thisVec(xCoords[*it], yCoords[*it], zCoords[*it]);
        myVec3 cross = crossProduct(vecSubtract(thisVec, firstVec), comp);
        if (mymagnitude(cross) > 1e-10) {
          areCollinear = false;
          break;
        }
      }
    }
    if (areCollinear) {
      // find the endpoints of segment describing all the points
      // compare x, if tie compare y, if tie compare z
      Iter min = aggNodes.begin();
      Iter max = aggNodes.begin();
      Iter it  = ++aggNodes.begin();
      for (; it != aggNodes.end(); it++) {
        if (xCoords[*it] < xCoords[*min])
          min = it;
        else if (xCoords[*it] == xCoords[*min]) {
          if (yCoords[*it] < yCoords[*min])
            min = it;
          else if (yCoords[*it] == yCoords[*min]) {
            if (zCoords[*it] < zCoords[*min]) min = it;
          }
        }
        if (xCoords[*it] > xCoords[*max])
          max = it;
        else if (xCoords[*it] == xCoords[*max]) {
          if (yCoords[*it] > yCoords[*max])
            max = it;
          else if (yCoords[*it] == yCoords[*max]) {
            if (zCoords[*it] > zCoords[*max])
              max = it;
          }
        }
      }
      vertices.push_back(*min);
      vertices.push_back(*max);
      geomSizes.push_back(2);
      continue;
    }
    bool areCoplanar = true;
    {
      // number of points is known to be >= 3
      Iter vert = aggNodes.begin();
      myVec3 v1(xCoords[*vert], yCoords[*vert], zCoords[*vert]);
      vert++;
      myVec3 v2(xCoords[*vert], yCoords[*vert], zCoords[*vert]);
      vert++;
      myVec3 v3(xCoords[*vert], yCoords[*vert], zCoords[*vert]);
      vert++;
      // Make sure the first three points aren't also collinear (need a non-degenerate triangle to get a normal)
      while (mymagnitude(crossProduct(vecSubtract(v1, v2), vecSubtract(v1, v3))) < 1e-10) {
        // Replace the third point with the next point
        v3 = myVec3(xCoords[*vert], yCoords[*vert], zCoords[*vert]);
        vert++;
      }
      for (; vert != aggNodes.end(); vert++) {
        myVec3 pt(xCoords[*vert], yCoords[*vert], zCoords[*vert]);
        if (fabs(pointDistFromTri(pt, v1, v2, v3)) > 1e-12) {
          areCoplanar = false;
          break;
        }
      }
      if (areCoplanar) {
        // do 2D convex hull
        myVec3 planeNorm = getNorm(v1, v2, v3);
        planeNorm.x      = fabs(planeNorm.x);
        planeNorm.y      = fabs(planeNorm.y);
        planeNorm.z      = fabs(planeNorm.z);
        std::vector<myVec2> points;
        std::vector<int> nodes;
        if (planeNorm.x >= planeNorm.y && planeNorm.x >= planeNorm.z) {
          // project points to yz plane to make hull
          for (Iter it = aggNodes.begin(); it != aggNodes.end(); it++) {
            nodes.push_back(*it);
            points.push_back(myVec2(yCoords[*it], zCoords[*it]));
          }
        }
        if (planeNorm.y >= planeNorm.x && planeNorm.y >= planeNorm.z) {
          // use xz
          for (Iter it = aggNodes.begin(); it != aggNodes.end(); it++) {
            nodes.push_back(*it);
            points.push_back(myVec2(xCoords[*it], zCoords[*it]));
          }
        }
        if (planeNorm.z >= planeNorm.x && planeNorm.z >= planeNorm.y) {
          for (Iter it = aggNodes.begin(); it != aggNodes.end(); it++) {
            nodes.push_back(*it);
            points.push_back(myVec2(xCoords[*it], yCoords[*it]));
          }
        }
        std::vector<int> convhull2d = giftWrap(points, nodes, xCoords, yCoords);
        geomSizes.push_back(convhull2d.size());
        vertices.reserve(vertices.size() + convhull2d.size());
        for (size_t i = 0; i < convhull2d.size(); i++)
          vertices.push_back(convhull2d[i]);
        continue;
      }
    }
    Iter exIt        = aggNodes.begin();                            // iterator to be used for searching for min/max x/y/z
    int extremeSix[] = {*exIt, *exIt, *exIt, *exIt, *exIt, *exIt};  // nodes with minimumX, maxX, minY ...
    exIt++;
    for (; exIt != aggNodes.end(); exIt++) {
      Iter it = exIt;
      if (xCoords[*it] < xCoords[extremeSix[0]] ||
          (xCoords[*it] == xCoords[extremeSix[0]] && yCoords[*it] < yCoords[extremeSix[0]]) ||
          (xCoords[*it] == xCoords[extremeSix[0]] && yCoords[*it] == yCoords[extremeSix[0]] && zCoords[*it] < zCoords[extremeSix[0]]))
        extremeSix[0] = *it;
      if (xCoords[*it] > xCoords[extremeSix[1]] ||
          (xCoords[*it] == xCoords[extremeSix[1]] && yCoords[*it] > yCoords[extremeSix[1]]) ||
          (xCoords[*it] == xCoords[extremeSix[1]] && yCoords[*it] == yCoords[extremeSix[1]] && zCoords[*it] > zCoords[extremeSix[1]]))
        extremeSix[1] = *it;
      if (yCoords[*it] < yCoords[extremeSix[2]] ||
          (yCoords[*it] == yCoords[extremeSix[2]] && zCoords[*it] < zCoords[extremeSix[2]]) ||
          (yCoords[*it] == yCoords[extremeSix[2]] && zCoords[*it] == zCoords[extremeSix[2]] && xCoords[*it] < xCoords[extremeSix[2]]))
        extremeSix[2] = *it;
      if (yCoords[*it] > yCoords[extremeSix[3]] ||
          (yCoords[*it] == yCoords[extremeSix[3]] && zCoords[*it] > zCoords[extremeSix[3]]) ||
          (yCoords[*it] == yCoords[extremeSix[3]] && zCoords[*it] == zCoords[extremeSix[3]] && xCoords[*it] > xCoords[extremeSix[3]]))
        extremeSix[3] = *it;
      if (zCoords[*it] < zCoords[extremeSix[4]] ||
          (zCoords[*it] == zCoords[extremeSix[4]] && xCoords[*it] < xCoords[extremeSix[4]]) ||
          (zCoords[*it] == zCoords[extremeSix[4]] && xCoords[*it] == xCoords[extremeSix[4]] && yCoords[*it] < yCoords[extremeSix[4]]))
        extremeSix[4] = *it;
      if (zCoords[*it] > zCoords[extremeSix[5]] ||
          (zCoords[*it] == zCoords[extremeSix[5]] && xCoords[*it] > xCoords[extremeSix[5]]) ||
          (zCoords[*it] == zCoords[extremeSix[5]] && xCoords[*it] == xCoords[extremeSix[5]] && yCoords[*it] > zCoords[extremeSix[5]]))
        extremeSix[5] = *it;
    }
    myVec3 extremeVectors[6];
    for (int i = 0; i < 6; i++) {
      myVec3 thisExtremeVec(xCoords[extremeSix[i]], yCoords[extremeSix[i]], zCoords[extremeSix[i]]);
      extremeVectors[i] = thisExtremeVec;
    }
    double maxDist = 0;
    int base1      = 0;  // ints from 0-5: which pair out of the 6 extreme points are the most distant? (indices in extremeSix and extremeVectors)
    int base2      = 0;
    for (int i = 0; i < 5; i++) {
      for (int j = i + 1; j < 6; j++) {
        double thisDist = distance(extremeVectors[i], extremeVectors[j]);
        // Want to make sure thisDist is significantly larger than maxDist to prevent roundoff errors from impacting node choice.
        if (thisDist > maxDist && thisDist - maxDist > 1e-12) {
          maxDist = thisDist;
          base1   = i;
          base2   = j;
        }
      }
    }
    std::list<myTriangle> hullBuilding;  // each Triangle is a triplet of nodes (int IDs) that form a triangle
    // remove base1 and base2 iters from aggNodes, they are known to be in the hull
    aggNodes.remove(extremeSix[base1]);
    aggNodes.remove(extremeSix[base2]);
    // extremeSix[base1] and [base2] still have the myVec3 representation
    myTriangle tri;
    tri.v1 = extremeSix[base1];
    tri.v2 = extremeSix[base2];
    // Now find the point that is furthest away from the line between base1 and base2
    maxDist = 0;
    // need the vectors to do "quadruple product" formula
    myVec3 b1 = extremeVectors[base1];
    myVec3 b2 = extremeVectors[base2];
    Iter thirdNode;
    for (Iter node = aggNodes.begin(); node != aggNodes.end(); node++) {
      myVec3 nodePos(xCoords[*node], yCoords[*node], zCoords[*node]);
      double dist = mymagnitude(crossProduct(vecSubtract(nodePos, b1), vecSubtract(nodePos, b2))) / mymagnitude(vecSubtract(b2, b1));
      // Want to make sure dist is significantly larger than maxDist to prevent roundoff errors from impacting node choice.
      if (dist > maxDist && dist - maxDist > 1e-12) {
        maxDist   = dist;
        thirdNode = node;
      }
    }
    // Now know the last node in the first triangle
    tri.v3 = *thirdNode;
    hullBuilding.push_back(tri);
    myVec3 b3(xCoords[*thirdNode], yCoords[*thirdNode], zCoords[*thirdNode]);
    aggNodes.erase(thirdNode);
    // Find the fourth node (most distant from triangle) to form tetrahedron
    maxDist          = 0;
    int fourthVertex = -1;
    for (Iter node = aggNodes.begin(); node != aggNodes.end(); node++) {
      myVec3 thisNode(xCoords[*node], yCoords[*node], zCoords[*node]);
      double nodeDist = pointDistFromTri(thisNode, b1, b2, b3);
      // Want to make sure nodeDist is significantly larger than maxDist to prevent roundoff errors from impacting node choice.
      if (nodeDist > maxDist && nodeDist - maxDist > 1e-12) {
        maxDist      = nodeDist;
        fourthVertex = *node;
      }
    }
    aggNodes.remove(fourthVertex);
    myVec3 b4(xCoords[fourthVertex], yCoords[fourthVertex], zCoords[fourthVertex]);
    // Add three new triangles to hullBuilding to form the first tetrahedron
    // use tri to hold the triangle info temporarily before being added to list
    tri    = hullBuilding.front();
    tri.v1 = fourthVertex;
    hullBuilding.push_back(tri);
    tri    = hullBuilding.front();
    tri.v2 = fourthVertex;
    hullBuilding.push_back(tri);
    tri    = hullBuilding.front();
    tri.v3 = fourthVertex;
    hullBuilding.push_back(tri);
    // now orient all four triangles so that the vertices are oriented clockwise (so getNorm_ points outward for each)
    myVec3 barycenter((b1.x + b2.x + b3.x + b4.x) / 4.0, (b1.y + b2.y + b3.y + b4.y) / 4.0, (b1.z + b2.z + b3.z + b4.z) / 4.0);
    for (std::list<myTriangle>::iterator tetTri = hullBuilding.begin(); tetTri != hullBuilding.end(); tetTri++) {
      myVec3 triVert1(xCoords[tetTri->v1], yCoords[tetTri->v1], zCoords[tetTri->v1]);
      myVec3 triVert2(xCoords[tetTri->v2], yCoords[tetTri->v2], zCoords[tetTri->v2]);
      myVec3 triVert3(xCoords[tetTri->v3], yCoords[tetTri->v3], zCoords[tetTri->v3]);
      myVec3 trinorm   = getNorm(triVert1, triVert2, triVert3);
      myVec3 ptInPlane = tetTri == hullBuilding.begin() ? b1 : b4;  // first triangle definitely has b1 in it, other three definitely have b4
      if (isInFront(barycenter, ptInPlane, trinorm)) {
        // don't want the faces of the tetrahedron to face inwards (towards barycenter) so reverse orientation
        // by swapping two vertices
        int temp   = tetTri->v1;
        tetTri->v1 = tetTri->v2;
        tetTri->v2 = temp;
      }
    }
    // now, have starting polyhedron in hullBuilding (all faces are facing outwards according to getNorm_) and remaining nodes to process are in aggNodes
    // recursively, for each triangle, make a list of the points that are 'in front' of the triangle. Find the point with the maximum distance from the triangle.
    // Add three new triangles, each sharing one edge with the original triangle but now with the most distant point as a vertex. Remove the most distant point from
    // the list of all points that need to be processed. Also from that list remove all points that are in front of the original triangle but not in front of all three
    // new triangles, since they are now enclosed in the hull.
    // Construct point lists for each face of the tetrahedron individually.
    myVec3 trinorms[4];  // normals to the four tetrahedron faces, now oriented outwards
    int index = 0;
    for (std::list<myTriangle>::iterator tetTri = hullBuilding.begin(); tetTri != hullBuilding.end(); tetTri++) {
      myVec3 triVert1(xCoords[tetTri->v1], yCoords[tetTri->v1], zCoords[tetTri->v1]);
      myVec3 triVert2(xCoords[tetTri->v2], yCoords[tetTri->v2], zCoords[tetTri->v2]);
      myVec3 triVert3(xCoords[tetTri->v3], yCoords[tetTri->v3], zCoords[tetTri->v3]);
      trinorms[index] = getNorm(triVert1, triVert2, triVert3);
      index++;
    }
    std::list<int> startPoints1;
    std::list<int> startPoints2;
    std::list<int> startPoints3;
    std::list<int> startPoints4;
    // scope this so that 'point' is not in function scope
    {
      Iter point = aggNodes.begin();
      while (!aggNodes.empty())  // this removes points one at a time as they are put in startPointsN or are already done
      {
        myVec3 pointVec(xCoords[*point], yCoords[*point], zCoords[*point]);
        // Note: Because of the way the tetrahedron faces are constructed above,
        // face 1 definitely contains b1 and faces 2-4 definitely contain b4.
        if (isInFront(pointVec, b1, trinorms[0])) {
          startPoints1.push_back(*point);
          point = aggNodes.erase(point);
        } else if (isInFront(pointVec, b4, trinorms[1])) {
          startPoints2.push_back(*point);
          point = aggNodes.erase(point);
        } else if (isInFront(pointVec, b4, trinorms[2])) {
          startPoints3.push_back(*point);
          point = aggNodes.erase(point);
        } else if (isInFront(pointVec, b4, trinorms[3])) {
          startPoints4.push_back(*point);
          point = aggNodes.erase(point);
        } else {
          point = aggNodes.erase(point);  // points here are already inside tetrahedron.
        }
      }
      // Call processTriangle for each triangle in the initial tetrahedron, one at a time.
    }
    typedef std::list<myTriangle>::iterator TriIter;
    TriIter firstTri  = hullBuilding.begin();
    myTriangle start1 = *firstTri;
    firstTri++;
    myTriangle start2 = *firstTri;
    firstTri++;
    myTriangle start3 = *firstTri;
    firstTri++;
    myTriangle start4 = *firstTri;
    // kick off depth-first recursive filling of hullBuilding list with all triangles in the convex hull
    if (!startPoints1.empty())
      processTriangle(hullBuilding, start1, startPoints1, barycenter, xCoords, yCoords, zCoords);
    if (!startPoints2.empty())
      processTriangle(hullBuilding, start2, startPoints2, barycenter, xCoords, yCoords, zCoords);
    if (!startPoints3.empty())
      processTriangle(hullBuilding, start3, startPoints3, barycenter, xCoords, yCoords, zCoords);
    if (!startPoints4.empty())
      processTriangle(hullBuilding, start4, startPoints4, barycenter, xCoords, yCoords, zCoords);
    // hullBuilding now has all triangles that make up this hull.
    // Dump hullBuilding info into the list of all triangles for the scene.
    vertices.reserve(vertices.size() + 3 * hullBuilding.size());
    for (TriIter hullTri = hullBuilding.begin(); hullTri != hullBuilding.end(); hullTri++) {
      vertices.push_back(hullTri->v1);
      vertices.push_back(hullTri->v2);
      vertices.push_back(hullTri->v3);
      geomSizes.push_back(3);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doGraphEdges(std::vector<int>& vertices, std::vector<int>& geomSizes, Teuchos::RCP<LWGraph>& G, Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& /* fx */, Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& /* fy */, Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& /* fz */) {
  ArrayView<const Scalar> values;

  std::vector<std::pair<int, int> > vert1;  // vertices (node indices)

  ArrayView<const LocalOrdinal> indices;
  for (LocalOrdinal locRow = 0; locRow < LocalOrdinal(G->GetNodeNumVertices()); locRow++) {
    auto neighbors = G->getNeighborVertices(locRow);
    // Add those local indices (columns) to the list of connections (which are pairs of the form (localM, localN))
    for (int gEdge = 0; gEdge < int(neighbors.length); ++gEdge) {
      vert1.push_back(std::pair<int, int>(locRow, neighbors(gEdge)));
    }
  }
  for (size_t i = 0; i < vert1.size(); i++) {
    if (vert1[i].first > vert1[i].second) {
      int temp        = vert1[i].first;
      vert1[i].first  = vert1[i].second;
      vert1[i].second = temp;
    }
  }
  std::sort(vert1.begin(), vert1.end());
  std::vector<std::pair<int, int> >::iterator newEnd = unique(vert1.begin(), vert1.end());  // remove duplicate edges
  vert1.erase(newEnd, vert1.end());
  // std::vector<int> points1;
  vertices.reserve(2 * vert1.size());
  geomSizes.reserve(vert1.size());
  for (size_t i = 0; i < vert1.size(); i++) {
    vertices.push_back(vert1[i].first);
    vertices.push_back(vert1[i].second);
    geomSizes.push_back(2);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ccw(const myVec2& a, const myVec2& b, const myVec2& c) {
  const double ccw_zero_threshold = 1e-8;
  double val                      = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
  return (val > ccw_zero_threshold) ? 1 : ((val < -ccw_zero_threshold) ? -1 : 0);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
myVec3 VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::crossProduct(myVec3 v1, myVec3 v2) {
  return myVec3(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dotProduct(myVec2 v1, myVec2 v2) {
  return v1.x * v2.x + v1.y * v2.y;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dotProduct(myVec3 v1, myVec3 v2) {
  return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isInFront(myVec3 point, myVec3 inPlane, myVec3 n) {
  myVec3 rel(point.x - inPlane.x, point.y - inPlane.y, point.z - inPlane.z);  // position of the point relative to the plane
  return dotProduct(rel, n) > 1e-12 ? true : false;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mymagnitude(myVec2 vec) {
  return sqrt(dotProduct(vec, vec));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mymagnitude(myVec3 vec) {
  return sqrt(dotProduct(vec, vec));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::distance(myVec2 p1, myVec2 p2) {
  return mymagnitude(vecSubtract(p1, p2));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::distance(myVec3 p1, myVec3 p2) {
  return mymagnitude(vecSubtract(p1, p2));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
myVec2 VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::vecSubtract(myVec2 v1, myVec2 v2) {
  return myVec2(v1.x - v2.x, v1.y - v2.y);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
myVec3 VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::vecSubtract(myVec3 v1, myVec3 v2) {
  return myVec3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
myVec2 VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNorm(myVec2 v)  //"normal" to a 2D vector - just rotate 90 degrees to left
{
  return myVec2(v.y, -v.x);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
myVec3 VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNorm(myVec3 v1, myVec3 v2, myVec3 v3)  // normal to face of triangle (will be outward rel. to polyhedron) (v1, v2, v3 are in CCW order when normal is toward viewpoint)
{
  return crossProduct(vecSubtract(v2, v1), vecSubtract(v3, v1));
}

// get minimum distance from 'point' to plane containing v1, v2, v3 (or the triangle with v1, v2, v3 as vertices)
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::pointDistFromTri(myVec3 point, myVec3 v1, myVec3 v2, myVec3 v3) {
  using namespace std;
  myVec3 norm = getNorm(v1, v2, v3);
  // must normalize the normal vector
  double normScl = mymagnitude(norm);
  double rv      = 0.0;
  if (normScl > 1e-8) {
    norm.x /= normScl;
    norm.y /= normScl;
    norm.z /= normScl;
    rv = fabs(dotProduct(norm, vecSubtract(point, v1)));
  } else {
    // triangle is degenerated
    myVec3 test1  = vecSubtract(v3, v1);
    myVec3 test2  = vecSubtract(v2, v1);
    bool useTest1 = mymagnitude(test1) > 0.0 ? true : false;
    bool useTest2 = mymagnitude(test2) > 0.0 ? true : false;
    if (useTest1 == true) {
      double normScl1 = mymagnitude(test1);
      test1.x /= normScl1;
      test1.y /= normScl1;
      test1.z /= normScl1;
      rv = fabs(dotProduct(test1, vecSubtract(point, v1)));
    } else if (useTest2 == true) {
      double normScl2 = mymagnitude(test2);
      test2.x /= normScl2;
      test2.y /= normScl2;
      test2.z /= normScl2;
      rv = fabs(dotProduct(test2, vecSubtract(point, v1)));
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError,
                                 "VisualizationHelpers::pointDistFromTri: Could not determine the distance of a point to a triangle.");
    }
  }
  return rv;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::vector<myTriangle> VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::processTriangle(std::list<myTriangle>& tris, myTriangle tri, std::list<int>& pointsInFront, myVec3& barycenter, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& xCoords, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& yCoords, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& zCoords) {
  //*tri is in the tris list, and is the triangle to process here. tris is a complete list of all triangles in the hull so far. pointsInFront is only a list of the nodes in front of tri. Need coords also.
  // precondition: each triangle is already oriented so that getNorm_(v1, v2, v3) points outward (away from interior of hull)
  // First find the point furthest from triangle.
  using namespace std;
  typedef std::list<int>::iterator Iter;
  typedef std::list<myTriangle>::iterator TriIter;
  typedef list<pair<int, int> >::iterator EdgeIter;
  double maxDist = 0;
  // Need vector representations of triangle's vertices
  myVec3 v1(xCoords[tri.v1], yCoords[tri.v1], zCoords[tri.v1]);
  myVec3 v2(xCoords[tri.v2], yCoords[tri.v2], zCoords[tri.v2]);
  myVec3 v3(xCoords[tri.v3], yCoords[tri.v3], zCoords[tri.v3]);
  myVec3 farPointVec;  // useful to have both the point's coordinates and it's position in the list
  Iter farPoint = pointsInFront.begin();
  for (Iter point = pointsInFront.begin(); point != pointsInFront.end(); point++) {
    myVec3 pointVec(xCoords[*point], yCoords[*point], zCoords[*point]);
    double dist = pointDistFromTri(pointVec, v1, v2, v3);
    // Want to make sure nodeDist is significantly larger than maxDist to prevent roundoff errors from impacting node choice.
    if (dist > maxDist && dist - maxDist > 1e-12) {
      dist        = maxDist;
      farPointVec = pointVec;
      farPoint    = point;
    }
  }
  // Find all the triangles that the point is in front of (can be more than 1)
  // At the same time, remove them from tris, as every one will be replaced later
  vector<myTriangle> visible;  // use a list of iterators so that the underlying object is still in tris
  for (TriIter it = tris.begin(); it != tris.end();) {
    myVec3 vec1(xCoords[it->v1], yCoords[it->v1], zCoords[it->v1]);
    myVec3 vec2(xCoords[it->v2], yCoords[it->v2], zCoords[it->v2]);
    myVec3 vec3(xCoords[it->v3], yCoords[it->v3], zCoords[it->v3]);
    myVec3 norm = getNorm(vec1, vec2, vec3);
    if (isInFront(farPointVec, vec1, norm)) {
      visible.push_back(*it);
      it = tris.erase(it);
    } else
      it++;
  }
  // Figure out what triangles need to be destroyed/created
  // First create a list of edges (as std::pair<int, int>, where the two ints are the node endpoints)
  list<pair<int, int> > horizon;
  // For each triangle, add edges to the list iff the edge only appears once in the set of all
  // Have members of horizon have the lower node # first, and the higher one second
  for (vector<myTriangle>::iterator it = visible.begin(); it != visible.end(); it++) {
    pair<int, int> e1(it->v1, it->v2);
    pair<int, int> e2(it->v2, it->v3);
    pair<int, int> e3(it->v1, it->v3);
    //"sort" the pair values
    if (e1.first > e1.second) {
      int temp  = e1.first;
      e1.first  = e1.second;
      e1.second = temp;
    }
    if (e2.first > e2.second) {
      int temp  = e2.first;
      e2.first  = e2.second;
      e2.second = temp;
    }
    if (e3.first > e3.second) {
      int temp  = e3.first;
      e3.first  = e3.second;
      e3.second = temp;
    }
    horizon.push_back(e1);
    horizon.push_back(e2);
    horizon.push_back(e3);
  }
  // sort based on lower node first, then higher node (lexicographical ordering built in to pair)
  horizon.sort();
  // Remove all edges from horizon, except those that appear exactly once
  {
    EdgeIter it = horizon.begin();
    while (it != horizon.end()) {
      int occur = count(horizon.begin(), horizon.end(), *it);
      if (occur > 1) {
        pair<int, int> removeVal = *it;
        while (removeVal == *it && !(it == horizon.end()))
          it = horizon.erase(it);
      } else
        it++;
    }
  }
  // Now make a list of new triangles being created, each of which take 2 vertices from an edge and one from farPoint
  list<myTriangle> newTris;
  for (EdgeIter it = horizon.begin(); it != horizon.end(); it++) {
    myTriangle t(it->first, it->second, *farPoint);
    newTris.push_back(t);
  }
  // Ensure every new triangle is oriented outwards, using the barycenter of the initial tetrahedron
  vector<myTriangle> trisToProcess;
  vector<list<int> > newFrontPoints;
  for (TriIter it = newTris.begin(); it != newTris.end(); it++) {
    myVec3 t1(xCoords[it->v1], yCoords[it->v1], zCoords[it->v1]);
    myVec3 t2(xCoords[it->v2], yCoords[it->v2], zCoords[it->v2]);
    myVec3 t3(xCoords[it->v3], yCoords[it->v3], zCoords[it->v3]);
    if (isInFront(barycenter, t1, getNorm(t1, t2, t3))) {
      // need to swap two vertices to flip orientation of triangle
      int temp       = it->v1;
      myVec3 tempVec = t1;
      it->v1         = it->v2;
      t1             = t2;
      it->v2         = temp;
      t2             = tempVec;
    }
    myVec3 outwardNorm = getNorm(t1, t2, t3);  // now definitely points outwards
    // Add the triangle to tris
    tris.push_back(*it);
    trisToProcess.push_back(tris.back());
    // Make a list of the points that are in front of nextToProcess, to be passed in for processing
    list<int> newInFront;
    for (Iter point = pointsInFront.begin(); point != pointsInFront.end();) {
      myVec3 pointVec(xCoords[*point], yCoords[*point], zCoords[*point]);
      if (isInFront(pointVec, t1, outwardNorm)) {
        newInFront.push_back(*point);
        point = pointsInFront.erase(point);
      } else
        point++;
    }
    newFrontPoints.push_back(newInFront);
  }
  vector<myTriangle> allRemoved;  // list of all invalid iterators that were erased by calls to processmyTriangle below
  for (int i = 0; i < int(trisToProcess.size()); i++) {
    if (!newFrontPoints[i].empty()) {
      // Comparing the 'triangle to process' to the one for this call prevents infinite recursion/stack overflow.
      // TODO: Why was it doing that? Rounding error? Make more robust fix. But this does work for the time being.
      if (find(allRemoved.begin(), allRemoved.end(), trisToProcess[i]) == allRemoved.end() && !(trisToProcess[i] == tri)) {
        vector<myTriangle> removedList = processTriangle(tris, trisToProcess[i], newFrontPoints[i], barycenter, xCoords, yCoords, zCoords);
        for (int j = 0; j < int(removedList.size()); j++)
          allRemoved.push_back(removedList[j]);
      }
    }
  }
  return visible;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::vector<int> VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::giftWrap(std::vector<myVec2>& points, std::vector<int>& nodes, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& xCoords, const Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& yCoords) {
  TEUCHOS_TEST_FOR_EXCEPTION(points.size() < 3, Exceptions::RuntimeError,
                             "CoarseningVisualization::giftWrap: Gift wrap algorithm input has to have at least 3 points!");

#if 1  // TAW's version to determine "minimal" node
  // determine minimal x and y coordinates
  double min_x = points[0].x;
  double min_y = points[0].y;
  for (std::vector<int>::iterator it = nodes.begin(); it != nodes.end(); it++) {
    int i = it - nodes.begin();
    if (points[i].x < min_x) min_x = points[i].x;
    if (points[i].y < min_y) min_y = points[i].y;
  }
  // create dummy min coordinates
  min_x -= 1.0;
  min_y -= 1.0;
  myVec2 dummy_min(min_x, min_y);

  // loop over all nodes and determine nodes with minimal distance to (min_x, min_y)
  std::vector<int> hull;
  myVec2 min                         = points[0];
  double mindist                     = distance(min, dummy_min);
  std::vector<int>::iterator minNode = nodes.begin();
  for (std::vector<int>::iterator it = nodes.begin(); it != nodes.end(); it++) {
    int i = it - nodes.begin();
    if (distance(points[i], dummy_min) < mindist) {
      mindist = distance(points[i], dummy_min);
      min     = points[i];
      minNode = it;
    }
  }
  hull.push_back(*minNode);
#else  // Brian's code
  std::vector<int> hull;
  std::vector<int>::iterator minNode = nodes.begin();
  myVec2 min                         = points[0];
  for (std::vector<int>::iterator it = nodes.begin(); it != nodes.end(); it++) {
    int i = it - nodes.begin();
    if (points[i].x < min.x || (fabs(points[i].x - min.x) < 1e-10 && points[i].y < min.y)) {
      min     = points[i];
      minNode = it;
    }
  }
  hull.push_back(*minNode);
#endif

  bool includeMin = false;
  // int debug_it = 0;
  while (1) {
    std::vector<int>::iterator leftMost = nodes.begin();
    if (!includeMin && leftMost == minNode) {
      leftMost++;
    }
    std::vector<int>::iterator it = leftMost;
    it++;
    for (; it != nodes.end(); it++) {
      if (it == minNode && !includeMin)  // don't compare to min on very first sweep
        continue;
      if (*it == hull.back())
        continue;
      // see if it is in front of line containing nodes thisHull.back() and leftMost
      // first get the left normal of leftMost - thisHull.back() (<dy, -dx>)
      myVec2 leftMostVec = points[leftMost - nodes.begin()];
      myVec2 lastVec(xCoords[hull.back()], yCoords[hull.back()]);
      myVec2 testNorm = getNorm(vecSubtract(leftMostVec, lastVec));
      // now dot testNorm with *it - leftMost. If dot is positive, leftMost becomes it. If dot is zero, take one further from thisHull.back().
      myVec2 itVec(xCoords[*it], yCoords[*it]);
      double dotProd = dotProduct(testNorm, vecSubtract(itVec, lastVec));
      if (-1e-8 < dotProd && dotProd < 1e-8) {
        // thisHull.back(), it and leftMost are collinear.
        // Just sum the differences in x and differences in y for each and compare to get further one, don't need distance formula
        myVec2 itDist       = vecSubtract(itVec, lastVec);
        myVec2 leftMostDist = vecSubtract(leftMostVec, lastVec);
        if (fabs(itDist.x) + fabs(itDist.y) > fabs(leftMostDist.x) + fabs(leftMostDist.y)) {
          leftMost = it;
        }
      } else if (dotProd > 0) {
        leftMost = it;
      }
    }
    // if leftMost is min, then the loop is complete.
    if (*leftMost == *minNode)
      break;
    hull.push_back(*leftMost);
    includeMin = true;  // have found second point (the one after min) so now include min in the searches
    // debug_it ++;
    // if(debug_it > 100) exit(0); //break;
  }
  return hull;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::vector<int> VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::makeUnique(std::vector<int>& vertices) const {
  using namespace std;
  vector<int> uniqueNodes = vertices;
  sort(uniqueNodes.begin(), uniqueNodes.end());
  vector<int>::iterator newUniqueFineEnd = unique(uniqueNodes.begin(), uniqueNodes.end());
  uniqueNodes.erase(newUniqueFineEnd, uniqueNodes.end());
  // uniqueNodes is now a sorted list of the nodes whose info actually goes in file
  // Now replace values in vertices with locations of the old values in uniqueFine
  for (int i = 0; i < int(vertices.size()); i++) {
    int lo     = 0;
    int hi     = uniqueNodes.size() - 1;
    int mid    = 0;
    int search = vertices[i];
    while (lo <= hi) {
      mid = lo + (hi - lo) / 2;
      if (uniqueNodes[mid] == search)
        break;
      else if (uniqueNodes[mid] > search)
        hi = mid - 1;
      else
        lo = mid + 1;
    }
    if (uniqueNodes[mid] != search)
      throw runtime_error("Issue in makeUnique_() - a point wasn't found in list.");
    vertices[i] = mid;
  }
  return uniqueNodes;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceAll(std::string result, const std::string& replaceWhat, const std::string& replaceWithWhat) const {
  while (1) {
    const int pos = result.find(replaceWhat);
    if (pos == -1)
      break;
    result.replace(pos, replaceWhat.size(), replaceWithWhat);
  }
  return result;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getFileName(int numProcs, int myRank, int level, const Teuchos::ParameterList& pL) const {
  std::string filenameToWrite = getBaseFileName(numProcs, level, pL);
  filenameToWrite             = this->replaceAll(filenameToWrite, "%PROCID", toString(myRank));
  return filenameToWrite;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getBaseFileName(int numProcs, int level, const Teuchos::ParameterList& pL) const {
  std::string filenameToWrite = pL.get<std::string>("visualization: output filename");
  int timeStep                = pL.get<int>("visualization: output file: time step");
  int iter                    = pL.get<int>("visualization: output file: iter");

  if (filenameToWrite.rfind(".vtu") == std::string::npos)
    filenameToWrite.append(".vtu");
  if (numProcs > 1 && filenameToWrite.rfind("%PROCID") == std::string::npos)  // filename can't be identical between processsors in parallel problem
    filenameToWrite.insert(filenameToWrite.rfind(".vtu"), "-proc%PROCID");

  filenameToWrite = this->replaceAll(filenameToWrite, "%LEVELID", toString(level));
  filenameToWrite = this->replaceAll(filenameToWrite, "%TIMESTEP", toString(timeStep));
  filenameToWrite = this->replaceAll(filenameToWrite, "%ITER", toString(iter));
  return filenameToWrite;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getPVTUFileName(int numProcs, int /* myRank */, int level, const Teuchos::ParameterList& pL) const {
  std::string filenameToWrite = getBaseFileName(numProcs, level, pL);
  std::string masterStem      = filenameToWrite.substr(0, filenameToWrite.rfind(".vtu"));
  masterStem                  = this->replaceAll(masterStem, "%PROCID", "");
  std::string pvtuFilename    = masterStem + "-master.pvtu";
  return pvtuFilename;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeFileVTKOpening(std::ofstream& fout, std::vector<int>& uniqueFine, std::vector<int>& geomSizesFine) const {
  std::string styleName = "PointCloud";  // TODO adapt this

  std::string indent = "      ";
  fout << "<!--" << styleName << " Aggregates Visualization-->" << std::endl;
  fout << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"LittleEndian\">" << std::endl;
  fout << "  <UnstructuredGrid>" << std::endl;
  fout << "    <Piece NumberOfPoints=\"" << uniqueFine.size() << "\" NumberOfCells=\"" << geomSizesFine.size() << "\">" << std::endl;
  fout << "      <PointData Scalars=\"Node Aggregate Processor\">" << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeFileVTKNodes(std::ofstream& fout, std::vector<int>& uniqueFine, Teuchos::RCP<const Map>& nodeMap) const {
  std::string indent = "      ";
  fout << "        <DataArray type=\"Int32\" Name=\"Node\" format=\"ascii\">" << std::endl;
  indent             = "          ";
  bool localIsGlobal = GlobalOrdinal(nodeMap->getGlobalNumElements()) == GlobalOrdinal(nodeMap->getLocalNumElements());
  for (size_t i = 0; i < uniqueFine.size(); i++) {
    if (localIsGlobal) {
      fout << uniqueFine[i] << " ";  // if all nodes are on this processor, do not map from local to global
    } else
      fout << nodeMap->getGlobalElement(uniqueFine[i]) << " ";
    if (i % 10 == 9)
      fout << std::endl
           << indent;
  }
  fout << std::endl;
  fout << "        </DataArray>" << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeFileVTKData(std::ofstream& fout, std::vector<int>& uniqueFine, LocalOrdinal myAggOffset, ArrayRCP<LocalOrdinal>& vertex2AggId, int myRank) const {
  std::string indent = "          ";
  fout << "        <DataArray type=\"Int32\" Name=\"Aggregate\" format=\"ascii\">" << std::endl;
  fout << indent;
  for (int i = 0; i < int(uniqueFine.size()); i++) {
    fout << myAggOffset + vertex2AggId[uniqueFine[i]] << " ";
    if (i % 10 == 9)
      fout << std::endl
           << indent;
  }
  fout << std::endl;
  fout << "        </DataArray>" << std::endl;
  fout << "        <DataArray type=\"Int32\" Name=\"Processor\" format=\"ascii\">" << std::endl;
  fout << indent;
  for (int i = 0; i < int(uniqueFine.size()); i++) {
    fout << myRank << " ";
    if (i % 20 == 19)
      fout << std::endl
           << indent;
  }
  fout << std::endl;
  fout << "        </DataArray>" << std::endl;
  fout << "      </PointData>" << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeFileVTKCoordinates(std::ofstream& fout, std::vector<int>& uniqueFine, Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& fx, Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& fy, Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>& fz, int dim) const {
  std::string indent = "      ";
  fout << "      <Points>" << std::endl;
  fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
  fout << indent;
  for (int i = 0; i < int(uniqueFine.size()); i++) {
    fout << fx[uniqueFine[i]] << " " << fy[uniqueFine[i]] << " ";
    if (dim == 2)
      fout << "0 ";
    else
      fout << fz[uniqueFine[i]] << " ";
    if (i % 3 == 2)
      fout << std::endl
           << indent;
  }
  fout << std::endl;
  fout << "        </DataArray>" << std::endl;
  fout << "      </Points>" << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeFileVTKCells(std::ofstream& fout, std::vector<int>& /* uniqueFine */, std::vector<LocalOrdinal>& vertices, std::vector<LocalOrdinal>& geomSize) const {
  std::string indent = "      ";
  fout << "      <Cells>" << std::endl;
  fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
  fout << indent;
  for (int i = 0; i < int(vertices.size()); i++) {
    fout << vertices[i] << " ";
    if (i % 10 == 9)
      fout << std::endl
           << indent;
  }
  fout << std::endl;
  fout << "        </DataArray>" << std::endl;
  fout << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
  fout << indent;
  int accum = 0;
  for (int i = 0; i < int(geomSize.size()); i++) {
    accum += geomSize[i];
    fout << accum << " ";
    if (i % 10 == 9)
      fout << std::endl
           << indent;
  }
  fout << std::endl;
  fout << "        </DataArray>" << std::endl;
  fout << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">" << std::endl;
  fout << indent;
  for (int i = 0; i < int(geomSize.size()); i++) {
    switch (geomSize[i]) {
      case 1:
        fout << "1 ";  // Point
        break;
      case 2:
        fout << "3 ";  // Line
        break;
      case 3:
        fout << "5 ";  // Triangle
        break;
      default:
        fout << "7 ";  // Polygon
    }
    if (i % 30 == 29)
      fout << std::endl
           << indent;
  }
  fout << std::endl;
  fout << "        </DataArray>" << std::endl;
  fout << "      </Cells>" << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeFileVTKClosing(std::ofstream& fout) const {
  fout << "    </Piece>" << std::endl;
  fout << "  </UnstructuredGrid>" << std::endl;
  fout << "</VTKFile>" << std::endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writePVTU(std::ofstream& pvtu, std::string baseFname, int numProcs, bool bFineEdges, bool /* bCoarseEdges */) const {
  // If using vtk, filenameToWrite now contains final, correct ***.vtu filename (for the current proc)
  // So the root proc will need to use its own filenameToWrite to make a list of the filenames of all other procs to put in
  // pvtu file.
  pvtu << "<VTKFile type=\"PUnstructuredGrid\" byte_order=\"LittleEndian\">" << std::endl;
  pvtu << "  <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
  pvtu << "    <PPointData Scalars=\"Node Aggregate Processor\">" << std::endl;
  pvtu << "      <PDataArray type=\"Int32\" Name=\"Node\"/>" << std::endl;
  pvtu << "      <PDataArray type=\"Int32\" Name=\"Aggregate\"/>" << std::endl;
  pvtu << "      <PDataArray type=\"Int32\" Name=\"Processor\"/>" << std::endl;
  pvtu << "    </PPointData>" << std::endl;
  pvtu << "    <PPoints>" << std::endl;
  pvtu << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>" << std::endl;
  pvtu << "    </PPoints>" << std::endl;
  for (int i = 0; i < numProcs; i++) {
    // specify the piece for each proc (the replaceAll expression matches what the filenames will be on other procs)
    pvtu << "    <Piece Source=\"" << replaceAll(baseFname, "%PROCID", toString(i)) << "\"/>" << std::endl;
  }
  // reference files with graph pieces, if applicable
  if (bFineEdges) {
    for (int i = 0; i < numProcs; i++) {
      std::string fn = replaceAll(baseFname, "%PROCID", toString(i));
      pvtu << "    <Piece Source=\"" << fn.insert(fn.rfind(".vtu"), "-finegraph") << "\"/>" << std::endl;
    }
  }
  /*if(doCoarseGraphEdges_)
  {
    for(int i = 0; i < numProcs; i++)
    {
      std::string fn = replaceAll(baseFname, "%PROCID", toString(i));
      pvtu << "    <Piece Source=\"" << fn.insert(fn.rfind(".vtu"), "-coarsegraph") << "\"/>" << std::endl;
    }
  }*/
  pvtu << "  </PUnstructuredGrid>" << std::endl;
  pvtu << "</VTKFile>" << std::endl;
  pvtu.close();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void VisualizationHelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>::buildColormap() const {
  try {
    std::ofstream color("random-colormap.xml");
    color << "<ColorMap name=\"MueLu-Random\" space=\"RGB\">" << std::endl;
    // Give -1, -2, -3 distinctive colors (so that the style functions can have constrasted geometry types)
    // Do red, orange, yellow to constrast with cool color spectrum for other types
    color << "  <Point x=\"" << -1 << "\" o=\"1\" r=\"1\" g=\"0\" b=\"0\"/>" << std::endl;
    color << "  <Point x=\"" << -2 << "\" o=\"1\" r=\"1\" g=\"0.6\" b=\"0\"/>" << std::endl;
    color << "  <Point x=\"" << -3 << "\" o=\"1\" r=\"1\" g=\"1\" b=\"0\"/>" << std::endl;
    srand(time(NULL));
    for (int i = 0; i < 5000; i += 4) {
      color << "  <Point x=\"" << i << "\" o=\"1\" r=\"" << (rand() % 50) / 256.0 << "\" g=\"" << (rand() % 256) / 256.0 << "\" b=\"" << (rand() % 256) / 256.0 << "\"/>" << std::endl;
    }
    color << "</ColorMap>" << std::endl;
    color.close();
  } catch (std::exception& e) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError,
                               "VisualizationHelpers::buildColormap: Error while building colormap file: " << e.what());
  }
}

}  // namespace MueLu

#endif /* MUELU_VISUALIZATIONHELPERS_DEF_HPP_ */
