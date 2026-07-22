// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef percept_FiniteVolumeMesh_h
#define percept_FiniteVolumeMesh_h

#include <iostream>
#include <vector>                       // for vector
#include <stddef.h>                     // for size_t
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/Field.hpp"      // for Field

#include "stk_mesh/base/Selector.hpp"   // for BulkData

namespace percept {


  // taken and modified from Sierra/Aero DualMesh
  /*
   * CellFaceAndEdgeNodes.h
   *
   *  Created on: Apr 10, 2014
   *      Author: swbova
   */

  const int QuadEdgeNodeOrder[4][2] = // [edge][edge_node]
    { {0,1}, {1,2}, {2,3}, {0, 3} };

  const int QuadSubcontrolNodeTable[4][4] = {
    {0, 4, 8, 7},
    {4, 1, 5, 8},
    {8, 5, 2, 6},
    {7, 8, 6, 3}
  };


  const int TriEdgeNodeOrder[3][2] = // [edge][edge_node]
    { {0,1}, {1,2}, {0, 2} };

  const int TriSubcontrolNodeTable[3][4] = {
    {0, 3, 6, 5},
    {3, 1, 4, 6},
    {5, 6, 4, 2}
  };

  // four nodes for each of six faces
  const int HexFaceTable[6][4] = {
    {0, 1, 2, 3},
    {4, 5, 6, 7},
    {0, 1, 5, 4},
    {3, 7, 6, 2},
    {1, 2, 6, 5},
    {0, 3, 7, 4} };

  //8 nodes for each of 8 subcontrol volumes
  const int HexSubcontrolNodeTable[8][8] = {
    {0, 8, 12, 11, 19, 20, 26, 25},
    {8, 1, 9, 12, 20, 18, 24, 26},
    {12, 9, 2, 10, 26, 24, 22, 23},
    {11, 12, 10, 3, 25, 26, 23, 21},
    {19, 20, 26, 25, 4, 13, 17, 16},
    {20, 18, 24, 26, 13, 5, 14, 17},
    {26, 24, 22, 23, 17, 14, 6, 15},
    {25, 26, 23, 21, 16, 17, 15, 7} };

  const int HexEdgeFacetTable[12][4] = {
    {20, 8, 12, 26},
    {24, 9, 12, 26},
    {10, 12, 26, 23},
    {11, 25, 26, 12},
    {13, 20, 26, 17},
    {17, 14, 24, 26},
    {17, 15, 23, 26},
    {16, 17, 26, 25},
    {19, 20, 26, 25},
    {20, 18, 24, 26},
    {22, 23, 26, 24},
    {21, 25, 26, 23}
  };
  const int HexEdgeNodeOrder[12][2] = // [edge][edge_node]
    { {0,1}, {1,2}, {2,3}, {3, 0},
      {4,5}, {5,6}, {6,7}, {7, 4},
      {0,4}, {1,5}, {2,6}, {3,7} };


  const int TetFaceTable[4][3] = { {0, 1, 2},
                                   {1, 2, 3},
                                   {0, 2, 3},
                                   {0, 1, 3} };

  // 8 nodes for each of four subcontrol hex volumes
  const int TetSubcontrolNodeTable[4][8] = {
    {0, 4, 7, 6, 11, 13, 14, 12},
    {1, 5, 7, 4, 9, 10, 14, 13},
    {2, 6, 7, 5, 8, 12, 14, 10},
    {3, 9, 13, 11, 8, 10, 14, 12} };

  const int TetEdgeFacetTable[6][4] = {  // [edge][subcontrol facet node]
    {4, 7, 14, 13},
    {7, 14, 10, 5},
    {6, 12, 14, 7},
    {11, 13, 14, 12},
    {13, 9, 10, 14},
    {10, 8, 12, 14}  };
  const int TetEdgeNodeOrder[6][2] = // [edge][edge_node]
    { {0,1}, {1,2}, {2,0}, {0,3}, {1,3}, {2,3} };


  const int WedgeFaceTable[5][4] = {
    {0, 2, 1, -1},
    {3, 5, 4, -1},
    {0, 1, 4, 3},
    {1, 4, 5, 2},
    {5, 3, 0, 2} };


  // 8 nodes for each of six subcontrol hex volumes
  const int WedgeSubcontrolNodeTable[6][8] = {
    {0, 15, 16, 6, 8, 19, 20, 9}, //node 0
    {9, 6, 1, 7, 20, 16, 14, 18}, //node 1
    {8, 9, 7, 2, 19, 20, 18, 17}, // node 2
    {19, 15, 16, 20, 12, 3, 10, 13}, //node 3
    {20, 16, 14, 18, 13, 10, 4, 11},  //node 4
    {19, 20, 18, 17, 12, 13, 11, 5} };  //node 5

  const int WedgeEdgeFacetTable[9][4] = {  //[edge][subcontrol facet node]
    {6, 9, 20, 16},
    {7, 9, 20, 18},
    {9, 8, 19, 20},
    {10, 16, 20, 13},
    {13, 11, 18, 20},
    {12, 13, 20, 19},
    {15, 16, 20, 19},
    {16, 14, 18, 20},
    {19, 20, 18, 17}  };

  const int WedgeEdgeNodeOrder[9][2] = // [edge][edge_node]
    {  {0,1},  {1,2}, {2,0},
       {3,4},  {4,5}, {5,3},
       {0,3},  {1,4}, {2,5}  };


  const int PyramidFaceTable[5][4] = {
    {0, 1, 2, 3},
    {0, 1, 4, -1},
    {1, 2, 4, -1},
    {3, 4, 2, -1},
    {0, 4, 3, -1} };


  // 8 nodes for each of four subcontrol hex volumes plus the apex octohedron
  const int PyramidSubcontrolNodeTable[5][10] = {
    {0, 5, 9, 8, 11, 12, 18, 17, -1, -1},  //node 0
    {5, 1, 6, 9, 12, 10, 14, 18, -1, -1}, //node 1
    {6, 2, 7, 9, 14, 13, 16, 18, -1, -1},// node 2
    {8, 9, 7, 3, 17, 18, 16, 15, -1, -1}, //node 3
    { 4, 18, 15 , 17, 11, 12, 10, 14, 13, 16 } };  //node 4


  const int PyramidEdgeFacetTable[8][4] = { //[edge][subcontrol face node]
    {5, 9, 18, 12}, //edge 0,1
    {6, 9, 18, 14}, //edge 1,2
    {7, 9, 18, 16}, //edge 2,3
    {8, 17, 18, 9}, //edge 3.0
    {11, 12, 18, 17}, //edge 0, 4
    {10, 14, 18, 12}, //edge 1,4
    {13, 16, 18, 14}, //edge 2,4
    {15, 17, 18, 16} };//edge 3,4

  const int PyramidEdgeNodeOrder[8][2] = // [edge][edge_node]
    {  {0,1}, {1,2}, {2,3}, {3,0},
       {0,4}, {1,4}, {2,4}, {3,4},
    };




class FiniteVolumeMesh {
public:

  std::ofstream myfile;

  FiniteVolumeMesh(const int d, stk::mesh::BulkData & r, stk::mesh::FieldBase *coordinatesField_ = 0,
                   stk::mesh::FieldBase *controlVolumeField_ = 0, stk::mesh::FieldBase *scVolumeField_ = 0);
  virtual ~FiniteVolumeMesh();

  virtual void computeAreaAndVolume() = 0;

  // defaults to universal_part
  stk::mesh::Selector active_part_selector() { return active_part_selector_ ; }
  void set_active_part_selector(stk::mesh::Selector s ) { active_part_selector_ = s; }

  //returns 0 for a volume that passes
  int check_subcontrol_volume(const double volume, const double * centroid,
                              const stk::mesh::Entity elem, const stk::mesh::Entity * nodes);

  void getMinAndMaxVolume();
  void zero_area_and_volume();

  //double check_area_closure();

  struct MeshStatistics {
    size_t globalNumNodes;
    size_t globalNumEdges;
    size_t globalNumElements;

    int localNumNodes;
    int localNumEdges;
    int localNumElements;

    int minLocalNumNodes;
    int minLocalNumEdges;
    int minLocalNumElements;

    int maxLocalNumNodes;
    int maxLocalNumEdges;
    int maxLocalNumElements;

    double minCellVolume;
    double maxCellVolume;
    double totalMeshVolume;
    double minEdgeLength;
    double minEdgeArea;
    double maxAreaClosureError;
    std::vector<double> minCoordinates;
    std::vector<double> maxCoordinates;

    void print(std::ostream & outStream);

    MeshStatistics(int nDim) :
      globalNumNodes(0),
      globalNumEdges(0),
      globalNumElements(0),
      localNumNodes(0),
      localNumEdges(0),
      localNumElements(0),
      totalMeshVolume(0){
      minCoordinates.resize(nDim);
      maxCoordinates.resize(nDim);
    }
  private:
    MeshStatistics();
  };

  MeshStatistics meshStats;

protected:
  int nDim_;
  stk::mesh::BulkData & bulkData_;

  stk::mesh::FieldBase * coordinatesField_;
  stk::mesh::FieldBase * controlVolumeField_;
  stk::mesh::FieldBase * scVolumeField_;

  stk::mesh::Selector active_part_selector_;

  typedef std::map<stk::mesh::Entity, double> NodeVolumeMap;
  typedef std::map<stk::mesh::Entity, std::vector<double> > ElementSCVolumeMap;
  NodeVolumeMap nodalVolume_;
  ElementSCVolumeMap scVolumeMap_;

  void collectCellsWithNegativeSubElements();

  struct BadCell {
    size_t global_id;
    size_t node_ids[8];
    double centroid[3];

    BadCell() : global_id(size_t(-1)) {
      for(int i = 0; i < 8; ++i)
        node_ids[i] = size_t(-1);
      for(int i = 0; i < 3; ++i)
        centroid[i] = 0;
    }
  };

  std::vector<BadCell> badCellList;

};

class FiniteVolumeMesh2D : public FiniteVolumeMesh{
public:

  FiniteVolumeMesh2D(stk::mesh::BulkData & r, stk::mesh::FieldBase *coordinatesField_ = 0,
                   stk::mesh::FieldBase *controlVolumeField_ = 0, stk::mesh::FieldBase *scVolumeField_ = 0);
  virtual ~FiniteVolumeMesh2D() {}

  virtual void computeAreaAndVolume() override;

  //virtual void computeBoundaryArea();

private:
  void getCentroidAndEdgeMidpoints(int numNodes,
                                   const stk::mesh::Entity * nodes,
                                   double centroid[2],
                                   double edge_midpoints[4][2]);

};
class FiniteVolumeMesh3D : public FiniteVolumeMesh{
public:

  FiniteVolumeMesh3D(stk::mesh::BulkData & r, stk::mesh::FieldBase *coordinatesField_ = 0,
                   stk::mesh::FieldBase *controlVolumeField_ = 0, stk::mesh::FieldBase *scVolumeField_ = 0);
  virtual ~FiniteVolumeMesh3D() {}

  virtual void computeAreaAndVolume() override;

  //virtual EdgeListType & buildEdgeList();

  //virtual void computeBoundaryArea();

  virtual bool elementVolume(stk::mesh::Entity elem, double *sc_volume);

private:
  void hexSubcontrolCoords(const double *centroid,
                           double subcontrol_coords[][3]);
  void tetSubcontrolCoords(const double *centroid,
                           double subcontrol_coords[][3]);
  void wedgeSubcontrolCoords(const double *centroid,
                             double subcontrol_coords[][3]);

  double hexVolume(const double coords[][3]);  //look below at HexFaceTable to see how the nodes are numbered

  void pyramidSubcontrolCoords(const double *centroid,
                               double subcontrol_coords[][3]);

  double pyramidTipVolume(const double coords[][3]);  //look below at HexFaceTable to see how the nodes are numbered

  double polyhedralVolume(const int numCoords,
                          const double x[][3], const int numTriangles, const int triangularFaceTable[][3]);

  void elemCoordsAndCentroid(stk::mesh::Entity elem,
                             const int numNodes,
                             double centroid[3],
                             double elem_coords[][3]);

  void accumulateSubcontrolAreaVectors(stk::mesh::Entity elem,
                                       const double subcontrol_elem_coords[][3],
                                       const int edgeFacetTable[][4],
                                       const int edgeNodeOrder[][2],
                                       const int numEdges);

  void quadFaceAreaVector(const double faceCoords[4][3], double * areaVector);

};



}

#endif
