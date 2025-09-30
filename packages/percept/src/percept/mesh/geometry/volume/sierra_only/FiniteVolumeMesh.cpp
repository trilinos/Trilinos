// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/mesh/geometry/volume/sierra_only/FiniteVolumeMesh.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>  // for stk::mesh::field_fill
#include <stddef.h>                     // for size_t
#include <algorithm>                    // for max, min, lower_bound, sort
#include <cmath>                        // for abs, sqrt
#include <fstream>                      // for operator<<, basic_ostream, etc
#include <iomanip>                      // for operator<<, setw
#include <map>                          // for _Rb_tree_iterator, etc
#include <set>                          // for set, set<>::iterator
#include <stk_mesh/base/FieldParallel.hpp>  // for parallel_sum
#include <stk_mesh/base/GetEntities.hpp>  // for get_selected_entities
#include <stk_topology/topology.hpp>    // for topology, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for all_reduce_max, etc
#include <string>                       // for basic_string, operator<<
#include <utility>                      // for pair
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, Bucket::iterator
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator&
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityVector, etc
#include "stk_topology/topology.hpp"    // for topology::num_vertices
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequire
#include "stk_util/util/PairIter.hpp"   // for PairIter

using std::ostream;

namespace percept {

  // taken and modified from Sierra/Aero DualMesh

  FiniteVolumeMesh::FiniteVolumeMesh(const int d, stk::mesh::BulkData & r, stk::mesh::FieldBase *coordinatesField ,
                   stk::mesh::FieldBase *controlVolumeField , stk::mesh::FieldBase *scVolumeField )
    :
    meshStats(d),
    nDim_(d), bulkData_(r), coordinatesField_(coordinatesField), controlVolumeField_(controlVolumeField), scVolumeField_(scVolumeField)
  {
    if (!coordinatesField)
      coordinatesField_ = const_cast<stk::mesh::FieldBase *>(bulkData_.mesh_meta_data().coordinate_field());

    active_part_selector_ = bulkData_.mesh_meta_data().universal_part();

    // const StringNames & sName = StringNames::self();

    // const stk::mesh::MetaData & meta_data = bulkData_.mesh_meta_data();
    // //STK Fields
    // coordinatesField_ = &bulkData_.current_stk_coordinates();
    // controlVolumeField_ = meta_data.get_field<GeneralField_type>(stk::topology::NODE_RANK, sName.control_volume);
    // scVolumeField_ = meta_data.get_field<GeneralField_type>(stk::topology::ELEMENT_RANK, sName.subcontrol_volume);

  }

  FiniteVolumeMesh::~FiniteVolumeMesh() {
  }

  int FiniteVolumeMesh::check_subcontrol_volume(const double volume, const double * centroid, const stk::mesh::Entity elem,
                                                const stk::mesh::Entity * nodes)
  {
    if (volume > 0.0) return 0;

    BadCell inverted;

    const stk::mesh::BulkData & bulk_data = bulkData_;

    for (unsigned k = 0; k < bulk_data.num_nodes(elem); ++k) {
      inverted.node_ids[k] = bulk_data.identifier(nodes[k]);
    }

    inverted.global_id =  bulk_data.identifier(elem);

    inverted.centroid[0] = centroid[0];
    inverted.centroid[1] = centroid[1];
    inverted.centroid[2] = centroid[2];

    badCellList.push_back(inverted);

    return 1;
  }
  void FiniteVolumeMesh::getMinAndMaxVolume() {
    //sanity check
    const stk::mesh::BulkData& bulk_data = bulkData_;

    stk::mesh::BucketVector const& nodeBuckets = bulk_data.get_buckets( stk::topology::NODE_RANK,
                                                                        active_part_selector());

    double minVolume = std::numeric_limits<double>::max();
    double maxVolume = 0;

    for (stk::mesh::BucketVector::const_iterator ib = nodeBuckets.begin(); ib != nodeBuckets.end(); ++ib) {

      const stk::mesh::Bucket & b = **ib;
      const size_t length = b.size();

      for (size_t n = 0; n < length; n++) {
        double * volume = (controlVolumeField_ ? (double *)stk::mesh::field_data(*controlVolumeField_, b[n]) : &nodalVolume_[b[n]]);

        minVolume = std::min(minVolume, *volume);
        maxVolume = std::max(maxVolume, *volume);

      }
    }

    double gMinVolume = 0;
    stk::all_reduce_min(bulk_data.parallel(), &minVolume, &gMinVolume, 1);
    double gMaxVolume = 0;
    stk::all_reduce_max(bulk_data.parallel(), &maxVolume, &gMaxVolume, 1);

    meshStats.minCellVolume = gMinVolume;
    meshStats.maxCellVolume = gMaxVolume;

  }

  void FiniteVolumeMesh::zero_area_and_volume() {
    if (controlVolumeField_)
      stk::mesh::field_fill(0.0, *controlVolumeField_);
    else if (nodalVolume_.size())
      {
        for (auto & it : nodalVolume_)
          {
            //stk::mesh::Entity node = it.first;
            double& vol = it.second;
            vol = 0.0;
          }
      }
    if (scVolumeField_)
      {
        //stk::mesh::field_fill(0.0, *controlVolumeField_);
      }
    else if (scVolumeMap_.size())
      {
        for (auto & it : scVolumeMap_)
          {
            std::vector<double>& vol = it.second;
            vol.assign(vol.size(), 0.0);
          }
      }
  }

  void FiniteVolumeMesh::MeshStatistics::print(std::ostream & outstream) {

    int ndim = minCoordinates.size();

    std::stringstream buffer;

    buffer.setf(std::ios_base::right,std::ios_base::adjustfield);
    const unsigned cellw = 32;

    buffer << std::setw(cellw) << "total number of elements:" << std::setw(cellw) << globalNumElements << "\n";
    buffer << std::setw(cellw) << "total number of edges:" << std::setw(cellw) << globalNumEdges << "\n";
    buffer << std::setw(cellw) << "total number of nodes:" << std::setw(cellw) << globalNumNodes << "\n";
    buffer << std::setw(cellw) << "minimum cell volume:" << std::setw(cellw) << minCellVolume << "\n";
    buffer << std::setw(cellw) << "maximum cell volume:" << std::setw(cellw) << maxCellVolume << "\n";
    buffer << std::setw(cellw) << "total mesh volume:" << std::setw(cellw) << totalMeshVolume << "\n";
    buffer << std::setw(cellw) << "minimum edge length:" << std::setw(cellw) << minEdgeLength << "\n";
    buffer << std::setw(cellw) << "minimum edge area:" << std::setw(cellw) << minEdgeArea << "\n";
    buffer << std::setw(cellw) << "maximum area closure error:" << std::setw(cellw) << maxAreaClosureError << "\n";
    std::stringstream bbLength;

    buffer << std::setw(cellw) << "mesh bounding box lengths:";
    bbLength.precision(2);
    bbLength << "(";
    for(int i = 0; i< ndim; ++i) {
      bbLength  << std::scientific << maxCoordinates[i] - minCoordinates[i];
      if (i < ndim-1 ) {
        bbLength << ", ";
      } else {
        bbLength << ")";
      }
    }
    buffer  << std::setw(cellw) << bbLength.str() << "\n";

    double loadBal = double(maxLocalNumNodes)/double(minLocalNumNodes);
    buffer  << std::setw(cellw) << "node load balance (max/min):" << std::setw(cellw) <<
      loadBal << "\n";
    loadBal = double(maxLocalNumEdges)/double(minLocalNumEdges);
    buffer  << std::setw(cellw) << "edge load balance (max/min):" << std::setw(cellw) <<
      loadBal << "\n";
    loadBal = double(maxLocalNumElements)/double(minLocalNumElements);
    buffer  << std::setw(cellw) << "element load balance (max/min):" << std::setw(cellw) <<
      loadBal << "\n";


    outstream << std::endl << std::endl;
    outstream << "==================================================================== \n";
    outstream << "                          Mesh Statistics \n\n";
    outstream << "==================================================================== \n";
    outstream << buffer.str();
    outstream << "==================================================================== \n";
    outstream << std::endl << std::endl;
  }


  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------
  FiniteVolumeMesh2D::FiniteVolumeMesh2D(stk::mesh::BulkData& r, stk::mesh::FieldBase *coordinatesField ,
                   stk::mesh::FieldBase *controlVolumeField , stk::mesh::FieldBase *scVolumeField )
    : FiniteVolumeMesh(2, r, coordinatesField, controlVolumeField, scVolumeField)
  {
  }


  void FiniteVolumeMesh2D::computeAreaAndVolume() {
    /* %TRACE% */


    zero_area_and_volume();

    const stk::mesh::BulkData& bulk_data = bulkData_;

    stk::mesh::BucketVector const& buckets = bulk_data.get_buckets(stk::topology::ELEM_RANK,
                                                                   active_part_selector());

    for (stk::mesh::BucketVector::const_iterator ib = buckets.begin(); ib != buckets.end(); ++ib) {

      const stk::mesh::Bucket& b = **ib;
      stk::mesh::Bucket::iterator elements = b.begin();
      const size_t nelem = b.size();

      meshStats.localNumElements += nelem;

      const stk::topology top = b.topology();
      const int numNodes = top.num_vertices();

      for (size_t e = 0; e < nelem; e++) {
        stk::mesh::Entity elem = elements[e];

        double edge_midpoints[4][2];
        double centroid[3] = {0, 0, 0};
        if (!scVolumeField_)
          {
            scVolumeMap_[elem].resize(numNodes);
            scVolumeMap_[elem].assign(numNodes, 0.0);
          }
        double* sc_volume = (scVolumeField_ ? (double *)stk::mesh::field_data(*scVolumeField_, elem) : &scVolumeMap_[elem][0]);
        const stk::mesh::Entity * nodes = bulk_data.begin_nodes(elem);

        getCentroidAndEdgeMidpoints(numNodes, nodes, centroid, edge_midpoints);

        // compute dual areas by breaking subcontrol quad into two subtriangles
        int last_n = numNodes - 1;

        for (int n = 0; n < numNodes; ++n) {
          stk::mesh::Entity  node = nodes[n];
          // if (!controlVolumeField_)
          //   {
          //     nodalVolume_[node]
          //   }
          double* volume = (controlVolumeField_? (double *)stk::mesh::field_data(*controlVolumeField_, node) : &nodalVolume_[node]);
          const double* x_n = (double *)stk::mesh::field_data(*coordinatesField_, node);
          // first triangle
          const double vol1 = 0.5
            * std::abs(
                       x_n[0] * (edge_midpoints[n][1] - centroid[1]) + edge_midpoints[n][0] * (centroid[1] - x_n[1])
                       + centroid[0] * (x_n[1] - edge_midpoints[n][1]));
          // second triangle
          const double vol2 = 0.5
            * std::abs(
                       x_n[0] * (centroid[1] - edge_midpoints[last_n][1]) + centroid[0] * (edge_midpoints[last_n][1] - x_n[1])
                       + edge_midpoints[last_n][0] * (x_n[1] - centroid[1]));
          last_n = n;
          const double scVol = vol1 + vol2;
          sc_volume[n] = scVol;
          *volume += scVol;

          const int badSub = check_subcontrol_volume(scVol, centroid, elem, nodes);
          if(badSub != 0) break;

        }
      }
    }

    //doParallelAndCheckForErrors();

  }


  void
  FiniteVolumeMesh2D::getCentroidAndEdgeMidpoints(int numNodes, const stk::mesh::Entity * nodes,
                                                  double centroid[2], double edge_midpoints[4][2]) {

    centroid[0] = 0;
    centroid[1] = 0;

    for (int n = 0; n < numNodes; ++n) {
      int next_n = n + 1 < numNodes ?
        n + 1 : 0;
      stk::mesh::Entity  node = nodes[n];
      stk::mesh::Entity  next_node = nodes[next_n];
      const double* x_n = (double *)stk::mesh::field_data(*coordinatesField_, node);
      const double* next_x = (double *)stk::mesh::field_data(*coordinatesField_, next_node);
      edge_midpoints[n][0] = 0.5 * (x_n[0] + next_x[0]);
      edge_midpoints[n][1] = 0.5 * (x_n[1] + next_x[1]);
      centroid[0] += x_n[0];
      centroid[1] += x_n[1];
    }

    centroid[0] = centroid[0] / double(numNodes);
    centroid[1] = centroid[1] / double(numNodes);

  }

  //------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------
  FiniteVolumeMesh3D::FiniteVolumeMesh3D(stk::mesh::BulkData& r, stk::mesh::FieldBase *coordinatesField ,
                   stk::mesh::FieldBase *controlVolumeField , stk::mesh::FieldBase *scVolumeField )
    : FiniteVolumeMesh(3, r, coordinatesField, controlVolumeField, scVolumeField)
  {
  }

  // return true if neg vol
  bool FiniteVolumeMesh3D::elementVolume(stk::mesh::Entity elem, double *sc_volume)
  {
    const stk::topology & top = bulkData_.bucket(elem).topology();
    const int numNodes = top.num_vertices();

    // all elements in a bucket are homogeneous
    //Hex 8 is the worst case.  Volumes are calculated by decomposing all
    // elements into hexes
    double subcontrol_elem_coords[27][3];

    double centroid[3] = {0, 0, 0};

    elemCoordsAndCentroid(elem, numNodes, centroid, subcontrol_elem_coords);

    if (!scVolumeField_)
      {
        scVolumeMap_[elem].resize(numNodes);
      }

    const stk::mesh::Entity * nodes = bulkData_.begin_nodes(elem);

    std::set<stk::mesh::Entity> nodeList;
    for (int n = 0; n < numNodes; ++n)
      {
        nodeList.insert(nodes[n]);
      }
    STK_ThrowRequire((int)nodeList.size() == numNodes);

    if (top.value() == stk::topology::HEX_8) {

      hexSubcontrolCoords(centroid, subcontrol_elem_coords);

      for (int n = 0; n < numNodes; ++n) {
        double hex_coords[8][3];
        for (int j = 0; j < 8; ++j) {
          for (int i = 0; i < 3; ++i) {
            hex_coords[j][i] = subcontrol_elem_coords[HexSubcontrolNodeTable[n][j]][i];
          }
        }
        const double scVol = hexVolume(hex_coords);
        sc_volume[n] = scVol;
        const int badSub = check_subcontrol_volume(scVol, centroid, elem, nodes);
        if(badSub != 0) return true;
      }

      // accumulateSubcontrolAreaVectors(elem, subcontrol_elem_coords, HexEdgeFacetTable,
      //                                 HexEdgeNodeOrder, 12);

    } else if (top.value() == stk::topology::TET_4) {

      tetSubcontrolCoords(centroid, subcontrol_elem_coords);

      for (int n = 0; n < numNodes; ++n) {
        //double * volume = (double *)stk::mesh::field_data(*controlVolumeField_, node);
        double hex_coords[8][3];
        for (int j = 0; j < 8; ++j) {
          for (int i = 0; i < 3; ++i) {
            hex_coords[j][i] = subcontrol_elem_coords[TetSubcontrolNodeTable[n][j]][i];
          }
        }

        const double scVol = hexVolume(hex_coords);
        sc_volume[n] = scVol;

        const int badSub = check_subcontrol_volume(scVol, centroid, elem, nodes);
        if(badSub != 0) return true;

      }

      // accumulateSubcontrolAreaVectors(elem, subcontrol_elem_coords, TetEdgeFacetTable,
      //                                 TetEdgeNodeOrder, 6);

    } else if (top.value() == stk::topology::WEDGE_6) {

      STK_ThrowRequire(6 == numNodes);
      wedgeSubcontrolCoords(centroid, subcontrol_elem_coords);

      for (int n = 0; n < numNodes; ++n) {
        double hex_coords[8][3];
        for (int j = 0; j < 8; ++j) {
          for (int i = 0; i < 3; ++i) {
            hex_coords[j][i] = subcontrol_elem_coords[WedgeSubcontrolNodeTable[n][j]][i];
          }
        }

        const double scVol = hexVolume(hex_coords);
        sc_volume[n] = scVol;
        const int badSub = check_subcontrol_volume(scVol, centroid, elem, nodes);
        if(badSub != 0) return true;

      }
      // accumulateSubcontrolAreaVectors(elem, subcontrol_elem_coords, WedgeEdgeFacetTable,
      //                                 WedgeEdgeNodeOrder, 9);

    } else if (top.value() == stk::topology::PYRAMID_5) {

      STK_ThrowRequire(5 == numNodes);
      pyramidSubcontrolCoords(centroid, subcontrol_elem_coords);

      for (int n = 0; n < numNodes - 1; ++n) //skip apex since it is an octohedron not a hex
        {
          //double * volume = (double *)stk::mesh::field_data(*controlVolumeField_, node);
          //double* volume = (controlVolumeField_? (double *)stk::mesh::field_data(*controlVolumeField_, node) : &nodalVolume_[node]);

          double hex_coords[8][3];
          for (int j = 0; j < 8; ++j) {
            for (int i = 0; i < 3; ++i) {
              hex_coords[j][i] = subcontrol_elem_coords[PyramidSubcontrolNodeTable[n][j]][i];
            }
          }

          const double scVol = hexVolume(hex_coords);
          sc_volume[n] = scVol;
          const int badSub = check_subcontrol_volume(scVol, centroid, elem, nodes);
          if(badSub != 0) return true;

        }

      // now do octohedron for apex
      double pyramidTipCoords[10][3];

      // stk::mesh::Entity node = nodes[numNodes-1];
      //double * volume = (double *)stk::mesh::field_data(*controlVolumeField_, node);
      //double* volume = (controlVolumeField_? (double *)stk::mesh::field_data(*controlVolumeField_, node) : &nodalVolume_[node]);
      for (int j = 0; j < 10; ++j) {
        for (int i = 0; i < 3; ++i) {
          pyramidTipCoords[j][i] = subcontrol_elem_coords[PyramidSubcontrolNodeTable[4][j]][i];
        }
      }

      const double scVol = pyramidTipVolume(pyramidTipCoords);
      sc_volume[numNodes-1] = scVol;

      const int badSub = check_subcontrol_volume(scVol, centroid, elem, nodes);
      if(badSub != 0) return true;

      // accumulateSubcontrolAreaVectors(elem, subcontrol_elem_coords, PyramidEdgeFacetTable,
      //                                 PyramidEdgeNodeOrder, 8);

    } else {

      throw std::runtime_error( "Unsupported element of type " + top.name() + " encountered.");
    }
    return false;
  }

  void FiniteVolumeMesh3D::computeAreaAndVolume() {

    zero_area_and_volume();

    const stk::mesh::BulkData& bulk_data = bulkData_;

    stk::mesh::BucketVector const& buckets = bulk_data.get_buckets(stk::topology::ELEM_RANK,
                                                                   active_part_selector());

    for (stk::mesh::BucketVector::const_iterator ib = buckets.begin(); ib != buckets.end(); ++ib) {

      const stk::mesh::Bucket & b = **ib;
      stk::mesh::Bucket::iterator elements = b.begin();

      const size_t nelem = b.size();
      meshStats.localNumElements += nelem;

      const stk::topology & top = b.topology();
      const int numNodes = top.num_vertices();

      for (size_t e = 0; e < nelem; e++) {

        stk::mesh::Entity elem = elements[e];
        double* sc_volume = (scVolumeField_ ? (double *)stk::mesh::field_data(*scVolumeField_, elem) : &scVolumeMap_[elem][0]);

        elementVolume(elem, sc_volume);

        const stk::mesh::Entity * nodes = bulk_data.begin_nodes(elem);
        for (int n = 0; n < numNodes; ++n)
          {
            stk::mesh::Entity node = nodes[n];
            double* volume = (controlVolumeField_? (double *)stk::mesh::field_data(*controlVolumeField_, node) : &nodalVolume_[node]);
            *volume += sc_volume[n];
          }
      }
    }
  }

  double FiniteVolumeMesh3D::hexVolume(const double coord[][3]) {
    /* This function works by defining an irregular hex to be a collection
     * of 24 triangular facets.  by inserting a node at the center of each quad face and subsequently
     *  breaking each quad face of the hex into four triangles.  (see "Efficient computation of volume of
     *  Hexahedral Cells", Jeffrey Grandy, LLNL, UCRL-ID-128886, October 30, 1997.
     */
    const int numCoords = 14;
    const int numTriangles = 24;
    double x[14][3];
    for (int n = 0; n < 8; ++n) {
      for (int i = 0; i < 3; ++i) {
        x[n][i] = coord[n][i];
      }
    }

    //compute the six face midpoints
    for (int i = 0; i < 3; ++i) {
      x[8][i] = 0.25 * (coord[0][i] + coord[1][i] + coord[2][i] + coord[3][i]);
      x[9][i] = 0.25 * (coord[4][i] + coord[5][i] + coord[6][i] + coord[7][i]);
      x[10][i] = 0.25 * (coord[0][i] + coord[1][i] + coord[5][i] + coord[4][i]);
      x[11][i] = 0.25 * (coord[3][i] + coord[2][i] + coord[6][i] + coord[7][i]);
      x[12][i] = 0.25 * (coord[1][i] + coord[2][i] + coord[6][i] + coord[5][i]);
      x[13][i] = 0.25 * (coord[0][i] + coord[3][i] + coord[7][i] + coord[4][i]);
    }
    int triangularFaceTable[24][3] = { {0, 8, 1}, {8, 2, 1}, {3, 2, 8}, {3, 8, 0}, {6, 9, 5}, {7, 9, 6}, {4, 9, 7}, {4, 5,
                                                                                                                     9}, {10, 0, 1}, {5, 10, 1}, {4, 10, 5}, {4, 0, 10}, {7, 6, 11}, {6, 2, 11}, {2, 3, 11}, {3, 7, 11}, {6, 12, 2}, {
        5, 12, 6}, {5, 1, 12}, {1, 2, 12}, {0, 4, 13}, {4, 7, 13}, {7, 3, 13}, {3, 0, 13}};
    return polyhedralVolume(numCoords, x, numTriangles, triangularFaceTable);
  }

  void FiniteVolumeMesh3D::hexSubcontrolCoords(const double* centroid, double element_coords[][3]) {
    //this node ordering must match the definitions in HexSubcontrolNodeTable
    const int numNodesPerElem = 8;
    const int numNodesPerFace = 4;
    const double numFaceNodesInv = 1. / double(numNodesPerFace);
    int counter = numNodesPerElem;
    double face_centroid[3] = {0, 0, 0};
    int next_node;
    // face 0
    for (int node = 0; node < 4; ++node) {
      next_node = node == 3 ?
        0 : node + 1;
      for (int i = 0; i < 3; ++i) {
        element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
        face_centroid[i] += element_coords[node][i] * numFaceNodesInv;
      }
      counter++;
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    // face 1
    for (int node = 4; node < 8; ++node) {
      next_node = (node == 7) ?
        4 : node + 1;
      for (int i = 0; i < 3; ++i) {
        element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
        face_centroid[i] += element_coords[node][i] * numFaceNodesInv;
      }
      counter++;
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //face 2
    int node = 1;
    next_node = 5;
    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
    }
    counter++;
    node = 0;
    next_node = 4;
    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
    }
    counter++;
    for (int n = 0; n < numNodesPerFace; ++n) {
      for (int i = 0; i < 3; ++i) {
        face_centroid[i] += element_coords[HexFaceTable[2][n]][i] * numFaceNodesInv;
      }
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //face 3
    node = 3;
    next_node = 7;
    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
    }
    counter++;
    node = 2;
    next_node = 6;
    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
    }
    counter++;
    for (int n = 0; n < numNodesPerFace; ++n) {
      for (int i = 0; i < 3; ++i) {
        face_centroid[i] += element_coords[HexFaceTable[3][n]][i] * numFaceNodesInv;
      }
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //face 4
    for (int n = 0; n < numNodesPerFace; ++n) {
      for (int i = 0; i < 3; ++i) {
        face_centroid[i] += element_coords[HexFaceTable[4][n]][i] * numFaceNodesInv;
      }
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //face 5
    for (int n = 0; n < numNodesPerFace; ++n) {
      for (int i = 0; i < 3; ++i) {
        face_centroid[i] += element_coords[HexFaceTable[5][n]][i] * numFaceNodesInv;
      }
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //element centroid
    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = centroid[i];
    }
    counter++;
    if (counter != 27)
      throw std::runtime_error( "counter should be 27");
  }

  void FiniteVolumeMesh3D::tetSubcontrolCoords(const double* centroid, double element_coords[][3]) {
    //this node ordering must match the definitions in TetSubcontrolNodeTable
    const int numNodesPerElem = 4;
    const int numNodesPerFace = 3;
    const double numFaceNodesInv = 1. / double(numNodesPerFace);
    int counter = numNodesPerElem;
    double face_centroid[3] = {0, 0, 0};
    // first face
    for (int node = 0; node < numNodesPerFace; ++node) {
      int next_node = node == 2 ?
        0 : node + 1;
      for (int i = 0; i < 3; ++i) {
        element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
        face_centroid[i] += element_coords[node][i] * numFaceNodesInv;
      }
      counter++;
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //second face
    for (int node = 2; node < 4; ++node) {
      int next_node = node == 3 ?
        1 : node + 1;
      for (int i = 0; i < 3; ++i) {
        element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
      }
      counter++;
    }

    for (int n = 0; n < numNodesPerFace; ++n) {
      for (int i = 0; i < 3; ++i) {
        face_centroid[i] += element_coords[TetFaceTable[1][n]][i] * numFaceNodesInv;
      }
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //third face
    int node = 0;
    int next_node = 3;
    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
    }
    counter++;
    for (int n = 0; n < numNodesPerFace; ++n) {
      for (int i = 0; i < 3; ++i) {
        face_centroid[i] += element_coords[TetFaceTable[2][n]][i] * numFaceNodesInv;
      }
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //fourth face
    for (int n = 0; n < numNodesPerFace; ++n) {
      for (int i = 0; i < 3; ++i) {
        face_centroid[i] += element_coords[TetFaceTable[3][n]][i] * numFaceNodesInv;
      }
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //element centroid
    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = centroid[i];
    }
    counter++;
    if (counter != 15)
      throw std::runtime_error( "counter should be 15");
  }

  void FiniteVolumeMesh3D::wedgeSubcontrolCoords(const double* centroid, double element_coords[][3]) {
    //this node ordering must match the definitions in WedgeubcontrolNodeTable
    const int numNodesPerElem = 6;
    const int numNodesPerTriFace = 3;
    const int numNodesPerQuadFace = 4;
    const double numTriFaceNodesInv = 1. / double(numNodesPerTriFace);
    const double numQuadFaceNodesInv = 1. / double(numNodesPerQuadFace);
    int counter = numNodesPerElem;
    double face_centroid[3] = {0, 0, 0};
    // first face
    for (int node = 0; node < numNodesPerTriFace; ++node) {
      int next_node = node == 2 ?
        0 : node + 1;
      for (int i = 0; i < 3; ++i) {
        element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
        face_centroid[i] += element_coords[node][i] * numTriFaceNodesInv;
      }
      counter++;
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //second face
    for (int node = 3; node < 6; ++node) {
      int next_node = node == 5 ?
        3 : node + 1;
      for (int i = 0; i < 3; ++i) {
        element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
      }
      counter++;
    }

    for (int n = 0; n < numNodesPerTriFace; ++n) {
      for (int i = 0; i < 3; ++i) {
        face_centroid[i] += element_coords[WedgeFaceTable[1][n]][i] * numTriFaceNodesInv;
      }
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //third face
    int node = 1;
    int next_node = 4;
    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
    }
    counter++;
    node = 3;
    next_node = 0;
    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
    }
    counter++;
    for (int n = 0; n < numNodesPerQuadFace; ++n) {
      for (int i = 0; i < 3; ++i) {
        face_centroid[i] += element_coords[WedgeFaceTable[2][n]][i] * numQuadFaceNodesInv;
      }
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //fourth face
    node = 5;
    next_node = 2;
    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
    }
    counter++;
    for (int n = 0; n < numNodesPerQuadFace; ++n) {
      for (int i = 0; i < 3; ++i) {
        face_centroid[i] += element_coords[WedgeFaceTable[3][n]][i] * numQuadFaceNodesInv;
      }
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    // fifth face
    for (int n = 0; n < numNodesPerQuadFace; ++n) {
      for (int i = 0; i < 3; ++i) {
        face_centroid[i] += element_coords[WedgeFaceTable[4][n]][i] * numQuadFaceNodesInv;
      }
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //element centroid
    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = centroid[i];
    }
    counter++;
    if (counter != 21)
      throw std::runtime_error( "counter should be 21");
  }

  void FiniteVolumeMesh3D::pyramidSubcontrolCoords(const double* centroid, double element_coords[][3]) {
    //this node ordering must match the definitions in TetSubcontrolNodeTable
    const int numNodesPerElem = 5;
    const int numNodesPerTriFace = 3;
    const int numNodesPerQuadFace = 4;
    const double numTriFaceNodesInv = 1. / double(numNodesPerTriFace);
    const double numQuadFaceNodesInv = 1. / double(numNodesPerQuadFace);
    int counter = numNodesPerElem;
    double face_centroid[3] = {0, 0, 0};
    // first face
    for (int node = 0; node < numNodesPerQuadFace; ++node) {
      int next_node = node == 3 ?
        0 : node + 1;
      for (int i = 0; i < 3; ++i) {
        element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
        face_centroid[i] += element_coords[node][i] * numQuadFaceNodesInv;
      }
      counter++;
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //second face
    int node = 1;
    int next_node = 4;
    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
    }
    counter++;
    node = 4;
    next_node = 0;
    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
    }
    counter++;
    for (int n = 0; n < numNodesPerTriFace; ++n) {
      for (int i = 0; i < 3; ++i) {
        face_centroid[i] += element_coords[PyramidFaceTable[1][n]][i] * numTriFaceNodesInv;
      }
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //third face
    node = 2;
    next_node = 4;
    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
    }
    counter++;
    for (int n = 0; n < numNodesPerTriFace; ++n) {
      for (int i = 0; i < 3; ++i) {
        face_centroid[i] += element_coords[PyramidFaceTable[2][n]][i] * numTriFaceNodesInv;
      }
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //fourth face
    node = 4;
    next_node = 3;
    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = 0.5 * (element_coords[node][i] + element_coords[next_node][i]);
    }
    counter++;
    for (int n = 0; n < numNodesPerTriFace; ++n) {
      for (int i = 0; i < 3; ++i) {
        face_centroid[i] += element_coords[PyramidFaceTable[3][n]][i] * numTriFaceNodesInv;
      }
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //fifth face
    for (int n = 0; n < numNodesPerTriFace; ++n) {
      for (int i = 0; i < 3; ++i) {
        face_centroid[i] += element_coords[PyramidFaceTable[4][n]][i] * numTriFaceNodesInv;
      }
    }

    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = face_centroid[i];
      face_centroid[i] = 0;
    }
    counter++;
    //element centroid
    for (int i = 0; i < 3; ++i) {
      element_coords[counter][i] = centroid[i];
    }
    counter++;
    if (counter != 19)
      throw std::runtime_error( "counter should be 19");
  }

  double FiniteVolumeMesh3D::pyramidTipVolume(const double coord[][3]) {
    /*  This function works by recognizing that  the subcontrol volume of the pyramid tip
     *  has four planar quad faces and four nonplanar qud faces.  It inserts a node at the center of
     *  each nonplanar quad face and subsequently
     *  breaking each nonplanar quad face of the hex into four triangles.  Note that the quads formed by
     *  traversing the vertex->edge midpoint->face midpoint->edge midpoint->vertex of an original triangular face are
     *  planar and do not need to have a node injected.
     */
    const int numCoords = 14;
    const int numTriangles = 24;
    double x[14][3];
    for (int n = 0; n < 10; ++n) {
      for (int i = 0; i < 3; ++i) {
        x[n][i] = coord[n][i];
      }
    }

    //compute the four face midpoints
    for (int i = 0; i < 3; ++i) {
      x[10][i] = 0.25 * (coord[1][i] + coord[3][i] + coord[2][i] + coord[9][i]);
      x[11][i] = 0.25 * (coord[1][i] + coord[5][i] + coord[4][i] + coord[3][i]);
      x[12][i] = 0.25 * (coord[1][i] + coord[7][i] + coord[6][i] + coord[5][i]);
      x[13][i] = 0.25 * (coord[1][i] + coord[9][i] + coord[8][i] + coord[7][i]);
    }
    static int triangularFaceTable[24][3] = { {1, 3, 10}, {2, 10, 3}, {2, 9, 10}, {10, 9, 1}, {4, 3, 11}, {3, 1, 11},
                                              {11, 1, 5}, {4, 11, 5}, {1, 12, 5}, {1, 7, 12}, {12, 7, 6}, {5, 12, 6}, {9, 8, 13}, {13, 8, 7}, {13, 7, 1}, {9, 13, 1},
                                              {4, 5, 0}, {5, 6, 0}, {6, 7, 0}, {7, 8, 0}, {0, 8, 9}, {0, 9, 2}, {0, 2, 3}, {0, 3, 4}};
    return polyhedralVolume(numCoords, x, numTriangles, triangularFaceTable);
  }

  double FiniteVolumeMesh3D::polyhedralVolume(
                                              const int /*numCoords*/,
                                              const double x[][3],
                                              const int numTriangles,
                                              const int triangularFaceTable[][3])
  {
    /* This function works by considering a polyhedron to be bonded by a collection of triangular facets.
     * .  Then it uses the gauss divergence theorem to compute the volume as 1/3 the
     * integral of the divergence of the coordinate vector :
     * 3V = integral{ (x0,x1,x2) dot A}.
     */
    double volume = 0;
    //iterate the triangles
    for (int t = 0; t < numTriangles; ++t) {
      const int p = triangularFaceTable[t][0];
      const int q = triangularFaceTable[t][1];
      const int r = triangularFaceTable[t][2];
      double xFace[3] = {0, 0, 0};
      for (int i = 0; i < 3; ++i) {
        xFace[i] = x[p][i] + x[q][i] + x[r][i];
      }
      volume += xFace[0] * ((x[q][1] - x[p][1]) * (x[r][2] - x[p][2]) - (x[r][1] - x[p][1]) * (x[q][2] - x[p][2]))
        - xFace[1] * ((x[q][0] - x[p][0]) * (x[r][2] - x[p][2]) - (x[r][0] - x[p][0]) * (x[q][2] - x[p][2]))
        + xFace[2] * ((x[q][0] - x[p][0]) * (x[r][1] - x[p][1]) - (x[r][0] - x[p][0]) * (x[q][1] - x[p][1]));
    }

    volume = volume / 18.0; // = 1/3 * 1/2 * 1/3 factor comes from green gauss, area vector and face average

    return volume;
  }
  void FiniteVolumeMesh3D::elemCoordsAndCentroid(
                                                 stk::mesh::Entity elem,
                                                 const int numNodes,
                                                 double centroid[3],
                                                 double element_coords[][3])
  {
    const stk::mesh::BulkData& bulk_data = bulkData_;
    centroid[0] = 0;
    centroid[1] = 0;
    centroid[2] = 0;
    const stk::mesh::Entity * nodes = bulk_data.begin_nodes(elem);
    for (int n = 0; n < numNodes; ++n) {
      stk::mesh::Entity node = nodes[n];
      const double* x = (double *)stk::mesh::field_data(*coordinatesField_, node);
      centroid[0] += x[0] / double(numNodes);
      centroid[1] += x[1] / double(numNodes);
      centroid[2] += x[2] / double(numNodes);
      element_coords[n][0] = x[0];
      element_coords[n][1] = x[1];
      element_coords[n][2] = x[2];
    }
  }

  void FiniteVolumeMesh3D::quadFaceAreaVector(const double coord[4][3], double* areaVector) {
    double x[5][3];
    x[4][0] = 0;
    x[4][1] = 0;
    x[4][2] = 0;
    for (int n = 0; n < 4; n++) {
      for (int i = 0; i < 3; ++i) {
        x[n][i] = coord[n][i];
        x[4][i] += coord[n][i];
      }
    }

    x[4][0] *= .25;
    x[4][1] *= .25;
    x[4][2] *= .25;
    static int triangularFaceTable[4][3] = { {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}};
    double twiceArea[3] = {0, 0, 0};
    for (int tri = 0; tri < 4; ++tri) {
      double r1[3];
      double r2[3];
      int* indx = &triangularFaceTable[tri][0];
      for (int i = 0; i < 3; ++i) {
        r1[i] = x[indx[1]][i] - x[indx[0]][i];
        r2[i] = x[indx[2]][i] - x[indx[0]][i];
      }
      twiceArea[0] += r1[1] * r2[2] - r2[1] * r1[2];
      twiceArea[1] += r1[2] * r2[0] - r1[0] * r2[2];
      twiceArea[2] += r1[0] * r2[1] - r2[0] * r1[1];
    }

    areaVector[0] = .5 * twiceArea[0];
    areaVector[1] = .5 * twiceArea[1];
    areaVector[2] = .5 * twiceArea[2];
  }



}

