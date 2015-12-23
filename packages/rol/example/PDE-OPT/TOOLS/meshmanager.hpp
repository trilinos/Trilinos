// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  meshmanager.hpp
    \brief Defines the MeshManger classes.
*/

#ifndef MESHMANAGER_HPP
#define MESHMANAGER_HPP

#include "Teuchos_ParameterList.hpp"
#include "Intrepid_FieldContainer.hpp"


/** \class  MeshManager
    \brief  This is the pure virtual parent class for mesh construction
            and management; it enables the generation of a few select
            meshes and the setup of data structures for computing
            cell/subcell adjacencies, for example, the cell-to-node
            adjacency map and the cell-to-edge adjacency map.
*/
template <class Real>
class MeshManager {

public:
  /** \brief  Destructor
   */
  virtual ~MeshManager() {}

  /** \brief Computes node coordinates.
             Format: number_of_nodes x 2 (Real)
                     (node_index)  x, y coordinates
  */
  virtual void getNodes(Intrepid::FieldContainer<Real> &nodes) = 0;

  /** \brief Computes cell to node adjacencies.
             Format: number_of_cells x number_of_nodes_per_cell (int)
                     (cell_index)  node_index1  node_index2  ...
  */
  virtual void getCellToNodeMap(Intrepid::FieldContainer<int> &ctn) = 0;

  /** \brief Computes cell to edge adjacencies.
             Format: number_of_cells x number_of_edges_per_cell (int)
                     (cell_index)  edge_index1  edge_index2  ...
  */
  virtual void getCellToEdgeMap(Intrepid::FieldContainer<int> &ctn) = 0;

  /** \brief Returns number of cells.
  */
  virtual int getNumCells() const = 0;

  /** \brief Returns number of nodes.
  */
  virtual int getNumNodes() const = 0;

  /** \brief Returns number of edges.
  */
  virtual int getNumEdges() const = 0;

}; // MeshManager



/** \class  MeshManager_BackwardFacingStepChannel
    \brief  Mesh construction and mesh management for the
            backward-facing step channel geometry, on
            quadrilateral grids.
*/
template <class Real>
class MeshManager_BackwardFacingStepChannel : public MeshManager<Real> {

/* Backward-facing step geometry.
 
  ***************************************************
  *        *           *                            *
  *   3    *     4     *             5              *
  *        *           *                            *
  ********** * * * * * * * * * * * * * * * * * * * **
           *     1     *             2              *
           *           *                            *
           ******************************************

*/

private:
  Real channelH_;   // channel height (height of regions 1+4)
  Real channelW_;   // channel width (width of regions 3+4+5)
  Real stepH_;      // step height (height of region 1)
  Real stepW_;      // step width (width of region 3)
  Real observeW_;   // width of observation region (width of region 1)
  
  int ref_;          // mesh refinement level
  bool store_;       // storage flag
  bool initialized_; // initialization flag

  int nx1_;
  int nx2_;
  int nx3_;
  int nx4_;
  int nx5_;
  int ny1_;
  int ny2_;
  int ny3_;
  int ny4_;
  int ny5_;

  int numCells_;
  int numNodes_;
  int numEdges_;
  int numMeshEntities_;

  Intrepid::FieldContainer<Real> meshNodes_;
  Intrepid::FieldContainer<int>  meshCellToNodeMap_;
  Intrepid::FieldContainer<int>  meshCellToEdgeMap_;

public:

  MeshManager_BackwardFacingStepChannel(Teuchos::ParameterList &parlist) {
    // Geometry data.
    channelH_ = parlist.sublist("Geometry").get(   "Channel height", 1.0);
    channelW_ = parlist.sublist("Geometry").get(    "Channel width", 8.0);
    stepH_    = parlist.sublist("Geometry").get(      "Step height", 0.5);
    stepW_    = parlist.sublist("Geometry").get(       "Step width", 1.0);
    observeW_ = parlist.sublist("Geometry").get("Observation width", 3.0);
    // Mesh data.
    ref_   = parlist.sublist("Geometry").get(      "Refinement level", 1);
    store_ = parlist.sublist("Geometry").get(       "Store mesh data", true);
    nx1_   = parlist.sublist("Geometry").get( "Observation region NX", 4*ref_);
    ny1_   = parlist.sublist("Geometry").get( "Observation region NY", 5*ref_);
    nx2_   = parlist.sublist("Geometry").get("Laminar flow region NX", 2*ref_);
    ny2_   = ny1_;
    nx3_   = parlist.sublist("Geometry").get(      "Inflow region NX", 1*ref_);
    ny3_   = parlist.sublist("Geometry").get(      "Inflow region NY", 2*ref_);
    nx4_   = nx1_;
    ny4_   = ny3_;
    nx5_   = nx2_;
    ny5_   = ny3_;
    initialized_ = false;
    numCells_ = (nx1_ + nx2_)*ny1_  +  (nx3_ + nx1_ + nx2_)*ny3_;
    numNodes_ = (nx1_ + nx2_ + 1)*ny1_  +  (nx3_ + nx1_ + nx2_ + 1)*(ny3_ + 1);
    numEdges_ = (2*(nx1_ + nx2_)+1)*ny1_ + (2*(nx3_ + nx1_ + nx2_)+1)*ny3_ + (nx3_ + nx1_ + nx2_);
    numMeshEntities_ = numCells_ + numNodes_ + numEdges_;
    // Compute and store mesh data structures.
    if (store_) {
      getMeshNodes(meshNodes_); 
      getMeshCellToNodeMap(meshCellToNodeMap_); 
      getMeshCellToEdgeMap(meshCellToEdgeMap_);
    }
    initialized_ = true;
  }


  void getNodes(Intrepid::FieldContainer<Real> &nodes) {

    if ((store_) && (initialized_)) {
      nodes = meshNodes_;
      return;
    }

    Real dy1 = stepH_ / ny1_;
    Real dy3 = (channelH_ - stepH_) / ny3_;
    Real dx1 = observeW_ / nx1_;
    Real dx2 = (channelW_ - stepW_ - observeW_) / nx2_;
    Real dx3 = stepW_ / nx3_;
    int nodeCt = 0;
    nodes.resize(numNodes_, 2);

    // bottom region
    for (int j=0; j<ny1_; ++j) {
      for (int i=0; i<=nx1_; ++i) {
        nodes(nodeCt, 0) = stepW_ + i*dx1;
        nodes(nodeCt, 1) = j*dy1; 
        ++nodeCt;
      }
      for (int i=0; i<nx2_; ++i) {
        nodes(nodeCt, 0) = stepW_ + observeW_ + (i+1)*dx2;
        nodes(nodeCt, 1) = j*dy1;
        ++nodeCt;
      }
    }

    // top region
    for (int j=0; j<=ny3_; ++j) {
      for (int i=0; i<=nx3_; ++i) {
        nodes(nodeCt, 0) = i*dx3;
        nodes(nodeCt, 1) = stepH_ + j*dy3; 
        ++nodeCt;
      }
      for (int i=0; i<nx1_; ++i) {
        nodes(nodeCt, 0) = stepW_ + (i+1)*dx1;
        nodes(nodeCt, 1) = stepH_ + j*dy3;
        ++nodeCt;
      }
      for (int i=0; i<nx2_; ++i) {
        nodes(nodeCt, 0) = stepW_ + observeW_ + (i+1)*dx2;
        nodes(nodeCt, 1) = stepH_ + j*dy3;
        ++nodeCt;
      }
    }

  }  // getMeshNodes


  void getMeshCellToNodeMap(Intrepid::FieldContainer<int> &ctn) {

    if ((store_) && (initialized_)) {
      ctn = meshCellToNodeMap_;
      return;
    }

    int cellCt = 0;
    ctn.resize(numCells_, 4);

    // bottom region
    for (int j=0; j<ny1_-1; ++j) {
      for (int i=0; i<nx1_+nx2_; ++i) {
        ctn(cellCt, 0) = j*(nx1_+nx2_+1) + i;
        ctn(cellCt, 1) = j*(nx1_+nx2_+1) + (i+1);
        ctn(cellCt, 2) = (j+1)*(nx1_+nx2_+1) + (i+1);
        ctn(cellCt, 3) = (j+1)*(nx1_+nx2_+1) + i;
        ++cellCt;
      }
    }

    // transition region
    for (int i=0; i<nx1_+nx2_; ++i) {
      ctn(cellCt, 0) = (ny1_-1)*(nx1_+nx2_+1) + i;
      ctn(cellCt, 1) = (ny1_-1)*(nx1_+nx2_+1) + (i+1);
      ctn(cellCt, 2) = ny1_*(nx1_+nx2_+1) + nx3_ + (i+1);
      ctn(cellCt, 3) = ny1_*(nx1_+nx2_+1) + nx3_ + i;
      ++cellCt;
    }

    // top region
    for (int j=0; j<ny3_; ++j) {
      for (int i=0; i<nx3_+nx1_+nx2_; ++i) {
        ctn(cellCt, 0) = ny1_*(nx1_+nx2_+1) + j*(nx3_+nx1_+nx2_+1) + i;
        ctn(cellCt, 1) = ny1_*(nx1_+nx2_+1) + j*(nx3_+nx1_+nx2_+1) + (i+1);
        ctn(cellCt, 2) = ny1_*(nx1_+nx2_+1) + (j+1)*(nx3_+nx1_+nx2_+1) + (i+1);
        ctn(cellCt, 3) = ny1_*(nx1_+nx2_+1) + (j+1)*(nx3_+nx1_+nx2_+1) + i;
        ++cellCt;
      }
    }

  } // getMeshCellToNodeMap


  void getMeshCellToEdgeMap(Intrepid::FieldContainer<int> &cte) {

    if ((store_) && (initialized_)) {
      cte = meshCellToEdgeMap_;
      return;
    }

    int cellCt = 0;
    cte.resize(numCells_, 4);

    // bottom region
    for (int j=0; j<ny1_-1; ++j) {
      for (int i=0; i<nx1_+nx2_; ++i) {
        cte(cellCt, 0) = j*(2*(nx1_+nx2_)+1) + i;
        cte(cellCt, 1) = j*(2*(nx1_+nx2_)+1) + (nx1_+nx2_) + (i+1);
        cte(cellCt, 2) = (j+1)*(2*(nx1_+nx2_)+1) + i;
        cte(cellCt, 3) = j*(2*(nx1_+nx2_)+1) + (nx1_+nx2_) + i;
        ++cellCt;
      }
    }

    // transition region
    for (int i=0; i<nx1_+nx2_; ++i) {
      cte(cellCt, 0) = (ny1_-1)*(2*(nx1_+nx2_)+1) + i;
      cte(cellCt, 1) = (ny1_-1)*(2*(nx1_+nx2_)+1) + (nx1_+nx2_) + (i+1); 
      cte(cellCt, 2) = ny1_*(2*(nx1_+nx2_)+1) + nx3_ + i;
      cte(cellCt, 3) = (ny1_-1)*(2*(nx1_+nx2_)+1) + (nx1_+nx2_) + i;
      ++cellCt;
    }

    // top region
    for (int j=0; j<ny3_; ++j) {
      for (int i=0; i<nx3_+nx1_+nx2_; ++i) {
        cte(cellCt, 0) = ny1_*(2*(nx1_+nx2_)+1) + j*(2*(nx3_+nx1_+nx2_)+1) + i;
        cte(cellCt, 1) = ny1_*(2*(nx1_+nx2_)+1) + j*(2*(nx3_+nx1_+nx2_)+1) + (nx3_+nx1_+nx2_) + (i+1);
        cte(cellCt, 2) = ny1_*(2*(nx1_+nx2_)+1) + (j+1)*(2*(nx3_+nx1_+nx2_)+1) + i;
        cte(cellCt, 3) = ny1_*(2*(nx1_+nx2_)+1) + j*(2*(nx3_+nx1_+nx2_)+1) + (nx3_+nx1_+nx2_) + i;
        ++cellCt;
      }
    }

  } // getMeshCellToEdgeMap


  int getNumCells() const {
    return numCells_;
  } // getNumCells


  int getNumNodes() const {
    return numNodes_;
  } // getNumNodes


  int getNumEdges() const {
    return numEdges_;
  } // getNumEdges


  void setRefinementLevel(const int &refLevel) {
    ref_ = refLevel;
    if (store_) {
      initialized_ = false;
      getMeshNodes(meshNodes_);
      getMeshCellToNodeMap(meshCellToNodeMap_);
      getMeshCellToEdgeMap(meshCellToEdgeMap_);
      initialized_ = true;
    }
  } // setRefinementLevel

}; // MeshManager_BackwardFacingStepChannel



/** \class  MeshManager_Rectangle
    \brief  Mesh construction and mesh management for the
            backward-facing step channel geometry, on
            quadrilateral grids.
*/
template <class Real>
class MeshManager_Rectangle : public MeshManager<Real> {

/* Rectangle geometry.
 
      ***********************
      *                     *   :
      *                     *   |
      *                     * height
      *                     *   |
      *                     *   :
      *                     *
      ***********************
  (X0,Y0)   :--width--:

*/

private:
  Real width_;   // rectangle height
  Real height_;  // rectangle width
  Real X0_;      // x coordinate of bottom left corner
  Real Y0_;      // y coordinate of bottom left corner

  int nx_;
  int ny_;

  bool store_;
  bool initialized_;

  int numCells_;
  int numNodes_;
  int numEdges_;
  int numMeshEntities_;

  Intrepid::FieldContainer<Real> meshNodes_;
  Intrepid::FieldContainer<int>  meshCellToNodeMap_;
  Intrepid::FieldContainer<int>  meshCellToEdgeMap_;

public:

  MeshManager_Rectangle(Teuchos::ParameterList &parlist) {
    // Geometry data.
    height_ = parlist.sublist("Geometry").get("Height", 1.0);
    width_  = parlist.sublist("Geometry").get( "Width", 3.0);
    X0_     = parlist.sublist("Geometry").get(    "X0", 0.0);
    Y0_     = parlist.sublist("Geometry").get(    "Y0", 0.0);
    // Mesh data.
    nx_ = parlist.sublist("Geometry").get("NX", 1);
    ny_ = parlist.sublist("Geometry").get("NY", 3);
    store_ = parlist.sublist("Geometry").get("Store mesh data", true);
    initialized_ = false;
    numCells_ = nx_ * ny_;
    numNodes_ = (nx_+1) * (ny_+1);
    numEdges_ = (nx_+1)*ny_ + (ny_+1)*nx_;
    numMeshEntities_ = numCells_ + numNodes_ + numEdges_;
    // Compute and store mesh data structures.
    if (store_) {
      getNodes(meshNodes_); 
      getCellToNodeMap(meshCellToNodeMap_); 
      getCellToEdgeMap(meshCellToEdgeMap_);
    }
    initialized_ = true;
  }


  void getNodes(Intrepid::FieldContainer<Real> &nodes) {

    if ((store_) && (initialized_)) {
      nodes = meshNodes_;
      return;
    }

    Real dx = width_ / nx_;
    Real dy = height_ / ny_;
    int nodeCt = 0;
    nodes.resize(numNodes_, 2);

    for (int j=0; j<=ny_; ++j) {
      Real ycoord = Y0_ + j*dy;
      for (int i=0; i<=nx_; ++i) {
        nodes(nodeCt, 0) = X0_ + i*dx;
        nodes(nodeCt, 1) = ycoord; 
        ++nodeCt;
      }
    }

  } // getMeshNodes


  void getCellToNodeMap(Intrepid::FieldContainer<int> &ctn) {

    if ((store_) && (initialized_)) {
      ctn = meshCellToNodeMap_;
      return;
    }

    int cellCt = 0;
    ctn.resize(numCells_, 4);

    for (int j=0; j<ny_; ++j) {
      for (int i=0; i<nx_; ++i) {
        ctn(cellCt, 0) = j*(nx_+1) + i;
        ctn(cellCt, 1) = j*(nx_+1) + (i+1);
        ctn(cellCt, 2) = (j+1)*(nx_+1) + (i+1);
        ctn(cellCt, 3) = (j+1)*(nx_+1) + i;
        ++cellCt;
      }
    }

  } // getMeshCellToNodeMap


  void getCellToEdgeMap(Intrepid::FieldContainer<int> &cte) {

    if ((store_) && (initialized_)) {
      cte = meshCellToEdgeMap_;
      return;
    }

    int cellCt = 0;
    cte.resize(numCells_, 4);

    for (int j=0; j<ny_; ++j) {
      for (int i=0; i<nx_; ++i) {
        cte(cellCt, 0) = j*(2*nx_+1) + i;
        cte(cellCt, 1) = j*(2*nx_+1) + nx_ + (i+1);
        cte(cellCt, 2) = (j+1)*(2*nx_+1) + i;
        cte(cellCt, 3) = j*(2*nx_+1) + nx_ + i;
        ++cellCt;
      }
    }

  } // getMeshCellToEdgeMap


  int getNumCells() const {
    return numCells_;
  } // getNumCells


  int getNumNodes() const {
    return numNodes_;
  } // getNumNodes


  int getNumEdges() const {
    return numEdges_;
  } // getNumEdges

}; // MeshManager_Rectangle

#endif
