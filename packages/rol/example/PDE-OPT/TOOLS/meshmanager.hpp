// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  meshmanager.hpp
    \brief Defines the MeshManager classes.
*/

#ifndef MESHMANAGER_HPP
#define MESHMANAGER_HPP

#include "Teuchos_ParameterList.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "ROL_Ptr.hpp"

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

  /** \brief Returns node coordinates.
             Format: number_of_nodes x 2 (Real)
                     (node_index)  x, y coordinates
  */
  virtual ROL::Ptr<Intrepid::FieldContainer<Real> > getNodes() const = 0;

  /** \brief Returns cell to node adjacencies.
             Format: number_of_cells x number_of_nodes_per_cell (int)
                     (cell_index)  node_index1  node_index2  ...
  */
  virtual ROL::Ptr<Intrepid::FieldContainer<int> > getCellToNodeMap() const = 0;

  /** \brief Returns cell to edge adjacencies.
             Format: number_of_cells x number_of_edges_per_cell (int)
                     (cell_index)  edge_index1  edge_index2  ...
  */
  virtual ROL::Ptr<Intrepid::FieldContainer<int> > getCellToEdgeMap() const {
    return ROL::makePtr<Intrepid::FieldContainer<int>>(); // default due to lack of edges in 1D
  }

  /** \brief Returns cell to face adjacencies.
             Format: number_of_cells x number_of_faces_per_cell (int)
                     (cell_index)  face_index1  face_index2  ...
  */
  virtual ROL::Ptr<Intrepid::FieldContainer<int> > getCellToFaceMap() const {
    return ROL::makePtr<Intrepid::FieldContainer<int>>(); // default due to lack of faces in 1D and 2D
  }

  /** \brief Returns cell IDs per processor.
             Format: number_of_procs vectors, each with number_of_cells cell IDs (ints)
  */
  virtual ROL::Ptr<std::vector<std::vector<int>>> getProcCellIds() const {
    return ROL::makePtr<std::vector<std::vector<int>>>();
    // default due to support for in-line partitioning,
    // where this functionality is bypassed
  }

  /** \brief Returns sideset information.
             Format: The std::vector components are indexed by the local side number (0, 1, 2, ...);
                     the FieldConTainer is a 1D array of cell indices.
             Input:  Sideset number.  Its meaning is context-dependent.
  */
  //virtual ROL::Ptr<std::vector<std::vector<Intrepid::FieldContainer<int> > > > getSideSets(
  //            std::ostream & outStream = std::cout,
  //            const bool verbose = false) const = 0;
  virtual ROL::Ptr<std::vector<std::vector<std::vector<int> > > > getSideSets(
              const bool verbose = false,
              std::ostream & outStream = std::cout) const = 0;

  /** \brief Returns number of cells.
  */
  virtual int getNumCells() const = 0;

  /** \brief Returns number of nodes.
  */
  virtual int getNumNodes() const = 0;

  /** \brief Returns number of edges.
  */
  virtual int getNumEdges() const {
    return 0; // default due to lack of edges in 1D
  }

  /** \brief Returns number of faces.
  */
  virtual int getNumFaces() const {
    return 0; // default due to lack of faces in 1D and 2D
  }

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

  ROL::Ptr<Intrepid::FieldContainer<Real> > meshNodes_;
  ROL::Ptr<Intrepid::FieldContainer<int> >  meshCellToNodeMap_;
  ROL::Ptr<Intrepid::FieldContainer<int> >  meshCellToEdgeMap_;

  ROL::Ptr<std::vector<std::vector<std::vector<int> > > >  meshSideSets_;

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
    /*nx1_   = parlist.sublist("Geometry").get( "Observation region NX", 3*2*ref_);
    ny1_   = parlist.sublist("Geometry").get( "Observation region NY", 1*ref_);
    nx2_   = parlist.sublist("Geometry").get("Laminar flow region NX", 4*2*ref_);
    ny2_   = ny1_;
    nx3_   = parlist.sublist("Geometry").get(      "Inflow region NX", 1*2*ref_);
    ny3_   = parlist.sublist("Geometry").get(      "Inflow region NY", 1*ref_);*/
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
    numCells_ = (nx1_ + nx2_)*ny1_  +  (nx3_ + nx1_ + nx2_)*ny3_;
    numNodes_ = (nx1_ + nx2_ + 1)*ny1_  +  (nx3_ + nx1_ + nx2_ + 1)*(ny3_ + 1);
    numEdges_ = (2*(nx1_ + nx2_)+1)*ny1_ + (2*(nx3_ + nx1_ + nx2_)+1)*ny3_ + (nx3_ + nx1_ + nx2_);
    // Compute mesh data structures.
    computeNodes();
    computeCellToNodeMap();
    computeCellToEdgeMap();
    computeSideSets();
  }


  ROL::Ptr<Intrepid::FieldContainer<Real> > getNodes() const {
    return meshNodes_;
  }


  ROL::Ptr<Intrepid::FieldContainer<int> > getCellToNodeMap() const {
    return meshCellToNodeMap_;
  }


  ROL::Ptr<Intrepid::FieldContainer<int> > getCellToEdgeMap() const {
    return meshCellToEdgeMap_;
  }


  ROL::Ptr<std::vector<std::vector<std::vector<int> > > > getSideSets(
              const bool verbose = false,
              std::ostream & outStream = std::cout) const { 
    return meshSideSets_;
  }


  int getNumCells() const {
    return numCells_;
  } // getNumCells


  int getNumNodes() const {
    return numNodes_;
  } // getNumNodes


  int getNumEdges() const {
    return numEdges_;
  } // getNumEdges


private:

  void computeNodes() {

    meshNodes_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(numNodes_, 2);
    Intrepid::FieldContainer<Real> &nodes = *meshNodes_;

    Real dy1 = stepH_ / ny1_;
    Real dy3 = (channelH_ - stepH_) / ny3_;
    Real dx1 = observeW_ / nx1_;
    Real dx2 = (channelW_ - stepW_ - observeW_) / nx2_;
    Real dx3 = stepW_ / nx3_;
    int nodeCt = 0;

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

  }  // computeNodes


  void computeCellToNodeMap() {

    meshCellToNodeMap_ = ROL::makePtr<Intrepid::FieldContainer<int>>(numCells_, 4);
    Intrepid::FieldContainer<int> &ctn = *meshCellToNodeMap_;

    int cellCt = 0;

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

  } // computeCellToNodeMap


  void computeCellToEdgeMap() {

    meshCellToEdgeMap_ = ROL::makePtr<Intrepid::FieldContainer<int>>(numCells_, 4);
    Intrepid::FieldContainer<int> &cte = *meshCellToEdgeMap_;

    int cellCt = 0;

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

  } // computeCellToEdgeMap


  void setRefinementLevel(const int &refLevel) {
    ref_ = refLevel;
    computeNodes();
    computeCellToNodeMap();
    computeCellToEdgeMap();
  } // setRefinementLevel


  virtual void computeSideSets() {

    meshSideSets_ = ROL::makePtr<std::vector<std::vector<std::vector<int> > >>(11);
    int numSides = 4;
    (*meshSideSets_)[0].resize(numSides); // bottom
    (*meshSideSets_)[1].resize(numSides); // right lower
    (*meshSideSets_)[2].resize(numSides); // right upper
    (*meshSideSets_)[3].resize(numSides); // top
    (*meshSideSets_)[4].resize(numSides); // left upper
    (*meshSideSets_)[5].resize(numSides); // middle
    (*meshSideSets_)[6].resize(numSides); // left lower
    (*meshSideSets_)[7].resize(1); // L-corner cell
    (*meshSideSets_)[8].resize(numSides); // top corner cell
    (*meshSideSets_)[9].resize(1); // between side 2 and 3
    (*meshSideSets_)[10].resize(numSides); // top corner cell
    (*meshSideSets_)[0][0].resize(nx1_+nx2_);
    (*meshSideSets_)[1][1].resize(ny2_);
    (*meshSideSets_)[2][1].resize(ny5_);
    (*meshSideSets_)[3][2].resize(nx3_+nx4_+nx5_);
    (*meshSideSets_)[4][3].resize(ny3_);
    (*meshSideSets_)[5][0].resize(nx3_);
    (*meshSideSets_)[6][3].resize(ny1_-1);
    (*meshSideSets_)[7][0].resize(1);
    (*meshSideSets_)[8][3].resize(1);
    (*meshSideSets_)[9][0].resize(1);
    (*meshSideSets_)[10][3].resize(ny1_);

    for (int i=0; i<nx1_+nx2_; ++i) {
      (*meshSideSets_)[0][0][i] = i;
    }
    for (int i=0; i<ny2_; ++i) {
      (*meshSideSets_)[1][1][i] = (i+1)*(nx1_+nx2_) - 1;
    }
    int offset = nx1_*ny1_+nx2_*ny2_;
    for (int i=0; i<ny5_; ++i) {
      (*meshSideSets_)[2][1][i] = offset + (i+1)*(nx3_+nx4_+nx5_) - 1;
    }
    for (int i=0; i<nx3_+nx4_+nx5_; ++i) {
      (*meshSideSets_)[3][2][i] = offset + (ny3_-1)*(nx3_+nx4_+nx5_) + i;
    }
    for (int i=0; i<ny3_; ++i) {
      (*meshSideSets_)[4][3][i] = offset + i*(nx3_+nx4_+nx5_);
    }
    for (int i=0; i<nx3_; ++i) {
      (*meshSideSets_)[5][0][i] = offset + i;
    }
    for (int i=0; i<ny1_-1; ++i) {
      (*meshSideSets_)[6][3][i] = i*(nx1_+nx2_);
    }
    (*meshSideSets_)[7][0][0] = offset + nx3_;
    (*meshSideSets_)[8][3][0] = (ny1_-1)*(nx1_+nx2_);
    (*meshSideSets_)[9][0][0] = offset + ny5_*(nx3_+nx4_+nx5_) - 1;
    for (int i=0; i<ny1_; ++i) {
      (*meshSideSets_)[10][3][i] = i*(nx1_+nx2_);
    }

  } // computeSideSets

}; // MeshManager_BackwardFacingStepChannel



/** \class  MeshManager_Rectangle
    \brief  Mesh construction and mesh management for the
            rectangle geometry, on quadrilateral grids.
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

  int numCells_;
  int numNodes_;
  int numEdges_;

  ROL::Ptr<Intrepid::FieldContainer<Real> > meshNodes_;
  ROL::Ptr<Intrepid::FieldContainer<int> >  meshCellToNodeMap_;
  ROL::Ptr<Intrepid::FieldContainer<int> >  meshCellToEdgeMap_;

  ROL::Ptr<std::vector<std::vector<std::vector<int> > > >  meshSideSets_;

public:

  MeshManager_Rectangle(Teuchos::ParameterList &parlist) {
    // Geometry data.
    width_  = parlist.sublist("Geometry").get( "Width", 3.0);
    height_ = parlist.sublist("Geometry").get("Height", 1.0);
    X0_     = parlist.sublist("Geometry").get(    "X0", 0.0);
    Y0_     = parlist.sublist("Geometry").get(    "Y0", 0.0);
    // Mesh data.
    nx_ = parlist.sublist("Geometry").get("NX", 3);
    ny_ = parlist.sublist("Geometry").get("NY", 1);
    numCells_ = nx_ * ny_;
    numNodes_ = (nx_+1) * (ny_+1);
    numEdges_ = (nx_+1)*ny_ + (ny_+1)*nx_;
    // Compute and store mesh data structures.
    computeNodes();
    computeCellToNodeMap();
    computeCellToEdgeMap();
    computeSideSets();
  }


  ROL::Ptr<Intrepid::FieldContainer<Real> > getNodes() const {
    return meshNodes_;
  }


  ROL::Ptr<Intrepid::FieldContainer<int> > getCellToNodeMap() const {
    return meshCellToNodeMap_;
  }


  ROL::Ptr<Intrepid::FieldContainer<int> > getCellToEdgeMap() const {
    return meshCellToEdgeMap_;
  }


  virtual ROL::Ptr<std::vector<std::vector<std::vector<int> > > > getSideSets(
              const bool verbose = false,
              std::ostream & outStream = std::cout) const { 
    return meshSideSets_;
  }


  int getNumCells() const {
    return numCells_;
  } // getNumCells


  int getNumNodes() const {
    return numNodes_;
  } // getNumNodes


  int getNumEdges() const {
    return numEdges_;
  } // getNumEdges

private:

  void computeNodes() {

    meshNodes_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(numNodes_, 2);
    Intrepid::FieldContainer<Real> &nodes = *meshNodes_;

    Real dx = width_ / nx_;
    Real dy = height_ / ny_;
    int nodeCt = 0;

    for (int j=0; j<=ny_; ++j) {
      Real ycoord = Y0_ + j*dy;
      for (int i=0; i<=nx_; ++i) {
        nodes(nodeCt, 0) = X0_ + i*dx;
        nodes(nodeCt, 1) = ycoord; 
        ++nodeCt;
      }
    }

  } // computeNodes


  void computeCellToNodeMap() {

    meshCellToNodeMap_ = ROL::makePtr<Intrepid::FieldContainer<int>>(numCells_, 4);
    Intrepid::FieldContainer<int> &ctn = *meshCellToNodeMap_;

    int cellCt = 0;

    for (int j=0; j<ny_; ++j) {
      for (int i=0; i<nx_; ++i) {
        ctn(cellCt, 0) = j*(nx_+1) + i;
        ctn(cellCt, 1) = j*(nx_+1) + (i+1);
        ctn(cellCt, 2) = (j+1)*(nx_+1) + (i+1);
        ctn(cellCt, 3) = (j+1)*(nx_+1) + i;
        ++cellCt;
      }
    }

  } // computeCellToNodeMap


  void computeCellToEdgeMap() {

    meshCellToEdgeMap_ = ROL::makePtr<Intrepid::FieldContainer<int>>(numCells_, 4);
    Intrepid::FieldContainer<int> &cte = *meshCellToEdgeMap_;

    int cellCt = 0;

    for (int j=0; j<ny_; ++j) {
      for (int i=0; i<nx_; ++i) {
        cte(cellCt, 0) = j*(2*nx_+1) + i;
        cte(cellCt, 1) = j*(2*nx_+1) + nx_ + (i+1);
        cte(cellCt, 2) = (j+1)*(2*nx_+1) + i;
        cte(cellCt, 3) = j*(2*nx_+1) + nx_ + i;
        ++cellCt;
      }
    }

  } // computeCellToEdgeMap


  virtual void computeSideSets() {

    meshSideSets_ = ROL::makePtr<std::vector<std::vector<std::vector<int> > >>(1);
    int numSides = 4;
    (*meshSideSets_)[0].resize(numSides);
    (*meshSideSets_)[0][0].resize(nx_);
    (*meshSideSets_)[0][1].resize(ny_);
    (*meshSideSets_)[0][2].resize(nx_);
    (*meshSideSets_)[0][3].resize(ny_);

    for (int i=0; i<nx_; ++i) {
      (*meshSideSets_)[0][0][i] = i;
    }
    for (int i=0; i<ny_; ++i) {
      (*meshSideSets_)[0][1][i] = (i+1)*nx_-1;
    }
    for (int i=0; i<nx_; ++i) {
      (*meshSideSets_)[0][2][i] = i + nx_*(ny_-1);
    }
    for (int i=0; i<ny_; ++i) {
      (*meshSideSets_)[0][3][i] = i*nx_;
    }

  } // computeSideSets


}; // MeshManager_Rectangle


template<class Real>
class MeshManager_Interval : public MeshManager<Real> {

/* Interval geometry [X0,X0+width] */

private:
  Real width_;    // Interval width
  Real X0_;       // x coordinate left corner

  int nx_;

  int numCells_;
  int numNodes_;

  ROL::Ptr<Intrepid::FieldContainer<Real> > meshNodes_;
  ROL::Ptr<Intrepid::FieldContainer<int> >  meshCellToNodeMap_;

  ROL::Ptr<std::vector<std::vector<std::vector<int> > > > meshSideSets_;

public:

  MeshManager_Interval(Teuchos::ParameterList &parlist) {

    // Geometry data
    width_    = parlist.sublist("Geometry").get("Width", 1.0);
    X0_       = parlist.sublist("Geometry").get("X0", 0.0);

    // Mesh data
    nx_       = parlist.sublist("Geometry").get("NX",10);

    numCells_ = nx_;
    numNodes_ = nx_+1;

    // Compute and store mesh data structures
    computeNodes();
    computeCellToNodeMap();
    computeSideSets();

  }

  ROL::Ptr<Intrepid::FieldContainer<Real> > getNodes() const {
    return meshNodes_;
  }

  ROL::Ptr<Intrepid::FieldContainer<int> > getCellToNodeMap() const {
    return meshCellToNodeMap_;
  }

  ROL::Ptr<std::vector<std::vector<std::vector<int> > > > getSideSets(
              const bool verbose = false,
              std::ostream & outStream = std::cout) const {
    return meshSideSets_;
  }

  int getNumCells() const {
    return numCells_;
  }

  int getNumNodes() const {
    return numNodes_;
  } // getNumNodes


private:

  void computeNodes() {

    meshNodes_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(numNodes_,1);
    Intrepid::FieldContainer<Real> &nodes = *meshNodes_;

    Real dx = width_ / nx_;

    for( int i=0; i<nx_+1; ++i ) {
      nodes(i, 0) = X0_ + i*dx;
    }
  } // computeNodes


  void computeCellToNodeMap() {

    meshCellToNodeMap_ = ROL::makePtr<Intrepid::FieldContainer<int>>(numCells_,2);
    Intrepid::FieldContainer<int> &ctn = *meshCellToNodeMap_;

    for( int i=0; i<nx_; ++i ) {
      ctn(i,0) = i;
      ctn(i,1) = i+1;
    }
  } // computeCellToNodeMap


  virtual void computeSideSets() {
    int numSideSets = 2;
    int numSides = 2;
    meshSideSets_
      = ROL::makePtr<std::vector<std::vector<std::vector<int> > >>(numSideSets);

    (*meshSideSets_)[0].resize(numSides);
    (*meshSideSets_)[0][0].resize(1);
    (*meshSideSets_)[0][0][0] = 0;
    (*meshSideSets_)[0][1].resize(0);

    (*meshSideSets_)[1].resize(numSides);
    (*meshSideSets_)[1][0].resize(0);
    (*meshSideSets_)[1][1].resize(1);
    (*meshSideSets_)[1][1][0] = numCells_-1;

  } // computeSideSets

}; // MeshManager_Interval


template<class Real>
class MeshManager_Fractional_Cylinder : public MeshManager<Real> {

/* Line geometry [X0,X0+width] */

private:
  Real width_;    // Interval width
  Real X0_;       // x coordinate left corner

  int nx_;

  Real gamma_;

  int numCells_;
  int numNodes_;
  int numEdges_;

  ROL::Ptr<Intrepid::FieldContainer<Real> > meshNodes_;
  ROL::Ptr<Intrepid::FieldContainer<int> >  meshCellToNodeMap_;
  ROL::Ptr<Intrepid::FieldContainer<int> >  meshCellToEdgeMap_;

  ROL::Ptr<std::vector<std::vector<std::vector<int> > > > meshSideSets_;

public:

  MeshManager_Fractional_Cylinder(Teuchos::ParameterList &parlist) : X0_(0) {
    // Mesh data
    gamma_ = parlist.sublist("Geometry").sublist("Cylinder").get("Grading Parameter",0.5);
    nx_    = parlist.sublist("Geometry").sublist("Cylinder").get("NI",10);
    width_ = parlist.sublist("Geometry").sublist("Cylinder").get("Height",2.0);

    std::cout << "GAMMA: " << gamma_ << "  NX: " << nx_ << "  WIDTH: " << width_ << std::endl;

    numCells_ = nx_;
    numNodes_ = nx_+1;
    numEdges_ = 2;

    std::cout << "NUMCELLS: " << numCells_ << "  NUMNODES: " << numNodes_ << "  NUMEDGES: " << numEdges_ << std::endl;

    // Compute and store mesh data structures
    computeNodes();
    computeCellToNodeMap();
    computeCellToEdgeMap();
    computeSideSets();
  }

  ROL::Ptr<Intrepid::FieldContainer<Real> > getNodes() const {
    return meshNodes_;
  }

  ROL::Ptr<Intrepid::FieldContainer<int> > getCellToNodeMap() const {
    return meshCellToNodeMap_; 
  }

  ROL::Ptr<Intrepid::FieldContainer<int> > getCellToEdgeMap() const {
    return meshCellToEdgeMap_;
  }

  ROL::Ptr<std::vector<std::vector<std::vector<int> > > > getSideSets(
              const bool verbose = false,
              std::ostream & outStream = std::cout) const { 
    return meshSideSets_;
  }

  int getNumCells() const {
    return numCells_;
  }

  int getNumNodes() const {
    return numNodes_;
  } // getNumNodes


  int getNumEdges() const {
    return numEdges_;
  } // getNumEdges


private:

  void computeNodes() {
    meshNodes_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(numNodes_,1);
    Intrepid::FieldContainer<Real> &nodes = *meshNodes_;

    for( int i=0; i<nx_+1; ++i ) {
      nodes(i, 0) = X0_ + std::pow(static_cast<Real>(i)/nx_,gamma_) * width_;
    }
  } // computeNodes


  void computeCellToNodeMap() {
    meshCellToNodeMap_ = ROL::makePtr<Intrepid::FieldContainer<int>>(numCells_,2);
    Intrepid::FieldContainer<int> &ctn = *meshCellToNodeMap_;

    for( int i=0; i<nx_; ++i ) {
      ctn(i,0) = i;
      ctn(i,1) = i+1;
    }
  } // computeCellToNodeMap


  void computeCellToEdgeMap() {
    meshCellToEdgeMap_ = ROL::makePtr<Intrepid::FieldContainer<int>>(numCells_,1);
    Intrepid::FieldContainer<int> &cte = *meshCellToEdgeMap_;

    for( int i=0; i<nx_; ++i ) {
      cte(i,0) = i;
    }
  } // computeCellToEdgeMap


  virtual void computeSideSets() {
    int numSideSets = 2;
    int numSides = 2;
    meshSideSets_
      = ROL::makePtr<std::vector<std::vector<std::vector<int> > >>(numSideSets);

    (*meshSideSets_)[0].resize(numSides);
    (*meshSideSets_)[0][0].resize(1);
    (*meshSideSets_)[0][0][0] = 0;
    (*meshSideSets_)[0][1].resize(0);

    (*meshSideSets_)[1].resize(numSides);
    (*meshSideSets_)[1][0].resize(0);
    (*meshSideSets_)[1][1].resize(1);
    (*meshSideSets_)[1][1][0] = numCells_-1;

  } // computeSideSets

}; // MeshManager_Fractional_Cylinder


/** \class  MeshManager_Brick
    \brief  Mesh construction and mesh management for the
            brick geometry, on hexahedral grids.
*/
template <class Real>
class MeshManager_Brick : public MeshManager<Real> {

/* Brick geometry.

                   ***********************
                 *                     * *
  :--depth--:  *                     *   *
             *                     *     *
           ***********************       *
           *                     *   :   *
           *                     *   |   *
           *                     * height
           *                     *   |
           *                     *   :
           *                     * *
           ***********************
      (X0,Y0,Z0) :--width--:

  (X0,Y0,Z0) denote the *smallest* values of each coordinate.

*/

private:
  Real width_;   // rectangle width
  Real depth_;   // rectangle depth
  Real height_;  // rectangle height
  Real X0_;      // x coordinate of bottom left front corner
  Real Y0_;      // y coordinate of bottom left front corner
  Real Z0_;      // z coordinate of bottom left front corner

  int nx_;
  int ny_;
  int nz_;

  int numCells_;
  int numNodes_;
  int numEdges_;

  ROL::Ptr<Intrepid::FieldContainer<Real> > meshNodes_;
  ROL::Ptr<Intrepid::FieldContainer<int> >  meshCellToNodeMap_;
  ROL::Ptr<Intrepid::FieldContainer<int> >  meshCellToEdgeMap_;

  ROL::Ptr<std::vector<std::vector<std::vector<int> > > >  meshSideSets_;

public:

  MeshManager_Brick(Teuchos::ParameterList &parlist) {
    // Geometry data.
    width_  = parlist.sublist("Geometry").get( "Width", 3.0);
    depth_  = parlist.sublist("Geometry").get( "Depth", 2.0);
    height_ = parlist.sublist("Geometry").get("Height", 1.0);
    X0_     = parlist.sublist("Geometry").get(    "X0", 0.0);
    Y0_     = parlist.sublist("Geometry").get(    "Y0", 0.0);
    Z0_     = parlist.sublist("Geometry").get(    "Z0", 0.0);
    // Mesh data.
    nx_ = parlist.sublist("Geometry").get("NX", 3);
    ny_ = parlist.sublist("Geometry").get("NY", 2);
    nz_ = parlist.sublist("Geometry").get("NZ", 1);
    numCells_ = nx_ * ny_ * nz_;
    numNodes_ = (nx_+1) * (ny_+1) * (nz_+1);
    numEdges_ = ((nx_+1)*ny_ + (ny_+1)*nx_)*(nz_+1) + (nx_+1)*(ny_+1)*nz_;
    // Compute and store mesh data structures.
    computeNodes(); 
    computeCellToNodeMap(); 
    computeCellToEdgeMap();
    computeSideSets();
  }


  ROL::Ptr<Intrepid::FieldContainer<Real> > getNodes() const {
    return meshNodes_;
  }


  ROL::Ptr<Intrepid::FieldContainer<int> > getCellToNodeMap() const {
    return meshCellToNodeMap_;
  }


  ROL::Ptr<Intrepid::FieldContainer<int> > getCellToEdgeMap() const {
    return meshCellToEdgeMap_;
  }


  ROL::Ptr<std::vector<std::vector<std::vector<int> > > > getSideSets(
              const bool verbose = false,
              std::ostream & outStream = std::cout) const { 
    return meshSideSets_;
  }


  int getNumCells() const {
    return numCells_;
  } // getNumCells


  int getNumNodes() const {
    return numNodes_;
  } // getNumNodes


  int getNumEdges() const {
    return numEdges_;
  } // getNumEdges

private:

  void computeNodes() {

    meshNodes_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(numNodes_, 3);
    Intrepid::FieldContainer<Real> &nodes = *meshNodes_;

    Real dx = width_ / nx_;
    Real dy = depth_ / ny_;
    Real dz = height_ / nz_;
    int nodeCt = 0;

    for (int k=0; k<=nz_; ++k) {
      Real zcoord = Z0_ + k*dz;
      for (int j=0; j<=ny_; ++j) {
        Real ycoord = Y0_ + j*dy;
        for (int i=0; i<=nx_; ++i) {
          nodes(nodeCt, 0) = X0_ + i*dx;
          nodes(nodeCt, 1) = ycoord; 
          nodes(nodeCt, 2) = zcoord; 
          ++nodeCt;
        }
      }
    }

  } // computeNodes


  void computeCellToNodeMap() {

    meshCellToNodeMap_ = ROL::makePtr<Intrepid::FieldContainer<int>>(numCells_, 8);
    Intrepid::FieldContainer<int> &ctn = *meshCellToNodeMap_;

    int cellCt = 0;

    int numNodesXY = (nx_+1)*(ny_+1);

    for (int k=0; k<nz_; ++k) {
      for (int j=0; j<ny_; ++j) {
        for (int i=0; i<nx_; ++i) {
          //
          ctn(cellCt, 0) = k*numNodesXY + j*(nx_+1) + i;
          ctn(cellCt, 1) = k*numNodesXY + j*(nx_+1) + (i+1);
          ctn(cellCt, 2) = k*numNodesXY + (j+1)*(nx_+1) + (i+1);
          ctn(cellCt, 3) = k*numNodesXY + (j+1)*(nx_+1) + i;
          //
          ctn(cellCt, 4) = (k+1)*numNodesXY + j*(nx_+1) + i;
          ctn(cellCt, 5) = (k+1)*numNodesXY + j*(nx_+1) + (i+1);
          ctn(cellCt, 6) = (k+1)*numNodesXY + (j+1)*(nx_+1) + (i+1);
          ctn(cellCt, 7) = (k+1)*numNodesXY + (j+1)*(nx_+1) + i;
          //
          ++cellCt;
        }
      }
    }

  } // computeCellToNodeMap


  void computeCellToEdgeMap() {

    meshCellToEdgeMap_ = ROL::makePtr<Intrepid::FieldContainer<int>>(numCells_, 12);
    Intrepid::FieldContainer<int> &cte = *meshCellToEdgeMap_;

    int cellCt = 0;

    int numEdgesXY  = ((nx_+1)*ny_ + (ny_+1)*nx_);
    int numEdgesZ   = (nx_+1)*(ny_+1);
    int numEdgesXYZ = numEdgesXY + numEdgesZ;

    for (int k=0; k<nz_; ++k) {
      for (int j=0; j<ny_; ++j) {
        for (int i=0; i<nx_; ++i) {
          // bottom
          cte(cellCt,  0) = k*numEdgesXYZ + j*(2*nx_+1) + i;
          cte(cellCt,  1) = k*numEdgesXYZ + j*(2*nx_+1) + nx_ + (i+1);
          cte(cellCt,  2) = k*numEdgesXYZ + (j+1)*(2*nx_+1) + i;
          cte(cellCt,  3) = k*numEdgesXYZ + j*(2*nx_+1) + nx_ + i;
          // verticals
          cte(cellCt,  8) = k*numEdgesXYZ + numEdgesXY + j*(nx_+1) + i;
          cte(cellCt,  9) = k*numEdgesXYZ + numEdgesXY + j*(nx_+1) + (i+1);
          cte(cellCt, 10) = k*numEdgesXYZ + numEdgesXY + (j+1)*(nx_+1) + (i+1);
          cte(cellCt, 11) = k*numEdgesXYZ + numEdgesXY + (j+1)*(nx_+1) + i;
          // top
          cte(cellCt,  4) = (k+1)*numEdgesXYZ + j*(2*nx_+1) + i;
          cte(cellCt,  5) = (k+1)*numEdgesXYZ + j*(2*nx_+1) + nx_ + (i+1);
          cte(cellCt,  6) = (k+1)*numEdgesXYZ + (j+1)*(2*nx_+1) + i;
          cte(cellCt,  7) = (k+1)*numEdgesXYZ + j*(2*nx_+1) + nx_ + i;
          ++cellCt;
        }
      }
    }

  } // computeCellToEdgeMap


  virtual void computeSideSets() {

    // single sideset (all of the boundary)
    meshSideSets_ = ROL::makePtr<std::vector<std::vector<std::vector<int> > >>(1);
    // the sideset has six sides with local side ids from 0 to 5
    int numSides = 6;
    (*meshSideSets_)[0].resize(numSides);
    (*meshSideSets_)[0][0].resize(nx_*nz_);
    (*meshSideSets_)[0][1].resize(ny_*nz_);
    (*meshSideSets_)[0][2].resize(nx_*nz_);
    (*meshSideSets_)[0][3].resize(ny_*nz_);
    (*meshSideSets_)[0][4].resize(nx_*ny_);
    (*meshSideSets_)[0][5].resize(nx_*ny_);

    int nxny = nx_*ny_;

    for (int j=0; j<nz_; ++j) {
      for (int i=0; i<nx_; ++i) {
        (*meshSideSets_)[0][0][i+nx_*j] = i+nxny*j;
      }
    }
    for (int j=0; j<nz_; ++j) {
      for (int i=0; i<ny_; ++i) {
        (*meshSideSets_)[0][1][i+ny_*j] = (i+1)*nx_-1 + nxny*j;
      }
    }
    for (int j=0; j<nz_; ++j) {
      for (int i=0; i<nx_; ++i) {
        (*meshSideSets_)[0][2][i+nx_*j] = nx_*(ny_-1) + i + nxny*j;
      }
    }
    for (int j=0; j<nz_; ++j) {
      for (int i=0; i<ny_; ++i) {
        (*meshSideSets_)[0][3][i+ny_*j] = i*nx_ + nxny*j;
      }
    }
    for (int j=0; j<ny_; ++j) {
      for (int i=0; i<nx_; ++i) {
        (*meshSideSets_)[0][4][i+nx_*j] = i + nx_*j;
      }
    }
    for (int j=0; j<ny_; ++j) {
      for (int i=0; i<nx_; ++i) {
        (*meshSideSets_)[0][5][i+nx_*j] = i + nx_*j + nxny*(nz_-1);
      }
    }

  } // computeSideSets


}; // MeshManager_Brick

#endif
