// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  DofManager.hpp
    \brief Sets up data structures for the indexing of the global degrees
           of freedom (Dof's).
*/

#include "Teuchos_ParameterList.hpp"
#include "Intrepid_FieldContainer.hpp"

template <class Real>
class DofManager {

typedef Intrepid::FieldContainer<Real> FCR;
typedef Intrepid::FieldContainer<int>  FCI;

enum FieldType {
  VELX, VELY, PRES
};

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

  int ref_;

  int numCells_;
  int numNodes_;
  int numEdges_;
  int numMeshEntities_;

  FCR meshNodes_;
  FCI meshCellToNodeMap_;
  FCI meshCellToEdgeMap_;

  int numDofs_;
  int numPres_;
  int numVelX_;
  int numVelY_;
  int numLocalDofs_;
  int numLocalPres_;
  int numLocalVelX_;
  int numLocalVelY_;

  FCI dofArray_;
  FCI cellDofs_;
  FCI cellPres_;
  FCI cellVelX_;
  FCI cellVelY_;

public:

  DofManager(Teuchos::ParameterList &parlist) {
    // Geometry data.
    channelH_ = parlist.sublist("Navier Stokes").get(   "Channel height", 1.0);
    channelW_ = parlist.sublist("Navier Stokes").get(    "Channel width", 8.0);
    stepH_    = parlist.sublist("Navier Stokes").get(      "Step height", 0.5);
    stepW_    = parlist.sublist("Navier Stokes").get(       "Step width", 1.0);
    observeW_ = parlist.sublist("Navier Stokes").get("Observation width", 3.0);
    // Mesh data.
    ref_ = parlist.sublist("Navier Stokes").get(      "Refinement level", 1);
    nx1_ = parlist.sublist("Navier Stokes").get( "Observation region NX", 4*ref_);
    ny1_ = parlist.sublist("Navier Stokes").get( "Observation region NY", 5*ref_);
    nx2_ = parlist.sublist("Navier Stokes").get("Laminar flow region NX", 2*ref_);
    ny2_ = ny1_;
    nx3_ = parlist.sublist("Navier Stokes").get(      "Inflow region NX", 1*ref_);
    ny3_ = parlist.sublist("Navier Stokes").get(      "Inflow region NY", 2*ref_);
    nx4_ = nx1_;
    ny4_ = ny3_;
    nx5_ = nx2_;
    ny5_ = ny3_;
    numCells_ = (nx1_ + nx2_)*ny1_  +  (nx3_ + nx1_ + nx2_)*ny3_;
    numNodes_ = (nx1_ + nx2_ + 1)*ny1_  +  (nx3_ + nx1_ + nx2_ + 1)*(ny3_ + 1);
    numEdges_ = (2*(nx1_ + nx2_)+1)*ny1_ + (2*(nx3_ + nx1_ + nx2_)+1)*ny3_ + (nx3_ + nx1_ + nx2_);
    numMeshEntities_ = numCells_ + numNodes_ + numEdges_;
    // Mesh and Dof data structures.
    numLocalPres_ = 4;
    numLocalVelX_ = 9;
    numLocalVelY_ = 9;
    numLocalDofs_ = numLocalPres_ + numLocalVelX_ + numLocalVelY_;
    computeMeshNodes(meshNodes_); 
    computeMeshCellToNodeMap(meshCellToNodeMap_); 
    computeMeshCellToEdgeMap(meshCellToEdgeMap_); 
    computeDofArray(dofArray_, numDofs_);
    computeCellDofs(cellDofs_, dofArray_, meshCellToNodeMap_, meshCellToEdgeMap_);
    computeCellPres(cellPres_, cellDofs_);
    computeCellVelX(cellVelX_, cellDofs_);
    computeCellVelY(cellVelY_, cellDofs_);
  }


  void computeMeshNodes(FCR &nodes) const {

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

  }  // computeMeshNodes


  void computeMeshCellToNodeMap(FCI &ctn) {

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

  } // computeMeshCellToNodeMap


  void computeMeshCellToEdgeMap(FCI &cte) {

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

  } // computeMeshCellToEdgeMap


  void computeDofArray(FCI &dof, int &numDofs, const int &globOrder=1) {

    dof.resize(numMeshEntities_, 3);
    int dofCt = -1;

    if (globOrder == 1) {
      //
      // This is the independent node id --> edge id --> cell id ordering:
      // 
      //     VelX VelY Pres
      //     -------------- 
      // n1  1    1    1
      // n2  1    1    1
      // ...
      // e1  1    1    0
      // e2  1    1    0
      // ...
      // c1  1    1    0
      // c2  1    1    0
      // ...
      //
      for (int i=0; i<numNodes_; ++i) {
        // dof(i, 1) --> 1 Velocity X
        // dof(i, 1) --> 1 Velocity Y
        // dof(i, 2) --> 1 Pressure
        dof(i, 0) = ++dofCt; // VelX
        dof(i, 1) = ++dofCt; // VelY
        dof(i, 2) = ++dofCt; // Pres
      }
      for (int i=numNodes_; i<numNodes_+numEdges_; ++i) {
        // dof(i, 1) --> 1 Velocity X
        // dof(i, 1) --> 1 Velocity Y
        // dof(i, 2) --> 0 Pressure
        dof(i, 0) = ++dofCt; // VelX
        dof(i, 1) = ++dofCt; // VelY
        dof(i, 2) = -1;
      }
      for (int i=numNodes_+numEdges_; i<numMeshEntities_; ++i) {
        // dof(i, 1) --> 1 Velocity X
        // dof(i, 1) --> 1 Velocity Y
        // dof(i, 2) --> 0 Pressure
        dof(i, 0) = ++dofCt; // VelX
        dof(i, 1) = ++dofCt; // VelY
        dof(i, 2) = -1;
      }
    }
    else if (globOrder==2) {
      //
      // This is the cell-driven cell node id --> cell edge id --> cell id ordering:
      // 
      //         VelX VelY Pres
      //         -------------- 
      // cell 1
      // n1      1    1    1
      // n2      1    1    1
      // ...
      // e1      1    1    0
      // e2      1    1    0
      // ...
      // c1      1    1    0
      // c2      1    1    0
      // ...
      // cell 2
      // ...
      // (unique ids only, of course)
      //
      // TO BE IMPLEMENTED!!!
    }
    numDofs = dofCt + 1;
  } // computeDofArray


  void computeCellDofs(FCI &cdofs, const FCI &dofs, const FCI &ctn, const FCI &cte) {

    cdofs.resize(numCells_, numLocalDofs_);

    for (int i=0; i<numCells_; ++i) {
      int ct = -1;
      int gid0 = ctn(i,0);
      int gid1 = ctn(i,1);
      int gid2 = ctn(i,2);
      int gid3 = ctn(i,3);
      int gid4 = numNodes_ + cte(i,0);
      int gid5 = numNodes_ + cte(i,1);
      int gid6 = numNodes_ + cte(i,2);
      int gid7 = numNodes_ + cte(i,3);
      int gid8 = numNodes_ + numEdges_ + i;
      cdofs(i,++ct) = dofs(gid0, 0);
      cdofs(i,++ct) = dofs(gid0, 1);
      cdofs(i,++ct) = dofs(gid0, 2);
      cdofs(i,++ct) = dofs(gid1, 0);
      cdofs(i,++ct) = dofs(gid1, 1);
      cdofs(i,++ct) = dofs(gid1, 2);
      cdofs(i,++ct) = dofs(gid2, 0);
      cdofs(i,++ct) = dofs(gid2, 1);
      cdofs(i,++ct) = dofs(gid2, 2);
      cdofs(i,++ct) = dofs(gid3, 0);
      cdofs(i,++ct) = dofs(gid3, 1);
      cdofs(i,++ct) = dofs(gid3, 2);
      cdofs(i,++ct) = dofs(gid4, 0);
      cdofs(i,++ct) = dofs(gid4, 1);
      cdofs(i,++ct) = dofs(gid5, 0);
      cdofs(i,++ct) = dofs(gid5, 1);
      cdofs(i,++ct) = dofs(gid6, 0);
      cdofs(i,++ct) = dofs(gid6, 1);
      cdofs(i,++ct) = dofs(gid7, 0);
      cdofs(i,++ct) = dofs(gid7, 1);
      cdofs(i,++ct) = dofs(gid8, 0);
      cdofs(i,++ct) = dofs(gid8, 1);
    }
  } // computeCellDofs


  void computeCellPres(FCI &cpres, const FCI &cdofs) {
    cpres.resize(numCells_, numLocalPres_);
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalPres_; ++j) {
        cpres(i,j) = cdofs(i, j*3+2);
      }
    }
  } // computeCellPres


  void computeCellVelX(FCI &cvelx, const FCI &cdofs) {
    cvelx.resize(numCells_, numLocalVelX_);
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalPres_; ++j) {
        cvelx(i,j) = cdofs(i, j*3);
      }
      for (int k=0; k<numLocalVelX_-numLocalPres_; ++k) {
        cvelx(i,k+numLocalPres_) = cdofs(i, 4*3 + k*2);
      }
    }
  } // computeCellVelX


  void computeCellVelY(FCI &cvely, const FCI &cdofs) {
    cvely.resize(numCells_, numLocalVelY_);
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalPres_; ++j) {
        cvely(i,j) = cdofs(i, j*3+1);
      }
      for (int k=0; k<numLocalVelY_-numLocalPres_; ++k) {
        cvely(i,k+numLocalPres_) = cdofs(i, 4*3 + k*2 + 1);
      }
    }
  } // computeCellVelY


  std::vector<int>& getCellDofs(const int& cid, const FieldType& ft) const {
    std::vector<int> dofVec;
    int nPres = 4;
    int nVelX = 9;
    int nVelY = 9;
    switch(ft) {
      case VELX:
      break;

      case VELY:
      break;

      case PRES:
        for (int i=0; i<nPres; ++i) {
          dofVec.push_back(cellDofs_(cid));
        }
      break;
    }
    return dofVec;
  } // getDofArray


  void setRefinementLevel(const int &refLevel) {
    ref_ = refLevel;
  } // setRefinementLevel


  int getNumCells() const {
    return numCells_;
  } // getNumCells


  int getNumNodes() const {
    return numNodes_;
  } // getNumNodes


  int getNumEdges() const {
    return numEdges_;
  } // getNumEdges


  const FCR& getMeshNodes() const {
    return meshNodes_;
  } // getMeshNodes


  const FCI& getMeshCellToNodeMap() const {
    return meshCellToNodeMap_;
  } // getMeshCellToNodeMap


  const FCI& getMeshCellToEdgeMap() const {
    return meshCellToEdgeMap_;
  } // getMeshCellToEdgeMap


  int getNumDofs() const {
    return numDofs_;
  } // getNumDofs


  const FCI& getDofArray() const {
    return dofArray_;
  } // getDofArray


  const FCI& getCellDofs() const {
    return cellDofs_;
  } // getCellDofs


  const FCI& getCellPres() const {
    return cellPres_;
  } // getCellPres


  const FCI& getCellVelX() const {
    return cellVelX_;
  } // getCellVelX


  const FCI& getCellVelY() const {
    return cellVelY_;
  } // getCellVelY

};
