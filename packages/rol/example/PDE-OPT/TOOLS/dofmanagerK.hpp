// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  dofmanager.hpp
    \brief Given a mesh with connectivity information, sets up data
           structures for the indexing of the global degrees of
           freedom (Dofs).
*/

#ifndef DOFMANAGERK_HPP
#define DOFMANAGERK_HPP

#include "Teuchos_ParameterList.hpp"
#include "Intrepid2_Basis.hpp"
#include "meshmanagerK.hpp"

template <class Real, class DeviceType>
class DofManager {

public:
  using scalar_view = Kokkos::DynRankView<Real,DeviceType>;
  using int_view = Kokkos::DynRankView<int,DeviceType>;
  using basis_ptr = Intrepid2::BasisPtr<DeviceType, Real, Real>;

private:
  ROL::Ptr<MeshManager<Real,DeviceType> > meshManager_;
  std::vector<basis_ptr> intrepidBasis_;

  int cellDim_;          // cell dimension

  int numCells_;         // number of mesh cells
  int numNodes_;         // number of mesh nodes
  int numEdges_;         // number of mesh edges
  int numFaces_;         // number of mesh faces

  int_view meshCellToNodeMap_;
  int_view meshCellToEdgeMap_;
  int_view meshCellToFaceMap_;

  // Local dof information.
  int numBases_;                                  // number of bases (fields)
  int numLocalNodes_;                             // number of nodes in the basic cell topology
  int numLocalEdges_;                             // number of edges in the basic cell topology
  int numLocalFaces_;                             // number of faces in the basic cell topology
  int numLocalVoids_;                             // number of voids in the basic cell topology, =1
  int numLocalNodeDofs_;                          // number of local (single-cell) node dofs for all bases
  int numLocalEdgeDofs_;                          // number of local (single-cell) edge dofs for all bases
  int numLocalFaceDofs_;                          // number of local (single-cell) face dofs for all bases
  int numLocalVoidDofs_;                          // number of local (single-cell) face dofs for all bases
  int numLocalDofs_;                              // total number of local (single-cell) face dofs for all bases
  std::vector<int> numDofsPerNode_;               // number of dofs per node in a cell (assumed constant), per basis
  std::vector<int> numDofsPerEdge_;               // number of dofs per edge in a cell (assumed constant), per basis
  std::vector<int> numDofsPerFace_;               // number of dofs per face in a cell (assumed constant), per basis
  std::vector<int> numDofsPerVoid_;               // number of dofs per void in a cell (assumed constant), per basis
  std::vector<std::vector<int> > fieldPattern_;   // local indexing of fields into the array [0,1,...,numLocalDofs-1];
                                                  // contains [number of bases] index vectors, where each index vector
                                                  // is of size [basis cardinality]
  // Global dof information.
  int numNodeDofs_;  // number of global node dofs
  int numEdgeDofs_;  // number of global edge dofs
  int numFaceDofs_;  // number of global face dofs
  int numVoidDofs_;  // number of global void dofs
  int numDofs_;      // total number of global dofs

  int_view nodeDofs_;  // global node dofs, of size [numNodes_ x numLocalNodeDofs_]
  int_view edgeDofs_;  // global edge dofs, of size [numEdges_ x numLocalEdgeDofs_]
  int_view faceDofs_;  // global face dofs, of size [numFaces_ x numLocalFaceDofs_]
  int_view voidDofs_;  // global face dofs, of size [numFaces_ x numLocalFaceDofs_]
  int_view cellDofs_;  // global cell dofs, of size [numCells_ x numLocalDofs_];
                                                           // ordered by subcell (node, then edge, then face) and basis index
  std::vector<int> mapToIntrepidPattern_;
  std::vector<int> mapToFieldPattern_;
  std::vector<int_view> fieldDofs_;  // global cell dofs, indexed by field/basis, of size [numCells_ x basis cardinality];
                                                                          // ordered by subcell (node, then edge, then face)

public:

  DofManager(ROL::Ptr<MeshManager<Real,DeviceType> > &meshManager,
             std::vector<basis_ptr> intrepidBasis) {

    meshManager_ = meshManager;
    intrepidBasis_ = intrepidBasis;
    cellDim_  = intrepidBasis_[0]->getBaseCellTopology().getDimension();
    numCells_ = meshManager_->getNumCells();
    numNodes_ = meshManager_->getNumNodes();
    numEdges_ = meshManager_->getNumEdges();
    numFaces_ = meshManager_->getNumFaces();

    // Mesh data structures.
    meshCellToNodeMap_ = meshManager_->getCellToNodeMap();
    meshCellToEdgeMap_ = meshManager_->getCellToEdgeMap();
    meshCellToFaceMap_ = meshManager_->getCellToFaceMap();

    // Local degree-of-freedom footprint.
    numBases_ = static_cast<int>(intrepidBasis_.size());
    numDofsPerNode_.resize(numBases_, 0);
    numDofsPerEdge_.resize(numBases_, 0);
    numDofsPerFace_.resize(numBases_, 0);
    numDofsPerVoid_.resize(numBases_, 0);
    numLocalDofs_ = 0;
    for (int i=0; i<numBases_; ++i) {
      auto dofTags = (intrepidBasis_[i])->getAllDofTags();
      for (int j=0; j<(intrepidBasis_[i])->getCardinality(); ++j) {
        if (dofTags(j,0) == 0) {
          numDofsPerNode_[i] = dofTags(j,3);
        }
        else if (dofTags(j,0) == 1) {
          if (cellDim_ == 1) { // 1D
            numDofsPerVoid_[i] = dofTags(j,3);
          }
          else { // 2D, 3D
            numDofsPerEdge_[i] = dofTags(j,3);
          }
        }
        else if (dofTags(j,0) == 2) {
          if (cellDim_ == 2) { // 2D
            numDofsPerVoid_[i] = dofTags(j,3);
          }
          else { // 3D
            numDofsPerFace_[i] = dofTags(j,3);
          }
        }
        else if (dofTags(j,0) == 3) {
          numDofsPerVoid_[i] = dofTags(j,3);
        }
      }
      numLocalDofs_ += (intrepidBasis_[i])->getCardinality();
    }
    numLocalNodeDofs_ = 0;
    numLocalEdgeDofs_ = 0;
    numLocalFaceDofs_ = 0;
    numLocalVoidDofs_ = 0;
    for (int i=0; i<numBases_; ++i) {
      numLocalNodeDofs_ += numDofsPerNode_[i];
      numLocalEdgeDofs_ += numDofsPerEdge_[i];
      numLocalFaceDofs_ += numDofsPerFace_[i];
      numLocalVoidDofs_ += numDofsPerVoid_[i];
    }
    numLocalNodes_ = static_cast<int>( (intrepidBasis_[0])->getBaseCellTopology().getVertexCount() );
    numLocalEdges_ = static_cast<int>( (intrepidBasis_[0])->getBaseCellTopology().getEdgeCount() );
    numLocalFaces_ = static_cast<int>( (intrepidBasis_[0])->getBaseCellTopology().getFaceCount() );
    numLocalVoids_ = 1;
    computeFieldPattern();

    // Global data structures.
    computeDofArrays();
    computeCellDofs();
    computeFieldDofs();
    computeDofTransforms();
  }


  int_view getNodeDofs() const {
    return nodeDofs_;
  }


  int_view getEdgeDofs() const {
    return edgeDofs_;
  }


  int_view getFaceDofs() const {
    return faceDofs_;
  }


  int_view getVoidDofs() const {
    return voidDofs_;
  }


  int_view getCellDofs() const {
    return cellDofs_;
  }


  int_view getFieldDofs() const {
    return fieldDofs_;
  }


  int_view getFieldDofs(const int & fieldNumber) const {
    return fieldDofs_[fieldNumber];
  }


  int getNumNodeDofs() const {
    return numNodeDofs_;
  }


  int getNumEdgeDofs() const {
    return numEdgeDofs_;
  }


  int getNumFaceDofs() const {
    return numFaceDofs_;
  }


  int getNumVoidDofs() const {
    return numVoidDofs_;
  }


  int getNumDofs() const {
    return numDofs_;
  }


  int getNumFields() const {
    return numBases_;
  }


  int getLocalFieldSize() const {
    return numLocalDofs_;
  }


  int getLocalFieldSize(const int & fieldNumber) const {
    return (intrepidBasis_[fieldNumber])->getCardinality();
  }


  std::vector<std::vector<int> > getFieldPattern() const {
    return fieldPattern_;
  }


  std::vector<int> getFieldPattern(const int & fieldNumber) const {
    return fieldPattern_[fieldNumber];
  }

  void transformToIntrepidPattern(const scalar_view array) const {
    if ( array.is_allocated() ) {
      int rank = array.rank();
      int nc   = array.extent(0);
      if ( rank == 2 ) {
        int nf = array.extent(1);
        scalar_view tmp("tmp2d", nc, nf);
        for (int c = 0; c < nc; ++c) {
          for (int f = 0; f < nf; ++f) {
            tmp(c, mapToIntrepidPattern_[f]) = array(c, f);
          }
        }
        Kokkos::deep_copy(array, tmp);
      }
      else if (rank == 3 ) {
        int nf1 = array.extent_int(1);
        int nf2 = array.extent_int(2);
        scalar_view tmp("tmp3d", nc, nf1, nf2);
        for (int c = 0; c < nc; ++c) {
          for (int f1 = 0; f1 < nf1; ++f1) {
            for (int f2 = 0; f2 < nf2; ++f2) {
              tmp(c, mapToIntrepidPattern_[f1], mapToIntrepidPattern_[f2]) = array(c, f1, f2);
            }
          }
        }
        Kokkos::deep_copy(array, tmp);
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
          ">>> PDE-OPT/TOOLS/dofmanager.hpp (transformToIntrepidPattern): Input array rank not 2 or 3!");
      }
    }
  }

  int mapToFieldPattern(int f) const {
    return mapToFieldPattern_[f];
  }

  void transformToFieldPattern(const scalar_view array) const {
    if ( array.is_allocated() ) {
      int rank = array.rank();
      int nc   = array.extent_int(0);
      if ( rank == 2 ) {
        int nf = array.extent_int(1);
        scalar_view tmp("tmp2d", nc, nf);
        for (int c = 0; c < nc; ++c) {
          for (int f = 0; f < nf; ++f) {
            tmp(c, mapToFieldPattern_[f]) = array(c, f);
          }
        }
        Kokkos::deep_copy(array, tmp);
      }
      else if (rank == 3 ) {
        int nf1 = array->dimension(1);
        int nf2 = array->dimension(2);
        scalar_view tmp("tmp3d", nc, nf1, nf2);
        for (int c = 0; c < nc; ++c) {
          for (int f1 = 0; f1 < nf1; ++f1) {
            for (int f2 = 0; f2 < nf2; ++f2) {
              tmp(c, mapToFieldPattern_[f1], mapToFieldPattern_[f2]) = array(c, f1, f2);
            }
          }
        }
        Kokkos::deep_copy(array, tmp);
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
          ">>> PDE-OPT/TOOLS/dofmanager.hpp (transformToFieldPattern): Input array rank not 2 or 3!");
      }
    }
  }

private:

  void computeDofArrays() {

    nodeDofs_ = int_view("nodeDofs_", numNodes_, numLocalNodeDofs_);
    edgeDofs_ = int_view("edgeDofs_", numEdges_, numLocalEdgeDofs_);
    faceDofs_ = int_view("faceDofs_", numFaces_, numLocalFaceDofs_);
    voidDofs_ = int_view("voidDofs_", numCells_, numLocalVoidDofs_);

    int dofCt = -1;

    //
    // This is the independent node id --> edge id --> cell id ordering.
    // For example, for a Q2-Q2-Q1 basis (Taylor-Hood), we would have:
    // 
    //     Q2   Q2   Q1
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

    // count node dofs
    for (int i=0; i<numNodes_; ++i) {
      int locNodeCt = -1;
      for (int j=0; j<numBases_; ++j) {
        for (int k=0; k<numDofsPerNode_[j]; ++k) {
          nodeDofs_(i, ++locNodeCt) = ++dofCt;
        }
      }
    }
    numNodeDofs_ = dofCt+1;

    // count edge dofs
    for (int i=0; i<numEdges_; ++i) {
      int locEdgeCt = -1;
      for (int j=0; j<numBases_; ++j) {
        for (int k=0; k<numDofsPerEdge_[j]; ++k) {
          edgeDofs_(i, ++locEdgeCt) = ++dofCt;
        }
      }
    }
    numEdgeDofs_ = dofCt+1-numNodeDofs_;

    // count face dofs
    for (int i=0; i<numFaces_; ++i) {
      int locFaceCt = -1;
      for (int j=0; j<numBases_; ++j) {
        for (int k=0; k<numDofsPerFace_[j]; ++k) {
          faceDofs_(i, ++locFaceCt) = ++dofCt;
        }
      }
    }
    numFaceDofs_ = dofCt+1-numNodeDofs_-numEdgeDofs_;

    // count void dofs
    for (int i=0; i<numCells_; ++i) {
      int locVoidCt = -1;
      for (int j=0; j<numBases_; ++j) {
        for (int k=0; k<numDofsPerVoid_[j]; ++k) {
          voidDofs_(i, ++locVoidCt) = ++dofCt;
        }
      }
    }
    numVoidDofs_ = dofCt+1-numNodeDofs_-numEdgeDofs_-numFaceDofs_;

    numDofs_ = numNodeDofs_+numEdgeDofs_+numFaceDofs_+numVoidDofs_;

  } // computeDofArrays


  void computeCellDofs() {

    cellDofs_ = int_view("cellDofs_", numCells_, numLocalDofs_);

    for (int i=0; i<numCells_; ++i) {
      int ct = -1;
      for (int j=0; j<numLocalNodes_; ++j) {
        for (int k=0; k<numLocalNodeDofs_; ++k) {
          cellDofs_(i,++ct) = nodeDofs_(meshCellToNodeMap_(i,j), k);
        }
      }
      for (int j=0; j<numLocalEdges_; ++j) {
        for (int k=0; k<numLocalEdgeDofs_; ++k) {
          cellDofs_(i,++ct) = edgeDofs_(meshCellToEdgeMap_(i,j), k);
        }
      }
      for (int j=0; j<numLocalFaces_; ++j) {
        for (int k=0; k<numLocalFaceDofs_; ++k) {
          cellDofs_(i,++ct) = faceDofs_(meshCellToFaceMap_(i,j), k);
        }
      }
      for (int j=0; j<numLocalVoids_; ++j) {
        for (int k=0; k<numLocalVoidDofs_; ++k) {
          cellDofs_(i,++ct) = voidDofs_(i, k);
        }
      }
    }
  } // computeCellDofs


  void computeFieldPattern() {

    fieldPattern_.resize(numBases_);

    int dofCt = -1;

    // count node dofs
    for (int i=0; i<numLocalNodes_; ++i) {
      for (int j=0; j<numBases_; ++j) {
        for (int k=0; k<numDofsPerNode_[j]; ++k) {
          fieldPattern_[j].push_back(++dofCt);
        }
      }
    }

    // count edge dofs
    for (int i=0; i<numLocalEdges_; ++i) {
      for (int j=0; j<numBases_; ++j) {
        for (int k=0; k<numDofsPerEdge_[j]; ++k) {
          fieldPattern_[j].push_back(++dofCt);
        }
      }
    }

    // count face dofs
    for (int i=0; i<numLocalFaces_; ++i) {
      for (int j=0; j<numBases_; ++j) {
        for (int k=0; k<numDofsPerFace_[j]; ++k) {
          fieldPattern_[j].push_back(++dofCt);
        }
      }
    }

    // count void dofs
    for (int i=0; i<numLocalVoids_; ++i) {
      for (int j=0; j<numBases_; ++j) {
        for (int k=0; k<numDofsPerVoid_[j]; ++k) {
          fieldPattern_[j].push_back(++dofCt);
        }
      }
    }

  } // computeFieldPattern


  void computeFieldDofs() {

    fieldDofs_.resize(numBases_);

    for (int fieldNum=0; fieldNum<numBases_; ++fieldNum) { 
      int basisCard = intrepidBasis_[fieldNum]->getCardinality();
      int_view fdofs("fieldDofs_", numCells_, basisCard);
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<basisCard; ++j) {
          fdofs(i,j) = cellDofs_(i, fieldPattern_[fieldNum][j]);
        }
      }
      fieldDofs_[fieldNum] = fdofs;
    }
  } // computeFieldDofs

  void computeDofTransforms(void) {
    // Build basis maps
    std::vector<std::vector<int> > map2IP(numBases_);
    std::vector<std::vector<int> > map2FP(numBases_);
    shards::CellTopology cellType = intrepidBasis_[0]->getBaseCellTopology();
    int nv = cellType.getVertexCount();
    for (int f=0; f<numBases_; ++f) { 
      int basisDeg = intrepidBasis_[f]->getDegree();
      if (cellDim_ == 1) {
        if (basisDeg == 0) {
          map2IP[f] = {0};
          map2FP[f] = {0};
        }
        else if (basisDeg == 1) {
          map2IP[f] = {0, 1};
          map2FP[f] = {0, 1};
        }
        else if (basisDeg == 2) {
          // C2 LINE basis does *not* follow the LINE<3> topology node order
          map2IP[f] = {0, 2, 1};
          map2FP[f] = {0, 2, 1};
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
            ">>> PDE-OPT/TOOLS/dofmanager.hpp (ComputeDofTransforms): basisDeg > 2");
        }
      }
      else if (cellDim_ == 2) {
        if (basisDeg == 0) {
          map2IP[f] = {0};
          map2FP[f] = {0};
        }
        else if (basisDeg == 1) {
          map2IP[f].resize(nv);
          map2FP[f].resize(nv);
          for (int i = 0; i < nv; ++i) {
            map2IP[f][i] = i;
            map2FP[f][i] = i;
          }
          //map2IP[f] = {0, 1, 2, 3}; // C1 QUAD
          //map2FP[f] = {0, 1, 2, 3}; // C1 QUAD
        }
        else if (basisDeg == 2) {
          int nbfs = (nv==3 ? 6 : 9);
          map2IP[f].resize(nbfs);
          map2FP[f].resize(nbfs);
          for (int i = 0; i < nbfs; ++i) {
            map2IP[f][i] = i;
            map2FP[f][i] = i;
          }
          //map2IP[f] = {0, 1, 2, 3, 4, 5, 6, 7, 8}; // C2 QUAD
          //map2FP[f] = {0, 1, 2, 3, 4, 5, 6, 7, 8}; // C2 QUAD
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
            ">>> PDE-OPT/TOOLS/dofmanager.hpp (ComputeDofTransforms): basisDeg > 2");
        }
      }
      else if (cellDim_ == 3) {
        if (basisDeg == 0) {
          map2IP[f] = {0};
          map2FP[f] = {0};
        }
        else if (basisDeg == 1) {
          map2IP[f].resize(nv);
          map2FP[f].resize(nv);
          for (int i = 0; i < nv; ++i) {
            map2IP[f][i] = i;
            map2FP[f][i] = i;
          }
          //map2IP[f] = {0, 1, 2, 3, 4, 5, 6, 7}; // C1 HEX
          //map2FP[f] = {0, 1, 2, 3, 4, 5, 6, 7}; // C1 HEX
        }
        else if (basisDeg == 2) {
          if (nv == 4) { // C2 TET basis follows the TET<10> topology node order
            int nbfs = 10;
            map2IP[f].resize(nbfs);
            map2FP[f].resize(nbfs);
            for (int i = 0; i < nbfs; ++i) {
              map2IP[f][i] = i;
              map2FP[f][i] = i;
            }
          }
          else if (nv == 8) { // C2 HEX basis does *not* follow the HEX<27> topology node order
            map2IP[f] = {0, 1, 2, 3, 4 ,5, 6, 7,                       // Nodes
                         8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15, // Edges
                         25, 24, 26, 23, 21, 22,                       // Faces
                         20};                                          // Voids
            map2FP[f] = {0, 1, 2, 3, 4 ,5, 6, 7,                       // Nodes
                         8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15, // Edges
                         26,                                           // Voids
                         24, 25, 23, 21, 20, 22};                      // Faces
          }
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
            ">>> PDE-OPT/TOOLS/dofmanager.hpp (ComputeDofTransforms): basisDeg > 2");
        }
      }
    }
    // Apply transformation to ind
    mapToIntrepidPattern_.clear(); mapToIntrepidPattern_.resize(numDofs_);
    mapToFieldPattern_.clear();    mapToFieldPattern_.resize(numDofs_);
    for (int i = 0; i < numBases_; ++i) {
      for (int j = 0; j < static_cast<int>(map2IP[i].size()); ++j) {
        mapToIntrepidPattern_[fieldPattern_[i][j]] = fieldPattern_[i][map2IP[i][j]];
        mapToFieldPattern_[fieldPattern_[i][j]]    = fieldPattern_[i][map2FP[i][j]];
      }
    }
  }

}; // DofManager

#endif
