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

/*! \file  dofmanager.hpp
    \brief Given a mesh with connectivity information, sets up data
           structures for the indexing of the global degrees of
           freedom (Dofs).
*/

#ifndef DOFMANAGER_HPP
#define DOFMANAGER_HPP

#include "Teuchos_ParameterList.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_Basis.hpp"
#include "meshmanager.hpp"

template <class Real>
class DofManager {

private:
  ROL::Ptr<MeshManager<Real> > meshManager_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > intrepidBasis_;

  int numCells_;         // number of mesh cells
  int numNodes_;         // number of mesh nodes
  int numEdges_;         // number of mesh edges

  ROL::Ptr<Intrepid::FieldContainer<int> > meshCellToNodeMap_;
  ROL::Ptr<Intrepid::FieldContainer<int> > meshCellToEdgeMap_;

  // Local dof information.
  int numBases_;                                  // number of bases (fields)
  int numLocalNodes_;                             // number of nodes in the basic cell topology
  int numLocalEdges_;                             // number of edges in the basic cell topology
  int numLocalFaces_;                             // number of faces in the basic cell topology
  int numLocalNodeDofs_;                          // number of local (single-cell) node dofs for all bases
  int numLocalEdgeDofs_;                          // number of local (single-cell) edge dofs for all bases
  int numLocalFaceDofs_;                          // number of local (single-cell) face dofs for all bases
  int numLocalDofs_;                              // total number of local (single-cell) face dofs for all bases
  std::vector<int> numDofsPerNode_;               // number of dofs per node in a cell (assumed constant), per basis
  std::vector<int> numDofsPerEdge_;               // number of dofs per edge in a cell (assumed constant), per basis
  std::vector<int> numDofsPerFace_;               // number of dofs per face in a cell (assumed constant), per basis
  std::vector<std::vector<int> > fieldPattern_;   // local indexing of fields into the array [0,1,...,numLocalDofs-1];
                                                  // contains [number of bases] index vectors, where each index vector
                                                  // is of size [basis cardinality]
  // Global dof information. 
  int numNodeDofs_;  // number of global node dofs
  int numEdgeDofs_;  // number of global edge dofs
  int numFaceDofs_;  // number of global face dofs
  int numDofs_;      // total number of global dofs

  ROL::Ptr<Intrepid::FieldContainer<int> > nodeDofs_;  // global node dofs, of size [numNodes_ x numLocalNodeDofs_]
  ROL::Ptr<Intrepid::FieldContainer<int> > edgeDofs_;  // global edge dofs, of size [numEdges_ x numLocalEdgeDofs_]
  ROL::Ptr<Intrepid::FieldContainer<int> > faceDofs_;  // global face dofs, of size [numFaces_ x numLocalFaceDofs_]
  ROL::Ptr<Intrepid::FieldContainer<int> > cellDofs_;  // global cell dofs, of size [numCells_ x numLocalDofs_];
                                                           // ordered by subcell (node, then edge, then face) and basis index

  std::vector<ROL::Ptr<Intrepid::FieldContainer<int> > > fieldDofs_;  // global cell dofs, indexed by field/basis, of size [numCells_ x basis cardinality];
                                                                          // ordered by subcell (node, then edge, then face)

public:

  DofManager(ROL::Ptr<MeshManager<Real> >                                                    &meshManager,
             std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > &intrepidBasis) {

    meshManager_ = meshManager;
    intrepidBasis_ = intrepidBasis;
    numCells_ = meshManager_->getNumCells();
    numNodes_ = meshManager_->getNumNodes();
    numEdges_ = meshManager_->getNumEdges();

    // Mesh data structures.
    meshCellToNodeMap_ = meshManager_->getCellToNodeMap();
    meshCellToEdgeMap_ = meshManager_->getCellToEdgeMap();

    // Local degree-of-freedom footprint.
    numBases_ = static_cast<int>(intrepidBasis_.size());
    numDofsPerNode_.resize(numBases_, 0);
    numDofsPerEdge_.resize(numBases_, 0);
    numDofsPerFace_.resize(numBases_, 0);
    numLocalDofs_ = 0;
    for (int i=0; i<numBases_; ++i) {
      std::vector<std::vector<int> > dofTags = (intrepidBasis_[i])->getAllDofTags();
      for (int j=0; j<(intrepidBasis_[i])->getCardinality(); j++) {
        if (dofTags[j][0] == 0) {
          numDofsPerNode_[i] = dofTags[j][3];
        }
        else if (dofTags[j][0] == 1) {
          numDofsPerEdge_[i] = dofTags[j][3];
        }
        else if (dofTags[j][0] == 2) {
          numDofsPerFace_[i] = dofTags[j][3];
        }
      }
      numLocalDofs_ += (intrepidBasis_[i])->getCardinality();
    }
    numLocalNodeDofs_ = 0;
    numLocalEdgeDofs_ = 0;
    numLocalFaceDofs_ = 0;
    for (int i=0; i<numBases_; ++i) {
      numLocalNodeDofs_ += numDofsPerNode_[i];
      numLocalEdgeDofs_ += numDofsPerEdge_[i];
      numLocalFaceDofs_ += numDofsPerFace_[i];
    }
    numLocalNodes_ = static_cast<int>( (intrepidBasis_[0])->getBaseCellTopology().getVertexCount() );
    numLocalEdges_ = static_cast<int>( (intrepidBasis_[0])->getBaseCellTopology().getEdgeCount() );
    numLocalFaces_ = 1;
    computeFieldPattern();

    // Global data structures.
    computeDofArrays();
    computeCellDofs();
    computeFieldDofs();
  }


  ROL::Ptr<Intrepid::FieldContainer<int> > getNodeDofs() const {
    return nodeDofs_;
  }


  ROL::Ptr<Intrepid::FieldContainer<int> > getEdgeDofs() const {
    return edgeDofs_;
  }


  ROL::Ptr<Intrepid::FieldContainer<int> > getFaceDofs() const {
    return faceDofs_;
  }


  ROL::Ptr<Intrepid::FieldContainer<int> > getCellDofs() const {
    return cellDofs_;
  }


  std::vector<ROL::Ptr<Intrepid::FieldContainer<int> > > getFieldDofs() const {
    return fieldDofs_;
  }


  ROL::Ptr<Intrepid::FieldContainer<int> > getFieldDofs(const int & fieldNumber) const {
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


private:

  void computeDofArrays() {

    nodeDofs_ = ROL::makePtr<Intrepid::FieldContainer<int>>(numNodes_, numLocalNodeDofs_);
    edgeDofs_ = ROL::makePtr<Intrepid::FieldContainer<int>>(numEdges_, numLocalEdgeDofs_);
    faceDofs_ = ROL::makePtr<Intrepid::FieldContainer<int>>(numCells_, numLocalFaceDofs_);

    Intrepid::FieldContainer<int> &nodeDofs = *nodeDofs_;
    Intrepid::FieldContainer<int> &edgeDofs = *edgeDofs_;
    Intrepid::FieldContainer<int> &faceDofs = *faceDofs_;

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
          nodeDofs(i, ++locNodeCt) = ++dofCt;
        }
      }
    }
    numNodeDofs_ = dofCt+1;

    // count edge dofs
    for (int i=0; i<numEdges_; ++i) {
      int locEdgeCt = -1;
      for (int j=0; j<numBases_; ++j) {
        for (int k=0; k<numDofsPerEdge_[j]; ++k) {
          edgeDofs(i, ++locEdgeCt) = ++dofCt;
        }
      }
    }
    numEdgeDofs_ = dofCt+1-numNodeDofs_;

    // count face dofs
    for (int i=0; i<numCells_; ++i) {
      int locFaceCt = -1;
      for (int j=0; j<numBases_; ++j) {
        for (int k=0; k<numDofsPerFace_[j]; ++k) {
          faceDofs(i, ++locFaceCt) = ++dofCt;
        }
      }
    }
    numFaceDofs_ = dofCt+1-numNodeDofs_-numEdgeDofs_;

    numDofs_ = numNodeDofs_+numEdgeDofs_+numFaceDofs_;

  } // computeDofArrays


  void computeCellDofs() {

    cellDofs_ = ROL::makePtr<Intrepid::FieldContainer<int>>(numCells_, numLocalDofs_);

    // Grab object references, for easier indexing.
    Intrepid::FieldContainer<int> &cdofs = *cellDofs_;
    Intrepid::FieldContainer<int> &nodeDofs = *nodeDofs_;
    Intrepid::FieldContainer<int> &edgeDofs = *edgeDofs_;
    Intrepid::FieldContainer<int> &faceDofs = *faceDofs_;
    Intrepid::FieldContainer<int> &ctn = *meshCellToNodeMap_;
    Intrepid::FieldContainer<int> &cte = *meshCellToEdgeMap_;

    for (int i=0; i<numCells_; ++i) {
      int ct = -1;
      for (int j=0; j<numLocalNodes_; ++j) {
        for (int k=0; k<numLocalNodeDofs_; ++k) {
          cdofs(i,++ct) = nodeDofs(ctn(i,j), k);
        }
      }
      for (int j=0; j<numLocalEdges_; ++j) {
        for (int k=0; k<numLocalEdgeDofs_; ++k) {
          cdofs(i,++ct) = edgeDofs(cte(i,j), k);
        }
      }
      for (int j=0; j<numLocalFaces_; ++j) {
        for (int k=0; k<numLocalFaceDofs_; ++k) {
          cdofs(i,++ct) = faceDofs(i, k);
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

  } // computeFieldPattern


  void computeFieldDofs() {

    fieldDofs_.resize(numBases_);

    Intrepid::FieldContainer<int> &cdofs = *cellDofs_;

    for (int fieldNum=0; fieldNum<numBases_; ++fieldNum) { 
      int basisCard = intrepidBasis_[fieldNum]->getCardinality();
      fieldDofs_[fieldNum] = ROL::makePtr<Intrepid::FieldContainer<int>>(numCells_, basisCard);
      Intrepid::FieldContainer<int> &fdofs = *(fieldDofs_[fieldNum]);
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<basisCard; ++j) {
          fdofs(i,j) = cdofs(i, fieldPattern_[fieldNum][j]);
        }
      }  
    }
  } // computeFieldDofs

}; // DofManager

#endif
