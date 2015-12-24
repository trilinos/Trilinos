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
  Teuchos::RCP<MeshManager<Real> > meshManager_;
  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > intrepidBasis_;

  int numCells_;         // number of mesh cells
  int numNodes_;         // number of mesh nodes
  int numEdges_;         // number of mesh edges
  int numMeshEntities_;  // number of all mesh entities

  Intrepid::FieldContainer<int> meshCellToNodeMap_;
  Intrepid::FieldContainer<int> meshCellToEdgeMap_;

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
                                                  // contains (number of bases) index vectors, where each index vector
                                                  // is of size (basis cardinality)
  // Global dof information. 
  int numNodeDofs_;  // number of global node dofs
  int numEdgeDofs_;  // number of global edge dofs
  int numFaceDofs_;  // number of global face dofs
  int numDofs_;      // total number of global dofs

  Intrepid::FieldContainer<int> nodeDofs_;  // global node dofs, of size numNodes_ x numLocalNodeDofs_
  Intrepid::FieldContainer<int> edgeDofs_;  // global edge dofs, of size numEdges_ x numLocalEdgeDofs_
  Intrepid::FieldContainer<int> faceDofs_;  // global face dofs, of size numFaces_ x numLocalFaceDofs_
  Intrepid::FieldContainer<int> cellDofs_;  // global cell dofs, of size numCells_ x numLocalDofs_;
                                            // ordered by subcell (node, then edge, then face) and basis index

public:

  DofManager(Teuchos::RCP<MeshManager<Real> >                                                    &meshManager,
             std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > &intrepidBasis) {
    meshManager_ = meshManager;
    intrepidBasis_ = intrepidBasis;
    numCells_ = meshManager_->getNumCells();
    numNodes_ = meshManager_->getNumNodes();
    numEdges_ = meshManager_->getNumEdges();
    numMeshEntities_ = numCells_ + numNodes_ + numEdges_;

    // Mesh data structures.
    meshManager_->getCellToNodeMap(meshCellToNodeMap_);
    meshManager_->getCellToEdgeMap(meshCellToEdgeMap_);

    // Local degree-of-freedom footprint.
    numBases_ = static_cast<int>(intrepidBasis.size());
    numDofsPerNode_.resize(numBases_, 0);
    numDofsPerEdge_.resize(numBases_, 0);
    numDofsPerFace_.resize(numBases_, 0);
    numLocalDofs_ = 0;
    for (int i=0; i<numBases_; ++i) {
      std::vector<std::vector<int> > dofTags = (intrepidBasis[i])->getAllDofTags();
      for (int j=0; j<(intrepidBasis[i])->getCardinality(); j++) {
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
      numLocalDofs_ += (intrepidBasis[i])->getCardinality();
    }
    numLocalNodeDofs_ = 0;
    numLocalEdgeDofs_ = 0;
    numLocalFaceDofs_ = 0;
    for (int i=0; i<numBases_; ++i) {
      numLocalNodeDofs_ += numDofsPerNode_[i];
      numLocalEdgeDofs_ += numDofsPerEdge_[i];
      numLocalFaceDofs_ += numDofsPerFace_[i];
    }
    std::cout << "Number of dofs on node, per basis: ";
    for (int i=0; i<numBases_; ++i) {
      std::cout << numDofsPerNode_[i] << "   ";
    }
    std::cout << "\nNumber of dofs on edge, per basis: ";
    for (int i=0; i<numBases_; ++i) {
      std::cout << numDofsPerEdge_[i] << "   ";
    }
    std::cout << "\nNumber of dofs on cell, per basis: ";
    for (int i=0; i<numBases_; ++i) {
      std::cout << numDofsPerFace_[i] << "   ";
    }
    std::cout << std::endl;
    numLocalNodes_ = static_cast<int>( (intrepidBasis[0])->getBaseCellTopology().getVertexCount() );
    numLocalEdges_ = static_cast<int>( (intrepidBasis[0])->getBaseCellTopology().getEdgeCount() );
    numLocalFaces_ = 1;
    getFieldPattern(fieldPattern_);
    for (int i=0; i<numBases_; ++i) {
      std::cout << "\nBasis " << i << ":   "; 
      for (int j=0; j<(intrepidBasis[i])->getCardinality(); ++j) {
        std::cout << fieldPattern_[i][j] << " ";
      }
    }

    getDofArrays(nodeDofs_, numNodeDofs_, edgeDofs_, numEdgeDofs_, faceDofs_, numFaceDofs_, numDofs_);

    std::cout << numNodeDofs_ << nodeDofs_;
    std::cout << numEdgeDofs_ << edgeDofs_;
    std::cout << numFaceDofs_ << faceDofs_;
    std::cout << numDofs_;

    getCellDofs(cellDofs_, nodeDofs_, edgeDofs_, faceDofs_, meshCellToNodeMap_, meshCellToEdgeMap_);

    std::cout << cellDofs_;

    Intrepid::FieldContainer<int> fdofs;
    getFieldDofs(fdofs, 0);
    std::cout << fdofs;
    getFieldDofs(fdofs, 1);
    std::cout << fdofs;
    getFieldDofs(fdofs, 2);
    std::cout << fdofs;

  }


  void getDofArrays(Intrepid::FieldContainer<int> &nodeDofs, int &numNodeDofs,
                    Intrepid::FieldContainer<int> &edgeDofs, int &numEdgeDofs,
                    Intrepid::FieldContainer<int> &faceDofs, int &numFaceDofs,
                    int &numDofs) {

    nodeDofs.resize(numNodes_, numLocalNodeDofs_);  nodeDofs.initialize(-1);
    edgeDofs.resize(numEdges_, numLocalEdgeDofs_);  edgeDofs.initialize(-1);
    faceDofs.resize(numCells_, numLocalFaceDofs_);  faceDofs.initialize(-1);

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
    numNodeDofs = dofCt+1;

    // count edge dofs
    for (int i=0; i<numEdges_; ++i) {
      int locEdgeCt = -1;
      for (int j=0; j<numBases_; ++j) {
        for (int k=0; k<numDofsPerEdge_[j]; ++k) {
          edgeDofs(i, ++locEdgeCt) = ++dofCt;
        }
      }
    }
    numEdgeDofs = dofCt+1-numNodeDofs;

    // count face dofs
    for (int i=0; i<numCells_; ++i) {
      int locFaceCt = -1;
      for (int j=0; j<numBases_; ++j) {
        for (int k=0; k<numDofsPerFace_[j]; ++k) {
          faceDofs(i, ++locFaceCt) = ++dofCt;
        }
      }
    }
    numFaceDofs = dofCt+1-numNodeDofs-numEdgeDofs;

    numDofs = numNodeDofs+numEdgeDofs+numFaceDofs;

  } // getDofArrays


  void getCellDofs(Intrepid::FieldContainer<int> &cdofs,
                   const Intrepid::FieldContainer<int> &nodeDofs,
                   const Intrepid::FieldContainer<int> &edgeDofs,
                   const Intrepid::FieldContainer<int> &faceDofs,
                   const Intrepid::FieldContainer<int> &ctn,
                   const Intrepid::FieldContainer<int> &cte) {

    cdofs.resize(numCells_, numLocalDofs_);

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


  void getFieldPattern(std::vector<std::vector<int> > &fieldPattern) {
    
    fieldPattern.resize(numBases_);

    int dofCt = -1;

    // count node dofs
    for (int i=0; i<numLocalNodes_; ++i) {
      for (int j=0; j<numBases_; ++j) {
        for (int k=0; k<numDofsPerNode_[j]; ++k) {
          fieldPattern[j].push_back(++dofCt);
        }
      }
    }

    // count edge dofs
    for (int i=0; i<numLocalEdges_; ++i) {
      for (int j=0; j<numBases_; ++j) {
        for (int k=0; k<numDofsPerEdge_[j]; ++k) {
          fieldPattern[j].push_back(++dofCt);
        }
      }
    }

    // count face dofs
    for (int i=0; i<numLocalFaces_; ++i) {
      for (int j=0; j<numBases_; ++j) {
        for (int k=0; k<numDofsPerFace_[j]; ++k) {
          fieldPattern[j].push_back(++dofCt);
        }
      }
    }

  } // getFieldPattern


  void getFieldDofs(Intrepid::FieldContainer<int> & fdofs, const int & fieldNumber) {
    int basisCard = intrepidBasis_[fieldNumber]->getCardinality();
    fdofs.resize(numCells_, basisCard);
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<basisCard; ++j) {
        fdofs(i,j) = cellDofs_(i, fieldPattern_[fieldNumber][j]);
      }
    }  
  }

/*

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
*/

};

#endif
