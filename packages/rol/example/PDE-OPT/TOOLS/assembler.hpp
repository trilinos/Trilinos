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

/*! \file  data.hpp
    \brief Generates and manages data for the Poisson example, including
           all mesh and discretization data, matrices, etc.
*/

#ifndef ROL_PDEOPT_ASSEMBLER_H
#define ROL_PDEOPT_ASSEMBLER_H

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "MatrixMarket_Tpetra.hpp"

#include "Amesos2.hpp"

#include "fe.hpp"
#include "pde.hpp"
#include "qoi.hpp"
#include "dofmanager.hpp"
#include "meshmanager.hpp"
#include "boundarycells.hpp"

//// Global Timers.
//#ifdef ROL_TIMERS
//Teuchos::RCP<Time> FactorTime_example_PDEOPT_TOOLS_PDEFEM_GLOB = TimeMonitor::getNewCounter("ROL: Factorization Time in PDEFEM");
//Teuchos::RCP<Time> LUSubstitutionTime_example_PDEOPT_TOOLS_PDEFEM_GLOB = TimeMonitor::getNewCounter("ROL: LU Substitution Time in PDEFEM");
//Teuchos::RCP<Time> SolverUpdateTime_example_PDEOPT_TOOLS_PDEFEM_GLOB = TimeMonitor::getNewCounter("ROL: Solver Update Time in PDEFEM");
//Teuchos::RCP<Time> LocalAssemblyTime_example_PDEOPT_TOOLS_PDEFEM_GLOB = TimeMonitor::getNewCounter("ROL: Local Assembly Time in PDEFEM");
//Teuchos::RCP<Time> ConstraintDerivativeTime_example_PDEOPT_TOOLS_PDEFEM_GLOB = TimeMonitor::getNewCounter("ROL: Constraint Derivative Application Time in PDEFEM");
//#endif

template<class Real>
class Assembler {
private:
  // Timers
//  Teuchos::RCP<Teuchos::Time::Time> timerSolverFactorization_;
//  Teuchos::RCP<Teuchos::Time::Time> timerSolverSubstitution_;
//  Teuchos::RCP<Teuchos::Time::Time> timerAssemblyNonlinear_;
//  Teuchos::RCP<Teuchos::Time::Time> timerSolverUpdate_;

  // Set in Constructor
  bool verbose_;
  bool isJ1Transposed_, isJ2Transposed_;
  bool isSolverConstructed_, isSolverTransConstructed_;

  // Set in SetCommunicator
  Teuchos::RCP<const Teuchos::Comm<int> > comm_;
  int myRank_, numProcs_;

  // Set in SetBasis
  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;

  // Set in SetDiscretization
  Teuchos::RCP<MeshManager<Real> > meshMgr_;
  Teuchos::RCP<DofManager<Real> >  dofMgr_;

  // Set in SetParallelStructure
  int numCells_;
  Teuchos::Array<int> myCellIds_;
  Teuchos::Array<int> cellOffsets_;
  Teuchos::RCP<const Tpetra::Map<> > myOverlapStateMap_;
  Teuchos::RCP<const Tpetra::Map<> > myUniqueStateMap_;
  Teuchos::RCP<const Tpetra::Map<> > myOverlapControlMap_;
  Teuchos::RCP<const Tpetra::Map<> > myUniqueControlMap_;
  Teuchos::RCP<const Tpetra::Map<> > myOverlapResidualMap_;
  Teuchos::RCP<const Tpetra::Map<> > myUniqueResidualMap_;
  Teuchos::RCP<Tpetra::CrsGraph<> >  matJ1Graph_;
  Teuchos::RCP<Tpetra::CrsGraph<> >  matJ2Graph_;
  Teuchos::RCP<Tpetra::CrsGraph<> >  matH11Graph_;
  Teuchos::RCP<Tpetra::CrsGraph<> >  matH12Graph_;
  Teuchos::RCP<Tpetra::CrsGraph<> >  matH21Graph_;
  Teuchos::RCP<Tpetra::CrsGraph<> >  matH22Graph_;

  // Set in SetCellNodes
  Teuchos::RCP<Intrepid::FieldContainer<Real> > volCellNodes_;
  Teuchos::RCP<std::vector<std::vector<Intrepid::FieldContainer<int> > > >  bdryCellIds_;
  //std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<int> > > >  bdryCellLocIds_;
  std::vector<std::vector<std::vector<int> > >  bdryCellLocIds_;
  std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;

  // Finite element vectors and matrices
  Teuchos::RCP<Tpetra::CrsMatrix<> >    matJ1_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >    matJ1_trans_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >    matJ2_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >    matJ2_trans_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >    matH11_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >    matH12_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >    matH21_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >    matH22_;
  Teuchos::RCP<Tpetra::MultiVector<> >  vecR_;
  Teuchos::RCP<Tpetra::MultiVector<> >  vecR_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> >  vecG1_;
  Teuchos::RCP<Tpetra::MultiVector<> >  vecG1_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> >  vecG2_;
  Teuchos::RCP<Tpetra::MultiVector<> >  vecG2_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> >  vecH11_;
  Teuchos::RCP<Tpetra::MultiVector<> >  vecH11_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> >  vecH12_;
  Teuchos::RCP<Tpetra::MultiVector<> >  vecH12_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> >  vecH21_;
  Teuchos::RCP<Tpetra::MultiVector<> >  vecH21_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> >  vecH22_;
  Teuchos::RCP<Tpetra::MultiVector<> >  vecH22_overlap_;

  // Linear solvers for Jacobian and adjoint Jacobian
  Teuchos::RCP<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > solver_;
  Teuchos::RCP<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > solver_trans_;

private:

  void setCommunicator(const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                       Teuchos::ParameterList &parlist,
                       std::ostream &outStream = std::cout) {
    comm_ = comm;
    // Get number of processors and my rank
    myRank_    = comm->getRank();
    numProcs_  = comm->getSize();
    // Parse parameter list
    verbose_ = parlist.sublist("PDE FEM").get("Verbose Output",false);
    if (verbose_ && myRank_==0 ) {
      outStream << "Initialized communicator. " << std::endl;
    }
    if (verbose_ && myRank_==0 ) {
      outStream << "Total number of processors: " << numProcs_ << std::endl;
    }
  }

  void setBasis(const std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > &basisPtrs,
              Teuchos::ParameterList &parlist,
              std::ostream &outStream = std::cout) {
    basisPtrs_ = basisPtrs;
    if (verbose_ && myRank_==0) {
      outStream << "Initialized PDE." << std::endl;
    }
  }

  void setDiscretization(Teuchos::ParameterList &parlist,
                         const Teuchos::RCP<MeshManager<Real> > &meshMgr = Teuchos::null,
                         std::ostream &outStream = std::cout) {
    if ( meshMgr != Teuchos::null ) {
      // Use MeshManager object if supplied
      meshMgr_ = meshMgr;
    }
    else {
      // Otherwise construct MeshManager objective from parameter list
    }
    dofMgr_ = Teuchos::rcp(new DofManager<Real>(meshMgr_,basisPtrs_));
    if (verbose_ && myRank_==0) {
      outStream << "Initialized discretization (MeshManager and DofManager)." << std::endl;
    }
  }

  void setParallelStructure(Teuchos::ParameterList &parlist,
                            std::ostream &outStream = std::cout) {
    int cellSplit = parlist.sublist("Geometry").get<int>("Partition type");
    /****************************************************/
    /*** Build parallel communication infrastructure. ***/
    /****************************************************/
    // Partition the cells in the mesh.  We use a basic quasi-equinumerous partitioning,
    // where the remainder, if any, is assigned to the last processor.
    Teuchos::Array<int> myGlobalIds;
    cellOffsets_.assign(numProcs_, 0);
    int totalNumCells = meshMgr_->getNumCells();
    int cellsPerProc  = totalNumCells / numProcs_;
    numCells_         = cellsPerProc;
    switch(cellSplit) {
      case 0:
        if (myRank_ == 0) {
          // remainder in the first
          numCells_ += totalNumCells % numProcs_;
        }
        for (int i=1; i<numProcs_; ++i) {
          cellOffsets_[i] = cellOffsets_[i-1] + cellsPerProc
                            + (static_cast<int>(i==1))*(totalNumCells % numProcs_);
        }
        break;
      case 1:
        if (myRank_ == numProcs_-1) {
          // remainder in the last
          numCells_ += totalNumCells % numProcs_;
        }
        for (int i=1; i<numProcs_; ++i) {
          cellOffsets_[i] = cellOffsets_[i-1] + cellsPerProc;
        }
        break;
      case 2:
        if (myRank_ < (totalNumCells%numProcs_)) {
          // spread remainder, starting from the first
          numCells_++;
        }
        for (int i=1; i<numProcs_; ++i) {
          cellOffsets_[i] = cellOffsets_[i-1] + cellsPerProc
                            + (static_cast<int>(i-1<(totalNumCells%numProcs_)));
        }
        break;
    }

    Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
    int numLocalDofs = cellDofs.dimension(1);
    if (verbose_) {
      outStream << "Cell offsets across processors: " << cellOffsets_
                << std::endl;
    }
    for (int i=0; i<numCells_; ++i) {
      myCellIds_.push_back(cellOffsets_[myRank_]+i);
      for (int j=0; j<numLocalDofs; ++j) {
        myGlobalIds.push_back( cellDofs(cellOffsets_[myRank_]+i,j) );
      }
    }
    std::sort(myGlobalIds.begin(), myGlobalIds.end());
    myGlobalIds.erase( std::unique(myGlobalIds.begin(),myGlobalIds.end()),myGlobalIds.end() );

    // Build maps.
    myOverlapStateMap_ = Teuchos::rcp(new Tpetra::Map<>(
                         Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                         myGlobalIds, 0, comm_));
    //std::cout << std::endl << myOverlapMap_->getNodeElementList()<<std::endl;
    /** One can also use the non-member function:
        myOverlapMap_ = Tpetra::createNonContigMap<int,int>(myGlobalIds_, comm_);
        to build the overlap map.
    **/
    myUniqueStateMap_ = Tpetra::createOneToOne<int,int>(myOverlapStateMap_);
    //std::cout << std::endl << myUniqueMap_->getNodeElementList() << std::endl;
    myOverlapControlMap_  = myOverlapStateMap_;
    myUniqueControlMap_   = myUniqueStateMap_;
    myOverlapResidualMap_ = myOverlapStateMap_;
    myUniqueResidualMap_  = myUniqueStateMap_;
//    myCellMap_ = Teuchos::rcp(new Tpetra::Map<>(
//                 Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
//                 myCellIds_, 0, comm_));

    /****************************************/
    /*** Assemble global graph structure. ***/
    /****************************************/
    matJ1Graph_ = Teuchos::rcp(new Tpetra::CrsGraph<>(myUniqueStateMap_, 0));
    Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matJ1Graph_->insertGlobalIndices(cellDofs(myCellIds_[i],j),
          cellDofsArrayRCP(myCellIds_[i]*numLocalDofs, numLocalDofs));
      }
    }
    matJ1Graph_->fillComplete();
    matJ2Graph_ = matJ1Graph_;
    matH11Graph_ = matJ1Graph_;
    matH12Graph_ = matJ1Graph_;
    matH21Graph_ = matJ2Graph_;
    matH22Graph_ = matJ2Graph_;
    // Initialize residual vectors
    vecR_            = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueResidualMap_, 1, true));
    vecR_overlap_    = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapResidualMap_, 1, true));
    // Initialize state gradient vectors
    vecG1_           = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueStateMap_, 1, true));
    vecG1_overlap_   = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapStateMap_, 1, true));
    // Initialize control gradient vectors
    vecG2_           = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueControlMap_, 1, true));
    vecG2_overlap_   = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapControlMap_, 1, true));
    // Initialize state-state HessVec vectors
    vecH11_          = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueStateMap_, 1, true));
    vecH11_overlap_  = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapStateMap_, 1, true));
    // Initialize state-control HessVec vectors
    vecH12_          = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueStateMap_, 1, true));
    vecH12_overlap_  = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapStateMap_, 1, true));
    // Initialize control-state HessVec vectors
    vecH12_          = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueControlMap_, 1, true));
    vecH12_overlap_  = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapControlMap_, 1, true));
    // Initialize control-control HessVec vectors
    vecH22_          = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueControlMap_, 1, true));
    vecH22_overlap_  = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapControlMap_, 1, true));
    // Initialize Jacobian matrices
    matJ1_ = Teuchos::rcp(new Tpetra::CrsMatrix<>(matJ1Graph_));
    matJ2_ = Teuchos::rcp(new Tpetra::CrsMatrix<>(matJ2Graph_));
    // Initialize Hessian matrices
    matH11_ = Teuchos::rcp(new Tpetra::CrsMatrix<>(matH11Graph_));
    matH12_ = Teuchos::rcp(new Tpetra::CrsMatrix<>(matH12Graph_));
    matH21_ = Teuchos::rcp(new Tpetra::CrsMatrix<>(matH21Graph_));
    matH22_ = Teuchos::rcp(new Tpetra::CrsMatrix<>(matH22Graph_));

    if (verbose_ && myRank_==0) {
      outStream << "Initialized parallel structures." << std::endl;
    }
  }

  void setCellNodes(void) {
    // Build volume cell nodes
    shards::CellTopology cellType = basisPtrs_[0]->getBaseCellTopology();
    int spaceDim = cellType.getDimension();
    int numNodesPerCell = cellType.getNodeCount();
    volCellNodes_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, numNodesPerCell, spaceDim));
    Intrepid::FieldContainer<Real> &nodes = *meshMgr_->getNodes();
    Intrepid::FieldContainer<int>  &ctn   = *meshMgr_->getCellToNodeMap();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numNodesPerCell; ++j) {
        for (int k=0; k<spaceDim; ++k) {
          (*volCellNodes_)(i, j, k) = nodes(ctn(myCellIds_[i],j), k);
        }
      }
    }
    // Build boundary cell nodes
    bdryCellIds_    = meshMgr_->getSideSets();
    int numSideSets = bdryCellIds_->size();
    if (numSideSets > 0) {
      bdryCellNodes_.resize(numSideSets);
      bdryCellLocIds_.resize(numSideSets);
      for (int i=0; i<numSideSets; ++i) {
        int numLocSides = (*bdryCellIds_)[i].size();
        bdryCellNodes_[i].resize(numLocSides);
        bdryCellLocIds_[i].resize(numLocSides);
        for (int j=0; j<numLocSides; ++j) {
          int numCellsSide = (*bdryCellIds_)[i][j].dimension(0);
          for (int k=0; k<numCellsSide; ++k) {
            int idx = mapGlobalToLocalCellId((*bdryCellIds_)[i][j](k));
            if (idx > -1) {
              bdryCellLocIds_[i][j].push_back(idx);
              //if (myRank_==1) {std::cout << "\nrank " << myRank_ << "   bcid " << i << "  " << j << "  " << k << "  " << myCellIds_[idx] << "  " << idx;}
            }
          }
          int myNumCellsSide = bdryCellLocIds_[i][j].size();
          if (myNumCellsSide > 0) {
            bdryCellNodes_[i][j] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(myNumCellsSide, numNodesPerCell, spaceDim));
          }
          for (int k=0; k<myNumCellsSide; ++k) {
            for (int l=0; l<numNodesPerCell; ++l) {
              for (int m=0; m<spaceDim; ++m) {
                (*bdryCellNodes_[i][j])(k, l, m) = nodes(ctn(myCellIds_[bdryCellLocIds_[i][j][k]],l), m);
              }
            }
          }
          //if ((myRank_==1) && (myNumCellsSide > 0)) {std::cout << (*bdryCellNodes_[i][j]);}
        }
      }
    }
    bdryCellNodes_.resize(numSideSets);
  }

  int mapGlobalToLocalCellId(const int & gid) {
    int minId = cellOffsets_[myRank_];
    int maxId = cellOffsets_[myRank_]+numCells_-1;
    if ((gid >= minId) && (gid <= maxId)) {
      return (gid - cellOffsets_[myRank_]);
    }
    else {
      return -1;
    }
  }

  void getCoeffFromStateVector(Teuchos::RCP<Intrepid::FieldContainer<Real> > &xcoeff,
                               const Tpetra::MultiVector<> &x) const {
    // Perform import onto myOverlapMap
    Teuchos::RCP<Tpetra::MultiVector<> > xshared =
      Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapStateMap_, 1, true));
    Tpetra::Import<> importer(myUniqueStateMap_, myOverlapStateMap_);
    xshared->doImport(x,importer,Tpetra::REPLACE);
    // Populate xcoeff
    Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
    int lfs = dofMgr_->getLocalFieldSize(0);
    xcoeff = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs));
    Teuchos::ArrayRCP<const Real> xdata = xshared->get1dView();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<lfs; ++j) {
        (*xcoeff)(i, j) = xdata[xshared->getMap()->getLocalElement(cellDofs(myCellIds_[i],j))];
      }
    }
  }

  void getCoeffFromControlVector(Teuchos::RCP<Intrepid::FieldContainer<Real> > &xcoeff,
                                 const Tpetra::MultiVector<> &x) const {
    // Perform import onto myOverlapMap
    Teuchos::RCP<Tpetra::MultiVector<> > xshared =
      Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapControlMap_, 1, true));
    Tpetra::Import<> importer(myUniqueControlMap_, myOverlapControlMap_);
    xshared->doImport(x,importer,Tpetra::REPLACE);
    // Populate xcoeff
    Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
    int lfs = dofMgr_->getLocalFieldSize(0);
    xcoeff = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs));
    Teuchos::ArrayRCP<const Real> xdata = xshared->get1dView();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<lfs; ++j) {
        (*xcoeff)(i, j) = xdata[xshared->getMap()->getLocalElement(cellDofs(myCellIds_[i],j))];
      }
    }
  }

  void constructSolver(const bool transpose = false) {
    // Construct solver using Amesos2 factory.
//    #ifdef ROL_TIMERS
//    Teuchos::TimeMonitor LocalTimer(*FactorTime_example_PDEOPT_TOOLS_PDEFEM_GLOB);
//    #endif
    if ( transpose ) {
      try{
        solver_trans_
          = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("KLU2", matJ1_trans_);
      } catch (std::invalid_argument e) {
        std::cout << e.what() << std::endl;
      }
      solver_trans_->numericFactorization();
    }
    else {
      try {
        solver_
          = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("KLU2", matJ1_);
      }
      catch (std::invalid_argument e) {
        std::cout << e.what() << std::endl;
      }
      solver_->numericFactorization();
    }
  }

public:
  // destructor
  virtual ~Assembler() {}

  // Constuctor: Discretization set from ParameterList
  Assembler(const std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > &basisPtrs,
          const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
          Teuchos::ParameterList &parlist,
          std::ostream &outStream = std::cout)
    : isJ1Transposed_(false),
      isJ2Transposed_(false),
      isSolverConstructed_(false),
      isSolverTransConstructed_(false) {
    setCommunicator(comm,parlist,outStream);
    setBasis(basisPtrs,parlist,outStream);
    setDiscretization(parlist,Teuchos::null,outStream);
    setParallelStructure(parlist,outStream);
    setCellNodes();
  }

  // Constructor: Discretization set from MeshManager input
  Assembler(const std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > &basisPtrs,
          const Teuchos::RCP<MeshManager<Real> > &meshMgr,
          const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
          Teuchos::ParameterList &parlist,
          std::ostream &outStream = std::cout)
    : isJ1Transposed_(false),
      isJ2Transposed_(false),
      isSolverConstructed_(false),
      isSolverTransConstructed_(false) {
    setCommunicator(comm,parlist,outStream);
    setBasis(basisPtrs,parlist,outStream);
    setDiscretization(parlist,meshMgr,outStream);
    setParallelStructure(parlist,outStream);
    setCellNodes();
  }

  void setCellNodes(PDE<Real> &pde) const {
    // Set PDE cell nodes
    pde.setCellNodes(volCellNodes_, bdryCellNodes_, bdryCellLocIds_);
  }

  void assemblePDEResidual(const Tpetra::MultiVector<> &u,
                           const Tpetra::MultiVector<> &z,
                           PDE<Real> &pde) {
    const Real zero(0);
    vecR_->scale(zero);
    vecR_overlap_->scale(zero);
    Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
    int numLocalDofs = cellDofs.dimension(1);
    // Initialize res
    Teuchos::RCP<Intrepid::FieldContainer<Real> > res;
    // Get u_coeff from u and z_coeff from z
    Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff, z_coeff;
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    // Compute PDE residual
    pde.residual(res,u_coeff,z_coeff);
    // assembly on the overlap map
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        vecR_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                          0,
                                          (*res)[i*numLocalDofs+j]);
      }
    }
    // change map
    Tpetra::Export<> exporter(vecR_overlap_->getMap(), vecR_->getMap()); // redistribution
    vecR_->doExport(*vecR_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
  }

  void assemblePDEJacobian1(const Tpetra::MultiVector<> &u,
                            const Tpetra::MultiVector<> &z,
                            PDE<Real> &pde) {
    const Real zero(0);
    matJ1_->resumeFill();
    matJ1_->setAllToScalar(zero);
    Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
    int numLocalDofs = cellDofs.dimension(1);
    int numLocalMatEntries = numLocalDofs * numLocalDofs;
    Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
    // Initialize jac
    Teuchos::RCP<Intrepid::FieldContainer<Real> > jac;
    // Get u_coeff from u and z_coeff from z
    Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff, z_coeff;
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    // Compute PDE Jacobian
    pde.Jacobian_1(jac,u_coeff,z_coeff);
    Teuchos::ArrayRCP<const Real> jacArrayRCP = jac->getData();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matJ1_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                    cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                    jacArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
      }
    }
    matJ1_->fillComplete();
    isJ1Transposed_ = false;
    isSolverConstructed_ = false;
    isSolverTransConstructed_ = false;
  }

  void assemblePDEJacobian2(const Tpetra::MultiVector<> &u,
                            const Tpetra::MultiVector<> &z,
                            PDE<Real> &pde) { 
    const Real zero(0);
    matJ2_->resumeFill();
    matJ2_->setAllToScalar(zero);
    Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
    int numLocalDofs = cellDofs.dimension(1);
    int numLocalMatEntries = numLocalDofs * numLocalDofs;
    Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
    Teuchos::RCP<Intrepid::FieldContainer<Real> > jac, u_coeff, z_coeff;
    // Get u_coeff from u and z_coeff from z
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    // Compute PDE Jacobian
    pde.Jacobian_2(jac,u_coeff,z_coeff);
    Teuchos::ArrayRCP<const Real> jacArrayRCP = jac->getData();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matJ2_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                    cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                    jacArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
      }
    }
    matJ2_->fillComplete();
    // Create matrix transposes.
    Tpetra::RowMatrixTransposer<> transposerJ2(matJ2_);
    matJ2_trans_ = transposerJ2.createTranspose();
    isJ2Transposed_ = false;
  }

  void assemblePDEHessian11(const Tpetra::MultiVector<> &u,
                            const Tpetra::MultiVector<> &z,
                            const Tpetra::MultiVector<> &l,
                            PDE<Real> &pde) { 
    try {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > hess, u_coeff, z_coeff, l_coeff;
      // Get u_coeff from u and z_coeff from z
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      getCoeffFromStateVector(l_coeff,l);
      // Compute PDE Hessian
      pde.Hessian_11(hess,u_coeff,z_coeff,l_coeff);
      // Zero Hessian
      const Real zero(0);
      matH11_->resumeFill();
      matH11_->setAllToScalar(zero);
      // Assemble PDE Hessian
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> hessArrayRCP = hess->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          matH11_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                       cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                       hessArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      matH11_->fillComplete();
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assemblePDEHessian11): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian11): Hessian not implemented.");
    }
  }

  void assemblePDEHessian12(const Tpetra::MultiVector<> &u,
                            const Tpetra::MultiVector<> &z,
                            const Tpetra::MultiVector<> &l,
                            PDE<Real> &pde) { 
    try {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > hess, u_coeff, z_coeff, l_coeff;
      // Get u_coeff from u and z_coeff from z
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      getCoeffFromStateVector(l_coeff,l);
      // Compute PDE Hessian
      pde.Hessian_12(hess,u_coeff,z_coeff,l_coeff);
      // Zero Hessian
      const Real zero(0);
      matH12_->resumeFill();
      matH12_->setAllToScalar(zero);
      // Assemble PDE Hessian
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> hessArrayRCP = hess->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          matH12_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                       cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                       hessArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      matH12_->fillComplete();
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assemblePDEHessian12): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian12): Hessian not implemented.");
    }
  }

  void assemblePDEHessian21(const Tpetra::MultiVector<> &u,
                            const Tpetra::MultiVector<> &z,
                            const Tpetra::MultiVector<> &l,
                            PDE<Real> &pde) { 
    try {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > hess, u_coeff, z_coeff, l_coeff;
      // Get u_coeff from u and z_coeff from z
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      getCoeffFromStateVector(l_coeff,l);
      // Compute PDE Hessian
      pde.Hessian_21(hess,u_coeff,z_coeff,l_coeff);
      // Zero Hessian
      const Real zero(0);
      matH21_->resumeFill();
      matH21_->setAllToScalar(zero);
      // Assemble PDE Hessian
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> hessArrayRCP = hess->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          matH21_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                       cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                       hessArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      matH21_->fillComplete();
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assemblePDEHessian21): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian21): Hessian not implemented.");
    }
  }

  void assemblePDEHessian22(const Tpetra::MultiVector<> &u,
                            const Tpetra::MultiVector<> &z,
                            const Tpetra::MultiVector<> &l,
                            PDE<Real> &pde) { 
    try {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > hess, u_coeff, z_coeff, l_coeff;
      // Get u_coeff from u and z_coeff from z
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      getCoeffFromStateVector(l_coeff,l);
      // Compute PDE Hessian
      pde.Hessian_22(hess,u_coeff,z_coeff,l_coeff);
      // Zero Hessian
      const Real zero(0);
      matH22_->resumeFill();
      matH22_->setAllToScalar(zero);
      // Assemble PDE Hessian
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> hessArrayRCP = hess->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          matH22_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                       cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                       hessArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      matH22_->fillComplete();
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assemblePDEHessian22): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian22): Hessian not implemented.");
    }
  }

  Real assembleQoIValue(const Tpetra::MultiVector<> &u,
                        const Tpetra::MultiVector<> &z,
                        QoI<Real> &qoi) {
    // Integrate obj object
    Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff, z_coeff, locVal;
    // Get u_coeff from u and z_coeff from z
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    // Get OBJ_CELL value
    qoi.value(locVal,u_coeff,z_coeff);
    // Assembly
    Real val(0), myval(0);
    for (int i=0; i<numCells_; ++i) {
      myval += (*locVal)(i);
    }
    Teuchos::reduceAll<int,Real>(*comm_,Teuchos::REDUCE_SUM,1,&myval,&val);
    return val;
  }

  void assembleQoIGradient1(const Tpetra::MultiVector<> &u,
                            const Tpetra::MultiVector<> &z,
                            QoI<Real> &qoi) {
    const Real zero(0);
    vecG1_->scale(zero);
    vecG1_overlap_->scale(zero);
    Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
    int numLocalDofs = cellDofs.dimension(1);
    // Get u_coeff from u and z_coeff from z
    Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff, z_coeff, locGrad;
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    // Compute local gradient
    qoi.gradient_1(locGrad,u_coeff,z_coeff);
    // Assembly in to the overlap gradient
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        vecG1_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                           0,
                                           (*locGrad)[i*numLocalDofs+j]);
      }
    }
    // Change to local map
    Tpetra::Export<> exporter(vecG1_overlap_->getMap(), vecG1_->getMap()); // redistribution
    vecG1_->doExport(*vecG1_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
  }

  void assembleQoIGradient2(const Tpetra::MultiVector<> &u,
                            const Tpetra::MultiVector<> &z,
                            QoI<Real> &qoi) {
    const Real zero(0);
    vecG2_->scale(zero);
    vecG2_overlap_->scale(zero);
    Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
    int numLocalDofs = cellDofs.dimension(1);
    // Get u_coeff from u and z_coeff from z
    Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff, z_coeff, locGrad;
    getCoeffFromStateVector(u_coeff,u);
    getCoeffFromControlVector(z_coeff,z);
    // Compute local gradient
    qoi.gradient_2(locGrad,u_coeff,z_coeff);
    // Assembly in to the overlap gradient
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        vecG2_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                           0,
                                           (*locGrad)[i*numLocalDofs+j]);
      }
    }
    // Change to local map
    Tpetra::Export<> exporter(vecG2_overlap_->getMap(), vecG2_->getMap()); // redistribution
    vecG2_->doExport(*vecG2_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
  }

  void assembleQoIHessVec11(const Tpetra::MultiVector<> &v,
                            const Tpetra::MultiVector<> &u,
                            const Tpetra::MultiVector<> &z,
                            QoI<Real> &qoi) {
    try {
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > v_coeff, u_coeff, z_coeff, locHess;
      getCoeffFromStateVector(v_coeff,v);
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi.HessVec_11(locHess, v_coeff, u_coeff,z_coeff);
      // Assembly in to the overlap gradient
      const Real zero(0);
      vecH11_->scale(zero);
      vecH11_overlap_->scale(zero);
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          vecH11_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                              0,
                                              (*locHess)[i*numLocalDofs+j]);
        }
      }
      // Change to local map
      Tpetra::Export<> exporter(vecH11_overlap_->getMap(), vecH11_->getMap()); // redistribution
      vecH11_->doExport(*vecH11_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assembleQoIHessVec11): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec11): Hessian not implemented.");
    }
  }

  void assembleQoIHessVec12(const Tpetra::MultiVector<> &v,
                            const Tpetra::MultiVector<> &u,
                            const Tpetra::MultiVector<> &z,
                            QoI<Real> &qoi) {
    try {
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > v_coeff, u_coeff, z_coeff, locHess;
      getCoeffFromControlVector(v_coeff,v);
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi.HessVec_12(locHess, v_coeff, u_coeff,z_coeff);
      // Assembly in to the overlap gradient
      const Real zero(0);
      vecH12_->scale(zero);
      vecH12_overlap_->scale(zero);
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          vecH12_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                              0,
                                              (*locHess)[i*numLocalDofs+j]);
        }
      }
      // Change to local map
      Tpetra::Export<> exporter(vecH12_overlap_->getMap(), vecH12_->getMap()); // redistribution
      vecH12_->doExport(*vecH12_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assembleQoIHessVec12): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec12): Hessian not implemented.");
    }
  }

  void assembleQoIHessVec21(const Tpetra::MultiVector<> &v,
                            const Tpetra::MultiVector<> &u,
                            const Tpetra::MultiVector<> &z,
                            QoI<Real> &qoi) {
    try {
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > v_coeff, u_coeff, z_coeff, locHess;
      getCoeffFromStateVector(v_coeff,v);
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi.HessVec_21(locHess, v_coeff, u_coeff,z_coeff);
      // Assembly in to the overlap gradient
      const Real zero(0);
      vecH21_->scale(zero);
      vecH21_overlap_->scale(zero);
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          vecH21_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                              0,
                                              (*locHess)[i*numLocalDofs+j]);
        }
      }
      // Change to local map
      Tpetra::Export<> exporter(vecH21_overlap_->getMap(), vecH21_->getMap()); // redistribution
      vecH21_->doExport(*vecH21_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assembleQoIHessVec21): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec21): Hessian not implemented.");
    }
  }

  void assembleQoIHessVec22(const Tpetra::MultiVector<> &v,
                            const Tpetra::MultiVector<> &u,
                            const Tpetra::MultiVector<> &z,
                            QoI<Real> &qoi) {
    try {
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > v_coeff, u_coeff, z_coeff, locHess;
      getCoeffFromControlVector(v_coeff,v);
      getCoeffFromStateVector(u_coeff,u);
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi.HessVec_22(locHess, v_coeff, u_coeff,z_coeff);
      // Assembly in to the overlap gradient
      const Real zero(0);
      vecH22_->scale(zero);
      vecH22_overlap_->scale(zero);
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          vecH22_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                              0,
                                              (*locHess)[i*numLocalDofs+j]);
        }
      }
      // Change to local map
      Tpetra::Export<> exporter(vecH22_overlap_->getMap(), vecH22_->getMap()); // redistribution
      vecH22_->doExport(*vecH22_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assembleQoIHessVec22): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec22): Hessian not implemented.");
    }
  }

  void linearPDEsolve(Teuchos::RCP<Tpetra::MultiVector<> > &u,
                      const Teuchos::RCP<const Tpetra::MultiVector<> > &r,
                      const bool transpose = false) {
    if ( transpose ) {
      if (!isJ1Transposed_) {
        // Create matrix transposes.
        Tpetra::RowMatrixTransposer<> transposerJ1(matJ1_);
        matJ1_trans_ = transposerJ1.createTranspose();
        isJ1Transposed_ = true;
      }
      if (!isSolverTransConstructed_) {
        constructSolver(true);
        isSolverTransConstructed_ = true;
      }
      solver_trans_->setX(u);
      solver_trans_->setB(r);
      solver_trans_->solve();
    }
    else {
      if (!isSolverConstructed_) {
        constructSolver(false);
        isSolverConstructed_ = true;
      }
      solver_->setX(u);
      solver_->setB(r);
      solver_->solve();
    }
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getPDEResidual(void) const {
    return vecR_;
  }

  Teuchos::RCP<Tpetra::CrsMatrix<> > getPDEJacobian1(const bool transpose = false) {
    if ( transpose ) {
      if (!isJ1Transposed_) {
        // Create matrix transposes.
        Tpetra::RowMatrixTransposer<> transposerJ1(matJ1_);
        matJ1_trans_ = transposerJ1.createTranspose();
        isJ1Transposed_ = true;
      }
      return matJ1_trans_;
    }
    else {
      return matJ1_;
    }
  }

  Teuchos::RCP<Tpetra::CrsMatrix<> > getPDEJacobian2(const bool transpose = false) {
    if ( transpose ) {
      if (!isJ2Transposed_) {
        // Create matrix transposes.
        Tpetra::RowMatrixTransposer<> transposerJ2(matJ2_);
        matJ2_trans_ = transposerJ2.createTranspose();
        isJ2Transposed_ = true;
      }
      return matJ2_trans_;
    }
    else {
      return matJ2_;
    }
  }

  Teuchos::RCP<Tpetra::CrsMatrix<> > getPDEHessian11(void) const {
    return matH11_;
  }

  Teuchos::RCP<Tpetra::CrsMatrix<> > getPDEHessian12(void) const {
    return matH12_;
  }

  Teuchos::RCP<Tpetra::CrsMatrix<> > getPDEHessian21(void) const {
    return matH21_;
  }

  Teuchos::RCP<Tpetra::CrsMatrix<> > getPDEHessian22(void) const {
    return matH22_;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getQoIGradient1(void) const {
    return vecG1_;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getQoIGradient2(void) const {
    return vecG2_;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getQoIHessVec11(void) const {
    return vecH11_;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getQoIHessVec12(void) const {
    return vecH12_;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getQoIHessVec21(void) const {
    return vecH21_;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getQoIHessVec22(void) const {
    return vecH22_;
  }

  void printMeshData(std::ostream &outStream) const {
    Teuchos::RCP<Intrepid::FieldContainer<Real> > nodesPtr = meshMgr_->getNodes();
    Teuchos::RCP<Intrepid::FieldContainer<int> >  cellToNodeMapPtr = meshMgr_->getCellToNodeMap();
    Intrepid::FieldContainer<Real>  &nodes = *nodesPtr;
    Intrepid::FieldContainer<int>   &cellToNodeMap = *cellToNodeMapPtr;
    outStream << std::endl;
    outStream << "Number of nodes = " << meshMgr_->getNumNodes() << std::endl;
    outStream << "Number of cells = " << meshMgr_->getNumCells() << std::endl;
    outStream << "Number of edges = " << meshMgr_->getNumEdges() << std::endl;
    // Print mesh to file.
    if ((myRank_ == 0)) {
      std::ofstream meshfile;
      meshfile.open("cell_to_node_quad.txt");
      for (int i=0; i<cellToNodeMap.dimension(0); ++i) {
        for (int j=0; j<cellToNodeMap.dimension(1); ++j) {
          meshfile << cellToNodeMap(i,j) << "  ";
        }
        meshfile << std::endl;
      }
      meshfile.close();
      
      meshfile.open("cell_to_node_tri.txt");
      for (int i=0; i<cellToNodeMap.dimension(0); ++i) {
        for (int j=0; j<3; ++j) {
          meshfile << cellToNodeMap(i,j) << "  ";
        }
        meshfile << std::endl;
        for (int j=2; j<5; ++j) {
          meshfile << cellToNodeMap(i,j%4) << "  ";
        }
        meshfile << std::endl;
      }
      meshfile.close();
     
      meshfile.open("nodes.txt");
      meshfile.precision(16);
      for (int i=0; i<nodes.dimension(0); ++i) {
        for (int j=0; j<nodes.dimension(1); ++j) {
          meshfile << std::scientific << nodes(i,j) << "  ";
        }
        meshfile << std::endl;
      }
      meshfile.close();
      /* This somewhat clunky output is for gnuplot.
      meshfile.open("mesh.txt");
      for (int i=0; i<cellToNodeMap.dimension(0); ++i) {
        meshfile << nodes(cellToNodeMap(i,0), 0) << "  " << nodes(cellToNodeMap(i,0), 1) << std::endl;
        meshfile << nodes(cellToNodeMap(i,1), 0) << "  " << nodes(cellToNodeMap(i,1), 1) << std::endl;
        meshfile << nodes(cellToNodeMap(i,2), 0) << "  " << nodes(cellToNodeMap(i,2), 1) << std::endl;
        meshfile << nodes(cellToNodeMap(i,3), 0) << "  " << nodes(cellToNodeMap(i,3), 1) << std::endl;
        meshfile << nodes(cellToNodeMap(i,0), 0) << "  " << nodes(cellToNodeMap(i,0), 1) << std::endl;
        meshfile << nodes(cellToNodeMap(i,1), 0) << "  " << nodes(cellToNodeMap(i,1), 1) << std::endl;
        meshfile << nodes(cellToNodeMap(i,2), 0) << "  " << nodes(cellToNodeMap(i,2), 1) << std::endl;
      }
      meshfile.close();
      */
    } //myRank 0 print
  } // prinf function end

  void outputTpetraData() const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> > matWriter;
    matWriter.writeSparseFile("jacobian1", matJ1_);
    matWriter.writeSparseFile("jacobian2", matJ2_);
    matWriter.writeDenseFile("residual", vecR_);
  }

  void outputTpetraVector(const Teuchos::RCP<const Tpetra::MultiVector<> > &vec,
                          const std::string &filename) const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> > vecWriter;
    vecWriter.writeDenseFile(filename, vec);
  }

  const Teuchos::RCP<const Tpetra::Map<> > getStateMap(void) const {
    return myUniqueStateMap_;
  }

  const Teuchos::RCP<const Tpetra::Map<> > getControlMap(void) const {
    return myUniqueControlMap_;
  }

  const Teuchos::RCP<const Tpetra::Map<> > getResidualMap(void) const {
    return myUniqueResidualMap_;
  }
 
  Teuchos::RCP<Tpetra::MultiVector<> > createStateVector(void) const {
    return Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueStateMap_, 1, true));
  }
 
  Teuchos::RCP<Tpetra::MultiVector<> > createControlVector(void) const {
    return Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueControlMap_, 1, true));
  }
 
  Teuchos::RCP<Tpetra::MultiVector<> > createResidualVector(void) const {
    return Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueResidualMap_, 1, true));
  }

}; // class Assembler

#endif
