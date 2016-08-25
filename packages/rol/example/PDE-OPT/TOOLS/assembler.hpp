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
  bool isSolverR1Constructed_, isSolverR2Constructed_;

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
  Teuchos::RCP<Tpetra::CrsGraph<> >  matR1Graph_;
  Teuchos::RCP<Tpetra::CrsGraph<> >  matR2Graph_;
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

  // Finite element vectors and matrices for PDE
  Teuchos::RCP<Tpetra::MultiVector<> >           pde_vecR_;
  Teuchos::RCP<Tpetra::MultiVector<> >           pde_vecR_overlap_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >             pde_matJ1_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >             pde_matJ1_trans_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >             pde_matJ2_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >             pde_matJ2_trans_;
  Teuchos::RCP<Tpetra::MultiVector<> >           pde_vecJ3_;
  Teuchos::RCP<Tpetra::MultiVector<> >           pde_vecJ3_overlap_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >             pde_matH11_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >             pde_matH12_;
  Teuchos::RCP<Tpetra::MultiVector<> >           pde_vecH13_;
  Teuchos::RCP<Tpetra::MultiVector<> >           pde_vecH13_overlap_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >             pde_matH21_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >             pde_matH22_;
  Teuchos::RCP<Tpetra::MultiVector<> >           pde_vecH23_;
  Teuchos::RCP<Tpetra::MultiVector<> >           pde_vecH23_overlap_;
  Teuchos::RCP<std::vector<std::vector<Real> > > pde_matH33_;

  // Finite element vectors and matrices for QoI
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecG1_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecG1_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecG2_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecG2_overlap_;
  Teuchos::RCP<std::vector<Real> >     qoi_vecG3_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH11_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH11_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH12_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH12_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH13_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH13_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH21_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH21_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH22_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH22_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH23_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH23_overlap_;
  Teuchos::RCP<std::vector<Real> >     qoi_vecH31_;
  Teuchos::RCP<std::vector<Real> >     qoi_vecH32_;
  Teuchos::RCP<std::vector<Real> >     qoi_vecH33_;

  // Finite element Riesz matrices
  Teuchos::RCP<Tpetra::CrsMatrix<> > matR1_;
  Teuchos::RCP<Tpetra::CrsMatrix<> > matR2_;

  // Linear solvers for Jacobian and adjoint Jacobian
  Teuchos::RCP<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > solver_;
  Teuchos::RCP<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > solver_trans_;

  // Linear solvers for Riesz maps
  Teuchos::RCP<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > solverR1_;
  Teuchos::RCP<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > solverR2_;

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

  void setBasis(
         const std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > &basisPtrs,
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
    matJ2Graph_  = matJ1Graph_;
    matR1Graph_  = matJ1Graph_;
    matR2Graph_  = matJ2Graph_;
    matH11Graph_ = matJ1Graph_;
    matH12Graph_ = matJ1Graph_;
    matH21Graph_ = matJ2Graph_;
    matH22Graph_ = matJ2Graph_;

    if (verbose_ && myRank_==0) {
      outStream << "Initialized parallel structures." << std::endl;
    }
  }

  void setCellNodes(std::ostream &outStream = std::cout) {
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
    bdryCellIds_    = meshMgr_->getSideSets(outStream, (verbose_ && myRank_==0));
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
                               const Teuchos::RCP<const Tpetra::MultiVector<> > &x) const {
    if ( x != Teuchos::null ) {
      // Perform import onto myOverlapMap
      Teuchos::RCP<Tpetra::MultiVector<> > xshared =
        Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapStateMap_, 1, true));
      Tpetra::Import<> importer(myUniqueStateMap_, myOverlapStateMap_);
      xshared->doImport(*x,importer,Tpetra::REPLACE);
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
    else {
      xcoeff = Teuchos::null;
    }
  }

  void getCoeffFromControlVector(Teuchos::RCP<Intrepid::FieldContainer<Real> > &xcoeff,
                                 const Teuchos::RCP<const Tpetra::MultiVector<> > &x) const {
    if ( x != Teuchos::null ) {
      // Perform import onto myOverlapMap
      Teuchos::RCP<Tpetra::MultiVector<> > xshared =
        Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapControlMap_, 1, true));
      Tpetra::Import<> importer(myUniqueControlMap_, myOverlapControlMap_);
      xshared->doImport(*x,importer,Tpetra::REPLACE);
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
    else {
      xcoeff = Teuchos::null;
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
          = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("KLU2", pde_matJ1_trans_);
      } catch (std::invalid_argument e) {
        std::cout << e.what() << std::endl;
      }
      solver_trans_->numericFactorization();
    }
    else {
      try {
        solver_
          = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("KLU2", pde_matJ1_);
      }
      catch (std::invalid_argument e) {
        std::cout << e.what() << std::endl;
      }
      solver_->numericFactorization();
    }
  }

  void constructSolverR1(void) {
    // Construct solver using Amesos2 factory.
//    #ifdef ROL_TIMERS
//    Teuchos::TimeMonitor LocalTimer(*FactorTime_example_PDEOPT_TOOLS_PDEFEM_GLOB);
//    #endif
    try {
      solverR1_
        = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("KLU2", matR1_);
    }
    catch (std::invalid_argument e) {
      std::cout << e.what() << std::endl;
    }
    solverR1_->numericFactorization();
  }

  void constructSolverR2(void) {
    // Construct solver using Amesos2 factory.
//    #ifdef ROL_TIMERS
//    Teuchos::TimeMonitor LocalTimer(*FactorTime_example_PDEOPT_TOOLS_PDEFEM_GLOB);
//    #endif
    try {
      solverR2_
        = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("KLU2", matR2_);
    }
    catch (std::invalid_argument e) {
      std::cout << e.what() << std::endl;
    }
    solverR2_->numericFactorization();
  }

public:
  // destructor
  virtual ~Assembler() {}

  // Constuctor: Discretization set from ParameterList
  Assembler(const std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > &basisPtrs,
          const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
          Teuchos::ParameterList &parlist,
          std::ostream &outStream = std::cout)
    : isJ1Transposed_(false),      isJ2Transposed_(false),
      isSolverConstructed_(false), isSolverTransConstructed_(false),
      pde_vecR_(Teuchos::null),    pde_vecR_overlap_(Teuchos::null),
      pde_matJ1_(Teuchos::null),   pde_matJ2_(Teuchos::null),
      pde_vecJ3_(Teuchos::null),   pde_vecJ3_overlap_(Teuchos::null),
      pde_matH11_(Teuchos::null),  pde_matH12_(Teuchos::null),
      pde_vecH13_(Teuchos::null),  pde_vecH13_overlap_(Teuchos::null),
      pde_matH21_(Teuchos::null),  pde_matH22_(Teuchos::null),
      pde_vecH23_(Teuchos::null),  pde_vecH23_overlap_(Teuchos::null),
      pde_matH33_(Teuchos::null),
      qoi_vecG1_(Teuchos::null),   qoi_vecG1_overlap_(Teuchos::null),
      qoi_vecG2_(Teuchos::null),   qoi_vecG2_overlap_(Teuchos::null),
      qoi_vecG3_(Teuchos::null),
      qoi_vecH11_(Teuchos::null),  qoi_vecH11_overlap_(Teuchos::null),
      qoi_vecH12_(Teuchos::null),  qoi_vecH12_overlap_(Teuchos::null),
      qoi_vecH13_(Teuchos::null),  qoi_vecH13_overlap_(Teuchos::null),
      qoi_vecH21_(Teuchos::null),  qoi_vecH21_overlap_(Teuchos::null),
      qoi_vecH22_(Teuchos::null),  qoi_vecH22_overlap_(Teuchos::null),
      qoi_vecH23_(Teuchos::null),  qoi_vecH23_overlap_(Teuchos::null),
      qoi_vecH31_(Teuchos::null),  qoi_vecH32_(Teuchos::null),
      qoi_vecH33_(Teuchos::null),
      matR1_(Teuchos::null),       matR2_(Teuchos::null) {
    setCommunicator(comm,parlist,outStream);
    setBasis(basisPtrs,parlist,outStream);
    setDiscretization(parlist,Teuchos::null,outStream);
    setParallelStructure(parlist,outStream);
    setCellNodes(outStream);
  }

  // Constructor: Discretization set from MeshManager input
  Assembler(const std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > &basisPtrs,
          const Teuchos::RCP<MeshManager<Real> > &meshMgr,
          const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
          Teuchos::ParameterList &parlist,
          std::ostream &outStream = std::cout)
    : isJ1Transposed_(false),      isJ2Transposed_(false),
      isSolverConstructed_(false), isSolverTransConstructed_(false),
      pde_vecR_(Teuchos::null),    pde_vecR_overlap_(Teuchos::null),
      pde_matJ1_(Teuchos::null),   pde_matJ2_(Teuchos::null),
      pde_vecJ3_(Teuchos::null),   pde_vecJ3_overlap_(Teuchos::null),
      pde_matH11_(Teuchos::null),  pde_matH12_(Teuchos::null),
      pde_vecH13_(Teuchos::null),  pde_vecH13_overlap_(Teuchos::null),
      pde_matH21_(Teuchos::null),  pde_matH22_(Teuchos::null),
      pde_vecH23_(Teuchos::null),  pde_vecH23_overlap_(Teuchos::null),
      pde_matH33_(Teuchos::null),
      qoi_vecG1_(Teuchos::null),   qoi_vecG1_overlap_(Teuchos::null),
      qoi_vecG2_(Teuchos::null),   qoi_vecG2_overlap_(Teuchos::null),
      qoi_vecG3_(Teuchos::null),
      qoi_vecH11_(Teuchos::null),  qoi_vecH11_overlap_(Teuchos::null),
      qoi_vecH12_(Teuchos::null),  qoi_vecH12_overlap_(Teuchos::null),
      qoi_vecH13_(Teuchos::null),  qoi_vecH13_overlap_(Teuchos::null),
      qoi_vecH21_(Teuchos::null),  qoi_vecH21_overlap_(Teuchos::null),
      qoi_vecH22_(Teuchos::null),  qoi_vecH22_overlap_(Teuchos::null),
      qoi_vecH23_(Teuchos::null),  qoi_vecH23_overlap_(Teuchos::null),
      qoi_vecH31_(Teuchos::null),  qoi_vecH32_(Teuchos::null),
      qoi_vecH33_(Teuchos::null),
      matR1_(Teuchos::null),       matR2_(Teuchos::null) {
    setCommunicator(comm,parlist,outStream);
    setBasis(basisPtrs,parlist,outStream);
    setDiscretization(parlist,meshMgr,outStream);
    setParallelStructure(parlist,outStream);
    setCellNodes(outStream);
  }

  void setCellNodes(PDE<Real> &pde) const {
    // Set PDE cell nodes
    pde.setCellNodes(volCellNodes_, bdryCellNodes_, bdryCellLocIds_);
  }

  /***************************************************************************/
  /* PDE assembly routines                                                   */
  /***************************************************************************/
  void assemblePDEResidual(const Teuchos::RCP<PDE<Real> > &pde,
                           const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                           const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                           const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    // Initialize residual vectors if not done so
    if ( pde_vecR_ == Teuchos::null ) { // Unique components of residual vector
      pde_vecR_ = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueResidualMap_, 1, true));
    }
    if ( pde_vecR_overlap_ == Teuchos::null ) { // Overlapping components of residual vector
      pde_vecR_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapResidualMap_, 1, true));
    }
    // Set residual vectors to zero
    pde_vecR_->scale(static_cast<Real>(0));
    pde_vecR_overlap_->scale(static_cast<Real>(0));
    // Get degrees of freedom
    Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
    int numLocalDofs = cellDofs.dimension(1);
    // Initialize res container
    Teuchos::RCP<Intrepid::FieldContainer<Real> > res;
    // Get u_coeff from u and z_coeff from z
    Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
    getCoeffFromStateVector(u_coeff,u);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
    getCoeffFromControlVector(z_coeff,z);
    // Compute PDE residual
    pde->residual(res,u_coeff,z_coeff,z_param);
    // assembly on the overlap map
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        pde_vecR_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                              0,
                                              (*res)[i*numLocalDofs+j]);
      }
    }
    // change map
    Tpetra::Export<> exporter(pde_vecR_overlap_->getMap(), pde_vecR_->getMap()); // redistribution
    pde_vecR_->doExport(*pde_vecR_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
  }

  void assemblePDEJacobian1(const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    try {
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
      getCoeffFromStateVector(u_coeff,u);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
      getCoeffFromControlVector(z_coeff,z);
      // Compute PDE Jacobian
      Teuchos::RCP<Intrepid::FieldContainer<Real> > jac;
      pde->Jacobian_1(jac,u_coeff,z_coeff,z_param); // Throw if not implemented or is zero
      // Initialize Jacobian matrices
      if ( pde_matJ1_ == Teuchos::null ) {
        pde_matJ1_ = Teuchos::rcp(new Tpetra::CrsMatrix<>(matJ1Graph_));
      }
      // Zero PDE Jacobian
      pde_matJ1_->resumeFill(); pde_matJ1_->setAllToScalar(static_cast<Real>(0));
      // Assemble PDE Jacobian
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> jacArrayRCP = jac->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          pde_matJ1_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                          cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                          jacArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      pde_matJ1_->fillComplete();
      isJ1Transposed_ = false;
      isSolverConstructed_ = false;
      isSolverTransConstructed_ = false;
    }
    catch ( Exception::Zero & ez ) {
      throw Exception::Zero(">>> (Assembler::assemblePDEJacobian1): Jacobian is zero.");
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEJacobian1): Jacobian not implemented.");
    }
  }

  void assemblePDEJacobian2(const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    try {
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
      getCoeffFromStateVector(u_coeff,u);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
      getCoeffFromControlVector(z_coeff,z);
      // Compute PDE Jacobian
      Teuchos::RCP<Intrepid::FieldContainer<Real> > jac;
      pde->Jacobian_2(jac,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      // Initialize Jacobian matrices
      if ( pde_matJ2_ == Teuchos::null ) {
        pde_matJ2_ = Teuchos::rcp(new Tpetra::CrsMatrix<>(matJ2Graph_));
      }
      // Zero PDE Jacobian
      pde_matJ2_->resumeFill(); pde_matJ2_->setAllToScalar(static_cast<Real>(0));
      // Assemble Jacobian
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> jacArrayRCP = jac->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          pde_matJ2_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                          cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                          jacArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      pde_matJ2_->fillComplete();
      isJ2Transposed_ = false;
    }
    catch ( Exception::Zero & ez ) {
      throw Exception::Zero(">>> (Assembler::assemblePDEJacobian2): Jacobian is zero.");
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEJacobian2): Jacobian not implemented.");
    }
  }

  void assemblePDEJacobian3(const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    if ( z_param != Teuchos::null ) {
      try {
        int size = static_cast<int>(z_param->size());
        // Initialize res
        std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > jac(size);
        // Get u_coeff from u and z_coeff from z
        Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
        getCoeffFromStateVector(u_coeff,u);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
        getCoeffFromControlVector(z_coeff,z);
        // Compute PDE local Jacobian wrt parametric controls
        pde->Jacobian_3(jac,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        // Initialize Jacobian storage if not done so already
        if (pde_vecJ3_ == Teuchos::null) {
          pde_vecJ3_ = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueResidualMap_, size, true));
        }
        if ( pde_vecJ3_overlap_ == Teuchos::null) {
          pde_vecJ3_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapResidualMap_, size, true));
        }
        // Assemble PDE Jacobian wrt parametric controls
        pde_vecJ3_->scale(static_cast<Real>(0));
        pde_vecJ3_overlap_->scale(static_cast<Real>(0));
        for (int k = 0; k < size; ++k) {
          Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
          int numLocalDofs = cellDofs.dimension(1);
          // assembly on the overlap map
          for (int i=0; i<numCells_; ++i) {
            for (int j=0; j<numLocalDofs; ++j) {
              pde_vecJ3_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                                     k,
                                                     (*jac[k])[i*numLocalDofs+j]);
            }
          }
          // change map
          Tpetra::Export<> exporter(pde_vecJ3_overlap_->getMap(), pde_vecJ3_->getMap()); // redistribution
          pde_vecJ3_->doExport(*pde_vecJ3_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
        }
      }
      catch ( Exception::Zero & ez ) {
        throw Exception::Zero(">>> (Assembler::assemblePDEJacobian3): Jacobian is zero.");
      }
      catch ( Exception::NotImplemented & eni ) {
        throw Exception::NotImplemented(">>> (Assembler::assemblePDEJacobian3): Jacobian not implemented.");
      }
    }
    else {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEJacobian3): Jacobian not implemented.");
    }
  }

  void assemblePDEHessian11(const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    try {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > hess; 
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
      getCoeffFromStateVector(u_coeff,u);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
      getCoeffFromControlVector(z_coeff,z);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > l_coeff = Teuchos::null;
      getCoeffFromStateVector(l_coeff,l);
      // Compute PDE Hessian
      pde->Hessian_11(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      // Initialize Hessian storage if not done so already
      if ( pde_matH11_ == Teuchos::null ) {
        pde_matH11_ = Teuchos::rcp(new Tpetra::CrsMatrix<>(matH11Graph_));
      }
      // Zero Hessian
      pde_matH11_->resumeFill(); pde_matH11_->setAllToScalar(static_cast<Real>(0));
      // Assemble PDE Hessian
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> hessArrayRCP = hess->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          pde_matH11_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                           cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                           hessArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      pde_matH11_->fillComplete();
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assemblePDEHessian11): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian11): Hessian not implemented.");
    }
  }

  void assemblePDEHessian12(const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    try {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > hess;
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
      getCoeffFromStateVector(u_coeff,u);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
      getCoeffFromControlVector(z_coeff,z);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > l_coeff = Teuchos::null;
      getCoeffFromStateVector(l_coeff,l);
      // Compute PDE Hessian
      pde->Hessian_12(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      // Initialize Hessian storage if not done so already
      if ( pde_matH12_ == Teuchos::null ) {
        pde_matH12_ = Teuchos::rcp(new Tpetra::CrsMatrix<>(matH12Graph_));
      }
      // Zero Hessian
      pde_matH12_->resumeFill(); pde_matH12_->setAllToScalar(static_cast<Real>(0));
      // Assemble PDE Hessian
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> hessArrayRCP = hess->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          pde_matH12_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                           cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                           hessArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      pde_matH12_->fillComplete();
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assemblePDEHessian12): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian12): Hessian not implemented.");
    }
  }

  void assemblePDEHessian13(const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    if ( z_param != Teuchos::null ) {
      try {
        int size = static_cast<int>(z_param->size());
        // Initialize local hessian
        std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > hess(size);
        // Get u_coeff from u, z_coeff from z and l_coeff from l
        Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
        getCoeffFromStateVector(u_coeff,u);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
        getCoeffFromControlVector(z_coeff,z);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > l_coeff = Teuchos::null;
        getCoeffFromStateVector(l_coeff,l);
        // Compute PDE local Jacobian wrt parametric controls
        pde->Hessian_13(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        // Initialize Jacobian storage if not done so already
        if (pde_vecH13_ == Teuchos::null) {
          pde_vecH13_ = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueStateMap_, size, true));
        }
        if ( pde_vecH13_overlap_ == Teuchos::null) {
          pde_vecH13_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapStateMap_, size, true));
        }
        // Assemble PDE Jacobian wrt parametric controls
        pde_vecH13_->scale(static_cast<Real>(0));
        pde_vecH13_overlap_->scale(static_cast<Real>(0));
        for (int k = 0; k < size; ++k) {
          Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
          int numLocalDofs = cellDofs.dimension(1);
          // assembly on the overlap map
          for (int i=0; i<numCells_; ++i) {
            for (int j=0; j<numLocalDofs; ++j) {
              pde_vecH13_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                                      k,
                                                      (*hess[k])[i*numLocalDofs+j]);
            }
          }
          // change map
          Tpetra::Export<> exporter(pde_vecH13_overlap_->getMap(), pde_vecH13_->getMap()); // redistribution
          pde_vecH13_->doExport(*pde_vecH13_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
        }
      }
      catch ( Exception::Zero & ez ) {
        throw Exception::Zero(">>> (Assembler::assemblePDEHessian13): Hessian is zero.");
      }
      catch ( Exception::NotImplemented & eni ) {
        throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian13): Hessian not implemented.");
      }
    }
    else {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian13): Hessian not implemented.");
    }
  }

  void assemblePDEHessian21(const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    try {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > hess;
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
      getCoeffFromStateVector(u_coeff,u);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
      getCoeffFromControlVector(z_coeff,z);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > l_coeff = Teuchos::null;
      getCoeffFromStateVector(l_coeff,l);
      // Compute PDE Hessian
      pde->Hessian_21(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      // Initialize Hessian storage if not done so already
      if ( pde_matH21_ == Teuchos::null ) {
        pde_matH21_ = Teuchos::rcp(new Tpetra::CrsMatrix<>(matH21Graph_));
      }
      // Zero Hessian
      pde_matH21_->resumeFill(); pde_matH21_->setAllToScalar(static_cast<Real>(0));
      // Assemble PDE Hessian
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> hessArrayRCP = hess->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          pde_matH21_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                           cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                           hessArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      pde_matH21_->fillComplete();
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assemblePDEHessian21): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian21): Hessian not implemented.");
    }
  }

  void assemblePDEHessian22(const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    try {
      Teuchos::RCP<Intrepid::FieldContainer<Real> > hess;
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
      getCoeffFromStateVector(u_coeff,u);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
      getCoeffFromControlVector(z_coeff,z);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > l_coeff = Teuchos::null;
      getCoeffFromStateVector(l_coeff,l);
      // Compute PDE Hessian
      pde->Hessian_22(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      // Initialize Hessian storage if not done so already
      if ( pde_matH22_ == Teuchos::null ) {
        pde_matH22_ = Teuchos::rcp(new Tpetra::CrsMatrix<>(matH22Graph_));
      }
      // Zero Hessian
      pde_matH22_->resumeFill(); pde_matH22_->setAllToScalar(static_cast<Real>(0));
      // Assemble PDE Hessian
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> hessArrayRCP = hess->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          pde_matH22_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                           cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                           hessArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      pde_matH22_->fillComplete();
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assemblePDEHessian22): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian22): Hessian not implemented.");
    }
  }

  void assemblePDEHessian23(const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    if ( z_param != Teuchos::null ) {
      try {
        int size = static_cast<int>(z_param->size());
        // Initialize local hessian
        std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > hess(size);
        // Get u_coeff from u, z_coeff from z and l_coeff from l
        Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
        getCoeffFromStateVector(u_coeff,u);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
        getCoeffFromControlVector(z_coeff,z);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > l_coeff = Teuchos::null;
        getCoeffFromStateVector(l_coeff,l);
        // Compute PDE local Jacobian wrt parametric controls
        pde->Hessian_23(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        // Initialize Jacobian storage if not done so already
        if (pde_vecH23_ == Teuchos::null) {
          pde_vecH23_ = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueControlMap_, size, true));
        }
        if ( pde_vecH23_overlap_ == Teuchos::null) {
          pde_vecH23_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapControlMap_, size, true));
        }
        // Assemble PDE Jacobian wrt parametric controls
        pde_vecH23_->scale(static_cast<Real>(0));
        pde_vecH23_overlap_->scale(static_cast<Real>(0));
        for (int k = 0; k < size; ++k) {
          Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
          int numLocalDofs = cellDofs.dimension(1);
          // assembly on the overlap map
          for (int i=0; i<numCells_; ++i) {
            for (int j=0; j<numLocalDofs; ++j) {
              pde_vecH23_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                                      k,
                                                      (*hess[k])[i*numLocalDofs+j]);
            }
          }
          // change map
          Tpetra::Export<> exporter(pde_vecH23_overlap_->getMap(), pde_vecH23_->getMap()); // redistribution
          pde_vecH23_->doExport(*pde_vecH23_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
        }
      }
      catch ( Exception::Zero & ez ) {
        throw Exception::Zero(">>> (Assembler::assemblePDEHessian23): Hessian is zero.");
      }
      catch ( Exception::NotImplemented & eni ) {
        throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian23): Hessian not implemented.");
      }
    }
    else {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian23): Hessian not implemented.");
    }
  }

  void assemblePDEHessian31(const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    assemblePDEHessian13(pde,l,u,z,z_param);
  }

  void assemblePDEHessian32(const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    assemblePDEHessian23(pde,l,u,z,z_param);
  }

  void assemblePDEHessian33(const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    if ( z_param != Teuchos::null ) {
      try {
        int size = static_cast<int>(z_param->size());
        // Initialize local hessian
        std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > tmp(size,Teuchos::null);
        std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > hess(size,tmp);
        // Get u_coeff from u, z_coeff from z and l_coeff from l
        Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
        getCoeffFromStateVector(u_coeff,u);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
        getCoeffFromControlVector(z_coeff,z);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > l_coeff = Teuchos::null;
        getCoeffFromStateVector(l_coeff,l);
        // Compute PDE local Jacobian wrt parametric controls
        pde->Hessian_33(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        // Initialize Jacobian storage if not done so already
        if (pde_matH33_ == Teuchos::null) {
          std::vector<Real> col(size,static_cast<Real>(0));
          pde_matH33_ = Teuchos::rcp(new std::vector<std::vector<Real> >(size,col));
        }
        // Assemble PDE Jacobian wrt parametric controls
        int cnt = 0, matSize = static_cast<int>(0.5*static_cast<Real>((size+1)*size));
        std::vector<Real> myHess(matSize,static_cast<Real>(0));
        std::vector<Real> globHess(matSize,static_cast<Real>(0));
        for (int k = 0; k < size; ++k) {
          for (int j = k; j < size; ++j) {
            for (int i=0; i<numCells_; ++i) {
              myHess[cnt] += (*(hess[k][j]))(i);
            }
            cnt++;
          }
        }
        Teuchos::reduceAll<int,Real>(*comm_,Teuchos::REDUCE_SUM,matSize,&myHess[0],&globHess[0]);
        cnt = 0;
        for (int k = 0; k < size; ++k) {
          for (int j = k; j < size; ++j) {
            (*pde_matH33_)[k][j] += globHess[cnt];
            if ( j != k ) { 
              (*pde_matH33_)[j][k] += globHess[cnt];
            }
            cnt++;
          }
        }
      }
      catch ( Exception::Zero & ez ) {
        throw Exception::Zero(">>> (Assembler::assemblePDEHessian33): Hessian is zero.");
      }
      catch ( Exception::NotImplemented & eni ) {
        throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian33): Hessian not implemented.");
      }
    }
    else {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian33): Hessian not implemented.");
    }
  }

  // Access routines
  Teuchos::RCP<Tpetra::MultiVector<> > getPDEResidual(void) const {
    return pde_vecR_;
  }

  Teuchos::RCP<Tpetra::CrsMatrix<> > getPDEJacobian1(const bool transpose = false) {
    if ( transpose ) {
      if (!isJ1Transposed_) {
        // Create matrix transposes.
        Tpetra::RowMatrixTransposer<> transposerJ1(pde_matJ1_);
        pde_matJ1_trans_ = transposerJ1.createTranspose();
        isJ1Transposed_ = true;
      }
      return pde_matJ1_trans_;
    }
    else {
      return pde_matJ1_;
    }
  }

  Teuchos::RCP<Tpetra::CrsMatrix<> > getPDEJacobian2(const bool transpose = false) {
    if ( transpose ) {
      if (!isJ2Transposed_) {
        // Create matrix transposes.
        Tpetra::RowMatrixTransposer<> transposerJ2(pde_matJ2_);
        pde_matJ2_trans_ = transposerJ2.createTranspose();
        isJ2Transposed_ = true;
      }
      return pde_matJ2_trans_;
    }
    else {
      return pde_matJ2_;
    }
  }

  Teuchos::RCP<Tpetra::CrsMatrix<> > getPDEHessian11(void) const {
    return pde_matH11_;
  }

  Teuchos::RCP<Tpetra::CrsMatrix<> > getPDEHessian12(void) const {
    return pde_matH12_;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getPDEHessian13(void) const {
    return pde_vecH13_;
  }

  Teuchos::RCP<Tpetra::CrsMatrix<> > getPDEHessian21(void) const {
    return pde_matH21_;
  }

  Teuchos::RCP<Tpetra::CrsMatrix<> > getPDEHessian22(void) const {
    return pde_matH22_;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getPDEHessian23(void) const {
    return pde_vecH23_;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getPDEHessian31(void) const {
    return pde_vecH13_;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getPDEHessian32(void) const {
    return pde_vecH23_;
  }

  Teuchos::RCP<std::vector<std::vector<Real> > > getPDEHessian33(void) const {
    return pde_matH33_;
  }

  // Application routines
  void applyPDEJacobian1(const Teuchos::RCP<Tpetra::MultiVector<> > &Jv,
                         const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                         const bool transpose = false) {
    if ( transpose ) {
      if (!isJ1Transposed_) {
        // Create matrix transposes.
        Tpetra::RowMatrixTransposer<> transposerJ1(pde_matJ1_);
        pde_matJ1_trans_ = transposerJ1.createTranspose();
        isJ1Transposed_ = true;
      }
      pde_matJ1_trans_->apply(*v,*Jv);
    }
    else {
      pde_matJ1_->apply(*v,*Jv);
    }
  }

  void applyInverseJacobian1(const Teuchos::RCP<Tpetra::MultiVector<> > &u,
                             const Teuchos::RCP<const Tpetra::MultiVector<> > &r,
                             const bool transpose = false) {
    if ( transpose ) {
      if (!isJ1Transposed_) {
        // Create matrix transposes.
        Tpetra::RowMatrixTransposer<> transposerJ1(pde_matJ1_);
        pde_matJ1_trans_ = transposerJ1.createTranspose();
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

  void applyPDEJacobian2(const Teuchos::RCP<Tpetra::MultiVector<> > &Jv,
                         const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                         const bool transpose = false) {
    if ( transpose ) {
      if (!isJ2Transposed_) {
        // Create matrix transposes.
        Tpetra::RowMatrixTransposer<> transposerJ2(pde_matJ2_);
        pde_matJ2_trans_ = transposerJ2.createTranspose();
        isJ2Transposed_ = true;
      }
      pde_matJ2_trans_->apply(*v,*Jv);
    }
    else {
      pde_matJ2_->apply(*v,*Jv);
    }
  }

  void applyPDEJacobian3(const Teuchos::RCP<Tpetra::MultiVector<> > &Jv,
                         const Teuchos::RCP<const std::vector<Real> > &v) const {
    Jv->putScalar(static_cast<Real>(0));
    const size_t size = static_cast<size_t>(v->size());
    for (size_t i = 0; i < size; ++i) {
      Teuchos::ArrayView<const size_t> col(&i,1);
      Jv->update((*v)[i],*(pde_vecJ3_->subView(col)),static_cast<Real>(1));
    }
  }

  void applyPDEAdjointJacobian3(const Teuchos::RCP<std::vector<Real> > &Jv,
                                const Teuchos::RCP<const Tpetra::MultiVector<> > &v) const {
    Teuchos::Array<Real> val(1,0);
    const size_t size = static_cast<size_t>(Jv->size());
    for (size_t i = 0; i < size; ++i) {
      Teuchos::ArrayView<const size_t> col(&i,1);
      pde_vecJ3_->subView(col)->dot(*v, val.view(0,1));
      (*Jv)[i] = val[0];
    }
  }

  void applyPDEHessian11(const Teuchos::RCP<Tpetra::MultiVector<> > &Hv,
                         const Teuchos::RCP<const Tpetra::MultiVector<> > &v) const {
    pde_matH11_->apply(*v,*Hv);
  }

  void applyPDEHessian12(const Teuchos::RCP<Tpetra::MultiVector<> > &Hv,
                         const Teuchos::RCP<const Tpetra::MultiVector<> > &v) const {
    pde_matH12_->apply(*v,*Hv);
  }

  void applyPDEHessian13(const Teuchos::RCP<std::vector<Real> > &Hv,
                         const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                         const bool zeroOut = true) const {
    const size_t size = static_cast<size_t>(Hv->size());
    if ( zeroOut ) {
      Hv->assign(size,static_cast<Real>(0));
    }
    Teuchos::Array<Real> val(1,0);
    for (size_t i = 0; i < size; ++i) {
      Teuchos::ArrayView<const size_t> col(&i,1);
      pde_vecH13_->subView(col)->dot(*v, val.view(0,1));
      (*Hv)[i] += val[0];
    }
  }

  void applyPDEHessian21(const Teuchos::RCP<Tpetra::MultiVector<> > &Hv,
                         const Teuchos::RCP<const Tpetra::MultiVector<> > &v) const {
    pde_matH21_->apply(*v,*Hv);
  }

  void applyPDEHessian22(const Teuchos::RCP<Tpetra::MultiVector<> > &Hv,
                         const Teuchos::RCP<const Tpetra::MultiVector<> > &v) const {
    pde_matH22_->apply(*v,*Hv);
  }

  void applyPDEHessian23(const Teuchos::RCP<std::vector<Real> > &Hv,
                         const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                         const bool zeroOut = true) const {
    const size_t size = static_cast<size_t>(Hv->size());
    if ( zeroOut ) {
      Hv->assign(size,static_cast<Real>(0));
    }
    Teuchos::Array<Real> val(1,0);
    for (size_t i = 0; i < size; ++i) {
      Teuchos::ArrayView<const size_t> col(&i,1);
      pde_vecH23_->subView(col)->dot(*v, val.view(0,1));
      (*Hv)[i] += val[0];
    }
  }

  void applyPDEHessian31(const Teuchos::RCP<Tpetra::MultiVector<> > &Hv,
                         const Teuchos::RCP<const std::vector<Real> > &v,
                         const bool zeroOut = true) const {
    if ( zeroOut ) {
      Hv->putScalar(static_cast<Real>(0));
    }
    const size_t size = static_cast<size_t>(v->size());
    for (size_t i = 0; i < size; ++i) {
      Teuchos::ArrayView<const size_t> col(&i,1);
      Hv->update((*v)[i],*(pde_vecH13_->subView(col)),static_cast<Real>(1));
    }
  }

  void applyPDEHessian32(const Teuchos::RCP<Tpetra::MultiVector<> > &Hv,
                         const Teuchos::RCP<const std::vector<Real> > &v,
                         const bool zeroOut = true) const {
    if ( zeroOut ) {
      Hv->putScalar(static_cast<Real>(0));
    }
    const size_t size = static_cast<size_t>(v->size());
    for (size_t i = 0; i < size; ++i) {
      Teuchos::ArrayView<const size_t> col(&i,1);
      Hv->update((*v)[i],*(pde_vecH23_->subView(col)),static_cast<Real>(1));
    }
  }

  void applyPDEHessian33(const Teuchos::RCP<std::vector<Real> > &Hv,
                         const Teuchos::RCP<const std::vector<Real> > &v,
                         const bool zeroOut = true ) const {
    const int size = Hv->size();
    if ( zeroOut ) {
      Hv->assign(size,static_cast<Real>(0));
    }
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
        (*Hv)[i] += (*pde_matH33_)[i][j]*(*v)[j];
      }
    }
  }
  /***************************************************************************/
  /* End of PDE assembly routines.                                           */
  /***************************************************************************/


  /***************************************************************************/
  /* QoI assembly routines                                                   */
  /***************************************************************************/
  Real assembleQoIValue(const Teuchos::RCP<QoI<Real> > &qoi,
                        const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                        const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                        const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    Real val(0);
    try {
      // Integrate obj object
      Teuchos::RCP<Intrepid::FieldContainer<Real> > locVal;
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
      getCoeffFromStateVector(u_coeff,u);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
      getCoeffFromControlVector(z_coeff,z);
      // Get OBJ_CELL value
      val = qoi->value(locVal,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      // Assembly
      if ( locVal != Teuchos::null ) {
        Real myval(0), gval(0);
        for (int i=0; i<numCells_; ++i) {
          myval += (*locVal)(i);
        }
        Teuchos::reduceAll<int,Real>(*comm_,Teuchos::REDUCE_SUM,1,&myval,&gval);
        val += gval;
      }
    }
    catch ( Exception::Zero & ez ) {
      val = static_cast<Real>(0);
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIValue): Value not implemented.");
    }
    return val;
  }

  void assembleQoIGradient1(const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    try {
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
      getCoeffFromStateVector(u_coeff,u);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      Teuchos::RCP<Intrepid::FieldContainer<Real> > locGrad;
      qoi->gradient_1(locGrad,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      // Initialize state QoI gradient vectors
      if ( qoi_vecG1_ == Teuchos::null ) {
        qoi_vecG1_ = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueStateMap_, 1, true));
      }
      if ( qoi_vecG1_overlap_ == Teuchos::null ) {
        qoi_vecG1_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapStateMap_, 1, true));
      }
      // Assembly in to the overlap gradient
      qoi_vecG1_->scale(static_cast<Real>(0));
      qoi_vecG1_overlap_->scale(static_cast<Real>(0));
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          qoi_vecG1_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                                 0,
                                                 (*locGrad)[i*numLocalDofs+j]);
        }
      }
      // Change to local map
      Tpetra::Export<> exporter(qoi_vecG1_overlap_->getMap(), qoi_vecG1_->getMap()); // redistribution
      qoi_vecG1_->doExport(*qoi_vecG1_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch ( Exception::Zero & ez ) {
      throw Exception::Zero(">>> (Assembler::assembleQoIGradient1): Gradient is zero.");
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIGradient1): Gradient not implemented.");
    }
  }

  void assembleQoIGradient2(const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    try {
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
      getCoeffFromStateVector(u_coeff,u);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      Teuchos::RCP<Intrepid::FieldContainer<Real> > locGrad;
      qoi->gradient_2(locGrad,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      // Initialize control gradient vectors
      if ( qoi_vecG2_ == Teuchos::null ) {
        qoi_vecG2_ = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueControlMap_, 1, true));
      }
      if ( qoi_vecG2_overlap_ == Teuchos::null ) {
        qoi_vecG2_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapControlMap_, 1, true));
      }
      // Assembly in to the overlap gradient
      qoi_vecG2_->scale(static_cast<Real>(0));
      qoi_vecG2_overlap_->scale(static_cast<Real>(0));
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          qoi_vecG2_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                                 0,
                                                 (*locGrad)[i*numLocalDofs+j]);
        }
      }
      // Change to local map
      Tpetra::Export<> exporter(qoi_vecG2_overlap_->getMap(), qoi_vecG2_->getMap()); // redistribution
      qoi_vecG2_->doExport(*qoi_vecG2_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch ( Exception::Zero & ez ) {
      throw Exception::Zero(">>> (Assembler::assembleQoIGradient2): Gradient is zero.");
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIGradient2): Gradient not implemented.");
    }
  }

  void assembleQoIGradient3(const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    if ( z_param != Teuchos::null ) {
      const int size = z_param->size();
      if ( qoi_vecG3_ == Teuchos::null ) {
        qoi_vecG3_ = Teuchos::rcp(new std::vector<Real>(size,0));
      }
      try {
        qoi_vecG3_->assign(size,0);
        // Initialize local gradient storage
        std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > locGrad(size,Teuchos::null);
        // Get u_coeff from u and z_coeff from z
        Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
        getCoeffFromStateVector(u_coeff,u);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
        getCoeffFromControlVector(z_coeff,z);
        // Compute gradient
        (*qoi_vecG3_) = qoi->gradient_3(locGrad,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        // Assembly
        std::vector<Real> myGrad(size,0), globGrad(size,0);
        for (int j = 0; j < size; ++j) {
          if ( locGrad[j] != Teuchos::null ) {
            for (int i=0; i<numCells_; ++i) {
              myGrad[j] += (*locGrad[j])(i);
            }
          }
        }
        Teuchos::reduceAll<int,Real>(*comm_,Teuchos::REDUCE_SUM,size,&myGrad[0],&globGrad[0]);
        for (int j = 0; j < size; ++j) {
          (*qoi_vecG3_)[j] += globGrad[j];
        }
      }
      catch ( Exception::Zero & ez ) {
        qoi_vecG3_->assign(size,0);
      }
      catch ( Exception::NotImplemented & eni ) {
        throw Exception::NotImplemented(">>> (Assembler::assembleQoIGradient3): Gradient not implemented.");
      }
    }
    else {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIGradient3): Gradient not implemented.");
    }
  }

  void assembleQoIHessVec11(const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    try {
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > locHess;
      Teuchos::RCP<Intrepid::FieldContainer<Real> > v_coeff = Teuchos::null;
      getCoeffFromStateVector(v_coeff,v);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
      getCoeffFromStateVector(u_coeff,u);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi->HessVec_11(locHess, v_coeff, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
      // Initialize state-state HessVec vectors
      if ( qoi_vecH11_ == Teuchos::null ) {
        qoi_vecH11_ = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueStateMap_, 1, true));
      }
      if ( qoi_vecH11_overlap_ == Teuchos::null ) {
        qoi_vecH11_overlap_  = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapStateMap_, 1, true));
      }
      // Assembly in to the overlap gradient
      qoi_vecH11_->scale(static_cast<Real>(0));
      qoi_vecH11_overlap_->scale(static_cast<Real>(0));
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          qoi_vecH11_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                                  0,
                                                  (*locHess)[i*numLocalDofs+j]);
        }
      }
      // Change to local map
      Tpetra::Export<> exporter(qoi_vecH11_overlap_->getMap(), qoi_vecH11_->getMap()); // redistribution
      qoi_vecH11_->doExport(*qoi_vecH11_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assembleQoIHessVec11): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec11): Hessian not implemented.");
    }
  }

  void assembleQoIHessVec12(const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    try {
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > locHess;
      Teuchos::RCP<Intrepid::FieldContainer<Real> > v_coeff = Teuchos::null;
      getCoeffFromControlVector(v_coeff,v);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
      getCoeffFromStateVector(u_coeff,u);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi->HessVec_12(locHess, v_coeff, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
      // Initialize state-control HessVec vectors
      if ( qoi_vecH12_ == Teuchos::null ) {
        qoi_vecH12_ = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueStateMap_, 1, true));
      }
      if ( qoi_vecH12_overlap_ == Teuchos::null ) {
        qoi_vecH12_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapStateMap_, 1, true));
      }
      // Assembly in to the overlap gradient
      qoi_vecH12_->scale(static_cast<Real>(0));
      qoi_vecH12_overlap_->scale(static_cast<Real>(0));
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          qoi_vecH12_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                                  0,
                                                  (*locHess)[i*numLocalDofs+j]);
        }
      }
      // Change to local map
      Tpetra::Export<> exporter(qoi_vecH12_overlap_->getMap(), qoi_vecH12_->getMap()); // redistribution
      qoi_vecH12_->doExport(*qoi_vecH12_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assembleQoIHessVec12): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec12): Hessian not implemented.");
    }
  }

  void assembleQoIHessVec13(const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const std::vector<Real> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > &z_param = Teuchos::null) {
    try {
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > locHess;
      Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
      getCoeffFromStateVector(u_coeff,u);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi->HessVec_13(locHess, v, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
      // Initialize state-control HessVec vectors
      if ( qoi_vecH13_ == Teuchos::null ) {
        qoi_vecH13_ = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueStateMap_, 1, true));
      }
      if ( qoi_vecH13_overlap_ == Teuchos::null ) {
        qoi_vecH13_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapStateMap_, 1, true));
      }
      // Assembly in to the overlap gradient
      qoi_vecH13_->scale(static_cast<Real>(0));
      qoi_vecH13_overlap_->scale(static_cast<Real>(0));
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          qoi_vecH13_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                                  0,
                                                  (*locHess)[i*numLocalDofs+j]);
        }
      }
      // Change to local map
      Tpetra::Export<> exporter(qoi_vecH13_overlap_->getMap(), qoi_vecH13_->getMap()); // redistribution
      qoi_vecH13_->doExport(*qoi_vecH13_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assembleQoIHessVec12): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec12): Hessian not implemented.");
    }
  }

  void assembleQoIHessVec21(const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    try {
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > locHess;
      Teuchos::RCP<Intrepid::FieldContainer<Real> > v_coeff = Teuchos::null;
      getCoeffFromStateVector(v_coeff,v);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
      getCoeffFromStateVector(u_coeff,u);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi->HessVec_21(locHess, v_coeff, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
      // Initialize control-state HessVec vectors
      if ( qoi_vecH21_ == Teuchos::null ) {
        qoi_vecH21_ = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueControlMap_, 1, true));
      }
      if ( qoi_vecH21_overlap_ == Teuchos::null ) {
        qoi_vecH21_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapControlMap_, 1, true));
      }
      // Assembly in to the overlap gradient
      qoi_vecH21_->scale(static_cast<Real>(0));
      qoi_vecH21_overlap_->scale(static_cast<Real>(0));
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          qoi_vecH21_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                                  0,
                                                  (*locHess)[i*numLocalDofs+j]);
        }
      }
      // Change to local map
      Tpetra::Export<> exporter(qoi_vecH21_overlap_->getMap(), qoi_vecH21_->getMap()); // redistribution
      qoi_vecH21_->doExport(*qoi_vecH21_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assembleQoIHessVec21): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec21): Hessian not implemented.");
    }
  }

  void assembleQoIHessVec22(const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    try {
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > locHess;
      Teuchos::RCP<Intrepid::FieldContainer<Real> > v_coeff = Teuchos::null;
      getCoeffFromControlVector(v_coeff,v);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
      getCoeffFromStateVector(u_coeff,u);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi->HessVec_22(locHess, v_coeff, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
      // Initialize control-control HessVec vectors
      if ( qoi_vecH22_ == Teuchos::null ) {
        qoi_vecH22_ = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueControlMap_, 1, true));
      }
      if ( qoi_vecH22_overlap_ == Teuchos::null ) {
        qoi_vecH22_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapControlMap_, 1, true));
      }
      // Assembly in to the overlap gradient
      qoi_vecH22_->scale(static_cast<Real>(0));
      qoi_vecH22_overlap_->scale(static_cast<Real>(0));
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          qoi_vecH22_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                                  0,
                                                  (*locHess)[i*numLocalDofs+j]);
        }
      }
      // Change to local map
      Tpetra::Export<> exporter(qoi_vecH22_overlap_->getMap(), qoi_vecH22_->getMap()); // redistribution
      qoi_vecH22_->doExport(*qoi_vecH22_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assembleQoIHessVec22): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec22): Hessian not implemented.");
    }
  }

  void assembleQoIHessVec23(const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const std::vector<Real> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    try {
      // Get u_coeff from u and z_coeff from z
      Teuchos::RCP<Intrepid::FieldContainer<Real> > locHess;
      Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
      getCoeffFromStateVector(u_coeff,u);
      Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi->HessVec_23(locHess, v, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
      // Initialize control-control HessVec vectors
      if ( qoi_vecH23_ == Teuchos::null ) {
        qoi_vecH23_ = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueControlMap_, 1, true));
      }
      if ( qoi_vecH23_overlap_ == Teuchos::null ) {
        qoi_vecH23_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapControlMap_, 1, true));
      }
      // Assembly in to the overlap gradient
      qoi_vecH23_->scale(static_cast<Real>(0));
      qoi_vecH23_overlap_->scale(static_cast<Real>(0));
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          qoi_vecH23_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                                  0,
                                                  (*locHess)[i*numLocalDofs+j]);
        }
      }
      // Change to local map
      Tpetra::Export<> exporter(qoi_vecH23_overlap_->getMap(), qoi_vecH23_->getMap()); // redistribution
      qoi_vecH23_->doExport(*qoi_vecH23_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assembleQoIHessVec23): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec23): Hessian not implemented.");
    }
  }

  void assembleQoIHessVec31(const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    if ( z_param != Teuchos::null ) {
      const int size = z_param->size();
      if ( qoi_vecH31_ == Teuchos::null ) {
        qoi_vecH31_ = Teuchos::rcp(new std::vector<Real>(size,0));
      }
      try {
        qoi_vecH31_->assign(size,0);
        // Initialize local gradient storage
        std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > locHess(size,Teuchos::null);
        // Get u_coeff from u and z_coeff from z
        Teuchos::RCP<Intrepid::FieldContainer<Real> > v_coeff = Teuchos::null;
        getCoeffFromStateVector(v_coeff,v);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
        getCoeffFromStateVector(u_coeff,u);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
        getCoeffFromControlVector(z_coeff,z);
        // Compute gradient
        (*qoi_vecH31_) = qoi->HessVec_31(locHess,v_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        // Assembly
        std::vector<Real> myHess(size,0), globHess(size,0);
        for (int j = 0; j < size; ++j) {
          if ( locHess[j] != Teuchos::null ) {
            for (int i=0; i<numCells_; ++i) {
              myHess[j] += (*locHess[j])(i);
            }
          }
        }
        Teuchos::reduceAll<int,Real>(*comm_,Teuchos::REDUCE_SUM,size,&myHess[0],&globHess[0]);
        for (int j = 0; j < size; ++j) {
          (*qoi_vecH31_)[j] += globHess[j];
        }
      }
      catch ( Exception::Zero & ez ) {
        qoi_vecH31_->assign(size,0);
      }
      catch ( Exception::NotImplemented & eni ) {
        throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec31): HessVec not implemented.");
      }
    }
    else {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec31): HessVec not implemented.");
    }
  }

  void assembleQoIHessVec32(const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    if ( z_param != Teuchos::null ) {
      const int size = z_param->size();
      if ( qoi_vecH32_ == Teuchos::null ) {
        qoi_vecH32_ = Teuchos::rcp(new std::vector<Real>(size,0));
      }
      try {
        qoi_vecH32_->assign(size,0);
        // Initialize local hessian times a vector storage
        std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > locHess(size,Teuchos::null);
        // Get u_coeff from u and z_coeff from z
        Teuchos::RCP<Intrepid::FieldContainer<Real> > v_coeff = Teuchos::null;
        getCoeffFromControlVector(v_coeff,v);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
        getCoeffFromStateVector(u_coeff,u);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
        getCoeffFromControlVector(z_coeff,z);
        // Compute local hessian times a vector
        (*qoi_vecH32_) = qoi->HessVec_32(locHess,v_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        // Assembly
        std::vector<Real> myHess(size,0), globHess(size,0);
        for (int j = 0; j < size; ++j) {
          if ( locHess[j] != Teuchos::null ) {
            for (int i=0; i<numCells_; ++i) {
              myHess[j] += (*locHess[j])(i);
            }
          }
        }
        Teuchos::reduceAll<int,Real>(*comm_,Teuchos::REDUCE_SUM,size,&myHess[0],&globHess[0]);
        for (int j = 0; j < size; ++j) {
          (*qoi_vecH32_)[j] += globHess[j];
        }
      }
      catch ( Exception::Zero & ez ) {
        qoi_vecH32_->assign(size,0);
      }
      catch ( Exception::NotImplemented & eni ) {
        throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec32): HessVec not implemented.");
      }
    }
    else {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec32): HessVec not implemented.");
    }
  }

  void assembleQoIHessVec33(const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const std::vector<Real> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > &z_param = Teuchos::null) {
    if ( z_param != Teuchos::null ) {
      const int size = z_param->size();
      if ( qoi_vecH33_ == Teuchos::null ) {
        qoi_vecH33_ = Teuchos::rcp(new std::vector<Real>(size,0));
      }
      try {
        qoi_vecH33_->assign(size,0);
        // Initialize local hessian times a vector storage
        std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > locHess(size,Teuchos::null);
        // Get u_coeff from u and z_coeff from z
        Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
        getCoeffFromStateVector(u_coeff,u);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
        getCoeffFromControlVector(z_coeff,z);
        // Compute local hessian times a vector
        (*qoi_vecH33_) = qoi->HessVec_33(locHess,v,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        // Assembly
        std::vector<Real> myHess(size,0), globHess(size,0);
        for (int j = 0; j < size; ++j) {
          if ( locHess[j] != Teuchos::null ) {
            for (int i=0; i<numCells_; ++i) {
              myHess[j] += (*locHess[j])(i);
            }
          }
        }
        Teuchos::reduceAll<int,Real>(*comm_,Teuchos::REDUCE_SUM,size,&myHess[0],&globHess[0]);
        for (int j = 0; j < size; ++j) {
          (*qoi_vecH33_)[j] += globHess[j];
        }
      }
      catch ( Exception::Zero & ez ) {
        qoi_vecH33_->assign(size,0);
      }
      catch ( Exception::NotImplemented & eni ) {
        throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec33): HessVec not implemented.");
      }
    }
    else {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec33): HessVec not implemented.");
    }
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getQoIGradient1(void) const {
    return qoi_vecG1_;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getQoIGradient2(void) const {
    return qoi_vecG2_;
  }

  Teuchos::RCP<std::vector<Real> > getQoIGradient3(void) const {
    return qoi_vecG3_;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getQoIHessVec11(void) const {
    return qoi_vecH11_;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getQoIHessVec12(void) const {
    return qoi_vecH12_;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getQoIHessVec13(void) const {
    return qoi_vecH13_;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getQoIHessVec21(void) const {
    return qoi_vecH21_;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getQoIHessVec22(void) const {
    return qoi_vecH22_;
  }

  Teuchos::RCP<Tpetra::MultiVector<> > getQoIHessVec23(void) const {
    return qoi_vecH23_;
  }

  Teuchos::RCP<std::vector<Real> > getQoIHessVec31(void) const {
    return qoi_vecH31_;
  }

  Teuchos::RCP<std::vector<Real> > getQoIHessVec32(void) const {
    return qoi_vecH32_;
  }

  Teuchos::RCP<std::vector<Real> > getQoIHessVec33(void) const {
    return qoi_vecH33_;
  }
  /***************************************************************************/
  /* End QoI assembly routines                                               */
  /***************************************************************************/


  /***************************************************************************/
  /* Assemble and apply Riesz operator corresponding to simulation variables */
  /***************************************************************************/
  void assemblePDERieszMap1(const Teuchos::RCP<PDE<Real> > &pde) {
    try {
      // Compute local state Riesz matrix
      Teuchos::RCP<Intrepid::FieldContainer<Real> > riesz;
      pde->RieszMap_1(riesz); // Throw if not implemented or zero
      // Initialize Riesz matrix if not done so already
      if ( matR1_ == Teuchos::null ) {
      matR1_ = Teuchos::rcp(new Tpetra::CrsMatrix<>(matR1Graph_));
      }
      matR1_->resumeFill(); matR1_->setAllToScalar(static_cast<Real>(0));
      // Assemble Riesz matrix
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> rieszArrayRCP = riesz->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          matR1_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                      cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                      rieszArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      matR1_->fillComplete();
      isSolverR1Constructed_ = false;
    }
    catch ( Exception::NotImplemented & eni ) {
      isSolverR1Constructed_ = true;
      throw Exception::NotImplemented(">>> (Assembler::assemblePDERieszMap1): Riesz map not implemented!");
    }
    catch ( Exception::Zero & ez ) {
      isSolverR1Constructed_ = true;
      throw Exception::Zero(">>> (Assembler::assemblePDERieszMap1): Riesz map is zero!");
    }
  }

  Teuchos::RCP<Tpetra::CrsMatrix<> > getPDERieszMap1(void) const {
    return matR1_;
  }

  void applyPDERieszMap1(const Teuchos::RCP<Tpetra::MultiVector<> > &Rv,
                         const Teuchos::RCP<const Tpetra::MultiVector<> > &v) const {
    if ( matR1_ == Teuchos::null ) {
      Rv->scale(static_cast<Real>(1),*v);
    }
    else {
      matR1_->apply(*v,*Rv);
    }
  }

  void applyPDEInverseRieszMap1(const Teuchos::RCP<Tpetra::MultiVector<> > &Rv,
                                const Teuchos::RCP<const Tpetra::MultiVector<> > &v) {
    if ( matR1_ == Teuchos::null ) {
      Rv->scale(static_cast<Real>(1),*v);
    }
    else {
      if (!isSolverR1Constructed_) {
        constructSolverR1();
        isSolverR1Constructed_ = true;
      }
      solverR1_->setX(Rv);
      solverR1_->setB(v);
      solverR1_->solve();
    }
  }
  /***************************************************************************/
  /* End of functions for Riesz operator of simulation variables.            */
  /***************************************************************************/


  /***************************************************************************/
  /* Assemble and apply Riesz operator corresponding to optimization         */
  /* variables                                                               */
  /***************************************************************************/
  void assemblePDERieszMap2(const Teuchos::RCP<PDE<Real> > &pde) {
    try {
      // Compute local control Riesz matrix
      Teuchos::RCP<Intrepid::FieldContainer<Real> > riesz;
      pde->RieszMap_2(riesz); // Throw if not implemented or zero
      // Initialize Riesz matrix if not done so already
      if ( matR2_ == Teuchos::null ) {
        matR2_ = Teuchos::rcp(new Tpetra::CrsMatrix<>(matR2Graph_));
      }
      // Assemble Riesz matrix
      matR2_->resumeFill(); matR2_->setAllToScalar(static_cast<Real>(0));
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> rieszArrayRCP = riesz->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          matR2_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                      cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                      rieszArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      matR2_->fillComplete();
      isSolverR2Constructed_ = false;
    }
    catch ( Exception::NotImplemented & eni ) {
      isSolverR2Constructed_ = true;
      throw Exception::NotImplemented(">>> (Assembler::assemblePDERieszMap2): Riesz map not implemented!");
    }
    catch ( Exception::Zero & ez ) {
      isSolverR2Constructed_ = true;
      throw Exception::Zero(">>> (Assembler::assemblePDERieszMap2): Riesz map is zero!");
    }
  }

  Teuchos::RCP<Tpetra::CrsMatrix<> > getPDERieszMap2(void) const {
    return matR2_;
  }

  void applyPDERieszMap2(const Teuchos::RCP<Tpetra::MultiVector<> > &Rv,
                         const Teuchos::RCP<const Tpetra::MultiVector<> > &v) const {
    if ( matR2_ == Teuchos::null ) {
      Rv->scale(static_cast<Real>(1),*v);
    }
    else {
      matR2_->apply(*v,*Rv);
    }
  }

  void applyPDEInverseRieszMap2(const Teuchos::RCP<Tpetra::MultiVector<> > &Rv,
                                const Teuchos::RCP<const Tpetra::MultiVector<> > &v) {
    if ( matR2_ == Teuchos::null ) {
      Rv->scale(static_cast<Real>(1),*v);
    }
    else {
      if (!isSolverR2Constructed_) {
        constructSolverR2();
        isSolverR2Constructed_ = true;
      }
      solverR2_->setX(Rv);
      solverR2_->setB(v);
      solverR2_->solve();
    }
  }
  /***************************************************************************/
  /* End of functions for Riesz operator of optimization variables.          */
  /***************************************************************************/


  /***************************************************************************/
  /* Output routines.                                                        */
  /***************************************************************************/
  void printMeshData(std::ostream &outStream) const {
    Teuchos::RCP<Intrepid::FieldContainer<Real> > nodesPtr = meshMgr_->getNodes();
    Teuchos::RCP<Intrepid::FieldContainer<int> >  cellToNodeMapPtr = meshMgr_->getCellToNodeMap();
    Intrepid::FieldContainer<Real>  &nodes = *nodesPtr;
    Intrepid::FieldContainer<int>   &cellToNodeMap = *cellToNodeMapPtr;
    if ( verbose_ && myRank_ == 0) {
      outStream << std::endl;
      outStream << "Number of nodes = " << meshMgr_->getNumNodes() << std::endl;
      outStream << "Number of cells = " << meshMgr_->getNumCells() << std::endl;
      outStream << "Number of edges = " << meshMgr_->getNumEdges() << std::endl;
    }
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
    matWriter.writeSparseFile("jacobian1", pde_matJ1_);
    matWriter.writeSparseFile("jacobian2", pde_matJ2_);
    matWriter.writeDenseFile("residual", pde_vecR_);
  }

  void outputTpetraVector(const Teuchos::RCP<const Tpetra::MultiVector<> > &vec,
                          const std::string &filename) const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> > vecWriter;
    vecWriter.writeDenseFile(filename, vec);
  }
  /***************************************************************************/
  /* End of output routines.                                                 */
  /***************************************************************************/


  /***************************************************************************/
  /* Vector generation routines.                                             */
  /***************************************************************************/
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
  /***************************************************************************/
  /* End of vector generation routines.                                      */
  /***************************************************************************/

}; // class Assembler

#endif
