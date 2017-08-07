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

/*! \file  assembler.hpp
    \brief Finite element assembly class.
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

#include "fe.hpp"
#include "fe_curl.hpp"
#include "pde.hpp"
#include "qoi.hpp"
#include "dofmanager.hpp"
#include "meshmanager.hpp"
#include "fieldhelper.hpp"

//// Global Timers.
#ifdef ROL_TIMERS
namespace ROL {
  namespace PDEOPT {
    Teuchos::RCP<Teuchos::Time> AssemblePDEResidual       = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Residual");
    Teuchos::RCP<Teuchos::Time> AssemblePDEJacobian1      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Jacobian1");
    Teuchos::RCP<Teuchos::Time> AssemblePDEJacobian2      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Jacobian2");
    Teuchos::RCP<Teuchos::Time> AssemblePDEJacobian3      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Jacobian3");
    Teuchos::RCP<Teuchos::Time> AssemblePDEHessian11      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian11");
    Teuchos::RCP<Teuchos::Time> AssemblePDEHessian12      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian12");
    Teuchos::RCP<Teuchos::Time> AssemblePDEHessian13      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian13");
    Teuchos::RCP<Teuchos::Time> AssemblePDEHessian21      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian21");
    Teuchos::RCP<Teuchos::Time> AssemblePDEHessian22      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian22");
    Teuchos::RCP<Teuchos::Time> AssemblePDEHessian23      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian23");
    Teuchos::RCP<Teuchos::Time> AssemblePDEHessian31      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian31");
    Teuchos::RCP<Teuchos::Time> AssemblePDEHessian32      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian32");
    Teuchos::RCP<Teuchos::Time> AssemblePDEHessian33      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian33");
    Teuchos::RCP<Teuchos::Time> AssembleQOIValue          = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Value");
    Teuchos::RCP<Teuchos::Time> AssembleQOIGradient1      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Gradient1");
    Teuchos::RCP<Teuchos::Time> AssembleQOIGradient2      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Gradient2");
    Teuchos::RCP<Teuchos::Time> AssembleQOIGradient3      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Gradient3");
    Teuchos::RCP<Teuchos::Time> AssembleQOIHessVec11      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec11");
    Teuchos::RCP<Teuchos::Time> AssembleQOIHessVec12      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec12");
    Teuchos::RCP<Teuchos::Time> AssembleQOIHessVec13      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec13");
    Teuchos::RCP<Teuchos::Time> AssembleQOIHessVec21      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec21");
    Teuchos::RCP<Teuchos::Time> AssembleQOIHessVec22      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec22");
    Teuchos::RCP<Teuchos::Time> AssembleQOIHessVec23      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec23");
    Teuchos::RCP<Teuchos::Time> AssembleQOIHessVec31      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec31");
    Teuchos::RCP<Teuchos::Time> AssembleQOIHessVec32      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec32");
    Teuchos::RCP<Teuchos::Time> AssembleQOIHessVec33      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec33");
  }
}
#endif

template<class Real>
class Assembler {

  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Map<>::node_type NO;
  typedef Tpetra::MultiVector<Real,LO,GO,NO> MV;
  typedef Tpetra::Operator<Real,LO,GO,NO> OP;

private:
  // Timers
//  Teuchos::RCP<Teuchos::Time::Time> timerSolverFactorization_;
//  Teuchos::RCP<Teuchos::Time::Time> timerSolverSubstitution_;
//  Teuchos::RCP<Teuchos::Time::Time> timerAssemblyNonlinear_;
//  Teuchos::RCP<Teuchos::Time::Time> timerSolverUpdate_;

  // Set in Constructor.
  bool verbose_;
  bool isJ1Transposed_, isJ2Transposed_;

  // Set in SetCommunicator.
  Teuchos::RCP<const Teuchos::Comm<int> > comm_;
  int myRank_, numProcs_;

  // Set in SetBasis.
  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;

  // Set in SetDiscretization.
  Teuchos::RCP<MeshManager<Real> > meshMgr_;
  Teuchos::RCP<DofManager<Real> >  dofMgr_;

  // Set in SetParallelStructure.
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

  // Set in SetCellNodes.
  Teuchos::RCP<Intrepid::FieldContainer<Real> > volCellNodes_;
  Teuchos::RCP<std::vector<std::vector<std::vector<int> > > >  bdryCellIds_;
  //std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<int> > > >  bdryCellLocIds_;
  std::vector<std::vector<std::vector<int> > >  bdryCellLocIds_;
  std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;

  // Finite element vectors and matrices for PDE.
  Teuchos::RCP<Tpetra::MultiVector<> > pde_vecR_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> > pde_vecJ3_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> > pde_vecH13_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> > pde_vecH23_overlap_;

  // Finite element vectors and matrices for QoI.
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecG1_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecG2_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH11_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH12_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH13_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH21_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH22_overlap_;
  Teuchos::RCP<Tpetra::MultiVector<> > qoi_vecH23_overlap_;

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
    bdryCellIds_    = meshMgr_->getSideSets((verbose_ && myRank_==0), outStream);
    int numSideSets = bdryCellIds_->size();
    if (numSideSets > 0) {
      bdryCellNodes_.resize(numSideSets);
      bdryCellLocIds_.resize(numSideSets);
      for (int i=0; i<numSideSets; ++i) {
        int numLocSides = (*bdryCellIds_)[i].size();
        bdryCellNodes_[i].resize(numLocSides);
        bdryCellLocIds_[i].resize(numLocSides);
        for (int j=0; j<numLocSides; ++j) {
          if ((*bdryCellIds_)[i][j].size() > 0) {
            int numCellsSide = (*bdryCellIds_)[i][j].size();
            for (int k=0; k<numCellsSide; ++k) {
              int idx = mapGlobalToLocalCellId((*bdryCellIds_)[i][j][k]);
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
      int lfs = dofMgr_->getLocalFieldSize();
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
      int lfs = dofMgr_->getLocalFieldSize();
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

public:
  // destructor
  virtual ~Assembler() {}

  // Constuctor: Discretization set from ParameterList
  Assembler(const std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > &basisPtrs,
          const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
          Teuchos::ParameterList &parlist,
          std::ostream &outStream = std::cout)
    : isJ1Transposed_(false), isJ2Transposed_(false) {
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
    : isJ1Transposed_(false), isJ2Transposed_(false) {
    setCommunicator(comm,parlist,outStream);
    setBasis(basisPtrs,parlist,outStream);
    setDiscretization(parlist,meshMgr,outStream);
    setParallelStructure(parlist,outStream);
    setCellNodes(outStream);
  }

  void setCellNodes(PDE<Real> &pde) const {
    // Set PDE cell nodes
    pde.setFieldPattern(dofMgr_->getFieldPattern());
    pde.setCellNodes(volCellNodes_, bdryCellNodes_, bdryCellLocIds_);
  }

  /***************************************************************************/
  /* PDE assembly routines                                                   */
  /***************************************************************************/
  void assemblePDEResidual(Teuchos::RCP<Tpetra::MultiVector<> > &r,
                           const Teuchos::RCP<PDE<Real> > &pde,
                           const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                           const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                           const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEResidual);
    #endif
    // Initialize residual vectors if not done so
    if ( r == Teuchos::null ) { // Unique components of residual vector
      r = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueResidualMap_, 1, true));
    }
    if ( pde_vecR_overlap_ == Teuchos::null ) { // Overlapping components of residual vector
      pde_vecR_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapResidualMap_, 1, true));
    }
    // Set residual vectors to zero
    r->scale(static_cast<Real>(0));
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
                                              (*res)(i,j));
      }
    }
    // change map
    Tpetra::Export<> exporter(pde_vecR_overlap_->getMap(), r->getMap()); // redistribution
    r->doExport(*pde_vecR_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
  }

  void assemblePDEJacobian1(Teuchos::RCP<Tpetra::CrsMatrix<> > &J1,
                            const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEJacobian1);
    #endif
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
      if ( J1 == Teuchos::null ) {
        J1 = Teuchos::rcp(new Tpetra::CrsMatrix<>(matJ1Graph_));
      }
      // Zero PDE Jacobian
      J1->resumeFill(); J1->setAllToScalar(static_cast<Real>(0));
      // Assemble PDE Jacobian
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> jacArrayRCP = jac->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          J1->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                  cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                  jacArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      J1->fillComplete();
      isJ1Transposed_ = false;
    }
    catch ( Exception::Zero & ez ) {
      throw Exception::Zero(">>> (Assembler::assemblePDEJacobian1): Jacobian is zero.");
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEJacobian1): Jacobian not implemented.");
    }
  }

  void assemblePDEJacobian2(Teuchos::RCP<Tpetra::CrsMatrix<> > &J2,
                            const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEJacobian2);
    #endif
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
      if ( J2 == Teuchos::null ) {
        J2 = Teuchos::rcp(new Tpetra::CrsMatrix<>(matJ2Graph_));
      }
      // Zero PDE Jacobian
      J2->resumeFill(); J2->setAllToScalar(static_cast<Real>(0));
      // Assemble Jacobian
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> jacArrayRCP = jac->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          J2->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                  cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                  jacArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      J2->fillComplete();
      isJ2Transposed_ = false;
    }
    catch ( Exception::Zero & ez ) {
      throw Exception::Zero(">>> (Assembler::assemblePDEJacobian2): Jacobian is zero.");
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEJacobian2): Jacobian not implemented.");
    }
  }

  void assemblePDEJacobian3(Teuchos::RCP<Tpetra::MultiVector<> > &J3,
                            const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEJacobian3);
    #endif
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
        if (J3 == Teuchos::null) {
          J3 = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueResidualMap_, size, true));
        }
        if ( pde_vecJ3_overlap_ == Teuchos::null) {
          pde_vecJ3_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapResidualMap_, size, true));
        }
        // Assemble PDE Jacobian wrt parametric controls
        J3->scale(static_cast<Real>(0));
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
          Tpetra::Export<> exporter(pde_vecJ3_overlap_->getMap(), J3->getMap()); // redistribution
          J3->doExport(*pde_vecJ3_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
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

  void assemblePDEHessian11(Teuchos::RCP<Tpetra::CrsMatrix<> > &H11,
                            const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian11);
    #endif
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
      if ( H11 == Teuchos::null ) {
        H11 = Teuchos::rcp(new Tpetra::CrsMatrix<>(matH11Graph_));
      }
      // Zero Hessian
      H11->resumeFill(); H11->setAllToScalar(static_cast<Real>(0));
      // Assemble PDE Hessian
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> hessArrayRCP = hess->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          H11->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                   cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                   hessArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      H11->fillComplete();
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assemblePDEHessian11): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian11): Hessian not implemented.");
    }
  }

  void assemblePDEHessian12(Teuchos::RCP<Tpetra::CrsMatrix<> > &H12,
                            const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian12);
    #endif
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
      if ( H12 == Teuchos::null ) {
        H12 = Teuchos::rcp(new Tpetra::CrsMatrix<>(matH12Graph_));
      }
      // Zero Hessian
      H12->resumeFill(); H12->setAllToScalar(static_cast<Real>(0));
      // Assemble PDE Hessian
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> hessArrayRCP = hess->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          H12->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                           cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                           hessArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      H12->fillComplete();
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assemblePDEHessian12): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian12): Hessian not implemented.");
    }
  }

  void assemblePDEHessian13(Teuchos::RCP<Tpetra::MultiVector<> > &H13,
                            const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian13);
    #endif
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
        if (H13 == Teuchos::null) {
          H13 = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueStateMap_, size, true));
        }
        if ( pde_vecH13_overlap_ == Teuchos::null) {
          pde_vecH13_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapStateMap_, size, true));
        }
        // Assemble PDE Jacobian wrt parametric controls
        H13->scale(static_cast<Real>(0));
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
          Tpetra::Export<> exporter(pde_vecH13_overlap_->getMap(), H13->getMap()); // redistribution
          H13->doExport(*pde_vecH13_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
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

  void assemblePDEHessian21(Teuchos::RCP<Tpetra::CrsMatrix<> > &H21,
                            const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian21);
    #endif
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
      if ( H21 == Teuchos::null ) {
        H21 = Teuchos::rcp(new Tpetra::CrsMatrix<>(matH21Graph_));
      }
      // Zero Hessian
      H21->resumeFill(); H21->setAllToScalar(static_cast<Real>(0));
      // Assemble PDE Hessian
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> hessArrayRCP = hess->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          H21->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                   cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                   hessArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      H21->fillComplete();
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assemblePDEHessian21): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian21): Hessian not implemented.");
    }
  }

  void assemblePDEHessian22(Teuchos::RCP<Tpetra::CrsMatrix<> > &H22,
                            const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian22);
    #endif
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
      if ( H22 == Teuchos::null ) {
        H22 = Teuchos::rcp(new Tpetra::CrsMatrix<>(matH22Graph_));
      }
      // Zero Hessian
      H22->resumeFill(); H22->setAllToScalar(static_cast<Real>(0));
      // Assemble PDE Hessian
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> hessArrayRCP = hess->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          H22->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                   cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                   hessArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      H22->fillComplete();
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assemblePDEHessian22): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assemblePDEHessian22): Hessian not implemented.");
    }
  }

  void assemblePDEHessian23(Teuchos::RCP<Tpetra::MultiVector<> > &H23,
                            const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian23);
    #endif
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
        if (H23 == Teuchos::null) {
          H23 = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueControlMap_, size, true));
        }
        if ( pde_vecH23_overlap_ == Teuchos::null) {
          pde_vecH23_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapControlMap_, size, true));
        }
        // Assemble PDE Jacobian wrt parametric controls
        H23->scale(static_cast<Real>(0));
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
          Tpetra::Export<> exporter(pde_vecH23_overlap_->getMap(), H23->getMap()); // redistribution
          H23->doExport(*pde_vecH23_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
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

  void assemblePDEHessian31(Teuchos::RCP<Tpetra::MultiVector<> > &H31,
                            const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian31);
    #endif
    assemblePDEHessian13(H31,pde,l,u,z,z_param);
  }

  void assemblePDEHessian32(Teuchos::RCP<Tpetra::MultiVector<> > &H32,
                            const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian32);
    #endif
    assemblePDEHessian23(H32,pde,l,u,z,z_param);
  }

  void assemblePDEHessian33(Teuchos::RCP<std::vector<std::vector<Real> > > &H33,
                            const Teuchos::RCP<PDE<Real> > &pde,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &l,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian33);
    #endif
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
        if (H33 == Teuchos::null) {
          std::vector<Real> col(size,static_cast<Real>(0));
          H33 = Teuchos::rcp(new std::vector<std::vector<Real> >(size,col));
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
            (*H33)[k][j] += globHess[cnt];
            if ( j != k ) { 
              (*H33)[j][k] += globHess[cnt];
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
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIValue);
    #endif
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

  void assembleQoIGradient1(Teuchos::RCP<Tpetra::MultiVector<> > &g1,
                            const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIGradient1);
    #endif
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
      if ( g1 == Teuchos::null ) {
        g1 = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueStateMap_, 1, true));
      }
      if ( qoi_vecG1_overlap_ == Teuchos::null ) {
        qoi_vecG1_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapStateMap_, 1, true));
      }
      // Assembly in to the overlap gradient
      g1->scale(static_cast<Real>(0));
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
      Tpetra::Export<> exporter(qoi_vecG1_overlap_->getMap(), g1->getMap()); // redistribution
      g1->doExport(*qoi_vecG1_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch ( Exception::Zero & ez ) {
      throw Exception::Zero(">>> (Assembler::assembleQoIGradient1): Gradient is zero.");
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIGradient1): Gradient not implemented.");
    }
  }

  void assembleQoIGradient2(Teuchos::RCP<Tpetra::MultiVector<> > &g2,
                            const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIGradient2);
    #endif
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
      if ( g2 == Teuchos::null ) {
        g2 = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueControlMap_, 1, true));
      }
      if ( qoi_vecG2_overlap_ == Teuchos::null ) {
        qoi_vecG2_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapControlMap_, 1, true));
      }
      // Assembly in to the overlap gradient
      g2->scale(static_cast<Real>(0));
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
      Tpetra::Export<> exporter(qoi_vecG2_overlap_->getMap(), g2->getMap()); // redistribution
      g2->doExport(*qoi_vecG2_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch ( Exception::Zero & ez ) {
      throw Exception::Zero(">>> (Assembler::assembleQoIGradient2): Gradient is zero.");
    }
    catch ( Exception::NotImplemented & eni ) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIGradient2): Gradient not implemented.");
    }
  }

  void assembleQoIGradient3(Teuchos::RCP<std::vector<Real> > &g3,
                            const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIGradient3);
    #endif
    if ( z_param != Teuchos::null ) {
      const int size = z_param->size();
      if ( g3 == Teuchos::null ) {
        g3 = Teuchos::rcp(new std::vector<Real>(size,0));
      }
      try {
        g3->assign(size,0);
        // Initialize local gradient storage
        std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > locGrad(size,Teuchos::null);
        // Get u_coeff from u and z_coeff from z
        Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
        getCoeffFromStateVector(u_coeff,u);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
        getCoeffFromControlVector(z_coeff,z);
        // Compute gradient
        (*g3) = qoi->gradient_3(locGrad,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
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
          (*g3)[j] += globGrad[j];
        }
      }
      catch ( Exception::Zero & ez ) {
        g3->assign(size,0);
      }
      catch ( Exception::NotImplemented & eni ) {
        throw Exception::NotImplemented(">>> (Assembler::assembleQoIGradient3): Gradient not implemented.");
      }
    }
    else {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIGradient3): Gradient not implemented.");
    }
  }

  void assembleQoIHessVec11(Teuchos::RCP<Tpetra::MultiVector<> > &H11,
                            const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec11);
    #endif
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
      if ( H11 == Teuchos::null ) {
        H11 = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueStateMap_, 1, true));
      }
      if ( qoi_vecH11_overlap_ == Teuchos::null ) {
        qoi_vecH11_overlap_  = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapStateMap_, 1, true));
      }
      // Assembly in to the overlap gradient
      H11->scale(static_cast<Real>(0));
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
      Tpetra::Export<> exporter(qoi_vecH11_overlap_->getMap(), H11->getMap()); // redistribution
      H11->doExport(*qoi_vecH11_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assembleQoIHessVec11): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec11): Hessian not implemented.");
    }
  }

  void assembleQoIHessVec12(Teuchos::RCP<Tpetra::MultiVector<> > &H12,
                            const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec12);
    #endif
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
      if ( H12 == Teuchos::null ) {
        H12 = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueStateMap_, 1, true));
      }
      if ( qoi_vecH12_overlap_ == Teuchos::null ) {
        qoi_vecH12_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapStateMap_, 1, true));
      }
      // Assembly in to the overlap gradient
      H12->scale(static_cast<Real>(0));
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
      Tpetra::Export<> exporter(qoi_vecH12_overlap_->getMap(), H12->getMap()); // redistribution
      H12->doExport(*qoi_vecH12_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assembleQoIHessVec12): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec12): Hessian not implemented.");
    }
  }

  void assembleQoIHessVec13(Teuchos::RCP<Tpetra::MultiVector<> > &H13,
                            const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const std::vector<Real> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > &z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec13);
    #endif
    if (z_param != Teuchos::null) {
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
        if ( H13 == Teuchos::null ) {
          H13 = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueStateMap_, 1, true));
        }
        if ( qoi_vecH13_overlap_ == Teuchos::null ) {
          qoi_vecH13_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapStateMap_, 1, true));
        }
        // Assembly in to the overlap gradient
        H13->scale(static_cast<Real>(0));
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
        Tpetra::Export<> exporter(qoi_vecH13_overlap_->getMap(), H13->getMap()); // redistribution
        H13->doExport(*qoi_vecH13_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
      }
      catch (Exception::Zero &ez) {
        throw Exception::Zero(">>> (Assembler::assembleQoIHessVec13): Hessian is zero.");
      }
      catch (Exception::NotImplemented &eni) {
        throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec13): Hessian not implemented.");
      }
    }
    else {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec13): Hessian not implemented.");
    }
  }

  void assembleQoIHessVec21(Teuchos::RCP<Tpetra::MultiVector<> > &H21,
                            const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec21);
    #endif
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
      if ( H21 == Teuchos::null ) {
        H21 = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueControlMap_, 1, true));
      }
      if ( qoi_vecH21_overlap_ == Teuchos::null ) {
        qoi_vecH21_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapControlMap_, 1, true));
      }
      // Assembly in to the overlap gradient
      H21->scale(static_cast<Real>(0));
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
      Tpetra::Export<> exporter(qoi_vecH21_overlap_->getMap(), H21->getMap()); // redistribution
      H21->doExport(*qoi_vecH21_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assembleQoIHessVec21): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec21): Hessian not implemented.");
    }
  }

  void assembleQoIHessVec22(Teuchos::RCP<Tpetra::MultiVector<> > &H22,
                            const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec22);
    #endif
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
      if ( H22 == Teuchos::null ) {
        H22 = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueControlMap_, 1, true));
      }
      if ( qoi_vecH22_overlap_ == Teuchos::null ) {
        qoi_vecH22_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapControlMap_, 1, true));
      }
      // Assembly in to the overlap gradient
      H22->scale(static_cast<Real>(0));
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
      Tpetra::Export<> exporter(qoi_vecH22_overlap_->getMap(), H22->getMap()); // redistribution
      H22->doExport(*qoi_vecH22_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    }
    catch (Exception::Zero &ez) {
      throw Exception::Zero(">>> (Assembler::assembleQoIHessVec22): Hessian is zero.");
    }
    catch (Exception::NotImplemented &eni) {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec22): Hessian not implemented.");
    }
  }

  void assembleQoIHessVec23(Teuchos::RCP<Tpetra::MultiVector<> > &H23,
                            const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const std::vector<Real> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec23);
    #endif
    if (z_param != Teuchos::null) {
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
        if ( H23 == Teuchos::null ) {
          H23 = Teuchos::rcp(new Tpetra::MultiVector<>(myUniqueControlMap_, 1, true));
        }
        if ( qoi_vecH23_overlap_ == Teuchos::null ) {
          qoi_vecH23_overlap_ = Teuchos::rcp(new Tpetra::MultiVector<>(myOverlapControlMap_, 1, true));
        }
        // Assembly in to the overlap gradient
        H23->scale(static_cast<Real>(0));
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
        Tpetra::Export<> exporter(qoi_vecH23_overlap_->getMap(), H23->getMap()); // redistribution
        H23->doExport(*qoi_vecH23_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
      }
      catch (Exception::Zero &ez) {
        throw Exception::Zero(">>> (Assembler::assembleQoIHessVec23): Hessian is zero.");
      }
      catch (Exception::NotImplemented &eni) {
        throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec23): Hessian not implemented.");
      }
    }
    else {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec23): Hessian not implemented.");
    }

  }

  void assembleQoIHessVec31(Teuchos::RCP<std::vector<Real> > &H31,
                            const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec31);
    #endif
    if ( z_param != Teuchos::null ) {
      const int size = z_param->size();
      if ( H31 == Teuchos::null ) {
        H31 = Teuchos::rcp(new std::vector<Real>(size,0));
      }
      try {
        H31->assign(size,0);
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
        (*H31) = qoi->HessVec_31(locHess,v_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
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
          (*H31)[j] += globHess[j];
        }
      }
      catch ( Exception::Zero & ez ) {
        H31->assign(size,0);
      }
      catch ( Exception::NotImplemented & eni ) {
        throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec31): HessVec not implemented.");
      }
    }
    else {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec31): HessVec not implemented.");
    }
  }

  void assembleQoIHessVec32(Teuchos::RCP<std::vector<Real> > &H32,
                            const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > & z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec32);
    #endif
    if ( z_param != Teuchos::null ) {
      const int size = z_param->size();
      if ( H32 == Teuchos::null ) {
        H32 = Teuchos::rcp(new std::vector<Real>(size,0));
      }
      try {
        H32->assign(size,0);
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
        (*H32) = qoi->HessVec_32(locHess,v_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
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
          (*H32)[j] += globHess[j];
        }
      }
      catch ( Exception::Zero & ez ) {
        H32->assign(size,0);
      }
      catch ( Exception::NotImplemented & eni ) {
        throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec32): HessVec not implemented.");
      }
    }
    else {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec32): HessVec not implemented.");
    }
  }

  void assembleQoIHessVec33(Teuchos::RCP<std::vector<Real> > &H33,
                            const Teuchos::RCP<QoI<Real> > &qoi,
                            const Teuchos::RCP<const std::vector<Real> > &v,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                            const Teuchos::RCP<const Tpetra::MultiVector<> > &z = Teuchos::null,
                            const Teuchos::RCP<const std::vector<Real> > &z_param = Teuchos::null) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec33);
    #endif
    if ( z_param != Teuchos::null ) {
      const int size = z_param->size();
      if ( H33 == Teuchos::null ) {
        H33 = Teuchos::rcp(new std::vector<Real>(size,0));
      }
      try {
        H33->assign(size,0);
        // Initialize local hessian times a vector storage
        std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > locHess(size,Teuchos::null);
        // Get u_coeff from u and z_coeff from z
        Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff = Teuchos::null;
        getCoeffFromStateVector(u_coeff,u);
        Teuchos::RCP<Intrepid::FieldContainer<Real> > z_coeff = Teuchos::null;
        getCoeffFromControlVector(z_coeff,z);
        // Compute local hessian times a vector
        (*H33) = qoi->HessVec_33(locHess,v,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
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
          (*H33)[j] += globHess[j];
        }
      }
      catch ( Exception::Zero & ez ) {
        H33->assign(size,0);
      }
      catch ( Exception::NotImplemented & eni ) {
        throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec33): HessVec not implemented.");
      }
    }
    else {
      throw Exception::NotImplemented(">>> (Assembler::assembleQoIHessVec33): HessVec not implemented.");
    }
  }
  /***************************************************************************/
  /* End QoI assembly routines                                               */
  /***************************************************************************/


  /***************************************************************************/
  /* Assemble and apply Riesz operator corresponding to simulation variables */
  /***************************************************************************/
  void assemblePDERieszMap1(Teuchos::RCP<Tpetra::CrsMatrix<> > &R1,
                            const Teuchos::RCP<PDE<Real> > &pde) {
    try {
      // Compute local state Riesz matrix
      Teuchos::RCP<Intrepid::FieldContainer<Real> > riesz;
      pde->RieszMap_1(riesz); // Throw if not implemented or zero
      // Initialize Riesz matrix if not done so already
      if ( R1 == Teuchos::null ) {
      R1 = Teuchos::rcp(new Tpetra::CrsMatrix<>(matR1Graph_));
      }
      R1->resumeFill(); R1->setAllToScalar(static_cast<Real>(0));
      // Assemble Riesz matrix
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> rieszArrayRCP = riesz->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          R1->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                      cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                      rieszArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      R1->fillComplete();
    }
    catch ( Exception::NotImplemented & eni ) {
      //throw Exception::NotImplemented(">>> (Assembler::assemblePDERieszMap1): Riesz map not implemented!");
    }
    catch ( Exception::Zero & ez ) {
      //throw Exception::Zero(">>> (Assembler::assemblePDERieszMap1): Riesz map is zero!");
    }
  }
  /***************************************************************************/
  /* End of functions for Riesz operator of simulation variables.            */
  /***************************************************************************/


  /***************************************************************************/
  /* Assemble and apply Riesz operator corresponding to optimization         */
  /* variables                                                               */
  /***************************************************************************/
  void assemblePDERieszMap2(Teuchos::RCP<Tpetra::CrsMatrix<> > &R2,
                            const Teuchos::RCP<PDE<Real> > &pde) {
    try {
      // Compute local control Riesz matrix
      Teuchos::RCP<Intrepid::FieldContainer<Real> > riesz;
      pde->RieszMap_2(riesz); // Throw if not implemented or zero
      // Initialize Riesz matrix if not done so already
      if ( R2 == Teuchos::null ) {
        R2 = Teuchos::rcp(new Tpetra::CrsMatrix<>(matR2Graph_));
      }
      // Assemble Riesz matrix
      R2->resumeFill(); R2->setAllToScalar(static_cast<Real>(0));
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int numLocalDofs = cellDofs.dimension(1);
      int numLocalMatEntries = numLocalDofs * numLocalDofs;
      Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
      Teuchos::ArrayRCP<const Real> rieszArrayRCP = riesz->getData();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<numLocalDofs; ++j) {
          R2->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                      cellDofsArrayRCP(myCellIds_[i] * numLocalDofs, numLocalDofs),
                                      rieszArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        }
      }
      R2->fillComplete();
    }
    catch ( Exception::NotImplemented & eni ) {
      //throw Exception::NotImplemented(">>> (Assembler::assemblePDERieszMap2): Riesz map not implemented!");
    }
    catch ( Exception::Zero & ez ) {
      //throw Exception::Zero(">>> (Assembler::assemblePDERieszMap2): Riesz map is zero!");
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
    }
  } // prinf function end

  void outputTpetraVector(const Teuchos::RCP<const Tpetra::MultiVector<> > &vec,
                          const std::string &filename) const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> > vecWriter;
    vecWriter.writeDenseFile(filename, vec);
  }

  void serialPrintStateEdgeField(const Teuchos::RCP<const Tpetra::MultiVector<> > &u,
                                 const Teuchos::RCP<FieldHelper<Real> > &fieldHelper,
                                 const std::string &filename,
                                 const Teuchos::RCP<FE_CURL<Real> > &fe) const {
    const int c = fe->curlN()->dimension(0);
    const int f = fe->curlN()->dimension(1);
    const int p = 1, d = 3;

    Teuchos::RCP<Intrepid::FieldContainer<Real> > u_coeff;
    getCoeffFromStateVector(u_coeff,u);

    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > U;
    fieldHelper->splitFieldCoeff(U, u_coeff);
    int numFields = U.size();

    // Transform cell center to physical
    Teuchos::RCP<Intrepid::FieldContainer<Real> > rx
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(p,d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > px
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p,d));
    fe->mapRefPointsToPhysical(px,rx);
    // Transform reference values into physical space.
    Teuchos::RCP<Intrepid::FieldContainer<Real> > cellJac
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p,d,d));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > cellJacInv
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p,d,d));
    Teuchos::RCP<shards::CellTopology> cellTopo
      = Teuchos::rcp(new shards::CellTopology(basisPtrs_[0]->getBaseCellTopology()));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valReference
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(f,p,d));
    basisPtrs_[0]->getValues(*valReference,*rx,Intrepid::OPERATOR_VALUE);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valPhysical
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,f,p,d));
    Intrepid::CellTools<Real>::setJacobian(*cellJac,*rx,*volCellNodes_,*cellTopo);
    Intrepid::CellTools<Real>::setJacobianInv(*cellJacInv, *cellJac);
    Intrepid::FunctionSpaceTools::HCURLtransformVALUE<Real>(*valPhysical,
                                                            *cellJacInv,
                                                            *valReference);

    std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > uval(numFields);
    for (int k = 0; k < numFields; ++k) {
      uval[k] = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c,p,d));
      Intrepid::FunctionSpaceTools::evaluate<Real>(*uval[k], *U[k], *valPhysical);
    }
    // Print
    std::ofstream file;
    file.open(filename);
    for (int i = 0; i < c; ++i) {
      file << std::scientific << std::setprecision(15);
      for (int j = 0; j < d; ++j) {
        file << std::setw(25) << (*px)(i,0,j);
      }
      for (int k = 0; k < numFields; ++k) {
        for (int j = 0; j < d; ++j) {
          file << std::setw(25) << (*uval[k])(i,0,j);
        }
      }
      file << std::endl;
    }
    file.close();
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

  /***************************************************************************/
  /* Accessor routines.                                                      */
  /***************************************************************************/
  const Teuchos::RCP<DofManager<Real> > getDofManager(void) const {
    return dofMgr_;
  }

  Teuchos::Array<int> getCellIds(void) const {
    return myCellIds_;
  }
  /***************************************************************************/
  /* End of accessor routines.                                               */
  /***************************************************************************/


}; // class Assembler

#endif
