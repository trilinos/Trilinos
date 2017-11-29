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

#include "Intrepid_DefaultCubatureFactory.hpp"

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
    ROL::SharedPointer<Teuchos::Time> AssemblePDEResidual       = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Residual");
    ROL::SharedPointer<Teuchos::Time> AssemblePDEJacobian1      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Jacobian1");
    ROL::SharedPointer<Teuchos::Time> AssemblePDEJacobian2      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Jacobian2");
    ROL::SharedPointer<Teuchos::Time> AssemblePDEJacobian3      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Jacobian3");
    ROL::SharedPointer<Teuchos::Time> AssemblePDEHessian11      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian11");
    ROL::SharedPointer<Teuchos::Time> AssemblePDEHessian12      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian12");
    ROL::SharedPointer<Teuchos::Time> AssemblePDEHessian13      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian13");
    ROL::SharedPointer<Teuchos::Time> AssemblePDEHessian21      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian21");
    ROL::SharedPointer<Teuchos::Time> AssemblePDEHessian22      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian22");
    ROL::SharedPointer<Teuchos::Time> AssemblePDEHessian23      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian23");
    ROL::SharedPointer<Teuchos::Time> AssemblePDEHessian31      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian31");
    ROL::SharedPointer<Teuchos::Time> AssemblePDEHessian32      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian32");
    ROL::SharedPointer<Teuchos::Time> AssemblePDEHessian33      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian33");
    ROL::SharedPointer<Teuchos::Time> AssembleQOIValue          = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Value");
    ROL::SharedPointer<Teuchos::Time> AssembleQOIGradient1      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Gradient1");
    ROL::SharedPointer<Teuchos::Time> AssembleQOIGradient2      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Gradient2");
    ROL::SharedPointer<Teuchos::Time> AssembleQOIGradient3      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Gradient3");
    ROL::SharedPointer<Teuchos::Time> AssembleQOIHessVec11      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec11");
    ROL::SharedPointer<Teuchos::Time> AssembleQOIHessVec12      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec12");
    ROL::SharedPointer<Teuchos::Time> AssembleQOIHessVec13      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec13");
    ROL::SharedPointer<Teuchos::Time> AssembleQOIHessVec21      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec21");
    ROL::SharedPointer<Teuchos::Time> AssembleQOIHessVec22      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec22");
    ROL::SharedPointer<Teuchos::Time> AssembleQOIHessVec23      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec23");
    ROL::SharedPointer<Teuchos::Time> AssembleQOIHessVec31      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec31");
    ROL::SharedPointer<Teuchos::Time> AssembleQOIHessVec32      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec32");
    ROL::SharedPointer<Teuchos::Time> AssembleQOIHessVec33      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec33");
  }
}
#endif

template<class Real>
class Solution {
public:
  virtual ~Solution() {}
  Solution(void) {}
  virtual Real evaluate(const std::vector<Real> &x, const int fieldNumber) const = 0;
};

template<class Real>
class Assembler {

  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Map<>::node_type NO;
  typedef Tpetra::MultiVector<Real,LO,GO,NO> MV;
  typedef Tpetra::Operator<Real,LO,GO,NO> OP;

private:
  // Timers
//  ROL::SharedPointer<Teuchos::Time::Time> timerSolverFactorization_;
//  ROL::SharedPointer<Teuchos::Time::Time> timerSolverSubstitution_;
//  ROL::SharedPointer<Teuchos::Time::Time> timerAssemblyNonlinear_;
//  ROL::SharedPointer<Teuchos::Time::Time> timerSolverUpdate_;

  // Set in Constructor.
  bool verbose_;
  bool isJ1Transposed_, isJ2Transposed_;

  // Set in SetCommunicator.
  ROL::SharedPointer<const Teuchos::Comm<int> > comm_;
  int myRank_, numProcs_;

  // Set in SetBasis.
  std::vector<ROL::SharedPointer<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;

  // Set in SetDiscretization.
  ROL::SharedPointer<MeshManager<Real> > meshMgr_;
  ROL::SharedPointer<DofManager<Real> >  dofMgr_;

  // Set in SetParallelStructure.
  int numCells_;
  Teuchos::Array<int> myCellIds_;
  Teuchos::Array<int> cellOffsets_;
  ROL::SharedPointer<const Tpetra::Map<> > myOverlapStateMap_;
  ROL::SharedPointer<const Tpetra::Map<> > myUniqueStateMap_;
  ROL::SharedPointer<const Tpetra::Map<> > myOverlapControlMap_;
  ROL::SharedPointer<const Tpetra::Map<> > myUniqueControlMap_;
  ROL::SharedPointer<const Tpetra::Map<> > myOverlapResidualMap_;
  ROL::SharedPointer<const Tpetra::Map<> > myUniqueResidualMap_;
  ROL::SharedPointer<Tpetra::CrsGraph<> >  matJ1Graph_;
  ROL::SharedPointer<Tpetra::CrsGraph<> >  matJ2Graph_;
  ROL::SharedPointer<Tpetra::CrsGraph<> >  matR1Graph_;
  ROL::SharedPointer<Tpetra::CrsGraph<> >  matR2Graph_;
  ROL::SharedPointer<Tpetra::CrsGraph<> >  matH11Graph_;
  ROL::SharedPointer<Tpetra::CrsGraph<> >  matH12Graph_;
  ROL::SharedPointer<Tpetra::CrsGraph<> >  matH21Graph_;
  ROL::SharedPointer<Tpetra::CrsGraph<> >  matH22Graph_;

  // Set in SetCellNodes.
  ROL::SharedPointer<Intrepid::FieldContainer<Real> > volCellNodes_;
  ROL::SharedPointer<std::vector<std::vector<std::vector<int> > > >  bdryCellIds_;
  //std::vector<std::vector<ROL::SharedPointer<Intrepid::FieldContainer<int> > > >  bdryCellLocIds_;
  std::vector<std::vector<std::vector<int> > >  bdryCellLocIds_;
  std::vector<std::vector<ROL::SharedPointer<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;

  // Finite element vectors and matrices for PDE.
  ROL::SharedPointer<Tpetra::MultiVector<> > pde_vecR_overlap_;
  ROL::SharedPointer<Tpetra::MultiVector<> > pde_vecJ3_overlap_;
  ROL::SharedPointer<Tpetra::MultiVector<> > pde_vecH13_overlap_;
  ROL::SharedPointer<Tpetra::MultiVector<> > pde_vecH23_overlap_;

  // Finite element vectors and matrices for QoI.
  ROL::SharedPointer<Tpetra::MultiVector<> > qoi_vecG1_overlap_;
  ROL::SharedPointer<Tpetra::MultiVector<> > qoi_vecG2_overlap_;
  ROL::SharedPointer<Tpetra::MultiVector<> > qoi_vecH11_overlap_;
  ROL::SharedPointer<Tpetra::MultiVector<> > qoi_vecH12_overlap_;
  ROL::SharedPointer<Tpetra::MultiVector<> > qoi_vecH13_overlap_;
  ROL::SharedPointer<Tpetra::MultiVector<> > qoi_vecH21_overlap_;
  ROL::SharedPointer<Tpetra::MultiVector<> > qoi_vecH22_overlap_;
  ROL::SharedPointer<Tpetra::MultiVector<> > qoi_vecH23_overlap_;

private:

  void setCommunicator(const ROL::SharedPointer<const Teuchos::Comm<int> > &comm,
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
         const std::vector<ROL::SharedPointer<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > &basisPtrs,
         Teuchos::ParameterList &parlist,
         std::ostream &outStream = std::cout) {
    basisPtrs_ = basisPtrs;
    if (verbose_ && myRank_==0) {
      outStream << "Initialized PDE." << std::endl;
    }
  }

  void setDiscretization(Teuchos::ParameterList &parlist,
                         const ROL::SharedPointer<MeshManager<Real> > &meshMgr = ROL::nullPointer,
                         std::ostream &outStream = std::cout) {
    if ( meshMgr != ROL::nullPointer ) {
      // Use MeshManager object if supplied
      meshMgr_ = meshMgr;
    }
    else {
      // Otherwise construct MeshManager objective from parameter list
    }
    dofMgr_ = ROL::makeShared<DofManager<Real>>(meshMgr_,basisPtrs_);
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
    myOverlapStateMap_ = ROL::makeShared<Tpetra::Map<>>(
                         Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                         myGlobalIds, 0, comm_);
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
//    myCellMap_ = ROL::makeShared<Tpetra::Map<>>(
//                 Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
//                 myCellIds_, 0, comm_);

    /****************************************/
    /*** Assemble global graph structure. ***/
    /****************************************/
    matJ1Graph_ = ROL::makeShared<Tpetra::CrsGraph<>>(myUniqueStateMap_, 0);
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
    volCellNodes_ = ROL::makeShared<Intrepid::FieldContainer<Real>>(numCells_, numNodesPerCell, spaceDim);
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
              bdryCellNodes_[i][j] = ROL::makeShared<Intrepid::FieldContainer<Real>>(myNumCellsSide, numNodesPerCell, spaceDim);
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

  void getCoeffFromStateVector(ROL::SharedPointer<Intrepid::FieldContainer<Real> > &xcoeff,
                               const ROL::SharedPointer<const Tpetra::MultiVector<> > &x) const {
    if ( x != ROL::nullPointer ) {
      // Perform import onto myOverlapMap
      ROL::SharedPointer<Tpetra::MultiVector<> > xshared =
        ROL::makeShared<Tpetra::MultiVector<>>(myOverlapStateMap_, 1, true);
      Tpetra::Import<> importer(myUniqueStateMap_, myOverlapStateMap_);
      xshared->doImport(*x,importer,Tpetra::REPLACE);
      // Populate xcoeff
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int lfs = dofMgr_->getLocalFieldSize();
      xcoeff = ROL::makeShared<Intrepid::FieldContainer<Real>>(numCells_, lfs);
      Teuchos::ArrayRCP<const Real> xdata = xshared->get1dView();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<lfs; ++j) {
          (*xcoeff)(i,j) = xdata[xshared->getMap()->getLocalElement(cellDofs(myCellIds_[i],j))];
        }
      }
      dofMgr_->transformToIntrepidPattern(xcoeff);
    }
    else {
      xcoeff = ROL::nullPointer;
    }
  }

  void getCoeffFromControlVector(ROL::SharedPointer<Intrepid::FieldContainer<Real> > &xcoeff,
                                 const ROL::SharedPointer<const Tpetra::MultiVector<> > &x) const {
    if ( x != ROL::nullPointer ) {
      // Perform import onto myOverlapMap
      ROL::SharedPointer<Tpetra::MultiVector<> > xshared =
        ROL::makeShared<Tpetra::MultiVector<>>(myOverlapControlMap_, 1, true);
      Tpetra::Import<> importer(myUniqueControlMap_, myOverlapControlMap_);
      xshared->doImport(*x,importer,Tpetra::REPLACE);
      // Populate xcoeff
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int lfs = dofMgr_->getLocalFieldSize();
      xcoeff = ROL::makeShared<Intrepid::FieldContainer<Real>>(numCells_, lfs);
      Teuchos::ArrayRCP<const Real> xdata = xshared->get1dView();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<lfs; ++j) {
          (*xcoeff)(i,j) = xdata[xshared->getMap()->getLocalElement(cellDofs(myCellIds_[i],j))];
        }
      }
      dofMgr_->transformToIntrepidPattern(xcoeff);
    }
    else {
      xcoeff = ROL::nullPointer;
    }
  }

public:
  // destructor
  virtual ~Assembler() {}

  // Constuctor: Discretization set from ParameterList
  Assembler(const std::vector<ROL::SharedPointer<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > &basisPtrs,
          const ROL::SharedPointer<const Teuchos::Comm<int> > &comm,
          Teuchos::ParameterList &parlist,
          std::ostream &outStream = std::cout)
    : isJ1Transposed_(false), isJ2Transposed_(false) {
    setCommunicator(comm,parlist,outStream);
    setBasis(basisPtrs,parlist,outStream);
    setDiscretization(parlist,ROL::nullPointer,outStream);
    setParallelStructure(parlist,outStream);
    setCellNodes(outStream);
  }

  // Constructor: Discretization set from MeshManager input
  Assembler(const std::vector<ROL::SharedPointer<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > &basisPtrs,
          const ROL::SharedPointer<MeshManager<Real> > &meshMgr,
          const ROL::SharedPointer<const Teuchos::Comm<int> > &comm,
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
  void assemblePDEResidual(ROL::SharedPointer<Tpetra::MultiVector<> > &r,
                           const ROL::SharedPointer<PDE<Real> > &pde,
                           const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                           const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                           const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEResidual);
    #endif
    // Initialize residual vectors if not done so
    if ( r == ROL::nullPointer ) { // Unique components of residual vector
      r = ROL::makeShared<Tpetra::MultiVector<>>(myUniqueResidualMap_, 1, true);
    }
    if ( pde_vecR_overlap_ == ROL::nullPointer ) { // Overlapping components of residual vector
      pde_vecR_overlap_ = ROL::makeShared<Tpetra::MultiVector<>>(myOverlapResidualMap_, 1, true);
    }
    // Set residual vectors to zero
    r->scale(static_cast<Real>(0));
    pde_vecR_overlap_->scale(static_cast<Real>(0));
    // Get degrees of freedom
    Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
    int numLocalDofs = cellDofs.dimension(1);
    // Initialize res container
    ROL::SharedPointer<Intrepid::FieldContainer<Real> > res;
    // Get u_coeff from u and z_coeff from z
    ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
    getCoeffFromStateVector(u_coeff,u);
    ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
    getCoeffFromControlVector(z_coeff,z);
    // Compute PDE residual
    pde->residual(res,u_coeff,z_coeff,z_param);
    dofMgr_->transformToFieldPattern(res);
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

  void assemblePDEJacobian1(ROL::SharedPointer<Tpetra::CrsMatrix<> > &J1,
                            const ROL::SharedPointer<PDE<Real> > &pde,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEJacobian1);
    #endif
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
      getCoeffFromStateVector(u_coeff,u);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
      getCoeffFromControlVector(z_coeff,z);
      // Compute PDE Jacobian
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > jac;
      pde->Jacobian_1(jac,u_coeff,z_coeff,z_param); // Throw if not implemented or is zero
      dofMgr_->transformToFieldPattern(jac);
      // Initialize Jacobian matrices
      if ( J1 == ROL::nullPointer ) {
        J1 = ROL::makeShared<Tpetra::CrsMatrix<>>(matJ1Graph_);
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

  void assemblePDEJacobian2(ROL::SharedPointer<Tpetra::CrsMatrix<> > &J2,
                            const ROL::SharedPointer<PDE<Real> > &pde,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEJacobian2);
    #endif
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
      getCoeffFromStateVector(u_coeff,u);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
      getCoeffFromControlVector(z_coeff,z);
      // Compute PDE Jacobian
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > jac;
      pde->Jacobian_2(jac,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(jac);
      // Initialize Jacobian matrices
      if ( J2 == ROL::nullPointer ) {
        J2 = ROL::makeShared<Tpetra::CrsMatrix<>>(matJ2Graph_);
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

  void assemblePDEJacobian3(ROL::SharedPointer<Tpetra::MultiVector<> > &J3,
                            const ROL::SharedPointer<PDE<Real> > &pde,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEJacobian3);
    #endif
    if ( z_param != ROL::nullPointer ) {
      try {
        int size = static_cast<int>(z_param->size());
        // Initialize res
        std::vector<ROL::SharedPointer<Intrepid::FieldContainer<Real> > > jac(size);
        // Get u_coeff from u and z_coeff from z
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
        getCoeffFromStateVector(u_coeff,u);
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
        getCoeffFromControlVector(z_coeff,z);
        // Compute PDE local Jacobian wrt parametric controls
        pde->Jacobian_3(jac,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        for (int i = 0; i < size; ++i) {
          dofMgr_->transformToFieldPattern(jac[i]);
        }
        // Initialize Jacobian storage if not done so already
        if (J3 == ROL::nullPointer) {
          J3 = ROL::makeShared<Tpetra::MultiVector<>>(myUniqueResidualMap_, size, true);
        }
        if ( pde_vecJ3_overlap_ == ROL::nullPointer) {
          pde_vecJ3_overlap_ = ROL::makeShared<Tpetra::MultiVector<>>(myOverlapResidualMap_, size, true);
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

  void assemblePDEHessian11(ROL::SharedPointer<Tpetra::CrsMatrix<> > &H11,
                            const ROL::SharedPointer<PDE<Real> > &pde,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &l,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian11);
    #endif
    try {
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > hess; 
      // Get u_coeff from u and z_coeff from z
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
      getCoeffFromStateVector(u_coeff,u);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
      getCoeffFromControlVector(z_coeff,z);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > l_coeff = ROL::nullPointer;
      getCoeffFromStateVector(l_coeff,l);
      // Compute PDE Hessian
      pde->Hessian_11(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(hess);
      // Initialize Hessian storage if not done so already
      if ( H11 == ROL::nullPointer ) {
        H11 = ROL::makeShared<Tpetra::CrsMatrix<>>(matH11Graph_);
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

  void assemblePDEHessian12(ROL::SharedPointer<Tpetra::CrsMatrix<> > &H12,
                            const ROL::SharedPointer<PDE<Real> > &pde,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &l,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian12);
    #endif
    try {
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > hess;
      // Get u_coeff from u and z_coeff from z
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
      getCoeffFromStateVector(u_coeff,u);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
      getCoeffFromControlVector(z_coeff,z);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > l_coeff = ROL::nullPointer;
      getCoeffFromStateVector(l_coeff,l);
      // Compute PDE Hessian
      pde->Hessian_12(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(hess);
      // Initialize Hessian storage if not done so already
      if ( H12 == ROL::nullPointer ) {
        H12 = ROL::makeShared<Tpetra::CrsMatrix<>>(matH12Graph_);
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

  void assemblePDEHessian13(ROL::SharedPointer<Tpetra::MultiVector<> > &H13,
                            const ROL::SharedPointer<PDE<Real> > &pde,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &l,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian13);
    #endif
    if ( z_param != ROL::nullPointer ) {
      try {
        int size = static_cast<int>(z_param->size());
        // Initialize local hessian
        std::vector<ROL::SharedPointer<Intrepid::FieldContainer<Real> > > hess(size);
        // Get u_coeff from u, z_coeff from z and l_coeff from l
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
        getCoeffFromStateVector(u_coeff,u);
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
        getCoeffFromControlVector(z_coeff,z);
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > l_coeff = ROL::nullPointer;
        getCoeffFromStateVector(l_coeff,l);
        // Compute PDE local Jacobian wrt parametric controls
        pde->Hessian_13(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        for (int i = 0; i < size; ++i) {
          dofMgr_->transformToFieldPattern(hess[i]);
        }
        // Initialize Jacobian storage if not done so already
        if (H13 == ROL::nullPointer) {
          H13 = ROL::makeShared<Tpetra::MultiVector<>>(myUniqueStateMap_, size, true);
        }
        if ( pde_vecH13_overlap_ == ROL::nullPointer) {
          pde_vecH13_overlap_ = ROL::makeShared<Tpetra::MultiVector<>>(myOverlapStateMap_, size, true);
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

  void assemblePDEHessian21(ROL::SharedPointer<Tpetra::CrsMatrix<> > &H21,
                            const ROL::SharedPointer<PDE<Real> > &pde,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &l,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian21);
    #endif
    try {
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > hess;
      // Get u_coeff from u and z_coeff from z
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
      getCoeffFromStateVector(u_coeff,u);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
      getCoeffFromControlVector(z_coeff,z);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > l_coeff = ROL::nullPointer;
      getCoeffFromStateVector(l_coeff,l);
      // Compute PDE Hessian
      pde->Hessian_21(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(hess);
      // Initialize Hessian storage if not done so already
      if ( H21 == ROL::nullPointer ) {
        H21 = ROL::makeShared<Tpetra::CrsMatrix<>>(matH21Graph_);
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

  void assemblePDEHessian22(ROL::SharedPointer<Tpetra::CrsMatrix<> > &H22,
                            const ROL::SharedPointer<PDE<Real> > &pde,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &l,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian22);
    #endif
    try {
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > hess;
      // Get u_coeff from u and z_coeff from z
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
      getCoeffFromStateVector(u_coeff,u);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
      getCoeffFromControlVector(z_coeff,z);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > l_coeff = ROL::nullPointer;
      getCoeffFromStateVector(l_coeff,l);
      // Compute PDE Hessian
      pde->Hessian_22(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(hess);
      // Initialize Hessian storage if not done so already
      if ( H22 == ROL::nullPointer ) {
        H22 = ROL::makeShared<Tpetra::CrsMatrix<>>(matH22Graph_);
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

  void assemblePDEHessian23(ROL::SharedPointer<Tpetra::MultiVector<> > &H23,
                            const ROL::SharedPointer<PDE<Real> > &pde,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &l,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian23);
    #endif
    if ( z_param != ROL::nullPointer ) {
      try {
        int size = static_cast<int>(z_param->size());
        // Initialize local hessian
        std::vector<ROL::SharedPointer<Intrepid::FieldContainer<Real> > > hess(size);
        // Get u_coeff from u, z_coeff from z and l_coeff from l
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
        getCoeffFromStateVector(u_coeff,u);
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
        getCoeffFromControlVector(z_coeff,z);
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > l_coeff = ROL::nullPointer;
        getCoeffFromStateVector(l_coeff,l);
        // Compute PDE local Jacobian wrt parametric controls
        pde->Hessian_23(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        for (int i = 0; i < size; ++i) {
          dofMgr_->transformToFieldPattern(hess[i]);
        }
        // Initialize Jacobian storage if not done so already
        if (H23 == ROL::nullPointer) {
          H23 = ROL::makeShared<Tpetra::MultiVector<>>(myUniqueControlMap_, size, true);
        }
        if ( pde_vecH23_overlap_ == ROL::nullPointer) {
          pde_vecH23_overlap_ = ROL::makeShared<Tpetra::MultiVector<>>(myOverlapControlMap_, size, true);
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

  void assemblePDEHessian31(ROL::SharedPointer<Tpetra::MultiVector<> > &H31,
                            const ROL::SharedPointer<PDE<Real> > &pde,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &l,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian31);
    #endif
    assemblePDEHessian13(H31,pde,l,u,z,z_param);
  }

  void assemblePDEHessian32(ROL::SharedPointer<Tpetra::MultiVector<> > &H32,
                            const ROL::SharedPointer<PDE<Real> > &pde,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &l,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian32);
    #endif
    assemblePDEHessian23(H32,pde,l,u,z,z_param);
  }

  void assemblePDEHessian33(ROL::SharedPointer<std::vector<std::vector<Real> > > &H33,
                            const ROL::SharedPointer<PDE<Real> > &pde,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &l,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian33);
    #endif
    if ( z_param != ROL::nullPointer ) {
      try {
        int size = static_cast<int>(z_param->size());
        // Initialize local hessian
        std::vector<ROL::SharedPointer<Intrepid::FieldContainer<Real> > > tmp(size,ROL::nullPointer);
        std::vector<std::vector<ROL::SharedPointer<Intrepid::FieldContainer<Real> > > > hess(size,tmp);
        // Get u_coeff from u, z_coeff from z and l_coeff from l
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
        getCoeffFromStateVector(u_coeff,u);
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
        getCoeffFromControlVector(z_coeff,z);
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > l_coeff = ROL::nullPointer;
        getCoeffFromStateVector(l_coeff,l);
        // Compute PDE local Jacobian wrt parametric controls
        pde->Hessian_33(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        for (int i = 0; i < size; ++i) {
          for (int j = 0; j < size; ++j) {
            dofMgr_->transformToFieldPattern(hess[i][j]);
          }
        }
        // Initialize Jacobian storage if not done so already
        if (H33 == ROL::nullPointer) {
          std::vector<Real> col(size,static_cast<Real>(0));
          H33 = ROL::makeShared<std::vector<std::vector<Real> >>(size,col);
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
  Real assembleQoIValue(const ROL::SharedPointer<QoI<Real> > &qoi,
                        const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                        const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                        const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIValue);
    #endif
    Real val(0);
    try {
      // Integrate obj object
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > locVal;
      // Get u_coeff from u and z_coeff from z
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
      getCoeffFromStateVector(u_coeff,u);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
      getCoeffFromControlVector(z_coeff,z);
      // Get OBJ_CELL value
      val = qoi->value(locVal,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      // Assembly
      if ( locVal != ROL::nullPointer ) {
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

  void assembleQoIGradient1(ROL::SharedPointer<Tpetra::MultiVector<> > &g1,
                            const ROL::SharedPointer<QoI<Real> > &qoi,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIGradient1);
    #endif
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
      getCoeffFromStateVector(u_coeff,u);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > locGrad;
      qoi->gradient_1(locGrad,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(locGrad);
      // Initialize state QoI gradient vectors
      if ( g1 == ROL::nullPointer ) {
        g1 = ROL::makeShared<Tpetra::MultiVector<>>(myUniqueStateMap_, 1, true);
      }
      if ( qoi_vecG1_overlap_ == ROL::nullPointer ) {
        qoi_vecG1_overlap_ = ROL::makeShared<Tpetra::MultiVector<>>(myOverlapStateMap_, 1, true);
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

  void assembleQoIGradient2(ROL::SharedPointer<Tpetra::MultiVector<> > &g2,
                            const ROL::SharedPointer<QoI<Real> > &qoi,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIGradient2);
    #endif
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
      getCoeffFromStateVector(u_coeff,u);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > locGrad;
      qoi->gradient_2(locGrad,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(locGrad);
      // Initialize control gradient vectors
      if ( g2 == ROL::nullPointer ) {
        g2 = ROL::makeShared<Tpetra::MultiVector<>>(myUniqueControlMap_, 1, true);
      }
      if ( qoi_vecG2_overlap_ == ROL::nullPointer ) {
        qoi_vecG2_overlap_ = ROL::makeShared<Tpetra::MultiVector<>>(myOverlapControlMap_, 1, true);
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

  void assembleQoIGradient3(ROL::SharedPointer<std::vector<Real> > &g3,
                            const ROL::SharedPointer<QoI<Real> > &qoi,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIGradient3);
    #endif
    if ( z_param != ROL::nullPointer ) {
      const int size = z_param->size();
      if ( g3 == ROL::nullPointer ) {
        g3 = ROL::makeShared<std::vector<Real>>(size,0);
      }
      try {
        g3->assign(size,0);
        // Initialize local gradient storage
        std::vector<ROL::SharedPointer<Intrepid::FieldContainer<Real> > > locGrad(size,ROL::nullPointer);
        // Get u_coeff from u and z_coeff from z
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
        getCoeffFromStateVector(u_coeff,u);
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
        getCoeffFromControlVector(z_coeff,z);
        // Compute gradient
        (*g3) = qoi->gradient_3(locGrad,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        for (int i = 0; i < size; ++i) {
          dofMgr_->transformToFieldPattern(locGrad[i]);
        }
        // Assembly
        std::vector<Real> myGrad(size,0), globGrad(size,0);
        for (int j = 0; j < size; ++j) {
          if ( locGrad[j] != ROL::nullPointer ) {
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

  void assembleQoIHessVec11(ROL::SharedPointer<Tpetra::MultiVector<> > &H11,
                            const ROL::SharedPointer<QoI<Real> > &qoi,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &v,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec11);
    #endif
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > locHess;
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > v_coeff = ROL::nullPointer;
      getCoeffFromStateVector(v_coeff,v);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
      getCoeffFromStateVector(u_coeff,u);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi->HessVec_11(locHess, v_coeff, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(locHess);
      // Initialize state-state HessVec vectors
      if ( H11 == ROL::nullPointer ) {
        H11 = ROL::makeShared<Tpetra::MultiVector<>>(myUniqueStateMap_, 1, true);
      }
      if ( qoi_vecH11_overlap_ == ROL::nullPointer ) {
        qoi_vecH11_overlap_  = ROL::makeShared<Tpetra::MultiVector<>>(myOverlapStateMap_, 1, true);
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

  void assembleQoIHessVec12(ROL::SharedPointer<Tpetra::MultiVector<> > &H12,
                            const ROL::SharedPointer<QoI<Real> > &qoi,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &v,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec12);
    #endif
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > locHess;
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > v_coeff = ROL::nullPointer;
      getCoeffFromControlVector(v_coeff,v);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
      getCoeffFromStateVector(u_coeff,u);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi->HessVec_12(locHess, v_coeff, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(locHess);
      // Initialize state-control HessVec vectors
      if ( H12 == ROL::nullPointer ) {
        H12 = ROL::makeShared<Tpetra::MultiVector<>>(myUniqueStateMap_, 1, true);
      }
      if ( qoi_vecH12_overlap_ == ROL::nullPointer ) {
        qoi_vecH12_overlap_ = ROL::makeShared<Tpetra::MultiVector<>>(myOverlapStateMap_, 1, true);
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

  void assembleQoIHessVec13(ROL::SharedPointer<Tpetra::MultiVector<> > &H13,
                            const ROL::SharedPointer<QoI<Real> > &qoi,
                            const ROL::SharedPointer<const std::vector<Real> > &v,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > &z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec13);
    #endif
    if (z_param != ROL::nullPointer) {
      try {
        // Get u_coeff from u and z_coeff from z
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > locHess;
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
        getCoeffFromStateVector(u_coeff,u);
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
        getCoeffFromControlVector(z_coeff,z);
        // Compute local gradient
        qoi->HessVec_13(locHess, v, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
        dofMgr_->transformToFieldPattern(locHess);
        // Initialize state-control HessVec vectors
        if ( H13 == ROL::nullPointer ) {
          H13 = ROL::makeShared<Tpetra::MultiVector<>>(myUniqueStateMap_, 1, true);
        }
        if ( qoi_vecH13_overlap_ == ROL::nullPointer ) {
          qoi_vecH13_overlap_ = ROL::makeShared<Tpetra::MultiVector<>>(myOverlapStateMap_, 1, true);
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

  void assembleQoIHessVec21(ROL::SharedPointer<Tpetra::MultiVector<> > &H21,
                            const ROL::SharedPointer<QoI<Real> > &qoi,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &v,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec21);
    #endif
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > locHess;
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > v_coeff = ROL::nullPointer;
      getCoeffFromStateVector(v_coeff,v);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
      getCoeffFromStateVector(u_coeff,u);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi->HessVec_21(locHess, v_coeff, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(locHess);
      // Initialize control-state HessVec vectors
      if ( H21 == ROL::nullPointer ) {
        H21 = ROL::makeShared<Tpetra::MultiVector<>>(myUniqueControlMap_, 1, true);
      }
      if ( qoi_vecH21_overlap_ == ROL::nullPointer ) {
        qoi_vecH21_overlap_ = ROL::makeShared<Tpetra::MultiVector<>>(myOverlapControlMap_, 1, true);
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

  void assembleQoIHessVec22(ROL::SharedPointer<Tpetra::MultiVector<> > &H22,
                            const ROL::SharedPointer<QoI<Real> > &qoi,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &v,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec22);
    #endif
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > locHess;
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > v_coeff = ROL::nullPointer;
      getCoeffFromControlVector(v_coeff,v);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
      getCoeffFromStateVector(u_coeff,u);
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi->HessVec_22(locHess, v_coeff, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(locHess);
      // Initialize control-control HessVec vectors
      if ( H22 == ROL::nullPointer ) {
        H22 = ROL::makeShared<Tpetra::MultiVector<>>(myUniqueControlMap_, 1, true);
      }
      if ( qoi_vecH22_overlap_ == ROL::nullPointer ) {
        qoi_vecH22_overlap_ = ROL::makeShared<Tpetra::MultiVector<>>(myOverlapControlMap_, 1, true);
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

  void assembleQoIHessVec23(ROL::SharedPointer<Tpetra::MultiVector<> > &H23,
                            const ROL::SharedPointer<QoI<Real> > &qoi,
                            const ROL::SharedPointer<const std::vector<Real> > &v,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec23);
    #endif
    if (z_param != ROL::nullPointer) {
      try {
        // Get u_coeff from u and z_coeff from z
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > locHess;
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
        getCoeffFromStateVector(u_coeff,u);
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
        getCoeffFromControlVector(z_coeff,z);
        // Compute local gradient
        qoi->HessVec_23(locHess, v, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
        dofMgr_->transformToFieldPattern(locHess);
        // Initialize control-control HessVec vectors
        if ( H23 == ROL::nullPointer ) {
          H23 = ROL::makeShared<Tpetra::MultiVector<>>(myUniqueControlMap_, 1, true);
        }
        if ( qoi_vecH23_overlap_ == ROL::nullPointer ) {
          qoi_vecH23_overlap_ = ROL::makeShared<Tpetra::MultiVector<>>(myOverlapControlMap_, 1, true);
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

  void assembleQoIHessVec31(ROL::SharedPointer<std::vector<Real> > &H31,
                            const ROL::SharedPointer<QoI<Real> > &qoi,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &v,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec31);
    #endif
    if ( z_param != ROL::nullPointer ) {
      const int size = z_param->size();
      if ( H31 == ROL::nullPointer ) {
        H31 = ROL::makeShared<std::vector<Real>>(size,0);
      }
      try {
        H31->assign(size,0);
        // Initialize local gradient storage
        std::vector<ROL::SharedPointer<Intrepid::FieldContainer<Real> > > locHess(size,ROL::nullPointer);
        // Get u_coeff from u and z_coeff from z
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > v_coeff = ROL::nullPointer;
        getCoeffFromStateVector(v_coeff,v);
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
        getCoeffFromStateVector(u_coeff,u);
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
        getCoeffFromControlVector(z_coeff,z);
        // Compute gradient
        (*H31) = qoi->HessVec_31(locHess,v_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        for (int i = 0; i < size; ++i) {
          dofMgr_->transformToFieldPattern(locHess[i]);
        }
        // Assembly
        std::vector<Real> myHess(size,0), globHess(size,0);
        for (int j = 0; j < size; ++j) {
          if ( locHess[j] != ROL::nullPointer ) {
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

  void assembleQoIHessVec32(ROL::SharedPointer<std::vector<Real> > &H32,
                            const ROL::SharedPointer<QoI<Real> > &qoi,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &v,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > & z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec32);
    #endif
    if ( z_param != ROL::nullPointer ) {
      const int size = z_param->size();
      if ( H32 == ROL::nullPointer ) {
        H32 = ROL::makeShared<std::vector<Real>>(size,0);
      }
      try {
        H32->assign(size,0);
        // Initialize local hessian times a vector storage
        std::vector<ROL::SharedPointer<Intrepid::FieldContainer<Real> > > locHess(size,ROL::nullPointer);
        // Get u_coeff from u and z_coeff from z
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > v_coeff = ROL::nullPointer;
        getCoeffFromControlVector(v_coeff,v);
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
        getCoeffFromStateVector(u_coeff,u);
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
        getCoeffFromControlVector(z_coeff,z);
        // Compute local hessian times a vector
        (*H32) = qoi->HessVec_32(locHess,v_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        for (int i = 0; i < size; ++i) {
          dofMgr_->transformToFieldPattern(locHess[i]);
        }
        // Assembly
        std::vector<Real> myHess(size,0), globHess(size,0);
        for (int j = 0; j < size; ++j) {
          if ( locHess[j] != ROL::nullPointer ) {
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

  void assembleQoIHessVec33(ROL::SharedPointer<std::vector<Real> > &H33,
                            const ROL::SharedPointer<QoI<Real> > &qoi,
                            const ROL::SharedPointer<const std::vector<Real> > &v,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                            const ROL::SharedPointer<const Tpetra::MultiVector<> > &z = ROL::nullPointer,
                            const ROL::SharedPointer<const std::vector<Real> > &z_param = ROL::nullPointer) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec33);
    #endif
    if ( z_param != ROL::nullPointer ) {
      const int size = z_param->size();
      if ( H33 == ROL::nullPointer ) {
        H33 = ROL::makeShared<std::vector<Real>>(size,0);
      }
      try {
        H33->assign(size,0);
        // Initialize local hessian times a vector storage
        std::vector<ROL::SharedPointer<Intrepid::FieldContainer<Real> > > locHess(size,ROL::nullPointer);
        // Get u_coeff from u and z_coeff from z
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPointer;
        getCoeffFromStateVector(u_coeff,u);
        ROL::SharedPointer<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPointer;
        getCoeffFromControlVector(z_coeff,z);
        // Compute local hessian times a vector
        (*H33) = qoi->HessVec_33(locHess,v,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        for (int i = 0; i < size; ++i) {
          dofMgr_->transformToFieldPattern(locHess[i]);
        }
        // Assembly
        std::vector<Real> myHess(size,0), globHess(size,0);
        for (int j = 0; j < size; ++j) {
          if ( locHess[j] != ROL::nullPointer ) {
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
  void assemblePDERieszMap1(ROL::SharedPointer<Tpetra::CrsMatrix<> > &R1,
                            const ROL::SharedPointer<PDE<Real> > &pde) {
    try {
      // Compute local state Riesz matrix
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > riesz;
      pde->RieszMap_1(riesz); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(riesz);
      // Initialize Riesz matrix if not done so already
      if ( R1 == ROL::nullPointer ) {
      R1 = ROL::makeShared<Tpetra::CrsMatrix<>>(matR1Graph_);
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
  void assemblePDERieszMap2(ROL::SharedPointer<Tpetra::CrsMatrix<> > &R2,
                            const ROL::SharedPointer<PDE<Real> > &pde) {
    try {
      // Compute local control Riesz matrix
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > riesz;
      pde->RieszMap_2(riesz); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(riesz);
      // Initialize Riesz matrix if not done so already
      if ( R2 == ROL::nullPointer ) {
        R2 = ROL::makeShared<Tpetra::CrsMatrix<>>(matR2Graph_);
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
  /* Compute error routines.                                                 */
  /***************************************************************************/
  Real computeStateError(const ROL::SharedPointer<const Tpetra::MultiVector<> > &soln,
                         const ROL::SharedPointer<Solution<Real> > &trueSoln,
                         const int cubDeg = 6,
                         const ROL::SharedPointer<FieldHelper<Real> > &fieldHelper = ROL::nullPointer) const {
    Real totalError(0);
    // populate inCoeffs
    ROL::SharedPointer<Intrepid::FieldContainer<Real> > inCoeffs0;
    getCoeffFromStateVector(inCoeffs0, soln);
    // split fields
    int numFields = 1;
    std::vector<ROL::SharedPointer<Intrepid::FieldContainer<Real> > > inCoeffs;
    if (fieldHelper != ROL::nullPointer) {
      numFields = fieldHelper->numFields();
      fieldHelper->splitFieldCoeff(inCoeffs,inCoeffs0);
    }
    else {
      inCoeffs.push_back(inCoeffs0);
    }
    // compute error
    for (int fn = 0; fn < numFields; ++fn) {
      // create fe object for error computation
      Intrepid::DefaultCubatureFactory<Real> cubFactory;
      shards::CellTopology cellType = basisPtrs_[fn]->getBaseCellTopology();
      ROL::SharedPointer<Intrepid::Cubature<Real> > cellCub = cubFactory.create(cellType, cubDeg);
      ROL::SharedPointer<FE<Real> > fe
        = ROL::makeShared<FE<Real>>(volCellNodes_,basisPtrs_[fn],cellCub);

      // get dimensions
      int c = fe->gradN()->dimension(0);
      int p = fe->gradN()->dimension(2);
      int d = fe->gradN()->dimension(3);

      // evaluate input coefficients on fe basis
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > funcVals
        = ROL::makeShared<Intrepid::FieldContainer<Real>>(c, p);
      fe->evaluateValue(funcVals, inCoeffs[fn]);

      // subtract off true solution
      std::vector<Real> x(d);
      for (int i=0; i<c; ++i) {
        for (int j=0; j<p; ++j) {
          for (int k=0; k<d; ++k) {
            x[k] = (*fe->cubPts())(i, j, k);
          }
          (*funcVals)(i, j) -= trueSoln->evaluate(x,fn);
        }
      }

      // compute norm squared of local error
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > normSquaredError
        = ROL::makeShared<Intrepid::FieldContainer<Real>>(c);
      fe->computeIntegral(normSquaredError, funcVals, funcVals, false);

      Real localErrorSum(0);
      Real globalErrorSum(0);
      for (int i=0; i<numCells_; ++i) {
        localErrorSum += (*normSquaredError)(i);
      }
      Teuchos::reduceAll<int, Real>(*comm_, Teuchos::REDUCE_SUM, 1, &localErrorSum, &globalErrorSum);
      totalError += globalErrorSum;
    }

    return std::sqrt(totalError);
  }

  Real computeControlError(const ROL::SharedPointer<const Tpetra::MultiVector<> > &soln,
                           const ROL::SharedPointer<Solution<Real> > &trueSoln,
                           const int cubDeg = 6,
                           const ROL::SharedPointer<FieldHelper<Real> > &fieldHelper = ROL::nullPointer) const {
    Real totalError(0);
    // populate inCoeffs
    ROL::SharedPointer<Intrepid::FieldContainer<Real> > inCoeffs0;
    getCoeffFromControlVector(inCoeffs0, soln);
    // split fields
    int numFields = 1;
    std::vector<ROL::SharedPointer<Intrepid::FieldContainer<Real> > > inCoeffs;
    if (fieldHelper != ROL::nullPointer) {
      numFields = fieldHelper->numFields();
      fieldHelper->splitFieldCoeff(inCoeffs,inCoeffs0);
    }
    else {
      inCoeffs.push_back(inCoeffs0);
    }
    // compute error
    for (int fn = 0; fn < numFields; ++fn) {
      // create fe object for error computation
      Intrepid::DefaultCubatureFactory<Real> cubFactory;
      shards::CellTopology cellType = basisPtrs_[fn]->getBaseCellTopology();
      ROL::SharedPointer<Intrepid::Cubature<Real> > cellCub = cubFactory.create(cellType, cubDeg);
      ROL::SharedPointer<FE<Real> > fe
        = ROL::makeShared<FE<Real>>(volCellNodes_,basisPtrs_[fn],cellCub);

      // get dimensions
      int c = fe->gradN()->dimension(0);
      int p = fe->gradN()->dimension(2);
      int d = fe->gradN()->dimension(3);

      // evaluate input coefficients on fe basis
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > funcVals
        = ROL::makeShared<Intrepid::FieldContainer<Real>>(c, p);
      fe->evaluateValue(funcVals, inCoeffs[fn]);

      // subtract off true solution
      std::vector<Real> x(d);
      for (int i=0; i<c; ++i) {
        for (int j=0; j<p; ++j) {
          for (int k=0; k<d; ++k) {
            x[k] = (*fe->cubPts())(i, j, k);
          }
          (*funcVals)(i, j) -= trueSoln->evaluate(x,fn);
        }
      }

      // compute norm squared of local error
      ROL::SharedPointer<Intrepid::FieldContainer<Real> > normSquaredError
        = ROL::makeShared<Intrepid::FieldContainer<Real>>(c);
      fe->computeIntegral(normSquaredError, funcVals, funcVals, false);

      Real localErrorSum(0);
      Real globalErrorSum(0);
      for (int i=0; i<numCells_; ++i) {
        localErrorSum += (*normSquaredError)(i);
      }
      Teuchos::reduceAll<int, Real>(*comm_, Teuchos::REDUCE_SUM, 1, &localErrorSum, &globalErrorSum);
      totalError += globalErrorSum;
    }

    return std::sqrt(totalError);
  }
  /***************************************************************************/
  /* End of compute solution routines.                                       */
  /***************************************************************************/

  /***************************************************************************/
  /* Output routines.                                                        */
  /***************************************************************************/
  void printMeshData(std::ostream &outStream) const {
    ROL::SharedPointer<Intrepid::FieldContainer<Real> > nodesPtr = meshMgr_->getNodes();
    ROL::SharedPointer<Intrepid::FieldContainer<int> >  cellToNodeMapPtr = meshMgr_->getCellToNodeMap();
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

  void outputTpetraVector(const ROL::SharedPointer<const Tpetra::MultiVector<> > &vec,
                          const std::string &filename) const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> > vecWriter;
    vecWriter.writeDenseFile(filename, vec);
  }

  void serialPrintStateEdgeField(const ROL::SharedPointer<const Tpetra::MultiVector<> > &u,
                                 const ROL::SharedPointer<FieldHelper<Real> > &fieldHelper,
                                 const std::string &filename,
                                 const ROL::SharedPointer<FE_CURL<Real> > &fe) const {
    const int c = fe->curlN()->dimension(0);
    const int f = fe->curlN()->dimension(1);
    const int p = 1, d = 3;

    ROL::SharedPointer<Intrepid::FieldContainer<Real> > u_coeff;
    getCoeffFromStateVector(u_coeff,u);

    std::vector<ROL::SharedPointer<Intrepid::FieldContainer<Real> > > U;
    fieldHelper->splitFieldCoeff(U, u_coeff);
    int numFields = U.size();

    // Transform cell center to physical
    ROL::SharedPointer<Intrepid::FieldContainer<Real> > rx
      = ROL::makeShared<Intrepid::FieldContainer<Real>>(p,d);
    ROL::SharedPointer<Intrepid::FieldContainer<Real> > px
      = ROL::makeShared<Intrepid::FieldContainer<Real>>(c,p,d);
    fe->mapRefPointsToPhysical(px,rx);
    // Transform reference values into physical space.
    ROL::SharedPointer<Intrepid::FieldContainer<Real> > cellJac
      = ROL::makeShared<Intrepid::FieldContainer<Real>>(c,p,d,d);
    ROL::SharedPointer<Intrepid::FieldContainer<Real> > cellJacInv
      = ROL::makeShared<Intrepid::FieldContainer<Real>>(c,p,d,d);
    ROL::SharedPointer<shards::CellTopology> cellTopo
      = ROL::makeShared<shards::CellTopology>(basisPtrs_[0]->getBaseCellTopology());
    ROL::SharedPointer<Intrepid::FieldContainer<Real> > valReference
      = ROL::makeShared<Intrepid::FieldContainer<Real>>(f,p,d);
    basisPtrs_[0]->getValues(*valReference,*rx,Intrepid::OPERATOR_VALUE);
    ROL::SharedPointer<Intrepid::FieldContainer<Real> > valPhysical
      = ROL::makeShared<Intrepid::FieldContainer<Real>>(c,f,p,d);
    Intrepid::CellTools<Real>::setJacobian(*cellJac,*rx,*volCellNodes_,*cellTopo);
    Intrepid::CellTools<Real>::setJacobianInv(*cellJacInv, *cellJac);
    Intrepid::FunctionSpaceTools::HCURLtransformVALUE<Real>(*valPhysical,
                                                            *cellJacInv,
                                                            *valReference);

    std::vector<ROL::SharedPointer<Intrepid::FieldContainer<Real> > > uval(numFields);
    for (int k = 0; k < numFields; ++k) {
      uval[k] = ROL::makeShared<Intrepid::FieldContainer<Real>>(c,p,d);
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
  const ROL::SharedPointer<const Tpetra::Map<> > getStateMap(void) const {
    return myUniqueStateMap_;
  }

  const ROL::SharedPointer<const Tpetra::Map<> > getControlMap(void) const {
    return myUniqueControlMap_;
  }

  const ROL::SharedPointer<const Tpetra::Map<> > getResidualMap(void) const {
    return myUniqueResidualMap_;
  }
 
  ROL::SharedPointer<Tpetra::MultiVector<> > createStateVector(void) const {
    return ROL::makeShared<Tpetra::MultiVector<>>(myUniqueStateMap_, 1, true);
  }
 
  ROL::SharedPointer<Tpetra::MultiVector<> > createControlVector(void) const {
    return ROL::makeShared<Tpetra::MultiVector<>>(myUniqueControlMap_, 1, true);
  }
 
  ROL::SharedPointer<Tpetra::MultiVector<> > createResidualVector(void) const {
    return ROL::makeShared<Tpetra::MultiVector<>>(myUniqueResidualMap_, 1, true);
  }
  /***************************************************************************/
  /* End of vector generation routines.                                      */
  /***************************************************************************/

  /***************************************************************************/
  /* Accessor routines.                                                      */
  /***************************************************************************/
  const ROL::SharedPointer<DofManager<Real> > getDofManager(void) const {
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
