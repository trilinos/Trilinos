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
    ROL::Ptr<Teuchos::Time> AssemblePDEResidual       = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Residual");
    ROL::Ptr<Teuchos::Time> AssemblePDEJacobian1      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Jacobian1");
    ROL::Ptr<Teuchos::Time> AssemblePDEJacobian2      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Jacobian2");
    ROL::Ptr<Teuchos::Time> AssemblePDEJacobian3      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Jacobian3");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian11      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian11");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian12      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian12");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian13      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian13");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian21      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian21");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian22      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian22");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian23      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian23");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian31      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian31");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian32      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian32");
    ROL::Ptr<Teuchos::Time> AssemblePDEHessian33      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble PDE Hessian33");
    ROL::Ptr<Teuchos::Time> AssembleQOIValue          = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Value");
    ROL::Ptr<Teuchos::Time> AssembleQOIGradient1      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Gradient1");
    ROL::Ptr<Teuchos::Time> AssembleQOIGradient2      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Gradient2");
    ROL::Ptr<Teuchos::Time> AssembleQOIGradient3      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI Gradient3");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec11      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec11");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec12      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec12");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec13      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec13");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec21      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec21");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec22      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec22");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec23      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec23");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec31      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec31");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec32      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec32");
    ROL::Ptr<Teuchos::Time> AssembleQOIHessVec33      = Teuchos::TimeMonitor::getNewCounter("ROL::PDEOPT: Assemble QOI HessVec33");
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
//  ROL::Ptr<Teuchos::Time::Time> timerSolverFactorization_;
//  ROL::Ptr<Teuchos::Time::Time> timerSolverSubstitution_;
//  ROL::Ptr<Teuchos::Time::Time> timerAssemblyNonlinear_;
//  ROL::Ptr<Teuchos::Time::Time> timerSolverUpdate_;

  // Set in Constructor.
  bool verbose_;
  bool isJ1Transposed_, isJ2Transposed_;

  // Set in SetCommunicator.
  ROL::Ptr<const Teuchos::Comm<int> > comm_;
  int myRank_, numProcs_;

  // Set in SetBasis.
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;

  // Set in SetDiscretization.
  ROL::Ptr<MeshManager<Real> > meshMgr_;
  ROL::Ptr<DofManager<Real> >  dofMgr_;

  // Set in SetParallelStructure.
  int numCells_;
  Teuchos::Array<int> myCellIds_;
  Teuchos::Array<int> cellOffsets_;
  ROL::Ptr<const Tpetra::Map<> > myOverlapStateMap_;
  ROL::Ptr<const Tpetra::Map<> > myUniqueStateMap_;
  ROL::Ptr<const Tpetra::Map<> > myOverlapControlMap_;
  ROL::Ptr<const Tpetra::Map<> > myUniqueControlMap_;
  ROL::Ptr<const Tpetra::Map<> > myOverlapResidualMap_;
  ROL::Ptr<const Tpetra::Map<> > myUniqueResidualMap_;
  ROL::Ptr<Tpetra::CrsGraph<> >  matJ1Graph_;
  ROL::Ptr<Tpetra::CrsGraph<> >  matJ2Graph_;
  ROL::Ptr<Tpetra::CrsGraph<> >  matR1Graph_;
  ROL::Ptr<Tpetra::CrsGraph<> >  matR2Graph_;
  ROL::Ptr<Tpetra::CrsGraph<> >  matH11Graph_;
  ROL::Ptr<Tpetra::CrsGraph<> >  matH12Graph_;
  ROL::Ptr<Tpetra::CrsGraph<> >  matH21Graph_;
  ROL::Ptr<Tpetra::CrsGraph<> >  matH22Graph_;

  // Set in SetCellNodes.
  ROL::Ptr<Intrepid::FieldContainer<Real> > volCellNodes_;
  ROL::Ptr<std::vector<std::vector<std::vector<int> > > >  bdryCellIds_;
  //std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<int> > > >  bdryCellLocIds_;
  std::vector<std::vector<std::vector<int> > >  bdryCellLocIds_;
  std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > bdryCellNodes_;

  // Finite element vectors and matrices for PDE.
  ROL::Ptr<Tpetra::MultiVector<> > pde_vecR_overlap_;
  ROL::Ptr<Tpetra::MultiVector<> > pde_vecJ3_overlap_;
  ROL::Ptr<Tpetra::MultiVector<> > pde_vecH13_overlap_;
  ROL::Ptr<Tpetra::MultiVector<> > pde_vecH23_overlap_;

  // Finite element vectors and matrices for QoI.
  ROL::Ptr<Tpetra::MultiVector<> > qoi_vecG1_overlap_;
  ROL::Ptr<Tpetra::MultiVector<> > qoi_vecG2_overlap_;
  ROL::Ptr<Tpetra::MultiVector<> > qoi_vecH11_overlap_;
  ROL::Ptr<Tpetra::MultiVector<> > qoi_vecH12_overlap_;
  ROL::Ptr<Tpetra::MultiVector<> > qoi_vecH13_overlap_;
  ROL::Ptr<Tpetra::MultiVector<> > qoi_vecH21_overlap_;
  ROL::Ptr<Tpetra::MultiVector<> > qoi_vecH22_overlap_;
  ROL::Ptr<Tpetra::MultiVector<> > qoi_vecH23_overlap_;

private:

  void setCommunicator(const ROL::Ptr<const Teuchos::Comm<int> > &comm,
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
         const std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > &basisPtrs,
         Teuchos::ParameterList &parlist,
         std::ostream &outStream = std::cout) {
    basisPtrs_ = basisPtrs;
    if (verbose_ && myRank_==0) {
      outStream << "Initialized PDE." << std::endl;
    }
  }

  void setDiscretization(Teuchos::ParameterList &parlist,
                         const ROL::Ptr<MeshManager<Real> > &meshMgr = ROL::nullPtr,
                         std::ostream &outStream = std::cout) {
    if ( meshMgr != ROL::nullPtr ) {
      // Use MeshManager object if supplied
      meshMgr_ = meshMgr;
    }
    else {
      // Otherwise construct MeshManager objective from parameter list
    }
    dofMgr_ = ROL::makePtr<DofManager<Real>>(meshMgr_,basisPtrs_);
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
    myOverlapStateMap_ = ROL::makePtr<Tpetra::Map<>>(
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
//    myCellMap_ = ROL::makePtr<Tpetra::Map<>>(
//                 Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
//                 myCellIds_, 0, comm_);

    /****************************************/
    /*** Assemble global graph structure. ***/
    /****************************************/
    matJ1Graph_ = ROL::makePtr<Tpetra::CrsGraph<>>(myUniqueStateMap_, 0);
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
    volCellNodes_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numNodesPerCell, spaceDim);
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
              bdryCellNodes_[i][j] = ROL::makePtr<Intrepid::FieldContainer<Real>>(myNumCellsSide, numNodesPerCell, spaceDim);
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

  void getCoeffFromStateVector(ROL::Ptr<Intrepid::FieldContainer<Real> > &xcoeff,
                               const ROL::Ptr<const Tpetra::MultiVector<> > &x) const {
    if ( x != ROL::nullPtr ) {
      // Perform import onto myOverlapMap
      ROL::Ptr<Tpetra::MultiVector<> > xshared =
        ROL::makePtr<Tpetra::MultiVector<>>(myOverlapStateMap_, 1, true);
      Tpetra::Import<> importer(myUniqueStateMap_, myOverlapStateMap_);
      xshared->doImport(*x,importer,Tpetra::REPLACE);
      // Populate xcoeff
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int lfs = dofMgr_->getLocalFieldSize();
      xcoeff = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs);
      Teuchos::ArrayRCP<const Real> xdata = xshared->get1dView();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<lfs; ++j) {
          (*xcoeff)(i,j) = xdata[xshared->getMap()->getLocalElement(cellDofs(myCellIds_[i],j))];
        }
      }
      dofMgr_->transformToIntrepidPattern(xcoeff);
    }
    else {
      xcoeff = ROL::nullPtr;
    }
  }

  void getCoeffFromControlVector(ROL::Ptr<Intrepid::FieldContainer<Real> > &xcoeff,
                                 const ROL::Ptr<const Tpetra::MultiVector<> > &x) const {
    if ( x != ROL::nullPtr ) {
      // Perform import onto myOverlapMap
      ROL::Ptr<Tpetra::MultiVector<> > xshared =
        ROL::makePtr<Tpetra::MultiVector<>>(myOverlapControlMap_, 1, true);
      Tpetra::Import<> importer(myUniqueControlMap_, myOverlapControlMap_);
      xshared->doImport(*x,importer,Tpetra::REPLACE);
      // Populate xcoeff
      Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
      int lfs = dofMgr_->getLocalFieldSize();
      xcoeff = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs);
      Teuchos::ArrayRCP<const Real> xdata = xshared->get1dView();
      for (int i=0; i<numCells_; ++i) {
        for (int j=0; j<lfs; ++j) {
          (*xcoeff)(i,j) = xdata[xshared->getMap()->getLocalElement(cellDofs(myCellIds_[i],j))];
        }
      }
      dofMgr_->transformToIntrepidPattern(xcoeff);
    }
    else {
      xcoeff = ROL::nullPtr;
    }
  }

public:
  // destructor
  virtual ~Assembler() {}

  // Constuctor: Discretization set from ParameterList
  Assembler(const std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > &basisPtrs,
          const ROL::Ptr<const Teuchos::Comm<int> > &comm,
          Teuchos::ParameterList &parlist,
          std::ostream &outStream = std::cout)
    : isJ1Transposed_(false), isJ2Transposed_(false) {
    setCommunicator(comm,parlist,outStream);
    setBasis(basisPtrs,parlist,outStream);
    setDiscretization(parlist,ROL::nullPtr,outStream);
    setParallelStructure(parlist,outStream);
    setCellNodes(outStream);
  }

  // Constructor: Discretization set from MeshManager input
  Assembler(const std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > &basisPtrs,
          const ROL::Ptr<MeshManager<Real> > &meshMgr,
          const ROL::Ptr<const Teuchos::Comm<int> > &comm,
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
  void assemblePDEResidual(ROL::Ptr<Tpetra::MultiVector<> > &r,
                           const ROL::Ptr<PDE<Real> > &pde,
                           const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                           const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                           const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEResidual);
    #endif
    // Initialize residual vectors if not done so
    if ( r == ROL::nullPtr ) { // Unique components of residual vector
      r = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueResidualMap_, 1, true);
    }
    if ( pde_vecR_overlap_ == ROL::nullPtr ) { // Overlapping components of residual vector
      pde_vecR_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapResidualMap_, 1, true);
    }
    // Set residual vectors to zero
    r->scale(static_cast<Real>(0));
    pde_vecR_overlap_->scale(static_cast<Real>(0));
    // Get degrees of freedom
    Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
    int numLocalDofs = cellDofs.dimension(1);
    // Initialize res container
    ROL::Ptr<Intrepid::FieldContainer<Real> > res;
    // Get u_coeff from u and z_coeff from z
    ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
    getCoeffFromStateVector(u_coeff,u);
    ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
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

  void assemblePDEJacobian1(ROL::Ptr<Tpetra::CrsMatrix<> > &J1,
                            const ROL::Ptr<PDE<Real> > &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEJacobian1);
    #endif
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
      getCoeffFromStateVector(u_coeff,u);
      ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
      getCoeffFromControlVector(z_coeff,z);
      // Compute PDE Jacobian
      ROL::Ptr<Intrepid::FieldContainer<Real> > jac;
      pde->Jacobian_1(jac,u_coeff,z_coeff,z_param); // Throw if not implemented or is zero
      dofMgr_->transformToFieldPattern(jac);
      // Initialize Jacobian matrices
      if ( J1 == ROL::nullPtr ) {
        J1 = ROL::makePtr<Tpetra::CrsMatrix<>>(matJ1Graph_);
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

  void assemblePDEJacobian2(ROL::Ptr<Tpetra::CrsMatrix<> > &J2,
                            const ROL::Ptr<PDE<Real> > &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEJacobian2);
    #endif
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
      getCoeffFromStateVector(u_coeff,u);
      ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
      getCoeffFromControlVector(z_coeff,z);
      // Compute PDE Jacobian
      ROL::Ptr<Intrepid::FieldContainer<Real> > jac;
      pde->Jacobian_2(jac,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(jac);
      // Initialize Jacobian matrices
      if ( J2 == ROL::nullPtr ) {
        J2 = ROL::makePtr<Tpetra::CrsMatrix<>>(matJ2Graph_);
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

  void assemblePDEJacobian3(ROL::Ptr<Tpetra::MultiVector<> > &J3,
                            const ROL::Ptr<PDE<Real> > &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEJacobian3);
    #endif
    if ( z_param != ROL::nullPtr ) {
      try {
        int size = static_cast<int>(z_param->size());
        // Initialize res
        std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > jac(size);
        // Get u_coeff from u and z_coeff from z
        ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
        getCoeffFromStateVector(u_coeff,u);
        ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
        getCoeffFromControlVector(z_coeff,z);
        // Compute PDE local Jacobian wrt parametric controls
        pde->Jacobian_3(jac,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        for (int i = 0; i < size; ++i) {
          dofMgr_->transformToFieldPattern(jac[i]);
        }
        // Initialize Jacobian storage if not done so already
        if (J3 == ROL::nullPtr) {
          J3 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueResidualMap_, size, true);
        }
        if ( pde_vecJ3_overlap_ == ROL::nullPtr) {
          pde_vecJ3_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapResidualMap_, size, true);
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

  void assemblePDEHessian11(ROL::Ptr<Tpetra::CrsMatrix<> > &H11,
                            const ROL::Ptr<PDE<Real> > &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &l,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian11);
    #endif
    try {
      ROL::Ptr<Intrepid::FieldContainer<Real> > hess; 
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
      getCoeffFromStateVector(u_coeff,u);
      ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
      getCoeffFromControlVector(z_coeff,z);
      ROL::Ptr<Intrepid::FieldContainer<Real> > l_coeff = ROL::nullPtr;
      getCoeffFromStateVector(l_coeff,l);
      // Compute PDE Hessian
      pde->Hessian_11(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(hess);
      // Initialize Hessian storage if not done so already
      if ( H11 == ROL::nullPtr ) {
        H11 = ROL::makePtr<Tpetra::CrsMatrix<>>(matH11Graph_);
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

  void assemblePDEHessian12(ROL::Ptr<Tpetra::CrsMatrix<> > &H12,
                            const ROL::Ptr<PDE<Real> > &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &l,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian12);
    #endif
    try {
      ROL::Ptr<Intrepid::FieldContainer<Real> > hess;
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
      getCoeffFromStateVector(u_coeff,u);
      ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
      getCoeffFromControlVector(z_coeff,z);
      ROL::Ptr<Intrepid::FieldContainer<Real> > l_coeff = ROL::nullPtr;
      getCoeffFromStateVector(l_coeff,l);
      // Compute PDE Hessian
      pde->Hessian_12(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(hess);
      // Initialize Hessian storage if not done so already
      if ( H12 == ROL::nullPtr ) {
        H12 = ROL::makePtr<Tpetra::CrsMatrix<>>(matH12Graph_);
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

  void assemblePDEHessian13(ROL::Ptr<Tpetra::MultiVector<> > &H13,
                            const ROL::Ptr<PDE<Real> > &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &l,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian13);
    #endif
    if ( z_param != ROL::nullPtr ) {
      try {
        int size = static_cast<int>(z_param->size());
        // Initialize local hessian
        std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > hess(size);
        // Get u_coeff from u, z_coeff from z and l_coeff from l
        ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
        getCoeffFromStateVector(u_coeff,u);
        ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
        getCoeffFromControlVector(z_coeff,z);
        ROL::Ptr<Intrepid::FieldContainer<Real> > l_coeff = ROL::nullPtr;
        getCoeffFromStateVector(l_coeff,l);
        // Compute PDE local Jacobian wrt parametric controls
        pde->Hessian_13(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        for (int i = 0; i < size; ++i) {
          dofMgr_->transformToFieldPattern(hess[i]);
        }
        // Initialize Jacobian storage if not done so already
        if (H13 == ROL::nullPtr) {
          H13 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueStateMap_, size, true);
        }
        if ( pde_vecH13_overlap_ == ROL::nullPtr) {
          pde_vecH13_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapStateMap_, size, true);
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

  void assemblePDEHessian21(ROL::Ptr<Tpetra::CrsMatrix<> > &H21,
                            const ROL::Ptr<PDE<Real> > &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &l,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian21);
    #endif
    try {
      ROL::Ptr<Intrepid::FieldContainer<Real> > hess;
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
      getCoeffFromStateVector(u_coeff,u);
      ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
      getCoeffFromControlVector(z_coeff,z);
      ROL::Ptr<Intrepid::FieldContainer<Real> > l_coeff = ROL::nullPtr;
      getCoeffFromStateVector(l_coeff,l);
      // Compute PDE Hessian
      pde->Hessian_21(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(hess);
      // Initialize Hessian storage if not done so already
      if ( H21 == ROL::nullPtr ) {
        H21 = ROL::makePtr<Tpetra::CrsMatrix<>>(matH21Graph_);
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

  void assemblePDEHessian22(ROL::Ptr<Tpetra::CrsMatrix<> > &H22,
                            const ROL::Ptr<PDE<Real> > &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &l,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian22);
    #endif
    try {
      ROL::Ptr<Intrepid::FieldContainer<Real> > hess;
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
      getCoeffFromStateVector(u_coeff,u);
      ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
      getCoeffFromControlVector(z_coeff,z);
      ROL::Ptr<Intrepid::FieldContainer<Real> > l_coeff = ROL::nullPtr;
      getCoeffFromStateVector(l_coeff,l);
      // Compute PDE Hessian
      pde->Hessian_22(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(hess);
      // Initialize Hessian storage if not done so already
      if ( H22 == ROL::nullPtr ) {
        H22 = ROL::makePtr<Tpetra::CrsMatrix<>>(matH22Graph_);
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

  void assemblePDEHessian23(ROL::Ptr<Tpetra::MultiVector<> > &H23,
                            const ROL::Ptr<PDE<Real> > &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &l,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian23);
    #endif
    if ( z_param != ROL::nullPtr ) {
      try {
        int size = static_cast<int>(z_param->size());
        // Initialize local hessian
        std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > hess(size);
        // Get u_coeff from u, z_coeff from z and l_coeff from l
        ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
        getCoeffFromStateVector(u_coeff,u);
        ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
        getCoeffFromControlVector(z_coeff,z);
        ROL::Ptr<Intrepid::FieldContainer<Real> > l_coeff = ROL::nullPtr;
        getCoeffFromStateVector(l_coeff,l);
        // Compute PDE local Jacobian wrt parametric controls
        pde->Hessian_23(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        for (int i = 0; i < size; ++i) {
          dofMgr_->transformToFieldPattern(hess[i]);
        }
        // Initialize Jacobian storage if not done so already
        if (H23 == ROL::nullPtr) {
          H23 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueControlMap_, size, true);
        }
        if ( pde_vecH23_overlap_ == ROL::nullPtr) {
          pde_vecH23_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapControlMap_, size, true);
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

  void assemblePDEHessian31(ROL::Ptr<Tpetra::MultiVector<> > &H31,
                            const ROL::Ptr<PDE<Real> > &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &l,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian31);
    #endif
    assemblePDEHessian13(H31,pde,l,u,z,z_param);
  }

  void assemblePDEHessian32(ROL::Ptr<Tpetra::MultiVector<> > &H32,
                            const ROL::Ptr<PDE<Real> > &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &l,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian32);
    #endif
    assemblePDEHessian23(H32,pde,l,u,z,z_param);
  }

  void assemblePDEHessian33(ROL::Ptr<std::vector<std::vector<Real> > > &H33,
                            const ROL::Ptr<PDE<Real> > &pde,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &l,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssemblePDEHessian33);
    #endif
    if ( z_param != ROL::nullPtr ) {
      try {
        int size = static_cast<int>(z_param->size());
        // Initialize local hessian
        std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > tmp(size,ROL::nullPtr);
        std::vector<std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > > hess(size,tmp);
        // Get u_coeff from u, z_coeff from z and l_coeff from l
        ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
        getCoeffFromStateVector(u_coeff,u);
        ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
        getCoeffFromControlVector(z_coeff,z);
        ROL::Ptr<Intrepid::FieldContainer<Real> > l_coeff = ROL::nullPtr;
        getCoeffFromStateVector(l_coeff,l);
        // Compute PDE local Jacobian wrt parametric controls
        pde->Hessian_33(hess,l_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        for (int i = 0; i < size; ++i) {
          for (int j = 0; j < size; ++j) {
            dofMgr_->transformToFieldPattern(hess[i][j]);
          }
        }
        // Initialize Jacobian storage if not done so already
        if (H33 == ROL::nullPtr) {
          std::vector<Real> col(size,static_cast<Real>(0));
          H33 = ROL::makePtr<std::vector<std::vector<Real> >>(size,col);
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
  Real assembleQoIValue(const ROL::Ptr<QoI<Real> > &qoi,
                        const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                        const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                        const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIValue);
    #endif
    Real val(0);
    try {
      // Integrate obj object
      ROL::Ptr<Intrepid::FieldContainer<Real> > locVal;
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
      getCoeffFromStateVector(u_coeff,u);
      ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
      getCoeffFromControlVector(z_coeff,z);
      // Get OBJ_CELL value
      val = qoi->value(locVal,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      // Assembly
      if ( locVal != ROL::nullPtr ) {
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

  void assembleQoIGradient1(ROL::Ptr<Tpetra::MultiVector<> > &g1,
                            const ROL::Ptr<QoI<Real> > &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIGradient1);
    #endif
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
      getCoeffFromStateVector(u_coeff,u);
      ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      ROL::Ptr<Intrepid::FieldContainer<Real> > locGrad;
      qoi->gradient_1(locGrad,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(locGrad);
      // Initialize state QoI gradient vectors
      if ( g1 == ROL::nullPtr ) {
        g1 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueStateMap_, 1, true);
      }
      if ( qoi_vecG1_overlap_ == ROL::nullPtr ) {
        qoi_vecG1_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapStateMap_, 1, true);
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

  void assembleQoIGradient2(ROL::Ptr<Tpetra::MultiVector<> > &g2,
                            const ROL::Ptr<QoI<Real> > &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIGradient2);
    #endif
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
      getCoeffFromStateVector(u_coeff,u);
      ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      ROL::Ptr<Intrepid::FieldContainer<Real> > locGrad;
      qoi->gradient_2(locGrad,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(locGrad);
      // Initialize control gradient vectors
      if ( g2 == ROL::nullPtr ) {
        g2 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueControlMap_, 1, true);
      }
      if ( qoi_vecG2_overlap_ == ROL::nullPtr ) {
        qoi_vecG2_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapControlMap_, 1, true);
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

  void assembleQoIGradient3(ROL::Ptr<std::vector<Real> > &g3,
                            const ROL::Ptr<QoI<Real> > &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIGradient3);
    #endif
    if ( z_param != ROL::nullPtr ) {
      const int size = z_param->size();
      if ( g3 == ROL::nullPtr ) {
        g3 = ROL::makePtr<std::vector<Real>>(size,0);
      }
      try {
        g3->assign(size,0);
        // Initialize local gradient storage
        std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > locGrad(size,ROL::nullPtr);
        // Get u_coeff from u and z_coeff from z
        ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
        getCoeffFromStateVector(u_coeff,u);
        ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
        getCoeffFromControlVector(z_coeff,z);
        // Compute gradient
        (*g3) = qoi->gradient_3(locGrad,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        for (int i = 0; i < size; ++i) {
          dofMgr_->transformToFieldPattern(locGrad[i]);
        }
        // Assembly
        std::vector<Real> myGrad(size,0), globGrad(size,0);
        for (int j = 0; j < size; ++j) {
          if ( locGrad[j] != ROL::nullPtr ) {
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

  void assembleQoIHessVec11(ROL::Ptr<Tpetra::MultiVector<> > &H11,
                            const ROL::Ptr<QoI<Real> > &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &v,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec11);
    #endif
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real> > locHess;
      ROL::Ptr<Intrepid::FieldContainer<Real> > v_coeff = ROL::nullPtr;
      getCoeffFromStateVector(v_coeff,v);
      ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
      getCoeffFromStateVector(u_coeff,u);
      ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi->HessVec_11(locHess, v_coeff, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(locHess);
      // Initialize state-state HessVec vectors
      if ( H11 == ROL::nullPtr ) {
        H11 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueStateMap_, 1, true);
      }
      if ( qoi_vecH11_overlap_ == ROL::nullPtr ) {
        qoi_vecH11_overlap_  = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapStateMap_, 1, true);
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

  void assembleQoIHessVec12(ROL::Ptr<Tpetra::MultiVector<> > &H12,
                            const ROL::Ptr<QoI<Real> > &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &v,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec12);
    #endif
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real> > locHess;
      ROL::Ptr<Intrepid::FieldContainer<Real> > v_coeff = ROL::nullPtr;
      getCoeffFromControlVector(v_coeff,v);
      ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
      getCoeffFromStateVector(u_coeff,u);
      ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi->HessVec_12(locHess, v_coeff, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(locHess);
      // Initialize state-control HessVec vectors
      if ( H12 == ROL::nullPtr ) {
        H12 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueStateMap_, 1, true);
      }
      if ( qoi_vecH12_overlap_ == ROL::nullPtr ) {
        qoi_vecH12_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapStateMap_, 1, true);
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

  void assembleQoIHessVec13(ROL::Ptr<Tpetra::MultiVector<> > &H13,
                            const ROL::Ptr<QoI<Real> > &qoi,
                            const ROL::Ptr<const std::vector<Real> > &v,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > &z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec13);
    #endif
    if (z_param != ROL::nullPtr) {
      try {
        // Get u_coeff from u and z_coeff from z
        ROL::Ptr<Intrepid::FieldContainer<Real> > locHess;
        ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
        getCoeffFromStateVector(u_coeff,u);
        ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
        getCoeffFromControlVector(z_coeff,z);
        // Compute local gradient
        qoi->HessVec_13(locHess, v, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
        dofMgr_->transformToFieldPattern(locHess);
        // Initialize state-control HessVec vectors
        if ( H13 == ROL::nullPtr ) {
          H13 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueStateMap_, 1, true);
        }
        if ( qoi_vecH13_overlap_ == ROL::nullPtr ) {
          qoi_vecH13_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapStateMap_, 1, true);
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

  void assembleQoIHessVec21(ROL::Ptr<Tpetra::MultiVector<> > &H21,
                            const ROL::Ptr<QoI<Real> > &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &v,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec21);
    #endif
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real> > locHess;
      ROL::Ptr<Intrepid::FieldContainer<Real> > v_coeff = ROL::nullPtr;
      getCoeffFromStateVector(v_coeff,v);
      ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
      getCoeffFromStateVector(u_coeff,u);
      ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi->HessVec_21(locHess, v_coeff, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(locHess);
      // Initialize control-state HessVec vectors
      if ( H21 == ROL::nullPtr ) {
        H21 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueControlMap_, 1, true);
      }
      if ( qoi_vecH21_overlap_ == ROL::nullPtr ) {
        qoi_vecH21_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapControlMap_, 1, true);
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

  void assembleQoIHessVec22(ROL::Ptr<Tpetra::MultiVector<> > &H22,
                            const ROL::Ptr<QoI<Real> > &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &v,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec22);
    #endif
    try {
      // Get u_coeff from u and z_coeff from z
      ROL::Ptr<Intrepid::FieldContainer<Real> > locHess;
      ROL::Ptr<Intrepid::FieldContainer<Real> > v_coeff = ROL::nullPtr;
      getCoeffFromControlVector(v_coeff,v);
      ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
      getCoeffFromStateVector(u_coeff,u);
      ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
      getCoeffFromControlVector(z_coeff,z);
      // Compute local gradient
      qoi->HessVec_22(locHess, v_coeff, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(locHess);
      // Initialize control-control HessVec vectors
      if ( H22 == ROL::nullPtr ) {
        H22 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueControlMap_, 1, true);
      }
      if ( qoi_vecH22_overlap_ == ROL::nullPtr ) {
        qoi_vecH22_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapControlMap_, 1, true);
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

  void assembleQoIHessVec23(ROL::Ptr<Tpetra::MultiVector<> > &H23,
                            const ROL::Ptr<QoI<Real> > &qoi,
                            const ROL::Ptr<const std::vector<Real> > &v,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec23);
    #endif
    if (z_param != ROL::nullPtr) {
      try {
        // Get u_coeff from u and z_coeff from z
        ROL::Ptr<Intrepid::FieldContainer<Real> > locHess;
        ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
        getCoeffFromStateVector(u_coeff,u);
        ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
        getCoeffFromControlVector(z_coeff,z);
        // Compute local gradient
        qoi->HessVec_23(locHess, v, u_coeff, z_coeff, z_param); // Throw if not implemented or zero
        dofMgr_->transformToFieldPattern(locHess);
        // Initialize control-control HessVec vectors
        if ( H23 == ROL::nullPtr ) {
          H23 = ROL::makePtr<Tpetra::MultiVector<>>(myUniqueControlMap_, 1, true);
        }
        if ( qoi_vecH23_overlap_ == ROL::nullPtr ) {
          qoi_vecH23_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapControlMap_, 1, true);
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

  void assembleQoIHessVec31(ROL::Ptr<std::vector<Real> > &H31,
                            const ROL::Ptr<QoI<Real> > &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &v,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec31);
    #endif
    if ( z_param != ROL::nullPtr ) {
      const int size = z_param->size();
      if ( H31 == ROL::nullPtr ) {
        H31 = ROL::makePtr<std::vector<Real>>(size,0);
      }
      try {
        H31->assign(size,0);
        // Initialize local gradient storage
        std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > locHess(size,ROL::nullPtr);
        // Get u_coeff from u and z_coeff from z
        ROL::Ptr<Intrepid::FieldContainer<Real> > v_coeff = ROL::nullPtr;
        getCoeffFromStateVector(v_coeff,v);
        ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
        getCoeffFromStateVector(u_coeff,u);
        ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
        getCoeffFromControlVector(z_coeff,z);
        // Compute gradient
        (*H31) = qoi->HessVec_31(locHess,v_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        for (int i = 0; i < size; ++i) {
          dofMgr_->transformToFieldPattern(locHess[i]);
        }
        // Assembly
        std::vector<Real> myHess(size,0), globHess(size,0);
        for (int j = 0; j < size; ++j) {
          if ( locHess[j] != ROL::nullPtr ) {
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

  void assembleQoIHessVec32(ROL::Ptr<std::vector<Real> > &H32,
                            const ROL::Ptr<QoI<Real> > &qoi,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &v,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > & z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec32);
    #endif
    if ( z_param != ROL::nullPtr ) {
      const int size = z_param->size();
      if ( H32 == ROL::nullPtr ) {
        H32 = ROL::makePtr<std::vector<Real>>(size,0);
      }
      try {
        H32->assign(size,0);
        // Initialize local hessian times a vector storage
        std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > locHess(size,ROL::nullPtr);
        // Get u_coeff from u and z_coeff from z
        ROL::Ptr<Intrepid::FieldContainer<Real> > v_coeff = ROL::nullPtr;
        getCoeffFromControlVector(v_coeff,v);
        ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
        getCoeffFromStateVector(u_coeff,u);
        ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
        getCoeffFromControlVector(z_coeff,z);
        // Compute local hessian times a vector
        (*H32) = qoi->HessVec_32(locHess,v_coeff,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        for (int i = 0; i < size; ++i) {
          dofMgr_->transformToFieldPattern(locHess[i]);
        }
        // Assembly
        std::vector<Real> myHess(size,0), globHess(size,0);
        for (int j = 0; j < size; ++j) {
          if ( locHess[j] != ROL::nullPtr ) {
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

  void assembleQoIHessVec33(ROL::Ptr<std::vector<Real> > &H33,
                            const ROL::Ptr<QoI<Real> > &qoi,
                            const ROL::Ptr<const std::vector<Real> > &v,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                            const ROL::Ptr<const Tpetra::MultiVector<> > &z = ROL::nullPtr,
                            const ROL::Ptr<const std::vector<Real> > &z_param = ROL::nullPtr) {
    #ifdef ROL_TIMERS
      Teuchos::TimeMonitor LocalTimer(*ROL::PDEOPT::AssembleQOIHessVec33);
    #endif
    if ( z_param != ROL::nullPtr ) {
      const int size = z_param->size();
      if ( H33 == ROL::nullPtr ) {
        H33 = ROL::makePtr<std::vector<Real>>(size,0);
      }
      try {
        H33->assign(size,0);
        // Initialize local hessian times a vector storage
        std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > locHess(size,ROL::nullPtr);
        // Get u_coeff from u and z_coeff from z
        ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff = ROL::nullPtr;
        getCoeffFromStateVector(u_coeff,u);
        ROL::Ptr<Intrepid::FieldContainer<Real> > z_coeff = ROL::nullPtr;
        getCoeffFromControlVector(z_coeff,z);
        // Compute local hessian times a vector
        (*H33) = qoi->HessVec_33(locHess,v,u_coeff,z_coeff,z_param); // Throw if not implemented or zero
        for (int i = 0; i < size; ++i) {
          dofMgr_->transformToFieldPattern(locHess[i]);
        }
        // Assembly
        std::vector<Real> myHess(size,0), globHess(size,0);
        for (int j = 0; j < size; ++j) {
          if ( locHess[j] != ROL::nullPtr ) {
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
  void assemblePDERieszMap1(ROL::Ptr<Tpetra::CrsMatrix<> > &R1,
                            const ROL::Ptr<PDE<Real> > &pde) {
    try {
      // Compute local state Riesz matrix
      ROL::Ptr<Intrepid::FieldContainer<Real> > riesz;
      pde->RieszMap_1(riesz); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(riesz);
      // Initialize Riesz matrix if not done so already
      if ( R1 == ROL::nullPtr ) {
      R1 = ROL::makePtr<Tpetra::CrsMatrix<>>(matR1Graph_);
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
  void assemblePDERieszMap2(ROL::Ptr<Tpetra::CrsMatrix<> > &R2,
                            const ROL::Ptr<PDE<Real> > &pde) {
    try {
      // Compute local control Riesz matrix
      ROL::Ptr<Intrepid::FieldContainer<Real> > riesz;
      pde->RieszMap_2(riesz); // Throw if not implemented or zero
      dofMgr_->transformToFieldPattern(riesz);
      // Initialize Riesz matrix if not done so already
      if ( R2 == ROL::nullPtr ) {
        R2 = ROL::makePtr<Tpetra::CrsMatrix<>>(matR2Graph_);
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
  Real computeStateError(const ROL::Ptr<const Tpetra::MultiVector<> > &soln,
                         const ROL::Ptr<Solution<Real> > &trueSoln,
                         const int cubDeg = 6,
                         const ROL::Ptr<FieldHelper<Real> > &fieldHelper = ROL::nullPtr) const {
    Real totalError(0);
    // populate inCoeffs
    ROL::Ptr<Intrepid::FieldContainer<Real> > inCoeffs0;
    getCoeffFromStateVector(inCoeffs0, soln);
    // split fields
    int numFields = 1;
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > inCoeffs;
    if (fieldHelper != ROL::nullPtr) {
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
      ROL::Ptr<Intrepid::Cubature<Real> > cellCub = cubFactory.create(cellType, cubDeg);
      ROL::Ptr<FE<Real> > fe
        = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrs_[fn],cellCub);

      // get dimensions
      int c = fe->gradN()->dimension(0);
      int p = fe->gradN()->dimension(2);
      int d = fe->gradN()->dimension(3);

      // evaluate input coefficients on fe basis
      ROL::Ptr<Intrepid::FieldContainer<Real> > funcVals
        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
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
      ROL::Ptr<Intrepid::FieldContainer<Real> > normSquaredError
        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
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

  Real computeControlError(const ROL::Ptr<const Tpetra::MultiVector<> > &soln,
                           const ROL::Ptr<Solution<Real> > &trueSoln,
                           const int cubDeg = 6,
                           const ROL::Ptr<FieldHelper<Real> > &fieldHelper = ROL::nullPtr) const {
    Real totalError(0);
    // populate inCoeffs
    ROL::Ptr<Intrepid::FieldContainer<Real> > inCoeffs0;
    getCoeffFromControlVector(inCoeffs0, soln);
    // split fields
    int numFields = 1;
    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > inCoeffs;
    if (fieldHelper != ROL::nullPtr) {
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
      ROL::Ptr<Intrepid::Cubature<Real> > cellCub = cubFactory.create(cellType, cubDeg);
      ROL::Ptr<FE<Real> > fe
        = ROL::makePtr<FE<Real>>(volCellNodes_,basisPtrs_[fn],cellCub);

      // get dimensions
      int c = fe->gradN()->dimension(0);
      int p = fe->gradN()->dimension(2);
      int d = fe->gradN()->dimension(3);

      // evaluate input coefficients on fe basis
      ROL::Ptr<Intrepid::FieldContainer<Real> > funcVals
        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c, p);
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
      ROL::Ptr<Intrepid::FieldContainer<Real> > normSquaredError
        = ROL::makePtr<Intrepid::FieldContainer<Real>>(c);
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
    ROL::Ptr<Intrepid::FieldContainer<Real> > nodesPtr = meshMgr_->getNodes();
    ROL::Ptr<Intrepid::FieldContainer<int> >  cellToNodeMapPtr = meshMgr_->getCellToNodeMap();
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

  void outputTpetraVector(const ROL::Ptr<const Tpetra::MultiVector<> > &vec,
                          const std::string &filename) const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> > vecWriter;
    vecWriter.writeDenseFile(filename, vec);
  }

  void serialPrintStateEdgeField(const ROL::Ptr<const Tpetra::MultiVector<> > &u,
                                 const ROL::Ptr<FieldHelper<Real> > &fieldHelper,
                                 const std::string &filename,
                                 const ROL::Ptr<FE_CURL<Real> > &fe) const {
    const int c = fe->curlN()->dimension(0);
    const int f = fe->curlN()->dimension(1);
    const int p = 1, d = 3;

    ROL::Ptr<Intrepid::FieldContainer<Real> > u_coeff;
    getCoeffFromStateVector(u_coeff,u);

    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > U;
    fieldHelper->splitFieldCoeff(U, u_coeff);
    int numFields = U.size();

    // Transform cell center to physical
    ROL::Ptr<Intrepid::FieldContainer<Real> > rx
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(p,d);
    ROL::Ptr<Intrepid::FieldContainer<Real> > px
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
    fe->mapRefPointsToPhysical(px,rx);
    // Transform reference values into physical space.
    ROL::Ptr<Intrepid::FieldContainer<Real> > cellJac
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d,d);
    ROL::Ptr<Intrepid::FieldContainer<Real> > cellJacInv
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d,d);
    ROL::Ptr<shards::CellTopology> cellTopo
      = ROL::makePtr<shards::CellTopology>(basisPtrs_[0]->getBaseCellTopology());
    ROL::Ptr<Intrepid::FieldContainer<Real> > valReference
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(f,p,d);
    basisPtrs_[0]->getValues(*valReference,*rx,Intrepid::OPERATOR_VALUE);
    ROL::Ptr<Intrepid::FieldContainer<Real> > valPhysical
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,p,d);
    Intrepid::CellTools<Real>::setJacobian(*cellJac,*rx,*volCellNodes_,*cellTopo);
    Intrepid::CellTools<Real>::setJacobianInv(*cellJacInv, *cellJac);
    Intrepid::FunctionSpaceTools::HCURLtransformVALUE<Real>(*valPhysical,
                                                            *cellJacInv,
                                                            *valReference);

    std::vector<ROL::Ptr<Intrepid::FieldContainer<Real> > > uval(numFields);
    for (int k = 0; k < numFields; ++k) {
      uval[k] = ROL::makePtr<Intrepid::FieldContainer<Real>>(c,p,d);
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
  const ROL::Ptr<const Tpetra::Map<> > getStateMap(void) const {
    return myUniqueStateMap_;
  }

  const ROL::Ptr<const Tpetra::Map<> > getControlMap(void) const {
    return myUniqueControlMap_;
  }

  const ROL::Ptr<const Tpetra::Map<> > getResidualMap(void) const {
    return myUniqueResidualMap_;
  }
 
  ROL::Ptr<Tpetra::MultiVector<> > createStateVector(void) const {
    return ROL::makePtr<Tpetra::MultiVector<>>(myUniqueStateMap_, 1, true);
  }
 
  ROL::Ptr<Tpetra::MultiVector<> > createControlVector(void) const {
    return ROL::makePtr<Tpetra::MultiVector<>>(myUniqueControlMap_, 1, true);
  }
 
  ROL::Ptr<Tpetra::MultiVector<> > createResidualVector(void) const {
    return ROL::makePtr<Tpetra::MultiVector<>>(myUniqueResidualMap_, 1, true);
  }
  /***************************************************************************/
  /* End of vector generation routines.                                      */
  /***************************************************************************/

  /***************************************************************************/
  /* Accessor routines.                                                      */
  /***************************************************************************/
  const ROL::Ptr<DofManager<Real> > getDofManager(void) const {
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
