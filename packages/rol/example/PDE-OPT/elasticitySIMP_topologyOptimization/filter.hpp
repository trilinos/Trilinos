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

/*! \file  filter.hpp
    \brief Density filtering based on solving a PDE.
*/

#ifndef ROL_DENSITY_FILTER_H
#define ROL_DENSITY_FILTER_H

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

#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "Amesos2.hpp"

#include "../TOOLS/dofmanager.hpp"

template<class Real>
class DensityFilter {

private:
  Teuchos::RCP<MeshManager<Real> > meshMgr_;
  Teuchos::RCP<DofManager<Real> >  dofMgr_;
  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;

  Teuchos::RCP<const Teuchos::Comm<int> > commPtr_;
  int myRank_;
  int numProcs_;

  int  basisOrder_;

  Teuchos::RCP<const Tpetra::Map<> >    myOverlapMap_;
  Teuchos::RCP<const Tpetra::Map<> >    myUniqueMap_;
  Teuchos::RCP<const Tpetra::Map<> >    myBColumnMap_;
  Teuchos::RCP<Tpetra::CrsGraph<> >     matAGraph_;
  Teuchos::RCP<Tpetra::CrsGraph<> >     matBGraph_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >    matA_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >    matB_;
  Teuchos::RCP<Tpetra::CrsMatrix<> >    matB_trans_;
  Teuchos::RCP<Tpetra::MultiVector<> >  vecCellVolumes_;

  Teuchos::Array<int> myCellIds_;
  Teuchos::Array<Real> myCellVolumes_;

  Teuchos::RCP<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > solverA_;

  shards::CellTopology cellType_;
  int spaceDim_;
  int numNodesPerCell_;
  int numCubPoints_;

  int totalNumCells_;
  int totalNumDofs_;
  int numCells_;

  Teuchos::RCP<Intrepid::FieldContainer<Real> > cubPoints_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cubWeights_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellNodes_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellJac_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellJacInv_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellJacDet_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellWeightedMeasure_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > valReference_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradReference_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > valPhysical_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradPhysical_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > kappaGradPhysical_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > valPhysicalWeighted_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradPhysicalWeighted_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradgradMats_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > valvalMats_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > valMats_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > onesVec_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cubPointsPhysical_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > kappa_;

  Real lengthScale_;
  bool enableFilter_;

public:

  DensityFilter(const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                const Teuchos::RCP<Teuchos::ParameterList> &parlist,
                const Teuchos::RCP<std::ostream> &outStream) {

    lengthScale_  = parlist->sublist("Density Filter").get<Real>("Length Scale");
    lengthScale_  = std::pow(lengthScale_/static_cast<Real>(2*std::sqrt(3)), 2);
    enableFilter_ = parlist->sublist("Density Filter").get<bool>("Enable");

    /************************************/
    /*** Retrieve communication data. ***/
    /************************************/
    commPtr_  = comm;
    myRank_   = commPtr_->getRank();
    numProcs_ = commPtr_->getSize();
    *outStream << "Total number of processors: " << numProcs_ << std::endl;
    /************************************/
    /************************************/

    /*************************************/
    /*** Retrieve parameter list data. ***/
    /*************************************/
    basisOrder_ = parlist->sublist("PDE FEM").get("Order of FE Discretization", 1);
    int cellSplit = parlist->sublist("Geometry").get("Partition type", 1);
    /*************************************/
    /*************************************/

    /****************************************************************************/
    /*** Initialize mesh / finite element fields / degree-of-freedom manager. ***/
    /****************************************************************************/

    // Mesh manager.
    meshMgr_ = Teuchos::rcp(new MeshManager_Rectangle<Real>(*parlist));
    // Finite element fields.
    Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr;
    if (basisOrder_ == 1) {
      basisPtr = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >);
    }
    else if (basisOrder_ == 2) {
      basisPtr = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >);
    }
    basisPtrs_.resize(1, Teuchos::null);
    basisPtrs_[0] = basisPtr;
    // DOF coordinate interface.
    Teuchos::RCP<Intrepid::DofCoordsInterface<Intrepid::FieldContainer<Real> > > coord_iface =
      Teuchos::rcp_dynamic_cast<Intrepid::DofCoordsInterface<Intrepid::FieldContainer<Real> > >(basisPtrs_[0]);
    // Degree-of-freedom manager.
    dofMgr_ = Teuchos::rcp(new DofManager<Real>(meshMgr_, basisPtrs_));
    // Retrieve total number of cells in the mesh.
    totalNumCells_ = meshMgr_->getNumCells();
    // Retrieve total number of degrees of freedom in the mesh.
    totalNumDofs_ = dofMgr_->getNumDofs();

    /****************************************************************************/
    /****************************************************************************/


    /****************************************************/
    /*** Build parallel communication infrastructure. ***/
    /****************************************************/

    // Partition the cells in the mesh.  We use a basic quasi-equinumerous partitioning,
    // where the remainder, if any, is assigned to the last processor.
    Teuchos::Array<int> myGlobIds_;
    Teuchos::Array<int> cellOffsets_(numProcs_, 0);
    int cellsPerProc = totalNumCells_ / numProcs_;
    numCells_ = cellsPerProc;
    switch(cellSplit) {
      case 0:
        if (myRank_ == 0) {  // remainder in the first
          numCells_ += totalNumCells_ % numProcs_;
        }
        for (int i=1; i<numProcs_; ++i) {
          cellOffsets_[i] = cellOffsets_[i-1] + cellsPerProc + (static_cast<int>(i==1))*(totalNumCells_ % numProcs_);
        }
        break;
      case 1:
        if (myRank_ == numProcs_-1) { // remainder in the last
          numCells_ += totalNumCells_ % numProcs_;
        }
        for (int i=1; i<numProcs_; ++i) {
          cellOffsets_[i] = cellOffsets_[i-1] + cellsPerProc;
        }
        break;
      case 2:
        if (myRank_ < (totalNumCells_%numProcs_)) { // spread remainder, starting from the first
          numCells_++;
        }
        for (int i=1; i<numProcs_; ++i) {
          cellOffsets_[i] = cellOffsets_[i-1] + cellsPerProc + (static_cast<int>(i-1<(totalNumCells_%numProcs_)));
        }
        break;
    }
    Intrepid::FieldContainer<int> &cellDofs = *(dofMgr_->getCellDofs());
    int numLocalDofs = cellDofs.dimension(1);
    *outStream << "Cell offsets across processors: " << cellOffsets_ << std::endl;
    for (int i=0; i<numCells_; ++i) {
      myCellIds_.push_back(cellOffsets_[myRank_]+i);
      for (int j=0; j<numLocalDofs; ++j) {
        myGlobIds_.push_back( cellDofs(cellOffsets_[myRank_]+i,j) );
      }
    }
    std::sort(myGlobIds_.begin(), myGlobIds_.end());
    myGlobIds_.erase( std::unique(myGlobIds_.begin(), myGlobIds_.end()), myGlobIds_.end() );

    // Build maps.
    myOverlapMap_ = Teuchos::rcp(new Tpetra::Map<>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                                                   myGlobIds_, 0, comm));
    //std::cout << std::endl << myOverlapMap_->getNodeElementList();
    /** One can also use the non-member function:
          myOverlapMap_ = Tpetra::createNonContigMap<int,int>(myGlobIds_, comm);
        to build the overlap map.
    **/
    myUniqueMap_ = Tpetra::createOneToOne<int,int>(myOverlapMap_);
    //std::cout << std::endl << myUniqueMap_->getNodeElementList() << std::endl;

    myBColumnMap_ = Teuchos::rcp(new Tpetra::Map<>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                                 myCellIds_, 0, comm));

    /****************************************************/
    /****************************************************/


    /****************************************************/
    /*** Set up local discretization data and arrays. ***/
    /****************************************************/

    // Retrieve some basic cell information.
    cellType_ = (basisPtrs_[0])->getBaseCellTopology();   // get the cell type from any basis
    spaceDim_ = cellType_.getDimension();                 // retrieve spatial dimension
    numNodesPerCell_ = cellType_.getNodeCount();          // retrieve number of nodes per cell

    // Cubature data.
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                                          // create cubature factory
    int cubDegree = 4;                                                                          // set cubature degree, e.g., 2
    Teuchos::RCP<Intrepid::Cubature<Real> > cellCub = cubFactory.create(cellType_, cubDegree);  // create default cubature
    numCubPoints_ = cellCub->getNumPoints();                                                    // retrieve number of cubature points

    int lfs = dofMgr_->getLocalFieldSize(0);

    // Discretization data. 
    cubPoints_            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCubPoints_, spaceDim_));
    cubWeights_           = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCubPoints_));
    cubPointsPhysical_    = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, numCubPoints_, spaceDim_));
    cellNodes_            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, numNodesPerCell_, spaceDim_));
    cellJac_              = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, numCubPoints_, spaceDim_, spaceDim_));
    cellJacInv_           = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, numCubPoints_, spaceDim_, spaceDim_));
    cellJacDet_           = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, numCubPoints_));
    cellWeightedMeasure_  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, numCubPoints_));
    valReference_         = Teuchos::rcp(new Intrepid::FieldContainer<Real>(lfs, numCubPoints_));  
    gradReference_        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(lfs, numCubPoints_, spaceDim_));  
    valPhysical_          = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs, numCubPoints_));
    gradPhysical_         = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs, numCubPoints_, spaceDim_));
    kappaGradPhysical_    = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs, numCubPoints_, spaceDim_));
    valPhysicalWeighted_  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs, numCubPoints_));
    gradPhysicalWeighted_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs, numCubPoints_, spaceDim_));
    gradgradMats_         = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs, lfs));
    valvalMats_           = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs, lfs));
    valMats_              = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, lfs, 1));
    onesVec_              = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, 1, numCubPoints_));
    kappa_                = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_, numCubPoints_));

    // Geometric definition of the cells in the mesh, based on the cell-to-node map and the domain partition.
    Intrepid::FieldContainer<Real> &nodes = *meshMgr_->getNodes();
    Intrepid::FieldContainer<int>  &ctn   = *meshMgr_->getCellToNodeMap();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numNodesPerCell_; ++j) {
        for (int k=0; k<spaceDim_; ++k) {
          (*cellNodes_)(i, j, k) = nodes(ctn(myCellIds_[i],j), k);
        }
      }
    }

    /****************************************************/
    /****************************************************/


    /****************************************************************/
    /*** Assemble cellwise contributions to vectors and matrices. ***/
    /****************************************************************/

    cellCub->getCubature(*cubPoints_, *cubWeights_);                                         // retrieve cubature points and weights
    (*basisPtrs_[0]).getValues(*gradReference_, *cubPoints_, Intrepid::OPERATOR_GRAD);       // evaluate grad operator at cubature points
    (*basisPtrs_[0]).getValues(*valReference_, *cubPoints_, Intrepid::OPERATOR_VALUE);       // evaluate value operator at cubature points

    Intrepid::CellTools<Real>::setJacobian(*cellJac_, *cubPoints_, *cellNodes_, cellType_);  // compute cell Jacobians
    Intrepid::CellTools<Real>::setJacobianInv(*cellJacInv_, *cellJac_);                      // compute inverses of cell Jacobians
    Intrepid::CellTools<Real>::setJacobianDet(*cellJacDet_, *cellJac_);                      // compute determinants of cell Jacobians

    Intrepid::FunctionSpaceTools::computeCellMeasure<Real>(*cellWeightedMeasure_,            // compute weighted cell measure
                                                           *cellJacDet_,
                                                           *cubWeights_);

    Intrepid::CellTools<Real>::mapToPhysicalFrame(*cubPointsPhysical_,                       // map reference cubature points to physical space
                                                  *cubPoints_,
                                                  *cellNodes_,
                                                  cellType_);

    Intrepid::FunctionSpaceTools::HGRADtransformGRAD<Real>(*gradPhysical_,                   // transform reference gradients into physical space
                                                           *cellJacInv_,
                                                           *gradReference_);
    Intrepid::FunctionSpaceTools::multiplyMeasure<Real>(*gradPhysicalWeighted_,              // multiply with weighted measure
                                                        *cellWeightedMeasure_,
                                                        *gradPhysical_);
    for (int i=0; i<numCells_; ++i) {                                                        // evaluate conductivity kappa at cubature points
      for (int j=0; j<numCubPoints_; ++j) {
        (*kappa_)(i, j) = funcKappa((*cubPointsPhysical_)(i, j, 0),
                                    (*cubPointsPhysical_)(i, j, 1));
      }
    }
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(*kappaGradPhysical_,         // multiply with conductivity kappa
                                                                *kappa_,
                                                                *gradPhysical_);
    Intrepid::FunctionSpaceTools::integrate<Real>(*gradgradMats_,                            // compute local grad.(kappa)grad (stiffness) matrices
                                                  *kappaGradPhysical_,
                                                  *gradPhysicalWeighted_,
                                                  Intrepid::COMP_CPP);

    Intrepid::FunctionSpaceTools::HGRADtransformVALUE<Real>(*valPhysical_,                   // transform reference values into physical space
                                                            *valReference_);
    Intrepid::FunctionSpaceTools::multiplyMeasure<Real>(*valPhysicalWeighted_,               // multiply with weighted measure
                                                        *cellWeightedMeasure_,
                                                        *valPhysical_);
    Intrepid::FunctionSpaceTools::integrate<Real>(*valvalMats_,                              // compute local val.val (mass) matrices
                                                  *valPhysical_,
                                                  *valPhysicalWeighted_,
                                                  Intrepid::COMP_CPP);

    for (int i=0; i<numCells_; ++i) {                                                        // fill vector of ones
      for (int j=0; j<numCubPoints_; ++j) {
        (*onesVec_)(i, 0, j) = static_cast<Real>(1);
      }
    }
    Intrepid::FunctionSpaceTools::integrate<Real>(*valMats_,                                 // compute local val.1 matrices
                                                  *valPhysicalWeighted_,
                                                  *onesVec_,
                                                  Intrepid::COMP_CPP);

    /****************************************************************/
    /****************************************************************/

    /****************************************/
    /*** Assemble global data structures. ***/
    /****************************************/

    // Assemble graphs.
    matAGraph_ = Teuchos::rcp(new Tpetra::CrsGraph<>(myUniqueMap_, 0));
    Teuchos::ArrayRCP<const int> cellDofsArrayRCP = cellDofs.getData();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matAGraph_->insertGlobalIndices(cellDofs(myCellIds_[i],j), cellDofsArrayRCP(myCellIds_[i]*numLocalDofs, numLocalDofs));
      }
    }
    matAGraph_->fillComplete();
    matBGraph_ = Teuchos::rcp(new Tpetra::CrsGraph<>(myUniqueMap_, myBColumnMap_, 0));
    Teuchos::ArrayRCP<const int> cellIdsArrayRCP = Teuchos::arcpFromArray(myCellIds_);
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matBGraph_->insertGlobalIndices(cellDofs(myCellIds_[i],j), cellIdsArrayRCP(i, 1));
      }
    }
    matBGraph_->fillComplete(myBColumnMap_, myUniqueMap_);

    // Assemble matrices.
    // Filter matrix = stiffness matrix plus mass matrix.
    matA_ = Tpetra::rcp(new Tpetra::CrsMatrix<>(matAGraph_));
    int numLocalMatEntries = numLocalDofs*numLocalDofs;
    Teuchos::ArrayRCP<const Real> gradgradArrayRCP = gradgradMats_->getData();
    Teuchos::ArrayRCP<const Real> valvalArrayRCP = valvalMats_->getData();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matA_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                   cellDofsArrayRCP(myCellIds_[i]*numLocalDofs, numLocalDofs),
                                   gradgradArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
        matA_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                   cellDofsArrayRCP(myCellIds_[i]*numLocalDofs, numLocalDofs),
                                   valvalArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
      }
    }
    matA_->fillComplete();
    // B matrix.
    matB_ = Tpetra::rcp(new Tpetra::CrsMatrix<>(matBGraph_));
    Teuchos::ArrayRCP<const Real> valArrayRCP = valMats_->getData();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matB_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                   cellIdsArrayRCP(i, 1),
                                   valArrayRCP(i*numLocalDofs+j, 1));
      }
    }
    matB_->fillComplete();

    // Compute cell colume vector.
    computeCellVolumes();

    // Create matrix transposes.
    Tpetra::RowMatrixTransposer<> transposerB(matB_);
    matB_trans_ = transposerB.createTranspose();

    /*********************************/
    /*** Construct solver objects. ***/
    /*********************************/

    // Construct solver using Amesos2 factory.
    try{
      solverA_ = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("KLU2", matA_);
    } catch (std::invalid_argument e) {
      std::cout << e.what() << std::endl;
    }
    solverA_->numericFactorization();

    /****************************************/
    /****************************************/

    //outputTpetraData();

  }

  void apply(Teuchos::RCP<Tpetra::MultiVector<> > & Fx, const Teuchos::RCP<const Tpetra::MultiVector<> > & x) {
    if (enableFilter_) {
      Teuchos::RCP<Tpetra::MultiVector<> > Bx = Teuchos::rcp(new Tpetra::MultiVector<>(matB_->getRangeMap(), 1));
      Teuchos::RCP<Tpetra::MultiVector<> > AiBx = Teuchos::rcp(new Tpetra::MultiVector<>(matA_->getDomainMap(), 1));
      Teuchos::RCP<Tpetra::MultiVector<> > Fx_unscaled = Teuchos::rcp(new Tpetra::MultiVector<>(matB_trans_->getRangeMap(), 1));
      matB_->apply(*x, *Bx);
      solverA_->setX(AiBx);
      solverA_->setB(Bx);
      solverA_->solve();
      //outputTpetraVector(AiBx, "density_nodal.txt");
      matB_trans_->apply(*AiBx, *Fx_unscaled);
      Teuchos::RCP<Tpetra::MultiVector<> > vecInvCellVolumes = Teuchos::rcp(new Tpetra::MultiVector<>(matB_->getDomainMap(), 1));
      vecInvCellVolumes->reciprocal(*vecCellVolumes_);
      Fx->elementWiseMultiply(1.0, *(vecInvCellVolumes->getVector(0)), *Fx_unscaled, 0.0);
    }
    else {
      Fx->update(1.0, *x, 0.0);
    }
  }

  Teuchos::RCP<Tpetra::CrsMatrix<> > getMatA() const {
    return matA_;
  }

  Teuchos::RCP<Tpetra::CrsMatrix<> > getMatB(const bool &transpose = false) const {
    if (transpose) {
      return matB_trans_;
    }
    else {
      return matB_;
    }
  }

  Teuchos::RCP<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > getSolver() const {
    return solverA_;
  }

  Real funcKappa(const Real &x1, const Real &x2) const {
    return lengthScale_;
  }

  void computeCellVolumes() {
    for (int i=0; i<numCells_; ++i){
      Real temp = 0.0;
      for(int j=0; j<numCubPoints_; ++j){
        temp += (*cellWeightedMeasure_)(i, j);
      }
      myCellVolumes_.push_back(temp);
    }

    vecCellVolumes_ = Teuchos::rcp(new Tpetra::MultiVector<>(matB_->getDomainMap(), 1, true));
    for (int i=0; i<numCells_; ++i){
        vecCellVolumes_->replaceGlobalValue(myCellIds_[i], 0, myCellVolumes_[i]);
    }
  }

  void outputTpetraData() const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> >   matWriter;
    matWriter.writeSparseFile("filter_mat", matA_, true);
    matWriter.writeSparseFile("projection_mat", matB_, true);
    matWriter.writeSparseFile("projection_trans_mat", matB_trans_, true);
  }

  void outputTpetraVector(const Teuchos::RCP<const Tpetra::MultiVector<> > &vec,
                          const std::string &filename) const {
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<> > vecWriter;
    vecWriter.writeDenseFile(filename, vec);
  }

}; // class StefanBoltzmannData

#endif
