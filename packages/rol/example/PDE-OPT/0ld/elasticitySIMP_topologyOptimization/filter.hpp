// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  filter.hpp
    \brief Density filtering based on solving a PDE.
*/

#ifndef ROL_DENSITY_FILTER_H
#define ROL_DENSITY_FILTER_H

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"

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
  using GO = typename Tpetra::Map<>::global_ordinal_type;

  ROL::Ptr<MeshManager<Real> > meshMgr_;
  ROL::Ptr<DofManager<Real> >  dofMgr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;

  ROL::Ptr<const Teuchos::Comm<int> > commPtr_;
  int myRank_;
  int numProcs_;

  int  basisOrder_;

  ROL::Ptr<const Tpetra::Map<> >    myOverlapMap_;
  ROL::Ptr<const Tpetra::Map<> >    myUniqueMap_;
  ROL::Ptr<const Tpetra::Map<> >    myBColumnMap_;
  ROL::Ptr<Tpetra::CrsGraph<> >     matAGraph_;
  ROL::Ptr<Tpetra::CrsGraph<> >     matBGraph_;
  ROL::Ptr<Tpetra::CrsMatrix<> >    matA_;
  ROL::Ptr<Tpetra::CrsMatrix<> >    matB_;
  ROL::Ptr<Tpetra::CrsMatrix<> >    matB_trans_;
  ROL::Ptr<Tpetra::MultiVector<> >  vecCellVolumes_;

  Teuchos::Array<GO> myCellIds_;
  Teuchos::Array<Real> myCellVolumes_;

  ROL::Ptr<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > solverA_;

  shards::CellTopology cellType_;
  int spaceDim_;
  int numNodesPerCell_;
  int numCubPoints_;

  GO totalNumCells_;
  GO totalNumDofs_;
  GO numCells_;

  ROL::Ptr<Intrepid::FieldContainer<Real> > cubPoints_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > cubWeights_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > cellNodes_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > cellJac_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > cellJacInv_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > cellJacDet_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > cellWeightedMeasure_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > valReference_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > gradReference_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > valPhysical_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > gradPhysical_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > kappaGradPhysical_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > valPhysicalWeighted_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > gradPhysicalWeighted_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > gradgradMats_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > valvalMats_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > valMats_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > onesVec_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > cubPointsPhysical_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > kappa_;

  Real lengthScale_;
  bool enableFilter_;

public:

  DensityFilter(const ROL::Ptr<const Teuchos::Comm<int> > &comm,
                const Teuchos::RCP<Teuchos::ParameterList> &parlist,
                const ROL::Ptr<std::ostream> &outStream) {

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
    meshMgr_ = ROL::makePtr<MeshManager_Rectangle<Real>>(*parlist);
    // Finite element fields.
    ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr;
    if (basisOrder_ == 1) {
      basisPtr = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
    }
    else if (basisOrder_ == 2) {
      basisPtr = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >>();
    }
    basisPtrs_.resize(1, ROL::nullPtr);
    basisPtrs_[0] = basisPtr;
    // DOF coordinate interface.
    ROL::Ptr<Intrepid::DofCoordsInterface<Intrepid::FieldContainer<Real> > > coord_iface =
      ROL::dynamicPtrCast<Intrepid::DofCoordsInterface<Intrepid::FieldContainer<Real> > >(basisPtrs_[0]);
    // Degree-of-freedom manager.
    dofMgr_ = ROL::makePtr<DofManager<Real>>(meshMgr_, basisPtrs_);
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
    Teuchos::Array<GO> myGlobIds_;
    Teuchos::Array<GO> cellOffsets_(numProcs_, 0);
    GO cellsPerProc = totalNumCells_ / numProcs_;
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
    Intrepid::FieldContainer<GO> &cellDofs = *(dofMgr_->getCellDofs());
    int numLocalDofs = cellDofs.dimension(1);
    *outStream << "Cell offsets across processors: " << cellOffsets_ << std::endl;
    for (GO i=0; i<numCells_; ++i) {
      myCellIds_.push_back(cellOffsets_[myRank_]+i);
      for (int j=0; j<numLocalDofs; ++j) {
        myGlobIds_.push_back( cellDofs(cellOffsets_[myRank_]+i,j) );
      }
    }
    std::sort(myGlobIds_.begin(), myGlobIds_.end());
    myGlobIds_.erase( std::unique(myGlobIds_.begin(), myGlobIds_.end()), myGlobIds_.end() );

    // Build maps.
    myOverlapMap_ = ROL::makePtr<Tpetra::Map<>>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                                                   myGlobIds_, 0, comm);
    //std::cout << std::endl << myOverlapMap_->getLocalElementList();
    /** One can also use the non-member function:
          myOverlapMap_ = Tpetra::createNonContigMap<int,int>(myGlobIds_, comm);
        to build the overlap map.
    **/
    myUniqueMap_ = Tpetra::createOneToOne(myOverlapMap_);
    //std::cout << std::endl << myUniqueMap_->getLocalElementList() << std::endl;

    myBColumnMap_ = ROL::makePtr<Tpetra::Map<>>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                                 myCellIds_, 0, comm);

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
    ROL::Ptr<Intrepid::Cubature<Real> > cellCub = cubFactory.create(cellType_, cubDegree);  // create default cubature
    numCubPoints_ = cellCub->getNumPoints();                                                    // retrieve number of cubature points

    int lfs = dofMgr_->getLocalFieldSize(0);

    // Discretization data. 
    cubPoints_            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCubPoints_, spaceDim_);
    cubWeights_           = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCubPoints_);
    cubPointsPhysical_    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_, spaceDim_);
    cellNodes_            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numNodesPerCell_, spaceDim_);
    cellJac_              = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_, spaceDim_, spaceDim_);
    cellJacInv_           = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_, spaceDim_, spaceDim_);
    cellJacDet_           = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_);
    cellWeightedMeasure_  = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_);
    valReference_         = ROL::makePtr<Intrepid::FieldContainer<Real>>(lfs, numCubPoints_);  
    gradReference_        = ROL::makePtr<Intrepid::FieldContainer<Real>>(lfs, numCubPoints_, spaceDim_);  
    valPhysical_          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, numCubPoints_);
    gradPhysical_         = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, numCubPoints_, spaceDim_);
    kappaGradPhysical_    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, numCubPoints_, spaceDim_);
    valPhysicalWeighted_  = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, numCubPoints_);
    gradPhysicalWeighted_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, numCubPoints_, spaceDim_);
    gradgradMats_         = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, lfs);
    valvalMats_           = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, lfs);
    valMats_              = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, 1);
    onesVec_              = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, 1, numCubPoints_);
    kappa_                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_);

    // Geometric definition of the cells in the mesh, based on the cell-to-node map and the domain partition.
    Intrepid::FieldContainer<Real> &nodes = *meshMgr_->getNodes();
    Intrepid::FieldContainer<GO>  &ctn   = *meshMgr_->getCellToNodeMap();
    for (GO i=0; i<numCells_; ++i) {
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
    matAGraph_ = ROL::makePtr<Tpetra::CrsGraph<>>(myUniqueMap_, 0);
    Teuchos::ArrayRCP<const GO> cellDofsArrayRCP = cellDofs.getData();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matAGraph_->insertGlobalIndices(cellDofs(myCellIds_[i],j), cellDofsArrayRCP(myCellIds_[i]*numLocalDofs, numLocalDofs));
      }
    }
    matAGraph_->fillComplete();
    matBGraph_ = ROL::makePtr<Tpetra::CrsGraph<>>(myUniqueMap_, myBColumnMap_, 0);
    Teuchos::ArrayRCP<const GO> cellIdsArrayRCP = Teuchos::arcpFromArray(myCellIds_);
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matBGraph_->insertGlobalIndices(cellDofs(myCellIds_[i],j), cellIdsArrayRCP(i, 1));
      }
    }
    matBGraph_->fillComplete(myBColumnMap_, myUniqueMap_);

    // Assemble matrices.
    // Filter matrix = stiffness matrix plus mass matrix.
    matA_ = ROL::makePtr<Tpetra::CrsMatrix<>>(matAGraph_);
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
    matB_ = ROL::makePtr<Tpetra::CrsMatrix<>>(matBGraph_);
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
    } catch (std::invalid_argument& e) {
      std::cout << e.what() << std::endl;
    }
    solverA_->numericFactorization();

    /****************************************/
    /****************************************/

    //outputTpetraData();

  }

  void apply(ROL::Ptr<Tpetra::MultiVector<> > & Fx, const ROL::Ptr<const Tpetra::MultiVector<> > & x) {
    if (enableFilter_) {
      ROL::Ptr<Tpetra::MultiVector<> > Bx = ROL::makePtr<Tpetra::MultiVector<>>(matB_->getRangeMap(), 1);
      ROL::Ptr<Tpetra::MultiVector<> > AiBx = ROL::makePtr<Tpetra::MultiVector<>>(matA_->getDomainMap(), 1);
      ROL::Ptr<Tpetra::MultiVector<> > Fx_unscaled = ROL::makePtr<Tpetra::MultiVector<>>(matB_trans_->getRangeMap(), 1);
      matB_->apply(*x, *Bx);
      solverA_->setX(AiBx);
      solverA_->setB(Bx);
      solverA_->solve();
      //outputTpetraVector(AiBx, "density_nodal.txt");
      matB_trans_->apply(*AiBx, *Fx_unscaled);
      ROL::Ptr<Tpetra::MultiVector<> > vecInvCellVolumes = ROL::makePtr<Tpetra::MultiVector<>>(matB_->getDomainMap(), 1);
      vecInvCellVolumes->reciprocal(*vecCellVolumes_);
      Fx->elementWiseMultiply(1.0, *(vecInvCellVolumes->getVector(0)), *Fx_unscaled, 0.0);
    }
    else {
      Fx->update(1.0, *x, 0.0);
    }
  }

  ROL::Ptr<Tpetra::CrsMatrix<> > getMatA() const {
    return matA_;
  }

  ROL::Ptr<Tpetra::CrsMatrix<> > getMatB(const bool &transpose = false) const {
    if (transpose) {
      return matB_trans_;
    }
    else {
      return matB_;
    }
  }

  ROL::Ptr<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > getSolver() const {
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

    vecCellVolumes_ = ROL::makePtr<Tpetra::MultiVector<>>(matB_->getDomainMap(), 1, true);
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

  void outputTpetraVector(const ROL::Ptr<const Tpetra::MultiVector<> > &vec,
                          const std::string &filename) const {
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<> > vecWriter;
    vecWriter.writeDenseFile(filename, vec);
  }

}; // class StefanBoltzmannData

#endif
