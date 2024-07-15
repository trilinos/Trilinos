// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  data.hpp
    \brief Generates and manages data for the Poisson example, including
           all mesh and discretization data, matrices, etc.
*/

#ifndef ROL_PDEOPT_STOCH_ADV_DIFF_DATA_H
#define ROL_PDEOPT_STOCH_ADV_DIFF_DATA_H

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
class PoissonData {

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
  ROL::Ptr<Tpetra::CrsGraph<> >     matGraph_;
  ROL::Ptr<Tpetra::CrsMatrix<> >    matA_;
  ROL::Ptr<Tpetra::CrsMatrix<> >    matA_dirichlet_;
  ROL::Ptr<Tpetra::CrsMatrix<> >    matA_dirichlet_trans_;
  ROL::Ptr<Tpetra::CrsMatrix<> >    matM_;
//  ROL::Ptr<Tpetra::CrsMatrix<> >    matB_;
//  ROL::Ptr<Tpetra::CrsMatrix<> >    matB_dirichlet_;
//  ROL::Ptr<Tpetra::CrsMatrix<> >    matB_dirichlet_trans_;
  ROL::Ptr<Tpetra::MultiVector<> >  vecUd_;
  ROL::Ptr<Tpetra::MultiVector<> >  vecF_;
  ROL::Ptr<Tpetra::MultiVector<> >  vecF_overlap_;
  ROL::Ptr<Tpetra::MultiVector<> >  vecF_dirichlet_;
  ROL::Ptr<Tpetra::MultiVector<> >  vecWeights_;

  ROL::Ptr<Tpetra::MultiVector<> > matB_;
  ROL::Ptr<Tpetra::MultiVector<> > matB_dirichlet_;
  ROL::Ptr<Tpetra::MultiVector<> > matB_overlap_;

  Teuchos::Array<GO> myCellIds_;
  Teuchos::Array<GO> myDirichletDofs_;

  ROL::Ptr<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > solverA_;
  ROL::Ptr<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > solverA_trans_;

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
  ROL::Ptr<Intrepid::FieldContainer<Real> > advGradPhysical_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > funcValPhysical_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > valPhysicalWeighted_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > gradPhysicalWeighted_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > gradgradMats_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > valvalMats_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > valfuncvalMats_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > advgradvalMats_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > pdeMats_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > cubPointsPhysical_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > kappa_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > cfunc_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > adv_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > dataF_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > datavalVecF_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > dofPoints_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > dofPointsPhysical_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > dataUd_;

  ROL::Ptr<Intrepid::FieldContainer<Real> > dataB_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > datavalB_;
  std::vector<Intrepid::FieldContainer<Real> > vectorDataValB_;

public:

  PoissonData(const ROL::Ptr<const Teuchos::Comm<int> > &comm,
              const Teuchos::RCP<Teuchos::ParameterList> &parlist,
              const ROL::Ptr<std::ostream> &outStream) {
    std::vector<Real> param(37,1);

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
    basisOrder_ = parlist->sublist("Problem").get("Order of FE discretization", 1);
    int cellSplit = parlist->sublist("Geometry").get("Partition type", 1);
    /*************************************/
    /*************************************/

    /****************************************************************************/
    /*** Initialize mesh / finite element fields / degree-of-freedom manager. ***/
    /****************************************************************************/

    // Mesh manager.
    meshMgr_ = ROL::makePtr<MeshManager_Rectangle<Real>>(*parlist);
    printMeshData(*outStream);
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
    for (int i=0; i<numCells_; ++i) {
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
    dofPoints_            = ROL::makePtr<Intrepid::FieldContainer<Real>>(lfs, spaceDim_);
    dofPointsPhysical_    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, spaceDim_);
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
    advGradPhysical_      = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, numCubPoints_);
    valPhysicalWeighted_  = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, numCubPoints_);
    funcValPhysical_      = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, numCubPoints_);
    gradPhysicalWeighted_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, numCubPoints_, spaceDim_);
    gradgradMats_         = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, lfs);
    valvalMats_           = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, lfs);
    valfuncvalMats_       = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, lfs);
    advgradvalMats_       = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, lfs);
    pdeMats_              = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, lfs);
    kappa_                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_);
    cfunc_                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_);
    adv_                  = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_, spaceDim_);
    dataF_                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_);
    datavalVecF_          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs);
    dataUd_               = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs);

    dataB_                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_);
    datavalB_             = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs);

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
    Intrepid::CellTools<Real>::mapToPhysicalFrame(*cubPointsPhysical_,                       // map reference cubature points to physical space
                                                  *cubPoints_,
                                                  *cellNodes_,
                                                  cellType_);
    (*basisPtrs_[0]).getValues(*gradReference_, *cubPoints_, Intrepid::OPERATOR_GRAD);       // evaluate grad operator at cubature points
    (*basisPtrs_[0]).getValues(*valReference_, *cubPoints_, Intrepid::OPERATOR_VALUE);       // evaluate value operator at cubature points

    Intrepid::CellTools<Real>::setJacobian(*cellJac_, *cubPoints_, *cellNodes_, cellType_);  // compute cell Jacobians
    Intrepid::CellTools<Real>::setJacobianInv(*cellJacInv_, *cellJac_);                      // compute inverses of cell Jacobians
    Intrepid::CellTools<Real>::setJacobianDet(*cellJacDet_, *cellJac_);                      // compute determinants of cell Jacobians

    Intrepid::FunctionSpaceTools::computeCellMeasure<Real>(*cellWeightedMeasure_,            // compute weighted cell measure
                                                           *cellJacDet_,
                                                           *cubWeights_);

    Intrepid::FunctionSpaceTools::HGRADtransformGRAD<Real>(*gradPhysical_,                   // transform reference gradients into physical space
                                                           *cellJacInv_,
                                                           *gradReference_);
    Intrepid::FunctionSpaceTools::multiplyMeasure<Real>(*gradPhysicalWeighted_,              // multiply with weighted measure
                                                        *cellWeightedMeasure_,
                                                        *gradPhysical_);
    for (GO i=0; i<numCells_; ++i) {                                                        // evaluate conductivity kappa at cubature points
      for (int j=0; j<numCubPoints_; ++j) {
        (*kappa_)(i, j) = funcKappa((*cubPointsPhysical_)(i, j, 0),
                                    (*cubPointsPhysical_)(i, j, 1),
                                    param);
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

//    for (int i=0; i<numCells_; ++i) {                                                        // evaluate conductivity kappa at cubature points
//      for (int j=0; j<numCubPoints_; ++j) {
//        (*cfunc_)(i, j) = funcControl((*cubPointsPhysical_)(i, j, 0),
//                                      (*cubPointsPhysical_)(i, j, 1));
//      }
//    }
//    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*funcValPhysical_,
//                                                                *cfunc_,
//                                                                *valPhysical_);
//    Intrepid::FunctionSpaceTools::integrate<Real>(*valfuncvalMats_,                          // compute local val.func.val (mass) matrices
//                                                  *funcValPhysical_,
//                                                  *valPhysicalWeighted_,
//                                                  Intrepid::COMP_CPP);
//
    for (GO i=0; i<numCells_; ++i) {                                                        // evaluate advection adv at cubature points
      for (int j=0; j<numCubPoints_; ++j) {
        for (int k=0; k<spaceDim_; ++k) {
          std::vector<Real> adv(spaceDim_, 0);
          funcAdv(adv, (*cubPointsPhysical_)(i, j, 0),
                       (*cubPointsPhysical_)(i, j, 1),
                       param);
          (*adv_)(i, j, k) = adv[k];
        }
      }
    }
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*advGradPhysical_,              // multiply with gradient in physical space
                                                             *adv_,
                                                             *gradPhysical_);
    // IMPORTANT: This is a nonsymmetric form; watch the field order in the integrand!!!
    Intrepid::FunctionSpaceTools::integrate<Real>(*advgradvalMats_,                          // compute local (adv.grad).val (advection) matrices
                                                  *valPhysicalWeighted_,
                                                  *advGradPhysical_,
                                                  Intrepid::COMP_CPP);

    Intrepid::RealSpaceTools<Real>::add(*pdeMats_, *gradgradMats_, *advgradvalMats_);        // add diffusion and advection matrices
    for (GO i=0; i<numCells_; ++i) {                                                        // evaluate functions at these points
      for (int j=0; j<numCubPoints_; ++j) {
        (*dataF_)(i, j) = funcRHS((*cubPointsPhysical_)(i, j, 0),
                                  (*cubPointsPhysical_)(i, j, 1),
                                  param);
      }
    }
    Intrepid::FunctionSpaceTools::integrate<Real>(*datavalVecF_,                             // compute local data.val vectors for RHS F
                                                  *dataF_,
                                                  *valPhysicalWeighted_,
                                                  Intrepid::COMP_CPP);

    coord_iface->getDofCoords(*dofPoints_);                                                  // get coordinates of DOFs in reference cell
    Intrepid::CellTools<Real>::mapToPhysicalFrame(*dofPointsPhysical_,                       // map reference DOF locations to physical space
                                                  *dofPoints_,
                                                  *cellNodes_,
                                                  cellType_);
    for (GO i=0; i<numCells_; ++i) {                                                        // evaluate functions at these points
      for (int j=0; j<lfs; ++j) {
        (*dataUd_)(i, j) = funcTarget((*dofPointsPhysical_)(i, j, 0), (*dofPointsPhysical_)(i, j, 1));
      }
    }

    // Control Operator
    const int nz = 9;
    const std::vector<Real> x = {0.25, 0.50, 0.75, 0.25, 0.50, 0.75, 0.25, 0.50, 0.75};
    const std::vector<Real> y = {0.25, 0.25, 0.25, 0.50, 0.50, 0.50, 0.75, 0.75, 0.75};
    for (int k=0; k<nz; ++k) {
      for (GO i=0; i<numCells_; ++i) {
        for (int j=0; j<numCubPoints_; ++j) {
          (*dataB_)(i, j) = funcControl((*cubPointsPhysical_)(i, j, 0),
                                        (*cubPointsPhysical_)(i, j, 1),
                                        x[k],y[k]);
        }
      }
      Intrepid::FunctionSpaceTools::integrate<Real>(*datavalB_,
                                                    *dataB_,
                                                    *valPhysicalWeighted_,
                                                    Intrepid::COMP_CPP);
      vectorDataValB_.push_back(*datavalB_);
    }


    /****************************************************************/
    /****************************************************************/

    /****************************************/
    /*** Assemble global data structures. ***/
    /****************************************/

    // Assemble graph.
    Teuchos::ArrayRCP<const GO> cellDofsArrayRCP = cellDofs.getData();
    Teuchos::Array<size_t> graphEntriesPerRow(myUniqueMap_->getLocalNumElements());
    for (GO i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        GO gid = cellDofs(myCellIds_[i],j);
        if(myUniqueMap_->isNodeGlobalElement(gid))
          graphEntriesPerRow[myUniqueMap_->getLocalElement(gid)] += numLocalDofs;
      }
    }
    matGraph_ = ROL::makePtr<Tpetra::CrsGraph<>>(myUniqueMap_, graphEntriesPerRow());
    for (GO i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matGraph_->insertGlobalIndices(cellDofs(myCellIds_[i],j), cellDofsArrayRCP(myCellIds_[i]*numLocalDofs, numLocalDofs));
      }
    }
    matGraph_->fillComplete();

    // Assemble matrices.
    // PDE matrix A.
    matA_ = ROL::makePtr<Tpetra::CrsMatrix<>>(matGraph_);
    int numLocalMatEntries = numLocalDofs*numLocalDofs;
    Teuchos::ArrayRCP<const Real> pdeArrayRCP = pdeMats_->getData();
    for (GO i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matA_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                   cellDofsArrayRCP(myCellIds_[i]*numLocalDofs, numLocalDofs),
                                   pdeArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
      }
    }
    matA_->fillComplete();

    // Mass matrix B.
    matB_ = ROL::makePtr<Tpetra::MultiVector<>>(matA_->getRangeMap(), nz, true);
    matB_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapMap_, nz, true);
    for (int k=0; k<nz; ++k) {// assembly over number of control sources
      for (GO i=0; i<numCells_; ++i) {// assembly on the overlap map
        for (int j=0; j<numLocalDofs; ++j) {
          matB_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                            k,
                                            vectorDataValB_[k][i*numLocalDofs+j]);
        }
      }
    }
    Tpetra::Export<> exporterB(matB_overlap_->getMap(), matB_->getMap()); // redistribution:
    matB_->doExport(*matB_overlap_, exporterB, Tpetra::ADD); // from the overlap map to the unique map
//    matB_ = ROL::makePtr<Tpetra::CrsMatrix<>>(matGraph_);
//    Teuchos::ArrayRCP<const Real> valfuncvalArrayRCP = valfuncvalMats_->getData();
//    for (int i=0; i<numCells_; ++i) {
//      for (int j=0; j<numLocalDofs; ++j) {
//        matB_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
//                                   cellDofsArrayRCP(myCellIds_[i]*numLocalDofs, numLocalDofs),
//                                   valfuncvalArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
//      }
//    }
//    matB_->fillComplete();

    // Mass matrix M.
    matM_ = ROL::makePtr<Tpetra::CrsMatrix<>>(matGraph_);
    Teuchos::ArrayRCP<const Real> valvalArrayRCP = valvalMats_->getData();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matM_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                   cellDofsArrayRCP(myCellIds_[i]*numLocalDofs, numLocalDofs),
                                   valvalArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
      }
    }
    matM_->fillComplete();

    // Assemble vectors.
    // vecF_ requires assembly using vecF_overlap_ and redistribution
    vecF_ = ROL::makePtr<Tpetra::MultiVector<>>(matA_->getRangeMap(), 1, true);
    vecF_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapMap_, 1, true);
    for (GO i=0; i<numCells_; ++i) {                                                 // assembly on the overlap map
      for (int j=0; j<numLocalDofs; ++j) {
        vecF_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                          0,
                                          (*datavalVecF_)[i*numLocalDofs+j]);
      }
    }
    Tpetra::Export<> exporter(vecF_overlap_->getMap(), vecF_->getMap());              // redistribution:
    vecF_->doExport(*vecF_overlap_, exporter, Tpetra::ADD);                           // from the overlap map to the unique map
    // vecUd_ does not require assembly
    vecUd_ = ROL::makePtr<Tpetra::MultiVector<>>(matA_->getDomainMap(), 1, true);
    for (GO i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        if (vecUd_->getMap()->isNodeGlobalElement(cellDofs(myCellIds_[i],j))) {
          vecUd_->replaceGlobalValue(cellDofs(myCellIds_[i],j),
                                     0,
                                     (*dataUd_)[i*numLocalDofs+j]);
        }
      }
    }
    // vecWeights
    vecWeights_ = ROL::makePtr<Tpetra::MultiVector<>>(matA_->getDomainMap(), 1, true);
    vecWeights_->putScalar(1.0);
    /*vecWeights_->putScalar(0.0);
    Real mask1_x1 = 0.30, mask1_x2 = 0.90, mask1_y1 = 0.20, mask1_y2 = 0.80;
    Real mask2_x1 = 0.15, mask2_x2 = 0.25, mask2_y1 = 0.55, mask2_y2 = 0.65;
    Real mask3_x1 = 0.45, mask3_x2 = 0.55, mask3_y1 = 0.05, mask3_y2 = 0.15;
    for (int i=0; i<vecWeights_->getMap()->getMaxLocalIndex(); ++i) {
      vecWeights_->replaceLocalValue(i, 0, !(i%5) ? ((static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX))) : static_cast<Real>(0) );
      //vecWeights_->replaceLocalValue(i, 0, !(i%5) ? 1 : static_cast<Real>(0) );
    }
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numCubPoints_; ++j) {
        if ( !( ((*cubPointsPhysical_)(i, j, 0) >= mask1_x1) &&
                ((*cubPointsPhysical_)(i, j, 0) <= mask1_x2) &&
                ((*cubPointsPhysical_)(i, j, 1) >= mask1_y1) &&
                ((*cubPointsPhysical_)(i, j, 1) <= mask1_y2) ) &&
             !( ((*cubPointsPhysical_)(i, j, 0) >= mask2_x1) &&
                ((*cubPointsPhysical_)(i, j, 0) <= mask2_x2) &&
                ((*cubPointsPhysical_)(i, j, 1) >= mask2_y1) &&
                ((*cubPointsPhysical_)(i, j, 1) <= mask2_y2) ) &&
             !( ((*cubPointsPhysical_)(i, j, 0) >= mask3_x1) &&
                ((*cubPointsPhysical_)(i, j, 0) <= mask3_x2) &&
                ((*cubPointsPhysical_)(i, j, 1) >= mask3_y1) &&
                ((*cubPointsPhysical_)(i, j, 1) <= mask3_y2) ) ) {
          vecWeights_->replaceGlobalValue(cellDofs(myCellIds_[i],j), 0, 0);
        }
      }
    }*/


    // Apply Dirichlet conditions.
    // Stiffness matrix with Dirichlet conditions:
    //  AD = [ A11  A12 ]  where A = [ A11 A12 ]
    //       [  0    I  ]            [ A21 A22 ]
    // Mass matrix with Dirichlet conditions:
    //  MD = [ M11  M12 ]  where M = [ M11 M12 ]
    //       [  0    0  ]            [ M21 M22 ]
    // Vector F with Dirichlet conditions G:
    //  FD = [ F1 ]  where F = [ F1 ]
    //       [ G  ]            [ F2 ]
    matA_dirichlet_ = ROL::makePtr<Tpetra::CrsMatrix<>>(*matA_, Teuchos::DataAccess::Copy);
    matB_dirichlet_ = ROL::makePtr<Tpetra::MultiVector<>>(matA_->getRangeMap(), nz, true);
    vecF_dirichlet_ = ROL::makePtr<Tpetra::MultiVector<>>(matA_->getRangeMap(), 1, true);
    Tpetra::deep_copy(*matB_dirichlet_, *matB_);
    Tpetra::deep_copy(*vecF_dirichlet_, *vecF_);
    ROL::Ptr<std::vector<std::vector<Intrepid::FieldContainer<GO> > > > dirichletSideSets = meshMgr_->getSideSets();
    std::vector<std::vector<Intrepid::FieldContainer<GO> > > &dss = *dirichletSideSets;
    Teuchos::Array<GO> mySortedCellIds_(myCellIds_);
    std::sort(mySortedCellIds_.begin(), mySortedCellIds_.end());
    mySortedCellIds_.erase( std::unique(mySortedCellIds_.begin(), mySortedCellIds_.end()), mySortedCellIds_.end() );
    std::vector<Teuchos::Array<GO> > myDirichletCellIds_(dss[0].size());
    for (int i=0; i<static_cast<int>(dss[0].size()); ++i) {
      if (i != 1) {  // exclude right boundary of the square domain (pure Neumann)
        for (int j=0; j<dss[0][i].dimension(0); ++j) {
          if (std::binary_search(mySortedCellIds_.begin(), mySortedCellIds_.end(), dss[0][i](j))) {
            myDirichletCellIds_[i].push_back(dss[0][i](j));
          }
        }
      }
    }
    Intrepid::FieldContainer<GO>  &cte      = *(meshMgr_->getCellToEdgeMap());
    Intrepid::FieldContainer<GO>  &nodeDofs = *(dofMgr_->getNodeDofs());
    Intrepid::FieldContainer<GO>  &edgeDofs = *(dofMgr_->getEdgeDofs());
    std::vector<std::vector<int> > dofTags = (basisPtrs_[0])->getAllDofTags();
    int numDofsPerNode = 0;
    int numDofsPerEdge = 0;
    for (int j=0; j<(basisPtrs_[0])->getCardinality(); ++j) {
      if (dofTags[j][0] == 0) {
        numDofsPerNode = dofTags[j][3];
      }
      if (dofTags[j][0] == 1) {
        numDofsPerEdge = dofTags[j][3];
      }
    }
    for (int i=0; i<static_cast<int>(myDirichletCellIds_.size()); ++i) {
      for (int j=0; j<myDirichletCellIds_[i].size(); ++j) {
        for (int k=0; k<numDofsPerNode; ++k) {
          const CellTopologyData * ctd = cellType_.getCellTopologyData();
          Teuchos::ArrayView<unsigned> locNodes(const_cast<unsigned *>(ctd->subcell[spaceDim_-1][i].node), cellType_.getVertexCount(spaceDim_-1, i));
          for (int l=0; l<static_cast<int>(cellType_.getVertexCount(spaceDim_-1, i)); ++l) {
            myDirichletDofs_.push_back(nodeDofs(ctn(myDirichletCellIds_[i][j], locNodes[l]), k));
          }
        }
        for (int k=0; k<numDofsPerEdge; ++k) {
           myDirichletDofs_.push_back(edgeDofs(cte(myDirichletCellIds_[i][j], i), k));
        }
      }
    }
    std::sort(myDirichletDofs_.begin(), myDirichletDofs_.end());
    myDirichletDofs_.erase( std::unique(myDirichletDofs_.begin(), myDirichletDofs_.end()), myDirichletDofs_.end() );
    matA_dirichlet_->resumeFill();
    for (int i=0; i<myDirichletDofs_.size(); ++i) {
      if (myUniqueMap_->isNodeGlobalElement(myDirichletDofs_[i])) {
        size_t numRowEntries = matA_dirichlet_->getNumEntriesInGlobalRow(myDirichletDofs_[i]);
        Teuchos::Array<GO> indices(numRowEntries, 0);    
        Teuchos::Array<Real> values(numRowEntries, 0);
        Teuchos::Array<Real> canonicalValues(numRowEntries, 0);    
        Teuchos::Array<Real> zeroValues(numRowEntries, 0);    
        matA_dirichlet_->getGlobalRowCopy(myDirichletDofs_[i], indices, values, numRowEntries);
        for (int j=0; j<indices.size(); ++j) {
          if (myDirichletDofs_[i] == indices[j]) {
            canonicalValues[j] = 1.0;
          }
        }
        matA_dirichlet_->replaceGlobalValues(myDirichletDofs_[i], indices, canonicalValues);
        for (int k=0; k<nz; k++) {
          matB_dirichlet_->replaceGlobalValue(myDirichletDofs_[i], k, 0);
        }
        vecF_dirichlet_->replaceGlobalValue(myDirichletDofs_[i], 0, 0);
      }
    }
    matA_dirichlet_->fillComplete();
//    matB_dirichlet_->scale(-1.0);
//    matB_dirichlet_->fillComplete();

    // Create matrix transposes.
    Tpetra::RowMatrixTransposer<> transposerA(matA_dirichlet_);
//    Tpetra::RowMatrixTransposer<> transposerB(matB_dirichlet_);
    matA_dirichlet_trans_ = transposerA.createTranspose();
//    matB_dirichlet_trans_ = transposerB.createTranspose();


    /*********************************/
    /*** Construct solver objects. ***/
    /*********************************/

    // Construct solver using Amesos2 factory.
    try{
      solverA_ = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("KLU2", matA_dirichlet_);
    } catch (std::invalid_argument& e) {
      std::cout << e.what() << std::endl;
    }
    try{
      solverA_trans_ = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("KLU2", matA_dirichlet_trans_);
    } catch (std::invalid_argument& e) {
      std::cout << e.what() << std::endl;
    }
    solverA_->numericFactorization();
    solverA_trans_->numericFactorization();


    /****************************************/
    /****************************************/

    //outputTpetraData();

  }


  void updateF(const std::vector<Real> &param) {
    Intrepid::FieldContainer<GO> &cellDofs = *(dofMgr_->getCellDofs());
    int numLocalDofs = cellDofs.dimension(1);
    for (GO i=0; i<numCells_; ++i) {                                                        // evaluate RHS function at cubature points
      for (int j=0; j<numCubPoints_; ++j) {
        (*dataF_)(i, j) = funcRHS((*cubPointsPhysical_)(i, j, 0),
                                  (*cubPointsPhysical_)(i, j, 1),
                                  param);
      }
    }
    Intrepid::FunctionSpaceTools::integrate<Real>(*datavalVecF_,                             // compute local data.val vectors for RHS F
                                                  *dataF_,
                                                  *valPhysicalWeighted_,
                                                  Intrepid::COMP_CPP);

    // vecF_ requires assembly using vecF_overlap_ and redistribution
    vecF_ = ROL::makePtr<Tpetra::MultiVector<>>(matA_->getRangeMap(), 1, true);
    vecF_overlap_ = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapMap_, 1, true);
    for (GO i=0; i<numCells_; ++i) {                                                 // assembly on the overlap map
      for (int j=0; j<numLocalDofs; ++j) {
        vecF_overlap_->sumIntoGlobalValue(cellDofs(myCellIds_[i],j),
                                          0,
                                          (*datavalVecF_)[i*numLocalDofs+j]);
      }
    }
    Tpetra::Export<> exporter(vecF_overlap_->getMap(), vecF_->getMap());              // redistribution:
    vecF_->doExport(*vecF_overlap_, exporter, Tpetra::ADD);                           // from the overlap map to the unique map

    Tpetra::deep_copy(*vecF_dirichlet_, *vecF_);
    for (int i=0; i < static_cast<int>(myDirichletDofs_.size()); ++i) {
      if (myUniqueMap_->isNodeGlobalElement(myDirichletDofs_[i])) {
        vecF_dirichlet_->replaceGlobalValue(myDirichletDofs_[i], 0, 0);
      }
    }
  }


  void updateA(const std::vector<Real> &param) {
    const Real one(1);

    Intrepid::FieldContainer<GO> &cellDofs = *(dofMgr_->getCellDofs());
    Teuchos::ArrayRCP<const GO> cellDofsArrayRCP = cellDofs.getData();
    int numLocalDofs = cellDofs.dimension(1);

    for (GO i=0; i<numCells_; ++i) {                                                        // evaluate conductivity kappa at cubature points
      for (int j=0; j<numCubPoints_; ++j) {
        (*kappa_)(i, j) = funcKappa((*cubPointsPhysical_)(i, j, 0),
                                    (*cubPointsPhysical_)(i, j, 1),
                                    param);
      }
    }
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(*kappaGradPhysical_,         // multiply with conductivity kappa
                                                                *kappa_,
                                                                *gradPhysical_);
    Intrepid::FunctionSpaceTools::integrate<Real>(*gradgradMats_,                            // compute local grad.(kappa)grad (stiffness) matrices
                                                  *kappaGradPhysical_,
                                                  *gradPhysicalWeighted_,
                                                  Intrepid::COMP_CPP);
    for (GO i=0; i<numCells_; ++i) {                                                        // evaluate advection adv at cubature points
      for (int j=0; j<numCubPoints_; ++j) {
        for (int k=0; k<spaceDim_; ++k) {
          std::vector<Real> adv(spaceDim_, 0);
          funcAdv(adv, (*cubPointsPhysical_)(i, j, 0),
                       (*cubPointsPhysical_)(i, j, 1),
                       param);
          (*adv_)(i, j, k) = adv[k];
        }
      }
    }
    Intrepid::FunctionSpaceTools::dotMultiplyDataField<Real>(*advGradPhysical_,              // multiply with gradient in physical space
                                                             *adv_,
                                                             *gradPhysical_);
    // IMPORTANT: This is a nonsymmetric form; watch the field order in the integrand!!!
    Intrepid::FunctionSpaceTools::integrate<Real>(*advgradvalMats_,                          // compute local (adv.grad).val (advection) matrices
                                                  *valPhysicalWeighted_,
                                                  *advGradPhysical_,
                                                  Intrepid::COMP_CPP);

    Intrepid::RealSpaceTools<Real>::add(*pdeMats_, *gradgradMats_, *advgradvalMats_);        // add diffusion and advection matrices

    int numLocalMatEntries = numLocalDofs*numLocalDofs;
    Teuchos::ArrayRCP<const Real> pdeArrayRCP = pdeMats_->getData();
    matA_->resumeFill();
    matA_->setAllToScalar(0);
    for (GO i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matA_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                   cellDofsArrayRCP(myCellIds_[i]*numLocalDofs, numLocalDofs),
                                   pdeArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
      }
    }
    matA_->fillComplete();

    matA_dirichlet_ = ROL::makePtr<Tpetra::CrsMatrix<>>(*matA_, Teuchos::DataAccess::Copy);
    matA_dirichlet_->resumeFill();
    for (int i=0; i < static_cast<int>(myDirichletDofs_.size()); ++i) {
      if (myUniqueMap_->isNodeGlobalElement(myDirichletDofs_[i])) {
        size_t numRowEntries = matA_dirichlet_->getNumEntriesInGlobalRow(myDirichletDofs_[i]);
        Teuchos::Array<GO> indices(numRowEntries, 0);
        Teuchos::Array<Real> values(numRowEntries, 0);
        Teuchos::Array<Real> canonicalValues(numRowEntries, 0);
        matA_dirichlet_->getGlobalRowCopy(myDirichletDofs_[i], indices, values, numRowEntries);
        for (int j=0; j<indices.size(); ++j) {
          if (myDirichletDofs_[i] == indices[j]) {
            canonicalValues[j] = one;
          }
	}
	matA_dirichlet_->replaceGlobalValues(myDirichletDofs_[i], indices, canonicalValues);
      }
    }
    matA_dirichlet_->fillComplete();

    Tpetra::RowMatrixTransposer<> transposerA(matA_dirichlet_);
    matA_dirichlet_trans_ = transposerA.createTranspose();
    
    // Construct solver using Amesos2 factory.
    try{
      solverA_ = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("KLU2", matA_dirichlet_);
    } catch (std::invalid_argument& e) {
      std::cout << e.what() << std::endl;
    }
    try{
      solverA_trans_ = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("KLU2", matA_dirichlet_trans_);
    } catch (std::invalid_argument& e) {
      std::cout << e.what() << std::endl;
    }
    solverA_->numericFactorization();
    solverA_trans_->numericFactorization();
  }


  ROL::Ptr<Tpetra::CrsMatrix<> > getMatA(const bool &transpose = false) const {
    if (transpose) {
      return matA_dirichlet_trans_;
    }
    else {
      return matA_dirichlet_;
    }
  }

  void applyMatB(ROL::Ptr<Tpetra::MultiVector<> > &out, const std::vector<Real> &in,
                 const bool sumInto = false) const {
    const size_t nz = matB_->getNumVectors();
    if (!sumInto) {
      out->scale(0.0);
    }
    for (size_t i = 0; i < nz; ++i) {
      Teuchos::ArrayView<const size_t> col(&i,1);
      out->update(in[i],*(matB_->subView(col)),1.0);
    }
  }

  void applyMatBtranspose(ROL::Ptr<std::vector<Real> > &out, const Tpetra::MultiVector<> &in) {
    const size_t nz = matB_->getNumVectors();
    out->assign(nz,0.0);
    std::vector<Real> val(1);
    for (size_t i = 0; i < nz; ++i) {
      val[0] = 0.0;
      Teuchos::ArrayView<const size_t> col(&i,1);
      in.dot(*(matB_->subView(col)),val);
      (*out)[i] = val[0];
    }
  }

//  ROL::Ptr<Tpetra::CrsMatrix<> > getMatB(const bool &transpose = false) const {
//    if (transpose) {
//      return matB_dirichlet_trans_;
//    }
//    else {
//      return matB_dirichlet_;
//    }
//  }


  ROL::Ptr<Tpetra::CrsMatrix<> > getMatM() const {
    return matM_;
  }


  ROL::Ptr<Tpetra::CrsMatrix<> > getMatR() const {
    return matM_;
  }


  ROL::Ptr<Tpetra::MultiVector<> > getVecUd() const {
    return vecUd_;
  }


  ROL::Ptr<Tpetra::MultiVector<> > getVecF() const {
    return vecF_dirichlet_;
  }


  ROL::Ptr<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > getSolver(const bool &transpose = false) const {
    if (transpose) {
      return solverA_trans_;
    }
    else {
      return solverA_;
    }
  }


  inline Real funcRHS(const Real &x1, const Real &x2,
                      const std::vector<Real> &param) const {
    const int ns = 5;
    const Real half(0.5), one(1), two(2);
    Real source(0), arg1(0), arg2(0), mag(0), x0(0), y0(0), sx(0), sy(0);
    // Upper and lower bounds on source magintudes
    const std::vector<Real> ml = {1.5, 1.2, 1.5, 1.2, 1.1};
    const std::vector<Real> mu = {2.5, 1.8, 1.9, 2.6, 1.5};
    // Upper and lower bounds on source locations
    const std::vector<Real> xl = {0.45, 0.75, 0.40, 0.05, 0.85};
    const std::vector<Real> xu = {0.55, 0.85, 0.60, 0.35, 0.95};
    const std::vector<Real> yl = {0.25, 0.55, 0.50, 0.45, 0.45};
    const std::vector<Real> yu = {0.35, 0.65, 0.70, 0.65, 0.55};
    // Upper and lower bounds on source widths
    const std::vector<Real> sxl = {0.03, 0.02, 0.01, 0.02, 0.015};
    const std::vector<Real> sxu = {0.07, 0.04, 0.05, 0.04, 0.025};
    const std::vector<Real> syl = {0.04, 0.01, 0.02, 0.02, 0.01};
    const std::vector<Real> syu = {0.12, 0.05, 0.04, 0.04, 0.03};
    for (int i=0; i<ns; ++i) {
      mag  = ml[i] + (mu[i]-ml[i])*half*(param[i]+one);
      x0   = xl[i] + (xu[i]-xl[i])*half*(param[i+1*ns]+one);
      y0   = yl[i] + (yu[i]-yl[i])*half*(param[i+3*ns]+one);
      sx   = sxl[i] + (sxu[i]-sxl[i])*half*(param[i+2*ns]+one);
      sy   = syl[i] + (syu[i]-syl[i])*half*(param[i+4*ns]+one);
      arg1 = std::pow((x1-x0)/sx, two);
      arg2 = std::pow((x2-y0)/sy, two);
      source += mag*std::exp(-half*(arg1+arg2));
    }
    return source;
  }


  inline Real funcControl(const Real &x1, const Real &x2, const Real &x, const Real &y) const {
      const Real sx(0.05), sy(0.05), half(0.5);
      return std::exp(- half*(x1-x)*(x1-x) / (sx*sx)
                      - half*(x2-y)*(x2-y) / (sy*sy));
//    const std::vector<Real> x0  = {0.2, 0.2, 0.8, 0.8, 0.5};
//    const std::vector<Real> y0  = {0.2, 0.8, 0.2, 0.8, 0.5};
//    const std::vector<Real> sx0 = {0.05, 0.05, 0.05, 0.05, 0.05};
//    const std::vector<Real> sy0 = {0.05, 0.05, 0.05, 0.05, 0.05};
//    const int ns = x0.size();
//    const Real half(0.5);
//    Real val(0);
//    for (int i=0; i<ns; ++i) {
//      val += std::exp(- half*(x1-x0[i])*(x1-x0[i]) / (sx0[i]*sx0[i])
//                      - half*(x2-y0[i])*(x2-y0[i]) / (sy0[i]*sy0[i]) );
//    }
//    return val;
    //return (x1 <= 0.1 ? static_cast<Real>(1) : static_cast<Real>(0));
  }


  inline Real funcTarget(const Real &x1, const Real &x2) const {
    return static_cast<Real>(0);
  }


  inline Real funcKappa(const Real &x1, const Real &x2,
                        const std::vector<Real> &param) const {
    // random diffusion coefficient from i. babuska, f. nobile, r. tempone 2010.
    // simplified model for random stratified media.
    const int ns = 10;
    const Real one(1), two(2), three(3), eight(8), sixteen(16), half(0.5);
    const Real lc = one/sixteen, sqrtpi = std::sqrt(M_PI);
    const Real xi = std::sqrt(sqrtpi*lc), sqrt3 = std::sqrt(three);
    Real f(0), phi(0);
    Real val = one + sqrt3*param[25]*std::sqrt(sqrtpi*lc*half);
    Real arg = one + sqrt3*std::sqrt(sqrtpi*lc*half);
    for (int i = 1; i < ns; ++i) {
      f = floor(half*static_cast<Real>(i+1));
      phi = ((i+1)%2 ? std::sin(f*M_PI*x1) : std::cos(f*M_PI*x1));
      val += xi*std::exp(-std::pow(f*M_PI*lc,two)/eight)*phi*sqrt3*param[i+25];
      arg += xi*std::exp(-std::pow(f*M_PI*lc,two)/eight)*std::abs(phi)*sqrt3;
    }
    return half + two*std::exp(val)/std::exp(arg);
  }


  inline void funcAdv(std::vector<Real> &adv,
                      const Real &x1, const Real &x2,
                      const std::vector<Real> &param) const {
    const Real half(0.5), one(1), five(5), ten(10);
    const Real a = five*half*(param[36]+one);
    const Real b = five + (ten-five)*half*(param[35]+one);
    adv[0] = b - a*x1;
    adv[1] = a*x2;
//    const Real half(0.5), one(1), four(4), ten(10);
//    const Real mag = one + (ten-one)*half*(param[35]+one);
//    const Real ang = M_PI/four * param[36];
//    adv[0] = mag*std::cos(ang);
//    adv[1] = mag*std::sin(ang);
  }


  void printMeshData(std::ostream &outStream) const {
    ROL::Ptr<Intrepid::FieldContainer<Real> > nodesPtr = meshMgr_->getNodes();
    ROL::Ptr<Intrepid::FieldContainer<GO> >  cellToNodeMapPtr = meshMgr_->getCellToNodeMap();
    Intrepid::FieldContainer<Real>  &nodes = *nodesPtr;
    Intrepid::FieldContainer<GO>   &cellToNodeMap = *cellToNodeMapPtr;
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
    }
  }


  void outputTpetraData() const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> >   matWriter;
    matWriter.writeSparseFile("stiffness_mat", matA_);
    matWriter.writeSparseFile("dirichlet_mat", matA_dirichlet_);
    matWriter.writeSparseFile("mass_mat", matM_);
    matWriter.writeDenseFile("Ud_vec", vecUd_);
  }


  void outputTpetraVector(const ROL::Ptr<const Tpetra::MultiVector<> > &vec,
                          const std::string &filename) const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> > vecWriter;
    vecWriter.writeDenseFile(filename, vec);
  }

}; // class PoissonData

#endif
