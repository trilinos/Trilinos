// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  data.hpp
    \brief Generates and manages data for the Stefan-Boltzmann example, including
           all mesh and discretization data, matrices, etc.
*/

#ifndef ROL_PDEOPT_STEFANBOLTZMANN_DATA_H
#define ROL_PDEOPT_STEFANBOLTZMANN_DATA_H

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
class StefanBoltzmannData {

private:
  using GO = typename Tpetra::Map<>::global_ordinal_type;

  ROL::Ptr<MeshManager<Real> > meshMgr_;
  ROL::Ptr<DofManager<Real> >  dofMgr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;

  ROL::Ptr<const Teuchos::Comm<int> > commPtr_;
  int myRank_;
  int numProcs_;

  Real alpha_;
  int  basisOrder_;

  ROL::Ptr<const Tpetra::Map<> >    myOverlapMap_;
  ROL::Ptr<const Tpetra::Map<> >    myUniqueMap_;
  ROL::Ptr<Tpetra::CrsGraph<> >     matGraph_;
  ROL::Ptr<Tpetra::CrsMatrix<> >    matA_;
  ROL::Ptr<Tpetra::CrsMatrix<> >    matA_dirichlet_;
  ROL::Ptr<Tpetra::CrsMatrix<> >    matA_dirichlet_trans_;
  ROL::Ptr<Tpetra::CrsMatrix<> >    matM_;
  ROL::Ptr<Tpetra::CrsMatrix<> >    matM_dirichlet_;
  ROL::Ptr<Tpetra::CrsMatrix<> >    matM_dirichlet_trans_;
  ROL::Ptr<Tpetra::MultiVector<> >  vecUd_;
  ROL::Ptr<Tpetra::MultiVector<> >  vecF_;
  ROL::Ptr<Tpetra::MultiVector<> >  vecF_overlap_;
  ROL::Ptr<Tpetra::MultiVector<> >  vecF_dirichlet_;

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
  ROL::Ptr<Intrepid::FieldContainer<Real> > valPhysicalWeighted_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > gradPhysicalWeighted_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > gradgradMats_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > valvalMats_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > cubPointsPhysical_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > kappa_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > dataF_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > datavalVecF_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > dofPoints_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > dofPointsPhysical_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > dataUd_;

  int sdim_;

public:

  StefanBoltzmannData(const ROL::Ptr<const Teuchos::Comm<int> > &comm,
                      const Teuchos::RCP<Teuchos::ParameterList> &parlist,
                      const ROL::Ptr<std::ostream> &outStream) {
    sdim_ = parlist->sublist("Problem").get("Stochastic Dimension",6);
    std::vector<Real> par(sdim_,1.0);

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
    alpha_ = parlist->sublist("Problem").get("Penalty", 1e-2);
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
    valPhysicalWeighted_  = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, numCubPoints_);
    gradPhysicalWeighted_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, numCubPoints_, spaceDim_);
    gradgradMats_         = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, lfs);
    valvalMats_           = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs, lfs);
    kappa_                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_);
    dataF_                = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_);
    datavalVecF_          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs);
    dataUd_               = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs);

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
    for (GO i=0; i<numCells_; ++i) {                                                        // evaluate conductivity kappa at cubature points
      for (int j=0; j<numCubPoints_; ++j) {
        (*kappa_)(i, j) = funcKappa((*cubPointsPhysical_)(i, j, 0),
                                    (*cubPointsPhysical_)(i, j, 1),
                                    par);
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

    for (GO i=0; i<numCells_; ++i) {                                                        // evaluate RHS function at cubature points
      for (int j=0; j<numCubPoints_; ++j) {
        (*dataF_)(i, j) = funcRHS((*cubPointsPhysical_)(i, j, 0),
                                  (*cubPointsPhysical_)(i, j, 1),
                                  par);
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
    for (GO i=0; i<numCells_; ++i) {                                                        // evaluate target function at these points
      for (int j=0; j<lfs; ++j) {
        (*dataUd_)(i, j) = funcTarget((*dofPointsPhysical_)(i, j, 0), (*dofPointsPhysical_)(i, j, 1));
      }
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
    // Stiffness matrix A.
    matA_ = ROL::makePtr<Tpetra::CrsMatrix<>>(matGraph_);
    int numLocalMatEntries = numLocalDofs*numLocalDofs;
    Teuchos::ArrayRCP<const Real> gradgradArrayRCP = gradgradMats_->getData();
    for (GO i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matA_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                   cellDofsArrayRCP(myCellIds_[i]*numLocalDofs, numLocalDofs),
                                   gradgradArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
      }
    }
    matA_->fillComplete();

    // Mass matrix M.
    matM_ = ROL::makePtr<Tpetra::CrsMatrix<>>(matGraph_);
    Teuchos::ArrayRCP<const Real> valvalArrayRCP = valvalMats_->getData();
    for (GO i=0; i<numCells_; ++i) {
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
    matM_dirichlet_ = ROL::makePtr<Tpetra::CrsMatrix<>>(*matA_, Teuchos::DataAccess::Copy);
    vecF_dirichlet_ = ROL::makePtr<Tpetra::MultiVector<>>(matA_->getRangeMap(), 1, true);
    Tpetra::deep_copy(*vecF_dirichlet_, *vecF_);
    ROL::Ptr<std::vector<std::vector<Intrepid::FieldContainer<GO> > > > dirichletSideSets = meshMgr_->getSideSets();
    std::vector<std::vector<Intrepid::FieldContainer<GO> > > &dss = *dirichletSideSets;
    Teuchos::Array<GO> mySortedCellIds_(myCellIds_);
    std::sort(mySortedCellIds_.begin(), mySortedCellIds_.end());
    mySortedCellIds_.erase( std::unique(mySortedCellIds_.begin(), mySortedCellIds_.end()), mySortedCellIds_.end() );
    std::vector<Teuchos::Array<GO> > myDirichletCellIds_(dss[0].size());
    for (int i=0; i<static_cast<int>(dss[0].size()); ++i) {
      for (int j=0; j<dss[0][i].dimension(0); ++j) {
        if (std::binary_search(mySortedCellIds_.begin(), mySortedCellIds_.end(), dss[0][i](j))) {
          myDirichletCellIds_[i].push_back(dss[0][i](j));
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
    matM_dirichlet_->resumeFill();
    for (int i=0; i<myDirichletDofs_.size(); ++i) {
      if (myUniqueMap_->isNodeGlobalElement(myDirichletDofs_[i])) {
        size_t numRowEntries = matA_dirichlet_->getNumEntriesInGlobalRow(myDirichletDofs_[i]);
        Teuchos::Array<GO> indices(numRowEntries, 0);    
        Teuchos::Array<Real> values(numRowEntries, 0);
        Teuchos::Array<Real> canonicalValues(numRowEntries, 0);    
        Teuchos::Array<Real> zeroValues(numRowEntries, 0);    
        matA_dirichlet_->getGlobalRowCopy(myDirichletDofs_[i], indices, values, numRowEntries);
        matM_dirichlet_->getGlobalRowCopy(myDirichletDofs_[i], indices, values, numRowEntries);
        for (int j=0; j<indices.size(); ++j) {
          if (myDirichletDofs_[i] == indices[j]) {
            canonicalValues[j] = 1.0;
          }
        }
        matA_dirichlet_->replaceGlobalValues(myDirichletDofs_[i], indices, canonicalValues);
        matM_dirichlet_->replaceGlobalValues(myDirichletDofs_[i], indices, zeroValues);
        vecF_dirichlet_->replaceGlobalValue(myDirichletDofs_[i], 0, 0);
      }
    }
    matA_dirichlet_->fillComplete();
    matM_dirichlet_->fillComplete();

    // Create matrix transposes.
    Tpetra::RowMatrixTransposer<> transposerA(matA_dirichlet_);
    Tpetra::RowMatrixTransposer<> transposerM(matM_dirichlet_);
    matA_dirichlet_trans_ = transposerA.createTranspose();
    matM_dirichlet_trans_ = transposerM.createTranspose();

    updateA(par);

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


  ROL::Ptr<Tpetra::CrsMatrix<> > getMatA(const bool &transpose = false) const {
    if (transpose) {
      return matA_dirichlet_trans_;
    }
    else {
      return matA_dirichlet_;
    }
  }


  ROL::Ptr<Tpetra::CrsMatrix<> > getMatB(const bool &transpose = false) const {
    if (transpose) {
      return matM_dirichlet_trans_;
    }
    else {
      return matM_dirichlet_;
    }
  }


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


  Real funcRHS(const Real &x1, const Real &x2, const std::vector<Real> &par) const {
    Real xi1 = 0.5*(par[sdim_-2]+1.0), xi2 = 0.5*(par[sdim_-1]+1.0);
    return std::exp(-1.e-2*(std::pow(x1-xi1,2.0)+std::pow(x2-xi2,2.0)));
    // return 2.0*M_PI*M_PI*std::sin(M_PI*x1)*std::sin(M_PI*x2) + (1.0/(alpha_*128.0*M_PI*M_PI))*std::sin(8.0*M_PI*x1)*std::sin(8.0*M_PI*x2);
  }


  Real funcTarget(const Real &x1, const Real &x2) const {
    return 1.0;
    // return std::sin(M_PI*x1)*std::sin(M_PI*x2) - std::sin(8.0*M_PI*x1)*std::sin(8.0*M_PI*x2);
  }


  Real funcStateSolution(const Real &x1, const Real &x2) const {
    return std::sin(M_PI*x1)*std::sin(M_PI*x2);
  }


  Real funcKappa(const Real &x1, const Real &x2, const std::vector<Real> &par) const {
    // Random diffusion coefficient from I. Babuska, F. Nobile, R. Tempone 2010.
    // Simplified model for random stratified media.
    Real Lc = 1.0/16.0, sqrtpi = std::sqrt(M_PI), xi = std::sqrt(sqrtpi*Lc);
    Real sqrt3 = std::sqrt(3.0), f = 0.0, phi = 0.0;
    Real val = 1.0 + sqrt3*par[0]*std::sqrt(sqrtpi*Lc*0.5);
    for (int i = 1; i < sdim_-2; i++) {
      f = floor((Real)(i+1)/2.0);
      phi = ((i+1)%2 ? std::sin(f*M_PI*x1) : std::cos(f*M_PI*x1));
      val += xi*std::exp(-std::pow(f*M_PI*Lc,2.0)/8.0)*phi*sqrt3*par[i];
    }
    return 0.5 + std::exp(val);
  }

  void updateF(const std::vector<Real> &par) {
    Intrepid::FieldContainer<GO> &cellDofs = *(dofMgr_->getCellDofs());
    int numLocalDofs = cellDofs.dimension(1);
    for (GO i=0; i<numCells_; ++i) {                                                        // evaluate RHS function at cubature points
      for (int j=0; j<numCubPoints_; ++j) {
        (*dataF_)(i, j) = funcRHS((*cubPointsPhysical_)(i, j, 0),
                                  (*cubPointsPhysical_)(i, j, 1),
                                  par);
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
    for (int i=0; i<myDirichletDofs_.size(); ++i) {
      if (myUniqueMap_->isNodeGlobalElement(myDirichletDofs_[i])) {
        vecF_dirichlet_->replaceGlobalValue(myDirichletDofs_[i], 0, 0);
      }
    }
  }

  void updateA(const std::vector<Real> &par) {

    Intrepid::FieldContainer<GO> &cellDofs = *(dofMgr_->getCellDofs());
    Teuchos::ArrayRCP<const GO> cellDofsArrayRCP = cellDofs.getData();
    int numLocalDofs = cellDofs.dimension(1);

    for (GO i=0; i<numCells_; ++i) {                                                        // evaluate conductivity kappa at cubature points
      for (int j=0; j<numCubPoints_; ++j) {
        (*kappa_)(i, j) = funcKappa((*cubPointsPhysical_)(i, j, 0),
                                    (*cubPointsPhysical_)(i, j, 1),
                                    par);
      }
    }
    Intrepid::FunctionSpaceTools::tensorMultiplyDataField<Real>(*kappaGradPhysical_,         // multiply with conductivity kappa
                                                                *kappa_,
                                                                *gradPhysical_);
    Intrepid::FunctionSpaceTools::integrate<Real>(*gradgradMats_,                            // compute local grad.(kappa)grad (stiffness) matrices
                                                  *kappaGradPhysical_,
                                                  *gradPhysicalWeighted_,
                                                  Intrepid::COMP_CPP);

    int numLocalMatEntries = numLocalDofs*numLocalDofs;
    Teuchos::ArrayRCP<const Real> gradgradArrayRCP = gradgradMats_->getData();
    matA_->resumeFill();
    matA_->setAllToScalar(0);
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs; ++j) {
        matA_->sumIntoGlobalValues(cellDofs(myCellIds_[i],j),
                                   cellDofsArrayRCP(myCellIds_[i]*numLocalDofs, numLocalDofs),
                                   gradgradArrayRCP(i*numLocalMatEntries+j*numLocalDofs, numLocalDofs));
      }
    }
    matA_->fillComplete();

    matA_dirichlet_ = ROL::makePtr<Tpetra::CrsMatrix<>>(*matA_, Teuchos::DataAccess::Copy);
    matA_dirichlet_->resumeFill();
    for (int i=0; i<myDirichletDofs_.size(); ++i) {
      if (myUniqueMap_->isNodeGlobalElement(myDirichletDofs_[i])) {
        size_t numRowEntries = matA_dirichlet_->getNumEntriesInGlobalRow(myDirichletDofs_[i]);
        Teuchos::Array<GO> indices(numRowEntries, 0);
        Teuchos::Array<Real> values(numRowEntries, 0);
        Teuchos::Array<Real> canonicalValues(numRowEntries, 0);
        matA_dirichlet_->getGlobalRowCopy(myDirichletDofs_[i], indices, values, numRowEntries);
        for (int j=0; j<indices.size(); ++j) {
          if (myDirichletDofs_[i] == indices[j]) {
            canonicalValues[j] = 1.0;
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

  Real computeStateError(const ROL::Ptr<const Tpetra::MultiVector<> > &soln) const {

    ROL::Ptr<Tpetra::MultiVector<> > soln_overlap =
      ROL::makePtr<Tpetra::MultiVector<>>(vecF_overlap_->getMap(), 1, true);
    Tpetra::Import<> importer(vecUd_->getMap(), soln_overlap->getMap());              // redistribution:
    soln_overlap->doImport(*soln, importer, Tpetra::REPLACE);                         // from the unique map to the overlap map

    Intrepid::DefaultCubatureFactory<Real> cubFactory;                                       // create cubature factory
    int cubDeg = 6;                                                                          // set cubature degree, e.g., 6
    ROL::Ptr<Intrepid::Cubature<Real> > cellCub = cubFactory.create(cellType_, cubDeg);  // create cubature for error computation
    int numCubPts = cellCub->getNumPoints();                                                 // retrieve number of cubature points
    int lfs = dofMgr_->getLocalFieldSize(0);
    Intrepid::FieldContainer<Real> cubPts(numCubPts, spaceDim_);
    Intrepid::FieldContainer<Real> cubWts(numCubPts);
    Intrepid::FieldContainer<Real> cubPtsPhys(numCells_, numCubPts, spaceDim_);
    Intrepid::FieldContainer<Real> jac(numCells_, numCubPts, spaceDim_, spaceDim_);
    Intrepid::FieldContainer<Real> jacDet(numCells_, numCubPts);
    Intrepid::FieldContainer<Real> valRef(lfs, numCubPts);
    Intrepid::FieldContainer<Real> valPhys(numCells_, lfs, numCubPts);
    Intrepid::FieldContainer<Real> wtMeas(numCells_, numCubPts);
    Intrepid::FieldContainer<Real> inCoeffs(numCells_, lfs);
    Intrepid::FieldContainer<Real> funcVals(numCells_, numCubPts);
    Intrepid::FieldContainer<Real> funcValsWt(numCells_, numCubPts);
    Intrepid::FieldContainer<Real> normSquaredError(numCells_);

    cellCub->getCubature(cubPts, cubWts);                                         // retrieve cubature points and weights
    (*basisPtrs_[0]).getValues(valRef, cubPts, Intrepid::OPERATOR_VALUE);         // evaluate value operator at cubature points

    Intrepid::CellTools<Real>::setJacobian(jac, cubPts, *cellNodes_, cellType_);  // compute cell Jacobians
    Intrepid::CellTools<Real>::setJacobianDet(jacDet, jac);                       // compute determinants of cell Jacobians

    Intrepid::FunctionSpaceTools::HGRADtransformVALUE<Real>(valPhys,              // transform reference values into physical space
                                                            valRef);

    Intrepid::FunctionSpaceTools::computeCellMeasure<Real>(wtMeas,                // compute weighted cell measure
                                                           jacDet,
                                                           cubWts);

    Intrepid::CellTools<Real>::mapToPhysicalFrame(cubPtsPhys,                     // map reference cubature points to physical space
                                                  cubPts,
                                                  *cellNodes_,
                                                  cellType_);

    Intrepid::FieldContainer<GO> &cellDofs = *(dofMgr_->getCellDofs());
    Teuchos::ArrayRCP<const Real> soln_data = soln_overlap->get1dView();                // populate inCoeffs
    for (GO i=0; i<numCells_; ++i) {
      for (int j=0; j<lfs; ++j) {
        inCoeffs(i, j) = soln_data[soln_overlap->getMap()->getLocalElement(cellDofs(myCellIds_[i],j))];
      }
    }

    Intrepid::FunctionSpaceTools::evaluate<Real>(funcVals, inCoeffs, valPhys);

    for (GO i=0; i<numCells_; ++i) {                                                   // compute error
      for (int j=0; j<numCubPts; ++j) {
        funcVals(i, j) -= funcStateSolution(cubPtsPhys(i, j, 0), cubPtsPhys(i, j, 1));
      }
    }

    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(funcValsWt,              // multiply with weighted measure
                                                               wtMeas,
                                                               funcVals);

    Intrepid::FunctionSpaceTools::integrate<Real>(normSquaredError,                     // compute norm squared of local error
                                                  funcVals,
                                                  funcValsWt,
                                                  Intrepid::COMP_CPP);

    Real localErrorSum(0);
    Real globalErrorSum(0);
    for (GO i=0; i<numCells_; ++i) {
      localErrorSum += normSquaredError(i);
    }
    ROL::Ptr<const Teuchos::Comm<int> > comm = soln_overlap->getMap()->getComm();
    Teuchos::reduceAll<int, Real>(*comm, Teuchos::REDUCE_SUM, 1, &localErrorSum, &globalErrorSum);

    return globalErrorSum;
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

}; // class StefanBoltzmannData

#endif
