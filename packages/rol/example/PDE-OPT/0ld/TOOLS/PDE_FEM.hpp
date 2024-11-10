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

#ifndef ROL_PDEOPT_PDE_FEM_H
#define ROL_PDEOPT_PDE_FEM_H

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
#include "./dofmanager.hpp"

// Global Timers.
#ifdef ROL_TIMERS
ROL::Ptr<Teuchos::Time> FactorTime_example_PDEOPT_TOOLS_PDEFEM_GLOB =
  Teuchos::TimeMonitor::getNewCounter("ROL: Factorization Time in PDEFEM");
ROL::Ptr<Teuchos::Time> LUSubstitutionTime_example_PDEOPT_TOOLS_PDEFEM_GLOB =
  Teuchos::TimeMonitor::getNewCounter("ROL: LU Substitution Time in PDEFEM");
ROL::Ptr<Teuchos::Time> SolverUpdateTime_example_PDEOPT_TOOLS_PDEFEM_GLOB =
  Teuchos::TimeMonitor::getNewCounter("ROL: Solver Update Time in PDEFEM");
ROL::Ptr<Teuchos::Time> LocalAssemblyTime_example_PDEOPT_TOOLS_PDEFEM_GLOB =
  Teuchos::TimeMonitor::getNewCounter("ROL: Local Assembly Time in PDEFEM");
ROL::Ptr<Teuchos::Time> ConstraintDerivativeTime_example_PDEOPT_TOOLS_PDEFEM_GLOB =
  Teuchos::TimeMonitor::getNewCounter("ROL: Constraint Derivative Application Time in PDEFEM");
#endif

template<class Real>
class PDE_FEM {
private:
  using GO = typename Tpetra::Map<>::global_ordinal_type;

protected:

  ROL::Ptr<MeshManager<Real> > meshMgr_;
  ROL::Ptr<DofManager<Real> >  dofMgr_;
  std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  ROL::Ptr<Intrepid::DofCoordsInterface<Intrepid::FieldContainer<Real> > > coord_iface_;
  Intrepid::FieldContainer<GO> cellDofs_;

  Teuchos::RCP<Teuchos::ParameterList> parlist_;
  ROL::Ptr<const Teuchos::Comm<int> > commPtr_;
  ROL::Ptr<std::ostream> outStream_;
    
  int numLocalDofs_;
  Intrepid::FieldContainer<GO> ctn_;
  Intrepid::FieldContainer<GO> cte_;

  int myRank_;
  int numProcs_;

  Real alpha_;
  int  basisOrder_;

  ROL::Ptr<const Tpetra::Map<> >    myCellMap_;
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
  
  Teuchos::Array<Real> myCellMeasure_;
  Teuchos::Array<GO> myCellIds_;
// Elements on Boundary  
  std::vector<Teuchos::Array<GO> > myBoundaryCellIds_;
// DBC
  Teuchos::Array<GO> myDirichletDofs_;
// BC Sides 
  std::vector<GO> my_dbc_;
  std::vector<GO> my_nbc_;

//Point load on Bundary!
  std::vector<GO> myPointLoadDofs_;
  std::vector<Real> myPointLoadVals_;

  ROL::Ptr<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > solverA_;
  ROL::Ptr<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > solverA_trans_;

  shards::CellTopology cellType_;
  
  int sideDim_;
  int spaceDim_;
  int numNodesPerCell_;
  int numCubPoints_;
  int lfs_;

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
  ROL::Ptr<Intrepid::FieldContainer<Real> > valPhysicalWeighted_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > gradPhysicalWeighted_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > gradgradMats_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > valvalMats_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > cubPointsPhysical_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > dataF_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > datavalVecF_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > dofPoints_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > dofPointsPhysical_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > dataUd_;

  bool verbose_;

protected:
public:

  // constructor
  //PDE_FEM() {}
  PDE_FEM() {}

  // destructor
  virtual ~PDE_FEM() {}
  
  virtual void Initialize(const ROL::Ptr<const Teuchos::Comm<int> > &comm,
                          const Teuchos::RCP<Teuchos::ParameterList> &parlist,
                          const ROL::Ptr<std::ostream> &outStream) {
    commPtr_   = comm;
    parlist_   = parlist;
    outStream_ = outStream;
    myRank_    = comm->getRank();
    numProcs_  = comm->getSize();

    verbose_ = parlist->sublist("PDE FEM").get("Verbose Output",false);
    if ( verbose_ ) {
      *outStream_ << "Total number of processors: " << numProcs_ << std::endl;
    }

    basisOrder_ = parlist->sublist("PDE FEM").get<int>("Order of FE Discretization");
  }

  void SetDiscretization(const ROL::Ptr<MeshManager<Real> > &meshMgr,
                         const std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > &basisPtrs) {
    meshMgr_ = meshMgr;
    totalNumCells_ = meshMgr_->getNumCells();

    basisPtrs_ = basisPtrs;

    // Retrieve some basic cell information.
    cellType_ = (basisPtrs_[0])->getBaseCellTopology(); // get the cell type from any basis
    spaceDim_ = cellType_.getDimension();        // retrieve spatial dimension
    numNodesPerCell_ = cellType_.getNodeCount(); // retrieve number of nodes per cell
    
    sideDim_ = spaceDim_ - 1;

    coord_iface_ = ROL::dynamicPtrCast<Intrepid::DofCoordsInterface<Intrepid::FieldContainer<Real> > >(basisPtrs_[0]);
    dofMgr_ = ROL::makePtr<DofManager<Real>>(meshMgr_, basisPtrs_);
  }

  virtual void SetParallelStructure(void) {
    int cellSplit = parlist_->sublist("Geometry").get<int>("Partition type");
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
      
    cellDofs_ = *(dofMgr_->getCellDofs());
    numLocalDofs_ = cellDofs_.dimension(1);
    if ( verbose_ ) {
      *outStream_ << "Cell offsets across processors: " << cellOffsets_
                  << std::endl;
    }
    for (int i=0; i<numCells_; ++i) {
      myCellIds_.push_back(cellOffsets_[myRank_]+i);
      for (int j=0; j<numLocalDofs_; ++j) {
        myGlobIds_.push_back( cellDofs_(cellOffsets_[myRank_]+i,j) );
      }
    }
    std::sort(myGlobIds_.begin(), myGlobIds_.end());
    myGlobIds_.erase( std::unique(myGlobIds_.begin(), myGlobIds_.end()), myGlobIds_.end() );

    // Build maps.
    myOverlapMap_ = ROL::makePtr<Tpetra::Map<>>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),                                            myGlobIds_, 0, commPtr_);
    //std::cout << std::endl << myOverlapMap_->getLocalElementList()<<std::endl;
    /** One can also use the non-member function:
        myOverlapMap_ = Tpetra::createNonContigMap<int,int>(myGlobIds_, commPtr_);
        to build the overlap map.
    **/
    myUniqueMap_ = Tpetra::createOneToOne(myOverlapMap_);
    //std::cout << std::endl << myUniqueMap_->getLocalElementList() << std::endl;
    myCellMap_ = ROL::makePtr<Tpetra::Map<>>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                              this->myCellIds_, 0, this->commPtr_);
  }
  
  virtual void SetUpLocalIntrepidArrays(void) {
    /****************************************************/
    /*** Set up local discretization data and arrays. ***/
    /****************************************************/
    // Cubature data.
    Intrepid::DefaultCubatureFactory<Real> cubFactory;                                          // create cubature factory
    int cubDegree = 4;                                                                          // set cubature degree, e.g., 2
    ROL::Ptr<Intrepid::Cubature<Real> > cellCub = cubFactory.create(cellType_, cubDegree);  // create default cubature
    numCubPoints_ = cellCub->getNumPoints();                                                    // retrieve number of cubature points

    lfs_ = dofMgr_->getLocalFieldSize(0);
    // Discretization data. 
    cubPoints_            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCubPoints_, spaceDim_);
    cubWeights_           = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCubPoints_);
    cubPointsPhysical_    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_, spaceDim_);
    dofPoints_            = ROL::makePtr<Intrepid::FieldContainer<Real>>(lfs_, spaceDim_);
    dofPointsPhysical_    = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs_, spaceDim_);
    cellNodes_            = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numNodesPerCell_, spaceDim_);
    cellJac_              = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_, spaceDim_, spaceDim_);
    cellJacInv_           = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_, spaceDim_, spaceDim_);
    cellJacDet_           = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_);
    cellWeightedMeasure_  = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, numCubPoints_);
    valReference_         = ROL::makePtr<Intrepid::FieldContainer<Real>>(lfs_, numCubPoints_);  
    gradReference_        = ROL::makePtr<Intrepid::FieldContainer<Real>>(lfs_, numCubPoints_, spaceDim_);  
    valPhysical_          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs_, numCubPoints_);
    gradPhysical_         = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs_, numCubPoints_, spaceDim_);
    valPhysicalWeighted_  = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs_, numCubPoints_);
    gradPhysicalWeighted_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells_, lfs_, numCubPoints_, spaceDim_);
    
    // Geometric definition of the cells in the mesh, based on the cell-to-node map and the domain partition.
    Intrepid::FieldContainer<Real> &nodes = *meshMgr_->getNodes();
    ctn_   = *meshMgr_->getCellToNodeMap();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numNodesPerCell_; ++j) {
        for (int k=0; k<spaceDim_; ++k) {
          (*cellNodes_)(i, j, k) = nodes(ctn_(myCellIds_[i],j), k);
        }
      }
    }

    //Compute Intrepid Arrays
    cellCub->getCubature(*cubPoints_, *cubWeights_);                                         // retrieve cubature points and weights
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

    Intrepid::FunctionSpaceTools::HGRADtransformVALUE<Real>(*valPhysical_,                   // transform reference values into physical space
                                                            *valReference_);
    Intrepid::FunctionSpaceTools::multiplyMeasure<Real>(*valPhysicalWeighted_,               // multiply with weighted measure
                                                        *cellWeightedMeasure_,
                                                        *valPhysical_);
    Intrepid::CellTools<Real>::mapToPhysicalFrame(*cubPointsPhysical_,                        // map reference cubature points to physical space
                                                  *cubPoints_,
                                                  *cellNodes_,
                                                  cellType_);
    
    ComputeCellTotalMeasures();
  } 
 
  virtual void ComputeCellTotalMeasures(void)
  {
	for (int i=0; i<numCells_; ++i){
		Real temp = 0.0;
		for(int j=0; j<numCubPoints_; ++j){	
    			temp += (*cellWeightedMeasure_)(i, j);
		}
		myCellMeasure_.push_back(temp);
	}
	std::cout<<"First cell total measure: "<<myCellMeasure_[0]<<std::endl;
  }


  virtual void ComputeLocalSystemMats(void) { }
   
  virtual void ComputeLocalForceVec(void) { }
  
  virtual void SetUpTrueDataOnNodes(void) { }
  
  virtual void AssembleSystemGraph(void) { 
    /****************************************/
    /*** Assemble global graph structure. ***/
    /****************************************/

    Teuchos::ArrayRCP<const GO> cellDofsArrayRCP = cellDofs_.getData();
    Teuchos::Array<size_t> graphEntriesPerRow(myUniqueMap_->getLocalNumElements());
    for (GO i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs_; ++j) {
        GO gid = cellDofs_(myCellIds_[i],j);
        if(myUniqueMap_->isNodeGlobalElement(gid))
          graphEntriesPerRow[myUniqueMap_->getLocalElement(gid)] += numLocalDofs_;
      }
    }
    matGraph_ = ROL::makePtr<Tpetra::CrsGraph<>>(myUniqueMap_, graphEntriesPerRow());
    for (GO i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs_; ++j) {
        matGraph_->insertGlobalIndices(cellDofs_(myCellIds_[i],j), cellDofsArrayRCP(myCellIds_[i]*numLocalDofs_, numLocalDofs_));
      }
    }
    matGraph_->fillComplete();
  }

  virtual void AssembleSystemMats(void) { 
    /****************************************/
    /*** Assemble global data structures. ***/
    /****************************************/
    // Assemble matrices.
    // Stiffness matrix A.
    matA_ = ROL::makePtr<Tpetra::CrsMatrix<>>(matGraph_);
    matA_dirichlet_ = ROL::makePtr<Tpetra::CrsMatrix<>>(matGraph_);
    int numLocalMatEntries = numLocalDofs_ * numLocalDofs_;
    Teuchos::ArrayRCP<const GO> cellDofsArrayRCP = cellDofs_.getData();
    Teuchos::ArrayRCP<const Real> gradgradArrayRCP = gradgradMats_->getData();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs_; ++j) {
        matA_->sumIntoGlobalValues(cellDofs_(myCellIds_[i],j),
                                   cellDofsArrayRCP(myCellIds_[i] * numLocalDofs_, numLocalDofs_),
                                   gradgradArrayRCP(i*numLocalMatEntries+j*numLocalDofs_, numLocalDofs_));
        matA_dirichlet_->sumIntoGlobalValues(cellDofs_(myCellIds_[i],j),
                                             cellDofsArrayRCP(myCellIds_[i] * numLocalDofs_, numLocalDofs_),
                                             gradgradArrayRCP(i*numLocalMatEntries+j*numLocalDofs_, numLocalDofs_));
      }
    }
    matA_->fillComplete();
    matA_dirichlet_->fillComplete();
    // Mass matrix M.
/*
    matM_ = ROL::makePtr<Tpetra::CrsMatrix<>>(matGraph_);
    matM_dirichlet_ = ROL::makePtr<Tpetra::CrsMatrix<>>(matGraph_);
    Teuchos::ArrayRCP<const Real> valvalArrayRCP = valvalMats_->getData();
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs_; ++j) {
        matM_->sumIntoGlobalValues(cellDofs_(myCellIds_[i],j),
                                   cellDofsArrayRCP(myCellIds_[i]*numLocalDofs_, numLocalDofs_),
                                   valvalArrayRCP(i*numLocalMatEntries+j*numLocalDofs_, numLocalDofs_));
        matM_dirichlet_->sumIntoGlobalValues(cellDofs_(myCellIds_[i],j),
                                             cellDofsArrayRCP(myCellIds_[i]*numLocalDofs_, numLocalDofs_),
                                             valvalArrayRCP(i*numLocalMatEntries+j*numLocalDofs_, numLocalDofs_));
      }
    }
    matM_->fillComplete();
    matM_dirichlet_->fillComplete();
*/
  }

  virtual void AssembleRHSVector(void) {
    // Assemble vectors.
    // vecF_ requires assembly using vecF_overlap_ and redistribution
    vecF_           = ROL::makePtr<Tpetra::MultiVector<>>(matA_->getRangeMap(), 1, true);
    vecF_dirichlet_ = ROL::makePtr<Tpetra::MultiVector<>>(matA_->getRangeMap(), 1, true);
    vecF_overlap_   = ROL::makePtr<Tpetra::MultiVector<>>(myOverlapMap_, 1, true);
    // assembly on the overlap map
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs_; ++j) {
        vecF_overlap_->sumIntoGlobalValue(cellDofs_(myCellIds_[i],j),
                                          0,
                                          (*datavalVecF_)[i*numLocalDofs_+j]);
      }
    }
 
    //Assemble the point loads!
    for (unsigned i=0; i<myPointLoadDofs_.size(); ++i) {
      vecF_overlap_->sumIntoGlobalValue(myPointLoadDofs_[i],
                                        0,
                                        myPointLoadVals_[i]);
    }

    // change map
    Tpetra::Export<> exporter(vecF_overlap_->getMap(), vecF_->getMap()); // redistribution
    vecF_->doExport(*vecF_overlap_, exporter, Tpetra::ADD);              // from the overlap map to the unique map
    vecF_dirichlet_->doExport(*vecF_overlap_, exporter, Tpetra::ADD);    // from the overlap map to the unique map
  }
  
  
  virtual void AssembleDataVector(void) {
    // vecUd_ does not require assembly
    vecUd_ = ROL::makePtr<Tpetra::MultiVector<>>(matA_->getDomainMap(), 1, true);
    for (int i=0; i<numCells_; ++i) {
      for (int j=0; j<numLocalDofs_; ++j) {
        if (vecUd_->getMap()->isNodeGlobalElement(cellDofs_(myCellIds_[i],j))) {
          vecUd_->replaceGlobalValue(cellDofs_(myCellIds_[i],j),
                                     0,
                                     (*dataUd_)[i*numLocalDofs_+j]);
        }
      }
    }
  }
 
  // find the local index of a global cell
  virtual int find_local_index(GO globalCell)
  {
	return myCellMap_->getLocalElement(globalCell);
/*
	for(int i=0; i<numCells_; i++)
	{
		if(myCellIds_[i] == globalCell)
			return i;
	}
	return -1;
*/
  }

  // check to see if a globaldof is on dirichlet boundary
  virtual bool check_myGlobalDof_on_boundary(GO globalDof) {
    if (std::binary_search(myDirichletDofs_.begin(), myDirichletDofs_.end(), globalDof)) {
      return true;
    }
    return false;
  }

  //create myBoundaryCellIds_ and myDirichletDofs_
  virtual void SetUpMyDBCInfo(bool ifUseCoordsToSpecifyBC, std::vector<GO> dbc_side) 
{
    if(ifUseCoordsToSpecifyBC){
	my_dbc_.resize(4);
	my_dbc_ = {0, 1, 2, 3};
    }else{
    	my_dbc_ = dbc_side;
    }
    //Print to check my BC info
    if ( verbose_ ) {
      *outStream_ << "My dbc numbers: ";
      for(unsigned i=0; i<my_dbc_.size(); ++i) {
        *outStream_ << my_dbc_[i];
      }
      *outStream_ << std::endl;
    }
    //
    ROL::Ptr<std::vector<std::vector<Intrepid::FieldContainer<GO> > > > dirichletSideSets = meshMgr_->getSideSets();
    std::vector<std::vector<Intrepid::FieldContainer<GO> > > &dss = *dirichletSideSets;
    Teuchos::Array<GO> mySortedCellIds_(myCellIds_);
    std::sort(mySortedCellIds_.begin(), mySortedCellIds_.end());
    mySortedCellIds_.erase( std::unique(mySortedCellIds_.begin(), mySortedCellIds_.end()), mySortedCellIds_.end() );
    
    myBoundaryCellIds_.resize(dss[0].size());
    for (int i=0; i<static_cast<int>(dss[0].size()); ++i) {
      for (int j=0; j<dss[0][i].dimension(0); ++j) {
        if (std::binary_search(mySortedCellIds_.begin(), mySortedCellIds_.end(), dss[0][i](j))) {
          myBoundaryCellIds_[i].push_back(dss[0][i](j));
        }
      }
    }
    
    cte_ = *(meshMgr_->getCellToEdgeMap());
    Intrepid::FieldContainer<GO>  &nodeDofs = *(dofMgr_->getNodeDofs());
    //Intrepid::FieldContainer<int>  &edgeDofs = *(dofMgr_->getEdgeDofs());
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
    // vector field
    int nfields = basisPtrs_.size();
    numDofsPerNode = numDofsPerNode * nfields;
    numDofsPerEdge = numDofsPerEdge * nfields;
    	
    Intrepid::FieldContainer<Real> &nodes = *meshMgr_->getNodes();
    int n_dbc = my_dbc_.size();
    for (int i=0; i<static_cast<int>(myBoundaryCellIds_.size()); ++i) {
      bool isdbc = false;
      for(int i_dbc = 0; i_dbc < n_dbc; i_dbc++) {
        if(i == my_dbc_[i_dbc]) {
          isdbc = true;
          break;
        }
      }
      if(!isdbc)
        continue;	
      
      for (int j=0; j<myBoundaryCellIds_[i].size(); ++j) {
/*	
	int ifDBCCell = check_DBC_By_Coords(myBoundaryCellIds_[i][j], i);
	if(ifDBCCell < 1)
		continue;
*/      
	std::vector<Real> x(spaceDim_);
        const CellTopologyData * ctd = cellType_.getCellTopologyData();
        Teuchos::ArrayView<unsigned> locNodes(const_cast<unsigned *>(ctd->subcell[spaceDim_-1][i].node), cellType_.getVertexCount(spaceDim_-1, i));
        for (int l=0; l<static_cast<int>(cellType_.getVertexCount(spaceDim_-1, i)); ++l) {
          x[0]=nodes(ctn_(myBoundaryCellIds_[i][j], locNodes[l]), 0);
          x[1]=nodes(ctn_(myBoundaryCellIds_[i][j], locNodes[l]), 1);
	  // use coordinates to check if a DOF is DBC DOF
	  int ifDBCNode = check_DBC_Coords_Range( x );
	  if(ifDBCNode < 0){
            continue;
          }
          else if(ifDBCNode == 0){
            if ( verbose_ ) {
              *outStream_ << "DBC node: "
                          << ctn_(myBoundaryCellIds_[i][j], locNodes[l])
                          << ", fixing X direction" << std::endl;
            }
            myDirichletDofs_.push_back(nodeDofs(ctn_(myBoundaryCellIds_[i][j], locNodes[l]), 0));
          }
          else if(ifDBCNode == 1 && numDofsPerNode >= 2){
            if ( verbose_ ) {
              *outStream_ << "DBC node: "
                          << ctn_(myBoundaryCellIds_[i][j], locNodes[l])
                          << ", fixing Y direction" << std::endl;
            }
            myDirichletDofs_.push_back(nodeDofs(ctn_(myBoundaryCellIds_[i][j], locNodes[l]), 1));
          }
          else {
            if ( verbose_ ) {
              *outStream_ << "DBC node: "
                          << ctn_(myBoundaryCellIds_[i][j], locNodes[l])
                          << ", fixing ALL direction" << std::endl;
            }
            for (int k=0; k<numDofsPerNode; ++k) {
              myDirichletDofs_.push_back(nodeDofs(ctn_(myBoundaryCellIds_[i][j], locNodes[l]), k));
            }
          }
        }
	// edge dofs are NOT in use	
	/*
        for (int k=0; k<numDofsPerEdge; ++k) {
           myDirichletDofs_.push_back(edgeDofs(cte_(myBoundaryCellIds_[i][j], i), k));
        }
	*/
      }
    }
    std::sort(myDirichletDofs_.begin(), myDirichletDofs_.end());
    myDirichletDofs_.erase( std::unique(myDirichletDofs_.begin(), myDirichletDofs_.end()), myDirichletDofs_.end() );
  }

  virtual int check_DBC_Coords_Range( const std::vector<Real> &x ) const {
    // return value :
    // -1: not a DBC node
    //  0: x direction fixed
    //  1: y direction fixed
    //  5: both x, y direction are fixed
    return 5;
  }
//
//
//
  virtual void process_loading_information(const Teuchos::RCP<Teuchos::ParameterList> &parlist) {}
//
//note that the point load is only allowed to applied on the boundary, not inside the body! 2D only
  virtual void check_and_Apply_PointLoad_By_Coords(int globalCellNum, int pointload_bc) {
    Intrepid::FieldContainer<Real> &nodes = *meshMgr_->getNodes();
    const CellTopologyData * ctd = cellType_.getCellTopologyData();
    Teuchos::ArrayView<unsigned> locNodes(const_cast<unsigned *>(ctd->subcell[spaceDim_-1][pointload_bc].node), cellType_.getVertexCount(spaceDim_-1, pointload_bc));
    std::vector<Real> x1(spaceDim_);
    std::vector<Real> x2(spaceDim_);
    std::vector<int > localNodeNum(2);
    x1[0]=nodes(ctn_(globalCellNum, locNodes[0]), 0);
    x1[1]=nodes(ctn_(globalCellNum, locNodes[0]), 1);
    x2[0]=nodes(ctn_(globalCellNum, locNodes[1]), 0);
    x2[1]=nodes(ctn_(globalCellNum, locNodes[1]), 1);
    ApplyPointLoad(pointload_bc, globalCellNum, localNodeNum, x1, x2);
  }

  virtual void ApplyPointLoad(const int pointload_bc,
                              const GO globalCellNum,
                              const std::vector<int> &localNodeNum,
                              const std::vector<Real> &coord1,
                              const std::vector<Real> &coord2) { 
    	//Intrepid::FieldContainer<int>  &nodeDofs = *(dofMgr_->getNodeDofs());
	//bool isLoadPosContainedInCurrentSegment = false;
	//int whichNodeIsCloserToPos = -1;
	return;
  }
//
//
//
  virtual void EnforceDBC(void) {
    // Apply Dirichlet conditions.
    // zero out row and column, make matrix symmetric
    //ROL::Ptr<Tpetra::Details::DefaultTypes::node_type> node = matA_->getNode();
    //matA_dirichlet_ = matA_->clone(node);
    //matM_dirichlet_ = matM_->clone(node);
    //vecF_dirichlet_ = ROL::makePtr<Tpetra::MultiVector<>>(matA_->getRangeMap(), 1, true);
    //Tpetra::deep_copy(*vecF_dirichlet_, *vecF_);
 
    matA_dirichlet_->resumeFill();
    //matM_dirichlet_->resumeFill();
 
    GO gDof;
    for(int i=0; i<numCells_; i++) {
      for(int j=0; j<numLocalDofs_; j++) {
        gDof = cellDofs_(myCellIds_[i], j);
        if (myUniqueMap_->isNodeGlobalElement(gDof) && check_myGlobalDof_on_boundary(gDof)) {
          size_t numRowEntries = matA_dirichlet_->getNumEntriesInGlobalRow(gDof);
          Teuchos::Array<GO> indices(numRowEntries, 0);
          Teuchos::Array<Real> values(numRowEntries, 0);
          Teuchos::Array<Real> canonicalValues(numRowEntries, 0);
          Teuchos::Array<Real> zeroValues(numRowEntries, 0);
          matA_dirichlet_->getGlobalRowCopy(gDof, indices, values, numRowEntries);
          //matM_dirichlet_->getGlobalRowCopy(gDof, indices, values, numRowEntries);
          for (int k=0; k<indices.size(); k++) {
            if (indices[k] == gDof) {
              canonicalValues[k] = 1.0;
            }
          }
          matA_dirichlet_->replaceGlobalValues(gDof, indices, canonicalValues);
          //matM_dirichlet_->replaceGlobalValues(gDof, indices, zeroValues);
          vecF_dirichlet_->replaceGlobalValue (gDof, 0, 0);
        }
 
        if (!check_myGlobalDof_on_boundary(gDof)) {
          size_t numDBCDofs = myDirichletDofs_.size();
          Teuchos::Array<GO> indices(numDBCDofs, 0);
          Teuchos::Array<Real> zeroValues(numDBCDofs, 0);
          for (int k=0; k<indices.size(); k++) {
            indices[k] = myDirichletDofs_[k];
          }
          matA_dirichlet_->replaceGlobalValues(gDof, indices, zeroValues);
          //matM_dirichlet_->replaceGlobalValues(gDof, indices, zeroValues);
        }
      }
    }
    matA_dirichlet_->fillComplete();
    //matM_dirichlet_->fillComplete();

//    GenerateTransposedMats();
//    ConstructSolvers();
//    ConstructAdjointSolvers();
  }

//  virtual void MatrixRemoveDBC(void) {
//    // Apply Dirichlet conditions.
//    // zero out row and column, make matrix symmetric
//    //ROL::Ptr<Tpetra::Details::DefaultTypes::node_type> node = matA_->getNode();
//    //matA_dirichlet_ = matA_->clone(node);
//    //matM_dirichlet_ = matM_->clone(node);
// 
//    matA_dirichlet_->resumeFill();
//    matM_dirichlet_->resumeFill();
// 
//    GO gDof;
//    for(int i=0; i<numCells_; i++) {
//      for(int j=0; j<numLocalDofs_; j++) {
//        gDof = cellDofs_(myCellIds_[i], j);
//        if (myUniqueMap_->isNodeGlobalElement(gDof) && check_myGlobalDof_on_boundary(gDof)) {
//          size_t numRowEntries = matA_dirichlet_->getNumEntriesInGlobalRow(gDof);
//          Teuchos::Array<GO> indices(numRowEntries, 0);
//          Teuchos::Array<Real> values(numRowEntries, 0);
//          Teuchos::Array<Real> canonicalValues(numRowEntries, 0);
//          Teuchos::Array<Real> zeroValues(numRowEntries, 0);
//          matA_dirichlet_->getGlobalRowCopy(gDof, indices, values, numRowEntries);
//          matM_dirichlet_->getGlobalRowCopy(gDof, indices, values, numRowEntries);
//          for (int k=0; k<indices.size(); k++) {
//            if (indices[k] == gDof) {
//              canonicalValues[k] = 1.0;
//            }
//          }
//          matA_dirichlet_->replaceGlobalValues(gDof, indices, canonicalValues);
//          matM_dirichlet_->replaceGlobalValues(gDof, indices, zeroValues);
//        }
// 
//        if (!check_myGlobalDof_on_boundary(gDof)) {
//          size_t numDBCDofs = myDirichletDofs_.size();
//          Teuchos::Array<GO> indices(numDBCDofs, 0);
//          Teuchos::Array<Real> zeroValues(numDBCDofs, 0);
//          for (int k=0; k<indices.size(); k++) {
//            indices[k] = myDirichletDofs_[k];
//          }
//          matA_dirichlet_->replaceGlobalValues(gDof, indices, zeroValues);
//          matM_dirichlet_->replaceGlobalValues(gDof, indices, zeroValues);
//        }
//      }
//    }
//    matA_dirichlet_->fillComplete();
//    matM_dirichlet_->fillComplete();
//  }

  virtual void MatrixRemoveDBC(void) {
    // Apply Dirichlet conditions.
    // zero out row and column, make matrix symmetric
    matA_dirichlet_->resumeFill();
    //matM_dirichlet_->resumeFill();
 
    for (GO i=0; i<myDirichletDofs_.size(); ++i) {
      if (myUniqueMap_->isNodeGlobalElement(myDirichletDofs_[i])) {
        size_t numRowEntries = matA_dirichlet_->getNumEntriesInGlobalRow(myDirichletDofs_[i]);
        Teuchos::Array<GO> indices(numRowEntries, 0);    
        Teuchos::Array<Real> values(numRowEntries, 0);
        Teuchos::Array<Real> canonicalValues(numRowEntries, 0);    
        Teuchos::Array<Real> zeroValues(numRowEntries, 0);    
        matA_dirichlet_->getGlobalRowCopy(myDirichletDofs_[i], indices, values, numRowEntries);
        //matM_dirichlet_->getGlobalRowCopy(myDirichletDofs_[i], indices, values, numRowEntries);
        for (int j=0; j < static_cast<int>(indices.size()); ++j) {
          if (myDirichletDofs_[i] == indices[j]) {
            canonicalValues[j] = 1.0;
          }
        }
        matA_dirichlet_->replaceGlobalValues(myDirichletDofs_[i], indices, canonicalValues);
        //matM_dirichlet_->replaceGlobalValues(myDirichletDofs_[i], indices, zeroValues);
        for (int j=0; j < static_cast<int>(indices.size()); ++j) {
          Teuchos::Array<GO> ind(1, myDirichletDofs_[i]);
          Teuchos::Array<Real> valA(1, canonicalValues[j]); 
          Teuchos::Array<Real> valM(1, zeroValues[j]); 
          matA_dirichlet_->replaceGlobalValues(indices[j], ind, valA);
          //matM_dirichlet_->replaceGlobalValues(indices[j], ind, valM);
        }
      }
    }
    matA_dirichlet_->fillComplete();
    //matM_dirichlet_->fillComplete();

//    GenerateTransposedMats();
//    ConstructSolvers();
//    ConstructAdjointSolvers();
  }

  virtual void VectorRemoveDBC(void) {
    // Apply Dirichlet conditions.
    //vecF_dirichlet_ = ROL::makePtr<Tpetra::MultiVector<>>(matA_->getRangeMap(), 1, true);
    //Tpetra::deep_copy(*vecF_dirichlet_, *vecF_);
 
    GO gDof(0);
    for(GO i=0; i<numCells_; i++) {
      for(int j=0; j<numLocalDofs_; j++) {
        gDof = cellDofs_(myCellIds_[i], j);
        if (myUniqueMap_->isNodeGlobalElement(gDof) && check_myGlobalDof_on_boundary(gDof)) {
          vecF_dirichlet_->replaceGlobalValue (gDof, 0, 0);
        }
      }
    }
  }

  virtual void GenerateTransposedMats() {
    // Create matrix transposes.
    Tpetra::RowMatrixTransposer<> transposerA(matA_dirichlet_);
    //Tpetra::RowMatrixTransposer<> transposerM(matM_dirichlet_);
    matA_dirichlet_trans_ = transposerA.createTranspose();
    //matM_dirichlet_trans_ = transposerM.createTranspose();
  }


  virtual void ConstructSolvers(void) {
    // Construct solver using Amesos2 factory.
    #ifdef ROL_TIMERS
    Teuchos::TimeMonitor LocalTimer(*FactorTime_example_PDEOPT_TOOLS_PDEFEM_GLOB);
    #endif
    try {
      solverA_ = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("KLU2", matA_dirichlet_);
      //solverA_ = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("PardisoMKL", matA_dirichlet_);
      //solverA_ = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("SuperLU", matA_dirichlet_);
    }
    catch (std::invalid_argument& e) {
      std::cout << e.what() << std::endl;
    }
    // Perform factorization.
    solverA_->numericFactorization();
  }

  virtual void ConstructAdjointSolvers() {
    // Construct solver using Amesos2 factory.
    try{
      solverA_trans_ = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("KLU2", matA_dirichlet_trans_);
      //solverA_trans_ = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("PardisoMKL", matA_dirichlet_trans_);
      //solverA_trans_ = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >("SuperLU", matA_dirichlet_trans_);
    } catch (std::invalid_argument& e) {
      std::cout << e.what() << std::endl;
    }
    solverA_trans_->numericFactorization();
  }

  ROL::Ptr<Tpetra::CrsMatrix<> > getMatA(const bool &transpose = false) const {
    if (transpose) {
      //return matA_dirichlet_trans_;
      return matA_dirichlet_;
    }
    else {
      return matA_dirichlet_;
    }
  }


  ROL::Ptr<Tpetra::CrsMatrix<> > getMatB(const bool &transpose = false) const {
    if (transpose) {
      //return matM_dirichlet_trans_;
      return matM_dirichlet_;
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
      //return solverA_trans_;
      return solverA_;
    }
    else {
      return solverA_;
    }
  }

 virtual Real funcTarget(const Real &x1, const Real &x2) const { return 0.0; }


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
      for (GO i=0; i<cellToNodeMap.dimension(0); ++i) {
        for (GO j=0; j<cellToNodeMap.dimension(1); ++j) {
          meshfile << cellToNodeMap(i,j) << "  ";
        }
        meshfile << std::endl;
      }
      meshfile.close();
      
      meshfile.open("cell_to_node_tri.txt");
      for (GO i=0; i<cellToNodeMap.dimension(0); ++i) {
        for (GO j=0; j<3; ++j) {
          meshfile << cellToNodeMap(i,j) << "  ";
        }
        meshfile << std::endl;
        for (GO j=2; j<5; ++j) {
          meshfile << cellToNodeMap(i,j%4) << "  ";
        }
        meshfile << std::endl;
      }
      meshfile.close();
     
      meshfile.open("nodes.txt");
      meshfile.precision(16);
      for (GO i=0; i<nodes.dimension(0); ++i) {
        for (GO j=0; j<nodes.dimension(1); ++j) {
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
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> >   matWriter;
    matWriter.writeSparseFile("stiffness_mat", matA_);
    matWriter.writeSparseFile("dirichlet_mat", matA_dirichlet_);
    //matWriter.writeSparseFile("mass_mat", matM_);
    matWriter.writeDenseFile("Ud_vec", vecUd_);
  }


  void outputTpetraVector(const ROL::Ptr<const Tpetra::MultiVector<> > &vec,
                          const std::string &filename) const {
    Tpetra::MatrixMarket::Writer< Tpetra::CrsMatrix<> > vecWriter;
    vecWriter.writeDenseFile(filename, vec);
  }

  // ACCESSOR FUNCTIONS
  ROL::Ptr<MeshManager<Real> >& GetMeshManager(void) {
    TEUCHOS_TEST_FOR_EXCEPTION(meshMgr_==ROL::nullPtr, std::logic_error,
      ">>> (PDE_FEM::GetMeshManager): MeshManager not initialized!");
    return meshMgr_;
  }

  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > >& GetBasisPtr(const int ind) {
    TEUCHOS_TEST_FOR_EXCEPTION(ind > spaceDim_-1 || ind < 0, std::logic_error,
      ">>> (PDE_FEM::GetBasisPtr): ind out of bounds!");
    TEUCHOS_TEST_FOR_EXCEPTION(basisPtrs_.size()==0, std::logic_error,
      ">>> (PDE_FEM::GetBasisPtr): BasisPntrs not initialized!");
    return basisPtrs_[ind];
  }

  int GetBasisOrder(void) const {
    return basisOrder_;
  }

  int GetSpaceDim(void) const {
    return spaceDim_;
  }

  GO GetNumCells(void) const {
    return numCells_;
  }

  int GetNumLocalDofs(void) const {
    return numLocalDofs_;
  }

  int GetNumCubPoints(void) const {
    return numCubPoints_;
  }

  int GetLocalFieldSize(void) const {
    return lfs_;
  }

}; // class PDE_FEM

#endif
