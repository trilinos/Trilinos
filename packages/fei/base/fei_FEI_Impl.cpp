/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <fei_utils.hpp>

#include <fei_FEI_Impl.hpp>
#include <fei_Record.hpp>
#include <fei_TemplateUtils.hpp>
#include <fei_ParameterSet.hpp>
#include <fei_base.hpp>

#include <fei_Pattern.hpp>
#include <fei_LibraryWrapper.hpp>
#include <fei_Data.hpp>
#include <fei_defs.h>

#include <stdexcept>
#include <cmath>

#undef fei_file
#define fei_file "fei_FEI_Impl.cpp"

#include <fei_ErrMacros.hpp>

fei::FEI_Impl::FEI_Impl(fei::SharedPtr<LibraryWrapper> wrapper,
		      MPI_Comm comm,
		      int masterRank)
  : wrapper_(1),
    nodeIDType_(0),
    elemIDType_(1),
    constraintIDType_(2),
    factory_(1),
    rowSpace_(NULL),
    matGraph_(),
    x_(),
    b_(),
    A_(),
    linSys_(),
    newData_(false),
    soln_fei_matrix_(NULL),
    soln_fei_vector_(NULL),
    comm_(comm),
    masterRank_(masterRank),
    localProc_(0),
    numProcs_(1),
    numParams_(0),
    paramStrings_(NULL),
    matrixIDs_(),
    rhsIDs_(),
    matScalars_(),
    matScalarsSet_(false),
    rhsScalars_(),
    rhsScalarsSet_(false),
    constraintID_(0),
    index_soln_(0),
    index_current_(0),
    index_current_rhs_row_(0),
    solveType_(0),
    iterations_(0),
    setSolveTypeCalled_(false),
    initPhaseIsComplete_(false),
    aggregateSystemFormed_(false),
    newMatrixDataLoaded_(0),
    solveCounter_(0),
    initTime_(0.0),
    loadTime_(0.0),
    solveTime_(0.0),
    solnReturnTime_(0.0),
    iwork_(),
    nodeset_(),
    nodeset_filled_(false),
    block_dof_per_elem_(),
    any_blocks_have_elem_dof_(false)
{
  wrapper_[0] = wrapper;

  factory_[0].reset(new snl_fei::Factory(comm, wrapper_[0]));
  createdFactory_ = true;

  basic_initializations();
}

fei::FEI_Impl::FEI_Impl(const fei::Factory* factory,
		      MPI_Comm comm,
		      int masterRank)
  : wrapper_(1),
    nodeIDType_(0),
    elemIDType_(1),
    constraintIDType_(2),
    factory_(1),
    createdFactory_(false),
    rowSpace_(NULL),
    matGraph_(),
    x_(),
    b_(),
    A_(),
    linSys_(),
    newData_(false),
    soln_fei_matrix_(NULL),
    soln_fei_vector_(NULL),
    comm_(comm),
    masterRank_(masterRank),
    localProc_(0),
    numProcs_(1),
    numParams_(0),
    paramStrings_(NULL),
    matrixIDs_(),
    rhsIDs_(),
    matScalars_(),
    matScalarsSet_(false),
    rhsScalars_(),
    rhsScalarsSet_(false),
    constraintID_(0),
    index_soln_(0),
    index_current_(0),
    index_current_rhs_row_(0),
    solveType_(0),
    iterations_(0),
    setSolveTypeCalled_(false),
    initPhaseIsComplete_(false),
    aggregateSystemFormed_(false),
    newMatrixDataLoaded_(0),
    solveCounter_(0),
    initTime_(0.0),
    loadTime_(0.0),
    solveTime_(0.0),
    solnReturnTime_(0.0),
    iwork_(),
    nodeset_(),
    nodeset_filled_(false),
    block_dof_per_elem_(),
    any_blocks_have_elem_dof_(false)
{
  wrapper_[0].reset(0);

  const snl_fei::Factory* snlfactory
    = dynamic_cast<const snl_fei::Factory*>(factory);
  if (snlfactory != NULL) {
    wrapper_[0] = snlfactory->get_LibraryWrapper();
  }

  factory_[0] = factory->clone();

  basic_initializations();
}

fei::FEI_Impl::~FEI_Impl()
{
  for(int k=0; k<numParams_; k++) {
    delete [] paramStrings_[k];
  }
  delete [] paramStrings_;

  if (soln_fei_matrix_ != NULL && wrapper_[0].get() != NULL) {
    fei::SharedPtr<LinearSystemCore> lsc = wrapper_[0]->getLinearSystemCore();
    if (lsc.get() != NULL) {
      lsc->destroyMatrixData(*soln_fei_matrix_);
      delete soln_fei_matrix_;
      soln_fei_matrix_ = NULL;
    }
  }

  if (soln_fei_vector_ != NULL && wrapper_[0].get() != NULL) {
    fei::SharedPtr<LinearSystemCore> lsc = wrapper_[0]->getLinearSystemCore();
    if (lsc.get() != NULL) {
      lsc->destroyVectorData(*soln_fei_vector_);
      delete soln_fei_vector_;
      soln_fei_vector_ = NULL;
    }
  }
}

void fei::FEI_Impl::basic_initializations()
{
  localProc_ = fei::localProc(comm_);
  numProcs_ = fei::numProcs(comm_);

  constraintID_ = localProc_*100000;

  matrixIDs_.resize(1);
  matrixIDs_[0] = 0;
  A_.resize(1);
  rhsIDs_.resize(1);
  rhsIDs_[0] = 0;
  b_.resize(1);

  rowSpace_ = factory_[0]->createVectorSpace(comm_, (const char*)NULL);

  rowSpace_->defineIDTypes(1, &nodeIDType_);
  rowSpace_->defineIDTypes(1, &elemIDType_);
  rowSpace_->defineIDTypes(1, &constraintIDType_);

  matGraph_ = factory_[0]->createMatrixGraph(rowSpace_, rowSpace_,(const char*)NULL);
  if (matGraph_.get() == NULL) {
    voidERReturn;
  }
}

int fei::FEI_Impl::setIDLists(int numMatrices,
			       const int* matrixIDs,
			       int numRHSs,
			       const int* rhsIDs)
{
  if (numMatrices < 1 || numRHSs < 1) {
    fei::console_out() << "fei::FEI_Impl::setIDLists ERROR, numMatrices and numRHSs "
	 << "must both be greater than 0."<<FEI_ENDL;
    ERReturn(-1);
  }

  matrixIDs_.resize(0);
  A_.resize(numMatrices);
  for(int i=0; i<numMatrices; ++i) {
    fei::sortedListInsert(matrixIDs[i], matrixIDs_);
  }
  if ((int)matrixIDs_.size() != numMatrices) {
    fei::console_out() << "fei::FEI_Impl::setIDLists ERROR creating matrixIDs_ list."<<FEI_ENDL;
    ERReturn(-1);
  }

  rhsIDs_.resize(0);
  b_.resize(numRHSs);
  for(int i=0; i<numRHSs; ++i) {
    fei::sortedListInsert(rhsIDs[i], rhsIDs_);
  }
  if ((int)rhsIDs_.size() != numRHSs) {
    fei::console_out() << "fei::FEI_Impl::setIDLists ERROR creating rhsIDs_ list."<<FEI_ENDL;
    ERReturn(-1);
  }

  if (wrapper_[0].get() != NULL) {
    fei::SharedPtr<LinearSystemCore> linSysCore = wrapper_[0]->getLinearSystemCore();
    if (linSysCore.get() != NULL) {
      linSysCore->setNumRHSVectors(rhsIDs_.size(), &rhsIDs_[0]);
    }
  }

  return(0);
}

fei::SharedPtr<fei::LinearSystem> fei::FEI_Impl::getLinearSystem()
{
  return( linSys_ );
}

int fei::FEI_Impl::parameters(int numParams, 
			       const char *const* paramStrings)
{
  // merge these parameters with any others we may have, for later use.
  snl_fei::mergeStringLists(paramStrings_, numParams_,
				   paramStrings, numParams);

  if (wrapper_[0].get() != NULL) {
    if (wrapper_[0]->haveLinearSystemCore()) {
      CHK_ERR( wrapper_[0]->getLinearSystemCore()->parameters(numParams, (char**)paramStrings) );
    }
    if (wrapper_[0]->haveFiniteElementData()) {
      CHK_ERR( wrapper_[0]->getFiniteElementData()->parameters(numParams, (char**)paramStrings) );
    }
  }

  std::vector<std::string> stdstrings;
  fei::utils::char_ptrs_to_strings(numParams, paramStrings, stdstrings);
  fei::ParameterSet paramset;
  fei::utils::parse_strings(stdstrings, " ", paramset);
  factory_[0]->parameters(paramset);

  return(0);
}

int fei::FEI_Impl::setSolveType(int solveType)
{
  solveType_ = solveType;
  setSolveTypeCalled_ = true;

  return(0);
}

int fei::FEI_Impl::initFields(int numFields, 
			       const int *fieldSizes, 
			       const int *fieldIDs,
             const int *fieldTypes)
{
  rowSpace_->defineFields(numFields, fieldIDs, fieldSizes, fieldTypes);

  return(0);
}

int fei::FEI_Impl::initElemBlock(GlobalID elemBlockID,
				  int numElements,
				  int numNodesPerElement,
				  const int* numFieldsPerNode,
				  const int* const* nodalFieldIDs,
				  int numElemDofFieldsPerElement,
				  const int* elemDOFFieldIDs,
				  int interleaveStrategy)
{
  //define pattern that describes the layout of fields for elements in
  //this element-block

  int numIDs = numNodesPerElement;
  if (numElemDofFieldsPerElement > 0) ++numIDs;
  std::vector<int> idTypes;
  std::vector<int> numFieldsPerID;
  std::vector<int> fieldIDs;

  int i, j;
  for(i=0; i<numNodesPerElement; ++i) {
    idTypes.push_back(nodeIDType_);
    numFieldsPerID.push_back(numFieldsPerNode[i]);

    for(j=0; j<numFieldsPerNode[i]; ++j) {
      fieldIDs.push_back(nodalFieldIDs[i][j]);
    }
  }

  if (numElemDofFieldsPerElement>0) {
    idTypes.push_back(elemIDType_);
    numFieldsPerID.push_back(numElemDofFieldsPerElement);
    for(i=0; i<numElemDofFieldsPerElement; ++i) {
      fieldIDs.push_back(elemDOFFieldIDs[i]);
    }

    block_dof_per_elem_.insert(std::pair<int,int>(elemBlockID, numElemDofFieldsPerElement));
    any_blocks_have_elem_dof_ = true;
  }

  int pattID = matGraph_->definePattern(numIDs,
			   &idTypes[0],
			   &numFieldsPerID[0],
			   &fieldIDs[0]);

  //initialize connectivity-block
  CHK_ERR( matGraph_->initConnectivityBlock(elemBlockID, numElements,
					    pattID) );

  return(0);
}

int fei::FEI_Impl::initElem(GlobalID elemBlockID,
			     GlobalID elemID,
			     const GlobalID* elemConn)
{
  bool elemdof = false;

  if (any_blocks_have_elem_dof_) {
    std::map<int,int>::const_iterator
      b_iter = block_dof_per_elem_.find(elemBlockID);
    if (b_iter != block_dof_per_elem_.end()) {
      int numIDs = matGraph_->getNumIDsPerConnectivityList(elemBlockID);
      iwork_.resize(numIDs);
      int* iworkPtr = &iwork_[0];
      for(int i=0; i<numIDs-1; ++i) {
	iworkPtr[i] = elemConn[i];
      }
      iworkPtr[numIDs-1] = elemID;

      CHK_ERR( matGraph_->initConnectivity(elemBlockID, elemID, iworkPtr) );
      elemdof = true;
    }
  }

  if (!elemdof) {
    CHK_ERR( matGraph_->initConnectivity(elemBlockID, elemID, elemConn) );
  }

  nodeset_filled_ = false;

  return(0);
}

int fei::FEI_Impl::initSlaveVariable(GlobalID slaveNodeID, 
				      int slaveFieldID,
				      int offsetIntoSlaveField,
				      int numMasterNodes,
				      const GlobalID* masterNodeIDs,
				      const int* masterFieldIDs,
				      const double* weights,
				      double rhsValue)
{
  throw std::runtime_error("FEI_Impl::initSlaveVariable not implemented.");
  return(0);
}

int fei::FEI_Impl::deleteMultCRs()
{
  throw std::runtime_error("FEI_Impl::deleteMultCRs not implemented.");
  return(0);
}

int fei::FEI_Impl::initSharedNodes(int numSharedNodes,
				    const GlobalID *sharedNodeIDs,  
				    const int* numProcsPerNode, 
				    const int *const *sharingProcIDs)
{
  CHK_ERR( rowSpace_->initSharedIDs(numSharedNodes,
				    nodeIDType_,
				    sharedNodeIDs,
				    numProcsPerNode,
				    sharingProcIDs) );

  return(0);
}

int fei::FEI_Impl::initCRMult(int numCRNodes,
			       const GlobalID* CRNodeIDs,
			       const int *CRFieldIDs,
			       int& CRID)
{
  iwork_.assign(numCRNodes,nodeIDType_);

  CRID = constraintID_++;

  CHK_ERR( matGraph_->initLagrangeConstraint(CRID,
					     constraintIDType_,
					     numCRNodes,
					     &iwork_[0],
					     CRNodeIDs,
					     CRFieldIDs) );

  return(0);
}

int fei::FEI_Impl::initCRPen(int numCRNodes,
			      const GlobalID* CRNodeIDs, 
			      const int *CRFieldIDs,
			      int& CRID)
{
  iwork_.assign(numCRNodes, nodeIDType_);

  CRID = constraintID_++;

  CHK_ERR( matGraph_->initPenaltyConstraint(CRID,
					    constraintIDType_,
					    numCRNodes,
					    &iwork_[0],
					    CRNodeIDs,
					    CRFieldIDs) );

  return(0);
}

int fei::FEI_Impl::initComplete()
{
  CHK_ERR( matGraph_->initComplete() );

  if (matrixIDs_.size() < 1 || factory_.size() < 1 ||
      A_.size() < 1 || b_.size() < 1) {
    ERReturn(-1);
  }

  A_[0] = factory_[0]->createMatrix(matGraph_);

  std::vector<std::string> stdstrings;
  fei::utils::char_ptrs_to_strings(numParams_, paramStrings_, stdstrings);
  fei::ParameterSet params;
  fei::utils::parse_strings(stdstrings, " ", params);

  CHK_ERR( A_[0]->parameters(params) );

  b_[0] = factory_[0]->createVector(matGraph_);

  if (matrixIDs_.size() > 1) {
    bool multiple_factories = false;

    if (wrapper_[0].get() != NULL) {
      multiple_factories = true;

      fei::SharedPtr<LinearSystemCore> linsyscore = wrapper_[0]->getLinearSystemCore();
      if (linsyscore.get() == NULL) {
	fei::console_out() << "fei::FEI_Impl ERROR, multiple matrix/rhs assembly not supported "
	     << "non-null LibraryWrapper holds null LinearSystemCore."<<FEI_ENDL;
	ERReturn(-1);
      }

      wrapper_.resize(matrixIDs_.size());
      factory_.resize(matrixIDs_.size());
      for(unsigned i=1; i<matrixIDs_.size(); ++i) {
	fei::SharedPtr<LinearSystemCore> lscclone(linsyscore->clone());
	wrapper_[i].reset(new LibraryWrapper(lscclone));
	factory_[i].reset(new snl_fei::Factory(comm_, wrapper_[i]));
      }
    }

    fei::SharedPtr<fei::Factory> factory;
    for(unsigned i=1; i<matrixIDs_.size(); ++i) {
      factory = multiple_factories ? factory_[i] : factory_[0];

      A_[i] = factory->createMatrix(matGraph_);
      CHK_ERR( A_[i]->parameters(params) );
    }

    for(unsigned i=1; i<rhsIDs_.size(); ++i) {
      factory = multiple_factories ? factory_[i] : factory_[0];

      b_[i] = factory->createVector(matGraph_);
    }

    if (wrapper_[0].get() != NULL) {
      fei::SharedPtr<LinearSystemCore> linsyscore
        = wrapper_[0]->getLinearSystemCore();

      linsyscore->setNumRHSVectors(1, &(rhsIDs_[0]));

      unsigned num = rhsIDs_.size();
      if (matrixIDs_.size() < num) num = matrixIDs_.size();
      for(unsigned i=1; i<num; ++i) {
	fei::SharedPtr<LinearSystemCore> lsc = wrapper_[i]->getLinearSystemCore();

	if (i==num-1 && rhsIDs_.size() > num) {
	  int numRHSs = rhsIDs_.size() - matrixIDs_.size() + 1;
	  lsc->setNumRHSVectors(numRHSs, &(rhsIDs_[i]));
	}
	else {
	  lsc->setNumRHSVectors(1, &(rhsIDs_[i]));
	}
      }

      if (rhsIDs_.size() < matrixIDs_.size()) {
	int dummyID = -1;
	for(unsigned i=rhsIDs_.size(); i<matrixIDs_.size(); ++i) {
	  wrapper_[i]->getLinearSystemCore()->setNumRHSVectors(1, &dummyID);
	}
      }
    }

    for(unsigned i=1; i<matrixIDs_.size(); ++i) {
      factory = multiple_factories ? factory_[i] : factory_[0];

      A_[i] = factory->createMatrix(matGraph_);
      CHK_ERR( A_[i]->parameters(params) );
    }

    for(unsigned i=1; i<rhsIDs_.size(); ++i) {
      factory = multiple_factories ? factory_[i] : factory_[0];

      b_[i] = factory->createVector(matGraph_);
    }
  }

  x_ = factory_[0]->createVector(matGraph_, true);

  linSys_ = factory_[0]->createLinearSystem(matGraph_);

  CHK_ERR( linSys_->parameters(numParams_, paramStrings_) );

  linSys_->setMatrix(A_[0]);
  linSys_->setRHS(b_[0]);
  linSys_->setSolutionVector(x_);

  return(0);
}

int fei::FEI_Impl::setCurrentMatrix(int matID)
{
  std::vector<int>::const_iterator
    iter = std::lower_bound(matrixIDs_.begin(), matrixIDs_.end(), matID);
  if (iter == matrixIDs_.end() || *iter != matID) {
    fei::console_out() << "fei::FEI_Impl::setCurrentMatrix: matID ("<<matID
       <<") not found." <<FEI_ENDL;
    return(-1);
  }

  index_current_ = iter-matrixIDs_.begin();

  return(0);
}

int fei::FEI_Impl::setCurrentRHS(int rhsID)
{
  std::vector<int>::const_iterator
    iter = std::lower_bound(rhsIDs_.begin(), rhsIDs_.end(), rhsID);
  if (iter == rhsIDs_.end() || *iter != rhsID) {
    fei::console_out() << "fei::FEI_Impl::setCurrentRHS: rhsID ("<<rhsID<<") not found."
	 << FEI_ENDL;
    return(-1);
  }

  index_current_rhs_row_ = iter - rhsIDs_.begin();

  return(0);
}

int fei::FEI_Impl::resetSystem(double s)
{
  int err = A_[index_current_]->putScalar(s);
  err += x_->putScalar(s);
  err += b_[index_current_rhs_row_]->putScalar(s);

  return(err);
}

int fei::FEI_Impl::resetMatrix(double s)
{
  return( A_[index_current_]->putScalar(s) );
}

int fei::FEI_Impl::resetRHSVector(double s)
{
  return( b_[index_current_rhs_row_]->putScalar(s) );
}

int fei::FEI_Impl::resetInitialGuess(double s)
{
  return( x_->putScalar(s) );
}

int fei::FEI_Impl::loadNodeBCs(int numNodes,
                                const GlobalID *nodeIDs,
                                int fieldID,
                                const int* offsetsIntoField,
                                const double* prescribedValues)
{
  CHK_ERR( linSys_->loadEssentialBCs(numNodes, nodeIDs,
                                     nodeIDType_, fieldID,
                                     offsetsIntoField, prescribedValues) );

  newData_ = true;

  return(0);
}

int fei::FEI_Impl::loadElemBCs(int numElems,
                    const GlobalID* elemIDs,  
                    int fieldID,
                    const double *const *alpha,  
                    const double *const *beta,  
                    const double *const *gamma)
{
  throw std::runtime_error("FEI_Impl::loadElemBCs not implemented.");
  return(0);
}

int fei::FEI_Impl::sumInElem(GlobalID elemBlockID,
                 GlobalID elemID,
                 const GlobalID* elemConn,
                 const double* const* elemStiffness,
                 const double* elemLoad,
                 int elemFormat)
{
  CHK_ERR( A_[index_current_]->sumIn(elemBlockID, elemID, elemStiffness, elemFormat) );

  int num = matGraph_->getConnectivityNumIndices(elemBlockID);
  std::vector<int> indices(num);
  CHK_ERR( matGraph_->getConnectivityIndices(elemBlockID, elemID, num,
                                             &indices[0], num) );
  CHK_ERR( b_[index_current_rhs_row_]->sumIn(num, &indices[0], elemLoad, 0) );

  newData_ = true;

  return(0);
}

int fei::FEI_Impl::sumInElemMatrix(GlobalID elemBlockID,
				    GlobalID elemID,
				    const GlobalID* elemConn,
				    const double* const* elemStiffness,
				    int elemFormat)
{
  CHK_ERR( A_[index_current_]->sumIn(elemBlockID, elemID, elemStiffness, elemFormat) );

  newData_ = true;

  return(0);
}

int fei::FEI_Impl::sumInElemRHS(GlobalID elemBlockID,
				 GlobalID elemID,
				 const GlobalID* elemConn,
				 const double* elemLoad)
{
  int num = matGraph_->getConnectivityNumIndices(elemBlockID);
  std::vector<int> indices(num);
  CHK_ERR( matGraph_->getConnectivityIndices(elemBlockID, elemID, num,
                                             &indices[0], num) );
  CHK_ERR( b_[index_current_rhs_row_]->sumIn(num, &indices[0], elemLoad, 0) );

  newData_ = true;

  return(0);
}

int fei::FEI_Impl::loadCRMult(int CRMultID,
			       int numCRNodes,
			       const GlobalID* CRNodeIDs,
			       const int* CRFieldIDs,
			       const double* CRWeights,
			       double CRValue)
{
  newData_ = true;

  CHK_ERR( linSys_->loadLagrangeConstraint(CRMultID, CRWeights, CRValue) );

  return(0);
}

int fei::FEI_Impl::loadCRPen(int CRPenID,
                 int numCRNodes,
                 const GlobalID* CRNodeIDs,
                 const int* CRFieldIDs,
                 const double* CRWeights,
                 double CRValue,
                 double penValue)
{
  newData_ = true;

  CHK_ERR( linSys_->loadPenaltyConstraint(CRPenID, CRWeights, penValue, CRValue) );

  return(0);
}

int fei::FEI_Impl::putIntoRHS(int IDType,
			       int fieldID,
			       int numIDs,
			       const GlobalID* IDs,
			       const double* coefficients)
{
  CHK_ERR( inputRHS(IDType, fieldID, numIDs, IDs, coefficients, false) );

  newData_ = true;

  return(0);
}

int fei::FEI_Impl::sumIntoRHS(int IDType,
			       int fieldID,
			       int numIDs,
			       const GlobalID* IDs,
			       const double* coefficients)
{
  CHK_ERR( inputRHS(IDType, fieldID, numIDs, IDs, coefficients, true) );

  newData_ = true;

  return(0);
}

int fei::FEI_Impl::inputRHS(int IDType,
			     int fieldID,
			     int numIDs,
			     const GlobalID* IDs,
			     const double* coefficients,
			     bool sumInto)
{
  int fieldSize = rowSpace_->getFieldSize(fieldID);

  int offset = 0;
  for(int i=0; i<numIDs; ++i) {
    int globalIndex = 0;
    CHK_ERR( rowSpace_->getGlobalIndex(IDType, IDs[i], globalIndex) );

    for(int j=0; j<fieldSize; ++j) {
      int eqn = globalIndex+j;
      if (sumInto) {
	CHK_ERR( b_[index_current_rhs_row_]->sumIn(1, &eqn, &(coefficients[offset++])) );
      }
      else {
	CHK_ERR( b_[index_current_rhs_row_]->copyIn(1, &eqn, &(coefficients[offset++])) );
      }
    }
  }

  return(0);
}

int fei::FEI_Impl::setMatScalars(int numScalars,
				  const int* IDs, 
				  const double* scalars)
{
  matScalars_.resize(matrixIDs_.size());

  for(int i=0; i<numScalars; ++i) {
    std::vector<int>::const_iterator
      iter = std::lower_bound(matrixIDs_.begin(), matrixIDs_.end(), IDs[i]);
    if (iter == matrixIDs_.end() || *iter != IDs[i]) {
      continue;
    }

    unsigned index = iter - matrixIDs_.begin();
    matScalars_[index] = scalars[i];
  }

  matScalarsSet_ = true;

  return(0);
}

int fei::FEI_Impl::setRHSScalars(int numScalars,
				  const int* IDs,
				  const double* scalars)
{
  rhsScalars_.resize(rhsIDs_.size());

  for(int i=0; i<numScalars; ++i) {
    std::vector<int>::const_iterator
      iter = std::lower_bound(rhsIDs_.begin(), rhsIDs_.end(), IDs[i]);
    if (iter == rhsIDs_.end() || *iter != IDs[i]) {
      continue;
    }

    unsigned index = iter - rhsIDs_.begin();
    rhsScalars_[index] = scalars[i];
  }

  rhsScalarsSet_ = true;

  return(0);
}

int fei::FEI_Impl::loadComplete(bool applyBCs,
                                bool globalAssemble)
{
  if (!newData_) {
    return(0);
  }

  if (linSys_.get() == NULL) {
    fei::console_out() << "fei::FEI_Impl::loadComplete: loadComplete can not be called"
	 << " until after initComplete has been called."<<FEI_ENDL;
    return(-1);
  }

  if (solveType_ == FEI_AGGREGATE_SUM) {
    for(unsigned i=0; i<A_.size(); ++i) {
	CHK_ERR( A_[i]->gatherFromOverlap() );
    }

    for(unsigned j=0; j<b_.size(); ++j) {
      CHK_ERR( b_[j]->gatherFromOverlap() );
    }
  
    CHK_ERR( aggregateSystem() );
  }

  CHK_ERR( linSys_->loadComplete(applyBCs, globalAssemble) );

  if (b_.size() > 1) {
    int rhs_counter = 0;
    for(unsigned i=0; i<b_.size(); ++i) {
      FEI_OSTRINGSTREAM osstr;
      osstr << "rhs_" << rhs_counter++;
      CHK_ERR( linSys_->putAttribute(osstr.str().c_str(), b_[i].get()) );
    }
  }

  newData_ = false;

  return(0);
}

int fei::FEI_Impl::aggregateSystem()
{
  if (wrapper_[0].get() == NULL) {
    ERReturn(-1);
  }

  if (wrapper_[0].get() != NULL) {
    CHK_ERR( aggregateSystem_LinSysCore() );
  }

  return(0);
}

int fei::FEI_Impl::aggregateSystem_LinSysCore()
{
  fei::SharedPtr<LinearSystemCore> lsc = wrapper_[0]->getLinearSystemCore();
  if (lsc.get() == NULL) return(-1);

  if (soln_fei_matrix_ == NULL) {
    soln_fei_matrix_ = new Data();

    CHK_ERR( lsc->copyOutMatrix(1.0, *soln_fei_matrix_) );
  }

  if (soln_fei_vector_ == NULL) {
    soln_fei_vector_ = new Data();

    CHK_ERR( lsc->setRHSID(rhsIDs_[0]) );
    CHK_ERR( lsc->copyOutRHSVector(1.0, *soln_fei_vector_) );
  }

  Data tmp, tmpv;
  for(unsigned i=0; i<matrixIDs_.size(); ++i) {
    fei::SharedPtr<LinearSystemCore> lsci = wrapper_[i]->getLinearSystemCore();
    if (lsci.get() == NULL) return(-1);

    if (i==0) {
      CHK_ERR( lsci->copyInMatrix(matScalars_[i], *soln_fei_matrix_) );
    }
    else {
      CHK_ERR( lsci->getMatrixPtr(tmp) );
      CHK_ERR( lsc->sumInMatrix(matScalars_[i], tmp) );
    }
  }

  int last_mat = matrixIDs_.size() - 1;
  for(unsigned i=0; i<rhsIDs_.size(); ++i) {
    fei::SharedPtr<LinearSystemCore> lsci = (int)i<last_mat ?
      wrapper_[i]->getLinearSystemCore() : wrapper_[last_mat]->getLinearSystemCore();

    if (i==0) {
      CHK_ERR( lsci->copyInRHSVector(rhsScalars_[i], *soln_fei_vector_) );
    }
    else {
      CHK_ERR( lsci->setRHSID(rhsIDs_[i]) );
      CHK_ERR( lsci->getRHSVectorPtr(tmpv) );
      CHK_ERR( lsc->sumInRHSVector(rhsScalars_[i], tmpv) );
    }
  }

  return(0);
}

int fei::FEI_Impl::residualNorm(int whichNorm,
				 int numFields,
				 int* fieldIDs,
				 double* norms)
{
  CHK_ERR( loadComplete() );

  fei::SharedPtr<fei::VectorSpace> rowspace =
    matGraph_->getRowSpace();
  int numLocalEqns = rowspace->getNumIndices_Owned();

  std::vector<double> residValues(numLocalEqns);
  double* residValuesPtr = &residValues[0];

  std::vector<int> globalEqnOffsets;
  rowspace->getGlobalIndexOffsets(globalEqnOffsets);
  int firstLocalOffset = globalEqnOffsets[localProc_];

  if (wrapper_[0].get() == NULL) {
    fei::SharedPtr<fei::Vector> r = factory_[0]->createVector(rowspace);

    //form r = A*x
    CHK_ERR( A_[index_current_]->multiply(x_.get(), r.get()) );

    //form r = b - A*x
    CHK_ERR( r->update(1.0, b_[index_current_rhs_row_].get(), -1.0) );

    //now put the values from r into the residValues array.
    for(int ii=0; ii<numLocalEqns; ++ii) {
      int index = firstLocalOffset+ii;
      CHK_ERR( r->copyOut(1, &index, &(residValuesPtr[ii]) ) );
    }
  }
  else {
    fei::SharedPtr<LinearSystemCore> linSysCore = wrapper_[0]->getLinearSystemCore();
    if (linSysCore.get() != NULL) {
      CHK_ERR( linSysCore->formResidual(residValuesPtr, numLocalEqns) );
    }
    else {
      fei::console_out() << "FEI_Impl::residualNorm: warning: residualNorm not implemented"
	   << " for FiniteElementData."<<FEI_ENDL;
      int offset = 0;
      for(int ii=0; ii<numFields; ++ii) {
	int fieldSize = rowspace->getFieldSize(fieldIDs[ii]);
	for(int jj=0; jj<fieldSize; ++jj) {
	  norms[offset++] = -99.9;
	}
      }
      return(0);
    }
  }

  int numLocalNodes = rowspace->getNumOwnedIDs(nodeIDType_);
  std::vector<int> nodeIDs(numLocalNodes);
  int* nodeIDsPtr = &nodeIDs[0];
  int check;
  CHK_ERR( rowspace->getOwnedIDs(nodeIDType_, numLocalNodes,
					nodeIDsPtr, check) );
  std::vector<int> indices(numLocalEqns);
  int* indicesPtr = &indices[0];

  std::vector<double> tmpNorms(numFields, 0.0);

  for(int i=0; i<numFields; ++i) {
    int fieldSize = rowspace->getFieldSize(fieldIDs[i]);

    CHK_ERR( rowspace->getGlobalIndices(numLocalNodes, nodeIDsPtr,
					nodeIDType_, fieldIDs[i],
					indicesPtr) );

    double tmp = 0.0;
    for(int j=0; j<fieldSize*numLocalNodes; ++j) {
      if (indicesPtr[j] < 0) {
	continue;
      }

      double val = std::abs(residValuesPtr[indicesPtr[j]-firstLocalOffset]);
      switch(whichNorm) {
      case 0:
	if (val > tmp) tmp = val;
	break;
      case 1:
	tmp += val;
	break;
      case 2:
	tmp += val*val;
	break;
      default:
	fei::console_out() << "FEI_Impl::residualNorm: whichNorm="<<whichNorm<<" not recognized"
	     << FEI_ENDL;
	return(-1);
      }
    }
    tmpNorms[i] = tmp;
  }

  std::vector<double> normsArray(numFields);
  for(int i=0; i<numFields; ++i) {
    normsArray[i] = norms[i];
  }

  switch(whichNorm) {
  case 0:
    CHK_ERR( fei::GlobalMax(comm_, tmpNorms, normsArray) );
    break;
  default:
    CHK_ERR( fei::GlobalSum(comm_, tmpNorms, normsArray) );
  }

  for(int i=0; i<numFields; ++i) {
    norms[i] = normsArray[i];
  }

  if (whichNorm == 2) {
    for(int ii=0; ii<numFields; ++ii) {
      norms[ii] = std::sqrt(norms[ii]);
    }
  }

  return(0);
}

int fei::FEI_Impl::solve(int& status)
{
  CHK_ERR( loadComplete() );

  fei::SharedPtr<fei::Solver> solver = factory_[0]->createSolver();

  std::vector<std::string> stdstrings;
  fei::utils::char_ptrs_to_strings(numParams_, paramStrings_, stdstrings);
  fei::ParameterSet params;
  fei::utils::parse_strings(stdstrings, " ", params);

  int err = solver->solve(linSys_.get(),
			  NULL,//preconditioningMatrix
			  params, iterations_, status);

  CHK_ERR( x_->scatterToOverlap() );

  return(err);
}

int fei::FEI_Impl::iterations(int& itersTaken) const
{
  itersTaken = iterations_;

  return(0);
}

int fei::FEI_Impl::version(const char*& versionString)
{
  versionString = fei::utils::version();

  return(0);
}

int fei::FEI_Impl::cumulative_cpu_times(double& initTime,
					 double& loadTime,
					 double& solveTime,
					 double& solnReturnTime)
{
  initTime = initTime_;
  loadTime = loadTime_;
  solveTime = solveTime_;
  solnReturnTime = solnReturnTime_;

  return(0);
}

int fei::FEI_Impl::getBlockNodeSolution(GlobalID elemBlockID,  
					 int numNodes, 
					 const GlobalID *nodeIDs, 
					 int *offsets,
					 double *results)
{
  (void)elemBlockID;
  return ( getNodalSolution(numNodes, nodeIDs, offsets, results) );
}

int fei::FEI_Impl::getNodalSolution(int numNodes,
				     const GlobalID* nodeIDs,
				     int* offsets,
				     double* results)
{
  int j;
  int offset = 0;
  std::vector<int> fieldIDs;
  for(int i=0; i<numNodes; ++i) {
    offsets[i] = offset;

    GlobalID nodeID = nodeIDs[i];

    int numFields = rowSpace_->getNumFields(nodeIDType_, nodeID);
    iwork_.resize( numFields*2 );
    int* fieldSizes = &iwork_[0];
    rowSpace_->getFields(nodeIDType_, nodeID, fieldIDs);

    int numDOF = 0;
    for(j=0; j<numFields; ++j) {
      fieldSizes[j] = rowSpace_->getFieldSize(fieldIDs[j]);
      numDOF += fieldSizes[j];
    }

    if (!rowSpace_->isLocal(nodeIDType_, nodeID)) {
      offset += numDOF;
      continue;
    }

    int results_offset = offset;
    for(j=0; j<numFields; ++j) {
      CHK_ERR( x_->copyOutFieldData(fieldIDs[j], nodeIDType_,
				    1, &nodeID,
				    &(results[results_offset])));

      results_offset += fieldSizes[j];
    }

    offset += numDOF;
  }

  offsets[numNodes] = offset;

  return(0);
}

int fei::FEI_Impl::getBlockFieldNodeSolution(GlobalID elemBlockID,
                                  int fieldID,
                                  int numNodes, 
                                  const GlobalID *nodeIDs, 
                                  double *results)
{
  throw std::runtime_error("FEI_Impl::getBlockFieldNodeSolution not implemented.");
  return(0);
}

int fei::FEI_Impl::getBlockElemSolution(GlobalID elemBlockID,  
					 int numElems, 
					 const GlobalID *elemIDs,
					 int& numElemDOFPerElement,
					 double *results)
{
  std::map<int,int>::const_iterator b_iter = block_dof_per_elem_.find(elemBlockID);
  if (b_iter == block_dof_per_elem_.end()) return(-1);
  numElemDOFPerElement = (*b_iter).second;

  fei::ConnectivityBlock* block = matGraph_->getConnectivityBlock(elemBlockID);
  if (block==NULL) {
    FEI_OSTRINGSTREAM osstr;
    osstr<< "fei::FEI_Impl::getBlockElemSolution ERROR, elemBlockID "
	 << elemBlockID << " not valid.";
    throw std::runtime_error(osstr.str());
  }

  fei::Pattern* pattern = block->getRowPattern();

  int numIDs = pattern->getNumIDs();
  const int* idTypes = pattern->getIDTypes();
  const int* fIDs= pattern->getFieldIDs();
  const int* numFieldsPerID = pattern->getNumFieldsPerID();

  std::vector<int> fieldIDs;
  int foffset = 0;
  int i;
  for(i=0; i<numIDs; ++i) {
    if (idTypes[i] == elemIDType_) {
      for(int j=0; j<numFieldsPerID[i]; ++j) {
	fieldIDs.push_back(fIDs[foffset++]);
      }
      break;
    }
    foffset += numFieldsPerID[i];
  }

  if (fieldIDs.size() < 1) {
    ERReturn(-1);
  }

  int offset = 0;
  for(i=0; i<numElems; ++i) {
    foffset = offset;
    for(size_t j=0; j<fieldIDs.size(); ++j) {
      int fieldSize;
      getFieldSize(fieldIDs[j], fieldSize);

      CHK_ERR( x_->copyOutFieldData(fieldIDs[j], elemIDType_, 1, &(elemIDs[i]),
				    &(results[foffset])) );
      foffset += fieldSize;
    }

    offset += numElemDOFPerElement;
  }

  return(0);
}

int fei::FEI_Impl::getNumCRMultipliers(int& numMultCRs)
{
  numMultCRs = rowSpace_->getNumOwnedAndSharedIDs(constraintIDType_);
  return(0);
}

int fei::FEI_Impl::getCRMultIDList(int numMultCRs, int* multIDs)
{
  int checkNum;
  CHK_ERR( rowSpace_->getOwnedAndSharedIDs(constraintIDType_, numMultCRs,
				  multIDs, checkNum) );
  return(0);
}

int fei::FEI_Impl::getCRMultipliers(int numCRs,
				     const int* CRIDs,
				     double *multipliers)
{
  iwork_.resize(numCRs);

  for(int i=0; i<numCRs; ++i) {
    CHK_ERR( rowSpace_->getGlobalIndex(constraintIDType_, CRIDs[i], iwork_[i]));
  }

  CHK_ERR( x_->copyOut(numCRs, &iwork_[0], multipliers) );

  return(0);
}

int fei::FEI_Impl::putBlockNodeSolution(GlobalID elemBlockID, 
                             int numNodes, 
                             const GlobalID *nodeIDs, 
                             const int *offsets,
                             const double *estimates)
{
  throw std::runtime_error("FEI_Impl::putBlockNodeSolution not implemented.");
  return(0);
}

int fei::FEI_Impl::putBlockFieldNodeSolution(GlobalID elemBlockID, 
                                  int fieldID, 
                                  int numNodes, 
                                  const GlobalID *nodeIDs, 
                                  const double *estimates)
{
  throw std::runtime_error("FEI_Impl::putBlockFieldNodeSolution not implemented.");
  return(0);
}

int fei::FEI_Impl::putBlockElemSolution(GlobalID elemBlockID,  
                             int numElems, 
                             const GlobalID *elemIDs, 
                             int dofPerElem,
                             const double *estimates)
{
  throw std::runtime_error("FEI_Impl::putBlockElemSolution not implemented.");
  return(0);
}

int fei::FEI_Impl::putCRMultipliers(int numMultCRs, 
                         const int* CRIDs,
                         const double* multEstimates)
{
  throw std::runtime_error("FEI_Impl::putCRMultipliers not implemented.");
  return(0);
}

int fei::FEI_Impl::getBlockNodeIDList(GlobalID elemBlockID,
                           int numNodes,
                           GlobalID *nodeIDs)
{
  if (!nodeset_filled_ || elemBlockID != nodeset_blockid_) {
    CHK_ERR( fillNodeset(elemBlockID) );
  }

  fei::copySetToArray(nodeset_, numNodes, nodeIDs);
  return(0);
}

int fei::FEI_Impl::fillNodeset(int blockID) const
{
  if (nodeset_filled_ && blockID == nodeset_blockid_) {
    return(0);
  }

  fei::ConnectivityBlock* block = matGraph_->getConnectivityBlock(blockID);
  if (block==NULL) {
    FEI_OSTRINGSTREAM osstr;
    osstr<< "fei::FEI_Impl::fillNodeset ERROR, blockID "
	 << blockID << " not valid.";
    throw std::runtime_error(osstr.str());
  }

  int nodeType = 0;
  snl_fei::RecordCollection* nodeRecords = NULL;
  matGraph_->getRowSpace()->getRecordCollection(nodeType, nodeRecords);
  std::vector<int>& nodes = block->getRowConnectivities();

  nodeset_.clear();

  for(unsigned i=0; i<nodes.size(); ++i) {
    nodeset_.insert(nodeRecords->getRecordWithLocalID(nodes[i])->getID());
  }

  nodeset_filled_ = true;
  nodeset_blockid_ = blockID;

  return(0);
}

int fei::FEI_Impl::getBlockElemIDList(GlobalID elemBlockID, 
                          int numElems, 
                          GlobalID* elemIDs)
{
  fei::ConnectivityBlock* block = matGraph_->getConnectivityBlock(elemBlockID);
  if (block==NULL) {
    FEI_OSTRINGSTREAM osstr;
    osstr<< "fei::FEI_Impl::getBlockElemIDList ERROR, elemBlockID "
	 << elemBlockID << " not valid.";
    throw std::runtime_error(osstr.str());
  }

  std::map<int,int>& elemIDSet = block->getConnectivityIDs();

  fei::copyKeysToArray(elemIDSet, numElems, elemIDs);

  return(0);
}

int fei::FEI_Impl::getNumSolnParams(GlobalID nodeID, int& numSolnParams) const
{
  numSolnParams = rowSpace_->getNumDegreesOfFreedom( nodeIDType_, nodeID);
  return(0);
}

int fei::FEI_Impl::getNumElemBlocks(int& numElemBlocks) const
{
  numElemBlocks = matGraph_->getNumConnectivityBlocks();
  return(0);
}

int fei::FEI_Impl::getNumBlockActNodes(GlobalID blockID, int& numNodes) const
{
  if (!nodeset_filled_ || blockID != nodeset_blockid_) {
    CHK_ERR( fillNodeset(blockID) );
  }

  numNodes = nodeset_.size();
  return(0);
}

int fei::FEI_Impl::getNumBlockActEqns(GlobalID blockID, int& numEqns) const
{
  throw std::runtime_error("fei::FEI_Impl::getNumBlockActEqns not implemented.");
  return(0);
}

int fei::FEI_Impl::getNumNodesPerElement(GlobalID blockID, int& nodesPerElem) const
{
  nodesPerElem = matGraph_->getNumIDsPerConnectivityList(blockID);
  return(0);
}

int fei::FEI_Impl::getNumEqnsPerElement(GlobalID blockID, int& numEqns) const
{
  numEqns = matGraph_->getConnectivityNumIndices(blockID);
  return(0);
}

int fei::FEI_Impl::getNumBlockElements(GlobalID blockID, int& numElems) const
{
  const fei::ConnectivityBlock* cblock =
    matGraph_->getConnectivityBlock(blockID);
  if (cblock != NULL) {
    numElems = cblock->getConnectivityIDs().size();
  }
  else {
    numElems = 0;
  }
  return(0);
}

int fei::FEI_Impl::getNumBlockElemDOF(GlobalID blockID, int& DOFPerElem) const
{
  std::map<int,int>::const_iterator b_iter = block_dof_per_elem_.find(blockID);
  if (b_iter == block_dof_per_elem_.end()) DOFPerElem = 0;
  else DOFPerElem = (*b_iter).second;

  return(0);
}

int fei::FEI_Impl::getParameters(int& numParams, char**& paramStrings)
{
  numParams = numParams_;
  paramStrings = paramStrings_;
  return(0);
}

int fei::FEI_Impl::getFieldSize(int fieldID, int& numScalars)
{
  numScalars = rowSpace_->getFieldSize(fieldID);
  return(0);
}

int fei::FEI_Impl::getEqnNumbers(GlobalID ID,
				  int idType, 
				  int fieldID,
				  int& numEqns,
				  int* eqnNumbers)
{
  numEqns = rowSpace_->getFieldSize(fieldID);
  CHK_ERR( rowSpace_->getGlobalIndices(1, &ID, idType, fieldID, eqnNumbers) );
  return(0);
}

int fei::FEI_Impl::getNodalFieldSolution(int fieldID,
					  int numNodes,
					  const GlobalID* nodeIDs,
					  double* results)
{
  CHK_ERR( x_->copyOutFieldData(fieldID, nodeIDType_, numNodes,
				nodeIDs, results) );
  return(0);
}

int fei::FEI_Impl::getNumLocalNodes(int& numNodes)
{
  numNodes = rowSpace_->getNumOwnedAndSharedIDs(nodeIDType_);
  return(0);
}

int fei::FEI_Impl::getLocalNodeIDList(int& numNodes,
				       GlobalID* nodeIDs,
				       int lenNodeIDs)
{
  CHK_ERR( rowSpace_->getOwnedAndSharedIDs(nodeIDType_, lenNodeIDs, nodeIDs, numNodes) );
  return(0);
}

int fei::FEI_Impl::putNodalFieldData(int fieldID,
				      int numNodes,
				      const GlobalID* nodeIDs,
				      const double* nodeData)
{
  int err = 0;
  if (fieldID < 0) {
    bool data_passed = false;
    if (wrapper_[0].get() != NULL) {
      std::vector<int> numbers(numNodes);
      for(int i=0; i<numNodes; ++i) {
        err = rowSpace_->getGlobalBlkIndex(nodeIDType_, nodeIDs[i], numbers[i]);
        if (err != 0) {
          fei::console_out() << "fei::FEI_Impl::putNodalFieldData ERROR, nodeID "
            << nodeIDs[i] << " not found."<<FEI_ENDL;
          ERReturn(-1);
        }
      }

      int fieldSize = 0;
      try {
        fieldSize = rowSpace_->getFieldSize(fieldID);
      }
      catch (std::runtime_error& exc) {
        fei::console_out() << "fei::FEI_Impl::putNodalFieldData ERROR: " <<exc.what()<<FEI_ENDL;
        ERReturn(-1);
      }

      fei::SharedPtr<LinearSystemCore> linSysCore = wrapper_[0]->getLinearSystemCore();
      if (linSysCore.get() != NULL) {
        linSysCore->putNodalFieldData(fieldID, fieldSize, 
            &numbers[0], numNodes, nodeData);
        data_passed = true;
      }
      else {
        //If we enter this block, we're probably dealing with a FiniteElementData
        //instance.
        fei::SharedPtr<FiniteElementData> fedata = wrapper_[0]->getFiniteElementData();
        if (fedata.get() != NULL) {
          fedata->putNodalFieldData(fieldID, fieldSize, numNodes,
              &numbers[0], nodeData);
          data_passed = true;
        }
      }
    }

    if (!data_passed) {
      //If we get to here and data_passed is false, wrapper_[0] is probably NULL.
      if (wrapper_[0].get() == NULL) {
        fei::SharedPtr<fei::Vector> dataVector =factory_[0]->createVector(matGraph_);

        CHK_ERR( dataVector->copyInFieldData(fieldID, nodeIDType_,
              numNodes, nodeIDs, nodeData) );
        if (fieldID == -3) {
          CHK_ERR( linSys_->putAttribute("coordinates", dataVector.get()) );
        }
        else {
          FEI_OSTRINGSTREAM osstr;
          osstr << "fieldID:" << fieldID;
          CHK_ERR( linSys_->putAttribute(osstr.str().c_str(), dataVector.get()) );
        }
      }
      else {
        fei::console_out() << "fei::FEI_Impl::putNodalFieldData ERROR, non-null LibraryWrapper"
          << " contains neither LinearSystemCore or FiniteElementData. " <<FEI_ENDL;
        ERReturn(-1);
      }
    }
    return(0);
  }

  CHK_ERR( x_->copyInFieldData(fieldID, nodeIDType_,
        numNodes, nodeIDs, nodeData));
  return(0);
}
