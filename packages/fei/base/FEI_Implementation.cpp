/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_CommUtils.hpp>
#include <fei_iostream.hpp>
#include <fei_fstream.hpp>
#include <fei_sstream.hpp>

#include <fei_utils.hpp>

#include <fei_Data.hpp>
#include <fei_LinearSystemCore.hpp>
#include <fei_FiniteElementData.hpp>

#include <fei_LibraryWrapper.hpp>
#include <fei_Filter.hpp>

#include <fei_LinSysCoreFilter.hpp>
#include <fei_FEDataFilter.hpp>

#include <SNL_FEI_Structure.hpp>
#include <fei_BlockDescriptor.hpp>
#include <fei_NodeDatabase.hpp>
#include <fei_ConnectivityTable.hpp>
#include <FEI_Implementation.hpp>
#include <snl_fei_Utils.hpp>

#undef fei_file
#define fei_file "FEI_Implementation.cpp"

#include <fei_ErrMacros.hpp>

//------------------------------------------------------------------------------
FEI_Implementation::FEI_Implementation(fei::SharedPtr<LibraryWrapper> libWrapper,
                                       MPI_Comm comm, int masterRank)
 : wrapper_(libWrapper),
   linSysCore_(NULL),
   lscArray_(),
   haveLinSysCore_(false),
   haveFEData_(false),
   problemStructure_(NULL),
   filter_(NULL),
   numInternalFEIs_(0),
   internalFEIsAllocated_(false),
   matrixIDs_(),
   numRHSIDs_(),
   rhsIDs_(),
   IDsAllocated_(false),
   matScalars_(),
   matScalarsSet_(false),
   rhsScalars_(),
   rhsScalarsSet_(false),
   index_soln_filter_(0),
   index_current_filter_(0),
   index_current_rhs_row_(0),
   solveType_(-1),
   setSolveTypeCalled_(false),
   initPhaseIsComplete_(false),
   aggregateSystemFormed_(false),
   newMatrixDataLoaded_(0),
   soln_fei_matrix_(NULL),
   soln_fei_vector_(NULL),
   comm_(comm),
   masterRank_(0),
   localRank_(0),
   numProcs_(1),
   outputLevel_(0),
   solveCounter_(1),
   debugOutput_(0),
   dbgOStreamPtr_(NULL),
   dbgFileOpened_(false),
   dbgFStreamPtr_(NULL),
   initTime_(0.0),
   loadTime_(0.0),
   solveTime_(0.0),
   solnReturnTime_(0.0),
   numParams_(0),
   paramStrings_(NULL)
{
  localRank_ = fei::localProc(comm);
  numProcs_ = fei::numProcs(comm);

   problemStructure_ = new SNL_FEI_Structure(comm_);

   if (problemStructure_ == NULL) {
     messageAbort("problemStructure allocation failed");
   }

   //If we have a FiniteElementData instance as the underlying
   //data receptacle and solver, then we'll set the shared-node-ownership rule
   //to make sure shared nodes are owned by a proc which contains them in
   //local elements.
   haveFEData_ = wrapper_->haveFiniteElementData();
   if (haveFEData_) {
     haveFEData_ = true;
     NodeCommMgr& nodeCommMgr = problemStructure_->getNodeCommMgr();
     nodeCommMgr.setSharedOwnershipRule(NodeCommMgr::PROC_WITH_LOCAL_ELEM);
   }

   haveLinSysCore_ = wrapper_->haveLinearSystemCore();
   if (haveLinSysCore_) {
     linSysCore_ = wrapper_->getLinearSystemCore();
     lscArray_.push_back(linSysCore_);
   }

   numInternalFEIs_ = 1;
   matrixIDs_.resize(1);
   matrixIDs_[0] = 0;
   numRHSIDs_.resize(1);
   numRHSIDs_[0] = 1;
   rhsIDs_.resize(1);
   rhsIDs_[0] = new int[1];
   rhsIDs_[0][0] = 0;
   rhsScalars_.resize(numInternalFEIs_);
   for(int ii=0; ii<numInternalFEIs_; ii++) rhsScalars_[ii] = NULL;

   //  and the time spent in the constructor is...

   return;
}

//------------------------------------------------------------------------------
FEI_Implementation::~FEI_Implementation()
{
  //
  //  Destructor function. Free allocated memory, etc.
  //

  if (debugOutput_) {
    (*dbgOStreamPtr_) << "FEI: destructor" << FEI_ENDL;
  }

  if (soln_fei_matrix_) {
      linSysCore_->destroyMatrixData(*soln_fei_matrix_);
      delete soln_fei_matrix_;
      soln_fei_matrix_ = NULL;
   }

   if (soln_fei_vector_) {
      linSysCore_->destroyVectorData(*soln_fei_vector_);
      delete soln_fei_vector_;
      soln_fei_vector_ = NULL;
   }

   deleteIDs();

   if (internalFEIsAllocated_) {
      for(int j=0; j<numInternalFEIs_; j++){
         delete filter_[j];
      }
      delete [] filter_;
   }

   deleteRHSScalars();

   internalFEIsAllocated_ = false;
   numInternalFEIs_ = 0;

   delete problemStructure_;

   for(int k=0; k<numParams_; k++) delete [] paramStrings_[k];
   delete [] paramStrings_;

   if (dbgFileOpened_ == true) {
     dbgFStreamPtr_->close(); delete dbgFStreamPtr_;
   }
   else delete dbgOStreamPtr_;
}


//------------------------------------------------------------------------------
void FEI_Implementation::deleteIDs()
{
  matrixIDs_.resize(0);
  for(size_t i=0; i<rhsIDs_.size(); i++) {
    delete [] rhsIDs_[i];
  }
  rhsIDs_.resize(0);
  numRHSIDs_.resize(0);
}

//------------------------------------------------------------------------------
void FEI_Implementation::deleteRHSScalars()
{
  for(size_t i=0; i<rhsScalars_.size(); i++) {
    delete [] rhsScalars_[i];
  }
  rhsScalars_.resize(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::setCurrentMatrix(int matID)
{
   if (debugOutput_) {
     (*dbgOStreamPtr_) << "FEI: setCurrentMatrix" << FEI_ENDL << "#matrix-id"
                       << FEI_ENDL<<matID<<FEI_ENDL;
   }

   index_current_filter_ = -1;

   for(int i=0; i<numInternalFEIs_; i++){
      if (matrixIDs_[i] == matID) index_current_filter_ = i;
   }

   if (debugOutput_) {
     (*dbgOStreamPtr_) << "#--- ID: " << matID
                       << ", ind: "<<index_current_filter_<<FEI_ENDL;
   }

   //if matID wasn't found, return non-zero (error)
   if (index_current_filter_ == -1) {
      fei::console_out() << "FEI_Implementation::setCurrentMatrix: ERROR, invalid matrix ID "
           << "supplied" << FEI_ENDL;
      return(-1);
   }

   debugOut("#FEI_Implementation leaving setCurrentMatrix");

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::getParameters(int& numParams, char**& paramStrings)
{
  numParams = numParams_;
  paramStrings = paramStrings_;
  return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::setCurrentRHS(int rhsID)
{
  if (debugOutput_) {
    (*dbgOStreamPtr_) << "FEI: setCurrentRHS" << FEI_ENDL << "#rhs-id"
                      << FEI_ENDL<<rhsID<<FEI_ENDL;
  }

  bool found = false;

  for(int j=0; j<numInternalFEIs_; j++){
    int index = fei::searchList(rhsID, rhsIDs_[j], numRHSIDs_[j]);
    if (index >= 0) {
      index_current_rhs_row_ = j;
      CHK_ERR( filter_[index_current_rhs_row_]->setCurrentRHS(rhsID) )
        found = true;
      break;
    }
  }

  if (!found) {
    fei::console_out() << "FEI_Implementation::setCurrentRHS: ERROR, invalid RHS ID" 
         << FEI_ENDL;
    ERReturn(-1);
  }

  debugOut("#FEI_Implementation leaving setCurrentRHS");

  return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::setSolveType(int solveType)
{
  if (debugOutput_) {
    (*dbgOStreamPtr_)<<"FEI: setSolveType"<<FEI_ENDL;
    (*dbgOStreamPtr_)<<solveType<<FEI_ENDL;
  }

  solveType_ = solveType;

  if (solveType_ == FEI_SINGLE_SYSTEM) {
    if (matrixIDs_.size() > 1) {
      messageAbort("setSolveType: solve-type is FEI_SINGLE_SYSTEM, but setIDLists() has been called with numMatrices > 1.");
    }
  }
  else if (solveType_ == FEI_EIGEN_SOLVE) {
  }
  else if (solveType_ == FEI_AGGREGATE_SUM) {
    //solving a linear-combination of separately
    //assembled matrices and rhs vectors
  }
  else if (solveType_ == FEI_AGGREGATE_PRODUCT) {
    //solving a product of separately assembled
    //matrices -- i.e., (C^T*M*C)x = rhs
  }
  else if (solveType_ == 4) {
    //4 means we'll be doing a multi-level solution
  }

  return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::setIDLists(int numMatrices, const int* matrixIDs,
                                   int numRHSs, const int* rhsIDs)
{
   if (debugOutput_) {
     FEI_OSTREAM& os = *dbgOStreamPtr_;
     os << "FEI: setIDLists" << FEI_ENDL
        << "#num-matrices" << FEI_ENDL << numMatrices << FEI_ENDL
        << "#matrixIDs" << FEI_ENDL;
     int i;
     for(i=0; i<numMatrices; ++i) os << matrixIDs[i] << " ";
     os << FEI_ENDL << "#num-rhs's" << FEI_ENDL;
     for(i=0; i<numRHSs; ++i) os << rhsIDs[i] << " ";
     os << FEI_ENDL;
   }

   deleteIDs();

   // We will try to assign the rhs's evenly over the matrices. i.e., give
   // roughly equal numbers of rhs's to each matrix.

   //first, let's make sure we have at least 1 matrixID to which we can assign
   //rhs's...
   int myNumMatrices = numMatrices;
   if (myNumMatrices == 0) myNumMatrices = 1;

   matrixIDs_.resize(myNumMatrices);

   if (rhsScalars_.size() != 0) deleteRHSScalars();

   numInternalFEIs_ = myNumMatrices;

   if (numMatrices == 0) {
      matrixIDs_[0] = 0;
   }
   else {
      for(int i=0; i<numMatrices; i++) matrixIDs_[i] = matrixIDs[i];
   }

   int quotient = numRHSs/myNumMatrices;
   int rem = numRHSs%numMatrices;

   //the allocateInternalFEIs function which will be called later from within
   //initComplete(), takes a list of matrixIDs, and a list
   //of numRHSsPerMatrix, and then a table of rhsIDs, where the table has a row
   //for each matrixID. Each of those rows is a list of the rhsIDs assigned to
   //the corresponding matrix. Is that clear???

   numRHSIDs_.resize(myNumMatrices);
   rhsIDs_.resize(myNumMatrices);

   int offset = 0;
   for(int i=0; i<myNumMatrices; i++) {
      numRHSIDs_[i] = quotient;
      if (i < rem) numRHSIDs_[i]++;

      rhsIDs_[i] = numRHSIDs_[i] > 0 ? new int[numRHSIDs_[i]] : NULL ;

      for(int j=0; j<numRHSIDs_[i]; j++) {
        rhsIDs_[i][j] = rhsIDs[offset+j];
      }

      offset += numRHSIDs_[i];
   }

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::initFields(int numFields,
                                   const int *fieldSizes,
                                   const int *fieldIDs,
                                   const int *fieldTypes)
{
    CHK_ERR( problemStructure_->initFields(numFields, fieldSizes, fieldIDs, fieldTypes) );

    return(0);
}
 
//------------------------------------------------------------------------------
int FEI_Implementation::initElemBlock(GlobalID elemBlockID,
                                      int numElements,
                                      int numNodesPerElement,
                                      const int* numFieldsPerNode,
                                      const int* const* nodalFieldIDs,
                                      int numElemDofFieldsPerElement,
                                      const int* elemDOFFieldIDs,
                                      int interleaveStrategy)
{
   CHK_ERR( problemStructure_->initElemBlock(elemBlockID,
                                             numElements,
                                             numNodesPerElement,
                                             numFieldsPerNode,
                                             nodalFieldIDs,
                                             numElemDofFieldsPerElement,
                                             elemDOFFieldIDs,
                                             interleaveStrategy) );

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::initElem(GlobalID elemBlockID,
                                 GlobalID elemID,
                                 const GlobalID* elemConn)
{
   CHK_ERR( problemStructure_->initElem(elemBlockID, elemID, elemConn) );

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::initSlaveVariable(GlobalID slaveNodeID, 
                                          int slaveFieldID,
                                          int offsetIntoSlaveField,
                                          int numMasterNodes,
                                          const GlobalID* masterNodeIDs,
                                          const int* masterFieldIDs,
                                          const double* weights,
                                          double rhsValue)
{
   CHK_ERR( problemStructure_->initSlaveVariable(slaveNodeID, slaveFieldID,
                                            offsetIntoSlaveField,
                                                   numMasterNodes, masterNodeIDs,
                                            masterFieldIDs, weights, rhsValue));

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::deleteMultCRs()
{
  debugOut("FEI: deleteMultCRs");

  CHK_ERR( problemStructure_->deleteMultCRs() );

  int err = -1;
  if (internalFEIsAllocated_) {
    err = filter_[index_current_filter_]->deleteMultCRs();
  }

  return(err);
}

//------------------------------------------------------------------------------
int FEI_Implementation::initSharedNodes(int numSharedNodes,
                                        const GlobalID *sharedNodeIDs,
                                        const int* numProcsPerNode,
                                        const int *const *sharingProcIDs)
{
  //
  //  In this function we simply accumulate the incoming data into
  // internal arrays in the problemStructure_ object.
  //
  CHK_ERR( problemStructure_->initSharedNodes(numSharedNodes,
                                              sharedNodeIDs,
                                              numProcsPerNode,
                                              sharingProcIDs));

  return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::initCRMult(int numCRNodes,
                                   const GlobalID* CRNodes,
                                   const int *CRFields,
                                   int& CRID)
{
//
// Store Lagrange Multiplier constraint data into internal structures, and
// return an identifier 'CRID' by which this constraint may be referred to
// later.
//

   CHK_ERR( problemStructure_->initCRMult(numCRNodes,
                                          CRNodes,
                                          CRFields,
                                          CRID));

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::initCRPen(int numCRNodes,
                                  const GlobalID* CRNodes,
                                  const int *CRFields,
                                  int& CRID)
{
//
// Store penalty constraint data and return an identifier 'CRID' by which the
// constraint may be referred to later.
//

   CHK_ERR( problemStructure_->initCRPen(numCRNodes,
                                         CRNodes,
                                         CRFields,
                                         CRID));

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::initComplete()
{
    bool generateGraph = !haveFEData_;

    CHK_ERR( problemStructure_->initComplete(generateGraph) );

    //now allocate one or more internal instances of Filter, depending on
    //whether the user has indicated that they're doing an aggregate solve
    //etc., via the functions setSolveType() and setIDLists().

    CHK_ERR( allocateInternalFEIs() );

    for(int i=0; i<numInternalFEIs_; ++i) {
      CHK_ERR( filter_[i]->initialize() );
    }

    problemStructure_->destroyMatIndices();

    initPhaseIsComplete_ = true;
    return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::resetSystem(double s)
{
//
//  This puts the value s throughout both the matrix and the vector.
//
   if (!internalFEIsAllocated_) return(0);

   CHK_ERR( filter_[index_current_filter_]->resetSystem(s))
 
   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::resetMatrix(double s)
{
   if (!internalFEIsAllocated_) return(0);

   CHK_ERR( filter_[index_current_filter_]->resetMatrix(s))

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::resetRHSVector(double s)
{
   if (!internalFEIsAllocated_) return(0);

   CHK_ERR( filter_[index_current_filter_]->resetRHSVector(s))

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::resetInitialGuess(double s)
{
   if (!internalFEIsAllocated_) return(0);

   CHK_ERR( filter_[index_current_filter_]->resetInitialGuess(s))

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::loadNodeBCs(int numNodes,
                                    const GlobalID *nodeIDs,
                                    int fieldID,
                                    const int* offsetsIntoField,
                                    const double* prescribedValues)
{
   if (!internalFEIsAllocated_)
      notAllocatedAbort("FEI_Implementation::loadNodeBCs");

   int index = index_current_filter_;
   if (solveType_ == 2) index = index_soln_filter_;

   CHK_ERR( filter_[index]->loadNodeBCs(numNodes,
                               nodeIDs, fieldID,
                               offsetsIntoField, prescribedValues));

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::loadElemBCs(int numElems,
                                    const GlobalID *elemIDs,
                                    int fieldID,
                                    const double *const *alpha,
                                    const double *const *beta,
                                    const double *const *gamma)
{
   if (!internalFEIsAllocated_)
      notAllocatedAbort("FEI_Implementation::loadElemBCs");

   int index = index_current_filter_;
   if (solveType_ == 2) index = index_soln_filter_;

   CHK_ERR( filter_[index]->loadElemBCs(numElems,
                               elemIDs, fieldID,
                               alpha, beta, gamma))

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::sumInElem(GlobalID elemBlockID,
                        GlobalID elemID,
                        const GlobalID* elemConn,
                        const double* const* elemStiffness,
                        const double* elemLoad,
                        int elemFormat)
{
  if (!internalFEIsAllocated_) {
    notAllocatedAbort("FEI_Implementation::sumInElem");
  }

  CHK_ERR( filter_[index_current_filter_]->sumInElem(elemBlockID, elemID,
                                                     elemConn, elemStiffness,
                                                     elemLoad, elemFormat));

  newMatrixDataLoaded_ = 1;

  return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::sumInElemMatrix(GlobalID elemBlockID,
                        GlobalID elemID,
                        const GlobalID* elemConn,
                        const double* const* elemStiffness,
                        int elemFormat)
{
  if (!internalFEIsAllocated_)
    notAllocatedAbort("FEI_Implementation::sumInElemMatrix");

  CHK_ERR( filter_[index_current_filter_]->sumInElemMatrix(elemBlockID,
                                          elemID, elemConn,
                                          elemStiffness, elemFormat))

  newMatrixDataLoaded_ = 1;


  return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::sumInElemRHS(GlobalID elemBlockID,
                        GlobalID elemID,
                        const GlobalID* elemConn,
                        const double* elemLoad)
{
   if (!internalFEIsAllocated_)
      notAllocatedAbort("FEI_Implementation::sumInElemRHS");

   CHK_ERR( filter_[index_current_rhs_row_]->sumInElemRHS(elemBlockID,
                                          elemID, elemConn, elemLoad))

   newMatrixDataLoaded_ = 1;

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::loadCRMult(int CRID,
                                   int numCRNodes,
                                   const GlobalID* CRNodes,
                                   const int* CRFields,
                                   const double* CRWeights,
                                   double CRValue)
{
   if (!internalFEIsAllocated_)
      notAllocatedAbort("FEI_Implementation::loadCRMult");

   newMatrixDataLoaded_ = 1;

   CHK_ERR( filter_[index_current_filter_]->loadCRMult(CRID,
                                           numCRNodes, CRNodes,
                                           CRFields, CRWeights, CRValue));

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::loadCRPen(int CRID,
                                  int numCRNodes,
                                  const GlobalID* CRNodes,
                                  const int* CRFields,
                                  const double* CRWeights,
                                  double CRValue,
                                  double penValue)
{
   if (!internalFEIsAllocated_)
      notAllocatedAbort("FEI_Implementation::loadCRPen");

   CHK_ERR( filter_[index_current_filter_]->loadCRPen(CRID,
                                          numCRNodes, CRNodes,
                                          CRFields, CRWeights,
                                          CRValue, penValue))

   newMatrixDataLoaded_ = 1;

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::sumIntoRHS(int IDType,
                                   int fieldID,
                                   int numIDs,
                                   const GlobalID* IDs,
                                   const double* rhsEntries)
{
  if (!internalFEIsAllocated_)
    notAllocatedAbort("FEI_Implementation::sumIntoRHS");

  CHK_ERR( filter_[index_current_rhs_row_]->sumIntoRHS(IDType, fieldID,
                                                       numIDs, IDs, rhsEntries) );
  newMatrixDataLoaded_ = 1;

  return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::sumIntoMatrixDiagonal(int IDType,
                             int fieldID,
                             int numIDs,
                             const GlobalID* IDs,
                             const double* coefficients)
{
  if (!internalFEIsAllocated_)
    notAllocatedAbort("FEI_Implementation::sumIntoMatrixDiagonal");

  CHK_ERR( filter_[index_current_filter_]->sumIntoMatrixDiagonal(IDType, fieldID,
                                                       numIDs, IDs, coefficients) );
  newMatrixDataLoaded_ = 1;

  return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::putIntoRHS(int IDType,
                                   int fieldID,
                                   int numIDs,
                                   const GlobalID* IDs,
                                   const double* rhsEntries)
{
  if (!internalFEIsAllocated_)
    notAllocatedAbort("FEI_Implementation::putIntoRHS");

  CHK_ERR( filter_[index_current_rhs_row_]->putIntoRHS(IDType, fieldID,
                                                       numIDs, IDs, rhsEntries) );
  newMatrixDataLoaded_ = 1;

  return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::setMatScalars(int numScalars,
                                      const int* IDs, const double* scalars)
{
   for(int i=0; i<numScalars; i++){
      std::vector<int>::iterator iter =
          std::find(matrixIDs_.begin(), matrixIDs_.end(), IDs[i]);
      if (iter != matrixIDs_.end()) {
         int index = iter - matrixIDs_.begin();
         matScalars_[index] = scalars[i];
      }
      else {
         fei::console_out() << "FEI_Implementation::setMatScalars: ERROR, invalid ID supplied"
              << FEI_ENDL;
         return(1);
      }
   }

   aggregateSystemFormed_ = false;
   matScalarsSet_ = true;

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::setRHSScalars(int numScalars,
                                      const int* IDs, const double* scalars)
{
  for(int i=0; i<numScalars; i++){
     bool found = false;

     for(int j=0; j<numInternalFEIs_; j++){
         int index = fei::searchList(IDs[i], rhsIDs_[j], numRHSIDs_[j]);
         if (index>=0) {
             rhsScalars_[j][index] = scalars[i];
             found = true;
             break;
         }
      }

      if (!found) {
         fei::console_out() << "FEI_Implementation::setRHSScalars: ERROR, invalid RHS ID supplied"
             << FEI_ENDL;
         return(1);
      }
   }

   aggregateSystemFormed_ = false;
   rhsScalarsSet_ = true;

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::parameters(int numParams, const char *const* paramStrings)
{
  // this function takes parameters and passes them to the internal
  // fei objects.

  if (numParams == 0 || paramStrings == NULL) {
    debugOut("#--- no parameters");
    return(0);
  }

  // merge these parameters with any others we may have, for later use.
  snl_fei::mergeStringLists(paramStrings_, numParams_,
                                   paramStrings, numParams);

  snl_fei::getIntParamValue("numMatrices", numParams,paramStrings, numInternalFEIs_);

  snl_fei::getIntParamValue("outputLevel", numParams,paramStrings, outputLevel_);

  const char* param = snl_fei::getParamValue("debugOutput",numParams,paramStrings);
  if (param != NULL) {
    setDebugOutput(param,"FEI_log");
  }

  if (debugOutput_) {
    (*dbgOStreamPtr_)<<"FEI: parameters"<<FEI_ENDL;
    (*dbgOStreamPtr_)<<"#FEI_Implementation, num-params "<<FEI_ENDL
                     <<numParams<<FEI_ENDL;
    (*dbgOStreamPtr_)<<"# "<<numParams<<" parameter lines follow:"<<FEI_ENDL;
    for(int i=0; i<numParams; i++){
      (*dbgOStreamPtr_)<<paramStrings[i]<<FEI_ENDL;
    }
  }

  if (haveLinSysCore_) {
    linSysCore_->parameters(numParams, (char**)paramStrings);
  }
  if (haveFEData_) {
    wrapper_->getFiniteElementData()->parameters(numParams, (char**)paramStrings);
  }

  problemStructure_->parameters(numParams, paramStrings);

  if (internalFEIsAllocated_){
    for(int i=0; i<numInternalFEIs_; i++){
      CHK_ERR( filter_[i]->parameters(numParams, paramStrings) );
    }
  }

  debugOut("#FEI_Implementation leaving parameters method");

  return(0);
}

//------------------------------------------------------------------------------
void FEI_Implementation::setDebugOutput(const char* path, const char* name)
{
  //
  //This function turns on debug output, and opens a file to put it in.
  //
  if (dbgFileOpened_) {
    dbgFStreamPtr_->close();
  }

  dbgFileOpened_ = false;
  delete dbgOStreamPtr_;

  FEI_OSTRINGSTREAM osstr;
  if (path != NULL) {
    osstr << path << "/";
  }
  osstr << name << "." << numProcs_ << "." << localRank_;

  debugOutput_ = 1;
  dbgFStreamPtr_ = new FEI_OFSTREAM(osstr.str().c_str(), IOS_APP);
  if (!dbgFStreamPtr_ || dbgFStreamPtr_->bad()){
    fei::console_out() << "couldn't open debug output file: " << osstr.str() << FEI_ENDL;
    debugOutput_ = 0;
  }

  if (debugOutput_) {
    const char* version_str = NULL;
    version(version_str);

    (*dbgFStreamPtr_) << version_str << FEI_ENDL;

    problemStructure_->setDbgOut(*dbgFStreamPtr_, path, "_0");
    dbgOStreamPtr_ = dbgFStreamPtr_;
    dbgOStreamPtr_->setf(IOS_SCIENTIFIC, IOS_FLOATFIELD);
    dbgFileOpened_ = true;

    if (internalFEIsAllocated_) {
      for(int i=0; i<numInternalFEIs_; ++i) {
        filter_[i]->setLogStream(dbgOStreamPtr_);
      }
    }
  }
}

//------------------------------------------------------------------------------
int FEI_Implementation::loadComplete(bool applyBCs,
                                     bool globalAssemble)
{
  (void)applyBCs;
  (void)globalAssemble;

  buildLinearSystem();

  return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::residualNorm(int whichNorm, int numFields,
                                     int* fieldIDs, double* norms)
{
   buildLinearSystem();

   double residTime = 0.0;

   int err = filter_[index_soln_filter_]->residualNorm(whichNorm, numFields,
                                                 fieldIDs, norms, residTime);

   solveTime_ += residTime;

   return(err);
}

//------------------------------------------------------------------------------
int FEI_Implementation::solve(int& status)
{
   buildLinearSystem();

   double sTime = 0.0;

   int err = filter_[index_soln_filter_]->solve(status, sTime);

   solveTime_ += sTime;

   return(err);
}

//------------------------------------------------------------------------------
int FEI_Implementation::iterations(int& itersTaken) const {
  itersTaken = filter_[index_soln_filter_]->iterations();
  return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::version(const char*& versionString)
{
  versionString = fei::utils::version();
  return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::cumulative_cpu_times(double& initTime,
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

//------------------------------------------------------------------------------
int FEI_Implementation::getBlockNodeSolution(GlobalID elemBlockID,
                                             int numNodes,
                                             const GlobalID *nodeIDs,
                                             int *offsets,
                                             double *results)
{
   CHK_ERR(filter_[index_soln_filter_]->getBlockNodeSolution(elemBlockID,
                                               numNodes,
                                               nodeIDs,
                                               offsets,
                                               results))
   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::getNodalSolution(int numNodes,
                                         const GlobalID *nodeIDs,
                                         int *offsets,
                                         double *results)
{
   CHK_ERR(filter_[index_soln_filter_]->getNodalSolution(numNodes,
                                                         nodeIDs,
                                                         offsets,
                                                         results))

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::getBlockFieldNodeSolution(GlobalID elemBlockID,
                                        int fieldID,
                                        int numNodes, 
                                        const GlobalID *nodeIDs, 
                                        double *results)
{
   CHK_ERR( filter_[index_soln_filter_]->getBlockFieldNodeSolution(elemBlockID,
                                                     fieldID, numNodes,
                                                     nodeIDs, results))

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::putBlockNodeSolution(GlobalID elemBlockID,
                                             int numNodes,
                                             const GlobalID *nodeIDs,
                                             const int *offsets,
                                             const double *estimates)
{
   CHK_ERR( filter_[index_soln_filter_]->putBlockNodeSolution(elemBlockID,
                                                numNodes, nodeIDs,
                                                offsets, estimates))
   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::putBlockFieldNodeSolution(GlobalID elemBlockID,
                                        int fieldID,
                                        int numNodes,
                                        const GlobalID *nodeIDs,
                                        const double *estimates)
{
   int err = filter_[index_soln_filter_]->putBlockFieldNodeSolution(elemBlockID,
                                                     fieldID, numNodes,
                                                     nodeIDs, estimates);
   return(err);
}

//------------------------------------------------------------------------------
int FEI_Implementation::getBlockElemSolution(GlobalID elemBlockID,  
                                   int numElems, 
                                   const GlobalID *elemIDs,
                                   int& numElemDOFPerElement,
                                   double *results)
{
   CHK_ERR( filter_[index_soln_filter_]->getBlockElemSolution(elemBlockID,
                                                       numElems, elemIDs,
                                                       numElemDOFPerElement,
                                                       results))
    return(0);
} 
      
//------------------------------------------------------------------------------
int FEI_Implementation::putBlockElemSolution(GlobalID elemBlockID,
                                   int numElems,
                                   const GlobalID *elemIDs,
                                   int dofPerElem,
                                   const double *estimates)
{
   CHK_ERR( filter_[index_soln_filter_]->putBlockElemSolution(elemBlockID,
                                                       numElems, elemIDs,
                                                       dofPerElem, estimates))
   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::getNumCRMultipliers(int& numMultCRs)
{
  numMultCRs = problemStructure_->getNumMultConstRecords();
  return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::getCRMultIDList(int numMultCRs,
                                        int* multIDs)
{
   if (numMultCRs > problemStructure_->getNumMultConstRecords()) return(-1);

   std::map<GlobalID,snl_fei::Constraint<GlobalID>*>::const_iterator
     cr_iter = problemStructure_->getMultConstRecords().begin(),
     cr_end = problemStructure_->getMultConstRecords().end();
   int i = 0;
   while(cr_iter != cr_end) {
      multIDs[i++] = (*cr_iter).first;
      ++cr_iter;
   }

   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::getCRMultipliers(int numMultCRs,
                                         const int* CRIDs,
                                         double* multipliers)
{
   CHK_ERR( filter_[index_soln_filter_]->getCRMultipliers(numMultCRs,
                                                     CRIDs, multipliers))
   return(0);
}

//------------------------------------------------------------------------------
int FEI_Implementation::putCRMultipliers(int numMultCRs,
                                         const int* CRIDs,
                                         const double *multEstimates)
{
    return(
    filter_[index_soln_filter_]->putCRMultipliers(numMultCRs,
                                            CRIDs,
                                            multEstimates)
    );
}

//-----------------------------------------------------------------------------
//  some utility functions to aid in using the "put" functions for passing
//  an initial guess to the solver
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
int FEI_Implementation::getBlockElemIDList(GlobalID elemBlockID,
                                           int numElems,
                                           GlobalID* elemIDs)
{
  //
  //  return the list of element IDs for a given block... the length parameter
  //  lenElemIDList should be used to check memory allocation for the calling
  //  method, as the calling method should have gotten a copy of this param
  //  from a call to getNumBlockElements before allocating memory for elemIDList
  //
  ConnectivityTable& connTable = problemStructure_->
    getBlockConnectivity(elemBlockID);

  std::map<GlobalID,int>& elemIDList = connTable.elemIDs;
  int len = elemIDList.size();
  if (len > numElems) len = numElems;

  fei::copyKeysToArray(elemIDList, len, elemIDs);

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEI_Implementation::getBlockNodeIDList(GlobalID elemBlockID,
                                           int numNodes,
                                           GlobalID *nodeIDs)
{
  //
  //  similar comments as for getBlockElemIDList(), except for returning the
  //  active node list
  //
  int numActiveNodes = problemStructure_->getNumActiveNodes();
  NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();

  int blk_idx = problemStructure_->getIndexOfBlock(elemBlockID);
  int offset = 0;
  for(int i=0; i<numActiveNodes; i++) {
    const NodeDescriptor* node = NULL;
    nodeDB.getNodeAtIndex(i, node);
    if (node==NULL) continue;
    if (node->hasBlockIndex(blk_idx))
      nodeIDs[offset++] = node->getGlobalNodeID();
    if (offset == numNodes) break;
  }

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEI_Implementation::getNumNodesPerElement(GlobalID blockID,
                                              int& nodesPerElem) const
{
  //
  //  return the number of nodes associated with elements of a given block ID
  //
  BlockDescriptor* block = NULL;
  CHK_ERR( problemStructure_->getBlockDescriptor(blockID, block) );

  nodesPerElem = block->getNumNodesPerElement();
  return(FEI_SUCCESS);
}
 
//------------------------------------------------------------------------------
int FEI_Implementation::getNumEqnsPerElement(GlobalID blockID, int& numEqns)
const
{
  //
  //  return the number of eqns associated with elements of a given block ID
  //
  BlockDescriptor* block = NULL;
  CHK_ERR( problemStructure_->getBlockDescriptor(blockID, block) );

  numEqns = block->getNumEqnsPerElement();
  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEI_Implementation::getNumSolnParams(GlobalID nodeID, int& numSolnParams)
const {
  //
  //  return the number of solution parameters at a given node
  //
  const NodeDescriptor* node = NULL;
  int err = problemStructure_->getNodeDatabase().getNodeWithID(nodeID, node);

  if (err != 0) {
    ERReturn(-1);
  }

  numSolnParams = node->getNumNodalDOF();
  return(0);
}
 
//------------------------------------------------------------------------------
int FEI_Implementation::getNumElemBlocks(int& numElemBlocks) const
{
  //
  //  the number of element blocks
  //
  numElemBlocks = problemStructure_->getNumElemBlocks();
  return( 0 );
}

//------------------------------------------------------------------------------
int FEI_Implementation::getNumBlockActNodes(GlobalID blockID, int& numNodes)
const {
  //
  //  return the number of active nodes associated with a given element block ID
  //
  BlockDescriptor* block = NULL;
  CHK_ERR( problemStructure_->getBlockDescriptor(blockID, block) );

  numNodes = block->getNumActiveNodes();
  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEI_Implementation::getNumBlockActEqns(GlobalID blockID, int& numEqns)
const {
//
// return the number of active equations associated with a given element
// block ID
//
  BlockDescriptor* block = NULL;
  CHK_ERR( problemStructure_->getBlockDescriptor(blockID, block) );

  numEqns = block->getTotalNumEqns();
  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEI_Implementation::getNumBlockElements(GlobalID blockID, int& numElems) const {
//
//  return the number of elements associated with a given elem blockID
//
  BlockDescriptor* block = NULL;
  CHK_ERR( problemStructure_->getBlockDescriptor(blockID, block) );

  numElems = block->getNumElements();
  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEI_Implementation::getNumBlockElemDOF(GlobalID blockID,
                                           int& DOFPerElem) const
{
  //
  //  return the number of elem equations associated with a given elem blockID
  //
  BlockDescriptor* block = NULL;
  CHK_ERR( problemStructure_->getBlockDescriptor(blockID, block) );

  DOFPerElem = block->getNumElemDOFPerElement();

  return(FEI_SUCCESS);
}

//------------------------------------------------------------------------------
int FEI_Implementation::getFieldSize(int fieldID,
                                     int& numScalars)
{
  //
  //  return the number of scalars associated with a given fieldID
  //

  numScalars = problemStructure_->getFieldSize(fieldID);
  return(0);  
}

//------------------------------------------------------------------------------
int FEI_Implementation::getEqnNumbers(GlobalID ID, int idType,
                                      int fieldID, int& numEqns,
                                      int* eqnNumbers)
{
  //
  // Translate from an ID/fieldID pair to a list of equation-numbers
  //

  return( problemStructure_->getEqnNumbers(ID, idType, fieldID,
                                           numEqns, eqnNumbers) );
}

//------------------------------------------------------------------------------
int FEI_Implementation::getNodalFieldSolution(int fieldID,
                                              int numNodes,
                                              const GlobalID* nodeIDs,
                                              double* results)
{
  return( filter_[index_soln_filter_]->getNodalFieldSolution(fieldID, numNodes,
                                                       nodeIDs, results) );
}

//------------------------------------------------------------------------------
int FEI_Implementation::getNumLocalNodes(int& numNodes)
{
  numNodes = problemStructure_->getNodeDatabase().getNodeIDs().size();
  return( 0 );
}

//------------------------------------------------------------------------------
int FEI_Implementation::getLocalNodeIDList(int& numNodes,
                                           GlobalID* nodeIDs,
                                           int lenNodeIDs)
{
  std::map<GlobalID,int>& nodes =
    problemStructure_->getNodeDatabase().getNodeIDs();
  numNodes = nodes.size();
  int len = numNodes;
  if (lenNodeIDs < len) len = lenNodeIDs;

  fei::copyKeysToArray(nodes, len, nodeIDs);

  return( 0 );
}

//------------------------------------------------------------------------------
int FEI_Implementation::putNodalFieldData(int fieldID,
                                          int numNodes,
                                          const GlobalID* nodeIDs,
                                          const double* nodeData)
{
  return( filter_[index_soln_filter_]->putNodalFieldData(fieldID, numNodes,
                                                   nodeIDs, nodeData) );
}

//------------------------------------------------------------------------------
void FEI_Implementation::buildLinearSystem()
{
  //
  //Private function.
  //
  //At the point when this function is called, all matrix assembly has
  //already taken place, with the data having been directed into the
  //appropriate Filter instance in the filter_ list. Now it's
  //time to finalize the linear system (get A,x,b ready to give to a
  //solver), performing this list of last-minute tasks:
  //
  //1. Have each Filter instance exchangeRemoteEquations.
  //2. Aggregate the system (form a linear combination of LHS's and
  //   RHS's) if solveType_==2.
  //3. call loadComplete to have the 'master' Filter instance
  //   (filter_[index_soln_filter_]) enforce any boundary conditions
  //   that have been loaded.
  //
  debugOut("#   buildLinearSystem");

  //
  //figure out if new matrix-data was loaded on any processor.
  int anyDataLoaded = newMatrixDataLoaded_;
#ifndef FEI_SER
  if (numProcs_ > 1) {
    if (MPI_Allreduce(&newMatrixDataLoaded_, &anyDataLoaded, 1, MPI_INT,
                      MPI_SUM, comm_) != MPI_SUCCESS) voidERReturn
  }
#endif

  if (anyDataLoaded) {
#ifndef FEI_SER
    for(int i=0; i<numInternalFEIs_; i++) {
      filter_[i]->exchangeRemoteEquations();
    }
#endif
    newMatrixDataLoaded_ = 0;
  }

  if (solveType_ == 2){
    //solveType_ == 2 means this is a linear-combination solve --
    //i.e., we're solving an aggregate system which is the sum of
    //several individual matrices and rhs's.

    if (!aggregateSystemFormed_) {
      if (!matScalarsSet_ || !rhsScalarsSet_) {
        FEI_COUT << "FEI_Implementation: WARNING: solveType_==2, aggregating system, but setMatScalars and/or setRHSScalars not yet called." << FEI_ENDL;
      }

      aggregateSystem();
    }
  }

  filter_[index_soln_filter_]->loadComplete();

  debugOut("#   leaving buildLinearSystem");
}

//------------------------------------------------------------------------------
int FEI_Implementation::aggregateSystem()
{
  debugOut("#   aggregateSystem");
   if (!haveLinSysCore_) ERReturn(-1);

   if (soln_fei_matrix_ == NULL) {
      soln_fei_matrix_ = new Data();

      CHK_ERR( lscArray_[index_soln_filter_]->
               copyOutMatrix(1.0, *soln_fei_matrix_) );
   }

   if (soln_fei_vector_ == NULL) {
      soln_fei_vector_ = new Data();

      CHK_ERR( lscArray_[index_soln_filter_]->
                     setRHSID(rhsIDs_[index_soln_filter_][0]) );

      CHK_ERR( lscArray_[index_soln_filter_]->
                    copyOutRHSVector(1.0, *soln_fei_vector_) );
   }

   Data tmp;
   Data tmpv;

   for(int i=0; i<numInternalFEIs_; i++){

      if (i == index_soln_filter_) {
         tmp.setTypeName(soln_fei_matrix_->getTypeName());
         tmp.setDataPtr(soln_fei_matrix_->getDataPtr());
         CHK_ERR( lscArray_[index_soln_filter_]->
                        copyInMatrix(matScalars_[i], tmp) );
      }
      else {
         CHK_ERR( lscArray_[i]->getMatrixPtr(tmp) );

         CHK_ERR( lscArray_[index_soln_filter_]->
                        sumInMatrix(matScalars_[i], tmp) );
      }

      for(int j=0; j<numRHSIDs_[i]; j++){
         if ((i == index_soln_filter_) && (j == 0)) {
            tmpv.setTypeName(soln_fei_vector_->getTypeName());
            tmpv.setDataPtr(soln_fei_vector_->getDataPtr());
         }
         else {
            CHK_ERR( lscArray_[i]->setRHSID(rhsIDs_[i][j]) );
            CHK_ERR( lscArray_[i]->getRHSVectorPtr(tmpv) );
         }

         if (i == index_soln_filter_) {
            CHK_ERR( lscArray_[index_soln_filter_]->
                         copyInRHSVector(rhsScalars_[i][j], tmpv) );
         }
         else {
            CHK_ERR( lscArray_[index_soln_filter_]->
                         sumInRHSVector(rhsScalars_[i][j], tmpv) );
         }
      }
   }

   aggregateSystemFormed_ = true;

   debugOut("#   leaving aggregateSystem");

   return(0);
}

//==============================================================================
int FEI_Implementation::allocateInternalFEIs(){
//
//This is a private FEI_Implementation function, to be called from within
//setSolveType or the other overloading of allocateInternalFEIs.
//Assumes that numInternalFEIs_ has already been set.
//It is also safe to assume that problemStructure_->initComplete() has already
//been called.
//

   if (internalFEIsAllocated_) return(0);

   matScalars_.resize(numInternalFEIs_);

   if (rhsScalars_.size() != 0) deleteRHSScalars();

   rhsScalars_.resize(numInternalFEIs_);

   for(int i=0; i<numInternalFEIs_; i++){
      matScalars_[i] = 1.0;

      rhsScalars_[i] = new double[numRHSIDs_[i]];

      for(int j=0; j<numRHSIDs_[i]; j++){
         rhsScalars_[i][j] = 1.0;
      }
   }

   IDsAllocated_ = true;

   if (numInternalFEIs_ > 0) {
      index_soln_filter_ = 0;
      index_current_filter_ = 0;
      filter_ = new Filter*[numInternalFEIs_];

      if (haveLinSysCore_) {
        if (numRHSIDs_[0] == 0) {
          int dummyID = -1;
          linSysCore_->setNumRHSVectors(1, &dummyID);
        }
        else {
          linSysCore_->setNumRHSVectors(numRHSIDs_[0], rhsIDs_[0]);
        }

        for(int i=1; i<numInternalFEIs_; i++) {
          fei::SharedPtr<LinearSystemCore> lsc(linSysCore_->clone());
          lsc->parameters(numParams_, paramStrings_);

          if (numRHSIDs_[i] == 0) {
            int dummyID = -1;
            lsc->setNumRHSVectors(1, &dummyID);
          }
          else {
            lsc->setNumRHSVectors(numRHSIDs_[i], rhsIDs_[i]);
          }

          lscArray_.push_back(lsc);
        }
      }

      for(int i=0; i<numInternalFEIs_; i++){

        if (haveLinSysCore_) {
          filter_[i] = new LinSysCoreFilter(this, comm_, problemStructure_,
                                            lscArray_[i].get(), masterRank_);
        }
        else if (haveFEData_) {
          filter_[i] = new FEDataFilter(this, comm_, problemStructure_,
                                        wrapper_.get(), masterRank_);
        }
        else {
          fei::console_out() << "FEI_Implementation: ERROR, don't have LinearSystemCore"
               << " or FiniteElementData implementation..." << FEI_ENDL;
          ERReturn(-1);
        }

        filter_[i]->setLogStream(dbgOStreamPtr_);

        FEI_OSTRINGSTREAM osstr;
        osstr<<"internalFei "<< i;
        std::string osstr_str = osstr.str();
        const char* param = osstr_str.c_str();
        filter_[i]->parameters(1, &param);

        if (debugOutput_) {
          (*dbgOStreamPtr_)<<"#-- fei["<<i<<"]->setNumRHSVectors "
                           <<numRHSIDs_[i]<<FEI_ENDL;
        }

        if (numRHSIDs_[i] == 0) {
          int dummyID = -1;
          filter_[i]->setNumRHSVectors(1, &dummyID);
        }
        else {
          filter_[i]->setNumRHSVectors(numRHSIDs_[i], rhsIDs_[i]);
        }
      }

      internalFEIsAllocated_ = true;
   }
   else {
     needParametersAbort("FEI_Implementation::allocateInternalFEIs");
   }

   return(0);
}

//==============================================================================
void FEI_Implementation::debugOut(const char* msg) {
  if (debugOutput_) { (*dbgOStreamPtr_)<<msg<<FEI_ENDL; }
}

//==============================================================================
void FEI_Implementation::debugOut(const char* msg, int whichFEI) {
   if (debugOutput_) {
     (*dbgOStreamPtr_)<<msg<<", -> fei["<<whichFEI<<"]"<<FEI_ENDL;
   }
}

//==============================================================================
void FEI_Implementation::messageAbort(const char* msg){

    fei::console_out() << "FEI_Implementation: ERROR " << msg << " Aborting." << FEI_ENDL;
    MPI_Abort(comm_, -1);
}

//==============================================================================
void FEI_Implementation::notAllocatedAbort(const char* name){

    fei::console_out() << name
         << FEI_ENDL << "ERROR, internal data structures not allocated."
         << FEI_ENDL << "'setIDLists' and/or 'setSolveType' must be called"
         << FEI_ENDL << "first to identify solveType and number of matrices"
         << FEI_ENDL << "to be assembled." << FEI_ENDL;
    MPI_Abort(comm_, -1);
}

//==============================================================================
void FEI_Implementation::needParametersAbort(const char* name){

   fei::console_out() << name
     << FEI_ENDL << "FEI_Implementation: ERROR, numMatrices has not been specified."
     << FEI_ENDL << "FEI_Implementation: 'parameters' must be called up front with"
     << FEI_ENDL << "FEI_Implementation: the string 'numMatrices n' to specify that"
     << FEI_ENDL << "FEI_Implementation: n matrices will be assembled." << FEI_ENDL;
    MPI_Abort(comm_, -1);
}

//==============================================================================
void FEI_Implementation::badParametersAbort(const char* name){

   fei::console_out() << name
        << FEI_ENDL << "FEI_Implementation: ERROR, inconsistent 'solveType' and"
        << FEI_ENDL << "FEI_Implementation: 'numMatrices' parameters specified."
        << FEI_ENDL << "FEI_Implementation: Aborting."
        << FEI_ENDL;
    MPI_Abort(comm_, -1);
}

