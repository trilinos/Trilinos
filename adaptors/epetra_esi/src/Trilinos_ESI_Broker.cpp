
//Define some architecture-specific symbols so that the compiler takes the right
//branch in Epetra_Object.h

#ifdef __PUMAGON__
#ifndef TFLOP
#define TFLOP
#endif
#endif

#ifdef sgi
#ifndef SGI
#define SGI
#endif
#endif

//We include Epetra_Object.h first to make sure we get the definitions of
//iostream, etc., from the same places that Epetra gets them.
#include "Epetra_Object.h"

#include "esi/ESI.h"

//We need to include the epetra-esi stuff that Trilinos_ESI_Broker is managing..

#ifndef EPETRA_ESI_INCLUDE_IMPLEMENTATION
//we're going to turn on the flag to include implementations, but if it isn't
//already on, we'll need to turn it back off afterwards...
#define Tril_ESI_Brok_UNDEF
#endif

#define EPETRA_ESI_INCLUDE_IMPLEMENTATION
#include "Epetra_ESI_Object.h"
#include "Epetra_ESI_Argv.h"
#include "Epetra_ESI_IndexSpace.h"
#include "Epetra_ESI_Vector.h"
#include "Epetra_ESI_CrsMatrix.h"
#include "Epetra_ESI_Operator.h"
#include "Aztec_ESI_Solver.h"

#ifdef Tril_ESI_Brok_UNDEF
#undef EPETRA_ESI_INCLUDE_IMPLEMENTATION
#undef Tril_ESI_Brok_UNDEF
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "Trilinos_ESI_Broker.h"

//Let's define some macros to make error-checking a little cleaner in the
//code down below...

#ifdef _fnName_
#undef _fnName_
#endif

//_fnName_ should be redefined at each function in the Trilinos_ESI_Broker
//implementation code below.
#define _fnName_ "<unknown function>"

#ifdef CHK_ERR
#undef CHK_ERR
#endif

#define CHK_ERR(a) { int chkErrorCode; if ((chkErrorCode = a) != 0) { \
                      cerr << _fnName_ << ", " << __FILE__ << ", line " \
                           << __LINE__ << " " << chkErrorCode << endl; \
                      return(chkErrorCode); \
                   } }

#ifdef ERReturn
#undef ERReturn
#endif

#define ERReturn(a) { cerr << _fnName_ << ", " << __FILE__ << ", line " \
                           << __LINE__ << endl; return(a); }

#ifdef voidERReturn
#undef voidERReturn
#endif

#define voidERReturn { cerr << _fnName_ << ", " << __FILE__ << ", line " \
			   << __LINE__ << endl; return; }

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::Trilinos_ESI_Broker"
Trilinos_ESI_Broker::Trilinos_ESI_Broker()
 : libName_(NULL),
   comm_(0),
   haveMPIComm_(false),
   petra_comm_(NULL),
   petra_graph_(NULL),
   ispaces_(0,1),
   ispaceNames_(0,1),
   matrices_(0,1),
   matrixNames_(0,1),
   vecs_(0,1),
   vecNames_(0,1),
   setGlobalOffsetsCalled_(false),
   setMatrixStructureCalled_(false),
   solvers_(0,1),
   solverNames_(0,1),
   localProc_(0),
   numProcs_(1),
   globalSize_(0),
   localSize_(0),
   localOffset_(0),
   rowLengths_(0, 256),
   params_(0,1),
   outputLevel_(0),
   debugOutput_(false),
   debugFile_(NULL)
{
   libName_ = new char[10];
   if (libName_ == NULL) voidERReturn;
   sprintf(libName_, "Trilinos");

}

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::Trilinos_ESI_Broker"
Trilinos_ESI_Broker::Trilinos_ESI_Broker(MPI_Comm comm)
 : libName_(NULL),
   comm_(comm),
   haveMPIComm_(true),
#ifdef EPETRA_MPI
   petra_comm_(new Epetra_MpiComm(comm)),
#else
   petra_comm_(new Epetra_SerialComm),
#endif
   petra_graph_(NULL),
   ispaces_(0,1),
   ispaceNames_(0,1),
   matrices_(0,1),
   matrixNames_(0,1),
   vecs_(0,1),
   vecNames_(0,1),
   setGlobalOffsetsCalled_(false),
   setMatrixStructureCalled_(false),
   solvers_(0,1),
   solverNames_(0,1),
   localProc_(0),
   numProcs_(1),
   globalSize_(0),
   localSize_(0),
   localOffset_(0),
   rowLengths_(0, 256),
   params_(0,1),
   outputLevel_(0),
   debugOutput_(false),
   debugFile_(NULL)
{
   libName_ = new char[10];
   if (libName_ == NULL) voidERReturn;
   sprintf(libName_, "Trilinos");

#ifdef EPETRA_MPI
   if (MPI_Comm_rank(comm_, &localProc_) != MPI_SUCCESS) voidERReturn;
   if (localProc_ < 0) voidERReturn;

   if (MPI_Comm_size(comm_, &numProcs_) != MPI_SUCCESS) voidERReturn;
   if (numProcs_ < 1) voidERReturn;
#endif
}

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::~Trilinos_ESI_Broker"
Trilinos_ESI_Broker::~Trilinos_ESI_Broker()
{
  delete[] libName_;

  //The extra curly brackets below are due to some compilers not liking the
  //re-declaration of 'i' in multiple for-loops in the same scope.

  { for(int i=0; i<matrices_.length(); i++)
    {delete matrices_[i]; delete[] matrixNames_[i];}
  }

  { for(int i=0; i<vecs_.length(); i++)
    {delete vecs_[i]; delete[] vecNames_[i];}
  }

  { for(int i=0; i<solvers_.length(); i++) 
    {delete solvers_[i]; delete[] solverNames_[i];}
  }

  { for(int i=0; i<ispaces_.length(); i++) 
    {delete ispaces_[i]; delete[] ispaceNames_[i];}
  }
  
  { for(int i=0; i<params_.length(); i++) delete[] params_[i]; }

  delete petra_graph_;
  delete petra_comm_;

  if (debugOutput_) fclose(debugFile_);
}

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::setMPIComm"
int Trilinos_ESI_Broker::setMPIComm(MPI_Comm comm)
{
  haveMPIComm_ = true;

#ifdef EPETRA_SERIAL
  petra_comm_ = new Epetra_SerialComm;
#else
  petra_comm_ = new Epetra_MpiComm(comm);

  if (MPI_Comm_rank(comm_, &localProc_) != MPI_SUCCESS) ERReturn(-1);
  if (localProc_ < 0) ERReturn(-1);

  if (MPI_Comm_size(comm_, &numProcs_) != MPI_SUCCESS) ERReturn(-1);
  if (numProcs_ < 1) ERReturn(-1);
#endif

  return(0);
}

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::getLibraryName"
int Trilinos_ESI_Broker::getLibraryName(char*& libName)
{
   libName = libName_;
   return(0);
}

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::parameters"
int Trilinos_ESI_Broker::parameters(int numParams, char** paramStrings)
{
   if (numParams <= 0) return(0);

   //We recognize only a couple parameters.

   char result[256];
   char* outputstring = NULL;

   if (getParam("debugOutput", numParams, paramStrings, result)) {
     char* dbgFileName = new char[64];
     char* dbgPath = new char[strlen(result)+1];
     strcpy(dbgPath, result);
     sprintf(dbgFileName, "PEBrkr_dbg.%d.%d", numProcs_, localProc_);

     CHK_ERR( openDebugOutput(dbgPath, dbgFileName) );
     delete [] dbgFileName;
     delete [] dbgPath;
   }

   if (getParam("outputLevel", numParams, paramStrings, result)) {
      sscanf(result, "%d", &outputLevel_);

      outputstring = new char[64];
      if (outputstring == NULL) ERReturn(-1);
      sprintf(outputstring, "AZ_output %d", outputLevel_);

      if (outputLevel_ == 1) {
	if (localProc_ !=  0) outputLevel_ = 0;
      }
   }

   //let's accumulate the parameters into our own copy, so we can give them
   //to solvers when we allocate one, etc.

   if (outputstring != NULL) {
     if (params_.append(outputstring) != 0) ERReturn(-1);
   }

   for(int i=0; i<numParams; i++) {
     char* param = new char[strlen(paramStrings[i])+1];
     if (param == NULL) return(-1);
     sprintf(param, paramStrings[i]);
     if (params_.append(param) != 0) ERReturn(-1);
   }

   return(0);
}

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::clone"
int Trilinos_ESI_Broker::clone(ESI_Broker*& esi_broker)
{
   esi_broker = new Trilinos_ESI_Broker(comm_);
   if (esi_broker == NULL) ERReturn(-1);

   return(0);
}

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::setGlobalOffsets"
int Trilinos_ESI_Broker::setGlobalOffsets(int len, int* nodeOffsets,
                        int* eqnOffsets, int* blkEqnOffsets)
{
   (void)nodeOffsets;
   (void)blkEqnOffsets;

   if (debugOutput_) {
     fprintf(debugFile_, "setGlobalOffsets\n");
     for(int i=0; i<len; i++) {
       fprintf(debugFile_, "   nodeOffsets[%d]: %d, eqnOffsets[%d]: %d, blkEqnOffsets[%d]: %d\n", i, nodeOffsets[i], i, eqnOffsets[i], i, blkEqnOffsets[i]);
       fflush(debugFile_);
     }
   }

   if (len <= localProc_+1) {
     cerr << "Trilinos_ESI_Broker::setGlobalOffsets: ERROR, len is " << len
          << " but localProc_ is " << localProc_ << ". len needs to be > "
          << "localProc_+1. Calling MPI_Abort..." << endl;
#ifdef EPETRA_SERIAL
     abort();
#else
     MPI_Abort(comm_, -1);
#endif
     ERReturn(-1);
   }

   localOffset_ = eqnOffsets[localProc_];
   localSize_   = eqnOffsets[localProc_+1] - localOffset_;
   globalSize_  = eqnOffsets[numProcs_];

   if (debugOutput_) {
     fprintf(debugFile_, "   localOffset_: %d, localSize_: %d\n",
	     localOffset_, localSize_);
     fflush(debugFile_);
   }

   if (localOffset_ < 0) {
     cerr << "Trilinos_ESI_Broker::setGlobalOffsets: localOffset_: " << localOffset_
          << ". Should be >= 0" << endl;
     ERReturn(-1);
   }

   if (localSize_ < 0) {
     cerr << "Trilinos_ESI_Broker::setGlobalOffsets: localSize_: " << localSize_
          << ". Should be >= 0" << endl;
     ERReturn(-1)
   }

   setGlobalOffsetsCalled_ = true;

   //Now let's get an indexspace established so that we'll be ready to service
   //matrix requests after setMatrixStructure is called.
   esi::IndexSpace<int>* ispace = NULL;
   CHK_ERR( getIndexSpaceInstance("ispace1", (void*&)ispace) );

   return(0);
}

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::setMultCREqns"
int Trilinos_ESI_Broker::setMultCREqns(int numCRs, 
                             const int* numNodesPerCR,
                             const int* const* nodeNumbers,
                             const int* const* eqnNumbers,
                             const int* fieldIDs,
                             const int* multiplierEqnNumbers)
{
   (void)numCRs; (void)numNodesPerCR; (void)nodeNumbers;
   (void)eqnNumbers; (void)fieldIDs; (void)multiplierEqnNumbers;
   return(0);
}

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::setMatrixStructure"
int Trilinos_ESI_Broker::setMatrixStructure(int** ptColIndices,
                                   int* ptRowLengths,
                                   int** blkColIndices,
                                   int* blkRowLengths,
                                   int* ptRowsPerBlkRow)
{
  if (!setGlobalOffsetsCalled_) {
    cerr << "Trilinos_ESI_Broker::setMatrixStructure: ERROR, setGlobalOffsets not"
         << " called yet." << endl;
    ERReturn(-1);
  }

  (void)ptColIndices;  (void)blkColIndices;
  (void)blkRowLengths;  (void)ptRowsPerBlkRow;

  rowLengths_.resize(localSize_);
  { for(int i=0; i<localSize_; i++) rowLengths_[i] = ptRowLengths[i]; }

  //
  //Now we can create and fill our Petra_CRS_Graph to use when we create
  //matrices later.
  //
  const char* ispcname = ispaceNames_[ispaces_.length()-1];
  esi::IndexSpace<int>* ispace = NULL;
  CHK_ERR( getIndexSpaceInstance(ispcname, (void*&)ispace) );

  epetra_esi::IndexSpace<int>* epetraIspace = NULL;
  CHK_ERR( ispace->getInterface("epetra_esi::IndexSpace",
				(void*&)epetraIspace) );

  petra_graph_ = new Epetra_CrsGraph(Copy, *epetraIspace, rowLengths_.dataPtr());
  if (petra_graph_ == NULL) ERReturn(-1);

  for(int i=0; i<localSize_; i++) {
    CHK_ERR( petra_graph_-> InsertGlobalIndices(i+localOffset_,
						ptRowLengths[i],
						ptColIndices[i]) );
  }

  CHK_ERR( petra_graph_->TransformToLocal() );

  setMatrixStructureCalled_ = true;
  return(0);
}

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::getInterfaceInstance"
int Trilinos_ESI_Broker::getInterfaceInstance(const char* instanceName,
					   const char* interfaceName,
					   void*& objectPtr)
{
  //We recognize the following interface requests:
  //
  //   Any of these matrix interfaces:
  //      esi::Operator
  //      esi::MatrixRowReadAccess
  //      esi::MatrixRowWriteAccess
  //      esi::MatrixData
  //      epetra_esi::CrsMatrix
  //
  //   Either of these vector interfaces:
  //      esi::Vector
  //      epetra_esi::Vector
  //
  //   Either of these solver interfaces:
  //      esi::Solver
  //      esi::SolverIterative
  //
  //   Either of these indexspace interfaces:
  //      esi::IndexSpace
  //      Epetra_ESI_IndexSpace

  if (stringsMatch(interfaceName, "esi::Operator") ||
      stringsMatch(interfaceName, "esi::MatrixRowReadAccess") ||
      stringsMatch(interfaceName, "esi::MatrixRowWriteAccess") ||
      stringsMatch(interfaceName, "esi::MatrixData") ||
      stringsMatch(interfaceName, "epetra_esi::CrsMatrix") ) {

    if (!setMatrixStructureCalled_) return(-1);

    esi::Object* mat = NULL;
    CHK_ERR( getMatrixInstance(ispaceNames_[ispaces_.length()-1], instanceName,
			       (void*&)mat) );

    CHK_ERR( mat->getInterface(interfaceName, objectPtr) );
    return(0);
  }

  if (stringsMatch(interfaceName, "esi::Vector") ||
      stringsMatch(interfaceName, "epetra_esi::Vector") ) {

    if (!setGlobalOffsetsCalled_) return(-1);

    esi::Vector<double,int>* vec = NULL;
    CHK_ERR( getVectorInstance(ispaceNames_[ispaces_.length()-1], instanceName,
			       (void*&)vec) );

    CHK_ERR( vec->getInterface(interfaceName, objectPtr) );
    return(0);
  }

  if (stringsMatch(interfaceName, "esi::Solver") ||
      stringsMatch(interfaceName, "esi::SolverIterative")) {

    esi::Object* solver = NULL;
    CHK_ERR( getSolverInstance(instanceName, (void*&)solver) );

    CHK_ERR( solver->getInterface(interfaceName, objectPtr) )
    return(0);
  }

  if (stringsMatch(interfaceName, "esi::IndexSpace") ||
      stringsMatch(interfaceName, "Epetra_ESI_IndexSpace") ) {

    if (!setGlobalOffsetsCalled_) return(-1);

    esi::IndexSpace<int>* ispace = NULL;
    CHK_ERR( getIndexSpaceInstance(instanceName, (void*&)ispace) );

    CHK_ERR( ispace->getInterface(interfaceName, objectPtr) )
    return(0);
  }

  return(-1);
}

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::setInterfaceInstance"
int Trilinos_ESI_Broker::setInterfaceInstance(const char* instanceName,
				       const char* interfaceName,
                                       void* objectPtr)
{
  //
  //We accept these interfaces:
  //
  // vector
  //    esi::Vector
  //
  // matrix
  //    esi::MatrixData
  //
  // lookup
  //    Lookup
  //
  // solver
  //    esi::Solver
  //

  if (stringsMatch(interfaceName, "esi::Vector")) {
    int err = vecs_.append((esi::Vector<double,int>*)objectPtr);
    if (err != 0) return(-1);

    char* name = new char[strlen(instanceName)+1];
    if (name == NULL) return(-1);

    sprintf(name, instanceName);
    err = vecNames_.append(name);
    return(err);
  }

  if (stringsMatch(interfaceName, "esi::MatrixData") ||
      stringsMatch(interfaceName, "esi::Operator") ||
      stringsMatch(interfaceName, "esi::MatrixRowReadAccess") ||
      stringsMatch(interfaceName, "esi::MatrixRowWriteAccess")) {

    int err = matrices_.append((esi::Object*)objectPtr);
    if (err != 0) return(-1);

    char* name = new char[strlen(instanceName)+1];
    if (name == NULL) return(-1);

    sprintf(name, instanceName);
    err = matrixNames_.append(name);
    return(err);
  }

  if (stringsMatch(interfaceName, "esi::Solver") ||
      stringsMatch(interfaceName, "esi::SolverIterative")) {

    int err = solvers_.append((esi::Object*)objectPtr);
    if (err != 0) return(-1);

    char* name = new char[strlen(instanceName)+1];
    if (name == NULL) return(-1);

    sprintf(name, instanceName);
    err = solverNames_.append(name);
    return(err);
  }

  return(-1);
}

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::openDebugOutput"
int Trilinos_ESI_Broker::openDebugOutput(const char* path, const char* fileName)
{
   char* dbgFileName = new char[strlen(path)+strlen(fileName)+5];
   sprintf(dbgFileName, "%s/%s", path, fileName);

   debugFile_ = fopen(dbgFileName, "w");
   if (debugFile_ == NULL) ERReturn(-1)

   delete [] dbgFileName;

   debugOutput_ = true;

   return(0);
}

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::getIndexSpaceInstance"
int Trilinos_ESI_Broker::getIndexSpaceInstance(const char* mapName, void*& mapPtr)
{
  if (!setGlobalOffsetsCalled_) return(-1);

  //If we already have a map with the specified name, return it.
  for(int i=0; i<ispaces_.length(); i++) {
    if (stringsMatch(mapName, ispaceNames_[i])) {
      mapPtr = ispaces_[i]; return(0);
    }
  }

  //Since we don't already have the specified map, create a new one and
  //return that.
  esi::IndexSpace<int>* newmap = 
    new epetra_esi::IndexSpace<int>(globalSize_, localSize_, 0, *petra_comm_);
  if (newmap == NULL) return(-1);
  int err = ispaces_.append(newmap);
  if (err != 0) return(err);


  char* newmapname = new char[strlen(mapName)+1];
  if (newmapname == NULL) return(-1);
  sprintf(newmapname,mapName);
  err = ispaceNames_.append(newmapname);
  if (err != 0) return(-1);

  mapPtr = newmap;
  return(0);
}

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::getVectorInstance"
int Trilinos_ESI_Broker::getVectorInstance(const char* mapName,
				       const char* vecName, void*& vecPtr)
{
  if (!setGlobalOffsetsCalled_) return(-1);

  //If we already have a vector with the specified name, return it.
  for(int i=0; i<vecs_.length(); i++) {
    if (stringsMatch(vecName, vecNames_[i])) {
      vecPtr = vecs_[i]; return(0);
    }
  }

  //Since we don't already have the specified vector, create a new one and
  //return that, after first obtaining a map to construct it with.
  esi::IndexSpace<int>* map = NULL;
  CHK_ERR( getIndexSpaceInstance(mapName, (void*&)map) );

  //If it isn't a epetra_esi::IndexSpace, then we can't use it.
  epetra_esi::IndexSpace<int>* epetraIspace = NULL;
  CHK_ERR( map->getInterface("epetra_esi::IndexSpace", (void*&)epetraIspace) );

  esi::Vector<double,int>* newvec = 
    new epetra_esi::Vector<double,int>(*epetraIspace);
  if (newvec == NULL) return(-1);

  int err = vecs_.append(newvec);
  if (err != 0) return(err);

  char* newvecname = new char[strlen(vecName)+1];
  if (newvecname == NULL) return(-1);

  sprintf(newvecname,vecName);
  err = vecNames_.append(newvecname);
  if (err != 0) return(-1);

  vecPtr = newvec;
  return(0);
}

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::getMatrixInstance"
int Trilinos_ESI_Broker::getMatrixInstance(const char* mapName,
				       const char* matName, void*& matPtr)
{
  if (!setMatrixStructureCalled_) return(-1);

  //If we already have a matrix with the specified name, return it.
  for(int i=0; i<matrices_.length(); i++) {
    if (stringsMatch(matName, matrixNames_[i])) {
      matPtr = matrices_[i]; return(0);
    }
  }

  //Since we don't already have the specified matrix, create a new one and
  //return that.

  epetra_esi::CrsMatrix<double,int>* newmat = 
    new epetra_esi::CrsMatrix<double,int>(Copy, *petra_graph_);
  if (newmat == NULL) return(-1);

  esi::Object* objPtr = NULL;
  int err = newmat->getInterface("esi::Object", (void*&)objPtr);
  if (err != 0) return(err);

  err = matrices_.append(objPtr);
  if (err != 0) return(err);

  char* newmatname = new char[strlen(matName)+1];
  if (newmatname == NULL) return(-1);

  sprintf(newmatname,matName);
  err = matrixNames_.append(newmatname);
  if (err != 0) return(-1);

  matPtr = objPtr;
  return(0);
}

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::getSolverInstance"
int Trilinos_ESI_Broker::getSolverInstance(const char* slvName, void*& slvPtr)
{
  //If we already have a solver with the specified name, return it.
  for(int i=0; i<solvers_.length(); i++) {
    if (stringsMatch(slvName, solverNames_[i])) {
      slvPtr = solvers_[i]; return(0);
    }
  }

  //Since we don't already have the specified solver, create a new one and
  //return that.
  esi::Object* newsolver = NULL;
  epetra_esi::CrsMatrix<double,int>* pmat = NULL;
  CHK_ERR( matrices_[matrices_.length()-1]->
	   getInterface("epetra_esi::CrsMatrix", (void*&)pmat) );

  newsolver = 
    new aztecoo_esi::Solver<double,int>(pmat);

  if (newsolver == NULL) return(-1);

  esi::Solver<double,int>* eslv = NULL;
  CHK_ERR( newsolver->getInterface("esi::Solver", (void*&)eslv));
  
  CHK_ERR( eslv->parameters(params_.length(), params_.dataPtr()) );

  esi::Object* objPtr = NULL;
  CHK_ERR( newsolver->getInterface("esi::Object", (void*&)objPtr));

  int err = solvers_.append(objPtr);
  if (err != 0) return(err);

  char* newname = new char[strlen(slvName)+1];
  if (newname == NULL) return(-1);
  sprintf(newname,slvName);
  err = solverNames_.append(newname);
  if (err != 0) return(-1);

  slvPtr = objPtr;
  return(0);
}

//------------------------------------------------------------------------------
#undef _fnName_
#define _fnName_ "Trilinos_ESI_Broker::getParam"
bool Trilinos_ESI_Broker::getParam(const char* flag, int numParams,
                 char** paramStrings, char* result)
{
//
// paramStrings is a collection of string pairs - each string in
// paramStrings consists of two strings separated by a space.
// This function looks through the strings in paramStrings, looking
// for one that contains flag in the first string. The second string
// is then returned in result.
// Assumes that result is allocated by the caller.
//

   char temp[128];

   if (flag == NULL || paramStrings == NULL) {
      return(false); // flag or paramStrings is the NULL pointer
   }

   for(int i = 0; i<numParams; i++) {
      if (paramStrings[i] != NULL)  { // check for NULL pointer
         if (strlen(flag) <= strlen(paramStrings[i])) {
         if (strncmp(flag, paramStrings[i], strlen(flag)) == 0) {
            /* flag found */
            sscanf(paramStrings[i], "%s %s", temp, result);
            return(true);
         }
         }
      }
   }

   return(false);  /* flag was not found in paramStrings */
}

