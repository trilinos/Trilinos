#ifndef _Aztec_ESI_Solver_cpp_
#define _Aztec_ESI_Solver_cpp_

#include "Aztec_ESI_Translator.h"
#include "Aztec_ESI_Solver.h"

#include "Epetra_ESI_utils.h"

//------------------------------------------------------------------------------
template<class Scalar, class Ordinal>
aztecoo_esi::Solver<Scalar,Ordinal>::
Solver(epetra_esi::CrsMatrix<Scalar,Ordinal>* A)
   : epetra_esi::Object(), petraA_(A->getEpetra_CrsMatrix()),
     maxIters_(300), tolerance_(1.e-8), whichConstructor_(0),
     petraAalloced_(false)
{
  aztecoo_ = new AztecOO;
  int err = addInterfaces();
  if (err) msgAbort("aztecoo_esi::Solver ctor addInterfaces() error");
}

//------------------------------------------------------------------------------
template<class Scalar, class Ordinal>
aztecoo_esi::Solver<Scalar,Ordinal>::
Solver(esi::Operator<Scalar, Ordinal>* esi_op)
   : epetra_esi::Object(),
     maxIters_(300), tolerance_(1.e-8), whichConstructor_(1),
     petraAalloced_(false)
{
  aztecoo_ = new AztecOO;
  petraA_ = new epetra_esi::Operator<double,int>(*esi_op);
  if (petraA_ == NULL) {
    msgAbort("aztecoo_esi::Solver ctor failed to alloc epetra_esi::Operator.");
  }

  int err = addInterfaces();
  if (err) msgAbort("aztecoo_esi::Solver ctor addInterfaces() error");
}

//------------------------------------------------------------------------------
template<class Scalar, class Ordinal>
int aztecoo_esi::Solver<Scalar,Ordinal>::addInterfaces()
{
  int err = addInterface("esi::Object", (void*)((esi::Object*)this) );
  if (err) {
    msgAbort("aztecoo_esi::Solver ctor addInterface(esi::Object) error");
  }

  err = addInterface("esi::Operator",
                     (void*)((esi::Operator<Scalar,Ordinal>*)this));
  if (err) {
    msgAbort("aztecoo_esi::Solver ctor addInterface(esi::Operator) error");
  }

  err = addInterface("esi::Solver",
                     (void*)((esi::Solver<Scalar,Ordinal>*)this));
  if (err) {
    msgAbort("aztecoo_esi::Solver ctor addInterface(esi::Solver) error");
  }

  err = addInterface("esi::SolverIterative",
                     (void*)((esi::SolverIterative<Scalar,Ordinal>*)this));
  if (err) {
    msgAbort("aztecoo_esi::Solver ctor addInterface(esi::SolverIterative) err");
  }

  err = addInterface("aztecoo_esi::Solver",
                     (void*)((aztecoo_esi::Solver<Scalar,Ordinal>*)this));
  if (err) {
    msgAbort("aztecoo_esi::Solver ctor addInterface(aztecoo_esi::Solver) err");
  }

  err = addInterface("AztecOO",
                     (void*)((AztecOO*)this));
  if (err) msgAbort("aztecoo_esi::Solver ctor addInterface(AztecOO) error");

  return(err);
}

//------------------------------------------------------------------------------
template<class Scalar, class Ordinal>
esi::ErrorCode aztecoo_esi::Solver<Scalar,Ordinal>::parameters(int numParams,
                                             char** paramStrings)
{
  //First, let's get those string-arrays in Epetra_Array objects...
  //(This is a very light-weight operation.)

  Epetra_Array<const char*>& azDefStrings = Translator::get_az_def_map();

  Epetra_Array<const char*>& azOptionStrings = Translator::get_az_option_strs();

  Epetra_Array<const char*>& azParamStrings = Translator::get_az_param_strs();

  //Now let's set up some work-strings...

  Epetra_Array<char> keyString(128), valString(128);
  char* keyStr = keyString.dataPtr();
  char* valStr = valString.dataPtr();

  //Now loop over the input parameters, and do the deciphering.
  for(int i=0; i<numParams; i++) {
    //
    //first see if this input parameter is a space-separated pair of strings.
    //if it isn't, we aren't interested in it.
    //
    int num = sscanf(paramStrings[i], "%s %s", keyStr, valStr);
    if (num < 2) continue;

    //Now we need to determine whether the key-string is an aztec-option
    //or an aztec-param. (Note the implicit assumption: it can't be both
    //an aztec-option AND an aztec-param. So, we jump out of this loop-
    //iteration if it's an aztec-option.)

    int optIndx = epetra_esi::findString(azOptionStrings, keyStr);
    if (optIndx >= 0) {
      handleAzOption(optIndx, keyStr, valStr, azDefStrings);
      continue;
    }

    int prmIndx = epetra_esi::findString(azParamStrings, keyStr);
    if (prmIndx >= 0) {
      handleAzParam(prmIndx, keyStr, valStr, azDefStrings);
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
template<class Scalar, class Ordinal>
int aztecoo_esi::Solver<Scalar,Ordinal>::handleAzOption(int optionIndx,
                                                     const char* keyStr,
                                                     const char* valStr,
                                         Epetra_Array<const char*>& azDefStrings
)
{
  static char tmpStr[128];

  int optionValue;
  int num = sscanf(valStr, "%d", &optionValue);

  if (num == 0) {
    int indx = epetra_esi::findHasSubString(azDefStrings, valStr);
    if (indx >= 0) {
      num = sscanf(azDefStrings[indx], "%s %d", tmpStr, &optionValue);
    }
  }

  if (num > 0) {
    aztecoo_->SetAztecOption(optionIndx, optionValue);
    if (optionIndx == AZ_max_iter) maxIters_ = optionValue;

    return(0);
  }

  cerr << "aztecoo_esi::Solver warning: handleAzOption("<<keyStr<<") failed to"
       << " convert '"<< valStr << "' to an integer." << endl;
  return(0);
}


//------------------------------------------------------------------------------
template<class Scalar, class Ordinal>
int aztecoo_esi::Solver<Scalar,Ordinal>::handleAzParam(int paramIndx,
                                                    const char* keyStr,
                                                    const char* valStr,
                                         Epetra_Array<const char*>& azDefStrings
)
{
  static char tmpStr[128];

  float paramValue;
  int num = sscanf(valStr, "%e", &paramValue);

  if (num == 0) {
    int indx = epetra_esi::findHasSubString(azDefStrings, valStr);
    if (indx >= 0) {
      num = sscanf(azDefStrings[indx], "%s %e", tmpStr, &paramValue);
    }
  }

  if (num > 0) {
    aztecoo_->SetAztecParam(paramIndx, (double)paramValue);
    if (paramIndx == AZ_tol) tolerance_ = paramValue;

    return(0);
  }

  cerr << "aztecoo_esi::Solver warning: handleAzParam("<<keyStr<<") failed to"
       << " convert '"<< valStr << "' to a double." << endl;
  return(0);
}

//------------------------------------------------------------------------------
template<class Scalar, class Ordinal>
int aztecoo_esi::Solver<Scalar,Ordinal>::createPetraVectorsThenSolve(
        esi::Vector<Scalar, Ordinal>& b, esi::Vector<Scalar, Ordinal>& x)
{
   int err, indexBase = 0;
   int globalSize, localSize;
   CHK_ERR( b.getGlobalSize(globalSize) );
   CHK_ERR( b.getLocalSize(localSize) );

   Epetra_Comm* comm = NULL;
#ifdef EPETRA_MPI
   //Since MPI is being used, let's try to get a communicator to create the
   //Epetra_Comm object with.
   MPI_Comm* mpicomm;
   err = getRunTimeModel("MPI", (void*&)mpicomm);
   if (err != 0) comm = new Epetra_SerialComm;
   else comm = new Epetra_MpiComm(*mpicomm);
#else
   comm = new Epetra_SerialComm;
#endif

   Epetra_Map pmap(globalSize, localSize, indexBase, *comm);

   double* xdata = NULL;
   double* bdata = NULL;
   CHK_ERR( x.getCoefPtrReadWriteLock(xdata) );
   CHK_ERR( b.getCoefPtrReadLock(bdata) );

   Epetra_Vector* px = new Epetra_Vector(Copy, pmap, xdata);
   Epetra_Vector* pb = new Epetra_Vector(Copy, pmap, bdata);
   if (px == NULL || pb == NULL) return(-1);

   err = aztecoo_->Iterate(petraA_, px, pb, maxIters_, tolerance_);

   delete px;
   delete pb;

   if (err != 0) return(err);

   CHK_ERR( x.releaseCoefPtrLock(xdata) );
   CHK_ERR( b.releaseCoefPtrLock(bdata) );

   return(err);
}

#endif

