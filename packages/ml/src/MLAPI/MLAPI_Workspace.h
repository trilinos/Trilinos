#ifndef MLAPI_WORKSPACE_H
#define MLAPI_WORKSPACE_H

#include "ml_include.h"
#include "Epetra_Comm.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#define ML_THROW(a) \
  { cerr << "*ML*ERR* (" << __FILE__ << "+:" << __LINE__ \
         << ") " << endl; \
    cerr << "*ML*ERR* " << a << endl; \
    throw(a); \
  }

namespace MLAPI {

static ML_Comm* MLAPIComm_ = 0;
static Epetra_Comm* EpetraComm_ = 0;

void Init() {
  ML_Comm_Create(&MLAPIComm_);
#ifdef HAVE_MPI
  EpetraComm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  EpetraComm_ = new Epetra_SerialComm;
#endif
  
}

void Finalize() {
  ML_Comm_Destroy(&MLAPIComm_);
}

ML_Comm* GetMLComm() {
  return(MLAPIComm_);
}

Epetra_Comm& GetEpetraComm() {
  return(*EpetraComm_);
}

int MyPID()
{
  return(EpetraComm_->MyPID());
}

int NumProc()
{
  return(EpetraComm_->NumProc());
}

int PrintLevel() 
{
  if (MyPID())
    return(0);
  else
    return(ML_Get_PrintLevel());
}
  
void SetPrintLevel(int Level) 
{
  ML_Set_PrintLevel(Level);
}

} // namespace MLAPI

#endif
