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

//! Returns a pointer to the ML_Comm object defined on MPI_COMM_WORLD.
ML_Comm* GetMLComm() {
  return(MLAPIComm_);
}

//! Returns a reference to the Epetra_Comm object defined on MPI_COMM_WORLD.
Epetra_Comm& GetEpetraComm() {
  return(*EpetraComm_);
}

//! Returns the ID of the calling process.
int MyPID()
{
  return(EpetraComm_->MyPID());
}

//! Returns the total number of processes in the computation.
int NumProc()
{
  return(EpetraComm_->NumProc());
}

//! Retutns the level of output (always 0 if MyPID() != 0).
int PrintLevel() 
{
  if (MyPID())
    return(0);
  else
    return(ML_Get_PrintLevel());
}
  
//! Sets the level of output (from 0 to 10, 0 being verbose).
void SetPrintLevel(int Level) 
{
  ML_Set_PrintLevel(Level);
}

//! Initialize the MLAPI workspace.
void Init() {
  ML_Comm_Create(&MLAPIComm_);
#ifdef HAVE_MPI
  EpetraComm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  EpetraComm_ = new Epetra_SerialComm;
#endif

    char * str = (char *) getenv("ML_BREAK_FOR_DEBUGGER");
  int i = 0, j = 0;
  char buf[80];
  char go = ' ';
  char hostname[80];
  if (str != NULL) i++;

  FILE * ML_capture_flag;
  ML_capture_flag = fopen("ML_debug_now","r");
  if(ML_capture_flag) {
    i++;
    fclose(ML_capture_flag);
  }

  GetEpetraComm().SumAll(&i, &j, 1);

  if (j != 0)
  {
    if (MyPID()  == 0) cout << "Host and Process Ids for tasks" << endl;
    for (i = 0; i <NumProc() ; i++) {
      if (i == MyPID() ) {
#ifdef COUGAR
        sprintf(buf, "Host: %s   PID: %d", "janus", getpid());
#else
        gethostname(hostname, sizeof(hostname));
        sprintf(buf, "Host: %s\tMyPID(): %d\tPID: %d",
                hostname, MyPID(), getpid());
#endif
        printf("%s\n",buf);
        fflush(stdout);
        sleep(1);
      }
    }
    if (MyPID() == 0) {
      printf("\n");
      printf("** Pausing because environment variable ML_BREAK_FOR_DEBUGGER has been set,\n");
      puts("** or file ML_debug_now has been created");
      printf("**\n");
      printf("** You may now attach debugger to the processes listed above.\n");
      printf( "**\n");
      printf( "** Enter a character to continue > "); fflush(stdout);
      scanf("%c",&go);
    }
  }

}

//! Destroys the MLAPI workspace.
void Finalize() {
  ML_Comm_Destroy(&MLAPIComm_);
}

} // namespace MLAPI

#endif
