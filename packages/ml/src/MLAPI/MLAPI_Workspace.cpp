/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_common.h"
#if defined(HAVE_ML_MLAPI)
#include "Epetra_Comm.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "MLAPI_Error.h"
#include "MLAPI_Workspace.h"
#ifdef _MSC_VER
#include "winprocess.h"
#endif

namespace MLAPI {

static ML_Comm* ML_Comm_ = 0;
static Epetra_Comm* Epetra_Comm_ = 0;

// ====================================================================== 
ML_Comm* GetML_Comm() 
{
  if (ML_Comm_ == 0) Init();

  return(ML_Comm_);
}

// ====================================================================== 
Epetra_Comm& GetEpetra_Comm() 
{
  if (Epetra_Comm_ == 0) Init();

  return(*Epetra_Comm_);
}

// ====================================================================== 
void Barrier()
{
  if (Epetra_Comm_ == 0) Init();

  Epetra_Comm_->Barrier();
}

// ====================================================================== 
int GetMyPID()
{
  if (Epetra_Comm_ == 0) Init();

  return(Epetra_Comm_->MyPID());
}

// ====================================================================== 
int GetNumProcs()
{
  if (Epetra_Comm_ == 0) Init();

  return(Epetra_Comm_->NumProc());
}

// ====================================================================== 
int GetPrintLevel() 
{
  if (GetMyPID())
    return(0);
  else
    return(ML_Get_PrintLevel());
}
  
// ====================================================================== 
void SetPrintLevel(int Level) 
{
  ML_Set_PrintLevel(Level);
}

// ====================================================================== 
void Init() 
{
  if (ML_Comm_ == 0) ML_Comm_Create(&ML_Comm_);
  if (Epetra_Comm_ == 0) {
#ifdef HAVE_MPI
    Epetra_Comm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
    Epetra_Comm_ = new Epetra_SerialComm;
#endif
  }

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

  GetEpetra_Comm().SumAll(&i, &j, 1);

  if (j != 0)
  {
    if (GetMyPID()  == 0) cout << "Host and Process Ids for tasks" << endl;
    for (i = 0; i < GetNumProcs() ; i++) {
      if (i == GetMyPID() ) {
#if defined(TFLOP) || defined(JANUS_STLPORT) || defined(COUGAR)
        sprintf(buf, "Host: %s   PID: %d", "janus", getpid());
#else
        gethostname(hostname, sizeof(hostname));
        sprintf(buf, "Host: %s\tMyPID(): %d\tPID: %d",
                hostname, GetMyPID(), getpid());
#endif
        printf("%s\n",buf);
        fflush(stdout);
#ifdef ICL
        Sleep(1);
#else
        sleep(1);
#endif
      }
    }
    if (GetMyPID() == 0) {
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

  ML_Set_PrintLevel(10);
}

// ====================================================================== 
void Finalize() 
{
  if (ML_Comm_) {
    ML_Comm_Destroy(&ML_Comm_);
    ML_Comm_ = 0;
  }

  if (Epetra_Comm_) {
    delete Epetra_Comm_;
    Epetra_Comm_ = 0;
  }
}

// ====================================================================== 
string GetString(const int& x) 
{
  char s[100];
  sprintf(s, "%d", x);
  return string(s);
}

// ====================================================================== 
string GetString(const double& x)
{
  char s[100];
  sprintf(s, "%g", x);
  return string(s);
}

// ====================================================================== 
int GetMatrixType() 
{
  return(ML_CSR_MATRIX);
}

} // namespace MLAPI
#endif
