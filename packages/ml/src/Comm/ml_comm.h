/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the ML communicator structure                         */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : September, 1998                                      */
/* ******************************************************************** */

#ifndef _MLCOMM_
#define _MLCOMM_

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include <stdio.h>
#include "ml_common.h"
#include "ml_defs.h"
#include "ml_memory.h"

/* ******************************************************************** */
/* communicator can be loaded by users or use the resident mpi          */
/* -------------------------------------------------------------------- */

#ifdef ML_MPI
#include <mpi.h>
#define USR_COMM MPI_Comm
#define USR_COMM_WORLD MPI_COMM_WORLD // CHECK: ALLOW MPI_COMM_WORLD
#define USR_REQ  MPI_Request

#define USR_ERRHANDLER MPI_Errhandler

/* two typedefs that are useful in MPI reduction calls */
typedef struct ml_IntLoc_struct {
  int value;
  int rank;
} ml_IntLoc;

typedef struct ml_DblLoc_struct {
  double value;
  int rank;
} ml_DblLoc;

#ifdef ML_CATCH_MPI_ERRORS_IN_DEBUGGER

/* the following might be necessary to avoid compiler warnings, since
   in newer versions of MPI (e.g., recent lam implementations), certain
   error handling functions may be deprecated. */
/*
#define USR_ERRHANDLER_FUNCTION MPI_Comm_errhandler_fn
#define USR_ERRHANDLER_SET      MPI_Comm_set_errhandler
#define USR_ERRHANDLER_CREATE   MPI_Comm_create_errhandler
*/
#define USR_ERRHANDLER_FUNCTION MPI_Handler_function
#define USR_ERRHANDLER_SET      MPI_Errhandler_set
#define USR_ERRHANDLER_CREATE   MPI_Errhandler_create

#endif /* ifdef ML_CATCH_MPI_ERRORS_IN_DEBUGGER */

#else /* ifdef ML_MPI */

#define USR_COMM int
#define USR_COMM_WORLD 0
#define USR_REQ  int

#define USR_ERRHANDLER int

#ifdef ML_CATCH_MPI_ERRORS_IN_DEBUGGER
#define USR_ERRHANDLER_FUNCTION void
#endif /* ifdef ML_CATCH_MPI_ERRORS_IN_DEBUGGER */

#endif /* ifdef ML_MPI */

/* ******************************************************************** */
/* components of the ML communicator                                    */
/* -------------------------------------------------------------------- */

typedef struct ML_Comm_Struct
{
   int      ML_id;
   int      ML_mypid;
   int      ML_nprocs;
   USR_COMM USR_comm;
   int      (*USR_sendbytes)(void*,unsigned int,int,int,USR_COMM);
   int      (*USR_irecvbytes)(void*,unsigned int,int*,int*,USR_COMM,USR_REQ*);
   int      (*USR_waitbytes)(void*,unsigned int,int*,int*,USR_COMM,USR_REQ*);
   void     (*USR_cheapwaitbytes)(void*,unsigned int,int*,int*,USR_COMM,USR_REQ*);
   USR_ERRHANDLER USR_errhandler;

} ML_Comm;
#if defined(__cplusplus) && defined(_MSC_VER)
extern "C" ML_Comm *global_comm; /* should be made obsolete */
#else
extern ML_Comm *global_comm; /* should be made obsolete */
#endif
/* ******************************************************************** */
/* functions for the ML communicator                                    */
/* -------------------------------------------------------------------- */

#ifndef ML_CPP
#ifdef __cplusplus
extern "C"
{
#endif
#endif

extern int    ML_Comm_Create( ML_Comm ** comm );
extern int    ML_Comm_Create2( ML_Comm ** comm, USR_COMM com );
extern int    ML_Comm_Destroy( ML_Comm ** comm );
extern int    ML_Comm_Check( ML_Comm *comm );

extern int    ML_Comm_Set_UsrComm( ML_Comm *comm, USR_COMM com );
extern int    ML_Comm_Set_Mypid( ML_Comm *comm, int mypid );
extern int    ML_Comm_Set_Nprocs( ML_Comm *comm, int nprocs);
extern int    ML_Comm_Set_SendFcn( ML_Comm *comm, int (*SendFcn)(void*,unsigned int,int,int,USR_COMM));
extern int    ML_Comm_Set_RecvFcn( ML_Comm *comm, int (*RecvFcn)(void*,unsigned int,int*,int*,USR_COMM,USR_REQ*));
extern int    ML_Comm_Set_WaitFcn( ML_Comm *comm, int (*WaitFcn)(void*,unsigned int,int*,int*,USR_COMM,USR_REQ*));

extern int    ML_Comm_GmaxInt( ML_Comm *comm, int intdata );
extern double ML_Comm_GmaxDouble(ML_Comm *comm, double ddata );
extern int    ML_Comm_GsumInt( ML_Comm *comm, int intdata );
extern double ML_Comm_GsumDouble(ML_Comm *comm, double ddata );
extern int    ML_Comm_GappendInt(ML_Comm *,int *invec,int *lleng, int tleng);
extern int    ML_Comm_GappendBigInt(ML_Comm *,ml_big_int *invec,int *lleng, int tleng);
extern int    ML_Comm_GappendDouble(ML_Comm*,double *dvec,int *lleng,int tleng);

extern int    ML_Comm_Irecv(void*,unsigned int,int *,int *,USR_COMM,USR_REQ*);
extern int    ML_Comm_Wait (void*,unsigned int,int *,int *,USR_COMM,USR_REQ*);
extern void   ML_Comm_CheapWait (void*,unsigned int,int *,int *,USR_COMM,USR_REQ*);
extern int    ML_Comm_Send (void*,unsigned int,int,  int,  USR_COMM );
extern int    ML_Comm_Barrier(ML_Comm *com_ptr);
extern int ML_gpartialsum_int(int val, ML_Comm *comm);


#ifdef ML_CATCH_MPI_ERRORS_IN_DEBUGGER
extern int    ML_Comm_ErrorHandlerSet(USR_COMM, USR_ERRHANDLER);
/*
extern int    ML_Comm_ErrorHandlerCreate(
                 void *(*HandlerFcn)(USR_COMM*,int*),
                 USR_ERRHANDLER*);
                 exit(1);
*/
extern int    ML_Comm_ErrorHandlerCreate(
                 USR_ERRHANDLER_FUNCTION (*HandlerFcn),
                 USR_ERRHANDLER*);
extern int    ML_Comm_ErrorHandlerDestroy(USR_ERRHANDLER*);
extern void   ML_Comm_ErrorHandler(USR_COMM*, int*,...);
#endif

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
