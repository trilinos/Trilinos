/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/******************************************************************/
/* FILE  ***********    mpi_implementation.h   ********************/
/******************************************************************/
/* Author : Lisa Alano June 18 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/

#ifndef _MPI_IMPL_H
#define _MPI_IMPL_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _MPI_NO_ASSERTIONS
#include <stdio.h>
#define _MPI_Assert( condition ) if ( !(condition) ) { fprintf(stderr,"Assertion failure on line %d in file %s: %s\n", __LINE__, __FILE__, #condition); exit(77); }
#endif

#define _MPI_COVERAGE()

#define _MPI_STANDARD_MPI_CHECK(comm) if ( ! _MPI_CHECK_STATUS() ) { /* Comm->handler(....) */; return _MPI_NOT_OK; }

#define _MPI_DEBUG	 	(1)
#define _MPI_TRUE		(1)
#define _MPI_FALSE		(0)

#define _MPI_RANK 		(0)
#define _MPI_SIZE 		(1)

#define _MPI_OK		        (MPI_SUCCESS)
#define _MPI_NOT_OK		(-1)

#define _MPI_NOT_STARTED	(1000)
#define _MPI_STARTED		(1111)
#define _MPI_ENDED		(2222)

#define _MPI_VALID		(500)
#define _MPI_NOT_VALID		(-500)
#define _MPI_NULL		(-1)

#define _MPI_OP_START       	(MPI_MAX)	
#define _MPI_OP_END       	(MPI_MAXLOC)
#define _MPI_OP_OFFSET		(_MPI_OP_END+1)
#define _MPI_TYPE_START         (MPI_CHAR)
#define _MPI_TYPE_END           (MPI_PACKED)
#define _MPI_TYPE_OFFSET        (_MPI_TYPE_END+1)

#define _MPI_REQUEST_NULL	(-90)
#define _MPI_INDEXED		(-111)
#define _MPI_VECTOR             (-222)
#define _MPI_HINDEXED		(-333)
#define _MPI_HVECTOR		(-444)
#define _MPI_DEFAULT		(-555)
#define _MPI_CONTIG		(-666)
#define _MPI_STRUCT		(-777)

#define _MPI_COMM_WORLD_GROUP   ( (MPI_Group) MPI_COMM_WORLD)

/*============================================================================*/
#ifndef _MPI_PREALLOCATION_SIZE
#define _MPI_PREALLOCATION_SIZE 1000
#endif

/* "Complex" datatypes */
typedef struct _MPI_FLOAT_INT {
  float afloat;
  int   aint;
} _MPI_FLOAT_INT;

typedef struct _MPI_LONG_INT {
  long along;
  int aint;
} _MPI_LONG_INT;

typedef struct _MPI_DOUBLE_INT {
  double adouble;
  int aint;
} _MPI_DOUBLE_INT;

typedef struct _MPI_SHORT_INT {
  short ashort;
  int aint;
} _MPI_SHORT_INT;

typedef struct _MPI_2INT {
  int aint;
  int bint;
} _MPI_2INT;

typedef struct _MPI_LONG_DOUBLE_INT {
  long along;
  double adouble;
  int aint;
} _MPI_LONG_DOUBLE_INT;

typedef struct _MPI_LONG_LONG_INT {
  long along;
  long blong;
  int aint;
} _MPI_LONG_LONG_INT;


/* FORTRAN ==========================================
#define MPI_REAL                ( (MPI_Datatype) 1019)
#define MPI_INTEGER             ( (MPI_Datatype) 1020)
#define MPI_LOGICAL             ( (MPI_Datatype) 1021)
#define MPI_DOUBLE_PRECISION    ( (MPI_Datatype) 1022)
#define MPI_COMPLEX             ( (MPI_Datatype) 1023)
#define MPI_DOUBLE_COMPLEX      ( (MPI_Datatype) 1024)
#define MPI_INTEGER1            ( (MPI_Datatype) 1025)
#define MPI_INTEGER2            ( (MPI_Datatype) 1026)
#define MPI_INTEGER4            ( (MPI_Datatype) 1027)
#define MPI_REAL4               ( (MPI_Datatype) 1028)
#define MPI_REAL8               ( (MPI_Datatype) 1029)
#define MPI_2INTEGER            ( (MPI_Datatype) 1030)
#define MPI_2REAL               ( (MPI_Datatype) 1031)
#define MPI_2DOUBLE_PRECISION   ( (MPI_Datatype) 1032)
#define MPI_2COMPLEX            ( (MPI_Datatype) 1033)
#define MPI_2DOUBLE_COMPLEX     ( (MPI_Datatype) 1034)
 ================================================= */


typedef struct _MPI_COMM_IMPL {
  MPI_Comm comm;
  int valid;
  MPI_Handler_function_pointer handler;
  /*MPI_Handler_function handler;*/
  MPI_Group group; 
} _MPI_COMM_IMPL;

typedef struct _MPI_TYPE_DES {
  int id;
  MPI_Aint size;
  MPI_Aint extent;
  MPI_Aint ub;
  MPI_Aint lb;
  int sendType;
  struct _MPI_TYPE_DES* next;
  struct _MPI_TYPE_INFO* info;
} _MPI_TYPE_DES;

typedef struct _MPI_TYPE_INFO {
  int count;
  int* blocklen;
  int* stride;
  int* types;
} _MPI_TYPE_INFO;

typedef struct _MPI_DATA_ENTRY {
  int valid;
  void* buffer;
  int count;
  MPI_Datatype type;
  int tag;
  MPI_Comm comm; 
  int user;      /*Lets you know if you should free the buffer. No, if it is the user's buffer*/
} _MPI_DATA_ENTRY;

typedef struct _MPI_OP_TYPE {
  int valid;
  void (*function) ( void * a, void * b, int * len, MPI_Datatype * ); 
  int commute;
} _MPI_OP_TYPE;


/* Variables */

extern _MPI_COMM_IMPL* _MPI_COMM_LIST;
extern _MPI_DATA_ENTRY* _MPI_DATA_BUFF;
extern _MPI_TYPE_DES* _MPI_TYPE_LIST;
extern _MPI_OP_TYPE* _MPI_OP_LIST;
extern _MPI_REQUEST_OBJECT** _MPI_REQ_LIST_OF_LISTS;

extern int _MPI_INIT_STATUS;
extern int _MPI_FINALIZED_FLAG;
extern int _MPI_INITIALIZED_FLAG;
extern int _MPI_TYPE_COUNT;
extern int _MPI_DATA_BUFF_COUNT;
extern int _MPI_COMM_COUNT;
extern int _MPI_OP_COUNT;
extern int _MPI_REQ_COUNT;

extern int _MPI_COMM_ARRAY_SIZE;
extern int _MPI_DATA_ARRAY_SIZE;
extern int _MPI_TYPE_ARRAY_SIZE;
extern int _MPI_OP_ARRAY_SIZE;
extern int _MPI_REQ_ARRAY_SIZE;

extern MPI_Request* _MPI_REQNULL;

/* Prototypes */

int _MPI_ERR_ROUTINE (MPI_Error_Class error, char* message);
int _MPI_CHECK_STATUS (MPI_Comm* comm);

/* MPI_Comm_Create.c  and MPI_Comm_Init.c */
int _MPI_Comm_Insert (int index);
int _MPI_Comm_Insert0 (int index, MPI_Comm set_comm, MPI_Handler_function_pointer function);
int _MPI_Comm_Invalid (int index);
int _MPI_Comm_check (int index);
int _MPI_Comm_check_legal (MPI_Comm comm, int *index);
int _MPI_Group_check (MPI_Group group);

/* MPI_Send and MPI_Recv Helper Functions - UTILITY*/
int _MPI_calculateSize(int count, MPI_Datatype type);
int _MPI_getSize(MPI_Datatype type);
int _MPI_calculateStructureSize (MPI_Datatype type);
int _MPI_checkCommunicator(MPI_Comm comm);
int _MPI_checks (void *message, int count, MPI_Datatype data, int dest, int tag, MPI_Comm comm);
int _MPI_checkSendType (MPI_Datatype type);
int _MPI_checkTag(int tag);
int _MPI_checkDestination(int dest);
int _MPI_checkDatatype(MPI_Datatype type);
int _MPI_checkCount(int count);
int _MPI_checkBuffer(void* message);
int _MPI_checkRequest (MPI_Request request);
int _MPI_COPY_UTIL(void *sendbuf,int sendcount, MPI_Datatype sendtype,void *recvbuf,int recvcount, MPI_Datatype recvtype);

/* Internal Data buffer */
int _MPI_Data_Invalid(int index);
int _MPI_Buff_Insert(void *message,int count,MPI_Datatype datatype,int tag,MPI_Comm comm);
int _MPI_Bbuff_Insert(void *message,int count,MPI_Datatype datatype,int tag,MPI_Comm comm);
int _MPI_Buff_Find(int tag, MPI_Comm comm);
int _MPI_Buff_Ifind(int tag, MPI_Comm comm);

int _MPI_Type_Invalid(int index);

/* Internal Operations */
int _MPI_checkOp(MPI_Op op);
void _MPI_Default_Op ( void* a, void* b, int* len, MPI_Datatype* type );
int _MPI_Op_insert (MPI_User_function function, int commute, MPI_Op* id);
int _MPI_Op_invalid(int index);
int _MPI_Function_check(MPI_User_function* function);
int _MPI_Commute_check(int commute);
int _MPI_check_overlap(void* sendbuf, void* recvbuf, int size);

/* Internal Requests Posted */
MPI_Request _MPI_New_Request(void* message, int count, MPI_Datatype datatype, int tag, MPI_Comm comm,int send);
int _MPI_Req_Invalid (MPI_Request request);
int _MPI_Check_Request_Array(int count, MPI_Request array[]);


/* Type operations */
int _MPI_BasicType (MPI_Datatype datatype);
int _MPI_FindType (MPI_Datatype datatype);
int _MPI_Free_datatype (MPI_Datatype datatype);
void _MPI_deleteAll (_MPI_TYPE_DES* parent);
int _MPI_Find_free (void);


#include "_MPI_UTILITY.h"

#ifdef __cplusplus
}
#endif

#endif
