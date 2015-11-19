/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: rrdrake $
 *    Date: 2009/07/17 15:14:49 $
 *    Revision: 1.4 $
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/******************************************************************/
/* FILE  ***********        mpi.h              ********************/
/******************************************************************/
/* Author : Lisa Alano June 18 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/

#ifndef _MPI_H
#define _MPI_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

#define MPI_VERSION             1
#define MPI_SUBVERSION          2

typedef int MPI_Datatype;

#define MPI_CHAR		( (MPI_Datatype) 1001)
#define MPI_BYTE		( (MPI_Datatype) 1002)
#define MPI_SHORT		( (MPI_Datatype) 1003)
#define MPI_INT			( (MPI_Datatype) 1004)
#define MPI_LONG		( (MPI_Datatype) 1005)
#define MPI_FLOAT		( (MPI_Datatype) 1006)
#define MPI_DOUBLE		( (MPI_Datatype) 1007)
#define MPI_UNSIGNED_CHAR	( (MPI_Datatype) 1008)
#define MPI_UNSIGNED_SHORT	( (MPI_Datatype) 1009)
#define MPI_UNSIGNED		( (MPI_Datatype) 1010)
#define MPI_UNSIGNED_LONG	( (MPI_Datatype) 1011)
#define MPI_LONG_DOUBLE		( (MPI_Datatype) 1012)
#define MPI_LONG_LONG           ( (MPI_Datatype) 1013)  /* KDD Added 11/19/15.*/
                                                        /* KDD Incremented all*/
                                                        /* KDD following types*/
                                                        /* KDD as code assumes*/
                                                        /* KDD MPI_PACKED is  */
                                                        /* KDD the last type  */
                                                        /* KDD specified.     */

/* Dataypes for the MPI functions MPI_MAXLOC and MPI_MINLOC */

#define MPI_FLOAT_INT		( (MPI_Datatype) 1014)
#define MPI_LONG_INT		( (MPI_Datatype) 1015)
#define MPI_DOUBLE_INT		( (MPI_Datatype) 1016)
#define MPI_SHORT_INT		( (MPI_Datatype) 1017)
#define MPI_2INT		( (MPI_Datatype) 1018)
#define MPI_LONG_DOUBLE_INT	( (MPI_Datatype) 1019)
#define MPI_LONG_LONG_INT	( (MPI_Datatype) 1020)

/* Special Dataypes */

#define MPI_REAL		( (MPI_Datatype) 1021)
#define MPI_INTEGER		( (MPI_Datatype) 1022)
#define MPI_LOGICAL		( (MPI_Datatype) 1023)
#define MPI_DOUBLE_PRECISION	( (MPI_Datatype) 1024)
#define MPI_COMPLEX		( (MPI_Datatype) 1025)
#define MPI_DOUBLE_COMPLEX	( (MPI_Datatype) 1026)

/* Optional Datatypes */

#define MPI_INTEGER1		( (MPI_Datatype) 1027)
#define MPI_INTEGER2		( (MPI_Datatype) 1028)
#define MPI_INTEGER4		( (MPI_Datatype) 1029)
#define MPI_REAL4		( (MPI_Datatype) 1030)
#define MPI_REAL8		( (MPI_Datatype) 1031)

/* Datatypes for the MPI functions MPI_MAXLOC and MPI_MINLOC */

#define MPI_2INTEGER		( (MPI_Datatype) 1032)
#define MPI_2REAL		( (MPI_Datatype) 1033)
#define MPI_2DOUBLE_PRECISION	( (MPI_Datatype) 1034)
#define MPI_2COMPLEX		( (MPI_Datatype) 1035)
#define MPI_2DOUBLE_COMPLEX	( (MPI_Datatype) 1036)

#define MPI_UB			( (MPI_Datatype) 1037)
#define MPI_LB			( (MPI_Datatype) 1038)
#define MPI_PACKED		( (MPI_Datatype) 1039)


typedef int MPI_Group;


#define MPI_GROUP_EMPTY		( (MPI_Group) 0)
#define MPI_GROUP_WORLD		( (MPI_Group) 1)

/* Results of the compare Operations */

#define MPI_IDENT		(111)
#define MPI_CONGRUENT		(222)
#define MPI_SIMILAR		(333)
#define MPI_UNEQUAL		(444)

/* Collective Operations */

typedef int MPI_Op;

#define MPI_MAX			( (MPI_Op) 49)
#define MPI_MIN			( (MPI_Op) 50)
#define MPI_SUM			( (MPI_Op) 51)
#define MPI_PROD		( (MPI_Op) 52)
#define MPI_LAND		( (MPI_Op) 53)
#define MPI_BAND		( (MPI_Op) 54)
#define MPI_LOR			( (MPI_Op) 55)
#define MPI_BOR			( (MPI_Op) 56)
#define MPI_LXOR		( (MPI_Op) 57)
#define MPI_BXOR		( (MPI_Op) 57)
#define MPI_MINLOC		( (MPI_Op) 58)
#define MPI_MAXLOC		( (MPI_Op) 59)

/* Permanent Key Vaues */

#define MPI_TAG_UB		(100)
#define MPI_HOST		(0)
#define MPI_IO			(0)
#define MPI_WTIME_IS_GLOBAL	(0)

/* Null Objects */

#define MPI_COMM_NULL	      	(-1)
#define MPI_OP_NULL		(-1)
#define MPI_GROUP_NULL		(-1)
#define MPI_DATATYPE_NULL	(-1)
#define MPI_ERRHANDLER_NULL	(-1)

/* Predefined Constants */

#define MPI_MAX_PROCESSOR_NAME	(30)
#define MPI_MAX_ERROR_STRING	(300)
#define MPI_UNDEFINED		(-100)
#define MPI_UNDEFINED_RANK	(-50)
#define MPI_KEYVAL_INVALID	(-150)
#define MPI_BSEND_OVERHEAD	(5)
#define MPI_PROC_NULL		(-200)
#define MPI_ANY_SOURCE		(-300)
#define MPI_ANY_TAG		(-400)
#define MPI_BOTTOM		(-500)

/* Toplogy Types */

#define MPI_GRAPH		(6000)
#define MPI_CART		(6001)

typedef int MPI_Comm;
 
typedef struct MPI_Status {
  int MPI_SOURCE;
  int MPI_TAG;
  int MPI_ERROR;
  int __count;
} MPI_Status;

#define MPI_STATUS_IGNORE ((MPI_Status *)0)
#define MPI_STATUSES_IGNORE ((MPI_Status *)0)

typedef struct _MPI_REQUEST_OBJECT {
  void* buffer; 
  int count; 
  int type;
  int tag;
  MPI_Comm comm;
  int send;
  int valid;
  int cancel;  /* true if this request was cancelled */
} _MPI_REQUEST_OBJECT;

typedef _MPI_REQUEST_OBJECT *MPI_Request;

/*MPI_Request MPI_REQUEST_NULL;*/
#define MPI_REQUEST_NULL	( (MPI_Request)NULL)

/* Special MPI types and functions */

/*
MPI_NULL_COPY_FN
MPI_NULL_DELETE_FN
MPI_DUP_FN
*/

#include "mpi_config.h"

typedef size_t MPI_Aint;
typedef int MPIO_Request;
typedef int MPI_Errhandler;
typedef int MPI_Fint;
typedef int MPI_File;
typedef int MPI_Info;
typedef int MPI_Offset;
typedef int MPI_Handle_enum;
typedef int MPI_Handle_type;

#define MPI_FILE_NULL ((MPI_File)(-1))
#define MPI_INFO_NULL ((MPI_Info)(-1))

#define MPI_MODE_RDONLY              2
#define MPI_MODE_RDWR                8
#define MPI_MODE_WRONLY              4
#define MPI_MODE_CREATE              1
#define MPI_MODE_EXCL               64
#define MPI_MODE_DELETE_ON_CLOSE    16
#define MPI_MODE_UNIQUE_OPEN        32
#define MPI_MODE_APPEND            128
#define MPI_MODE_SEQUENTIAL        256

#define MPI_MAX_INFO_KEY       255
#define MPI_MAX_INFO_VAL      1024

typedef int MPI_Error_Class;

typedef void MPI_Handler_function (MPI_Comm* comm, int* error_code, ...);
typedef void (*MPI_Handler_function_pointer) (MPI_Comm* comm, int* error_code, ...);

void MPI_ERRORS_ARE_FATAL (MPI_Comm* comm, int* error_code, ...);
void MPI_ERRORS_RETURN (MPI_Comm* comm, int* error_code, ...);
typedef void MPI_User_function(void* invec, void* inoutvec, int* len, MPI_Datatype* dataype);

#include "mpi_implementation.h"


#define MPI_COMM_WORLD ( (MPI_Comm) -100)		
#define MPI_COMM_SELF	( (MPI_Comm) -100)

typedef int MPI_Copy_function (MPI_Comm oldcomm, int keyval, 
	void* extra_arg, void* attribute_val_in, void* attribute_val_out, int* flag);
typedef int MPI_Delete_function (MPI_Comm comm, int keyval, void* attribute_val, void* extra_arg);



 
/* MPI Error Classes */

/* Success=0 is safest as some codes wrongly assume zero=good non-zero=bad */
#define MPI_SUCCESS 		( (MPI_Error_Class) 0) 
#define MPI_ERR_BUFFER		( (MPI_Error_Class) -5001)
#define MPI_ERR_COUNT		( (MPI_Error_Class) -5002)
#define MPI_ERR_TYPE		( (MPI_Error_Class) -5003)
#define MPI_ERR_TAG		( (MPI_Error_Class) -5004)
#define MPI_ERR_COMM		( (MPI_Error_Class) -5005)
#define MPI_ERR_RANK		( (MPI_Error_Class) -5006)
#define MPI_ERR_ROOT		( (MPI_Error_Class) -5007)
#define MPI_ERR_GROUP		( (MPI_Error_Class) -5008)
#define MPI_ERR_OP		( (MPI_Error_Class) -5009)
#define MPI_ERR_TOPOLOGY	( (MPI_Error_Class) -5010)
#define MPI_ERR_DIMS		( (MPI_Error_Class) -5011)
#define MPI_ERR_ARG		( (MPI_Error_Class) -5012)
#define MPI_ERR_UNKNOWN		( (MPI_Error_Class) -5013)
#define MPI_ERR_TRUNCATE	( (MPI_Error_Class) -5014)
#define MPI_ERR_OTHER		( (MPI_Error_Class) -5015)
#define MPI_ERR_IN_STATUS	( (MPI_Error_Class) -5016)
#define MPI_ERR_PENDING		( (MPI_Error_Class) -5017)
#define MPI_ERR_REQUEST		( (MPI_Error_Class) -5018)
#define MPI_ERR_LASTCODE	( (MPI_Error_Class) -5019)
#define MPI_ERR_INTERN		( (MPI_Error_Class) -5020)


/* Function Prototypes */

#include "mpi_prototypes.h"

#include "mpi_profile.h"

/* Additional variables */



#if 0
FORTRAN: DOUBLE PRECISION MPI_WTIME
FORTRAN: EXTERNAL MPI_WTIME

FORTRAN: DOUBLE PRECISION MPI_WTICK
FORTRAN: EXTERNAL MPI_WTICK

FORTRAN: DOUBLE PRECISION PMPI_WTIME
FORTRAN: EXTERNAL PMPI_WTIME

FORTRAN: DOUBLE PRECISION PMPI_WTICK
FORTRAN: EXTERNAL PMPI_WTICK

SKIP FOR NOW FORTRAN: EXTERNAL MPI_NULL_COPY_FN
SKIP FOR NOW FORTRAN: EXTERNAL MPI_NULL_DELETE_FN
SKIP FOR NOW FORTRAN: EXTERNAL MPI_DUP_FN
#endif

	
#ifdef __cplusplus
}
#endif

#endif
