/*
// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

/* */
/* MPI declarations used by RTOp example program. */
/* These where taken from mpich for Windows NT. */
/* */

#ifndef RTOP_MPI_H
#define RTOP_MPI_H

#include "RTOp_ConfigDefs.hpp" /* This C++ file has a proper C mode */

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////// */
/* MPI declarations */

#define MPI_Aint int
typedef int MPI_Datatype;
#define MPI_CHAR           ((MPI_Datatype)1)
#define MPI_INT            ((MPI_Datatype)6)
#define MPI_FLOAT          ((MPI_Datatype)10)
#define MPI_DOUBLE         ((MPI_Datatype)11)
typedef int MPI_Comm;
#define MPI_COMM_WORLD 91
#define MPI_COMM_NULL      ((MPI_Comm)0)
typedef int MPI_Op;
#define MPI_OP_NULL        ((MPI_Op)0)
#define MPI_MAX            (MPI_Op)(100)
#define MPI_MIN            (MPI_Op)(101)
#define MPI_SUM            (MPI_Op)(102)
#define MPI_DATATYPE_NULL  ((MPI_Datatype)0)
typedef struct { int MPI_SOURCE; int MPI_TAG; int MPI_ERROR; } MPI_Status;
typedef void (MPI_User_function) ( void *, void *, int *, MPI_Datatype * ); 

/* ////////////////////////// */
/* MPI functions */

#define EXPORT_MPI_API
EXPORT_MPI_API int MPI_Init(int *, char ***);
EXPORT_MPI_API int MPI_Finalize(void);
EXPORT_MPI_API int MPI_Comm_size(MPI_Comm, int *);
EXPORT_MPI_API int MPI_Comm_rank(MPI_Comm, int *);
EXPORT_MPI_API int MPI_Type_struct(int, int *, MPI_Aint *, MPI_Datatype *, MPI_Datatype *);
EXPORT_MPI_API int MPI_Type_commit(MPI_Datatype *);
EXPORT_MPI_API int MPI_Type_free(MPI_Datatype *);
EXPORT_MPI_API int MPI_Op_create(MPI_User_function *, int, MPI_Op *);
EXPORT_MPI_API int MPI_Op_free( MPI_Op *);
EXPORT_MPI_API int MPI_Send(void*, int, MPI_Datatype, int, int, MPI_Comm);
EXPORT_MPI_API int MPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status*);
EXPORT_MPI_API int MPI_Sendrecv_replace(void*, int, MPI_Datatype, int, int, int, int, MPI_Comm, MPI_Status*);
EXPORT_MPI_API int MPI_Reduce(void* , void*, int, MPI_Datatype, MPI_Op, int, MPI_Comm);
EXPORT_MPI_API int MPI_Allreduce(void* , void*, int, MPI_Datatype, MPI_Op, MPI_Comm);
EXPORT_MPI_API int MPI_Barrier(MPI_Comm);
EXPORT_MPI_API int MPI_Bcast(void*, int, MPI_Datatype, int, MPI_Comm );
EXPORT_MPI_API int MPI_Gather(void* , int, MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm); 

#ifdef __cplusplus
}
#endif

#endif /* RTOP_MPI_H */
