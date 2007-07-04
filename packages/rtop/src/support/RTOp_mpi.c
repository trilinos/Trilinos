/*
// @HEADER
// ***********************************************************************
// 
//      Thyra: Interfaces and Support Code for the Interoperability of Abstract Numerical Algorithms 
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

/*
 * This file includes selected hollow MPI function definitions for a
 * sinlge process implementation.
 */

#include <assert.h>

/*
 * RAB: 2004/01/22: This file is included because it includes
 * Thyra_Config.h which then defines RTOp_USE_MPI or not.  If
 * RTOp_USE_MPI is defined then this header file will also include
 * RTOp_mpi.h for these delcarations.
 */
#include "RTOp_MPI_config.h"

#ifndef RTOp_USE_MPI

int MPI_Init(int *argc, char ***argv)
{
  return 0;
}

int MPI_Finalize(void)
{
  return 0;
}

int MPI_Comm_size(MPI_Comm comm, int *size)
{
  *size = 1;
  return 0;
}

int MPI_Comm_rank(MPI_Comm comm, int *rank)
{
  *rank = 0;
  return 0;
}

int MPI_Type_struct(int count , int *array_of_blocklengths, MPI_Aint *array_of_displacements
  , MPI_Datatype *array_of_types, MPI_Datatype *data_type)
{
  /* Make the mpi datatype just the extent (needed latter!) */
  int len = 0, extent = 0, k = 0;
  for( k = 0; k < count; ++k ) {
    switch( array_of_types[k] ) {
      case MPI_CHAR:
        len = sizeof(char);
        break;
      case MPI_INT:
        len = sizeof(int);
        break;
      case MPI_FLOAT:
        len = sizeof(float);
        break;
      case MPI_DOUBLE:
        len = sizeof(double);
        break;
      default:
        assert(0);
    }
    len = array_of_displacements[k] + array_of_blocklengths[k] * len;
    if( len > extent )
      extent = len;
  }
  *data_type = extent;
  return 0;
}

int MPI_Type_commit(MPI_Datatype *datatype)
{
  return 0;
}

int MPI_Type_free(MPI_Datatype *op)
{
  *op = MPI_DATATYPE_NULL;
  return 0;
}

int MPI_Op_create(MPI_User_function *func, int communitive, MPI_Op *op)
{
  *op = (MPI_Op)*func;
  return 0;
}
int MPI_Op_free( MPI_Op *op)
{
  *op = MPI_OP_NULL;
  return 0;
}

int MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
{
  assert(0); /* Should never be called in serial mode */
  return 0;
}
int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status)
{
  assert(0); /* Should never be called in serial mode */
  return 0;
}

int MPI_Sendrecv_replace(void* buff, int count, MPI_Datatype datatype, int dest, int sendtag, int source, int recvtag, MPI_Comm comm, MPI_Status* status)
{
  assert(0); /* Should never be called in serial mode */
  return 0;
}

int MPI_Reduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op
  , int root, MPI_Comm comm)
{
  char
    *_sendbuf = sendbuf,
    *_recvbuf = recvbuf;
  int k;
  for( k = 0; k < count * datatype; ++k )
    _recvbuf[k] =_sendbuf[k]; /* just copy bit for bit */
  return 0;
}

int MPI_Allreduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype
  , MPI_Op op, MPI_Comm comm)
{
  char
    *_sendbuf = sendbuf,
    *_recvbuf = recvbuf;
  int k;
  for( k = 0; k < count * datatype; ++k )
    _recvbuf[k] =_sendbuf[k]; /* just copy bit for bit */
  return 0;
}

int MPI_Barrier(MPI_Comm comm)
{
  return 0;
}

int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm )
{
  return 0;
}

int MPI_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype
         , void* recvbuf, int recvcount, MPI_Datatype recvtype, int root , MPI_Comm comm )
{
  char
    *_sendbuf = sendbuf,
    *_recvbuf = recvbuf;
  int k;
  assert(sendtype == recvtype);
  assert(sendcount == recvcount);
  for( k = 0; k < sendcount * sendtype; ++k )
    _recvbuf[k] =_sendbuf[k]; /* just copy bit for bit */
  return 0;
}

#endif /* RTOp_USE_MPI */
