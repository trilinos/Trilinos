
//@HEADER
/*
************************************************************************

              Epetra: Linear Algebra Services Package 
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

#define EPETRA_FORTRAN
#ifdef EPETRA_ADDRESS64BIT
#define EPETRA_OBJECT_PTR long integer
#else
#define EPETRA_OBJECT_PTR integer
#endif
#define EPETRA_ERROR_FLAG integer

      EPETRA_OBJECT_PTR       epetra_mpicomm_create1
      EPETRA_OBJECT_PTR       epetra_mpicomm_create2
      EPETRA_OBJECT_PTR       epetra_serialcomm_create

      EPETRA_ERROR_FLAG       epetra_comm_mypid
      EPETRA_ERROR_FLAG       epetra_comm_numproc

      EPETRA_OBJECT_PTR       epetra_map_create1
      EPETRA_OBJECT_PTR       epetra_map_create2
      EPETRA_OBJECT_PTR       epetra_map_create3
      EPETRA_ERROR_FLAG       epetra_map_numMyelements
      EPETRA_ERROR_FLAG       epetra_map_numGlobalelements

      EPETRA_OBJECT_PTR       epetra_vector_create1
      EPETRA_OBJECT_PTR       epetra_vector_create2
      EPETRA_ERROR_FLAG       epetra_vector_update
      EPETRA_ERROR_FLAG       epetra_vector_norm1
      EPETRA_ERROR_FLAG       epetra_vector_norm2
      EPETRA_ERROR_FLAG       epetra_vector_random
      EPETRA_ERROR_FLAG       epetra_vector_putscalar


      EPETRA_OBJECT_PTR       epetra_localmap_create
      EPETRA_OBJECT_PTR       epetra_blockmap_create1
      EPETRA_OBJECT_PTR       epetra_blockmap_create2
      EPETRA_OBJECT_PTR       epetra_localblockmap_create1
      EPETRA_OBJECT_PTR       epetra_localblockmap_create2

      EPETRA_OBJECT_PTR       epetra_crsmatrix_create
      EPETRA_OBJECT_PTR       epetra_vbrmatrix_create

      EPETRA_OBJECT_PTR       epetra_multivector_create


      EPETRA_ERROR_FLAG  epetra_vbrmatrix_allocate
      EPETRA_ERROR_FLAG  epetra_vbrmatrix_putblockrow
      EPETRA_ERROR_FLAG  epetra_vbrmatrix_fillcomplete
      EPETRA_ERROR_FLAG  epetra_vbrmatrix_matvec

      EPETRA_ERROR_FLAG  epetra_crsmatrix_allocate
      EPETRA_ERROR_FLAG  epetra_crsmatrix_putrow
      EPETRA_ERROR_FLAG  epetra_crsmatrix_sumintodiagonal
      EPETRA_ERROR_FLAG  epetra_crsmatrix_fillcomplete
      EPETRA_ERROR_FLAG  epetra_crsmatrix_matvec

