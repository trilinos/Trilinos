/*
//@HEADER
// ************************************************************************
//
//               Epetra: Linear Algebra Services Package
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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


#if defined(Epetra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Epetra package is deprecated"
#endif
#endif
