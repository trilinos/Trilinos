
C * Copyright (2001) Sandia Corportation. Under the terms of Contract 
C * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
C * work by or on behalf of the U.S. Government.  Export of this program
C * may require a license from the United States Government.


C * NOTICE:  The United States Government is granted for itself and others
C * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
C * license in ths data to reproduce, prepare derivative works, and
C * perform publicly and display publicly.  Beginning five (5) years from
C * July 25, 2001, the United States Government is granted for itself and
C * others acting on its behalf a paid-up, nonexclusive, irrevocable
C * worldwide license in this data to reproduce, prepare derivative works,
C * distribute copies to the public, perform publicly and display
C * publicly, and to permit others to do so.
C * 
C * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
C * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
C * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
C * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
C * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
C * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.


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

