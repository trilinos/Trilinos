
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */


#define EPETRA_FORTRAN

#ifdef EPETRA_ADDRESS64BIT

#define EPETRA_MATRIX long integer
#define EPETRA_VECTOR long integer
#define EPETRA_MULTIVECTOR long integer
#define EPETRA_COMM long integer
#define EPETRA_MAP long integer
#define EPETRA_LOCALMAP long integer
#define EPETRA_BLOCKMAP long integer
#define EPETRA_LOCALBLOCKMAP long integer

#else

#define EPETRA_MATRIX integer
#define EPETRA_VECTOR integer
#define EPETRA_MULTIVECTOR integer
#define EPETRA_COMM integer
#define EPETRA_MAP integer
#define EPETRA_LOCALMAP integer
#define EPETRA_BLOCKMAP integer
#define EPETRA_LOCALBLOCKMAP integer

#endif

#define EPETRA_ERROR_FLAG integer

      EPETRA_COMM             petra_comm_create
      EPETRA_COMM             petra_comm_create_serial

      EPETRA_MAP              petra_map_create
      EPETRA_LOCALMAP         petra_localmap_create
      EPETRA_BLOCKMAP         petra_blockmap_create1
      EPETRA_BLOCKMAP         petra_blockmap_create2
      EPETRA_LOCALBLOCKMAP    petra_localblockmap_create1
      EPETRA_LOCALBLOCKMAP    petra_localblockmap_create2

      EPETRA_MATRIX           petra_rdp_dcrs_matrix_create
      EPETRA_MATRIX           petra_rdp_dvbr_matrix_create

      EPETRA_VECTOR           petra_rdp_vector_create
      EPETRA_MULTIVECTOR     petra_rdp_multivector_create

      EPETRA_ERROR_FLAG  petra_comm_getmypid
      EPETRA_ERROR_FLAG  petra_comm_getnumproc

      EPETRA_ERROR_FLAG  petra_rdp_dvbr_matrix_allocate
      EPETRA_ERROR_FLAG  petra_rdp_dvbr_matrix_putblockrow
      EPETRA_ERROR_FLAG  petra_rdp_dvbr_matrix_fillcomplete
      EPETRA_ERROR_FLAG  petra_rdp_dvbr_matrix_matvec

      EPETRA_ERROR_FLAG  petra_rdp_dcrs_matrix_allocate
      EPETRA_ERROR_FLAG  petra_rdp_dcrs_matrix_putrow
      EPETRA_ERROR_FLAG  petra_rdp_dcrs_matrix_sumintodiagonal
      EPETRA_ERROR_FLAG  petra_rdp_dcrs_matrix_fillcomplete
      EPETRA_ERROR_FLAG  petra_rdp_dcrs_matrix_matvec

      EPETRA_ERROR_FLAG  petra_rdp_vector_putvector
      EPETRA_ERROR_FLAG  petra_rdp_vector_lincomb
      EPETRA_ERROR_FLAG  petra_rdp_vector_norm2
