
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

#ifndef EPETRA_C_WRAPPERS_H
#define EPETRA_C_WRAPPERS_H

#ifdef EPETRA_FORTRAN

typedef double * EPETRA_DOUBLE;
typedef int    * EPETRA_INT;
#define EPETRA_DEREF(a) *a

#ifdef EPETRA_ADDRESS64BIT

typedef long int EPETRA_OBJECT_PTR;
typedef long int & EPETRA_OBJECT_REF;

#else

typedef int EPETRA_OBJECT_PTR;
typedef int & EPETRA_OBJECT_REF;

#endif
#else

/* These typedefs act as new types for the Epetra C interface */

typedef double EPETRA_DOUBLE;
typedef int    EPETRA_INT;
#define EPETRA_DEREF(a) a

typedef void * EPETRA_OBJECT_PTR;
typedef void * EPETRA_OBJECT_REF;


#endif
 
#ifdef EPETRA_FORTRAN
#if defined(TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE)
#define MANGLE(x) x
#else
#define MANGLE(x) x ## __
#endif
#else
#define MANGLE(x) x
#endif



#ifdef __cplusplus
extern "C" {
#endif

  /*****************************************************/
  /**                  Epetra_Comm                    **/
  /***************************************************/

#ifdef EPETRA_MPI
  EPETRA_OBJECT_PTR MANGLE(epetra_mpicomm_create1)();
  EPETRA_OBJECT_PTR MANGLE(epetra_mpicomm_create2)(MPI_Comm * comm);
#endif
  EPETRA_OBJECT_PTR MANGLE(epetra_serialcomm_create)();

  int MANGLE(epetra_comm_mypid)(EPETRA_OBJECT_REF communicator);

  int MANGLE(epetra_comm_numproc)(EPETRA_OBJECT_REF communicator);

  void MANGLE(epetra_comm_barrier)(EPETRA_OBJECT_REF communicator);

  void MANGLE(epetra_comm_destroy)(EPETRA_OBJECT_REF communicator);

  /*****************************************************/
  /**                  Epetra_Map                     **/
  /***************************************************/

  EPETRA_OBJECT_PTR MANGLE(epetra_map_create1)(EPETRA_INT numGlobalEquations, 
					       EPETRA_INT indexBase,
					       EPETRA_OBJECT_REF comm);

  EPETRA_OBJECT_PTR MANGLE(epetra_map_create2)(EPETRA_INT numGlobalEquations, 
					       EPETRA_INT numMyElements,
					       EPETRA_INT indexBase,
					       EPETRA_OBJECT_REF comm);

  EPETRA_OBJECT_PTR MANGLE(epetra_map_create3)(EPETRA_INT numGlobalEquations, 
					       EPETRA_INT numlocalEquations,
					       int *updateList, EPETRA_INT indexBase,
					       EPETRA_OBJECT_REF comm);
  int MANGLE(epetra_map_nummyelements)(EPETRA_OBJECT_REF map);
  int MANGLE(epetra_map_numglobalelements)(EPETRA_OBJECT_REF map);

#ifndef EPETRA_FORTRAN  /* Fortran cannot receive a pointer to int */
  int * MANGLE(epetra_map_myglobalelements)(EPETRA_OBJECT_REF map);
#endif

  EPETRA_OBJECT_PTR MANGLE(epetra_map_comm)(EPETRA_OBJECT_REF map);

  void MANGLE(epetra_map_destroy)(EPETRA_OBJECT_REF map);

  /*****************************************************/
  /**              Epetra_Vector                  **/
  /***************************************************/

  EPETRA_OBJECT_PTR MANGLE(epetra_vector_create1)(EPETRA_OBJECT_REF map);
  EPETRA_OBJECT_PTR MANGLE(epetra_vector_create2)(EPETRA_INT Copy, EPETRA_OBJECT_REF map,
						  double * V);

  int MANGLE(epetra_vector_putscalar)(EPETRA_OBJECT_REF x, EPETRA_DOUBLE scalar);

  int MANGLE(epetra_vector_update)(EPETRA_OBJECT_REF x, EPETRA_DOUBLE scalara,
				   EPETRA_OBJECT_REF a,
				   EPETRA_DOUBLE scalarb, EPETRA_OBJECT_REF b,
				   EPETRA_DOUBLE scalarx);

  int MANGLE(epetra_vector_norm1)(EPETRA_OBJECT_REF x, double *result);

  int MANGLE(epetra_vector_norm2)(EPETRA_OBJECT_REF x, double *result);

  int MANGLE(epetra_vector_random)(EPETRA_OBJECT_REF x);

  void MANGLE(epetra_vector_print)(EPETRA_OBJECT_REF x);

  void MANGLE(epetra_vector_destroy)(EPETRA_OBJECT_REF x);

#ifdef SKIP4NOW
  /*****************************************************/
  /**              petra_dvbr_matrix             **/
  /***************************************************/

  EPETRA_OBJECT_PTR MANGLE(epetra_dvbr_matrix_create)(EPETRA_MAP rowmap);

  int MANGLE(epetra_dvbr_matrix_allocate)(EPETRA_MATRIX A, int* numNzBlks, int* blkColInds);

  int MANGLE(epetra_dvbr_matrix_putblockrow)(EPETRA_MATRIX A, EPETRA_INT
					     blk_row, EPETRA_INT num_nz_blocks,
					     double* vals, int* blk_col_inds);

  int MANGLE(epetra_dvbr_matrix_fillcomplete)(EPETRA_MATRIX A);

  int  MANGLE(epetra_dvbr_matrix_matvec)(EPETRA_MATRIX A, EPETRA_VECTOR x, EPETRA_VECTOR y);

  int MANGLE(epetra_dvbr_matrix_matmultivec)(EPETRA_MATRIX A,
					     EPETRA_MULTIVECTOR x,
					     EPETRA_MULTIVECTOR y);

  void MANGLE(epetra_dvbr_matrix_destroy)(EPETRA_MATRIX A);

  /*****************************************************/
  /**                  petra_dcrs_matrix         **/
  /***************************************************/

  EPETRA_OBJECT_PTR MANGLE(epetra_dcrs_matrix_create)(EPETRA_MAP rowmap);

  int  MANGLE(epetra_dcrs_matrix_allocate)(EPETRA_MATRIX A, int* rowLengths);

  int MANGLE(epetra_dcrs_matrix_putrow)(EPETRA_MATRIX A, EPETRA_INT row, 
					EPETRA_INT num_nz,
					double* vals, int* col_inds);

  int MANGLE(epetra_dcrs_matrix_sumintodiagonal)(EPETRA_MATRIX A, 
						 double* diagonal);

  int MANGLE(epetra_dcrs_matrix_fillcomplete)(EPETRA_MATRIX A);

  int MANGLE(epetra_dcrs_matrix_matvec)(EPETRA_MATRIX A, EPETRA_VECTOR x, 
					EPETRA_VECTOR y);

  int MANGLE(epetra_dcrs_matrix_matmultivec)(EPETRA_MATRIX A,
					     EPETRA_MULTIVECTOR x,
					     EPETRA_MULTIVECTOR y);

  void MANGLE(epetra_dcrs_matrix_destroy)(EPETRA_MATRIX A);

  /*****************************************************/
  /**               petra_multivector            **/
  /***************************************************/

  /* create empty shell WITHOUT float storage, fill later with put functions */
  EPETRA_OBJECT_PTR MANGLE(epetra_multivector_create)();

  /* create empty shell WITH float storage, fill later with put functions */
  EPETRA_OBJECT_PTR MANGLE(epetra_multivector_create1)(EPETRA_MAP map, EPETRA_INT numvectors);

  /* Build multivector from a Fortran-style 2D array
     NOTE: User storage is not copied, user must keep A intact!! */
  EPETRA_OBJECT_PTR MANGLE(epetra_multivector_create2)(EPETRA_MAP map,
								double *A, EPETRA_INT lda, EPETRA_INT numvectors);

  /* Build multivector from a double **
     NOTE: User storage is not copied, user must keep A intact!! */
  EPETRA_OBJECT_PTR MANGLE(epetra_multivector_create3)(EPETRA_MAP map,
								double **in_multivector, EPETRA_INT numvectors);

  /* Copy constructor */
  EPETRA_OBJECT_PTR MANGLE(epetra_multivector_create4)(EPETRA_MULTIVECTOR
								in_multivector);

  /* creates a new multivector from numvector number of vectors of an existing
   * multivector where the vectors to be copied are listed in
   * vecIndices.  
   */
  EPETRA_OBJECT_PTR MANGLE(epetra_multivector_create5)(EPETRA_MULTIVECTOR 
								in_multivector, EPETRA_INT numvectors, int *vecIndices);

  EPETRA_OBJECT_PTR MANGLE(epetra_multivector_create6)(EPETRA_MULTIVECTOR
								in_multiVector, EPETRA_INT startindex, EPETRA_INT numvectors);

  int MANGLE(epetra_multivector_putmultivector)(EPETRA_MULTIVECTOR multivector, 
						double **in_multivector);

  /* Allocates space for a multivector created by the default
   * constructor */
  int MANGLE(epetra_multivector_allocate)(EPETRA_MULTIVECTOR multivector, 
					  EPETRA_MAP map, EPETRA_INT numvectors);

  int MANGLE(epetra_multivector_putscalar)(EPETRA_MULTIVECTOR multivector, EPETRA_DOUBLE scalar);

  int MANGLE(epetra_multivector_scale)
       (EPETRA_MULTIVECTOR multiVector, EPETRA_DOUBLE scalar);

  int MANGLE(epetra_multivector_scalecopy)
       (EPETRA_MULTIVECTOR multiVector, EPETRA_MULTIVECTOR multiVector_in,
	EPETRA_DOUBLE scalar);

  int MANGLE(epetra_multivector_dotprod)
       (EPETRA_MULTIVECTOR multiVector, EPETRA_MULTIVECTOR multiVector_in,
	double *scalar);

  int MANGLE(epetra_multivector_addvec)
       (EPETRA_MULTIVECTOR multiVector, EPETRA_DOUBLE scalar, 
	EPETRA_MULTIVECTOR multiVector_in);

  int MANGLE(epetra_multivector_norm1)(EPETRA_MULTIVECTOR multivector, double *result);

  int MANGLE(epetra_multivector_norm2)(EPETRA_MULTIVECTOR multivector, double *result);

  int MANGLE(epetra_multivector_lincomb)(EPETRA_MULTIVECTOR multivector,
					 EPETRA_MULTIVECTOR b, 
					 EPETRA_DOUBLE scalar, EPETRA_MULTIVECTOR c);

  int MANGLE(epetra_multivector_random)(EPETRA_MULTIVECTOR multivector);

  /* Note: The return value for this function is the number of vectors
     in the multivector */
  int MANGLE(epetra_multivector_numvectors)(EPETRA_MULTIVECTOR multivector);

  int MANGLE(epetra_multivector_reduce)(EPETRA_MULTIVECTOR multivector);

  int MANGLE(eepetra_multivector_gemm)(EPETRA_MULTIVECTOR multivector,
				       EPETRA_INT transa, EPETRA_INT transb, EPETRA_DOUBLE alpha,
				       EPETRA_MULTIVECTOR A, EPETRA_MULTIVECTOR B,
				       EPETRA_DOUBLE beta );

  void MANGLE(epetra_multivector_destroy)(EPETRA_MULTIVECTOR multivector);

  /*****************************************************/
  /**                  petra_blockmap                **/
  /***************************************************/

  EPETRA_OBJECT_PTR MANGLE(epetra_blockmap_create1)(
							  EPETRA_INT numGlobalEquations, EPETRA_INT numlocalEquations, int *updateList,
							  EPETRA_INT numGlobalblocks, EPETRA_INT numlocalblocks, int *blockUpdateList,
							  int* blockSizes, EPETRA_INT indexBase, EPETRA_COMM comm);

  EPETRA_OBJECT_PTR MANGLE(epetra_blockmap_create2)(
							  EPETRA_INT numGlobalblocks, EPETRA_INT numlocalblocks, int *blockUpdateList,
							  int* blockSizes, EPETRA_INT indexBase, EPETRA_COMM comm);

  void MANGLE(epetra_blockmap_destroy)(EPETRA_BLOCKMAP blockmap);

  /*****************************************************/
  /**                  petra_localmap                **/
  /***************************************************/

  EPETRA_OBJECT_PTR MANGLE(epetra_localmap_create)(EPETRA_INT numlocalEquations,
							 EPETRA_INT indexBase, EPETRA_COMM comm);

  void MANGLE(epetra_localmap_destroy)(EPETRA_LOCALMAP localmap);

  /*****************************************************/
  /**                  petra_localblockmap           **/
  /***************************************************/

  EPETRA_OBJECT_PTR MANGLE(epetra_localblockmap_create1)(
								    EPETRA_INT numlocalEquations,
								    EPETRA_INT numlocalblocks,
								    int* blockSizes,
								    EPETRA_INT indexBase, EPETRA_COMM comm);

  EPETRA_OBJECT_PTR MANGLE(epetra_localblockmap_create2)(
								    EPETRA_INT numlocalblocks,
								    int* blockSizes,
								    EPETRA_INT indexBase, EPETRA_COMM comm);

  void MANGLE(epetra_localblockmap_destroy)(EPETRA_LOCALBLOCKMAP localblockmap);
#endif /* 0 */

#ifdef __cplusplus
}
#endif

#endif /* EPETRA_C_WRAPPERS_H */
