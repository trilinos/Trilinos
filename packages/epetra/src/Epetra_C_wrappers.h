
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

#ifndef _EPETRA_C_WRAPPERS_H_
#define _EPETRA_C_WRAPPERS_H_

#ifdef EPETRA_FORTRAN

typedef double * EPETRA_DOUBLE;
typedef int    * EPETRA_INT;
#define DEREF_ *

#ifdef EPETRA_ADDRESS64BIT

typedef long int EPETRA_MATRIX_POINTER;
typedef long int EPETRA_VECTOR_POINTER;
typedef long int EPETRA_MULTIVECTOR_POINTER;
typedef long int EPETRA_COMM_POINTER;
typedef long int EPETRA_MAP_POINTER;
typedef long int EPETRA_LOCALMAP_POINTER;
typedef long int EPETRA_BLOCKMAP_POINTER;
typedef long int EPETRA_LOCALBLOCKMAP_POINTER;

typedef long int & EPETRA_MATRIX;
typedef long int & EPETRA_VECTOR;
typedef long int & EPETRA_MULTIVECTOR;
typedef long int & EPETRA_COMM;
typedef long int & EPETRA_MAP;
typedef long int & EPETRA_LOCALMAP;
typedef long int & EPETRA_BLOCKMAP;
typedef long int & EPETRA_LOCALBLOCKMAP;


#else

typedef int EPETRA_MATRIX_POINTER;
typedef int EPETRA_VECTOR_POINTER;
typedef int EPETRA_MULTIVECTOR_POINTER;
typedef int EPETRA_COMM_POINTER;
typedef int EPETRA_MAP_POINTER;
typedef int EPETRA_LOCALMAP_POINTER;
typedef int EPETRA_BLOCKMAP_POINTER;
typedef int EPETRA_LOCALBLOCKMAP_POINTER;

typedef int & EPETRA_MATRIX;
typedef int & EPETRA_VECTOR;
typedef int & EPETRA_MULTIVECTOR;
typedef int & EPETRA_COMM;
typedef int & EPETRA_MAP;
typedef int & EPETRA_LOCALMAP;
typedef int & EPETRA_BLOCKMAP;
typedef int & EPETRA_LOCALBLOCKMAP;

#endif
#else

/* These typedefs act as new types for the petra C interface */

typedef double EPETRA_DOUBLE;
typedef int    EPETRA_INT;
#define DEREF_ 

typedef void * EPETRA_MATRIX_POINTER;
typedef void * EPETRA_VECTOR_POINTER;
typedef void * EPETRA_MULTIVECTOR_POINTER;
typedef void * EPETRA_COMM_POINTER;
typedef void * EPETRA_MAP_POINTER;
typedef void * EPETRA_LOCALMAP_POINTER;
typedef void * EPETRA_BLOCKMAP_POINTER;
typedef void * EPETRA_LOCALBLOCKMAP_POINTER;

typedef void * EPETRA_MATRIX;
typedef void * EPETRA_VECTOR;
typedef void * EPETRA_MULTIVECTOR;
typedef void * EPETRA_COMM;
typedef void * EPETRA_MAP;
typedef void * EPETRA_LOCALMAP;
typedef void * EPETRA_BLOCKMAP;
typedef void * EPETRA_LOCALBLOCKMAP;


#endif
 
#if defined(EPETRA_FORTRAN) && defined(TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE)
#define MANGLE(x) x
#else
#define MANGLE(x) x ## __
#endif


#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************/
/**              petra_rdp_dvbr_matrix             **/
/***************************************************/

EPETRA_MATRIX_POINTER MANGLE(petra_rdp_dvbr_matrix_create)(EPETRA_MAP rowmap);

int MANGLE(petra_rdp_dvbr_matrix_allocate)(EPETRA_MATRIX A, int* numNzBlks, int* blkColInds);

int MANGLE(petra_rdp_dvbr_matrix_putblockrow)(EPETRA_MATRIX A, EPETRA_INT
blk_row, EPETRA_INT num_nz_blocks,
				  double* vals, int* blk_col_inds);

int MANGLE(petra_rdp_dvbr_matrix_fillcomplete)(EPETRA_MATRIX A);

int  MANGLE(petra_rdp_dvbr_matrix_matvec)(EPETRA_MATRIX A, EPETRA_VECTOR x, EPETRA_VECTOR y);

int MANGLE(petra_rdp_dvbr_matrix_matmultivec)(EPETRA_MATRIX A,
                                          EPETRA_MULTIVECTOR x,
                                          EPETRA_MULTIVECTOR y);

void MANGLE(petra_rdp_dvbr_matrix_destroy)(EPETRA_MATRIX A);

/*****************************************************/
/**                  petra_rdp_dcrs_matrix         **/
/***************************************************/

EPETRA_MATRIX_POINTER MANGLE(petra_rdp_dcrs_matrix_create)(EPETRA_MAP rowmap);

int  MANGLE(petra_rdp_dcrs_matrix_allocate)(EPETRA_MATRIX A, int* rowLengths);

int MANGLE(petra_rdp_dcrs_matrix_putrow)(EPETRA_MATRIX A, EPETRA_INT row, 
                                         EPETRA_INT num_nz,
				                     double* vals, int* col_inds);

int MANGLE(petra_rdp_dcrs_matrix_sumintodiagonal)(EPETRA_MATRIX A, 
                                                  double* diagonal);

int MANGLE(petra_rdp_dcrs_matrix_fillcomplete)(EPETRA_MATRIX A);

int MANGLE(petra_rdp_dcrs_matrix_matvec)(EPETRA_MATRIX A, EPETRA_VECTOR x, 
                                         EPETRA_VECTOR y);

int MANGLE(petra_rdp_dcrs_matrix_matmultivec)(EPETRA_MATRIX A,
                                          EPETRA_MULTIVECTOR x,
                                          EPETRA_MULTIVECTOR y);

void MANGLE(petra_rdp_dcrs_matrix_destroy)(EPETRA_MATRIX A);

/*****************************************************/
/**              petra_rdp_vector                  **/
/***************************************************/

EPETRA_VECTOR_POINTER MANGLE(petra_rdp_vector_create)(EPETRA_MAP map);

int MANGLE(petra_rdp_vector_putvector)(EPETRA_VECTOR x, double *vector);

int MANGLE(petra_rdp_vector_putscalar)(EPETRA_VECTOR x, EPETRA_DOUBLE scalar);

int MANGLE(petra_rdp_vector_lincomb)(EPETRA_VECTOR x, EPETRA_VECTOR b,
           EPETRA_DOUBLE scalar, EPETRA_VECTOR c);

int MANGLE(petra_rdp_vector_norm1)(EPETRA_VECTOR x, double *result);

int MANGLE(petra_rdp_vector_norm2)(EPETRA_VECTOR x, double *result);

int MANGLE(petra_rdp_vector_random)(EPETRA_VECTOR x);

void MANGLE(petra_rdp_vector_destroy)(EPETRA_VECTOR x);

/*****************************************************/
/**               petra_rdp_multivector            **/
/***************************************************/

  /* create empty shell WITHOUT float storage, fill later with put functions */
EPETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create)();

  /* create empty shell WITH float storage, fill later with put functions */
EPETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create1)(EPETRA_MAP map, EPETRA_INT numvectors);

  /* Build multivector from a Fortran-style 2D array
     NOTE: User storage is not copied, user must keep A intact!! */
EPETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create2)(EPETRA_MAP map,
                                double *A, EPETRA_INT lda, EPETRA_INT numvectors);

  /* Build multivector from a double **
     NOTE: User storage is not copied, user must keep A intact!! */
EPETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create3)(EPETRA_MAP map,
                                double **in_multivector, EPETRA_INT numvectors);

  /* Copy constructor */
EPETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create4)(EPETRA_MULTIVECTOR
                                               in_multivector);

  /* creates a new multivector from numvector number of vectors of an existing
   * multivector where the vectors to be copied are listed in
   * vecIndices.  
   */
EPETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create5)(EPETRA_MULTIVECTOR 
                        in_multivector, EPETRA_INT numvectors, int *vecIndices);

EPETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create6)(EPETRA_MULTIVECTOR
                        in_multiVector, EPETRA_INT startindex, EPETRA_INT numvectors);

int MANGLE(petra_rdp_multivector_putmultivector)(EPETRA_MULTIVECTOR multivector, 
                                       double **in_multivector);

  /* Allocates space for a multivector created by the default
   * constructor */
int MANGLE(petra_rdp_multivector_allocate)(EPETRA_MULTIVECTOR multivector, 
                                EPETRA_MAP map, EPETRA_INT numvectors);

int MANGLE(petra_rdp_multivector_putscalar)(EPETRA_MULTIVECTOR multivector, EPETRA_DOUBLE scalar);

int MANGLE(petra_rdp_multivector_scale)
          (EPETRA_MULTIVECTOR multiVector, EPETRA_DOUBLE scalar);

int MANGLE(petra_rdp_multivector_scalecopy)
          (EPETRA_MULTIVECTOR multiVector, EPETRA_MULTIVECTOR multiVector_in,
           EPETRA_DOUBLE scalar);

int MANGLE(petra_rdp_multivector_dotprod)
          (EPETRA_MULTIVECTOR multiVector, EPETRA_MULTIVECTOR multiVector_in,
           double *scalar);

int MANGLE(petra_rdp_multivector_addvec)
          (EPETRA_MULTIVECTOR multiVector, EPETRA_DOUBLE scalar, 
           EPETRA_MULTIVECTOR multiVector_in);

int MANGLE(petra_rdp_multivector_norm1)(EPETRA_MULTIVECTOR multivector, double *result);

int MANGLE(petra_rdp_multivector_norm2)(EPETRA_MULTIVECTOR multivector, double *result);

int MANGLE(petra_rdp_multivector_lincomb)(EPETRA_MULTIVECTOR multivector,
                               EPETRA_MULTIVECTOR b, 
                               EPETRA_DOUBLE scalar, EPETRA_MULTIVECTOR c);

int MANGLE(petra_rdp_multivector_random)(EPETRA_MULTIVECTOR multivector);

  /* Note: The return value for this function is the number of vectors
     in the multivector */
int MANGLE(petra_rdp_multivector_numvectors)(EPETRA_MULTIVECTOR multivector);

int MANGLE(petra_rdp_multivector_reduce)(EPETRA_MULTIVECTOR multivector);

int MANGLE(petra_rdp_multivector_gemm)(EPETRA_MULTIVECTOR multivector,
			    EPETRA_INT transa, EPETRA_INT transb, EPETRA_DOUBLE alpha,
			    EPETRA_MULTIVECTOR A, EPETRA_MULTIVECTOR B,
			    EPETRA_DOUBLE beta );

void MANGLE(petra_rdp_multivector_destroy)(EPETRA_MULTIVECTOR multivector);

/*****************************************************/
/**                  Epetra_Map                     **/
/***************************************************/

EPETRA_MAP_POINTER MANGLE(petra_map_create1)(EPETRA_INT numGlobalEquations, 
                               EPETRA_COMM comm);

EPETRA_MAP_POINTER MANGLE(petra_map_create2)(EPETRA_INT numGlobalEquations, 
                               EPETRA_INT numlocalEquations,
                               int *updateList, EPETRA_INT indexBase,
                               EPETRA_COMM comm);
int MANGLE(petra_map_numlocalequations)(EPETRA_MAP map);

#ifndef EPETRA_FORTRAN  /* Fortran cannot receive a pointer to int */
int * MANGLE(petra_map_getupdatelist)(EPETRA_MAP map);
#endif

EPETRA_COMM_POINTER MANGLE(petra_map_getcommunicator)(EPETRA_MAP map);

void MANGLE(petra_map_destroy)(EPETRA_MAP map);

/*****************************************************/
/**                  petra_blockmap                **/
/***************************************************/

EPETRA_BLOCKMAP_POINTER MANGLE(petra_blockmap_create1)(
                EPETRA_INT numGlobalEquations, EPETRA_INT numlocalEquations, int *updateList,
                EPETRA_INT numGlobalblocks, EPETRA_INT numlocalblocks, int *blockUpdateList,
                int* blockSizes, EPETRA_INT indexBase, EPETRA_COMM comm);

EPETRA_BLOCKMAP_POINTER MANGLE(petra_blockmap_create2)(
                EPETRA_INT numGlobalblocks, EPETRA_INT numlocalblocks, int *blockUpdateList,
                int* blockSizes, EPETRA_INT indexBase, EPETRA_COMM comm);

void MANGLE(petra_blockmap_destroy)(EPETRA_BLOCKMAP blockmap);

/*****************************************************/
/**                  petra_localmap                **/
/***************************************************/

EPETRA_LOCALMAP_POINTER MANGLE(petra_localmap_create)(EPETRA_INT numlocalEquations,
                               EPETRA_INT indexBase, EPETRA_COMM comm);

void MANGLE(petra_localmap_destroy)(EPETRA_LOCALMAP localmap);

/*****************************************************/
/**                  petra_localblockmap           **/
/***************************************************/

EPETRA_LOCALBLOCKMAP_POINTER MANGLE(petra_localblockmap_create1)(
                EPETRA_INT numlocalEquations,
                EPETRA_INT numlocalblocks,
                int* blockSizes,
                EPETRA_INT indexBase, EPETRA_COMM comm);

EPETRA_LOCALBLOCKMAP_POINTER MANGLE(petra_localblockmap_create2)(
                EPETRA_INT numlocalblocks,
                int* blockSizes,
                EPETRA_INT indexBase, EPETRA_COMM comm);

void MANGLE(petra_localblockmap_destroy)(EPETRA_LOCALBLOCKMAP localblockmap);

/*****************************************************/
/**                  petra_comm                    **/
/***************************************************/

#ifdef EPETRA_MPI
EPETRA_COMM_POINTER MANGLE(petra_comm_create)(MPI_Comm MPI_communicator);
#else
EPETRA_COMM_POINTER MANGLE(petra_comm_create_serial)();
#endif

int MANGLE(petra_comm_getmypid)(EPETRA_COMM communicator);

int MANGLE(petra_comm_getnumproc)(EPETRA_COMM communicator);

void MANGLE(petra_comm_barrier)(EPETRA_COMM communicator);

void MANGLE(petra_comm_destroy)(EPETRA_COMM communicator);

#ifdef __cplusplus
}
#endif

#endif /* _EPETRA_C_WRAPPERS_H_ */
