#ifndef _PETRA_C_WRAPPERS_H_
#define _PETRA_C_WRAPPERS_H_

#ifdef PETRA_FORTRAN

typedef double * PETRA_DOUBLE;
typedef int    * PETRA_INT;
#define DEREF_ *

#ifdef PETRA_ADDRESS64BIT

typedef long int PETRA_MATRIX_POINTER;
typedef long int PETRA_VECTOR_POINTER;
typedef long int PETRA_MULTIVECTOR_POINTER;
typedef long int PETRA_COMM_POINTER;
typedef long int PETRA_MAP_POINTER;
typedef long int PETRA_LOCALMAP_POINTER;
typedef long int PETRA_BLOCKMAP_POINTER;
typedef long int PETRA_LOCALBLOCKMAP_POINTER;

typedef long int & PETRA_MATRIX;
typedef long int & PETRA_VECTOR;
typedef long int & PETRA_MULTIVECTOR;
typedef long int & PETRA_COMM;
typedef long int & PETRA_MAP;
typedef long int & PETRA_LOCALMAP;
typedef long int & PETRA_BLOCKMAP;
typedef long int & PETRA_LOCALBLOCKMAP;


#else

typedef int PETRA_MATRIX_POINTER;
typedef int PETRA_VECTOR_POINTER;
typedef int PETRA_MULTIVECTOR_POINTER;
typedef int PETRA_COMM_POINTER;
typedef int PETRA_MAP_POINTER;
typedef int PETRA_LOCALMAP_POINTER;
typedef int PETRA_BLOCKMAP_POINTER;
typedef int PETRA_LOCALBLOCKMAP_POINTER;

typedef int & PETRA_MATRIX;
typedef int & PETRA_VECTOR;
typedef int & PETRA_MULTIVECTOR;
typedef int & PETRA_COMM;
typedef int & PETRA_MAP;
typedef int & PETRA_LOCALMAP;
typedef int & PETRA_BLOCKMAP;
typedef int & PETRA_LOCALBLOCKMAP;

#endif
#else

/* These typedefs act as new types for the petra C interface */

typedef double PETRA_DOUBLE;
typedef int    PETRA_INT;
#define DEREF_ 

typedef void * PETRA_MATRIX_POINTER;
typedef void * PETRA_VECTOR_POINTER;
typedef void * PETRA_MULTIVECTOR_POINTER;
typedef void * PETRA_COMM_POINTER;
typedef void * PETRA_MAP_POINTER;
typedef void * PETRA_LOCALMAP_POINTER;
typedef void * PETRA_BLOCKMAP_POINTER;
typedef void * PETRA_LOCALBLOCKMAP_POINTER;

typedef void * PETRA_MATRIX;
typedef void * PETRA_VECTOR;
typedef void * PETRA_MULTIVECTOR;
typedef void * PETRA_COMM;
typedef void * PETRA_MAP;
typedef void * PETRA_LOCALMAP;
typedef void * PETRA_BLOCKMAP;
typedef void * PETRA_LOCALBLOCKMAP;


#endif
 
#if defined(PETRA_FORTRAN) && defined(TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE)
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

PETRA_MATRIX_POINTER MANGLE(petra_rdp_dvbr_matrix_create)(PETRA_MAP rowmap);

int MANGLE(petra_rdp_dvbr_matrix_allocate)(PETRA_MATRIX A, int* numNzBlks, int* blkColInds);

int MANGLE(petra_rdp_dvbr_matrix_putblockrow)(PETRA_MATRIX A, PETRA_INT
blk_row, PETRA_INT num_nz_blocks,
				  double* vals, int* blk_col_inds);

int MANGLE(petra_rdp_dvbr_matrix_fillcomplete)(PETRA_MATRIX A);

int  MANGLE(petra_rdp_dvbr_matrix_matvec)(PETRA_MATRIX A, PETRA_VECTOR x, PETRA_VECTOR y);

int MANGLE(petra_rdp_dvbr_matrix_matmultivec)(PETRA_MATRIX A,
                                          PETRA_MULTIVECTOR x,
                                          PETRA_MULTIVECTOR y);

void MANGLE(petra_rdp_dvbr_matrix_destroy)(PETRA_MATRIX A);

/*****************************************************/
/**                  petra_rdp_dcrs_matrix         **/
/***************************************************/

PETRA_MATRIX_POINTER MANGLE(petra_rdp_dcrs_matrix_create)(PETRA_MAP rowmap);

int  MANGLE(petra_rdp_dcrs_matrix_allocate)(PETRA_MATRIX A, int* rowLengths);

int MANGLE(petra_rdp_dcrs_matrix_putrow)(PETRA_MATRIX A, PETRA_INT row, 
                                         PETRA_INT num_nz,
				                     double* vals, int* col_inds);

int MANGLE(petra_rdp_dcrs_matrix_sumintodiagonal)(PETRA_MATRIX A, 
                                                  double* diagonal);

int MANGLE(petra_rdp_dcrs_matrix_fillcomplete)(PETRA_MATRIX A);

int MANGLE(petra_rdp_dcrs_matrix_matvec)(PETRA_MATRIX A, PETRA_VECTOR x, 
                                         PETRA_VECTOR y);

int MANGLE(petra_rdp_dcrs_matrix_matmultivec)(PETRA_MATRIX A,
                                          PETRA_MULTIVECTOR x,
                                          PETRA_MULTIVECTOR y);

void MANGLE(petra_rdp_dcrs_matrix_destroy)(PETRA_MATRIX A);

/*****************************************************/
/**              petra_rdp_vector                  **/
/***************************************************/

PETRA_VECTOR_POINTER MANGLE(petra_rdp_vector_create)(PETRA_MAP map);

int MANGLE(petra_rdp_vector_putvector)(PETRA_VECTOR x, double *vector);

int MANGLE(petra_rdp_vector_putscalar)(PETRA_VECTOR x, PETRA_DOUBLE scalar);

int MANGLE(petra_rdp_vector_lincomb)(PETRA_VECTOR x, PETRA_VECTOR b,
           PETRA_DOUBLE scalar, PETRA_VECTOR c);

int MANGLE(petra_rdp_vector_norm1)(PETRA_VECTOR x, double *result);

int MANGLE(petra_rdp_vector_norm2)(PETRA_VECTOR x, double *result);

int MANGLE(petra_rdp_vector_random)(PETRA_VECTOR x);

void MANGLE(petra_rdp_vector_destroy)(PETRA_VECTOR x);

/*****************************************************/
/**               petra_rdp_multivector            **/
/***************************************************/

  /* create empty shell WITHOUT float storage, fill later with put functions */
PETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create)();

  /* create empty shell WITH float storage, fill later with put functions */
PETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create1)(PETRA_MAP map, PETRA_INT numvectors);

  /* Build multivector from a Fortran-style 2D array
     NOTE: User storage is not copied, user must keep A intact!! */
PETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create2)(PETRA_MAP map,
                                double *A, PETRA_INT lda, PETRA_INT numvectors);

  /* Build multivector from a double **
     NOTE: User storage is not copied, user must keep A intact!! */
PETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create3)(PETRA_MAP map,
                                double **in_multivector, PETRA_INT numvectors);

  /* Copy constructor */
PETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create4)(PETRA_MULTIVECTOR
                                               in_multivector);

  /* creates a new multivector from numvector number of vectors of an existing
   * multivector where the vectors to be copied are listed in
   * vecIndices.  
   */
PETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create5)(PETRA_MULTIVECTOR 
                        in_multivector, PETRA_INT numvectors, int *vecIndices);

PETRA_MULTIVECTOR_POINTER MANGLE(petra_rdp_multivector_create6)(PETRA_MULTIVECTOR
                        in_multiVector, PETRA_INT startindex, PETRA_INT numvectors);

int MANGLE(petra_rdp_multivector_putmultivector)(PETRA_MULTIVECTOR multivector, 
                                       double **in_multivector);

  /* Allocates space for a multivector created by the default
   * constructor */
int MANGLE(petra_rdp_multivector_allocate)(PETRA_MULTIVECTOR multivector, 
                                PETRA_MAP map, PETRA_INT numvectors);

int MANGLE(petra_rdp_multivector_putscalar)(PETRA_MULTIVECTOR multivector, PETRA_DOUBLE scalar);

int MANGLE(petra_rdp_multivector_scale)
          (PETRA_MULTIVECTOR multiVector, PETRA_DOUBLE scalar);

int MANGLE(petra_rdp_multivector_scalecopy)
          (PETRA_MULTIVECTOR multiVector, PETRA_MULTIVECTOR multiVector_in,
           PETRA_DOUBLE scalar);

int MANGLE(petra_rdp_multivector_dotprod)
          (PETRA_MULTIVECTOR multiVector, PETRA_MULTIVECTOR multiVector_in,
           double *scalar);

int MANGLE(petra_rdp_multivector_addvec)
          (PETRA_MULTIVECTOR multiVector, PETRA_DOUBLE scalar, 
           PETRA_MULTIVECTOR multiVector_in);

int MANGLE(petra_rdp_multivector_norm1)(PETRA_MULTIVECTOR multivector, double *result);

int MANGLE(petra_rdp_multivector_norm2)(PETRA_MULTIVECTOR multivector, double *result);

int MANGLE(petra_rdp_multivector_lincomb)(PETRA_MULTIVECTOR multivector,
                               PETRA_MULTIVECTOR b, 
                               PETRA_DOUBLE scalar, PETRA_MULTIVECTOR c);

int MANGLE(petra_rdp_multivector_random)(PETRA_MULTIVECTOR multivector);

  /* Note: The return value for this function is the number of vectors
     in the multivector */
int MANGLE(petra_rdp_multivector_numvectors)(PETRA_MULTIVECTOR multivector);

int MANGLE(petra_rdp_multivector_reduce)(PETRA_MULTIVECTOR multivector);

int MANGLE(petra_rdp_multivector_gemm)(PETRA_MULTIVECTOR multivector,
			    PETRA_INT transa, PETRA_INT transb, PETRA_DOUBLE alpha,
			    PETRA_MULTIVECTOR A, PETRA_MULTIVECTOR B,
			    PETRA_DOUBLE beta );

void MANGLE(petra_rdp_multivector_destroy)(PETRA_MULTIVECTOR multivector);

/*****************************************************/
/**                  Petra_Map                     **/
/***************************************************/

PETRA_MAP_POINTER MANGLE(petra_map_create1)(PETRA_INT numGlobalEquations, 
                               PETRA_COMM comm);

PETRA_MAP_POINTER MANGLE(petra_map_create2)(PETRA_INT numGlobalEquations, 
                               PETRA_INT numlocalEquations,
                               int *updateList, PETRA_INT indexBase,
                               PETRA_COMM comm);
int MANGLE(petra_map_numlocalequations)(PETRA_MAP map);

#ifndef PETRA_FORTRAN  /* Fortran cannot receive a pointer to int */
int * MANGLE(petra_map_getupdatelist)(PETRA_MAP map);
#endif

PETRA_COMM_POINTER MANGLE(petra_map_getcommunicator)(PETRA_MAP map);

void MANGLE(petra_map_destroy)(PETRA_MAP map);

/*****************************************************/
/**                  petra_blockmap                **/
/***************************************************/

PETRA_BLOCKMAP_POINTER MANGLE(petra_blockmap_create1)(
                PETRA_INT numGlobalEquations, PETRA_INT numlocalEquations, int *updateList,
                PETRA_INT numGlobalblocks, PETRA_INT numlocalblocks, int *blockUpdateList,
                int* blockSizes, PETRA_INT indexBase, PETRA_COMM comm);

PETRA_BLOCKMAP_POINTER MANGLE(petra_blockmap_create2)(
                PETRA_INT numGlobalblocks, PETRA_INT numlocalblocks, int *blockUpdateList,
                int* blockSizes, PETRA_INT indexBase, PETRA_COMM comm);

void MANGLE(petra_blockmap_destroy)(PETRA_BLOCKMAP blockmap);

/*****************************************************/
/**                  petra_localmap                **/
/***************************************************/

PETRA_LOCALMAP_POINTER MANGLE(petra_localmap_create)(PETRA_INT numlocalEquations,
                               PETRA_INT indexBase, PETRA_COMM comm);

void MANGLE(petra_localmap_destroy)(PETRA_LOCALMAP localmap);

/*****************************************************/
/**                  petra_localblockmap           **/
/***************************************************/

PETRA_LOCALBLOCKMAP_POINTER MANGLE(petra_localblockmap_create1)(
                PETRA_INT numlocalEquations,
                PETRA_INT numlocalblocks,
                int* blockSizes,
                PETRA_INT indexBase, PETRA_COMM comm);

PETRA_LOCALBLOCKMAP_POINTER MANGLE(petra_localblockmap_create2)(
                PETRA_INT numlocalblocks,
                int* blockSizes,
                PETRA_INT indexBase, PETRA_COMM comm);

void MANGLE(petra_localblockmap_destroy)(PETRA_LOCALBLOCKMAP localblockmap);

/*****************************************************/
/**                  petra_comm                    **/
/***************************************************/

#ifdef PETRA_MPI
PETRA_COMM_POINTER MANGLE(petra_comm_create)(MPI_Comm MPI_communicator);
#else
PETRA_COMM_POINTER MANGLE(petra_comm_create_serial)();
#endif

int MANGLE(petra_comm_getmypid)(PETRA_COMM communicator);

int MANGLE(petra_comm_getnumproc)(PETRA_COMM communicator);

void MANGLE(petra_comm_barrier)(PETRA_COMM communicator);

void MANGLE(petra_comm_destroy)(PETRA_COMM communicator);

#ifdef __cplusplus
}
#endif

#endif /* _PETRA_C_WRAPPERS_H_ */
