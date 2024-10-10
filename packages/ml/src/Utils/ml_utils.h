/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Miscellaneous functions for efficient searching and sorting          */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : September, 1998                                      */
/* ******************************************************************** */

#ifndef __MLUTILH__
#define __MLUTILH__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#ifndef __cplusplus
#if defined(ICL) || defined(_MSC_VER)
#include <time.h>
#else
#include <sys/time.h>
#endif
#endif

#ifndef ICL
#ifndef _WIN32
#include <unistd.h>
#endif
#endif

/*#include "ml_struct.h"*/
#include "ml_config.h"
#include "ml_common.h"
#include "ml_defs.h"
#include "ml_comm.h"
#include "ml_operator.h"

#define ML_dabs(x) (((x) > 0.) ? x : (-(x)))
#define ML_abs(x) (((x) > 0) ? x : (-(x)))
#define ML_min(a,b) (((a) <= (b)) ? (a) : (b))
#define ML_max(x,y) (((x) > (y)) ? (x) : (y))
#define ML_avoid_unused_param(x) ML_use_param(x,0)

/* A fast integer hash function written by Bob Jenkins. */
#define ml_hash_function(a) \
{ \
  a = (a+0x7ed55d16) + (a<<12); \
  a = (a^0xc761c23c) ^ (a>>19); \
  a = (a+0x165667b1) + (a<<5); \
  a = (a+0xd3a2646c) ^ (a<<9); \
  a = (a+0xfd7046c5) + (a<<3); \
  a = (a^0xb55a4f09) ^ (a>>16); \
}

#define ML_UseInlinedHashFunction

#ifdef ML_UseInlinedHashFunction
/* work around for compiling on qt (ax_create_stdint_h.m4 didn't work) */
#ifdef _MSC_VER
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef __int64 int64_t;
typedef unsigned __int64 uint64_t;
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#else
typedef unsigned int uint32_t;
#endif

extern uint32_t ml_unew_val;
/* Important: If you want to use ML_fast_hash, the table size must be 2^k for a
   positive integer k. */
#define ML_fast_hash(new_val, hash_table, hlm1, used, hash_index) \
{ \
  ml_unew_val = new_val; \
  ml_hash_function(ml_unew_val); \
  *hash_index = ((int) ml_unew_val) & hlm1; \
  while ( hash_table[*hash_index] != new_val) { \
    if (hash_table[*hash_index] == -1) { (*used)++; break;} \
    (*hash_index)++; \
    *hash_index = (*hash_index) & hlm1; \
  } \
}

#define ML_hash_it(new_val,hash_table,hash_length,used, hash_index) \
{ \
  *hash_index = new_val<<1; \
  if (*hash_index < 0) *hash_index = new_val; \
  *hash_index = (*hash_index) % hash_length; \
  while ( hash_table[*hash_index] != new_val) { \
    if (hash_table[*hash_index] == -1) { (*used)++; break;} \
    (*hash_index)++; \
    *hash_index = (*hash_index) % hash_length; \
  } \
}
#endif /*ifdef ML_UseInlinedHashFunction */

/* JJH FIXME
#ifdef __GNUC__

#define ML_Enter() \
 if(ML_DebugLocation()) { if( NEVADA::comm.rank() ==
NEVADA::comm.IO_processor()) printf("+++ Entering: %s\n",__PRETTY_FUNCTION__);}

 #define ML_Leave() \
 if(ML_DebugLocation()) { if( NEVADA::comm.rank() ==
NEVADA::comm.IO_processor()) printf("--- Leaving: %s\n",__PRETTY_FUNCTION__);}

 #else
 #define MLEnter()
 #define MLLeave()
 #endif
 */

 /*
 #define MLEnter() \
 if(Debug_Location()){ if( NEVADA::comm.rank() == NEVADA::comm.IO_processor())
std::cout << "+++
 Entering: " << __FILE__ << ":" << __LINE__  << std::endl;}

 #define MLLeave() \
 if(Debug_Location()){ if( NEVADA::comm.rank() == NEVADA::comm.IO_processor())
std::cout << "--- Leaving:
 " << __FILE__ << ":" << __LINE__ << std::endl;}
*/


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* The following is included in the siesta SUN4 environment to solve    */
/* the random number generation problem                                 */
/* -------------------------------------------------------------------- */
#if defined(SUN4) || defined(SUN5)

#ifndef ML_CPP
#ifdef __cplusplus
   extern "C" {
#endif
#endif

     double drand48();
     void srand48(long seedval);

#ifndef ML_CPP
#ifdef __cplusplus
   }
#endif
#endif

#endif

#define million  0.1e7

#ifndef ML_CPP
#ifdef __cplusplus
   extern "C" {
#endif
#endif

   void   pr_error(const char *fmt,  ... );
   double GetClock(void);
   void StartTimer(double* t0);
   void StopTimer(double* t0, double* delta);
   void ReportTimer(double t0, const char *msgString, ML_Comm *comm);

   int    ML_crude_search( int, int, int * );
   int    ML_sorted_search( int, int, int * );
   int    ML_fastsorted_search( int, int, int * , int);
   int    ML_sorted_search2( int, int, int *, int, int ** );
   int    ML_search_insert_sort( int, int *, int *, int * );
   int    ML_split_dsort(double list[], int N, int *, int);
   int    ML_selection_dsort(double list[], int N, int *, int);
   int    ML_sort( int, int * );
   void   ML_dsort2(double *, int N, int *);

   int    ML_Check_Context( void * );
   int    ML_randomize( int , int * );
   int    ML_random_init(void);
   int    ML_get_random_seed();
   void   ML_set_random_seed(int seed);
   void   ML_random_vec(double u[], int N, ML_Comm *comm);
   double ML_srandom1(int *seed);

   void   ML_serial_start(ML_Comm *comm);
   void   ML_serial_end(ML_Comm *comm);
   int    ML_Coord2RBM(int Nnodes, double x[], double y[], double z[], double rbm[], int Ndof, int NscalarDof);
   void ML_az_dsort2(double dlist[], int N, int list2[]);

   /* these are functions used by Ray in his RAP thing */

   void   ML_az_sort(int list[], int N, int list2[], double list3[]);
   void   ML_az_dsort(double list[], int N);
   void   ML_gsum_scalar_int(int vals[], int vals2[], ML_Comm *comm);
   void   ML_gsum_vec_int(int *vals[], int *vals2[], int, ML_Comm *comm);
   void   ML_rm_duplicates(int array[], int *N);
   void   ML_splitup_big_msg(int, char *, char *, unsigned int, int *,
                             int *, int *, int *, int , int *, ML_Comm *);
   double ML_gdot(int N, double r[], double z[], ML_Comm *comm);
   double ML_gsum_double(double val, ML_Comm *comm);
   void   ML_gsum_vec_double(double *vals[], double *vals2[], int, ML_Comm *comm);
   double ML_gmax_double(double val, ML_Comm *comm);
   int    ML_gmax_int(int val, ML_Comm *comm);
   int    ML_find_index(int key, int list[], int length);
   void   ML_use_param(void *data, int junk);
   void   ML_BreakForDebugger(ML_Comm *comm);
   void ML_Pause(ML_Comm *comm);
   void ML_print_line (const char *charstr, int ntimes);

   /*MS*/
   int ML_gsum_int(int val, ML_Comm *comm);
   int ML_gmin_int(int val, ML_Comm *comm);
   double ML_gmin_double(double val, ML_Comm *comm);
   /*ms*/
   extern int ML_Operator_Print_UsingGlobalOrdering( ML_Operator *matrix,
                                           const char label[],
                                           int *, int *);
   extern int ML_build_global_numbering( ML_Operator *Amat,
                                         int **pglobal_numbering,
                                         const char *rowsOrCols );


   int ML_Operator_Lump(ML_Operator *A, ML_Operator **B);
   double ML_Global_Standard_Deviation(double sample, int n,
                                       int activeflag, ML_Comm *comm);

   int ML_SetupCoordinates(ML *ml_ptr, int level, int NumPDEEqns,
                        double *in_x_coord, double *in_y_coord,
                        double *in_z_coord);
   int ML_hash_init(int hash_list[], int hash_length, int *hash_used);
#ifndef ML_UseInlinedHashFunction
   void ML_hash_it(int value, int table[], int tableLength,int *spaceUsed,
                   int *hashKey);
   void ML_fast_hash(int value, int table[], int tableLengthMinusOne,
                     int *spaceUsed, int *hashKey);
#endif
   void ML_print_align(int int2match, char *space, int pad);

   int ML_estimate_avg_nz_per_row(ML_Operator * matrix, double * avg_nz);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
