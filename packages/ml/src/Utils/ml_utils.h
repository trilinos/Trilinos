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

#include <sys/time.h>
/*#include "ml_struct.h"*/
#include "ml_defs.h"
#include "ml_comm.h"

#define ML_dabs(x) (((x) > 0.) ? x : (-(x)))
#define ML_abs(x) (((x) > 0) ? x : (-(x)))
#define ML_min(a,b) (((a) <= (b)) ? (a) : (b))
#define ML_max(x,y) (((x) > (y)) ? (x) : (y))



/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* The following is included in the siesta SUN4 environment to solve    */
/* the random number generation problem                                 */
/* -------------------------------------------------------------------- */
#if defined(SUN4) || defined(SUN5)

#ifdef __cplusplus
   extern "C" {
#endif

     double drand48();
     void srand48(long seedval);

#ifdef __cplusplus
   }
#endif

#endif




#define million  0.1e7

#ifdef __cplusplus
   extern "C" {
#endif

   int    pr_error(char *fmt,  ... );
   double GetClock(void);
   void   StartTimer(void);
   void   StopTimer(void);
   double GetElapsedTime(void);

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
   void   ML_random_vec(double u[], int N, ML_Comm *comm);
   double ML_srandom1(int *seed);

   void   ML_serial_start(ML_Comm *comm);
   void   ML_serial_end(ML_Comm *comm);
   int    ML_Coord2RBM(int Nnodes, double x[], double y[], double z[],
                       double rbm[], int Ndof);
   void ML_az_dsort2(double dlist[], int N, int list2[]);

   /* these are functions used by Ray in his RAP thing */

   void   ML_az_sort(int list[], int N, int list2[], double list3[]);
   void   ML_az_dsort(double list[], int N);
   void   ML_gsum_vec_int(int vals[], int vals2[], int, ML_Comm *comm);
   void   ML_rm_duplicates(int array[], int *N);
   void   ML_splitup_big_msg(int, char *, char *, unsigned int, int *,
                             int *, int *, int *, int , int *, ML_Comm *);
   double ML_gdot(int N, double r[], double z[], ML_Comm *comm);
   double ML_gsum_double(double val, ML_Comm *comm);
   double ML_gmax_double(double val, ML_Comm *comm);
   int    ML_gmax_int(int val, ML_Comm *comm);
   int    ML_find_index(int key, int list[], int length);

#ifdef __cplusplus
}
#endif

#endif

