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

/*
Current revision: $Revision$
Branch:           $Branch$
Last modified:    $Date$
Modified by:      $Author$
*/

#include <math.h>
#include <stdlib.h>
#include "ml_utils.h"
#include "ml_comm.h"
#include "ml_lapack.h"
#include <time.h>
#include "ml_viz_stats.h"
#ifdef _MSC_VER
#pragma comment(lib, "Ws2_32.lib") 
# include <Winsock2.h>
# include <process.h>
void sleep(int sec)
{
  Sleep(sec * 1000);
}
#endif

/* For hashing macro defined in ml_utils.h. */
uint32_t ml_unew_val;

/* ******************************************************************** */
/* Timer subroutine                                                     */
/* ******************************************************************** */

/* Generic timer ... I hope it works on lots of different machines */
#ifdef AZTEC

#ifdef __cplusplus
extern "C" 
{
#endif
extern double machine_dependent_second(void);
#ifdef __cplusplus
}
#endif

#endif
double GetClock(void)
{
#ifdef AZTEC
return( machine_dependent_second());
#else
#ifdef SMOS
  double dclock();
  return (dclock());
#else
  static clock_t last_num_ticks = 0;
  static double  inv_clocks_per_sec = 1./(double)CLOCKS_PER_SEC;
  static double  clock_width =
    (double)(1L<<((int)sizeof(clock_t)*8-2))*4./(double)CLOCKS_PER_SEC;
  static int     clock_rollovers = 0;
  double value;
  clock_t num_ticks = clock();
  if(num_ticks < last_num_ticks) clock_rollovers++;
  value = num_ticks * inv_clocks_per_sec;
  if(clock_rollovers) value += clock_rollovers * clock_width;
  last_num_ticks = num_ticks;
  return(value);
#endif
#endif
}


/* Random Numbers */
static int ml_random_seed = 0;
static int ml_random_start = 1;



/* ******************************************************************** */
/* StartTimer                                                           */
/*   t0   (out) start time                                              */ 
/* ******************************************************************** */
void StartTimer(double* t0)
{
#ifdef ML_TIMING
  *t0 = GetClock();
#else
  return;
#endif
}

/* ******************************************************************** */
/* StopTimer                                                            */
/*   t0   (in/out) On entry, the start time recorded by a previous call */
/*                 to either StartTimer() or StopTimer().               */
/*                 On exit, the stop time.                              */ 
/*   delta (out)   Difference between current time and that passed      */
/*                 in via t0.                                           */
/* ******************************************************************** */

void StopTimer(double* t0, double *delta)
{
#ifdef ML_TIMING
  double now;

  now = GetClock();
  *delta = now - *t0;
  *t0 = now;
#else
  return;
#endif
}

/* ******************************************************************** */
/* ReportTimer                                                          */
/*   Reports average time over all processors & standard deviation.     */
/*   Nothing prints if the ML print level is 0.                         */
/*                                                                      */
/*   t0        (in) The time to report.  This will be averaged over all */
/*                  processors.                                         */
/*   msgString (in) Pointer to message to print with time reported.     */ 
/*   comm      (in) Pointer to an ML communicator.                      */
/* ******************************************************************** */

void ReportTimer(double t0, char *msgString, ML_Comm *comm)
{
#ifdef ML_TIMING
  double t1;
#ifdef ML_MPI
  ml_DblLoc srct,mint,maxt;
#endif
  if (ML_Get_PrintLevel() == 0)
    return;
  t1 = ML_gsum_double(t0, comm);
  t1 = t1/((double) comm->ML_nprocs);
# ifdef ML_MPI
  srct.value = t0; 
  srct.rank = comm->ML_mypid;
  MPI_Reduce(&srct,&maxt,1,MPI_DOUBLE_INT,MPI_MAXLOC,0,comm->USR_comm);
  MPI_Reduce(&srct,&mint,1,MPI_DOUBLE_INT,MPI_MINLOC,0,comm->USR_comm);
  if (comm->ML_mypid==0 && comm->ML_nprocs > 1)
    printf("%s \t avg = %1.3e seconds, min = %1.3e (%d), max = %1.3e (%d)\n",
           msgString, t1, mint.value, mint.rank, maxt.value, maxt.rank);
# endif /*ifdef ML_MPI*/
  if (comm->ML_nprocs == 1)
    printf("%s  = %1.3e seconds\n",msgString,t1);
#else
  return;
#endif /*ifdef ML_TIMING*/
}

/* ******************************************************************** */
/* Given an unsorted list of indices and the key, see whether the key   */
/* is in the list                                                       */
/* -------------------------------------------------------------------- */

int ML_crude_search(int key, int nlist, int *list) 
{
   int  i, found=0;

   for ( i = 0; i < nlist; i++ )
   {
      if ( list[i] == key ) 
      {
         found = 1;
         break;
      }
   }
   if ( found == 1 ) return 0;
   else              return -1;
}

/* ******************************************************************** */
/* Given a sorted list of indices and the key, find the position of the */
/* key in the list.  If not found, return the index of the position     */
/* corresponding to where it would have been stored.                    */
/* -------------------------------------------------------------------- */

int ML_fastsorted_search(int key, int nlist, int *list, int init_guess) 
{
   int  nfirst, nlast, nmid, found, index = 0;

   if (nlist <= 0) return -1;
   if (list[init_guess] == key) return init_guess;
   if (list[init_guess] > key ) {
     nlast = init_guess;
     nfirst = init_guess - 5;
     if (nfirst < 0) nfirst = 0;
     while ( list[nfirst] > key) {
       nfirst -= 5;
       if (nfirst < 0) nfirst = 0;
     }
   }
   else {
     nfirst = init_guess;
     nlast = init_guess + 5;
     if (nlast > nlist-1) nlast = nlist-1;
     while ( list[nlast] < key) {
       nlast += 5;
       if (nlast > nlist-1) nlast = nlist-1;
     }
   }
   if (key > list[nlast])  return -(nlast+1);
   if (key < list[nfirst]) return -(nfirst+1);
   found = 0;
   while ((found == 0) && ((nlast-nfirst)>1)) {
      nmid = (nfirst + nlast) / 2;
      if (key == list[nmid])     {index  = nmid; found = 1;}
      else if (key > list[nmid])  nfirst = nmid;
      else                        nlast  = nmid;
   }
   if (found == 1)               return index;
   else if (key == list[nfirst]) return nfirst;
   else if (key == list[nlast])  return nlast;
   else                          return -(nfirst+1);
}
/* ******************************************************************** */
/* Given a sorted list of indices and the key, find the position of the */
/* key in the list.  If not found, return the index of the position     */
/* corresponding to where it would have been stored.                    */
/* -------------------------------------------------------------------- */

int ML_sorted_search(int key, int nlist, int *list) 
{
   int  nfirst, nlast, nmid, found, index = 0;

   if (nlist <= 0) return -1;
   nfirst = 0;  
   nlast  = nlist-1;
   if (key > list[nlast])  return -(nlast+1);
   if (key < list[nfirst]) return -(nfirst+1);
   found = 0;
   while ((found == 0) && ((nlast-nfirst)>1)) {
      nmid = (nfirst + nlast) / 2;
      if (key == list[nmid])     {index  = nmid; found = 1;}
      else if (key > list[nmid])  nfirst = nmid;
      else                        nlast  = nmid;
   }
   if (found == 1)               return index;
   else if (key == list[nfirst]) return nfirst;
   else if (key == list[nlast])  return nlast;
   else                          return -(nfirst+1);
}

/* ******************************************************************** */
/* This is another search subroutine which saves space by using a bit   */
/* map to store where a given key has been accessed before.  It only    */
/* returns a success if the key is found and it has never been accessed */
/* before.                                                              */
/* -------------------------------------------------------------------- */

int ML_sorted_search2(int key,int nlist,int *list,int cnum,int **map) 
{
   int  nfirst, nlast, nmid, found, index = 0, retdata = 0, col, digit, mask;
   int  nbit_int=sizeof(int)*8;

   nfirst = 0;  
   nlast  = nlist-1;
   found = 0;
   while ((found == 0) && ((nlast-nfirst)>1)) {
      nmid = (nfirst + nlast) / 2;
      if (key == list[nmid])     {index  = nmid; found = 1;}
      else if (key > list[nmid])  nfirst = nmid;
      else                        nlast  = nmid;
   }
   if (found == 1)               retdata = index;
   else if (key == list[nfirst]) retdata = nfirst;
   else if (key == list[nlast])  retdata = nlast;
   col   = cnum / nbit_int;
   digit = cnum % nbit_int;
   mask  = 1 << digit;
   if ((map[retdata][col] & mask) == 0) {
      map[retdata][col] = map[retdata][col] | mask;
      return retdata;
   } else return -1;   
}

/* ******************************************************************** */
/* Given a sorted list of indices and the key, this subroutine finds    */
/* the the position of the key in the list.  If not found, insert the   */
/* key at the appropriate place. (This will be slow.)                   */
/* -------------------------------------------------------------------- */

int ML_search_insert_sort(int key, int *list, int *nlist, int *cnt_list)
{
   int  i, index, n;

   n = (*nlist);
   index = ML_sorted_search(key, *nlist, list);
   if (index < 0) {
      index = - (index + 1);
      if ( n == 0 ) { 
         index = 0; 
         list[0] = key; 
         if ( cnt_list != 0 ) cnt_list[0] = 1;
      } else {
         for ( i = n-1; i > index; i-- ) list[i+1] = list[i];
         if ( cnt_list != 0 ) 
            for ( i = n-1; i > index; i-- ) cnt_list[i+1] = cnt_list[i];
         if (key > list[index]) {
            index++;
            list[index] = key;
            if ( cnt_list != 0 ) cnt_list[index] = 1; 
         } else {
            list[index+1] = list[index]; 
            list[index] = key;
            if ( cnt_list != 0 ) {
               cnt_list[index+1] = cnt_list[index];
               cnt_list[index] = 1;
            }
         }
      }
      (*nlist)++;
   } else {
      if ( cnt_list != 0 ) cnt_list[index]++;
   }
   return index;
}

/* ******************************************************************** */
/* Check a given object to see what it is.                              */
/* -------------------------------------------------------------------- */

int ML_Check_Context( void *obj )
{
   int     id;
   ML_Comm *comm;

   comm = (ML_Comm *) obj;
   
   id = comm->ML_id;
   switch ( id ) {
      case ML_ID_ML : printf("It is a ML object.\n");
                      break;
      case ML_ID_SL : printf("It is a SL object.\n");
                      break;
      default :       printf("Object not recognized. \n");
   }
   return id;
}

/* ******************************************************************** */
/* sort a given list in increasing order                                */
/* -------------------------------------------------------------------- */

int ML_sort(int nlist, int *list) 
{
   int  i, key, *cnt1_array, *cnt2_array, begin, count1, count2;

   if ( nlist <= 1 ) return 0;
   if ( nlist == 2 ) 
   {
      if ( list[0] > list[1] ) 
      {
         key = list[0];
         list[0] = list[1];
         list[1] = key;
      }
      return 0;
   }
   key = list[0];
   count1 = 0;
   count2 = 0;
   cnt1_array = (int*) ML_allocate( nlist * sizeof(int) );
   cnt2_array = (int*) ML_allocate( nlist * sizeof(int) );
   for ( i = 1; i < nlist; i++ ) 
   {
      if ( list[i] <  key ) cnt1_array[count1++] = list[i];
      if ( list[i] >= key ) cnt2_array[count2++] = list[i];
   }
   for ( i = 0; i < count1; i++ ) list[i] = cnt1_array[i];
   list[count1] = key;
   for ( i = 0; i < count2; i++ ) list[count1+1+i] = cnt2_array[i];
   ML_free( cnt1_array );
   ML_free( cnt2_array );
   ML_sort( count1, list );
   begin = count1+1;
   for ( i = count1+1; i < nlist; i++ ) 
   {
      if ( list[i] != key ) break;
      else                  { begin++; count2--; }
   }
   ML_sort( count2, &list[begin] );
   return 0;
}

/* ******************************************************************** */
/* sort a given list in increasing order                                */
/* -------------------------------------------------------------------- */

int ML_split_dsort(double *dlist, int nlist, int *ilist, int limit) 
{
   int    itemp, *iarray1, *iarray2, count1, count2, i;
   double dtemp, *darray1, *darray2;

   if ( nlist <= 1 ) return 0;
   if ( nlist == 2 ) 
   {
      if ( dlist[0] < dlist[1] ) 
      {
         dtemp = dlist[0]; dlist[0] = dlist[1]; dlist[1] = dtemp;
         itemp = ilist[0]; ilist[0] = ilist[1]; ilist[1] = itemp;
      }
      return 0;
   }
   count1 = 0;
   count2 = 0;
   iarray1 = (int *)   ML_allocate( 2 * nlist * sizeof(int) );
   iarray2 = iarray1 + nlist;
   darray1 = (double*) ML_allocate( 2 * nlist * sizeof(double) );
   darray2 = darray1 + nlist;

   if ( darray2 == NULL )
   {
      printf("ERROR : malloc\n");
      exit(1);
   }
   dtemp  = dlist[0];
   itemp  = ilist[0];
   for ( i = 1; i < nlist; i++ ) 
   {
      if (dlist[i] >= dtemp  ) 
      {
         darray1[count1] = dlist[i];
         iarray1[count1++] = ilist[i];
      } 
      else if (dlist[i] <  dtemp) 
      {
         darray2[count2] = dlist[i];
         iarray2[count2++] = ilist[i];
      }
   }
   dlist[count1] = dtemp;
   ilist[count1] = itemp;
   for ( i = 0; i < count1; i++ ) 
   {
      dlist[i] = darray1[i];
      ilist[i] = iarray1[i];
   }
   for ( i = 0; i < count2; i++ ) 
   {
      dlist[count1+1+i] = darray2[i];
      ilist[count1+1+i] = iarray2[i];
   }
   ML_free( darray1 );
   ML_free( iarray1 );
   if ( count1+1 == limit ) return 0;
   else if ( count1+1 < limit )
      ML_split_dsort(&(dlist[count1+1]),count2,&(ilist[count1+1]),limit-count1-1);
   else
      ML_split_dsort( dlist, count1, ilist, limit );
   return 0;
}

/* ******************************************************************** */
/* selection sort                                                       */
/* -------------------------------------------------------------------- */

int ML_selection_dsort(double *vals, int length, int *cols, int limit)
{
   int    i, k, ind1, ind2, ind3, loopcnt, base, level, newcount;
   double *darray, *darray1, *darray2, **treeArray;
   int    expLeng, tmpLeng, *treeLengs, **treeIArray, *iarray, *iarray1, *iarray2;
   int    col1, col2;
   double val1, val2;

   level = (int) (log( 2.0 * length ) / log( 2.0 ));
   printf("level = %d\n", level);

   /* set up data structure */

   expLeng    = (int) pow(2.0, (double) (level+1));
   iarray     = (int    *)  ML_allocate(expLeng   * sizeof(int));
   darray     = (double *)  ML_allocate(expLeng   * sizeof(double));
   treeLengs  = (int *)     ML_allocate((level+1) * sizeof(int));
   treeArray  = (double **) ML_allocate((level+1) * sizeof(double*));
   treeIArray = (int **)    ML_allocate((level+1) * sizeof(int*));
   treeLengs[level]  = length;
   base              = expLeng >> 1;
   treeArray[level]  = &(darray[base]);
   treeIArray[level] = &(iarray[base]);
   for ( k = level-1; k >= 0; k-- )
   {
      base = base >> 1;
      treeArray[k] = &(darray[base]);
      treeIArray[k] = &(iarray[base]);
      treeLengs[k] = ( treeLengs[k+1] + 1 ) / 2;
   }
   darray1 = treeArray[level];
   iarray1 = treeIArray[level];
   for ( i = 0; i < length; i++ )
   {
      darray1[i] = vals[i];
      iarray1[i] = cols[i];
   }
   if ( length < expLeng )
   {
      darray1[length] = 0.0;
      iarray1[length] = 0;
   }

   /* pre-sort */

   for ( k = level; k > 0; k-- )
   {
      darray1 = treeArray[k];
      iarray1 = treeIArray[k];
      darray2 = treeArray[k-1];
      iarray2 = treeIArray[k-1];
      tmpLeng = treeLengs[k];
      for ( i = 0; i < tmpLeng; i+=2 )
      {
         ind1 = i + 1;
         ind2 = i >> 1;
         if ( darray1[i] > darray1[ind1] )
         {
            iarray2[ind2] = iarray1[i];
            darray2[ind2] = darray1[i];
         }
         else
         {
            iarray2[ind2] = iarray1[ind1];
            darray2[ind2] = darray1[ind1];
         }
      }
      if ( tmpLeng % 2 == 1 )
      {
         iarray2[treeLengs[k-1]-1] = iarray1[tmpLeng-1];
         darray2[treeLengs[k-1]-1] = darray1[tmpLeng-1];
      }
   }

   /* post-sort */

   loopcnt = limit;
   newcount = 0;
   while ( loopcnt > 0 )
   {
      if ( loopcnt % 100000 == 0 ) printf("loopcnt = %d\n", loopcnt);
      vals[newcount] = darray[1];
      cols[newcount++] = iarray[1];
      darray1 = treeArray[level];
      darray1[iarray[1]] = 0.0;
      ind3 = iarray[1] >> 1;
      ind1 = ind3 << 1;
      ind2 = ind1 + 1;
      for ( k = level; k > 0; k-- )
      {
         iarray2 = treeIArray[k-1];
         darray2 = treeArray[k-1];
         col1 = treeIArray[k][ind1];
         col2 = treeIArray[k][ind2];
         val1 = treeArray[k][ind1];
         val2 = treeArray[k][ind2];
         if ( val1 > val2 )
         {
            iarray2[ind3] = col1;
            darray2[ind3] = val1;
         }
         else
         {
            iarray2[ind3] = col2;
            darray2[ind3] = val2;
         }
         ind3 = ind3 >> 1;
         ind1 = ind3 << 1;
         ind2 = ind1 + 1;
      }
      loopcnt--;
   }

   ML_free(darray);
   ML_free(treeArray);
   ML_free(treeLengs);
   return 0;
}

/* ******************************************************************** */
/* randomize an integer array                                           */
/* -------------------------------------------------------------------- */

int ML_random_init() 
{
/*
   double stime;
   stime = GetClock();
*/
   return 0;
}

/* ******************************************************************** */
/* randomize an integer array                                           */
/* -------------------------------------------------------------------- */
int ML_randomize(int nlist, int *list) 
{
   int    i, nm1, iran1, iran2, itmp;

   nm1    = nlist - 1;
   for ( i = 0; i < 3*nlist; i++ )
   {
      iran1 = (int) ( nm1 * ML_srandom1(&ml_random_seed) );
      iran2 = (int) ( nm1 * ML_srandom1(&ml_random_seed) );
      if ( iran1 != iran2 )
      {
         itmp        = list[iran2];
         list[iran2] = list[iran1];
         list[iran1] = itmp;
      }
   }
   return 0;
}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

/**************************************************************************

  This routine was taken from Knuth: Sorting and Searching. It puts the input
  data list into a heap and then sorts it.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     double, maximum value in vector 'vec'
  ============

  Parameter list:
  ===============

  list:            On input, values to be sorted. On output, sorted values
                   (i.e., list[i] <= list[i+1]).

  N:               length of vector 'vec'.

  list2:           If on input,
                   a) list2 = NULL: it is unchanged on output,
                   b) list2 is a list associated with 'list':
                   on output, if list[k] on input is now element 'j' on output,
                   list2[j] on output is list2[k].

  list3:           If on input,
                   a) list3 = NULL: it is unchanged on output,
                   b) list3 is a list associated with 'list':
                   on output, list3[j] is assigned the input value of list3[k],
                   if list[j] has been assigned the input value of list[k].

**************************************************************************/

void ML_az_dsort(double list[], int N)
{

  /* local variables */

  int    l, r, j, i, flag;
  double RR, K;

  /*********************** execution begins ******************************/

  if (N <= 1) return;

  l   = N / 2 + 1;
  r   = N - 1;
  l   = l - 1;
  RR  = list[l - 1];
  K   = list[l - 1];

    while (r != 0) {
      j = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (list[j] > list[j - 1]) j = j + 1;

          if (list[j - 1] > K) {
            list[ i - 1] = list[ j - 1];
          }
          else {
            flag = 0;
          }
        }
      }

      list[ i - 1] = RR;

      if (l == 1) {
        RR  = list [r];

        K = list[r];
        list[r ] = list[0];
        r = r - 1;
      }
      else {
        l   = l - 1;
        RR  = list[ l - 1];
        K   = list[l - 1];
      }
    }

    list[ 0] = RR;

} 

/* ******************************************************************** */
/* sort an integer array                                                */
/* -------------------------------------------------------------------- */

void ML_az_sort(int list[], int N, int list2[], double list3[])
{

  /* local variables */

  int    l, r, RR, K, j, i, flag;
  int    RR2;
  double RR3;

  /*********************** execution begins ******************************/

  if (N <= 1) return;

  l   = N / 2 + 1;
  r   = N - 1;
  l   = l - 1;
  RR  = list[l - 1];
  K   = list[l - 1];

  if ((list2 != NULL) && (list3 != NULL)) {
    RR2 = list2[l - 1];
    RR3 = list3[l - 1];
    while (r != 0) {
      j = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (list[j] > list[j - 1]) j = j + 1;

          if (list[j - 1] > K) {
            list[ i - 1] = list[ j - 1];
            list2[i - 1] = list2[j - 1];
            list3[i - 1] = list3[j - 1];
          }
          else {
            flag = 0;
          }
        }
      }

      list[ i - 1] = RR;
      list2[i - 1] = RR2;
      list3[i - 1] = RR3;

      if (l == 1) {
        RR  = list [r];
        RR2 = list2[r];
        RR3 = list3[r];

        K = list[r];
        list[r ] = list[0];
        list2[r] = list2[0];
        list3[r] = list3[0];
        r = r - 1;
      }
      else {
        l   = l - 1;
        RR  = list[ l - 1];
        RR2 = list2[l - 1];
        RR3 = list3[l - 1];
        K   = list[l - 1];
      }
    }

    list[ 0] = RR;
    list2[0] = RR2;
    list3[0] = RR3;
  }
  else if (list2 != NULL) {
    RR2 = list2[l - 1];
    while (r != 0) {
      j = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (list[j] > list[j - 1]) j = j + 1;

          if (list[j - 1] > K) {
            list[ i - 1] = list[ j - 1];
            list2[i - 1] = list2[j - 1];
          }
          else {
            flag = 0;
          }
        }
      }

      list[ i - 1] = RR;
      list2[i - 1] = RR2;

      if (l == 1) {
        RR  = list [r];
        RR2 = list2[r];

        K = list[r];
        list[r ] = list[0];
        list2[r] = list2[0];
        r = r - 1;
      }
      else {
        l   = l - 1;
        RR  = list[ l - 1];
        RR2 = list2[l - 1];
        K   = list[l - 1];
      }
    }

    list[ 0] = RR;
    list2[0] = RR2;
  }
  else if (list3 != NULL) {
    RR3 = list3[l - 1];
    while (r != 0) {
      j = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (list[j] > list[j - 1]) j = j + 1;

          if (list[j - 1] > K) {
            list[ i - 1] = list[ j - 1];
            list3[i - 1] = list3[j - 1];
          }
          else {
            flag = 0;
          }
        }
      }

      list[ i - 1] = RR;
      list3[i - 1] = RR3;

      if (l == 1) {
        RR  = list [r];
        RR3 = list3[r];

        K = list[r];
        list[r ] = list[0];
        list3[r] = list3[0];
        r = r - 1;
      }
      else {
        l   = l - 1;
        RR  = list[ l - 1];
        RR3 = list3[l - 1];
        K   = list[l - 1];
      }
    }

    list[ 0] = RR;
    list3[0] = RR3;

  }
  else {
    while (r != 0) {
      j = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (list[j] > list[j - 1]) j = j + 1;

          if (list[j - 1] > K) {
            list[ i - 1] = list[ j - 1];
          }
          else {
            flag = 0;
          }
        }
      }

      list[ i - 1] = RR;

      if (l == 1) {
        RR  = list [r];

        K = list[r];
        list[r ] = list[0];
        r = r - 1;
      }
      else {
        l   = l - 1;
        RR  = list[ l - 1];
        K   = list[l - 1];
      }
    }

    list[ 0] = RR;
  }

} 

/* ******************************************************************** */
/* sort a double array                                                  */
/* -------------------------------------------------------------------- */

void ML_az_dsort2(double dlist[], int N, int list2[])
{
  int    l, r, j, i, flag;
  int    RR2;
  double dRR, dK;

  if (N <= 1) return;

  l    = N / 2 + 1;
  r    = N - 1;
  l    = l - 1;
  dRR  = dlist[l - 1];
  dK   = dlist[l - 1];

  if (list2 != NULL) {
     RR2 = list2[l - 1];
     while (r != 0) {
        j = l;
        flag = 1;

        while (flag == 1) {
           i = j;
           j = j + j;

           if (j > r + 1)
              flag = 0;
           else {
              if (j < r + 1)
                 if (dlist[j] > dlist[j - 1]) j = j + 1;

              if (dlist[j - 1] > dK) {
                 dlist[ i - 1] = dlist[ j - 1];
                 list2[i - 1] = list2[j - 1];
              }
              else {
                 flag = 0;
              }
           }
        }
        dlist[ i - 1] = dRR;
        list2[i - 1] = RR2;

        if (l == 1) {
           dRR  = dlist [r];
           RR2 = list2[r];
           dK = dlist[r];
           dlist[r ] = dlist[0];
           list2[r] = list2[0];
           r = r - 1;
         }
         else {
            l   = l - 1;
            dRR  = dlist[ l - 1];
            RR2 = list2[l - 1];
            dK   = dlist[l - 1];
         }
      }
      dlist[ 0] = dRR;
      list2[0] = RR2;
   }
   else {
      while (r != 0) {
         j = l;
         flag = 1;
         while (flag == 1) {
            i = j;
            j = j + j;
            if (j > r + 1)
               flag = 0;
            else {
               if (j < r + 1)
                  if (dlist[j] > dlist[j - 1]) j = j + 1;
               if (dlist[j - 1] > dK) {
                  dlist[ i - 1] = dlist[ j - 1];
               }
               else {
                  flag = 0;
               }
            }
         }
         dlist[ i - 1] = dRR;
         if (l == 1) {
            dRR  = dlist [r];
            dK = dlist[r];
            dlist[r ] = dlist[0];
            r = r - 1;
         }
         else {
            l   = l - 1;
            dRR  = dlist[ l - 1];
            dK   = dlist[l - 1];
         }
      }
      dlist[ 0] = dRR;
   }

}

/******************************************************************************/

void ML_gsum_scalar_int(int vals[], int vals2[], ML_Comm *comm)
{
#ifdef ML_MPI
  MPI_Allreduce((void *) vals,(void *) vals2, 1, MPI_INT, MPI_SUM,
                comm->USR_comm);
  *vals = *vals2;
#endif
  return;
} /* ML_gsum_scalar_int */

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

/**************************************************************************

  For each element in vals[], perform a global sum with the other processors.
  That is, on output vals[i] is equal to the sum of the input values in vals[i]
  on all the processors.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  tvals:           On input, vals[i] on this processor is to be summed with
                   vals[i] on all the other processors.
                   On output, vals[i] is the sum of the input values in val[i]
                   defined on all processors.

  tvals2:          Work space of size 'length'.

  node:            Current processor number.

  nprocs:          Number of processors in the current machine configuration.

  length:          Number of values in 'vals' (i.e. number of global sums).

**************************************************************************/

void ML_gsum_vec_int(int **tvals, int **tvals2, int length, ML_Comm *comm)
{
#ifdef ML_MPI
  int *tmpptr;
  MPI_Allreduce((void *) *tvals,(void *) *tvals2, length, MPI_INT, MPI_SUM,
                comm->USR_comm);
  tmpptr = *tvals;
  *tvals = *tvals2;
  *tvals2 = tmpptr;
#endif
  return;
} /* ML_gsum_vec_int */



/* Just like the  ML_gsum_vec_int but for double vectors */

void ML_gsum_vec_double(double **tvals, double **tvals2, int length, ML_Comm *comm)
{
#ifdef ML_MPI
  double *tmpptr;

  MPI_Allreduce((void *) *tvals,(void *) *tvals2, length, MPI_DOUBLE, MPI_SUM,
                comm->USR_comm);
  tmpptr = *tvals;
  *tvals = *tvals2;
  *tvals2 = tmpptr;
#endif
  return;
} /* ML_gsum_vec_double */

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
void ML_rm_duplicates(int array[], int *N)
{
/*
 * remove any duplicates that might appear in the SORTED
 * array 'array'.
 *
 */
  int k, kk;

  kk = 0;
  for (k = 1; k < *N; k++) {
    if (array[kk] != array[k]) {
      kk++;
      array[kk] = array[k];
    }
  }
  if (*N != 0) kk++;

  *N= kk;
}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************

   Author:        John Shadid, SNL, 1421 Date: 8/1/94
  =======

  Return code:    (void).
  ============

  Parameter list:
  ===============

  num_neighbors:      total number of neighbors to communicate with
  buffer:             on input  - buffer that holds the information to be sent
                      on output - buffer contains the recieve information
  start_send_proc:    contains index of start of information that is to be
                      sent to each processor
  actual_send_length: contains number of double precision entries to be sent to
                      each processor
  actual_recv_length: number of entries to be recieve from each processor
  proc_num_neighbor:  contains the processor number for each neighbor
  type:               the message type to be used in the mesage
  total_num_recv:     on output - total number of actual recvieved entried

**************************************************************************/
void ML_splitup_big_msg(int num_neighbors, char *ibuffer, char *obuffer,
			unsigned int element_size, int *start_send_proc,
                        int *actual_send_length, int *actual_recv_length, int
                        *proc_num_neighbor, int type, int *total_num_recv,
                        ML_Comm *comm)
{

  /*
   * This function handshakes big messages between all the neighbors.  The
   * length of the messages are calculated conservatively to not allow overflow
   * of the message buffers.
   */

  int     m, n, rtype, j, dummy_int;
  int     max_neighbors, messg_size_doubles, doubles_sent;
  int     total_doubles_to_send, dest, messg_from, messg_type;
  int     total_doubles_to_recv, total_send_size;
  int     *finished_send_messg;
  int     *finished_recv_messg;
  int     number_of_messages, *start_recv_proc;
  int     allowed_buff_size, num_recv;
  int     max_buffer_size = 0, max_messg_size;
  char    *send_buffer;
  char   *char_ptr;
  char   *yo = "ML_splitup_big_msg ";
  int     split_up = ML_FALSE;
  int     dummy_add;
  int     debug = ML_FALSE;
  unsigned int length, size;
  
  
  USR_REQ     *request;  /* Message handle */

  /**************************** execution begins ****************************/

  finished_send_messg = (int *) ML_allocate( (num_neighbors+10)*sizeof(int));
  finished_recv_messg = (int *) ML_allocate( (num_neighbors+10)*sizeof(int));
  start_recv_proc     = (int *) ML_allocate( (num_neighbors+10)*sizeof(int));
  request             = (USR_REQ *) ML_allocate( (num_neighbors+10)*sizeof(USR_REQ));
  if ( (request == NULL) || (start_recv_proc == NULL))
     pr_error("ML_splitup_big_msg: out of space\n");

  /* Compute the global maximum message buffer size needed */

  
  for (n = 0; n < num_neighbors; n++) {
    max_buffer_size += actual_recv_length[n];
  }
  max_buffer_size = ML_gmax_int(max_buffer_size, comm);

  /* Determine if splitting of messages is necessary */

  /* Too big for message buffers */
  /*  10/31/03  This code causes setup scaling problems.  You should only
       uncomment it if you're encountering MPI buffer overflow problems. */
  /*
  if (max_buffer_size > (int) (ML_MAX_MSG_BUFF_SIZE / (2 * element_size))) {
     split_up = ML_TRUE;
  }
  */

  if (split_up == ML_TRUE) {

    /*
     * Compute maximum total message size in bytes that any processor will
     * recieve and the maximum number of neighbors that any proc must
     * communicate with. Also initalize some logical arrays.
     */

    max_messg_size = 0;
    for (n = 0; n < num_neighbors; n++) {
      max_messg_size = ML_max(max_messg_size, actual_recv_length[n]);
      finished_send_messg[n] = finished_recv_messg[n] = ML_FALSE;
    }
    max_messg_size = ML_gmax_int(max_messg_size, comm);
    max_neighbors  = ML_gmax_int(num_neighbors,  comm);

    /*
     * Total received nonzeros and starting location for each processors
     * message that will be received
     */

    num_recv = 0;
    for (n = 0; n < num_neighbors; n++) {
      start_recv_proc[n] = num_recv;
      num_recv          += actual_recv_length[n];
    }
    *total_num_recv = num_recv;

    /*
     * Compute the global maximum allowed message size and the maximum number of
     * messages to send all the required information.
     */

    allowed_buff_size  = (int) floor(((double) ML_MAX_MSG_BUFF_SIZE /
                                      (double) (3*element_size)));

    messg_size_doubles = (int) floor((double) allowed_buff_size / 
                                      (double) max_neighbors);

    number_of_messages = (int) ceil((double) max_messg_size /
                                    (double) (messg_size_doubles));

    if (comm->ML_mypid == 0 && debug == ML_TRUE) {
      (void) printf("\n\t\tSplitting up messages in splitup_big_msg\n");
      (void) printf("\t\tmax_buffer_size required  (bytes): %d\n",
                    max_buffer_size*element_size);
      (void) printf("\t\tmax_buffer_size allocated (bytes): %d\n",
                    allowed_buff_size*element_size);
      (void) printf("\t\tindividual message size   (bytes): %d\n",
                    messg_size_doubles*element_size);
      (void) printf("\t\ttotal number of split messages to be sent: %d\n\n",
                    number_of_messages);
    }

    if (ibuffer == obuffer) {
       /*
        * The input and output buffers are the same. Allocate a temporary 
        * send buffer that can hold all out going messages.
        * Then copy all info to this buffer.
        */

        total_send_size = 0;
        for (n = 0; n < num_neighbors; n++) {
           total_send_size += actual_send_length[n];
        }

        send_buffer =(char *) ML_allocate((total_send_size+1)*element_size);
        if (send_buffer == NULL) {
           (void) fprintf(stderr,
                          "no space in ML_splitup_big_msg: send_buffer \n");
           exit(-1);
        }
        for (n = 0; n < (int) (total_send_size*element_size) ; n++)
          send_buffer[n] = ibuffer[n];
    }
    else send_buffer = ibuffer;

    /*
     * Send and receive messages in a series of communications. Each set of
     * exchanges is followed by a syncronization to not allow message buffers to
     * overflow.
     */

    doubles_sent = 0;

    for (m = 0; m < number_of_messages; m++) {
      type++;

      /* post recieves for split messages */

      for (n = 0; n < num_neighbors; n++) {

        total_doubles_to_recv = actual_recv_length[n];
        messg_from            = proc_num_neighbor[n];
        dummy_int             = type;

        if (doubles_sent + messg_size_doubles < total_doubles_to_recv ) {

          /* read messg_size_doubles bytes */

          length = messg_size_doubles*element_size;

          char_ptr = (char *) (&obuffer[start_recv_proc[n]*element_size] + 
                                       doubles_sent*element_size);
          comm->USR_irecvbytes((void *) char_ptr, length, &messg_from, 
                               &dummy_int,  comm->USR_comm, request+n);
        }
        else if (doubles_sent+messg_size_doubles >= total_doubles_to_recv &&
                 finished_recv_messg[n] == ML_FALSE) {

          /* read actual_recv_length[n] - doubles_sent bytes */

          length = (total_doubles_to_recv - doubles_sent)*element_size;

          char_ptr = (char *) (&obuffer[start_recv_proc[n]*element_size] + 
                               doubles_sent*element_size);
          comm->USR_irecvbytes((void *) char_ptr, length, &messg_from, 
                               &dummy_int,  comm->USR_comm, request+n);
        }
        else if (finished_recv_messg[n] == ML_TRUE) {

          /* read integer dummy message */

          length = sizeof(int);
          comm->USR_irecvbytes((void *) &dummy_add, length, &messg_from, 
                                &dummy_int,  comm->USR_comm, request+n);
        }
      }

      /* write split messages */

      for (n = 0; n < num_neighbors; n++) {
        total_doubles_to_send = actual_send_length[n];
        dest                  = proc_num_neighbor[n];

        if (doubles_sent + messg_size_doubles < total_doubles_to_send) {

          /* send out messg_size_doubles bytes */

          length = messg_size_doubles*element_size;
          char_ptr = (char *) (&send_buffer[element_size*start_send_proc[n]] + 
                               doubles_sent*element_size);
          (void) comm->USR_sendbytes((void *) char_ptr, length, dest, type,
				     comm->USR_comm);
        }
        else if (doubles_sent + messg_size_doubles >= total_doubles_to_send &&
                 finished_send_messg[n] == ML_FALSE) {

          /* send out actual_send_length[n] - doubles_sent bytes */

          length = (total_doubles_to_send - doubles_sent)*element_size;

          char_ptr = (char *) (&send_buffer[start_send_proc[n]*element_size] + 
                               doubles_sent*element_size);
          (void) comm->USR_sendbytes((void *) char_ptr, length, dest, type, 
				     comm->USR_comm);

          finished_send_messg[n] = ML_TRUE;
        }
        else if (finished_send_messg[n] == ML_TRUE) {

          /* send out integer dummy message */

          length = sizeof(int);
          (void) comm->USR_sendbytes((void *) &dummy_add, length, dest, type, 
				     comm->USR_comm);
        }
      }

      /* read split messages */

      for (n = 0; n < num_neighbors; n++) {
        total_doubles_to_recv = actual_recv_length[n];
        messg_from            = proc_num_neighbor[n];
        messg_type            = type;

        if (doubles_sent + messg_size_doubles < total_doubles_to_recv ) {

          /* read messg_size_doubles bytes */

          length = messg_size_doubles*element_size;
          char_ptr = (char *) (&obuffer[start_recv_proc[n]*element_size] + 
                               doubles_sent*element_size);
          size =  comm->USR_waitbytes((void *) char_ptr, length, &messg_from,
                             &messg_type, comm->USR_comm, request+n); 

          if (length > size) {
           (void) fprintf(stderr,"%sE4ROR on node %d\nwait failed, message type = %d    %d %d (%d)\n", yo, comm->ML_mypid, 
                          messg_type,length,size,messg_from);
           exit(-1);
          }
        }
        else if (doubles_sent+messg_size_doubles >= total_doubles_to_recv &&
                 finished_recv_messg[n] == ML_FALSE) {

          /* read actual_recv_length[n] - doubles_sent bytes */

          length = (total_doubles_to_recv - doubles_sent)*element_size;
          char_ptr = (char *) (&obuffer[start_recv_proc[n]*element_size] + 
                               doubles_sent*element_size);
          size =  comm->USR_waitbytes((void *) char_ptr, length, &messg_from,
                                      &messg_type, comm->USR_comm, request+n); 

          if (length > size) {
           (void) fprintf(stderr,"%sE3ROR on node %d\nwait failed, message type = %d   %d %d  (%d)\n", yo, comm->ML_mypid, 
                          messg_type,length,size,messg_from);
           exit(-1);
          }

          finished_recv_messg[n] = ML_TRUE;
        }
        else if (finished_recv_messg[n] == ML_TRUE) {

          /* read integer dummy message */

          length = sizeof(int);
          size =  comm->USR_waitbytes((void *) &dummy_add, length, &messg_from,
                             &messg_type, comm->USR_comm, request+n); 

          if (length > size) {
           (void) fprintf(stderr,"%sE2ROR on node %d\nwait failed, message type = %d %d %d (%d)\n", yo, comm->ML_mypid, 
                          messg_type,length,size,messg_from);
           exit(-1);
          }

        }
      }

      doubles_sent += messg_size_doubles;


      j = 3; j = ML_gmax_int(j, comm); /* synchronize procs */
    }

    if (ibuffer == obuffer) {ML_free(send_buffer); send_buffer = NULL;}
    ML_free(request);
    ML_free(start_recv_proc);
    ML_free(finished_recv_messg);
    ML_free(finished_send_messg);
    return;
  }

  else {
     type++;
     
     if (ibuffer == obuffer ) {
        /* Allocate a send buffer, if the input */
        /* and output buffers are the same.     */
        total_send_size = 0;
        for (n = 0; n < num_neighbors; n++) {
           total_send_size += actual_send_length[n];
        }
        send_buffer = (char *) ML_allocate((total_send_size+1)*element_size);
        if (send_buffer == NULL) {
           (void) fprintf(stderr,"no space ML_splitup_big_msg: send_buffer \n");
           exit(-1);
        }
   
        for (n = 0; n < (int) (total_send_size*element_size) ; n++) 
           send_buffer[n] = ibuffer[n];
     }
     else send_buffer = ibuffer;
     
     /* post receives for message */
     
     j = 0;
     for (n = 0; n < num_neighbors; n++) {
        messg_from = proc_num_neighbor[n];
        dummy_int = type;
        size      = actual_recv_length[n]*element_size;

        comm->USR_irecvbytes((void *) &obuffer[j], size, &messg_from, 
                             &dummy_int, comm->USR_comm, request+n);
        j += actual_recv_length[n]*element_size;
     }

     /* send messages to each neighbor */

     for (n = 0; n < num_neighbors; n++) {
        size = actual_send_length[n]*element_size;
        (void) comm->USR_sendbytes((void *) &send_buffer[start_send_proc[n]*
                             element_size], size, proc_num_neighbor[n], type, 
			     comm->USR_comm);
     }             

     /* wait for all messages */

     j = 0;
     for (n = 0; n < num_neighbors; n++) {
        messg_from = proc_num_neighbor[n];
        rtype     = type;
        size      = actual_recv_length[n]*element_size;
        length =  comm->USR_waitbytes((void *) &obuffer[j], size, &messg_from,
                                      &rtype, comm->USR_comm, request+n); 
        if ((length != size) && (size !=0) ) {
           (void) fprintf(stderr, "%sERROR on node %d\nwait failed, message type = %d    %d %d (%d)\n", yo, comm->ML_mypid, 
                          rtype,length,size,messg_from);
           exit(-1);
        }
        j += length;
     }
     *total_num_recv = j/element_size;
     if (ibuffer == obuffer) {ML_free(send_buffer); send_buffer = NULL;}
  }
  ML_free(request);
  ML_free(start_recv_proc);
  ML_free(finished_recv_messg);
  ML_free(finished_send_messg);

} /* ML_splitup_big_msg */

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/**************************************************************************

  Routine to perform dot product of r and z with unit stride. This routine call
  the BLAS routine ddot to do the local vector dot product and then uses the
  global summation routine ML_gsum_double to obtain the reguired global result.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     double, dot product of vectors 'r' and 'z'
  ============

  Parameter list:
  ===============

  N:               Length of vector 'vec'.

  r, z:            Vectors of length 'N'.

**************************************************************************/

double ML_gdot(int N, double r[], double z[], ML_Comm *comm)
{

  static int one = 1;
  int        add_N;

  add_N = N;

  return ML_gsum_double(DDOT_F77(&add_N, r, &one, z, &one), comm);

} /* dot */

/**************************************************************************

  Global double sum.

  Author:
  =======

  Return code:     double, result of global sum.
  ============

  Parameter list:
  ===============

  val:             Individual processor value to be summed.

**************************************************************************/

double ML_gsum_double(double val, ML_Comm *comm)
{
#ifdef ML_MPI
  double val2;          /* arriving value to add */
  MPI_Allreduce((void *) &val,(void *) &val2, 1, MPI_DOUBLE, MPI_SUM,
                comm->USR_comm);
  return val2;
#else
  return val;
#endif
} /* ML_gsum_double */

/**************************************************************************

  Global double max.

  Author:
  =======

  Return code:     double, result of global sum.
  ============

  Parameter list:
  ===============

  val:             Individual processor value to be summed.

**************************************************************************/

double ML_gmax_double(double val, ML_Comm *comm)
{
#ifdef ML_MPI
  double val2;          /* arriving value to add */
  MPI_Allreduce((void *) &val,(void *) &val2, 1, MPI_DOUBLE, MPI_MAX,
                comm->USR_comm);
  return val2;
#else
  return val;
#endif
} /* ML_gmax_double */

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/**************************************************************************

  Global max of type int.

  Author:
  =======

  Return code:     int, maximum value across all processors.
  ============

  Parameter list:
  ===============

  val:             Individual processor value.


**************************************************************************/

int ML_gmax_int(int val, ML_Comm *comm)
{
#ifdef ML_MPI
  int   val2;                     /* arriving value to add */
  MPI_Allreduce((void *) &val,(void *) &val2, 1, MPI_INT, MPI_MAX,
                comm->USR_comm);
  return val2;
#else
  return val;
#endif
} /* ML_gmax_int */

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/**************************************************************************

  Find 'key' in 'list' and return the index number.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     int, -1 = key not found, i = list[i] = key
  ============

  Parameter list:
  ===============

  key:             Element to be search for in list.

  list:            List to be searched.

  length:          Length of list.

**************************************************************************/

int ML_find_index(int key, int list[], int length)
{

  /* local variables */

  int start, end;
  int mid;

  /*********************** execution begins ******************************/
  if (length == 0) return -1;

  start = 0;
  end   = length - 1;

  while (end - start > 1) {
    mid = (start + end) / 2;
    if (list[mid] < key) start = mid;
    else end = mid;
  }

  if (list[start] == key) return start;
  if (list[end] == key)   return end;
  return -1;

} /* ML_find_index */

/******************************************************************************/


int ML_get_random_seed(){
/*******************************************************************************

  Returns the random seed

  Author:          Chris Siefert
  =======
*******************************************************************************/


  return ml_random_seed;
}


void ML_set_random_seed(int seed){
/*******************************************************************************

  Sets the random seed for the ML

  Author:          Chris Siefert
  =======

  Parameter list:
  ===============

  seed:            New seed to use
*******************************************************************************/

  ml_random_seed=seed;
  ml_random_start=0;
}



void ML_random_vec(double u[], int N, ML_Comm *comm) 

/*******************************************************************************

  Set u to a random vector.

  Author:          Ray Tuminaro
  =======

  Parameter list:
  ===============

  u:               On output, vector is initialized to random numbers.

  N:               On input, length of 'u'.

*******************************************************************************/
{

  /* local variables */
  int        i;
  int maxint = 2147483647; /* 2^31 -1 */

  /*********************** BEGIN EXECUTION *********************************/

  /* Distribute the seeds evenly in [1,maxint-1].  This guarantees nothing
   * about where in random number stream we are, but avoids overflow situations
   * in parallel when multiplying by a PID.  It would be better to use
   * a good parallel random number generator. */

  if (ml_random_start) {
    ml_random_seed = (int)((maxint-1) * (1.0 -(comm->ML_mypid+1)/(comm->ML_nprocs+1.0)) );
    ml_random_start = 0;
  }
  if (ml_random_seed < 1 || ml_random_seed == maxint)
    pr_error("ML*ERR* Problem detected in ML_random_vec with seed = %d.\nML*ERR* It should be in the interval [1,2^31-2].\n",ml_random_seed);

  for (i = 0; i < N; i++) u[i] = ML_srandom1(&ml_random_seed);

} /* ML_random_vec */

/******************************************************************************/
double ML_srandom1(int *seed)

/*******************************************************************************
  Random number generator.
  From Park, S. K. and Miller, K. W. 1988. "Random number generators: good ones
  are hard to find." Commun. ACM 31, 10 (Oct. 1988), 1192-1201.

  Parameter list:
  ===============
  seed:            Random number seed.

  Return value:
  ===============
  scalar double in interval (0,1)

*******************************************************************************/
{
  int    a = 16807, m = 2147483647, q = 127773, r = 2836;
  int    lo, hi, test;
  double rand_num;

  /**************************** execution begins ******************************/

  hi   = *seed / q;
  lo   = *seed % q;
  test = a * lo - r * hi;

  if (test > 0) *seed = test;
  else *seed = test + m;

  rand_num = (double) *seed / (double) m;
  return rand_num;

} /* ML_srandom1 */


/* Essentially, like printf , but we exit at the end */
#include <stdarg.h>
void pr_error(char *fmt,  ... ) 
{
  char ml_message_string[800];
  va_list ap;
  va_start(ap, fmt);
  vsprintf(ml_message_string,fmt, ap);
  /*
  char *p;
  int  ival;
  double dval;

  va_start(ap, fmt);
  for (p = fmt; *p; p++) {
     if (*p != '%') { putchar(*p); continue; }
     switch (*++p) {
     case 'd':
        ival = va_arg(ap, int);
        printf("%d",ival);
        break;
     case 'e':
        dval = va_arg(ap,double);
        printf("%e",dval);
        break;
     default:
        putchar(*p);
        break;
     }
  }
  */
  va_end(ap);
  fprintf(stderr,"\n%sn",ml_message_string);
# ifdef ML_MPI
  MPI_Abort(MPI_COMM_WORLD,1);
# endif
  exit(EXIT_FAILURE);
}

#define SERIAL_TYPE 790331
void ML_serial_start(ML_Comm *comm)
{
   int data = 0, type = SERIAL_TYPE, neighbor;
   USR_REQ   request;

   if (comm->ML_mypid != 0) {
       neighbor = comm->ML_mypid-1;
       comm->USR_irecvbytes((void *) &data, sizeof(int), &neighbor, &type,
                comm->USR_comm, &request);
       comm->USR_cheapwaitbytes((void *) &data, sizeof(int), &neighbor, &type,
                        comm->USR_comm, &request);
   }
}
void ML_serial_end(ML_Comm *comm) 
{
   int data = 0, type = SERIAL_TYPE, neighbor;

   neighbor = comm->ML_mypid + 1;
   if (neighbor != comm->ML_nprocs) {
      comm->USR_sendbytes((void *) &data, sizeof(int), neighbor, type,
                          comm->USR_comm);
   }
}

/* ******************************************************************** */
/* Routine that pauses execution, prints out process id's, and allows   */
/* the developer to attach a debugger if desired.                       */
/* (Based on code from ALEGRA).                                         */
/* ******************************************************************** */
#ifndef ICL
#include <unistd.h>
#endif
void ML_BreakForDebugger(ML_Comm *comm)
{
  int i,j;
  int mypid = comm->ML_mypid;
  int nproc = comm->ML_nprocs;
  char buf[80];
  char hostname[80];
  char go = ' ';
  char *str;
  FILE * ML_capture_flag;

  str = (char *) getenv("ML_BREAK_FOR_DEBUGGER");
  i = 0;
  if (str != NULL) i++;

  ML_capture_flag = fopen("ML_debug_now","r");
  if(ML_capture_flag) {
    i++;
    fclose(ML_capture_flag);
  }

  ML_gsum_scalar_int(&i, &j, comm);
  if (i != 0)
  {
    if (mypid == 0) printf("Host and Process Ids for tasks\n");
    for (i = 0; i < nproc; i++) {
      if (i == mypid) {
#if defined(TFLOP) || defined(JANUS_STLPORT) || defined(COUGAR)
        sprintf(buf, "Host: %s   PID: %d", "janus", getpid());
#else
        gethostname(hostname, sizeof(hostname));
        sprintf(buf, "Host: %s   PID: %d (mpi task %d)", hostname, getpid(),mypid);
#endif
        printf("%s\n",buf);
        fflush(stdout);
#ifdef ICL
        Sleep(1);
#else
        sleep(1);
#endif
      }
    }
    if(mypid == 0) {
      printf("\n");
      printf("** Pausing because environment variable ML_BREAK_FOR_DEBUGGER\n");
      printf("** has been set, or file ML_debug_now exists.\n");
      printf("**\n");
      printf("** You may now attach debugger to the processes listed above.\n");
      printf( "**\n");
      printf( "** Enter a character to continue > "); fflush(stdout);
	  scanf("%c",&go);
    }
  }
}

/* Function to sync up processors and execute one line at a time. */
                                                                                
void ML_Pause(ML_Comm *comm)
{
  char go = ' ';

  ML_Comm_Barrier(comm);
                                                                                
  if (comm->ML_mypid == 0) {
      printf( "** Press enter to continue > "); fflush(stdout);
      scanf("%c",&go);
  }
  ML_Comm_Barrier(comm);
}


void ML_use_param(void *data, int junk)
{
  if ( (junk == -365) && (data == NULL)) printf("ML_avoid_unused_param\n");
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void ML_print_line (const char *charstr, int ntimes)

{
  int i;
  for (i = 0; i < ntimes; i++) printf("%c", *charstr);
  printf("\n");
}

/*MS*/
int ML_gsum_int(int val, ML_Comm *comm)
{

  int i;
  ML_gsum_scalar_int(&val, &i,comm);
  return val;
}

int ML_gmin_int(int val, ML_Comm *comm)
{
  int min;
#ifdef HAVE_MPI
  MPI_Allreduce(&val, &min, 1, MPI_INT, MPI_MIN, comm->USR_comm );
#else
  min = val;
#endif
  return min;
}

double ML_gmin_double(double val, ML_Comm *comm)
{
  double min;
#ifdef HAVE_MPI
  MPI_Allreduce(&val, &min, 1, MPI_DOUBLE, MPI_MIN, comm->USR_comm );
#else
  min = val;
#endif
  return min;
}

#include "ml_operator.h"

/* ******************************************************************** */
/* print a ML_Operator into MATLAB format. Only one file is generated   */
/* using global ordering.                                               */
/*                                                                      */
/* Parameter list:                                                      */
/* matrix :             ML_Operator, distributed among processes.       *
 *                      The matrix may be square or rectangular.        *
 *                                                                      *
 * label :              matrix will be written in MATLAB (i,j,k) format *
 *                      to file "label.m". Note that only ONE file will *
 *                      be created that contains the ENTIRE operator.   *
 *                                                                      *
 * global_row_ordering: optional global row numbering.                  *
 * global_col_ordering: optional global column numbering.               *
 *                      If either of these are null, they will be cal-  *
 *                      culated in this function.                       *
 *                                                                      */
/* Albuquerque, 30-Oct-03                                               */
/* ******************************************************************** */

int ML_Operator_Print_UsingGlobalOrdering( ML_Operator *matrix, 
                                           const char label[],
                                           int *global_row_ordering,
                                           int *global_col_ordering)
{

   int    i, j, iproc;
   int    *bindx;
   int    MyPID, NumProc;
   double *val;
   int    allocated, row_length;
   char   filename[80];
   FILE   *fid;
   int    Nrows, NglobalRows, NglobalCols=0;
   int    is_globalRows_allocated = 0;
   int    is_globalCols_allocated = 0;
   int    length=0;
   ML_Comm * comm = matrix->comm;
   
   NglobalRows = matrix->outvec_leng;
   ML_gsum_scalar_int(&NglobalRows,&i, comm);
   NglobalCols = matrix->invec_leng;
   ML_gsum_scalar_int(&NglobalCols,&i, comm);

   if( global_row_ordering == NULL ) {
     length = ML_build_global_numbering(matrix, &global_row_ordering,"rows");
     is_globalRows_allocated = 1;
   }

   if( global_col_ordering == NULL ) {
     if (NglobalRows == NglobalCols)
       global_col_ordering = global_row_ordering; 
     else {
       length = ML_build_global_numbering(matrix, &global_col_ordering,"cols");
       is_globalCols_allocated = 1;
     }
   }
 
   if ( matrix->getrow == NULL) return(1);

   MyPID = comm->ML_mypid;
   NumProc = comm->ML_nprocs;

   allocated = matrix->max_nz_per_row;
   bindx = (int    *)  ML_allocate( allocated*sizeof(int   ));
   val   = (double *)  ML_allocate( allocated*sizeof(double));

   Nrows = matrix->outvec_leng;

   if( label != NULL ) {
     sprintf( filename, "%s.m", label );
     if( MyPID == 0 ) printf("Writing matrix to file %s...\n",filename);
   } else {
     if( MyPID == 0 ) printf("Writing matrix to stdout...\n");
   }

   for ( iproc=0 ; iproc<NumProc ; iproc++ ) {

     if ( MyPID == iproc )
     {
       if( label != NULL ) {
         if( MyPID == 0 ) fid = fopen(filename,"w");
         else             fid = fopen(filename,"a");
       } else {
         fid = stdout;
       }
 
       if( MyPID == 0 ) {
         fprintf(fid,"%%N_global_rows = %d\n", NglobalRows );
         fprintf(fid,"%%N_global_cols = %d\n", NglobalCols );
         fprintf(fid,"%%Number of processors = %d\n", NumProc );
         fprintf(fid,"%% To load this data into Matlab:\n");
         fprintf(fid,"%%    load(filename); A = spconvert(filename);\n");
         fprintf(fid,"%% This ordering may be different than the application's matrix ordering!\n\n");
         fprintf(fid,"%% NOTE: If there are no entries in column %d or row %d,\n",NglobalCols,NglobalRows);
         fprintf(fid,"%% Matlab may get the matrix dimensions wrong.\n");
       }
       
       fprintf( fid,
		"%%Writing data for processor %d\n%%N_rows = %d\n%%outvec_leng = %d\n%%invec_leng = %d\n",
		iproc,
		Nrows,
        matrix->outvec_leng,
        matrix->invec_leng );
       fprintf(fid,"%%\n");
       for (i=0; i<length;i++)
         fprintf(fid,"%%global_row_ordering[%d] = %d\n",i,global_row_ordering[i]);
       fprintf(fid,"%%\n");

       for (i = 0 ; i < Nrows; i++) {
         ML_get_matrix_row(matrix, 1, &i, &allocated, &bindx, &val,
                 &row_length, 0);
         for  (j = 0; j < row_length; j++)
           fprintf(fid,"%d  %d  %20.13e\n",
                   global_row_ordering[i]+1,
                   global_col_ordering[bindx[j]]+1,
                   val[j]);
       }
       if( label != NULL ) fclose(fid);
     } /*if ( MyPID == iproc ) */
#ifdef ML_MPI
     ML_Comm_Barrier( matrix->comm);
#endif
   } /*for ( iproc=0 ; iproc<NumProc ; iproc++ ) */

   /* free memory and return */
   
   fflush(stdout);
   ML_free(val);
   ML_free(bindx);

   if( is_globalRows_allocated == 1 ) ML_free( global_row_ordering );
   if( is_globalCols_allocated == 1 ) ML_free( global_col_ordering );

   return 0;
}

/* ******************************************************************** *
 * Create global numbering for a ML_Operator. I suppose that ML uses a  *
 * linear decomposition among the processes (that is, proc 0 is assigned*
 * the first Nrows elements, and so on). This is enough to define the   *
 * global numbering of local nodes. For the ghost nodes (columns), I use*
 * ML_exchange_bdry.                                                    *
 *
 * This now supports global numbering for either rows or columns, which *
 * may be of interest for rectangular matrices.  If the 3rd input arg   *
 * is either "columns" or "cols", global column numbers are found.      *
 * Otherwise, global row numbers are found.                             *
 *                                                                      *
 * Albuquerque, 30-Oct-03                                               *
 * ******************************************************************** */

int ML_build_global_numbering( ML_Operator *Amat,
			       int **pglobal_numbering,
                   const char *rowsOrCols )
{

  int    i;
  int    Nloc, Nghosts, offset;
  double * dtemp = NULL;
  int * global_numbering;
  ML_Comm *comm = Amat->comm;
  int NglobalRows, NglobalCols;
  int findingRows;

  NglobalRows = Amat->outvec_leng;
  ML_gsum_scalar_int(&NglobalRows,&i, comm);
  NglobalCols = Amat->invec_leng;
  ML_gsum_scalar_int(&NglobalCols,&i, comm);
  
  if (strcmp(rowsOrCols,"cols") == 0 ||  strcmp(rowsOrCols,"columns") == 0) {
    Nloc = Amat->invec_leng;
    findingRows = 0;
  }
  else {
    Nloc = Amat->outvec_leng;
    findingRows = 1;
  }

  if (findingRows && NglobalRows != NglobalCols) {
#ifdef ML_MPI
    MPI_Scan ( &Nloc, &offset, 1, MPI_INT, MPI_SUM,
	       Amat->comm->USR_comm );
    offset -= Nloc;
#else
    offset = 0;
#endif
    global_numbering = (int *)ML_allocate(sizeof(int)*(Nloc+1));
    for (i=0; i<Nloc; i++) global_numbering[i] = i+offset;
    *pglobal_numbering = global_numbering;
    return Nloc;
  }

  if (Amat->getrow->pre_comm == NULL) Nghosts = 0;
  else {
    if (Amat->getrow->pre_comm->total_rcv_length <= 0)
      ML_CommInfoOP_Compute_TotalRcvLength(Amat->getrow->pre_comm);
    Nghosts = Amat->getrow->pre_comm->total_rcv_length;
  }
  
  /* allocate +1 because it is possible that some procs will have
     no rows at all (with ParMETIS) */

  dtemp = (double *) ML_allocate( sizeof(double) * (Nloc+Nghosts+1));
  if( dtemp == NULL )
    pr_error("*ML*ERR* not enough memory to allocated %d bytes\n*ML*ERR* (file %s, line %d)\n",
	     (Nloc+Nghosts) * (int)sizeof(double),
	     __FILE__,
	     __LINE__ );

#ifdef ML_MPI
  MPI_Scan ( &Nloc, &offset, 1, MPI_INT, MPI_SUM,
	     Amat->comm->USR_comm );
  offset -= Nloc;
#else
  offset = 0;
#endif

  /* global numbering for local nodes. Note that ML always
     supposes to have contiguous local nodes (that is, the
     global set of nodes has been subdivided into contiguous
     chunks). This may not be true for the first level
     (ML uses an order which is not the physical one). So, don't
     be surprised that a tridiagonal matrix (before AZ_transform)
     is no longer tridiagonal, for instance... */
    
  for( i=0 ; i<Nloc ; i++ ) dtemp[i] = 1.0*(i+offset);
  for (i=0 ; i<Nghosts; i++) dtemp[i+Nloc] = -1;

  /* I exchange this information using ML_exchange_bdry,
     which is coded for double vectors. */

  ML_exchange_bdry(dtemp,Amat->getrow->pre_comm, Nloc,
		   comm, ML_OVERWRITE,NULL);

  /* allocates memory for global_ordering (+1 as before) */
  
  global_numbering = (int *)ML_allocate(sizeof(int)*(Nloc+Nghosts+1));
       
  if( global_numbering == NULL )
    pr_error("*ML*ERR* not enough memory to allocated %d bytes\n*ML*ERR* (file %s, line %d)\n",
	     (Nloc+Nghosts) * (int)sizeof(int),
	     __FILE__,
	     __LINE__ );

  /* put the received double vectors in the integer vector */

  for( i=0 ; i<Nloc+Nghosts ; i++ )
    global_numbering[i] = (int)dtemp[i];

  *pglobal_numbering = global_numbering;
  
  ML_free( dtemp ); dtemp = NULL;

  return Nloc+Nghosts;
    
} /* ML_build_global_numbering() */

/*******************************************************************************
 * ML_Operator_Lump is intended to create a lumped matrix, e.g., a lumped
 * mass matrix.   It does the obvious thing. Let v = [1,1,...1]' and
 * d = A*v.  Then A'=diag(d) is the lumped version of A.
 *
 * To create the lumped matrix, I do a half clone, then populate the getrow,
 * matvec, and data fields of the new matrix.
 *
 * This works in serial, but hasn't been tested in parallel.
 *
 * JJH, 3/5/2004
 ******************************************************************************/

int ML_Operator_Lump(ML_Operator *A, ML_Operator **B)
{
  double *vin,*vout;
  int mm,nn,i;
  struct ML_CSR_MSRdata *csr_data;

  mm = A->invec_leng;
  nn = A->outvec_leng;
  vin = (double *) ML_allocate(mm * sizeof(double) );
  vout = (double *) ML_allocate((nn+1) * sizeof(double) );
  for (i=0;i<mm;i++) vin[i] = 1.0;
  ML_Operator_Apply(A,mm,vin,nn,vout);

  *B = ML_Operator_halfClone(A);
  (*B)->halfclone = ML_FALSE;
  (*B)->N_nonzeros = nn;
  ML_Operator_Set_Getrow(*B, nn, MSR_getrows);

  csr_data = (struct ML_CSR_MSRdata *)
               ML_allocate(sizeof(struct ML_CSR_MSRdata));

  csr_data->rowptr = NULL;
  csr_data->values = vout;
  csr_data->columns = (int *) ML_allocate((nn+1) * sizeof(int));
  for (i=0; i<nn+1; i++) csr_data->columns[i] = nn+1;

  ML_Operator_Set_ApplyFuncData( *B, mm, nn,
                  csr_data, nn, MSR_matvec, 0);

  ML_free(vin);
  return 0;
}

/*******************************************************************************
 Calculate standard deviation based on the formula

        sigma = ( 1 / (n-1) * sum( (a_i - a)^2 ) )^(0.5)

 where a_i are the samples & a is the average of the a_i's.

     sample     -- value to be analyzed
     n          -- number of samples
     activeflag -- nonzero if this processor is participating
     comm       -- ML communicator

 Note: there are more efficient formulas to calculate this.
*******************************************************************************/


double ML_Global_Standard_Deviation(double sample, int n,
                                    int activeflag, ML_Comm *comm)
{
  double avg = 0.0;
  double sum = 0.0;;

  if (n <= 0) return -999.0;
  if (n == 1) return 0.0;

/*  printf("(%d) sample = %e, n = %d active = %d\n",comm->ML_mypid, sample, n,
   activeflag);*/
  if (activeflag == 0)
     sample = 0.0;
  avg = ML_gsum_double(sample, comm) / n;
  /* printf("(%d) avg = %e\n",comm->ML_mypid, avg); */
  if (activeflag)
    sample -= avg;

  sample = sample*sample;
  sum = ML_gsum_double(sample, comm);
 
  return sqrt(sum / (n-1));
}

/* ************************************************************************ * */

int ML_SetupCoordinates(ML *ml_ptr, int level, int NumPDEEqns,
                        double *in_x_coord, double *in_y_coord,
                        double *in_z_coord)
{
  int    NumDimensions, n, i, Nghost;
  double *x_coord, *y_coord, *z_coord, *tmp;
  ML_Operator* AAA;
  ML_Aggregate_Viz_Stats *grid_info;
  
  NumDimensions  = 0;

  if (!(in_x_coord == 0 && in_y_coord == 0 && in_z_coord == 0))
  {
    grid_info = (ML_Aggregate_Viz_Stats *) ml_ptr->Grid[level].Grid;
    AAA = &(ml_ptr->Amat[level]);

    n = AAA->invec_leng;
    Nghost = 0;

    if (AAA->getrow->pre_comm) 
    {
      if (AAA->getrow->pre_comm->total_rcv_length <= 0)
        ML_CommInfoOP_Compute_TotalRcvLength(AAA->getrow->pre_comm);
      Nghost = AAA->getrow->pre_comm->total_rcv_length;
    }

    tmp = (double *) ML_allocate(sizeof(double) * (Nghost+n));
    for (i = 0 ; i < Nghost + n ; ++i)
      tmp[i] = 0.0;

    n /= NumPDEEqns;
    Nghost /= NumPDEEqns;

    if (in_x_coord) 
    {
      NumDimensions++;
      x_coord = (double *) ML_allocate(sizeof(double) * (Nghost+n));

      for (i = 0 ; i < n ; ++i)
        tmp[i * NumPDEEqns] = in_x_coord[i];

      ML_exchange_bdry(&tmp[0],AAA->getrow->pre_comm, NumPDEEqns * n, 
                       AAA->comm, ML_OVERWRITE,NULL);

      for (i = 0 ; i < n + Nghost ; ++i)
        x_coord[i] = tmp[i * NumPDEEqns];

      grid_info->x = x_coord;
    } /* if (in_x_coord) */

    if (in_y_coord) 
    {
      NumDimensions++;
      y_coord = (double *) ML_allocate(sizeof(double) * (Nghost+n));

      for (i = 0 ; i < n ; ++i)
        tmp[i * NumPDEEqns] = in_y_coord[i];

      ML_exchange_bdry(&tmp[0],AAA->getrow->pre_comm, NumPDEEqns * n, 
                       AAA->comm, ML_OVERWRITE,NULL);

      for (i = 0 ; i < n + Nghost ; ++i)
        y_coord[i] = tmp[i * NumPDEEqns];

      grid_info->y = y_coord;
    } /* if (in_y_coord) */

    if (in_z_coord) 
    {
      NumDimensions++;
      z_coord = (double *) ML_allocate(sizeof(double) * (Nghost+n));

      for (i = 0 ; i < n ; ++i)
        tmp[i * NumPDEEqns] = in_z_coord[i];

      ML_exchange_bdry(&tmp[0],AAA->getrow->pre_comm, NumPDEEqns * n, 
                       AAA->comm, ML_OVERWRITE,NULL);

      for (i = 0 ; i < n + Nghost ; ++i)
        z_coord[i] = tmp[i * NumPDEEqns];

      grid_info->z = z_coord;
    } /* if (in_z_coord) */

    grid_info->Ndim = NumDimensions;
    ML_free(tmp);

  } /* if (!(in_x_coord == 0 && in_y_coord == 0 && in_z_coord == 0)) */

  return(0);
} /* ML_SetupCoordinates() */

/* ************************************************************************ * */

int ML_hash_init(int hash_list[], int hash_length, int *hash_used)
{
  int i;

  for (i = 0; i < hash_length; i++) hash_list[i] = -1;
  *hash_used = 0;
  return 0;
}

/* ************************************************************************ * */

#ifndef ML_UseInlinedHashFunction
void ML_hash_it( int new_val, int hash_list[], int hash_length,int *hash_used,
                 int *index)
{
  int ui;

  ui = new_val<<1;
  if (ui < 0) ui = new_val;
  *index = ui % hash_length;
  while ( hash_list[*index] != new_val) {
    if (hash_list[*index] == -1) { (*hash_used)++; break;}
    (*index)++;
    *index = (*index) % hash_length;
  }
}

/* ************************************************************************ * */
/*
   Important: If you want to use ML_fast_hash, the table size must be 2^k for
   a positive integer k.
*/

void ML_fast_hash(int new_val, int hash_list[], int hlm1,int *hash_used,
                  int *index)
{

  uint32_t ui;

  ui = new_val;
  ml_hash_function(ui);
  /* fast mod because hash_length = 2^k */
  *index = ((int) ui) & hlm1;

  while ( hash_list[*index] != new_val) {
    if (hash_list[*index] == -1) { (*hash_used)++; break;}
    (*index)++;
    *index = (*index) & hlm1;
  }
}
#endif /*ifndef ML_UseInlinedHashFunction*/

/* ************************************************************************* */
/*
  Takes an integer as input, and returns a character array containing a space that
  is the same print length (i.e., the same number of characters on the screen.)
  For example, '-12345' would return '      ' (6 spaces).  Useful for aligning
  prints.
*/
void ML_print_align(int int2match, char *space, int pad)
{
  int jj;
  double dtemp;
  if (int2match < 0) pad++;
  dtemp=fabs((double)int2match);
  for (jj=0; jj<pad; jj++) space[jj]=' ';
  jj=pad;
  while (dtemp>=1) {space[jj++]=' '; dtemp /= 10;}
  if (int2match == 0) space[jj++]=' ';
  space[jj] = '\0';
}
