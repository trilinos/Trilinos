/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions to manipulate the Int_lists data structure.                */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : April, 1998                                          */
/* ******************************************************************** */

#include <stdlib.h>
#include "ml_intlist.h"

/* ******************************************************************** */
/* initialization                                                       */
/* -------------------------------------------------------------------- */

int ML_IntList_Create(ML_IntList **ilistp, int n, int ncnt)
{
   int        nbytes;
   ML_IntList *ilist;

   ML_memory_alloc( (void**) ilistp, sizeof(ML_IntList), "IL1" );
   ilist = (*ilistp);
   ilist->ML_id = ML_ID_ILIST;
   ilist->length = 0;

   if ( n > 0 ) 
   {
      nbytes = ( n + 1 ) * sizeof(int);
      ML_memory_alloc( (void**) &(ilist->start), (unsigned int) nbytes, "IL2" );
      ilist->start[0] = 0;
   } else 
      ilist->start = 0;

   if (ncnt > 0) 
   {
      ilist->cur_mem_leng = ncnt;
      nbytes = ncnt * sizeof(int);
      ML_memory_alloc( (void**) &(ilist->members), (unsigned int) nbytes, "IL3" );
   } 
   else if (n > 0) 
   {
      ilist->cur_mem_leng = n;
      nbytes = n * sizeof(int);
      ML_memory_alloc( (void**) &(ilist->members), (unsigned int) nbytes, "IL4" );
   } 
   else 
   {
      ilist->cur_mem_leng = 0;
      ilist->members      = 0;
   }
   return 0;
}

/* ******************************************************************** */
/* destructor for the list structure                                    */
/* -------------------------------------------------------------------- */

int ML_IntList_Destroy(ML_IntList **ilistp)
{
   ML_IntList *ilist;

   ilist = (*ilistp);
   if ( ilist->ML_id != ML_ID_ILIST ) 
   {
      printf("ML_IntList_Destroy : wrong object. \n");
      exit(1);
   }
   if (ilist->start != 0) ML_memory_free( (void **) &(ilist->start) ); 
   ilist->length = 0;
   if (ilist->members != 0) ML_memory_free( (void **) &(ilist->members) ); 
   ilist->cur_mem_leng = 0;
   ilist->ML_id = -1;
   ML_memory_free( (void **) ilistp ); 
   return 0;
}

/* ******************************************************************** */
/* load a sublist to the structure                                      */
/* -------------------------------------------------------------------- */

int ML_IntList_Load_Sublist(ML_IntList *ilist, int ncnt, int *list)
{
   int   i, begin, end, icnt=0, leng, *itmp, cur_leng;

   if ( ilist->ML_id != ML_ID_ILIST ) 
   {
      printf("ML_IntList_Load_Sublist : wrong object. \n");
      exit(1);
   }
   cur_leng = ilist->cur_mem_leng;
   leng     = ilist->length;
   begin    = ilist->start[leng];
   end      = begin + ncnt;
   if ( end >= cur_leng ) 
   {
      cur_leng = end + (10 * ncnt);
      ilist->cur_mem_leng = cur_leng;
      ML_memory_alloc( (void**) &itmp, cur_leng * sizeof(int), "IL5");
      for ( i = 0; i < begin; i++ ) itmp[i] = ilist->members[i];
      ML_memory_free( (void **) &(ilist->members) ); 
      ilist->members = itmp;
   }
   for ( i = begin; i < end; i++ ) ilist->members[i] = list[icnt++];
   ilist->start[leng+1] = ilist->start[leng] + ncnt;
   ilist->length++;
   return 0;
}

/* ******************************************************************** */
/* get a sublist from the structure                                     */
/* -------------------------------------------------------------------- */

int ML_IntList_Get_Sublist(ML_IntList *ilist, int index, int *ncnt,
                            int *list)
{
   int   i, begin, end;
  
   if ( ilist->ML_id != ML_ID_ILIST ) 
   {
      printf("ML_IntList_Get_Sublist : wrong object. \n");
      exit(1);
   }
   begin = ilist->start[index];
   end   = ilist->start[index+1];
   (*ncnt) = end - begin;
   for ( i = begin; i < end; i++) list[i-begin] = ilist->members[i];
   return 0;
}

/* ******************************************************************** */
/* print the data in the list structure                                 */
/* -------------------------------------------------------------------- */

int ML_IntList_Print(ML_IntList *ilist)
{
   int  i, j;

   if ( ilist->ML_id != ML_ID_ILIST ) 
   {
      printf("ML_IntList_Print : wrong object. \n");
      exit(1);
   }
   printf("int_lists : length = %d \n",ilist->length);
   for ( i = 0; i < ilist->length; i++ ) 
   {
      printf("int_lists : sublist %d \n", i);
      for ( j = ilist->start[i]; j < ilist->start[i+1]; j++ ) 
         printf("int_lists : member %d = %d \n",i,ilist->members[j]);
   }
   return 0;
}

