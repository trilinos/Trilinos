/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* local memory management functions                                    */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : January, 1998                                        */
/* ******************************************************************** */

#include "ml_comm.h"
#include "ml_memory.h"

static long malloc_initialized=-1;
static long malloc_leng_log[MAX_MALLOC_LOG];
static long malloc_addr_log[MAX_MALLOC_LOG];
static char malloc_name_log[MAX_MALLOC_LOG][3];

/* ******************************************************************** */
/* memory allocation function                                           */
/* -------------------------------------------------------------------- */

int ML_memory_alloc( void **memptr, unsigned int leng, char *name )
{
   int  i, *int_ptr, nchunks, ndouble=sizeof(double);
   char *var_ptr;

   /* ----------------------------------------------------------------- */
   /* if the first time, initialize the bookkeeping array               */
   /* ----------------------------------------------------------------- */

   if (malloc_initialized == -1)
   { 
      for (i=0; i<MAX_MALLOC_LOG; i++) malloc_leng_log[i] = -1;
      malloc_leng_log[0] = 0;
      malloc_initialized = 0; 
   }

   /* ----------------------------------------------------------------- */
   /* check for error before allocating any storage space               */
   /* ----------------------------------------------------------------- */

   if (leng < 0)
   {
      printf("ML_malloc error : %s - leng < 0 (%d).\n",name,leng);
      exit(-1);
   }
/*
TAKING THIS OUT TO HANDLE CASE WHEN THERE ARE NO POINTS ON PROC
   else if (leng == 0)
   { 
      printf("ML_malloc warning : %s - leng = 0 (%d).\n",name,leng);
      (*memptr) = NULL;
      return 0;
   }
*/
   else
   {

      /* -------------------------------------------------------------- */
      /* allocate more space than requested for error checking          */
      /* -------------------------------------------------------------- */

      nchunks = leng / ndouble;
      if ((nchunks * ndouble) < leng) nchunks = nchunks + 3;
      else                            nchunks = nchunks + 2;
      var_ptr = (char*) calloc( nchunks, ndouble );

      /* -------------------------------------------------------------- */
      /* if an error is returned from malloc, terminate execution       */
      /* -------------------------------------------------------------- */

      if (var_ptr == NULL)
      {
         printf("ML_malloc error : unable to allocate ");
         printf("%d bytes to %s.\n", leng, name );
         exit(-1);
      }

      /* -------------------------------------------------------------- */
      /* check if there is any available slot in the  bookkeeping array */
      /* -------------------------------------------------------------- */

      for ( i = 1; i < MAX_MALLOC_LOG; i++ ) 
         if (malloc_leng_log[i] == -1) break;

      /* -------------------------------------------------------------- */
      /* if the reply is positive, register the allocation information  */ 
      /* -------------------------------------------------------------- */

      if ( i < MAX_MALLOC_LOG )
      {
         int_ptr    = (int *) var_ptr;
         (*int_ptr) = i + 1;
         int_ptr    = (int*) ((long) var_ptr + nchunks*ndouble - ndouble);
         (*int_ptr) = i + 1;
         malloc_addr_log[i] = (long) memptr;
         malloc_leng_log[i] = nchunks * ndouble;
         malloc_name_log[i][0] = name[0];
         malloc_name_log[i][1] = name[1];
         malloc_name_log[i][2] = name[2];
         var_ptr = (char*) ((long) var_ptr + ndouble);
         (*memptr) = (void *) var_ptr; 
         return i;

      /* -------------------------------------------------------------- */
      /* otherwise, signal the inability to track error for future      */
      /* investigation of this allocation                               */
      /* -------------------------------------------------------------- */

      }
      else
      {
         /*
         if (malloc_initialized != 1) 
            printf("ML_malloc message : unable to track any more %s.\n",
                    name);
         */
         malloc_initialized = 1;
         int_ptr    = (int*) var_ptr;
         (*int_ptr) = -1;
         int_ptr    = (int*) ((long) var_ptr + nchunks*ndouble - ndouble);
         (*int_ptr) = -1;
         var_ptr = (char*) ((long) var_ptr + ndouble);
         (*memptr) = (void *) var_ptr; 
         return 0;
      }
   }
   return 0;
}

/* ******************************************************************** */
/* memory deallocation function                                         */
/* -------------------------------------------------------------------- */

int ML_memory_free(void ** var_ptr)
{
   int  ndouble=sizeof(double), index, index2, *int_ptr;
   char *char_ptr;

   /* ------------------------------------------------------------------ */
   /* Extract head and tail information and check to see if they match.  */
   /* If not, flag error.                                                */
   /* ------------------------------------------------------------------ */

   char_ptr = (char *) (*var_ptr);

   if (char_ptr != NULL)
   {
      int_ptr = (int *) ((long) char_ptr - ndouble);
      index   = (*int_ptr) - 1;
      if ( index >= 0 )
      {
         if (index > MAX_MALLOC_LOG)
         {
            if ( global_comm != NULL )
               printf("%d : ML_memory_free error : header invalid(%d).\n",
                            global_comm->ML_mypid, index);
            else
               printf("ML_memory_free error : header invalid(%d).\n",index);
            exit(-1);
         }
         int_ptr = (int *) ((long) char_ptr + malloc_leng_log[index] - 
                                   2 * ndouble);
         index2   = (*int_ptr);
         if (index != index2-1)
         {
            if ( global_comm == NULL )
               printf("ML_memory_free warning : header/tail mismatch - %d\n",
                       index);
            else
               printf("%d : ML_memory_free warning : header/tail mismatch - %d\n",
                         global_comm->ML_mypid, index);
            printf("   (1) : header,tail indices = %d %d \n",index,index2);
            printf("   (2) : %.3s length = %ld \n", malloc_name_log[index],
                                                   malloc_leng_log[index]);
         }
         /* ########## This check may be messed up by pointer exchanges
         if ( ((long) var_ptr) != malloc_addr_log[index])
         { 
            printf("ML_memory_free warning : \n");
            printf("    %.3s - header and log mismatch.\n",
                               malloc_name_log[index]);
            printf("MEM LOG %d : %d \n", index, (int) malloc_addr_log[index]); 
         }
         ############## */
         malloc_leng_log[index] = -1;
      } 
/*
      else
         printf("ML_memory_free : variable not found.\n");
*/

      int_ptr = (int *) ((long) char_ptr - ndouble);
      free(int_ptr);
   }
   (*var_ptr) = NULL;
   return 0;
}

/* ******************************************************************** */
/* inquire about a variable                                             */
/* -------------------------------------------------------------------- */

int ML_memory_check_var(void* var_ptr)
{
   int  ndouble=sizeof(double), index, index2, *int_ptr;
   char *char_ptr;

   /* ------------------------------------------------------------------ */
   /* Extract head and tail information and check to see if they match.  */
   /* If not, flag error.                                                */
   /* ------------------------------------------------------------------ */

   char_ptr = (char *) var_ptr;

   if (char_ptr != NULL)
   {
      if ( global_comm != NULL )
         printf("%d : ML_memory_check_var : %ld\n",global_comm->ML_mypid,
                (long)var_ptr);
      else
         printf("ML_memory_check_var : %ld\n", (long)var_ptr);

      int_ptr = (int *) ((long) char_ptr - ndouble);
      index   = (*int_ptr) - 1;
      if ( index >= 0 && index < MAX_MALLOC_LOG)
      {
         if ( global_comm != NULL )
         {
            printf("%d : ML_memory_check_var : index, length = %d %d\n",
                   global_comm->ML_mypid, index, (int) malloc_leng_log[index]);
         }
         else
         {
            printf("ML_memory_check_var : index, length = %d %d\n",
                   index, (int) malloc_leng_log[index]);
         }
         if (index > MAX_MALLOC_LOG)
         {
            if ( global_comm != NULL )
               printf("%d : ML_memory_check_var error : header invalid(%d).\n",
                            global_comm->ML_mypid, index);
            else
               printf("ML_memory_check_var error : header invalid(%d)\n",index);
            exit(-1);
         }
         int_ptr = (int *) ((long) char_ptr + malloc_leng_log[index] - 
                                   2 * ndouble);
         index2   = (*int_ptr);
         if (index != index2-1)
         {
            if ( global_comm != NULL )
              printf("%d : ML_memory_check_var warning : header,tail mismatch-%d\n",
                     global_comm->ML_mypid, index);
            else
              printf("ML_memory_check_var warning : header,tail mismatch-%d\n",
                     index);
            printf("   (1) : header,tail indices = %d %d \n",index,index2);
            printf("   (2) : %.3s length = %ld \n", malloc_name_log[index],
                                                   malloc_leng_log[index]);
         }
         /* ########## This check may be messed up by pointer exchanges
         if ( ((long) var_ptr) != malloc_addr_log[index])
         { 
            printf("ML_memory_check_var warning : \n");
            printf("    %.3s - header and log mismatch.\n",
                               malloc_name_log[index]);
            printf("MEM LOG %d : %d \n", index, (int) malloc_addr_log[index]); 
         }
         ############## */
      } 
      else
      {
         if ( global_comm != NULL )
            printf("%d : ML_memory_check_var : invalid index = %d\n", 
                   global_comm->ML_mypid, index);
         else
            printf("ML_memory_check_var : invalid index = %d\n", index);
      }
   }
   return 0;
}

/* ******************************************************************** */
/* memory enquiry function                                              */
/* -------------------------------------------------------------------- */

int ML_memory_inquire()
{
   int  i, isum, icnt;

   if (malloc_initialized == 1)
   {
      printf("ML_memory_inquire : memory usage not available. \n"); 
      return 0;
   }
   else
   {
      isum = icnt = 0;
      for ( i = 0; i < MAX_MALLOC_LOG; i++ ) 
      {
         if ( malloc_leng_log[i] > 0 )
         {
            icnt++; 
            isum += malloc_leng_log[i];
            printf("ML_memory_inquire : %d - %.3s (%ld bytes) is nonempty.\n",
                   i, malloc_name_log[i], malloc_leng_log[i]);
         }
      }
      if ( isum > 0 )
      {
         printf("ML_memory_inquire Summary : \n");
         printf("ML_memory_inquire : %d bytes allocated. \n", isum);
         printf("ML_memory_inquire : %d slots allocated. \n", icnt);
      }
/*
      else 
         printf("ML_memory_inquire Summary : nothing allocated.\n");
*/
      return isum;
   }
}

/* ******************************************************************** */
/* short form of the memory enquiry function                            */
/* -------------------------------------------------------------------- */

int ML_memory_inquire_short(int id)
{
   int  i, isum, icnt;

   if (malloc_initialized == 1)
   {
      printf("ML_memory_inquire : memory usage not available. \n"); 
      return 0;
   }
   else
   {
      isum = icnt = 0;
      for ( i = 0; i < MAX_MALLOC_LOG; i++ ) 
      {
         if (malloc_leng_log[i] > 0)
         {
            icnt++; 
            isum += malloc_leng_log[i];
         }
      }
      printf("%d : ML_memory_inquire : %d bytes allocated.\n", id, isum);
      return isum;
   }
}

/* ******************************************************************** */
/* clean memory of a special variable name                              */
/* -------------------------------------------------------------------- */

int ML_memory_clean( char *name, int inlen )
{
   int i, j, clean_flag, leng;

   leng = inlen;
   if ( inlen > 3 ) leng = 3;
   if ( inlen < 0 ) leng = 0;

   for ( i = 0; i < MAX_MALLOC_LOG; i++ )
   {
      if (malloc_leng_log[i] != -1)
      {
         clean_flag = 0;
         for ( j = 0; j < leng; j++ )
         {
            if ( malloc_name_log[i][j] != name[j] )
            {
               clean_flag = 1; 
               break;
            }
         }
         if ( clean_flag == 0 )
         {
            free( (int *) malloc_addr_log[i] );
            malloc_leng_log[i] = -1;
         }
      }
   }
   return 0;
}

