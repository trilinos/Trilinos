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

#ifdef ML_CPP
#include <stdlib.h>
#endif

#include "ml_comm.h"
#include "ml_memory.h"

/*
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
*/

static long malloc_initialized=-1;
static long malloc_leng_log[MAX_MALLOC_LOG];
static long malloc_addr_log[MAX_MALLOC_LOG];
static char malloc_name_log[MAX_MALLOC_LOG][3];
void *ml_void_mem_ptr;

/* ******************************************************************** */
/* memory allocation function                                           */
/* -------------------------------------------------------------------- */

#define ML_FUNCTION_NAME ML_memory_alloc
int ML_memory_alloc( void **memptr, unsigned int leng, char *name )
{
   int  i, *int_ptr, nchunks, ndouble=sizeof(double);
   char *var_ptr;
   double *dptr;

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

/*
TAKING THIS OUT TO HANDLE CASE WHEN THERE ARE NO POINTS ON PROC
   else if (leng == 0)
   { 
      printf("ML_malloc warning : %s - leng = 0 (%d).\n",name,leng);
      (*memptr) = NULL;
      return 0;
   }
*/

   {

      /* -------------------------------------------------------------- */
      /* allocate more space than requested for error checking          */
      /* -------------------------------------------------------------- */

      nchunks = leng / ndouble;
      if ((nchunks * ndouble) < leng) nchunks = nchunks + 3;
      else                            nchunks = nchunks + 2;
      var_ptr = (char *) ML_allocate(nchunks*ndouble);
      dptr = (double *) var_ptr;
      for (i = 0; i < nchunks; i++) dptr[i] = 0.;

      /* -------------------------------------------------------------- */
      /* if an error is returned from malloc, terminate execution       */
      /* -------------------------------------------------------------- */

      if (var_ptr == NULL)
      {
         int mypid=0;
#ifdef ML_MPI
         MPI_Comm_rank(MPI_COMM_WORLD,&mypid);
#endif
         pr_error("(%d) %s: unable to allocate %d bytes to %s.\n",mypid, ML_FUNCTION_NAME,leng, name );
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
}
#ifdef ML_FUNCTION_NAME
#undef ML_FUNCTION_NAME
#endif

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
      ML_free(int_ptr);
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
   void *mem_ptr;

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
            mem_ptr = (void *) malloc_addr_log[i];
            ML_free( mem_ptr );
            malloc_leng_log[i] = -1;
         }
      }
   }
   return 0;
}


int ml_allo_count = 0, ml_free_count = 0;

#ifdef ML_MEM_CHECK
/* sophisticated wrappers for allocating memory */

struct ml_widget {                  /* ml_widget is used to maintain a linked   */
   int order;                    /* list of all memory that was allocated */
   int size;
   char *address;
   struct ml_widget *next;
};
struct ml_widget *ml_widget_head =  NULL;  /* points to first element of allocated */
                                     /* memory linked list.                  */

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

void ML_print_it() {
/*
 * Print out the allocated memory linked list
 */

   struct ml_widget *current;

   current = ml_widget_head;
   while (current != NULL) {
      printf("(%d,%d,%u)\n",current->order, current->size, current->address);
      current = current->next;
   }
}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

char *ML_allocate(unsigned int isize) {

/* 
 *  Allocate memory and record the event by placing an entry in the 
 *  ml_widget_head list. Also recored the size of this entry as well.
 *
 *  Note: we actually allocate more memory than is requested (7 doubles more).
 *  This additional memory is used to record the 'size' and to mark the
 *  memory with a header and trailer which we can later check to see if
 *  they were overwritten.
 *
 */

    char *ptr, *header_start, *header_end;
    struct ml_widget *ml_widget;
    int *size_ptr, i, size;
    double *dptr;

    size = (int) isize;

    size = size + 7*sizeof(double);
    ml_widget = (struct ml_widget *) ML_allocate(sizeof(struct ml_widget));
    if (ml_widget == NULL) return(NULL);
    ptr = (char *) ML_allocate(size);
    if (ptr == NULL) {
       ML_free(ml_widget);
       return(NULL);
    }
    ml_allo_count++;

    /* put trash in the space to make sure nobody is expecting zeros */
    for (i = 0 ; i < size/sizeof(char) ; i++ ) 
       ptr[i] = 'f';


    /* record the entry */

    ml_widget->order = ml_allo_count;
if (size == -7*sizeof(double) ) {
printf("allocating 0 space %u (%d)\n",ptr,size);
 i = 0;
 size = 1/i;
 ml_widget = NULL;
}
    ml_widget->size  = size - 7*sizeof(double);
    ml_widget->next  = ml_widget_head;
    ml_widget_head   = ml_widget;
    ml_widget->address = ptr;

    size_ptr = (int *) ptr;
    size_ptr[0] = size - 7*sizeof(double);
    dptr     = (double *) ptr;

    /* mark the header */

    header_start = (char *) &(dptr[1]);

    for (i = 0 ; i < 3*sizeof(double)/sizeof(char) ; i++ )
       header_start[i] = 'x';

    /* mark the trailer */

    header_end = &(ptr[ (size/sizeof(char)) - 1]);
    header_start = (char *) &(dptr[4]);
    header_start = & (header_start[(size-7*sizeof(double))/sizeof(char)]);

    while (header_start <= header_end) {
       *header_start = 'x';
       header_start++;
    }

    return( (char *) &(dptr[4]) );
}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

void ML_myfree(void *vptr) {
/*
 * Free memory and remove the corresponding entry from the ml_widget_head
 * list. Additionally, check that the size stored in the header is correct
 * and that there has not been any memory overwritten in the header or
 * the trailer.
 *
 */

   struct ml_widget *current, *prev;
   double *dptr;
   int *iptr, size, i;
   char *header_start, *header_end, *ptr;

    ptr = (char *) vptr;
    ml_free_count++;
    if (ptr == NULL) {
       printf("Trying to free a NULL ptr\n");
i = 0;
size = 1/i;
ml_widget_head = NULL;
    }
    else {
       current = ml_widget_head;
       prev    = NULL;
       dptr = (double *) ptr;
       --dptr;
       --dptr;
       --dptr;
       --dptr;
       ptr = (char *) dptr;
       while (current != NULL) {
          if (current->address == ptr) break;
          else { prev = current; current = current->next; }
       }
       if (current == NULL) {
          printf("the pointer %u was not found and thus can not be freed.\n",
                  ptr);
	  size = 0;
	  size = 1/size;
	  sqrt(-23.);
	  while (1 == 1) ;
          exit(1);
       }
       else {
           /* check to see if the header is corrupted */
           iptr = (int *) ptr;
           header_start = (char *) &(dptr[1]);

           for (i = 0 ; i < 3*sizeof(double)/sizeof(char) ; i++ ) {
              if (header_start[i] != 'x') {
                 printf("header is corrupted for %u (%d,%d)\n",ptr,
                         current->size,current->order);
                 size =  0;
                 size = 1/size;
              }
           }
           size = iptr[0];

           /* check to see if the sizes are different */

           if (current->size != size) {
              printf("Freeing %u whose size has changed (%d,%d)\n",
                     current->address,current->size,size);
              exit(1);
           }

           /* check to see if the trailer is corrupted */

           header_end = &(ptr[ ((size+7*sizeof(double))/sizeof(char)) - 1]);
           header_start = (char *) &(dptr[4]);
           header_start = &(header_start[size/sizeof(char)]);

           while (header_start <= header_end) {
              if ( *header_start != 'x') {
                 printf("trailer is corrupted for %u (%d,%d)\n",
                         ptr, size,
                         current->order);
                 size =  0;
                 size = 1/size;
              }
              header_start++;
           }

           /* free the space and the ml_widget */

           ML_free(ptr);
           if (ml_widget_head == current) ml_widget_head = current->next;
           else prev->next = current->next;
           ML_free(current);

       }
   }

}

char *ML_realloc(void *vptr, unsigned int new_size) {

   struct ml_widget *current, *prev;
   int i, *iptr, size, *new_size_ptr;
   char *header_start, *header_end, *ptr;
   char *data1, *data2, *new_ptr, *new_header_start, *new_header_end;
   int newmsize, smaller;
   double *dptr, *new_dptr;

    ptr = (char *) vptr;
    if (ptr == NULL) {
       printf("Trying to realloc a NULL ptr\n");
       exit(1);
    }
    else {
       current = ml_widget_head;
       prev    = NULL;
data1 = ptr;
       dptr = (double *) ptr;
       --dptr;
       --dptr;
       --dptr;
       --dptr;
       ptr = (char *) dptr;
       while (current != NULL) {
          if (current->address == ptr) break;
          else { prev = current; current = current->next; }
       }
       if (current == NULL) {
          printf("the pointer %u was not found and thus can not be realloc.\n",
                  ptr);
          exit(1);
       }
       else {
           /* check to see if the header is corrupted */
           iptr = (int *) ptr;
           header_start = (char *) &(dptr[1]);

           for (i = 0 ; i < 3*sizeof(double)/sizeof(char) ; i++ ) {
              if (header_start[i] != 'x') {
                 printf("realloc header is corrupted for %u (%d,%d)\n",ptr,
                         current->size,current->order);
                 size =  0;
                 size = 1/size;
              }
/* DO WE CHECK THE TRAILER ???? */
           }
           size = iptr[0];


    newmsize = new_size + 7*sizeof(double);
    new_ptr = (char *) ML_allocate(newmsize);
    if (new_ptr == NULL) return(NULL);


    new_size_ptr = (int *) new_ptr;
    new_size_ptr[0] = new_size;
    new_dptr     = (double *) new_ptr;
data2 = (char *) &(new_dptr[4]);

    new_header_start = (char *) &(new_dptr[1]);

    for (i = 0 ; i < 3*sizeof(double)/sizeof(char) ; i++ )
       new_header_start[i] = 'x';

    new_header_end = &(new_ptr[ (newmsize/sizeof(char)) - 1]);
    new_header_start = (char *) &(new_dptr[4]);
    new_header_start= &(new_header_start[new_size/sizeof(char)]);

    while (new_header_start <= new_header_end) {
       *new_header_start = 'x';
       new_header_start++;
    }

    smaller = current->size;
    if (smaller > new_size ) smaller = new_size;
    for (i = 0 ; i < smaller; i++) data2[i] = data1[i];

    ML_free(dptr);
    current->size  = new_size;
    current->address = (char *) new_dptr;
    return( (char *) &(new_dptr[4]));

       }
    }
}
   




void ML_spit_it_out()
{
printf("malloc/free %d %d\n",ml_allo_count,ml_free_count);
if (ml_allo_count != ml_free_count) 
printf("WHOA XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
ML_print_it();
}
#endif
 

/*****************************************************************************
 *  ML_memory_check: 
 *    Print out memory usage information. ML_memory_check's arguments
 *    are used to print an identifying string along with the memory information
 *    and are identical to printf's arguments.
 *
 *    On return, ML_memory_check() returns a pointer to the last string
 *    that ML_memory_check() used when print. 
 *
 *    ML_memory_check(NULL) does not print and only returns a pointer to
 *    the last string used to print. This enables low level routines to
 *    take previous string to ML_memory_check() and append to it.
 *
******************************************************************************/
#define ML_NIntStats 6
#define ML_NDblStats 1
#include "malloc.h"
#include "ml_utils.h"
#include <stdarg.h>
#ifdef ML_TFLOP
int heap_info(long *a, long *b, long *c, long *d);
#endif
char * ML_memory_check(char *fmt, ... )
{
#ifdef ML_MEMORY_CHK
   long int fragments, total_free, largest_free, total_used;
   static double start_time = -1.;
   double elapsed_time;
   int id, nnodes, i;
   long isrcvec[ML_NIntStats],imaxvec[ML_NIntStats],
        iminvec[ML_NIntStats],iavgvec[ML_NIntStats];
   double dsrcvec[ML_NDblStats],dmaxvec[ML_NDblStats],
        dminvec[ML_NDblStats],davgvec[ML_NDblStats];
   static char *ml_memory_label = NULL;
   va_list ap;
#ifndef  ML_TFLOP
   struct mallinfo M;
   int *junk;
   static long unsigned int ml_total_mem = 0;
#endif

  /* allocate space for string that is printed with memory information */

  if (ml_memory_label == NULL) {
    ml_memory_label = (char *) malloc(sizeof(char)*200);
    ml_memory_label[0] = '\0';
  }

  /* if fmt is NULL just return the current string associated with */
  /* the memory printing. The idea is that an low level function   */
  /* can use this to get the string, append any additional info    */
  /* and use this when it invokes this routine a second time.      */
  if (fmt == NULL) return(ml_memory_label);

  /* Take variable argument and transform it to string that will   */
  /* is printed with memory statistics.                            */

  va_start(ap, fmt);
  vsprintf(ml_memory_label,fmt, ap);
  va_end(ap);


  elapsed_time = GetClock();
  if (start_time == -1.) start_time = elapsed_time;
  elapsed_time = elapsed_time - start_time;

#ifdef ML_TFLOP
   heap_info(&fragments, &total_free, &largest_free, &total_used);  
#else
   /* Try to estimate the amount of total space in the system */
   /* by doing mallocs and frees on increasing size chunks.   */
   /* Ideally 'ML_memory_check' should be called in the beginning of   */
   /* the program to calibrate the total amount of avail. mem.*/

   if (ml_total_mem == 0) {
     junk = NULL;
     ml_total_mem = 10000;
     while ( (junk = malloc(sizeof(int)*ml_total_mem)) != NULL) {
       ml_total_mem *= 2;
       free(junk);
     }
     ml_total_mem = ml_total_mem/2;
     ml_total_mem = (int) (((double)ml_total_mem)*1.1);
     while ( (junk = malloc(sizeof(int)*ml_total_mem)) != NULL) {
       ml_total_mem = (int) (((double)ml_total_mem)*1.1);
       free(junk);
     }
     ml_total_mem = (int) (((double)ml_total_mem)/1.1);
     while ( (junk = malloc(sizeof(int)*ml_total_mem)) != NULL) {
       ml_total_mem = ml_total_mem + 1000;
       free(junk);
     }
     ml_total_mem = ml_total_mem - 1000;
     ml_total_mem = ml_total_mem*sizeof(int);
   }

   /* now use system call to get memory used information */

   M = mallinfo();
   fragments = M.ordblks + M.smblks + M.hblks;
   total_free = M.fsmblks + M.fordblks;
   total_used = M.hblkhd + M.usmblks + M.uordblks;
   total_free = ml_total_mem - total_used;
   largest_free = -1024;
#endif

   /* convert to Kbytes */

   total_free   = total_free/1024;
   largest_free = largest_free/1024;
   total_used   = total_used/1024;

   /* Only print if fmt string is not empty */
   /* This allows for an intialization of   */
   /* ml_total_mem without any printing     */
   if (strlen(fmt) == 0)    return(ml_memory_label);


   isrcvec[0] = fragments;
   isrcvec[1] = total_free;
   isrcvec[2] = largest_free;
   isrcvec[3] = total_used;
   isrcvec[4] = total_free+total_used;
   isrcvec[5] = (int) (((double)total_used*1000)/
                          ((double)(total_free+total_used)));
   dsrcvec[0] = elapsed_time;
   nnodes = 1;
   id = 0;
#ifdef ML_MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&id);
   MPI_Comm_size(MPI_COMM_WORLD,&nnodes);
   MPI_Reduce(isrcvec,imaxvec,ML_NIntStats,MPI_LONG,MPI_MAX,0,MPI_COMM_WORLD); 
   MPI_Reduce(isrcvec,iminvec,ML_NIntStats,MPI_LONG,MPI_MIN,0,MPI_COMM_WORLD);
   MPI_Reduce(isrcvec,iavgvec,ML_NIntStats,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
   MPI_Reduce(dsrcvec,dmaxvec,ML_NDblStats,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD); 
   MPI_Reduce(dsrcvec,dminvec,ML_NDblStats,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
   MPI_Reduce(dsrcvec,davgvec,ML_NDblStats,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
   for (i =0; i < ML_NIntStats; i++) {
     imaxvec[i] = isrcvec[i];
     iminvec[i] = isrcvec[i];
     iavgvec[i] = isrcvec[i];
   }
   for (i =0; i < ML_NDblStats; i++) {
     dmaxvec[i] = dsrcvec[i];
     dminvec[i] = dsrcvec[i];
     davgvec[i] = dsrcvec[i];
   }
#endif
   printf("%s(%d): blks = %d, free = %d, max free = %d, used = %d, total = %d, %% used = %e, time = %e\n",
	  ml_memory_label,id,fragments, total_free, largest_free, total_used,
	  total_free+total_used, 
          ((double)total_used)/((double)(total_free+total_used)),elapsed_time);


   if (id == 0)  {
    for (i =0; i < ML_NIntStats; i++) iavgvec[i] = iavgvec[i]/((double) nnodes);
    printf("Summary Heap data at %s\n",ml_memory_label);
    printf("                        avg           min             max\n");
    printf("---------------------------------------------------------\n");
    printf("    blks     %14d %14d %14d %s\n",iavgvec[0], iminvec[0],imaxvec[0],
	   ml_memory_label);
    printf("    free     %14d %14d %14d %s\n",iavgvec[1], iminvec[1],imaxvec[1],
	   ml_memory_label);
    if (iavgvec[2] != -1)
      printf("    max free %14d %14d %14d %s\n",iavgvec[2], iminvec[2],imaxvec[2],
	   ml_memory_label);
    printf("    used     %14d %14d %14d %s\n",iavgvec[3], iminvec[3],imaxvec[3],
	   ml_memory_label);
    printf("    total    %14d %14d %14d %s\n",iavgvec[4], iminvec[4],imaxvec[4],
	   ml_memory_label);
    printf("    %% used   %14.1f %14.1f %14.1f %s\n",
	   ((double)iavgvec[5])/10.,((double)iminvec[5])/10.,((double)imaxvec[5])/10.,
	   ml_memory_label);
    printf("    time   %14.1f %14.1f %14.1f %s\n",
	   davgvec[0],dminvec[0],dmaxvec[0], ml_memory_label);
   }
   return(ml_memory_label);
#else
   return(NULL);
#endif
}

/* returns the maximum allocatable contiguous memory, in Mbytes, using malloc() */
int ML_MaxAllocatableSize()
{
  
  size_t left_size = 1024;
  size_t right_size = 1024*1024*1024;
  size_t max_size;
  void * ptr;
  
  do {

    max_size = (right_size + left_size)/2;
    
    ptr = malloc( max_size );
    if( ptr == 0 )  right_size = max_size;
    else {
      left_size = max_size;
      free( ptr );
    }

    
  } while( right_size - left_size > 16*1024 );

  return (int)(max_size/(1024*1024));
  
}

    
#include "ml_utils.h" 
/* returns the maximum allocatable memory, in Mbytes, using mallinfo() */
int ML_MaxMemorySize()
{ 
  long int fragments, total_free, largest_free, total_used; 
  int percent;
  
#ifndef  ML_TFLOP 
  struct mallinfo M; 
#endif 
  
#ifdef ML_TFLOP 
  heap_info(&fragments, &total_free, &largest_free, &total_used);  
#else
  /* use system call to get memory used information */ 
  
  M = mallinfo(); 
  fragments = M.ordblks + M.smblks + M.hblks; 
  total_free = M.fsmblks + M.fordblks; 
  total_used = M.hblkhd + M.usmblks + M.uordblks; 
  /*  total_free = ml_total_mem - total_used;  */
  largest_free = -1024; 
#endif 
  /* convert to Kbytes */ 
  
  return( (int)(total_used/(1024*1024)) ); 
  
} 
