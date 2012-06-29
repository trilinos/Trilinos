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

#include "ml_utils.h"
#include "ml_comm.h"
#include "ml_memory.h"
#include "string.h"
#include "limits.h"

/*
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
*/

static ml_size_t malloc_initialized=-1;
static ml_size_t malloc_leng_log[MAX_MALLOC_LOG];
static ml_size_t malloc_addr_log[MAX_MALLOC_LOG];
static char malloc_name_log[MAX_MALLOC_LOG][3];
void *ml_void_mem_ptr;

/* ******************************************************************** */
/* memory allocation function                                           */
/* -------------------------------------------------------------------- */

#define ML_FUNCTION_NAME "ML_memory_alloc"
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
      if ((nchunks * ndouble) < (int) leng) nchunks = nchunks + 3;
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
	 pr_error("(%d) %s: unable to allocate %d bytes to %s.\n",mypid, ML_FUNCTION_NAME, leng, name );
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
         int_ptr    = (int*) ((ml_size_t) var_ptr + nchunks*ndouble - ndouble);
         (*int_ptr) = i + 1;
         malloc_addr_log[i] = (ml_size_t) memptr;
         malloc_leng_log[i] = nchunks * ndouble;
         malloc_name_log[i][0] = name[0];
         malloc_name_log[i][1] = name[1];
         malloc_name_log[i][2] = name[2];
         var_ptr = (char*) ((ml_size_t) var_ptr + ndouble);
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
         int_ptr    = (int*) ((ml_size_t) var_ptr + nchunks*ndouble - ndouble);
         (*int_ptr) = -1;
         var_ptr = (char*) ((ml_size_t) var_ptr + ndouble);
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
      int_ptr = (int *) ((ml_size_t) char_ptr - ndouble);
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
         int_ptr = (int *) ((ml_size_t) char_ptr + malloc_leng_log[index] - 
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
         if ( ((ml_size_t) var_ptr) != malloc_addr_log[index])
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

      int_ptr = (int *) ((ml_size_t) char_ptr - ndouble);

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
                (ml_size_t)var_ptr);
      else
         printf("ML_memory_check_var : %ld\n", (ml_size_t)var_ptr);

      int_ptr = (int *) ((ml_size_t) char_ptr - ndouble);
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
         int_ptr = (int *) ((ml_size_t) char_ptr + malloc_leng_log[index] - 
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
         if ( ((ml_size_t) var_ptr) != malloc_addr_log[index])
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
    ml_widget = (struct ml_widget *) malloc(sizeof(struct ml_widget));
                                     /* THIS MALLOC NEEDS TO STAY A */
                                     /* MALLOC AND NOT AN ML_ALLOCATE */
    if (ml_widget == NULL) return(NULL);
    ptr = (char *) malloc(size);     /* THIS MALLOC NEEDS TO STAY A */
                                     /* MALLOC AND NOT AN ML_ALLOCATE */

    if (ptr == NULL) {
       ML_free(ml_widget);           /* THIS FREE() NEEDS TO STAY A */
                                     /* FREE AND NOT AN ML_FREE     */
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
       while (1 == 1) ;
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

           free(ptr);  /* THIS FREE SHOULD NOT BE CHANGE TO ML_FREE */
           if (ml_widget_head == current) ml_widget_head = current->next;
           else prev->next = current->next;
           free(current);  /* THIS FREE SHOULD NOT BE CHANGE TO ML_FREE */

       }
   }

}

char *ML_myrealloc(void *vptr, unsigned int new_size) {

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
    new_ptr = (char *) malloc(newmsize);  /* THIS MALLOC NEEDS TO STAY A */
                                          /* MALLOC AND NOT AN ML_ALLOCATE */
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
 *    Print out memory usage information in megabytes. ML_memory_check's
 *    arguments are used to print an identifying string along with the
 *    memory information and are identical to printf's arguments.
 *
 *    On return, ML_memory_check() returns a pointer to the last string
 *    that ML_memory_check() used when print. 
 *
 *    ML_memory_check(NULL) does not print and only returns a pointer to
 *    the last string used to print. This enables low level routines to
 *    take previous string to ML_memory_check() and append to it.
 *
 *    Ideally 'ML_memory_check' should be called in the beginning of
 *    the program to calibrate the total amount of available memory.
 *
 *    Implementation notes:  All calculations are done in megabytes.
 *    I'm also not sure what happens on a shared memory machine if multiple
 *    processors are reading /proc/meminfo.
 *
******************************************************************************/
#define ML_NIntStats 10
#define ML_NDblStats 2
/* OS X defines malloc in a non-standard place */
#ifdef HAVE_MALLOC_H
#ifndef __APPLE__
#  include <malloc.h>
#else
#  include <sys/malloc.h>
#endif
#else
#  include <stdlib.h>
#endif
#include "ml_utils.h"
#include <stdarg.h>

char * ML_memory_check(char *fmt, ... )
{
#ifdef ML_MEMORY_CHK
   size_t fragments=0;
   int total_free=0, largest_free=0, total_used=0;
   int total_swap=0, total_swap_free=0, total_swap_used=0;
   static double start_time = -1.;
   double elapsed_time;
   int id, nnodes, i;
   ml_IntLoc isrcvec[ML_NIntStats],imaxvec[ML_NIntStats],
        iminvec[ML_NIntStats];
   int isrcvec_copy[ML_NIntStats];
   int iavgvec[ML_NIntStats];
   ml_DblLoc dsrcvec[ML_NDblStats],dmaxvec[ML_NDblStats],
        dminvec[ML_NDblStats];
   double dsrcvec_copy[ML_NDblStats];
   double davgvec[ML_NDblStats];
   static char *ml_memory_label = NULL;
   va_list ap;
#ifdef  ML_TFLOP
   unsigned long ultotal_free=0, ullargest_free=0, ultotal_used=0;
#else
   struct mallinfo M;
   static int ml_total_mem = 0;
#endif
   FILE *fid;
#  define ml_meminfo_size 23
   int haveMemInfo=0, overflowDetected = 0;
   char method[80];
   int mypid=0;

  /* allocate space for string that is printed with memory information */

  if (ml_memory_label == NULL) {
    ml_memory_label = (char *) malloc(sizeof(char)*200);
                                    /* THIS MALLOC NEEDS TO STAY A */
                                    /* MALLOC AND NOT AN ML_ALLOCATE */
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
   /* Memory statistics for Red Storm.  FYI, heap_info returns bytes. */
#ifndef NO_HEAPINFO
   heap_info(&fragments, &ultotal_free, &ullargest_free, &ultotal_used);  
#ifdef ML_MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&mypid);
#endif
   total_free=(int) (ultotal_free / (1024*1024));
   largest_free= (int) (ullargest_free / (1024*1024));
   total_used = (int) (ultotal_used / (1024*1024));
   sprintf(method,"Using heap_info()");
#else
   total_free=0;
   largest_free=0;
   total_used=0;
#endif
#else
   /*
      Memory statistics for all other platforms, via the system call mallinfo()
      and reading file /proc/meminfo, which is available under most Linux OS's.
   */

   M = mallinfo();

   fid = fopen("/proc/meminfo","r");
   if (fid != NULL) {
     char str[80], units[10];
     int k;
     for (i=0; i< ml_meminfo_size; i++) {
       if (fscanf(fid,"%s%d%s", str, &k,units) == 3) {
         if (strcmp(str,"MemTotal:") == 0 && (ml_total_mem==0))
            ml_total_mem = k/1024;
         if (strcmp(str,"MemFree:") == 0)  {total_free = k/1024; }
         if (strcmp(str,"SwapTotal:") == 0)  {total_swap = k/1024; }
         if (strcmp(str,"SwapFree:") == 0)  {total_swap_free = k/1024; }
       }
     }
     fclose(fid);
     total_used = ml_total_mem - total_free;
     total_swap_used = total_swap - total_swap_free;
     sprintf(method,"Using /proc/meminfo");
     haveMemInfo = 1;
   }

   /* If /proc/meminfo doesn't exist, use mallinfo() instead. */
   if ( !haveMemInfo )
   {
     if (ml_total_mem == 0) ml_total_mem = ML_MaxAllocatableSize();
     if (M.hblkhd < 0) { /* try to fix overflow */
       double delta = fabs(((double) INT_MIN) - ((double) M.hblkhd)) + 1;
       total_used = (int) ( (((double) INT_MAX) + delta) / (1024*1024) );
       overflowDetected = 1;
     }
     /*Ignore this field upon overflow because I'm don't know how to handle it*/
     if (M.uordblks > 0) total_used += M.uordblks / (1024*1024);
     total_free = ml_total_mem - total_used;
     sprintf(method,"Using mallinfo()");
   }
   fragments = M.ordblks + M.hblks;

   largest_free = -1;
#endif /*ifdef ML_TFLOP*/

   /* Only print if fmt string is not empty */
   /* This allows for an intialization of   */
   /* ml_total_mem without any printing     */
   if (strlen(fmt) == 0)    return(ml_memory_label);


   /*isrcvec[0].value = fragments; */
   isrcvec[0].value = 0;
   isrcvec[1].value = total_free;
   isrcvec[2].value = largest_free;
   isrcvec[3].value = total_used;
   isrcvec[4].value = total_free + total_used;
   /*TODO could this overflow?*/
   isrcvec[5].value = (int) ( ((double)total_used*1000) /
                        ((double)(total_free+total_used)) );
   isrcvec[6].value = total_swap_free;
   isrcvec[7].value = total_swap_used;
   isrcvec[8].value = total_swap;
   /*TODO could this overflow?*/
   isrcvec[9].value = (int) ( ((double)total_swap_used*1000) /
                        ((double)(total_swap)) );
   dsrcvec[0].value = elapsed_time;
   dsrcvec[1].value = fragments;

#ifdef ML_MPI
   for (i =0; i < ML_NIntStats; i++)
      MPI_Comm_rank(MPI_COMM_WORLD,&(isrcvec[i].rank));
   for (i =0; i < ML_NDblStats; i++)
      MPI_Comm_rank(MPI_COMM_WORLD,&(dsrcvec[i].rank));
#endif

   for (i =0; i < ML_NIntStats; i++) isrcvec_copy[i] = isrcvec[i].value;
   for (i =0; i < ML_NDblStats; i++) dsrcvec_copy[i] = dsrcvec[i].value;

   nnodes = 1;
   id = 0;
#ifdef ML_MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&id);
   MPI_Comm_size(MPI_COMM_WORLD,&nnodes);
   MPI_Reduce(isrcvec,imaxvec,ML_NIntStats,MPI_2INT,MPI_MAXLOC,0,MPI_COMM_WORLD); 
   MPI_Reduce(isrcvec,iminvec,ML_NIntStats,MPI_2INT,MPI_MINLOC,0,MPI_COMM_WORLD);
   MPI_Reduce(isrcvec_copy,iavgvec,ML_NIntStats,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
   MPI_Reduce(dsrcvec,dmaxvec,ML_NDblStats,MPI_DOUBLE_INT,MPI_MAXLOC,0,MPI_COMM_WORLD); 
   MPI_Reduce(dsrcvec,dminvec,ML_NDblStats,MPI_DOUBLE_INT,MPI_MINLOC,0,MPI_COMM_WORLD);
   MPI_Reduce(dsrcvec_copy,davgvec,ML_NDblStats,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   MPI_Reduce(&overflowDetected,&i,1,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
   overflowDetected = i;
#else
   for (i =0; i < ML_NIntStats; i++) {
     imaxvec[i].value = isrcvec[i].value;
     iminvec[i].value = isrcvec[i].value;
     iavgvec[i] = isrcvec[i].value;
   }
   for (i =0; i < ML_NDblStats; i++) {
     dmaxvec[i].value = dsrcvec[i].value;
     dminvec[i].value = dsrcvec[i].value;
     davgvec[i] = dsrcvec[i].value;
   }
#endif
/* uncomment lines below if you want individual processor information */
/*
   printf("%s(%d): blks = %ld, free = %ld, max free = %ld, used = %ld, total = %ld, %% used = %e, time = %e\n",
	  ml_memory_label,id,fragments, total_free, largest_free, total_used,
	  total_free+total_used, 
          ((double)total_used)/((double)(total_free+total_used)),elapsed_time);
*/

   if (id == 0 && ML_Get_PrintLevel() > 0) {
     for (i =0; i < ML_NIntStats; i++)
       iavgvec[i] = (int) (iavgvec[i]/((double) nnodes));
     for (i =0; i < ML_NDblStats; i++)
       davgvec[i] = davgvec[i] / nnodes;
     printf("-------------------------------------------------------------\n");
     printf("Summary Heap data (Mbytes) at %s\n",ml_memory_label);
     printf("%s\n",method);
     if (overflowDetected)
        printf("*WARNING* mallinfo() counter overflow detected\n");
     printf("                       avg           min             max\n");
     printf("-------------------------------------------------------------\n");
     printf(" blks       %11d %11d (%5d) %11d (%5d) %s\n",
            (int) davgvec[1], (int) dminvec[1].value, dminvec[1].rank,
            (int) dmaxvec[1].value, dmaxvec[1].rank, ml_memory_label);
     printf(" free       %11d %11d (%5d) %11d (%5d) %s\n",
            iavgvec[1], iminvec[1].value, iminvec[1].rank,
            imaxvec[1].value, imaxvec[1].rank, ml_memory_label);
     if (iavgvec[2] != -1)
       printf(" max free   %11d %11d (%5d) %11d (%5d) %s\n",
              iavgvec[2], iminvec[2].value, iminvec[2].rank,
              imaxvec[2].value, imaxvec[2].rank, ml_memory_label);
     printf(" used       %11d %11d (%5d) %11d (%5d) %s\n",
              iavgvec[3], iminvec[3].value, iminvec[3].rank,
              imaxvec[3].value, imaxvec[3].rank, ml_memory_label);
     printf(" total      %11d %11d (%5d) %11d (%5d) %s\n",
              iavgvec[4], iminvec[4].value, iminvec[4].rank,
              imaxvec[4].value, imaxvec[4].rank, ml_memory_label);
     printf(" %% used       %9.1f   %9.1f (%5d)   %9.1f (%5d) %s\n",
            ((double)iavgvec[5])/10., ((double)iminvec[5].value)/10.,
            iminvec[5].rank,
            ((double)imaxvec[5].value)/10., imaxvec[5].rank, ml_memory_label);
     printf(" time         %9.1f   %9.1f (%5d)   %9.1f (%5d) %s\n",
            davgvec[0],dminvec[0].value,dminvec[0].rank,
            dmaxvec[0].value, dmaxvec[0].rank, ml_memory_label);
     if (haveMemInfo) {
       printf(" swap free  %11d %11d (%5d) %11d (%5d) %s\n",
              iavgvec[6], iminvec[6].value,iminvec[6].rank,
              imaxvec[6].value, iminvec[6].rank, ml_memory_label);
       printf(" swap used  %11d %11d (%5d) %11d (%5d) %s\n",
                iavgvec[7], iminvec[7].value, iminvec[7].rank,
                imaxvec[7].value, imaxvec[7].rank, ml_memory_label);
       printf(" total swap %11d %11d (%5d) %11d (%5d) %s\n",
                iavgvec[8], iminvec[8].value, iminvec[8].rank,
                imaxvec[8].value, imaxvec[8].rank, ml_memory_label);
       printf(" %% swap used  %9.1f   %9.1f (%5d)   %9.1f (%5d) %s\n",
              ((double)iavgvec[9])/10., ((double)iminvec[9].value)/10.,
              iminvec[9].rank,
              ((double)imaxvec[9].value)/10., imaxvec[9].rank, ml_memory_label);
     }
   } /*if (id == 0 ... */
   return(ml_memory_label);
#else
   return(NULL);
#endif
} /*ML_memory_check*/

/** ********************************************************************** **/

/* Returns an estimate (in megabytes) of the amount of total memory */
/* available on the system by using malloc.                         */
/*                                                                  */
/* A major caveat is that there is no way to detect swap space.     */
/* Experience has shown that this approach is unreliable, as malloc */
/* may very well allocate memory that it cannot back.               */
/*                                                                  */
/* Start by doubling memory requests until reaching a size that     */
/* cannot be malloc'd. Then use a binary search to get close to     */
/* the upper limit.  Finally, step by 4 byte increments to get the  */
/* maximum size that can be malloc'd.                               */
/*                                                                  */
int ML_MaxAllocatableSize()
{
  int* junk = NULL;
  long long int upper, mid, lower;
  long long int ml_total_mem;

  ml_total_mem = 10000;
  while ( (junk = (int*) malloc(sizeof(int)*ml_total_mem)) != NULL) {
    ml_total_mem *= 2;
    free(junk);
  }
  ml_total_mem = ml_total_mem/2;

  lower = ml_total_mem;
  upper = ml_total_mem * 2;
  while (upper > (lower+10)) {
    mid   = (lower+upper)/2;
    junk = (int*) malloc(sizeof(int)*(mid));
    if (junk == NULL)
      upper = mid;
    else {
      lower = mid; 
      free(junk);
    }
  }
  ml_total_mem = lower;

  /* step by 4 bytes the rest of the way */
  while ( (junk = (int*) malloc(sizeof(int)*ml_total_mem)) != NULL) {
    ml_total_mem = ml_total_mem + 1;
    free(junk);
  }
  ml_total_mem = ml_total_mem - 1;
  return (int)( ml_total_mem * sizeof(int) / (1024*1024) );

} /*ML_MaxAllocatableSize()*/

/** ********************************************************************** **/
    
#include "ml_utils.h" 
#ifdef ML_MALLINFO
#ifndef __APPLE__
#include <malloc.h>
#else
#include <sys/malloc.h>
#endif
#endif
/* returns the maximum allocatable memory, in Mbytes, using mallinfo() */
int ML_MaxMemorySize()
{
#ifdef ML_MALLINFO
  size_t fragments;
  long long int total_free, largest_free, total_used; 
  /* int percent; */
  
#ifndef  ML_TFLOP 
  struct mallinfo M; 
#endif 
#ifdef  ML_TFLOP
   unsigned long ultotal_free=0, ullargest_free=0, ultotal_used=0;
#endif 
  
#ifdef ML_TFLOP 
#ifndef NO_HEAPINFO
  heap_info(&fragments, &ultotal_free, &ullargest_free, &ultotal_used);  
  total_free=(int) (ultotal_free / 1024);
  largest_free= (int) (ullargest_free / 1024);
  total_used = (int) (ultotal_used /1024);
#else
  total_free=0;
  largest_free=0;
  total_used=0;
#endif
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
#else
  return -1;
#endif
}
