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
 

/******************************************************************************
   To use this routine, invoke with a call like

       ML_mchk(__FILE__,__LINE__);

*******************************************************************************/

#ifdef ML_JANUS


int heap_info(long *a, long *b, long *c, long *d);

void ML_mchk(char *ptr,int line)
{
 long frags, tfree, lfree, tused;
 long nnodes,me;
 long sourcevec[4],maxvec[4],minvec[4],avgvec[4];
 int i;
 double inverse;
 
 MPI_Comm_size(MPI_COMM_WORLD,&nnodes);
 MPI_Comm_rank(MPI_COMM_WORLD,&me);
 if (heap_info(&frags, &tfree, &lfree, &tused) == -1) {
   printf("get_heap_info failed in mchk\n");
   exit(1);
 }
 
 sourcevec[0] = frags;
/*convert to Kbytes on the fly*/
 sourcevec[1] = tfree/1024;
 sourcevec[2] = lfree/1024;
 sourcevec[3] = tused/1024;
 
 
 MPI_Reduce(sourcevec,maxvec, 4, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD); 
 if (me == 0) {
 
        inverse = 1.0/nnodes;
        for (i =0; i < 4; i++) avgvec[i] = avgvec[i]*inverse;
 
        printf("Heap data at %s line %d \n",ptr, line);
        printf("fragments         : (%8ld -- %8ld -- %8ld)\n",minvec[0],avgvec[0],maxvec[0]);
        printf("free (KB)         : (%8ld -- %8ld -- %8ld)\n",minvec[1],avgvec[1],maxvec[1]);
        printf("largest free (KB) : (%8ld -- %8ld -- %8ld)\n",minvec[2],avgvec[2],maxvec[2]);
        printf("total used  (KB)  : (%8ld -- %8ld -- %8ld)\n",minvec[3],avgvec[3],maxvec[3]);
        printf("\n");
 }
}
 MPI_Reduce(sourcevec,minvec, 4, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
 MPI_Reduce(sourcevec,avgvec, 4, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
