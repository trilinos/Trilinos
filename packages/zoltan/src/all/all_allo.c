/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 * $Name$
 *====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "mpi.h"
#include "lb_const.h"
#include "all_allo_const.h"
#include "params_const.h"

#ifdef __STDC__
#include <stdarg.h>
#else
#include <varargs.h>
#endif

static int DEBUG_MEMORY = 0;	/* Flag for detecting memory leaks */
static int bytes_used = 0;	/* Sum of active allocations */
static int bytes_max = 0;	/* Largest total of active allocations */

static int nmalloc = 0;         /* number of calls to malloc */
static int nfree = 0;           /* number of calls to free */

static struct malloc_debug_data {
  int       order;			/* which malloc call is it? */
  int       size;			/* size of malloc invocation */
  double   *ptr;			/* memory location returned */
  char  file[MAX_PARAM_STRING_LEN+1];   /* file name */
  int       line;                       /* line number */
  struct malloc_debug_data *next;	/* pointer to next element */
} *top = NULL;

/* Fortran memory allocation callback functions */

static LB_FORT_MALLOC_INT_FN *LB_Fort_Malloc_int;
static LB_FORT_MALLOC_GID_FN *LB_Fort_Malloc_GID;
static LB_FORT_MALLOC_LID_FN *LB_Fort_Malloc_LID;
static LB_FORT_FREE_INT_FN *LB_Fort_Free_int;
static LB_FORT_FREE_GID_FN *LB_Fort_Free_GID;
static LB_FORT_FREE_LID_FN *LB_Fort_Free_LID;


int LB_Set_Malloc_Param(
char *name,			/* name of variable */
char *val)			/* value of variable */
{
    int status;
    PARAM_UTYPE result;		/* value returned from Check_Param */
    int index;			/* index returned from Check_Param */
    PARAM_VARS malloc_params[] = {
	{ "DEBUG_MEMORY", &DEBUG_MEMORY, "INT" },
	{ NULL, NULL, NULL }
    };

    status = LB_Check_Param(name, val, malloc_params, &result, &index);
    if (status == 0 && index == 0) {
	DEBUG_MEMORY = result.ival;
	status = 3;
    }

    return(status);
}


/******************************************************************************
 *
 *                    Dynamic Allocation of Multidimensional Arrays
 *-----------------------------------------------------------------------------
 *
 * Example Usage:
 *
 *     typedef  struct
 *       {      int     bus1;
 *              int     bus2;
 *              int     dest;
 *      }       POINT;
 *
 *      POINT    **points, corner;
 *
 *      points = (POINT **)LB_Array_Alloc(file, lineno, 2, x, y, sizeof(POINT));
 *                                        ^     ^       ^  ^  ^
 *                                        |     |       |  |  |
 *                 name of calling file---*     |       |  |  |
 *                                              |       |  |  |
 *                  line number of call---------*       |  |  |
 *                                                      |  |  |
 *                 number of dimensions-----------------+  |  |
 *                                                         |  |
 *                  first dimension max--------------------+  |
 *                                                            |
 *                 second dimension max-----------------------+
 *
 *         (points may be now be used as if it were declared
 *          POINT points[x][y])
 *
 *          This particular version is limited to dimensions of 4 or less.
 *
 *      corner = points[2][3]; (refer to the structure as you would any array)
 *
 *      free (points); (frees the entire structure in one fell swoop)
 *
 *****************************************************************************/
/******************************************************************************
*       The following section is a commented section containing
*       an example main code:
*******************************************************************************
*double *LB_Array_Alloc();
*main()
*{
*  int ***temp;
*   int *temp2;
*   int i, j, k;
*   int il, jl, kl;
*
*   malloc_debug(2);
*   il = 2;
*   jl = 3;
*   kl = 3;
*   temp = (int ***) Array_Alloc(__FILE__,__LINE__,3,il,jl,kl,sizeof(int));
*   for (i=0; i<il; i++) {
*      for (j=0; j<jl; j++) {
*         for (k=0; k<kl; k++) temp[i][j][k] = 1;
*      }
*   }
*
*   temp2 = (int *) malloc(10*sizeof(int));
*   for (i=0; i<10; i++) temp2[i] = 0;
*
*   for (i=0; i<il; i++) {
*      for (j=0; j<jl; j++) {
*         for (k=0; k<kl; k++) (void) printf(" %d\n", temp[i][j][k]);
*      }
*   }
*   malloc_verify();
*}
******************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#ifdef __STDC__

double *LB_Array_Alloc(char *file, int lineno, int numdim, ...)

#else

double *LB_Array_Alloc(va_alist)
va_dcl

#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

{
  char *yo = "LB_Array_Alloc";
  int i, j;
  struct dimension {
    long index;  /* Number of elements in the dimension  */
    long total;  /* Total number of elements             */
    long size;   /* Size of a single element in bytes    */
    long off;    /* offset from beginning of array       */
  } dim[4];      /* Info about each dimension            */

#ifndef __STDC__
  char *file;           /* Filename of source code from call.   */
  int lineno;           /* Line number of call.                 */
  int numdim;           /* Number of dimensions                 */
#endif

  long    total;        /* Total size of the array              */
  double *dfield;       /* ptr to avoid lint complaints         */
  char   *field;        /* The multi-dimensional array          */
  char  **ptr;          /* Pointer offset                       */
  char   *data;         /* Data offset                          */
  va_list va;           /* Current pointer in the argument list */

#ifdef __STDC__
  va_start(va, numdim);
#else
  va_start(va);
  file = va_arg(va, char *);
  lineno = va_arg(va, int);
  numdim = va_arg(va, int);
#endif

  if (numdim <= 0) {
    fprintf(stderr, "%s (%s: %d) ERROR: number of dimensions, %d, is <=0\n",
            yo, file, lineno, numdim);
    return((double *) NULL);
  }
  else if (numdim > 4) {
    fprintf(stderr, "%s (%s: %d) ERROR: number of dimensions, %d, is > 4\n",
            yo, file, lineno, numdim);
    return((double *) NULL);
  }

  dim[0].index = va_arg(va, int);

  if (dim[0].index <= 0) {
    if (DEBUG_MEMORY > 0) {
      fprintf(stderr, "WARNING, %s (%s: %d) called with first "
            "dimension <= 0, %ld; will return NULL\n", 
            yo, file, lineno, dim[0].index);
    }
    return((double *) NULL);
  }

  dim[0].total = dim[0].index;
  dim[0].size  = sizeof(void *);
  dim[0].off   = 0;
  for (i = 1; i < numdim; i++) {
    dim[i].index = va_arg(va, int);
    if (dim[i].index <= 0) {
      fprintf(stderr, "WARNING: %s (%s: %d) called with dimension %d <= 0, "
              "%ld; will return NULL\n",
              yo, file, lineno, i+1, dim[i].index);
      return((double *) NULL);
    }
    dim[i].total = dim[i-1].total * dim[i].index;
    dim[i].size  = sizeof(void *);
    dim[i].off   = dim[i-1].off + dim[i-1].total * dim[i-1].size;
  }

  dim[numdim-1].size = va_arg(va, int);
  va_end(va);

  /* Round up the last offset value so data is properly aligned. */

  dim[numdim-1].off = dim[numdim-1].size *
    ((dim[numdim-1].off+dim[numdim-1].size-1)/dim[numdim-1].size);

  total = dim[numdim-1].off + dim[numdim-1].total * dim[numdim-1].size;

  dfield = (double *) LB_Malloc((int) total, file, lineno); 

  if (dfield != NULL) {
    field  = (char *) dfield;
    for (i = 0; i < numdim - 1; i++) {
      ptr  = (char **) (field + dim[i].off);
      data = (char *) (field + dim[i+1].off);
      for (j = 0; j < dim[i].total; j++) {
        ptr[j] = data + j * dim[i+1].size * dim[i+1].index;
      }
    }
  }

  return dfield;

} /* LB_Array_Alloc */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/


/* Safe version of malloc.  Does not initialize memory .*/

double *LB_Malloc(int n, char *filename, int lineno)

{
  char *yo = "LB_Malloc";
  struct malloc_debug_data *new_ptr;     /* data structure for malloc data */
  int       proc;             /* processor ID for debugging msg */
  double *pntr;           /* return value */

  if (n > 0) {
    pntr = (double *) malloc((unsigned) n);
    if (pntr == NULL) {
      MPI_Comm_rank(MPI_COMM_WORLD, &proc);
      fprintf(stderr, "%s (from %s,%d) No space on proc %d - number of bytes "
              "requested = %d\n", yo, filename, lineno, proc, n);
      return ((double *) NULL);
    }
    nmalloc++;
  }
  else if (n == 0)
    pntr = NULL;
  else {		/* n < 0 */
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    fprintf(stderr, "%s (from %s,%d) ERROR on proc %d: "
	    "Negative malloc argument. (%d)\n", yo, filename, lineno, proc, n);
    return ((double *) NULL);
  }

  if (DEBUG_MEMORY > 1 && pntr != NULL) {
    /* Add allocation to linked list of active allocations. */
    new_ptr = (struct malloc_debug_data *)
      malloc(sizeof(struct malloc_debug_data));

    if (new_ptr == NULL) {
      MPI_Comm_rank(MPI_COMM_WORLD, &proc);
      fprintf(stderr, "WARNING: No space on proc %d for malloc_debug %d.\n",
	proc, n);
      return (pntr);
    }

    new_ptr->order = nmalloc;
    new_ptr->size = n;
    new_ptr->ptr = pntr;
    strncpy(new_ptr->file, filename, MAX_PARAM_STRING_LEN);
    new_ptr->line = lineno;
    new_ptr->next = top;
    top = new_ptr;
    bytes_used += n;
    if (bytes_used > bytes_max) {
      bytes_max = bytes_used;
    }
  }

  if (DEBUG_MEMORY > 2) {
    /* Print out details of allocation. */
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    fprintf(stderr, "Proc %d: order=%d, size=%d, location=0x%lx, file=%s, line=%d\n",
      proc, nmalloc, n, (long) pntr, filename, lineno);
  }


  return pntr;

} /* LB_Malloc */

/* Safe version of realloc. Does not initialize memory. */

double *LB_Realloc(void *ptr, int n, char *filename, int lineno)
{
  char *yo = "LB_Realloc";
  struct malloc_debug_data *dbptr;   /* loops through debug list */
  int       proc;             /* processor ID */
  double   *p;                /* returned pointer */

  if (ptr == NULL) {	/* Previous allocation not important */
    if (n == 0) {
      p = NULL;
    }
    else {
      p = LB_Malloc(n, filename, lineno);
    }
  }
  else {
    if (n == 0) {
      LB_Free((void **) &ptr, filename, lineno);
      p = NULL;
    }
    else {
      p = (double *) realloc((char *) ptr, n);

      if (DEBUG_MEMORY > 1) {
        /* Need to replace item in allocation list */
        for (dbptr = top; dbptr != NULL && (char *) dbptr->ptr != ptr;
	   dbptr = dbptr->next);
	if (dbptr == NULL) {	/* previous allocation not found in list. */
	   MPI_Comm_rank(MPI_COMM_WORLD, &proc);
	   fprintf(stderr, "Proc %d: Memory error: "
	     "In realloc, address not found in debug list (0x%lx)\n",
	     proc, (long) ptr);
	}
	else {	/* Update entry in allocation list */
	  bytes_used += n - dbptr->size;
	  dbptr->size = n;
	  dbptr->ptr = p;
	  if (bytes_used > bytes_max) {
	    bytes_max = bytes_used;
	  }
	}
      }

      if (p == NULL) {
        MPI_Comm_rank(MPI_COMM_WORLD, &proc);
        fprintf(stderr, "%s (from %s,%d) No space on proc %d - "
		"number of bytes requested = %d\n",
		yo, filename, lineno, proc, n);
      }
    }
  }

  return (p);
} /* LB_Realloc */


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_Free (void **ptr, char *filename, int lineno)
{
  struct malloc_debug_data *dbptr;   /* loops through debug list */
  struct malloc_debug_data **prev;   /* holds previous pointer */
  int       proc;             /* processor ID */

/*
 *  This version of free calls the system's free function.  It doesn't call
 *  free if ptr is the NULL pointer.
 */
 
  if (ptr == NULL || *ptr == NULL) 
    return;

  nfree++;

  if (DEBUG_MEMORY > 1) {
    /* Remove allocation of list of active allocations */
    prev = &top;
    for (dbptr = top; dbptr != NULL && (char *) dbptr->ptr != *ptr;
      dbptr = dbptr->next) {
      prev = &(dbptr->next);
    }
    if (dbptr == NULL) {
      MPI_Comm_rank(MPI_COMM_WORLD, &proc);
      fprintf(stderr, "Proc %d: Memory error: In free, address (0x%lx) "
	"not found in debug list. File=%s, line=%d.\n", proc, 
        (long) *ptr, filename, lineno);
   }
   else {
       *prev = dbptr->next;
       bytes_used = -dbptr->size;
       free((char *) dbptr);
       }
   }
 
  free(*ptr);
 
  /* Set value of ptr to NULL, to flag further references to it. */
  *ptr = NULL;

}  /* LB_Free */


#if 0 /* Not sure we need this feature */

/* Free n pointers. Variable number of arguments is allowed.  */

#ifdef __STDC__

void LB_Multifree(char *filename, int lineno, int n, ...)
{
  int i;
  va_list va;
  
  va_start(va, n);
  for (i=0; i<n; i++){
    LB_Free(va_arg(va, void **), filename, lineno);
  }
  va_end(va);
}

#else

void LB_Multifree(va_alist)
va_decl
{
  int i, n, lineno;
  char *filename;
  va_list va;
   
  va_start(va);
  filename = va_arg(va, char *);
  lineno = va_arg(va, int);
  n = va_arg(va, int);
  for (i=0; i<n; i++){
    LB_Free(va_arg(va, void **), filename, lineno);
  }
  va_end(va);
}

#endif

#endif /* if 0 */

/* Print out status of malloc/free calls.  Flag any memory leaks. */

void      LB_Memory_Stats()
{
    struct malloc_debug_data *dbptr;	/* loops through debug list */
    int       proc;		/* processor ID */


    if (DEBUG_MEMORY == 1) {
        MPI_Comm_rank(MPI_COMM_WORLD, &proc);
	fprintf(stderr, "Proc %d: Calls to malloc = %d,  Calls to free = %d\n",
		proc, nmalloc, nfree);
    }
    if (DEBUG_MEMORY > 1) {
        MPI_Comm_rank(MPI_COMM_WORLD, &proc);
	fprintf(stderr, "Proc %d: Calls to malloc = %d,  Calls to free = %d, maximum bytes = %d\n",
		proc, nmalloc, nfree, bytes_max);
	if (top != NULL) {
	    fprintf(stderr, "Proc %d: Remaining allocations:\n", proc);
	    for (dbptr = top; dbptr != NULL; dbptr = dbptr->next) {
		fprintf(stderr, " order=%d, size=%d, location=0x%lx, file=%s, line=%d\n", 
                  dbptr->order, dbptr->size, (long) dbptr->ptr,
                  dbptr->file, dbptr->line);
	    }
	}
    }
} /* LB_Memory_Stats */


/* Return number associated with the most recent malloc call. */
int       LB_Malloc_Num()
{
  return (nmalloc);
} /* LB_Memory_Num */


/******************************************************************************
 * Special allocation for routines that allocate an array and return pointer.
 *
 * LB_Special_Malloc allows the allocation to be done from either C or Fortran.
 *
 * LB_Special_Free frees memory allocated by LB_Special_Malloc
 *
 * LB_Register_Fort_Malloc is called by the wrappers for the Fortran
 * interface to provide pointers to the Fortran allocation/free routines.
 *
 * int LB_Special_Malloc(struct LB_Struct *lb, void **array, int size,
 *                       LB_SPECIAL_MALLOC_TYPE type)
 *
 *   lb    -- the load balancing structure in use
 *   array -- int**, struct LB_GID**, or struct LB_LID**; returned as a
 *            pointer to the allocated space
 *   size  -- number of elements to be allocated in the array
 *   type  -- the type of array; LB_SPECIAL_MALLOC_INT, LB_SPECIAL_MALLOC_GID,
 *            or LB_SPECIAL_MALLOC_LID
 *
 * The return value is 1 if the allocation succeeded, 0 if it failed.
 *
 * int LB_Special_Free(struct LB_Struct *lb, void **array,
                       LB_SPECIAL_MALLOC_TYPE type)
 *
 *****************************************************************************/

void LB_Register_Fort_Malloc(LB_FORT_MALLOC_INT_FN *fort_malloc_int,
                             LB_FORT_MALLOC_GID_FN *fort_malloc_GID,
                             LB_FORT_MALLOC_LID_FN *fort_malloc_LID,
                             LB_FORT_FREE_INT_FN *fort_free_int,
                             LB_FORT_FREE_GID_FN *fort_free_GID,
                             LB_FORT_FREE_LID_FN *fort_free_LID)
{
   LB_Fort_Malloc_int = fort_malloc_int;
   LB_Fort_Malloc_GID = fort_malloc_GID;
   LB_Fort_Malloc_LID = fort_malloc_LID;
   LB_Fort_Free_int = fort_free_int;
   LB_Fort_Free_GID = fort_free_GID;
   LB_Fort_Free_LID = fort_free_LID;
}

int LB_Special_Malloc(struct LB_Struct *lb, void **array, int size,
                      LB_SPECIAL_MALLOC_TYPE type)
{
   int *ret_addr, success;

   success = 1;
   if (lb->Fortran) {

/* allocation from Fortran */

      switch(type) {
      case LB_SPECIAL_MALLOC_INT:
#ifdef PGI /* special case for PGI Fortran compiler */
         LB_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr,array[2]);
#else
#ifdef FUJITSU /* special case for Fujitsu and Lahey Fortran compilers */
         LB_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr,array[2],0,0);
#else
         LB_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr);
#endif
#endif
         if (ret_addr==0) success=0;
         break;
      case LB_SPECIAL_MALLOC_GID:
         if (LB_GID_IS_INT) {
#ifdef PGI
            LB_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr,array[2]);
#else
#ifdef FUJITSU
            LB_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr,array[2],0,0);
#else
            LB_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr);
#endif
#endif
         }else{
#ifdef PGI
            LB_Fort_Malloc_GID((LB_GID *)(array[1]),&size,&ret_addr,array[2]);
#else
#ifdef FUJITSU
            LB_Fort_Malloc_GID((LB_GID *)(array[1]),&size,&ret_addr,array[2],0,0);
#else
            LB_Fort_Malloc_GID((LB_GID *)(array[1]),&size,&ret_addr);
#endif
#endif
         }
         if (ret_addr==0) success=0;
         break;
      case LB_SPECIAL_MALLOC_LID:
         if (LB_LID_IS_INT) {
#ifdef PGI
            LB_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr,array[2]);
#else
#ifdef FUJITSU
            LB_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr,array[2],0,0);
#else
            LB_Fort_Malloc_int((int *)(array[1]),&size,&ret_addr);
#endif
#endif
         }else{
#ifdef PGI
            LB_Fort_Malloc_LID((LB_LID *)(array[1]),&size,&ret_addr,array[2]);
#else
#ifdef FUJITSU
            LB_Fort_Malloc_LID((LB_LID *)(array[1]),&size,&ret_addr,array[2],0,0);
#else
            LB_Fort_Malloc_LID((LB_LID *)(array[1]),&size,&ret_addr);
#endif
#endif
         }
         if (ret_addr==0) success=0;
         break;
      default:
         fprintf(stderr, "Error: illegal value passed for type in LB_Special_Malloc\n");
         success = 0;
      }
      if (success) {
         array[0] = ret_addr;
      }else{
         array[0] = NULL;
      }

   }else{

/* allocation from C */

      switch(type) {
      case LB_SPECIAL_MALLOC_INT:
         *array = (int *) LB_Malloc(size*sizeof(int),__FILE__,__LINE__);
         break;
      case LB_SPECIAL_MALLOC_GID:
         *array = (LB_GID *) LB_Malloc(size*sizeof(LB_GID),__FILE__,__LINE__);
         break;
      case LB_SPECIAL_MALLOC_LID:
         *array = (LB_LID *) LB_Malloc(size*sizeof(LB_LID),__FILE__,__LINE__);
         break;
      default:
         fprintf(stderr, "Error: illegal value passed for type in LB_Special_Malloc\n");
         *array = NULL;
      }
      if (*array==NULL) success=0;
   }
   return success;
}

int LB_Special_Free(struct LB_Struct *lb, void **array,
                    LB_SPECIAL_MALLOC_TYPE type)
{
   int success;

   success = 1;
   if (lb->Fortran) {

/* deallocation from Fortran */

      switch(type) {
      case LB_SPECIAL_MALLOC_INT:
#ifdef PGI /* special case for PGI Fortran compiler */
         LB_Fort_Free_int((int *)(array[1]),array[2]);
#else
#ifdef FUJITSU /* special case for Fujitsu and Lahey Fortran compilers */
         LB_Fort_Free_int((int *)(array[1]),array[2]);
#else
         LB_Fort_Free_int((int *)(array[1]));
#endif
#endif
         break;
      case LB_SPECIAL_MALLOC_GID:
         if (LB_GID_IS_INT) {
#ifdef PGI
            LB_Fort_Free_int((int *)(array[1]),array[2]);
#else
#ifdef FUJITSU
            LB_Fort_Free_int((int *)(array[1]),array[2]);
#else
            LB_Fort_Free_int((int *)(array[1]));
#endif
#endif
         }else{
#ifdef PGI
            LB_Fort_Free_GID((LB_GID *)(array[1]),array[2]);
#else
#ifdef FUJITSU
            LB_Fort_Free_GID((LB_GID *)(array[1]),array[2]);
#else
            LB_Fort_Free_GID((LB_GID *)(array[1]));
#endif
#endif
         }
         break;
      case LB_SPECIAL_MALLOC_LID:
         if (LB_LID_IS_INT) {
#ifdef PGI
            LB_Fort_Free_int((int *)(array[1]),array[2]);
#else
#ifdef FUJITSU
            LB_Fort_Free_int((int *)(array[1]),array[2]);
#else
            LB_Fort_Free_int((int *)(array[1]));
#endif
#endif
         }else{
#ifdef PGI
            LB_Fort_Free_LID((LB_LID *)(array[1]),array[2]);
#else
#ifdef FUJITSU
            LB_Fort_Free_LID((LB_LID *)(array[1]),array[2]);
#else
            LB_Fort_Free_LID((LB_LID *)(array[1]));
#endif
#endif
         }
         break;
      default:
         fprintf(stderr, "Error: illegal value passed for type in LB_Special_Free\n");
         success = 0;
      }

   }else{

/* deallocation from C */

      LB_FREE(array);
   }
   return success;
}

/*****************************************************************************/
/*                      END of all_allo.c                                     */
/*****************************************************************************/
