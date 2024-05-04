/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <stddef.h> // for size_t
#include <stdio.h>  // for NULL, fprintf, printf, FILE, etc
#include <stdlib.h> // for malloc, free, realloc

static int    nmalloc    = 0; /* number of calls to malloc */
static int    nfree      = 0; /* number of calls to free */
static size_t bytes_used = 0; /* current dynamic memory usage */
static size_t bytes_max  = 0; /* largest dynamic memory usage */

static struct smalloc_debug_data
{
  int                        order; /* which smalloc call is it? */
  size_t                     size;  /* size of malloc invocation */
  double                    *ptr;   /* memory location returned */
  struct smalloc_debug_data *next;  /* pointer to next element */
} *top = NULL;

void message(char *msg, size_t n, FILE *ofile)
{
  fprintf(stderr, msg, n);
  if (ofile != NULL) {
    fprintf(ofile, msg, n);
  }
}

void sfree(void *ptr);

/* Safe version of malloc.  Does not initialize memory .*/
void *smalloc(size_t n)
{
  extern FILE *Output_File;       /* output file or null */
  extern int   DEBUG_MEMORY;      /* use debug memory allocator? */
  void        *ptr;               /* return value */
  struct smalloc_debug_data *new; /* data structure for malloc data */
  void bail(char *msg, int status);

  nmalloc++;
  if (n == 0) {
    message("ERROR: Non-positive argument to smalloc (%zu).\n", n, Output_File);
    bail(NULL, 1);
  }

  ptr = malloc(n);

  if (ptr == NULL) {
    message("Program out of space while attempting to allocate %zu.\n", n, Output_File);
    bail(NULL, 1);
  }

  if (DEBUG_MEMORY > 1) {
    new = (struct smalloc_debug_data *)malloc(sizeof(struct smalloc_debug_data));

    if (new == NULL) {
      message("WARNING: No space for malloc_debug %zu.\n", n, Output_File);
      return (ptr);
    }

    new->order = nmalloc;
    new->size  = n;
    new->ptr   = ptr;
    new->next  = top;
    top        = new;
    bytes_used += n;
    if (bytes_used > bytes_max) {
      bytes_max = bytes_used;
    }
  }

  if (DEBUG_MEMORY > 2) {
    printf(" order=%d, size=%zu, location=%p\n", nmalloc, n, ptr);
  }

  return (ptr);
}

/* Safe version of malloc.  Does not initialize memory .*/
/* Returns instead of dying if it fails. */
void *smalloc_ret(size_t n)
{
  extern FILE *Output_File;       /* output file or null */
  extern int   DEBUG_MEMORY;      /* use debug memory allocator? */
  void        *ptr;               /* return value */
  struct smalloc_debug_data *new; /* data structure for malloc data */
  void bail(char *msg, int status);

  ptr = NULL;
  if (n == 0) {
    message("ERROR: Non-positive argument to smalloc_ret (%zu).\n", n, Output_File);
    bail(NULL, 1);
  }

  else {
    nmalloc++;
    ptr = malloc(n);

    if (ptr == NULL) {
      nmalloc--;
      message("\nERROR: Unable to allocate %zu bytes of memory.\n", n, Output_File);
    }
    else {

      if (DEBUG_MEMORY > 1) {
        new = (struct smalloc_debug_data *)malloc(sizeof(struct smalloc_debug_data));

        if (new == NULL) {
          message("WARNING: No space for malloc_debug %zu.\n", n, Output_File);
          return (ptr);
        }

        new->order = nmalloc;
        new->size  = n;
        new->ptr   = ptr;
        new->next  = top;
        top        = new;
        bytes_used += n;
        if (bytes_used > bytes_max) {
          bytes_max = bytes_used;
        }
      }

      if (DEBUG_MEMORY > 2) {
        printf(" order=%d, size=%zu, location=%p\n", nmalloc, n, ptr);
      }
    }
  }

  return (ptr);
}

/* Safe version of realloc */
void *srealloc(void *ptr, size_t n)
{
  extern FILE               *Output_File;  /* output file or null */
  void                      *p;            /* returned pointer */
  extern int                 DEBUG_MEMORY; /* use debug memory allocator? */
  struct smalloc_debug_data *dbptr;        /* loops through debug list */
  void                       bail(char *msg, int status);

  if (ptr == NULL) {
    if (n == 0) {
      return (NULL);
    }

    p = smalloc(n);
  }
  else {
    if (n == 0) {
      sfree(ptr);
      return (NULL);
    }

    p = realloc(ptr, n);
    if (DEBUG_MEMORY > 1) {
      for (dbptr = top; dbptr != NULL && dbptr->ptr != ptr; dbptr = dbptr->next) {
        ;
      }
      if (dbptr == NULL) {
        fprintf(stderr, "Memory error: In srealloc, address not found in debug list (%p)\n", ptr);
      }
      else {
        dbptr->size = n;
        dbptr->ptr  = p;
        bytes_used += n;
        if (bytes_used > bytes_max) {
          bytes_max = bytes_used;
        }
      }
    }
  }

  if (p == NULL) {
    message("Program out of space while attempting to reallocate %zu.\n", n, Output_File);
    bail(NULL, 1);
  }
  return (p);
}

/* Safe version of realloc */
/* Returns instead of dying if it fails. */
void *srealloc_ret(void *ptr, size_t n)
{
  extern FILE               *Output_File;  /* output file or null */
  void                      *p;            /* returned pointer */
  extern int                 DEBUG_MEMORY; /* use debug memory allocator? */
  struct smalloc_debug_data *dbptr;        /* loops through debug list */

  if (ptr == NULL) {
    if (n == 0) {
      return (NULL);
    }

    p = smalloc(n);
  }
  else {
    if (n == 0) {
      sfree(ptr);
      return (NULL);
    }

    p = realloc(ptr, n);
    if (DEBUG_MEMORY > 1) {
      for (dbptr = top; dbptr != NULL && dbptr->ptr != ptr; dbptr = dbptr->next) {
        ;
      }
      if (dbptr == NULL) {
        fprintf(stderr, "Memory error: In srealloc_ret, address not found in debug list (%p)\n",
                ptr);
      }
      else {
        dbptr->size = n;
        dbptr->ptr  = p;
        bytes_used += n;
        if (bytes_used > bytes_max) {
          bytes_max = bytes_used;
        }
      }
    }
  }

  if (p == NULL && DEBUG_MEMORY > 0) {
    message("WARNING: Program out of space while attempting to reallocate %zu.\n", n, Output_File);
  }
  return (p);
}

/* Safe version of free. */
void sfree(void *ptr)
{
  extern FILE                *Output_File;  /* output file or null */
  extern int                  DEBUG_MEMORY; /* use debug memory allocator? */
  struct smalloc_debug_data  *dbptr;        /* loops through debug list */
  struct smalloc_debug_data **prev;         /* holds previous pointer */

  if (DEBUG_MEMORY > 1) {
    if (ptr != NULL) { /* search through debug list for it */
      prev = &top;
      for (dbptr = top; dbptr != NULL && dbptr->ptr != ptr; dbptr = dbptr->next) {
        prev = &(dbptr->next);
      }
      if (dbptr == NULL) {
        fprintf(stderr, "Memory error: In sfree, address not found in debug list (%p)\n", ptr);
        if (Output_File != NULL) {
          fprintf(Output_File, "Memory error: In sfree, address not found in debug list (%p)\n",
                  ptr);
        }
      }
      else {
        *prev = dbptr->next;
        bytes_used -= dbptr->size;
        free(dbptr);
      }
    }
  }

  if (ptr != NULL) {
    nfree++;
    free(ptr);
  }
}

void smalloc_stats(void)
{
  extern int                 DEBUG_MEMORY; /* use debug memory allocator? */
  struct smalloc_debug_data *dbptr;        /* loops through debug list */

  if (DEBUG_MEMORY == 1) {
    printf("Calls to smalloc = %d,  Calls to sfree = %d\n", nmalloc, nfree);
  }
  if (DEBUG_MEMORY > 1) {
    printf("Calls to smalloc = %d,  Calls to sfree = %d, maximum bytes = %zu\n", nmalloc, nfree,
           bytes_max);
    if (top != NULL) {
      printf("Remaining allocations:\n");
      for (dbptr = top; dbptr != NULL; dbptr = dbptr->next) {
        printf(" order=%d, size=%zu, location=%p\n", dbptr->order, dbptr->size, (void *)dbptr->ptr);
      }
    }
  }
}

int smalloc_num(void) { return (nmalloc); }
