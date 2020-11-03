/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#include "fmt/ostream.h"
#include <cstdarg> // for va_arg, va_end, va_list, etc
#include <cstddef> // for size_t
#include <cstdio>  // for nullptr, stderr
#include <cstdlib> // for exit, free, malloc

/*#include "rf_allo.h"*/

#ifdef DEBUG
extern int Proc;
#endif

static void *smalloc(size_t n, char *filename, int lineno);

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
 *      points = (POINT **) array_alloc (2, x, y, sizeof(POINT));
 *                               ^ ^ ^
 *                               | | |
 *         number of dimensions--+ | |
 *                                 | |
 *          first dimension max----+ |
 *                                   |
 *         second dimension max------+
 *
 *         (points may be now be used as if it were declared
 *          POINT points[x][y])
 *
 *          This particular version is limited to dimensions of 3 or less.
 *
 *      corner = points[2][3]; (refer to the structure as you would any array)
 *
 *      free (points); (frees the entire structure in one fell swoop)
 *
 *****************************************************************************/
/******************************************************************************
 *      The following section is a commented section containing
 *      an example main code:
 *******************************************************************************
 *double *array_alloc();
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
 *   temp = (int ***) array_alloc(3,il,jl,kl,sizeof(int));
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
 *         for (k=0; k<kl; k++) fmt::print(" {}\n", temp[i][j][k]);
 *      }
 *   }
 *   malloc_verify();
 *}
 ******************************************************************************/

/*****************************************************************************/
/*****************************************************************************/

void *array_alloc(const char *file, int lineno, int numdim, ...)
{
  struct dim
  {
    size_t index; /* Number of elements in the dimension        */
    size_t total; /* Total number of elements           */
    size_t size;  /* Size of a single element in bytes  */
    size_t off;   /* offset from beginning of array     */
  } dim[3];       /* Info about each dimension          */
  size_t  total;  /* Total size of the array            */
  void *  dfield; /* ptr to avoid lint complaints               */
  char *  field;  /* The multi-dimensional array                */
  char ** ptr;    /* Pointer offset                     */
  char *  data;   /* Data offset                                */
  va_list va;     /* Current pointer in the argument list       */

  va_start(va, numdim);

  if (numdim <= 0) {
    fmt::print(stderr, "{} ({}: {}) ERROR: number of dimensions, {}, is <=0\n", __func__, file,
               lineno, numdim);
    exit(1);
  }
  else if (numdim > 3) {
    fmt::print(stderr, "{} ({}: {}) ERROR: number of dimensions, {}, is > 3\n", __func__, file,
               lineno, numdim);
    exit(1);
  }

  dim[0].index = va_arg(va, size_t);

  if (dim[0].index == 0) {
#ifdef DEBUG
    fmt::print(stderr,
               "WARNING, {} ({}: {}) called with first "
               "dimension == 0; will return nullptr\n",
               __func__, file, lineno);
#endif
    va_end(va);
    return (nullptr);
  }

  dim[0].total = dim[0].index;
  dim[0].size  = sizeof(void *);
  dim[0].off   = 0;
  for (int i = 1; i < numdim; i++) {
    dim[i].index = va_arg(va, size_t);
    if (dim[i].index == 0) {
      fmt::print(stderr,
                 "WARNING: {} ({}: {}) called with dimension {} == 0, "
                 " will return nullptr\n",
                 __func__, file, lineno, i + 1);
      va_end(va);
      return (nullptr);
    }
    dim[i].total = dim[i - 1].total * dim[i].index;
    dim[i].size  = sizeof(void *);
    dim[i].off   = dim[i - 1].off + dim[i - 1].total * dim[i - 1].size;
  }
  dim[numdim - 1].size = va_arg(va, size_t);
  va_end(va);

  /* Round up the last offset value so data is properly aligned. */
  dim[numdim - 1].off = dim[numdim - 1].size *
                        ((dim[numdim - 1].off + dim[numdim - 1].size - 1) / dim[numdim - 1].size);

  total = dim[numdim - 1].off + dim[numdim - 1].total * dim[numdim - 1].size;

  dfield = smalloc(total, const_cast<char *>(file), lineno);
  field  = reinterpret_cast<char *>(dfield);

  for (int i = 0; i < numdim - 1; i++) {
    ptr  = reinterpret_cast<char **>(field + dim[i].off);
    data = (field + dim[i + 1].off);
    for (size_t j = 0; j < dim[i].total; j++) {
      ptr[j] = data + j * dim[i + 1].size * dim[i + 1].index;
    }
  }

  return (dfield);
}

/* Safe version of malloc.  Does not initialize memory .*/

/* Modified by Scott Hutchinson (1421) 20 January 1993 */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void *smalloc(size_t n, char *filename, int lineno)

{
  void *pntr = nullptr; /* return value */

  if (n != 0) {
    pntr = malloc(n);
  }

  if (pntr == nullptr && n != 0) {
    fmt::print(stderr,
               "{} (from {},{}) Out of space - number of bytes "
               "requested = {:n}\n",
               __func__, filename, lineno, n);
    exit(0);
  }

  return pntr;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void safe_free(void **ptr)
{
  /*
   *  This version of free calls the system's free function
   *  with maximum error checking. It also doesn't call free if ptr is
   *  the nullptr pointer.
   */

  if (*ptr != nullptr) {

    free(*ptr);

    /*
     *  Set the value of ptr to nullptr, so that further references
     *  to it will be flagged.
     */

    *ptr = nullptr;
  }
} /* safe_free */

/*****************************************************************************/
/*                       END of rf_allo.c                                    */
/*****************************************************************************/
