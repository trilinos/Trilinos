#ifndef CHACO_UTIL_SMALLOC_H
#define CHACO_UTIL_SMALLOC_H

/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#include <stddef.h>

/* Safe version of malloc.  Does not initialize memory .*/
extern void *smalloc(size_t n);

/* Safe version of malloc.  Does not initialize memory .*/
/* Returns instead of dying if it fails. */
extern void *smalloc_ret(size_t n);

/* Safe version of realloc */
extern void *srealloc(void *ptr, size_t n);

/* Safe version of realloc */
/* Returns instead of dying if it fails. */
extern void *srealloc_ret(void *ptr, size_t n);

/* Safe version of free. */
extern void sfree(void *ptr);

#endif
