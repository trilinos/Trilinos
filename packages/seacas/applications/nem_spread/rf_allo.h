/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/* function declarations for dynamic array allocation */

extern void *array_alloc(const char *file, int lineno, int numdim, ...);

extern void safe_free(void **ptr);
