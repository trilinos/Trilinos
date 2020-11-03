/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef SUPLIB_C_COPY_STRING
#define SUPLIB_C_COPY_STRING
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
char *copy_string(char *dest, char const *source, size_t elements);
#ifdef __cplusplus
} /* close brackets on extern "C" declaration */
#endif
#endif
