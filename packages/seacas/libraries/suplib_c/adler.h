/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef ADLER_H
#define ADLER_H

#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
size_t adler(size_t adler, const void *vbuf, size_t len);
#ifdef __cplusplus
} /* close brackets on extern "C" declaration */
#endif
#endif
