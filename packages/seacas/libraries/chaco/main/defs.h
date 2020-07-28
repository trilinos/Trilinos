#ifndef CHACO_MAIN_DEFS_H
#define CHACO_MAIN_DEFS_H

/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#include <math.h>
#include <stdlib.h>

#ifndef max
#define max(A, B) ((A) > (B) ? (A) : (B))
#endif
#ifndef min
#define min(A, B) ((A) < (B) ? (A) : (B))
#endif
#ifndef sign
#define sign(A) ((A) < 0 ? -1 : 1)
#endif
#ifndef absval
#define absval(A) ((A) < 0 ? -(A) : (A))
#endif
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/* Define constants that are needed in various places */
#if defined(M_PI)
#define PI M_PI
#define TWOPI (2.0 * M_PI)
#define HALFPI (0.5 * M_PI)
#else
#define PI 3.141592653589793
#define TWOPI 6.283185307179586
#define HALFPI 1.570796326794896
#endif

#endif
