/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef CHACO_UTIL_CHACO_RANDOM_H
#define CHACO_UTIL_CHACO_RANDOM_H

#ifdef __cplusplus
extern "C" {
#endif

long   init_rand_port(long seed);
long   get_init_rand_port(void);
long   genr_rand_port(long init_rand);
long   rand_port(void);
double rand_rect_port(void);

#ifdef __cplusplus
} /* close brackets on extern "C" declaration */
#endif

#endif
