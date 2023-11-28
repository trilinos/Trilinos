/*
 * Copyright(C) 1999-2020, 2022, 2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

unsigned long init_rand_port(unsigned long seed);
unsigned long get_init_rand_port(void);
unsigned long genr_rand_port(unsigned long init_rand);
unsigned long rand_port(void);
double        rand_rect_port(void);

#ifdef __cplusplus
} /* close brackets on extern "C" declaration */
#endif
