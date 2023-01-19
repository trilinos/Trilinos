// Copyright(C) 1999-2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#ifdef __cplusplus
extern "C" {
#endif

char *ap_getline_int(char *); /* read a line of input */
void  ap_gl_setwidth(int);    /* specify width of screen */
void  ap_gl_histadd(char *);  /* adds entries to hist */

#ifdef __cplusplus
} /* close brackets on extern "C" declaration */
#endif
