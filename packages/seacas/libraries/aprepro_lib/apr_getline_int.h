// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef APREPRO_GETLINE_H
#define APREPRO_GETLINE_H

#ifdef __cplusplus
extern "C" {
#endif

char *ap_getline_int(char *); /* read a line of input */
void  ap_gl_setwidth(int);    /* specify width of screen */
void  ap_gl_histadd(char *);  /* adds entries to hist */

#ifdef __cplusplus
} /* close brackets on extern "C" declaration */
#endif

#endif /* APREPRO_GETLINE_H */
