/*
 * Copyright(C) 1999-2021 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#ifndef GETLINE_H
#define GETLINE_H

#define GL_BUF_SIZE 1024

char *getline_int(char *); /* read a line of input */
void  gl_setwidth(int);    /* specify width of screen */
void  gl_histadd(char *);  /* adds entries to hist */

#endif /* GETLINE_H */
