/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifndef CH_INIT_DIST_CONST_H
#define CH_INIT_DIST_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "dr_const.h"
#include "dr_input_const.h"

/* define the Chaco initial distribution types */
#define INITIAL_FILE   0
#define INITIAL_LINEAR 1
#define INITIAL_CYCLIC 2
#define INITIAL_OWNER  3

extern void ch_dist_init(int, int, PARIO_INFO_PTR, short **, int, MPI_Comm);
extern int ch_dist_num_vtx(int, short *);
extern int ch_dist_max_num_vtx(short *);
extern void ch_dist_vtx_list(int*, int*, int, short *);
extern int ch_dist_proc(int, short *, int);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif  /* CH_INIT_DIST_CONST_H */
