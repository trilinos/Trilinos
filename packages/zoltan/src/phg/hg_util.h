/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __HG_UTIL_H
#define __HG_UTIL_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#ifdef HGEXEC
#include "hypergraph.h"
#else
#include "zz_const.h"
#endif


#include "hg_hypergraph.h"



/* Hypergraph utilities */
extern void Zoltan_HG_HGraph_Init (HGraph*);
extern int Zoltan_HG_HGraph_Free  (HGraph*);
extern int Zoltan_HG_Create_Mirror(ZZ*, HGraph*);

extern void Zoltan_HG_Graph_Init  (Graph*);
extern int Zoltan_HG_Graph_Free   (Graph*);
extern int Zoltan_HG_Info         (ZZ*, HGraph*);
extern int Zoltan_HG_Check        (ZZ*, HGraph*);
extern int Zoltan_HG_Graph_to_HGraph(ZZ*, Graph*,  HGraph*);
extern void Zoltan_HG_Print(ZZ*, HGraph*, FILE*);

extern unsigned int Zoltan_HG_Rand (void);
extern void         Zoltan_HG_Srand (unsigned int);
extern void         Zoltan_HG_Rand_Perm_Int (int*, int);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __HG_UTIL_H */
