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

#ifndef __PHG_UTIL_H
#define __PHG_UTIL_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdarg.h>
#include "zz_const.h"
#include "phg.h"



/* Hypergraph utilities */
extern void Zoltan_PHG_HGraph_Init (PHGraph*);
extern int Zoltan_PHG_HGraph_Free  (PHGraph*);
extern int Zoltan_PHG_Create_Mirror(ZZ*, PHGraph*);

extern void Zoltan_PHG_Graph_Init  (PGraph*);
extern int Zoltan_PHG_Graph_Free   (PGraph*);
extern int Zoltan_PHG_Info         (ZZ*, PHGraph*);
extern int Zoltan_PHG_Check        (ZZ*, PHGraph*);
extern int Zoltan_PHG_Graph_to_HGraph (ZZ*, PGraph*, PHGraph*);
extern void Zoltan_PHG_Print (ZZ*, PHGraph*, FILE*);

extern unsigned int Zoltan_PHG_Rand (void);
extern void         Zoltan_PHG_Srand (unsigned int);
extern void         Zoltan_PHG_Rand_Perm_Int (int*, int);

    /* UVC: some utility functions not particularly related to hypergraph */
extern char *uMe(PHGComm *);
extern void uprintf(PHGComm *, char *,...);
extern void errexit(char *,...);
    
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __HG_UTIL_H */
