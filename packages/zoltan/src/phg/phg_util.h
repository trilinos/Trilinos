// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __PHG_UTIL_H
#define __PHG_UTIL_H

#include <stdarg.h>
#include "phg_comm.h"
#include "phg_hypergraph.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Definitions to allow simplicity in PHG code 
 * while protecting application namespace. 
 */
#define uMe Zoltan_PHG_uMe
#define uprintf Zoltan_PHG_uprintf
#define errexit Zoltan_PHG_errexit

/* UVC: some utility functions not particularly related to hypergraph */
extern char *Zoltan_PHG_uMe(PHGComm *);
extern void Zoltan_PHG_uprintf(PHGComm *, char *,...);
extern void Zoltan_PHG_errexit(char *,...);

extern int Zoltan_PHG_isPrime(int);

extern void Zoltan_PHG_Find_Root(int, int, MPI_Comm, int *, int *);

extern int Zoltan_PHG_LoadBalStat(ZZ *zz, HGraph *);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __HG_UTIL_H */
