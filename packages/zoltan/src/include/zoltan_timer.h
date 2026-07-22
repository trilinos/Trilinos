// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __ZOLTANTIMER_H
#define __ZOLTANTIMER_H

#include <stdio.h>
#include <mpi.h>


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Constants used in Zoltan timer routines */
#define ZOLTAN_TIME_WALL 1
#define ZOLTAN_TIME_CPU  2
#define ZOLTAN_TIME_USER 3

/* Macros to add line/file info */

#define ZOLTAN_TIMER_START(a, b, c) \
  Zoltan_Timer_Start(a, b, c, __FILE__, __LINE__)
#define ZOLTAN_TIMER_STOP(a, b, c) \
  Zoltan_Timer_Stop(a, b, c, __FILE__, __LINE__)

/* Function prototypes */

struct Zoltan_Timer;

struct Zoltan_Timer *Zoltan_Timer_Create(int);
int Zoltan_Timer_Init(struct Zoltan_Timer *, int, const char *);
struct Zoltan_Timer *Zoltan_Timer_Copy(struct Zoltan_Timer *zt);
int Zoltan_Timer_Copy_To(struct Zoltan_Timer **to, struct Zoltan_Timer *from);
int Zoltan_Timer_Reset(struct Zoltan_Timer *, int, int, const char*);
int Zoltan_Timer_ChangeFlag(struct Zoltan_Timer *, int);
int Zoltan_Timer_Start(struct Zoltan_Timer *, int, MPI_Comm, char *, int);
int Zoltan_Timer_Stop(struct Zoltan_Timer *, int, MPI_Comm, char *, int);
int Zoltan_Timer_Print(struct Zoltan_Timer *, int, int, MPI_Comm, FILE *);
int Zoltan_Timer_PrintAll(struct Zoltan_Timer *, int, MPI_Comm, FILE *);
void Zoltan_Timer_Destroy(struct Zoltan_Timer **);

extern double Zoltan_Time(int);
extern double Zoltan_Time_Resolution(int);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
