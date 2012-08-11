/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */

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
