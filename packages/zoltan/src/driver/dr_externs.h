// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _DR_EXTERNS_H
#define _DR_EXTERNS_H

#include "dr_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Global variables for driver */
extern int Debug_Driver;
extern int Debug_Chaco_Input;
extern int Number_Iterations;
extern int Driver_Action;
extern int Chaco_In_Assign_Inv;
extern struct Test_Flags Test;
extern struct Output_Flags Output;
extern double Total_Partition_Time;

#define DEBUG_TRACE_START(proc,yo) \
  if (((proc) == 0 && Debug_Driver > 1) || (Debug_Driver > 2))  \
    printf("%d DRIVER ENTERING %s\n", (proc), yo);
#define DEBUG_TRACE_END(proc,yo) \
  if (((proc) == 0 && Debug_Driver > 1) || (Debug_Driver > 2))  \
    printf("%d DRIVER LEAVING %s\n", (proc), yo);
#define DEBUG_TRACE_DETAIL(proc,yo,str) \
  if (Debug_Driver > 2) \
    printf("%d DRIVER %s: %s\n", proc,yo, str);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
#endif /* _DR_EXTERNS_H */
