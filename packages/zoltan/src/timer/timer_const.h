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

#ifndef __TIMER_CONST_H
#define __TIMER_CONST_H

/* Constants used in Zoltan timer routines */
#define ZOLTAN_TIME_WALL 1
#define ZOLTAN_TIME_CPU  2
#define ZOLTAN_TIME_USER 3
#define ZOLTAN_TIME_USERSYS 4

/* Function prototypes */
extern double Zoltan_Time(int);
extern double Zoltan_Time_Resolution(int);
extern int Zoltan_Set_Timer_Param(char *, char *, int *);

#endif
