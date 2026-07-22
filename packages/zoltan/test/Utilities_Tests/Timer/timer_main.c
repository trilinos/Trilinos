// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "zoltan_timer.h"
#include <math.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Test program to exercise Timer capabilities. */

void first_test(struct Zoltan_Timer*);
void second_test(struct Zoltan_Timer*);
void third_test(struct Zoltan_Timer*);

/****************************************************************************/

int main(int argc, char *argv[])
{
struct Zoltan_Timer *zt1, *zt2, *zt3, *zt4;
int i, me;
const int MAINLOOP=20;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  zt1 = Zoltan_Timer_Create(ZOLTAN_TIME_WALL);
  zt2 = Zoltan_Timer_Create(ZOLTAN_TIME_USER);
  zt3 = Zoltan_Timer_Create(ZOLTAN_TIME_WALL);

  for (i = 0; i < MAINLOOP; i++) {
    if (me == 0) printf("\n\n\t****Beginning first test****\n");
    first_test(zt1);

    if (me == 0) printf("\n\n\t****Beginning second test****\n");
    second_test(zt2);
  }

  if (me == 0) printf("\n\nFINAL RESULTS -- FIRST TEST:\n");
  Zoltan_Timer_PrintAll(zt1, 0, MPI_COMM_WORLD, stdout);
  if (me == 0) printf("\n\nFINAL RESULTS -- SECOND TEST:\n");
  Zoltan_Timer_PrintAll(zt2, 0, MPI_COMM_WORLD, stdout);

  /* Copy tests */
  Zoltan_Timer_Copy_To(&zt3, zt1);
  zt4 = Zoltan_Timer_Copy(zt2);
  for (i = 0; i < MAINLOOP; i++) {
    if (me == 0) printf("\n\n\t****Beginning first copy test****\n");
    first_test(zt3);

    if (me == 0) printf("\n\n\t****Beginning second copy test****\n");
    second_test(zt4);
  }
  if (me == 0) printf("\n\nFINAL RESULTS -- FIRST COPY TEST:\n");
  Zoltan_Timer_PrintAll(zt3, 0, MPI_COMM_WORLD, stdout);
  if (me == 0) printf("\n\nFINAL RESULTS -- SECOND COPY TEST:\n");
  Zoltan_Timer_PrintAll(zt4, 0, MPI_COMM_WORLD, stdout);

  /* Test printing while timer is still running.  */
  if (me == 0) printf("\n\n\t****Intermediate print test****\n");
  third_test(zt1);
  if (me == 0) printf("\n\nFINAL RESULTS -- INTERMEDIATE PRINT TEST:\n");
  Zoltan_Timer_PrintAll(zt1, 0, MPI_COMM_WORLD, stdout);

  if (me == 0) printf("\n\nTHE END\n");

  Zoltan_Timer_Destroy(&zt1);
  Zoltan_Timer_Destroy(&zt2);

  MPI_Finalize();

  return 0;
}

/****************************************************************************/

void first_test(struct Zoltan_Timer *zt)
{
/* First test of Timer:  This test accrues times through
 * separate calls to first_test.  
 * The time for timer two should be roughly twice that of timer one.
 * The time for timer three should be roughly four times that of timer one.
 */
int i, j, me; 
static int firsttime=1;
const int LOOP1=1000,
          LOOP2=2000,
          LOOP3=4000;
const int MAINLOOP=100;
const int USE_BARRIER=1;
static int t1=-1, t2=-1, t3=-1;

  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  for (i = 0; i < MAINLOOP; i++) {

    if (firsttime)
      t1 = Zoltan_Timer_Init(zt, USE_BARRIER, "Loop 1");

    ZOLTAN_TIMER_START(zt, t1, MPI_COMM_WORLD);
    for (j = 0; j < LOOP1; j++) {
      double a;
      a = sqrt((double) (j * LOOP1));
    }
    ZOLTAN_TIMER_STOP(zt, t1, MPI_COMM_WORLD);

    if (firsttime)
      t2 = Zoltan_Timer_Init(zt, USE_BARRIER, "Loop 2");

    ZOLTAN_TIMER_START(zt, t2, MPI_COMM_WORLD);
    for (j = 0; j < LOOP2; j++) {
      double a;
      a = sqrt((double) (j * LOOP2));
    }
    ZOLTAN_TIMER_STOP(zt, t2, MPI_COMM_WORLD);

    if (firsttime)
      t3 = Zoltan_Timer_Init(zt, USE_BARRIER, "Loop 3");

    ZOLTAN_TIMER_START(zt, t3, MPI_COMM_WORLD);
    for (j = 0; j < LOOP3; j++) {
      double a;
      a = sqrt((double) (j * LOOP3));
    }
    ZOLTAN_TIMER_STOP(zt, t3, MPI_COMM_WORLD);

    firsttime=0;
  }

  Zoltan_Timer_PrintAll(zt, 0, MPI_COMM_WORLD, stdout);

}

/****************************************************************************/

void second_test(struct Zoltan_Timer *zt)
{
/* Second test of Timer:  This test does not accrue times through
 * separate function calls.  It exercises the REALLOC when more timers
 * than INITLENGTH are requested. 
 * Computation is similar to first test.  Timer_flag used is different.
 */
int i, j, me; 
const int LOOP1=1000,
          LOOP2=2000,
          LOOP3=4000;
const int MAINLOOP=100;
const int USE_BARRIER=1;
int t1=-1, t2=-1, t3=-1;
char str[200];
static int cnt = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  sprintf(str, "STLoop 1 %d", cnt);
  t1 = Zoltan_Timer_Init(zt, USE_BARRIER, str);
  sprintf(str, "STLoop 2 %d", cnt);
  t2 = Zoltan_Timer_Init(zt, USE_BARRIER, str);
  sprintf(str, "STLoop 3 %d", cnt);
  t3 = Zoltan_Timer_Init(zt, USE_BARRIER, str);
  cnt++;

  for (i = 0; i < MAINLOOP; i++) {
    ZOLTAN_TIMER_START(zt, t1, MPI_COMM_WORLD);
    for (j = 0; j < LOOP1; j++) {
      double a;
      a = sqrt((double) (j * LOOP1));
    }
    ZOLTAN_TIMER_STOP(zt, t1, MPI_COMM_WORLD);

    ZOLTAN_TIMER_START(zt, t2, MPI_COMM_WORLD);
    for (j = 0; j < LOOP2; j++) {
      double a;
      a = sqrt((double) (j * LOOP2));
    }
    ZOLTAN_TIMER_STOP(zt, t2, MPI_COMM_WORLD);

    ZOLTAN_TIMER_START(zt, t3, MPI_COMM_WORLD);
    for (j = 0; j < LOOP3; j++) {
      double a;
      a = sqrt((double) (j * LOOP3));
    }
    ZOLTAN_TIMER_STOP(zt, t3, MPI_COMM_WORLD);
  }

  Zoltan_Timer_PrintAll(zt, 0, MPI_COMM_WORLD, stdout);
}

void third_test(struct Zoltan_Timer *zt)
{
/* Third test of Timer:  This test accrues times through
 * separate calls to third_test.  
 * The time for timer two should be roughly twice that of timer one.
 * The time for timer three should be roughly four times that of timer one.
 * Intermediate print statements are included (i.e., prints while the timer
 * is still running).
 */
int i, j, me; 
const int LOOP1=1000,
          LOOP2=2000,
          LOOP3=4000;
const int MAINLOOP=100;
const int t1=0, t2=1, t3=2;

  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  for (i = 0; i < MAINLOOP; i++) {

    ZOLTAN_TIMER_START(zt, t1, MPI_COMM_WORLD);
    for (j = 0; j < LOOP1; j++) {
      double a;
      a = sqrt((double) (j * LOOP1));
      if (!(j%1000)) {
        if (me == 0) printf("LOOP1 %d: \n", j);
        Zoltan_Timer_Print(zt, t1, 0, MPI_COMM_WORLD, stdout);
      }
    }
    ZOLTAN_TIMER_STOP(zt, t1, MPI_COMM_WORLD);

    ZOLTAN_TIMER_START(zt, t2, MPI_COMM_WORLD);
    for (j = 0; j < LOOP2; j++) {
      double a;
      a = sqrt((double) (j * LOOP2));
      if (!(j%1000)) {
        if (me == 0) printf("LOOP2 %d: \n", j);
        Zoltan_Timer_Print(zt, t2, 0, MPI_COMM_WORLD, stdout);
      }
    }
    ZOLTAN_TIMER_STOP(zt, t2, MPI_COMM_WORLD);

    ZOLTAN_TIMER_START(zt, t3, MPI_COMM_WORLD);
    for (j = 0; j < LOOP3; j++) {
      double a;
      a = sqrt((double) (j * LOOP3));
      if (!(j%1000)) {
        if (me == 0) printf("LOOP3 %d: \n", j);
        Zoltan_Timer_Print(zt, t3, 0, MPI_COMM_WORLD, stdout);
      }
    }
    ZOLTAN_TIMER_STOP(zt, t3, MPI_COMM_WORLD);
  }

  Zoltan_Timer_PrintAll(zt, 0, MPI_COMM_WORLD, stdout);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
