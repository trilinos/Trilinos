/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "Timer_dh.h"
#include "Mem_dh.h"

#undef __FUNC__
#define __FUNC__ "Timer_dhCreate"
void
Timer_dhCreate (Timer_dh * t)
{
  START_FUNC_DH
    struct _timer_dh *tmp =
    (struct _timer_dh *) MALLOC_DH (sizeof (struct _timer_dh));
  CHECK_V_ERROR;
  *t = tmp;

  tmp->isRunning = false;
  tmp->begin_wall = 0.0;
  tmp->end_wall = 0.0;
#ifdef WIN32
  tmp->sc_clk_tck = CLOCKS_PER_SEC;
#else
  tmp->sc_clk_tck = sysconf (128);
#endif

#if defined(EUCLID_TIMING)
  sprintf (msgBuf_dh, "using EUCLID_TIMING; _SC_CLK_TCK = %i",
	   (int) tmp->sc_clk_tck);
  SET_INFO (msgBuf_dh);
#elif defined(MPI_TIMING)
  SET_INFO ("using MPI timing")
#else
  SET_INFO ("using JUNK timing")
#endif
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "Timer_dhDestroy"
void
Timer_dhDestroy (Timer_dh t)
{
  START_FUNC_DH FREE_DH (t);
END_FUNC_DH}

/*-------------------------------------------------------------------------------
 * EUCLID_TIMING timing methods; these use times() to record 
 * both wall and cpu time.
 *-------------------------------------------------------------------------------*/

#ifdef EUCLID_TIMING

#undef __FUNC__
#define __FUNC__ "Timer_dhStart"
void
Timer_dhStart (Timer_dh t)
{
  START_FUNC_DH t->begin_wall = times (&(t->begin_cpu));
  t->isRunning = true;
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "Timer_dhStop"
void
Timer_dhStop (Timer_dh t)
{
  START_FUNC_DH t->end_wall = times (&(t->end_cpu));
  t->isRunning = false;
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "Timer_dhReadWall"
double
Timer_dhReadWall (Timer_dh t)
{
  START_FUNC_DH double retval = 0.0;
  long int sc_clk_tck = t->sc_clk_tck;
  if (t->isRunning)
    t->end_wall = times (&(t->end_cpu));
  retval = (double) (t->end_wall - t->begin_wall) / (double) sc_clk_tck;
END_FUNC_VAL (retval)}

#undef __FUNC__
#define __FUNC__ "Timer_dhReadCPU"
double
Timer_dhReadCPU (Timer_dh t)
{
  START_FUNC_DH double retval;
  long int sc_clk_tck = t->sc_clk_tck;
  if (t->isRunning)
    t->end_wall = times (&(t->end_cpu));
  retval = (double) (t->end_cpu.tms_utime - t->begin_cpu.tms_utime
		     + t->end_cpu.tms_stime - t->begin_cpu.tms_stime
		     + t->end_cpu.tms_cutime - t->begin_cpu.tms_cutime
		     + t->end_cpu.tms_cstime - t->begin_cpu.tms_cstime)
    / (double) sc_clk_tck;
END_FUNC_VAL (retval)}

#undef __FUNC__
#define __FUNC__ "Timer_dhReadUsage"
double
Timer_dhReadUsage (Timer_dh t)
{
  START_FUNC_DH double cpu = Timer_dhReadCPU (t);
  double wall = Timer_dhReadWall (t);
  double retval = 100.0 * cpu / wall;
  END_FUNC_VAL (retval);
}

/*-------------------------------------------------------------------------------
 * Parallel timing functions; these use MPI_Wtime() to record 
 * wall-clock time only.
 *-------------------------------------------------------------------------------*/

#elif defined(MPI_TIMING)

#undef __FUNC__
#define __FUNC__ "Timer_dhStart"
void
Timer_dhStart (Timer_dh t)
{
  START_FUNC_DH t->begin_wall = MPI_Wtime ();
  t->isRunning = true;
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "Timer_dhStop"
void
Timer_dhStop (Timer_dh t)
{
  START_FUNC_DH t->end_wall = MPI_Wtime ();
  t->isRunning = false;
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "Timer_dhReadWall"
double
Timer_dhReadWall (Timer_dh t)
{
  START_FUNC_DH double retval;
  if (t->isRunning)
    t->end_wall = MPI_Wtime ();
  retval = t->end_wall - t->begin_wall;
END_FUNC_VAL (retval)}

#undef __FUNC__
#define __FUNC__ "Timer_dhReadCPU"
double
Timer_dhReadCPU (Timer_dh t)
{
START_FUNC_DH END_FUNC_VAL (-1.0)}

#undef __FUNC__
#define __FUNC__ "Timer_dhReadUsage"
double
Timer_dhReadUsage (Timer_dh t)
{
  START_FUNC_DH END_FUNC_VAL (-1.0);
}


/*-------------------------------------------------------------------------------
 * junk timing methods -- these do nothing!
 *-------------------------------------------------------------------------------*/

#else

#undef __FUNC__
#define __FUNC__ "Timer_dhStart"
void
Timer_dhStart (Timer_dh t)
{
START_FUNC_DH END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "Timer_dhStop"
void
Timer_dhStop (Timer_dh t)
{
START_FUNC_DH END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "Timer_dhReadWall"
double
Timer_dhReadWall (Timer_dh t)
{
START_FUNC_DH END_FUNC_VAL (-1.0)}

#undef __FUNC__
#define __FUNC__ "Timer_dhReadCPU"
double
Timer_dhReadCPU (Timer_dh t)
{
START_FUNC_DH END_FUNC_VAL (-1.0)}

#undef __FUNC__
#define __FUNC__ "Timer_dhReadUsage"
double
Timer_dhReadUsage (Timer_dh t)
{
  START_FUNC_DH END_FUNC_VAL (-1.0);
}

#endif
