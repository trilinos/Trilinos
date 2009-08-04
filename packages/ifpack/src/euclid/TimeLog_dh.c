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

#include "TimeLog_dh.h"
#include "Timer_dh.h"
#include "Mem_dh.h"

#define MAX_TIME_MARKS  100
#define MAX_DESC_LENGTH 60

struct _timeLog_dh
{
  int first;
  int last;
  double time[MAX_TIME_MARKS];
  char desc[MAX_TIME_MARKS][MAX_DESC_LENGTH];
  Timer_dh timer;
};

#undef __FUNC__
#define __FUNC__ "TimeLog_dhCreate"
void
TimeLog_dhCreate (TimeLog_dh * t)
{
  START_FUNC_DH int i;
  struct _timeLog_dh *tmp =
    (struct _timeLog_dh *) MALLOC_DH (sizeof (struct _timeLog_dh));
  CHECK_V_ERROR;
  *t = tmp;
  tmp->first = tmp->last = 0;
  Timer_dhCreate (&tmp->timer);
  for (i = 0; i < MAX_TIME_MARKS; ++i)
    strcpy (tmp->desc[i], "X");
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "TimeLog_dhDestroy"
void
TimeLog_dhDestroy (TimeLog_dh t)
{
  START_FUNC_DH Timer_dhDestroy (t->timer);
  FREE_DH (t);
END_FUNC_DH}


#undef __FUNC__
#define __FUNC__ "TimeLog_dhStart"
void
TimeLog_dhStart (TimeLog_dh t)
{
  START_FUNC_DH Timer_dhStart (t->timer);
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "TimeLog_dhStop"
void
TimeLog_dhStop (TimeLog_dh t)
{
  START_FUNC_DH Timer_dhStop (t->timer);
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "TimeLog_dhMark"
void
TimeLog_dhMark (TimeLog_dh t, char *desc)
{
  START_FUNC_DH if (t->last < MAX_TIME_MARKS - 3)
    {
/*     SET_V_ERROR("overflow; please increase MAX_TIME_MARKS and recompile"); */
      Timer_dhStop (t->timer);
      t->time[t->last] = Timer_dhReadWall (t->timer);
      Timer_dhStart (t->timer);
      sprintf (t->desc[t->last], desc);
      t->last += 1;
    }
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "TimeLog_dhReset"
void
TimeLog_dhReset (TimeLog_dh t)
{
  START_FUNC_DH if (t->last < MAX_TIME_MARKS - 2)
    {
      double total = 0.0;
      int i, first = t->first, last = t->last;
      for (i = first; i < last; ++i)
	total += t->time[i];
      t->time[last] = total;
      sprintf (t->desc[last], "========== totals, and reset ==========\n");
      t->last += 1;
      t->first = t->last;
      Timer_dhStart (t->timer);
    }
END_FUNC_DH}


#undef __FUNC__
#define __FUNC__ "TimeLog_dhPrint"
void
TimeLog_dhPrint (TimeLog_dh t, FILE * fp, bool allPrint)
{
  START_FUNC_DH int i;
  double total = 0.0;
  double timeMax[MAX_TIME_MARKS];
  double timeMin[MAX_TIME_MARKS];
  static bool wasSummed = false;


  if (!wasSummed)
    {
      for (i = t->first; i < t->last; ++i)
	total += t->time[i];
      t->time[t->last] = total;
      sprintf (t->desc[t->last], "========== totals, and reset ==========\n");
      t->last += 1;

      MPI_Allreduce (t->time, timeMax, t->last, MPI_DOUBLE, MPI_MAX, comm_dh);
      MPI_Allreduce (t->time, timeMin, t->last, MPI_DOUBLE, MPI_MIN, comm_dh);
      wasSummed = true;
    }

  if (fp != NULL)
    {
      if (myid_dh == 0 || allPrint)
	{
	  fprintf (fp,
		   "\n----------------------------------------- timing report\n");
	  fprintf (fp, "\n   self     max     min\n");
	  for (i = 0; i < t->last; ++i)
	    {
	      fprintf (fp, "%7.3f %7.3f %7.3f   #%s\n", t->time[i],
		       timeMax[i], timeMin[i], t->desc[i]);
	    }
	  fflush (fp);
	}
    }				/* if (fp != NULL) */
END_FUNC_DH}
