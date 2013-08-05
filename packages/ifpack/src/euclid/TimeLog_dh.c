/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
