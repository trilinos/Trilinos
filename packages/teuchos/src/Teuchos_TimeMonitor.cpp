// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_ConfigDefs.hpp"

using namespace Teuchos;

Array<RefCountPtr<Time> > TimeMonitor::timers_;

RefCountPtr<Time> TimeMonitor::getNewTimer(const string& name)
{
  RefCountPtr<Time> rtn = rcp(new Time(name), true);
  timers_.append(rtn);
  return rtn;
}


void TimeMonitor::summarize(ostream &out)
{
  Array<string> names(timers_.length());
  Array<double> timings(timers_.length());

  for (int i=0; i<timers_.length(); i++)
    {
      names[i] = timers_[i]->name();
      timings[i] = timers_[i]->totalElapsedTime();
    }

  int np=1;
  int rank=0;
#ifdef HAVE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if (np==1)
    {
      for (int i=0; i<names.length(); i++)
        {
          out << std::left << std::setw(40) << names[i]
              << ": " << timings[i] << endl;
        }
    }
  else
    {
      Array<double> minTime(timers_.length());
      Array<double> maxTime(timers_.length());
      Array<double> avgTime(timers_.length());
      gatherTimings(timings, minTime, avgTime, maxTime);
      if (rank==0)
        {
          for (int i=0; i<names.length(); i++)
            {
              out << std::left << std::setw(30) << names[i]
                  << ": " << std::left << std::setw(12) << minTime[i]
                  << std::left << std::setw(12) << avgTime[i]
                  << std::left << std::setw(12) << maxTime[i] << endl;
            }
        }
    }

}

void TimeMonitor::gatherTimings(const Array<double>& timings,
                                Array<double>& minTime,
                                Array<double>& avgTime,
                                Array<double>& maxTime)
{
#ifdef HAVE_MPI
  int np;
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  void* tPtr = (void*) &(timings[0]);
  void* minPtr = (void*) &(minTime[0]);
  void* avgPtr = (void*) &(avgTime[0]);
  void* maxPtr = (void*) &(maxTime[0]);

  int count = (int) timings.length();

  MPI_Allreduce(tPtr, minPtr, count, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(tPtr, avgPtr, count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(tPtr, maxPtr, count, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  for (int i=0; i<avgTime.length(); i++)
    {
      avgTime[i] = avgTime[i]/((double) np);
    }
#endif
}

