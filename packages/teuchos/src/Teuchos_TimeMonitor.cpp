// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER


#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_TableColumn.hpp"
#include "Teuchos_TableFormat.hpp"
#include "Teuchos_MPIContainerComm.hpp"


namespace Teuchos {


void TimeMonitor::zeroOutTimers()
{
  
  const Array<RCP<Time> > timers = counters();
  
  const int numTimers = timers.size();
  
  for( int i = 0; i < numTimers; ++i ) {
    Time &timer = *timers[i];
#ifdef TEUCHOS_DEBUG
    TEST_FOR_EXCEPTION(
      timer.isRunning(), std::logic_error,
      "Teuchos::TimeMonitor::zeroOutTimers():\n\n"
      "Error, the timer i = " << i << " with name \"" << timer.name() << "\""
      " is current running and not not be set to zero!"
      );
#endif
    timer.reset();
  }
  
}


void TimeMonitor::summarize(
  std::ostream &out,
  const bool alwaysWriteLocal,
  const bool writeGlobalStats,
  const bool writeZeroTimers
  )
{

  Array<std::string> localNames(counters().length());
  Array<double> localTimings(counters().length());
  Array<double> localCallCounts(counters().length());

  for (int i=0; i<counters().length(); i++)
  {
    localNames[i] = counters()[i]->name();
    localTimings[i] = counters()[i]->totalElapsedTime();
    localCallCounts[i] = counters()[i]->numCalls();
  }
  
  // Gather timings from all procs, in case some timers have been activated on
  // other processors but not on this one.
  Array<std::string> names;
  Array<Array<double> > data(2);
  PerformanceMonitorUtils::synchValues(MPIComm::world(), localNames, 
    tuple(localTimings, localCallCounts),
    names, data);
  
  Array<double>& timings = data[0];
  Array<double>& calls = data[1];

  // Remove zero rows
  if (!writeZeroTimers) {
    Teuchos::Array<std::string>::iterator nms=names.end();
    Teuchos::Array<double>::iterator tms=timings.end();
    Teuchos::Array<double>::iterator cls=calls.end();
    for (int i=names.length()-1; i>=0; i--)  {
      cls--; nms--; tms--;
      if (calls[i]<0.1) {
        cls = calls.erase(cls);
        nms = names.erase(nms);
        tms = timings.erase(tms);
      }
    }
  }
  
  // Form the table data
  MPIComm comm = MPIComm::world();
  int np = comm.getNProc();

  int precision = format().precision();

  Array<TableColumn> columnsToWrite;

  TableColumn nameCol(names);
  Array<std::string> titles;
  titles.append("Timer Name");

  columnsToWrite.append(nameCol);

  Array<int> columnWidths;
  columnWidths.append(format().computeRequiredColumnWidth(titles[titles.size()-1], 
      nameCol));

	  if (np==1 || alwaysWriteLocal)
	  {
	    TableColumn timeAndCalls(timings, calls, precision, true);
	    titles.append("Local time (num calls)");
	    columnsToWrite.append(timeAndCalls);
	    columnWidths.append(format().computeRequiredColumnWidth(titles[titles.size()-1], 
		timeAndCalls));
	  }
	  
	  if (np > 1 && writeGlobalStats)
	  {

	    titles.append("Min over procs");
	      
	    Array<double> minTimings;
	    PerformanceMonitorUtils::reduce(comm, EMin, timings, minTimings);
	      
	    Array<double> minCalls;
	    PerformanceMonitorUtils::reduce(comm, EMin, calls, minCalls);

	    TableColumn timeAndCalls(minTimings, minCalls, precision, true);
	    columnsToWrite.append(timeAndCalls);
	      
	    columnWidths.append(format().computeRequiredColumnWidth(titles[titles.size()-1], 
		timeAndCalls));

	  }
	  
	  if (np > 1 && writeGlobalStats)
	  {

	    titles.append("Avg over procs");
	      
	    Array<double> avgTimings;
	    PerformanceMonitorUtils::reduce(comm, EAvg, timings, avgTimings);
	      
	    Array<double> avgCalls;
	    PerformanceMonitorUtils::reduce(comm, EAvg, calls, avgCalls);

	    TableColumn timeAndCalls(avgTimings, avgCalls, precision, true);
	    columnsToWrite.append(timeAndCalls);
	      
	    columnWidths.append(format().computeRequiredColumnWidth(titles[titles.size()-1], 
		timeAndCalls));

	  }
	  
	  if (np > 1 && writeGlobalStats)
	  {

	    titles.append("Max over procs");
	      
	    Array<double> maxTimings;
	    PerformanceMonitorUtils::reduce(comm, EMax, timings, maxTimings);
	      
	    Array<double> maxCalls;
	    PerformanceMonitorUtils::reduce(comm, EMax, calls, maxCalls);

	    TableColumn timeAndCalls(maxTimings, maxCalls, precision, true);
	    columnsToWrite.append(timeAndCalls);
	      
	    columnWidths.append(format().computeRequiredColumnWidth(titles[titles.size()-1], 
		timeAndCalls));

	  }

	  format().setColumnWidths(columnWidths);

	  const bool writeOnThisProcessor = ( comm.getRank()==0 || alwaysWriteLocal );
	  if (writeOnThisProcessor)
	  {
	    format().writeWholeTable(out, "TimeMonitor Results",
	      titles, columnsToWrite);
	  }

	}


} // namespace Teuchos
