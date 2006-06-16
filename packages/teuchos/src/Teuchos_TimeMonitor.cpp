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
#include "Teuchos_TableColumn.hpp"
#include "Teuchos_TableFormat.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_MPIContainerComm.hpp"
#include "Teuchos_ConfigDefs.hpp"
using namespace Teuchos;



void TimeMonitor::summarize(ostream &out, bool alwaysWriteLocal,
                            bool writeGlobalStats)
{
  Array<string> localNames(counters().length());
  Array<double> localTimings(counters().length());
  Array<double> localCallCounts(counters().length());

  for (int i=0; i<counters().length(); i++)
    {
      localNames[i] = counters()[i]->name();
      localTimings[i] = counters()[i]->totalElapsedTime();
      localCallCounts[i] = counters()[i]->numCalls();
    }
  
  /* Gather timings from all procs, in case some timers have been activated
   * on other processors but not on this one.  */
  Array<string> names;
  Array<Array<double> > data(2);
  PerformanceMonitorUtils::synchValues(MPIComm::world(), localNames, 
                                       tuple(localTimings, localCallCounts),
                                       names, data);
  
  const Array<double>& timings = data[0];
  const Array<double>& calls = data[1];


  
  /* form the table data */
  MPIComm comm = MPIComm::world();
  int np = comm.getNProc();

  int precision = format().precision();

  Array<TableColumn> columnsToWrite;

  TableColumn nameCol(names);
  Array<string> titles;
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

  RefCountPtr<std::ostream> outPtr = rcp(&out, false);
  bool writeOnThisProcessor = ( comm.getRank()==0 || alwaysWriteLocal );
  if (writeOnThisProcessor)
    {
      format().writeWholeTable(outPtr, "TimeMonitor Results",
                               titles, columnsToWrite);
    }
}
