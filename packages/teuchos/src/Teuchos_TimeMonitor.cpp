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
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_TableColumn.hpp"
#include "Teuchos_TableFormat.hpp"
//#include "Teuchos_MPIContainerComm.hpp"

namespace Teuchos {

  TimeMonitor::TimeMonitor (Time& timer, bool reset) 
    : PerformanceMonitorBase<Time>(timer, reset)
  {
    if (!isRecursiveCall()) counter().start(reset);
  }

  TimeMonitor::~TimeMonitor() {
    if (!isRecursiveCall()) counter().stop();
  }

  void 
  TimeMonitor::zeroOutTimers()
  {
    const Array<RCP<Time> > timers = counters();
  
    typedef Array<RCP<Time> >::size_type size_type;
    const size_type numTimers = timers.size();
    for (size_type i = 0; i < numTimers; ++i) 
      {
	Time &timer = *timers[i];
#ifdef TEUCHOS_DEBUG
	TEST_FOR_EXCEPTION(timer.isRunning(), std::logic_error,
			   "The timer i = " << i << " with name \"" 
			   << timer.name() << "\" is currently running and may not "
			   "be reset.");
#endif // TEUCHOS_DEBUG
	timer.reset();
      }
  }

  void
  TimeMonitor::filterZeroData (Array<timer_datum_t>& timerData)
  {
    Array<timer_datum_t> newTimerData;
    for (Array<timer_datum_t>::const_iterator it = timerData.begin(); 
	 it != timerData.end(); ++it)
      {
	if (it->second.second > 0)
	  newTimerData.push_back (*it);
      }
    timerData.swap (newTimerData);
  }

  void
  TimeMonitor::mergeTimers (const Comm<int>& comm, 
			    const Array<timer_datum_t>& localTimerData,
			    Array<timer_datum_t>& globalTimerData,
			    const bool intersect)
  {
    const int myRank = comm.getRank ();
    const int left = 0;
    const int right = comm.getSize() - 1;
    Array<timer_datum_t> theGlobalTimerData;
    mergeTimersHelper (comm, myRank, left, right, 
		       localTimerData, theGlobalTimerData, intersect);
    // "Transactional" semantics ensure strong exception safety for
    // output.
    globalTimerData.swap (theGlobalTimerData);
  }

  void
  TimeMonitor::mergeTimersHelper (const Comm<int>& comm, 
				  const int myRank,
				  const int left,
				  const int right, // inclusive range [left, right]
				  const Array<timer_datum_t>& localTimerData,
				  Array<timer_datum_t>& globalTimerData,
				  const bool intersect)
  {
    // Correctness proof:
    //
    // 1. Both set intersection and set union are associative (and
    //    indeed even commutative) operations.
    // 2. mergeTimersHelper() is just a reduction by binary tree.
    // 3. Reductions may use any tree shape as long as the binary
    //    operation is associative.
    //
    // Recursive "reduction" algorithm:
    //
    // If the (intersection, union) of [left, mid-1] and the
    // (intersection, union) of [mid, right] are both computed
    // correctly, then the (intersection, union) of these two sets is
    // the union of [left, right].
    //
    // The first base case is left == right: the (intersection, union)
    // of one set is simply that set.  For safety, we include a second
    // base case left > right, which means an empty interval: the
    // (intersection, union) of an empty set of sets is the empty set.
    if (left > right)
      return;
    else if (left == right)
      {
	globalTimerData.resize (localTimerData.length());
	std::copy (localTimerData.begin(), localTimerData.end(), 
		   globalTimerData.begin());
      }
    else
      { // You're sending messages across the network, so don't bother
	// to optimize away a few branches here.
	//
	// Recurse on [left, mid-1] or [mid, right], depending on myRank.
	const int mid = (right - left + 1) / 2;
	if (myRank >= left && myRank <= mid-1)
	  mergeTimersHelper (comm, myRank, left, mid-1, 
			     localTimerData, globalTimerData, intersect);
	else if (myRank >= mid && myRank <= right)
	  mergeTimersHelper (comm, myRank, mid, right,
			     localTimerData, globalTimerData, intersect);
	// Combine the results of the recursive step.
	if (myRank == left || myRank == mid)
	  mergeTimersPair (comm, myRank, left, mid, 
			   localTimerData, globalTimerData, intersect);
      }
  }

  namespace {
    //
    // Local utility functions for operating on timer_datum_t objects.
    //
    inline bool
    lessTimerData (const std::pair<std::string, std::pair<double, int> >& x,
		   const std::pair<std::string, std::pair<double, int> >& y)
    {
      return x.first < y.first;
    }

    inline bool
    eqTimerData (const std::pair<std::string, std::pair<double, int> >& x,
		 const std::pair<std::string, std::pair<double, int> >& y)
    {
      return x.first == y.first;
    }

    inline void
    mergeTimerData (std::pair<std::string, std::pair<double, int> >& x,
		    const std::pair<std::string, std::pair<double, int> >& y)
    {
      x.second.first += y.second.first;
      x.second.second += y.second.second;
    }

    inline const std::string& 
    timerDatumName (const std::pair<std::string, std::pair<double, int> >& x)
    {
      return x.first;
    }
  } // namespace (anonymous)

  void
  TimeMonitor::mergeTimersPair (const Comm<int>& comm, 
				const int myRank,
				const int left,
				const int mid,
				const Array<timer_datum_t>& localTimerData,
				Array<timer_datum_t>& globalTimerData,
				const bool intersect)
  {
    using std::string;

    if (myRank == left)
      { // Receive timer data from the other process, and merge its
	// timers with the timers on this process.
	Array<timer_datum_t> otherTimerData;
	(void) receive (comm, mid, &otherTimerData);

	// Assume that the other process' timer data is (jointly)
	// sorted, in alphabetical order by timer names.  Assume also
	// that our timer data (as stored in (names, timings, calls)) is
	// sorted in the same way.  Compute the set intersection / union
	// as specified.
	Array<timer_datum_t> newTimerData;
	if (intersect)
	  std::set_intersection (globalTimerData.begin(), globalTimerData.end(),
				 otherTimerData.begin(), otherTimerData.end(),
				 std::back_inserter (newTimerData));
	else
	  std::set_union (globalTimerData.begin(), globalTimerData.end(),
			  otherTimerData.begin(), otherTimerData.end(),
			  std::back_inserter (newTimerData));
	globalTimerData.swap (newTimerData);
      }
    else if (myRank == mid)
      send (comm, localTimerData, left);
    else
      TEST_FOR_EXCEPTION(myRank != left && myRank != mid, 
			 std::logic_error,
			 "Should never get here! myRank=" << myRank 
			 << " is neither left=" << left << " nor mid=" 
			 << mid << ".");
  }

  void
  TimeMonitor::collectLocalTimerData (Array<timer_datum_t>& localData)
  {
    using std::make_pair;
    using std::pair;
    using std::string;

    Array<RCP<Time> > theLocalCounters = counters();
    const Array<RCP<Time> >::size_type numLocalCounters = theLocalCounters.length();
  
    Array<timer_datum_t> theLocalData;
    theLocalData.reserve (numLocalCounters);
    for (Array<RCP<Time> >::const_iterator it = theLocalCounters.begin();
	 it != theLocalCounters.end(); ++it)
      {
	const string& name = (*it)->name();
	const double timing = (*it)->totalElapsedTime();
	const int numCalls = (*it)->numCalls();
	theLocalData.push_back (make_pair (name, make_pair (timing, numCalls)));
      }

    // Sort the local timer data in place, using timer name as the key.
    std::sort (theLocalData.begin(), theLocalData.end(), lessTimerData);

    // Merge timers with duplicate labels by summing their total elapsed
    // times and call counts.
    Array<timer_datum_t>::iterator first = theLocalData.begin();
    Array<timer_datum_t>::iterator second = first;
    while (first != theLocalData.end() && second != theLocalData.end())
      {
	++second;
	if (second != theLocalData.end())
	  {
	    if (eqTimerData (*first, *second))
	      mergeTimerData (*first, *second);
	    else
	      ++first;
	  }
	else 
	  // "second" is at the (old) end of the data.
	  // Make "first" mark the new end.
	  ++first;
      }
    Array<timer_datum_t> finalLocalData;
    finalLocalData.reserve (theLocalData.length());
    std::copy (theLocalData.begin(), first, std::back_inserter (finalLocalData));

    // This avoids copying the array yet another time, and also gives
    // a nice "transactional" semantics to this method (for strong
    // exception safety).
    localData.swap (finalLocalData);
  }

  void 
  TimeMonitor::summarize (std::ostream &out,
			  const bool alwaysWriteLocal,
			  const bool writeGlobalStats,
			  const bool writeZeroTimers,
			  const bool globalUnionOfTimers)
  {
    using std::string;

    // The default "MPI_COMM_WORLD" communicator.
    RCP<const Comm<int> > pComm = DefaultComm<int>::getComm ();
    const int numProcs = pComm->getSize();
    const int myRank = pComm->getRank();

    // This is subject to revision, since we can't print both global
    // and local stats in the same table if they don't share the same
    // set of timer labels.  In any case, only MPI Rank 0 prints local
    // data.
    bool writeLocalStats = alwaysWriteLocal && myRank == 0;

    // Collect and sort local timer data by timer names.
    Array<timer_datum_t> localTimerData;
    collectLocalTimerData (localTimerData);

    // globalTimerData is only valid if writeGlobalStats is true.
    Array<timer_datum_t> globalTimerData; 
    if (writeGlobalStats)
      { // This does the correct and inexpensive thing if numProcs==1.
	mergeTimers (*pComm, localTimerData, globalTimerData, 
		     ! globalUnionOfTimers);
	if (! writeZeroTimers)
	  filterZeroData (globalTimerData);
      }
    // Wait until after merging local into global timer data before
    // filtering out local timers with zero call counts.  This ensures
    // that if at least one MPI process had a nonzero call count for a
    // particular timer, it won't be filtered out.  We only need to do
    // this work if we writing local stats.  Note also that the
    // filtering step may make the set of local timers different than
    // the set of global timers, which would prevent us from printing
    // the local stats.
    //
    // A smarter thing to do would be to let the filtering of local
    // timers be guided by the set of filtered global timers, but that
    // would be more complicated and would not help the common case
    // (in which the local and global sets of timers, as well as their
    // call counts, are the same).
    if ((alwaysWriteLocal || ! writeGlobalStats) && ! writeZeroTimers)
      filterZeroData (localTimerData);

    // Extract the timer names (global or local) into a single array
    // of strings, representing the column labels in the table.
    Array<string> timerNames;
    timerNames.reserve (globalTimerData.length());
    if (writeGlobalStats)
      {
	std::transform (globalTimerData.begin(), globalTimerData.end(), 
			std::back_inserter (timerNames), timerDatumName);
	// If the set of local timer names doesn't match the set of
	// global timer names, we can't print both in the same table.
	// Check this first!  In this case, rather than throwing an
	// exception, we ignore the "alwaysWriteLocal" option and only
	// print global timer data.
	if (alwaysWriteLocal && 
	    ! std::equal (localTimerData.begin(), localTimerData.end(),
			  globalTimerData.begin(), eqTimerData))
	  writeLocalStats = false;
      }
    else // Use local timer names as the column labels.
      std::transform (localTimerData.begin(), localTimerData.end(), 
		      std::back_inserter (timerNames), timerDatumName);

    const int precision = format().precision();

    // All columns of the table, in order.
    Array<TableColumn> tableColumns;

    // Labels of all the columns of the table.
    // We will append to this when we add each column.
    Array<string> titles;

    // Widths (in number of characters) of each column.
    // We will append to this when we add each column.
    Array<int> columnWidths;

    // Table column containing all timer labels.
    {
      titles.append ("Timer Name");
      TableColumn nameCol (timerNames);
      tableColumns.append (nameCol);

      // Each column is as wide as it needs to be to hold both its title
      // and all of the column data.  This column's title is the current
      // last entry of the titles array.
      columnWidths.append (format().computeRequiredColumnWidth (titles.back(), 
								nameCol));
    }

    // Table column containing local timer stats, if applicable.  We
    // only write local stats for MPI Proc 0, only if asked, and only if
    // the list of local timer labels matches the list of global timer
    // labels.  (The latter is a common case, but not the only
    // possibility.)
    if (writeLocalStats)
      {
	titles.append ("Local time (num calls)");

	// Copy local timer data out of the array-of-structs into
	// separate arrays, for display in the table.
	Array<double> localTimings;
	Array<double> localNumCalls;
	localTimings.reserve (localTimerData.length());
	localNumCalls.reserve (localTimerData.length());
	for (Array<timer_datum_t>::const_iterator it = localTimerData.begin();
	     it != localTimerData.end(); ++it)
	  {
	    localTimings.push_back (it->second.first);
	    localNumCalls.push_back (static_cast<double> (it->second.second));
	  }
	TableColumn timeAndCalls (localTimings, localNumCalls, precision, true);
	tableColumns.append (timeAndCalls);
	columnWidths.append (format().computeRequiredColumnWidth (titles.back(), 
								  timeAndCalls));
      }

    if (writeGlobalStats)
      {
	const Array<timer_datum_t>::size_type numGlobalTimers = 
	  globalTimerData.length();

	// Copy global timer data out of the array-of-structs into
	// separate arrays, for display in the table and/or for
	// computing statistics.
	Array<double> globalTimings;
	Array<double> globalNumCalls;
	globalTimings.reserve (globalTimerData.length());
	globalNumCalls.reserve (globalTimerData.length());
	for (Array<timer_datum_t>::const_iterator it = globalTimerData.begin();
	     it != globalTimerData.end(); ++it)
	  {
	    globalTimings.push_back (it->second.first);
	    globalNumCalls.push_back (static_cast<double> (it->second.second));
	  }

	if (numProcs == 1)
	  { // Don't display statistics in the case of only 1 MPI process.
	    titles.append ("Global time (num calls)");
	    TableColumn timeAndCalls (globalTimings, globalNumCalls, precision, true);
	    tableColumns.append (timeAndCalls);
	    columnWidths.append (format().computeRequiredColumnWidth (titles.back(), 
								      timeAndCalls));
	  }
	else // numProcs > 1
	  {
	    // Table column containing min of global timer stats, if
	    // applicable.
	    {
	      titles.append ("Min over procs");
	      Array<double> minGlobalTimings (numGlobalTimers);
	      Array<double> minGlobalNumCalls (numGlobalTimers);

	      // Teuchos_CommHelpers.hpp doesn't currently have a reduce(); it
	      // only has a reduceAll().  It would be better just to reduce to
	      // Proc 0, but that has to wait until reduce() is implemented.
	      if (numGlobalTimers > 0)
		{
		  reduceAll (*pComm, REDUCE_MIN, 
			     static_cast<int> (numGlobalTimers), 
			     &globalTimings[0], &minGlobalTimings[0]);
		  reduceAll (*pComm, REDUCE_MIN, 
			     static_cast<int> (numGlobalTimers),
			     &globalNumCalls[0], &minGlobalNumCalls[0]);
		}
	      TableColumn timeAndCalls (minGlobalTimings, minGlobalNumCalls, 
					precision, true);
	      tableColumns.append (timeAndCalls);
	      columnWidths.append (format().computeRequiredColumnWidth (titles.back(), 
									timeAndCalls));
	    }
	    // Table column containing arithmetic mean of global timer
	    // stats, if applicable.
	    {
	      titles.append ("Mean over procs");

	      // Scale first, so that the reduction can sum.  This avoids
	      // unnecessary overflow, in case the sum is large but the number
	      // of processors is also large.
	      Array<double> scaledGlobalTimings (numGlobalTimers);
	      Array<double> scaledGlobalNumCalls (numGlobalTimers);
	      std::transform (globalTimings.begin(), globalTimings.end(), 
			      scaledGlobalTimings.begin(), 
			      std::bind2nd (std::divides<double>(), 
					    static_cast<double> (numProcs)));
	      std::transform (globalNumCalls.begin(), globalNumCalls.end(), 
			      scaledGlobalNumCalls.begin(), 
			      std::bind2nd (std::divides<double>(), 
					    static_cast<double> (numProcs)));
	      Array<double> avgGlobalTimings (numGlobalTimers);
	      Array<double> avgGlobalNumCalls (numGlobalTimers);
	      if (numGlobalTimers > 0)
		{
		  reduceAll (*pComm, REDUCE_SUM, 
			     static_cast<int> (numGlobalTimers),
			     &scaledGlobalTimings[0], &avgGlobalTimings[0]);
		  reduceAll (*pComm, REDUCE_SUM, 
			     static_cast<int> (numGlobalTimers),
			     &scaledGlobalNumCalls[0], &avgGlobalNumCalls[0]);
		}
	      TableColumn timeAndCalls (avgGlobalTimings, avgGlobalNumCalls, 
					precision, true);
	      tableColumns.append (timeAndCalls);
	      columnWidths.append (format().computeRequiredColumnWidth (titles.back(), 
									timeAndCalls));
	    }

	    // Table column containing max of global timer stats, if
	    // applicable.
	    {
	      titles.append("Max over procs");
	      Array<double> maxGlobalTimings (numGlobalTimers);
	      Array<double> maxGlobalNumCalls (numGlobalTimers);

	      // Teuchos_CommHelpers.hpp doesn't currently have a reduce(); it
	      // only has a reduceAll().  It would be better just to reduce to
	      // Proc 0, but that has to wait until reduce() is implemented.
	      if (numGlobalTimers > 0)
		{
		  reduceAll (*pComm, REDUCE_MAX, 
			     static_cast<int> (numGlobalTimers),
			     &globalTimings[0], &maxGlobalTimings[0]);
		  reduceAll (*pComm, REDUCE_MAX, 
			     static_cast<int> (numGlobalTimers),
			     &globalNumCalls[0], &maxGlobalNumCalls[0]);
		}
	      TableColumn timeAndCalls (maxGlobalTimings, maxGlobalNumCalls, 
					precision, true);
	      tableColumns.append (timeAndCalls);
	      columnWidths.append (format ().computeRequiredColumnWidth (titles.back(), 
									 timeAndCalls));
	    }
	  }
      }

    // Print the whole table to the given output stream on MPI Rank 0.
    format ().setColumnWidths (columnWidths);
    if (myRank == 0)
      format().writeWholeTable (out, "TimeMonitor Results", titles, tableColumns);
  }


} // namespace Teuchos
