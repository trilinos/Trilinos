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

  typedef std::map<std::string, std::pair<double, int> > timer_map_t;

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
  
    // In debug mode, loop first to check whether any of the timers
    // are running, before resetting them.  This ensures that this
    // method satisfies the strong exception guarantee (either it
    // completes normally, or there are no side effects).
#ifdef TEUCHOS_DEBUG
    typedef Array<RCP<Time> >::size_type size_type;
    const size_type numTimers = timers.size();
    for (size_type i = 0; i < numTimers; ++i) 
      {
	Time &timer = *timers[i];
	TEST_FOR_EXCEPTION(timer.isRunning(), std::logic_error,
			   "The timer i = " << i << " with name \"" 
			   << timer.name() << "\" is currently running and may not "
			   "be reset.");
      }
#endif // TEUCHOS_DEBUG

    for (Array<RCP<Time> >::const_iterator it = timers.begin(); 
	 it != timers.end(); ++it)
      (*it)->reset ();
  }

  namespace {

    // \brief Helper function for \c mergeTimersHelper().
    //
    // mergeTimersHelper() implements the set union resp. intersection
    // (depending on the \c setOp argument) of the MPI process' sets
    // of timer names as a parallel reduction.  The \c
    // mergeTimersPair() function implements the associative operator
    // which computes the set union resp. intersection of two sets:
    // the "left" process' intermediate reduction result (global timer
    // names), and the "mid" process' local timer names.
    // 
    // \param comm [in] Communicator for which mergeTimersHelper() was
    //   called.
    //
    // \param myRank [in] Rank of the calling MPI process; must be
    //   either == left or == mid.
    //
    // \param left [in] The "left" input argument of
    //   mergeTimersHelper().
    //
    // \param mid [in] The value of "mid" in the implementation of
    //   mergeTimersHelper().
    //
    // \param localTimerNames [in] List of timer names belonging to
    //   the calling MPI process.
    //
    // \param globalTimerNames [in/out] Only accessed if myRank ==
    //   left.  If so, on input: the intermediate reduction result of
    //   the union resp. intersection (depending on \c setOp).  On
    //   output: the union resp. intersection of the input value of
    //   globalTimerNames with the other MPI process' localTimerNames.
    //
    // \param setOp [in] If Intersection, compute the set intersection
    //   of timer names, else if Union, compute the set union of timer
    //   names.
    void
    mergeTimersPair (const Comm<int>& comm, 
		     const int myRank,
		     const int left,
		     const int mid,
		     const Array<std::string>& localTimerNames,
		     Array<std::string>& globalTimerNames,
		     const ETimerSetOp setOp)
    {
      using std::string;

      if (myRank == left)
	{ // Receive timer names from the other process, and merge its
	  // names with the timer names on this process.
	  Array<string> otherTimerNames;
	  (void) receive (comm, mid, &otherTimerNames);

	  // Assume that both globalTimerNames and otherTimerNames are
	  // sorted.  Compute the set intersection / union as
	  // specified by the enum.
	  Array<string> newTimerNames;
	  if (setOp == Intersection)
	    std::set_intersection (globalTimerNames.begin(), globalTimerNames.end(),
				   otherTimerNames.begin(), otherTimerNames.end(),
				   std::back_inserter (newTimerNames));
	  else if (setOp == Union)
	    std::set_union (globalTimerNames.begin(), globalTimerNames.end(),
			    otherTimerNames.begin(), otherTimerNames.end(),
			    std::back_inserter (newTimerNames));
	  else
	    TEST_FOR_EXCEPTION(setOp != Intersection && setOp != Union,
			       std::logic_error,
			       "Invalid set operation enum value.  Please "
			       "report this bug to the Teuchos developers.");
	  globalTimerNames.swap (newTimerNames);
	}
      else if (myRank == mid)
	send (comm, localTimerNames, left);
      else
	TEST_FOR_EXCEPTION(myRank != left && myRank != mid, 
			   std::logic_error,
			   "myRank=" << myRank << " is neither left=" << left
			   << " nor mid=" << mid << ".  Please report this "
			   "bug to the Teuchos developers.");
    }

    // Recursive helper function for \c mergeTimers().
    //
    // mergeTimersHelper() implements the set union resp. intersection
    // (depending on the \c setOp argument) of the MPI process' sets
    // of timer names as a parallel reduction. (Since the Teuchos comm
    // wrappers as of 11 July 2011 lack a wrapper for MPI_Reduce(), we
    // hand-roll the reduction using a binary tree via recursion.  We
    // don't need an all-reduce in this case.)
    void
    mergeTimersHelper (const Comm<int>& comm, 
		       const int myRank,
		       const int left,
		       const int right, // inclusive range [left, right]
		       const Array<std::string>& localTimerNames,
		       Array<std::string>& globalTimerNames,
		       const ETimerSetOp setOp)
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
	  globalTimerNames.resize (localTimerNames.length());
	  std::copy (localTimerNames.begin(), localTimerNames.end(), 
		     globalTimerNames.begin());
	}
      else
	{ // You're sending messages across the network, so don't bother
	  // to optimize away a few branches here.
	  //
	  // Recurse on [left, mid-1] or [mid, right], depending on myRank.
	  const int mid = (right - left + 1) / 2;
	  if (myRank >= left && myRank <= mid-1)
	    mergeTimersHelper (comm, myRank, left, mid-1, 
			       localTimerNames, globalTimerNames, setOp);
	  else if (myRank >= mid && myRank <= right)
	    mergeTimersHelper (comm, myRank, mid, right,
			       localTimerNames, globalTimerNames, setOp);
	  // Combine the results of the recursive step.
	  if (myRank == left || myRank == mid)
	    mergeTimersPair (comm, myRank, left, mid, 
			     localTimerNames, globalTimerNames, setOp);
	}
    }

    /// \brief Merge timer data over all processors.
    //
    // \param comm [in] Communicator over which to merge.
    // \param localTimerData [in] Each processor's timer data.
    // \param globalTimerData [out] On output, on MPI Proc 0: the
    //   results of merging the timer data.
    // \param setOp [in] If Intersection, globalTimerData on output
    //   contains the intersection of all timers.  If Union,
    //   globalTimerData on output contains the union of all timers.
    void
    mergeTimers (const Comm<int>& comm, 
		 const Array<std::string>& localTimerNames,
		 Array<std::string>& globalTimerNames,
		 const ETimerSetOp setOp)
    {
      const int myRank = comm.getRank ();
      const int left = 0;
      const int right = comm.getSize() - 1;
      Array<std::string> theGlobalTimerNames;
      mergeTimersHelper (comm, myRank, left, right, 
			 localTimerNames, theGlobalTimerNames, setOp);
      // "Transactional" semantics ensure strong exception safety for
      // output.
      globalTimerNames.swap (theGlobalTimerNames);
    }

    //
    // \brief Locally filter out timer data with zero call counts.
    //
    void
    filterZeroData (std::map<std::string, std::pair<double, int> >& timerData)
    {
      timer_map_t newTimerData;
      for (timer_map_t::const_iterator it = timerData.begin(); 
	   it != timerData.end(); ++it)
	{
	  if (it->second.second > 0)
	    newTimerData[it->first] = it->second;
	}
      timerData.swap (newTimerData);
    }

    //
    // \brief Collect and sort local timer data by timer names.
    //
    void
    collectLocalTimerData (std::map<std::string, std::pair<double, int> >& localData,
			   const Array<RCP<Time> >& localCounters)
    {
      using std::make_pair;
      typedef timer_map_t::const_iterator const_iter_t;
      typedef timer_map_t::iterator iter_t;

      timer_map_t theLocalData;
      for (Array<RCP<Time> >::const_iterator it = localCounters.begin();
	   it != localCounters.end(); ++it)
	{
	  const std::string& name = (*it)->name();
	  const double timing = (*it)->totalElapsedTime();
	  const int numCalls = (*it)->numCalls();

	  // Merge timers with duplicate labels, by summing their
	  // total elapsed times and call counts.
	  iter_t loc = theLocalData.find (name);
	  if (loc == theLocalData.end())
	    // Use loc as an insertion location hint.
	    theLocalData.insert (loc, make_pair (name, make_pair (timing, numCalls)));
	  else
	    {
	      loc->second.first += timing;
	      loc->second.second += numCalls;
	    }
	}
      // This avoids copying the map, and also makes this method
      // satisfy the strong exception guarantee.
      localData.swap (theLocalData);
    }
  } // namespace (anonymous)

  namespace {
    //
    // Local utility functions for operating on timer_datum_t objects.
    //

    inline bool
    eqTimerData (const std::pair<std::string, std::pair<double, int> >& x,
		 const std::pair<std::string, std::pair<double, int> >& y)
    {
      return x.first == y.first;
    }

    inline const std::string& 
    timerDatumName (const std::pair<std::string, std::pair<double, int> >& x)
    {
      return x.first;
    }
  } // namespace (anonymous)

  void 
  TimeMonitor::summarize (std::ostream &out,
			  const bool alwaysWriteLocal,
			  const bool writeGlobalStats,
			  const bool writeZeroTimers,
			  const ETimerSetOp setOp)
  {
    using std::make_pair;
    using std::string;

    // The default "MPI_COMM_WORLD" communicator.
    RCP<const Comm<int> > pComm = DefaultComm<int>::getComm ();
    const int numProcs = pComm->getSize();
    const int myRank = pComm->getRank();

    // Only MPI Rank 0 prints local timer stats, if asked always to
    // print local data.
    bool writeLocalStats = alwaysWriteLocal && myRank == 0;

    // Collect and sort local timer data by timer names.
    timer_map_t localTimerData;
    collectLocalTimerData (localTimerData, counters());

    // Extract the set of local timer names.  The std::map keeps them
    // sorted alphabetically.
    Array<string> localTimerNames;
    std::transform (localTimerData.begin(), localTimerData.end(), 
		    std::back_inserter (localTimerNames), 
		    std::mem_fun (&timer_map_t::value_type::first));
    // for (timer_map_t::const_iterator it = localTimerData.begin(); 
    // 	 it != localTimerData.end(); ++it)
    //   localTimerNames.push_back (it->first);
      
    // globalTimerData and globalTimerNames are only valid if
    // writeGlobalStats is true.
    Array<string> globalTimerNames;
    timer_map_t globalTimerData;

    // If setOp == Union, there may be some global timers that are not
    // local timers.  In that case, if we want to print local timers,
    // we forbid filtering out zero data, for reasons discussed below.
    bool foundGlobalNotLocal = false;
    if (writeGlobalStats)
      { // This does the correct and inexpensive thing (just copies
	// the timer data) if numProcs == 1.
	mergeTimers (*pComm, localTimerNames, globalTimerNames, setOp);

	// mergeTimers() just merges timer names, not their actual
	// data.  Now we need to fill globalTimerData with this MPI
	// process' timer data for those timers in globalTimerNames.
	// If setOp == Union, there may be some global timers that are
	// not local timers.  In that case, if we want to print local
	// timers, we insert a local timer datum with zero elapsed
	// time and zero call count, and forbid filtering out zero
	// data.  Inserting the local timer datum makes sure that both
	// global and local timer columns in the output table have the
	// same number of rows.
	//
	// Insertion optimization; if the iterator given to
	// map::insert points right before where we want to insert,
	// insertion is O(1).  globalTimerNames is sorted, so feeding
	// the iterator output of map::insert into the next
	// invocation's input should make the whole insertion O(N)
	// where N is the number of entries in globalTimerNames.
	timer_map_t::iterator globalMapIter = globalTimerData.begin();
	for (Array<string>::const_iterator it = globalTimerNames.begin(); 
	     it != globalTimerNames.end(); ++it)
	  {
	    const std::string& globalName = *it;
	    localMapIter = localTimerData.find (globalName);
	    if (localMapIter == localTimerData.end())
	      {
		foundGlobalNotLocal = true;

		// We really only need to do the following on MPI Rank
		// 0, which is the only process that will print local
		// or global output.  However, we do it on all MPI
		// processes just in case someone later wants to
		// modify this function to print out local timer data
		// for some other MPI process other than Rank 0.
		std::pair<timer_map_t::iterator, bool> insertResult;
		if (alwaysWriteLocal)
		  {
		    insertResult = localTimerData.insert (localMapIter, make_pair (globalName, make_pair (double(0), int(0))));
		    localMapIter = insertResult.first;
		  }
	      }
	    else 
	      // We know that globalTimerNames is unique, so we don't
	      // have to worry about inserting an element whose key
	      // already exists in the global timer data map.
	      globalTimerData.insert (globalMapIter, make_pair (globalName, localMapIter->second));
	
	    // Don't filter out zero data if we had to insert a local
	    // timer datum corresponding to a global timer name with no
	    // local timer name.
	    if (! writeZeroTimers && ! foundGlobalNotLocal)
	      filterZeroData (globalTimerData);
	  }
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
    //
    // If writing local stats or not writing global stats, and if not
    // writing zero timers, and if (if writing global stats, we
    // didn't find any global-not-local timer names), filter zero data
    // from local timers.
    if ((alwaysWriteLocal || ! writeGlobalStats) && ! writeZeroTimers && (! writeGlobalStatus || foundGlobalNotLocal))
      filterZeroData (localTimerData);

    // Extract the timer names (global or local) into a single array
    // of strings, representing the column labels in the table.
    Array<string> timerNames (globalTimerData.size());
    if (writeGlobalStats)
      {
	std::copy (globalTimerNames.begin(), globalTimerNames.end(), timerNames.begin());

	// FIXME (mfh 15 Jul 2011) We've already dealt with this
	// above.  localTimerNames might be different than
	// globalTimerNames, but if writeGlobalStats==true and
	// alwaysWriteLocal==true, then localTimerData has had "dummy"
	// data inserted to pad out the table.
	if (false)
	  {
	    // If the set of local timer names doesn't match the set of
	    // global timer names, we can't print both in the same table.
	    // Check this first!  In this case, rather than throwing an
	    // exception, we ignore the "alwaysWriteLocal" option and only
	    // print global timer data.
	    if (alwaysWriteLocal && 
		! std::equal (localTimerNames.begin(), localTimerNames.end(), 
			      globalTimerNames.begin()))
	      writeLocalStats = false;
	  }
      }
    else // Use local timer names as the column labels.
      std::copy (localTimerNames.begin(), localTimerNames.end(), timerNames.begin());

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
	for (timer_map_t::const_iterator it = localTimerData.begin();
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
	const timer_map_t::size_type numGlobalTimers = globalTimerData.size();

	// Copy global timer data out of the array-of-structs into
	// separate arrays, for display in the table and/or for
	// computing statistics.
	Array<double> globalTimings;
	Array<double> globalNumCalls;
	globalTimings.reserve (globalTimerData.length());
	globalNumCalls.reserve (globalTimerData.length());
	for (timer_map_t::const_iterator it = globalTimerData.begin();
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
