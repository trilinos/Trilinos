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
#include <functional>

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
	// We throw a runtime_error rather than a logic_error, because
	// logic_error suggests a bug in the implementation of
	// TimeMonitor.  Calling zeroOutTimers() when a timer is
	// running is not TimeMonitor's fault.
	TEST_FOR_EXCEPTION(timer.isRunning(), std::runtime_error,
			   "The timer i = " << i << " with name \"" 
			   << timer.name() << "\" is currently running and may not "
			   "be reset.");
      }
#endif // TEUCHOS_DEBUG

    for (Array<RCP<Time> >::const_iterator it = timers.begin(); 
	 it != timers.end(); ++it)
      (*it)->reset ();
  }

  // An anonymous namespace is the standard way of limiting linkage of
  // its contained routines to file scope.
  namespace {
    // \brief Return an "empty" local timer datum.
    // 
    // "Empty" means the datum has zero elapsed time and zero call
    // count.  This function does not actually create a timer.
    //
    // \param name The timer's name.
    std::pair<std::string, std::pair<double, int> >
    makeEmptyTimerDatum (const std::string& name)
    {
      return std::make_pair (name, std::make_pair (double(0), int(0)));
    }

    // \brief Locally filter out timer data with zero call counts.
    //
    // \param timerData [in/out]
    void
    filterZeroData (timer_map_t& timerData)
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
    collectLocalTimerData (timer_map_t& localData,
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

  void 
  TimeMonitor::summarize (std::ostream &out,
			  const bool alwaysWriteLocal,
			  const bool writeGlobalStats,
			  const bool writeZeroTimers,
			  const ECounterSetOp setOp)
  {
    using std::cerr;
    using std::endl;
    using std::make_pair;
    using std::string;
    
    const bool debug = false;

    // The default communicator.  If Trilinos was built with MPI
    // enabled, this should be MPI_COMM_WORLD.  Otherwise, this should
    // be a "serial" (no MPI, one "process") communicator.
    RCP<const Comm<int> > pComm = DefaultComm<int>::getComm ();

    // Callers may or may not have initialized MPI before calling
    // summarize().  Just because they built with MPI, doesn't mean
    // they want to use MPI.  It's not my responsibility to initialize
    // MPI for them, and I don't have the context I need in order to
    // do so anyway.  Thus, if Trilinos was built with MPI and MPI has
    // not yet been initialized, make pComm a "serial" communicator.
#ifdef HAVE_MPI
    {
      int mpiHasBeenStarted = 0;
      MPI_Initialized (&mpiHasBeenStarted);
      if (! mpiHasBeenStarted)
        {
	  // mfh 19 Jul 2011
	  //
	  // The first line commented out below compiles and runs
	  // correctly with GCC 4.5.1, but gives a compiler error with
	  // Intel's C++ compiler (version 11.1).  Intel's compiler
	  // also rejects the two commented-out lines that follow.  It
	  // seems that as of 19 July 2011, this code is the only code
	  // in Trilinos that uses getDefaultSerialComm().  Intel's
	  // compiler claims that the "Ordinal" template parameter of
	  // DefaultComm is really the Teuchos::Ordinal typedef; the
	  // Ordinal template parameter _should_ shadow the typedef in
	  // Teuchos, and does with GCC, but does not in Intel's
	  // compiler.  This may be the case with other compilers as
	  // well, but I haven't tested them yet.
	  //
	  //pComm = DefaultComm<int>::getDefaultSerialComm (null);
	  //
	  //RCP<const Comm<int> > nullComm; // is null.
	  //pComm = DefaultComm<int>::getDefaultSerialComm (nullComm);

	  pComm = rcp_implicit_cast<const Comm<int> > (rcp (new SerialComm<int> ()));
	}
    }
#endif // HAVE_MPI

    const int numProcs = pComm->getSize();
    const int myRank = pComm->getRank();

    if (debug && myRank == 0)
      {
	cerr << "summarize (out, "
	     << "alwaysWriteLocal=" << alwaysWriteLocal 
	     << ", writeGlobalStats=" << writeGlobalStats 
	     << ", writeZeroTimers=" << writeZeroTimers 
	     << ", setOp=" << (setOp==Union ? "Union" : "Intersection") 
	     << ")" << endl;
      }

    // Collect and sort local timer data by timer names.
    timer_map_t localTimerData;
    collectLocalTimerData (localTimerData, counters());

    // In debug mode, print out local timer data on each process,
    // before possibly filtering out data with zero call counts.
    if (debug)
      {
	for (int p = 0; p < numProcs; ++p)
	  {
	    if (myRank == p)
	      {
		cerr << "Proc " << myRank << ": Local timer data:" << endl;
		for (timer_map_t::const_iterator it = localTimerData.begin(); 
		     it != localTimerData.end(); ++it)
		  cerr << "-- " << it->first << ", " << it->second.first 
		       << ", " << it->second.second << endl;
	      }
	    // Two barriers generally synchronize output, at least
	    // when debugging with multiple MPI processes on one node.
	    barrier (*pComm);
	    barrier (*pComm);
	  }
      }

    // Filter out zero data locally first.  This ensures that if we
    // are writing global stats, and if a timer name exists in the set
    // of global names, then that timer has a nonzero call count on at
    // least one MPI process.
    if (! writeZeroTimers)
      {
	filterZeroData (localTimerData);

	// In debug mode, print out local timer data on each process,
	// after possibly filtering out data with zero call counts.
	if (debug)
	  {
	    for (int p = 0; p < numProcs; ++p)
	      {
		if (myRank == p)
		  {
		    cerr << "Proc " << myRank << ": Local timer data, "
		      "after filtering zero call counts:" << endl;
		    for (timer_map_t::const_iterator it = localTimerData.begin(); 
			 it != localTimerData.end(); ++it)
		      cerr << "-- " << it->first << ", " << it->second.first 
			   << ", " << it->second.second << endl;
		  }
		// Two barriers generally synchronize output, at least
		// when debugging with multiple MPI processes on one node.
		barrier (*pComm);
		barrier (*pComm);
	      }
	  }
      }

    // Extract the set of local timer names.  The std::map keeps them
    // sorted alphabetically.
    Array<string> localTimerNames;
    localTimerNames.reserve (localTimerData.size());
    for (timer_map_t::const_iterator it = localTimerData.begin(); 
     	 it != localTimerData.end(); ++it)
      localTimerNames.push_back (it->first);

    if (debug)
      {
	for (int p = 0; p < numProcs; ++p)
	  {
	    if (myRank == p)
	      {
		cerr << "Proc " << myRank << ": Local timer names:" << endl;
		for (Array<string>::const_iterator it = localTimerNames.begin(); 
		     it != localTimerNames.end(); ++it)
		  cerr << "-- " << *it << endl;
	      }
	    barrier (*pComm);
	    barrier (*pComm);
	  }
      }

    // globalTimerData and globalTimerNames are only valid if
    // writeGlobalStats is true.
    Array<string> globalTimerNames;
    timer_map_t globalTimerData;

    // If writeGlobalStats is true (i.e., if we are computing global
    // stats), there may be some global timers that are not local
    // timers on the calling MPI process(es).  In that case, if
    // alwaysWriteLocal is true, then we need to fill in the "missing"
    // local timers.  That will ensure that both global and local
    // timer columns in the output table have the same number of rows.
    // If writeZeroTimers==false, we already filtered out local timers
    // with zero call counts above.  Thus, any inserted local timer
    // data won't get filtered out again.
    //
    // Inserting new local data with zero call counts will result in
    // local timers with zero call counts, which violates the
    // expectation of writeZeroTimers == false for the local call
    // counts.  However, writeZeroTimers == false will still do what
    // it says for the global call counts.
    if (writeGlobalStats)
      { 
	// This does the correct and inexpensive thing (just copies
	// the timer data) if numProcs == 1.  Otherwise, it initiates
	// a communication with \f$O(\log P)\f$ messages along the
	// critical path.
	mergeCounterNames (*pComm, localTimerNames, globalTimerNames, setOp);

	if (debug)
	  {
	    // Sanity check that all MPI procs have the name number of
	    // global timer names.
	    const timer_map_t::size_type myNumGlobalNames = globalTimerNames.size();
	    timer_map_t::size_type minNumGlobalNames = 0;
	    timer_map_t::size_type maxNumGlobalNames = 0;
	    reduceAll (*pComm, REDUCE_MIN, myNumGlobalNames, 
		       outArg (minNumGlobalNames));
	    reduceAll (*pComm, REDUCE_MAX, myNumGlobalNames, 
		       outArg (maxNumGlobalNames));
	    TEST_FOR_EXCEPTION(minNumGlobalNames != maxNumGlobalNames,
			       std::logic_error,
			       "Min # global timer names = " << minNumGlobalNames 
			       << " != max # global timer names = " << maxNumGlobalNames
			       << ".  Please report this bug to the Teuchos developers.");
	    TEST_FOR_EXCEPTION(myNumGlobalNames != minNumGlobalNames,
			       std::logic_error,
			       "My # global timer names = " << myNumGlobalNames 
			       << " != min # global timer names = " << minNumGlobalNames
			       << ".  Please report this bug to the Teuchos developers.");
	  }

	// mergeTimers() just merges timer names, not their actual
	// data.  Now we need to fill globalTimerData with this MPI
	// process' timer data for those timers in globalTimerNames.
	//
	// All processes need the full list of global timers, since
	// there may be some global timers that are not local timers.
	// That's why mergeCounterNames() has to be an all-reduce, not
	// just a reduction to Proc 0.
	//
	// If there are some global timers that are not local timers,
	// and if we want to print local timers, we insert a local
	// timer datum with zero elapsed time and zero call count into
	// localTimerData as well, for the reason mentioned above.
	//
	// Insertion optimization; if the iterator given to
	// map::insert points right before where we want to insert,
	// insertion is O(1).  globalTimerNames is sorted, so feeding
	// the iterator output of map::insert into the next
	// invocation's input should make the whole insertion O(N)
	// where N is the number of entries in globalTimerNames.
	timer_map_t::iterator globalMapIter = globalTimerData.begin();
	timer_map_t::iterator localMapIter;
	for (Array<string>::const_iterator it = globalTimerNames.begin(); 
	     it != globalTimerNames.end(); ++it)
	  {
	    const std::string& globalName = *it;
	    localMapIter = localTimerData.find (globalName);

	    if (localMapIter == localTimerData.end())
	      {
		if (alwaysWriteLocal)
		  {
		    // We really only need to do this on MPI Proc 0, which is
		    // the only process that currently may print local timers.
		    // However, we do it on all MPI processes, just in case
		    // someone later wants to modify this function to print
		    // out local timer data for some MPI process other than
		    // Rank 0.  This extra computation won't affect the cost
		    // along the critical path, for future computations in
		    // which Proc 0 participates.
		    localMapIter = 
		      localTimerData.insert (localMapIter, 
					     makeEmptyTimerDatum (globalName));

		    // Make sure the missing global name gets added to
		    // the list of local names.  We'll resort it below.
		    localTimerNames.push_back (globalName);
		  }
		globalMapIter = 
		  globalTimerData.insert (globalMapIter, 
					  makeEmptyTimerDatum (globalName));
	      }
	    else
	      globalMapIter = 
		globalTimerData.insert (globalMapIter, 
					make_pair (globalName, 
						   localMapIter->second));
	  }

	if (alwaysWriteLocal)
	  // Resort the list of local timer names, since we may have
	  // inserted "missing" names above.
	  std::sort (localTimerNames.begin(), localTimerNames.end());

	if (debug)
	  {
	    // Sanity check that all MPI procs have the name number of
	    // global timers.
	    const timer_map_t::size_type myNumGlobalTimers = globalTimerData.size();
	    timer_map_t::size_type minNumGlobalTimers = 0;
	    timer_map_t::size_type maxNumGlobalTimers = 0;
	    reduceAll (*pComm, REDUCE_MIN, myNumGlobalTimers, 
		       outArg (minNumGlobalTimers));
	    reduceAll (*pComm, REDUCE_MAX, myNumGlobalTimers, 
		       outArg (maxNumGlobalTimers));
	    TEST_FOR_EXCEPTION(minNumGlobalTimers != maxNumGlobalTimers,
			       std::logic_error,
			       "Min # global timers = " << minNumGlobalTimers 
			       << " != max # global timers = " << maxNumGlobalTimers
			       << ".  Please report this bug to the Teuchos developers.");
	    TEST_FOR_EXCEPTION(myNumGlobalTimers != minNumGlobalTimers,
			       std::logic_error,
			       "My # global timers = " << myNumGlobalTimers 
			       << " != min # global timers = " << minNumGlobalTimers
			       << ".  Please report this bug to the Teuchos developers.");
	  }
      } // if (writeGlobalStats)

    // Extract the timer names (global or local) into a single array
    // of strings, representing the column labels in the table.
    Array<string> timerNames;
    timerNames.reserve (globalTimerData.size());
    if (writeGlobalStats) 
      // Use global timer names as the column names.
      std::copy (globalTimerNames.begin(), globalTimerNames.end(), 
		 std::back_inserter (timerNames));
    else 
      // Use local timer names as the column labels.
      std::copy (localTimerNames.begin(), localTimerNames.end(), 
		 std::back_inserter (timerNames));

    if (debug)
      {
	for (int p = 0; p < numProcs; ++p)
	  {
	    if (myRank == p)
	      {
		cerr << "Proc " << myRank << ": Global timer names:" << endl;
		for (Array<std::string>::const_iterator it = globalTimerNames.begin(); 
		     it != globalTimerNames.end(); ++it)
		  cerr << "-- " << *it << endl;
	      }
	    barrier (*pComm);
	    barrier (*pComm);
	  }
	for (int p = 0; p < numProcs; ++p)
	  {
	    if (myRank == p)
	      {
		cerr << "Proc " << myRank << ": Global timer data:" << endl;
		for (timer_map_t::const_iterator it = globalTimerData.begin(); 
		     it != globalTimerData.end(); ++it)
		  cerr << "-- " << it->first << ", " << it->second.first 
		       << ", " << it->second.second << endl;
	      }
	    barrier (*pComm);
	    barrier (*pComm);
	  }
      }

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
    // only write local stats if asked, only on MPI Proc 0, and only
    // if there is more than one MPI process in the communicator
    // (otherwise local stats == global stats, so we just print the
    // global stats).
    if (alwaysWriteLocal && numProcs > 1 && myRank == 0)
      {
	titles.append ("Local time (num calls)");

	// Copy local timer data out of the array-of-structs into
	// separate arrays, for display in the table.
	Array<double> localTimings;
	Array<double> localNumCalls;
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
	for (timer_map_t::const_iterator it = globalTimerData.begin();
	     it != globalTimerData.end(); ++it)
	  {
	    globalTimings.push_back (it->second.first);
	    globalNumCalls.push_back (static_cast<double> (it->second.second));
	  }

	if (numProcs == 1)
	  { 
	    // Don't display statistics in the case of only 1 MPI process.
	    // Just display the elapsed times and the call counts.
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
	    //
	    // NOTE (mfh 18 Jul 2011) The report minimum global number
	    // of calls may be for a different MPI process than the
	    // minimum global timing.  Ditto for the other statistics.
	    // What this means is that you should not divide the
	    // reported (minimum, mean, maximum) elapsed time by the
	    // reported (minimum, mean, maximum) call count to get an
	    // "average," since those quantities are not necessarily
	    // comparable.
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
      {
	std::ostringstream theTitle;
	theTitle << "TimeMonitor results over " << numProcs << " processor" 
		 << (numProcs > 1 ? "s" : "");
	format().writeWholeTable (out, theTitle.str(), titles, tableColumns);
      }
  }


} // namespace Teuchos
