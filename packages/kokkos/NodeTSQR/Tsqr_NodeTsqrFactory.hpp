//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef __TSQR_NodeTsqrFactory_hpp
#define __TSQR_NodeTsqrFactory_hpp

#include <Kokkos_ConfigDefs.hpp> // HAVE_KOKKOS_TBB

#ifdef HAVE_KOKKOS_TBB
#  include <Kokkos_TBBNode.hpp>
#  include <TbbTsqr.hpp>
#  include <tbb/task_scheduler_init.h>
#endif // HAVE_KOKKOS_TBB

#include <Kokkos_SerialNode.hpp>
#include <Tsqr_SequentialTsqr.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListExceptions.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class NodeTsqrFactory
  /// \brief Common interface for intranode TSQR.
  template<class Node, class Scalar, class LocalOrdinal>
  class NodeTsqrFactory {
  public:
    typedef Node node_type;
    typedef Teuchos::RCP<node_type> node_ptr;
    // Just a default
    typedef SequentialTsqr<LocalOrdinal, Scalar> node_tsqr_type;

    /// \brief Default parameter list for intranode TSQR.
    ///
    /// \note The default implementation returns an empty (not null)
    ///   parameter list.  Each specialization for a specific Node
    ///   type redefines this method to return a parameter list
    ///   appropriate for that Node type's TSQR implementation.
    static Teuchos::RCP<const Teuchos::ParameterList>
    getDefaultParameters ()
    {
      using Teuchos::ParameterList;
      using Teuchos::RCP;

      // FIXME (mfh 28 Oct 2010) The use of persistent data means that
      // this routine is NOT reentrant.  That means that it is not
      // thread-safe, for example.  One way to make this routine
      // reentrant would be to use a construct like pthread_once().
      static RCP< const ParameterList > defaultParams_;
      if (defaultParams_.is_null())
	{
	  RCP<ParameterList> params = Teuchos::parameterList();
	  defaultParams_ = params;
	}
      return defaultParams_;
    }

    /// \brief Return a pointer to the intranode TSQR implementation
    ///
    /// In a proper implementation, this would teturn a pointer to the
    /// intranode TSQR implementation.  This class method is
    /// implemented with correct behavior for those Kokkos Node types
    /// for which we have an intranode TSQR implementation.
    static Teuchos::RCP< node_tsqr_type >
    makeNodeTsqr (const Teuchos::RCP<const Teuchos::ParameterList>& plist)
    {
      throw std::logic_error("TSQR is not supported on your Kokkos Node type");
    }
  };

  template< class T >
  static void
  getParamValue (const Teuchos::RCP<const Teuchos::ParameterList>& plist,
		 const char name[],
		 T& value, // set to default before calling
		 bool& gotValue)
  {
    if (plist.is_null())
      {
	gotValue = false;
	return;
      }
    // All this try/catch stuff is because the C++ compiler can't
    // deduce the right two-argument get() function (second argument
    // would be the default).
    T retrievedValue;
    try {
      const std::string paramName (name);
      // We know from above that plist is not null.
      retrievedValue = plist->get< T > (paramName);
      gotValue = true;
    } catch (Teuchos::Exceptions::InvalidParameter&) { 
      // Could be wrong type (InvalidParameterType) or wrong name
      // (InvalidParameterName).  In either case, we just say that the
      // value doesn't exist.  
      //
      // For now, we just ignore the parameter if the type is wrong.
      // This is because we want TSQR "just to work" even if the
      // parameters are wrong; the parameters are "options" and should
      // be "optional."
      gotValue = false;
    }
    // Only write to the output argument if we got a value out of the
    // parameter list.  (This means that getParamValue() has the
    // strong exception guarantee (no destructive updates if it
    // throws), as long as T's assignment operator also has the strong
    // exception guarantee.)
    if (gotValue)
      value = retrievedValue;
  }

  static size_t
  getCacheSizeHint (const Teuchos::RCP<const Teuchos::ParameterList>& params, 
		    const Teuchos::RCP<const Teuchos::ParameterList>& defaultParams)
  {
    // We try to guess among some reasonable names.  The first in the
    // list is the canonical name.  "cacheBlockSize" and related names
    // are retained only for backwards compatibility and may be
    // removed at any time.
    const char* possibleNames[] = {"cacheSizeHint",
				   "cache_size_hint",
				   "cacheSize",
				   "cache_size",
				   "cacheBlockSize",
				   "cache_block_size"};
    const int numPossibleNames = 6;
    size_t cacheSizeHint = 0;
    bool gotCacheSizeHint = false;
      
    for (int trial = 0; trial < numPossibleNames && ! gotCacheSizeHint; ++trial)
      getParamValue< size_t > (params, possibleNames[trial], 
			       cacheSizeHint, gotCacheSizeHint);
    if (! gotCacheSizeHint)
      {
	// Default parameters had better have the value, so we don't
	// try to catch any exceptions here.
	const std::string canonicalName (possibleNames[0]);
	cacheSizeHint = defaultParams->get< size_t > (canonicalName);
      }
    return cacheSizeHint;
  }

#ifdef HAVE_KOKKOS_TBB
  static int
  getNumCores (const Teuchos::RCP<const Teuchos::ParameterList>& params,
	       const Teuchos::RCP<const Teuchos::ParameterList>& defaultParams)
  {
    // We try to guess among some reasonable names.
    // The first in the list is the canonical name.
    const char* possibleNames[] = {"numCores",
				   "ncores",
				   "numThreads",
				   "nthreads", 
				   "numTasks",
				   "ntasks"};
    const int numPossibleNames = 6;
    int numCores = 1; // will reset this below
    bool gotNumCores = false;

    for (int trial = 0; trial < numPossibleNames && ! gotNumCores; ++trial)
      getParamValue< int > (params, possibleNames[trial], 
			    numCores, gotNumCores);
    if (! gotNumCores)
      {
	// Default parameters had better have the value, so we don't
	// try to catch any exceptions here.
	const std::string canonicalName (possibleNames[0]);
	numCores = defaultParams->get< int > (canonicalName);
      }
    return numCores;
  }

  template<class Scalar, class LocalOrdinal>
  class NodeTsqrFactory<Kokkos::TBBNode, Scalar, LocalOrdinal> {
  public:
    typedef Kokkos::TBBNode node_type;
    typedef Teuchos::RCP<node_type> node_ptr;
    typedef TBB::TbbTsqr<LocalOrdinal, Scalar> node_tsqr_type;

    /// \brief Default parameters
    ///
    /// Return default parameters and their values for setting up the
    /// Intel Threading Building Blocks (TBB) version of intranode
    /// TSQR.
    ///
    /// Currently, the only parameters that matter for this
    /// implementation are "cacheSizeHint" (a size_t) and "numCores"
    /// (an int).  Not setting "cacheSizeHint" or setting it to zero
    /// will ask TSQR to pick a reasonable default.  Not setting
    /// "numCores" will ask TSQR to pick a reasonable default, namely,
    /// tbb::task_scheduler_init::default_num_threads().  Mild
    /// oversubscription will not decrease performance, but
    /// undersubscription may decrease performance.
    static Teuchos::RCP<const Teuchos::ParameterList>
    getDefaultParameters ()
    {
      using Teuchos::ParameterList;
      using Teuchos::RCP;

      // FIXME (mfh 28 Oct 2010) The use of persistent data means that
      // this routine is NOT reentrant.  That means that it is not
      // thread-safe, for example.  One way to make this routine
      // reentrant would be to use a construct like pthread_once().
      static RCP<const ParameterList> defaultParams_;
      if (defaultParams_.is_null())
	{
	  RCP<ParameterList> params = Teuchos::parameterList();
	  // The TSQR implementation sets a reasonable default value
	  // if you tell it that the cache block size is zero.
	  const size_t defaultCacheSizeHint = 0;
	  params->set ("cacheSizeHint", defaultCacheSizeHint,
		       "Cache size hint in bytes (as a size_t) to use for intr"
		       "anode TSQR.  If zero, the intranode TSQR implementation"
		       " will pick a reasonable default.");
	  // TSQR uses a recursive division of the tall skinny matrix
	  // into TBB tasks.  Each task works on a block row.  The TBB
	  // task scheduler ensures that oversubscribing TBB tasks
	  // won't oversubscribe cores, so it's OK if
	  // default_num_threads() is too many.  For example, TBB
	  // might say default_num_threads() is the number of cores on
	  // the node, but the TBB task scheduler might have been
	  // initialized with the number of cores per NUMA region, for
	  // hybrid MPI + TBB parallelism.
	  const int defaultNumCores = 
	    tbb::task_scheduler_init::default_num_threads();
	  params->set ("numCores", defaultNumCores, 
		       "Number of tasks (\"threads\") to use in the intranode "
		       "parallel part of TSQR.  There is little/no performance "
		       "penalty for mild oversubscription, but a potential "
		       "penalty for undersubscription.");
	  defaultParams_ = params;
	}
      return defaultParams_;
    }

    static Teuchos::RCP<node_tsqr_type>
    makeNodeTsqr (const Teuchos::RCP<const Teuchos::ParameterList>& params)
    {
      Teuchos::RCP<const Teuchos::ParameterList> defaultParams = 
	getDefaultParameters ();
      const size_t cacheSizeHint = getCacheSizeHint (params, defaultParams);
      const int numCores = getNumCores (params, defaultParams);
      return Teuchos::rcp (new node_tsqr_type (numCores, cacheSizeHint));
    }
  };
#endif // HAVE_KOKKOS_TBB

  template< class Scalar, class LocalOrdinal >
  class NodeTsqrFactory<Kokkos::SerialNode, Scalar, LocalOrdinal> {
  public:
    typedef Kokkos::SerialNode node_type;
    typedef Teuchos::RCP<node_type> node_ptr;
    typedef SequentialTsqr<LocalOrdinal, Scalar> node_tsqr_type;

    /// \brief Default parameters
    ///
    /// Return default parameters and their values for setting up the
    /// sequential implementation of the intranode part of TSQR.
    ///
    /// Currently, the only parameter that matters for this
    /// implementation is "cacheSizeHint" (a size_t).  Not setting
    /// that parameter or setting it to zero will pick a reasonable
    /// default.
    static Teuchos::RCP< const Teuchos::ParameterList >
    getDefaultParameters ()
    {
      using Teuchos::ParameterList;
      using Teuchos::RCP;

      // FIXME (mfh 28 Oct 2010) The use of persistent data means that
      // this routine is NOT reentrant.  That means that it is not
      // thread-safe, for example.  One way to make this routine
      // reentrant would be to use a construct like pthread_once().
      static RCP<const ParameterList> defaultParams_;
      if (defaultParams_.is_null())
	{
	  RCP<ParameterList> params = Teuchos::parameterList();
	  // The TSQR implementation sets a reasonable default value
	  // if you give it a cache size hint of zero.
	  const size_t defaultCacheSizeHint = 0;
	  params->set ("cacheSizeHint", defaultCacheSizeHint,
		       "Cache size hint in bytes (as a size_t) to use for intra"
		       "node TSQR.  If zero, the intranode TSQR implementation "
		       "will pick a reasonable default.");
	  defaultParams_ = params;
	}
      return defaultParams_;
    }

    static Teuchos::RCP<node_tsqr_type>
    makeNodeTsqr (const Teuchos::RCP<const Teuchos::ParameterList>& params)
    {
      Teuchos::RCP<const Teuchos::ParameterList> defaultParams = 
	getDefaultParameters ();
      const size_t cacheSizeHint = getCacheSizeHint (params, defaultParams);
      return Teuchos::rcp (new node_tsqr_type (cacheSizeHint));
    }
  };

} // namespace TSQR

#endif // __TSQR_NodeTsqrFactory_hpp
