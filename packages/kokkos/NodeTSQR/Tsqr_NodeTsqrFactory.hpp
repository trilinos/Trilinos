// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2010) Sandia Corporation
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

#ifndef __TSQR_NodeTsqrFactory_hpp
#define __TSQR_NodeTsqrFactory_hpp

#include <Kokkos_ConfigDefs.hpp> // HAVE_KOKKOS_TBB

#ifdef HAVE_KOKKOS_TBB
#  include <Kokkos_TBBNode.hpp>
#  include <TbbTsqr.hpp>
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
  /// \brief Common interface for intranode TSQR
  template< class Node, class Scalar, class LocalOrdinal >
  class NodeTsqrFactory {
  public:
    typedef Node node_type;
    typedef Teuchos::RCP< node_type > node_ptr;
    // Just a default
    typedef SequentialTsqr< LocalOrdinal, Scalar > node_tsqr_type;
    
    static Teuchos::RCP< node_tsqr_type >
    makeNodeTsqr (const Teuchos::ParameterList& plist)
    {
      throw std::logic_error("TSQR is not supported on your Kokkos Node type");
    }
  };

  template< class T >
  static void
  getParamValue (const Teuchos::ParameterList& plist,
		 const char name[],
		 T& value, // set to default before calling
		 bool& gotValue)
  {
    // All this try/catch stuff is because the C++ compiler can't
    // deduce the right two-argument get() function (second argument
    // would be the default).
    try {
      const std::string paramName (name);
      value = plist.get< T > (paramName);
      gotValue = true;
    } catch (InvalidParameter&) {
      gotValue = false;
    }
  }

  static size_t
  getCacheBlockSize (const Teuchos::ParameterList& plist)
  {
    using Teuchos::Exceptions::InvalidParameter;
    const char* possibleNames[] = {"cacheBlockSize",
				   "cacheSize",
				   "cache_block_size",
				   "cache_size"};
    const int numPossibleNames = 4;
    size_t cacheBlockSize = 0;
    bool gotCacheBlockSize = false;

    for (int trial = 0; trial < numPossibleNames && ! gotCacheBlockSize; ++trial)
      getParamValue< size_t > (plist, possibleNames[trial], 
			       cacheBlockSize, gotCacheBlockSize);
    return cacheBlockSize;
  }

  static int
  getNumCores (const Teuchos::ParameterList& plist)
  {
    using Teuchos::Exceptions::InvalidParameter;
    const char* possibleNames[] = {"numCores",
				   "num_cores",
				   "ncores",
				   "numThreads",
				   "Num Threads",
				   "num_threads",
				   "nthreads"};
    const int numPossibleNames = 7;
    int numCores = 1;
    bool gotNumCores = false;

    for (int trial = 0; trial < numPossibleNames && ! gotNumCores; ++trial)
      getParamValue< int > (plist, possibleNames[trial], numCores, gotNumCores);
    return numCores;
  }

#ifdef HAVE_KOKKOS_TBB
  template< class Scalar, class LocalOrdinal >
  class NodeTsqrFactory< Kokkos::TBBNode, Scalar, LocalOrdinal > {
  public:
    typedef Kokkos::TBBNode node_type;
    typedef Teuchos::RCP< node_type > node_ptr;
    typedef TbbTsqr< LocalOrdinal, Scalar > node_tsqr_type;

    static Teuchos::RCP< node_tsqr_type >
    makeNodeTsqr (const Teuchos::ParameterList& plist)
    {
      const size_t cacheBlockSize = getCacheBlockSize (plist);
      const int numCores = getNumCores (plist);
      return Teuchos::rcp (new node_tsqr_type (numCores, cacheBlockSize));
    }
  };
#endif // HAVE_KOKKOS_TBB

  template< class Scalar, class LocalOrdinal >
  class NodeTsqrFactory< Kokkos::SerialNode, Scalar, LocalOrdinal > {
    typedef Scalar scalar_type;
    typedef LocalOrdinal ordinal_type;

    typedef Kokkos::SerialNode node_type;
    typedef Teuchos::RCP< node_type > node_ptr;
    typedef SequentialTsqr< LocalOrdinal, Scalar > node_tsqr_type;

    static Teuchos::RCP< node_tsqr_type >
    makeNodeTsqr (const Teuchos::ParameterList& plist)
    {
      const size_t cacheBlockSize = getCacheBlockSize (plist);
      return Teuchos::rcp (new node_tsqr_type (cacheBlockSize));
    }
  };

} // namespace TSQR

#endif // __TSQR_NodeTsqrFactory_hpp
