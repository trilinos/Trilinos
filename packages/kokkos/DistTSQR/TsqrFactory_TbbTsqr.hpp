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

#ifndef __TSQR_Trilinos_TsqrFactory_TbbTsqr_hpp
#define __TSQR_Trilinos_TsqrFactory_TbbTsqr_hpp

/// \file TsqrFactory_TbbTsqr.hpp
///
/// \warning Anasazi users should _not_ include this file directly.

#include "Kokkos_ConfigDefs.hpp" // HAVE_KOKKOS_TBB
#include "Teuchos_ParameterListExceptions.hpp"

#ifdef HAVE_KOKKOS_TBB
#  include "TbbTsqr.hpp"
#endif // HAVE_KOKKOS_TBB

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

#ifdef HAVE_KOKKOS_TBB
    template< class LO, class S >
    class TbbTsqrFactory :
      public TsqrFactory< LO, S, TSQR::TBB::TbbTsqr< LO, S >, DistTsqr< LO, S > >
    {
    public:
      // Help C++ pull in the typedefs from the base class.  C++ needs
      // help when both the base and the derived classes are
      // templated.
      typedef TsqrFactory< LO, S, TSQR::TBB::TbbTsqr< LO, S >, DistTsqr< LO, S > > base_type;
      typedef typename base_type::scalar_messenger_ptr scalar_messenger_ptr;
      typedef typename base_type::node_tsqr_type node_tsqr_type;
      typedef typename base_type::node_tsqr_ptr  node_tsqr_ptr;
      typedef typename base_type::dist_tsqr_type dist_tsqr_type;
      typedef typename base_type::dist_tsqr_ptr  dist_tsqr_ptr;
      typedef typename base_type::tsqr_type tsqr_type;
      typedef typename base_type::tsqr_ptr  tsqr_ptr;

      TbbTsqrFactory () {}
      virtual ~TbbTsqrFactory () {}

    private:
      /// \param [in] plist Parameter list with the following keys:
      ///   \li "cacheSizeHint": Cache size hint in bytes.  Default is
      ///       zero, which means that TSQR will guess a reasonable
      ///       cache size.  The size should correspond to that of the
      ///       largest cache that is private to each CPU core, if
      ///       such a private cache exists; alternately, it should
      ///       correspond to the amount of shared cache, divided by
      ///       the number of cores sharing that cache.
      ///   \li "numCores": If node_tsqr_type requires this (TbbTsqr
      ///       does, for Intel Threading Building Blocks intranode
      ///       parallelism), it should be set to the number of CPU
      ///       cores that are to participate in the intranode
      ///       parallel part of TSQR.  If node_tsqr_type does not
      ///       require this, then the parameter is ignored.
      virtual node_tsqr_ptr
      makeNodeTsqr (const Teuchos::ParameterList& plist) const
      {
	using Teuchos::Exceptions::InvalidParameter;
	int numCores = 1;
	size_t cacheSizeHint = 0;

	try {
	  const std::string numCoresParamName ("numCores");
	  numCores = plist.get<int>(numCoresParamName);
	} catch (InvalidParameter&) {
	  numCores = 1;
	}

	try {
	  const std::string cacheSizeHintParamName ("cacheSizeHint");
	  cacheSizeHint = plist.get<size_t>(cacheSizeHintParamName);
	} catch (InvalidParameter&) {
	  // Intranode TSQR interprets cacheSizeHint==0 as "set cache
	  // block size to a reasonable positive default value."
	  cacheSizeHint = 0;
	}

	node_tsqr_ptr nodeTsqr (new node_tsqr_type (numCores, cacheSizeHint));
	return nodeTsqr;
      }

      virtual dist_tsqr_ptr
      makeDistTsqr (const scalar_messenger_ptr& messenger,
		    const Teuchos::ParameterList& plist) const
      {
	return Teuchos::rcp (new dist_tsqr_type (messenger));
      }
    };
#endif // HAVE_KOKKOS_TBB

  } // namespace Trilinos
} // namespace TSQR


#endif // __TSQR_Trilinos_TsqrFactory_TbbTsqr_hpp
