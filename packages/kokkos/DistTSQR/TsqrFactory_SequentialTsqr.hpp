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

#ifndef __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp
#define __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp

/// \file TsqrFactory_SequentialTsqr.hpp
/// \brief Declaration and definition of SequentialTsqrFactory.
///

#include "Tsqr_SequentialTsqr.hpp"
#include "Tsqr.hpp"
#include "Teuchos_ParameterListExceptions.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    /// \class SequentialTsqrFactory
    ///
    /// Subclass of \c TsqrFactory that knows how to instantiate
    /// \c SequentialTsqr as the intranode TSQR implementation.
    template<class LO, class S>
    class SequentialTsqrFactory :
      public TsqrFactory<LO, S, SequentialTsqr<LO, S>, DistTsqr<LO, S> >
    {
    public:
      // This class' parent class.
      typedef TsqrFactory<LO, S, SequentialTsqr<LO, S>, DistTsqr<LO, S> > base_type;

      // Pull in the typedefs from the base class.  C++ doesn't do
      // this when both the base and the derived classes are
      // templated.
      typedef typename base_type::scalar_messenger_ptr scalar_messenger_ptr;
      typedef typename base_type::node_tsqr_type node_tsqr_type;
      typedef typename base_type::node_tsqr_ptr  node_tsqr_ptr;
      typedef typename base_type::dist_tsqr_type dist_tsqr_type;
      typedef typename base_type::dist_tsqr_ptr  dist_tsqr_ptr;
      typedef typename base_type::tsqr_type tsqr_type;
      typedef typename base_type::tsqr_ptr  tsqr_ptr;

      SequentialTsqrFactory () {}
      virtual ~SequentialTsqrFactory () {}

    private:

      /// \brief Implementation of the parent class' pure virtual method.
      ///
      /// \param [in] plist Parameter list with the following key:
      ///   \li "cacheSizeHint": Cache size hint in bytes.  Default is
      ///       zero, which means that TSQR will guess a reasonable
      ///       cache size.  The size should correspond to that of the
      ///       largest cache that is private to each CPU core, if
      ///       such a private cache exists; alternately, it should
      ///       correspond to the amount of shared cache, divided by
      ///       the number of cores sharing that cache.
      virtual node_tsqr_ptr
      makeNodeTsqr (const Teuchos::ParameterList& plist) const
      {
	using Teuchos::Exceptions::InvalidParameter;
	size_t cacheBlockHint = 0;

	// All this try/catch stuff is because the C++ compiler can't
	// deduce the right two-argument get() function (second argument
	// would be the default).
	try {
	  const std::string cacheBlockHintParamName ("cacheBlockHint");
	  cacheBlockHint = plist.get< size_t > (cacheBlockHintParamName);
	} catch (InvalidParameter&) {
	  cacheBlockHint = 0;
	}
	
	node_tsqr_ptr node_tsqr (new node_tsqr_type (cacheBlockHint));
	return node_tsqr;
      }

      //! Implementation of the parent class' pure virtual method.
      virtual dist_tsqr_ptr
      makeDistTsqr (const scalar_messenger_ptr& messenger,
		    const Teuchos::ParameterList& plist) const
      {
	return Teuchos::rcp (new dist_tsqr_type (messenger));
      }
    };

  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrFactory_SequentialTsqr_hpp
