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

#ifndef __TSQR_Trilinos_TsqrAdaptor_Tpetra_MultiVector_SerialNode_hpp
#define __TSQR_Trilinos_TsqrAdaptor_Tpetra_MultiVector_SerialNode_hpp

/// \file TsqrAdaptor_Tpetra_MultiVector.hpp
///
/// \warning Anasazi users should _not_ include this file directly.
///   Include "AnasaziTpetraAdaptor.hpp" instead.

#include "TsqrTypeAdaptor_Tpetra_MultiVector_SerialNode.hpp"
#ifdef HAVE_KOKKOS_TBB
#  include "TsqrTypeAdaptor_Tpetra_MultiVector_TBBNode.hpp"
#endif // HAVE_KOKKOS_TBB
#include "TsqrCommFactory_Tpetra.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    template< class S, class LO, class GO, class Node >
    class TsqrTpetraAdaptor :
      public TsqrAdaptor< S, LO, GO, Tpetra::MultiVector< S, LO, GO, Node > >
    {
    public:
      typedef Node node_type;
      // The MV typedef is just shorthand to make the base_type
      // typedef fit on one line.  It's not a circular definition,
      // because C++ doesn't inherit typedefs from base to derived
      // class, when both the base and the derived classes are
      // templated.
      typedef Tpetra::MultiVector< S, LO, GO, Node > MV;
      typedef TsqrAdaptor< S, LO, GO, MV > base_type;
      typedef typename base_type::comm_ptr comm_ptr;
      typedef typename base_type::multivector_type multivector_type;
      typedef typename base_type::scalar_type scalar_type;
      typedef typename base_type::local_ordinal_type local_ordinal_type;
      typedef typename base_type::scalar_messenger_ptr scalar_messenger_ptr;
      typedef typename base_type::ordinal_messenger_ptr ordinal_messenger_ptr;

      /// \brief Constructor
      ///
      /// \param mv [in] Multivector object, used only to access the
      ///   underlying communicator object (in this case,
      ///   Teuchos::Comm<int>, accessed via the Tpetra::Map belonging
      ///   to the multivector).  All multivector objects with which
      ///   this Adaptor works must use the same map and communicator.
      ///
      /// \param plist [in] List of parameters for configuring TSQR.
      ///   The specific parameter keys that are read depend on the
      ///   TSQR implementation.  "cacheBlockSize" (cache block size
      ///   per core, in bytes) tends to be defined for all of the
      ///   non-GPU implementations.  For details, check the specific
      ///   TsqrFactory implementation.
      TsqrTpetraAdaptor (const multivector_type& mv,
			 const Teuchos::ParameterList& plist)
      {
	init (mv, plist);
      }

    private:
      virtual void
      fetchDims (const multivector_type& A,
		 local_ordinal_type& nrowsLocal, 
		 local_ordinal_type& ncols, 
		 local_ordinal_type& LDA) const
      {
	nrowsLocal = A.getLocalLength();
	ncols = A.getNumVectors();
	if (nrowsLocal < ncols)
	  {
	    std::ostringstream os;
	    os << "TSQR: The local component of the input matrix has fewer row"
	      "s (" << nrowsLocal << ") than columns (" << ncols << ").  TSQR "
	      "does not support this case.";
	    throw std::runtime_error (os.str());
	  }
	if (! A.isConstantStride())
	  {
	    // FIXME (mfh 14 June 2010) Storage of A uses nonconstant
	    // stride internally, but that doesn't necessarily mean we
	    // can't run TSQR.  It depends on what get1dViewNonConst()
	    // returns.  If it's copied and packed into a matrix with
	    // constant stride, then we are free to run TSQR.
	    std::ostringstream os;
	    os << "TSQR does not support Tpetra::MultiVector inputs that do no"
	      "t have constant stride.";
	    throw std::runtime_error (os.str());
	  }
	LDA = A.getStride();
      }

      virtual Teuchos::ArrayRCP< scalar_type > 
      fetchNonConstView (multivector_type& A) const
      {
	return A.get1dViewNonConst();
      }

      virtual Teuchos::ArrayRCP< const scalar_type > 
      fetchConstView (const multivector_type& A) const
      {
	return A.get1dView();
      }

      virtual void
      fetchMessengers (const multivector_type& mv,
		       scalar_messenger_ptr& pScalarMessenger,
		       ordinal_messenger_ptr& pOrdinalMessenger) const
      {
	typedef TpetraCommFactory< S, LO, GO, Node > comm_factory_type;
	comm_ptr pComm = mv.getMap()->getComm();
	comm_factory_type factory;
	factory.makeMessengers (pComm, pScalarMessenger, pOrdinalMessenger);
      }
    };

  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrAdaptor_Tpetra_MultiVector_SerialNode_hpp
