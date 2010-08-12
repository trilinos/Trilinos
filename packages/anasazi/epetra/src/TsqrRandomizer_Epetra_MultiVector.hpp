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

#ifndef __TSQR_Trilinos_Randomizer_Epetra_MultiVector_hpp
#define __TSQR_Trilinos_Randomizer_Epetra_MultiVector_hpp

/// \file TsqrRandomizer_Epetra_MultiVector.hpp
///
/// \warning Anasazi users should _not_ include this file directly.

#include "TsqrTypeAdaptor_Epetra_Multivector.hpp"
#include "TsqrCommFactory_Epetra.hpp"
#include "Teuchos_ArrayView.hpp"
#include <limits>
#include <stdexcept>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    template< class Gen >
    class EpetraRandomizer :
      public Randomizer< double, int, int, Epetra_MultiVector, Gen >
    {
    public:
      // The C++ compiler needs help inheriting typedefs, when both
      // the base class and the derived class are templated.
      typedef Randomizer< double, int, int, Epetra_MultiVector, Gen > base_type;
      typedef typename base_type::multivector_type multivector_type;
      typedef typename base_type::local_ordinal_type local_ordinal_type;
      typedef typename base_type::scalar_type scalar_type;

      typedef typename base_type::normalgen_ptr normalgen_ptr;
      typedef typename base_type::ordinal_messenger_ptr ordinal_messenger_ptr;
      typedef typename base_type::scalar_messenger_ptr scalar_messenger_ptr;

      /// \brief Constructor
      ///
      /// \param mv [in] Multivector object, used only to access the
      ///   underlying map and its underlying communicator object (in
      ///   this case, Epetra_(Block)Map resp. Epetra_Comm).  All
      ///   multivector objects with which EpetraRandomizer works must
      ///   use the same map and communicator.
      /// \param pGen [in/out] Pointer to generator of pseudorandom
      ///   normal(0,1) sequence.
      EpetraRandomizer (const multivector_type& mv,
			const normalgen_ptr& pGen)
      {
	init (mv, pGen);
      }

    private:
      // mfh 14 Jul 2010: For who knows what reason, this refuses to
      // compile if you replace "double" with "scalar_type".  And no,
      // adding "typename" doesn't work either.
      typedef Teuchos::ArrayView< double >::size_type size_type;

      virtual void
      fetchDims (const multivector_type& A,
		 local_ordinal_type& nrowsLocal,
		 local_ordinal_type& ncols,
		 local_ordinal_type& LDA) const
      {
	nrowsLocal = A.MyLength();
	ncols = A.NumVectors();
	if (nrowsLocal < ncols)
	  {
	    std::ostringstream os;
	    os << "TSQR: The local component of the input matrix has fewer row"
	      "s (" << nrowsLocal << ") than columns (" << ncols << ").  TSQR "
	      "does not support this case.";
	    throw std::runtime_error (os.str());
	  }
	if (! A.ConstantStride())
	  {
	    std::ostringstream os;
	    os << "TSQR does not support Epetra_MultiVector inputs that do not"
	      " have constant stride.";
	    throw std::runtime_error (os.str());
	  }
	LDA = A.Stride();
      }

      size_type
      fetchArrayLength (const multivector_type& A) const
      {
	//typedef Teuchos::ArrayView< scalar_type >::size_type size_type;

	// Compute length (nelts) of the A.Values() array.  This only
	// makes sense if A has constant stride.  We convert from int
	// (the local_ordinal_type) to size_type, hoping that the
	// latter is no smaller than local_ordinal_type.  This gives
	// us a better chance that LDA*ncols fits in a size_type.  We
	// check for overflow both when converting each of LDA and
	// ncols from local_ordinal_type to size_type, and when
	// computing LDA*ncols (in size_type arithmetic).

	// Fetch dimensions of local part of A, as local_ordinal_type.
	local_ordinal_type nrowsLocal_LO, ncols_LO, LDA_LO;
	fetchDims (A, nrowsLocal_LO, ncols_LO, LDA_LO);

	// Convert dimensions from local_ordinal_type to size_type,
	// ensuring that the conversions don't overflow.
	const size_type LDA = static_cast< size_type > (LDA_LO);
	if (static_cast< local_ordinal_type > (LDA) != LDA_LO)
	  throw std::runtime_error ("TSQR::Trilinos::EpetraRandomizer: "
				    "Leading dimension of local part of "
				    "Epetra_MultiVector overflowed on co"
				    "nversion from local_ordinal_type to"
				    " ArrayView::size_type");
	const size_type ncols = static_cast< size_type > (ncols_LO);
	if (static_cast< local_ordinal_type > (ncols) != ncols_LO)
	  throw std::runtime_error ("TSQR::Trilinos::EpetraRandomizer: "
				    "Number of columns of the given "
				    "Epetra_MultiVector overflowed on co"
				    "nversion from local_ordinal_type to"
				    " ArrayView::size_type");
	// Make sure that the length of A.Values(), which is the
	// product of the leading dimension and the number of columns
	// of A, fits in a size_type.
	const size_type nelts = LDA * ncols;
	if (nelts / LDA != ncols)
	  throw std::runtime_error ("TSQR::Trilinos::EpetraRandomizer: "
				    " Length of A.Values() does not fit "
				    "in an ArrayView::size_type");
	return nelts;
      }

      virtual Teuchos::ArrayRCP< scalar_type > 
      fetchNonConstView (multivector_type& A) const
      {
	using Teuchos::arcpFromArrayView;
	using Teuchos::arrayView;

	const size_type nelts = fetchArrayLength (A);
	// The returned ArrayRCP does NOT own A.Values().
	return arcpFromArrayView (arrayView (A.Values(), nelts));
      }

      virtual void
      fetchMessengers (const multivector_type& mv,
		       scalar_messenger_ptr& pScalarMessenger,
		       ordinal_messenger_ptr& pOrdinalMessenger) const
      {
	typedef EpetraCommFactory comm_factory_type;
	// mv.Comm() returns a "const Epetra_Comm&".  rcpFromRef()
	// makes a non-owning (weak) RCP out of it.  This requires
	// that the mv's Epetra_Comm object not fall out of scope,
	// which it shouldn't as long as we are computing with
	// multivectors distributed according to that communicator.
	comm_ptr pComm = Teuchos::rcpFromRef (mv.Comm());
	comm_factory_type factory;
	factory.makeMessengers (pComm, pScalarMessenger, pOrdinalMessenger);
      }
    };

  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_Randomizer_Epetra_MultiVector_hpp
