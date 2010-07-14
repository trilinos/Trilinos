#ifndef __TSQR_Trilinos_Randomizer_Epetra_MultiVector_hpp
#define __TSQR_Trilinos_Randomizer_Epetra_MultiVector_hpp

/// \file TsqrRandomizer_Epetra_MultiVector.hpp
///
/// \warning Users should _not_ include this file directly.  Include
///   "TsqrRandomizer.hpp" instead.  If HAVE_ANASAZI_EPETRA is
///   defined, then this file will be included automatically.  If for
///   some reason you need to include this file directly, be sure to
///   include "TsqrRandomizer.hpp" first.

#include "Epetra_MultiVector.h" // sic (not .hpp)
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
      EpetraRandomizer (const normalgen_ptr& pGen,
			const ordinal_messenger_ptr& pOrdinalMess,
			const scalar_messenger_ptr& pScalarMess) : 
	base_type (pGen, pOrdinalMess, pScalarMess) 
      {}

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
    };

  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_Randomizer_Epetra_MultiVector_hpp
