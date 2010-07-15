#ifndef __TSQR_Trilinos_TsqrAdaptor_Epetra_MultiVector_hpp
#define __TSQR_Trilinos_TsqrAdaptor_Epetra_MultiVector_hpp

/// \file TsqrAdaptor_Epetra_MultiVector.hpp
///
/// \warning Users should _not_ include this file directly.  Include
///   "TsqrAdaptor.hpp" instead.  If HAVE_ANASAZI_EPETRA is defined,
///   then this file will be included automatically.  If for some
///   reason you need to include this file directly, be sure to
///   include "TsqrAdaptor.hpp" first.

#include "TsqrTypeAdaptor_Epetra_Multivector.hpp"
#include "Epetra_MultiVector.h" // sic (not .hpp)
#include "TsqrCommFactory.hpp"
#include <stdexcept>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    class TsqrEpetraAdaptor : 
      public TsqrAdaptor< double, int, int, Epetra_MultiVector >
    {
    public:
      // C++ doesn't inherit typedefs from base to derived class, when
      // both the base and the derived classes are templated.  To
      // avoid code duplication, we pull in typedefs from base_type.
      typedef TsqrAdaptor< double, int, int, Epetra_MultiVector > base_type;
      typedef base_type::comm_ptr comm_ptr;
      typedef base_type::multivector_type multivector_type;
      typedef base_type::scalar_type scalar_type;
      typedef base_type::local_ordinal_type local_ordinal_type;
      typedef typename base_type::scalar_messenger_ptr scalar_messenger_ptr;
      typedef typename base_type::ordinal_messenger_ptr ordinal_messenger_ptr;

      /// \brief Constructor
      ///
      /// \param mv [in] Multivector object, used only to access the
      ///   underlying map and its underlying communicator object (in
      ///   this case, Epetra_(Block)Map resp. Epetra_Comm).  All
      ///   multivector objects with which this Adaptor works must use
      ///   the same map and communicator.
      /// \param plist [in] List of parameters for configuring TSQR.
      ///   The specific parameter keys that are read depend on the
      ///   TSQR implementation.  "cacheBlockSize" (cache block size
      ///   per core, in bytes) tends to be defined for all of the
      ///   non-GPU implementations.  For details, check the specific
      ///   TsqrFactory implementation.
      TsqrEpetraAdaptor (const multivector_type& mv,
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

      Teuchos::ArrayView< scalar_type >::size_type 
      fetchArrayLength (const multivector_type& A) const
      {
	typedef Teuchos::ArrayView< scalar_type >::size_type size_type;

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
	  throw std::runtime_error ("TSQR::Trilinos::TsqrEpetraAdaptor: "
				    "Leading dimension of local part of "
				    "Epetra_MultiVector overflowed on co"
				    "nversion from local_ordinal_type to"
				    " ArrayView::size_type");
	const size_type ncols = static_cast< size_type > (ncols_LO);
	if (static_cast< local_ordinal_type > (ncols) != ncols_LO)
	  throw std::runtime_error ("TSQR::Trilinos::TsqrEpetraAdaptor: "
				    "Number of columns of the given "
				    "Epetra_MultiVector overflowed on co"
				    "nversion from local_ordinal_type to"
				    " ArrayView::size_type");
	// Make sure that the length of A.Values(), which is the
	// product of the leading dimension and the number of columns
	// of A, fits in a size_type.
	const size_type nelts = LDA * ncols;
	if (nelts / LDA != ncols)
	  throw std::runtime_error ("TSQR::Trilinos::TsqrEpetraAdaptor: "
				    " Length of A.Values() does not fit "
				    "in an ArrayView::size_type");
	return nelts;
      }

      virtual Teuchos::ArrayRCP< scalar_type > 
      fetchNonConstView (multivector_type& A) const
      {
	using Teuchos::arcpFromArrayView;
	using Teuchos::arrayView;
	typedef Teuchos::ArrayView< scalar_type >::size_type size_type;

	const size_type nelts = fetchArrayLength (A);
	// The returned ArrayRCP does NOT own A.Values().
	return arcpFromArrayView (arrayView (A.Values(), nelts));
      }

      virtual Teuchos::ArrayRCP< const scalar_type > 
      fetchConstView (const multivector_type& A) const
      {
	using Teuchos::arcpFromArrayView;
	using Teuchos::arrayView;
	using Teuchos::ArrayView;
	typedef ArrayView< scalar_type >::size_type size_type;

	const size_type nelts = fetchArrayLength (A);
	const scalar_type* A_ptr = A.Values();
	ArrayView< const scalar_type > A_view = arrayView (A_ptr, nelts);

	// The returned ArrayRCP does NOT own A.Values().
	return arcpFromArrayView (A_view);
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

#endif // __TSQR_Trilinos_TsqrAdaptor_Epetra_MultiVector_hpp
