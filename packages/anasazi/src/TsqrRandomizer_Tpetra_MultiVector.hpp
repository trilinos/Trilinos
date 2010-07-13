#ifndef __TSQR_Trilinos_Randomizer_Tpetra_MultiVector_hpp
#define __TSQR_Trilinos_Randomizer_Tpetra_MultiVector_hpp

/// \file TsqrRandomizer_Tpetra_MultiVector.hpp
///
/// \warning Users should _not_ include this file directly.  Include
///   "TsqrRandomizer.hpp" instead.  If HAVE_ANASAZI_TPETRA is defined,
///   then this file will be included automatically.  If for some
///   reason you need to include this file directly, be sure to
///   include "TsqrRandomizer.hpp" first.

#include "Tpetra_MultiVector.hpp"
#ifdef HAVE_KOKKOS_TBB
#  include "Kokkos_TBBNode.hpp"
#endif // HAVE_KOKKOS_TBB

#include <stdexcept>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    template< class S, class LO, class GO, class NodeType, class Gen >
    class Randomizer< S, LO, GO, 
		      Tpetra::MultiVector< S, LO, GO, NodeType >,
		      Gen >
    {
    public:
      typedef Tpetra::MultiVector< S, LO, GO, NodeType > MV;

      typedef S   scalar_type;
      typedef LO  local_ordinal_type;
      typedef GO  global_ordinal_type;
      typedef MV  multivector_type;
      typedef Gen normalgen_type;

      typedef typename TSQR::ScalarTraits< S >::magnitude_type magnitude_type;
      typedef Teuchos::RCP< Gen > normalgen_ptr;
      typedef TSQR::Random::MatrixGenerator< S, LO, Gen > matgen_type;
      typedef Teuchos::RCP< MessengerBase< LO > > ordinal_messenger_ptr;
      typedef Teuchos::RCP< MessengerBase< S > > scalar_messenger_ptr;

      Randomizer (const normalgen_ptr& pGen,
		  const ordinal_messenger_ptr& pOrdinalMess,
		  const scalar_messenger_ptr& pScalarMess) : 
	pGen_ (pGen),
	pOrdinalMess_ (pOrdinalMess),
	pScalarMess_ (pScalarMess)
      {}

      void
      randomMultiVector (MV& A, const magnitude_type singularValues[])
      {
	using TSQR::Random::randomGlobalMatrix;
	using Teuchos::ArrayRCP;
	typedef MatView< local_ordinal_type, scalar_type > mat_view_type;

	local_ordinal_type nrowsLocal, ncols, LDA;
	fetchDims (A, nrowsLocal, ncols, LDA);
	ArrayRCP< scalar_type > A_ptr = pEntries (A);
	mat_view_type A_view (nrowsLocal, ncols, A_ptr.get(), LDA);

	randomGlobalMatrix (pGen_.get(), A_view, singularValues,
			    pOrdinalMess_.get(), pScalarMess_.get());
      }

    private:
      void
      fetchDims (const multivector_type& A,
		 local_ordinal_type& nrowsLocal, 
		 local_ordinal_type& ncols, 
		 local_ordinal_type& LDA)
      {
	nrowsLocal = A.getLocalLength();
	ncols = A.getNumVectors();
	if (nrowsLocal < ncols)
	  {
	    std::ostringstream os;
	    os << "The local component of the input matrix has fewer row"
	      "s (" << nrowsLocal << ") than columns (" << ncols << ").  "
	      "TSQR::Trilinos::Randomizer does not support this case.";
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
	    os << "TSQR::Trilinos::Randomizer does not support Tpetra::"
	      "MultiVector inputs that do not have constant stride.";
	    throw std::runtime_error (os.str());
	  }
	LDA = A.getStride();
      }

      Teuchos::ArrayRCP< scalar_type > 
      pEntries (multivector_type& A)
      {
	// This won't be efficient if the entries of A live in a
	// different memory space (e.g., GPU node), but it will be
	// correct for any Tpetra::MultiVector type.
	return A.get1dViewNonConst();
      }

      normalgen_ptr pGen_;
      ordinal_messenger_ptr pOrdinalMess_;
      scalar_messenger_ptr pScalarMess_;
    };
  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_Randomizer_Tpetra_MultiVector_hpp
