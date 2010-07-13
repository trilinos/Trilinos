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
#include <limits>
#include <stdexcept>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    template< class Gen >
    class Randomizer< double, int, int, Epetra_MultiVector, Gen > {
    public:
      typedef double             scalar_type;
      typedef int                local_ordinal_type;
      typedef int                global_ordinal_type;
      typedef Epetra_MultiVector multivector_type;
      typedef Gen                normalgen_type;

      typedef double magnitude_type;
      typedef Teuchos::RCP< Gen > normalgen_ptr;
      typedef TSQR::Random::MatrixGenerator< double, int, Gen > matgen_type;
      typedef Teuchos::RCP< MessengerBase< int > > ordinal_messenger_ptr;
      typedef Teuchos::RCP< MessengerBase< double > > scalar_messenger_ptr;

      /// \brief Constructor
      ///
      /// TsqrRandomizer constructor.
      Randomizer (const normalgen_ptr& pGen,
		  const ordinal_messenger_ptr& pOrdinalMess,
		  const scalar_messenger_ptr& pScalarMess) : 
	pGen_ (pGen),
	pOrdinalMess_ (pOrdinalMess),
	pScalarMess_ (pScalarMess)
      {}

      /// \brief Fill A with a random (distributed) matrix
      ///
      /// Fill the MultiVector A with a random (distributed) matrix
      /// with the given singular values.  Given the same singular
      /// values and the same pseudorandom number sequence (produced
      /// by the Gen object), this function will always produce the
      /// same matrix, no matter the number of processors.  It
      /// achieves this at the cost of scalability; only Proc 0
      /// invokes the pseudorandom number generator.
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
      /// \brief Return dimensions of a multivector object
      ///
      /// For a given multivector A, return the number of rows stored
      /// locally on this process, the number of columns (multivectors
      /// are stored in a block row layout, so all columns of this
      /// process' row block are stored on this process), and the
      /// leading dimension of this process' row block (>= # rows on
      /// this process).
      ///
      /// \param A [in] The multivector object
      /// \param nrowsLocal [out] Number of rows of A stored locally 
      ///   on this process
      /// \param ncols [out] Number of columns of A
      /// \param LDA [out] Leading dimension of this process' row 
      ///   block of A
      void 
      fetchDims (const multivector_type& A, 
		 local_ordinal_type& nrowsLocal, 
		 local_ordinal_type& ncols, 
		 local_ordinal_type& LDA)
      {
	nrowsLocal = A.MyLength();
	ncols = A.NumVectors();
	if (nrowsLocal < ncols)
	  {
	    std::ostringstream os;
	    os << "The local component of the input matrix has fewer row"
	      "s (" << nrowsLocal << ") than columns (" << ncols << ").  "
	      "TSQR::Trilinos::Randomizer does not support this case.";
	    throw std::runtime_error (os.str());
	  }
	if (! A.ConstantStride())
	  {
	    std::ostringstream os;
	    os << "TSQR::Trilinos::Randomizer does not support Epetra_MultiVec"
	      "tor inputs that do not have constant stride.";
	    throw std::runtime_error (os.str());
	  }
	LDA = A.Stride();
      }

      Teuchos::ArrayRCP< scalar_type > 
      pEntries (multivector_type& A)
      {
	using Teuchos::arcpFromArrayView;
	using Teuchos::arrayView;
	typedef ArrayView< scalar_type >::size_type size_type;

	//
	// Compute length (nelts) of the A.Values() array.  This only
	// makes sense if A has constant stride.  We convert from int
	// (the local_ordinal_type) to size_type, hoping that the
	// latter is no smaller than local_ordinal_type.  This gives
	// us a better chance that LDA*ncols fits in a size_type.  We
	// check for overflow both when converting each of LDA and
	// ncols from local_ordinal_type to size_type, and when
	// computing LDA*ncols (in size_type arithmetic).
	//
	if (! A.ConstantStride())
	  {
	    std::ostringstream os;
	    os << "TSQR::Trilinos::Randomizer does not support Epetra_MultiVec"
	      "tor inputs that do not have constant stride.";
	    throw std::runtime_error (os.str());
	  }
	// Stride of A.  We assume column-major order and LDA >= # rows.
	const size_type LDA = static_cast< size_type > (A.Stride());
	// Number of columns of A
	const size_type ncols = static_cast< size_type > (A.NumVectors());
	// Make sure that the above conversions didn't overflow.
	if (static_cast< local_ordinal_type > (LDA) != A.Stride() ||
	    static_cast< local_ordinal_type > (ncols) != A.NumVectors())
	  throw std::runtime_error ("Overflow of local_ordinal_type into "
				    "ArrayView::size_type");
	// Make sure that the product of the above conversions doesn't
	// overflow.
	const size_type nelts = LDA * ncols;
	if (nelts / LDA != ncols)
	  throw std::runtime_error ("Overflow of LDA*ncols in "
				    "ArrayView::size_type");
	// The returned ArrayRCP does NOT own A.Values().
	return arcpFromArrayView (arrayView (A.Values(), nelts));
      }

      normalgen_ptr pGen_;
      ordinal_messenger_ptr pOrdinalMess_;
      scalar_messenger_ptr pScalarMess_;
    };

  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_Randomizer_Epetra_MultiVector_hpp
