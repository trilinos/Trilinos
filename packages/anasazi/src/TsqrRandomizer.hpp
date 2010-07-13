#ifndef __TSQR_Trilinos_Randomizer_hpp
#define __TSQR_Trilinos_Randomizer_hpp

#include "AnasaziConfigDefs.hpp"
#include "Tsqr_ScalarTraits.hpp"
#include "Tsqr_Random_GlobalMatrix.hpp"
#include "TsqrTrilinosMessenger.hpp"
#include <string>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {
    /// \class Randomizer
    /// \brief Generates random test problems for TSQR
    ///
    /// Randomizer knows how to fill in an instance of the given
    /// MultiVector class MV with a (pseudo)random test problem,
    /// using a generator of type Gen.
    /// \li S type of the matrix entries
    /// \li LO local ordinal type
    /// \li GO global ordinal type
    /// \li MV MultiVector type
    /// \li Gen normal(0,1) pseudorandom generator type
    template< class S, class LO, class GO, class MV, class Gen >
    class Randomizer {
    public:
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

      /// \brief Constructor
      ///
      /// Constructor takes RCP smart pointers of the three objects
      /// needed to generate random test problems for TSQR.  RCPs of
      /// these objects are held by Randomizer for the lifetime of the
      /// latter.
      ///
      /// \param pGen [in/out] Pointer to generator of pseudorandom
      ///   normal(0,1) sequence.
      /// \param pOrdinalMess [in] Handles communication of 
      ///   local_ordinal_type numbers
      /// \param pScalarMess [in] Handles communication of 
      ///   scalar_type numbers
      Randomizer (const normalgen_ptr& pGen,
		  const ordinal_messenger_ptr& pOrdinalMess,
		  const scalar_messenger_ptr& pScalarMess) : 
	pGen_ (pGen),
	pOrdinalMess_ (pOrdinalMess),
	pScalarMess_ (pScalarMess)
      {}

      /// \brief Fill A with a (pseudo)random (distributed) matrix
      ///
      /// Fill the MultiVector A with a (pseudo)random (distributed)
      /// matrix with the given singular values.  Given the same
      /// singular values and the same pseudorandom number sequence
      /// (produced by the Gen object), this function will always
      /// produce the same matrix, no matter the number of processors.
      /// It achieves this at the cost of scalability; only Proc 0
      /// invokes the pseudorandom number generator.
      void
      randomMultiVector (multivector_type& A, 
			 const magnitude_type singularValues[])
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
		 local_ordinal_type& LDA);

      Teuchos::ArrayRCP< scalar_type > 
      pEntries (multivector_type& A);

      normalgen_ptr pGen_;
      ordinal_messenger_ptr pOrdinalMess_;
      scalar_messenger_ptr pScalarMess_;
    };

  } // namespace Trilinos
} // namespace TSQR

#ifdef HAVE_ANASAZI_EPETRA
#  include "TsqrRandomizer_Epetra_MultiVector.hpp"
#endif // HAVE_ANASAZI_EPETRA

#ifdef HAVE_ANASAZI_TPETRA
#  include "TsqrRandomizer_Tpetra_MultiVector.hpp"
#endif // HAVE_ANASAZI_TPETRA

#endif // __TSQR_Trilinos_Randomizer_hpp
