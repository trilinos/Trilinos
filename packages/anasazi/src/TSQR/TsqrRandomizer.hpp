// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_Trilinos_Randomizer_hpp
#define __TSQR_Trilinos_Randomizer_hpp

#include "AnasaziConfigDefs.hpp"
#include "TsqrTypeAdaptor.hpp"
#include "TsqrCommFactory.hpp"

#include "Tsqr_ScalarTraits.hpp"
#include "Tsqr_Random_GlobalMatrix.hpp"

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
      typedef Teuchos::RCP< Gen > normalgen_ptr;
      typedef TSQR::Random::MatrixGenerator< S, LO, Gen > matgen_type;

      typedef typename TSQR::ScalarTraits< S >::magnitude_type magnitude_type;

      typedef TsqrTypeAdaptor< S, LO, GO, MV >      type_adaptor;
      typedef typename type_adaptor::comm_type      comm_type;
      typedef typename type_adaptor::comm_ptr       comm_ptr;
      typedef Teuchos::RCP< MessengerBase< LO > > ordinal_messenger_ptr;
      typedef Teuchos::RCP< MessengerBase< S > >  scalar_messenger_ptr;

      virtual ~Randomizer() {}

      /// \brief Fill A with a (pseudo)random (distributed) matrix
      ///
      /// Fill the MultiVector A with a (pseudo)random (distributed)
      /// matrix with the given singular values.  Given the same
      /// singular values and the same pseudorandom number sequence
      /// (produced by the Gen object), this function will always
      /// produce the same matrix, no matter the number of processors.
      /// It achieves this at the cost of scalability; only Proc 0
      /// invokes the pseudorandom number generator.
      virtual void
      randomMultiVector (multivector_type& A, 
			 const magnitude_type singularValues[])
      {
	using TSQR::Random::randomGlobalMatrix;
	using Teuchos::ArrayRCP;
	typedef MatView< local_ordinal_type, scalar_type > matview_type;

	local_ordinal_type nrowsLocal, ncols, LDA;
	fetchDims (A, nrowsLocal, ncols, LDA);
	ArrayRCP< scalar_type > A_ptr = fetchNonConstView (A);
	matview_type A_view (nrowsLocal, ncols, A_ptr.get(), LDA);

	randomGlobalMatrix (pGen_.get(), A_view, singularValues,
			    pOrdinalMessenger_.get(), pScalarMessenger_.get());
      }

    protected:
      /// Like the constructor, except you're not supposed to call the
      /// constructor of a pure virtual class.
      ///
      /// \param mv [in] Only used to extract the underlying
      ///   communication object (e.g., Epetra_Comm or
      ///   Teuchos::Comm<int>).
      /// \param pGen [in/out] Pointer to generator of pseudorandom
      ///   normal(0,1) sequence.
      void 
      init (const multivector_type& mv,
	    const normalgen_ptr& pGen)
      {
	pGen_ = pGen;
	// This is done in a multivector type - dependent way.
	fetchMessengers (mv, pScalarMessenger_, pOrdinalMessenger_);
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
      virtual void 
      fetchDims (const multivector_type& A, 
		 local_ordinal_type& nrowsLocal, 
		 local_ordinal_type& ncols, 
		 local_ordinal_type& LDA) const = 0;

      /// \return Non-const smart pointer to the node-local data in A
      ///
      /// \note Child classes should implement this in such a way as
      /// to make the above public methods always correct (though not
      /// necessarily efficient) for all multivector types.  (It may
      /// not be efficient if the ArrayRCP copies between different
      /// memory spaces.)
      virtual Teuchos::ArrayRCP< scalar_type > 
      fetchNonConstView (multivector_type& A) const = 0;

      /// Maps from multivector_type object to (scalar_messenger_ptr,
      /// ordinal_messenger_ptr).
      virtual void
      fetchMessengers (const multivector_type& mv,
		       scalar_messenger_ptr& pScalarMessenger,
		       ordinal_messenger_ptr& pOrdinalMessenger) const = 0;

      normalgen_ptr pGen_;
      ordinal_messenger_ptr pOrdinalMessenger_;
      scalar_messenger_ptr pScalarMessenger_;
    };

  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_Randomizer_hpp
