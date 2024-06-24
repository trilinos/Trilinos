// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Anasazi_StubTsqrAdapter_hpp
#define __Anasazi_StubTsqrAdapter_hpp

#include <AnasaziConfigDefs.hpp>
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <stdexcept>

/// \file AnasaziStubTsqrAdapter.hpp
/// \brief "Stub" TSQR adapter for unsupported multivector types.

namespace Anasazi {
namespace details {

  /// \class StubTsqrAdapter
  /// \brief "Stub" TSQR adaptor for unsupported multivector types.
  ///
  /// TSQR (Tall Skinny QR factorization) is an orthogonalization
  /// kernel that is as accurate as Householder QR, yet requires only
  /// \f$2 \log P\f$ messages between $P$ MPI processes, independently
  /// of the number of columns in the multivector.  
  ///
  /// TSQR works independently of the particular multivector
  /// implementation, and interfaces to the latter via an adapter
  /// class.  Each multivector type MV needs its own adapter class.
  /// The specialization of MultiVecTraits for MV refers to its
  /// corresponding adapter class as its \c tsqr_adaptor_type [sic;
  /// sorry about the lack of standard spelling of "adapter"] typedef.
  /// For examples, please refer to the Epetra_MultiVector and
  /// Tpetra::MultiVector specializations of Anasazi::MultiVecTraits.
  ///
  /// Nevertheless, there may be multivector types for which a TSQR
  /// adapter has not yet been written.  This "stub" adapter
  /// implements the interface that TSQR adapters must implement, but
  /// all of its methods throw std::logic_error to indicate that this
  /// is a stub.  Thus, it allows Anasazi classes like
  /// TsqrOrthoManagerImpl to compile successfully for unsupported MV
  /// types.  This in turn allows OrthoManagerFactory to be templated
  /// on the MV type.
  template<class MultiVectorType>
  class StubTsqrAdapter : public Teuchos::ParameterListAcceptorDefaultBase {
  public:
    typedef MultiVectorType MV;
    typedef double scalar_type; // This doesn't really matter
    typedef int ordinal_type; // This doesn't matter either
    typedef int node_type; // Nor does this
    typedef Teuchos::SerialDenseMatrix<ordinal_type, scalar_type> dense_matrix_type;
    typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

    /// \brief Constructor (that accepts a parameter list).
    ///
    /// \param plist [in] List of parameters for configuring TSQR.
    ///   The specific parameter keys that are read depend on the TSQR
    ///   implementation.  For details, call \c getValidParameters()
    ///   and examine the documentation embedded therein.
    StubTsqrAdapter (const Teuchos::RCP<Teuchos::ParameterList>& plist) 
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "TSQR adapter for "
        "multivector type \"" << Teuchos::TypeNameTraits<MV>::name() 
        << " is not yet implemented.");
    }

    //! Default constructor (stub; throws std::logic_error).
    StubTsqrAdapter ()
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "TSQR adapter for "
        "multivector type \"" << Teuchos::TypeNameTraits<MV>::name() 
        << " is not yet implemented.");
    }

    //! Copy constructor (throws std::logic_error).
    StubTsqrAdapter (const StubTsqrAdapter& rhs)
    {
      (void) rhs;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "TSQR adapter for "
        "multivector type \"" << Teuchos::TypeNameTraits<MV>::name() 
        << " is not yet implemented.");
    }

    //! Get list of valid default parameters (stub; throws std::logic_error).
    Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters () const
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "TSQR adapter for "
        "multivector type \"" << Teuchos::TypeNameTraits<MV>::name() 
        << " is not yet implemented.");
    }

    //! Set parameters (stub; throws std::logic_error).
    void 
    setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "TSQR adapter for "
        "multivector type \"" << Teuchos::TypeNameTraits<MV>::name() 
        << " is not yet implemented.");
    }

    //! Compute QR factorization [Q,R] = qr(A,0) (stub; throws std::logic_error).
    void
    factorExplicit (MV& A,
		    MV& Q,
		    dense_matrix_type& R,
		    const bool forceNonnegativeDiagonal=false)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "TSQR adapter for "
        "multivector type \"" << Teuchos::TypeNameTraits<MV>::name() 
        << " is not yet implemented.");
    }

    //! Rank-revealing decomposition (stub; does nothing).
    int
    revealRank (MV& Q,
		dense_matrix_type& R,
		const magnitude_type& tol)
    {
      // mfh 08 Sep 2012: In order to prevent compiler warnings on
      // some platforms, we simply return some value.  This code can
      // never execute anyway, since it is in an instance method and
      // all of the constructors throw exceptions.  (We've overridden
      // the default and copy constructors to throw exceptions.)
      return 0; 
    }
  };

} // namespace details
} // namespace Anasazi

#endif // __Anasazi_StubTsqrAdapter_hpp
