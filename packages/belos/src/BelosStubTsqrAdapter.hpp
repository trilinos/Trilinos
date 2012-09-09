// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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

#ifndef __Belos_StubTsqrAdapter_hpp
#define __Belos_StubTsqrAdapter_hpp

#include <BelosConfigDefs.hpp>
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <stdexcept>

/// \file BelosStubTsqrAdapter.hpp
/// \brief "Stub" TSQR adapter for unsupported multivector types.

namespace Belos {
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
  /// Tpetra::MultiVector specializations of Belos::MultiVecTraits.
  ///
  /// Nevertheless, there may be multivector types for which a TSQR
  /// adapter has not yet been written.  This "stub" adapter
  /// implements the interface that TSQR adapters must implement, but
  /// all of its methods throw std::logic_error to indicate that this
  /// is a stub.  Thus, it allows Belos classes like
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
      // mfh 07 Sep 2012: In order to prevent compiler warnings on
      // some platforms, we simply return some value.  This code can
      // never execute anyway, since it is in an instance method and
      // all of the constructors throw exceptions.  (We've overridden
      // the default and copy constructors to throw exceptions.)
      return 0; 
    }
  };

} // namespace details
} // namespace Belos

#endif // __Belos_StubTsqrAdapter_hpp

