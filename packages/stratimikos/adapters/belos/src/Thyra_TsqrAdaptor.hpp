// @HEADER
// ***********************************************************************
//
//                 Belos: Block Eigensolvers Package
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

#ifndef __Thyra_TsqrAdaptor_hpp
#define __Thyra_TsqrAdaptor_hpp

#include <BelosConfigDefs.hpp>

#include <Thyra_DetachedMultiVectorView.hpp>
#include <Thyra_MultiVectorBase.hpp>
#include <Thyra_MultiVectorStdOps.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace Thyra {

  /// \class TsqrAdaptor
  /// \brief Stub adaptor from Thyra::MultiVectorBase to TSQR
  ///
  /// TSQR (Tall Skinny QR factorization) is an orthogonalization
  /// kernel that is as accurate as Householder QR, yet requires only
  /// \f$2 \log P\f$ messages between $P$ MPI processes, independently
  /// of the number of columns in the multivector.  
  ///
  /// TSQR works independently of the particular multivector
  /// implementation, and interfaces to the latter via an adaptor
  /// class.  Thyra::TsqrAdaptor is the adaptor class for
  /// Thyra::MultiVectorBase.  It templates on the MultiVector (MV)
  /// type so that it can pick up that class' typedefs.  In
  /// particular, TSQR chooses its intranode implementation based on
  /// the Kokkos Node type of the multivector.
  ///
  /// \warning This is a stub adaptor that just placates the compiler
  ///   and does nothing.
  ///
  template< class Scalar >
  class TsqrAdaptor {
  public:
    typedef Thyra::MultiVectorBase< Scalar > MV;
    typedef Scalar scalar_type;
    typedef int ordinal_type; // MultiVectorBase really does use int for this
    typedef int node_type; // FIXME (mfh 26 Oct 2010) stub for now
    typedef Teuchos::SerialDenseMatrix< ordinal_type, scalar_type > dense_matrix_type;
    typedef typename Teuchos::ScalarTraits< scalar_type >::magnitudeType magnitude_type;

  public:
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
    ///   NodeTsqrFactory implementation.
    TsqrAdaptor (const MV& mv,
		 const Teuchos::ParameterList& plist) 
    {}

    /// \brief Compute QR factorization [Q,R] = qr(A,0)
    ///
    void
    factorExplicit (MV& A,
		    MV& Q,
		    dense_matrix_type& R)
    {}

    /// \brief Rank-revealing decomposition
    ///
    /// Using the R factor and explicit Q factor from
    /// factorExplicit(), compute the singular value decomposition
    /// (SVD) of R (\f$R = U \Sigma V^*\f$).  If R is full rank (with
    /// respect to the given relative tolerance tol), don't change Q
    /// or R.  Otherwise, compute \f$Q := Q \cdot U\f$ and \f$R :=
    /// \Sigma V^*\f$ in place (the latter may be no longer upper
    /// triangular).
    ///
    /// \param Q [in/out] On input: explicit Q factor computed by
    ///   factorExplicit().  (Must be an orthogonal resp. unitary
    ///   matrix.)  On output: If R is of full numerical rank with
    ///   respect to the tolerance tol, Q is unmodified.  Otherwise, Q
    ///   is updated so that the first rank columns of Q are a basis
    ///   for the column space of A (the original matrix whose QR
    ///   factorization was computed by factorExplicit()).  The
    ///   remaining columns of Q are a basis for the null space of A.
    ///
    /// \param R [in/out] On input: ncols by ncols upper triangular
    ///   matrix with leading dimension ldr >= ncols.  On output: if
    ///   input is full rank, R is unchanged on output.  Otherwise, if
    ///   \f$R = U \Sigma V^*\f$ is the SVD of R, on output R is
    ///   overwritten with $\Sigma \cdot V^*$.  This is also an ncols by
    ///   ncols matrix, but may not necessarily be upper triangular.
    ///
    /// \param tol [in] Relative tolerance for computing the numerical
    ///   rank of the matrix R.
    ///
    /// \return Rank \f$r\f$ of R: \f$ 0 \leq r \leq ncols\f$.
    ///
    int
    revealRank (MV& Q,
		dense_matrix_type& R,
		const magnitude_type& tol)
    {}
  };

} // namespace Tpetra

#endif // HAVE_KOKKOS_TSQR

#endif // __Thyra_TsqrAdaptor_hpp

