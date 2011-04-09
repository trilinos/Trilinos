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

#ifndef __TSQR_Tsqr_NodeTsqr_hpp
#define __TSQR_Tsqr_NodeTsqr_hpp

#include <Tsqr_ApplyType.hpp>
#include <Tsqr_NodeTsqrBase.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  class NodeTsqrFactorOutput {
  public:
    NodeTsqrFactorOutput () {}
    virtual ~NodeTsqrFactorOutput() = 0;
  };

  /// \class NodeTsqr
  /// \brief Common functionality and interface for intranode TSQR.
  ///
  /// NodeTsqr provides a generic interface for TSQR operations within
  /// a node ("intranode"), as well as basic rank-revealing
  /// functionality.
  ///
  /// NodeTsqrFactorOutput is the FactorOutput type for the specific
  /// TSQR implementation.  When you inherit from this base class,
  /// fill in the specific NodeTsqrFactorOutput type.  It's awkward,
  /// but it gives us some flexibility in the interface.
  template<class LocalOrdinal, class Scalar>
  class NodeTsqr : public NodeTsqrBase<LocalOrdinal, Scalar> {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType magnitude_type;
    typedef LocalOrdinal ordinal_type;
    typedef NodeTsqrFactorOutput factor_output_type;

    //! Constructor
    NodeTsqr() {}

    //! Virtual destructor ensures safe polymorphic destruction.
    virtual ~NodeTsqr() {}

    /// \brief Reveal rank of TSQR's R factor.
    ///
    /// Compute SVD \f$R = U \Sigma V^*\f$, not in place.  Use the
    /// resulting singular values to compute the numerical rank of R,
    /// with respect to the relative tolerance tol.  If R is full
    /// rank, return without modifying R.  If R is not full rank,
    /// overwrite R with \f$\Sigma \cdot V^*\f$.
    ///
    /// \return Numerical rank of R: 0 <= rank <= ncols.
    virtual LocalOrdinal
    reveal_R_rank (const LocalOrdinal ncols,
		   Scalar R[],
		   const LocalOrdinal ldr,
		   Scalar U[],
		   const LocalOrdinal ldu,
		   const magnitude_type tol) const
    {
      if (tol < 0)
	{
	  std::ostringstream os;
	  os << "reveal_R_rank: negative tolerance tol = "
	     << tol << " is not allowed.";
	  throw std::logic_error (os.str());
	}
      else if (ncols < 0)
	{
	  std::ostringstream os;
	  os << "reveal_R_rank: number of columns is negative (= " << ncols << ")";
	  throw std::logic_error (os.str());
	}
      else if (ldr < ncols)
	{
	  std::ostringstream os;
	  os << "reveal_R_rank: stride of R (ldr = " << ldr 
	     << ") is less than the number of columns (ncols = " 
	     << ncols << ")";
	  throw std::logic_error (os.str());
	}
      else if (ldu < ncols)
	{
	  std::ostringstream os;
	  os << "reveal_R_rank: stride of U (ldu = " << ldu
	     << ") is less than the number of columns (ncols = " 
	     << ncols << ")";
	  throw std::logic_error (os.str());
	}

      // Take the easy exit if available.
      if (ncols == 0)
	return 0;

      LAPACK< LocalOrdinal, Scalar > lapack;
      MatView< LocalOrdinal, Scalar > R_view (ncols, ncols, R, ldr);
      Matrix< LocalOrdinal, Scalar > B (R_view); // B := R (deep copy)
      MatView< LocalOrdinal, Scalar > U_view (ncols, ncols, U, ldu);
      Matrix< LocalOrdinal, Scalar > VT (ncols, ncols, Scalar(0));

      std::vector< magnitude_type > svd_rwork (5*ncols);
      std::vector< magnitude_type > singular_values (ncols);
      LocalOrdinal svd_lwork = -1; // -1 for LWORK query; will be changed
      int svd_info = 0;

      // LAPACK LWORK query for singular value decomposition.  WORK
      // array is always of ScalarType, even in the complex case.
      {
	Scalar svd_lwork_scalar = Scalar(0);
	lapack.GESVD ("A", "A", ncols, ncols, 
		       B.get(), B.lda(), &singular_values[0], 
		       U_view.get(), U_view.lda(), VT.get(), VT.lda(),
		       &svd_lwork_scalar, svd_lwork, &svd_rwork[0], &svd_info);
	if (svd_info != 0)
	  {
	    std::ostringstream os;
	    os << "reveal_R_rank: GESVD LWORK query returned nonzero INFO = "
	       << svd_info;
	    throw std::logic_error (os.str());
	  }
	svd_lwork = static_cast< LocalOrdinal > (svd_lwork_scalar);
	// Check the LWORK cast.  LAPACK shouldn't ever return LWORK
	// that won't fit in an OrdinalType, but it's not bad to make
	// sure.
	if (static_cast< Scalar > (svd_lwork) != svd_lwork_scalar)
	  {
	    std::ostringstream os;
	    os << "reveal_R_rank: GESVD LWORK query returned LWORK that "
	      "doesn\'t fit in LocalOrdinal.  As a Scalar, LWORK = "
	       << svd_lwork_scalar << ", but cast to LocalOrdinal, "
	      "LWORK = " << svd_lwork << ".";
	    throw std::logic_error (os.str());
	  }
	// Make sure svd_lwork >= 0.
	if (std::numeric_limits< LocalOrdinal >::is_signed && svd_lwork < 0)
	  {
	    std::ostringstream os;
	    os << "reveal_R_rank: GESVD LWORK query returned negative "
	      "LWORK = " << svd_lwork;
	    throw std::logic_error (os.str());
	  }
      }
      // Allocate workspace for SVD.
      std::vector< Scalar > svd_work (svd_lwork);

      // Compute SVD $B := U \Sigma V^*$.  B is overwritten, which is
      // why we copied R into B (so that we don't overwrite R if R is
      // full rank).
      lapack.GESVD ("A", "A", ncols, ncols, 
		    B.get(), B.lda(), &singular_values[0], 
		    U_view.get(), U_view.lda(), VT.get(), VT.lda(),
		    &svd_work[0], svd_lwork, &svd_rwork[0], &svd_info);

      // GESVD computes singular values in decreasing order and they
      // are all nonnegative.  We know by now that ncols > 0.  "tol"
      // is a relative tolerance: relative to the largest singular
      // value, which is the 2-norm of the matrix.
      const magnitude_type absolute_tolerance = tol * singular_values[0];

      // Determine rank of B, using singular values.  
      LocalOrdinal rank = ncols;
      for (LocalOrdinal k = 1; k < ncols; ++k)
	// "<=" in case singular_values[0] == 0.
	if (singular_values[k] <= absolute_tolerance)
	  {
	    rank = k;
	    break;
	  }

      if (rank == ncols)
	return rank; // Don't modify Q or R, if R is full rank.

      //
      // R is not full rank.  
      //
      // 1. Compute \f$R := \Sigma V^*\f$.
      // 2. Return rank (0 <= rank < ncols).
      //

      // Compute R := \Sigma VT.  \Sigma is diagonal so we apply it
      // column by column (normally one would think of this as row by
      // row, but this "Hadamard product" formulation iterates more
      // efficiently over VT).  
      //
      // After this computation, R may no longer be upper triangular.
      // R may be zero if all the singular values are zero, but we
      // don't need to check for this case; it's rare in practice, and
      // the computations below will be correct regardless.
      for (LocalOrdinal j = 0; j < ncols; ++j)
	{
	  const Scalar* const VT_j = &VT(0,j);
	  Scalar* const R_j = &R_view(0,j);

	  for (LocalOrdinal i = 0; i < ncols; ++i)
	    R_j[i] = singular_values[i] * VT_j[i];
	}

      return rank;
    }

    /// \brief Compute rank-revealing decomposition.
    ///
    /// Using the R factor from factor() and the explicit Q factor
    /// from explicit_Q(), compute the SVD of R (\f$R = U \Sigma
    /// V^*\f$).  R.  If R is full rank (with respect to the given
    /// relative tolerance tol), don't change Q or R.  Otherwise,
    /// compute \f$Q := Q \cdot U\f$ and \f$R := \Sigma V^*\f$ in
    /// place (the latter may be no longer upper triangular).
    ///
    /// \return Rank \f$r\f$ of R: \f$ 0 \leq r \leq ncols\f$.
    ///
    virtual LocalOrdinal
    reveal_rank (const LocalOrdinal nrows,
		 const LocalOrdinal ncols,
		 Scalar Q[],
		 const LocalOrdinal ldq,
		 Scalar R[],
		 const LocalOrdinal ldr,
		 const magnitude_type tol,
		 const bool contiguousCacheBlocks) const
    {
      // Take the easy exit if available.
      if (ncols == 0)
	return 0;
      Matrix< LocalOrdinal, Scalar > U (ncols, ncols, Scalar(0));
      const LocalOrdinal rank = 
	reveal_R_rank (ncols, R, ldr, U.get(), U.ldu(), tol);
      
      if (rank < ncols)
	{
	  // If R is not full rank: reveal_R_rank() already computed
	  // the SVD \f$R = U \Sigma V^*\f$ of (the input) R, and
	  // overwrote R with \f$\Sigma V^*\f$.  Now, we compute \f$Q
	  // := Q \cdot U\f$, respecting cache blocks of Q.
	  Q_times_B (nrows, ncols, Q, ldq, U.get(), U.lda(), 
		     contiguousCacheBlocks);
	}
      return rank;
    }
  };
} // namespace TSQR


#endif __TSQR_Tsqr_NodeTsqr_hpp
