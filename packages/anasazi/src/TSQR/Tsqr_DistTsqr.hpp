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

#ifndef __TSQR_Tsqr_DistTsqr_hpp
#define __TSQR_Tsqr_DistTsqr_hpp

#include <Tsqr_ApplyType.hpp>
#include <Tsqr_DistTsqrHelper.hpp>
#include <Tsqr_ScalarTraits.hpp>

#include <Teuchos_RCP.hpp>

#include <algorithm> // std::min, std::max
#include <stdexcept>
#include <string>
#include <utility> // std::pair

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class DistTsqr
  /// \brief Internode part of TSQR
  ///
  /// This class combines the ncols by ncols R factors computed by the
  /// node TSQR factorization on individual MPI processes.
  template< class LocalOrdinal, class Scalar >
  class DistTsqr {
  public:
    typedef Scalar scalar_type;
    typedef LocalOrdinal ordinal_type;
    typedef std::vector< std::vector< Scalar > > VecVec;
    typedef std::pair< VecVec, VecVec > FactorOutput;
    typedef int rank_type;

    /// Constructor
    ///
    /// \param messenger [in/out] Encapsulation of communication
    ///   operations.  Not owned by *this (delete NOT called in
    ///   ~DistTsqr()).
    DistTsqr (const Teuchos::RCP< MessengerBase< Scalar > >& messenger) :
      messenger_ (messenger) 
    {}

    rank_type rank() const { return messenger_->rank(); }
    rank_type size() const { return messenger_->size(); }

    ///
    /// Destructor (does nothing)
    ///
    ~DistTsqr () {}

    /// Whether or not all diagonal entries of the R factor computed
    /// by the QR factorization are guaranteed to be nonnegative.
    bool QR_produces_R_factor_with_nonnegative_diagonal () const {
      return Combine< LocalOrdinal, Scalar >::QR_produces_R_factor_with_nonnegative_diagonal();
    }

    /// \brief Compute QR factorization of R factors, one per MPI process
    ///
    /// Compute the QR factorization of the P*ncols by ncols matrix
    /// consisting of all P nodes' R_mine upper triangular matrices
    /// stacked on top of each other.  Generally these upper triangular
    /// matrices should come from the QR factorization (perhaps computed
    /// by sequential or node-parallel TSQR) of a general matrix on each
    /// node.
    ///
    /// \note "ncols" below is the number of columns in the matrix to
    ///   factor.  Should be the same on all nodes.
    /// 
    /// \param R_mine [in,out] On input, an ncols by ncols upper triangular 
    ///   matrix with leading dimension ncols, stored unpacked (as a general 
    ///   matrix).  Elements below the diagonal are ignored.  On output, the
    ///   final R factor of the QR factorization of all nodes' different
    ///   R_mine inputs.  The final R factor is replicated over all nodes.
    ///
    /// \return Two arrays with the same number of elements: first, an
    ///   array of "local Q factors," and second, an array of "local tau
    ///   arrays."  These together form an implicit representation of
    ///   the Q factor.  They should be passed into the apply() and
    ///   explicit_Q() functions as the "factorOutput" parameter.
    FactorOutput
    factor (MatView< LocalOrdinal, Scalar > R_mine)
    {
      VecVec Q_factors, tau_arrays;
      DistTsqrHelper< LocalOrdinal, Scalar > helper;
      const LocalOrdinal ncols = R_mine.ncols();

      std::vector< Scalar > R_local (ncols*ncols);
      copy_matrix (ncols, ncols, &R_local[0], ncols, R_mine.get(), R_mine.lda());

      const int P = messenger_->size();
      const int my_rank = messenger_->rank();
      const int first_tag = 0;
      std::vector< Scalar > work (ncols);
      helper.factor_helper (ncols, R_local, my_rank, 0, P-1, first_tag, 
			    messenger_.get(), Q_factors, tau_arrays, work);
      copy_matrix (ncols, ncols, R_mine.get(), R_mine.lda(), &R_local[0], ncols);
      return std::make_pair (Q_factors, tau_arrays);
    }

    void
    apply (const ApplyType& apply_type,
	   const LocalOrdinal ncols_C,
	   const LocalOrdinal ncols_Q,
	   Scalar C_mine[],
	   const LocalOrdinal ldc_mine,
	   const FactorOutput& factor_output)
    {
      const bool transposed = apply_type.transposed();

      if (transposed)
	throw std::logic_error("DistTsqr: Applying Q^T or Q^H "
			       "not yet implemented");

      const int P = messenger_->size();
      const int my_rank = messenger_->rank();
      const int first_tag = 0;
      std::vector< Scalar > C_other (ncols_C * ncols_C);
      std::vector< Scalar > work (ncols_C);
  
      const VecVec& Q_factors = factor_output.first;
      const VecVec& tau_arrays = factor_output.second;

      // assert (Q_factors.size() == tau_arrays.size());
      const int cur_pos = Q_factors.size() - 1;
      DistTsqrHelper< ordinal_type, scalar_type > helper;
      helper.apply_helper (apply_type, ncols_C, ncols_Q, C_mine, ldc_mine, 
			   &C_other[0], my_rank, 0, P-1, first_tag, 
			   messenger_.get(), Q_factors, tau_arrays, cur_pos, 
			   work);
    }

    void
    explicit_Q (const LocalOrdinal ncols_Q,
		Scalar Q_mine[],
		const LocalOrdinal ldq_mine,
		const FactorOutput& factor_output)
    {
      const int my_rank = messenger_->rank ();
      fill_matrix (ncols_Q, ncols_Q, Q_mine, ldq_mine, Scalar(0));
      if (my_rank == 0)
	{
	  for (LocalOrdinal j = 0; j < ncols_Q; ++j)
	    Q_mine[j + j*ldq_mine] = Scalar (1);
	}
      apply (ApplyType::NoTranspose, ncols_Q, ncols_Q, 
	     Q_mine, ldq_mine, factor_output);
    }

  private:
    Teuchos::RCP< MessengerBase< Scalar > > messenger_;
  };

} // namespace TSQR

#endif // __TSQR_Tsqr_DistTsqr_hpp
