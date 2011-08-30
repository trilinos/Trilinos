//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
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
// ************************************************************************
//@HEADER

/// \file Tsqr_DistTsqr.hpp
/// \brief Internode part of TSQR.
///
#ifndef __TSQR_Tsqr_DistTsqr_hpp
#define __TSQR_Tsqr_DistTsqr_hpp

#include <Tsqr_DistTsqrHelper.hpp>
#include <Tsqr_DistTsqrRB.hpp>

#include <utility> // std::pair


namespace TSQR {

  /// \class DistTsqr
  /// \brief Internode part of TSQR.
  /// \author Mark Hoemmen
  ///
  /// \tparam LocalOrdinal Index type for dense matrices of Scalar. 
  /// \tparam Scalar Value type for matrices to factor. 
  ///
  /// This class combines the square R factors computed by the
  /// intranode TSQR factorization (\c NodeTsqr subclass) on
  /// individual MPI processes.
  ///
  /// It should be possible to instantiate
  /// DistTsqr<LocalOrdinal,Scalar> for any LocalOrdinal and Scalar
  /// types for which \c Combine<LocalOrdinal, Scalar> and \c
  /// LAPACK<LocalOrdinal, Scalar> can be instantiated.
  template<class LocalOrdinal, class Scalar>
  class DistTsqr {
  public:
    typedef Scalar scalar_type;
    typedef LocalOrdinal ordinal_type;
    typedef MatView<ordinal_type, scalar_type > matview_type;
    typedef std::vector<std::vector<scalar_type> > VecVec;
    typedef std::pair<VecVec, VecVec> FactorOutput;
    typedef int rank_type;

    /// \brief Constructor.
    ///
    /// \param messenger [in/out] Wrapper of communication operations
    ///   between MPI processes.
    DistTsqr (const Teuchos::RCP<MessengerBase<scalar_type> >& messenger) :
      messenger_ (messenger),
      reduceBroadcastImpl_ (messenger)
    {}

    /// \brief Rank of this (MPI) process.
    /// 
    /// Rank is computed via MPI_Comm_rank() on the underlying
    /// communicator, if the latter is an MPI communicator.  If it's a
    /// serial "communicator," the rank is always zero.
    rank_type rank() const { return messenger_->rank(); }

    /// \brief Total number of MPI processes in this communicator.
    ///
    /// The size is communicated via MPI_Comm_size() on the underlying
    /// communicator, if the latter is an MPI communicator.  If it's a
    /// serial "communicator," the size is always one.
    rank_type size() const { return messenger_->size(); }

    /// \brief Destructor.
    ///
    /// The destructor doesn't need to do anything, thanks to smart
    /// pointers.
    ~DistTsqr () {}

    /// \brief Does the R factor have a nonnegative diagonal?
    ///
    /// DistTsqr implements a QR factorization (of a distributed
    /// matrix with a special structure).  Some, but not all, QR
    /// factorizations produce an R factor whose diagonal may include
    /// negative entries.  This Boolean tells you whether DistTsqr
    /// promises to compute an R factor whose diagonal entries are all
    /// nonnegative.
    bool QR_produces_R_factor_with_nonnegative_diagonal () const {
      typedef Combine<ordinal_type, scalar_type> combine_type;
      return combine_type::QR_produces_R_factor_with_nonnegative_diagonal() &&
	reduceBroadcastImpl_.QR_produces_R_factor_with_nonnegative_diagonal();
    }

    /// \brief Internode TSQR with explicit Q factor.
    ///
    /// Call this routine, instead of \c factor() and \c explicit_Q(),
    /// if you want to compute the QR factorization and only want the
    /// Q factor in explicit form (i.e., as a matrix).
    ///
    /// \param R_mine [in/out] View of a matrix with at least as many
    ///   rows as columns.  On input: upper triangular matrix (R
    ///   factor from intranode TSQR); different on each process..  On
    ///   output: R factor from intranode QR factorization; bitwise
    ///   identical on all processes, since it is effectively
    ///   broadcast from Proc 0.
    ///
    /// \param Q_mine [out] View of a matrix with the same number of
    ///   rows as R_mine has columns.  On output: this process'
    ///   component of the internode Q factor.  (Write into the top
    ///   block of this process' entire Q factor, fill the rest of Q
    ///   with zeros, and call intranode TSQR's apply() on it, to get
    ///   the final explicit Q factor.)
    void
    factorExplicit (matview_type R_mine, matview_type Q_mine)
    {
      reduceBroadcastImpl_.factorExplicit (R_mine, Q_mine);
    }

    /// \brief Get cumulative timings for \c factorExplicit().
    ///
    /// Fill in the timings vector with cumulative timings from
    /// factorExplicit().  The vector gets resized to fit all the
    /// timings.
    void 
    getFactorExplicitTimings (std::vector<TimeStats>& stats) const
    {
      reduceBroadcastImpl_.getStats (stats);
    }

    /// \brief Get labels for timings for \c factorExplicit().
    ///
    /// Fill in the labels vector with the string labels for the
    /// timings from factorExplicit().  The vector gets resized to fit
    /// all the labels.
    void
    getFactorExplicitTimingLabels (std::vector<std::string>& labels) const
    {
      reduceBroadcastImpl_.getStatsLabels (labels);
    }

    /// \brief Compute QR factorization of R factors, one per MPI process.
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
    factor (matview_type R_mine)
    {
      VecVec Q_factors, tau_arrays;
      DistTsqrHelper< ordinal_type, scalar_type > helper;
      const ordinal_type ncols = R_mine.ncols();

      std::vector< scalar_type > R_local (ncols*ncols);
      copy_matrix (ncols, ncols, &R_local[0], ncols, R_mine.get(), R_mine.lda());

      const int P = messenger_->size();
      const int my_rank = messenger_->rank();
      const int first_tag = 0;
      std::vector<scalar_type> work (ncols);
      helper.factor_helper (ncols, R_local, my_rank, 0, P-1, first_tag, 
			    messenger_.get(), Q_factors, tau_arrays, work);
      copy_matrix (ncols, ncols, R_mine.get(), R_mine.lda(), &R_local[0], ncols);
      return std::make_pair (Q_factors, tau_arrays);
    }

    //! Apply the result of \c factor() to the distributed matrix C.
    void
    apply (const ApplyType& apply_type,
	   const ordinal_type ncols_C,
	   const ordinal_type ncols_Q,
	   scalar_type C_mine[],
	   const ordinal_type ldc_mine,
	   const FactorOutput& factor_output)
    {
      const bool transposed = apply_type.transposed();

      if (transposed)
	throw std::logic_error("DistTsqr: Applying Q^T or Q^H "
			       "not yet implemented");

      const int P = messenger_->size();
      const int my_rank = messenger_->rank();
      const int first_tag = 0;
      std::vector< scalar_type > C_other (ncols_C * ncols_C);
      std::vector< scalar_type > work (ncols_C);
  
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

    //! Apply the result of \c factor() to compute the explicit Q factor.
    void
    explicit_Q (const ordinal_type ncols_Q,
		scalar_type Q_mine[],
		const ordinal_type ldq_mine,
		const FactorOutput& factor_output)
    {
      const int my_rank = messenger_->rank ();
      fill_matrix (ncols_Q, ncols_Q, Q_mine, ldq_mine, Scalar(0));
      if (my_rank == 0)
	{
	  for (ordinal_type j = 0; j < ncols_Q; ++j)
	    Q_mine[j + j*ldq_mine] = Scalar (1);
	}
      apply (ApplyType::NoTranspose, ncols_Q, ncols_Q, 
	     Q_mine, ldq_mine, factor_output);
    }

  private:
    Teuchos::RCP<MessengerBase<scalar_type> > messenger_;
    DistTsqrRB<ordinal_type, scalar_type> reduceBroadcastImpl_;
  };

} // namespace TSQR

#endif // __TSQR_Tsqr_DistTsqr_hpp
