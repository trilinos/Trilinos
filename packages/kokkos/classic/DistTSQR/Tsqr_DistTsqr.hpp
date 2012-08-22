//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>
#include <Teuchos_ScalarTraits.hpp>
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
  class DistTsqr : public Teuchos::ParameterListAcceptorDefaultBase {
  public:
    typedef Scalar scalar_type;
    typedef LocalOrdinal ordinal_type;
    typedef MatView<ordinal_type, scalar_type > matview_type;
    typedef std::vector<std::vector<scalar_type> > VecVec;
    typedef std::pair<VecVec, VecVec> FactorOutput;
    typedef int rank_type;

  private:
    typedef Teuchos::ScalarTraits<Scalar> STS;

  public:

    /// \brief Constructor (that accepts a parameter list).
    ///
    /// \param plist [in/out] List of parameters for configuring TSQR.
    ///   The specific parameter keys that are read depend on the TSQR
    ///   implementation.  For details, call \c getValidParameters()
    ///   and examine the documentation embedded therein.
    DistTsqr (const Teuchos::RCP<Teuchos::ParameterList>& plist)
    {
      setParameterList (plist);
    }

    //! Constructor (that uses default parameters).
    DistTsqr ()
    {
      setParameterList (Teuchos::null);
    }

    void 
    setParameterList (const Teuchos::RCP<Teuchos::ParameterList>& plist)
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;

      RCP<ParameterList> params = plist.is_null() ? 
	parameterList (*getValidParameters()) : plist;

      // Do nothing for now, other than store the list.
      this->setMyParamList (params);
    }

    Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters() const
    {
      return Teuchos::parameterList ("DistTsqr"); // Empty list for now.
    }
    
    /// \brief Finish initialization using the messenger object.
    ///
    /// \param messenger [in/out] An object handling communication
    ///   between (MPI) process(es).
    void init (const Teuchos::RCP<MessengerBase<scalar_type> >& messenger)
    {
      messenger_ = messenger;
      reduceBroadcastImpl_ = 
	Teuchos::rcp (new DistTsqrRB<ordinal_type, scalar_type> (messenger_));
    }

    /// \brief Rank of this (MPI) process.
    /// 
    /// Rank is computed via MPI_Comm_rank() on the underlying
    /// communicator, if the latter is an MPI communicator.  If it's a
    /// serial "communicator," the rank is always zero.
    rank_type rank() const { 
      TEUCHOS_TEST_FOR_EXCEPTION(! ready(), std::logic_error,
				 "Before using DistTsqr computational methods, "
				 "you must first call init() with a valid "
				 "MessengerBase instance.");
      return messenger_->rank(); 
    }

    /// \brief Total number of MPI processes in this communicator.
    ///
    /// The size is communicated via MPI_Comm_size() on the underlying
    /// communicator, if the latter is an MPI communicator.  If it's a
    /// serial "communicator," the size is always one.
    rank_type size() const { 
      TEUCHOS_TEST_FOR_EXCEPTION(! ready(), std::logic_error,
				 "Before using DistTsqr computational methods, "
				 "you must first call init() with a valid "
				 "MessengerBase instance.");
      return messenger_->size(); 
    }

    /// \brief Destructor.
    ///
    /// The destructor doesn't need to do anything, thanks to smart
    /// pointers.
    virtual ~DistTsqr () {}

    /// \brief Does the R factor have a nonnegative diagonal?
    ///
    /// DistTsqr implements a QR factorization (of a distributed
    /// matrix with a special structure).  Some, but not all, QR
    /// factorizations produce an R factor whose diagonal may include
    /// negative entries.  This Boolean tells you whether DistTsqr
    /// promises to compute an R factor whose diagonal entries are all
    /// nonnegative.
    bool QR_produces_R_factor_with_nonnegative_diagonal () const {
      TEUCHOS_TEST_FOR_EXCEPTION(! ready(), std::logic_error,
				 "Before using DistTsqr computational methods, "
				 "you must first call init() with a valid "
				 "MessengerBase instance.");
      typedef Combine<ordinal_type, scalar_type> combine_type;
      return combine_type::QR_produces_R_factor_with_nonnegative_diagonal() &&
	reduceBroadcastImpl_->QR_produces_R_factor_with_nonnegative_diagonal();
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
    ///
    /// \param forceNonnegativeDiagonal [in] If true, then (if
    ///   necessary) do extra work (modifying both the Q and R
    ///   factors) in order to force the R factor to have a
    ///   nonnegative diagonal.
    void
    factorExplicit (matview_type R_mine, 
		    matview_type Q_mine,
		    const bool forceNonnegativeDiagonal=false)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(! ready(), std::logic_error,
				 "Before using DistTsqr computational methods, "
				 "you must first call init() with a valid "
				 "MessengerBase instance.");
      reduceBroadcastImpl_->factorExplicit (R_mine, Q_mine, 
					    forceNonnegativeDiagonal);
    }

    /// \brief Get cumulative timings for \c factorExplicit().
    ///
    /// Fill in the timings vector with cumulative timings from
    /// factorExplicit().  The vector gets resized to fit all the
    /// timings.
    void 
    getFactorExplicitTimings (std::vector<TimeStats>& stats) const
    {
      TEUCHOS_TEST_FOR_EXCEPTION(! ready(), std::logic_error,
				 "Before using DistTsqr computational methods, "
				 "you must first call init() with a valid "
				 "MessengerBase instance.");
      reduceBroadcastImpl_->getStats (stats);
    }

    /// \brief Get labels for timings for \c factorExplicit().
    ///
    /// Fill in the labels vector with the string labels for the
    /// timings from factorExplicit().  The vector gets resized to fit
    /// all the labels.
    void
    getFactorExplicitTimingLabels (std::vector<std::string>& labels) const
    {
      TEUCHOS_TEST_FOR_EXCEPTION(! ready(), std::logic_error,
				 "Before using DistTsqr computational methods, "
				 "you must first call init() with a valid "
				 "MessengerBase instance.");
      reduceBroadcastImpl_->getStatsLabels (labels);
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
      TEUCHOS_TEST_FOR_EXCEPTION(! ready(), std::logic_error,
				 "Before using DistTsqr computational methods, "
				 "you must first call init() with a valid "
				 "MessengerBase instance.");
      VecVec Q_factors, tau_arrays;
      DistTsqrHelper<ordinal_type, scalar_type> helper;
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
      TEUCHOS_TEST_FOR_EXCEPTION(! ready(), std::logic_error,
				 "Before using DistTsqr computational methods, "
				 "you must first call init() with a valid "
				 "MessengerBase instance.");
      const bool transposed = apply_type.transposed();
      TEUCHOS_TEST_FOR_EXCEPTION(transposed, std::logic_error,
				 "DistTsqr: Applying Q^T or Q^H has not yet "
				 "been implemented.");
      const int P = messenger_->size();
      const int my_rank = messenger_->rank();
      const int first_tag = 0;
      std::vector<scalar_type> C_other (ncols_C * ncols_C);
      std::vector<scalar_type> work (ncols_C);
  
      const VecVec& Q_factors = factor_output.first;
      const VecVec& tau_arrays = factor_output.second;

      // assert (Q_factors.size() == tau_arrays.size());
      const int cur_pos = Q_factors.size() - 1;
      DistTsqrHelper<ordinal_type, scalar_type> helper;
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
      TEUCHOS_TEST_FOR_EXCEPTION(! ready(), std::logic_error,
				 "Before using DistTsqr computational methods, "
				 "you must first call init() with a valid "
				 "MessengerBase instance.");
      const int myRank = messenger_->rank ();
      fill_matrix (ncols_Q, ncols_Q, Q_mine, ldq_mine, STS::zero());
      if (myRank == 0) {
	for (ordinal_type j = 0; j < ncols_Q; ++j)
	  Q_mine[j + j*ldq_mine] = STS::one();
      }
      apply (ApplyType::NoTranspose, ncols_Q, ncols_Q, 
	     Q_mine, ldq_mine, factor_output);
    }

  private:
    Teuchos::RCP<MessengerBase<scalar_type> > messenger_;
    Teuchos::RCP<DistTsqrRB<ordinal_type, scalar_type> > reduceBroadcastImpl_;

    /// \brief Whether this object is ready to perform computations.
    /// 
    /// It is <i>not</i> ready until after \c init() has been called.
    bool ready() const {
      return ! messenger_.is_null() && ! reduceBroadcastImpl_.is_null();
    }
  };

} // namespace TSQR

#endif // __TSQR_Tsqr_DistTsqr_hpp
