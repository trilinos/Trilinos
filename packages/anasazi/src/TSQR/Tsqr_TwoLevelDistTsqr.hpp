// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// ***********************************************************************
// @HEADER

#ifndef __TSQR_Tsqr_TwoLevelDistTsqr_hpp
#define __TSQR_Tsqr_TwoLevelDistTsqr_hpp

#include <Tsqr_DistTsqr.hpp>
#include <Tsqr_MpiCommFactory.hpp>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class TwoLevelDistTsqr
  /// \brief Interprocess part of TSQR
  ///
  /// Interprocess part of TSQR, composed of an internode part and an
  /// intranode part (but only working between MPI processes, not
  /// within a process).
  template< class LocalOrdinal, 
	    class Scalar, 
	    class DistTsqrType = DistTsqr< LocalOrdinal, Scalar > >
  class TwoLevelDistTsqr {
  public:
    typedef Scalar scalar_type;
    typedef LocalOrdinal ordinal_type;
    typedef DistTsqrType dist_tsqr_type;
    typedef typename dist_tsqr_type::rank_type rank_type;
    typedef typename Teuchos::RCP< dist_tsqr_type > dist_tsqr_ptr;
    typedef typename dist_tsqr_type::FactorOutput DistTsqrFactorOutput;
    typedef std::pair< DistTsqrFactorOutput, DistTsqrFactorOutput > FactorOutput;

    /// \brief Constructor
    ///
    TwoLevelDistTsqr () :
      worldMess_ (TSQR::MPI::makeMpiCommWorld()),
      nodeDistTsqr_ (TSQR::MPI::makeMpiCommNode()),
      networkDistTsqr_ (TSQR::MPI::makeMpiCommNetwork())
    {}

    ///
    /// \brief Destructor
    ///
    ~TwoLevelDistTsqr () {}

    /// Whether or not all diagonal entries of the R factor computed
    /// by the QR factorization are guaranteed to be nonnegative.
    bool QR_produces_R_factor_with_nonnegative_diagonal () const {
      return nodeDistTsqr_->QR_produces_R_factor_with_nonnegative_diagonal() &&
	networkDistTsqr_->QR_produces_R_factor_with_nonnegative_diagonal();
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
    FactorOutput
    factor (MatView< LocalOrdinal, Scalar > R_mine)
    {
      DistTsqrFactorOutput nodeOutput = nodeDistTsqr_->factor (R_mine);
      DistTsqrFactorOutput networkOutput = networkDistTsqr_->factor (R_mine);
      return std::make_pair (nodeOutput, networkOutput);
    }

    void
    apply (const ApplyType& applyType,
	   const LocalOrdinal ncols_C,
	   const LocalOrdinal ncols_Q,
	   Scalar C_mine[],
	   const LocalOrdinal ldc_mine,
	   const FactorOutput& factorOutput)
    {
      if (applyType.transposed())
	{
	  nodeDistTsqr_->apply (applyType, ncols_C, ncols_Q, 
				C_mine, ldc_mine, factorOutput.first);
	  networkDistTsqr_->apply (applyType, ncols_C, ncols_Q, 
				   C_mine, ldc_mine, factorOutput.second);
	}
      else 
	{
	  networkDistTsqr_->apply (applyType, ncols_C, ncols_Q, 
				   C_mine, ldc_mine, factorOutput.second);
	  nodeDistTsqr_->apply (applyType, ncols_C, ncols_Q, 
				C_mine, ldc_mine, factorOutput.first);
	}
    }

    void
    explicit_Q (const LocalOrdinal ncols_Q,
		Scalar Q_mine[],
		const LocalOrdinal ldq_mine,
		const FactorOutput& factorOutput)
    {
      typedef MatView< LocalOrdinal, Scalar > matview_type;
      matview_type Q_view (ncols_Q, ncols_Q, Q_mine, ldq_mine, Scalar(0));
      Q_view.fill (Scalar(0));
	
      const rank_type myRank = worldMess_->rank();
      if (myRank == 0)
	{
	  if (networkMess_->rank() != 0)
	    {
	      std::ostringstream os;
	      os << "My rank with respect to MPI_COMM_WORLD is 0, but my rank "
		"with respect to MPI_COMM_NETWORK is nonzero (= " 
		 << networkMess_->rank() << ").  We could deal with this by "
		"swapping data between those two ranks, but we haven\'t "
		"implemented that fix yet.";
	      throw std::logic_error (os.str());
	    }
	  for (LocalOrdinal j = 0; j < ncols_Q; ++j)
	    Q_view(j, j) = Scalar (1);
	}
      apply (ApplyType::NoTranspose, ncols_Q, ncols_Q, 
	     Q_mine, ldq_mine, factorOutput);
    }

  private:
    Teuchos::RCP< MessengerBase< Scalar > > worldMess_;
    dist_tsqr_ptr nodeDistTsqr_, networkDistTsqr_;
  };

} // namespace TSQR

#endif // __TSQR_Tsqr_TwoLevelDistTsqr_hpp
