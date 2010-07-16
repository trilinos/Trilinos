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


/// \file AnasaziTsqrOrthoManager.hpp
/// \brief Orthogonalization manager based on Tall Skinny QR (TSQR)

#ifndef __AnasaziTsqrOrthoManager_hpp
#define __AnasaziTsqrOrthoManager_hpp

#include "AnasaziConfigDefs.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziTsqrAdaptor.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziMatOrthoManager.hpp"
#include "Teuchos_LAPACK.hpp"

namespace Anasazi {

  template< class ScalarType, class MV, class OP >
  class TsqrOrthoManager : 
    public MatOrthoManager< ScalarType, MV, OP > 
  {
  private:
    typedef typename Teuchos::ScalarTraits< ScalarType >::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits< ScalarType >    SCT;
    typedef Teuchos::ScalarTraits< MagnitudeType > SCTM;
    typedef MultiVecTraits< ScalarType, MV >       MVT;
    typedef OperatorTraits< ScalarType, MV, OP >   OPT;
    std::string dbgstr;

    typedef typename Anasazi::TsqrAdaptor< ScalarType, MV >::adaptor_type tsqr_adaptor_type;
    typedef Teuchos::RCP< adaptor_type > tsqr_adaptor_ptr;

  public:

    /// \brief Constructor
    ///
    /// Constructor specifying re-orthogonalization tolerance.
    ///
    /// \param tsqrParams [in] Configuration parameters for TSQR.
    /// \param Op [in] Inner product.  Don't set to anything not 
    ///   Teuchos::null.
    TsqrOrthoManager (const Teuchos::ParameterList& tsqrParams,
		      Teuchos::RCP<const OP> Op = Teuchos::null)
      : MatOrthoManager< ScalarType, MV, OP >(Op), 
	tsqrParams_ (tsqrParams),
	tsqrAdaptor_ (Teuchos::null),
	Q_ (Teuchos::null)
    {
      Teuchos::LAPACK< int, MagnitudeType > lapack;
      eps_ = lapack.LAMCH('E');
      if (Op != Teuchos::null)
	throw OrthoError("Sorry, TsqrOrthoManager doesn\'t work with "
			 "non-null Op");
    }

    // FIXME (mfh 16 Jul 2010) setOp() in the base class
    // (MatOrthoManager) is non-virtual, but should be...
    void 
    setOp (Teuchos::RCP< const OP > Op)
    {
      if (Op != Teuchos::null)
	throw OrthoError("Sorry, TsqrOrthoManager doesn\'t work with "
			 "non-null Op");
      else
	_Op = Op;
    }

    ~TsqrOrthoManager() {};

    typedef Teuchos::RCP< MV >       mv_ptr;
    typedef Teuchos::RCP< const MV > const_mv_ptr;
    typedef Teuchos::Array< const_mv_ptr >                       const_prev_mvs_type;
    typedef Teuchos::SerialDenseMatrix< int, ScalarType >        serial_matrix_type;
    typedef Teuchos::RCP< serial_matrix_type >                   serial_matrix_ptr;
    typedef Teuchos::Array< Teuchos::RCP< serial_matrix_type > > prev_coeffs_type;

    void 
    projectMat (MV& X, 
		const_prev_mvs_type Q,
		prev_coeffs_type C = Teuchos::tuple (serial_matrix_ptr (Teuchos::null)),
		mv_ptr MX = Teuchos::null,
		const_prev_mvs_type MQ = Teuchos::tuple (const_mv_ptr (Teuchos::null))) const
    {
      if (getOp() != Teuchos::null)
	throw OrthoError("Sorry, TsqrOrthoManager doesn\'t work with "
			 "non-null Op");
      else if (MX != Teuchos::null)
	throw OrthoError("Sorry, TsqrOrthoManager::normalizeMat() "
			 "doesn\'t work with MX non-null yet");

      // Test for quick exit
      const int num_Q_blocks = Q.length();
      if (num_Q_blocks == 0)
	return;
      const int nrows_X = MVT::GetVecLength (X);
      const int ncols_X = MVT::GetNumberVecs (X);
      if (nrows_X == 0 || ncols_X == 0)
	return;

      //
      // Make sure that for each i, the dimensions of X and Q[i] are
      // compatible.
      //
      int ncols_Q_total = 0; // total over all Q blocks
      for (int i = 0; i < num_Q_blocks; ++i)
	{
	  const int nrows_Q = MVT::GetVecLength (*Q[i]);
	  TEST_FOR_EXCEPTION( (nrows_Q != nrows_X), 
			      std::invalid_argument,
			      "Anasazi::TsqrOrthoManager::projectMat(): "
			      "Size of X not consistant with size of Q" );
	  ncols_Q_total += MVT::GetNumberVecs (*Q[i]);
	}
      if (ncols_Q_total == 0)
	return; // quick exit

      // If we don't have enough C, expanding it creates null references.
      // If we have too many, resizing just throws away the later ones.
      // If we have exactly as many as we have Q, this call has no effect.
      C.resize (num_Q_blocks);
      prev_coeffs_type newC (num_Q_blocks);
      for (int i = 0; i < num_Q_blocks; ++i) 
	{
	  const int ncols_Q_i = MVT::GetNumberVecs (*Q[i]);
	  // Create a new C[i] if necessary.
	  if (C[i] == Teuchos::null)
	    C[i] = Teuchos::rcp (new serial_matrix_type (ncols_Q_i, ncols_X));
	  // Fill C[i] with zeros.
	  C[i]->putScalar (ScalarType(0));
	  // Create newC[i] as a clone of C[i].  (All that really
	  // matters is that is has the same dimensions and is filled
	  // with zeros; cloning C[i] accomplishes both in this case.)
	  newC[i] = Teucos::rcp (new serial_matrix_type (*C[i]));
	}

      // Norms, and "inverse norms" (set to zero if norm is zero) of
      // the columns of X.
      std::vector< ScalarType > normX (ncols_X), invnormX (ncols_X);

      // Whether (another / a) round of (block) Gram-Schmidt
      // orthogonalization should be done.
      bool doGramSchmidt = true;

      // Orthogonalization tolerance
      MagnitudeType tolerance = MagnitudeType(1) / SCTM::squareroot(eps_);

      // Block Gram-Schmidt (re)orthogonalization loop: one loop
      // iteration corresponds to one projection step.
      while (doGramSchmidt)
	{
	  // "Modified Gram-Schmidt" version of Block Gram-Schmidt:
	  //
	  // $C_i := Q_i^* \cdot X$, and $X := X - Q_i \cdot C_i$.
	  for (int i = 0; i < num_Q_blocks; ++i)
	    {
	      innerProdMat (*Q[i], X, *newC[i], Teuchos::null);
	      MVT::MvTimesMatAddMv (ScalarType(-1), *Q[i], *newC[i], ScalarType(1), X);
	    }

	  // Decide whether we need to reorthogonalize.  We begin by
	  // computing the largest column norm of 
	  // \[
	  // C = 
	  // \begin{pmatrix}
	  //   C_0 \\
          //   C_1 \\
	  //   \vdots \\
	  //   C_{\text{num\_Q\_blocks} - 1} \\
          // \end{pmatrix}
	  //
	  MagnitudeType maxCNorm (0);

	  for (int j = 0; j < ncols_X; ++j)
	    {
	      MagnitudeType total (0);
	      for (int k = 0; i < num_Q_blocks; ++i)
		{
		  const int nrows_C_k = MVT::GetNumberVecs (*Q[k]);
		  scalar_matrix_type& C_k = *newC[k];
		  const ScalarType* const C_k_col_j = &C_k(0,j);

		  for (int i = 0; i < nrows_C_k; ++i)
		    total += SCT::magnitude (C_k_col_j[i]) * 
		      SCT::magnitude (C_k_col_j[i]);
		}
	      maxNorm = (total > maxNorm) ? total : maxNorm;
	    }

	  // FIXME (mfh 16 Jul 2010)
	  //
	  // 0.36 appears as a magic number in
	  // AnasaziSVQBOrthoManager.hpp.  Should justify it there or
	  // at least name it and document it as a magic number.
	  const MagnitudeType magicThreshold (0.36);

	  // FIXME (mfh 16 Jul 2010)
	  //
	  // Does this "magic threshold" test work, if we don't scale
	  // the columns of X?
	  //
	  // Test whether we should reorthogonalize.
	  if (maxNorm < magicThreshold)
	    doGramSchmidt = false;

	  // Accumulate results into C.
	  for (int i = 0; i < num_Q_blocks; ++i)
	    (*C[i]) += *newC[i];

	} // end while (doGramSchmidt)
    }

    int 
    normalizeMat (MV& X,
		  serial_matrix_ptr B = Teuchos::null,
		  mv_ptr MX = Teuchos::null) const
    {
      if (getOp() != Teuchos::null)
	throw OrthoError("Sorry, TsqrOrthoManager doesn\'t work with "
			 "non-null Op");
      else if (MX != Teuchos::null)
	throw OrthoError("Sorry, TsqrOrthoManager::normalizeMat() "
			 "doesn\'t work with MX non-null yet");

      // The TSQR adaptor object requires a specific MV object for
      // initialization.  As long as subsequent MV objects use the
      // same map and communicator (e.g., the same Tpetra::Map
      // resp. Teuchos::Comm<int> objects), we don't need to
      // reinitialize the adaptor.  
      //
      // FIXME (mfh 15 Jul 2010) If tsqrAdaptor_ has already been
      // initialized, check to make sure that X has the same map and
      // communicator as the multivector previously used to initialize
      // tsqrAdaptor_.
      if (tsqrAdaptor_ == Teuchos::null)
	tsqrAdaptor_ = Teuchos::rcp (X, tsqrParams_);
      
      // MVT returns int for these two quantities, even though
      // local_ordinal_type of the MV may be some other type.
      const int nrows = MVT::GetVecLength (X);
      const int ncols = MVT::GetNumberVecs (X);

      // TSQR's rank-revealing part doesn't work unless B is provided.
      const bool B_is_null_on_input = (B == Teuchos::null);
      if (B_is_null_on_input)
	B = Teuchos::rcp (new serial_matrix_type (ncols, ncols));

      // Q_ is temporary workspace.  It must have the same dimensions
      // as X.  If not, we have to reallocate.  We also have to
      // allocate (not "re-") if we haven't allocated Q_ before.  (We
      // can't allocate Q_ until we have some X, so we need a
      // multivector as the "prototype.")
      if (Q_ == Teuchos::null || 
	  nrows != MVT::GetVecLength(Q_) || 
	  ncols != MVT::GetVecLength(Q_))
	// Q is allocated to have the same dimensions as X.  Contents
	// of X are not copied.  We will fill in Q with the (explicit)
	// Q factor of X, suitably modified if X is not full rank.
	Q_ = MVT::Clone (X, ncols);

      // Compute rank-revealing decomposition (in this case, TSQR of X
      // followed by SVD of the R factor and appropriate updating of
      // the resulting Q factor) of X.  X is modified in place, and Q
      // contains the results.
      typedef typename tsqr_adaptor_type::factor_output_type factor_output_type;

      // Compute TSQR and SVD of X.  Resulting orthogonal vectors go
      // into Q_, and coefficients (not necessarily upper triangular)
      // go into B.
      OrdinalType rank;
      try {
	factor_output_type factorOutput = tsqrAdaptor_->factor (X, *B);
	tsqrAdaptor_->explicitQ (X, factorOutput, Q_);
	rank = tsqrAdaptor_->revealRank (Q_, B, relativeTolerance);
      } catch (std::exception& e) {
	throw OrthoError (e.what());
      }
      // Now we should copy Q_ back into X, but don't do it yet: if we
      // want to fill the last ncols-rank columns with random data, we
      // should do so in Q_, because it's fresh in the cache.
      if (fillWhenNotFullRank_)
	{
	  // If X did not have full (numerical rank), augment the last
	  // ncols-rank columns of X with random data.
	  const int ncolsToFill = ncols - rank;
	  if (ncolsToFill > 0)
	    {
	      // ind: Indices of columns of X to fill with random data.
	      std::vector< int > fillIndices (ncolsToFill);
	      for (int j = 0; j < ncolsToFill; ++j)
		fillIndices[j] = j + rank;

	      mv_ptr Q_null = MVT::CloneViewNonConst (Q_, fillIndices);
	      MVT::MvRandom (*Q_null);
	      Q_null = Teuchos::null;
	    }
	}
      // TODO (mfh 16 Jul 2010) Copy Q_ back into X.
      MVT::DeepCopy (X, Q_);

      // Get rid of B if the user didn't want it as output.
      if (B_is_null_on_input)
	B = Teuchos::null;

      return rank;
    }


    int 
    projectAndNormalizeMat (MV &X,
			    const_prev_mvs_type Q,
			    prev_coeffs_type C = Teuchos::tuple (serial_matrix_ptr (Teuchos::null)),
			    serial_matrix_ptr B = Teuchos::null, 
			    mv_ptr MX = Teuchos::null,
			    const_prev_mvs_type MQ = Teuchos::tuple (const_mv_ptr (Teuchos::null))) const 
    {
      projectMat (X, Q, C, MX, MQ);
      return normalizeMat (X, B, MX);
    }


    MagnitudeType 
    orthonormErrorMat (const MV &X, 
		       const_mv_ptr MX = Teuchos::null) const
    {
      const ScalarType ONE = SCT::one();
      const int rank = MVT::GetNumberVecs(X);
      Teuchos::SerialDenseMatrix< int, ScalarType > xTx(rank, rank);
      innerProdMat (X, X, xTx, MX, MX);
      for (int i = 0; i < rank; ++i) {
	xTx(i,i) -= ONE;
      }
      return xTx.normFrobenius();
    }

    MagnitudeType 
    orthogErrorMat (const MV &X, 
		    const MV &Y,
		    mv_ptr MX = Teuchos::null, 
		    mv_ptr MY = Teuchos::null) const
    {
      const int r1 = MVT::GetNumberVecs (X);
      const int r2 = MVT::GetNumberVecs (Y);
      serial_matrix_type xTx (r1, r2);
      innerProdMat (X, Y, xTx, MX, MY);
      return xTx.normFrobenius();
    }

  private:
    tsqr_adaptor_ptr tsqrAdaptor_;
    mv_ptr Q_;
    MagnitudeType eps_;
  };

} // namespace Anasazi

#endif // __AnasaziTsqrOrthoManager_hpp

