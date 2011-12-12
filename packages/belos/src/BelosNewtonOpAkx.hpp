#ifndef __Belos_NewtonOpAkx_hpp
#define __Belos_NewtonOpAkx_hpp

//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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

#include <BelosAkx.hpp>
#include <BelosNewtonBasis.hpp>
#include <Teuchos_BLAS.hpp>

// #include <numeric>

namespace Belos {

  /// \class NewtonOpAkx
  /// \brief Newton-basis Akx implementation using abstract operators.
  /// \author Mark Hoemmen
  ///
  /// \warning This is EXPERIMENTAL CODE.  DO NOT RELY ON THIS CODE.
  ///   The interface or implementation may change at any time.
  ///
  template<class Scalar, class MV, class OP>
  class NewtonOpAkx : public OpAkx<Scalar, MV, OP> {
  public:
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    NewtonOpAkx (const Teuchos::RCP<const OP>& A, 
		 const Teuchos::RCP<const OP>& M_left,
		 const Teuchos::RCP<const OP>& M_right,
		 const bool modified) :
      OpAkx (A, M_left, M_right), modified_ (modified) {}

    void 
    computeBasis (const MV& q_last, MV& V_cur, const int s) 
    {
      using Teuchos::Range1D;
      using Teuchos::RCP;
      using Teuchos::rcp_const_cast;
      typedef MultiVecTraits<Scalar, MV> MVT;
      typedef Teuchos::ScalarTraits<Scalar> STS;

      TEUCHOS_TEST_FOR_EXCEPTION(s < 0, std::invalid_argument, 
			 "Number of basis vectors to compute (s = " 
			 << s << ") is invalid; it must be positive.");
      if (s == 0)
	return; // Nothing to do
      TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(V_cur) < s, std::invalid_argument,
			 "You're asking to compute s = " << s << " basis "
			 "vectors, but only have room for " 
			 << MVT::GetNumberVecs(V_cur) << ".");
      RCP<const MV> v_prv = Teuchos::rcpFromRef (q_last);
      int k = 0;
      while (k < s)
	{
	  RCP<MV> v_cur = MVT::CloneViewNonConst (V_cur, Range1D(k,k));
	  if (imagParts_[k] == STS::zero())
	    { // The current shift is real-valued.
	      applyOp (v_prv, v_cur);
	      MVT::MvAddMv (-realParts_[k], v_prv, STS::one(), v_cur, v_cur);
	      // Shift the v_prv view over to the next column of
	      // V_cur.  v_cur gets created on demand, so we don't
	      // have to reassign it.
	      v_prv = rcp_const_cast<const MV> (v_cur);
	      ++k;
	    }
	  else // The current shift has nonzero imaginary part.
	    {  
	      // "Modified" Newton basis trick to ensure real arithmetic,
	      // even if the shifts are complex (as long as any complex
	      // shifts occur in consecutive complex conjugate pairs with
	      // the positive-real-part member of the pair first).
	      //
	      // v_cur := (A - \Re \theta I) v_prv
	      applyOp (v_prv, v_cur);
	      MVT::MvAddMv (-realParts_[k], v_prv, STS::one(), v_cur, v_cur);
	      // v_nxt := (A - \Re \theta I) v_cur + (\Im \theta)^2 v_prv
	      RCP<MV> v_nxt = MVT::CloneViewNonConst (V_cur, Range1D(k+1, k+1));
	      applyOp (v_nxt, v_cur);
	      MVT::MvAddMv (-realParts_[k], v_cur, STS::one(), v_nxt, v_nxt);
	      // FIXME (mfh 14 Feb 2011)
	      //
	      // You might think that we should be concerned that
	      // squaring the imaginary part of the current shift
	      // might overflow.  This is indeed possible, but the
	      // nonmodified Newton basis doesn't really promise any
	      // better (two steps more or less multiply by the square
	      // of the magnitude).  This probably needs some more
	      // analysis, but I don't think it's really broken.
	      {
		const magnitude_type imSqu = imagParts_[k] * imagParts_[k];
		MVT::MvAddMv (imSqu, v_prv, STS::one(), v_nxt, v_nxt);
	      }
	      // Shift the v_prv view over to the next column of
	      // V_cur.  v_cur and v_nxt get created on demand, so we
	      // don't have to reassign them.
	      v_prv = rcp_const_cast<const MV> (v_cur);
	      k += 2;
	    }
	}
    }

    void 
    computeFlexibleBasis (const MV& q_last, 
			  MV& Z_cur, 
			  MV& V_cur, 
			  const int s)
    {
      using Teuchos::Range1D;
      using Teuchos::RCP;
      using Teuchos::rcp_const_cast;
      typedef MultiVecTraits<Scalar, MV> MVT;

      TEUCHOS_TEST_FOR_EXCEPTION(s < 0, std::invalid_argument, 
			 "Number of basis vectors to compute (s = " 
			 << s << ") is invalid; it must be positive.");
      if (s == 0)
	return; // Nothing to do
      TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(V_cur) < s, std::invalid_argument,
			 "You're asking to compute s = " << s << " basis "
			 "vectors, but only have room for " 
			 << MVT::GetNumberVecs(V_cur) << " in V_cur.");
      TEUCHOS_TEST_FOR_EXCEPTION(MVT::GetNumberVecs(Z_cur) < s, std::invalid_argument,
			 "You're asking to compute s = " << s << " basis "
			 "vectors, but only have room for " 
			 << MVT::GetNumberVecs(Z_cur) << " in Z_cur.");
      RCP<const MV> v_prv = Teuchos::rcpFromRef (q_last);
      int k = 0;
      while (k < s)
	{
	  RCP<MV> z_cur = MVT::CloneViewNonConst (Z_cur, Range1D(k,k));
	  RCP<MV> v_cur = MVT::CloneViewNonConst (V_cur, Range1D(k,k));
	  applyFlexibleOp (v_prv, z_cur, v_cur);
	  if (imagParts_[k] == STS::zero())
	    {
	      // We've already applied the full operator to v_prv,
	      // producing v_cur.  Now we just need to apply the shift
	      // to v_cur.
	      MVT::MvAddMv (-realParts_[k], v_prv, STS::one(), v_cur, v_cur);
	      // Shift the v_prv view over to the next column of
	      // V_cur.  z_cur and v_cur get created on demand, so we
	      // don't have to reassign them.
	      v_prv = rcp_const_cast<const MV> (v_cur);
	      ++k;
	    }
	  else
	    {
	      // Apply the first shift (realParts_[k]) to v_cur.
	      MVT::MvAddMv (-realParts_[k], v_prv, STS::one(), v_cur, v_cur);
	      RCP<MV> z_nxt = MVT::CloneViewNonConst (Z_cur, Range1D(k+1, k+1));
	      RCP<MV> v_nxt = MVT::CloneViewNonConst (V_cur, Range1D(k+1, k+1));
	      // After the next line, z_nxt will be just the result of
	      // applying the right preconditioner to v_cur (to which
	      // the first shift has already been applied).  Thus,
	      // it's only necessary to apply the second shift to
	      // v_nxt.
	      applyFlexibleOp (v_cur, z_nxt, v_nxt);
	      MVT::MvAddMv (-realParts_[k], v_cur, STS::one(), v_nxt, v_nxt);
	      {
		const magnitude_type imSqu = imagParts_[k] * imagParts_[k];
		MVT::MvAddMv (imSqu, v_prv, STS::one(), v_nxt, v_nxt);
	      }
	      // Shift the v_prv view over to the next column of
	      // V_cur.  z_cur, v_cur, z_nxt, and v_nxt get created on
	      // demand, so we don't have to reassign them.
	      v_prv = rcp_const_cast<const MV> (v_cur);
	      k += 2;
	    }
	}
    }

    ///
    /// \note s may become s+1 if the s-th shift is the first member
    ///   of a complex conjugate pair.
    std::pair<Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,Scalar>, int>
    changeOfBasisMatrix (int s)
    {
      using Teuchos::Array;
      using Teuchos::ArrayView;
      typedef Teuchos::SerialDenseMatrix<int,Scalar> mat_type;
      typedef Teuchos::ScalarTraits<Scalar> STS;
      typedef Teuchos::ScalarTraits<magnitude_type> STM;

      // Only keep the first s shifts, counting multiplicities.
      int numShifts = 0, numUniqueShifts = 0;
      Array<int> mults;
      while (numShifts < s)
	{
	  const int k = numUniqueShifts;

	  // If the current shift is the first member of a complex
	  // conjugate pair (i.e., has positive imaginary part), see
	  // if there is room to add all multiplicities for both
	  // members of the complex conjugate pair.
	  if (imagParts_[k] > STM::zero())
	    { 
	      int curMult = 1;
	      // Can we add at least two more shifts (counting multiplicities)?
	      if (numShifts + 1 >= s)
		{
		  // We either have to increase the room so we can add
		  // both, or not add either member of the pair.  We
		  // decide what to do based on how many shifts we've
		  // added thus far: if <= 1, we add the pair and
		  // increase s accordingly, else we don't add the
		  // pair and stop.
		  if (numShifts > 1)
		    break; // Don't add any more shifts
		}
	      else 
		{ // We can add both members of the pair, with
		  // multiplicity at least 1.  Let's see how many
		  // pairs we can add.
		  const int leftover = numShifts + 2*multiplicities_[k] - s;
		  int curMult;
		  if (leftover > 0)
		    // If even, accept, else round up one.
		    curMult = (leftover % 2 == 0) ? (leftover/2) : (leftover/2 + 1);
		  else // We can add both shifts with all their multiplicities.
		    curMult = multiplicities_[k];
		}
	      numUniqueShifts += 2;
	      numShifts += 2*curMult;
	      mults.push_back (curMult);
	      mults.push_back (curMult);
	    }
	  else
	    {
	      const int curMult = 
		(numShifts + multiplicities_[k] > s) ? 
		(numShifts + multiplicities_[k] - s) : 
		multiplicities_[k];
	      numUniqueShifts++;
	      numShifts += curMult;
	      mults.push_back (curMult);
	    }
	}
      // Views of the shifts and their multiplicities of interest.
      ArrayView<const magnitude_type> realParts = realParts_.view (0, numUniqueShifts);
      ArrayView<const magnitude_type> imagParts = imagParts_.view (0, numUniqueShifts);

      // FIXME (mfh 09 Feb 2011) We must also recompute the matrix
      // when the shifts have changed.
      if (B_.is_null() || 
	  (B_->numRows() != s+1 || B_->numCols() != s) ||
	  shiftsChanged_)
	B_ = makeChangeOfBasisMatrix (realParts, imagParts, mults, modified_);
      return B_;
    }


    void
    updateBasis (const Teuchos::ArrayView<const magnitude_type>& realParts,
		 const Teuchos::ArrayView<const magnitude_type>& imagParts,
		 const Teuchos::ArrayView<const int>& multiplicities,
		 const int numValues)
    {
      using Teuchos::Array;
      using Teuchos::ArrayView;
      typedef Teuchos::SerialDenseMatrix<int,Scalar> mat_type;
      typedef Teuchos::ScalarTraits<Scalar> STS;
      typedef Teuchos::ScalarTraits<magnitude_type> STM;

      checkModifiedNewtonShifts (realParts, imagParts, multiplicities);

      // Only take in new shifts if there are more new shifts than old
      // ones.
      if (numValues < realParts_.size())
	return;
      else
	{
	  realParts_.resize (numValues);
	  imagParts_.resize (numValues);
	  mults_.resize (numValues);
	  
	  realParts_.assign (realParts.begin(), realParts.end());
	  imagParts_.assign (imagParts.begin(), imagParts.end());
	  mults_.assign (multiplicities.begin(), multiplicities.end());

	  // Recompute the change-of-basis matrix, since the shifts
	  // have changed.
	  typedef Teuchos::ScalarTraits<Scalar> STS;
	  if (B_.is_null() || B_->numRows() != numValues+1 || 
	      B_->numCols() != numValues)
	    B_ = makeNewtonChangeOfBasisMatrix<Scalar, STS::isComplex> (realParts_, imagParts_, mults_, modified_);
	}
    }

  private:

    //! Monomial change-of-basis matrix, cached to avoid allocations.
    Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,Scalar> > B_;

    //! Real parts of the (possibly complex-valued) shifts.
    Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> realParts_;
    /// Imaginary parts of the (possibly complex-valued) shifts.  If
    /// the shifts are all real, this array contains all zero entries,
    /// but always has the same number of entries as realParts_.
    Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> imagParts_;

    //! Whether or not this is the modified Newton basis.
    bool modified_;
  };

} // namespace Belos

#endif // __Belos_NewtonOpAkx_hpp
