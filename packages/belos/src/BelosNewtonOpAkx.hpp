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
// #include <numeric>

namespace Belos {

  /// \fn makeNewtonChangeOfBasisMatrix 
  /// \author Mark Hoemmen
  ///
  /// Return the change-of-matrix matrix for the Newton-basis matrix
  /// powers kernel, with as many shifts as there are entries of
  /// realParts (imagParts must have the same number of entries),
  /// where (realParts[k], imagParts[k]) has multiplicity
  /// multiplicities[k].
  ///
  /// \note The implementation depends syntactically on whether Scalar
  ///   is a complex type or a real (! isComplex) type; hence, the
  ///   second template parameter.
  /// 
  /// \param realParts [in] Real parts of the shifts.
  /// \param imagParts [in] Imaginary parts of the shifts.
  /// \param multiplicities [in] Multiplicities of the shifts.
  /// \param modified [in] Whether to compute the change-of-basis
  ///   matrix for the modified Newton basis, or the regular Newton
  ///   basis.
  template<class Scalar, bool isComplex>
  Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,Scalar> >
  makeNewtonChangeOfBasisMatrix (const Teuchos::ArrayView<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& realParts,
				 const Teuchos::ArrayView<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& imagParts,
				 const Teuchos::ArrayView<const int>& multiplicities,
				 const bool modified);

  /// \fn makeModifiedNewtonChangeOfBasisMatrix
  /// \author Mark Hoemmen
  ///
  /// Return the change-of-matrix matrix for the (modified Newton)
  /// basis matrix powers kernel, with as many shifts as there are
  /// entries of realParts (imagParts must have the same number of
  /// entries), where (realParts[k], imagParts[k]) has multiplicity
  /// multiplicities[k].
  //
  /// \note The implementation does not depend on whether Scalar is a
  ///   complex type.  This is because complex shifts may occur for
  ///   both real and complex matrices, and complex shift k is stored 
  ///   as the pair (realParts[k], imagParts[k]).
  template<class Scalar>
  Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,Scalar> >
  makeModifiedNewtonChangeOfBasisMatrix (const Teuchos::ArrayView<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& realParts,
					 const Teuchos::ArrayView<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& imagParts,
					 const Teuchos::ArrayView<const int>& multiplicities)
  {
    using Teuchos::RCP;
    typedef Teuchos::SerialDenseMatrix<int,Scalar> mat_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType magnitude_type;

    RCP<mat_type> B (new mat_type (s+1, s));

    int j = 0; // Current column of the B matrix
    int kk = 0;
    while (kk < multiplicities.size())
      {
	const magnitude_type curRealPart = realParts[kk];
	const magnitude_type curImagPart = imagParts[kk];
	const int curMult = multiplicities[kk];
	
	// If the current shift (kk) is the first element of a complex
	// conjugate pair, then the following must hold, else there is
	// an error:
	//
	// 1. This is at least one more shift left
	// 2. lambda(kk) == conj(lambda(kk+1))
	// 3. multiplicities[kk] == multiplicities[kk+1]
	if (imagParts[kk] > 0)
	  {
	    // Check for any of the possible error conditions
	    if (kk == multiplicities.size() - 1 || 
		curRealPart != realParts[kk+1] ||
		curImagPart != -imagParts[kk] ||
		curMult != multiplicities[kk+1])
	      {
		std::ostringstream os;
		os << "With the modified Newton basis, any complex shifts "
		  "must occur in consecutive complex conjugate pairs, with "
		  "the first member of each pair having positive imaginary "
		  "part.  However, shift number " << kk << " = (" 
		   << curRealPart << ", " << curImagPart << ") ";
		if (kk == multiplicities.size() - 1)
		  {
		    os << "is the last shift in the given order.";
		  }
		else if (curRealPart != realParts[kk+1] ||
			 curImagPart != -imagParts[kk])
		  {
		    os << " is followed by shift number " << (kk+1) 
		       << " = (" << realParts[kk+1] << ", " 
		       << imagParts[kk+1] << "), that is not the complex "
		      "conjugate of the shift before it.";
		  }
		else if (curMult != multiplicities[kk+1])
		  {
		    os << " has multiplicity " << curMult << ", but is "
		      "followed by shift " << (kk+1) << " = (" 
		       << realParts[kk+1] << ", " << imagParts[kk+1] << ") "
		      "with multiplicity " << multiplicities[kk+1] << ".";
		  }
		throw std::invalid_argument(os.str());
	      }
	    // Handle both elements of the complex conjugate pair
	    // at the same time.
	    for (int kkk = 0; kkk < mults[kk]; ++kkk)
	      {
		(*B)(j, j) = curRealPart;
		(*B)(j+1, j) = STS::one();
		++j;
		(*B)(j-1, j) = -(curImagPart * curImagPart);
		(*B)(j, j) = curRealPart;
		(*B)(j+1, j) = STS::one();
		++j;
	      }
	    // We've already used lambda(kk+1) implicitly, so skip ahead.
	    kk = kk + 2;
	  }
	else if (imagParts[kk] < 0)
	  {
	    std::ostringstream os;
	    os << "With the modified Newton basis, any complex shifts must "
	      "occur in consecutive complex conjugate pairs, with the first "
	      "member of each pair having positive imaginary part.  However, "
	      "shift " << kk << " = (" << realParts[kk] << ", " 
	       << imagParts[kk] << ") begins a complex conjugate pair, but "
	      "has negative imaginary part.";
	    throw std::invalid_argument(os.str());
	  }
	else // imagParts[kk] == 0 
	  {
	    for (int kkk = 0; kkk < curMult; ++kkk)
	      { // This is a real shift.  Hopefully Scalar has an operator=
		// that does the right thing for magnitude_type inputs.
		(*B)(j, j) = curRealPart;
		(*B)(j+1, j) = STS::one();
		++j;
	      }
	    kk = kk + 1;
	  }
      }
    return B;
  }


  // Specialization for real Scalar
  template<class Scalar>
  Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,Scalar> >
  makeNewtonChangeOfBasisMatrix<Scalar, false> (const Teuchos::ArrayView<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& realParts,
						const Teuchos::ArrayView<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& imagParts,
						const Teuchos::ArrayView<const int>& multiplicities,
						const bool modified)
  {
    using Teuchos::RCP;
    typedef Teuchos::SerialDenseMatrix<int,Scalar> mat_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType magnitude_type;

    RCP<mat_type> B (new mat_type (s+1, s));
    if (modified)
      makeModifiedNewtonChangeOfBasisMatrix<Scalar> (realParts, imagParts, multiplicities);
    else
      {
	typedef typename Teuchos::ArrayView<const magnitude_type>::size_type size_type;
	int j = 0; // Current column index of B
	// For each (multiple) shift value
	for (size_type k = 0; k < multiplicities.size(); ++k)
	  {
	    TEST_FOR_EXCEPTION(imagParts[k] != 0, std::invalid_argument,
			       "The (non-modified) Newton basis does not work when "
			       "Scalar is a real type, but the shifts are complex. "
			       " Shift number " << k << " = (" << realParts[k] 
			       << ", " << imagParts[k] << ").");
	    // For each multiplicity of that shift value
	    for (size_type kk = 0; kk < multiplicities[k]; ++kk)
	      {
		(*B)(j, j) = realParts[k];
		(*B)(j+1, j) = STS::one();
		++j;
	      }
	  }
      }
    return B;
  }


  // Specialization for complex Scalar
  template<class Scalar>
  Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,Scalar> >
  makeNewtonChangeOfBasisMatrix<Scalar, true> (const Teuchos::ArrayView<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& realParts,
					       const Teuchos::ArrayView<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& imagParts,
					       const Teuchos::ArrayView<const int>& multiplicities,
					       const bool modified)
  {
    using Teuchos::RCP;
    typedef Teuchos::SerialDenseMatrix<int,Scalar> mat_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType magnitude_type;

    RCP<mat_type> B (new mat_type (s+1, s));

    if (modified)
      return makeModifiedNewtonChangeOfBasisMatrix<Scalar> (realParts, imagParts, multiplicities);
    else
      {
	typedef typename Teuchos::ArrayView<const magnitude_type>::size_type size_type;
	int j = 0; // Current column index of B
	for (size_type k = 0; k < multiplicities.size(); ++k)
	  for (size_type kk = 0; kk < multiplicities[k]; ++kk)
	    {
	      // We have to template this function on whether or not
	      // Scalar is a complex-valued type, because the
	      // two-argument constructor for Scalar only works
	      // (syntactically) if Scalar is complex.
	      (*B)(j, j) = Scalar (realParts[k], imagParts[k]);
	      (*B)(j+1, j) = STS::one();
	      ++j;
	    }
      }
    return B;
  }

  /// \class MonomialOpAkx
  /// \brief Newton-basis Akx implementation using abstract operators.
  /// \author Mark Hoemmen
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

      TEST_FOR_EXCEPTION(s < 0, std::invalid_argument, 
			 "Number of basis vectors to compute (s = " 
			 << s << ") is invalid; it must be positive.");
      if (s == 0)
	return; // Nothing to do
      TEST_FOR_EXCEPTION(MVT::GetNumberVecs(V_cur) < s, std::invalid_argument,
			 "You're asking to compute s = " << s << " basis "
			 "vectors, but only have room for " 
			 << MVT::GetNumberVecs(V_cur) << ".");
      RCP<const MV> v_prv = Teuchos::rcpFromRef (q_last);
      int k = 0;
      while (k < s)
	{
	  RCP<MV> v_cur = MVT::CloneViewNonConst (V_cur, Range1D(k,k));
	  if (imagParts_[k] == STS::zero())
	    {
	      applyOp (v_prv, v_cur);
	      MVT::MvAddMv (-realParts_[k], v_prv, STS::one(), v_cur, v_cur);
	      v_prv = rcp_const_cast<const MV> (v_cur);
	      ++k;
	    }
	  // "Modified" Newton basis trick to ensure real arithmetic,
	  // even if the shifts are complex (as long as any complex
	  // shifts occur in consecutive complex conjugate pairs with
	  // the positive-real-part member of the pair first).
	  //
	  // (A - \overline{\theta} I) (A - \theta) x
	  // = A^2 x - (\overline{\theta} + \theta) A x + |\theta|^2 x
	  // = A (A x - (\overline{\theta} + \theta) x) + |\theta|^2 x.
	  else
	    {
	      applyOp (v_prv, v_cur);
	      MVT::MvAddMv (-realParts_[k], v_prv, STS::one(), v_cur, v_cur);
	      v_prv = rcp_const_cast<const MV> (v_cur);
	      v_cur = MVT::CloneViewNonConst (V_cur, Range1D(k+1,k+1));
	      applyOp (v_prv, v_cur);
	      {
		const magnitude_type theMag = 
		  carefulMagnitude (realParts[k], imagParts[k]);
		// FIXME (mfh 08 Feb 2011) What if squaring the
		// magnitude overflows?  Then we might want to use a
		// skipahead strategy.  If we want reproducible
		// results, we should then remember whether the
		// magnitude overflowed.
		MVT::MvAddMv (theMag*theMag, v_prv, STS::one(), v_cur, v_cur);
	      }
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

      TEST_FOR_EXCEPTION(s < 0, std::invalid_argument, 
			 "Number of basis vectors to compute (s = " 
			 << s << ") is invalid; it must be positive.");
      if (s == 0)
	return; // Nothing to do
      TEST_FOR_EXCEPTION(MVT::GetNumberVecs(V_cur) < s, std::invalid_argument,
			 "You're asking to compute s = " << s << " basis "
			 "vectors, but only have room for " 
			 << MVT::GetNumberVecs(V_cur) << " in V_cur.");
      TEST_FOR_EXCEPTION(MVT::GetNumberVecs(Z_cur) < s, std::invalid_argument,
			 "You're asking to compute s = " << s << " basis "
			 "vectors, but only have room for " 
			 << MVT::GetNumberVecs(Z_cur) << " in Z_cur.");
      RCP<const MV> v_prv = Teuchos::rcpFromRef (q_last);
      for (int k = 0; k < s; ++k, v_prv)
	{
	  RCP<MV> z_cur = MVT::CloneViewNonConst (Z_cur, Range1D(k,k));
	  RCP<MV> v_cur = MVT::CloneViewNonConst (V_cur, Range1D(k,k));
	  applyFlexibleOp (v_prv, z_cur, v_cur);
	  v_prv = rcp_const_cast<const MV> (v_cur);
	}
    }

    static void
    checkShifts (const Teuchos::ArrayView<const magnitude_type>& realParts,
		 const Teuchos::ArrayView<const magnitude_type>& imagParts,
		 const Teuchos::ArrayView<const int>& mults)
    {
      using Teuchos::ArrayView;
      typedef Teuchos::ScalarTraits<Scalar> STS;
      typedef Teuchos::ScalarTraits<magnitude_type> STM;
      
      const typename Teuchos::ArrayView<const int>::size_type size_type;
      TEST_FOR_EXCEPTION(mults.size() != realParts.size() || mults.size() != imagParts.size(),
			 std::logic_error, "Number of multiplicities is not "
			 "commensurate with the number of real or imaginary "
			 "parts of the shifts.");
      TEST_FOR_EXCEPTION(realParts.size() != imagParts.size(), std::logic_error,
			 "Number of real parts of the shifts is different "
			 "than the number of imaginary parts.");
      int k = 0;
      while (k < mults.size())
	{
	  if (imagParts_[k] > STM::zero())
	    {
	      TEST_FOR_EXCEPTION(k == multiplicities_.size() - 1, std::logic_error, "When computing the modified Newton basis' change-of-basis matrix: the last shift is the first member of a complex conjugate pair, with no second member.");
	      TEST_FOR_EXCEPTION(realParts_[k] != realParts_[k+1], std::logic_error, "When computing the modified Newton basis' change-of-basis matrix: complex shifts may only occur in complex conjugate pairs.");
	      TEST_FOR_EXCEPTION(imagParts_[k] != -imagParts_[k+1], std::logic_error, "When computing the modified Newton basis' change-of-basis matrix: complex shifts may only occur in complex conjugate pairs.");
	      TEST_FOR_EXCEPTION(multiplicities[k] != multiplicities_[k+1], std::logic_error, "When computing the modified Newton basis' change-of-basis matrix: in a complex conjugate pair, the two members of the pair have different multiplicities.");
	      k = k + 2;
	    }
	  else if (imagParts_[k] < STM::zero())
	    {
	      throw std::logic_error("Complex shift occurred first with negative imaginary part.");
	    }
	  else // imagParts_[k] == 0
	    ++k;
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

      // const typename Teuchos::ArrayView<const int>::size_type size_type;
      // const int totalNumShifts = 
      // 	std::accumulate (multiplicities.begin(), multiplicities.end(), int(0));
      checkShifts (realParts_, imagParts_, multiplicities_);

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

      // const typename Teuchos::ArrayView<const int>::size_type size_type;
      // const int totalNumShifts = 
      // 	std::accumulate (multiplicities.begin(), multiplicities.end(), int(0));
      checkShifts (realParts_, imagParts_, multiplicities_);

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


  private:

    static Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,Scalar> > 
    makeChangeOfBasisMatrix (const Teuchos::ArrayView<magnitude_type>& realParts,
			     const Teuchos::ArrayView<magnitude_type>& imagParts,
			     const Teuchos::ArrayView<const int>& multiplicities,
			     const bool modified)
    {
      TEST_FOR_EXCEPTION(realParts.size() != imagParts.size(), 
			 std::logic_error, 
			 "The real parts and imaginary parts arrays storing the"
			 " Newton basis shifts should always be the same length"
			 ", but here realParts.size() = " << realParts.size() 
			 << " and imagParts.size() = " << imagParts.size() 
			 << ".");
      TEST_FOR_EXCEPTION(realParts.size() != multiplicities.size(), 
			 std::logic_error, 
			 "The number of Newton basis shifts (modulo multiplici"
			 "ties) " << realParts.size() << " does not match the "
			 "number of multiplicities " << multiplicities.size() 
			 << ".");
      typedef Teuchos::ScalarTraits<Scalar> STS;
      return makeNewtonChangeOfBasisMatrix<Scalar, STS::isComplex> (realParts, imagParts, multiplicities, modified);
    }

  template<class Scalar, bool isComplex>
  Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,Scalar> >
  makeNewtonChangeOfBasisMatrix (const Teuchos::ArrayView<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& realParts,
				 const Teuchos::ArrayView<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& imagParts,
				 const Teuchos::ArrayView<const int>& multiplicities,
				 const bool modified);


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
