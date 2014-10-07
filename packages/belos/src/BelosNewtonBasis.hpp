#ifndef __Belos_NewtonBasis_hpp
#define __Belos_NewtonBasis_hpp

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

namespace Belos {

  /// Returning the ostream& lets you "chain" outputs:
  ///
  /// out << "foo: " << printShift(a,b) << ", bar" << endl;
  template<class Scalar>
  static std::ostream&
  printShift (std::ostream& out, 
	      const typename Teuchos::ScalarTraits<Scalar>::magnitudeType& realPart,
	      const typename Teuchos::ScalarTraits<Scalar>::magnitudeType& imagPart)
  {
    out << "(" << realPart << ", " << imagPart << ")";
    return out;
  }

  template<class Scalar, class ExceptionType>
  void
  checkModifiedNewtonShifts (const Teuchos::ArrayView<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& realParts,
	       const Teuchos::ArrayView<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& imagParts,
	       const Teuchos::ArrayView<const int>& mults)
  {
    using Teuchos::ArrayView;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    typedef Teuchos::ScalarTraits<magnitude_type> STM;
    typedef Teuchos::ArrayView<const int>::size_type size_type;

    TEUCHOS_TEST_FOR_EXCEPTION(realParts.size() != imagParts.size(), 
		       ExceptionType, 
		       "Newton basis shifts: There are " << realParts.size() 
		       << " real values, but " << imagParts.size() 
		       << " corresponding imaginary values.  realParts[k] + "
		       "i*imagParts[k] comprises a single complex-valued shift,"
		       " so the realParts and imagParts arrays should always "
		       "have the same length.";
    TEUCHOS_TEST_FOR_EXCEPTION(mults.size() != realParts.size() || 
		       mults.size() != imagParts.size(),
		       ExceptionType, 
		       "Newton basis shifts: There are " << realParts.size()
		       << " unique shift values, but " << mults.size()
		       << " multiplicities.  The realParts, imagParts, and "
		       "mults arrays should all have the same length.";
      int k = 0;
      while (k < mults.size())
	{
	  if (mults[k] <= 0)
	    {
	      std::ostringstream os;
	      os << "Error in shifts for the modified Newton basis: "
		 << "Shift number " << (k+1) << " of " << mults.size()
		 << " has nonpositive multiplicity " << mults[k] 
		 << ".  All shifts for the Newton basis must have "
		"positive multiplicity.";
	      throw ExceptionType (os.str());
	    }
	  else if (imagParts[k] > STM::zero())
	    {
	      const bool bad = k == multiplicities_.size() - 1 ||
		realParts[k] != realParts[k+1] ||
		imagParts[k] != -imagParts[k+1] ||
		multiplicities[k] != multiplicities_[k+1];
	      if (bad)
		{
		  std::ostringstream os;
		  os << "Error in shifts for the modified Newton basis: ";
		  if (multiplicities_.size() - 1)
		    os << "The last (of " << mults.size() << ") shift " 
		       << printShift(realParts[k], imagParts[k])
		       << " has nonzero imaginary part.  "
		      "This is forbidden, since complex shifts must occur in "
		      "consecutive complex conjugate pairs in the modified "
		      "Newton basis, but this complex shift occurs on its own.";
		  else if (realParts[k] != realParts[k+1] ||
			   imagParts[k] != -imagParts[k+1])
		    os << "The current (" << (k+1) << " of " << mults.size() 
		       << ") shift "
		       << printShift(realParts[k], imagParts[k])
		       << " is complex, but the shift that follows "
		       << printShift(realParts[k+1], imagParts[k+1])
		       << " is not its complex conjugate.  This is forbidden, "
		      "since shifts must occur in consecutive complex "
		      "conjugate pairs.";
		  else if (multiplicities[k] != multiplicities_[k+1])
		    os << "Shifts number " << (k+1) << " and " << (k+2) 
		       << " are consecutive shifts that form a complex "
		      "conjugate pair, but they have different "
		      "multiplicities.  This is forbidden in the modified "
		      "Newton basis, since complex shifts must be the "
		      "eigenvalues of a real matrix, and thus may only occur "
		      "in complex conjugate pairs with the same multiplicity.";
		  throw ExceptionType (os.str());
		}
	      k = k + 2;
	    }
	  else if (imagParts[k] < STM::zero())
	    {
	      std::ostringstream os;
	      os << "Error in shifts for the modified Newton basis: "
		 << "Shift number " << (k+1) << " of " << mults.size()
		 << " is complex with negative imaginary part.  Complex "
		"shifts are allowed, but must occur in complex conjugate "
		"pairs, with the first member of each pair having positive "
		"imaginary part.  In this case, the other member of the pair "
		"should have occurred first.";
	      throw ExceptionType (os.str());
	    }
	  else // imagParts[k] == 0
	    ++k;
	}
    }

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
	//
	// We assume that the given shifts have been run through
	// checkModifiedNewtonShifts(), which will detect any of the possible error
	// conditions.
	if (imagParts[kk] > 0)
	  {
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
	else if (imagParts[kk] == 0)
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
	else
	  throw std::logic_error("Should never get here!");
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
	    TEUCHOS_TEST_FOR_EXCEPTION(imagParts[k] != 0, std::invalid_argument,
			       "The (non-modified) Newton basis does not work when "
			       "Scalar is a real type, but the shifts are complex. "
			       " Shift number " << k << " = " 
			       << printShift(realParts[k], imagParts[k]) << ".";
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

} // namespace Belos

#endif // __Belos_NewtonBasis_hpp
