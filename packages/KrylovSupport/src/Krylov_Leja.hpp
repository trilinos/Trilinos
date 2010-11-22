#ifndef __Krylov_Leja_hpp
#define __Krylov_Leja_hpp

#include <Krylov_ProductIsZero.hpp>
#include <Krylov_ShiftsOutOfOrder.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace Krylov {

  /// \brief Return absolute value of x
  ///
  /// Return the absolute value of x, taking care that the result
  /// makes sense if x is +/-Inf or NaN.
  ///
  /// \note If x is -Inf, this will return +Inf.  If x is NaN, this
  ///   will return -NaN which is just NaN.
  template< class MagnitudeType >
  MagnitudeType carefulAbs (const MagnitudeType& x)
  {
    if (x < Teuchos::ScalarTraits< MagnitudeType >::zero())
      return -x;
    else
      return x;
  }

  /// \brief \fn$\sqrt(x^2 + y^2)\fn$ with careful scaling
  ///
  /// Compute \fn$\sqrt(x^2 + y^2)\fn$, carefully scaling to avoid
  /// unnecessary underflow.  The result is correct even if one or
  /// both of (x,y) are +/-Inf or NaN.
  ///
  /// \note If either of (x, y) is NaN, then this method will return
  ///   NaN, as it should.  Also, if either of (x, y) is +/-Inf but
  ///   neither is NaN, then this methods will return +Inf, as it
  ///   should.
  ///
  /// \note This method compares to LAPACK's _LAPY2 routine.
  ///
  template< class MagnitudeType >
  MagnitudeType LAPY2 (const MagnitudeType& x, const MagnitudeType& y)
  {
    using Teuchos::ScalarTraits;

    if (ScalarTraits< MagnitudeType >::isnaninf (x) ||
	ScalarTraits< MagnitudeType >::isnaninf (y))
      // This makes sense if you use a case-based analysis.
      //
      // Case 1: Either of x or y is NaN: then return NaN.  (NaN +
      //   anything is NaN.)  carefulAbs will carry NaN inputs
      //   through.
      //
      // Case 2: Either of x or y is +/-Inf, but either is NaN: then
      //   return +Inf.  (+Inf + (non-NaN) == +Inf.)
      //   carefulAbs(+/-Inf) == +Inf.
      return carefulAbs (x) + carefulAbs (y);

    // If x or y is NaN, std::max() will be incorrect.  It may be
    // wrong in different ways, depending on its implementation.
    // That's why we check for NaN above.
    const MagnitudeType x_abs = ScalarTraits< MagnitudeType >::magnitude (x);
    const MagnitudeType y_abs = ScalarTraits< MagnitudeType >::magnitude (y);
    const MagnitudeType maxAbs = std::max (x_abs, y_abs);

    if (maxAbs == ScalarTraits< MagnitudeType >::zero())
      return maxAbs; // Don't divide by zero
    else
      {
	const MagnitudeType x_scaled = x_abs / maxAbs;
	const MagnitudeType y_scaled = y_abs / maxAbs;
	return maxAbs * ScalarTraits< MagnitudeType >::squareroot(x_scaled * x_scaled + y_scaled * y_scaled);
      }
  }

  /// \class Leja
  /// \brief Compute the (modified) Leja ordering
  ///
  template< class Scalar >
  class Leja {
  public:
    typedef Scalar scalar_type;
    typedef magnitude_type Teuchos::ScalarTraits< scalar_type >::magnitudeType magnitude_type;
    typedef std::pair< magnitude_type, magnitude_type > shift_type;

    /// Compute the (modified) Leja ordering of the given shifts
    ///
    /// \param inputShifts [in] Array of unique possibly
    ///   complex-valued input shifts, each represented as an (real,
    ///   imaginary) pair.  If the modified Leja ordering is to be
    ///   computed, the shifts must be arranged such that complex
    ///   conjugate pairs appear consecutively with the eigenvalue
    ///   having the positive imaginary part occurring first in the
    ///   pair.  (LAPACK's SGEEV and DGEEV nonsymmetric dense
    ///   eigenvalue routines order their output in this way, so you
    ///   don't have to worry about arranging them that way yourself.)
    ///
    /// \param inputMultiplicities [in] Array with the same length as
    ///   inputShifts, such that inputShifts[k] is a shift that occurs
    ///   inputMultiplicities[k] times in the input set of shifts.
    ///   (The elements of inputShifts must be unique, but we allow
    ///   repeated shifts via the multiplicities mechanism.)
    ///
    /// \param outputShifts [out] Array of the same shifts (FIXME
    ///   modulo a tiny bit of rounding error due to scaling -- fix by
    ///   using outputIndices) ordered using the Leja scheme (or the
    ///   modified Leja scheme, if modified is true).
    /// 
    /// \param outputMultiplicities [out] The input multipliticies,
    ///   reordered using the same ordering as outputShifts, for your
    ///   convenience.
    ///
    /// \param numShifts [in] Total number of unique shifts (not
    ///   counting multiplicities).  Only this many of the first
    ///   consecutive elements of inputShifts and inputMultiplicities
    ///   are read.
    ///
    /// \param modified [in] Whether to compute the modified Leja
    /// ordering.  This affects the required ordering of the shifts on
    /// input.
    static void
    orderShifts (const std::vector< shift_type >& inputShifts, 
		 const std::vector< int >& inputMultiplicities,
		 std::vector< shift_type >& outputShifts, 
		 std::vector< int >& outputIndices,
		 std::vector< int >& outputMultiplicities,
		 const int numShifts, 
		 const bool modified);

    static bool
    shiftsValid (const std::vector< shift_type >& inputShifts,
		 const std::vector< int >& inputMultiplicities,
		 const int numShifts, 
		 const bool modified);

  private:
    static magnitude_type magnitude (const shift_type& shift) {
      // LAPACK routine for computing sqrt(x^2 + y^2), with careful
      // scaling to avoid unnecessary overflow.
      return LAPY2< magnitude_type > (shift.first, shift.second);
    }
    static magnitude_type diff (const shift_type& x, const shift_type& y) {
      return shift_type (x.first - y.first, x.second - y.second);
    }

    static int
    addOneShift (const std::vector< shift_type >& inputShifts,
		 const std::vector< int >& inputMultiplicies,
		 std::vector< shift_type >& outputShifts,
		 std::vector< int >& outputIndices,
		 std::vector< int >& outputMultiplicies,
		 const int numShifts,
		 const int maxIndex,
		 const int curPos);

    static int
    addShift (const std::vector< shift_type >& inputShifts,
	      const std::vector< int >& inputMultiplicies,
	      std::vector< shift_type >& outputShifts,
	      std::vector< int >& outputIndices,
	      std::vector< int >& outputMultiplicies,
	      const int numShifts,
	      const int maxIndex,
	      const int curPos, 
	      const bool modified);

    static int
    tableUpdate (std::vector< magnitude_type >& products, 
		 const std::vector< shift_type >& inputShifts,
		 const std::vector< int >& multiplicities,
		 const std::vector< shift_type >& outputShifts,
		 const int curPos);

    static int
    firstOutputShift (const std::vector< shift_type >& inputShifts,
		      const std::vector< int >& inputMultiplicities);

    static magnitude_type 
    scaleAndUpdateCapacityEstimate (const magnitude_type& C, 
				    std::vector< shift_type >& scaledinputShifts,
				    std::vector< shift_type >& outputShifts,
				    const std::vector< int >& outputMultiplicities,
				    std::vector< magnitude_type >& products,
				    const int curPos)
    {
      const magnitude_type C_new = magnitude_type(1);
      for (int j = 0; j < curPos; ++j)
	{
	  const magnitude_type curMag = magnitude (diff (outputShifts[curPos], outputShifts[j]));
	  const magnitude_type pwr = magnitude_type (outputMultiplicities[j]) / magnitude_type (curPos);
	  C_new = C_new * Teuchos::ScalarTraits< magnitude_type >::pow (curMag, pwr);
	}
      if (C_new == Teuchos::ScalarTraits< magnitude_type >::zero())
	throw CapacityEstimateIsZero (outputShifts, outputMultiplicities, C, curPos);

      //
      // Rescale the input and output shifts.
      //
      const magnitude_type divisor = C_new / C;
      for (int i = 0; i < numShifts; ++i)
	scaledInputShifts[i] = scaledInputShifts[i] / divisor;
      for (int j = 0; j < curPos; ++j)
	outputShifts[j] = outputShifts[j] / divisor;
      for (int j = 0; i < curPos; ++j)
	products[j] = products[j] / divisor;

      return C_new;
    }
  };

  template< class Scalar >
  int
  Leja< Scalar >::
  addOneShift (const std::vector< typename Leja< Scalar >::shift_type >& inputShifts,
	       const std::vector< int >& inputMultiplicies,
	       std::vector< typename Leja< Scalar >::shift_type >& outputShifts,
	       std::vector< int >& outputIndices,
	       std::vector< int >& outputMultiplicies,
	       const int numShifts,
	       const int maxIndex,
	       const int curPos)
  {
    if (maxIndex < 0 || maxIndex >= numShifts)
      throw ShiftsOutOfOrder (maxIndex, curPos);
    else if (curPos >= numShifts) // not enough room to add the new shift
      throw ShiftsOutOfOrder (maxIndex, curPos);
    outputShifts[curPos] = inputShifts[maxIndex];
    outputIndices[curPos] = maxIndex;
    outputMultiplicies[maxIndex] = inputMultiplicies[curPos];
    return curPos + 1;
  }


  template< class Scalar >
  int
  Leja< Scalar >::
  addShift (const std::vector< typename Leja< Scalar >::shift_type >& inputShifts,
	    const std::vector< int >& inputMultiplicies,
	    std::vector< typename Leja< Scalar >::shift_type >& outputShifts,
	    std::vector< int >& outputIndices,
	    std::vector< int >& outputMultiplicies,
	    const int numShifts,
	    const int maxIndex,
	    const int curPos, 
	    const bool modified)
  {
    // Always add at least one shift.
    const int nextPos = addOneShift (inputShifts, inputMultiplicities, outputShifts, outputIndices, 
				     outputMultiplicities, numShifts, maxIndex, curPos);
    if (! modified)
      return nextPos;
    else
      {
	const magnitude_type ZERO = Teuchos::ScalarTraits< magnitude_type >::zero();
	const magnitude_type imagPart = inputShifts[maxIndex].second;
	if (imagPart == ZERO)
	  return nextPos;
	else
	  {
	    // The imaginary part of the first shift in the complex
	    // conjugate pair must be positive.
	    if (imagPart < ZERO) 
	      throw ShiftsOutOfOrder (maxIndex, nextPos);
	    else
	      return addOneShift (inputShifts, inputMultiplicities, outputShifts, outputIndices, 
				  outputMultiplicities, numShifts, maxIndex+1, nextPos);
	  }
      }
  }

  template< class Scalar >
  bool
  Leja< Scalar >::
  shiftsValid (const std::vector< typename Leja< Scalar >::shift_type >& inputShifts,
	       const std::vector< int >& inputMultiplicities,
	       const int numShifts, 
	       const bool modified)
  {
    int k = 0; 
    while (k < numShifts)
      {
	// Each shift has to occur with positive multiplicity.
	if (inputMultiplicities[k] <= 0)
	  return false;

	const MagnitudeType ZERO = Teuchos::ScalarTraits< MagnitudeType >::zero();
	const MagnitudeType imagPart = inputShifts[k].second;
	// The non-modified Newton ordering doesn't care about the
	// value of the shift.  For the modified Newton ordering, real
	// shifts are allowed at any position, as long as the ordering
	// rules about complex shifts aren't violated.
	if (! modified || imagPart == ZERO)
	  k = k + 1;
	// For the modified Newton ordering, complex-valued shifts
	// must occur in consecutive complex conjugate pairs, with the
	// first of the two shifts in the pair having positive
	// imaginary part.
	else if (imagPart < ZERO)
	  return false; 
	else // imagPart > ZERO
	  {
	    // Check if the other shift in the conjugate pair is
	    // missing.
	    if (k >= numShifts - 1)
	      return false; 
	    // Check that the other shift in the conjugate pair has
	    // negative imaginary part.  In exact arithmetic, its
	    // imaginary part should be -imagPart.  However, it may
	    // have been computed subject to rounding error, so don't
	    // check for exact equality.  The caller is responsible
	    // for ensuring that the rounding error is sufficiently
	    // small.
	    else if (inputShifts[k+1].second >= 0)
	      return false;
	    // The multiplicity of the other shift in the conjugate
	    // pair should be the same as the multiplicity of the
	    // first shift in the pair.
	    else if (inputMultiplicities[k+1] != inputMultiplicities[k])
	      return false;
	    else
	      k = k + 2;
	  }
      }
    if (k == numShifts)
      return true;
    else
      return false;
  }


  template< class Scalar >
  int
  Leja< Scalar >::
  tableUpdate (std::vector< typename Leja< Scalar >::magnitude_type >& products,
	       const std::vector< typename Leja< Scalar >::shift_type >& inputShifts,
	       const std::vector< int >& inputMultiplicities,
	       const std::vector< typename Leja< Scalar >::shift_type >& outputShifts,
	       const int curPos)
  {
    using Teuchos::ScalarTraits;

    const magnitude_type ZERO = ScalarTraits< magnitude_type >::zero();
    magnitude_type maxMag (-1);
    int maxIndex = -1;
    
    for (int i = 0; i < numShifts; ++i)
      {
	// We use magnitude() and diff() in order to avoid complex
	// arithmetic, when we are not building Trilinos with complex
	// arithmetic support.  Note that this magnitude() has a
	// different interface than the magnitude() in Teuchos::
	// ScalarTraits.
	const magnitude_type curMag = magnitude (diff (inputShifts[i], outputShifts[curPos]));
	if (curMag > ZERO)
	  {
	    const magnitude_type pwr = magnitude_type (inputMultiplicities[i]);
	    products[i] = products[i] * ScalarTraits< magnitude_type >::pow (curMag, pwr);
	    // Use a strict inequality, so that we pick the minimum index
	    // (in case some values repeat, which they shouldn't, since
	    // multiplicities should take care of that).
	    if (products[i] > maxMag)
	      {
		maxMag = products[i];
		maxIndex = i;
	      }
	  }
      }
    // If all results are zero, we return -1 as a flag.
    return maxIndex;
  }

  template< class Scalar >
  int
  Leja< Scalar >::
  firstOutputShift (const std::vector< typename Leja< Scalar >::shift_type >& inputShifts,
		    const std::vector< int >& inputMultiplicities)
  {
    int maxIndex = -1;
    magnitude_type maxVal (-1);
    for (int i = 0; i < numShifts; ++i)
      {
	const magnitude_type curMag = magnitude (inputShifts[i]);
	const magnitude_type pwr = magnitude_type (inputMultiplicities[i]);
	const magnitude_type curVal = 
	  Teuchos::ScalarTraits< magnitude_type >::pow (curMag, pwr);

	// Use a strict inequality, so that we pick the minimum index
	// (in case some values repeat, which they shouldn't, since
	// multiplicities should take care of that).
	if (curVal > maxVal)
	  {
	    maxVal = curVal;
	    maxIndex = i;
	  }
      }
    // If all results are zero, we return -1 as a flag.
    return maxIndex;
  }

  template< class Scalar >
  void
  Leja< Scalar >::
  orderShifts (const std::vector< typename Leja< Scalar >::shift_type >& inputShifts, 
	       const std::vector< int >& inputMultiplicities,
	       std::vector< typename Leja< Scalar >::shift_type >& outputShifts, 
	       std::vector< int >& outputIndices,
	       std::vector< int >& outputMultiplicities,
	       const int numShifts, 
	       const bool modified)
  {
    using Teuchos::ScalarTraits;

    if (numShifts < 0)
      throw std::invalid_argument("numShifts < 0");
    else if (inputShifts.size() < numShifts)
      throw std::invalid_argument("inputShifts is too short");
    else if (inputMultiplicities.size() < numShifts)
      throw std::invalid_argument("inputMultiplicities is too short");

    outputShifts.resize (numShifts);
    outputIndices.resize (numShifts);
    outputMultiplicities.resize (numShifts);
    if (numShifts == 0) return; // Nothing to do

    // List of partial products of magnitudes of differences, used for
    // a dynamic programming algorithm for computing the Leja ordering.
    std::vector< magnitude_type > products (numShifts, magnitude_type(1));

    // Fill with flags to denote not-yet-initialized values
    std::fill (outputIndices.begin(), outputIndices.end(), -1);
    std::fill (outputMultiplicities.begin(), outputMultiplicities.end(), -1);

    // Copy input shifts, since our algorithm will rescale them in order
    // to avoid {under-,over-}flow when maximizing the product for the
    // (modified) Leja ordering.
    std::vector< shift_type > scaledInputShifts (numShifts);
    std::copy (inputShifts.begin(), inputShifts.end(), scaledInputShifts.begin());

    // Current capacity estimate
    magnitude_type C = ScalarTraits< magnitude_type >::one();
  
    // The first output shift is the shift with largest magnitude.
    int curPos = 0;
    maxIndex = firstOutputShift (scaledInputShifts, inputMultiplicities);
    if (maxIndex == -1) // flag for all zero results
      throw ProductIsZero (j, false);
    curPos = addShift (scaledInputShifts, inputMultiplicities, outputShifts, outputIndices, 
		       outputMultiplicities, numShifts, maxIndex, curPos, modified);
    //
    // Classic dynamic programming algorithm: O(numShifts^2) operations,
    // O(numShifts) storage.  
    //
    C = ScalarTraits< magnitude_type >::one();
    while (curPos < numShifts)
      {
	// Update capacity estimate.
	C = scaleAndUpdateCapacityEstimate (C, scaledInputShifts, outputShifts, outputMultiplicities, products, curPos);
	// Index of best input shift
	maxIndex = tableUpdate (products, scaledInputShifts, multiplicities, outputShifts, curPos, modified);
	if (maxIndex == -1) // flag for all zero results
	  throw ProductIsZero (curPos+1, false);
	curPos = addShift (scaledInputShifts, inputMultiplicities, outputShifts, outputIndices, 
			   outputMultiplicities, numShifts, maxIndex, curPos, modified);
      }

    // Undo the scaling of the output shifts.  Input shifts don't really
    // matter, since the caller should have made a copy of them anyway!
    for (int i = 0; i < numShifts; ++i)
      outputShifts[i] = outputShifts[i] * C;
  }

} // namespace Krylov
#endif // __Krylov_Leja_hpp
