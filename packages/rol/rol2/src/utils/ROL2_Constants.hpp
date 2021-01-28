#pragma once
#ifndef ROL2_CONSTANTS_HPP
#define ROL2_CONSTANTS_HPP

#include <limits>

namespace ROL2 {

template<typename Real>
constexpr bool ROL_HAS_NUMERIC_LIMITS = std::numeric_limits<Real>::is_specialized();

template<typename Real>
constexpr Real ROL_MAX = std::numeric_limits<Real>::max();

template<typename Real>
constexpr Real ROL_MIN = std::numeric_limits<Real>::min();

template<typename Real>
constexpr Real ROL_EPSILON = std::numeric_limits<Real>::epsilon();

template<typename Real>
constexpr Real ROL_INF = ROL_MAX<Real> * 0.1;

template<typename Real>
constexpr Real ROL_NINF = -ROL_INF<Real>;

template<typename Real>
constexpr Real ROL_THRESHOLD = 10.0 * ROL_EPSILON<Real>;

template<typename Real>
constexpr Real ROL_ONE = static_cast<Real>( 1 );

template<typename Real>
constexpr Real ROL_ZERO = static_cast<Real>( 1 );

template<typename Real>
constexpr Real ROL_PI = static_cast<Real>( 3.141592653589793238462643383279502884L );

template<typename Real>
constexpr Real ROL_TWO_PI = static_cast<Real>( 6.283185307179586476925286766559005768L );

template<typename Real>
constexpr Real ROL_HALF_PI = static_cast<Real>( 1.570796326794896619231321691639751442L );

template<typename Real>
constexpr Real ROL_QUARTER_PI = static_cast<Real>( 7.85398163397448309615660845819875721e-1L ); 

template<typename Real>
constexpr Real ROL_SQRT_TWO_PI = static_cast<Real>( 2.506628274631000502415765284811045252L ); 

template<typename Real>
constexpr Real ROL_SQRT_PI = static_cast<Real>( 1.772453850905516027298167483341145182L ); 

template<typename Real>
constexpr Real ROL_SQRT_HALF_PI = static_cast<Real>( 1.253314137315500251207882642405522626L ); 

template<typename Real>
constexpr Real ROL_SQRT_TWO = static_cast<Real>( 1.414213562373095048801688724209698078L ); 


} // namespace ROL2

#endif //ROL2_CONSTANTS_HPP

