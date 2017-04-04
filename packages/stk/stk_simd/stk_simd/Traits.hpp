// Copyright 2013 Sandia Corporation, Albuquerque, NM.

#ifndef STK_SIMD_TRAITS_H
#define STK_SIMD_TRAITS_H

#include <cmath>

namespace stk {

namespace {
    
// length

template <typename FloatType>
struct SimdLength {
  static constexpr int length = simd::ndoubles;
};

template <>
struct SimdLength<double> {
  static constexpr int length = 1;
};

template <>
struct SimdLength<float> {
  static constexpr int length = 1;
};

#if defined(STK_SIMD)
template <>
struct SimdLength<simd::Float> {
  static constexpr int length = simd::nfloats;
};
#endif

// bools (get bool type from real type)

template <typename T>
struct BoolT {
  typedef simd::Bool type;
};

template <>
struct BoolT<double> {
  typedef bool type;
};

template <>
struct BoolT<float> {
  typedef bool type;
};

template <>
struct BoolT<bool> {
  typedef bool type;
};

#if defined(STK_SIMD)
template <>
struct BoolT<simd::Float> {
  typedef simd::Boolf type;
};

template <>
struct BoolT<simd::Boolf> {
  typedef simd::Boolf type;
};
#endif

// // reals (get real type from bool type)

template <typename T>
struct RealT {
  typedef simd::Double type;
};

template <>
struct RealT<float> {
  typedef float type;
};

template <>
struct RealT<bool> {
  typedef double type;
};

#if defined(STK_SIMD)
template <>
struct RealT<simd::Float> {
  typedef simd::Float type;
};

template <>
struct RealT<simd::Boolf> {
  typedef simd::Float type;
};

template <>
struct RealT<double> {
  typedef double type;
};
#endif

// simd type from base type

template <typename T>
struct SimdT {
  typedef T type;
};
   
template <>
struct SimdT<double> {
  typedef simd::Double type;
};
 
template <>
struct SimdT<float> {
  typedef simd::Float type;
};

template <>
struct SimdT<bool> {
  typedef simd::Bool type;
};

template <>
struct SimdT<double&> {
  typedef simd::Double& type;
};
 
template <>
struct SimdT<float&> {
  typedef simd::Float& type;
};

template <>
struct SimdT<bool&> {
  typedef simd::Bool& type;
};

template <>
struct SimdT<const double&> {
  typedef const simd::Double& type;
};
 
template <>
struct SimdT<const float&> {
  typedef const simd::Float& type;
};

template <>
struct SimdT<const bool&> {
  typedef const simd::Bool& type;
};
   
// base type from simd type

template <typename T>
struct BaseT {
  typedef T type;
};
   
#if defined(STK_SIMD)
template <>
struct BaseT<simd::Double> {
  typedef double type;
};

template <>
struct BaseT<simd::Float> {
  typedef float type;
};

template <>
struct BaseT<simd::Boolf> {
  typedef bool type;
};

template <>
struct BaseT<simd::Bool> {
  typedef bool type;
};
#endif

}

template <typename Type>
struct Traits {
  // convert between real and bool types
  typedef typename RealT<Type>::type real_type;
  typedef typename BoolT<Type>::type bool_type;

  // convert between simd and base types
  typedef typename SimdT<Type>::type simd_type;
  typedef typename BaseT<Type>::type base_type;
 
  static const real_type ZERO;
  static const real_type SIGN_MASK;
  static const real_type ONE;
  static const real_type NEG_ONE;
  static const real_type TWO;
  static const real_type THREE;
  static const real_type HALF;
  static const real_type THIRD;
  static const real_type SQRT_TWO;
  static const real_type SQRT_THREE;

  static const bool_type FALSE_VAL;
  static const bool_type TRUE_VAL;

  static constexpr int length = SimdLength<Type>::length;
};

template<typename Type> const typename Traits<Type>::real_type Traits<Type>::ZERO(0.0);
template<typename Type> const typename Traits<Type>::real_type Traits<Type>::SIGN_MASK(-0.0);
template<typename Type> const typename Traits<Type>::real_type Traits<Type>::ONE(1.0);
template<typename Type> const typename Traits<Type>::real_type Traits<Type>::NEG_ONE(-1.0);
template<typename Type> const typename Traits<Type>::real_type Traits<Type>::TWO(2.0);
template<typename Type> const typename Traits<Type>::real_type Traits<Type>::THREE(3.0);
template<typename Type> const typename Traits<Type>::real_type Traits<Type>::HALF(0.5);
template<typename Type> const typename Traits<Type>::real_type Traits<Type>::THIRD(1.0/3.0);
template<typename Type> const typename Traits<Type>::real_type Traits<Type>::SQRT_TWO(std::sqrt(2.0));
template<typename Type> const typename Traits<Type>::real_type Traits<Type>::SQRT_THREE(std::sqrt(3.0));

template<typename Type> const typename Traits<Type>::bool_type Traits<Type>::FALSE_VAL(false);
template<typename Type> const typename Traits<Type>::bool_type Traits<Type>::TRUE_VAL(true);

} // namespace stk

#endif
