// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
//
// ***********************************************************************
//
// The forward-mode AD classes in Sacado are a derivative work of the
// expression template classes in the Fad package by Nicolas Di Cesare.
// The following banner is included in the original Fad source code:
//
// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//            CEMRACS 98 : C++ courses,
//         templates : new C++ techniques
//            for scientific computing
//
//********************************************************
//
//  NumericalTraits class to illustrate TRAITS
//
//********************************************************
// @HEADER

#ifndef SACADO_TRAITS_HPP
#define SACADO_TRAITS_HPP

#include <string>
#include <type_traits>

#include "Sacado_ConfigDefs.h"
#include "Sacado_dummy_arg.hpp"
#include "Sacado_mpl_enable_if.hpp"
#include "Sacado_mpl_disable_if.hpp"

#ifdef HAVE_SACADO_COMPLEX
#include <complex>
#endif

namespace Sacado {

  /*!
   * \brief Enum use to signal whether the derivative array should be
   * initialized in AD object constructors.
   */
  enum DerivInit {
    NoInitDerivArray = 0, //!< Do not initialize the derivative array
    InitDerivArray        //!< Initialize the derivative array
  };

  //! Is a type an expression
  template <typename T>
  struct IsExpr {
    static const bool value = false;
  };

  //! Determine whether a given type is a view
  template <typename T>
  struct IsView {
    static const bool value = false;
  };

  //! Get the base Fad type from a view/expression
  template <typename T>
  struct BaseExprType {
    typedef T type;
  };

  //! Get view type for any Fad type
  template <typename T,unsigned,unsigned> struct ViewFadType {};

  //! Specialize this for a given type T to disable default Promote rules
  template <typename T> struct OverrideDefaultPromote {
    static const bool value = false;
  };

  //! Base template specification for %Promote
  /*!
   * The %Promote classes provide a mechanism for computing the
   * promoted type of a binary operation.
   */
  template <typename A, typename B, typename Enabled = void> struct Promote {};

  //! Specialization of %Promote for a single type
  template <typename A>
  struct Promote< A, A,
                  typename mpl::enable_if_c< !OverrideDefaultPromote<A>::value >::type > {
    typedef typename BaseExprType<A>::type type;
  };

  //! Specialization of %Promote when A is convertible to B but not vice-versa
  template <typename A, typename B>
  struct Promote< A, B,
                  typename mpl::enable_if_c< std::is_convertible<A,B>::value &&
                                            !std::is_convertible<B,A>::value &&
                                            !OverrideDefaultPromote<A>::value &&
                                            !OverrideDefaultPromote<B>::value
                                           >::type > {
    typedef typename BaseExprType<B>::type type;
  };

  //! Specialization of %Promote when B is convertible to A but not vice-versa
  template <typename A, typename B>
  struct Promote< A, B,
                  typename mpl::enable_if_c< std::is_convertible<B,A>::value &&
                                            !std::is_convertible<A,B>::value &&
                                            !OverrideDefaultPromote<A>::value &&
                                            !OverrideDefaultPromote<B>::value
                                           >::type > {
    typedef typename BaseExprType<A>::type type;
  };

 /*!
  * \brief Specialization of Promote when A and B are convertible to each
  * other, and one of them is an expression.
  */
  template <typename A, typename B>
  struct Promote< A, B,
                  typename mpl::enable_if_c< std::is_convertible<A,B>::value &&
                                             std::is_convertible<B,A>::value &&
                                             !std::is_same<A,B>::value &&
                                             ( IsExpr<A>::value ||
                                               IsExpr<B>::value ) >::type >
  {
    typedef typename BaseExprType<A>::type A_base_fad_type;
    typedef typename BaseExprType<B>::type B_base_fad_type;
    typedef typename Promote< A_base_fad_type, B_base_fad_type >::type type;
  };

  /*!
   * \brief Specialization of Promote when A is an expression and B is
   * convertible to its value-type, e.g., Promote< fad-expression, double >
   * (using BaseExprType to remove ViewFad)
   */
  template <typename A, typename B>
  struct Promote< A, B,
                  typename mpl::enable_if_c< !std::is_convertible<A,B>::value &&
                                             !std::is_convertible<B,A>::value &&
                                             IsExpr<A>::value &&
                                             std::is_convertible< B, typename BaseExprType< typename A::value_type >::type >::value
                                             >::type >
  {
    typedef typename BaseExprType<A>::type type;
  };

  /*!
   * \brief Specialization of Promote when B is an expression and A is
   * convertible to its value-type, e.g., Promote< double, fad-expression >
   * (using BaseExprType to remove ViewFad)
   */
  template <typename A, typename B>
  struct Promote< A, B,
                  typename mpl::enable_if_c< !std::is_convertible<A,B>::value &&
                                             !std::is_convertible<B,A>::value &&
                                             IsExpr<B>::value &&
                                              std::is_convertible< A, typename BaseExprType< typename B::value_type >::type >::value
                                             >::type >
  {
    typedef typename BaseExprType<B>::type type;
  };

  /*!
   * \brief Specialization of Promote when A and B are (different) expressions,
   * with the same value type, e.g, Promote< fad-expr1, fad-expr2 >
   * (using BaseExprType to remove ViewFad)
   */
  template <typename A, typename B>
  struct Promote< A, B,
                  typename mpl::enable_if_c< !std::is_convertible<A,B>::value &&
                                             !std::is_convertible<B,A>::value &&
                                             IsExpr<A>::value &&
                                             IsExpr<B>::value &&
                                             std::is_same< typename BaseExprType< typename A::value_type >::type,
                                                           typename BaseExprType< typename B::value_type >::type >::value
                                             >::type >
  {
    typedef typename BaseExprType<A>::type A_base_expr_type;
    typedef typename BaseExprType<B>::type B_base_expr_type;
    typedef typename Promote< A_base_expr_type, B_base_expr_type >::type type;
  };

  //! Specialization of %Promote to builtin types
#define SACADO_PROMOTE_SPECIALIZATION(type1,type2,type3)   \
  template <> struct Promote< type1, type2, void > {       \
    typedef type3 type;                                    \
  };                                                       \
  template <> struct Promote< type2, type1, void > {       \
    typedef type3 type;                                    \
  };

  SACADO_PROMOTE_SPECIALIZATION(double,float,double)
  SACADO_PROMOTE_SPECIALIZATION(double,long,double)
  SACADO_PROMOTE_SPECIALIZATION(double,int,double)
  SACADO_PROMOTE_SPECIALIZATION(float,long,float)
  SACADO_PROMOTE_SPECIALIZATION(float,int,float)
#ifdef HAVE_SACADO_COMPLEX
  SACADO_PROMOTE_SPECIALIZATION(std::complex<double>,std::complex<float>,std::complex<double>)
  SACADO_PROMOTE_SPECIALIZATION(std::complex<double>,double,std::complex<double>)
  SACADO_PROMOTE_SPECIALIZATION(std::complex<double>,float,std::complex<double>)
  SACADO_PROMOTE_SPECIALIZATION(std::complex<double>,long,std::complex<double>)
  SACADO_PROMOTE_SPECIALIZATION(std::complex<double>,int,std::complex<double>)
  SACADO_PROMOTE_SPECIALIZATION(std::complex<float>,float,std::complex<float>)
  SACADO_PROMOTE_SPECIALIZATION(std::complex<float>,long,std::complex<float>)
  SACADO_PROMOTE_SPECIALIZATION(std::complex<float>,int,std::complex<float>)
#endif // HAVE_SACADO_COMPLEX

#undef SACADO_PROMOTE_SPECIALIZATION

   // Macros for building proper Promote specialization for AD types

#define SACADO_AD_PROMOTE_SPEC(NS, AD) /* */

#define SACADO_AD_PROMOTE_SPEC2(NS, AD) /* */

#define SACADO_FAD_PROMOTE_SPEC(NS, FAD) /* */

#define SACADO_SFAD_PROMOTE_SPEC(NS, FAD) /* */

#define SACADO_EXPR_PROMOTE_SPEC(NS) /* */

#define SACADO_VFAD_PROMOTE_SPEC(NS) /* */

#define SACADO_RAD_PROMOTE_SPEC(NS)                                     \
  namespace NS {                                                        \
    template <typename> class ADvar;                                    \
    template <typename> class ADvari;                                   \
  }                                                                     \
  template <typename T>                                                 \
  struct OverrideDefaultPromote< NS :: ADvari <T>& > {                  \
    static const bool value = true;                                     \
  };                                                                    \
  template <typename T>                                                 \
  struct Promote< NS :: ADvar <T>,                                      \
                  NS :: ADvari <T>& > {                                 \
    typedef NS :: ADvar <T> type;                                       \
  };                                                                    \
  template <typename T>                                                 \
  struct Promote< NS :: ADvari <T>&,                                    \
                  NS :: ADvar <T> > {                                   \
    typedef NS :: ADvar <T> type;                                       \
  };                                                                    \
  template <typename T>                                                 \
  struct Promote< NS :: ADvari <T>&,                                    \
                  typename NS :: ADvari <T>::value_type > {             \
    typedef NS :: ADvar <T> type;                                       \
  };                                                                    \
  template <typename T>                                                 \
  struct Promote< typename NS :: ADvari <T>::value_type,                \
                  NS :: ADvari <T>& > {                                 \
    typedef NS :: ADvar <T> type;                                       \
  };                                                                    \
  template <typename T>                                                 \
  struct Promote< NS :: ADvari <T>&,                                    \
                  typename dummy< typename NS :: ADvari <T>::value_type, \
                                  typename NS :: ADvari <T>::scalar_type \
                                  >::type > {                           \
    typedef NS :: ADvar <T> type;                                       \
  };                                                                    \
  template <typename T>                                                 \
  struct Promote< typename dummy< typename NS :: ADvari <T>::value_type, \
                                  typename NS :: ADvari <T>::scalar_type \
                                  >::type,                              \
                  NS :: ADvari <T>& > {                                 \
    typedef NS :: ADvar <T> type;                                       \
  };

  //
  // We define defaults for all of the traits to make Sacado easier to use.
  // The default choices are based on what appears to be the "safest" choice
  // for any scalar type.  They may not work in all cases, in which case a
  // specialization should be provided.
  //

  //! Base template specification for %ScalarType
  /*!
   * The %ScalarType classes provide a mechanism for computing the
   * base underlying type of nested AD classes
   */
  template <typename T> struct ScalarType {
    typedef T type;
  };

  //! Specialization of %ScalarType for const types
  /*!
   * This should work for most types
   */
  template <typename T> struct ScalarType<const T> {
    typedef const typename ScalarType<T>::type type;
  };

  //! Base template specification for %ValueType
  /*!
   * The %ValueType classes provide a mechanism for computing the
   * the type stored in AD classes
   */
  template <typename T> struct ValueType {
    typedef T type;
  };

  //! Specialization of %ValueType for const types
  /*!
   * This should work for most types
   */
  template <typename T> struct ValueType<const T> {
    typedef const typename ValueType<T>::type type;
  };

  //! Base template specification for %IsADType
  /*!
   * The %IsADType classes provide a mechanism for
   * determining whether a type is an AD type
   */
  template <typename T> struct IsADType {
    static const bool value = false;
  };

  //! Base template specification for %IsScalarType
  /*!
   * The %IsScalarType classes provide a mechanism for
   * determining whether a type is a scalar type (float, double, etc...)
   */
  template <typename T> struct IsScalarType {
    static const bool value = false;
  };

  //! Base template specification for %IsSimdType
  /*!
   * The %IsSimdType classes provide a mechanism for computing the
   * determining whether a type is a Simd type (Doubles, MP::Vector, ...)
   */
  template <typename T> struct IsSimdType {
    static const bool value = false;
  };

  //! Base template specification for %Value
  /*!
   * The %Value functor returns the value of an AD type.
   */
  template <typename T> struct Value {
    SACADO_INLINE_FUNCTION
    static const T& eval(const T& x) { return x; }
  };

  //! Specialization of Value for const types
  template <typename T> struct Value<const T> {
    typedef typename ValueType<T>::type value_type;
    SACADO_INLINE_FUNCTION
    static const value_type& eval(const T& x) {
      return Value<T>::eval(x);
    }
  };

  //! Base template specification for %ScalarValue
  /*!
   * The %ScalarValue functor returns the base scalar value of an AD type,
   * i.e., something that isn't an AD type.
   */
  template <typename T> struct ScalarValue {
    SACADO_INLINE_FUNCTION
    static const T& eval(const T& x) { return x; }
  };

  //! Specialization of ScalarValue for const types
  template <typename T> struct ScalarValue<const T> {
    typedef typename ScalarType<T>::type scalar_type;
    SACADO_INLINE_FUNCTION
    static const scalar_type& eval(const T& x) {
      return ScalarValue<T>::eval(x);
    }
  };

  //! A simple template function for invoking ScalarValue<>
  template <typename T>
  SACADO_INLINE_FUNCTION
  typename ScalarType<T>::type scalarValue(const T& x) {
    return ScalarValue<T>::eval(x);
  }

  //! Base template specification for marking constants
  template <typename T> struct MarkConstant {
    SACADO_INLINE_FUNCTION
    static void eval(T& x) {}
  };

  //! Base template specification for string names of types
  template <typename T> struct StringName {
    static std::string eval() { return ""; }
  };

  //! Base template specification for testing equivalence
  template <typename T> struct IsEqual {
    SACADO_INLINE_FUNCTION
    static bool eval(const T& x, const T& y) { return x == y; }
  };

  //! Base template specification for testing whether type is statically sized
  template <typename T> struct IsStaticallySized {
    static const bool value = false;
  };

  //! Specialization of %IsStaticallySized for const types
  /*!
   * This should work for most types
   */
  template <typename T> struct IsStaticallySized<const T> {
    static const bool value = IsStaticallySized<T>::value;
  };

  //! Base template specification for static size
  template <typename T> struct StaticSize {
    static const unsigned value = 0;
  };

  //! Base template specification for whether a type is a Fad type
  template <typename T> struct IsFad {
    static const bool value = false;
  };

  //! Base template specification for whether a type is a Fad type
  template <typename T> struct IsFad< const T >
  {
    static const bool value = IsFad<T>::value;
  };

  //! Remove const from a type
  template <typename T>
  struct RemoveConst {
    typedef T type;
  };

  //! Remove const from a type
  template <typename T>
  struct RemoveConst< const T > {
    typedef T type;
  };

  //! Specialization of above classes to builtin types
#define SACADO_BUILTIN_SPECIALIZATION(t,NAME)             \
  template <> struct ScalarType< t > {                    \
    typedef t type;                                       \
  };                                                      \
  template <> struct ValueType< t > {                     \
    typedef t type;                                       \
  };                                                      \
  template <> struct IsADType< t > {                      \
    static const bool value = false;                      \
  };                                                      \
  template <> struct IsScalarType< t > {                  \
    static const bool value = true;                       \
  };                                                      \
  template <> struct Value< t > {                         \
    SACADO_INLINE_FUNCTION                                \
    static const t& eval(const t& x) { return x; }        \
  };                                                      \
  template <> struct ScalarValue< t > {                   \
    SACADO_INLINE_FUNCTION                                \
    static const t& eval(const t& x) { return x; }        \
  };                                                      \
  template <> struct StringName< t > {                    \
    static std::string eval() { return NAME; }            \
  };                                                      \
  template <> struct IsEqual< t > {                       \
    SACADO_INLINE_FUNCTION                                \
    static bool eval(const t& x, const t& y) {            \
      return x == y; }                                    \
  };                                                      \
  template <> struct IsStaticallySized< t > {             \
    static const bool value = true;                       \
  };

#define SACADO_BUILTIN_SPECIALIZATION_COMPLEX(t,NAME)     \
  template <> struct ScalarType< t > {                    \
    typedef t type;                                       \
  };                                                      \
  template <> struct ValueType< t > {                     \
    typedef t type;                                       \
  };                                                      \
  template <> struct IsADType< t > {                      \
    static const bool value = false;                      \
  };                                                      \
  template <> struct IsScalarType< t > {                  \
    static const bool value = true;                       \
  };                                                      \
  template <> struct Value< t > {                         \
    static const t& eval(const t& x) { return x; }        \
  };                                                      \
  template <> struct ScalarValue< t > {                   \
    static const t& eval(const t& x) { return x; }        \
  };                                                      \
  template <> struct StringName< t > {                    \
    static std::string eval() { return NAME; }            \
  };                                                      \
  template <> struct IsEqual< t > {                       \
    static bool eval(const t& x, const t& y) {            \
      return x == y; }                                    \
  };                                                      \
  template <> struct IsStaticallySized< t > {             \
    static const bool value = true;                       \
  };

  SACADO_BUILTIN_SPECIALIZATION(char,"char")
  SACADO_BUILTIN_SPECIALIZATION(float,"float")
  SACADO_BUILTIN_SPECIALIZATION(double,"double")
  SACADO_BUILTIN_SPECIALIZATION(int,"int")
  SACADO_BUILTIN_SPECIALIZATION(unsigned int,"unsigned int")
  SACADO_BUILTIN_SPECIALIZATION(long,"long")
  SACADO_BUILTIN_SPECIALIZATION(unsigned long,"unsigned long")
  SACADO_BUILTIN_SPECIALIZATION(bool,"bool")
#ifdef HAVE_SACADO_COMPLEX
  SACADO_BUILTIN_SPECIALIZATION_COMPLEX(std::complex<double>,"std::complex<double>")
  SACADO_BUILTIN_SPECIALIZATION_COMPLEX(std::complex<float>,"std::complex<float>")
#endif

#undef SACADO_BUILTIN_SPECIALIZATION
#undef SACADO_BUILTIN_SPECIALIZATION_COMPLEX

template< typename T , T v , bool NonZero = ( v != T(0) ) >
struct integral_nonzero
{
  // Declaration of 'static const' causes an unresolved linker symbol in debug
  // static const T value = v ;
  enum { value = T(v) };
  typedef T value_type ;
  typedef integral_nonzero<T,v> type ;
  SACADO_INLINE_FUNCTION integral_nonzero() {}
  SACADO_INLINE_FUNCTION integral_nonzero( const T & ) {}
  SACADO_INLINE_FUNCTION integral_nonzero( const integral_nonzero & ) {}
  SACADO_INLINE_FUNCTION integral_nonzero& operator=(const integral_nonzero &) {return *this;}
  SACADO_INLINE_FUNCTION integral_nonzero& operator=(const T &) {return *this;}
};

template< typename T , T zero >
struct integral_nonzero<T,zero,false>
{
  T value ;
  typedef T value_type ;
  typedef integral_nonzero<T,0> type ;
  SACADO_INLINE_FUNCTION integral_nonzero() : value() {}
  SACADO_INLINE_FUNCTION integral_nonzero( const T & v ) : value(v) {}
  SACADO_INLINE_FUNCTION integral_nonzero( const integral_nonzero & v) : value(v.value) {}
  SACADO_INLINE_FUNCTION integral_nonzero& operator=(const integral_nonzero & v) { value = v.value; return *this; }
  SACADO_INLINE_FUNCTION integral_nonzero& operator=(const T & v) { value = v; return *this; }
};

} // namespace Sacado

#endif // SACADO_TRAITS_HPP
