// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_EXP_MP_VECTOR_HPP
#define SACADO_FAD_EXP_MP_VECTOR_HPP

#include "Sacado_MP_Vector.hpp"

namespace Sacado {

  namespace Fad {
  namespace Exp {

    //! Expression template specialization tag for Fad< MP::Vector >
    class ExprSpecMPVector {};

    //! Specialization of extender for MP::Vector scalar types
    /*!
     * Extends interface to add val(), dx(), fastAccessDx() functions that
     * take an ensemble component argument, for flattening nested scalar type
     * expression templates.
     */
    template <typename T>
    class Extender<
      T,
      typename std::enable_if<
        Sacado::is_mp_vector<typename T::value_type>::value >::type
      > : public T
    {
    public:

      typedef typename T::value_type value_type;
      typedef typename value_type::value_type val_type;

      // Define expression template specialization
      typedef ExprSpecMPVector expr_spec_type;

      // Bring in constructors
      using T::T;

      // Bring in methods we are overloading
      using T::val;
      using T::dx;
      using T::fastAccessDx;

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      const val_type& val(int j) const { return T::val().fastAccessCoeff(j); }

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      val_type& val(int j) { return T::val().fastAccessCoeff(j); }

      //! Returns derivative component \c i with bounds checking
      KOKKOS_INLINE_FUNCTION
      val_type dx(int i, int j) const {
        return this->size() ? this->dx_[i].fastAccessCoeff(j) : val_type(0.0);
      }

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      val_type& fastAccessDx(int i, int j) {
        return this->dx_[i].fastAccessCoeff(j);
      }

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      const val_type& fastAccessDx(int i, int j) const {
        return this->dx_[i].fastAccessCoeff(j);
      }

    };

    //! Specialization of ExprAssign for MP::Vector scalar types
    template <typename DstType>
    class ExprAssign<
      DstType,
      typename std::enable_if<
        std::is_same< typename DstType::expr_spec_type, ExprSpecMPVector >::value
        >::type
      > {
    public:

      //! Typename of values
      typedef typename DstType::value_type value_type;

      // MP::Vector size (assuming static, because that's all we care about)
      static const int VecNum = Sacado::StaticSize<value_type>::value;

      //! Implementation of dst = x
      template <typename SrcType>
      KOKKOS_INLINE_FUNCTION
      static void assign_equal(DstType& dst, const SrcType& x)
      {
        const int xsz = x.size();

        if (xsz != dst.size())
          dst.resizeAndZero(xsz);

        const int sz = dst.size();

        // For ViewStorage, the resize above may not in fact resize the
        // derivative array, so it is possible that sz != xsz at this point.
        // The only valid use case here is sz > xsz == 0, so we use sz in the
        // assignment below

        if (sz) {
          if (x.hasFastAccess()) {
            SACADO_FAD_DERIV_LOOP(i,sz)
              for (int j=0; j<VecNum; ++j)
                dst.fastAccessDx(i,j) = x.fastAccessDx(i,j);
          }
          else
            SACADO_FAD_DERIV_LOOP(i,sz)
              for (int j=0; j<VecNum; ++j)
                dst.fastAccessDx(i,j) = x.dx(i,j);
        }

        for (int j=0; j<VecNum; ++j)
          dst.val(j) = x.val(j);
      }

      //! Implementation of dst += x
      template <typename SrcType>
      KOKKOS_INLINE_FUNCTION
      static void assign_plus_equal(DstType& dst, const SrcType& x)
      {
        const int xsz = x.size(), sz = dst.size();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            if (x.hasFastAccess())
              SACADO_FAD_DERIV_LOOP(i,sz)
                for (int j=0; j<VecNum; ++j)
                  dst.fastAccessDx(i,j) += x.fastAccessDx(i,j);
            else
              for (int i=0; i<sz; ++i)
                for (int j=0; j<VecNum; ++j)
                  dst.fastAccessDx(i,j) += x.dx(i,j);
          }
          else {
            dst.resizeAndZero(xsz);
            if (x.hasFastAccess())
              SACADO_FAD_DERIV_LOOP(i,xsz)
                for (int j=0; j<VecNum; ++j)
                  dst.fastAccessDx(i,j) = x.fastAccessDx(i,j);
            else
              SACADO_FAD_DERIV_LOOP(i,xsz)
                for (int j=0; j<VecNum; ++j)
                  dst.fastAccessDx(i,j) = x.dx(i,j);
          }
        }

        for (int j=0; j<VecNum; ++j)
          dst.val(j) += x.val(j);
      }

      //! Implementation of dst -= x
      template <typename SrcType>
      KOKKOS_INLINE_FUNCTION
      static void assign_minus_equal(DstType& dst, const SrcType& x)
      {
        const int xsz = x.size(), sz = dst.size();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            if (x.hasFastAccess())
              SACADO_FAD_DERIV_LOOP(i,sz)
                for (int j=0; j<VecNum; ++j)
                  dst.fastAccessDx(i,j) -= x.fastAccessDx(i,j);
            else
              SACADO_FAD_DERIV_LOOP(i,sz)
                for (int j=0; j<VecNum; ++j)
                  dst.fastAccessDx(i,j) -= x.dx(i,j);
          }
          else {
            dst.resizeAndZero(xsz);
            if (x.hasFastAccess())
              SACADO_FAD_DERIV_LOOP(i,xsz)
                for (int j=0; j<VecNum; ++j)
                  dst.fastAccessDx(i,j) = -x.fastAccessDx(i,j);
            else
              SACADO_FAD_DERIV_LOOP(i,xsz)
                for (int j=0; j<VecNum; ++j)
                  dst.fastAccessDx(i,j) = -x.dx(i,j);
          }
        }

        for (int j=0; j<VecNum; ++j)
          dst.val(j) -= x.val(j);
      }

      //! Implementation of dst *= x
      template <typename SrcType>
      KOKKOS_INLINE_FUNCTION
      static void assign_times_equal(DstType& dst, const SrcType& x)
      {
        const int xsz = x.size(), sz = dst.size();
        const value_type xval = x.val();
        const value_type v = dst.val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            if (x.hasFastAccess())
              SACADO_FAD_DERIV_LOOP(i,sz)
                for (int j=0; j<VecNum; ++j)
                  dst.fastAccessDx(i) = v.fastAccessCoeff(j)*x.fastAccessDx(i,j) + dst.fastAccessDx(i,j)*xval.fastAccessCoeff(j);
            else
              SACADO_FAD_DERIV_LOOP(i,sz)
                for (int j=0; j<VecNum; ++j)
                  dst.fastAccessDx(i) = v.fastAccessCoeff(j)*x.dx(i,j) + dst.fastAccessDx(i,j)*xval.fastAccessCoeff(j);
          }
          else {
            dst.resizeAndZero(xsz);
            if (x.hasFastAccess())
              SACADO_FAD_DERIV_LOOP(i,xsz)
                for (int j=0; j<VecNum; ++j)
                  dst.fastAccessDx(i,j) = v.fastAccessCoeff(j)*x.fastAccessDx(i,j);
            else
              SACADO_FAD_DERIV_LOOP(i,xsz)
                for (int j=0; j<VecNum; ++j)
                  dst.fastAccessDx(i,j) = v.fastAccessCoeff(j)*x.dx(i,j);
          }
        }
        else {
          if (sz) {
            SACADO_FAD_DERIV_LOOP(i,sz)
              for (int j=0; j<VecNum; ++j)
                dst.fastAccessDx(i,j) *= xval.fastAccessCoeff(j);
          }
        }

        for (int j=0; j<VecNum; ++j)
          dst.val(j) *= xval.fastAccessCoeff(j);
      }

      //! Implementation of dst /= x
      template <typename SrcType>
      KOKKOS_INLINE_FUNCTION
      static void assign_divide_equal(DstType& dst, const SrcType& x)
      {
        const int xsz = x.size(), sz = dst.size();
        const value_type xval = x.val();
        const value_type v = dst.val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          const value_type xval2 = xval*xval;
          if (sz) {
            if (x.hasFastAccess())
              SACADO_FAD_DERIV_LOOP(i,sz)
                for (int j=0; j<VecNum; ++j)
                  dst.fastAccessDx(i,j) =
                    ( dst.fastAccessDx(i,j)*xval.fastAccessCoeff(j) - v.fastAccessCoeff(j)*x.fastAccessDx(i,j) ) / xval2.fastAccessCoeff(j);
            else
              SACADO_FAD_DERIV_LOOP(i,sz)
                for (int j=0; j<VecNum; ++j)
                  dst.fastAccessDx(i,j) =
                    ( dst.fastAccessDx(i,j)*xval.fastAccessCoeff(j) - v.fastAccessCoeff(j)*x.dx(i,j) ) / xval2.fastAccessCoeff(j);
          }
          else {
            dst.resizeAndZero(xsz);
            if (x.hasFastAccess())
              SACADO_FAD_DERIV_LOOP(i,xsz)
                for (int j=0; j<VecNum; ++j)
                  dst.fastAccessDx(i,j) = - v.fastAccessCoeff(j)*x.fastAccessDx(i,j) / xval2.fastAccessCoeff(j);
            else
              SACADO_FAD_DERIV_LOOP(i,xsz)
                for (int j=0; j<VecNum; ++j)
                  dst.fastAccessDx(i,j) = -v.fastAccessCoeff(j)*x.dx(i,j) / xval2.fastAccessCoeff(j);
          }
        }
        else {
          if (sz) {
            SACADO_FAD_DERIV_LOOP(i,sz)
              for (int j=0; j<VecNum; ++j)
                dst.fastAccessDx(i,j) /= xval.fastAccessCoeff(j);
          }
        }

        for (int j=0; j<VecNum; ++j)
          dst.val(j) /= xval.fastAccessCoeff(j);
      }

    };

    /*!
     * \brief Specialization of ExprAssign for statically sized Fad types
     * and MP::Vector types
     */
    template <typename DstType>
    class ExprAssign<
      DstType,
      typename std::enable_if<
        Sacado::IsStaticallySized<DstType>::value &&
        std::is_same< typename DstType::expr_spec_type, ExprSpecMPVector >::value
        >::type
      > {
    public:

      //! Typename of values
      typedef typename DstType::value_type value_type;

      // MP::Vector size (assuming static, because that's all we care about)
      static const int VecNum = Sacado::StaticSize<value_type>::value;

      //! Implementation of dst = x
      template <typename SrcType>
      KOKKOS_INLINE_FUNCTION
      static void assign_equal(DstType& dst, const SrcType& x)
      {
        const int sz = dst.size();
        SACADO_FAD_DERIV_LOOP(i,sz)
          for (int j=0; j<VecNum; ++j)
            dst.fastAccessDx(i,j) = x.fastAccessDx(i,j);
        for (int j=0; j<VecNum; ++j)
          dst.val(j) = x.val(j);
      }

      //! Implementation of dst += x
      template <typename SrcType>
      KOKKOS_INLINE_FUNCTION
      static void assign_plus_equal(DstType& dst, const SrcType& x)
      {
        const int sz = dst.size();
        SACADO_FAD_DERIV_LOOP(i,sz)
          for (int j=0; j<VecNum; ++j)
            dst.fastAccessDx(i,j) += x.fastAccessDx(i,j);
        for (int j=0; j<VecNum; ++j)
          dst.val(j) += x.val(j);
      }

      //! Implementation of dst -= x
      template <typename SrcType>
      KOKKOS_INLINE_FUNCTION
      static void assign_minus_equal(DstType& dst, const SrcType& x)
      {
        const int sz = dst.size();
        SACADO_FAD_DERIV_LOOP(i,sz)
          for (int j=0; j<VecNum; ++j)
            dst.fastAccessDx(i,j) -= x.fastAccessDx(i,j);
        for (int j=0; j<VecNum; ++j)
          dst.val(j) -= x.val(j);
      }

      //! Implementation of dst *= x
      template <typename SrcType>
      KOKKOS_INLINE_FUNCTION
      static void assign_times_equal(DstType& dst, const SrcType& x)
      {
        const int sz = dst.size();
        const value_type xval = x.val();
        const value_type v = dst.val();
        SACADO_FAD_DERIV_LOOP(i,sz)
          for (int j=0; j<VecNum; ++j)
            dst.fastAccessDx(i,j) = v.fastAccessCoeff(j)*x.fastAccessDx(i,j) + dst.fastAccessDx(i,j)*xval.fastAccessCoeff(j);
        for (int j=0; j<VecNum; ++j)
          dst.val(j) *= xval.fastAccessCoeff(j);
      }

      //! Implementation of dst /= x
      template <typename SrcType>
      KOKKOS_INLINE_FUNCTION
      static void assign_divide_equal(DstType& dst, const SrcType& x)
      {
        const int sz = dst.size();
        const value_type xval = x.val();
        const value_type xval2 = xval*xval;
        const value_type v = dst.val();
        SACADO_FAD_DERIV_LOOP(i,sz)
          for (int j=0; j<VecNum; ++j)
            dst.fastAccessDx(i,j) =
              ( dst.fastAccessDx(i,j)*xval.fastAccessCoeff(j) - v.fastAccessCoeff(j)*x.fastAccessDx(i,j) )/ xval2.fastAccessCoeff(j);
        for (int j=0; j<VecNum; ++j)
          dst.val(j) /= xval.fastAccessCoeff(j);
      }

    };

  } // namespace Exp
  } // namespace Fad

} // namespace Sacado

// Specialize expression template operators to add similar extensions
#include "Sacado_Fad_Exp_Ops.hpp"

#define FAD_UNARYOP_MACRO(OPNAME,OP,USING,MPVALUE,VALUE,DX,FASTACCESSDX) \
namespace Sacado {                                                      \
  namespace Fad {                                                       \
  namespace Exp {                                                       \
                                                                        \
    template <typename T>                                               \
    class OP< T,ExprSpecMPVector > :                                    \
      public Expr< OP< T,ExprSpecMPVector > > {                         \
    public:                                                             \
                                                                        \
      typedef typename std::remove_cv<T>::type ExprT;                   \
      typedef typename ExprT::value_type value_type;                    \
      typedef typename ExprT::scalar_type scalar_type;                  \
                                                                        \
      typedef typename value_type::value_type val_type;                 \
                                                                        \
      typedef ExprSpecMPVector expr_spec_type;                          \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      OP(const T& expr_) : expr(expr_)  {}                              \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      int size() const { return expr.size(); }                          \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool hasFastAccess() const {                                      \
        return expr.hasFastAccess();                                    \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        USING                                                           \
        return MPVALUE;                                                 \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      val_type val(int j) const {                                       \
        USING                                                           \
        return VALUE;                                                   \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      val_type dx(int i, int j) const {                                 \
        USING                                                           \
        return DX;                                                      \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      val_type fastAccessDx(int i, int j) const {                       \
        USING                                                           \
        return FASTACCESSDX;                                            \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      const T& expr;                                                    \
    };                                                                  \
                                                                        \
  }                                                                     \
  }                                                                     \
                                                                        \
}

FAD_UNARYOP_MACRO(operator+,
                  UnaryPlusOp,
                  ;,
                  expr.val(),
                  expr.val(j),
                  expr.dx(i,j),
                  expr.fastAccessDx(i,j))
FAD_UNARYOP_MACRO(operator-,
                  UnaryMinusOp,
                  ;,
                  -expr.val(),
                  -expr.val(j),
                  -expr.dx(i,j),
                  -expr.fastAccessDx(i,j))
FAD_UNARYOP_MACRO(exp,
                  ExpOp,
                  using std::exp;,
                  exp(expr.val()),
                  exp(expr.val(j)),
                  exp(expr.val(j))*expr.dx(i,j),
                  exp(expr.val(j))*expr.fastAccessDx(i,j))
FAD_UNARYOP_MACRO(log,
                  LogOp,
                  using std::log;,
                  log(expr.val()),
                  log(expr.val(j)),
                  expr.dx(i,j)/expr.val(j),
                  expr.fastAccessDx(i,j)/expr.val(j))
FAD_UNARYOP_MACRO(log10,
                  Log10Op,
                  using std::log10; using std::log;,
                  log10(expr.val()),
                  log10(expr.val(j)),
                  expr.dx(i,j)/( log(value_type(10))*expr.val()),
                  expr.fastAccessDx(i,j) / ( log(value_type(10))*expr.val()))
FAD_UNARYOP_MACRO(sqrt,
                  SqrtOp,
                  using std::sqrt;,
                  sqrt(expr.val()),
                  sqrt(expr.val(j)),
                  expr.dx(i,j)/(value_type(2)* sqrt(expr.val())),
                  expr.fastAccessDx(i,j)/(value_type(2)* sqrt(expr.val())))
FAD_UNARYOP_MACRO(cos,
                  CosOp,
                  using std::cos; using std::sin;,
                  cos(expr.val()),
                  cos(expr.val(j)),
                  -expr.dx(i,j)* sin(expr.val()),
                  -expr.fastAccessDx(i,j)* sin(expr.val()))
FAD_UNARYOP_MACRO(sin,
                  SinOp,
                  using std::cos; using std::sin;,
                  sin(expr.val()),
                  sin(expr.val(j)),
                  expr.dx(i,j)* cos(expr.val()),
                  expr.fastAccessDx(i,j)* cos(expr.val()))
FAD_UNARYOP_MACRO(tan,
                  TanOp,
                  using std::tan;,
                  tan(expr.val()),
                  tan(expr.val(j)),
                  expr.dx(i,j)*
                    (value_type(1)+ tan(expr.val())* tan(expr.val())),
                  expr.fastAccessDx(i,j)*
                    (value_type(1)+ tan(expr.val())* tan(expr.val())))
FAD_UNARYOP_MACRO(acos,
                  ACosOp,
                  using std::acos; using std::sqrt;,
                  acos(expr.val()),
                  acos(expr.val(j)),
                  -expr.dx(i,j)/ sqrt(value_type(1)-expr.val()*expr.val()),
                  -expr.fastAccessDx(i,j) /
                    sqrt(value_type(1)-expr.val()*expr.val()))
FAD_UNARYOP_MACRO(asin,
                  ASinOp,
                  using std::asin; using std::sqrt;,
                  asin(expr.val()),
                  asin(expr.val(j)),
                  expr.dx(i,j)/ sqrt(value_type(1)-expr.val()*expr.val()),
                  expr.fastAccessDx(i,j) /
                    sqrt(value_type(1)-expr.val()*expr.val()))
FAD_UNARYOP_MACRO(atan,
                  ATanOp,
                  using std::atan;,
                  atan(expr.val()),
                  atan(expr.val(j)),
                  expr.dx(i,j)/(value_type(1)+expr.val()*expr.val()),
                  expr.fastAccessDx(i,j)/(value_type(1)+expr.val()*expr.val()))
FAD_UNARYOP_MACRO(cosh,
                  CoshOp,
                  using std::cosh; using std::sinh;,
                  cosh(expr.val()),
                  cosh(expr.val(j)),
                  expr.dx(i,j)* sinh(expr.val()),
                  expr.fastAccessDx(i,j)* sinh(expr.val()))
FAD_UNARYOP_MACRO(sinh,
                  SinhOp,
                  using std::cosh; using std::sinh;,
                  sinh(expr.val()),
                  sinh(expr.val(j)),
                  expr.dx(i,j)* cosh(expr.val()),
                  expr.fastAccessDx(i,j)* cosh(expr.val()))
FAD_UNARYOP_MACRO(tanh,
                  TanhOp,
                  using std::tanh; using std::cosh;,
                  tanh(expr.val()),
                  tanh(expr.val(j)),
                  expr.dx(i,j)/( cosh(expr.val())* cosh(expr.val())),
                  expr.fastAccessDx(i,j) /
                    ( cosh(expr.val())* cosh(expr.val())))
FAD_UNARYOP_MACRO(acosh,
                  ACoshOp,
                  using std::acosh; using std::sqrt;,
                  acosh(expr.val()),
                  acosh(expr.val(j)),
                  expr.dx(i,j)/ sqrt((expr.val()-value_type(1)) *
                                       (expr.val()+value_type(1))),
                  expr.fastAccessDx(i,j)/ sqrt((expr.val()-value_type(1)) *
                                                 (expr.val()+value_type(1))))
FAD_UNARYOP_MACRO(asinh,
                  ASinhOp,
                  using std::asinh; using std::sqrt;,
                  asinh(expr.val()),
                  asinh(expr.val(j)),
                  expr.dx(i,j)/ sqrt(value_type(1)+expr.val()*expr.val()),
                  expr.fastAccessDx(i,j)/ sqrt(value_type(1)+
                                                 expr.val()*expr.val()))
FAD_UNARYOP_MACRO(atanh,
                  ATanhOp,
                  using std::atanh;,
                  atanh(expr.val()),
                  atanh(expr.val(j)),
                  expr.dx(i,j)/(value_type(1)-expr.val()*expr.val()),
                  expr.fastAccessDx(i,j)/(value_type(1)-
                                                 expr.val()*expr.val()))
FAD_UNARYOP_MACRO(abs,
                  AbsOp,
                  using std::abs;,
                  abs(expr.val()),
                  abs(expr.val(j)),
                  if_then_else( expr.val() >= 0, expr.dx(i,j), value_type(-expr.dx(i,j)) ),
                  if_then_else( expr.val() >= 0, expr.fastAccessDx(i,j), value_type(-expr.fastAccessDx(i,j)) ) )
FAD_UNARYOP_MACRO(fabs,
                  FAbsOp,
                  using std::fabs;,
                  fabs(expr.val()),
                  fabs(expr.val(j)),
                  if_then_else( expr.val() >= 0, expr.dx(i,j), value_type(-expr.dx(i,j)) ),
                  if_then_else( expr.val() >= 0, expr.fastAccessDx(i,j), value_type(-expr.fastAccessDx(i,j)) ) )
FAD_UNARYOP_MACRO(cbrt,
                  CbrtOp,
                  using std::cbrt;,
                  cbrt(expr.val()),
                  cbrt(expr.val(j)),
                  expr.dx(i,j)/(value_type(3)*cbrt(expr.val()*expr.val())),
                  expr.fastAccessDx(i,j)/(value_type(3)*cbrt(expr.val()*expr.val())))

#undef FAD_UNARYOP_MACRO

namespace Sacado {
  namespace Fad {
  namespace Exp {

    // For MP::Vector scalar type, promote constants up to expression's value
    // type.  If the constant type is the same as the value type, we can store
    // the constant as a reference.  If it isn't, we must copy it into a new
    // value type object.  We do this so that we can always access the constant
    // as a value type.
    template <typename ConstType, typename ValueType>
    struct ConstTypeRef {
      typedef ValueType type;
    };

    template <typename ValueType>
    struct ConstTypeRef<ValueType, ValueType> {
      typedef ValueType& type;
    };
  }
  }
}

#define FAD_BINARYOP_MACRO(OPNAME,OP,USING,MPVALUE,VALUE,DX,CDX1,CDX2,FASTACCESSDX,MPVAL_CONST_DX_1,MPVAL_CONST_DX_2,VAL_CONST_DX_1,VAL_CONST_DX_2,CONST_DX_1,CONST_DX_2,CONST_FASTACCESSDX_1,CONST_FASTACCESSDX_2) \
namespace Sacado {                                                      \
  namespace Fad {                                                       \
  namespace Exp {                                                       \
                                                                        \
    template <typename T1, typename T2 >                                \
    class OP< T1, T2, false, false, ExprSpecMPVector > :                \
      public Expr< OP< T1, T2, false, false, ExprSpecMPVector > > {     \
    public:                                                             \
                                                                        \
      typedef typename std::remove_cv<T1>::type ExprT1;                 \
      typedef typename std::remove_cv<T2>::type ExprT2;                 \
      typedef typename ExprT1::value_type value_type_1;                 \
      typedef typename ExprT2::value_type value_type_2;                 \
      typedef typename Sacado::Promote<value_type_1,                    \
                                       value_type_2>::type value_type;  \
                                                                        \
      typedef typename ExprT1::scalar_type scalar_type_1;               \
      typedef typename ExprT2::scalar_type scalar_type_2;               \
      typedef typename Sacado::Promote<scalar_type_1,                   \
                                       scalar_type_2>::type scalar_type; \
                                                                        \
      typedef typename value_type::value_type val_type;                 \
                                                                        \
      typedef ExprSpecMPVector expr_spec_type;                          \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      OP(const T1& expr1_, const T2& expr2_) :                          \
        expr1(expr1_), expr2(expr2_) {}                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      int size() const {                                                \
        const int sz1 = expr1.size(), sz2 = expr2.size();               \
        return sz1 > sz2 ? sz1 : sz2;                                   \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool hasFastAccess() const {                                      \
        return expr1.hasFastAccess() && expr2.hasFastAccess();          \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        USING                                                           \
        return MPVALUE;                                                 \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      val_type val(int j) const {                                       \
        USING                                                           \
        return VALUE;                                                   \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      val_type dx(int i, int j) const {                                 \
        USING                                                           \
        const int sz1 = expr1.size(), sz2 = expr2.size();               \
        if (sz1 > 0 && sz2 > 0)                                         \
          return DX;                                                    \
        else if (sz1 > 0)                                               \
          return CDX2;                                                  \
        else                                                            \
          return CDX1;                                                  \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      val_type fastAccessDx(int i, int j) const {                       \
        USING                                                           \
        return FASTACCESSDX;                                            \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      const T1& expr1;                                                  \
      const T2& expr2;                                                  \
                                                                        \
    };                                                                  \
                                                                        \
    template <typename T1, typename T2>                                 \
    class OP< T1, T2, false, true, ExprSpecMPVector > :                 \
      public Expr< OP< T1, T2, false, true, ExprSpecMPVector > > {      \
    public:                                                             \
                                                                        \
      typedef typename std::remove_cv<T1>::type ExprT1;                 \
      typedef T2 ConstT;                                                \
      typedef typename ExprT1::value_type value_type;                   \
      typedef typename ExprT1::scalar_type scalar_type;                 \
                                                                        \
      typedef typename value_type::value_type val_type;                 \
                                                                        \
      typedef ExprSpecMPVector expr_spec_type;                          \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      OP(const T1& expr1_, const ConstT& c_) :                          \
        expr1(expr1_), c(c_) {}                                         \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      int size() const {                                                \
        return expr1.size();                                            \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool hasFastAccess() const {                                      \
        return expr1.hasFastAccess();                                   \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        USING                                                           \
        return MPVAL_CONST_DX_2;                                        \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      val_type val(int j) const {                                       \
        USING                                                           \
        return VAL_CONST_DX_2;                                          \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      val_type dx(int i, int j) const {                                 \
        USING                                                           \
        return CONST_DX_2;                                              \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      val_type fastAccessDx(int i, int j) const {                       \
        USING                                                           \
        return CONST_FASTACCESSDX_2;                                    \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      const T1& expr1;                                                  \
      const typename ConstTypeRef<ConstT,value_type>::type c;           \
    };                                                                  \
                                                                        \
    template <typename T1, typename T2>                                 \
    class OP< T1, T2, true, false,ExprSpecMPVector > :                  \
      public Expr< OP< T1, T2, true, false, ExprSpecMPVector > > {      \
    public:                                                             \
                                                                        \
      typedef typename std::remove_cv<T2>::type ExprT2;                 \
      typedef T1 ConstT;                                                \
      typedef typename ExprT2::value_type value_type;                   \
      typedef typename ExprT2::scalar_type scalar_type;                 \
                                                                        \
      typedef typename value_type::value_type val_type;                 \
                                                                        \
      typedef ExprSpecMPVector expr_spec_type;                          \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      OP(const ConstT& c_, const T2& expr2_) :                          \
        c(c_), expr2(expr2_) {}                                         \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      int size() const {                                                \
        return expr2.size();                                            \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      bool hasFastAccess() const {                                      \
        return expr2.hasFastAccess();                                   \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      value_type val() const {                                          \
        USING                                                           \
        return MPVAL_CONST_DX_1;                                        \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      val_type val(int j) const {                                       \
        USING                                                           \
        return VAL_CONST_DX_1;                                          \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      val_type dx(int i, int j) const {                                 \
        USING                                                           \
        return CONST_DX_1;                                              \
      }                                                                 \
                                                                        \
      KOKKOS_INLINE_FUNCTION                                            \
      val_type fastAccessDx(int i, int j) const {                       \
        USING                                                           \
        return CONST_FASTACCESSDX_1;                                    \
      }                                                                 \
                                                                        \
    protected:                                                          \
                                                                        \
      const typename ConstTypeRef<ConstT,value_type>::type c;           \
      const T2& expr2;                                                  \
    };                                                                  \
                                                                        \
  }                                                                     \
  }                                                                     \
                                                                        \
}


FAD_BINARYOP_MACRO(operator+,
                   AdditionOp,
                   ;,
                   expr1.val() + expr2.val(),
                   expr1.val(j) + expr2.val(j),
                   expr1.dx(i,j) + expr2.dx(i,j),
                   expr2.dx(i,j),
                   expr1.dx(i,j),
                   expr1.fastAccessDx(i,j) + expr2.fastAccessDx(i,j),
                   c + expr2.val(),
                   expr1.val() + c,
                   c.fastAccessCoeff(j) + expr2.val(j),
                   expr1.val(j) + c.fastAccessCoeff(j),
                   expr2.dx(i,j),
                   expr1.dx(i,j),
                   expr2.fastAccessDx(i,j),
                   expr1.fastAccessDx(i,j))
FAD_BINARYOP_MACRO(operator-,
                   SubtractionOp,
                   ;,
                   expr1.val() - expr2.val(),
                   expr1.val(j) - expr2.val(j),
                   expr1.dx(i,j) - expr2.dx(i,j),
                   -expr2.dx(i,j),
                   expr1.dx(i,j),
                   expr1.fastAccessDx(i,j) - expr2.fastAccessDx(i,j),
                   c - expr2.val(),
                   expr1.val() - c,
                   c.fastAccessCoeff(j) - expr2.val(j),
                   expr1.val(j) - c.fastAccessCoeff(j),
                   -expr2.dx(i,j),
                   expr1.dx(i,j),
                   -expr2.fastAccessDx(i,j),
                   expr1.fastAccessDx(i,j))
FAD_BINARYOP_MACRO(operator*,
                   MultiplicationOp,
                   ;,
                   expr1.val() * expr2.val(),
                   expr1.val(j) * expr2.val(j),
                   expr1.val(j)*expr2.dx(i,j) + expr1.dx(i,j)*expr2.val(j),
                   expr1.val(j)*expr2.dx(i,j),
                   expr1.dx(i,j)*expr2.val(j),
                   expr1.val(j)*expr2.fastAccessDx(i,j) +
                     expr1.fastAccessDx(i,j)*expr2.val(j),
                   c * expr2.val(),
                   expr1.val() * c,
                   c.fastAccessCoeff(j) * expr2.val(j),
                   expr1.val(j) * c.fastAccessCoeff(j),
                   c.fastAccessCoeff(j)*expr2.dx(i,j),
                   expr1.dx(i,j)*c.fastAccessCoeff(j),
                   c.fastAccessCoeff(j)*expr2.fastAccessDx(i,j),
                   expr1.fastAccessDx(i,j)*c.fastAccessCoeff(j))
FAD_BINARYOP_MACRO(operator/,
                   DivisionOp,
                   ;,
                   expr1.val() / expr2.val(),
                   expr1.val(j) / expr2.val(j),
                   (expr1.dx(i,j)*expr2.val(j) - expr2.dx(i,j)*expr1.val(j)) /
                     (expr2.val(j)*expr2.val(j)),
                   -expr2.dx(i,j)*expr1.val(j) / (expr2.val(j)*expr2.val(j)),
                   expr1.dx(i,j)/expr2.val(j),
                   (expr1.fastAccessDx(i,j)*expr2.val(j) -
                      expr2.fastAccessDx(i,j)*expr1.val(j)) /
                      (expr2.val(j)*expr2.val(j)),
                   c / expr2.val(),
                   expr1.val() / c,
                   c.fastAccessCoeff(j) / expr2.val(j),
                   expr1.val(j) / c.fastAccessCoeff(j),
                   -expr2.dx(i,j)*c.fastAccessCoeff(j) / (expr2.val(j)*expr2.val(j)),
                   expr1.dx(i,j)/c.fastAccessCoeff(j),
                   -expr2.fastAccessDx(i,j)*c.fastAccessCoeff(j) / (expr2.val(j)*expr2.val(j)),
                   expr1.fastAccessDx(i,j)/c.fastAccessCoeff(j))
FAD_BINARYOP_MACRO(atan2,
                   Atan2Op,
                   using std::atan2;,
                   atan2(expr1.val(), expr2.val()),
                   atan2(expr1.val(j), expr2.val(j)),
                   (expr2.val(j)*expr1.dx(i,j) - expr1.val(j)*expr2.dx(i,j))/
                        (expr1.val(j)*expr1.val(j) + expr2.val(j)*expr2.val(j)),
                   -expr1.val(j)*expr2.dx(i,j)/
                        (expr1.val(j)*expr1.val(j) + expr2.val(j)*expr2.val(j)),
                   expr2.val(j)*expr1.dx(i,j)/
                        (expr1.val(j)*expr1.val(j) + expr2.val(j)*expr2.val(j)),
                   (expr2.val(j)*expr1.fastAccessDx(i,j) - expr1.val(j)*expr2.fastAccessDx(i,j))/
                        (expr1.val(j)*expr1.val(j) + expr2.val(j)*expr2.val(j)),
                   atan2(c, expr2.val()),
                   atan2(expr1.val(), c),
                   atan2(c.fastAccessCoeff(j), expr2.val(j)),
                   atan2(expr1.val(j), c.fastAccessCoeff(j)),
                   (-c.fastAccessCoeff(j)*expr2.dx(i,j)) / (c.fastAccessCoeff(j)*c.fastAccessCoeff(j) + expr2.val(j)*expr2.val(j)),
                   (c.fastAccessCoeff(j)*expr1.dx(i,j))/ (expr1.val(j)*expr1.val(j) + c.fastAccessCoeff(j)*c.fastAccessCoeff(j)),
                   (-c.fastAccessCoeff(j)*expr2.fastAccessDx(i,j))/ (c.fastAccessCoeff(j)*c.fastAccessCoeff(j) + expr2.val(j)*expr2.val(j)),
                   (c.fastAccessCoeff(j)*expr1.fastAccessDx(i,j))/ (expr1.val(j)*expr1.val(j) + c.fastAccessCoeff(j)*c.fastAccessCoeff(j)))
// FAD_BINARYOP_MACRO(pow,
//                    PowerOp,
//                    using std::pow; using std::log;,
//                    pow(expr1.val(), expr2.val()),
//                    pow(expr1.val(j), expr2.val(j)),
//                    if_then_else( expr1.val(j) == val_type(0.0), val_type(0.0), val_type((expr2.dx(i,j)*log(expr1.val(j))+expr2.val(j)*expr1.dx(i,j)/expr1.val(j))*pow(expr1.val(j),expr2.val(j))) ),
//                    if_then_else( expr1.val(j) == val_type(0.0), val_type(0.0), val_type(expr2.dx(i,j)*log(expr1.val(j))*pow(expr1.val(j),expr2.val(j))) ),
//                    if_then_else( expr1.val(j) == val_type(0.0), val_type(0.0), val_type(expr2.val(j)*expr1.dx(i,j)/expr1.val(j)*pow(expr1.val(j),expr2.val(j))) ),
//                    if_then_else( expr1.val(j) == val_type(0.0), val_type(0.0), val_type((expr2.fastAccessDx(i,j)*log(expr1.val(j))+expr2.val(j)*expr1.fastAccessDx(i,j)/expr1.val(j))*pow(expr1.val(j),expr2.val(j))) ),
//                    pow(c, expr2.val()),
//                    pow(expr1.val(), c),
//                    pow(c.fastAccessCoeff(j), expr2.val(j)),
//                    pow(expr1.val(j), c.fastAccessCoeff(j)),
//                    if_then_else( c.fastAccessCoeff(j) == val_type(0.0), val_type(0.0), val_type(expr2.dx(i,j)*log(c.fastAccessCoeff(j))*pow(c.fastAccessCoeff(j),expr2.val(j))) ),
//                    if_then_else( expr1.val(j) == val_type(0.0), val_type(0.0), val_type(c.fastAccessCoeff(j)*expr1.dx(i,j)/expr1.val(j)*pow(expr1.val(j),c.fastAccessCoeff(j))) ),
//                    if_then_else( c.fastAccessCoeff(j) == val_type(0.0), val_type(0.0), val_type(expr2.fastAccessDx(i,j)*log(c.fastAccessCoeff(j))*pow(c.fastAccessCoeff(j),expr2.val(j))) ),
//                    if_then_else( expr1.val(j) == val_type(0.0), val_type(0.0), val_type(c.fastAccessCoeff(j)*expr1.fastAccessDx(i,j)/expr1.val(j)*pow(expr1.val(j),c.fastAccessCoeff(j)))) )
FAD_BINARYOP_MACRO(max,
                   MaxOp,
                   ;,
                   if_then_else( expr1.val() >= expr2.val(),  expr1.val(), expr2.val() ),
                   if_then_else( expr1.val(j) >= expr2.val(j),  expr1.val(j), expr2.val(j) ),
                   if_then_else( expr1.val(j) >= expr2.val(j), expr1.dx(i,j), expr2.dx(i,j) ),
                   if_then_else( expr1.val(j) >= expr2.val(j), val_type(0.0), expr2.dx(i,j) ),
                   if_then_else( expr1.val(j) >= expr2.val(j), expr1.dx(i,j), val_type(0.0) ),
                   if_then_else( expr1.val(j) >= expr2.val(j), expr1.fastAccessDx(i,j), expr2.fastAccessDx(i,j) ),
                   if_then_else( c >= expr2.val(), c,  expr2.val() ),
                   if_then_else( expr1.val() >= c, expr1.val(), c ),
                   if_then_else( c.fastAccessCoeff(j) >= expr2.val(j), c.fastAccessCoeff(j),  expr2.val(j) ),
                   if_then_else( expr1.val(j) >= c.fastAccessCoeff(j), expr1.val(j), c.fastAccessCoeff(j) ),
                   if_then_else( c.fastAccessCoeff(j) >= expr2.val(j), val_type(0.0),  expr2.dx(i,j) ),
                   if_then_else( expr1.val(j) >= c.fastAccessCoeff(j), expr1.dx(i,j), val_type(0.0) ),
                   if_then_else( c.fastAccessCoeff(j) >= expr2.val(j), val_type(0.0), expr2.fastAccessDx(i,j) ),
                   if_then_else( expr1.val(j) >= c.fastAccessCoeff(j), expr1.fastAccessDx(i,j), val_type(0.0) ) )
FAD_BINARYOP_MACRO(min,
                   MinOp,
                   ;,
                   if_then_else( expr1.val() <= expr2.val(), expr1.val(), expr2.val() ),
                   if_then_else( expr1.val(j) <= expr2.val(j), expr1.val(j), expr2.val(j) ),
                   if_then_else( expr1.val(j) <= expr2.val(j), expr1.dx(i,j), expr2.dx(i,j) ),
                   if_then_else( expr1.val(j) <= expr2.val(j), val_type(0.0), expr2.dx(i,j) ),
                   if_then_else( expr1.val(j) <= expr2.val(j), expr1.dx(i,j), val_type(0.0) ),
                   if_then_else( expr1.val(j) <= expr2.val(j), expr1.fastAccessDx(i,j), expr2.fastAccessDx(i,j) ),
                   if_then_else( c <= expr2.val(), c, expr2.val() ),
                   if_then_else( expr1.val() <= c, expr1.val(), c ),
                   if_then_else( c.fastAccessCoeff(j) <= expr2.val(j), c.fastAccessCoeff(j), expr2.val(j) ),
                   if_then_else( expr1.val(j) <= c.fastAccessCoeff(j), expr1.val(j), c.fastAccessCoeff(j) ),
                   if_then_else( c.fastAccessCoeff(j) <= expr2.val(j), val_type(0), expr2.dx(i,j) ),
                   if_then_else( expr1.val(j) <= c.fastAccessCoeff(j), expr1.dx(i,j), val_type(0) ),
                   if_then_else( c.fastAccessCoeff(j) <= expr2.val(j), val_type(0), expr2.fastAccessDx(i,j) ),
                   if_then_else( expr1.val(j) <= c.fastAccessCoeff(j), expr1.fastAccessDx(i,j), val_type(0) ) )

// Special handling for std::pow() to provide specializations of PowerOp for
// "simd" value types that use if_then_else(). The only reason for not using
// if_then_else() always is to avoid evaluating the derivative if the value is
// zero to avoid throwing FPEs.
namespace Sacado {
  namespace Fad {
  namespace Exp {

    //
    // Implementation for simd type using if_then_else()
    //
    template <typename T1, typename T2>
    class PowerOp< T1, T2, false, false, ExprSpecMPVector, PowerImpl::Simd > :
      public Expr< PowerOp< T1, T2, false, false, ExprSpecMPVector,
                            PowerImpl::Simd > > {
    public:

      typedef typename std::remove_cv<T1>::type ExprT1;
      typedef typename std::remove_cv<T2>::type ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename value_type::value_type val_type;

      typedef ExprSpecMPVector expr_spec_type;

      KOKKOS_INLINE_FUNCTION
      PowerOp(const T1& expr1_, const T2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      KOKKOS_INLINE_FUNCTION
      int size() const {
        const int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      KOKKOS_INLINE_FUNCTION
      value_type val() const {
        using std::pow;
        return pow(expr1.val(), expr2.val());
      }

      KOKKOS_INLINE_FUNCTION
      val_type val(int j) const {
        using std::pow;
        return pow(expr1.val(j), expr2.val(j));
      }

      KOKKOS_INLINE_FUNCTION
      val_type dx(int i, int j) const {
        using std::pow; using std::log;
        const int sz1 = expr1.size(), sz2 = expr2.size();
        if (sz1 > 0 && sz2 > 0)
          return if_then_else( expr1.val(j) == val_type(0.0), val_type(0.0), val_type((expr2.dx(i,j)*log(expr1.val(j))+expr2.val(j)*expr1.dx(i,j)/expr1.val(j))*pow(expr1.val(j),expr2.val(j))) );
        else if (sz1 > 0)
          // Don't use formula (a(x)^b)' = b*a(x)^{b-1}*a'(x)
          // It seems less accurate and caused convergence problems in some codes
          return if_then_else( expr2.val(j) == scalar_type(1.0), expr1.dx(i,j), if_then_else( expr1.val(j) == val_type(0.0), val_type(0.0), val_type(expr2.val(j)*expr1.dx(i,j)/expr1.val(j)*pow(expr1.val(j),expr2.val(j))) ));
        else
          return if_then_else( expr1.val(j) == val_type(0.0), val_type(0.0), val_type(expr2.dx(i,j)*log(expr1.val(j))*pow(expr1.val(j),expr2.val(j))) );
      }

      KOKKOS_INLINE_FUNCTION
      val_type fastAccessDx(int i, int j) const {
        using std::pow; using std::log;
        return if_then_else( expr1.val(j) == val_type(0.0), val_type(0.0), val_type((expr2.fastAccessDx(i,j)*log(expr1.val(j))+expr2.val(j)*expr1.fastAccessDx(i,j)/expr1.val(j))*pow(expr1.val(j),expr2.val(j))) );
      }

    protected:

      const T1& expr1;
      const T2& expr2;

    };

    template <typename T1, typename T2>
    class PowerOp< T1, T2, false, true, ExprSpecMPVector, PowerImpl::Simd >
      : public Expr< PowerOp< T1, T2, false, true, ExprSpecMPVector,
                              PowerImpl::Simd > > {
    public:

      typedef typename std::remove_cv<T1>::type ExprT1;
      typedef T2 ConstT;
      typedef typename ExprT1::value_type value_type;
      typedef typename ExprT1::scalar_type scalar_type;

      typedef typename value_type::value_type val_type;

      typedef ExprSpecMPVector expr_spec_type;

      KOKKOS_INLINE_FUNCTION
      PowerOp(const T1& expr1_, const ConstT& c_) :
        expr1(expr1_), c(c_) {}

      KOKKOS_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      KOKKOS_INLINE_FUNCTION
      value_type val() const {
        using std::pow;
        return pow(expr1.val(), c);
      }

      KOKKOS_INLINE_FUNCTION
      val_type val(int j) const {
        using std::pow;
        return pow(expr1.val(j), c.fastAccessCoeff(j));
      }

      KOKKOS_INLINE_FUNCTION
      val_type dx(int i, int j) const {
        using std::pow;
        // Don't use formula (a(x)^b)' = b*a(x)^{b-1}*a'(x)
        // It seems less accurate and caused convergence problems in some codes
        return if_then_else( c.fastAccessCoeff(j) == scalar_type(1.0), expr1.dx(i,j), if_then_else( expr1.val(j) == val_type(0.0), val_type(0.0), val_type(c.fastAccessCoeff(j)*expr1.dx(i,j)/expr1.val(j)*pow(expr1.val(j),c.fastAccessCoeff(j))) ));
      }

      KOKKOS_INLINE_FUNCTION
      val_type fastAccessDx(int i, int j) const {
        using std::pow;
        // Don't use formula (a(x)^b)' = b*a(x)^{b-1}*a'(x)
        // It seems less accurate and caused convergence problems in some codes
        return if_then_else( c.fastAccessCoeff(j) == scalar_type(1.0), expr1.fastAccessDx(i,j), if_then_else( expr1.val(j) == val_type(0.0), val_type(0.0), val_type(c.fastAccessCoeff(j)*expr1.fastAccessDx(i,j)/expr1.val(j)*pow(expr1.val(j),c.fastAccessCoeff(j))) ));
      }

    protected:

      const T1& expr1;
      const ConstT& c;
    };

    template <typename T1, typename T2>
    class PowerOp< T1, T2, true, false, ExprSpecMPVector, PowerImpl::Simd >
      : public Expr< PowerOp< T1, T2, true, false, ExprSpecMPVector,
                              PowerImpl::Simd> > {
    public:

      typedef typename std::remove_cv<T2>::type ExprT2;
      typedef T1 ConstT;
      typedef typename ExprT2::value_type value_type;
      typedef typename ExprT2::scalar_type scalar_type;

      typedef typename value_type::value_type val_type;

      typedef ExprSpecMPVector expr_spec_type;

      KOKKOS_INLINE_FUNCTION
      PowerOp(const ConstT& c_, const T2& expr2_) :
        c(c_), expr2(expr2_) {}

      KOKKOS_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      KOKKOS_INLINE_FUNCTION
      value_type val() const {
        using std::pow;
        return pow(c, expr2.val());
      }

      KOKKOS_INLINE_FUNCTION
      val_type val(int j) const {
        using std::pow;
        return pow(c.fastAccessCoeff(j), expr2.val(j));
      }

      KOKKOS_INLINE_FUNCTION
      val_type dx(int i, int j) const {
        using std::pow; using std::log;
        return if_then_else( c.fastAccessCoeff(j) == val_type(0.0), val_type(0.0), val_type(expr2.dx(i,j)*log(c.fastAccessCoeff(j))*pow(c.fastAccessCoeff(j),expr2.val(j))) );
      }

      KOKKOS_INLINE_FUNCTION
      val_type fastAccessDx(int i, int j) const {
        using std::pow; using std::log;
        return if_then_else( c.fastAccessCoeff(j) == val_type(0.0), val_type(0.0), val_type(expr2.fastAccessDx(i,j)*log(c.fastAccessCoeff(j))*pow(c.fastAccessCoeff(j),expr2.val(j))) );
      }

    protected:

      const ConstT& c;
      const T2& expr2;
    };

    //
    // Specialization for nested derivatives.  This version does not use
    // if_then_else/ternary-operator on the base so that nested derivatives
    // are correct.
    //
    template <typename T1, typename T2>
    class PowerOp< T1, T2, false, false, ExprSpecMPVector,
                   PowerImpl::NestedSimd > :
      public Expr< PowerOp< T1, T2, false, false, ExprSpecMPVector,
                   PowerImpl::NestedSimd > > {
    public:

      typedef typename std::remove_cv<T1>::type ExprT1;
      typedef typename std::remove_cv<T2>::type ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename value_type::value_type val_type;

      typedef ExprSpecMPVector expr_spec_type;

      KOKKOS_INLINE_FUNCTION
      PowerOp(const T1& expr1_, const T2& expr2_) :
        expr1(expr1_), expr2(expr2_) {}

      KOKKOS_INLINE_FUNCTION
      int size() const {
        const int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      KOKKOS_INLINE_FUNCTION
      value_type val() const {
        using std::pow;
        return pow(expr1.val(), expr2.val());
      }

      KOKKOS_INLINE_FUNCTION
      val_type val(int j) const {
        using std::pow;
        return pow(expr1.val(j), expr2.val(j));
      }

      KOKKOS_INLINE_FUNCTION
      value_type dx(int i, int j) const {
        using std::pow; using std::log;
        const int sz1 = expr1.size(), sz2 = expr2.size();
        if (sz1 > 0 && sz2 > 0)
          return (expr2.dx(i,j)*log(expr1.val(j))+expr2.val(j)*expr1.dx(i,j)/expr1.val(j))*pow(expr1.val(j),expr2.val(j));
        else if (sz1 > 0)
          return if_then_else( expr2.val(j) == scalar_type(0.0), value_type(0.0), value_type((expr2.val(j)*expr1.dx(i,j))*pow(expr1.val(j),expr2.val(j)-scalar_type(1.0))));
        else
          return expr2.dx(i,j)*log(expr1.val(j))*pow(expr1.val(j),expr2.val(j));
      }

      KOKKOS_INLINE_FUNCTION
      value_type fastAccessDx(int i, int j) const {
        using std::pow; using std::log;
        return (expr2.fastAccessDx(i,j)*log(expr1.val(j))+expr2.val(j)*expr1.fastAccessDx(i,j)/expr1.val(j))*pow(expr1.val(j),expr2.val(j));
      }

    protected:

      const T1& expr1;
      const T2& expr2;

    };

    template <typename T1, typename T2>
    class PowerOp< T1, T2, false, true, ExprSpecMPVector,
                   PowerImpl::NestedSimd > :
      public Expr< PowerOp< T1, T2, false, true, ExprSpecMPVector,
                   PowerImpl::NestedSimd > > {
    public:

      typedef typename std::remove_cv<T1>::type ExprT1;
      typedef T2 ConstT;
      typedef typename ExprT1::value_type value_type;
      typedef typename ExprT1::scalar_type scalar_type;

      typedef typename value_type::value_type val_type;

      typedef ExprSpecMPVector expr_spec_type;

      KOKKOS_INLINE_FUNCTION
      PowerOp(const T1& expr1_, const ConstT& c_) :
        expr1(expr1_), c(c_) {}

      KOKKOS_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      KOKKOS_INLINE_FUNCTION
      value_type val() const {
        using std::pow;
        return pow(expr1.val(), c);
      }

      KOKKOS_INLINE_FUNCTION
      val_type val(int j) const {
        using std::pow;
        return pow(expr1.val(j), c.fastAccessCoeff(j));
      }

      KOKKOS_INLINE_FUNCTION
      value_type dx(int i, int j) const {
        using std::pow;
        return if_then_else( c.fastAccessCoeff(j) == scalar_type(0.0), value_type(0.0), value_type(c.fastAccessCoeff(j)*expr1.dx(i,j)*pow(expr1.val(j),c.fastAccessCoeff(j)-scalar_type(1.0))));
      }

      KOKKOS_INLINE_FUNCTION
      value_type fastAccessDx(int i, int j) const {
        using std::pow;
        return if_then_else( c.fastAccessCoeff(j) == scalar_type(0.0), value_type(0.0), value_type(c.fastAccessCoeff(j)*expr1.fastAccessDx(i,j)*pow(expr1.val(j),c.fastAccessCoeff(j)-scalar_type(1.0))));
      }

    protected:

      const T1& expr1;
      const ConstT& c;
    };

    template <typename T1, typename T2>
    class PowerOp<T1, T2, true, false, ExprSpecMPVector,
                   PowerImpl::NestedSimd > :
      public Expr< PowerOp< T1, T2, true, false, ExprSpecMPVector,
                   PowerImpl::NestedSimd > > {
    public:

      typedef typename std::remove_cv<T2>::type ExprT2;
      typedef T1 ConstT;
      typedef typename ExprT2::value_type value_type;
      typedef typename ExprT2::scalar_type scalar_type;

      typedef typename value_type::value_type val_type;

      typedef ExprSpecMPVector expr_spec_type;

      KOKKOS_INLINE_FUNCTION
      PowerOp(const ConstT& c_, const T2& expr2_) :
        c(c_), expr2(expr2_) {}

      KOKKOS_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      KOKKOS_INLINE_FUNCTION
      value_type val() const {
        using std::pow;
        return pow(c, expr2.val());
      }

      KOKKOS_INLINE_FUNCTION
      val_type val(int j) const {
        using std::pow;
        return pow(c.fastAccessCoeff(j), expr2.val(j));
      }

      KOKKOS_INLINE_FUNCTION
      value_type dx(int i, int j) const {
        using std::pow; using std::log;
        return expr2.dx(i,j)*log(c.fastAccessCoeff(j))*pow(c.fastAccessCoeff(j),expr2.val(j));
      }

      KOKKOS_INLINE_FUNCTION
      value_type fastAccessDx(int i, int j) const {
        using std::pow; using std::log;
        return expr2.fastAccessDx(i,j)*log(c.fastAccessCoeff(j))*pow(c.fastAccessCoeff(j),expr2.val(j));
      }

    protected:

      const ConstT& c;
      const T2& expr2;
    };

  }
  }
}

//--------------------------if_then_else operator -----------------------
// Can't use the above macros because it is a ternary operator (sort of).
// Also, relies on C++11

namespace Sacado {
  namespace Fad {
  namespace Exp {

    template <typename CondT, typename T1, typename T2>
    class IfThenElseOp< CondT,T1,T2,false,false,ExprSpecMPVector > :
      public Expr< IfThenElseOp< CondT, T1, T2, false, false, ExprSpecMPVector > > {

    public:

      typedef typename std::remove_cv<T1>::type ExprT1;
      typedef typename std::remove_cv<T2>::type ExprT2;
      typedef typename ExprT1::value_type value_type_1;
      typedef typename ExprT2::value_type value_type_2;
      typedef typename Sacado::Promote<value_type_1,
                                       value_type_2>::type value_type;

      typedef typename ExprT1::scalar_type scalar_type_1;
      typedef typename ExprT2::scalar_type scalar_type_2;
      typedef typename Sacado::Promote<scalar_type_1,
                                       scalar_type_2>::type scalar_type;

      typedef typename value_type::value_type val_type;

      typedef ExprSpecMPVector expr_spec_type;

      KOKKOS_INLINE_FUNCTION
      IfThenElseOp(const CondT& cond_, const T1& expr1_, const T2& expr2_) :
        cond(cond_), expr1(expr1_), expr2(expr2_) {}

      KOKKOS_INLINE_FUNCTION
      int size() const {
        int sz1 = expr1.size(), sz2 = expr2.size();
        return sz1 > sz2 ? sz1 : sz2;
      }

      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess() && expr2.hasFastAccess();
      }

      KOKKOS_INLINE_FUNCTION
      value_type val() const {
        return if_then_else( cond, expr1.val(), expr2.val() );
      }

      KOKKOS_INLINE_FUNCTION
      val_type val(int j) const {
        return if_then_else( cond, expr1.val(j), expr2.val(j) );
      }

      KOKKOS_INLINE_FUNCTION
      val_type dx(int i, int j) const {
        return if_then_else( cond, expr1.dx(i,j), expr2.dx(i,j) );
      }

      KOKKOS_INLINE_FUNCTION
      val_type fastAccessDx(int i, int j) const {
        return if_then_else( cond, expr1.fastAccessDx(i,j), expr2.fastAccessDx(i,j) );
      }

    protected:

      const CondT&  cond;
      const T1& expr1;
      const T2& expr2;

    };

    template <typename CondT, typename T1, typename T2>
    class IfThenElseOp< CondT, T1, T2, false, true, ExprSpecMPVector> :
      public Expr< IfThenElseOp< CondT, T1, T2, false, true, ExprSpecMPVector > > {

    public:

      typedef typename std::remove_cv<T1>::type ExprT1;
      typedef T2 ConstT;
      typedef typename ExprT1::value_type value_type;
      typedef typename ExprT1::scalar_type scalar_type;

      typedef typename value_type::value_type val_type;

      KOKKOS_INLINE_FUNCTION
      IfThenElseOp(const CondT& cond_, const T1& expr1_, const ConstT& c_) :
        cond(cond_), expr1(expr1_), c(c_) {}

      KOKKOS_INLINE_FUNCTION
      int size() const {
        return expr1.size();
      }

      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr1.hasFastAccess();
      }

      KOKKOS_INLINE_FUNCTION
      value_type val() const {
        return if_then_else( cond, expr1.val(), c );
      }

      KOKKOS_INLINE_FUNCTION
      val_type val(int j) const {
        return if_then_else( cond, expr1.val(j), c.fastAccessCoeff(j) );
      }

      KOKKOS_INLINE_FUNCTION
      val_type dx(int i, int j) const {
        return if_then_else( cond, expr1.dx(i,j), val_type(0.0) );
      }

      KOKKOS_INLINE_FUNCTION
      val_type fastAccessDx(int i, int j) const {
        return if_then_else( cond, expr1.fastAccessDx(i,j), val_type(0.0) );
      }

    protected:

      const CondT&  cond;
      const T1& expr1;
      const typename ConstTypeRef<ConstT,value_type>::type c;
    };

    template <typename CondT, typename T1, typename T2>
    class IfThenElseOp< CondT, T1, T2, true, false, ExprSpecMPVector> :
      public Expr< IfThenElseOp< CondT, T1, T2, true, false, ExprSpecMPVector > > {

    public:

      typedef typename std::remove_cv<T2>::type ExprT2;
      typedef T1 ConstT;
      typedef typename ExprT2::value_type value_type;
      typedef typename ExprT2::scalar_type scalar_type;

      typedef typename value_type::value_type val_type;

      typedef ExprSpecMPVector expr_spec_type;

      KOKKOS_INLINE_FUNCTION
      IfThenElseOp(const CondT& cond_, const ConstT& c_, const T2& expr2_) :
        cond(cond_), c(c_), expr2(expr2_) {}

      KOKKOS_INLINE_FUNCTION
      int size() const {
        return expr2.size();
      }

      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess() const {
        return expr2.hasFastAccess();
      }

      KOKKOS_INLINE_FUNCTION
      value_type val() const {
        return if_then_else( cond, c, expr2.val() );
      }

      KOKKOS_INLINE_FUNCTION
      val_type val(int j) const {
        return if_then_else( cond, c.fastAccessCoeff(j), expr2.val(j) );
      }

      KOKKOS_INLINE_FUNCTION
      val_type dx(int i, int j) const {
        return if_then_else( cond, val_type(0.0), expr2.dx(i,j) );
      }

      KOKKOS_INLINE_FUNCTION
      val_type fastAccessDx(int i, int j) const {
        return if_then_else( cond, val_type(0.0), expr2.fastAccessDx(i,j) );
      }

    protected:

      const CondT&  cond;
      const typename ConstTypeRef<ConstT,value_type>::type c;
      const T2& expr2;
    };

  }
  }
}

#undef FAD_BINARYOP_MACRO

#endif // SACADO_FAD_EXP_MP_VECTOR_HPP
