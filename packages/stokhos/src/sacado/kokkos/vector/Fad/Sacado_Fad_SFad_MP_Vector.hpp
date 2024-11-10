// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
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
//  A short implementation ( not all operators and
//  functions are overloaded ) of 1st order Automatic
//  Differentiation in forward mode (FAD) using
//  EXPRESSION TEMPLATES.
//
//********************************************************
// @HEADER

#ifndef SACADO_FAD_SFAD_MP_VECTOR_HPP
#define SACADO_FAD_SFAD_MP_VECTOR_HPP

#include "Sacado_Fad_SFad.hpp"
#include "Sacado_Fad_ExprSpec_MP_Vector.hpp"

namespace Stokhos {
  template <typename Ord, typename Val, int Num, typename Dev>
  class StaticFixedStorage;
}

namespace Sacado {

  namespace MP {
    template <typename S> class Vector;
  }

  namespace Fad {

    template <typename Ord, typename Val, int VecNum, typename Dev, int Num>
    struct ExprSpec< SFadExprTag< Sacado::MP::Vector< Stokhos::StaticFixedStorage<Ord,Val,VecNum,Dev> >, Num > > {
      typedef ExprSpecMPVector type;
    };
    template <typename Ord, typename Val, int VecNum, typename Dev, int Num>
    struct ExprSpec< SFad< Sacado::MP::Vector< Stokhos::StaticFixedStorage<Ord,Val,VecNum,Dev> >, Num > > {
      typedef ExprSpecMPVector type;
    };

    /*!
     * \brief Expression template forward-mode AD class with static memory
     * allocation.
     */
    /*!
     * This classes specializes Expr to SFad expressions.
     */
    template <typename Ord, typename Val, int VecNum, typename Dev, int Num>
    class Expr< SFadExprTag< Sacado::MP::Vector< Stokhos::StaticFixedStorage<Ord,Val,VecNum,Dev> >,Num >,ExprSpecMPVector > {

    public:

      typedef Sacado::MP::Vector< Stokhos::StaticFixedStorage<Ord,Val,VecNum,Dev> > T;

      //! Typename of values
      typedef typename RemoveConst<T>::type value_type;

      //! Typename of scalar's (which may be different from T)
      typedef typename ScalarType<value_type>::type scalar_type;

      //! Typename of base-expressions
      typedef SFad<value_type,Num> base_expr_type;

      typedef typename value_type::value_type val_type;

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor
      KOKKOS_INLINE_FUNCTION
      Expr() : val_( T(0.)) {
        ss_array<T>::zero(dx_, Num); }

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Expr(const S & x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        val_(x) {
        ss_array<T>::zero(dx_, Num);
      }

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      KOKKOS_INLINE_FUNCTION
      Expr(const int sz, const T & x, const DerivInit zero_out = InitDerivArray)  : val_(x) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (sz != Num)
          throw "SFad::SFad() Error:  Supplied derivative dimension does not match compile time length.";
#endif
        if (zero_out == InitDerivArray)
          ss_array<T>::zero(dx_, Num);
      }

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      KOKKOS_INLINE_FUNCTION
      Expr(const int sz, const int i, const T & x) :
        val_(x) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (sz != Num)
          throw "SFad::SFad() Error:  Supplied derivative dimension does not match compile time length.";
        if (i >= Num)
          throw "SFad::SFad() Error:  Invalid derivative index.";
#endif

        ss_array<T>::zero(dx_, Num);
        dx_[i]=1.;
      }

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      Expr(const Expr& x) :
        val_(x.val()) {
        for (int i=0; i<Num; i++)
          dx_[i] = x.dx_[i];
      }

      //! Copy constructor from any Expression object
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Expr(const Expr<S>& x, SACADO_ENABLE_EXPR_CTOR_DECL) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (x.size() != Num)
          throw "SFad::SFad() Error:  Attempt to assign with incompatible sizes";
#endif

        for(int i=0; i<Num; ++i) {
          for (int j=0; j<VecNum; ++j)
            dx_[i].fastAccessCoeff(j) = x.fastAccessDx(i,j);
        }

        for (int j=0; j<VecNum; ++j)
          this->val(j) = x.val(j);
      }

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~Expr() {}

      //! Set %Fad object as the \c ith independent variable
      /*!
       * Sets the derivative array of length \c n to the \c ith row of the
       * identity matrix and has the same affect as the
       * Implementation(const int sz, const int i, const T & x)
       * constructor.
       */
      KOKKOS_INLINE_FUNCTION
      void diff(const int ith, const int n) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (n != Num)
          throw "SFad::diff() Error:  Supplied derivative dimension does not match compile time length.";
#endif

        ss_array<T>::zero(dx_, Num);
        dx_[ith] = T(1.);
      }

      //! Resize derivative array to length \c sz
      /*!
       * Since the derivative array length is not dynamic, this method
       * throws an error if compiled with SACADO_DEBUG defined.
       */
      KOKKOS_INLINE_FUNCTION
      void resize(int sz) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (sz != Num)
          throw "SFad::resize() Error:  Cannot resize fixed derivative array dimension";
#endif
      }

      //! Expand derivative array to size sz
      /*!
       * Since the derivative array length is not dynamic, this method
       * throws an error if compiled with SACADO_DEBUG defined.
       */
      KOKKOS_INLINE_FUNCTION
      void expand(int sz) { resize(sz); }

      //! Zero out the derivative array
      KOKKOS_INLINE_FUNCTION
      void zero() { ss_array<T>::zero(dx_, Num); }

      //! Set whether this Fad object should update values
      KOKKOS_INLINE_FUNCTION
      void setUpdateValue(bool update_val) {  }

      //! Return whether this Fad object has an updated value
      KOKKOS_INLINE_FUNCTION
      bool updateValue() const { return true; }

      //! Returns whether two Fad objects have the same values
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(bool) isEqualTo(const Expr<S>& x) const {
        typedef IsEqual<value_type> IE;
        if (x.size() != this->size()) return false;
        bool eq = IE::eval(x.val(), this->val());
        for (int i=0; i<this->size(); i++)
          eq = eq && IE::eval(x.dx(i), this->dx(i));
        return eq;
      }

      //@}

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      const T& val() const { return val_;}

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      T& val() { return val_;}

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      const val_type& val(int j) const { return val_.fastAccessCoeff(j);}

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      val_type& val(int j) { return val_.fastAccessCoeff(j);}

      //@}

      /*!
       * @name Derivative accessor methods
       */
      //@{

      //! Returns number of derivative components
      KOKKOS_INLINE_FUNCTION
      int size() const { return Num;}

      /*!
       * \brief Returns number of derivative components that can be stored
       * without reallocation
       */
      KOKKOS_INLINE_FUNCTION
      int availableSize() const { return Num; }

      //! Returns true if derivative array is not empty
      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess() const { return true; }

      //! Returns true if derivative array is empty
      KOKKOS_INLINE_FUNCTION
      bool isPassive() const { return false; }

      //! Set whether variable is constant
      KOKKOS_INLINE_FUNCTION
      void setIsConstant(bool is_const) {}

      //! Returns derivative array
      KOKKOS_INLINE_FUNCTION
      const T* dx() const { return &(dx_[0]);}

      //! Returns derivative component \c i with bounds checking
      KOKKOS_INLINE_FUNCTION
      const T& dx(int i) const { return dx_[i]; }

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      T& fastAccessDx(int i) { return dx_[i];}

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      const T& fastAccessDx(int i) const { return dx_[i];}

      //! Returns derivative component \c i with bounds checking
      KOKKOS_INLINE_FUNCTION
      const val_type& dx(int i, int j) const { return dx_[i].fastAccessCoeff(j); }

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      val_type& fastAccessDx(int i, int j) { return dx_[i].fastAccessCoeff(j);}

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      const val_type& fastAccessDx(int i, int j) const { return dx_[i].fastAccessCoeff(j);}

      //@}

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(Expr&) operator=(const S& v) {
        val_ = v;
        ss_array<T>::zero(dx_, Num);
        return *this;
      }

      //! Assignment with Expr right-hand-side
      KOKKOS_INLINE_FUNCTION
      Expr& operator=(const Expr& x) {
        if (this != &x) {
          // Copy value
          val_ = x.val_;

          // Copy dx_
          for (int i=0; i<Num; i++)
            dx_[i] = x.dx_[i];
        }
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(Expr&) operator=(const Expr<S>& x) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (x.size() != Num)
          throw "SFad::operator=() Error:  Attempt to assign with incompatible sizes";
#endif

       for(int i=0; i<Num; ++i) {
          for (int j=0; j<VecNum; ++j)
            dx_[i].fastAccessCoeff(j) = x.fastAccessDx(i,j);
        }

        for (int j=0; j<VecNum; ++j)
          this->val(j) = x.val(j);

        return *this;
      }

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Addition-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(Expr&) operator += (const S& v) {
        this->val() += v;
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(Expr&) operator -= (const S& v) {
        this->val() -= v;
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(Expr&) operator *= (const S& v) {
        this->val() *= v;
        for (int i=0; i<Num; ++i)
          dx_[i] *= v;
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(Expr&) operator /= (const S& v) {
        this->val() /= v;
        for (int i=0; i<Num; ++i)
          dx_[i] /= v;
        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(Expr&) operator += (const Expr<S>& x) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (x.size() != Num)
          throw "SFad::operator+=() Error:  Attempt to assign with incompatible sizes";
#endif

        for(int i=0; i<Num; ++i) {
          for (int j=0; j<VecNum; ++j)
            dx_[i].fastAccessCoeff(j) += x.fastAccessDx(i,j);
        }

        for (int j=0; j<VecNum; ++j)
          this->val(j) += x.val(j);

        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(Expr&) operator -= (const Expr<S>& x) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (x.size() != Num)
          throw "SFad::operator-=() Error:  Attempt to assign with incompatible sizes";
#endif

        for(int i=0; i<Num; ++i) {
          for (int j=0; j<VecNum; ++j)
            dx_[i].fastAccessCoeff(j) -= x.fastAccessDx(i,j);
        }

        for (int j=0; j<VecNum; ++j)
          this->val(j) -= x.val(j);

        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(Expr&) operator *= (const Expr<S>& x) {
        T xval;
        for (int j=0; j<VecNum; ++j)
          xval.fastAccessCoeff(j) = x.val(j);

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (x.size() != Num)
          throw "SFad::operator*=() Error:  Attempt to assign with incompatible sizes";
#endif

        for(int i=0; i<Num; ++i)
          for (int j=0; j<VecNum; ++j)
            dx_[i].fastAccessCoeff(j) = val_.fastAccessCoeff(j) * x.fastAccessDx(i,j) + dx_[i] * xval.fastAccessCoeff(j);

        val_ *= xval;

        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(Expr&) operator /= (const Expr<S>& x) {
        T xval;
        for (int j=0; j<VecNum; ++j)
          xval.fastAccessCoeff(j) = x.val(j);
        T xval2 = xval*xval;

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (x.size() != Num)
          throw "SFad::operator/=() Error:  Attempt to assign with incompatible sizes";
#endif

        for(int i=0; i<Num; ++i)
          for (int j=0; j<VecNum; ++j)
            dx_[i].fastAccessCoeff(j) = ( dx_[i].fastAccessCoeff(j)*xval.fastAccessCoeff(j) - val_.fastAccessCoeff(j)*x.fastAccessDx(i,j) )/ (xval2.fastAccessCoeff(j));

        val_ /= xval;

        return *this;
      }

      //@}

    protected:

      //! Derivatives
      T dx_[Num];

      //! Value
      T val_;

    }; // class Expr<SFadExprTag>

  } // namespace Fad

} // namespace Sacado

#include "Sacado_Fad_Ops_MP_Vector.hpp"

#endif // SACADO_FAD_SFAD_MP_VECTOR_HPP
