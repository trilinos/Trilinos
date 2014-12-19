// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
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
//  A short implementation ( not all operators and
//  functions are overloaded ) of 1st order Automatic
//  Differentiation in forward mode (FAD) using
//  EXPRESSION TEMPLATES.
//
//********************************************************
// @HEADER

#ifndef SACADO_CACHEFAD_SFAD_HPP
#define SACADO_CACHEFAD_SFAD_HPP

#include "Sacado_CacheFad_SFadTraits.hpp"
#include "Sacado_CacheFad_Expression.hpp"
#include "Sacado_StaticArrayTraits.hpp"

namespace Sacado {

  //! Namespace for forward-mode AD classes
  namespace CacheFad {

    //! A tag for specializing Expr for SFad expressions
    template <typename T, int Num>
    struct SFadExprTag {};

    // Forward declaration
    template <typename T, int Num> class SFad;

    /*!
     * \brief Expression template forward-mode AD class with static memory
     * allocation.
     */
    /*!
     * This classes specializes Expr to SFad expressions.
     */
    template <typename T, int Num>
    class Expr< SFadExprTag<T,Num> > {

    public:

      //! Typename of values
      typedef typename RemoveConst<T>::type value_type;

      //! Typename of scalar's (which may be different from T)
      typedef typename ScalarType<value_type>::type scalar_type;

      //! Typename of base-expressions
      typedef SFad<T,Num> base_expr_type;

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor
      KOKKOS_INLINE_FUNCTION
      Expr() : val_( T(0.)), update_val_(true) { ss_array<T>::zero(dx_, Num); }

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Expr(const S & x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        val_(x), update_val_(true) {
        ss_array<T>::zero(dx_, Num);
      }

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      KOKKOS_INLINE_FUNCTION
      Expr(const int sz, const T & x)  : val_(x), update_val_(true) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (sz != Num)
          throw "CacheFad::SFad() Error:  Supplied derivative dimension does not match compile time length.";
#endif

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
        val_(x), update_val_(true) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (sz != Num)
          throw "CacheFad::SFad() Error:  Supplied derivative dimension does not match compile time length.";
        if (i >= Num)
          throw "CacheFad::SFad() Error:  Invalid derivative index.";
#endif

        ss_array<T>::zero(dx_, Num);
        dx_[i]=1.;
      }

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      Expr(const Expr& x)  :
        val_(x.val()), update_val_(x.update_val_) {
        for (int i=0; i<Num; i++)
          dx_[i] = x.dx_[i];
      }

      //! Copy constructor from any Expression object
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      Expr(const Expr<S>& x, SACADO_ENABLE_EXPR_CTOR_DECL)  :
        update_val_(x.updateValue()) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (x.size() != Num)
          throw "CacheFad::SFad() Error:  Attempt to assign with incompatible sizes";
#endif

        x.cache();

        if (update_val_)
          this->val() = x.val();
        else
          this->val() = T(0.);

        for(int i=0; i<Num; ++i)
          dx_[i] = x.fastAccessDx(i);
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
          throw "CacheFad::diff() Error:  Supplied derivative dimension does not match compile time length.";
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
          throw "CacheFad::resize() Error:  Cannot resize fixed derivative array dimension";
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
      void setUpdateValue(bool update_val) { update_val_ = update_val; }

      //! Return whether this Fad object has an updated value
      KOKKOS_INLINE_FUNCTION
      bool updateValue() const { return update_val_; }

      //! Cache values
      KOKKOS_INLINE_FUNCTION
      void cache() const {}

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
        // Copy value
        val_ = x.val_;

        // Copy dx_
        for (int i=0; i<Num; i++)
          dx_[i] = x.dx_[i];

        update_val_ = x.update_val_;

        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(Expr&) operator=(const Expr<S>& x) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (x.size() != Num)
          throw "CacheFad::operator=() Error:  Attempt to assign with incompatible sizes";
#endif

        x.cache();

        for(int i=0; i<Num; ++i)
          dx_[i] = x.fastAccessDx(i);

        update_val_ = x.updateValue();
        if (update_val_)
          val_ = x.val();

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
        if (update_val_) this->val() += v;
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(Expr&) operator -= (const S& v) {
        if (update_val_) this->val() -= v;
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(Expr&) operator *= (const S& v) {
        if (update_val_) this->val() *= v;
        for (int i=0; i<Num; ++i)
          dx_[i] *= v;
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(Expr&) operator /= (const S& v) {
        if (update_val_) this->val() /= v;
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
          throw "CacheFad::operator+=() Error:  Attempt to assign with incompatible sizes";
#endif

        x.cache();

        for (int i=0; i<Num; ++i)
          dx_[i] += x.fastAccessDx(i);

        update_val_ = x.updateValue();
        if (update_val_)
          val_ += x.val();

        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(Expr&) operator -= (const Expr<S>& x) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (x.size() != Num)
          throw "CacheFad::operator-=() Error:  Attempt to assign with incompatible sizes";
#endif

        x.cache();

        for(int i=0; i<Num; ++i)
          dx_[i] -= x.fastAccessDx(i);

        update_val_ = x.updateValue();
        if (update_val_)
          val_ -= x.val();

        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(Expr&) operator *= (const Expr<S>& x) {
        x.cache();

        T xval = x.val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (x.size() != Num)
          throw "CacheFad::operator*=() Error:  Attempt to assign with incompatible sizes";
#endif

        for(int i=0; i<Num; ++i)
          dx_[i] = val_ * x.fastAccessDx(i) + dx_[i] * xval;

        update_val_ = x.updateValue();
        if (update_val_)
          val_ *= xval;

        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(Expr&) operator /= (const Expr<S>& x) {
        x.cache();

        T xval = x.val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (x.size() != Num)
          throw "CacheFad::operator/=() Error:  Attempt to assign with incompatible sizes";
#endif

        for(int i=0; i<Num; ++i)
          dx_[i] = ( dx_[i]*xval - val_*x.fastAccessDx(i) )/ (xval*xval);

        update_val_ = x.updateValue();
        if (update_val_)
          val_ /= xval;

        return *this;
      }

      //@}

    protected:

      //! Value
      T val_;

      //! Derivatives
      T dx_[Num];

      //! Update value
      bool update_val_;

    }; // class Expr<SFadExprTag>

  } // namespace CacheFad

} // namespace Sacado

#define FAD_NS CacheFad
#include "Sacado_Fad_SFad_tmpl.hpp"
#undef FAD_NS

#include "Sacado_CacheFad_ViewFad.hpp"
#include "Sacado_CacheFad_Ops.hpp"

#endif // SACADO_CACHEFAD_SFAD_HPP
