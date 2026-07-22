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
//  A short implementation ( not all operators and
//  functions are overloaded ) of 1st order Automatic
//  Differentiation in forward mode (FAD) using
//  EXPRESSION TEMPLATES.
//
//********************************************************
// @HEADER

#ifndef SACADO_CACHEFAD_GENERALFAD_HPP
#define SACADO_CACHEFAD_GENERALFAD_HPP

#include "Sacado_CacheFad_Expression.hpp"

namespace Sacado {

  //! Namespace for forward-mode AD classes w/caching
  namespace CacheFad {

    //! Forward-mode AD class templated on the storage for the derivative array
    /*!
     * This class provides a general forward mode AD implementation for any
     * type of derivative array storage.  It does not incorporate expression
     * templates.
     *
     * This version of GeneralFad is virtually identical to
     * Sacado::Fad::GeneralFad, with small differences in the implementations
     * of some of the methods to support caching of "values" of intermediate
     * operations in an expression, recognizing that the "this" object might
     * also be on the right-hand-side of an expression.
     */
    template <typename T, typename Storage>
    class GeneralFad : public Storage {

    public:

      //! Typename of values
      typedef typename RemoveConst<T>::type value_type;

      //! Typename of scalar's (which may be different from T)
      typedef typename ScalarType<value_type>::type scalar_type;

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor
      SACADO_INLINE_FUNCTION
      GeneralFad() : Storage(T(0.)) {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      template <typename S>
      SACADO_INLINE_FUNCTION
      GeneralFad(const S& x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        Storage(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      SACADO_INLINE_FUNCTION
      GeneralFad(const int sz, const T & x, const DerivInit zero_out = InitDerivArray) :
        Storage(sz, x, zero_out) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      SACADO_INLINE_FUNCTION
      GeneralFad(const int sz, const int i, const T & x) :
        Storage(sz, x, InitDerivArray) {
        this->fastAccessDx(i)=1.;
      }

      //! Constructor with supplied storage \c s
      SACADO_INLINE_FUNCTION
      GeneralFad(const Storage& s) : Storage(s) {}

      //! Copy constructor
      SACADO_INLINE_FUNCTION
      GeneralFad(const GeneralFad& x) :
        Storage(x) {}

      //! Copy constructor from any Expression object
      template <typename S>
      SACADO_INLINE_FUNCTION
      GeneralFad(const Expr<S>& x, SACADO_ENABLE_EXPR_CTOR_DECL)  :
        Storage(x.size(), T(0.), NoInitDerivArray) {
        x.cache();

        const int sz = x.size();

        this->val() = x.val();

        if (sz) {
          if (x.hasFastAccess())
            for(int i=0; i<sz; ++i)
              this->fastAccessDx(i) = x.fastAccessDx(i);
          else
            for(int i=0; i<sz; ++i)
              this->fastAccessDx(i) = x.dx(i);
        }
      }

      //! Destructor
      SACADO_INLINE_FUNCTION
      ~GeneralFad() {}

      //! Set %GeneralFad object as the \c ith independent variable
      /*!
       * Sets the derivative array of length \c n to the \c ith row of the
       * identity matrix and has the same affect as the
       * Implementation(const int sz, const int i, const T & x)
       * constructor.
       */
      SACADO_INLINE_FUNCTION
      void diff(const int ith, const int n) {
        if (this->size() != n)
          this->resize(n);

        this->zero();
        this->fastAccessDx(ith) = T(1.);
      }

      //! Set whether this Fad object should update values
      SACADO_INLINE_FUNCTION
      void setUpdateValue(bool update_val) {  }

      //! Return whether this Fad object has an updated value
      SACADO_INLINE_FUNCTION
      bool updateValue() const { return true; }

      //! Cache values
      SACADO_INLINE_FUNCTION
      void cache() const {}

      //! Returns whether two Fad objects have the same values
      template <typename S>
      SACADO_INLINE_FUNCTION
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
       * @name Derivative accessor methods
       */
      //@{

      /*!
       * \brief Returns number of derivative components that can be stored
       * without reallocation
       */
      SACADO_INLINE_FUNCTION
      int availableSize() const { return this->length(); }

      //! Returns true if derivative array is not empty
      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const { return this->size()!=0;}

      //! Returns true if derivative array is empty
      SACADO_INLINE_FUNCTION
      bool isPassive() const { return this->size()==0;}

      //! Set whether variable is constant
      SACADO_INLINE_FUNCTION
      void setIsConstant(bool is_const) {
        if (is_const && this->size()!=0)
          this->resize(0);
      }

      //@}

      /*!
       * @name Assignment operators
       */
      //@{

     //! Assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(GeneralFad&) operator=(const S& v) {
        this->val() = v;
        if (this->size()) this->resize(0);
        return *this;
      }

      //! Assignment with Expr right-hand-side
      SACADO_INLINE_FUNCTION
      GeneralFad&
      operator=(const GeneralFad& x) {
        // Copy val_ and dx_
        Storage::operator=(x);
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(GeneralFad&) operator=(const Expr<S>& x) {
        x.cache();

        const int xsz = x.size();

        if (xsz != this->size())
          this->resizeAndZero(xsz);

        const int sz = this->size();

        // For ViewStorage, the resize above may not in fact resize the
        // derivative array, so it is possible that sz != xsz at this point.
        // The only valid use case here is sz > xsz == 0, so we use sz in the
        // assignment below

        if (sz) {
          if (x.hasFastAccess())
            for(int i=0; i<sz; ++i)
              this->fastAccessDx(i) = x.fastAccessDx(i);
          else
            for(int i=0; i<sz; ++i)
              this->fastAccessDx(i) = x.dx(i);
        }

        this->val() = x.val();

        return *this;
      }

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Addition-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(GeneralFad&) operator += (const S& v) {
        this->val() += v;
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(GeneralFad&) operator -= (const S& v) {
        this->val() -= v;
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(GeneralFad&) operator *= (const S& v) {
        const int sz = this->size();
        this->val() *= v;
        for (int i=0; i<sz; ++i)
          this->fastAccessDx(i) *= v;
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(GeneralFad&) operator /= (const S& v) {
        const int sz = this->size();
        this->val() /= v;
        for (int i=0; i<sz; ++i)
          this->fastAccessDx(i) /= v;
        return *this;
      }

      //! Addition-assignment operator with GeneralFad right-hand-side
      SACADO_INLINE_FUNCTION
      GeneralFad& operator += (const GeneralFad& x) {
        const int xsz = x.size(), sz = this->size();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            for (int i=0; i<sz; ++i)
              this->fastAccessDx(i) += x.fastAccessDx(i);
          }
          else {
            this->resizeAndZero(xsz);
            for (int i=0; i<xsz; ++i)
              this->fastAccessDx(i) = x.fastAccessDx(i);
          }
        }

        this->val() += x.val();

        return *this;
      }

      //! Subtraction-assignment operator with GeneralFad right-hand-side
      SACADO_INLINE_FUNCTION
      GeneralFad& operator -= (const GeneralFad& x) {
        const int xsz = x.size(), sz = this->size();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            for(int i=0; i<sz; ++i)
              this->fastAccessDx(i) -= x.fastAccessDx(i);
          }
          else {
            this->resizeAndZero(xsz);
            for(int i=0; i<xsz; ++i)
              this->fastAccessDx(i) = -x.fastAccessDx(i);
          }
        }

        this->val() -= x.val();


        return *this;
      }

      //! Multiplication-assignment operator with GeneralFad right-hand-side
      SACADO_INLINE_FUNCTION
      GeneralFad& operator *= (const GeneralFad& x) {
        const int xsz = x.size(), sz = this->size();
        T xval = x.val();
        T v = this->val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            for(int i=0; i<sz; ++i)
              this->fastAccessDx(i) = v*x.fastAccessDx(i) + this->fastAccessDx(i)*xval;
          }
          else {
            this->resizeAndZero(xsz);
            for(int i=0; i<xsz; ++i)
              this->fastAccessDx(i) = v*x.fastAccessDx(i);
          }
        }
        else {
          if (sz) {
            for (int i=0; i<sz; ++i)
              this->fastAccessDx(i) *= xval;
          }
        }

        this->val() *= xval;

        return *this;
      }

      //! Division-assignment operator with GeneralFad right-hand-side
      SACADO_INLINE_FUNCTION
      GeneralFad& operator /= (const GeneralFad& x) {
        const int xsz = x.size(), sz = this->size();
        T xval = x.val();
        T v = this->val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            for(int i=0; i<sz; ++i)
              this->fastAccessDx(i) =
                ( this->fastAccessDx(i)*xval - v*x.fastAccessDx(i) )/ (xval*xval);
          }
          else {
            this->resizeAndZero(xsz);
            for(int i=0; i<xsz; ++i)
              this->fastAccessDx(i) = - v*x.fastAccessDx(i) / (xval*xval);
          }
        }
        else {
          if (sz) {
            for (int i=0; i<sz; ++i)
              this->fastAccessDx(i) /= xval;
          }
        }

        this->val() /= xval;

        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(GeneralFad&) operator += (const Expr<S>& x) {
        x.cache();

        const int xsz = x.size(), sz = this->size();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            if (x.hasFastAccess())
              for (int i=0; i<sz; ++i)
                this->fastAccessDx(i) += x.fastAccessDx(i);
            else
              for (int i=0; i<sz; ++i)
                this->fastAccessDx(i) += x.dx(i);
          }
          else {
            this->resizeAndZero(xsz);
            if (x.hasFastAccess())
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = x.fastAccessDx(i);
            else
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = x.dx(i);
          }
        }

        this->val() += x.val();

        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(GeneralFad&) operator -= (const Expr<S>& x) {
        x.cache();

        const int xsz = x.size(), sz = this->size();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            if (x.hasFastAccess())
              for(int i=0; i<sz; ++i)
                this->fastAccessDx(i) -= x.fastAccessDx(i);
            else
              for (int i=0; i<sz; ++i)
                this->fastAccessDx(i) -= x.dx(i);
          }
          else {
            this->resizeAndZero(xsz);
            if (x.hasFastAccess())
              for(int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = -x.fastAccessDx(i);
            else
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = -x.dx(i);
          }
        }

        this->val() -= x.val();


        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(GeneralFad&) operator *= (const Expr<S>& x) {
        x.cache();

        const int xsz = x.size(), sz = this->size();
        T xval = x.val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            if (x.hasFastAccess())
              for(int i=0; i<sz; ++i)
                this->fastAccessDx(i) = this->val() * x.fastAccessDx(i) + this->fastAccessDx(i) * xval;
            else
              for (int i=0; i<sz; ++i)
                this->fastAccessDx(i) = this->val() * x.dx(i) + this->fastAccessDx(i) * xval;
          }
          else {
            this->resizeAndZero(xsz);
            if (x.hasFastAccess())
              for(int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = this->val() * x.fastAccessDx(i);
            else
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = this->val() * x.dx(i);
          }
        }
        else {
          if (sz) {
            for (int i=0; i<sz; ++i)
              this->fastAccessDx(i) *= xval;
          }
        }

        this->val() *= xval;

        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(GeneralFad&) operator /= (const Expr<S>& x) {
        x.cache();

        const int xsz = x.size(), sz = this->size();
        T xval = x.val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            if (x.hasFastAccess())
              for(int i=0; i<sz; ++i)
                this->fastAccessDx(i) = ( this->fastAccessDx(i)*xval - this->val()*x.fastAccessDx(i) )/ (xval*xval);
            else
              for (int i=0; i<sz; ++i)
                this->fastAccessDx(i) = ( this->fastAccessDx(i)*xval - this->val()*x.dx(i) )/ (xval*xval);
          }
          else {
            this->resizeAndZero(xsz);
            if (x.hasFastAccess())
              for(int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = - this->val()*x.fastAccessDx(i) / (xval*xval);
            else
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = -this->val() * x.dx(i) / (xval*xval);
          }
        }
        else {
          if (sz) {
            for (int i=0; i<sz; ++i)
              this->fastAccessDx(i) /= xval;
          }
        }

        this->val() /= xval;

        return *this;
      }

      //@}

    }; // class GeneralFad

  } // namespace CacheFad

} // namespace Sacado

#endif // SACADO_CACHEFAD_GENERALFAD_HPP
