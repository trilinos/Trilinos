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

#ifndef SACADO_FAD_GENERALFAD_MP_VECTOR_HPP
#define SACADO_FAD_GENERALFAD_MP_VECTOR_HPP

#include "Sacado_Fad_GeneralFad.hpp"
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

    template <typename Ord, typename Val, int VecNum, typename Dev,
              typename Storage>
    struct ExprSpec< GeneralFad< Sacado::MP::Vector< Stokhos::StaticFixedStorage<Ord,Val,VecNum,Dev> >, Storage > > {
      typedef ExprSpecMPVector type;
    };

    //! Forward-mode AD class templated on the storage for the derivative array
    /*!
     * This class provides a general forward mode AD implementation for any
     * type of derivative array storage.  It does not incorporate expression
     * templates.
     */
    template <typename Ord, typename Val, int VecNum, typename Dev,
              typename Storage>
    class GeneralFad< Sacado::MP::Vector< Stokhos::StaticFixedStorage<Ord,Val,VecNum,Dev> >, Storage> : public Storage {

    public:

      typedef Sacado::MP::Vector< Stokhos::StaticFixedStorage<Ord,Val,VecNum,Dev> > T;

      //! Typename of values
      typedef typename RemoveConst<T>::type value_type;

      //! Typename of scalar's (which may be different from T)
      typedef typename ScalarType<value_type>::type scalar_type;

      typedef typename value_type::value_type val_type;

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor
      KOKKOS_INLINE_FUNCTION
      GeneralFad() : Storage(T(0.)) {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      GeneralFad(const S& x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        Storage(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      KOKKOS_INLINE_FUNCTION
      GeneralFad(const int sz, const T & x, const DerivInit zero_out = InitDerivArray) :
        Storage(sz, x, zero_out) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      KOKKOS_INLINE_FUNCTION
      GeneralFad(const int sz, const int i, const T & x) :
        Storage(sz, x) {
        this->fastAccessDx(i)=1.;
      }

      //! Constructor with supplied storage \c s
      KOKKOS_INLINE_FUNCTION
      GeneralFad(const Storage& s) : Storage(s) {}

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      GeneralFad(const GeneralFad& x) :
        Storage(x) {}

      //! Copy constructor from any Expression object
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      GeneralFad(const Expr<S>& x, SACADO_ENABLE_EXPR_CTOR_DECL)  :
        Storage(x.size(), T(0.))
      {
        const int sz = x.size();

        if (sz) {
          if (x.hasFastAccess())
            for(int i=0; i<sz; ++i)
              for (int j=0; j<VecNum; ++j)
                fastAccessDx(i,j) = x.fastAccessDx(i,j);
          else
            for(int i=0; i<sz; ++i)
              for (int j=0; j<VecNum; ++j)
                fastAccessDx(i,j) = x.dx(i,j);
        }

        for (int j=0; j<VecNum; ++j)
          val(j) = x.val(j);
      }

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~GeneralFad() {}

      //! Set %GeneralFad object as the \c ith independent variable
      /*!
       * Sets the derivative array of length \c n to the \c ith row of the
       * identity matrix and has the same affect as the
       * Implementation(const int sz, const int i, const T & x)
       * constructor.
       */
      KOKKOS_INLINE_FUNCTION
      void diff(const int ith, const int n) {
        if (this->size() != n)
          this->resize(n);

        this->zero();
        this->fastAccessDx(ith) = T(1.);
      }

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
      const T& val() const { return Storage::val();}

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      T& val() { return Storage::val();}

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      const val_type& val(int j) const { return Storage::val().fastAccessCoeff(j);}

      //! Returns value
      KOKKOS_INLINE_FUNCTION
      val_type& val(int j) { return Storage::val().fastAccessCoeff(j);}

      //@}

      /*!
       * @name Derivative accessor methods
       */
      //@{

      /*!
       * \brief Returns number of derivative components that can be stored
       * without reallocation
       */
      KOKKOS_INLINE_FUNCTION
      int availableSize() const { return this->length(); }

      //! Returns true if derivative array is not empty
      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess() const { return this->size()!=0; }

      //! Returns true if derivative array is empty
      KOKKOS_INLINE_FUNCTION
      bool isPassive() const { return this->size()==0; }

      //! Set whether variable is constant
      KOKKOS_INLINE_FUNCTION
      void setIsConstant(bool is_const) {
        if (is_const && this->size()!=0)
          this->resize(0);
      }

     //! Returns derivative array
      KOKKOS_INLINE_FUNCTION
      const T* dx() const { return this->dx_; }

      //! Returns derivative component \c i with bounds checking
      KOKKOS_INLINE_FUNCTION
      T dx(int i) const { return this->size() ? this->dx_[i] : T(0.); }

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      T& fastAccessDx(int i) { return this->dx_[i];}

      //! Returns derivative component \c i without bounds checking
      KOKKOS_INLINE_FUNCTION
      const T& fastAccessDx(int i) const { return this->dx_[i];}

      //! Returns derivative component \c i with bounds checking
      KOKKOS_INLINE_FUNCTION
      val_type dx(int i, int j) const { return this->size() ? this->dx_[i].fastAccessCoeff(j) : val_type(0.0); }

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

      //@}

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(GeneralFad&) operator=(const S& v) {
        this->val() = v;
        if (this->size()) this->resize(0);
        return *this;
      }

      //! Assignment with GeneralFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      GeneralFad&
      operator=(const GeneralFad& x) {
        // Copy value & dx_
        Storage::operator=(x);
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(GeneralFad&) operator=(const Expr<S>& x) {
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
              for (int j=0; j<VecNum; ++j)
                fastAccessDx(i,j) = x.fastAccessDx(i,j);
          else
            for(int i=0; i<sz; ++i)
              for (int j=0; j<VecNum; ++j)
                fastAccessDx(i,j) = x.dx(i,j);
        }

        for (int j=0; j<VecNum; ++j)
          val(j) = x.val(j);

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
      SACADO_ENABLE_VALUE_FUNC(GeneralFad&) operator += (const S& v) {
        this->val() += v;
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(GeneralFad&) operator -= (const S& v) {
        this->val() -= v;
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(GeneralFad&) operator *= (const S& v) {
        const int sz = this->size();
        this->val() *= v;
        for (int i=0; i<sz; ++i)
          this->fastAccessDx(i) *= v;
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(GeneralFad&) operator /= (const S& v) {
        const int sz = this->size();
        this->val() /= v;
        for (int i=0; i<sz; ++i)
          this->fastAccessDx(i) /= v;
        return *this;
      }

      //! Addition-assignment operator with GeneralFad right-hand-side
      KOKKOS_INLINE_FUNCTION
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
      KOKKOS_INLINE_FUNCTION
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
      KOKKOS_INLINE_FUNCTION
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
      KOKKOS_INLINE_FUNCTION
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
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(GeneralFad&) operator += (const Expr<S>& x) {
        const int xsz = x.size();
        int sz = this->size();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz > sz) {
          this->resizeAndZero(xsz);
          sz = this->size();
        }

        if (sz) {
          if (x.hasFastAccess())
            for(int i=0; i<sz; ++i)
              for (int j=0; j<VecNum; ++j)
                fastAccessDx(i,j) += x.fastAccessDx(i,j);
          else
            for(int i=0; i<sz; ++i)
              for (int j=0; j<VecNum; ++j)
                fastAccessDx(i,j) += x.dx(i,j);
        }

        for (int j=0; j<VecNum; ++j)
          val(j) += x.val(j);

        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(GeneralFad&) operator -= (const Expr<S>& x) {
        const int xsz = x.size(), sz = this->size();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            if (x.hasFastAccess())
              for(int i=0; i<sz; ++i)
                for (int j=0; j<VecNum; ++j)
                  fastAccessDx(i,j) -= x.fastAccessDx(i,j);
            else
              for (int i=0; i<sz; ++i)
                for (int j=0; j<VecNum; ++j)
                  fastAccessDx(i,j) -= x.dx(i,j);
          }
          else {
            this->resizeAndZero(xsz);
            if (x.hasFastAccess())
              for(int i=0; i<xsz; ++i)
                for (int j=0; j<VecNum; ++j)
                  fastAccessDx(i,j) = -x.fastAccessDx(i,j);
            else
              for (int i=0; i<xsz; ++i)
                for (int j=0; j<VecNum; ++j)
                  fastAccessDx(i,j) = -x.dx(i,j);
          }
        }

        for (int j=0; j<VecNum; ++j)
          val(j) -= x.val(j);

        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(GeneralFad&) operator *= (const Expr<S>& x) {
        const int xsz = x.size(), sz = this->size();
        T xval = x.val();
        T v = this->val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            if (x.hasFastAccess())
              for(int i=0; i<sz; ++i)
                for (int j=0; j<VecNum; ++j)
                  fastAccessDx(i,j) = v.fastAccessCoeff(j)*x.fastAccessDx(i,j) + fastAccessDx(i,j)*xval.fastAccessCoeff(j);
            else
              for (int i=0; i<sz; ++i)
                for (int j=0; j<VecNum; ++j)
                  fastAccessDx(i,j) = v.fastAccessCoeff(j)*x.dx(i,j) + fastAccessDx(i,j)*xval.fastAccessCoeff(j);
          }
          else {
            this->resizeAndZero(xsz);
            if (x.hasFastAccess())
              for(int i=0; i<xsz; ++i)
                for (int j=0; j<VecNum; ++j)
                  fastAccessDx(i,j) = v.fastAccessCoeff(j)*x.fastAccessDx(i,j);
            else
              for (int i=0; i<xsz; ++i)
                for (int j=0; j<VecNum; ++j)
                  fastAccessDx(i,j) = v.fastAccessCoeff(j)*x.dx(i,j);
          }
        }
        else {
          if (sz) {
            for (int i=0; i<sz; ++i)
              for (int j=0; j<VecNum; ++j)
                fastAccessDx(i,j) *= xval.fastAccessCoeff(j);
          }
        }

        for (int j=0; j<VecNum; ++j)
          val(j) *= xval.fastAccessCoeff(j);

        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(GeneralFad&) operator /= (const Expr<S>& x) {
        const int xsz = x.size(), sz = this->size();
        T xval = x.val();
        T v = this->val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          T xval2 = xval*xval;
          if (sz) {
            if (x.hasFastAccess())
              for(int i=0; i<sz; ++i)
                for (int j=0; j<VecNum; ++j)
                  fastAccessDx(i,j) = ( fastAccessDx(i,j)*xval.fastAccessCoeff(j) - v.fastAccessCoeff(j)*x.fastAccessDx(i,j) )/ (xval2.fastAccessCoeff(j));
            else
              for (int i=0; i<sz; ++i)
                for (int j=0; j<VecNum; ++j)
                  fastAccessDx(i,j) = ( fastAccessDx(i,j)*xval.fastAccessCoeff(j) - v.fastAccessCoeff(j)*x.dx(i,j) )/ (xval2.fastAccessCoeff(j));
          }
          else {
            this->resizeAndZero(xsz);
            if (x.hasFastAccess())
              for(int i=0; i<xsz; ++i)
                for (int j=0; j<VecNum; ++j)
                  fastAccessDx(i,j) = - v.fastAccessCoeff(j)*x.fastAccessDx(i,j) / (xval2.fastAccessCoeff(j));
            else
              for (int i=0; i<xsz; ++i)
                for (int j=0; j<VecNum; ++j)
                  fastAccessDx(i,j) = -v.fastAccessCoeff(j)*x.dx(i,j) / (xval2.fastAccessCoeff(j));
          }
        }
        else {
          if (sz) {
            for (int i=0; i<sz; ++i)
              for (int j=0; j<VecNum; ++j)
                this->fastAccessDx(i,j) /= xval.fastAccessCoeff(j);
          }
        }

        for (int j=0; j<VecNum; ++j)
          val(j) /= xval.fastAccessCoeff(j);

        return *this;
      }

      //@}

    private:

      template <typename S>
      KOKKOS_INLINE_FUNCTION
      void
      fastAssign(const Expr<S>& x) {
        const int sz = this->size();
        for(int i=0; i<sz; ++i)
          for (int j=0; j<VecNum; ++j)
            fastAccessDx(i,j) = x.fastAccessDx(i,j);

        for (int j=0; j<VecNum; ++j)
          val(j) = x.val(j);
      }

      template <typename S>
      KOKKOS_INLINE_FUNCTION
      void
      slowAssign(const Expr<S>& x) {
        const int sz = this->size();
        for(int i=0; i<sz; ++i)
          for (int j=0; j<VecNum; ++j)
            fastAccessDx(i,j) = x.dx(i,j);

        for (int j=0; j<VecNum; ++j)
          val(j) = x.val(j);
      }

    }; // class GeneralFad

  } // namespace Fad

} // namespace Sacado

#endif // SACADO_FAD_GENERALFAD_HPP
