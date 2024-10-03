// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_TAY_CACHETAYLOR_HPP
#define SACADO_TAY_CACHETAYLOR_HPP

#include <valarray>

#include "Sacado_Tay_CacheTaylorExpr.hpp"

// forward decalarations
namespace Sacado {
  namespace Tay {
    template <class ExprT> class UnaryPlusOp;
    template <class ExprT> class UnaryMinusOp;
  }
}

namespace Sacado {

  namespace Tay {

    // Forward declaration
    template <typename T> class CacheTaylor;

    //! Taylor polynomial class using caching expression templates
    /*!
     * This class provides the implementation of the Taylor object required
     * for expression templating.  Class CacheTaylor provides the complete
     * user inteface.
     */
    template <typename T>
    class CacheTaylorImplementation {

    public:

      //! Typename of values
      typedef T value_type;

      //! Typename of scalar's (which may be different from ValueT)
      typedef typename ScalarType<T>::type scalar_type;

      //! Default constructor
      CacheTaylorImplementation() : coeff_(T(0.),1) {}

      //! Constructor with supplied value \c x
      /*!
       * Sets the first coefficient to x
       */
      CacheTaylorImplementation(const T& x) : coeff_(x,1) {}

      //! Constructor with degree d and value \c x
      /*!
       * Initializes first coeffienct to \c x and of a polynomial of degree d
       */
      CacheTaylorImplementation(int d, const T & x) :
        coeff_(T(0.),d+1) {
        coeff_[0] = x;
      }

      //! Copy constructor
      CacheTaylorImplementation(const CacheTaylorImplementation& x) :
        coeff_(x.coeff_) {}

      //! Destructor
      ~CacheTaylorImplementation() {}

       //! Resize polynomial to degree d
      /*!
       * Coefficients are preserved if \c keep_coeffs is \c true, otherwise
       * all coefficients are reset to zero.
       */
      void resize(int d, bool keep_coeffs) {
        if (keep_coeffs)
          resizeCoeffs(d);
        else
          coeff_.resize(d+1, T(0.));
      }

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      const T& val() const { return coeff_[0];}

      //! Returns value
      T& val() { return coeff_[0];}

      //@}

      /*!
       * @name Taylor coefficient accessor methods
       */
      //@{

      //! Returns degree of polynomial
      int degree() const { return coeff_size()-1;}

      //! Returns true if polynomial has degree >= d
      bool hasFastAccess(int d) const { return coeff_size()>=d+1;}

      //! Returns Taylor coefficient array
      const std::valarray<T>& coeff() const { return coeff_;}

      //! Returns degree \c i term with bounds checking
      const T coeff(int i) const {
        T tmp= i<coeff_size() ? coeff_[i]:T(0.); return tmp;}

      //! Returns degree \c i term with bounds checking
      T coeff(int i) {
        T tmp= i<coeff_size() ? coeff_[i]:T(0.); return tmp;}

      //! Returns degree \c i term without bounds checking
      T& fastAccessCoeff(int i) { return coeff_[i];}

      //! Returns degree \c i term without bounds checking
      const T& fastAccessCoeff(int i) const { return coeff_[i];}

      //! Allocate coefficient cache
      void allocateCache(int d) const {}

      //! Returns whether two Taylor objects have the same values
      template <typename S>
      bool isEqualTo(const Expr<S>& x) const {
        typedef IsEqual<value_type> IE;
        if (x.degree() != this->degree()) return false;
        bool eq = true;
        for (int i=0; i<=this->degree(); i++)
          eq = eq && IE::eval(x.coeff(i), this->coeff(i));
        return eq;
      }

      //@}

    protected:

      //! Resize coefficient array to new size
      void resizeCoeffs(int dnew) {
        std::valarray<T> tmp = coeff_;
        int sz = coeff_size();
        coeff_.resize(dnew+1,T(0.));
        if (sz > dnew+1) {
          std::slice s(0,dnew+1,1);
          coeff_ = tmp[s];
        }
        else {
          std::slice s(0,sz,1);
          coeff_[s] = tmp;
        }
      }

      int coeff_size() const { return coeff_.size(); }

    protected:

      //! Taylor polynomial coefficients
      std::valarray<T> coeff_;

    }; // class CacheTaylorImplementation

    //! CacheTaylor expression template specialization
    /*!
     * This template class represents a simple CacheTaylor expression.
     */
    template <typename T>
    class Expr< CacheTaylorImplementation<T> > :
      public CacheTaylorImplementation<T> {

    public:

      //! Typename of base-expressions
      typedef CacheTaylor<T> base_expr_type;

      //! Default constructor
      Expr() : CacheTaylorImplementation<T>() {}

      //! Constructor with supplied value \c x
      /*!
       * Sets the first coefficient to x
       */
      Expr(const T & x) : CacheTaylorImplementation<T>(x) {}

      //! Constructor with degree d and value \c x
      /*!
       * Initializes first coeffienct to \c x and of a polynomial of degree d
       */
      Expr(int d, const T & x) : CacheTaylorImplementation<T>(d,x) {}

      //! Copy constructor
      Expr(const Expr& x) : CacheTaylorImplementation<T>(x) {}

    }; // class Expr< CacheTaylorImplementation<T> >

    //! Forward-mode AD class using dynamic memory allocation
    /*!
     * This class provides the user interface of the Taylor object.  Class
     * CacheTaylorImplementation provides the implementation.
     */
    template <typename T>
    class CacheTaylor : public Expr< CacheTaylorImplementation<T> > {

    public:

      //! Typename of values
      typedef T value_type;

      //! Typename of scalar's (which may be different from ValueT)
      typedef typename ScalarType<T>::type scalar_type;

      //! Turn CacheTaylor into a meta-function class usable with mpl::apply
      template <typename U>
      struct apply {
        typedef CacheTaylor<U> type;
      };

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      CacheTaylor() : Expr< CacheTaylorImplementation<T> >() {}

      //! Constructor with supplied value \c x
      /*!
       * Sets the first coefficient to x
       */
      CacheTaylor(const T & x) : Expr< CacheTaylorImplementation<T> >(x) {}

      //! Constructor with supplied value \c x
      /*!
       * Sets the first coefficient to x.
       * Creates a dummy overload when ValueT and ScalarT are the same type.
       */
      CacheTaylor(const typename dummy<value_type,scalar_type>::type& x) :
        Expr< CacheTaylorImplementation<T> >(value_type(x)) {}

      //! Constructor with degree d and value \c x
      /*!
       * Initializes first coeffienct to \c x and of a polynomial of degree d
       */
      CacheTaylor(int d, const T & x) :
        Expr< CacheTaylorImplementation<T> >(d,x) {}

      //! Copy constructor
      CacheTaylor(const CacheTaylor& x) : Expr< CacheTaylorImplementation<T> >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> CacheTaylor(const Expr<S>& x);

      //@}

      //! Destructor
      ~CacheTaylor() {}

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      CacheTaylor<T>& operator=(const T& v);

      //! Assignment operator with constant right-hand-side
      /*!
       * Creates a dummy overload when value_type and scalar_type are
       * the same type.
       */
      CacheTaylor<T>&
      operator=(const typename dummy<value_type,scalar_type>::type& val) {
        return operator=(value_type(val));
      }

      //! Assignment with CacheTaylor right-hand-side
      CacheTaylor<T>& operator=(const CacheTaylor<T>& x);

      //! Assignment operator with any expression right-hand-side
      template <typename S> CacheTaylor<T>& operator=(const Expr<S>& x);

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Unary-plus operator
      inline Expr< UnaryExpr< CacheTaylor<T>, UnaryPlusOp > >
      operator + () const {
        typedef UnaryExpr< CacheTaylor<T>, UnaryPlusOp > expr_t;
        return Expr<expr_t>(expr_t(*this));
      }

      //! Unary-minus operator
      inline Expr< UnaryExpr< CacheTaylor<T>, UnaryMinusOp > >
      operator - () const {
        typedef UnaryExpr< CacheTaylor<T>, UnaryMinusOp > expr_t;
        return Expr<expr_t>(expr_t(*this));
      }

      //! Addition-assignment operator with constant right-hand-side
      CacheTaylor<T>& operator += (const T& x);

      //! Subtraction-assignment operator with constant right-hand-side
      CacheTaylor<T>& operator -= (const T& x);

      //! Multiplication-assignment operator with constant right-hand-side
      CacheTaylor<T>& operator *= (const T& x);

      //! Division-assignment operator with constant right-hand-side
      CacheTaylor<T>& operator /= (const T& x);

      //! Addition-assignment operator with Taylor right-hand-side
      template <typename S> CacheTaylor<T>& operator += (const S& x);

      //! Subtraction-assignment operator with Taylor right-hand-side
      template <typename S> CacheTaylor<T>& operator -= (const S& x);

      //! Multiplication-assignment operator with Taylor right-hand-side
      template <typename S> CacheTaylor<T>& operator *= (const S& x);

      //! Division-assignment operator with Taylor right-hand-side
      template <typename S> CacheTaylor<T>& operator /= (const S& x);

      //@}

    }; // class CacheTaylor<T>

  } // namespace Tay

} // namespace Sacado

#include "Sacado_Tay_CacheTaylorTraits.hpp"
#include "Sacado_Tay_CacheTaylorImp.hpp"
#include "Sacado_Tay_CacheTaylorOps.hpp"

#endif // SACADO_TAYLOR_CACHETAYLOR_HPP
