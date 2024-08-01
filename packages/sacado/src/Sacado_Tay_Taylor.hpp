// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_TAY_TAYLOR_HPP
#define SACADO_TAY_TAYLOR_HPP

#include "Sacado_ConfigDefs.h"
#include "Sacado_Base.hpp"
#include "Sacado_Handle.hpp"
#include <cmath>
#include <algorithm>    // for std::min and std::max
#include <ostream>      // for std::ostream
#include "Sacado_dummy_arg.hpp"

namespace Sacado {

  //! Namespace for Taylor polynomial AD classes
  namespace Tay {

    //! Taylor polynomial class
    /*!
     * Uses a handle and a "copy-on-write" strategy for efficient copying, but
     * no expression templating.
     */
    template <typename T>
    class Taylor : public Base< Taylor<T> > {
    public:

      //! Turn Taylor into a meta-function class usable with mpl::apply
      template <typename U>
      struct apply {
        typedef Taylor<U> type;
      };

      //! Typename of values
      typedef T value_type;

      //! Typename of scalar's (which may be different from value_type)
      typedef typename ScalarType<value_type>::type scalar_type;

      //! Default constructor
      Taylor();

      //! Constructor with supplied value \c x
      /*!
       * Sets the first coefficient to x
       */
      Taylor(const T& x);

      //! Constructor with supplied value \c x
      /*!
       * Sets the first coefficient to x.
       * Creates a dummy overload when ValueT and ScalarT are the same type.
       */
      Taylor(const typename dummy<value_type,scalar_type>::type& x);

      //! Constructor with degree d and value \c x
      /*!
       * Initializes first coeffienct to \c x and of a polynomial of degree d
       */
      Taylor(int d, const T & x);

      //! Constructor with degree d
      /*!
       * Initializes all components to zero
       */
      explicit Taylor(int d);

      //! Copy constructor
      Taylor(const Taylor& x);

      //! Destructor
      ~Taylor();

      //! Resize polynomial to degree d
      /*!
       * Coefficients are preserved if \c keep_coeffs is \c true, otherwise
       * all coefficients are reset to zero.
       */
      void resize(int d, bool keep_coeffs);

      //! Reserve space for a degree d polynomial
      /*!
       * Coefficients are preserved.
       */
      void reserve(int d);

      //! Prepare polynomial for writing
      /*!
       * This method prepares the polynomial for writing through coeff() and
       * fastAccessCoeff() member functions.  It ensures the handle for the
       * %Taylor coefficients is not shared among any other %Taylor polynomial
       * objects.  If the handle is not shared it does nothing, so there
       * is no cost in calling this method in this case.  If the handle is
       * shared and this method is not called, any changes to the coefficients
       * by coeff() or fastAccessCoeff() may change other polynomial objects.
       */
      void copyForWrite() { th.makeOwnCopy(); }

      //! Returns whether two Taylor objects have the same values
      bool isEqualTo(const Taylor& x) const {
        typedef IsEqual<value_type> IE;
        if (x.degree() != this->degree()) return false;
        bool eq = true;
        for (int i=0; i<=this->degree(); i++)
          eq = eq && IE::eval(x.coeff(i), this->coeff(i));
        return eq;
      }

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      Taylor<T>& operator=(const T& val);

      //! Assignment operator with constant right-hand-side
      /*!
       * Creates a dummy overload when value_type and scalar_type are
       * the same type.
       */
      Taylor<T>& operator=(const typename dummy<value_type,scalar_type>::type& val);

      //! Assignment with Taylor right-hand-side
      Taylor<T>& operator=(const Taylor<T>& x);

      //@}

      /*!
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      const T& val() const { return th->coeff_[0];}

      //! Returns value
      T& val() { return th->coeff_[0];}

      //@}

      /*!
       * @name Taylor coefficient accessor methods
       */
      //@{

      //! Returns degree of polynomial
      int degree() const { return th->deg_;}

      //! Returns true if polynomial has degree >= d
      bool hasFastAccess(int d) const { return th->deg_>=d;}

      //! Returns Taylor coefficient array
      const T* coeff() const { return th->coeff_;}

      //! Returns Taylor coefficient array
      T* coeff() { return th->coeff_;}

      //! Returns degree \c i term with bounds checking
      T coeff(int i) const {
        T tmp= i<=th->deg_ ? th->coeff_[i]:T(0.); return tmp;}

      //! Returns degree \c i term without bounds checking
      T& fastAccessCoeff(int i) { return th->coeff_[i];}

      //! Returns degree \c i term without bounds checking
      const T& fastAccessCoeff(int i) const { return th->coeff_[i];}

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Unary-plus operator
      Taylor<T> operator + () const;

      //! Unary-minus operator
      Taylor<T> operator - () const;

      //! Addition-assignment operator with constant right-hand-side
      Taylor<T>& operator += (const T& x);

      //! Subtraction-assignment operator with constant right-hand-side
      Taylor<T>& operator -= (const T& x);

      //! Multiplication-assignment operator with constant right-hand-side
      Taylor<T>& operator *= (const T& x);

      //! Division-assignment operator with constant right-hand-side
      Taylor<T>& operator /= (const T& x);

      //! Addition-assignment operator with Taylor right-hand-side
      Taylor<T>& operator += (const Taylor<T>& x);

      //! Subtraction-assignment operator with Taylor right-hand-side
      Taylor<T>& operator -= (const Taylor<T>& x);

      //! Multiplication-assignment operator with Taylor right-hand-side
      Taylor<T>& operator *= (const Taylor<T>& x);

      //! Division-assignment operator with Taylor right-hand-side
      Taylor<T>& operator /= (const Taylor<T>& x);

      //@}

    protected:

      //! Return length of array
      int length() const { return th->len_; }

      //! Resize coefficient array to new size
      void resizeCoeffs(int len);

    protected:

      struct TaylorData {

        //! Taylor polynomial coefficients
        T* coeff_;

        //! Degree of polynomial
        int deg_;

        //! Length of allocated polynomial array
        int len_;

        //! Default constructor
        TaylorData();

        //! Constructor with supplied value \c x
        TaylorData(const T& x);

        //! Constructor with degree d and value \c x
        TaylorData(int d, const T & x);

        //! Constructor with degree d
        TaylorData(int d);

        //! Constructor with degree d and length l
        TaylorData(int d, int l);

        //! Copy constructor
        TaylorData(const TaylorData& x);

        //! Destructor
        ~TaylorData();

        //! Assignment operator
        TaylorData& operator=(const TaylorData& x);

      };

      Sacado::Handle<TaylorData> th;

    }; // class Taylor

    //! Compute Taylor series of n-th derivative of x
    template <typename T>
    Taylor<T> diff(const Taylor<T>& x, int n = 1) {
      const int d = x.degree();
      if (n <= 0)
        return x;
      else if (n > d) {
        Taylor<T> y(0);
        return y;
      }
      Taylor<T> y(d-n);
      int c = 1;
      for (int i=1; i<=n; ++i)
        c *= i;
      for (int i=n; i<=d; ++i) {
        y.fastAccessCoeff(i-n) = x.fastAccessCoeff(i) * T(c);
        c = (c / (i-n+1)) * (i+1);
      }
      return y;
    }

    // Operations
    template <typename T> Taylor<T> operator+(const Base< Taylor<T> >& a,
                                              const Base< Taylor<T> >& b);
    template <typename T> Taylor<T> operator+(const typename Taylor<T>::value_type& a,
                                              const Base< Taylor<T> >& b);
    template <typename T> Taylor<T> operator+(const Base< Taylor<T> >& a,
                                              const typename Taylor<T>::value_type& b);
    template <typename T> Taylor<T> operator-(const Base< Taylor<T> >& a,
                                              const Base< Taylor<T> >& b);
    template <typename T> Taylor<T> operator-(const typename Taylor<T>::value_type& a,
                                              const Base< Taylor<T> >& b);
    template <typename T> Taylor<T> operator-(const Base< Taylor<T> >& a,
                                              const typename Taylor<T>::value_type& b);
    template <typename T> Taylor<T> operator*(const Base< Taylor<T> >& a,
                                              const Base< Taylor<T> >& b);
    template <typename T> Taylor<T> operator*(const typename Taylor<T>::value_type& a,
                                              const Base< Taylor<T> >& b);
    template <typename T> Taylor<T> operator*(const Base< Taylor<T> >& a,
                                              const typename Taylor<T>::value_type& b);
    template <typename T> Taylor<T> operator/(const Base< Taylor<T> >& a,
                                              const Base< Taylor<T> >& b);
    template <typename T> Taylor<T> operator/(const typename Taylor<T>::value_type& a,
                                              const Base< Taylor<T> >& b);
    template <typename T> Taylor<T> operator/(const Base< Taylor<T> >& a,
                                              const typename Taylor<T>::value_type& b);
    template <typename T> Taylor<T> exp(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> log(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> log10(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> sqrt(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> cbrt(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> pow(const Base< Taylor<T> >& a,
                                        const Base< Taylor<T> >& b);
    template <typename T> Taylor<T> pow(const typename Taylor<T>::value_type& a,
                                        const Base< Taylor<T> >& b);
    template <typename T> Taylor<T> pow(const Base< Taylor<T> >& a,
                                        const typename Taylor<T>::value_type& b);
    template <typename T> void sincos(const Base< Taylor<T> >& a,
                                      Taylor<T>& s, Taylor<T>& c);
    template <typename T> Taylor<T> cos(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> sin(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> tan(const Base< Taylor<T> >& a);
    template <typename T> void sinhcosh(const Base< Taylor<T> >& a,
                                        Taylor<T>& s, Taylor<T>& c);
    template <typename T> Taylor<T> cosh(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> sinh(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> tanh(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> quad(const typename Taylor<T>::value_type& c0,
                                         const Base< Taylor<T> >& a,
                                         const Base< Taylor<T> >& b);
    template <typename T> Taylor<T> acos(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> asin(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> atan(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> atan2(const Base< Taylor<T> >& a,
                                          const Base< Taylor<T> >& b);
    template <typename T> Taylor<T> atan2(const typename Taylor<T>::value_type& a,
                                          const Base< Taylor<T> >& b);
    template <typename T> Taylor<T> atan2(const Base< Taylor<T> >& a,
                                          const typename Taylor<T>::value_type& b);
    template <typename T> Taylor<T> acosh(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> asinh(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> atanh(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> abs(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> fabs(const Base< Taylor<T> >& a);
    template <typename T> Taylor<T> max(const Base< Taylor<T> >& a,
                                        const Base< Taylor<T> >& b);
    template <typename T> Taylor<T> max(const typename Taylor<T>::value_type& a,
                                        const Base< Taylor<T> >& b);
    template <typename T> Taylor<T> max(const Base< Taylor<T> >& a,
                                        const typename Taylor<T>::value_type& b);
    template <typename T> Taylor<T> min(const Base< Taylor<T> >& a,
                                        const Base< Taylor<T> >& b);
    template <typename T> Taylor<T> min(const typename Taylor<T>::value_type& a,
                                        const Base< Taylor<T> >& b);
    template <typename T> Taylor<T> min(const Base< Taylor<T> >& a,
                                        const typename Taylor<T>::value_type& b);
    template <typename T> bool operator==(const Base< Taylor<T> >& a,
                                          const Base< Taylor<T> >& b);
    template <typename T> bool operator==(const typename Taylor<T>::value_type& a,
                                          const Base< Taylor<T> >& b);
    template <typename T> bool operator==(const Base< Taylor<T> >& a,
                                          const typename Taylor<T>::value_type& b);
    template <typename T> bool operator!=(const Base< Taylor<T> >& a,
                                          const Base< Taylor<T> >& b);
    template <typename T> bool operator!=(const typename Taylor<T>::value_type& a,
                                          const Base< Taylor<T> >& b);
    template <typename T> bool operator!=(const Base< Taylor<T> >& a,
                                          const typename Taylor<T>::value_type& b);
    template <typename T> bool operator<=(const Base< Taylor<T> >& a,
                                          const Base< Taylor<T> >& b);
    template <typename T> bool operator<=(const typename Taylor<T>::value_type& a,
                                          const Base< Taylor<T> >& b);
    template <typename T> bool operator<=(const Base< Taylor<T> >& a,
                                          const typename Taylor<T>::value_type& b);
    template <typename T> bool operator>=(const Base< Taylor<T> >& a,
                                          const Base< Taylor<T> >& b);
    template <typename T> bool operator>=(const typename Taylor<T>::value_type& a,
                                          const Base< Taylor<T> >& b);
    template <typename T> bool operator>=(const Base< Taylor<T> >& a,
                                          const typename Taylor<T>::value_type& b);
    template <typename T> bool operator<(const Base< Taylor<T> >& a,
                                         const Base< Taylor<T> >& b);
    template <typename T> bool operator<(const typename Taylor<T>::value_type& a,
                                         const Base< Taylor<T> >& b);
    template <typename T> bool operator<(const Base< Taylor<T> >& a,
                                         const typename Taylor<T>::value_type& b);
    template <typename T> bool operator>(const Base< Taylor<T> >& a,
                                         const Base< Taylor<T> >& b);
    template <typename T> bool operator>(const typename Taylor<T>::value_type& a,
                                         const Base< Taylor<T> >& b);
    template <typename T> bool operator>(const Base< Taylor<T> >& a,
                                         const typename Taylor<T>::value_type& b);
    template <typename T> std::ostream& operator << (std::ostream& os,
                                                     const Base< Taylor<T> >& a);

  } // namespace Tay

} // namespace Sacado

#include "Sacado_Tay_TaylorTraits.hpp"
#include "Sacado_Tay_TaylorImp.hpp"

#endif // SACADO_TAY_TAYLOR_HPP
