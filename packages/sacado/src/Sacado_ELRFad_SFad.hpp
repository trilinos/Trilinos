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

#ifndef SACADO_ELRFAD_SFAD_HPP
#define SACADO_ELRFAD_SFAD_HPP

#include "Sacado_ELRFad_SFadTraits.hpp"
#include "Sacado_ELRFad_Expression.hpp"
#include "Sacado_StaticArrayTraits.hpp"
#include "Sacado_mpl_range_c.hpp"
#include "Sacado_mpl_for_each.hpp"

namespace Sacado {

  namespace ELRFad {

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
      typedef SFad<value_type,Num> base_expr_type;

      //! Number of arguments
      static const int num_args = 1;

      //! Is expression linear
      static const bool is_linear = true;

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor
      SACADO_INLINE_FUNCTION
      Expr() : val_( T(0.)) { ss_array<T>::zero(dx_, Num); }

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      template <typename S>
      SACADO_INLINE_FUNCTION
      Expr(const S & x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        val_(x) {
        ss_array<T>::zero(dx_, Num);
      }

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      SACADO_INLINE_FUNCTION
      Expr(const int sz, const T & x, const DerivInit zero_out = InitDerivArray) : val_(x) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (sz != Num)
          throw "SELRFad::SFad() Error:  Supplied derivative dimension does not match compile time length.";
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
      SACADO_INLINE_FUNCTION
      Expr(const int sz, const int i, const T & x) :
        val_(x) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (sz != Num)
          throw "SELRFad::SFad() Error:  Supplied derivative dimension does not match compile time length.";
        if (i >= Num)
          throw "SELRFad::SFad() Error:  Invalid derivative index.";
#endif

        ss_array<T>::zero(dx_, Num);
        dx_[i]=1.;
      }

      //! Copy constructor
      SACADO_INLINE_FUNCTION
      Expr(const Expr& x) : val_(x.val_) {
        //ss_array<T>::copy(x.dx_, dx_, Num);
        for (int i=0; i<Num; i++)
          dx_[i] = x.dx_[i];
      }

      //! Copy constructor from any Expression object
      template <typename S>
      SACADO_INLINE_FUNCTION
      Expr(const Expr<S>& x, SACADO_ENABLE_EXPR_CTOR_DECL) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (x.size() != Num)
          throw "SELRFad::SFad() Error:  Attempt to assign with incompatible sizes";
#endif

        // Compute partials
        LocalAccumOp< Expr<S> > op(x);

        // Compute each tangent direction
        for(int i=0; i<Num; ++i) {
          op.t = T(0.);
          op.i = i;

          // Automatically unrolled loop that computes
          // for (int j=0; j<N; j++)
          //   op.t += op.partials[j] * x.getTangent<j>(i);
          Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

          dx_[i] = op.t;

        }

        this->val() = x.val();
      }

      //! Destructor
      SACADO_INLINE_FUNCTION
      ~Expr() {}

      //! Set %Fad object as the \c ith independent variable
      /*!
       * Sets the derivative array of length \c n to the \c ith row of the
       * identity matrix and has the same affect as the
       * Implementation(const int sz, const int i, const T & x)
       * constructor.
       */
      SACADO_INLINE_FUNCTION
      void diff(const int ith, const int n) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (n != Num)
          throw "SELRFad::diff() Error:  Supplied derivative dimension does not match compile time length.";
#endif

        ss_array<T>::zero(dx_, Num);
        dx_[ith] = T(1.);
      }

      //! Resize derivative array to length \c sz
      /*!
       * Since the derivative array length is not dynamic, this method
       * throws an error if compiled with SACADO_DEBUG defined.
       */
      SACADO_INLINE_FUNCTION
      void resize(int sz) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (sz != Num)
          throw "SELRFad::resize() Error:  Cannot resize fixed derivative array dimension";
#endif
      }

      //! Expand derivative array to size sz
      /*!
       * Since the derivative array length is not dynamic, this method
       * throws an error if compiled with SACADO_DEBUG defined.
       */
      SACADO_INLINE_FUNCTION
      void expand(int sz) { resize(sz); }

      //! Zero out the derivative array
      SACADO_INLINE_FUNCTION
      void zero() { ss_array<T>::zero(dx_, Num); }

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
       * @name Value accessor methods
       */
      //@{

      //! Returns value
      SACADO_INLINE_FUNCTION
      const T& val() const { return val_;}

      //! Returns value
      SACADO_INLINE_FUNCTION
      T& val() { return val_;}

      //@}

      /*!
       * @name Derivative accessor methods
       */
      //@{

      //! Returns number of derivative components
      SACADO_INLINE_FUNCTION
      int size() const { return Num;}

      /*!
       * \brief Returns number of derivative components that can be stored
       * without reallocation
       */
      SACADO_INLINE_FUNCTION
      int availableSize() const { return Num; }

      //! Returns true if derivative array is not empty
      SACADO_INLINE_FUNCTION
      bool hasFastAccess() const { return true; }

      //! Returns true if derivative array is empty
      SACADO_INLINE_FUNCTION
      bool isPassive() const { return false; }

      //! Set whether variable is constant
      SACADO_INLINE_FUNCTION
      void setIsConstant(bool is_const) {}

      //! Returns derivative array
      SACADO_INLINE_FUNCTION
      const T* dx() const { return &(dx_[0]);}

      //! Returns derivative component \c i with bounds checking
      SACADO_INLINE_FUNCTION
      const T& dx(int i) const { return dx_[i]; }

      //! Returns derivative component \c i without bounds checking
      SACADO_INLINE_FUNCTION
      T& fastAccessDx(int i) { return dx_[i];}

      //! Returns derivative component \c i without bounds checking
      SACADO_INLINE_FUNCTION
      const T& fastAccessDx(int i) const { return dx_[i];}

      //! Return partials w.r.t. arguments
      SACADO_INLINE_FUNCTION
      void computePartials(const T& bar, value_type partials[]) const {
        partials[0] = bar;
      }

      //! Return tangent component \c i of arguments
      SACADO_INLINE_FUNCTION
      void getTangents(int i, value_type dots[]) const {
        dots[0] = this->dx_[i];
      }

      //! Return whether argument is active
      template <int Arg>
      SACADO_INLINE_FUNCTION
      bool isActive() const { return true; }

      //! Return tangent component \c i of argument \c Arg
      template <int Arg>
      SACADO_INLINE_FUNCTION
      const T& getTangent(int i) const { return this->dx_[i]; }

      //! Get dx array
      SACADO_INLINE_FUNCTION
      const value_type* getDx(int j) const { return this->dx(); }

      //@}

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(Expr&) operator=(const S& v) {
        val_ = v;
        ss_array<T>::zero(dx_, Num);
        return *this;
      }

      //! Assignment with Expr right-hand-side
      SACADO_INLINE_FUNCTION
      Expr& operator=(const Expr& x) {
        if (this != &x) {
          // Copy value
          val_ = x.val_;

          // Copy dx_
          //ss_array<T>::copy(x.dx_, dx_, Num);
          for (int i=0; i<Num; i++)
            dx_[i] = x.dx_[i];
        }
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(Expr&) operator=(const Expr<S>& x) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (x.size() != Num)
          throw "SELRFad::operator=() Error:  Attempt to assign with incompatible sizes";
#endif

        // Compute partials
        LocalAccumOp< Expr<S> > op(x);

        // Compute each tangent direction
        for(int i=0; i<Num; ++i) {
          op.t = T(0.);
          op.i = i;

          // Automatically unrolled loop that computes
          // for (int j=0; j<N; j++)
          //   op.t += op.partials[j] * x.getTangent<j>(i);
          Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

          dx_[i] = op.t;
        }

        // Compute value
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
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(Expr&) operator += (const S& v) {
        this->val() += v;
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(Expr&) operator -= (const S& v) {
        this->val() -= v;
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(Expr&) operator *= (const S& v) {
        this->val() *= v;
        for (int i=0; i<Num; ++i)
          dx_[i] *= v;
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(Expr&) operator /= (const S& v) {
        this->val() /= v;
        for (int i=0; i<Num; ++i)
          dx_[i] /= v;
        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(Expr&) operator += (const Expr<S>& x) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (x.size() != Num)
          throw "SELRFad::operator+=() Error:  Attempt to assign with incompatible sizes";
#endif

        // Compute partials
        LocalAccumOp< Expr<S> > op(x);

        // Compute each tangent direction
        for(int i=0; i<Num; ++i) {
          op.t = T(0.);
          op.i = i;

          // Automatically unrolled loop that computes
          // for (int j=0; j<N; j++)
          //   op.t += op.partials[j] * x.getTangent<j>(i);
          Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

          dx_[i] += op.t;
        }

        // Compute value
        val_ += x.val();

        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(Expr&) operator -= (const Expr<S>& x) {
#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (x.size() != Num)
          throw "SELRFad::operator-=() Error:  Attempt to assign with incompatible sizes";
#endif

        // Compute partials
        LocalAccumOp< Expr<S> > op(x);

        // Compute each tangent direction
        for(int i=0; i<Num; ++i) {
          op.t = T(0.);
          op.i = i;

          // Automatically unrolled loop that computes
          // for (int j=0; j<N; j++)
          //   op.t += op.partials[j] * x.getTangent<j>(i);
          Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

          dx_[i] -= op.t;
        }

        // Compute value
        val_ -= x.val();

        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(Expr&) operator *= (const Expr<S>& x) {
        T xval = x.val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (x.size() != Num)
          throw "SELRFad::operator*=() Error:  Attempt to assign with incompatible sizes";
#endif

        // Compute partials
        LocalAccumOp< Expr<S> > op(x);

        // Compute each tangent direction
        for(int i=0; i<Num; ++i) {
          op.t = T(0.);
          op.i = i;

          // Automatically unrolled loop that computes
          // for (int j=0; j<N; j++)
          //   op.t += op.partials[j] * x.getTangent<j>(i);
          Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

          dx_[i] = val_ * op.t + dx_[i] * xval;
        }

        // Compute value
        val_ *= xval;

        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(Expr&) operator /= (const Expr<S>& x) {
        T xval = x.val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if (x.size() != Num)
          throw "SELRFad::operator/=() Error:  Attempt to assign with incompatible sizes";
#endif

        // Compute partials
        LocalAccumOp< Expr<S> > op(x);

        T xval2 = xval*xval;

        // Compute each tangent direction
        for(int i=0; i<Num; ++i) {
          op.t = T(0.);
          op.i = i;

          // Automatically unrolled loop that computes
          // for (int j=0; j<N; j++)
          //   op.t += op.partials[j] * x.getTangent<j>(i);
          Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

          dx_[i] = (dx_[i] * xval - val_ * op.t) / xval2;
        }

        // Compute value
        val_ /= xval;

        return *this;
      }

      //@}

    protected:

      //! Derivatives
      T dx_[Num];

      //! Value
      T val_;

      // Functor for mpl::for_each to compute the local accumulation
      // of a tangent derivative
      //
      // We use getTangent<>() to get dx components from expression
      // arguments instead of getting the argument directly or extracting
      // the dx array due to striding in ViewFad (or could use striding
      // directly here if we need to get dx array).
      template <typename ExprT>
      struct LocalAccumOp {
        typedef typename ExprT::value_type value_type;
        static const int N = ExprT::num_args;
        const ExprT& x;
        mutable value_type t;
        value_type partials[N];
        int i;
        SACADO_INLINE_FUNCTION
        LocalAccumOp(const ExprT& x_) : x(x_) {
          x.computePartials(value_type(1.), partials);
        }
        SACADO_INLINE_FUNCTION
        void getTangents(int i_) { i = i_; }
        template <typename ArgT>
        SACADO_INLINE_FUNCTION
        void operator () (ArgT arg) const {
          const int Arg = ArgT::value;
          if (x.template isActive<Arg>())
            t += partials[Arg] * x.template getTangent<Arg>(i);
        }
      };

    }; // class Expr<SFadExprTag>

  } // namespace ELRFad

} // namespace Sacado

#define FAD_NS ELRFad
#include "Sacado_Fad_SFad_tmpl.hpp"
#undef FAD_NS

#include "Sacado_ELRFad_ViewFad.hpp"
#include "Sacado_ELRFad_Ops.hpp"

#endif // SACADO_ELRFAD_SFAD_HPP
