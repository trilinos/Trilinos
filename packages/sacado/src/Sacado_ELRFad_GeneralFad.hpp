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

#ifndef SACADO_ELRFAD_GENERALFAD_HPP
#define SACADO_ELRFAD_GENERALFAD_HPP

#include "Sacado_ELRFad_Expression.hpp"
#include "Sacado_dummy_arg.hpp"
#include "Sacado_mpl_range_c.hpp"
#include "Sacado_mpl_for_each.hpp"
#include<ostream>

namespace Sacado {

  //! Namespace for expression-level reverse forward-mode AD classes
  namespace ELRFad {

    //! Forward-mode AD class templated on the storage for the derivative array
    /*!
     * This class provides a general forward mode AD implementation for any
     * type of derivative array storage.  It does not incorporate expression
     * templates.
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
        const int sz = x.size();
        if (sz) {

          if (Expr<S>::is_linear) {
            if (x.hasFastAccess())
              for(int i=0; i<sz; ++i)
                this->fastAccessDx(i) = x.fastAccessDx(i);
            else
              for(int i=0; i<sz; ++i)
                this->fastAccessDx(i) = x.dx(i);
          }
          else {

            if (x.hasFastAccess()) {
              // Compute partials
              FastLocalAccumOp< Expr<S> > op(x);

              // Compute each tangent direction
              for(op.i=0; op.i<sz; ++op.i) {
                op.t = T(0.);

                // Automatically unrolled loop that computes
                // for (int j=0; j<N; j++)
                //   op.t += op.partials[j] * x.getTangent<j>(i);
                Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

                this->fastAccessDx(op.i) = op.t;
              }
            }
            else {
              // Compute partials
              SlowLocalAccumOp< Expr<S> > op(x);

              // Compute each tangent direction
              for(op.i=0; op.i<sz; ++op.i) {
                op.t = T(0.);

                // Automatically unrolled loop that computes
                // for (int j=0; j<N; j++)
                //   op.t += op.partials[j] * x.getTangent<j>(i);
                Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

                this->fastAccessDx(op.i) = op.t;
              }
            }

          }

        }

        // Compute value
        this->val() = x.val();
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

      //! Assignment with GeneralFad right-hand-side
      SACADO_INLINE_FUNCTION
      GeneralFad&
      operator=(const GeneralFad& x) {
        // Copy value & dx_
        Storage::operator=(x);
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
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

          if (Expr<S>::is_linear) {
            if (x.hasFastAccess())
              for(int i=0; i<sz; ++i)
                this->fastAccessDx(i) = x.fastAccessDx(i);
            else
              for(int i=0; i<sz; ++i)
                this->fastAccessDx(i) = x.dx(i);
          }
          else {

            if (x.hasFastAccess()) {
              // Compute partials
              FastLocalAccumOp< Expr<S> > op(x);

              // Compute each tangent direction
              for(op.i=0; op.i<sz; ++op.i) {
                op.t = T(0.);

                // Automatically unrolled loop that computes
                // for (int j=0; j<N; j++)
                //   op.t += op.partials[j] * x.getTangent<j>(i);
                Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

                this->fastAccessDx(op.i) = op.t;
              }
            }
            else {
              // Compute partials
              SlowLocalAccumOp< Expr<S> > op(x);

              // Compute each tangent direction
              for(op.i=0; op.i<sz; ++op.i) {
                op.t = T(0.);

                // Automatically unrolled loop that computes
                // for (int j=0; j<N; j++)
                //   op.t += op.partials[j] * x.getTangent<j>(i);
                Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

                this->fastAccessDx(op.i) = op.t;
              }
            }
          }
        }

        // Compute value
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
        const int xsz = x.size(), sz = this->size();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (Expr<S>::is_linear) {
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
        }
        else {

          if (xsz) {

            if (sz != xsz)
              this->resizeAndZero(xsz);

            if (x.hasFastAccess()) {
              // Compute partials
              FastLocalAccumOp< Expr<S> > op(x);

              // Compute each tangent direction
              for(op.i=0; op.i<xsz; ++op.i) {
                op.t = T(0.);

                // Automatically unrolled loop that computes
                // for (int j=0; j<N; j++)
                //   op.t += op.partials[j] * x.getTangent<j>(i);
                Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

                this->fastAccessDx(op.i) += op.t;
              }
            }
            else {
              // Compute partials
              SlowLocalAccumOp< Expr<S> > op(x);

              // Compute each tangent direction
              for(op.i=0; op.i<xsz; ++op.i) {
                op.t = T(0.);

                // Automatically unrolled loop that computes
                // for (int j=0; j<N; j++)
                //   op.t += op.partials[j] * x.getTangent<j>(i);
                Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

                this->fastAccessDx(op.i) += op.t;
              }
            }

          }

        }

        // Compute value
        this->val() += x.val();

        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(GeneralFad&) operator -= (const Expr<S>& x) {
        const int xsz = x.size(), sz = this->size();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (Expr<S>::is_linear) {
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
        }
        else {

          if (xsz) {

            if (sz != xsz)
              this->resizeAndZero(xsz);

            if (x.hasFastAccess()) {
              // Compute partials
              FastLocalAccumOp< Expr<S> > op(x);

              // Compute each tangent direction
              for(op.i=0; op.i<xsz; ++op.i) {
                op.t = T(0.);

                // Automatically unrolled loop that computes
                // for (int j=0; j<N; j++)
                //   op.t += op.partials[j] * x.getTangent<j>(i);
                Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

                this->fastAccessDx(op.i) -= op.t;
              }
            }
            else {
              // Compute partials
              SlowLocalAccumOp< Expr<S> > op(x);

              // Compute each tangent direction
              for(op.i=0; op.i<xsz; ++op.i) {
                op.t = T(0.);

                // Automatically unrolled loop that computes
                // for (int j=0; j<N; j++)
                //   op.t += op.partials[j] * x.getTangent<j>(i);
                Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

                this->fastAccessDx(op.i) -= op.t;
              }
            }
          }

        }

        this->val() -= x.val();

        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(GeneralFad&) operator *= (const Expr<S>& x) {
        const int xsz = x.size(), sz = this->size();
        T xval = x.val();
        T v = this->val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (Expr<S>::is_linear) {
          if (xsz) {
            if (sz) {
              if (x.hasFastAccess())
                for(int i=0; i<sz; ++i)
                  this->fastAccessDx(i) = v*x.fastAccessDx(i) + this->fastAccessDx(i)*xval;
              else
                for (int i=0; i<sz; ++i)
                  this->fastAccessDx(i) = v*x.dx(i) + this->fastAccessDx(i)*xval;
            }
            else {
              this->resizeAndZero(xsz);
              if (x.hasFastAccess())
                for(int i=0; i<xsz; ++i)
                  this->fastAccessDx(i) = v*x.fastAccessDx(i);
              else
                for (int i=0; i<xsz; ++i)
                  this->fastAccessDx(i) = v*x.dx(i);
            }
          }
          else {
            if (sz) {
              for (int i=0; i<sz; ++i)
                this->fastAccessDx(i) *= xval;
            }
          }
        }
        else {

          if (xsz) {

            if (sz) {

              if (x.hasFastAccess()) {
                // Compute partials
                FastLocalAccumOp< Expr<S> > op(x);

                // Compute each tangent direction
                for(op.i=0; op.i<xsz; ++op.i) {
                  op.t = T(0.);

                  // Automatically unrolled loop that computes
                  // for (int j=0; j<N; j++)
                  //   op.t += op.partials[j] * x.getTangent<j>(i);
                  Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

                  this->fastAccessDx(op.i) =
                    v * op.t + this->fastAccessDx(op.i) * xval;
                }
              }
              else {
                // Compute partials
                SlowLocalAccumOp< Expr<S> > op(x);

                // Compute each tangent direction
                for(op.i=0; op.i<xsz; ++op.i) {
                  op.t = T(0.);

                  // Automatically unrolled loop that computes
                  // for (int j=0; j<N; j++)
                  //   op.t += op.partials[j] * x.getTangent<j>(i);
                  Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

                  this->fastAccessDx(op.i) =
                    v * op.t + this->fastAccessDx(op.i) * xval;
                }
              }

            }

            else {

              this->resizeAndZero(xsz);

              if (x.hasFastAccess()) {
                // Compute partials
                FastLocalAccumOp< Expr<S> > op(x);

                // Compute each tangent direction
                for(op.i=0; op.i<xsz; ++op.i) {
                  op.t = T(0.);

                  // Automatically unrolled loop that computes
                  // for (int j=0; j<N; j++)
                  //   op.t += op.partials[j] * x.getTangent<j>(i);
                  Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

                  this->fastAccessDx(op.i) = v * op.t;
                }
              }
              else {
                // Compute partials
                SlowLocalAccumOp< Expr<S> > op(x);

                // Compute each tangent direction
                for(op.i=0; op.i<xsz; ++op.i) {
                  op.t = T(0.);

                  // Automatically unrolled loop that computes
                  // for (int j=0; j<N; j++)
                  //   op.t += op.partials[j] * x.getTangent<j>(i);
                  Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

                  this->fastAccessDx(op.i) = v * op.t;
                }
              }

            }

          }

          else {

            if (sz) {
              for (int i=0; i<sz; ++i)
                this->fastAccessDx(i) *= xval;
            }

          }

        }

        this->val() *= xval;

        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(GeneralFad&) operator /= (const Expr<S>& x) {
        const int xsz = x.size(), sz = this->size();
        T xval = x.val();
        T v = this->val();

#if defined(SACADO_DEBUG) && !defined(__CUDA_ARCH__ )
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

        if (Expr<S>::is_linear) {
          if (xsz) {
            if (sz) {
              if (x.hasFastAccess())
                for(int i=0; i<sz; ++i)
                  this->fastAccessDx(i) = ( this->fastAccessDx(i)*xval - v*x.fastAccessDx(i) )/ (xval*xval);
              else
                for (int i=0; i<sz; ++i)
                  this->fastAccessDx(i) = ( this->fastAccessDx(i)*xval - v*x.dx(i) )/ (xval*xval);
            }
            else {
              this->resizeAndZero(xsz);
              if (x.hasFastAccess())
                for(int i=0; i<xsz; ++i)
                  this->fastAccessDx(i) = - v*x.fastAccessDx(i) / (xval*xval);
              else
                for (int i=0; i<xsz; ++i)
                  this->fastAccessDx(i) = -v*x.dx(i) / (xval*xval);
            }
          }
          else {
            if (sz) {
              for (int i=0; i<sz; ++i)
                this->fastAccessDx(i) /= xval;
            }
          }
        }
        else {

          if (xsz) {

            T xval2 = xval*xval;

            if (sz) {

              if (x.hasFastAccess()) {
                // Compute partials
                FastLocalAccumOp< Expr<S> > op(x);

                // Compute each tangent direction
                for(op.i=0; op.i<xsz; ++op.i) {
                  op.t = T(0.);

                  // Automatically unrolled loop that computes
                  // for (int j=0; j<N; j++)
                  //   op.t += op.partials[j] * x.getTangent<j>(i);
                  Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

                  this->fastAccessDx(op.i) =
                    (this->fastAccessDx(op.i) * xval - v * op.t) / xval2;
                }
              }
              else {
                // Compute partials
                SlowLocalAccumOp< Expr<S> > op(x);

                // Compute each tangent direction
                for(op.i=0; op.i<xsz; ++op.i) {
                  op.t = T(0.);

                  // Automatically unrolled loop that computes
                  // for (int j=0; j<N; j++)
                  //   op.t += op.partials[j] * x.getTangent<j>(i);
                  Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

                  this->fastAccessDx(op.i) =
                    (this->fastAccessDx(op.i) * xval - v * op.t) / xval2;
                }
              }

            }

            else {

              this->resizeAndZero(xsz);

              if (x.hasFastAccess()) {
                // Compute partials
                FastLocalAccumOp< Expr<S> > op(x);

                // Compute each tangent direction
                for(op.i=0; op.i<xsz; ++op.i) {
                  op.t = T(0.);

                  // Automatically unrolled loop that computes
                  // for (int j=0; j<N; j++)
                  //   op.t += op.partials[j] * x.getTangent<j>(i);
                  Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

                  this->fastAccessDx(op.i) = (-v * op.t) / xval2;
                }
              }
              else {
                // Compute partials
                SlowLocalAccumOp< Expr<S> > op(x);

                // Compute each tangent direction
                for(op.i=0; op.i<xsz; ++op.i) {
                  op.t = T(0.);

                  // Automatically unrolled loop that computes
                  // for (int j=0; j<N; j++)
                  //   op.t += op.partials[j] * x.getTangent<j>(i);
                  Sacado::mpl::for_each< mpl::range_c< int, 0, Expr<S>::num_args > > f(op);

                  this->fastAccessDx(op.i) = (-v * op.t) / xval2;
                }
              }

            }

          }

          else {

            if (sz) {
              for (int i=0; i<sz; ++i)
                this->fastAccessDx(i) /= xval;
            }

          }

        }

        this->val() /= xval;

        return *this;
      }

      //@}

    protected:

      // Functor for mpl::for_each to compute the local accumulation
      // of a tangent derivative
      //
      // We use getTangent<>() to get dx components from expression
      // arguments instead of getting the argument directly or extracting
      // the dx array due to striding in ViewFad (or could use striding
      // directly here if we need to get dx array).
      template <typename ExprT>
      struct FastLocalAccumOp {
        typedef typename ExprT::value_type value_type;
        static const int N = ExprT::num_args;
        const ExprT& x;
        mutable value_type t;
        value_type partials[N];
        int i;
        SACADO_INLINE_FUNCTION
        FastLocalAccumOp(const ExprT& x_) : x(x_) {
          x.computePartials(value_type(1.), partials);
        }
        template <typename ArgT>
        SACADO_INLINE_FUNCTION
        void operator () (ArgT arg) const {
          const int Arg = ArgT::value;
          t += partials[Arg] * x.template getTangent<Arg>(i);
        }
      };

      template <typename ExprT>
      struct SlowLocalAccumOp : FastLocalAccumOp<ExprT> {
        SACADO_INLINE_FUNCTION
        SlowLocalAccumOp(const ExprT& x_) :
          FastLocalAccumOp<ExprT>(x_) {}
        template <typename ArgT>
        SACADO_INLINE_FUNCTION
        void operator () (ArgT arg) const {
          const int Arg = ArgT::value;
          if (this->x.template isActive<Arg>())
            this->t += this->partials[Arg] * this->x.template getTangent<Arg>(this->i);
        }
      };

    }; // class GeneralFad


    template <typename T, typename Storage>
    std::ostream& operator << (std::ostream& os,
                               const GeneralFad<T,Storage>& x) {
      os << x.val() << " [";

      for (int i=0; i< x.size(); i++) {
        os << " " << x.dx(i);
      }

      os << " ]";
      return os;
    }

  } // namespace ELRFad

} // namespace Sacado

#endif // SACADO_ELRFAD_GENERALFAD_HPP
