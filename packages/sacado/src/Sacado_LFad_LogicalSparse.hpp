// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_LFAD_LOGICALSPARSE_HPP
#define SACADO_LFAD_LOGICALSPARSE_HPP

#include "Sacado_LFad_LogicalSparseTraits.hpp"
#include "Sacado_LFad_ExpressionTraits.hpp"
#include "Sacado_Fad_DynamicStorage.hpp"

namespace Sacado {

  //! Namespace for logical forward-mode AD classes
  namespace LFad {

    //! Wrapper for a generic expression template
    /*!
     * This template class serves as a wrapper for all Fad expression
     * template classes.
     */
    template <typename ExprT> class Expr {};

    //! Meta-function for determining nesting with an expression
    /*!
     * This determines the level of nesting within nested Fad types.
     * The default implementation works for any type that isn't a Fad type
     * or an expression of Fad types.
     */
    template <typename T>
    struct ExprLevel {
      static const unsigned value = 0;
    };

    template <typename T>
    struct ExprLevel< Expr<T> > {
      static const unsigned value =
        ExprLevel< typename Expr<T>::value_type >::value + 1;
    };

    //! Determine whether a given type is an expression
    template <typename T>
    struct IsFadExpr {
      static const bool value = false;
    };

    template <typename T>
    struct IsFadExpr< Expr<T> > {
      static const bool value = true;
    };

    // Forward declaration
    template <typename ValT, typename LogT> class LogicalSparse;

    /*!
     * \brief Implementation class for computing the logical sparsity of a
     * derivative using forward-mode AD.
     */
    template <typename ValT, typename LogT>
    class LogicalSparseImp :
      public Fad::DynamicStorage<ValT,LogT> {

      typedef Fad::DynamicStorage<ValT,LogT> Storage;

    public:

      //! Typename of values (e.g., double)
      typedef ValT value_type;

      //! Typename of scalar's (which may be different from ValT)
      typedef typename ScalarType<value_type>::type scalar_type;

      //! Logical type (i.e., type for derivative array components (e.g., bool)
      typedef LogT logical_type;

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor
      LogicalSparseImp() : Storage(value_type(0)) {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      template <typename S>
      LogicalSparseImp(const S& x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        Storage(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      LogicalSparseImp(const int sz, const value_type & x) :
        Storage(sz, x, InitDerivArray) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      LogicalSparseImp(const int sz, const int i, const value_type & x) :
        Storage(sz, x, InitDerivArray) {
        this->fastAccessDx(i)=logical_type(1);
      }

      //! Copy constructor
      LogicalSparseImp(const LogicalSparseImp& x) :
        Storage(x) {}

      //! Copy constructor from any Expression object
      template <typename S>
      LogicalSparseImp(const Expr<S>& x, SACADO_ENABLE_EXPR_CTOR_DECL)  :
        Storage(value_type(0)) {
        int sz = x.size();

        if (sz != this->size())
          this->resize(sz);

        if (sz) {
          if (x.hasFastAccess())
            for(int i=0; i<sz; ++i)
              this->fastAccessDx(i) = x.fastAccessDx(i);
          else
            for(int i=0; i<sz; ++i)
              this->fastAccessDx(i) = x.dx(i);
        }

        this->val() = x.val();
      }

      //! Destructor
      ~LogicalSparseImp() {}

      //! Set %LogicalSparseImp object as the \c ith independent variable
      /*!
       * Sets the derivative array of length \c n to the \c ith row of the
       * identity matrix and has the same affect as the
       * Implementation(const int sz, const int i, const T & x)
       * constructor.
       */
      void diff(const int ith, const int n) {
        if (this->size() != n)
          this->resize(n);

        this->zero();
        this->fastAccessDx(ith) = logical_type(1);
      }

      //! Returns whether two Fad objects have the same values
      template <typename S>
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

      //! Returns true if derivative array is not empty
      bool hasFastAccess() const { return this->size()!=0;}

      //! Returns true if derivative array is empty
      bool isPassive() const { return this->size()!=0;}

      //! Set whether variable is constant
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
      SACADO_ENABLE_VALUE_FUNC(LogicalSparseImp&) operator=(const S& v)  {
        this->val() = v;
        if (this->size()) {
          this->zero();
          this->resize(0);
        }
        return *this;
      }

      //! Assignment with Expr right-hand-side
      LogicalSparseImp& operator=(const LogicalSparseImp& x) {
        // Copy value & dx_
        Storage::operator=(x);

        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      SACADO_ENABLE_EXPR_FUNC(LogicalSparseImp&) operator=(const Expr<S>& x) {
        int sz = x.size();

        if (sz != this->size())
          this->resize(sz);

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
      SACADO_ENABLE_VALUE_FUNC(LogicalSparseImp&) operator += (const S& v) {
        this->val() += v;
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_ENABLE_VALUE_FUNC(LogicalSparseImp&) operator -= (const S& v) {
        this->val() -= v;
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_ENABLE_VALUE_FUNC(LogicalSparseImp&) operator *= (const S& v) {
        this->val() *= v;
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_ENABLE_VALUE_FUNC(LogicalSparseImp&) operator /= (const S& v) {
        this->val() /= v;
        return *this;
      }

      //! Addition-assignment operator with LogicalSparseImp right-hand-side
      LogicalSparseImp& operator += (const LogicalSparseImp& x) {
        int xsz = x.size(), sz = this->size();

#ifdef SACADO_DEBUG
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "LFad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            for (int i=0; i<xsz; ++i)
              this->fastAccessDx(i) = this->fastAccessDx(i) || x.fastAccessDx(i);
          }
          else {
            this->resize(xsz);
            for (int i=0; i<xsz; ++i)
              this->fastAccessDx(i) = x.fastAccessDx(i);
          }
        }

        this->val() += x.val();

        return *this;
      }

      //! Subtraction-assignment operator with LogicalSparseImp right-hand-side
      LogicalSparseImp& operator -= (const LogicalSparseImp& x) {
        int xsz = x.size(), sz = this->size();

#ifdef SACADO_DEBUG
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "LFad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            for (int i=0; i<xsz; ++i)
              this->fastAccessDx(i) = this->fastAccessDx(i) || x.fastAccessDx(i);
          }
          else {
            this->resize(xsz);
            for (int i=0; i<xsz; ++i)
              this->fastAccessDx(i) = x.fastAccessDx(i);
          }
        }

        this->val() -= x.val();


        return *this;
      }

      //! Multiplication-assignment operator with LogicalSparseImp right-hand-side
      LogicalSparseImp& operator *= (const LogicalSparseImp& x) {
        int xsz = x.size(), sz = this->size();

#ifdef SACADO_DEBUG
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "LFad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            for (int i=0; i<xsz; ++i)
              this->fastAccessDx(i) = this->fastAccessDx(i) || x.fastAccessDx(i);
          }
          else {
            this->resize(xsz);
            for (int i=0; i<xsz; ++i)
              this->fastAccessDx(i) = x.fastAccessDx(i);
          }
        }

        this->val() *= x.val();

        return *this;
      }

      //! Division-assignment operator with LogicalSparseImp right-hand-side
      LogicalSparseImp& operator /= (const LogicalSparseImp& x) {
        int xsz = x.size(), sz = this->size();

#ifdef SACADO_DEBUG
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "LFad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            for (int i=0; i<xsz; ++i)
              this->fastAccessDx(i) = this->fastAccessDx(i) || x.fastAccessDx(i);
          }
          else {
            this->resize(xsz);
            for (int i=0; i<xsz; ++i)
              this->fastAccessDx(i) = x.fastAccessDx(i);
          }
        }

        this->val() /= x.val();

        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_ENABLE_EXPR_FUNC(LogicalSparseImp&) operator += (const Expr<S>& x){
        int xsz = x.size(), sz = this->size();

#ifdef SACADO_DEBUG
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "LFad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            if (x.hasFastAccess())
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = this->fastAccessDx(i) || x.fastAccessDx(i);
            else
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = this->fastAccessDx(i) || x.dx(i);
          }
          else {
            this->resize(xsz);
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
      SACADO_ENABLE_EXPR_FUNC(LogicalSparseImp&) operator -= (const Expr<S>& x){
        int xsz = x.size(), sz = this->size();

#ifdef SACADO_DEBUG
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "LFad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            if (x.hasFastAccess())
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = this->fastAccessDx(i) || x.fastAccessDx(i);
            else
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = this->fastAccessDx(i) || x.dx(i);
          }
          else {
            this->resize(xsz);
            if (x.hasFastAccess())
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = x.fastAccessDx(i);
            else
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = x.dx(i);
          }
        }

        this->val() -= x.val();


        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_ENABLE_EXPR_FUNC(LogicalSparseImp&) operator *= (const Expr<S>& x){
        int xsz = x.size(), sz = this->size();

#ifdef SACADO_DEBUG
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "LFad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            if (x.hasFastAccess())
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = this->fastAccessDx(i) || x.fastAccessDx(i);
            else
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = this->fastAccessDx(i) || x.dx(i);
          }
          else {
            this->resize(xsz);
            if (x.hasFastAccess())
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = x.fastAccessDx(i);
            else
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = x.dx(i);
          }
        }

        this->val() *= x.val();

        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_ENABLE_EXPR_FUNC(LogicalSparseImp&) operator /= (const Expr<S>& x){
        int xsz = x.size(), sz = this->size();

#ifdef SACADO_DEBUG
        if ((xsz != sz) && (xsz != 0) && (sz != 0))
          throw "LFad Error:  Attempt to assign with incompatible sizes";
#endif

        if (xsz) {
          if (sz) {
            if (x.hasFastAccess())
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = this->fastAccessDx(i) || x.fastAccessDx(i);
            else
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = this->fastAccessDx(i) || x.dx(i);
          }
          else {
            this->resize(xsz);
            if (x.hasFastAccess())
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = x.fastAccessDx(i);
            else
              for (int i=0; i<xsz; ++i)
                this->fastAccessDx(i) = x.dx(i);
          }
        }

        this->val() /= x.val();

        return *this;
      }

      //@}

    }; // class LogicalSparseImp

    //! Expression template specialization for LogicalSparse
    template <typename ValT, typename LogT>
    class Expr< LogicalSparseImp<ValT,LogT> > :
      public LogicalSparseImp<ValT,LogT> {

    public:

      //! Typename of values
      typedef typename LogicalSparseImp<ValT,LogT>::value_type value_type;

      //! Typename of scalar's (which may be different from value_type)
      typedef typename LogicalSparseImp<ValT,LogT>::scalar_type scalar_type;

      //! Typename of base-expressions
      typedef LogicalSparse<ValT,LogT> base_expr_type;

      //! Default constructor
      Expr() :
        LogicalSparseImp<ValT,LogT>() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      template <typename S>
      Expr(const S & x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        LogicalSparseImp<ValT,LogT>(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      Expr(const int sz, const ValT & x) :
        LogicalSparseImp<ValT,LogT>(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      Expr(const int sz, const int i, const ValT & x) :
        LogicalSparseImp<ValT,LogT>(sz,i,x) {}

      //! Copy constructor
      Expr(const Expr& x) :
        LogicalSparseImp<ValT,LogT>(x) {}

      //! Copy constructor from any Expression object
      template <typename S>
      Expr(const Expr<S>& x, SACADO_ENABLE_EXPR_CTOR_DECL) :
        LogicalSparseImp<ValT,LogT>(x) {}

      //! Destructor
      ~Expr() {}

    }; // class Expr<LogicalSparseImp>

    /*!
     * \brief User inteface class for computing the logical sparsity pattern
     * of a derivative via forward-mode AD.
     */
    template <typename ValT, typename LogT >
    class LogicalSparse : public Expr< LogicalSparseImp<ValT,LogT > > {

    public:

      //! Base classes
      typedef LogicalSparseImp< ValT,LogT > ImplType;
      typedef Expr<ImplType> ExprType;

      //! Typename of values
      typedef typename ExprType::value_type value_type;

      //! Typename of scalar's (which may be different from value_type)
      typedef typename ExprType::scalar_type scalar_type;

      //! Turn LogicalSparse into a meta-function class usable with mpl::apply
      template <typename T, typename U = LogT>
      struct apply {
        typedef LogicalSparse<T,U> type;
      };

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      LogicalSparse() :
        ExprType() {}

      //! Constructor with supplied value \c x of type ValueT
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      template <typename S>
      LogicalSparse(const S& x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        ExprType(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      LogicalSparse(const int sz, const ValT& x) :
        ExprType(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      LogicalSparse(const int sz, const int i, const ValT & x) :
        ExprType(sz,i,x) {}

      //! Copy constructor
      LogicalSparse(const LogicalSparse& x) :
        ExprType(x) {}

      //! Copy constructor from any Expression object
      template <typename S>
      LogicalSparse(const Expr<S>& x, SACADO_ENABLE_EXPR_CTOR_DECL) :
        ExprType(x) {}

      //@}

      //! Destructor
      ~LogicalSparse() {}

      //! Assignment operator with constant right-hand-side
      template <typename S>
      SACADO_ENABLE_VALUE_FUNC(LogicalSparse&) operator=(const S& v) {
        ImplType::operator=(v);
        return *this;
      }

      //! Assignment operator with LogicalSparse right-hand-side
      LogicalSparse& operator=(const LogicalSparse& x) {
        ImplType::operator=(static_cast<const ImplType&>(x));
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      SACADO_ENABLE_EXPR_FUNC(LogicalSparse&) operator=(const Expr<S>& x)
        {
          ImplType::operator=(x);
          return *this;
        }

      //! Addition-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_ENABLE_VALUE_FUNC(LogicalSparse&) operator += (const S& x) {
        ImplType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_ENABLE_VALUE_FUNC(LogicalSparse&) operator -= (const S& x) {
        ImplType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_ENABLE_VALUE_FUNC(LogicalSparse&) operator *= (const S& x) {
        ImplType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      SACADO_ENABLE_VALUE_FUNC(LogicalSparse&) operator /= (const S& x) {
        ImplType::operator/=(x);
        return *this;
      }

      //! Addition-assignment operator with LogicalSparse right-hand-side
      LogicalSparse& operator += (const LogicalSparse& x) {
        ImplType::operator+=(static_cast<const ImplType&>(x));
        return *this;
      }

      //! Subtraction-assignment operator with LogicalSparse right-hand-side
      LogicalSparse& operator -= (const LogicalSparse& x) {
        ImplType::operator-=(static_cast<const ImplType&>(x));
        return *this;
      }

      //! Multiplication-assignment operator with LogicalSparse right-hand-side
      LogicalSparse& operator *= (const LogicalSparse& x) {
        ImplType::operator*=(static_cast<const ImplType&>(x));
        return *this;
      }

      //! Division-assignment operator with LogicalSparse right-hand-side
      LogicalSparse& operator /= (const LogicalSparse& x) {
        ImplType::operator/=(static_cast<const ImplType&>(x));
        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_ENABLE_EXPR_FUNC(LogicalSparse&) operator += (const Expr<S>& x) {
        ImplType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_ENABLE_EXPR_FUNC(LogicalSparse&) operator -= (const Expr<S>& x) {
        ImplType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_ENABLE_EXPR_FUNC(LogicalSparse&) operator *= (const Expr<S>& x) {
        ImplType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      SACADO_ENABLE_EXPR_FUNC(LogicalSparse&) operator /= (const Expr<S>& x) {
        ImplType::operator/=(x);
        return *this;
      }

    }; // class LogicalSparse<ValT,LogT>

    template <typename T, typename L>
    struct ExprLevel< LogicalSparse<T,L> > {
      static const unsigned value =
        ExprLevel< typename LogicalSparse<T,L>::value_type >::value + 1;
    };

  } // namespace LFad

  template <typename T>
  struct IsExpr< LFad::Expr<T> > {
    static const bool value = true;
  };

  template <typename T>
  struct BaseExprType< LFad::Expr<T> > {
    typedef typename LFad::Expr<T>::base_expr_type type;
  };

  template <typename T, typename L>
  struct IsExpr< LFad::LogicalSparse<T,L> > {
    static const bool value = true;
  };

  template <typename T, typename L>
  struct BaseExprType< LFad::LogicalSparse<T,L> > {
    typedef typename LFad::LogicalSparse<T,L>::base_expr_type type;
  };

} // namespace Sacado

#include "Sacado_LFad_LogicalSparseOps.hpp"

#endif // SACADO_LFAD_LOGICALSPARSE_HPP
