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
// @HEADER

#ifndef SACADO_FAD_EXP_GENERALFAD_HPP
#define SACADO_FAD_EXP_GENERALFAD_HPP

#include "Sacado_Fad_Exp_GeneralFadTraits.hpp"
#include "Sacado_Fad_Exp_Expression.hpp"
#include "Sacado_Fad_Exp_Extender.hpp"
#include "Sacado_Fad_Exp_ExprAssign.hpp"

namespace Sacado {

  //! Namespace for forward-mode AD classes
  namespace Fad {
  namespace Exp {

    //! Forward-mode AD class templated on the storage for the derivative array
    /*!
     * This class provides a general forward mode AD implementation for any
     * type of derivative array storage.
     */
    template <typename Storage>
    class GeneralFad :
      public Expr< GeneralFad<Storage> >, // Brings in expression interface
      public Extender<Storage> // Brings in interface extensions & storage
    {
    public:

      //! Storage type
      typedef Storage StorageType;

      //! Expression type
      typedef Expr< GeneralFad<Storage> > ExprType;

      //! Extender type
      typedef Extender<Storage> ExtenderType;

      //! Typename of values
      using typename StorageType::value_type;

      //! Typename of scalar's (which may be different from T)
      typedef typename ScalarType<value_type>::type scalar_type;

      //! Whether we are a view
      static constexpr bool is_view = Storage::is_view;

      //! Turn GeneralFad into a meta-function class usable with mpl::apply
      template <typename T>
      struct apply {
        typedef typename Storage::template apply<T>::type S;
        typedef GeneralFad<S> type;
      };

      //! Replace static derivative length
      template <int N>
      struct apply_N {
        typedef typename Storage::template apply_N<N>::type S;
        typedef GeneralFad<S> type;
      };

      /*!
       * @name Initialization methods
       */
      //@{

       //! Inherit Storage's and Extender's constructors
      using ExtenderType::ExtenderType;

      //! Default constructor
      GeneralFad() = default;

      //! Copy constructor
      GeneralFad(const GeneralFad& x) = default;

      //! Move constructor
      GeneralFad(GeneralFad&& x) = default;

      //! Constructor with value (disabled for ViewFad)
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      GeneralFad(const S & x, SACADO_EXP_ENABLE_VALUE_CTOR_DECL) :
        ExtenderType(x) {}

      //! Copy constructor from any Expression object (disabled for ViewFad)
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      GeneralFad(const Expr<S>& x, SACADO_EXP_ENABLE_EXPR_CTOR_DECL) :
        ExtenderType(x.derived().size(), value_type(0.), NoInitDerivArray)
      {
        ExprAssign<GeneralFad>::assign_equal(*this, x.derived());
      }

      //! Destructor
      ~GeneralFad() = default;

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
        this->fastAccessDx(ith) = value_type(1.);
      }

      //! Set whether this Fad object should update values
      /*! Retained for backward compatibility.
       */
      KOKKOS_INLINE_FUNCTION
      void setUpdateValue(bool update_val) {}

      //! Return whether this Fad object has an updated value
      /*! Retained for backward compatibility.
       */
      KOKKOS_INLINE_FUNCTION
      bool updateValue() const { return true; }

      //! Cache values
      /*! Retained for backward compatibility.
       */
      KOKKOS_INLINE_FUNCTION
      void cache() const {}

      //! Returns whether two Fad objects have the same values
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_EXP_ENABLE_EXPR_FUNC(bool) isEqualTo(const Expr<S>& xx) const {
        typedef typename Expr<S>::derived_type expr_type;
        const expr_type& x = xx.derived();

        typedef IsEqual<value_type> IE;
        if (x.size() != this->size()) return false;
        bool eq = IE::eval(x.val(), this->val());
        const int sz = this->size();
        SACADO_FAD_DERIV_LOOP(i,sz)
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
      KOKKOS_INLINE_FUNCTION
      int availableSize() const { return this->length(); }

      //! Returns true if derivative array is not empty
      KOKKOS_INLINE_FUNCTION
      bool hasFastAccess() const { return this->size()!=0; }

      //! Set whether variable is constant
      KOKKOS_INLINE_FUNCTION
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
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(GeneralFad&) operator=(const S& v) {
        this->val() = v;
        if (this->size()) this->resize(0);
        return *this;
      }

      //! Assignment with GeneralFad right-hand-side
      GeneralFad&
      operator=(const GeneralFad& x) = default;

      //! Move assignment with GeneralFad right-hand-side
      GeneralFad&
      operator=(GeneralFad&& x) = default;

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_EXP_ENABLE_EXPR_FUNC(GeneralFad&) operator=(const Expr<S>& x) {
        ExprAssign<GeneralFad>::assign_equal(*this, x.derived());
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
        SACADO_FAD_DERIV_LOOP(i,sz)
          this->fastAccessDx(i) *= v;
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(GeneralFad&) operator /= (const S& v) {
        const int sz = this->size();
        this->val() /= v;
        SACADO_FAD_DERIV_LOOP(i,sz)
          this->fastAccessDx(i) /= v;
        return *this;
      }

      //! Addition-assignment operator with GeneralFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      GeneralFad& operator += (const GeneralFad& x) {
        ExprAssign<GeneralFad>::assign_plus_equal(*this, x);
        return *this;
      }

      //! Subtraction-assignment operator with GeneralFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      GeneralFad& operator -= (const GeneralFad& x) {
        ExprAssign<GeneralFad>::assign_minus_equal(*this, x);
        return *this;
      }

      //! Multiplication-assignment operator with GeneralFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      GeneralFad& operator *= (const GeneralFad& x) {
        ExprAssign<GeneralFad>::assign_times_equal(*this, x);
        return *this;
      }

      //! Division-assignment operator with GeneralFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      GeneralFad& operator /= (const GeneralFad& x) {
        ExprAssign<GeneralFad>::assign_divide_equal(*this, x);
        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_EXP_ENABLE_EXPR_FUNC(GeneralFad&) operator += (const Expr<S>& x) {
        ExprAssign<GeneralFad>::assign_plus_equal(*this, x.derived());
        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_EXP_ENABLE_EXPR_FUNC(GeneralFad&) operator -= (const Expr<S>& x) {
        ExprAssign<GeneralFad>::assign_minus_equal(*this, x.derived());
        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_EXP_ENABLE_EXPR_FUNC(GeneralFad&) operator *= (const Expr<S>& x) {
        ExprAssign<GeneralFad>::assign_times_equal(*this, x.derived());
        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_EXP_ENABLE_EXPR_FUNC(GeneralFad&) operator /= (const Expr<S>& x) {
        ExprAssign<GeneralFad>::assign_divide_equal(*this, x.derived());
        return *this;
      }

      //@}

    }; // class GeneralFad

    template <typename S>
    struct ExprLevel< GeneralFad<S> > {
      static constexpr unsigned value =
        ExprLevel< typename GeneralFad<S>::value_type >::value + 1;
    };

    template <typename S>
    struct IsFadExpr< GeneralFad<S> > {
      static constexpr bool value = true;
    };

  } // namespace Exp
  } // namespace Fad

  template <typename S>
  struct IsView< Fad::Exp::GeneralFad<S> > {
    static constexpr bool value =  Fad::Exp::GeneralFad<S>::is_view;
  };

  template <typename S>
  struct IsFad< Fad::Exp::GeneralFad<S> > {
    static constexpr bool value = true;
  };

  template <typename S>
  struct IsExpr< Fad::Exp::GeneralFad<S> > {
    static constexpr bool value = true;
  };

  template <typename S>
  struct BaseExprType< Fad::Exp::GeneralFad<S> > {
    typedef Fad::Exp::GeneralFad<S> type;
  };

} // namespace Sacado

#include "Sacado_Fad_Exp_Ops.hpp"

#endif // SACADO_FAD_EXP_GENERALFAD_HPP
