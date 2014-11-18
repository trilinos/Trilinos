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

#ifndef SACADO_ELRCACHEFAD_VIEWFAD_HPP
#define SACADO_ELRCACHEFAD_VIEWFAD_HPP

#include "Sacado_ELRCacheFad_GeneralFadExpr.hpp"
#include "Sacado_ELRCacheFad_ViewFadTraits.hpp"
#include "Sacado_Fad_ViewStorage.hpp"
#include "Sacado_dummy_arg.hpp"

namespace Sacado {

  namespace ELRCacheFad {

    /*!
     * \brief Forward-mode AD class using dynamic memory allocation and
     * expression templates.
     */
    /*!
     * This is the user-level class for forward mode AD with dynamic
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is not known at compile time.  The user
     * interface is provided by Sacado::Fad::GeneralFad.
     */
    template <typename ValueT, unsigned length, unsigned stride>
    class ViewFad :
      public Expr< GeneralFad<ValueT,Fad::ViewStorage<ValueT,length,stride> > > {

    public:

      //! Base classes
      typedef Fad::ViewStorage<ValueT,length,stride> StorageType;
      typedef GeneralFad<ValueT,StorageType> GeneralFadType;
      typedef Expr<GeneralFadType> ExprType;

      //! Typename of scalar's (which may be different from ValueT)
      typedef typename ScalarType<ValueT>::type ScalarT;

      //! Turn ViewFad into a meta-function class usable with mpl::apply
      template <typename T>
      struct apply {
        typedef ViewFad<T,length,stride> type;
      };

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      KOKKOS_INLINE_FUNCTION
      ViewFad() :
        ExprType() {}

      //! Constructor with supplied value \c x of type ValueT
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      KOKKOS_INLINE_FUNCTION
      ViewFad(const ValueT& x) :
        ExprType(x) {}

      //! Constructor with supplied value \c x of type ScalarT
      /*!
       * Initializes value to \c ValueT(x) and derivative array is empty.
       * Creates a dummy overload when ValueT and ScalarT are the same type.
       */
      KOKKOS_INLINE_FUNCTION
      ViewFad(const typename dummy<ValueT,ScalarT>::type& x) :
        ExprType(ValueT(x)) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      KOKKOS_INLINE_FUNCTION
      ViewFad(const int sz, const ValueT& x) :
        ExprType(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      KOKKOS_INLINE_FUNCTION
      ViewFad(const int sz, const int i, const ValueT & x) :
        ExprType(sz,i,x) {}

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      ViewFad(const ViewFad& x) :
        ExprType(x) {}

      //! Copy constructor from any Expression object
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      ViewFad(const Expr<S>& x) :
        ExprType(x) {}

      //! Constructor with supplied storage \c s
      KOKKOS_INLINE_FUNCTION
      ViewFad(const StorageType& s) :
        ExprType(s) {}

      //! View-specific constructor
      KOKKOS_INLINE_FUNCTION
      ViewFad(ValueT* v, const int arg_size = 0, const int arg_stride = 0) :
        ExprType( StorageType(v,arg_size,arg_stride) ) {}

      //@}

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~ViewFad() {}

      //! Assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator=(const ValueT& v) {
        GeneralFadType::operator=(v);
        return *this;
      }

      //! Assignment operator with constant right-hand-side
      /*!
       * Creates a dummy overload when ValueT and ScalarT are the same type.
       */
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator=(const typename dummy<ValueT,ScalarT>::type& v) {
        GeneralFadType::operator=(ValueT(v));
        return *this;
      }

      //! Assignment operator with ViewFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator=(const ViewFad& x) {
        GeneralFadType::operator=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator=(const Expr<S>& x)
      {
        GeneralFadType::operator=(x);
        return *this;
      }

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Addition-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator += (const ValueT& x) {
        GeneralFadType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator -= (const ValueT& x) {
        GeneralFadType::operator-=(x);
        return *this;
      };

      //! Multiplication-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator *= (const ValueT& x) {
        GeneralFadType::operator*=(x);
        return *this;
      };

      //! Division-assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator /= (const ValueT& x) {
        GeneralFadType::operator/=(x);
        return *this;
      };

      //! Addition-assignment operator with constant right-hand-side
      /*!
       * Creates a dummy overload when value_type and scalar_type are the
       * same type.
       */
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator += (const typename dummy<ValueT,ScalarT>::type& x) {
        GeneralFadType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      /*!
       * Creates a dummy overload when value_type and scalar_type are the
       * same type.
       */
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator -= (const typename dummy<ValueT,ScalarT>::type& x) {
        GeneralFadType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      /*!
       * Creates a dummy overload when value_type and scalar_type are the
       * same type.
       */
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator *= (const typename dummy<ValueT,ScalarT>::type& x) {
        GeneralFadType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      /*!
       * Creates a dummy overload when value_type and scalar_type are the
       * same type.
       */
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator /= (const typename dummy<ValueT,ScalarT>::type& x) {
        GeneralFadType::operator/=(x);
        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator += (const Expr<S>& x) {
        GeneralFadType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator -= (const Expr<S>& x) {
        GeneralFadType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator *= (const Expr<S>& x) {
        GeneralFadType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      ViewFad& operator /= (const Expr<S>& x) {
        GeneralFadType::operator/=(x);
        return *this;
      }

      //@}

    }; // class ViewFad<ValueT>

    template <typename T, unsigned l, unsigned s>
    struct BaseExpr< GeneralFad<T,Fad::ViewStorage<T,l,s> > > {
      typedef ViewFad<T,l,s> type;
    };

  } // namespace ELRCacheFad

} // namespace Sacado

#endif // SACADO_ELRCACHEFAD_VIEWFAD_HPP
