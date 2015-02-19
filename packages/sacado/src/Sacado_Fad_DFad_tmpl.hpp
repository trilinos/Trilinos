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

namespace Sacado {

  namespace FAD_NS {

    /*!
     * \brief Forward-mode AD class using dynamic memory allocation and
     * expression templates.
     */
    /*!
     * This is the user-level class for forward mode AD with dynamic
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is not known at compile time.  The user
     * interface is provided by Sacado::FAD_NS::GeneralFad.
     */
    template <typename ValueT>
    class DFad :
      public Expr< GeneralFad<ValueT,Fad::DynamicStorage<ValueT> > > {

    public:

      //! Base classes
      typedef Fad::DynamicStorage<ValueT> StorageType;
      typedef GeneralFad<ValueT,StorageType> GeneralFadType;
      typedef Expr<GeneralFadType> ExprType;

      //! Typename of values
      typedef typename ExprType::value_type value_type;

      //! Typename of scalar's (which may be different from value_type)
      typedef typename ExprType::scalar_type scalar_type;

      //! Typename of scalar's (which may be different from ValueT)
      typedef typename ScalarType<ValueT>::type ScalarT;

      //! Turn DFad into a meta-function class usable with mpl::apply
      template <typename T>
      struct apply {
        typedef DFad<T> type;
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
      DFad() :
        ExprType() {}

      //! Constructor with supplied value \c x convertible to ValueT
      /*!
       * Initializes value to \c ValueT(x) and derivative array is empty.
       */
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      DFad(const S& x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        ExprType(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      KOKKOS_INLINE_FUNCTION
      DFad(const int sz, const ValueT& x, const DerivInit zero_out = InitDerivArray) :
        ExprType(sz,x,zero_out) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      KOKKOS_INLINE_FUNCTION
      DFad(const int sz, const int i, const ValueT & x) :
        ExprType(sz,i,x) {}

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      DFad(const DFad& x) :
        ExprType(x) {}

      //! Copy constructor from any Expression object
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      DFad(const Expr<S>& x, SACADO_ENABLE_EXPR_CTOR_DECL) :
        ExprType(x) {}

      //@}

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~DFad() {}

      //! Assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(DFad&) operator=(const S& v) {
        GeneralFadType::operator=(v);
        return *this;
      }

      //! Assignment operator with DFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      DFad& operator=(const DFad& x) {
        GeneralFadType::operator=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(DFad&) operator=(const Expr<S>& x)
      {
        GeneralFadType::operator=(x);
        return *this;
      }

      //! Addition-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(DFad&) operator += (const S& x) {
        GeneralFadType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(DFad&) operator -= (const S& x) {
        GeneralFadType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(DFad&) operator *= (const S& x) {
        GeneralFadType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(DFad&) operator /= (const S& x) {
        GeneralFadType::operator/=(x);
        return *this;
      }

      //! Addition-assignment operator with DFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      DFad& operator += (const DFad& x) {
        GeneralFadType::operator+=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Subtraction-assignment operator with DFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      DFad& operator -= (const DFad& x) {
        GeneralFadType::operator-=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Multiplication-assignment operator with DFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      DFad& operator *= (const DFad& x) {
        GeneralFadType::operator*=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Division-assignment operator with DFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      DFad& operator /= (const DFad& x) {
        GeneralFadType::operator/=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(DFad&) operator += (const Expr<S>& x) {
        GeneralFadType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(DFad&) operator -= (const Expr<S>& x) {
        GeneralFadType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(DFad&) operator *= (const Expr<S>& x) {
        GeneralFadType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(DFad&) operator /= (const Expr<S>& x) {
        GeneralFadType::operator/=(x);
        return *this;
      }

    }; // class DFad<ValueT>

    template <typename T>
    struct BaseExpr< GeneralFad<T,Fad::DynamicStorage<T> > > {
      typedef DFad< typename GeneralFad<T,Fad::DynamicStorage<T> >::value_type > type;
    };

    template <typename T>
    struct ExprLevel< DFad<T> > {
      static const unsigned value =
        ExprLevel< typename DFad<T>::value_type >::value + 1;
    };

    template <typename T>
    struct IsFadExpr< DFad<T> > {
      static const bool value = true;
    };

  } // namespace Fad

  template <typename T>
  struct IsFad< FAD_NS::DFad<T> > {
    static const bool value = true;
  };

  template <typename T>
  struct IsExpr< FAD_NS::DFad<T> > {
    static const bool value = true;
  };

  template <typename T>
  struct BaseExprType< FAD_NS::DFad<T> > {
    typedef typename FAD_NS::DFad<T>::base_expr_type type;
  };

  template <typename,unsigned,unsigned> struct ViewFadType;
  namespace FAD_NS {
    template <typename,unsigned,unsigned,typename> class ViewFad;
  }

  //! The View Fad type associated with this type
  template< class ValueType, unsigned length, unsigned stride >
  struct ViewFadType< Sacado::FAD_NS::DFad< ValueType >, length, stride > {
    typedef Sacado::FAD_NS::ViewFad< ValueType,length,stride,Sacado::FAD_NS::DFad<ValueType> > type;
  };

  //! The View Fad type associated with this type
  /*!
   * Do not include the const in the base expr type.
   */
  template< class ValueType, unsigned length, unsigned stride >
  struct ViewFadType< const Sacado::FAD_NS::DFad< ValueType >, length, stride > {
    typedef Sacado::FAD_NS::ViewFad< const ValueType,length,stride,Sacado::FAD_NS::DFad<ValueType> > type;
  };

} // namespace Sacado

//
// Classes needed for Kokkos::View< DFad<...> ... > specializations
//
// Users can disable these view specializations either at configure time or
// by defining SACADO_DISABLE_FAD_VIEW_SPEC in their code.
//

#if defined(HAVE_SACADO_KOKKOSCORE) && defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#include "impl/Kokkos_AnalyzeShape.hpp"
#include "Kokkos_AnalyzeSacadoShape.hpp"

// Forward declarations
namespace Sacado {
  namespace FAD_NS {
    template <typename,unsigned,unsigned,typename> class ViewFad;
  }
}

namespace Kokkos {
namespace Impl {

// Forward declarations
struct ViewSpecializeSacadoFad;

/** \brief  Analyze the array shape of a Sacado::FAD_NS::DFad<T>.
 *
 *  This specialization is required so that the array shape of
 *  Kokkos::View< Sacado::FAD_NS::DFad<T>, ... >
 *  can be determined at compile-time.
 *
 *  For View purposes, DFad is treated as a dynamic dimension.
 */
template< class ValueType >
struct AnalyzeShape< Sacado::FAD_NS::DFad< ValueType > >
  : Shape< sizeof(Sacado::FAD_NS::DFad< ValueType >) , 0 > // Treat as a scalar
{
public:

  typedef ViewSpecializeSacadoFad specialize ;

  typedef Shape< sizeof(Sacado::FAD_NS::DFad< ValueType >) , 0 > shape ;

  typedef       Sacado::FAD_NS::DFad< ValueType >  array_intrinsic_type ;
  typedef const Sacado::FAD_NS::DFad< ValueType >  const_array_intrinsic_type ;
  typedef array_intrinsic_type non_const_array_intrinsic_type ;

  typedef       Sacado::FAD_NS::DFad< ValueType >  type ;
  typedef const Sacado::FAD_NS::DFad< ValueType >  const_type ;
  typedef       Sacado::FAD_NS::DFad< ValueType >  non_const_type ;

  typedef       Sacado::FAD_NS::DFad< ValueType >  value_type ;
  typedef const Sacado::FAD_NS::DFad< ValueType >  const_value_type ;
  typedef       Sacado::FAD_NS::DFad< ValueType >  non_const_value_type ;
};

/** \brief  Analyze the array shape of a Sacado::FAD_NS::DFad<T>.
 *
 *  This specialization is required so that the array shape of
 *  Kokkos::View< Sacado::FAD_NS::DFad<T>, ... >
 *  can be determined at compile-time.
 *
 *  For View purposes, DFad is treated as a dynamic dimension.
 */
template< class ValueType, class Layout >
struct AnalyzeSacadoShape< Sacado::FAD_NS::DFad< ValueType >, Layout >
  : ShapeInsert< typename AnalyzeSacadoShape< ValueType, Layout >::shape , 0 >::type
{
private:

  typedef AnalyzeSacadoShape< ValueType, Layout > nested ;

public:

  typedef ViewSpecializeSacadoFad specialize ;

  typedef typename ShapeInsert< typename nested::shape , 0 >::type shape ;

  typedef typename nested::array_intrinsic_type *        array_intrinsic_type ;
  typedef typename nested::const_array_intrinsic_type *  const_array_intrinsic_type ;
  typedef array_intrinsic_type non_const_array_intrinsic_type ;

  typedef       Sacado::FAD_NS::DFad< ValueType >  type ;
  typedef const Sacado::FAD_NS::DFad< ValueType >  const_type ;
  typedef       Sacado::FAD_NS::DFad< ValueType >  non_const_type ;

  typedef       Sacado::FAD_NS::DFad< ValueType >  value_type ;
  typedef const Sacado::FAD_NS::DFad< ValueType >  const_value_type ;
  typedef       Sacado::FAD_NS::DFad< ValueType >  non_const_value_type ;

  typedef typename nested::type           flat_array_type ;
  typedef typename nested::const_type     const_flat_array_type ;
  typedef typename nested::non_const_type non_const_flat_array_type ;
};

} // namespace Impl
} // namespace Kokkos

#endif
