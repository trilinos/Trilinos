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

    // Forward declaration
    template <typename T, int Num>
    class StaticStorage;

    /*!
     * \brief Forward-mode AD class using static memory allocation
     * with long arrays and expression templates.
     */
    /*!
     * This is the user-level class for forward mode AD with static
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is known at compile time.  The largest size
     * of the derivative array is fixed by the template parameter \c Num
     * while the actual size used is set by the \c sz argument to the
     * constructor or the \c n argument to diff().  The user
     * interface is provided by Sacado::FAD_NS::GeneralFad.
     */
    template <typename ValueT, int Num>
    class SLFad :
      public Expr< GeneralFad<ValueT,Fad::StaticStorage<ValueT,Num> > > {

    public:

      //! Base classes
      typedef Fad::StaticStorage<ValueT,Num> StorageType;
      typedef GeneralFad<ValueT,StorageType> GeneralFadType;
      typedef Expr<GeneralFadType> ExprType;

      //! Typename of values
      typedef typename ExprType::value_type value_type;

      //! Typename of scalar's (which may be different from value_type)
      typedef typename ExprType::scalar_type scalar_type;

      //! Typename of scalar's (which may be different from ValueT)
      typedef typename ScalarType<ValueT>::type ScalarT;

      //! Turn SLFad into a meta-function class usable with mpl::apply
      template <typename T>
      struct apply {
        typedef SLFad<T,Num> type;
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
      SLFad() :
        ExprType() {}

      //! Constructor with supplied value \c x convertible to ValueT
      /*!
       * Initializes value to \c ValueT(x) and derivative array is empty.
       */
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SLFad(const S& x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        ExprType(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      KOKKOS_INLINE_FUNCTION
      SLFad(const int sz, const ValueT & x, const bool zero_out = true) :
        ExprType(sz,x,zero_out) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      KOKKOS_INLINE_FUNCTION
      SLFad(const int sz, const int i, const ValueT & x) :
        ExprType(sz,i,x) {}

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      SLFad(const SLFad& x) :
        ExprType(x) {}

      //! Copy constructor from any Expression object
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SLFad(const Expr<S>& x, SACADO_ENABLE_EXPR_CTOR_DECL) :
        ExprType(x) {}

      //@}

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~SLFad() {}

      //! Assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SLFad&) operator=(const S& v) {
        GeneralFadType::operator=(v);
        return *this;
      }

      //! Assignment operator with SLFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      SLFad& operator=(const SLFad& x) {
        GeneralFadType::operator=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(SLFad&) operator=(const Expr<S>& x)
      {
        GeneralFadType::operator=(x);
        return *this;
      }

      //! Addition-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SLFad&) operator += (const S& x) {
        GeneralFadType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SLFad&) operator -= (const S& x) {
        GeneralFadType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SLFad&) operator *= (const S& x) {
        GeneralFadType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SLFad&) operator /= (const S& x) {
        GeneralFadType::operator/=(x);
        return *this;
      }

      //! Addition-assignment operator with SLFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      SLFad& operator += (const SLFad& x) {
        GeneralFadType::operator+=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Subtraction-assignment operator with SLFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      SLFad& operator -= (const SLFad& x) {
        GeneralFadType::operator-=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Multiplication-assignment operator with SLFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      SLFad& operator *= (const SLFad& x) {
        GeneralFadType::operator*=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Division-assignment operator with SLFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      SLFad& operator /= (const SLFad& x) {
        GeneralFadType::operator/=(static_cast<const GeneralFadType&>(x));
        return *this;
      }

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(SLFad&) operator += (const Expr<S>& x) {
        GeneralFadType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(SLFad&) operator -= (const Expr<S>& x) {
        GeneralFadType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(SLFad&) operator *= (const Expr<S>& x) {
        GeneralFadType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(SLFad&) operator /= (const Expr<S>& x) {
        GeneralFadType::operator/=(x);
        return *this;
      }

    }; // class SLFad<ValueT,Num>

    template <typename T, int N>
    struct BaseExpr< GeneralFad<T,Fad::StaticStorage<T,N> > > {
      typedef SLFad<typename GeneralFad<T,Fad::StaticStorage<T,N> >::value_type,N> type;
    };

    template <typename T, int N>
    struct ExprLevel< SLFad<T,N> > {
      static const unsigned value =
        ExprLevel< typename SLFad<T,N>::value_type >::value + 1;
    };

    template <typename T, int N>
    struct IsFadExpr< SLFad<T,N> > {
      static const bool value = true;
    };

  } // namespace Fad

  template <typename T, int N>
  struct IsFad< FAD_NS::SLFad<T,N> > {
    static const bool value = true;
  };

  template <typename T, int N>
  struct IsExpr< FAD_NS::SLFad<T,N> > {
    static const bool value = true;
  };

  template <typename T, int N>
  struct BaseExprType< FAD_NS::SLFad<T,N> > {
    typedef typename FAD_NS::SLFad<T,N>::base_expr_type type;
  };

  template <typename T,unsigned,unsigned> struct ViewFadType;
  namespace FAD_NS {
    template <typename,unsigned,unsigned,typename> class ViewFad;
  }

  //! The View Fad type associated with this type
  template< class ValueType, int N, unsigned length, unsigned stride >
  struct ViewFadType< Sacado::FAD_NS::SLFad< ValueType, N >, length, stride > {
    typedef Sacado::FAD_NS::ViewFad< ValueType,length,stride,Sacado::FAD_NS::SLFad<ValueType,N> > type;
  };

  //! The View Fad type associated with this type
  /*!
   * Do not include the const in the base expr type.
   */
  template< class ValueType, int N, unsigned length, unsigned stride >
  struct ViewFadType< const Sacado::FAD_NS::SLFad< ValueType, N >, length, stride > {
    typedef Sacado::FAD_NS::ViewFad< const ValueType,length,stride,Sacado::FAD_NS::SLFad<ValueType,N > > type;
  };

} // namespace Sacado

//
// Classes needed for Kokkos::View< SLFad<...> ... > specializations
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

/** \brief  Analyze the array shape of a Sacado::FAD_NS::SLFad<T,N>.
 *
 *  This specialization is required so that the array shape of
 *  Kokkos::View< Sacado::FAD_NS::SLFad<T,N>, ... >
 *  can be determined at compile-time.
 *
 *  For View purposes, SLFad is treated as a dynamic dimension.
 */
template< class ValueType, int N >
struct AnalyzeShape< Sacado::FAD_NS::SLFad< ValueType, N > >
  : Shape< sizeof(Sacado::FAD_NS::SLFad< ValueType, N >) , 0 > // Treat as a scalar
{
public:

  typedef ViewSpecializeSacadoFad specialize ;

  typedef Shape< sizeof(Sacado::FAD_NS::SLFad< ValueType, N >) , 0 > shape ;

  typedef       Sacado::FAD_NS::SLFad< ValueType, N >  array_intrinsic_type ;
  typedef const Sacado::FAD_NS::SLFad< ValueType, N >  const_array_intrinsic_type ;
  typedef array_intrinsic_type non_const_array_intrinsic_type ;

  typedef       Sacado::FAD_NS::SLFad< ValueType, N >  type ;
  typedef const Sacado::FAD_NS::SLFad< ValueType, N >  const_type ;
  typedef       Sacado::FAD_NS::SLFad< ValueType, N >  non_const_type ;

  typedef       Sacado::FAD_NS::SLFad< ValueType, N >  value_type ;
  typedef const Sacado::FAD_NS::SLFad< ValueType, N >  const_value_type ;
  typedef       Sacado::FAD_NS::SLFad< ValueType, N >  non_const_value_type ;
};

/** \brief  Analyze the array shape of a Sacado::FAD_NS::SLFad<T,N>.
 *
 *  This specialization is required so that the array shape of
 *  Kokkos::View< Sacado::FAD_NS::SLFad<T,N>, ... >
 *  can be determined at compile-time.
 *
 *  For View purposes, SLFad is treated as a dynamic dimension.
 */
template< class ValueType, class Layout, int N >
struct AnalyzeSacadoShape< Sacado::FAD_NS::SLFad< ValueType, N >, Layout >
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

  typedef       Sacado::FAD_NS::SLFad< ValueType, N >  type ;
  typedef const Sacado::FAD_NS::SLFad< ValueType, N >  const_type ;
  typedef       Sacado::FAD_NS::SLFad< ValueType, N >  non_const_type ;

  typedef       Sacado::FAD_NS::SLFad< ValueType, N >  value_type ;
  typedef const Sacado::FAD_NS::SLFad< ValueType, N >  const_value_type ;
  typedef       Sacado::FAD_NS::SLFad< ValueType, N >  non_const_value_type ;

  typedef typename nested::type           flat_array_type ;
  typedef typename nested::const_type     const_flat_array_type ;
  typedef typename nested::non_const_type non_const_flat_array_type ;
};

} // namespace Impl
} // namespace Kokkos

#endif
