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

    //! Forward-mode AD class using static memory allocation
    /*!
     * This is the user-level class for forward mode AD with static
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is known at compile time.  The size
     * of the derivative array is fixed by the template parameter \c Num.
     */
    template <typename ValueT, int Num>
    class SFad :
      public Expr< SFadExprTag<ValueT,Num > > {

    public:

      //! Base classes
      typedef Expr< SFadExprTag<ValueT,Num > > ExprType;

      //! Typename of values
      typedef typename ExprType::value_type value_type;

      //! Typename of scalar's (which may be different from value_type)
      typedef typename ExprType::scalar_type scalar_type;

      //! Typename of scalar's (which may be different from ValueT)
      typedef typename ScalarType<ValueT>::type ScalarT;

      //! Turn SFad into a meta-function class usable with mpl::apply
      template <typename T>
      struct apply {
        typedef SFad<T,Num> type;
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
      SFad() :
        ExprType() {}

      //! Constructor with supplied value \c x convertible to ValueT
      /*!
       * Initializes value to \c ValueT(x) and derivative array is empty.
       */
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SFad(const S& x, SACADO_ENABLE_VALUE_CTOR_DECL) :
        ExprType(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      KOKKOS_INLINE_FUNCTION
      SFad(const int sz, const ValueT & x) :
        ExprType(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      KOKKOS_INLINE_FUNCTION
      SFad(const int sz, const int i, const ValueT & x) :
        ExprType(sz,i,x) {}

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      SFad(const SFad& x) :
        ExprType(static_cast<const ExprType&>(x)) {}

      //! Copy constructor from any Expression object
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SFad(const Expr<S>& x, SACADO_ENABLE_EXPR_CTOR_DECL) :
        ExprType(x) {}

      //@}

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~SFad() {}

      //! Assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator=(const S& v) {
        ExprType::operator=(v);
        return *this;
      }

      //! Assignment operator with SFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      SFad& operator=(const SFad& x) {
        ExprType::operator=(static_cast<const ExprType&>(x));
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_EXPR_FUNC(SFad&) operator=(const Expr<S>& x)
      {
        ExprType::operator=(x);
        return *this;
      }

      //! Addition-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator += (const S& x) {
        ExprType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator -= (const S& x) {
        ExprType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator *= (const S& x) {
        ExprType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with constant right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator /= (const S& x) {
        ExprType::operator/=(x);
        return *this;
      }

      //! Addition-assignment operator with SFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      SFad& operator += (const SFad& x) {
        ExprType::operator+=(static_cast<const ExprType&>(x));
        return *this;
      }

      //! Subtraction-assignment operator with SFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      SFad& operator -= (const SFad& x) {
        ExprType::operator-=(static_cast<const ExprType&>(x));
        return *this;
      }

      //! Multiplication-assignment operator with SFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      SFad& operator *= (const SFad& x) {
        ExprType::operator*=(static_cast<const ExprType&>(x));
        return *this;
      }

      //! Division-assignment operator with SFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      SFad& operator /= (const SFad& x) {
        ExprType::operator/=(static_cast<const ExprType&>(x));
        return *this;
      }

       //! Addition-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator += (const Expr<S>& x) {
        ExprType::operator+=(x);
        return *this;
      }

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator -= (const Expr<S>& x) {
        ExprType::operator-=(x);
        return *this;
      }

      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator *= (const Expr<S>& x) {
        ExprType::operator*=(x);
        return *this;
      }

      //! Division-assignment operator with Expr right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SACADO_ENABLE_VALUE_FUNC(SFad&) operator /= (const Expr<S>& x) {
        ExprType::operator/=(x);
        return *this;
      }

    }; // class SFad<ValueT,Num>

    template <typename T, int Num>
    std::ostream& operator << (std::ostream& os,
                               const Expr< SFadExprTag<T,Num> >& x) {
      os << x.val() << " [";

      for (int i=0; i< x.size(); i++) {
        os << " " << x.dx(i);
      }

      os << " ]";
      return os;
    }

    template <typename T, int N>
    struct ExprLevel< SFad<T,N> > {
      static const unsigned value =
        ExprLevel< typename SFad<T,N>::value_type >::value + 1;
    };

  } // namespace Fad

  template <typename T,unsigned,unsigned> struct ViewFadType;
  namespace FAD_NS {
    template <typename,unsigned,unsigned,typename> class ViewFad;
  }

  //! The View Fad type associated with this type
  template< class ValueType, int N, unsigned length, unsigned stride >
  struct ViewFadType< Sacado::FAD_NS::SFad< ValueType, N >, length, stride > {
    typedef Sacado::FAD_NS::ViewFad< ValueType,length,stride,Sacado::FAD_NS::SFad< ValueType,N> > type;
};

   //! The View Fad type associated with this type
  template< class ValueType, int N, unsigned length, unsigned stride >
  struct ViewFadType< const Sacado::FAD_NS::SFad< ValueType, N >, length, stride > {
    typedef Sacado::FAD_NS::ViewFad< const ValueType,length,stride,Sacado::FAD_NS::SFad<const ValueType,N> > type;
  };

} // namespace Sacado

//
// Classes needed for Kokkos::View< SFad<...> ... > specializations
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

/** \brief  Analyze the array shape of a Sacado::FAD_NS::SFad<T,N>.
 *
 *  This specialization is required so that the array shape of
 *  Kokkos::View< Sacado::FAD_NS::SFad<T,N>, ... >
 *  can be determined at compile-time.
 *
 *  We add one to the SFad dimension (N) to store the value component.
 */
template< class ValueType, int N >
struct AnalyzeShape< Sacado::FAD_NS::SFad< ValueType, N > >
  : Shape< sizeof(Sacado::FAD_NS::SFad< ValueType, N >) , 0 > // Treat as a scalar
{
public:

  typedef ViewSpecializeSacadoFad specialize ;

  typedef Shape< sizeof(Sacado::FAD_NS::SFad< ValueType, N >) , 0 > shape ;

  typedef       Sacado::FAD_NS::SFad< ValueType, N >  array_intrinsic_type ;
  typedef const Sacado::FAD_NS::SFad< ValueType, N >  const_array_intrinsic_type ;
  typedef array_intrinsic_type non_const_array_intrinsic_type ;

  typedef       Sacado::FAD_NS::SFad< ValueType, N >  type ;
  typedef const Sacado::FAD_NS::SFad< ValueType, N >  const_type ;
  typedef       Sacado::FAD_NS::SFad< ValueType, N >  non_const_type ;

  typedef       Sacado::FAD_NS::SFad< ValueType, N >  value_type ;
  typedef const Sacado::FAD_NS::SFad< ValueType, N >  const_value_type ;
  typedef       Sacado::FAD_NS::SFad< ValueType, N >  non_const_value_type ;
};

/** \brief  Analyze the array shape of a Sacado::FAD_NS::SFad<T,N>.
 *
 *  This specialization is required so that the array shape of
 *  Kokkos::View< Sacado::FAD_NS::SFad<T,N>, ... >
 *  can be determined at compile-time.
 *
 *  We add one to the SFad dimension (N) to store the value component.
 */
template< class ValueType, class Layout, int N >
struct AnalyzeSacadoShape< Sacado::FAD_NS::SFad< ValueType, N >, Layout >
  : ShapeInsert< typename AnalyzeSacadoShape< ValueType, Layout >::shape , N+1 >::type
{
private:

  typedef AnalyzeSacadoShape< ValueType, Layout > nested ;

public:

  typedef ViewSpecializeSacadoFad specialize ;

  typedef typename ShapeInsert< typename nested::shape , N+1 >::type shape ;

  typedef typename nested::array_intrinsic_type         array_intrinsic_type [N+1];
  typedef typename nested::const_array_intrinsic_type   const_array_intrinsic_type [N+1] ;
  typedef array_intrinsic_type non_const_array_intrinsic_type ;

  typedef       Sacado::FAD_NS::SFad< ValueType, N >  type ;
  typedef const Sacado::FAD_NS::SFad< ValueType, N >  const_type ;
  typedef       Sacado::FAD_NS::SFad< ValueType, N >  non_const_type ;

  typedef       Sacado::FAD_NS::SFad< ValueType, N >  value_type ;
  typedef const Sacado::FAD_NS::SFad< ValueType, N >  const_value_type ;
  typedef       Sacado::FAD_NS::SFad< ValueType, N >  non_const_value_type ;

  typedef typename nested::type           flat_array_type ;
  typedef typename nested::const_type     const_flat_array_type ;
  typedef typename nested::non_const_type non_const_flat_array_type ;
};

} // namespace Impl
} // namespace Kokkos

#endif
