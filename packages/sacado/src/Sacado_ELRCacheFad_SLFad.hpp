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

#ifndef SACADO_ELRCACHEFAD_SLFAD_HPP
#define SACADO_ELRCACHEFAD_SLFAD_HPP

#include "Sacado_ELRCacheFad_GeneralFadExpr.hpp"
#include "Sacado_ELRCacheFad_SLFadTraits.hpp"
#include "Sacado_Fad_StaticStorage.hpp"
#include "Sacado_dummy_arg.hpp"

namespace Sacado {

  namespace ELRCacheFad {

    /*!
     * \brief Forward-mode AD class using static memory allocation
     * with long arrays, caching expression templates, and expression
     * level reverse mode.
     */
    /*!
     * This is the user-level class for forward mode AD with static
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is known at compile time.  The largest size
     * of the derivative array is fixed by the template parameter \c Num
     * while the actual size used is set by the \c sz argument to the
     * constructor or the \c n argument to diff().  The user
     * interface is provided by Sacado::Fad::GeneralFad.  It is similar
     * to Sacado::ELRFad::SLFad, except it uses the caching expression templates
     * that cache the results of val() calculations for later dx()
     * calculations.
     */
    template <typename ValueT, int Num>
    class SLFad :
      public Expr< GeneralFad<ValueT,Fad::StaticStorage<ValueT,Num> > > {

    public:

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
        Expr< GeneralFad< ValueT,Fad::StaticStorage<ValueT,Num> > >() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      KOKKOS_INLINE_FUNCTION
      SLFad(const ValueT& x) :
        Expr< GeneralFad< ValueT,Fad::StaticStorage<ValueT,Num> > >(x) {}

      //! Constructor with supplied value \c x of type ScalarT
      /*!
       * Initializes value to \c ValueT(x) and derivative array is empty.
       * Creates a dummy overload when ValueT and ScalarT are the same type.
       */
      KOKKOS_INLINE_FUNCTION
      SLFad(const typename dummy<ValueT,ScalarT>::type& x) :
        Expr< GeneralFad< ValueT,Fad::StaticStorage<ValueT,Num> > >(ValueT(x)) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      KOKKOS_INLINE_FUNCTION
      SLFad(const int sz, const ValueT & x) :
        Expr< GeneralFad< ValueT,Fad::StaticStorage<ValueT,Num> > >(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      KOKKOS_INLINE_FUNCTION
      SLFad(const int sz, const int i, const ValueT & x) :
        Expr< GeneralFad< ValueT,Fad::StaticStorage<ValueT,Num> > >(sz,i,x) {}

      //! Copy constructor
      KOKKOS_INLINE_FUNCTION
      SLFad(const SLFad& x) :
        Expr< GeneralFad< ValueT,Fad::StaticStorage<ValueT,Num> > >(x) {}

      //! Copy constructor from any Expression object
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SLFad(const Expr<S>& x) :
        Expr< GeneralFad< ValueT,Fad::StaticStorage<ValueT,Num> > >(x) {}

      //@}

      //! Destructor
      KOKKOS_INLINE_FUNCTION
      ~SLFad() {}

      //! Assignment operator with constant right-hand-side
      KOKKOS_INLINE_FUNCTION
      SLFad& operator=(const ValueT& v) {
        GeneralFad< ValueT,Fad::StaticStorage<ValueT,Num> >::operator=(v);
        return *this;
      }

      //! Assignment operator with constant right-hand-side
      /*!
       * Creates a dummy overload when ValueT and ScalarT are the same type.
       */
      KOKKOS_INLINE_FUNCTION
      SLFad& operator=(const typename dummy<ValueT,ScalarT>::type& v) {
        GeneralFad< ValueT,Fad::StaticStorage<ValueT,Num> >::operator=(ValueT(v));
        return *this;
      }

      //! Assignment operator with DFad right-hand-side
      KOKKOS_INLINE_FUNCTION
      SLFad& operator=(const SLFad& x) {
        GeneralFad< ValueT,Fad::StaticStorage<ValueT,Num> >::operator=(static_cast<const GeneralFad< ValueT,Fad::StaticStorage<ValueT,Num> >&>(x));
        return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S>
      KOKKOS_INLINE_FUNCTION
      SLFad& operator=(const Expr<S>& x)
      {
        GeneralFad< ValueT,Fad::StaticStorage<ValueT,Num> >::operator=(x);
        return *this;
      }

    }; // class SLFad<ValueT,Num>

  } // namespace ELRCacheFad

} // namespace Sacado

//
// Classes needed for Kokkos::View< SLFad<...> ... > specializations
//
// Users can disable these view specializations either at configure time or
// by defining SACADO_DISABLE_FAD_VIEW_SPEC in their code.
//

#if defined(HAVE_SACADO_KOKKOSCORE) && defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#include "impl/Kokkos_AnalyzeShape.hpp"
#include "Sacado_ELRCacheFad_ViewFad.hpp"

namespace Kokkos {
namespace Impl {

// Forward declarations
struct ViewSpecializeSacadoFad;
template <typename T,unsigned,unsigned> struct ViewFadType;

//! The View Fad type associated with this type
template< class ValueType, int N, unsigned length, unsigned stride >
struct ViewFadType< Sacado::ELRCacheFad::SLFad< ValueType, N >, length, stride > {
  typedef Sacado::ELRCacheFad::ViewFad<ValueType,length,stride> type;
};

//! The View Fad type associated with this type
template< class ValueType, int N, unsigned length, unsigned stride >
struct ViewFadType< const Sacado::ELRCacheFad::SLFad< ValueType, N >, length, stride > {
  typedef Sacado::ELRCacheFad::ViewFad<const ValueType,length,stride> type;
};

/** \brief  Analyze the array shape of a Sacado::ELRCacheFad::SLFad<T,N>.
 *
 *  This specialization is required so that the array shape of
 *  Kokkos::View< Sacado::ELRCacheFad::SLFad<T,N>, ... >
 *  can be determined at compile-time.
 *
 *  For View purposes, SLFad is treated as a dynamic dimension.
 */
template< class ValueType, int N >
struct AnalyzeShape< Sacado::ELRCacheFad::SLFad< ValueType, N > >
  : ShapeInsert< typename AnalyzeShape< ValueType >::shape , 0 >::type
{
private:

  typedef AnalyzeShape< ValueType > nested ;

public:

  typedef ViewSpecializeSacadoFad specialize ;

  typedef typename ShapeInsert< typename nested::shape , 0 >::type shape ;

  typedef typename nested::array_type *        array_type ;
  typedef typename nested::const_array_type *  const_array_type ;

  typedef array_type non_const_array_type ;

  typedef       Sacado::ELRCacheFad::SLFad< ValueType, N >  type ;
  typedef const Sacado::ELRCacheFad::SLFad< ValueType, N >  const_type ;
  typedef       Sacado::ELRCacheFad::SLFad< ValueType, N >  non_const_type ;

  typedef       Sacado::ELRCacheFad::SLFad< ValueType, N >  value_type ;
  typedef const Sacado::ELRCacheFad::SLFad< ValueType, N >  const_value_type ;
  typedef       Sacado::ELRCacheFad::SLFad< ValueType, N >  non_const_value_type ;
};

} // namespace Impl
} // namespace Kokkos

#endif

#endif // SACADO_ELRCACHEFAD_SLFAD_HPP
