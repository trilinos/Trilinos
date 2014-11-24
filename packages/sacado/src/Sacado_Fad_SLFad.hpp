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

#ifndef SACADO_FAD_SLFAD_HPP
#define SACADO_FAD_SLFAD_HPP

#include "Sacado_Fad_GeneralFadExpr.hpp"
#include "Sacado_Fad_SLFadTraits.hpp"
#include "Sacado_Fad_StaticStorage.hpp"

namespace Sacado {

  namespace Fad {

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
     * interface is provided by Sacado::Fad::GeneralFad.
     */
#include "Sacado_Fad_SLFad_tmpl.hpp"

  } // namespace Fad

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
#include "Sacado_Fad_ViewFad.hpp"

namespace Kokkos {
namespace Impl {

// Forward declarations
struct ViewSpecializeSacadoFad;
template <typename T,unsigned,unsigned> struct ViewFadType;

//! The View Fad type associated with this type
template< class ValueType, int N, unsigned length, unsigned stride >
struct ViewFadType< Sacado::Fad::SLFad< ValueType, N >, length, stride > {
  typedef Sacado::Fad::ViewFad<ValueType,length,stride> type;
};

//! The View Fad type associated with this type
template< class ValueType, int N, unsigned length, unsigned stride >
struct ViewFadType< const Sacado::Fad::SLFad< ValueType, N >, length, stride > {
  typedef Sacado::Fad::ViewFad<const ValueType,length,stride> type;
};

/** \brief  Analyze the array shape of a Sacado::Fad::SLFad<T,N>.
 *
 *  This specialization is required so that the array shape of
 *  Kokkos::View< Sacado::Fad::SLFad<T,N>, ... >
 *  can be determined at compile-time.
 *
 *  For View purposes, SLFad is treated as a dynamic dimension.
 */
template< class ValueType, int N >
struct AnalyzeShape< Sacado::Fad::SLFad< ValueType, N > >
  : Shape< sizeof(Sacado::Fad::SLFad< ValueType, N >) , 0 > // Treat as a scalar
{
public:

  typedef ViewSpecializeSacadoFad specialize ;

  typedef Shape< sizeof(Sacado::Fad::SLFad< ValueType, N >) , 0 > shape ;

  typedef       Sacado::Fad::SLFad< ValueType, N >        array_intrinsic_type ;
  typedef const Sacado::Fad::SLFad< ValueType, N >  const_array_intrinsic_type ;
  typedef array_intrinsic_type non_const_array_intrinsic_type ;

  typedef       Sacado::Fad::SLFad< ValueType, N >  type ;
  typedef const Sacado::Fad::SLFad< ValueType, N >  const_type ;
  typedef       Sacado::Fad::SLFad< ValueType, N >  non_const_type ;

  typedef       Sacado::Fad::SLFad< ValueType, N >  value_type ;
  typedef const Sacado::Fad::SLFad< ValueType, N >  const_value_type ;
  typedef       Sacado::Fad::SLFad< ValueType, N >  non_const_value_type ;
};

/** \brief  Analyze the array shape of a Sacado::Fad::SLFad<T,N>.
 *
 *  This specialization is required so that the array shape of
 *  Kokkos::View< Sacado::Fad::SLFad<T,N>, ... >
 *  can be determined at compile-time.
 *
 *  For View purposes, SLFad is treated as a dynamic dimension.
 */
template< class ValueType, class Layout, int N >
struct AnalyzeSacadoShape< Sacado::Fad::SLFad< ValueType, N >, Layout >
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

  typedef       Sacado::Fad::SLFad< ValueType, N >  type ;
  typedef const Sacado::Fad::SLFad< ValueType, N >  const_type ;
  typedef       Sacado::Fad::SLFad< ValueType, N >  non_const_type ;

  typedef       Sacado::Fad::SLFad< ValueType, N >  value_type ;
  typedef const Sacado::Fad::SLFad< ValueType, N >  const_value_type ;
  typedef       Sacado::Fad::SLFad< ValueType, N >  non_const_value_type ;

  typedef typename nested::type           flat_array_type ;
  typedef typename nested::const_type     const_flat_array_type ;
  typedef typename nested::non_const_type non_const_flat_array_type ;
};

} // namespace Impl
} // namespace Kokkos

#endif

#endif // SACADO_FAD_SLFAD_HPP
