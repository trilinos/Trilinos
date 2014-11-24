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

#ifndef SACADO_ELRFAD_DFAD_HPP
#define SACADO_ELRFAD_DFAD_HPP

#include "Sacado_ELRFad_GeneralFadExpr.hpp"
#include "Sacado_ELRFad_DFadTraits.hpp"
#include "Sacado_Fad_DynamicStorage.hpp"

namespace Sacado {

  namespace ELRFad {

    /*!
     * \brief Forward-mode AD class using dynamic memory allocation and
     * expression-level reverse mode expression templates.
     */
    /*!
     * This is the user-level class for forward mode AD with dynamic
     * memory allocation, and is appropriate for whenever the number
     * of derivative components is not known at compile time.  The user
     * interface is provided by Sacado::ELRFad::GeneralFad.
     */
#include "Sacado_Fad_DFad_tmpl.hpp"

  } // namespace ELRFad

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
#include "Sacado_ELRFad_ViewFad.hpp"

namespace Kokkos {
namespace Impl {

// Forward declarations
struct ViewSpecializeSacadoFad;
template <typename T,unsigned,unsigned> struct ViewFadType;

//! The View Fad type associated with this type
template< class ValueType, unsigned length, unsigned stride >
struct ViewFadType< Sacado::ELRFad::DFad< ValueType >, length, stride > {
  typedef Sacado::ELRFad::ViewFad<ValueType,length,stride> type;
};

//! The View Fad type associated with this type
template< class ValueType, unsigned length, unsigned stride >
struct ViewFadType< const Sacado::ELRFad::DFad< ValueType >, length, stride > {
  typedef Sacado::ELRFad::ViewFad<const ValueType,length,stride> type;
};

/** \brief  Analyze the array shape of a Sacado::ELRFad::DFad<T>.
 *
 *  This specialization is required so that the array shape of
 *  Kokkos::View< Sacado::ELRFad::DFad<T>, ... >
 *  can be determined at compile-time.
 *
 *  For View purposes, DFad is treated as a dynamic dimension.
 */
template< class ValueType >
struct AnalyzeShape< Sacado::ELRFad::DFad< ValueType > >
  : Shape< sizeof(Sacado::ELRFad::DFad< ValueType >) , 0 > // Treat as a scalar
{
public:

  typedef ViewSpecializeSacadoFad specialize ;

  typedef Shape< sizeof(Sacado::ELRFad::DFad< ValueType >) , 0 > shape ;

  typedef       Sacado::ELRFad::DFad< ValueType >        array_intrinsic_type ;
  typedef const Sacado::ELRFad::DFad< ValueType >  const_array_intrinsic_type ;
  typedef array_intrinsic_type non_const_array_intrinsic_type ;

  typedef       Sacado::ELRFad::DFad< ValueType >  type ;
  typedef const Sacado::ELRFad::DFad< ValueType >  const_type ;
  typedef       Sacado::ELRFad::DFad< ValueType >  non_const_type ;

  typedef       Sacado::ELRFad::DFad< ValueType >  value_type ;
  typedef const Sacado::ELRFad::DFad< ValueType >  const_value_type ;
  typedef       Sacado::ELRFad::DFad< ValueType >  non_const_value_type ;
};

/** \brief  Analyze the array shape of a Sacado::ELRFad::DFad<T>.
 *
 *  This specialization is required so that the array shape of
 *  Kokkos::View< Sacado::ELRFad::DFad<T>, ... >
 *  can be determined at compile-time.
 *
 *  For View purposes, DFad is treated as a dynamic dimension.
 */
template< class ValueType, class Layout >
struct AnalyzeSacadoShape< Sacado::ELRFad::DFad< ValueType >, Layout >
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

  typedef       Sacado::ELRFad::DFad< ValueType >  type ;
  typedef const Sacado::ELRFad::DFad< ValueType >  const_type ;
  typedef       Sacado::ELRFad::DFad< ValueType >  non_const_type ;

  typedef       Sacado::ELRFad::DFad< ValueType >  value_type ;
  typedef const Sacado::ELRFad::DFad< ValueType >  const_value_type ;
  typedef       Sacado::ELRFad::DFad< ValueType >  non_const_value_type ;

  typedef typename nested::type           flat_array_type ;
  typedef typename nested::const_type     const_flat_array_type ;
  typedef typename nested::non_const_type non_const_flat_array_type ;
};

} // namespace Impl
} // namespace Kokkos

#endif

#endif // SACADO_ELRFAD_DFAD_HPP
