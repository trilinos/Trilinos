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

#ifndef SACADO_FAD_EXP_VIEWFAD_HPP
#define SACADO_FAD_EXP_VIEWFAD_HPP

#include "Sacado_Fad_Exp_GeneralFad.hpp"
#include "Sacado_Fad_Exp_ViewStorage.hpp"

namespace Sacado {

  namespace Fad {
  namespace Exp {

    template <typename T, unsigned static_length, unsigned static_stride, typename U>
    using  ViewFad = GeneralFad< ViewStorage<T,static_length,static_stride,U> >;

  } // namespace Exp
  } // namespace Fad

  template <typename,unsigned,unsigned> struct ViewFadType;

  //! The View Fad type associated with this type
  template< class S, unsigned length, unsigned stride >
  struct ViewFadType< Fad::Exp::GeneralFad<S>, length, stride > {
    typedef Fad::Exp::ViewFad< typename S::value_type,length,stride,Fad::Exp::GeneralFad<S> > type;
  };

  //! The View Fad type associated with this type
  /*!
   * Do not include the const in the base expr type.
   */
  template< class S, unsigned length, unsigned stride >
  struct ViewFadType< const Fad::Exp::GeneralFad<S>, length, stride > {
    typedef Fad::Exp::ViewFad< const typename S::value_type,length,stride,Fad::Exp::GeneralFad<S> > type;
  };

  // Specialization of BaseExprType for ViewFad, to use the base fad type
  template <typename T, unsigned static_length, unsigned static_stride, typename U>
  struct BaseExprType< Fad::Exp::GeneralFad< Fad::Exp::ViewStorage<T,static_length,static_stride,U> > > {
    typedef U type;
  };

} // namespace Sacado

#endif // SACADO_FAD_EXP_VIEWFAD_HPP
