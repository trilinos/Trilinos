// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef FEAPP_TEMPLATETYPES_HPP
#define FEAPP_TEMPLATETYPES_HPP

// Includ the MPL vector
#include "Sacado_mpl_vector.hpp"

// Include all of our AD types
#include "Sacado_Fad_DFad.hpp"

// Typedef AD types to standard names
typedef double RealType;
typedef Sacado::Fad::DFad<double> FadType;

// Define which types we are using
#define REAL_ACTIVE 1
#define FAD_ACTIVE 1

// Turn on/off explicit template instantiation
#define SACADO_ETI

// Build the MPL vector containing our valid types
typedef Sacado::mpl::vector<> ValidTypes0;
#if REAL_ACTIVE
typedef Sacado::mpl::push_back<ValidTypes0, RealType>::type ValidTypes1;
#else
typedef ValidTypes0 ValidTypes1;
#endif
#if FAD_ACTIVE
typedef Sacado::mpl::push_back<ValidTypes1, FadType>::type ValidTypes2;
#else
typedef ValidTypes1 ValidTypes2;
#endif
typedef ValidTypes2 ValidTypes;

// Define macro for explicit template instantiation
#if REAL_ACTIVE
#define INSTANTIATE_TEMPLATE_CLASS_REAL(name) template class name<double>;
#else
#define INSTANTIATE_TEMPLATE_CLASS_REAL(name)
#endif

#if FAD_ACTIVE
#define INSTANTIATE_TEMPLATE_CLASS_FAD(name) template class name<FadType>;
#else
#define INSTANTIATE_TEMPLATE_CLASS_FAD(name)
#endif

#define INSTANTIATE_TEMPLATE_CLASS(name) \
  INSTANTIATE_TEMPLATE_CLASS_REAL(name)	 \
  INSTANTIATE_TEMPLATE_CLASS_FAD(name)

#endif // FEAPP_TEMPLATETYPES_HPP
