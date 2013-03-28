// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_SACADO_HPP
#define STOKHOS_SACADO_HPP

// We need to define math functions first for nested AD types
#include "Sacado_MathFunctions.hpp"
#include "Stokhos_Sacado_MathFunctions.hpp"

// Stokhos headers
#include "Stokhos.hpp"

// Traits classes
#include "Sacado_ETPCE_ExpressionTraits.hpp"
#include "Sacado_ETPCE_OrthogPolyTraits.hpp"
#include "Sacado_ETV_ExpressionTraits.hpp"
#include "Sacado_ETV_VectorTraits.hpp"
#include "Sacado_ETV_Vector2Traits.hpp"

// Sacado overloaded operators for Stokhos
#include "Stokhos_StandardStorage.hpp"
#include "Stokhos_StaticStandardStorage.hpp"
#include "Stokhos_StaticFixedStandardStorage.hpp"
#include "Sacado_PCE_OrthogPoly.hpp"
#include "Sacado_ETPCE_OrthogPoly.hpp"
#include "Sacado_ETV_Vector.hpp"
#include "Sacado_ETV_Vector2.hpp"

#endif // STOKHOS_SACADO_HPP
