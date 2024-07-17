// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

// Sacado overloaded operators for Stokhos
#include "Stokhos_StandardStorage.hpp"
#include "Stokhos_StaticStandardStorage.hpp"
#include "Stokhos_StaticFixedStandardStorage.hpp"
#include "Sacado_PCE_OrthogPoly.hpp"
#include "Sacado_ETPCE_OrthogPoly.hpp"

#endif // STOKHOS_SACADO_HPP
