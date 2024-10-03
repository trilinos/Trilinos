// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_HPP
#define STOKHOS_HPP

// Version string
#include "Stokhos_Version.hpp"

// Bases
#include "Stokhos_LegendreBasis.hpp"
#include "Stokhos_ClenshawCurtisLegendreBasis.hpp"
#include "Stokhos_GaussPattersonLegendreBasis.hpp"
#include "Stokhos_HermiteBasis.hpp"
#include "Stokhos_JacobiBasis.hpp"
#include "Stokhos_RysBasis.hpp"
#include "Stokhos_DiscretizedStieltjesBasis.hpp"
#include "Stokhos_PecosOneDOrthogPolyBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_TensorProductBasis.hpp"
#include "Stokhos_TotalOrderBasis.hpp"
#include "Stokhos_SmolyakBasis.hpp"
#include "Stokhos_BasisFactory.hpp"

// Expansion methods
#include "Stokhos_OrthogPolyApprox.hpp"
#include "Stokhos_ConstantOrthogPolyExpansion.hpp"
#include "Stokhos_AlgebraicOrthogPolyExpansion.hpp"
#include "Stokhos_QuadOrthogPolyExpansion.hpp"
#include "Stokhos_PseudoSpectralOrthogPolyExpansion.hpp"
//#include "Stokhos_DerivOrthogPolyExpansion.hpp"
#include "Stokhos_ForUQTKOrthogPolyExpansion.hpp"
#include "Stokhos_ExpansionFactory.hpp"

// Quadrature methods
#include "Stokhos_TensorProductQuadrature.hpp"
#include "Stokhos_SparseGridQuadrature.hpp"
#include "Stokhos_SmolyakSparseGridQuadrature.hpp"
#include "Stokhos_AnisoSparseGridQuadrature.hpp"
#include "Stokhos_QuadratureFactory.hpp"

// Pseudospectral methods
#include "Stokhos_TensorProductPseudoSpectralOperator.hpp"
#include "Stokhos_QuadraturePseudoSpectralOperator.hpp"
#include "Stokhos_SmolyakPseudoSpectralOperator.hpp"

// Tensors
#include "Stokhos_Sparse3Tensor.hpp"
#include "Stokhos_Dense3Tensor.hpp"

#endif // STOKHOS_HPP
