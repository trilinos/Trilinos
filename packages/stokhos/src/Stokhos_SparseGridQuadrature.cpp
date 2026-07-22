// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_SparseGridQuadrature.hpp"

#ifdef HAVE_STOKHOS_DAKOTA

// Instantiate the static sparse grid quadrature pointer for int/double
// (anything else and the user has to do it)
template <>
Stokhos::SparseGridQuadrature<int,double>*
Stokhos::SparseGridQuadrature<int,double>::sgq(NULL);

#endif // HAVE_STOKHOS_DAKOTA
