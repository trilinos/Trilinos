// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_BasisValues2.hpp"
#include "Panzer_BasisValues2_impl.hpp"

namespace panzer {

#define BASIS_VALUES_INSTANTIATION(SCALAR) \
template class BasisValues2<SCALAR>;

BASIS_VALUES_INSTANTIATION(panzer::Traits::RealType)

// Disabled due to long build times on cuda (30+ minutes for this
// instantiaiton alone. If we need sensitivities wrt coordinates, we
// can reenable.
//BASIS_VALUES_INSTANTIATION(panzer::Traits::FadType)

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
BASIS_VALUES_INSTANTIATION(panzer::Traits::HessianType)
#endif

}
