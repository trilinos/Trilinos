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

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_DotProduct_decl.hpp"
#include "Panzer_DotProduct_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(panzer::DotProduct)

#define DOT_PRODUCT_NON_MEMBER_CONST_ETI(EVALT,TRAITST) \
template  \
Teuchos::RCP<panzer::DotProduct<EVALT,TRAITST> > \
panzer::buildEvaluator_DotProduct<EVALT,TRAITST>(const std::string &, \
                                                 const panzer::PointRule &, \
                                                 const std::string &, \
                                                 const std::string &, \
                                                 double multiplier, \
                                                 const std::string &);

DOT_PRODUCT_NON_MEMBER_CONST_ETI(panzer::Traits::Residual,panzer::Traits)
DOT_PRODUCT_NON_MEMBER_CONST_ETI(panzer::Traits::Tangent,panzer::Traits)
DOT_PRODUCT_NON_MEMBER_CONST_ETI(panzer::Traits::Jacobian,panzer::Traits)
#ifdef Panzer_BUILD_HESSIAN_SUPPORT
DOT_PRODUCT_NON_MEMBER_CONST_ETI(panzer::Traits::Hessian,panzer::Traits)
#endif

#endif
