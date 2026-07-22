// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_IntegrationDescriptor.hpp"

#include "Panzer_HashUtils.hpp"

#include "Teuchos_Assert.hpp"

namespace panzer
{

IntegrationDescriptor::IntegrationDescriptor()
{
  setup(-1, NONE);
}

IntegrationDescriptor::IntegrationDescriptor(const int cubature_order, const int integration_type, const int side)
{
  setup(cubature_order, integration_type, side);
}

void
IntegrationDescriptor::setup(const int cubature_order, const int integration_type, const int side)
{
  _integration_type = integration_type;
  _cubature_order = cubature_order;
  _side = side;

  if((_integration_type == SIDE) or (_integration_type == CV_BOUNDARY)){
    TEUCHOS_ASSERT(side >= 0);
  } else {
    TEUCHOS_ASSERT(side == -1);
  }
  _key = std::hash<IntegrationDescriptor>()(*this);
}

}

std::size_t
std::hash<panzer::IntegrationDescriptor>::operator()(const panzer::IntegrationDescriptor& desc) const
{
  std::size_t seed = 0;

  panzer::hash_combine(seed,desc.getType());
  panzer::hash_combine(seed,desc.getOrder());
  panzer::hash_combine(seed,desc.getSide());

  return seed;
}
