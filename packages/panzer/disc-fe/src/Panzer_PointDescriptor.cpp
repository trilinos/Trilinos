// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_PointDescriptor.hpp"

#include "Panzer_HashUtils.hpp"

#include "Teuchos_Assert.hpp"

namespace panzer
{

PointDescriptor::
PointDescriptor(const std::string & type,const Teuchos::RCP<PointGenerator> & generator)
{
  setup(type,generator);
}

void
PointDescriptor::
setup(const std::string & type,const Teuchos::RCP<PointGenerator> & generator)
{
  _type = type;
  _generator = generator;

  // sanity check
  TEUCHOS_ASSERT(_generator!=Teuchos::null);

  _key = std::hash<PointDescriptor>()(*this);
}

}

std::size_t
std::hash<panzer::PointDescriptor>::operator()(const panzer::PointDescriptor& desc) const
{
  std::size_t seed = 0;
  const std::string prefix = "point_desc";

  panzer::hash_combine(seed,prefix);   // prevent collisions with other descriptors
  panzer::hash_combine(seed,desc.getType());

  return seed;
}
