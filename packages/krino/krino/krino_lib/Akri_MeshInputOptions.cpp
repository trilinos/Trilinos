// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MeshInputOptions.hpp>
#include <Akri_DiagWriter.hpp>
#include <stk_util/environment/RuntimeWarning.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace krino {

MeshInputOptions::Registry MeshInputOptions::the_registry;

std::shared_ptr<MeshInputOptions>
MeshInputOptions::get_or_create(const std::string & model_name)
{
  std::shared_ptr<MeshInputOptions> options = std::make_shared<MeshInputOptions>(model_name);
  std::pair<Registry::iterator,bool> result = the_registry.insert(options);

  if ( ! result.second ) {
    stk::RuntimeWarningAdHoc() << "A mesh database named '" << model_name
        << "' has already been defined and will be reused";
  }

  return *result.first;
}

MeshInputOptions *
MeshInputOptions::get(const std::string &model_name)
{
  MeshInputOptions *options = nullptr;
  std::shared_ptr<MeshInputOptions> tmp = std::make_shared<MeshInputOptions>(model_name);

  Registry::iterator iter = the_registry.find(tmp);

  if ( iter != the_registry.end() )
    options = iter->get() ;

  return options;
}

int MeshInputOptions::get_generated_mesh_spatial_dimension() const
{
  STK_ThrowRequire(use_generated_mesh());
  return (my_generated_mesh_domain_type == GENERATED_2D_MESH_FOR_INTERFACE_BOUNDING_BOX || my_generated_mesh_domain.size() == 4) ? 2 : 3;
}

} // namespace krino
