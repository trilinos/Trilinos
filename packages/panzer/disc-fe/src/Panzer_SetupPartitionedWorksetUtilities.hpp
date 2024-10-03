// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_SETUP_PARTITIONED_WORKSET_UTILITIES_HPP
#define PANZER_SETUP_PARTITIONED_WORKSET_UTILITIES_HPP

#include "Panzer_Workset.hpp"

#include "Teuchos_RCP.hpp"

#include <vector>

namespace panzer
{
  class WorksetDescriptor;
  struct WorksetNeeds;
  struct LocalMeshInfo;
}

namespace panzer
{

  /** Build worksets for a partitioned mesh
   *
   * \param[in] mesh_info Mesh info object
   * \param[in] description Description of workset
   * \param[in] orientations Orientations object for correcting basis information
   *
   * \returns vector of worksets for the corresponding element block.
   */
  Teuchos::RCP<std::vector<panzer::Workset> >
  buildPartitionedWorksets(const panzer::LocalMeshInfo & mesh_info,
                           const panzer::WorksetDescriptor & description,
                           const Teuchos::RCP<const OrientationsInterface> & orientations = Teuchos::null);

}

#endif
