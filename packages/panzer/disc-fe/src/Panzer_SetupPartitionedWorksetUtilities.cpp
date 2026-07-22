// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_SetupPartitionedWorksetUtilities.hpp"

#include "Panzer_LocalPartitioningUtilities.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_WorksetNeeds.hpp"
#include "Panzer_WorksetDescriptor.hpp"

namespace panzer
{

namespace
{
void
convertMeshPartitionToWorkset(const panzer::LocalMeshPartition & partition,
                              const Teuchos::RCP<const OrientationsInterface> & orientations,
                              panzer::Workset & workset)
{
  WorksetOptions options;
  options.side_assembly_ = false;
  options.align_side_points_ = false;
  options.orientations_ = orientations;

  // Construct the workset from the partition
  workset.setup(partition, options);

}
}

Teuchos::RCP<std::vector<panzer::Workset> >  
buildPartitionedWorksets(const panzer::LocalMeshInfo & mesh_info,
                         const panzer::WorksetDescriptor & description,
                         const Teuchos::RCP<const OrientationsInterface> & orientations)
{
  Teuchos::RCP<std::vector<panzer::Workset> > worksets = Teuchos::rcp(new std::vector<panzer::Workset>());

  // Make sure it makes sense to partition
  TEUCHOS_ASSERT(description.requiresPartitioning());

  // Each partition represents a chunk of the mesh
  std::vector<panzer::LocalMeshPartition> partitions;
  panzer::generateLocalMeshPartitions(mesh_info, description, partitions);

  int i=0;
  for(const auto & partition : partitions){
    worksets->push_back(panzer::Workset());
    convertMeshPartitionToWorkset(partition, orientations, worksets->back());

    // We hash in a unique id the the given workset
    size_t id = std::hash<WorksetDescriptor>()(description);
    panzer::hash_combine(id, i++);
    worksets->back().setIdentifier(id);
  }

  return worksets;
}

}
