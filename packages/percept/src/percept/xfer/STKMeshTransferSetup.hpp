// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_BuildTransfer_hpp
#define percept_BuildTransfer_hpp

#include <percept/xfer/FromMesh.hpp>
#include <percept/xfer/ToMesh.hpp>
#include <percept/xfer/LinInterp.hpp>

#include <stk_transfer/GeometricTransfer.hpp>
#include <stk_mesh/base/MeshUtils.hpp>

namespace percept
{

typedef FromMesh<> FMesh;
typedef ToMesh<> TMesh;

typedef stk::transfer::GeometricTransfer< class LinInterp<  FMesh, TMesh > > STKMeshTransfer;


template<class SMT>
inline
std::shared_ptr<SMT>
buildSTKMeshTransfer(stk::mesh::BulkData       &bulkData_from,
         stk::mesh::FieldBase * coordinates_from,
		     stk::mesh::FieldBase * field_from,
		     stk::mesh::BulkData       &bulkData_to,
         stk::mesh::FieldBase * coordinates_to,
		     stk::mesh::FieldBase * field_to,
		     const std::string &transfer_name,
                     SrcFieldType srcFieldType=SRC_FIELD,
                     const double expansion_factor = 1.5)
{
  TransferType transferType = THREED_TO_THREED;

  if (2 == bulkData_from.mesh_meta_data().spatial_dimension() &&
      2 == bulkData_to.mesh_meta_data().spatial_dimension()) {
    transferType = TWOD_TO_TWOD;
  }
  else if (2 == bulkData_from.mesh_meta_data().spatial_dimension() &&
	   3 == bulkData_to.mesh_meta_data().spatial_dimension()) {
    transferType = TWOD_AXI_TO_THREED;
  }

  std::shared_ptr<FMesh >
    from_mesh (new FMesh(bulkData_from,
                         coordinates_from,
                         field_from,
                         bulkData_from.parallel()));

  std::shared_ptr<TMesh >
    to_mesh (new TMesh(bulkData_to,
                       coordinates_to,
                       field_to,
                       bulkData_to.parallel(),
                       transferType,
                       srcFieldType));

  std::shared_ptr<SMT>
    mesh_transfer(new SMT(from_mesh, to_mesh, transfer_name, expansion_factor));

  return mesh_transfer;
}

template<class SMT>
inline
void
initializeSTKMeshTransfer(SMT * mesh_transfer)
{
  mesh_transfer->coarse_search();

  mesh_transfer->meshA()->fromBulkData_.modification_begin();

  // based on Transfer::change_ghosting() in conchas2
  //  which calls Transfer::ghost_from_elements();
  STKMeshTransfer::MeshA::EntityProcVec entity_keys;
  mesh_transfer->determine_entities_to_copy(entity_keys);
  mesh_transfer->meshA()->update_ghosting(entity_keys);

  stk::mesh::fixup_ghosted_to_shared_nodes(mesh_transfer->meshA()->fromBulkData_);
  mesh_transfer->meshA()->fromBulkData_.modification_end();

  mesh_transfer->local_search();
}

} //namespace percept

#endif
