// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_GPBuildTransfer_hpp
#define percept_GPBuildTransfer_hpp

#include "GPFromMesh.hpp"
#include "GPToMesh.hpp"

#include <percept/xfer/LinInterp.hpp>
#include <percept/xfer/STKMeshTransferSetup.hpp>

#include <stk_transfer/GeometricTransfer.hpp>
#include <stk_mesh/base/MeshUtils.hpp>

namespace percept
{

  //typedef stk::transfer::GeometricTransfer< class GPProjectToSurface< class GPFromMesh, class GPToMesh > > GPSTKMeshTransfer;
  typedef stk::transfer::GeometricTransfer< class LinInterp< class GPFromMesh, class GPToMesh > > GPSTKMeshTransfer;

  inline
  std::shared_ptr<GPSTKMeshTransfer>
  buildGPSTKMeshTransfer(PerceptMesh& fromMesh,
                       const stk::mesh::PartVector& fromParts,
                       const std::vector<stk::mesh::FieldBase *>& fromFields,
                       PerceptMesh& toMesh,
                       const stk::mesh::PartVector& toParts,
                       const std::vector<stk::mesh::FieldBase *>& toFields,
                       const std::string &transfer_name)
  {
    TransferType transferType = THREED_TO_THREED;

    if (2 == fromMesh.get_spatial_dim() &&
        2 == toMesh.get_spatial_dim()) {
      transferType = TWOD_TO_TWOD;
    }
    else if (2 == fromMesh.get_spatial_dim() &&
             3 == toMesh.get_spatial_dim()) {
      transferType = TWOD_AXI_TO_THREED;
    }

    std::shared_ptr<GPFromMesh >
      from_mesh (new GPFromMesh(fromMesh, fromParts, fromFields));

    std::shared_ptr<GPToMesh >
      to_mesh (new GPToMesh(toMesh, toParts, toFields, transferType));

#if 0
    std::shared_ptr<GPSTKMeshTransfer>
      mesh_transfer(new GPSTKMeshTransfer(from_mesh, to_mesh, transfer_name, 1.5, stk::search::OCTREE));
#else
    std::shared_ptr<GPSTKMeshTransfer>
      mesh_transfer(new GPSTKMeshTransfer(from_mesh, to_mesh, transfer_name));
#endif

    return mesh_transfer;
  }


} //namespace percept

#endif
