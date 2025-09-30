// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_transfer_util/spmd/GeometricTransferUtils.hpp"
#include "stk_transfer_util/spmd/ElementRecvMesh.hpp"
#include "stk_search_util/MeshUtility.hpp"
#include "stk_search_util/ObjectCoordinates.hpp"
#include "stk_io/IossBridge.hpp"
#include <stk_mesh/base/Entity.hpp>         // for Entity
#include <stk_mesh/base/FieldBase.hpp>      // for FieldBase
#include "stk_mesh/base/MetaData.hpp"       // for MetaData
#include "stk_mesh/base/Selector.hpp"       // for Selector
#include "stk_mesh/base/Bucket.hpp"         // for Bucket
#include <stk_mesh/base/BulkData.hpp>       // for BulkData
#include "stk_mesh/base/CompositeRank.hpp"                 // for CompositeRank
#include "stk_topology/topology.hpp"        // for topology, operator<<, top...
#include "stk_util/util/ReportHandler.hpp"  // for eval_test_condition, STK_...

#include <algorithm>                        // for max, sort, fill
#include <cstddef>                          // for size_t
#include <memory>                           // for shared_ptr, __shared_ptr_...
#include <type_traits>
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer {
namespace spmd {

stk::mesh::EntityRank
get_transfer_rank(const stk::mesh::MetaData& meta,
                  const std::vector<std::string>& partNames,
                  const std::string& transferName)
{
  STK_ThrowRequireMsg(!partNames.empty(), " Error:: no parts specified.\n"
                                          << " Check Specification of transfer named: " << transferName
                                          << "\n");

  std::vector<std::string>::const_iterator partNameIter = partNames.begin();

  const std::string partName(*partNameIter);

  stk::mesh::EntityRank transferRank = stk::topology::INVALID_RANK;

  stk::mesh::Part* part = meta.get_part(partName);

  if(part != nullptr) {
    transferRank = stk::mesh::CompositeRank::get_rank(part);
    STK_ThrowRequireMsg(transferRank != stk::topology::INVALID_RANK, "Found part with invalid composite rank\n");
  }
  else {
    std::ostringstream oss;
    oss << "Part named " << partName << " in transfer named: " << transferName << " not found. ";
    throw std::runtime_error(oss.str());
  }

  if(transferRank != stk::topology::FACE_RANK && transferRank != stk::topology::ELEMENT_RANK &&
     transferRank != stk::topology::NODE_RANK && transferRank != stk::topology::EDGE_RANK) {
    stk::RuntimeWarning() << "Part named " << partName << " in transfer named: " << transferName
                          << " has an invalid entity rank: " << (transferRank) << "\n"
                          << " This is based on first checking the part type of : "
                          << (part->primary_entity_rank()) << "\n"
                          << " Expected one of: " << (stk::topology::NODE_RANK) << ", " << (stk::topology::FACE_RANK)
                          << ", " << (stk::topology::EDGE_RANK) << ", " << (stk::topology::ELEMENT_RANK) << "\n"
                          << " Check Specification of transfer with name " << transferName << "\n"
                          << StackTrace;
    return stk::topology::INVALID_RANK;
  }

  ++partNameIter;
  for(; partNameIter != partNames.end(); ++partNameIter) {
    const std::string otherPartName(*partNameIter);
    stk::mesh::Part* otherPart = meta.get_part(otherPartName);
    if(otherPart != nullptr) {
      stk::mesh::EntityRank otherEntityRank = stk::mesh::CompositeRank::get_rank(otherPart);
      if(transferRank != otherEntityRank) {
        std::ostringstream oss;
        oss << "Part named " << otherPartName << " in transfer named: " << transferName
            << " has an invalid Mesh Object type: " << (otherEntityRank) << "\n"
            << " Expected " << (transferRank) << " Because first part had this type and all "
            << " parts must be of consistent type \n"
            << StackTrace;
        throw std::runtime_error(oss.str());
      }
    }
  }
  return (transferRank);
}

std::vector<std::string> prune_empty_assemblies(const stk::mesh::MetaData& meta,
                                                const std::vector<std::string>& partNames,
                                                const std::string& transferName)
{
  std::vector<std::string> prunedPartNames;

  for(const std::string& partName : partNames) {
    stk::mesh::Part* meshPart = meta.get_part(partName);
    STK_ThrowRequireMsg(nullptr != meshPart,
                    "Could not find part named: '" << partName << "' in transfer named: " << transferName);

    if(stk::io::is_part_assembly_io_part(*meshPart)) {
      if(!stk::io::get_unique_leaf_parts(*meshPart).empty()) {
        prunedPartNames.push_back(partName);
      }
      else {
        stk::RuntimeWarning() << "Attempting to use empty assembly: '" << meshPart->name()
                              << "' for transfer named: '" << transferName;
      }
    }
    else {
      prunedPartNames.push_back(partName);
    }
  }

  return prunedPartNames;
}

stk::mesh::PartVector get_parts(const stk::mesh::MetaData& meta,
                                const std::vector<std::string>& partNames,
                                const std::string& transferName,
                                const stk::mesh::Part* defaultPart)
{
  stk::mesh::PartVector parts;

  std::vector<std::string> prunedPartNames = prune_empty_assemblies(meta, partNames, transferName);

  if(!prunedPartNames.empty()) {
    (void)get_transfer_rank(meta, prunedPartNames, transferName);

    for(const std::string& partName : prunedPartNames) {
      stk::mesh::Part* meshPart = meta.get_part(partName);
      STK_ThrowRequireMsg(nullptr != meshPart,
                          "Could not find part named: " << partName << " in transfer named: " << transferName);
      parts.push_back(meshPart);
    }
  }
  else if(nullptr != defaultPart) {
    stk::mesh::Part* nonConstDefaultPart = const_cast<stk::mesh::Part*>(defaultPart);
    parts.push_back(nonConstDefaultPart);
  }

  return parts;
}

double mesh_bounding_box_geometric_tolerance(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* coordinateField, double scaleFactor)
{
  int nDim = bulk.mesh_meta_data().spatial_dimension();
  auto mesh_box = stk::search::get_mesh_bounding_box(coordinateField, bulk);
  double diag = 0.0;
  for(int d = 0; d < nDim; ++d) {
    const double dd = mesh_box.max_corner()[d] - mesh_box.min_corner()[d];
    diag += dd * dd;
  }
  diag = std::sqrt(diag);
  return scaleFactor * diag;
}

} // namespace spmd
} // namespace transfer
} // namespace stk
