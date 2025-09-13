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
#include "stk_search_util/ExternalPointHandler.hpp"
#include "stk_search_util/PointProjection.hpp"// for ProjectionResult, Proje...
#include "stk_search_util/MeshUtility.hpp"    // for compute_parametric_dist...
#include "stk_search/DistanceComparison.hpp"  // for distance_sq
#include <stk_mesh/base/Entity.hpp>           // for Entity
#include <stk_mesh/base/BulkData.hpp>         // for BulkData
#include <stk_mesh/base/MetaData.hpp>         // for MetaData
#include <limits>                             // for numeric_limits

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {

ExternalPointNoOpHandler::ExternalPointNoOpHandler(double /*tol*/) {}

bool ExternalPointNoOpHandler::handle_point(const stk::search::spmd::EntityKeyPair& /*k*/,
                                            const std::vector<double>& /*toCoords*/,
                                            std::vector<double>& /*parametricCoords*/,
                                            double& geometricDistanceSquared,
                                            bool& isWithinGeometricTolerance) const
{
  isWithinGeometricTolerance = false;
  geometricDistanceSquared = std::numeric_limits<double>::max();
  return false;
}

MasterElementExternalPointProjection::MasterElementExternalPointProjection(
    stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* coords,
    std::shared_ptr<MasterElementProviderInterface> masterElemProvider,
    const double parametricTol, const double geometricTol)
  : m_bulk(bulk)
  , m_coords(coords)
  , m_masterElemProvider(masterElemProvider)
  , m_parametricTol(parametricTol)
  , m_geometricTol(geometricTol)
{
}

bool MasterElementExternalPointProjection::handle_point(const stk::search::spmd::EntityKeyPair& k, const std::vector<double>& toCoords,
                                                        std::vector<double>& parametricCoords,
                                                        double& geometricDistanceSquared,
                                                        bool& isWithinGeometricTolerance) const
{
  stk::mesh::Entity elem = k;
  project_point_to_element_boundary(elem, toCoords, parametricCoords, geometricDistanceSquared);

  isWithinGeometricTolerance = (geometricDistanceSquared <= m_geometricTol * m_geometricTol);
  return true;
}

void MasterElementExternalPointProjection::project_point_to_element_boundary(stk::mesh::Entity elem, const std::vector<double>& tocoords,
                                                                             std::vector<double>& parametricCoords,
                                                                             double& nearestDistance) const
{
  unsigned nDim = m_bulk.mesh_meta_data().spatial_dimension();
  std::vector<double> interpolatedLoc(nDim);

  stk::search::ProjectionData data(m_bulk, m_masterElemProvider, tocoords, *m_coords);

  nearestDistance = stk::search::compute_parametric_distance(data, elem, parametricCoords, interpolatedLoc);
  if(nearestDistance <= (1 + m_parametricTol)) {
    nearestDistance = stk::search::distance_sq(nDim, interpolatedLoc.data(), tocoords.data());
    return;
  }

  stk::search::ProjectionResult result;
  stk::search::project_to_closest_face(data, elem, result);

  if(result.doneProjection) {
    nearestDistance = result.geometricDistanceSquared;
    parametricCoords.swap(result.parametricCoords);
  }
}

MasterElementExternalPointTruncation::MasterElementExternalPointTruncation(
    stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* coords,
    std::shared_ptr<MasterElementProviderInterface> masterElemProvider, const double geometricTol)
  : m_bulk(bulk)
  , m_coords(coords)
  , m_masterElemProvider(masterElemProvider)
  , m_geometricTol(geometricTol)
{
}

bool MasterElementExternalPointTruncation::handle_point(const stk::search::spmd::EntityKeyPair& k, const std::vector<double>& toCoords,
                                                        std::vector<double>& parametricCoords,
                                                        double& geometricDistanceSquared,
                                                        bool& isWithinGeometricTolerance) const
{
  stk::mesh::Entity elem = k;
  truncate_point_to_element_boundary(elem, toCoords, parametricCoords, geometricDistanceSquared);

  isWithinGeometricTolerance = (geometricDistanceSquared <= m_geometricTol * m_geometricTol);
  return true;
}

void MasterElementExternalPointTruncation::truncate_point_to_element_boundary(stk::mesh::Entity elem, const std::vector<double>& tocoords,
                                                                              std::vector<double>& parametricCoords,
                                                                              double& nearestDistance) const
{
  stk::search::ProjectionData data(m_bulk, m_masterElemProvider, tocoords, *m_coords);
  stk::search::ProjectionResult result;
  stk::search::truncate_to_entity(data, elem, result);

  nearestDistance = result.geometricDistanceSquared;
  parametricCoords.swap(result.parametricCoords);
}

} // namespace search
} // namespace stk
