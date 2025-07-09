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
#include "stk_search_util/FindParametricCoordinates.hpp"
#include "stk_mesh/base/Bucket.hpp"           // for Bucket
#include <stk_mesh/base/BulkData.hpp>         // for BulkData
#include <stk_mesh/base/Entity.hpp>           // for Entity
#include "stk_mesh/base/MetaData.hpp"         // for MetaData
#include "stk_search/DistanceComparison.hpp"  // for distance
#include "stk_util/util/ReportHandler.hpp"
#include "stk_mesh/base/EntityKey.hpp"        // for operator<<, EntityKey
#include "stk_topology/topology.hpp"          // for operator<<, topology
#include <algorithm>                          // for copy, fill, max
#include <memory>                             // for shared_ptr, __shared_pt...
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {

MasterElementParametricCoordsFinder::MasterElementParametricCoordsFinder(
    stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* coords,
    std::shared_ptr<MasterElementProviderInterface> masterElemProvider,
    const double parametricTolerance)
  : m_bulk(bulk)
  , m_meta(bulk.mesh_meta_data())
  , m_coords(coords)
  , m_masterElemProvider(masterElemProvider)
  , m_parametricTolerance(parametricTolerance)
  , m_spatialDimension(m_meta.spatial_dimension())
{
}

void MasterElementParametricCoordsFinder::gather_nodal_coordinates(const spmd::EntityKeyPair& key,
                                                                   [[maybe_unused]] const stk::topology topo) const
{
  unsigned numFieldComponents;
  unsigned numNodes;

  m_masterElemProvider->nodal_field_data(key, m_coords, numFieldComponents, numNodes, m_elementCoords);

  STK_ThrowAssertMsg(numFieldComponents == m_spatialDimension, "Invalid coordinate field: " << m_coords.name());
  STK_ThrowAssertMsg(numNodes == topo.num_nodes(), "Mismatch between key: " << key << " and topology: " << topo.name());
}

void
MasterElementParametricCoordsFinder::find_parametric_coords(const spmd::EntityKeyPair& k, const std::vector<double>& evalPoint,
                                                            std::vector<double>& paramCoords, double& paramDistance,
                                                            bool& isWithinParametricTolerance) const
{
  stk::mesh::Entity entity = k;

  const stk::mesh::Bucket& bucket = m_bulk.bucket(entity);
  const stk::topology topo = bucket.topology();

  gather_nodal_coordinates(k, topo);

  auto meTopo = SearchTopology(topo, k, &bucket);
  m_masterElemProvider->find_parametric_coordinates(meTopo, m_spatialDimension, m_elementCoords, evalPoint,
                                                    paramCoords, paramDistance);

  isWithinParametricTolerance = paramDistance <= (1 + m_parametricTolerance);
}

void
MasterElementParametricCoordsFinder::evaluate_parametric_coords(const spmd::EntityKeyPair& k,
                                                                const std::vector<double>& paramCoords,
                                                                std::vector<double>& evalPoint) const
{
  stk::mesh::Entity elem = k;

  const stk::mesh::Bucket& bucket = m_bulk.bucket(elem);
  const stk::topology topo = bucket.topology();

  gather_nodal_coordinates(k, topo);

  auto meTopo = SearchTopology(topo, k, &bucket);
  m_masterElemProvider->evaluate_field(meTopo, paramCoords, m_spatialDimension, m_elementCoords, evalPoint);
}

NodeParametricCoordsFinder::NodeParametricCoordsFinder(stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* coords)
  : m_meta(bulk.mesh_meta_data())
  , m_coords(coords)
{
}

void NodeParametricCoordsFinder::find_parametric_coords(const spmd::EntityKeyPair& k, const std::vector<double>& toCoords,
                                                        std::vector<double>& paramCoords, double& paramDistance,
                                                        bool& isWithinParametricTolerance) const
{
  STK_ThrowRequire(k.rank() == stk::topology::NODE_RANK);

  stk::mesh::Entity theNode = k;

  const double* fromCoords = static_cast<const double *>(stk::mesh::field_data(*m_coords, theNode));
  unsigned nDim = m_meta.spatial_dimension();

  paramCoords.assign(nDim, 0.0);

  paramDistance = distance(nDim, fromCoords, toCoords.data());

  isWithinParametricTolerance = true;
}

void NodeParametricCoordsFinder::evaluate_parametric_coords(const spmd::EntityKeyPair& k,
                                                            const std::vector<double>& /*paramCoords*/,
                                                            std::vector<double>& evalPoint) const
{
  STK_ThrowRequire(k.rank() == stk::topology::NODE_RANK);

  stk::mesh::Entity theNode = k;

  const double* fromCoords = static_cast<const double *>(stk::mesh::field_data(*m_coords, theNode));
  unsigned nDim = m_meta.spatial_dimension();

  evalPoint.assign(fromCoords, fromCoords + nDim);
}

} // namespace search
} // namespace stk
