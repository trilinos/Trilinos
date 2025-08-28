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
#include "stk_search_util/PointEvaluator.hpp"
#include "stk_search_util/ObjectCoordinates.hpp"                    // for stk_determine_c...
#include <stk_mesh/base/Entity.hpp>                   // for Entity
#include "stk_mesh/base/EntityKey.hpp"                // for EntityKey, oper...
#include "stk_mesh/base/MetaData.hpp"                 // for MetaData
#include "stk_util/util/ReportHandler.hpp"            // for eval_test_condi...
#include "stk_mesh/base/Bucket.hpp"                   // for Bucket
#include <stk_mesh/base/BulkData.hpp>                 // for BulkData
#include <stk_mesh/base/FieldBase.hpp>                // for FieldBase, fiel...
#include <algorithm>                                  // for copy, fill, max
#include <memory>                                     // for shared_ptr, __s...

// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {

MasterElementGaussPointEvaluator::MasterElementGaussPointEvaluator(
    stk::mesh::BulkData& bulk,
    const stk::mesh::FieldBase* coords,
    std::shared_ptr<MasterElementProviderInterface> masterElemProvider)
  : m_bulk(bulk)
  , m_coords(coords)
  , m_masterElemProvider(masterElemProvider)
  , m_spatialDimension(m_bulk.mesh_meta_data().spatial_dimension())
{
}

size_t MasterElementGaussPointEvaluator::num_points(const spmd::EntityKeyPair& k, const stk::topology& topo)
{
  auto meTopo = SearchTopology(topo, k);
  return m_masterElemProvider->num_integration_points(meTopo);
}

void MasterElementGaussPointEvaluator::gather_nodal_coordinates(const stk::mesh::Entity entity,
                                                                [[maybe_unused]] const stk::topology topo) const
{
  unsigned numFieldComponents;
  unsigned numNodes;

  spmd::EntityKeyPair key(entity, m_bulk.entity_key(entity));
  m_masterElemProvider->nodal_field_data(key, m_coords, numFieldComponents, numNodes, m_elementCoords);

  STK_ThrowAssertMsg(numFieldComponents == m_spatialDimension, "Invalid coordinate field: " << m_coords.name());
  STK_ThrowAssertMsg(numNodes == topo.num_nodes(), "Mismatch between key: " << m_bulk.entity_key(entity) << " and topology: " << topo.name());
}

void MasterElementGaussPointEvaluator::coordinates(const spmd::EntityKeyPair& k, size_t pointIndex, std::vector<double>& coords)
{
  stk::mesh::Entity elem = k;
  stk::mesh::Bucket& bucket = m_bulk.bucket(elem);
  stk::topology topo = bucket.topology();
  auto meTopo = SearchTopology(topo, k, &bucket);
  size_t numParCoor = m_masterElemProvider->num_parametric_coordinates(meTopo);

  m_masterElemProvider->integration_points(meTopo, m_coordVector);
  const double* gpLoc = &(m_coordVector[pointIndex * numParCoor]);

  m_paramCoordVector.assign(gpLoc, gpLoc+numParCoor);

  unsigned nDim = m_bulk.mesh_meta_data().spatial_dimension();

  gather_nodal_coordinates(elem, topo);

  m_masterElemProvider->evaluate_field(meTopo, m_paramCoordVector, nDim, m_elementCoords, coords);
}

CentroidEvaluator::CentroidEvaluator(stk::mesh::BulkData& /*bulk*/, const stk::mesh::FieldBase* coords)
  : m_coords(coords)
{
  STK_ThrowRequireMsg(m_coords->entity_rank() == stk::topology::NODE_RANK ||
                      m_coords->entity_rank() == stk::topology::ELEM_RANK,
                      "Centroid evaluator coordinate field must be either NODE_RANK or ELEM_RANK");
}

size_t CentroidEvaluator::num_points(const spmd::EntityKeyPair& /*k*/, const stk::topology& /*topo*/)
{
  return 1u;
}

void CentroidEvaluator::coordinates(const spmd::EntityKeyPair& k, size_t /*pointIndex*/, std::vector<double>& coords)
{
  stk::mesh::Entity elem = k;
  const unsigned nDim = m_coords->mesh_meta_data().spatial_dimension();

  if(m_coords->entity_rank() == stk::topology::NODE_RANK) {
    stk::search::determine_centroid(nDim, elem, *m_coords, coords);
  }
  else if(m_coords->entity_rank() == stk::topology::ELEM_RANK) {
    const double* coor = static_cast<double*>(stk::mesh::field_data(*m_coords, elem));
    coords.assign(coor, coor+nDim);
  }
}

NodeEvaluator::NodeEvaluator(stk::mesh::BulkData& /*bulk*/, const stk::mesh::FieldBase* coords)
  : m_coords(coords)
{
  STK_ThrowRequireMsg(m_coords->entity_rank() == stk::topology::NODE_RANK,
                      "Centroid evaluator coordinate field must be NODE_RANK");
}

size_t NodeEvaluator::num_points(const spmd::EntityKeyPair& /*k*/, const stk::topology& /*topo*/)
{
  return 1u;
}

void NodeEvaluator::coordinates(const spmd::EntityKeyPair& k, size_t /*pointIndex*/, std::vector<double>& coords)
{
  stk::mesh::Entity node = k;
  const unsigned nDim = m_coords->mesh_meta_data().spatial_dimension();

  const double* coor = static_cast<double*>(stk::mesh::field_data(*m_coords, node));
  coords.assign(coor, coor+nDim);
}

} // namespace search
} // namespace stk
