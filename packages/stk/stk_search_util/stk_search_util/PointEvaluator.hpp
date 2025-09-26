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

#ifndef STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_POINTEVALUATOR_HPP_
#define STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_POINTEVALUATOR_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_search/SearchInterface.hpp"             // for EvaluatePointsI...
#include "stk_search_util/MasterElementProvider.hpp"  // for StkProvideMaste...
#include "stk_search_util/spmd/EntityKeyPair.hpp"
#include "stk_topology/topology.hpp"                  // for topology
#include <cstddef>                                    // for size_t
#include <memory>                                     // for shared_ptr
#include <vector>                                     // for vector
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class MetaData; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {

using PointEvaluatorInterface = EvaluatePointsInterface<stk::topology, spmd::EntityKeyPair>;

class MasterElementGaussPointEvaluator : public PointEvaluatorInterface {
 public:
  MasterElementGaussPointEvaluator(stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* coords,
                                   std::shared_ptr<MasterElementProviderInterface> masterElemProvider);

  size_t num_points(const spmd::EntityKeyPair&, const stk::topology&) override;
  void coordinates(const spmd::EntityKeyPair&, size_t, std::vector<double>& coords) override;

  ~MasterElementGaussPointEvaluator() = default;

 private:
  stk::mesh::BulkData& m_bulk;
  const SearchField m_coords;
  std::shared_ptr<MasterElementProviderInterface> m_masterElemProvider;
  const unsigned m_spatialDimension{0};

  mutable std::vector<double> m_coordVector;
  mutable std::vector<double> m_paramCoordVector;
  mutable std::vector<double> m_elementCoords;

  void gather_nodal_coordinates(const stk::mesh::Entity entity, const stk::topology topo) const;
};

class CentroidEvaluator : public PointEvaluatorInterface {
 public:
  CentroidEvaluator(stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* coords);

  size_t num_points(const spmd::EntityKeyPair&, const stk::topology&) override;
  void coordinates(const spmd::EntityKeyPair&, size_t, std::vector<double>& coords) override;

  ~CentroidEvaluator() = default;

 private:
  const stk::mesh::FieldBase* m_coords{nullptr};

  mutable std::vector<double> m_coordVector;
};

class NodeEvaluator : public PointEvaluatorInterface {
 public:
  NodeEvaluator(stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* coords);

  size_t num_points(const spmd::EntityKeyPair&, const stk::topology&) override;
  void coordinates(const spmd::EntityKeyPair&, size_t, std::vector<double>& coords) override;

  ~NodeEvaluator() = default;

 private:
  const stk::mesh::FieldBase* m_coords{nullptr};

  mutable std::vector<double> m_coordVector;
};

} // namespace stk
} // namespace search

#endif /* STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_POINTEVALUATOR_HPP_ */
