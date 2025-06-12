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

#ifndef STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_POINTPROJECTION_HPP_
#define STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_POINTPROJECTION_HPP_

#include <string>
#include <vector>
#include <memory>

#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/EntityKey.hpp"
#include "stk_mesh/base/FieldState.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "stk_mesh/base/Types.hpp"
#include "stk_search/Box.hpp"
#include "stk_search/DistanceComparison.hpp"
#include "stk_search/IdentProc.hpp"
#include "stk_search/Point.hpp"
#include "stk_search/Sphere.hpp"

#include "stk_search_util/MasterElementProvider.hpp"

#include <Kokkos_Core_fwd.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { struct Entity; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class MetaData; } }

namespace stk {
namespace search {

struct ProjectionData {
  const stk::mesh::BulkData& bulk;
  const std::shared_ptr<MasterElementProviderInterface> masterElemProvider;
  const std::vector<double>& evalPoint;
  const stk::mesh::FieldBase& coordinateField;
  const SearchField projectedCoordinateField;

  mutable std::vector<spmd::EntityKeyPair> scratchKeySpace;

  ProjectionData(const stk::mesh::BulkData& bulk_,
                 const std::shared_ptr<MasterElementProviderInterface> masterElemProvider_,
                 const std::vector<double>& targetCoordinates_,
                 const stk::mesh::FieldBase& sourceCoordinateField_)
                 : bulk(bulk_),
                   masterElemProvider(masterElemProvider_),
                   evalPoint(targetCoordinates_),
                   coordinateField(sourceCoordinateField_),
                   projectedCoordinateField(&sourceCoordinateField_)
                 {}

  ProjectionData(const stk::mesh::BulkData& bulk_,
                 const std::shared_ptr<MasterElementProviderInterface> masterElemProvider_,
                 const std::vector<double>& targetCoordinates_,
                 const stk::mesh::FieldBase& sourceCoordinateField_,
                 double (*transformFunc_)(double),
                 double defaultValue_)
                 : bulk(bulk_),
                   masterElemProvider(masterElemProvider_),
                   evalPoint(targetCoordinates_),
                   coordinateField(sourceCoordinateField_),
                   projectedCoordinateField(&sourceCoordinateField_, transformFunc_, defaultValue_)
                 {}
};

struct ProjectionResult {
  double geometricDistanceSquared;
  double parametricDistance;
  std::vector<double> parametricCoords;
  bool doneProjection = false;
};

void project_to_entity(const ProjectionData& data, stk::mesh::Entity entity, ProjectionResult& result);

void project_to_closest_side(const ProjectionData& data, stk::mesh::Entity entity, ProjectionResult& result);

void project_to_closest_face(const ProjectionData& data, stk::mesh::Entity entity, ProjectionResult& result);

void truncate_to_entity(const ProjectionData& data, stk::mesh::Entity entity, ProjectionResult& result);

double fill_projected_object_location(const ProjectionData& data, stk::mesh::Entity entity,
                                      const std::vector<double>& location, std::vector<double>& parametricCoordinates);

double compute_parametric_distance(const ProjectionData& data, stk::mesh::Entity entity,
                                   const std::vector<double>& entityParCoord, std::vector<double>& evaluatedLocation);

} // end namespace search
} // end namespace stk

#endif /* STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_POINTPROJECTION_HPP_ */
