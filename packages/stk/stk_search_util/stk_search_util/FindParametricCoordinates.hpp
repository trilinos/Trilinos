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

#ifndef STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_FINDPARAMETRICCOORDINATES_HPP_
#define STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_FINDPARAMETRICCOORDINATES_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_search/SearchInterface.hpp"             // for FindParametricC...
#include "stk_search_util/MasterElementProvider.hpp"  // for ProvideMaste...
#include "stk_search_util/spmd/EntityKeyPair.hpp"
#include <memory>                                     // for shared_ptr
#include <vector>                                     // for vector
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class MetaData; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace search {

using FindParametricCoordsInterface = FindParametricCoordinatesInterface<spmd::EntityKeyPair>;

class MasterElementParametricCoordsFinder : public FindParametricCoordsInterface {
 public:
  MasterElementParametricCoordsFinder(stk::mesh::BulkData& bulk,
                                      const stk::mesh::FieldBase* coords,
                                      std::shared_ptr<MasterElementProviderInterface> masterElemProvider,
                                      const double parametricTolerance);

  // TODO: define an interface for the return type (this currently assumes Fmwk-type arguments)
  void find_parametric_coords(const spmd::EntityKeyPair& k, const std::vector<double>& evalPoint,
                              std::vector<double>& paramCoords, double& paramDistance,
                              bool& isWithinParametricTolerance) const override;
  void evaluate_parametric_coords(const spmd::EntityKeyPair& k, const std::vector<double>& paramCoords,
                                  std::vector<double>& evalPoint) const override;

  ~MasterElementParametricCoordsFinder() = default;

 private:
  stk::mesh::BulkData& m_bulk;
  stk::mesh::MetaData& m_meta;
  const SearchField m_coords;
  std::shared_ptr<MasterElementProviderInterface> m_masterElemProvider;
  double m_parametricTolerance;
  const unsigned m_spatialDimension{0};
  mutable std::vector<double> m_elementCoords;

  void gather_nodal_coordinates(const spmd::EntityKeyPair& key, const stk::topology topo) const;
};

class NodeParametricCoordsFinder : public FindParametricCoordsInterface {
 public:
  NodeParametricCoordsFinder(stk::mesh::BulkData& bulk, const stk::mesh::FieldBase* coords);

  void find_parametric_coords(const spmd::EntityKeyPair& k, const std::vector<double>& toCoords,
                              std::vector<double>& paramCoords, double& paramDistance,
                              bool& isWithinParametricTolerance) const override;
  void evaluate_parametric_coords(const spmd::EntityKeyPair& k, const std::vector<double>& paramCoords,
                                  std::vector<double>& evalPoint) const override;

  ~NodeParametricCoordsFinder() = default;

 private:
  const stk::mesh::MetaData& m_meta;
  const stk::mesh::FieldBase* m_coords{nullptr};
};

} // namespace search
} // namespace stk

#endif /* STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_FINDPARAMETRICCOORDINATES_HPP_ */
