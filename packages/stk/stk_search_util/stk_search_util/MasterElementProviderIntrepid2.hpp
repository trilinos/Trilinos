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

#ifndef STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_MASTERELEMENTPROVIDERINTREPID2_HPP_
#define STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_MASTERELEMENTPROVIDERINTREPID2_HPP_

#include "MasterElementProvider.hpp"
#include "Intrepid2BasisFactory.hpp"
#include <Intrepid2_DefaultCubatureFactory.hpp>
#include "Intrepid2BasisFactory.hpp"

namespace stk::search {

class MasterElementProviderIntrepid2 : public MasterElementProviderInterface
{
  public:
  using Topology  = MasterElementProviderInterface::Topology;
  using EntityKey = MasterElementProviderInterface::EntityKey;
  using Field     = MasterElementProviderInterface::Field;

  MasterElementProviderIntrepid2(bool useCompositeTet10 = false) :
    m_lastTopology(stk::topology::INVALID_TOPOLOGY),
    m_numEvalPoints(0),
    m_useCompositeTet10(useCompositeTet10)
  {}

  virtual ~MasterElementProviderIntrepid2() {}

  void evaluate_field(const Topology& topo,
      const std::vector<double>& paramCoords,     // (numParamCoords)
      const unsigned numFieldComponents,
      const std::vector<double>& fieldData,       // (numFieldComponents x numTopologyNodes)
      std::vector<double>& result) const override;     // (numFieldComponents)

  void evaluate_field(const Topology& topo,
      const unsigned numEvalPoints,
      const std::vector<double>& paramCoords,     // (numParamCoords x numEvalPoints)
      const unsigned numFieldComponents,
      const std::vector<double>& fieldData,       // (numFieldComponents x numTopologyNodes)
      std::vector<double>& result) const override;     // (numFieldComponents x numEvalPoints)


  void nodal_field_data(const EntityKey& key,
      const Field& field,
      unsigned& numFieldComponents,
      unsigned& numNodes,
      std::vector<double>& fieldData) const override;  // (numFieldComponents x numTopologyNodes)

  void nodal_field_data(const std::vector<EntityKey>& nodeKeys,
      const Field& field,
      unsigned& numFieldComponents,
      std::vector<double>& fieldData) const override;  // (numFieldComponents x numTopologyNodes)

  void find_parametric_coordinates(const Topology& topo,
      const unsigned numCoordComponents,
      const std::vector<double>& elemNodeCoords,  // (numCoordComponents x numTopologyNodes)
      const std::vector<double>& inputCoords,     // (numCoordComponents)
      std::vector<double>& paramCoords,
      double& paramDistance) const override;

  void coordinate_center(const Topology& topo, std::vector<double>& coords) const override;

  unsigned num_parametric_coordinates(const Topology& topo) const override;

  unsigned num_integration_points(const Topology& topo) const override;


  void integration_points(const Topology& topo, std::vector<double>& gaussPoints) const override;

  private:
    double compute_parametric_distance(const stk::topology& topo, const std::vector<double>& paramCoords) const;

    using PointViewType = Intrepid2::Cubature<BasisExecSpace>::PointViewTypeAllocatable;
    using WeightViewType = Intrepid2::Cubature<BasisExecSpace>::WeightViewTypeAllocatable;
    using Cubature = Intrepid2::Cubature<BasisExecSpace, double, double>;

    Teuchos::RCP<Cubature> get_default_cubature(const Topology& topo) const;

    void check_consistent_topology(const Topology& topo) const;

    mutable Teuchos::RCP<Intrepid2::Basis<I2DeviceType,double,double>> m_basis;
    mutable stk::topology m_lastTopology;
    mutable BasisViewType m_basisVals;
    mutable PointViewType m_paramCoords;
    mutable unsigned m_numEvalPoints;
    bool m_useCompositeTet10;
};
}

#endif
