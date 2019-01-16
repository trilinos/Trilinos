// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include "Panzer_IntegrationDescriptor.hpp"

#include "Panzer_HashUtils.hpp"

#include "Teuchos_Assert.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Shards_CellTopology.hpp"
#include "Panzer_PointDescriptor.hpp"
#include "Panzer_PointGenerator.hpp"   // includes Kokkos::DynRankView
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_RealSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_CubatureControlVolume.hpp"
#include "Intrepid2_CubatureControlVolumeSide.hpp"
#include "Intrepid2_CubatureControlVolumeBoundary.hpp"


namespace panzer
{

// Anonymous namespace that holds the coordinate generator,
// this hides any details about the generator from an external
// user and hopefully hides Kokkos from any file that doesn't need it.
namespace {

/** A generator that builds integration points
 *
 * \note: This is largely untested, we don't explicitly use this for anything other than volume points
 *
 * \note: The order of side points with be wrong w.r.t. workset points due to the need to align points on surfaces/sides
 *
  */
class IntegrationCoordsGenerator :
    public PointGenerator {
public:

  IntegrationCoordsGenerator() = delete;
  IntegrationCoordsGenerator(const IntegrationCoordsGenerator &) = delete;

  IntegrationCoordsGenerator(const int type,
                             const int order,
                             const int side = -1)
    : type_(type), order_(order), side_(side) {}

  virtual ~IntegrationCoordsGenerator() = default;

  Kokkos::DynRankView<double> getPoints(const shards::CellTopology & topo) const override
  {

    auto intrepid_cubature = getIntrepidCubature(type_, order_, topo, side_);

    const int num_points = intrepid_cubature->getNumPoints();

    Kokkos::DynRankView<double> points("cubature_ref_coords",num_points,topo.getDimension());

    // HACK: Need to implement this, but it isn't clear how (need to sort points per face based on physics position which can't be done without vertices)
    if(type_ == IntegrationDescriptor::SURFACE)
      return points;

    Kokkos::DynRankView<double> weights("cubature_weights",num_points);

    intrepid_cubature->getCubature(points, weights);

    return points;
  }

  virtual bool hasPoints(const shards::CellTopology & topo) const {return type_ != IntegrationDescriptor::SURFACE;}

  int numPoints(const shards::CellTopology & topo) const override
  {
    if(type_ == IntegrationDescriptor::SURFACE)
      return 0;
    return getIntrepidCubature(type_, order_, topo, side_)->getNumPoints();
  }

protected:

  Teuchos::RCP<Intrepid2::Cubature<PHX::Device::execution_space,double,double>>
  getIntrepidCubature(const int type,
                      const int order,
                      const shards::CellTopology & topo,
                      const int side = -1) const
  {
    typedef panzer::IntegrationDescriptor ID;

    if(type == ID::CV_SIDE){
      return Teuchos::rcp(new Intrepid2::CubatureControlVolumeSide<PHX::Device::execution_space,double,double>(topo));
    } else if(type == ID::CV_VOLUME){
      return Teuchos::rcp(new Intrepid2::CubatureControlVolume<PHX::Device::execution_space,double,double>(topo));
    } else if(type == ID::CV_BOUNDARY){
      TEUCHOS_ASSERT(side >= 0);
      return Teuchos::rcp(new Intrepid2::CubatureControlVolumeBoundary<PHX::Device::execution_space,double,double>(topo,side));
    } else {

      // Surface integration is not easily supported by Intrepid2
      // TODO: This can still be done, we just need to do it
      TEUCHOS_ASSERT(type != ID::SURFACE);

      Intrepid2::DefaultCubatureFactory cubature_factory;

      if(side >= 0){
        const auto side_topology = getSideTopology(topo,side);
        // Special case of side integration
        return cubature_factory.create<PHX::Device::execution_space,double,double>(*side_topology,order);

      }
      return cubature_factory.create<PHX::Device::execution_space,double,double>(topo,order);
    }
  }

  Teuchos::RCP<shards::CellTopology>
  getSideTopology(const shards::CellTopology & topo,
                  const int side) const
  {
    TEUCHOS_ASSERT(side >= 0);

    Teuchos::RCP<shards::CellTopology> sideTopo;
    const int spatial_dimension = topo.getDimension();

    if(spatial_dimension == 1)
      return Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData<shards::Node>()));;

    TEUCHOS_TEST_FOR_EXCEPTION( (side >= static_cast<int>(topo.getSideCount())),
                                std::runtime_error, "Error - local side "
                                << side << " is not in range (0->" << topo.getSideCount()-1
                                << ") of topologic entity!");

    return Teuchos::rcp(new shards::CellTopology(topo.getCellTopologyData(spatial_dimension-1,side)));
  }

  int type_;
  int order_;
  int side_;

};

} // end namespace <anonymous>

IntegrationDescriptor::IntegrationDescriptor()
{
  setup(-1, NONE);
}

IntegrationDescriptor::IntegrationDescriptor(const int cubature_order, const int integration_type, const int side)
{
  setup(cubature_order, integration_type, side);
}

void
IntegrationDescriptor::setup(const int cubature_order, const int integration_type, const int side)
{
  _integration_type = integration_type;
  _cubature_order = cubature_order;
  _side = side;

  if(_integration_type == SIDE or _integration_type == CV_BOUNDARY){
    TEUCHOS_ASSERT(side >= 0);
  } else {
    TEUCHOS_ASSERT(side == -1);
  }
  _key = std::hash<IntegrationDescriptor>()(*this);
}


PointDescriptor
IntegrationDescriptor::getPointDescriptor() const
{
  std::stringstream ss;
  ss << "Integration Points: Order " << _cubature_order << ", Type ";
  if(_integration_type == VOLUME)
    ss << "Volume";
  else if(_integration_type == SURFACE)
    ss << "Surface";
  else if(_integration_type == SIDE)
    ss << "Side " << _side;
  else if(_integration_type == CV_VOLUME)
    ss << "Control Volume";
  else if(_integration_type == CV_BOUNDARY)
    ss << "Control Volume Boundary";
  else if(_integration_type == CV_SIDE)
    ss << "Control Volume Side";

  return PointDescriptor(ss.str(), Teuchos::rcp(new IntegrationCoordsGenerator(_integration_type, _cubature_order, _side)));

}

}

std::size_t
std::hash<panzer::IntegrationDescriptor>::operator()(const panzer::IntegrationDescriptor& desc) const
{
  std::size_t seed = 0;

  panzer::hash_combine(seed,desc.getType());
  panzer::hash_combine(seed,desc.getOrder());
  panzer::hash_combine(seed,desc.getSide());

  return seed;
}
