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


#include "Panzer_IntegrationRule.hpp"

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_CubatureControlVolume.hpp"
#include "Intrepid2_CubatureControlVolumeSide.hpp"
#include "Intrepid2_CubatureControlVolumeBoundary.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_CellData.hpp"

panzer::IntegrationRule::
IntegrationRule(int in_cubature_degree, const panzer::CellData& cell_data)
   : PointRule() 
{
  setup(in_cubature_degree,cell_data);
}

panzer::IntegrationRule::
IntegrationRule(const panzer::CellData& cell_data, std::string in_cv_type)
   : PointRule()
{
  setup_cv(cell_data,in_cv_type);
}

void panzer::IntegrationRule::setup(int in_cubature_degree, const panzer::CellData& cell_data)
{
  cubature_degree = in_cubature_degree;
  cv_type = "none";
  int spatialDimension = cell_data.baseCellDimension();

  std::stringstream ss;
  ss << "CubaturePoints (Degree=" << cubature_degree;
  
  // Intrepid2 does not support a quadrature on a 0-dimensional object
  // (which doesn't make much sense anyway) to work around this we
  // will adjust the integration rule manually
  if(cell_data.isSide() && spatialDimension==1) {
     ss << ",side)";
     PointRule::setup(ss.str(),1,cell_data);

     return;
  }

  Teuchos::RCP<const shards::CellTopology> topo = cell_data.getCellTopology();
  Teuchos::RCP<shards::CellTopology> sideTopo = getSideTopology(cell_data);

  Intrepid2::DefaultCubatureFactory<double,Kokkos::DynRankView<double,PHX::Device> > cubature_factory;
  Teuchos::RCP<Intrepid2::Cubature<double,Kokkos::DynRankView<double,PHX::Device>  > > intrepid_cubature;

  // get side topology
  if (Teuchos::is_null(sideTopo)) {
    ss << ",volume)";
    intrepid_cubature = cubature_factory.create(*topo, cubature_degree);
  }
  else {
    ss << ",side)";
    intrepid_cubature = cubature_factory.create(*sideTopo, cubature_degree);
  }

  PointRule::setup(ss.str(),intrepid_cubature->getNumPoints(),cell_data);
}

void panzer::IntegrationRule::setup_cv(const panzer::CellData& cell_data, std::string in_cv_type)
{
  // set cubature degree to arbitrary constant for indexing
  cubature_degree = 1;
  cv_type = in_cv_type;
  if (cv_type == "volume") {
     cubature_degree = 75;
  }
  if (cv_type == "side") {
     cubature_degree = 85;
  }
  if (cv_type == "boundary") {
     cubature_degree = 95;
  }

  //int spatialDimension = cell_data.baseCellDimension();

  std::stringstream ss;
  ss << "CubaturePoints ControlVol (Index=" << cubature_degree;

  Teuchos::RCP<const shards::CellTopology> topo = cell_data.getCellTopology();

  Teuchos::RCP<Intrepid2::Cubature<double,Kokkos::DynRankView<double,PHX::Device>  > > intrepid_cubature;

  int num_points(0);
  if (cv_type == "volume") {
    ss << ",volume)";
    intrepid_cubature  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,Kokkos::DynRankView<double,PHX::Device>,Kokkos::DynRankView<double,PHX::Device> >(topo));
    num_points = intrepid_cubature->getNumPoints();
  }
  else if (cv_type == "side") {
    ss << ",side)";
    intrepid_cubature  = Teuchos::rcp(new Intrepid2::CubatureControlVolumeSide<double,Kokkos::DynRankView<double,PHX::Device>,Kokkos::DynRankView<double,PHX::Device> >(topo));
    num_points = intrepid_cubature->getNumPoints();
  }
  else if (cv_type == "boundary") {
    ss << ",boundary)";
    intrepid_cubature  = Teuchos::rcp(new 
           Intrepid2::CubatureControlVolumeBoundary<double,Kokkos::DynRankView<double,PHX::Device>,Kokkos::DynRankView<double,PHX::Device> >(topo,cell_data.side()));
    num_points = intrepid_cubature->getNumPoints();
  }

  PointRule::setup(ss.str(),num_points,cell_data);
}

int panzer::IntegrationRule::order() const
{ return cubature_degree; }

void panzer::IntegrationRule::print(std::ostream & os)
{
   os << "IntegrationRule ( "
      << "Name = " << getName()
      << ", Degree = " << cubature_degree 
      << ", Dimension = " << spatial_dimension 
      << ", Workset Size = " << workset_size
      << ", Num Points = " << num_points 
      << ", Side = " << side
      << " )";
}

void panzer::IntegrationRule::referenceCoordinates(Kokkos::DynRankView<double,PHX::Device> & cub_points)
{
    // build an interpid cubature rule
    Teuchos::RCP< Intrepid2::Cubature<double,Kokkos::DynRankView<double,PHX::Device> > > intrepid_cubature;
    Intrepid2::DefaultCubatureFactory<double,Kokkos::DynRankView<double,PHX::Device> > cubature_factory;
    
    if (!isSide())
      intrepid_cubature = cubature_factory.create(*(topology),cubature_degree);
    else
      intrepid_cubature = cubature_factory.create(*(side_topology),cubature_degree);

    int num_ip = intrepid_cubature->getNumPoints();
    Kokkos::DynRankView<double,PHX::Device> cub_weights("cub_weights",num_ip);

    // now compute weights (and throw them out) as well as reference points
    if (!isSide()) {
      cub_points = Kokkos::DynRankView<double,PHX::Device>("cub_points", num_ip, topology->getDimension());
      intrepid_cubature->getCubature(cub_points, cub_weights);
    }
    else {
      cub_points = Kokkos::DynRankView<double,PHX::Device>("cub_points", num_ip, side_topology->getDimension());
      intrepid_cubature->getCubature(cub_points, cub_weights);
    }
}
