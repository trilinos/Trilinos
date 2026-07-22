// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include "Shards_CellTopology.hpp"
#include "Shards_BasicTopologies.hpp"

#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Panzer_PointDescriptor.hpp"
#include "Panzer_BasisDescriptor.hpp"
#include "Panzer_PointGenerator.hpp"

#include "Panzer_IntrepidBasisFactory.hpp"

namespace panzer {

  TEUCHOS_UNIT_TEST(point_descriptor, construction)
  {
    TEST_THROW(PointDescriptor("none",Teuchos::null),std::logic_error);

    BasisDescriptor desc_a(4,"HDiv");
    BasisDescriptor desc_b(2,"HGrad");
    PointDescriptor pd_a = desc_a.getPointDescriptor();
    PointDescriptor pd_b = desc_b.getPointDescriptor();

    TEST_ASSERT(pd_a.getType()!=pd_b.getType());
    TEST_ASSERT(pd_a.getKey()!=pd_b.getKey());

    const PointGenerator & gen_a = pd_a.getGenerator();
    const PointGenerator & gen_b = pd_b.getGenerator();

    shards::CellTopology myHex(shards::getCellTopologyData<shards::Hexahedron<8> >());
    shards::CellTopology myTet(shards::getCellTopologyData<shards::Tetrahedron<4> >());
   
    // check gen_a
    {
      Teuchos::RCP<Intrepid2::Basis<PHX::Device::execution_space,double,double> > 
          intrepid_basis = createIntrepid2Basis<PHX::Device::execution_space,double,double>("HDiv",4,myHex);

      auto hex_pts = gen_a.getPoints(myHex);

      TEST_EQUALITY(Teuchos::as<int>(hex_pts.extent(0)),intrepid_basis->getCardinality());
      TEST_EQUALITY(Teuchos::as<int>(hex_pts.extent(1)),3);
    }

    // check gen_b
    {
      Teuchos::RCP<Intrepid2::Basis<PHX::Device::execution_space,double,double> > 
          intrepid_basis_hex = createIntrepid2Basis<PHX::Device::execution_space,double,double>("HGrad",2,myHex);
      Teuchos::RCP<Intrepid2::Basis<PHX::Device::execution_space,double,double> > 
          intrepid_basis_tet = createIntrepid2Basis<PHX::Device::execution_space,double,double>("HGrad",2,myTet);

      auto tet_pts = gen_b.getPoints(myTet);
      auto hex_pts = gen_b.getPoints(myHex);

      TEST_EQUALITY(Teuchos::as<int>(hex_pts.extent(0)),intrepid_basis_hex->getCardinality());
      TEST_EQUALITY(hex_pts.extent(1),3);

      TEST_EQUALITY(Teuchos::as<int>(tet_pts.extent(0)),intrepid_basis_tet->getCardinality());
      TEST_EQUALITY(tet_pts.extent(1),3);
    }

  }

}
