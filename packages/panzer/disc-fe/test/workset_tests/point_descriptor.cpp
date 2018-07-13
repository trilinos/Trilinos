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
