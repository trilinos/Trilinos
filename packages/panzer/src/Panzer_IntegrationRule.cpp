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
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_CellData.hpp"

panzer::IntegrationRule::
IntegrationRule(int in_cubature_degree, const panzer::CellData& cell_data)
   : PointRule() 
{
  setup(in_cubature_degree,cell_data);
}

void panzer::IntegrationRule::setup(int in_cubature_degree, const panzer::CellData& cell_data)
{
  cubature_degree = in_cubature_degree ;
  int spatialDimension = cell_data.baseCellDimension();

  std::stringstream ss;
  ss << "CubaturePoints (Degree=" << cubature_degree;
  
  // Intrepid does not support a quadrature on a 0-dimensional object
  // (which doesn't make much sense anyway) to work around this we
  // will adjust the integration rule manually
  if(cell_data.isSide() && spatialDimension==1) {
     ss << ",side)";
     PointRule::setup(ss.str(),1,cell_data);

     return;
  }

  Teuchos::RCP<const shards::CellTopology> topo = cell_data.getCellTopology();
  Teuchos::RCP<shards::CellTopology> sideTopo = getSideTopology(cell_data);

  Intrepid::DefaultCubatureFactory<double,Intrepid::FieldContainer<double> > cubature_factory;
  Teuchos::RCP<Intrepid::Cubature<double,Intrepid::FieldContainer<double>  > > intrepid_cubature;

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
