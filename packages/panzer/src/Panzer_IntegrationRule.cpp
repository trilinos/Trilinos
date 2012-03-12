
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
  
  // Intrepid does support a quadrature on a 0-dimensional object
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
