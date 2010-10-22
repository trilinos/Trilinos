
#include "Panzer_IntegrationRule.hpp"

#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_CellData.hpp"

panzer::IntegrationRule::
IntegrationRule(int in_cubature_degree, const panzer::CellData& cell_data) :
  side(-1)
{
  cubature_degree = in_cubature_degree ;
  spatial_dimension = cell_data.baseCellDimension();
  workset_size = cell_data.numCells();
  
  if (spatial_dimension == 3)
    topology = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<8> >()));
  else if (spatial_dimension == 2)
    topology = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<4> >()));
  
  TEST_FOR_EXCEPTION(Teuchos::is_null(topology), std::runtime_error,
		     "Failed to allocate cell topology!");
  
  Intrepid::DefaultCubatureFactory<double,Intrepid::FieldContainer<double> > 
    cubature_factory;

  if (cell_data.isSide()) {
    side = cell_data.side();

    TEST_FOR_EXCEPTION( (side >= static_cast<int>(topology->getSideCount())), 
			std::runtime_error, "Error - local side " 
			<< side << " is not in range (0->" << topology->getSideCount()-1 
			<< ") of topologic entity!");
    
    side_topology = Teuchos::rcp(new shards::CellTopology(topology->getTopology(topology->getDimension()-1,side)));
  }

  TEST_FOR_EXCEPTION(side >= 0 && Teuchos::is_null(side_topology), 
		     std::runtime_error,
		     "Failed to allocate side topology!");
  
  Teuchos::RCP<Intrepid::Cubature<double,Intrepid::FieldContainer<double>  > > 
    intrepid_cubature;

  if (Teuchos::is_null(side_topology))
    intrepid_cubature = cubature_factory.create(*topology, cubature_degree);
  else
    intrepid_cubature = cubature_factory.create(*side_topology, cubature_degree);

  num_points = intrepid_cubature->getNumPoints();

  using Teuchos::rcp;
  using PHX::MDALayout;
  
  dl_scalar = 
    rcp(new MDALayout<Cell,IP>(workset_size,num_points));
  
  dl_vector = 
    rcp(new MDALayout<Cell,IP,Dim>(workset_size, num_points,
				   spatial_dimension));
  
  dl_tensor = 
    rcp(new MDALayout<Cell,IP,Dim,Dim>(workset_size, num_points,
				       spatial_dimension,
				       spatial_dimension));

}

bool panzer::IntegrationRule::isSide()
{
  return (!Teuchos::is_null(side_topology));
}
