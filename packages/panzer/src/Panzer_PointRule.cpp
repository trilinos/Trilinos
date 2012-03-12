
#include "Panzer_PointRule.hpp"

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_CellData.hpp"

panzer::PointRule::
PointRule(const std::string & ptName,
          int np,
          const panzer::CellData& cell_data) :
  side(-1)
{
  setup(ptName,np,cell_data);
}

void panzer::PointRule::
setup(const std::string & ptName,
      int np, 
      const panzer::CellData& cell_data)
{
  point_name = ptName;
  num_points = np;
  spatial_dimension = cell_data.baseCellDimension();
  workset_size = cell_data.numCells();
  
  topology = cell_data.getCellTopology();
  TEUCHOS_TEST_FOR_EXCEPTION(topology==Teuchos::null,std::runtime_error,
                     "PointRule::setup - Base topology from cell_data cannot be null!");
  TEUCHOS_TEST_FOR_EXCEPTION(spatial_dimension!=(int) topology->getDimension(), std::runtime_error,
		     "PointRule::setup - Spatial dimension from cell_data does not match the cell topology.");
  
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(topology), std::runtime_error,
		     "PointRule::setup - Failed to allocate cell topology!");
  
  // handle side issues : first veryify the 0D side follows the rules
  if(cell_data.isSide() && spatial_dimension==1) {
     TEUCHOS_ASSERT(num_points==1); // only one point on a node 
  }

  // now extract side topology
  side_topology = getSideTopology(cell_data);
  if (side_topology!=Teuchos::null)
    side = cell_data.side();

  TEUCHOS_TEST_FOR_EXCEPTION(side >= 0 && Teuchos::is_null(side_topology), 
		     std::runtime_error,
		     "Failed to allocate side topology!");

  // allocate data layout objects
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

const std::string & panzer::PointRule::getName() const
{
   return point_name;
}

bool panzer::PointRule::isSide()
{
  return (!Teuchos::is_null(side_topology));
}

Teuchos::RCP<shards::CellTopology> panzer::PointRule::getSideTopology(const CellData & cell_data)
{
  Teuchos::RCP<shards::CellTopology> sideTopo;
  int spatial_dimension = cell_data.baseCellDimension();
  Teuchos::RCP<const shards::CellTopology> topology = cell_data.getCellTopology();

  if (cell_data.isSide() && spatial_dimension>1) {
    int side = cell_data.side();

    TEUCHOS_TEST_FOR_EXCEPTION( (side >= static_cast<int>(topology->getSideCount())), 
			std::runtime_error, "Error - local side " 
			<< side << " is not in range (0->" << topology->getSideCount()-1 
   			<< ") of topologic entity!");
    
    sideTopo = Teuchos::rcp(new shards::CellTopology(topology->getCellTopologyData(topology->getDimension()-1,side)));
  }
  else if(cell_data.isSide() && spatial_dimension==1) {
    sideTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData<shards::Node>()));
  }
  
  return sideTopo;
}

void panzer::PointRule::print(std::ostream & os)
{
   os << "panzer::PointRule ( "
      << "Name = " << getName()
      << ", Dimension = " << spatial_dimension 
      << ", Workset Size = " << workset_size
      << ", Num Points = " << num_points 
      << ", Side = " << side
      << " )";
}
