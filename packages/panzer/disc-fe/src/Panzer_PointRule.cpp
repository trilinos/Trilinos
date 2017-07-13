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


#include "Panzer_PointRule.hpp"

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"
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

panzer::PointRule::
PointRule(const std::string & point_rule_name,
          const int num_cells,
          const int num_points_per_cell,
          const int num_faces,
          const int num_points_per_face,
          const Teuchos::RCP<const shards::CellTopology> & cell_topology)
{
  setup(point_rule_name, num_cells, num_points_per_cell, num_faces, num_points_per_face, cell_topology);
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
  _num_faces = -1;
  _num_points_per_face = -1;
  
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

  dl_vector3 =
      rcp(new MDALayout<Cell,IP,Dim>(workset_size, num_points,3));
  dl_tensor3x3 =
      rcp(new MDALayout<Cell,IP,Dim,Dim>(workset_size, num_points,3,3));

}


void panzer::PointRule::
setup(const std::string & point_rule_name,
      const int num_cells,
      const int num_points_per_cell,
      const int num_faces,
      const int num_points_per_face,
      const Teuchos::RCP<const shards::CellTopology> & cell_topology)
{

  topology = cell_topology;
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(topology),std::runtime_error,
                       "PointRule::setup - Cell topology cannot be null.");

  point_name = point_rule_name;
  num_points = num_points_per_cell;
  spatial_dimension = cell_topology->getDimension();
  workset_size = num_cells;
  _num_faces = num_faces;
  _num_points_per_face = num_points_per_face;
  side = -1;

  // allocate data layout objects
  dl_scalar = Teuchos::rcp(new PHX::MDALayout<Cell,IP>(workset_size,num_points));
  dl_vector = Teuchos::rcp(new PHX::MDALayout<Cell,IP,Dim>(workset_size, num_points,spatial_dimension));
  dl_tensor = Teuchos::rcp(new PHX::MDALayout<Cell,IP,Dim,Dim>(workset_size, num_points,spatial_dimension,spatial_dimension));

  dl_vector3 = Teuchos::rcp(new PHX::MDALayout<Cell,IP,Dim>(workset_size, num_points,3));
  dl_tensor3x3 = Teuchos::rcp(new PHX::MDALayout<Cell,IP,Dim,Dim>(workset_size, num_points,3,3));

}

const std::string & panzer::PointRule::getName() const
{
   return point_name;
}

bool panzer::PointRule::isSide() const
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

Teuchos::RCP<PHX::DataLayout>
panzer::PointRule::getCellDataLayout() const
{return Teuchos::rcp(new PHX::MDALayout<Cell>(workset_size));}

Teuchos::RCP<PHX::DataLayout>
panzer::PointRule::getCellDataLayout(const int dim0) const
{return Teuchos::rcp(new PHX::MDALayout<Cell,Dim>(workset_size, dim0));}

Teuchos::RCP<PHX::DataLayout>
panzer::PointRule::getCellDataLayout(const int dim0, const int dim1) const
{return Teuchos::rcp(new PHX::MDALayout<Cell,Dim,Dim>(workset_size, dim0, dim1));}


Teuchos::RCP<PHX::DataLayout>
panzer::PointRule::getCellPointDataLayout() const
{return Teuchos::rcp(new PHX::MDALayout<Cell,IP>(workset_size, num_points));}

Teuchos::RCP<PHX::DataLayout>
panzer::PointRule::getCellPointDataLayout(const int dim0) const
{return Teuchos::rcp(new PHX::MDALayout<Cell,IP,Dim>(workset_size, num_points, dim0));}

Teuchos::RCP<PHX::DataLayout>
panzer::PointRule::getCellPointDataLayout(const int dim0, const int dim1) const
{return Teuchos::rcp(new PHX::MDALayout<Cell,IP,Dim,Dim>(workset_size, num_points, dim0, dim1));}


Teuchos::RCP<PHX::DataLayout>
panzer::PointRule::getFaceDataLayout() const
{return Teuchos::rcp(new PHX::MDALayout<Face>(_num_faces));}

Teuchos::RCP<PHX::DataLayout>
panzer::PointRule::getFaceDataLayout(const int dim0) const
{return Teuchos::rcp(new PHX::MDALayout<Face,Dim>(_num_faces, dim0));}

Teuchos::RCP<PHX::DataLayout>
panzer::PointRule::getFaceDataLayout(const int dim0, const int dim1) const
{return Teuchos::rcp(new PHX::MDALayout<Face,Dim,Dim>(_num_faces, dim0, dim1));}


Teuchos::RCP<PHX::DataLayout>
panzer::PointRule::getFacePointDataLayout() const
{return Teuchos::rcp(new PHX::MDALayout<Face,IP>(_num_faces, _num_points_per_face));}

Teuchos::RCP<PHX::DataLayout>
panzer::PointRule::getFacePointDataLayout(const int dim0) const
{return Teuchos::rcp(new PHX::MDALayout<Face,IP,Dim>(_num_faces, _num_points_per_face, dim0));}

Teuchos::RCP<PHX::DataLayout>
panzer::PointRule::getFacePointDataLayout(const int dim0, const int dim1) const
{return Teuchos::rcp(new PHX::MDALayout<Face,IP,Dim,Dim>(_num_faces, _num_points_per_face, dim0, dim1));}



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
