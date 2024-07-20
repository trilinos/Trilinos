// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_IntegrationRule.hpp"

#include "Teuchos_ArrayRCP.hpp"
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
   : PointRule(), IntegrationDescriptor()
{
  if(cell_data.isSide()){
    IntegrationDescriptor::setup(in_cubature_degree, IntegrationDescriptor::SIDE, cell_data.side());
  } else {
    IntegrationDescriptor::setup(in_cubature_degree, IntegrationDescriptor::VOLUME);
  }
  setup(in_cubature_degree,cell_data);
}

panzer::IntegrationRule::
IntegrationRule(const panzer::CellData& cell_data, const std::string & in_cv_type)
   : PointRule(), IntegrationDescriptor()
{
  // Cubature orders are only used for indexing so we make them large enough not to interfere with other rules.
  // TODO: This requirement (on arbitrary cubature order) will be dropped with the new Workset design (using descriptors to find integration rules)
  if(in_cv_type == "volume"){
    TEUCHOS_TEST_FOR_EXCEPT_MSG(cell_data.isSide(),
                                "IntegrationRule::IntegrationRule : Control Volume 'volume' type requested, but CellData is setup for sides.");
    IntegrationDescriptor::setup(75, IntegrationDescriptor::CV_VOLUME);
  } else if(in_cv_type == "side"){
    IntegrationDescriptor::setup(85, IntegrationDescriptor::CV_SIDE);
  } else if(in_cv_type == "boundary"){
    TEUCHOS_TEST_FOR_EXCEPT_MSG(not cell_data.isSide(),
                                "IntegrationRule::IntegrationRule : Control Volume 'boundary' type requested, but CellData is not setup for sides.");
    IntegrationDescriptor::setup(95, IntegrationDescriptor::CV_BOUNDARY, cell_data.side());
  } else {
    TEUCHOS_ASSERT(false);
  }
  setup_cv(cell_data,in_cv_type);
}

panzer::IntegrationRule::
IntegrationRule(const panzer::IntegrationDescriptor& description,
                const Teuchos::RCP<const shards::CellTopology> & cell_topology,
                const int num_cells,
                const int num_faces)
  : PointRule(), IntegrationDescriptor(description)
{

  TEUCHOS_ASSERT(description.getType() != panzer::IntegrationDescriptor::NONE);

  cubature_degree = description.getOrder();
  cv_type = "none";

  panzer::CellData cell_data;
  if(isSide()){
    cell_data = panzer::CellData(num_cells, getSide(), cell_topology);
  } else {
    cell_data = panzer::CellData(num_cells, cell_topology);
  }

  if(getType() == panzer::IntegrationDescriptor::VOLUME){
    setup(getOrder(), cell_data);
  } else if(description.getType() == panzer::IntegrationDescriptor::SIDE){
    setup(getOrder(), cell_data);
  } else if(description.getType() == panzer::IntegrationDescriptor::SURFACE){
    TEUCHOS_ASSERT(num_faces!=-1);
    setup_surface(cell_topology, num_cells, num_faces);
  } else if(description.getType() == panzer::IntegrationDescriptor::CV_VOLUME){
    setup_cv(cell_data, "volume");
  } else if(description.getType() == panzer::IntegrationDescriptor::CV_SIDE){
    setup_cv(cell_data, "side");
  } else if(description.getType() == panzer::IntegrationDescriptor::CV_BOUNDARY){
    setup_cv(cell_data, "boundary");
  } else {
    TEUCHOS_ASSERT(false);
  }

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

  const shards::CellTopology & topo = *cell_data.getCellTopology();
  Teuchos::RCP<shards::CellTopology> sideTopo = getSideTopology(cell_data);

  Intrepid2::DefaultCubatureFactory cubature_factory;
  Teuchos::RCP<Intrepid2::Cubature<PHX::Device::execution_space,double,double>> intrepid_cubature;

  // get side topology
  if (Teuchos::is_null(sideTopo)) {
    ss << ",volume)";
    intrepid_cubature = cubature_factory.create<PHX::Device::execution_space,double,double>(topo, cubature_degree);
  }
  else {
    ss << ",side)";
    intrepid_cubature = cubature_factory.create<PHX::Device::execution_space,double,double>(*sideTopo, cubature_degree);
  }

  PointRule::setup(ss.str(),intrepid_cubature->getNumPoints(),cell_data);
}

void panzer::IntegrationRule::setup_surface(const Teuchos::RCP<const shards::CellTopology> & cell_topology, const int num_cells, const int num_faces)
{

  const int cell_dim = cell_topology->getDimension();
  const int subcell_dim = cell_dim-1;
  const int num_faces_per_cell = cell_topology->getSubcellCount(subcell_dim);

  panzer::CellData cell_data(num_cells, cell_topology);

  std::string point_rule_name;
  {
    std::stringstream ss;
    ss << "CubaturePoints (Degree=" << getOrder() << ",surface)";
    point_rule_name = ss.str();
  }

  // We can skip some steps for 1D
  if(cell_dim == 1){
    const int num_points_per_cell = num_faces_per_cell;
    const int num_points_per_face = 1;
    PointRule::setup(point_rule_name, num_cells, num_points_per_cell, num_faces, num_points_per_face, cell_topology);
    _point_offsets.resize(3,0);
    _point_offsets[0] = 0;
    _point_offsets[1] = num_points_per_face;
    _point_offsets[2] = _point_offsets[1]+num_points_per_face;
    return;
  }

  Intrepid2::DefaultCubatureFactory cubature_factory;

  _point_offsets.resize(num_faces_per_cell+1,0);
  int test_face_size = -1;
  for(int subcell_index=0; subcell_index<num_faces_per_cell; ++subcell_index){
    Teuchos::RCP<shards::CellTopology> face_topology = Teuchos::rcp(new shards::CellTopology(cell_topology->getCellTopologyData(subcell_dim,subcell_index)));
    const auto & intrepid_cubature = cubature_factory.create<PHX::Device::execution_space,double,double>(*face_topology, getOrder());
    const int num_face_points = intrepid_cubature->getNumPoints();
    _point_offsets[subcell_index+1] = _point_offsets[subcell_index] + num_face_points;

    // Right now we only support each face having the same number of points
    if(test_face_size==-1){
      test_face_size = num_face_points;
    } else {
      TEUCHOS_ASSERT(num_face_points == test_face_size);
    }
  }

  const int num_points_per_cell = _point_offsets.back();
  const int num_points_per_face = _point_offsets[1];

  PointRule::setup(point_rule_name, num_cells, num_points_per_cell, num_faces, num_points_per_face, cell_topology);

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

  const shards::CellTopology & topo = *cell_data.getCellTopology();

  Teuchos::RCP<Intrepid2::Cubature<PHX::Device::execution_space,double,double> > intrepid_cubature;

  int tmp_num_points = 0;
  if (cv_type == "volume") {
    ss << ",volume)";
    intrepid_cubature  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<PHX::Device::execution_space,double,double>(topo));
    tmp_num_points = intrepid_cubature->getNumPoints();
  }
  else if (cv_type == "side") {
    ss << ",side)";
    intrepid_cubature  = Teuchos::rcp(new Intrepid2::CubatureControlVolumeSide<PHX::Device::execution_space,double,double>(topo));
    tmp_num_points = intrepid_cubature->getNumPoints();
  }
  else if (cv_type == "boundary") {
    ss << ",boundary)";
    intrepid_cubature  = Teuchos::rcp(new Intrepid2::CubatureControlVolumeBoundary<PHX::Device::execution_space,double,double>(topo,cell_data.side()));
    tmp_num_points = intrepid_cubature->getNumPoints();
  }

  PointRule::setup(ss.str(),tmp_num_points,cell_data);
}

int panzer::IntegrationRule::order() const
{ return cubature_degree; }


int panzer::IntegrationRule::getPointOffset(const int subcell_index) const
{
  // Need to make sure this is a surface integrator
  TEUCHOS_ASSERT(getType() == SURFACE);
  return _point_offsets[subcell_index];
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

void panzer::IntegrationRule::referenceCoordinates(Kokkos::DynRankView<double,PHX::Device> & cub_points)
{
    // build an interpid cubature rule
  Teuchos::RCP< Intrepid2::Cubature<PHX::Device::execution_space,double,double> > intrepid_cubature;
    Intrepid2::DefaultCubatureFactory cubature_factory;
    
    if (!isSide())
      intrepid_cubature = cubature_factory.create<PHX::Device::execution_space,double,double>(*(topology),cubature_degree);
    else
      intrepid_cubature = cubature_factory.create<PHX::Device::execution_space,double,double>(*(side_topology),cubature_degree);

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
