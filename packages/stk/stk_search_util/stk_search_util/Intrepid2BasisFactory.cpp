#include <stk_mesh/base/MetaData.hpp>
#include "Intrepid2BasisFactory.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"



namespace stk::search {

BasisPtr
lagrangeBasisFactory( stk::topology stkTopo, bool useCompositeTet10 )
{
  using outputValueType = double;
  using pointValueType = double;

  shards::CellTopology cellTopo = stk::mesh::get_cell_topology(stkTopo);
  BasisPtr r_val;

  //This essentially mirrors the code in Intrepid2_CellTools.hpp in the
  //function 'createHGradBasis'. The main difference is that we are
  //returning std::shared_ptr whereas Intrepid2 returns Teuchos::RCP.

  switch (cellTopo.getKey()) {
  case shards::Line<2>::key:          r_val = std::make_shared<Intrepid2::Basis_HGRAD_LINE_C1_FEM   <I2DeviceType,outputValueType,pointValueType>>(); break;
  case shards::Line<3>::key:          r_val = std::make_shared<Intrepid2::Basis_HGRAD_LINE_C2_FEM   <I2DeviceType,outputValueType,pointValueType>>(); break;
  case shards::Triangle<3>::key:      r_val = std::make_shared<Intrepid2::Basis_HGRAD_TRI_C1_FEM    <I2DeviceType,outputValueType,pointValueType>>(); break;
  case shards::Quadrilateral<4>::key: r_val = std::make_shared<Intrepid2::Basis_HGRAD_QUAD_C1_FEM   <I2DeviceType,outputValueType,pointValueType>>(); break;
  case shards::Tetrahedron<4>::key:   r_val = std::make_shared<Intrepid2::Basis_HGRAD_TET_C1_FEM    <I2DeviceType,outputValueType,pointValueType>>(); break;
  case shards::Hexahedron<8>::key:    r_val = std::make_shared<Intrepid2::Basis_HGRAD_HEX_C1_FEM    <I2DeviceType,outputValueType,pointValueType>>(); break;
  case shards::Wedge<6>::key:         r_val = std::make_shared<Intrepid2::Basis_HGRAD_WEDGE_C1_FEM  <I2DeviceType,outputValueType,pointValueType>>(); break;
  case shards::Pyramid<5>::key:       r_val = std::make_shared<Intrepid2::Basis_HGRAD_PYR_C1_FEM    <I2DeviceType,outputValueType,pointValueType>>(); break;

  case shards::Triangle<6>::key:      r_val = std::make_shared<Intrepid2::Basis_HGRAD_TRI_C2_FEM    <I2DeviceType,outputValueType,pointValueType>>(); break;
  case shards::Quadrilateral<8>::key: r_val = std::make_shared<Intrepid2::Basis_HGRAD_QUAD_I2_FEM   <I2DeviceType,outputValueType,pointValueType>>(); break;
  case shards::Quadrilateral<9>::key: r_val = std::make_shared<Intrepid2::Basis_HGRAD_QUAD_C2_FEM   <I2DeviceType,outputValueType,pointValueType>>(); break;
  case shards::Tetrahedron<10>::key:
    if (useCompositeTet10) {
      r_val = std::make_shared<Intrepid2::Basis_HGRAD_TET_COMP12_FEM<I2DeviceType,outputValueType,pointValueType>>();
    }
    else {
      r_val = std::make_shared<Intrepid2::Basis_HGRAD_TET_C2_FEM<I2DeviceType,outputValueType,pointValueType>>();
    }
    break;
  case shards::Tetrahedron<11>::key:  r_val = std::make_shared<Intrepid2::Basis_HGRAD_TET_COMP12_FEM<I2DeviceType,outputValueType,pointValueType>>(); break;
  case shards::Hexahedron<20>::key:   r_val = std::make_shared<Intrepid2::Basis_HGRAD_HEX_I2_FEM    <I2DeviceType,outputValueType,pointValueType>>(); break;
  case shards::Hexahedron<27>::key:   r_val = std::make_shared<Intrepid2::Basis_HGRAD_HEX_C2_FEM    <I2DeviceType,outputValueType,pointValueType>>(); break;
  case shards::Wedge<15>::key:        r_val = std::make_shared<Intrepid2::Basis_HGRAD_WEDGE_I2_FEM  <I2DeviceType,outputValueType,pointValueType>>(); break;
  case shards::Wedge<18>::key:        r_val = std::make_shared<Intrepid2::Basis_HGRAD_WEDGE_C2_FEM  <I2DeviceType,outputValueType,pointValueType>>(); break;
  case shards::Pyramid<13>::key:      r_val = std::make_shared<Intrepid2::Basis_HGRAD_PYR_I2_FEM    <I2DeviceType,outputValueType,pointValueType>>(); break;

  case shards::ShellQuadrilateral<4>::key: r_val = std::make_shared<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<I2DeviceType,outputValueType,pointValueType>>(); break;
  case shards::ShellTriangle<3>::key: r_val = std::make_shared<Intrepid2::Basis_HGRAD_TRI_C1_FEM<I2DeviceType,outputValueType,pointValueType>>(); break;

  case shards::Beam<2>::key:
  case shards::Beam<3>::key:
  case shards::ShellLine<2>::key:
  case shards::ShellLine<3>::key:
  case shards::ShellTriangle<6>::key:
  case shards::ShellQuadrilateral<8>::key:
  case shards::ShellQuadrilateral<9>::key:
  default:
    STK_ThrowRequireMsg(false, "Cell topology " + stkTopo.name() + " not supported.");
  }
  return r_val;
}


}
