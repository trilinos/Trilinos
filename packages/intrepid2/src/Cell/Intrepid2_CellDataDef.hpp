// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CellDataDef.hpp
    \brief  Definition file for the classes:
    Intrepid2::RefSubcellParametrization,
    Intrepid2::RefCellNodes,
    Intrepid2::RefCellCenter.
    \author Kyungjoo Kim
    \author Mauro Perego
 */

#ifndef __INTREPID2_CELLDATA_DEF_HPP__
#define __INTREPID2_CELLDATA_DEF_HPP__

namespace Intrepid2 {

template<typename DeviceType>
inline bool
RefSubcellParametrization<DeviceType>::
isSupported( const unsigned cellTopoKey ) {
  switch ( cellTopoKey ) {
  case shards::Line<2>::key:
  case shards::Line<3>::key:
  case shards::ShellLine<2>::key:
  case shards::ShellLine<3>::key:
  case shards::Beam<2>::key:
  case shards::Beam<3>::key:
  case shards::Triangle<3>::key:
  case shards::Triangle<4>::key:
  case shards::Triangle<6>::key:
  // case shards::ShellTriangle<3>::key:
  // case shards::ShellTriangle<6>::key:
  case shards::Quadrilateral<4>::key:
  case shards::Quadrilateral<8>::key:
  case shards::Quadrilateral<9>::key:
  // case shards::ShellQuadrilateral<4>::key:
  // case shards::ShellQuadrilateral<8>::key:
  // case shards::ShellQuadrilateral<9>::key:
  case shards::Tetrahedron<4>::key:
  case shards::Tetrahedron<8>::key:
  case shards::Tetrahedron<10>::key:
  case shards::Tetrahedron<11>::key:
  case shards::Hexahedron<8>::key:
  case shards::Hexahedron<20>::key:
  case shards::Hexahedron<27>::key:
  case shards::Pyramid<5>::key:
  case shards::Pyramid<13>::key:
  case shards::Pyramid<14>::key:
  case shards::Wedge<6>::key:
  case shards::Wedge<15>::key:
  case shards::Wedge<18>::key:
  return true;
  default:
    return false;
  }
}

template<typename DeviceType>
inline
typename RefSubcellParametrization<DeviceType>::ConstViewType
RefSubcellParametrization<DeviceType>::
get( const ordinal_type          subcellDim,
    const unsigned              parentCellKey ) {

  if(!isSubcellParametrizationSet_)
    set();

  ViewType subcellParam;

  switch (parentCellKey ) {
  case shards::Tetrahedron<4>::key:
  case shards::Tetrahedron<8>::key:
  case shards::Tetrahedron<10>::key:
  case shards::Tetrahedron<11>::key:       subcellParam = ( subcellDim == 2 ? tetFacesParam : tetEdgesParam ); break;

  case shards::Hexahedron<8>::key:
  case shards::Hexahedron<20>::key:
  case shards::Hexahedron<27>::key:        subcellParam = ( subcellDim == 2 ? hexFacesParam : hexEdgesParam ); break;

  case shards::Pyramid<5>::key:
  case shards::Pyramid<13>::key:
  case shards::Pyramid<14>::key:           subcellParam = ( subcellDim == 2 ? pyrFacesParam : pyrEdgesParam ); break;

  case shards::Wedge<6>::key:
  case shards::Wedge<15>::key:
  case shards::Wedge<18>::key:             subcellParam = ( subcellDim == 2 ? wedgeFacesParam : wedgeEdgesParam ); break;

  case shards::Triangle<3>::key:
  case shards::Triangle<4>::key:
  case shards::Triangle<6>::key:           subcellParam = triEdgesParam; break;

  case shards::Quadrilateral<4>::key:
  case shards::Quadrilateral<8>::key:
  case shards::Quadrilateral<9>::key:      subcellParam = quadEdgesParam; break;

  // case shards::ShellTriangle<3>::key:
  // case shards::ShellTriangle<6>::key:      subcellParam = ( subcellDim == 2 ? shellTriFacesParam : shellTriEdgesParam ); break;

  // case shards::ShellQuadrilateral<4>::key:
  // case shards::ShellQuadrilateral<8>::key:
  // case shards::ShellQuadrilateral<9>::key: subcellParam = ( subcellDim == 2 ? shellQuadFacesParam : shellQuadEdgesParam ); break;

  case shards::ShellLine<2>::key:
  case shards::ShellLine<3>::key:
  case shards::Beam<2>::key:
  case shards::Beam<3>::key:               subcellParam = lineEdgesParam; break;
  default: {
    INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
        ">>> ERROR (Intrepid2::RefSubcellParametrization::get): invalid cell topology.");
  }
  }
  return subcellParam;
}

template<typename DeviceType>
void
RefSubcellParametrization<DeviceType>::set() {

  if(isSubcellParametrizationSet_)
    return;

  ordinal_type subcellDim;
  {
    const auto tet = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4> >());

    subcellDim = 2;
    tetFacesParam = ViewType("CellTools::SubcellParametrization::tetFaces", tet.getSubcellCount(subcellDim), tet.getDimension(), subcellDim+1);
    auto subcell2dParamHost = Kokkos::create_mirror_view(tetFacesParam);
    set( subcell2dParamHost, subcellDim, tet );
    deep_copy(tetFacesParam,subcell2dParamHost);

    subcellDim = 1;
    tetEdgesParam = ViewType("CellTools::SubcellParametrization::tetEdges", tet.getSubcellCount(subcellDim), tet.getDimension(), subcellDim+1);
    auto subcellParamHost = Kokkos::create_mirror_view(tetEdgesParam);
    set( subcellParamHost, subcellDim, tet );
    deep_copy(tetEdgesParam,subcellParamHost);
  }
  {
    const auto hex = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >());

    subcellDim = 2;
    hexFacesParam = ViewType("CellTools::SubcellParametrization::hexFaces", hex.getSubcellCount(subcellDim), hex.getDimension(), subcellDim+1);
    auto subcell2dParamHost = Kokkos::create_mirror_view(hexFacesParam);
    set( subcell2dParamHost, subcellDim, hex );
    deep_copy(hexFacesParam,subcell2dParamHost);

    subcellDim = 1;
    hexEdgesParam = ViewType("CellTools::SubcellParametrization::hexEdges", hex.getSubcellCount(subcellDim), hex.getDimension(), subcellDim+1);
    auto subcellParamHost = Kokkos::create_mirror_view(hexEdgesParam);
    set( subcellParamHost, subcellDim, hex );
    deep_copy(hexEdgesParam,subcellParamHost);
  }
  {
    const auto pyr = shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<5> >());

    subcellDim = 2;
    pyrFacesParam = ViewType("CellTools::SubcellParametrization::pyrFaces", pyr.getSubcellCount(subcellDim), pyr.getDimension(), subcellDim+1);
    auto subcell2dParamHost = Kokkos::create_mirror_view(pyrFacesParam);
    set( subcell2dParamHost, subcellDim, pyr );
    deep_copy(pyrFacesParam,subcell2dParamHost);

    subcellDim = 1;
    pyrEdgesParam = ViewType("CellTools::SubcellParametrization::pyrEdges", pyr.getSubcellCount(subcellDim), pyr.getDimension(), subcellDim+1);
    auto subcellParamHost = Kokkos::create_mirror_view(pyrEdgesParam);
    set( subcellParamHost, subcellDim, pyr );
    deep_copy(pyrEdgesParam,subcellParamHost);
  }
  {
    const auto wedge = shards::CellTopology(shards::getCellTopologyData<shards::Wedge<6> >());

    subcellDim = 2;
    wedgeFacesParam = ViewType("CellTools::SubcellParametrization::wedgeFaces", wedge.getSubcellCount(subcellDim), wedge.getDimension(), subcellDim+1);
    auto subcell2dParamHost = Kokkos::create_mirror_view(wedgeFacesParam);
    set( subcell2dParamHost, subcellDim, wedge );
    deep_copy(wedgeFacesParam,subcell2dParamHost);

    subcellDim = 1;
    wedgeEdgesParam = ViewType("CellTools::SubcellParametrization::wedgeEdges", wedge.getSubcellCount(subcellDim), wedge.getDimension(), subcellDim+1);
    auto subcellParamHost = Kokkos::create_mirror_view(wedgeEdgesParam);
    set( subcellParamHost, subcellDim, wedge );
    deep_copy(wedgeEdgesParam,subcellParamHost);
  }
  {
    const auto tri = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >());

    subcellDim = 1;
    triEdgesParam = ViewType("CellTools::SubcellParametrization::triEdges", tri.getSubcellCount(subcellDim), tri.getDimension(), subcellDim+1);
    auto subcellParamHost = Kokkos::create_mirror_view(triEdgesParam);
    set( subcellParamHost, subcellDim, tri );
    deep_copy(triEdgesParam,subcellParamHost);
  }
  {
    const auto quad = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >());

    subcellDim = 1;
    quadEdgesParam = ViewType("CellTools::SubcellParametrization::quadEdges", quad.getSubcellCount(subcellDim), quad.getDimension(), subcellDim+1);
    auto subcellParamHost = Kokkos::create_mirror_view(quadEdgesParam);
    set( subcellParamHost, subcellDim, quad );
    deep_copy(quadEdgesParam,subcellParamHost);

  }
  {
    const auto line = shards::CellTopology(shards::getCellTopologyData<shards::ShellLine<2> >());

    subcellDim = 1;
    lineEdgesParam = ViewType("CellTools::SubcellParametrization::lineEdges", line.getSubcellCount(subcellDim), line.getDimension(), subcellDim+1);
    auto subcellParamHost = Kokkos::create_mirror_view(lineEdgesParam);
    set( subcellParamHost, subcellDim, line );
    deep_copy(lineEdgesParam,subcellParamHost);
  }

  Kokkos::push_finalize_hook( [=] {
    lineEdgesParam = ViewType();
    triEdgesParam = ViewType();
    quadEdgesParam = ViewType();
    shellTriEdgesParam = ViewType();
    shellQuadEdgesParam = ViewType();
    tetEdgesParam = ViewType();
    hexEdgesParam = ViewType();
    pyrEdgesParam = ViewType();
    wedgeEdgesParam = ViewType();
    shellTriFacesParam = ViewType();
    shellQuadFacesParam = ViewType();
    tetFacesParam = ViewType();
    hexFacesParam = ViewType();
    pyrFacesParam = ViewType();
    wedgeFacesParam = ViewType();
  });

  isSubcellParametrizationSet_= true;
}

template<typename DeviceType>
template <typename HostViewType>
void
RefSubcellParametrization<DeviceType>::
set( HostViewType               subcellParam,
    const ordinal_type         subcellDim,
    const shards::CellTopology parentCell ) {
  // subcellParametrization is rank-3 FieldContainer with dimensions (SC, PCD, COEF) where:
  //  - SC    is the subcell count of subcells with the specified dimension in the parent cell
  //  - PCD   is Parent Cell Dimension, which gives the number of coordinate functions in the map
  //          PCD = 2 for standard 2D cells and non-standard 2D cells: shell line and beam
  //          PCD = 3 for standard 3D cells and non-standard 3D cells: shell Tri and Quad
  //  - COEF  is number of coefficients needed to specify a coordinate function:
  //          COEFF = 2 for edge parametrizations
  //          COEFF = 3 for both Quad and Tri face parametrizations. Because all Quad reference faces
  //          are affine, the coefficient of the bilinear term u*v is zero and is not stored, i.e.,
  //          3 coefficients are sufficient to store Quad face parameterization maps.
  //
  // Edge parametrization maps [-1,1] to edge defined by (v0, v1)
  // Face parametrization maps [-1,1]^2 to quadrilateral face (v0, v1, v2, v3), or
  // standard 2-simplex  {(0,0),(1,0),(0,1)} to traingle face (v0, v1, v2).
  // This defines orientation-preserving parametrizations with respect to reference edge and
  // face orientations induced by their vertex order.

  // get subcellParametrization dimensions: (sc, pcd, coeff)
  const auto sc    = parentCell.getSubcellCount(subcellDim);
  const auto pcd   = parentCell.getDimension();

  INTREPID2_TEST_FOR_EXCEPTION( subcellDim < 1 || subcellDim > static_cast<ordinal_type>(pcd-1), std::invalid_argument,
      ">>> ERROR (Intrepid2::RefSubcellParametrization::set): Parametrizations defined in a range between 1 and (dim-1)");

  const auto refNodes = RefCellNodes<Kokkos::HostSpace>::get(parentCell.getKey());

  if (subcellDim == 1) {
    // Edge parametrizations of 2D and 3D cells (shell lines and beams are 2D cells with edges)
    for (size_type subcellOrd=0;subcellOrd<sc;++subcellOrd) {
      // vertexK[0] = x_k; vertexK[1] = y_k; vertexK[2] = z_k; z_k = 0 for 2D cells
      // Note that ShellLine and Beam are 2D cells!
      const auto v0ord = parentCell.getNodeMap(subcellDim, subcellOrd, 0);
      const auto v1ord = parentCell.getNodeMap(subcellDim, subcellOrd, 1);

      const auto v0 = Kokkos::subview(refNodes, v0ord, Kokkos::ALL());
      const auto v1 = Kokkos::subview(refNodes, v1ord, Kokkos::ALL());

      // x(t) = (x0 + x1)/2 + t*(x1 - x0)/2
      subcellParam(subcellOrd, 0, 0) = (v0(0) + v1(0))/2.0;
      subcellParam(subcellOrd, 0, 1) = (v1(0) - v0(0))/2.0;

      // y(t) = (y0 + y1)/2 + t*(y1 - y0)/2
      subcellParam(subcellOrd, 1, 0) = (v0(1) + v1(1))/2.0;
      subcellParam(subcellOrd, 1, 1) = (v1(1) - v0(1))/2.0;

      if( pcd == 3 ) {
        // z(t) = (z0 + z1)/2 + t*(z1 - z0)/2
        subcellParam(subcellOrd, 2, 0) = (v0(2) + v1(2))/2.0;
        subcellParam(subcellOrd, 2, 1) = (v1(2) - v0(2))/2.0;
      }
    }
  }
  else if (subcellDim == 2) {
    // Face parametrizations of 3D cells: (shell Tri and Quad are 3D cells with faces)
    // A 3D cell can have both Tri and Quad faces, but because they are affine images of the
    // parametrization domain, 3 coefficients are enough to store them in both cases.
    for (size_type subcellOrd=0;subcellOrd<sc;++subcellOrd) {

      switch (parentCell.getKey(subcellDim,subcellOrd)) {

      case shards::Triangle<3>::key:
      case shards::Triangle<4>::key:
      case shards::Triangle<6>::key: {
        const auto v0ord = parentCell.getNodeMap(subcellDim, subcellOrd, 0);
        const auto v1ord = parentCell.getNodeMap(subcellDim, subcellOrd, 1);
        const auto v2ord = parentCell.getNodeMap(subcellDim, subcellOrd, 2);

        const auto v0 = Kokkos::subview(refNodes, v0ord, Kokkos::ALL());
        const auto v1 = Kokkos::subview(refNodes, v1ord, Kokkos::ALL());
        const auto v2 = Kokkos::subview(refNodes, v2ord, Kokkos::ALL());

        // x(u,v) = x0 + (x1 - x0)*u + (x2 - x0)*v
        subcellParam(subcellOrd, 0, 0) = v0(0);
        subcellParam(subcellOrd, 0, 1) = v1(0) - v0(0);
        subcellParam(subcellOrd, 0, 2) = v2(0) - v0(0);

        // y(u,v) = y0 + (y1 - y0)*u + (y2 - y0)*v
        subcellParam(subcellOrd, 1, 0) = v0(1);
        subcellParam(subcellOrd, 1, 1) = v1(1) - v0(1);
        subcellParam(subcellOrd, 1, 2) = v2(1) - v0(1);

        // z(u,v) = z0 + (z1 - z0)*u + (z2 - z0)*v
        subcellParam(subcellOrd, 2, 0) = v0(2);
        subcellParam(subcellOrd, 2, 1) = v1(2) - v0(2);
        subcellParam(subcellOrd, 2, 2) = v2(2) - v0(2);
        break;
      }
      case shards::Quadrilateral<4>::key:
      case shards::Quadrilateral<8>::key:
      case shards::Quadrilateral<9>::key: {
        const auto v0ord = parentCell.getNodeMap(subcellDim, subcellOrd, 0);
        const auto v1ord = parentCell.getNodeMap(subcellDim, subcellOrd, 1);
        const auto v2ord = parentCell.getNodeMap(subcellDim, subcellOrd, 2);
        const auto v3ord = parentCell.getNodeMap(subcellDim, subcellOrd, 3);

        const auto v0 = Kokkos::subview(refNodes, v0ord, Kokkos::ALL());
        const auto v1 = Kokkos::subview(refNodes, v1ord, Kokkos::ALL());
        const auto v2 = Kokkos::subview(refNodes, v2ord, Kokkos::ALL());
        const auto v3 = Kokkos::subview(refNodes, v3ord, Kokkos::ALL());

        // x(u,v) = (x0+x1+x2+x3)/4+u*(-x0+x1+x2-x3)/4+v*(-x0-x1+x2+x3)/4+uv*(0=x0-x1+x2-x3)/4
        subcellParam(subcellOrd, 0, 0) = ( v0(0) + v1(0) + v2(0) + v3(0))/4.0;
        subcellParam(subcellOrd, 0, 1) = (-v0(0) + v1(0) + v2(0) - v3(0))/4.0;
        subcellParam(subcellOrd, 0, 2) = (-v0(0) - v1(0) + v2(0) + v3(0))/4.0;

        // y(u,v) = (y0+y1+y2+y3)/4+u*(-y0+y1+y2-y3)/4+v*(-y0-y1+y2+y3)/4+uv*(0=y0-y1+y2-y3)/4
        subcellParam(subcellOrd, 1, 0) = ( v0(1) + v1(1) + v2(1) + v3(1))/4.0;
        subcellParam(subcellOrd, 1, 1) = (-v0(1) + v1(1) + v2(1) - v3(1))/4.0;
        subcellParam(subcellOrd, 1, 2) = (-v0(1) - v1(1) + v2(1) + v3(1))/4.0;

        // z(u,v) = (z0+z1+z2+z3)/4+u*(-z0+z1+z2-z3)/4+v*(-z0-z1+z2+z3)/4+uv*(0=z0-z1+z2-z3)/4
        subcellParam(subcellOrd, 2, 0) = ( v0(2) + v1(2) + v2(2) + v3(2))/4.0;
        subcellParam(subcellOrd, 2, 1) = (-v0(2) + v1(2) + v2(2) - v3(2))/4.0;
        subcellParam(subcellOrd, 2, 2) = (-v0(2) - v1(2) + v2(2) + v3(2))/4.0;
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
            ">>> ERROR (Intrepid2::RefSubcellParametrization::set): parametrization not defined for the specified face topology.");
      }
      }
    }
  }
}



template<typename DeviceType>
bool
RefSubcellParametrization<DeviceType>::
isSubcellParametrizationSet_ = false;

#define DefineStaticRefParametrization(obj) template<typename DeviceType> \
    typename RefSubcellParametrization<DeviceType>::ViewType \
    RefSubcellParametrization<DeviceType>:: \
    obj = typename RefSubcellParametrization<DeviceType>::ViewType();

DefineStaticRefParametrization(lineEdgesParam)
DefineStaticRefParametrization(triEdgesParam)
DefineStaticRefParametrization(quadEdgesParam)
DefineStaticRefParametrization(shellTriEdgesParam)
DefineStaticRefParametrization(shellQuadEdgesParam)
DefineStaticRefParametrization(tetEdgesParam)
DefineStaticRefParametrization(hexEdgesParam)
DefineStaticRefParametrization(pyrEdgesParam)
DefineStaticRefParametrization(wedgeEdgesParam)
DefineStaticRefParametrization(shellTriFacesParam)
DefineStaticRefParametrization(shellQuadFacesParam)
DefineStaticRefParametrization(tetFacesParam)
DefineStaticRefParametrization(hexFacesParam)
DefineStaticRefParametrization(pyrFacesParam)
DefineStaticRefParametrization(wedgeFacesParam)


template<typename DeviceType>
void
RefCellNodes<DeviceType>::set() {

  if(isReferenceNodeDataSet_)
    return;

  auto createDataViewFromHostArray = [](const std::string& view_name, double const * source_array, ordinal_type dim){
    ViewType dest_view(view_name, dim, 3);
    auto host_view = Kokkos::create_mirror_view(dest_view);
    for(ordinal_type i=0; i<dim; ++i)
      for(ordinal_type j=0; j<3; ++j)
        host_view(i,j) = source_array[3*i+j];
    Kokkos::deep_copy(dest_view,host_view);
    return dest_view;
  };

  {
    // create memory on devices
    lineNodes           = createDataViewFromHostArray("CellTools::ReferenceNodeData::line", &refNodeDataStatic_.line[0][0], 2);
    line3Nodes          = createDataViewFromHostArray("CellTools::ReferenceNodeData::line_3", &refNodeDataStatic_.line_3[0][0], 3);

    triangleNodes       = createDataViewFromHostArray("CellTools::ReferenceNodeData::triangle", &refNodeDataStatic_.triangle[0][0], 3);
    triangle4Nodes      = createDataViewFromHostArray("CellTools::ReferenceNodeData::triangle_4", &refNodeDataStatic_.triangle_4[0][0], 4);
    triangle6Nodes      = createDataViewFromHostArray("CellTools::ReferenceNodeData::triangle_6", &refNodeDataStatic_.triangle_6[0][0], 6);

    quadrilateralNodes  = createDataViewFromHostArray("CellTools::ReferenceNodeData::quad", &refNodeDataStatic_.quadrilateral[0][0], 4);
    quadrilateral8Nodes = createDataViewFromHostArray("CellTools::ReferenceNodeData::quad_8", &refNodeDataStatic_.quadrilateral_8[0][0], 8);
    quadrilateral9Nodes = createDataViewFromHostArray("CellTools::ReferenceNodeData::quad_9", &refNodeDataStatic_.quadrilateral_9[0][0], 9);

    tetrahedronNodes    = createDataViewFromHostArray("CellTools::ReferenceNodeData::tet", &refNodeDataStatic_.tetrahedron[0][0], 4);
    tetrahedron8Nodes   = createDataViewFromHostArray("CellTools::ReferenceNodeData::tet_8", &refNodeDataStatic_.tetrahedron_8[0][0], 8);
    tetrahedron10Nodes  = createDataViewFromHostArray("CellTools::ReferenceNodeData::tet_10", &refNodeDataStatic_.tetrahedron_10[0][0], 10);
    tetrahedron11Nodes  = createDataViewFromHostArray("CellTools::ReferenceNodeData::tet_11", &refNodeDataStatic_.tetrahedron_11[0][0], 11);

    hexahedronNodes     = createDataViewFromHostArray("CellTools::ReferenceNodeData::hex", &refNodeDataStatic_.hexahedron[0][0], 8);
    hexahedron20Nodes   = createDataViewFromHostArray("CellTools::ReferenceNodeData::hex_20", &refNodeDataStatic_.hexahedron_20[0][0], 20);
    hexahedron27Nodes   = createDataViewFromHostArray("CellTools::ReferenceNodeData::hex_27", &refNodeDataStatic_.hexahedron_27[0][0], 27);

    pyramidNodes        = createDataViewFromHostArray("CellTools::ReferenceNodeData::pyr", &refNodeDataStatic_.pyramid[0][0], 5);
    pyramid13Nodes      = createDataViewFromHostArray("CellTools::ReferenceNodeData::pyr_13", &refNodeDataStatic_.pyramid_13[0][0], 13);
    pyramid14Nodes      = createDataViewFromHostArray("CellTools::ReferenceNodeData::pyr_14", &refNodeDataStatic_.pyramid_14[0][0], 14);

    wedgeNodes          = createDataViewFromHostArray("CellTools::ReferenceNodeData::wedge", &refNodeDataStatic_.wedge[0][0], 6);
    wedge15Nodes        = createDataViewFromHostArray("CellTools::ReferenceNodeData::wedge_15", &refNodeDataStatic_.wedge_15[0][0], 15);
    wedge18Nodes        = createDataViewFromHostArray("CellTools::ReferenceNodeData::wedge_18", &refNodeDataStatic_.wedge_18[0][0], 18);
  }

  Kokkos::push_finalize_hook( [=] {

    lineNodes           = ViewType();
    line3Nodes          = ViewType();

    triangleNodes       = ViewType();
    triangle4Nodes      = ViewType();
    triangle6Nodes      = ViewType();

    quadrilateralNodes  = ViewType();
    quadrilateral8Nodes = ViewType();
    quadrilateral9Nodes = ViewType();

    tetrahedronNodes    = ViewType();
    tetrahedron8Nodes   = ViewType();
    tetrahedron10Nodes  = ViewType();
    tetrahedron11Nodes  = ViewType();

    hexahedronNodes     = ViewType();
    hexahedron20Nodes   = ViewType();
    hexahedron27Nodes   = ViewType();

    pyramidNodes        = ViewType();
    pyramid13Nodes      = ViewType();
    pyramid14Nodes      = ViewType();

    wedgeNodes          = ViewType();
    wedge15Nodes        = ViewType();
    wedge18Nodes        = ViewType();
  } );

  isReferenceNodeDataSet_ = true;
}

template<typename DeviceType>
inline
typename RefCellNodes<DeviceType>::ConstViewType
RefCellNodes<DeviceType>::get(const unsigned cellTopoKey){

  if(!isReferenceNodeDataSet_)
    set();

  ViewType refNodes;

  switch (cellTopoKey ) {
  case shards::Line<2>::key:
  case shards::ShellLine<2>::key:
  case shards::Beam<2>::key:               refNodes = lineNodes; break;
  case shards::Line<3>::key:
  case shards::ShellLine<3>::key:
  case shards::Beam<3>::key:               refNodes = line3Nodes; break;

  case shards::Triangle<3>::key:
  case shards::ShellTriangle<3>::key:      refNodes = triangleNodes; break;
  case shards::Triangle<4>::key:           refNodes = triangle4Nodes; break;
  case shards::Triangle<6>::key:
  case shards::ShellTriangle<6>::key:      refNodes = triangle6Nodes; break;

  case shards::Quadrilateral<4>::key:
  case shards::ShellQuadrilateral<4>::key: refNodes = quadrilateralNodes; break;
  case shards::Quadrilateral<8>::key:
  case shards::ShellQuadrilateral<8>::key: refNodes = quadrilateral8Nodes; break;
  case shards::Quadrilateral<9>::key:
  case shards::ShellQuadrilateral<9>::key: refNodes = quadrilateral9Nodes; break;

  case shards::Tetrahedron<4>::key:        refNodes = tetrahedronNodes; break;
  case shards::Tetrahedron<8>::key:        refNodes = tetrahedron8Nodes; break;
  case shards::Tetrahedron<10>::key:       refNodes = tetrahedron10Nodes; break;
  case shards::Tetrahedron<11>::key:       refNodes = tetrahedron11Nodes; break;

  case shards::Hexahedron<8>::key:         refNodes = hexahedronNodes; break;
  case shards::Hexahedron<20>::key:        refNodes = hexahedron20Nodes; break;
  case shards::Hexahedron<27>::key:        refNodes = hexahedron27Nodes; break;

  case shards::Pyramid<5>::key:            refNodes = pyramidNodes; break;
  case shards::Pyramid<13>::key:           refNodes = pyramid13Nodes; break;
  case shards::Pyramid<14>::key:           refNodes = pyramid14Nodes; break;

  case shards::Wedge<6>::key:              refNodes = wedgeNodes; break;
  case shards::Wedge<15>::key:             refNodes = wedge15Nodes; break;
  case shards::Wedge<18>::key:             refNodes = wedge18Nodes; break;

  default: {
    INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE( true, std::invalid_argument,
        ">>> ERROR (Intrepid2::CellTools::getReferenceNode): invalid cell topology.");
  }
  }
  return refNodes;
}

template<typename DeviceType>
bool
RefCellNodes<DeviceType>::
isReferenceNodeDataSet_ = false;

#define DefineStaticRefNodes(obj) template<typename DeviceType> \
    typename RefCellNodes<DeviceType>::ViewType \
    RefCellNodes<DeviceType>:: \
    obj = typename RefCellNodes<DeviceType>::ViewType();

DefineStaticRefNodes(lineNodes)
DefineStaticRefNodes(line3Nodes)

DefineStaticRefNodes(triangleNodes)
DefineStaticRefNodes(triangle4Nodes)
DefineStaticRefNodes(triangle6Nodes)

DefineStaticRefNodes(quadrilateralNodes)
DefineStaticRefNodes(quadrilateral8Nodes)
DefineStaticRefNodes(quadrilateral9Nodes)

DefineStaticRefNodes(tetrahedronNodes)
DefineStaticRefNodes(tetrahedron8Nodes)
DefineStaticRefNodes(tetrahedron10Nodes)
DefineStaticRefNodes(tetrahedron11Nodes)

DefineStaticRefNodes(hexahedronNodes)
DefineStaticRefNodes(hexahedron20Nodes)
DefineStaticRefNodes(hexahedron27Nodes)

DefineStaticRefNodes(pyramidNodes)
DefineStaticRefNodes(pyramid13Nodes)
DefineStaticRefNodes(pyramid14Nodes)

DefineStaticRefNodes(wedgeNodes)
DefineStaticRefNodes(wedge15Nodes)
DefineStaticRefNodes(wedge18Nodes)

template<typename DeviceType>
const typename RefCellNodes<DeviceType>::ReferenceNodeDataStatic
RefCellNodes<DeviceType>::
refNodeDataStatic_ = {
    // line
    { // 2
        {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}
    },
    { // 3
        {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 0.0, 0.0}
    },
    // triangle
    { // 3
        { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}
    },
    { // 4
        { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 1.0/3.0, 1.0/3.0, 0.0}
    },
    { // 6
        { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
        { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}
    },
    // quad
    { // 4
        {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}
    },
    { // 8
        {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
        { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}
    },
    { // 9
        {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
        { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0}, { 0.0, 0.0, 0.0}
    },
    // tet
    { // 4
        { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0}
    },
    { // 8
        { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
        { 1.0/3.0, 0.0, 1.0/3.0}, { 1.0/3.0, 1.0/3.0, 1.0/3.0}, { 1.0/3.0, 1.0/3.0, 0.0}, { 0.0, 1.0/3.0, 1.0/3.0}
    },
    { // 10
        { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
        { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}, { 0.0, 0.0, 0.5}, { 0.5, 0.0, 0.5}, { 0.0, 0.5, 0.5}
    },
    { // 11
        { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
        { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}, { 0.0, 0.0, 0.5}, { 0.5, 0.0, 0.5}, { 0.0, 0.5, 0.5}
    },
    // hex
    { // 8
        {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
        {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0}
    },
    { // 20
        {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
        {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0},
        { 0.0,-1.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, {-1.0, 0.0,-1.0},
        {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
        { 0.0,-1.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0}, {-1.0, 0.0, 1.0}
    },
    { // 27
        {-1.0,-1.0,-1.0}, { 1.0,-1.0,-1.0}, { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0},
        {-1.0,-1.0, 1.0}, { 1.0,-1.0, 1.0}, { 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0},
        { 0.0,-1.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, {-1.0, 0.0,-1.0},
        {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0},
        { 0.0,-1.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0}, {-1.0, 0.0, 1.0},
        { 0.0, 0.0, 0.0},
        { 0.0, 0.0,-1.0}, { 0.0, 0.0, 1.0}, {-1.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, {0.0,-1.0, 0.0}, {0.0, 1.0, 0.0}
    },
    // pyramid
    { // 5
        {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0}
    },
    { // 13
        {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
        { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0},
        {-0.5,-0.5, 0.5}, { 0.5,-0.5, 0.5}, { 0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}
    },
    { // 14
        {-1.0,-1.0, 0.0}, { 1.0,-1.0, 0.0}, { 1.0, 1.0, 0.0}, {-1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0},
        { 0.0,-1.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0}, {-1.0, 0.0, 0.0},
        {-0.5,-0.5, 0.5}, { 0.5,-0.5, 0.5}, { 0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}, { 0.0, 0.0, 0.0}
    },
    // wedge
    { // 6
        { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0}
    },
    { // 15
        { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0},
        { 0.5, 0.0,-1.0}, { 0.5, 0.5,-1.0}, { 0.0, 0.5,-1.0}, { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
        { 0.5, 0.0, 1.0}, { 0.5, 0.5, 1.0}, { 0.0, 0.5, 1.0}
    },
    { // 18
        { 0.0, 0.0,-1.0}, { 1.0, 0.0,-1.0}, { 0.0, 1.0,-1.0}, { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 0.0, 1.0, 1.0},
        { 0.5, 0.0,-1.0}, { 0.5, 0.5,-1.0}, { 0.0, 0.5,-1.0}, { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 0.0, 1.0, 0.0},
        { 0.5, 0.0, 1.0}, { 0.5, 0.5, 1.0}, { 0.0, 0.5, 1.0},
        { 0.5, 0.0, 0.0}, { 0.5, 0.5, 0.0}, { 0.0, 0.5, 0.0}
    }
};

template<typename DeviceType>
void
RefCellCenter<DeviceType>::set() {

  if(isReferenceCellCenterDataSet_)
    return;

  auto createDataViewFromHostArray = [](const std::string& view_name, double const * source_array){
    ViewType dest_view(view_name, 3);
    auto host_view = Kokkos::create_mirror_view(dest_view);
    for(ordinal_type i=0; i<3; ++i) host_view[i] = source_array[i];
    Kokkos::deep_copy(dest_view, host_view);
    return dest_view;
  };

  {
    // create memory on devices
    lineCenter            = createDataViewFromHostArray("CellTools::ReferenceCenterData::line", &refCenterDataStatic_.line[0]);

    triangleCenter        = createDataViewFromHostArray("CellTools::ReferenceCenterData::triangle", &refCenterDataStatic_.triangle[0]);

    quadrilateralCenter   = createDataViewFromHostArray("CellTools::ReferenceCenterData::quad", &refCenterDataStatic_.quadrilateral[0]);

    tetrahedronCenter     = createDataViewFromHostArray("CellTools::ReferenceCenterData::tet", &refCenterDataStatic_.tetrahedron[0]);

    hexahedronCenter      = createDataViewFromHostArray("CellTools::ReferenceCenterData::hex", &refCenterDataStatic_.hexahedron[0]);

    pyramidCenter         = createDataViewFromHostArray("CellTools::ReferenceCenterData::pyr", &refCenterDataStatic_.pyramid[0]);

    wedgeCenter           = createDataViewFromHostArray("CellTools::ReferenceCenterData::wedge", &refCenterDataStatic_.wedge[0]);
  }

  Kokkos::push_finalize_hook( [=] {

    lineCenter            = ViewType();

    triangleCenter        = ViewType();

    quadrilateralCenter   = ViewType();

    tetrahedronCenter     = ViewType();

    hexahedronCenter      = ViewType();

    pyramidCenter         = ViewType();

    wedgeCenter           = ViewType();
  } );

  isReferenceCellCenterDataSet_ = true;
}

template<typename DeviceType>
inline
typename RefCellCenter<DeviceType>::ConstViewType
RefCellCenter<DeviceType>::get(const unsigned cellTopoKey){

  if(!isReferenceCellCenterDataSet_)
    set();

  ViewType cellCenter;

  switch (cellTopoKey ) {
  case shards::Line<2>::key:
  case shards::ShellLine<2>::key:
  case shards::Beam<2>::key:
  case shards::Line<3>::key:
  case shards::ShellLine<3>::key:
  case shards::Beam<3>::key:               cellCenter = lineCenter; break;

  case shards::Triangle<3>::key:
  case shards::ShellTriangle<3>::key:
  case shards::Triangle<4>::key:
  case shards::Triangle<6>::key:
  case shards::ShellTriangle<6>::key:      cellCenter = triangleCenter; break;

  case shards::Quadrilateral<4>::key:
  case shards::ShellQuadrilateral<4>::key:
  case shards::Quadrilateral<8>::key:
  case shards::ShellQuadrilateral<8>::key:
  case shards::Quadrilateral<9>::key:
  case shards::ShellQuadrilateral<9>::key: cellCenter = quadrilateralCenter; break;

  case shards::Tetrahedron<4>::key:
  case shards::Tetrahedron<8>::key:
  case shards::Tetrahedron<10>::key:
  case shards::Tetrahedron<11>::key:       cellCenter = tetrahedronCenter; break;

  case shards::Hexahedron<8>::key:
  case shards::Hexahedron<20>::key:
  case shards::Hexahedron<27>::key:        cellCenter = hexahedronCenter; break;

  case shards::Pyramid<5>::key:
  case shards::Pyramid<13>::key:
  case shards::Pyramid<14>::key:           cellCenter = pyramidCenter; break;

  case shards::Wedge<6>::key:
  case shards::Wedge<15>::key:
  case shards::Wedge<18>::key:             cellCenter = wedgeCenter; break;

  default: {
    INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE( true, std::invalid_argument,
        ">>> ERROR (Intrepid2::CellTools::getReferenceCellCenter): invalid cell topology.");
  }
  }
  return cellCenter;
}

template<typename DeviceType>
bool
RefCellCenter<DeviceType>::
isReferenceCellCenterDataSet_ = false;

#define DefineStaticRefCenter(obj) template<typename DeviceType> \
    typename RefCellCenter<DeviceType>::ViewType \
    RefCellCenter<DeviceType>:: \
    obj = typename RefCellCenter<DeviceType>::ViewType();

DefineStaticRefCenter(lineCenter)
DefineStaticRefCenter(triangleCenter)
DefineStaticRefCenter(quadrilateralCenter)
DefineStaticRefCenter(tetrahedronCenter)
DefineStaticRefCenter(hexahedronCenter)
DefineStaticRefCenter(pyramidCenter)
DefineStaticRefCenter(wedgeCenter)

template<typename DeviceType>
const typename RefCellCenter<DeviceType>::ReferenceCenterDataStatic
RefCellCenter<DeviceType>::
refCenterDataStatic_ = {
    // line
    {0.0, 0.0, 0.0},
    // triangle
    { 1.0/3.0, 1.0/3.0, 0.0},
    // quad
    {0.0, 0.0, 0.0},
    // tet
    { 0.25, 0.25, 0.25},
    // hex
    { 0.0, 0.0, 0.0},
    // pyramid
    { 0.0, 0.0, 0.25},
    // wedge
    { 1.0/3.0, 1.0/3.0, 0.0},
};


// Point Inclusion 


  template<typename PointViewType, typename ScalarType>
  KOKKOS_INLINE_FUNCTION
  bool
  PointInclusion<shards::Line<>::key>::
  check(const PointViewType &point, const ScalarType threshold) {
    //this implementation should work when PointType is a Sacado Fad<ScalarType>
    const ScalarType minus_one = -1.0 - threshold, plus_one = 1.0 + threshold;
    return (minus_one <= point(0) && point(0) <= plus_one);
  }  

  template<typename PointViewType, typename ScalarType>
  KOKKOS_INLINE_FUNCTION
  bool
  PointInclusion<shards::Triangle<>::key>::
  check(const PointViewType &point, const ScalarType threshold) {
    //this implementation should work when PointType is a Sacado Fad<ScalarType>
    using PointType = typename PointViewType::value_type; 
    const PointType one = 1.0;
    const PointType distance = max( max( -point(0), -point(1) ), point(0) + point(1) - one );
    return distance < threshold;
  }
  
  template<typename PointViewType, typename ScalarType>
  KOKKOS_INLINE_FUNCTION
  bool
  PointInclusion<shards::Quadrilateral<>::key>::
  check(const PointViewType &point, const ScalarType threshold) {
    //this implementation should work when PointType is a Sacado Fad<ScalarType>
    const ScalarType minus_one = -1.0 - threshold, plus_one = 1.0 + threshold;
    return ((minus_one <= point(0) && point(0) <= plus_one) &&
            (minus_one <= point(1) && point(1) <= plus_one));
  }  

  template<typename PointViewType, typename ScalarType>
  KOKKOS_INLINE_FUNCTION
  bool
  PointInclusion<shards::Tetrahedron<>::key>::
  check(const PointViewType &point, const ScalarType threshold) {
    //this implementation should work when PointType is a Sacado Fad<ScalarType>
    using PointType = typename PointViewType::value_type; 
    const PointType one = 1.0;
    const PointType distance = max( max(-point(0),-point(1)),
                                  max(-point(2), point(0) + point(1) + point(2) - one) );
    return distance < threshold;
  }

  template<typename PointViewType, typename ScalarType>
  KOKKOS_INLINE_FUNCTION
  bool
  PointInclusion<shards::Hexahedron<>::key>::
  check(const PointViewType &point, const ScalarType threshold) {
    //this implementation should work when PointType is a Sacado Fad<ScalarType>
    const ScalarType minus_one = -1.0 - threshold, plus_one = 1.0 + threshold;
    return ((minus_one <= point(0) && point(0) <= plus_one) &&
            (minus_one <= point(1) && point(1) <= plus_one) &&
            (minus_one <= point(2) && point(2) <= plus_one));
  }
  
  template<typename PointViewType, typename ScalarType>
  KOKKOS_INLINE_FUNCTION
  bool
  PointInclusion<shards::Pyramid<>::key>::
  check(const PointViewType &point, const ScalarType threshold) {
    //this implementation should work when PointType is a Sacado Fad<ScalarType>
    using PointType = typename PointViewType::value_type; 
    const PointType minus_one = -1.0 - threshold, plus_one = 1.0 + threshold, minus_zero = -threshold;
    const PointType left = minus_one + point(2);
    const PointType right =  plus_one - point(2);
    return ((left       <= point(0) && point(0) <= right) &&
            (left       <= point(1) && point(1) <= right) &&
            (minus_zero <= point(2) && point(2) <= plus_one));
  }

  template<typename PointViewType, typename ScalarType>
  KOKKOS_INLINE_FUNCTION
  bool
  PointInclusion<shards::Wedge<>::key>::
  check(const PointViewType &point, const ScalarType threshold) {
    //this implementation should work when PointType is a Sacado Fad<ScalarType>
    using PointType = typename PointViewType::value_type; 
    const ScalarType minus_one = -1.0 - threshold, plus_one = 1.0 + threshold;
    const PointType one = 1.0;
    const PointType distance = max( max( -point(0), -point(1) ), point(0) + point(1) - one );
    return (distance < threshold && (minus_one <= point(2) && point(2) <= plus_one));
  }

}  // Intrepid2 namespace

#endif

