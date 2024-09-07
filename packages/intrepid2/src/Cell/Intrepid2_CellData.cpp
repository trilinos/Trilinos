// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CellData.cpp
    \brief  Definition file for CellData
    \author Mauro Perego
 */

#ifndef __INTREPID2_CELLDATA_CPP__
#define __INTREPID2_CELLDATA_CPP__

#include "Intrepid2_CellData.hpp"

const CellTopologyData*
Intrepid2::getCellTopologyData(const unsigned& cellTopologyKey){
    const CellTopologyData* cellTopologyData;
    switch (cellTopologyKey) {
      case shards::Line<2>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Line<2>>();
      break;
      case shards::Triangle<3>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Triangle<3>>();
      break;
      case shards::Quadrilateral<4>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Quadrilateral<4>>();
      break;
      case shards::Tetrahedron<4>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Tetrahedron<4>>();
      break;
      case shards::Hexahedron<8>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Hexahedron<8>>();
      break;
      case shards::Wedge<6>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Wedge<6>>();
      break;
      case shards::Pyramid<5>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Pyramid<5>>();
      break;
      
      case shards::Line<3>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Line<3>>();
      break;
      case shards::Triangle<6>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Triangle<6>>();
      break;
      case shards::Quadrilateral<8>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Quadrilateral<8>>();
      break;
      case shards::Quadrilateral<9>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Quadrilateral<9>>();
      break;
      case shards::Tetrahedron<10>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Tetrahedron<10>>();
      break;
      case shards::Tetrahedron<11>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Tetrahedron<11>>();
      break;
      case shards::Hexahedron<20>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Hexahedron<20>>();
      break;
      case shards::Hexahedron<27>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Hexahedron<27>>();
      break;
      case shards::Wedge<15>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Wedge<15>>();
      break;
      case shards::Wedge<18>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Wedge<18>>();
      break;
      case shards::Pyramid<13>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Pyramid<13>>();
      break;        
      case shards::Pyramid<14>::key:
      cellTopologyData = shards::getCellTopologyData<shards::Pyramid<14>>();
      break;
      default: {
          INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
          ">>> ERROR (Intrepid2::getBaseCellTopology): invalid cell topology.");
      }
    }
    return cellTopologyData;
}

#endif

