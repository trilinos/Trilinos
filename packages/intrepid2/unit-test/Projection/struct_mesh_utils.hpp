// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Functions for creating simple structured meshes out of the hypercube [-1,1]^d, d=2,3.

    The functions take the integers NX, NY, NZ indicating the number of cells in the x, y and z directions.
    It returns the rank-2 view cellVertices, containig the vertex id given a cell and the local vertex numbering on the cell.
    It also returns the rank-2 view vertexCoords, containing the vertex coordinates for a given vertex id. 

    \author Created by Mauro Perego
 */
#ifndef __INTREPID2_TEST_STRUCTMESHUTILS_HPP__
#define __INTREPID2_TEST_STRUCTMESHUTILS_HPP__

#include "Intrepid2_config.h"
#include "Intrepid2_CellTools.hpp"
#include <random>

namespace Intrepid2 {

namespace Test {


  template<typename ViewType, typename IntViewType>
  void createStructHexMesh(ViewType& vertexCoords, IntViewType& cellVertices, ordinal_type NX, ordinal_type NY, ordinal_type NZ, bool perturb = false, std::ostream & outStream = std::cout){

    // *********************************** CELL TOPOLOGY **********************************

    // Get cell topology for base hexahedron
    using CellTopology = shards::CellTopology;
    CellTopology cellTopo(shards::getCellTopologyData<shards::Hexahedron<8> >() );

    // Get dimensions
    const ordinal_type dim = cellTopo.getDimension();
    ordinal_type numVerticesPerCell = cellTopo.getNodeCount();

    // *********************************** GENERATE MESH ************************************

    outStream << "Generating Hex mesh ... \n\n";

    outStream << "    NX" << "   NY" << "   NZ\n";
    outStream << std::setw(5) << NX <<
        std::setw(5) << NY <<
        std::setw(5) << NZ << "\n\n";

    // Print mesh information
    ordinal_type numCells = NX*NY*NZ;
    ordinal_type numVertices = (NX+1)*(NY+1)*(NZ+1);
    //ordinal_type numEdges = (NX)*(NY + 1)*(NZ + 1) + (NX + 1)*(NY)*(NZ + 1) + (NX + 1)*(NY + 1)*(NZ);
    //ordinal_type numFaces = (NX)*(NY)*(NZ + 1) + (NX)*(NY + 1)*(NZ) + (NX + 1)*(NY)*(NZ);
    outStream << " Number of Cells: " << numCells << " \n";
    outStream << " Number of Vertices: " << numVertices << " \n";

    // Cube
    double leftX = -1.0, rightX = 1.0;
    double leftY = -1.0, rightY = 1.0;
    double leftZ = -1.0, rightZ = 1.0;

    // Mesh spacing
    double hx = (rightX-leftX)/((double)NX);
    double hy = (rightY-leftY)/((double)NY);
    double hz = (rightZ-leftZ)/((double)NZ);

    // Get nodal coordinates
    vertexCoords = ViewType("vertexCoords", numVertices, dim);
    auto hVertexCoord = Kokkos::create_mirror_view(vertexCoords);
    for (ordinal_type ivertex=0, k=0; k<NZ+1; k++) {
      for (ordinal_type j=0; j<NY+1; j++) {
        for (ordinal_type i=0; i<NX+1; i++) {
          hVertexCoord(ivertex,0) = leftX + (double)i*hx;
          hVertexCoord(ivertex,1) = leftY + (double)j*hy;
          hVertexCoord(ivertex,2) = leftZ + (double)k*hz;
          ivertex++;
        }
      }
    }

    // Perturb mesh coordinates (only interior vertexs)
    if (perturb){
      for (ordinal_type k=1; k<NZ; k++) {
        for (ordinal_type j=1; j<NY; j++) {
          for (ordinal_type i=1; i<NX; i++) {
            ordinal_type ivertex = i + j * (NX + 1) + k * (NX + 1) * (NY + 1);
            // random numbers between -1.0 and 1.0
            double rx = 2.0 * (double)rand()/RAND_MAX - 1.0;
            double ry = 2.0 * (double)rand()/RAND_MAX - 1.0;
            double rz = 2.0 * (double)rand()/RAND_MAX - 1.0;
            // limit variation to 1/4 edge length
            hVertexCoord(ivertex,0) = hVertexCoord(ivertex,0) + 0.125 * hx * rx;
            hVertexCoord(ivertex,1) = hVertexCoord(ivertex,1) + 0.125 * hy * ry;
            hVertexCoord(ivertex,2) = hVertexCoord(ivertex,2) + 0.125 * hz * rz;
          }
        }
      }
    }
    deep_copy(vertexCoords,hVertexCoord);

    // Cell to Vertex map
    cellVertices = IntViewType("cellVertices", numCells, numVerticesPerCell);
    auto hCellVertices = Kokkos::create_mirror_view(cellVertices);
    ordinal_type icell = 0;
    for (ordinal_type k=0; k<NZ; k++) {
      for (ordinal_type j=0; j<NY; j++) {
        for (ordinal_type i=0; i<NX; i++) {
          hCellVertices(icell,0) = (NY + 1)*(NX + 1)*k + (NX + 1)*j + i;
          hCellVertices(icell,1) = (NY + 1)*(NX + 1)*k + (NX + 1)*j + i + 1;
          hCellVertices(icell,2) = (NY + 1)*(NX + 1)*k + (NX + 1)*(j + 1) + i + 1;
          hCellVertices(icell,3) = (NY + 1)*(NX + 1)*k + (NX + 1)*(j + 1) + i;
          hCellVertices(icell,4) = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*j + i;
          hCellVertices(icell,5) = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*j + i + 1;
          hCellVertices(icell,6) = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*(j + 1) + i + 1;
          hCellVertices(icell,7) = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*(j + 1) + i;
          icell++;
        }
      }
    }
    deep_copy(cellVertices,hCellVertices);
  }


  template<typename ViewType, typename IntViewType>
  void createStructTetMesh(ViewType& vertexCoords, IntViewType& cellVertices, ordinal_type NX, ordinal_type NY, ordinal_type NZ, bool perturb = false, std::ostream & outStream = std::cout){

    // *********************************** CELL TOPOLOGY **********************************

    // Get cell topology for base hexahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology cellTopo(shards::getCellTopologyData<shards::Tetrahedron<4> >() );

    // Get dimensions
    const ordinal_type dim = cellTopo.getDimension();
    ordinal_type numVerticesPerCell = cellTopo.getNodeCount();

    // *********************************** GENERATE MESH ************************************

    outStream << "Generating Tet mesh ... \n\n";

    outStream << "    NX" << "   NY" << "   NZ\n";
    outStream << std::setw(5) << NX <<
        std::setw(5) << NY <<
        std::setw(5) << NZ << "\n\n";

    // Print mesh information
    ordinal_type numCells = NX*NY*NZ*5;
    ordinal_type numVertices = (NX+1)*(NY+1)*(NZ+1);
    outStream << " Number of Cells: " << numCells << " \n";
    outStream << " Number of Vertices: " << numVertices << " \n";

    // Cube
    double leftX = -1.0, rightX = 1.0;
    double leftY = -1.0, rightY = 1.0;
    double leftZ = -1.0, rightZ = 1.0;

    // Mesh spacing
    double hx = (rightX-leftX)/((double)NX);
    double hy = (rightY-leftY)/((double)NY);
    double hz = (rightZ-leftZ)/((double)NZ);

    // Get nodal coordinates
    vertexCoords = ViewType("vertexCoords", numVertices, dim);
    auto hVertexCoord = Kokkos::create_mirror_view(vertexCoords);
    for (ordinal_type ivertex=0, k=0; k<NZ+1; k++) {
      for (ordinal_type j=0; j<NY+1; j++) {
        for (ordinal_type i=0; i<NX+1; i++) {
          hVertexCoord(ivertex,0) = leftX + (double)i*hx;
          hVertexCoord(ivertex,1) = leftY + (double)j*hy;
          hVertexCoord(ivertex,2) = leftZ + (double)k*hz;
          ivertex++;
        }
      }
    }

    // Perturb mesh coordinates (only interior vertexs)
    if (perturb){
      for (ordinal_type k=1; k<NZ; k++) {
        for (ordinal_type j=1; j<NY; j++) {
          for (ordinal_type i=1; i<NX; i++) {
            ordinal_type ivertex = i + j * (NX + 1) + k * (NX + 1) * (NY + 1);
            // random numbers between -1.0 and 1.0
            double rx = 2.0 * (double)rand()/RAND_MAX - 1.0;
            double ry = 2.0 * (double)rand()/RAND_MAX - 1.0;
            double rz = 2.0 * (double)rand()/RAND_MAX - 1.0;
            // limit variation to 1/4 edge length
            hVertexCoord(ivertex,0) = hVertexCoord(ivertex,0) + 0.125 * hx * rx;
            hVertexCoord(ivertex,1) = hVertexCoord(ivertex,1) + 0.125 * hy * ry;
            hVertexCoord(ivertex,2) = hVertexCoord(ivertex,2) + 0.125 * hz * rz;
          }
        }
      }
    }
    deep_copy(vertexCoords,hVertexCoord);

    // Cell to Vertex map
    cellVertices = IntViewType("cellVertices", numCells, numVerticesPerCell);
    auto hCellVertices = Kokkos::create_mirror_view(cellVertices);
    ordinal_type icell = 0;
    for (ordinal_type k=0; k<NZ; k++) {
      for (ordinal_type j=0; j<NY; j++) {
        for (ordinal_type i=0; i<NX; i++) {
          auto v0 = (NY + 1)*(NX + 1)*k + (NX + 1)*j + i;
          auto v1 = (NY + 1)*(NX + 1)*k + (NX + 1)*j + i + 1;
          auto v2 = (NY + 1)*(NX + 1)*k + (NX + 1)*(j + 1) + i + 1;
          auto v3 = (NY + 1)*(NX + 1)*k + (NX + 1)*(j + 1) + i;
          auto v4 = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*j + i;
          auto v5 = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*j + i + 1;
          auto v6 = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*(j + 1) + i + 1;
          auto v7 = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*(j + 1) + i;

          hCellVertices(icell,0) = v0;
          hCellVertices(icell,1) = v1;
          hCellVertices(icell,2) = v2;
          hCellVertices(icell,3) = v5;
          icell++;

          hCellVertices(icell,0) = v0;
          hCellVertices(icell,1) = v2;
          hCellVertices(icell,2) = v7;
          hCellVertices(icell,3) = v5;
          icell++;

          hCellVertices(icell,0) = v0;
          hCellVertices(icell,1) = v2;
          hCellVertices(icell,2) = v3;
          hCellVertices(icell,3) = v7;
          icell++;

          hCellVertices(icell,0) = v0;
          hCellVertices(icell,1) = v5;
          hCellVertices(icell,2) = v7;
          hCellVertices(icell,3) = v4;
          icell++;

          hCellVertices(icell,0) = v2;
          hCellVertices(icell,1) = v7;
          hCellVertices(icell,2) = v5;
          hCellVertices(icell,3) = v6;
          icell++;
        }
      }
    }
    deep_copy(cellVertices,hCellVertices);
  }

  template<typename ViewType, typename IntViewType>
  void createStructWedgeMesh(ViewType& vertexCoords, IntViewType& cellVertices, ordinal_type NX, ordinal_type NY, ordinal_type NZ, bool perturb = false, std::ostream & outStream = std::cout){

    // *********************************** CELL TOPOLOGY **********************************

    // Get cell topology for base hexahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology cellTopo(shards::getCellTopologyData<shards::Wedge<6> >() );

    // Get dimensions
    const ordinal_type dim = cellTopo.getDimension();
    ordinal_type numVerticesPerCell = cellTopo.getNodeCount();

    // *********************************** GENERATE MESH ************************************

    outStream << "Generating Tet mesh ... \n\n";

    outStream << "    NX" << "   NY" << "   NZ\n";
    outStream << std::setw(5) << NX <<
        std::setw(5) << NY <<
        std::setw(5) << NZ << "\n\n";

    // Print mesh information
    ordinal_type numCells = NX*NY*NZ*2;
    ordinal_type numVertices = (NX+1)*(NY+1)*(NZ+1);
    outStream << " Number of Cells: " << numCells << " \n";
    outStream << " Number of Vertices: " << numVertices << " \n";

    // Cube
    double leftX = -1.0, rightX = 1.0;
    double leftY = -1.0, rightY = 1.0;
    double leftZ = -1.0, rightZ = 1.0;

    // Mesh spacing
    double hx = (rightX-leftX)/((double)NX);
    double hy = (rightY-leftY)/((double)NY);
    double hz = (rightZ-leftZ)/((double)NZ);

    // Get nodal coordinates
    vertexCoords = ViewType("vertexCoords", numVertices, dim);
    auto hVertexCoord = Kokkos::create_mirror_view(vertexCoords);
    for (ordinal_type ivertex=0, k=0; k<NZ+1; k++) {
      for (ordinal_type j=0; j<NY+1; j++) {
        for (ordinal_type i=0; i<NX+1; i++) {
          hVertexCoord(ivertex,0) = leftX + (double)i*hx;
          hVertexCoord(ivertex,1) = leftY + (double)j*hy;
          hVertexCoord(ivertex,2) = leftZ + (double)k*hz;
          ivertex++;
        }
      }
    }

    // Perturb mesh coordinates (only interior vertexs)
    if (perturb){
      for (ordinal_type k=1; k<NZ; k++) {
        for (ordinal_type j=1; j<NY; j++) {
          for (ordinal_type i=1; i<NX; i++) {
            ordinal_type ivertex = i + j * (NX + 1) + k * (NX + 1) * (NY + 1);
            // random numbers between -1.0 and 1.0
            double rx = 2.0 * (double)rand()/RAND_MAX - 1.0;
            double ry = 2.0 * (double)rand()/RAND_MAX - 1.0;
            double rz = 2.0 * (double)rand()/RAND_MAX - 1.0;
            // limit variation to 1/4 edge length
            hVertexCoord(ivertex,0) = hVertexCoord(ivertex,0) + 0.125 * hx * rx;
            hVertexCoord(ivertex,1) = hVertexCoord(ivertex,1) + 0.125 * hy * ry;
            hVertexCoord(ivertex,2) = hVertexCoord(ivertex,2) + 0.125 * hz * rz;
          }
        }
      }
    }
    deep_copy(vertexCoords,hVertexCoord);

    // Cell to Vertex map
    cellVertices = IntViewType("cellVertices", numCells, numVerticesPerCell);
    auto hCellVertices = Kokkos::create_mirror_view(cellVertices);
    ordinal_type icell = 0;
    for (ordinal_type k=0; k<NZ; k++) {
      for (ordinal_type j=0; j<NY; j++) {
        for (ordinal_type i=0; i<NX; i++) {
          auto v0 = (NY + 1)*(NX + 1)*k + (NX + 1)*j + i;
          auto v1 = (NY + 1)*(NX + 1)*k + (NX + 1)*j + i + 1;
          auto v2 = (NY + 1)*(NX + 1)*k + (NX + 1)*(j + 1) + i + 1;
          auto v3 = (NY + 1)*(NX + 1)*k + (NX + 1)*(j + 1) + i;
          auto v4 = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*j + i;
          auto v5 = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*j + i + 1;
          auto v6 = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*(j + 1) + i + 1;
          auto v7 = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*(j + 1) + i;

          hCellVertices(icell,0) = v0;
          hCellVertices(icell,1) = v1;
          hCellVertices(icell,2) = v2;
          hCellVertices(icell,3) = v4;
          hCellVertices(icell,4) = v5;
          hCellVertices(icell,5) = v6;
          icell++;

          hCellVertices(icell,0) = v0;
          hCellVertices(icell,1) = v2;
          hCellVertices(icell,2) = v3;
          hCellVertices(icell,3) = v4;
          hCellVertices(icell,4) = v6;
          hCellVertices(icell,5) = v7;
          icell++;
        }
      }
    }
    deep_copy(cellVertices,hCellVertices);
  }


  template<typename ViewType, typename IntViewType>
  void createStructQuadMesh(ViewType& vertexCoords, IntViewType& cellVertices, ordinal_type NX, ordinal_type NY, bool perturb = false, std::ostream & outStream = std::cout){

    // *********************************** CELL TOPOLOGY **********************************

    // Get cell topology for base hexahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology cellTopo(shards::getCellTopologyData<shards::Quadrilateral<4> >() );

    // Get dimensions
    const ordinal_type dim = cellTopo.getDimension();
    ordinal_type numVerticesPerCell = cellTopo.getNodeCount();

    // *********************************** GENERATE MESH ************************************

  // *********************************** GENERATE MESH ************************************

    outStream << "Generating Quad mesh ... \n\n";

    outStream << "    NX" << "   NY\n";
    outStream << std::setw(5) << NX <<
        std::setw(5) << NY << "\n\n";

    // Print mesh information
    ordinal_type numCells = NX*NY;
    ordinal_type numVertices = (NX+1)*(NY+1);
    outStream << " Number of Cells: " << numCells << " \n";
    outStream << " Number of Vertices: " << numVertices << " \n";

    // Cube
    double leftX = -1.0, rightX = 1.0;
    double leftY = -1.0, rightY = 1.0;
    // Mesh spacing
    double hx = (rightX-leftX)/((double)NX);
    double hy = (rightY-leftY)/((double)NY);

    // Get nodal coordinates
    vertexCoords = ViewType("vertexCoords", numVertices, dim);
    auto hVertexCoord = Kokkos::create_mirror_view(vertexCoords);
    for (ordinal_type ivertex = 0, j=0; j<NY+1; j++) {
      for (ordinal_type i=0; i<NX+1; i++) {
        hVertexCoord(ivertex,0) = leftX + (double)i*hx;
        hVertexCoord(ivertex,1) = leftY + (double)j*hy;
        ivertex++;
      }
    }

    // Perturb mesh coordinates (only interior vertexs)
    if (perturb){
      for (ordinal_type j=1; j<NY; j++) {
        for (ordinal_type i=1; i<NX; i++) {
          ordinal_type ivertex = i + j * (NX + 1);
          // random numbers between -1.0 and 1.0
          double rx = 2.0 * (double)rand()/RAND_MAX - 1.0;
          double ry = 2.0 * (double)rand()/RAND_MAX - 1.0;
          // limit variation to 1/4 edge length
          hVertexCoord(ivertex,0) = hVertexCoord(ivertex,0) + 0.125 * hx * rx;
          hVertexCoord(ivertex,1) = hVertexCoord(ivertex,1) + 0.125 * hy * ry;
        }
      }
    }
    deep_copy(vertexCoords,hVertexCoord);

    // Cell to Vertex map
    cellVertices = IntViewType("cellVertices", numCells, numVerticesPerCell);
    auto hCellVertices = Kokkos::create_mirror_view(cellVertices);
    ordinal_type icell = 0;

    for (ordinal_type j=0; j<NY; j++) {
      for (ordinal_type i=0; i<NX; i++) {
        hCellVertices(icell,0) = (NX + 1)*j + i;
        hCellVertices(icell,1) = (NX + 1)*j + i + 1;
        hCellVertices(icell,2) = (NX + 1)*(j + 1) + i + 1;
        hCellVertices(icell,3) = (NX + 1)*(j + 1) + i;
        icell++;
      }
    }
    deep_copy(cellVertices,hCellVertices);
  }


  template<typename ViewType, typename IntViewType>
  void createStructTriMesh(ViewType& vertexCoords, IntViewType& cellVertices, ordinal_type NX, ordinal_type NY, bool perturb = false, std::ostream & outStream = std::cout){

    // *********************************** CELL TOPOLOGY **********************************

    // Get cell topology for base hexahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology cellTopo(shards::getCellTopologyData<shards::Quadrilateral<4> >() );

    // Get dimensions
    const ordinal_type dim = cellTopo.getDimension();
    ordinal_type numVerticesPerCell = cellTopo.getNodeCount();

    // *********************************** GENERATE MESH ************************************

  // *********************************** GENERATE MESH ************************************

    outStream << "Generating Tri mesh ... \n\n";

    outStream << "    NX" << "   NY\n";
    outStream << std::setw(5) << NX <<
        std::setw(5) << NY << "\n\n";

    // Print mesh information
    ordinal_type numCells = NX*NY*2;
    ordinal_type numVertices = (NX+1)*(NY+1);
    outStream << " Number of Cells: " << numCells << " \n";
    outStream << " Number of Vertices: " << numVertices << " \n";

    // Cube
    double leftX = -1.0, rightX = 1.0;
    double leftY = -1.0, rightY = 1.0;
    // Mesh spacing
    double hx = (rightX-leftX)/((double)NX);
    double hy = (rightY-leftY)/((double)NY);

    // Get nodal coordinates
    vertexCoords = ViewType("vertexCoords", numVertices, dim);
    auto hVertexCoord = Kokkos::create_mirror_view(vertexCoords);
    for (ordinal_type ivertex = 0, j=0; j<NY+1; j++) {
      for (ordinal_type i=0; i<NX+1; i++) {
        hVertexCoord(ivertex,0) = leftX + (double)i*hx;
        hVertexCoord(ivertex,1) = leftY + (double)j*hy;
        ivertex++;
      }
    }

    // Perturb mesh coordinates (only interior vertexs)
    if (perturb){
      for (ordinal_type j=1; j<NY; j++) {
        for (ordinal_type i=1; i<NX; i++) {
          ordinal_type ivertex = i + j * (NX + 1);
          // random numbers between -1.0 and 1.0
          double rx = 2.0 * (double)rand()/RAND_MAX - 1.0;
          double ry = 2.0 * (double)rand()/RAND_MAX - 1.0;
          // limit variation to 1/4 edge length
          hVertexCoord(ivertex,0) = hVertexCoord(ivertex,0) + 0.125 * hx * rx;
          hVertexCoord(ivertex,1) = hVertexCoord(ivertex,1) + 0.125 * hy * ry;
        }
      }
    }
    deep_copy(vertexCoords,hVertexCoord);

    // Cell to Vertex map
    cellVertices = IntViewType("cellVertices", numCells, numVerticesPerCell);
    auto hCellVertices = Kokkos::create_mirror_view(cellVertices);
    ordinal_type icell = 0;

    for (ordinal_type j=0; j<NY; j++) {
      for (ordinal_type i=0; i<NX; i++) {
        auto v0 = (NX + 1)*j + i;
        auto v1 = (NX + 1)*j + i + 1;
        auto v2 = (NX + 1)*(j + 1) + i + 1;
        auto v3 = (NX + 1)*(j + 1) + i;

        hCellVertices(icell,0) = v0;
        hCellVertices(icell,1) = v1;
        hCellVertices(icell,2) = v3;
        icell++;

        hCellVertices(icell,0) = v1;
        hCellVertices(icell,1) = v2;
        hCellVertices(icell,2) = v3;
        icell++;
      }
    }
    deep_copy(cellVertices,hCellVertices);
  }

  template<typename ViewType, typename IntViewType>
  void createStructMesh(ViewType& vertexCoords, IntViewType& cellVertices, const shards::CellTopology& cellTopo, ordinal_type NX, ordinal_type NY, ordinal_type NZ=-1, bool perturb = false, std::ostream & outStream = std::cout) {
    switch(cellTopo.getKey()) {
      case shards::Hexahedron<8>::key:
        createStructHexMesh(vertexCoords, cellVertices, NX, NY, NZ, perturb, outStream);
        break;
      case shards::Tetrahedron<4>::key:
        createStructTetMesh(vertexCoords, cellVertices, NX, NY, NZ, perturb, outStream);
        break;
      case shards::Wedge<6>::key:
        createStructWedgeMesh(vertexCoords, cellVertices, NX, NY, NZ, perturb, outStream);
        break;
      case shards::Quadrilateral<4>::key:
        createStructQuadMesh(vertexCoords, cellVertices, NX, NY, perturb, outStream);
        break;
      case shards::Triangle<3>::key:
        createStructTriMesh(vertexCoords, cellVertices, NX, NY, perturb, outStream);
        break;
      default:
        INTREPID2_TEST_FOR_EXCEPTION
         (true, std::runtime_error, "Intrepid2::Test::createStructMesh: Topology not supported");
    }
  }

}
}

#endif

