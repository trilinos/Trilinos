// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file
    \brief  Unit test (ControlVolume, Side, Boundary): correctness of values.
    \author Created by K. Peterson, P. Bochev, D. Ridzal and Kyungjoo Kim
*/


#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_CubatureControlVolume.hpp"
//#include "Intrepid2_CubatureControlVolumeSide.hpp"
//#include "Intrepid2_CubatureControlVolumeBoundary.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {

  namespace Test {
#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    } catch (std::logic_error err) {                                    \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    };      

#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)
    
    template<typename ValueType,
             typename DeviceSpaceType,
             typename cellCoordViewType,
             typename exactValueType,
             typename outStreamType>
    int Integration_Test25_CubatureControlVolume(const shards::CellTopology cellTopo,
                                                 const cellCoordViewType cellCoords,
                                                 const exactValueType *exactCubPoints,
                                                 const exactValueType *exactCubWeights,
                                                 outStreamType outStream) {
      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

      typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;
      typedef Kokkos::DynRankView<ValueType,HostSpaceType>   DynRankViewHost;
      
      typedef ValueType pointValueType;
      typedef ValueType weightValueType;
      typedef CubatureControlVolume<DeviceSpaceType,pointValueType,weightValueType> CubatureControlVolumeType;

      const auto tol = 10.0 * tolerence();

      *outStream << " --- CubatureControlVolume testing\n";
      int errorFlag = 0;
      try {
        const auto numCells = cellCoords.dimension(0);
        const auto numVerts = cellCoords.dimension(1);
        const auto spaceDim = cellCoords.dimension(2);

        CubatureControlVolumeType CubatureControlVolume(cellTopo);

        const auto numPoints = CubatureControlVolume.getNumPoints();
        if (numPoints != numVerts) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "Number of cubature points (" << numPoints 
                     << ") does not match to the correct # of points (4)\n";
        }
        
        // get cubature points and weights for volume integration over control volume
        DynRankViewHost ConstructWithLabel(cubPoints,  numCells, numPoints, spaceDim);
        DynRankViewHost ConstructWithLabel(cubWeights, numCells, numPoints);

        CubatureControlVolume.getCubature(cubPoints, cubWeights, cellCoords);

        auto cubPointsHost  = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), cubPoints);
        auto cubWeightsHost = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), cubWeights);

        Kokkos::deep_copy(cubPointsHost,  cubPoints);
        Kokkos::deep_copy(cubWeightsHost, cubWeights);
        
        // validation
        for (ordinal_type i=0,p=0,w=0;i<numCells;++i) 
          for (ordinal_type j=0;j<numPoints;++j,++w) {
            for (ordinal_type k=0;k<spaceDim;++k,++p) {

              // check values of cubature points
              if (std::isnan(cubPointsHost(i,j,k)) || 
                  std::abs(cubPointsHost(i,j,k) - exactCubPoints[p]) > tol) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed point: " << cubPointsHost(i,j,k)
                           << " but reference value: " << exactCubPoints[p] << "\n";
              }
            }
            
            // check values of cubature weights
            if (std::isnan(cubWeightsHost(i,j)) ||
                std::abs(cubWeightsHost(i,j) - exactCubWeights[w]) > tol) {
              errorFlag++;
              *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
              *outStream << " At multi-index { ";
              *outStream << i << " ";*outStream << j << " ";
              *outStream << "}  computed weight: " << cubWeightsHost(i,j)
                         << " but reference value: " << exactCubWeights[w] << "\n";
            }
          }
      } catch (std::exception err) {
        *outStream << err.what() << "\n";
        errorFlag = -10000;
      }  
      return errorFlag;
    }


    template<typename ValueType, typename DeviceSpaceType>
    int Integration_Test25(const bool verbose) {

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs, false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

      *outStream << "DeviceSpace::  "; DeviceSpaceType::print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType::print_configuration(*outStream, false);

      *outStream                                                        
        << "===============================================================================\n" 
        << "|                                                                             |\n" 
        << "|                         Unit Test (CubatureControlVolume)                   |\n" 
        << "|                                   (CubatureControlVolumeSide)               |\n" 
        << "|                                   (CubatureControlVolumeBoundary)           |\n" 
        << "|                                                                             |\n" 
        << "|     1) Correctness of cubature points and weights                           |\n" 
        << "|     2) Comparison of sub-control volume weights and primary cell volume     |\n" 
        << "|     3) Control volume integration                                           |\n" 
        << "|                                                                             |\n" 
        << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" 
        << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" 
        << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" 
        << "|                                                                             |\n" 
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" 
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" 
        << "|                                                                             |\n" 
        << "===============================================================================\n" 
        << "| TEST 1: correctness of cubature points and weights                          |\n" 
        << "===============================================================================\n";
      
      typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;
      typedef Kokkos::DynRankView<ValueType,HostSpaceType>   DynRankViewHost;

      int errorFlag = 0;

      try {
        *outStream << " - Quadrilateral cell testing \n";
        {
          // quadrilateral primary cells
          const ordinal_type numCells = 2, numVerts = 4, spaceDim = 2; 
          DynRankViewHost ConstructWithLabel(cellCoordsHost, numCells, numVerts, spaceDim);
          
          cellCoordsHost(0,0,0) = 0.0; cellCoordsHost(0,0,1) = 0.0;
          cellCoordsHost(0,1,0) = 0.5; cellCoordsHost(0,1,1) = 0.0;
          cellCoordsHost(0,2,0) = 0.5; cellCoordsHost(0,2,1) = 1.0;
          cellCoordsHost(0,3,0) = 0.0; cellCoordsHost(0,3,1) = 1.0;
          
          cellCoordsHost(1,0,0) = 0.5; cellCoordsHost(1,0,1) = 0.0;
          cellCoordsHost(1,1,0) = 1.0; cellCoordsHost(1,1,1) = 0.0;
          cellCoordsHost(1,2,0) = 1.0; cellCoordsHost(1,2,1) = 1.0;
          cellCoordsHost(1,3,0) = 0.5; cellCoordsHost(1,3,1) = 1.0;
          
          auto cellCoords = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), cellCoordsHost);
          Kokkos::deep_copy(cellCoords, cellCoordsHost);
          
          // points 
          const ValueType exactCubPoints[] = {
            0.125, 0.25, 0.375, 0.25,
            0.375, 0.75, 0.125, 0.75,
            0.625, 0.25, 0.875, 0.25,
            0.875, 0.75, 0.625, 0.75
          };
          
          // weights
          const ValueType exactCubWeights[] = {
            0.125, 0.125, 0.125, 0.125,
            0.125, 0.125, 0.125, 0.125
          };
          
          // // side points 
          // const ValueType exactSideCubPoints[] = {
          //   0.25, 0.25, 0.375, 0.5,
          //   0.25, 0.75, 0.125, 0.5,
          //   0.75, 0.25, 0.875, 0.5,
          //   0.75, 0.75, 0.625, 0.5
          // };
          
          // // side weights (these are weighted normals!)
          // const ValueType exactSideCubWeights[] = {
          //   0.5, 0.0, 0.0, 0.25,
          //   -0.5, 0.0, 0.0,-0.25,
          //   0.5, 0.0, 0.0, 0.25,
          //   -0.5, 0.0, 0.0,-0.25
          // };
          
          // // boundary points 
          // const ValueType exactBCCubPoints[] = {
          //   0.375, 1.0, 0.125, 1.0,
          //   0.875, 1.0, 0.625, 1.0
          // };
          
          // // boundary weights 
          // const ValueType exactBCCubWeights[] = {
          //   0.25, 0.25, 0.25, 0.25
          // };
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<> >());
          errorFlag += Integration_Test25_CubatureControlVolume<ValueType,DeviceSpaceType>(cellTopo,
                                                                                           cellCoords,
                                                                                           exactCubPoints,
                                                                                           exactCubWeights,
                                                                                           outStream);
        }
        
        *outStream << " - Triangle cell testing \n";
        {
          // triangle primary cells      
          const ordinal_type numCells = 2, numVerts = 3, spaceDim = 2; 
          DynRankViewHost ConstructWithLabel(cellCoordsHost, numCells, numVerts, spaceDim);
          
          cellCoordsHost(0,0,0) = 0.0; cellCoordsHost(0,0,1) = 0.0;
          cellCoordsHost(0,1,0) = 0.5; cellCoordsHost(0,1,1) = 0.0;
          cellCoordsHost(0,2,0) = 0.0; cellCoordsHost(0,2,1) = 0.5;
          cellCoordsHost(1,0,0) = 0.5; cellCoordsHost(1,0,1) = 0.0;
          cellCoordsHost(1,1,0) = 0.5; cellCoordsHost(1,1,1) = 0.5;
          cellCoordsHost(1,2,0) = 0.0; cellCoordsHost(1,2,1) = 0.5;
          
          auto cellCoords = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), cellCoordsHost);
          Kokkos::deep_copy(cellCoords, cellCoordsHost);
          
          // points 
          const ValueType exactCubPoints[] = {
            0.1041666666666667, 0.1041666666666667, 0.2916666666666667, 0.1041666666666667,
            0.1041666666666667, 0.2916666666666667, 0.3958333333333333, 0.2083333333333333,
            0.3958333333333333, 0.3958333333333333, 0.2083333333333333, 0.3958333333333333
          };
          
          // weights
          const ValueType exactCubWeights[] = {
            0.0416666666666667, 0.0416666666666667, 0.0416666666666667, 
            0.0416666666666667, 0.0416666666666667, 0.0416666666666667 
          };
          
          // // side points 
          // const ValueType exactSideCubPoints[] = {
          //   0.2083333333333333, 0.0833333333333333, 0.2083333333333333, 0.2083333333333333,
          //   0.0833333333333333, 0.2083333333333333, 0.4166666666666667, 0.2916666666666667,
          //   0.2916666666666667, 0.4166666666666667, 0.2916666666666667, 0.2916666666666667
          // };
          
          //  // side weights (these are weighted normals!)
          // const ValueType exactSideCubWeights[] = {
          //   0.1666666666666667, 0.0833333333333333,-0.0833333333333333, 0.0833333333333333,
          //  -0.0833333333333333,-0.1666666666666667, 0.0833333333333333, 0.1666666666666667,
          //  -0.1666666666666667,-0.0833333333333333, 0.0833333333333333,-0.0833333333333333
          // };
          
          // // boundary points 
          // const ValueType exactBCCubPoints[] = {
          //   0.375, 0.125, 0.125, 0.375,
          //   0.375, 0.5, 0.125, 0.5
          // };
          
          // // boundary weights 
          // const ValueType exactBCCubWeights[] = {
          //   0.353553390593274, 0.353553390593274, 0.25, 0.25
          // };
          
          // triangle primary cells      
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData< shards::Triangle<> >());
          errorFlag += Integration_Test25_CubatureControlVolume<ValueType,DeviceSpaceType>(cellTopo,
                                                                                           cellCoords,
                                                                                           exactCubPoints,
                                                                                           exactCubWeights,
                                                                                           outStream);
        }
        *outStream << " - Tetrahedron cell testing \n";
        {
          // tetrahedron primary cells      
          const ordinal_type numCells = 1, numVerts = 4, spaceDim = 3; 
          DynRankViewHost ConstructWithLabel(cellCoordsHost, numCells, numVerts, spaceDim);

          cellCoordsHost(0,0,0) = 0.0; cellCoordsHost(0,0,1) = 0.0; cellCoordsHost(0,0,2) = 0.0;
          cellCoordsHost(0,1,0) = 1.0; cellCoordsHost(0,1,1) = 0.0; cellCoordsHost(0,1,2) = 0.0;
          cellCoordsHost(0,2,0) = 0.0; cellCoordsHost(0,2,1) = 1.0; cellCoordsHost(0,2,2) = 0.0;
          cellCoordsHost(0,3,0) = 0.0; cellCoordsHost(0,3,1) = 0.0; cellCoordsHost(0,3,2) = 1.0;

          auto cellCoords = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), cellCoordsHost);
          Kokkos::deep_copy(cellCoords, cellCoordsHost);

          // points 
          const ValueType exactCubPoints[] = {
            0.17708333333333333,  0.17708333333333333,  0.17708333333333333, 
            0.46875000000000000,  0.17708333333333333,  0.17708333333333333, 
            0.17708333333333333,  0.46875000000000000,  0.17708333333333333, 
            0.17708333333333333,  0.17708333333333333,  0.46875000000000000
          };

          // weights
          const ValueType exactCubWeights[] = {
            0.0416666666666667, 0.0416666666666667, 
            0.0416666666666667, 0.0416666666666667 
          };

          // // side points 
          // const ValueType exactSideCubPoints[] = {
          //   0.3541666666666667, 0.1458333333333333, 0.1458333333333333,
          //   0.3541666666666667, 0.3541666666666667, 0.1458333333333333,
          //   0.1458333333333333, 0.3541666666666667, 0.1458333333333333,
          //   0.1458333333333333, 0.1458333333333333, 0.3541666666666667,
          //   0.3541666666666667, 0.1458333333333333, 0.3541666666666667,
          //   0.1458333333333333, 0.3541666666666667, 0.3541666666666667
          // };

          // // side weights (these are weighted normals!)
          // const ValueType exactSideCubWeights[] = {
          //   0.0833333333333333, 0.0416666666666667, 0.041666666666667,
          //   -0.0416666666666667, 0.0416666666666667, 0.000000000000000,
          //   -0.0416666666666667,-0.0833333333333333,-0.041666666666667,
          //   0.0416666666666667, 0.0416666666666667, 0.083333333333333,
          //   -0.0416666666666667, 0.0000000000000000, 0.041666666666667,
          //   0.0000000000000000,-0.0416666666666667, 0.041666666666667
          // };

          // // boundary points 
          // const ValueType exactBCCubPoints[] = {
          //   0.208333333333333, 0.00, 0.208333333333333,
          //   0.583333333333333, 0.00, 0.208333333333333, 
          //   0.208333333333333, 0.00, 0.583333333333333, 
          // };

          // // boundary weights 
          // const ValueType exactBCCubWeights[] = {
          //   0.166666666666667, 0.166666666666667, 0.166666666666667
          // };

          // tetrahedral primary cells
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<> >());
          errorFlag += Integration_Test25_CubatureControlVolume<ValueType,DeviceSpaceType>(cellTopo,
                                                                                           cellCoords,
                                                                                           exactCubPoints,
                                                                                           exactCubWeights,
                                                                                           outStream);
        }
        *outStream << " - Hexahedron cell testing \n";
        {
          const ordinal_type numCells = 1, numVerts = 8, spaceDim = 3; 
          DynRankViewHost ConstructWithLabel(cellCoordsHost, numCells, numVerts, spaceDim);

          cellCoordsHost(0,0,0) = 0.0; cellCoordsHost(0,0,1) = 0.0; cellCoordsHost(0,0,2) = 0.0;
          cellCoordsHost(0,1,0) = 2.0; cellCoordsHost(0,1,1) = 0.0; cellCoordsHost(0,1,2) = 0.0;
          cellCoordsHost(0,2,0) = 2.0; cellCoordsHost(0,2,1) = 1.5; cellCoordsHost(0,2,2) = 0.0;
          cellCoordsHost(0,3,0) = 0.0; cellCoordsHost(0,3,1) = 1.5; cellCoordsHost(0,3,2) = 0.0;
          cellCoordsHost(0,4,0) = 0.0; cellCoordsHost(0,4,1) = 0.0; cellCoordsHost(0,4,2) = 1.0;
          cellCoordsHost(0,5,0) = 2.0; cellCoordsHost(0,5,1) = 0.0; cellCoordsHost(0,5,2) = 1.0;
          cellCoordsHost(0,6,0) = 2.0; cellCoordsHost(0,6,1) = 1.5; cellCoordsHost(0,6,2) = 1.0;
          cellCoordsHost(0,7,0) = 0.0; cellCoordsHost(0,7,1) = 1.5; cellCoordsHost(0,7,2) = 1.0;

          auto cellCoords = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), cellCoordsHost);
          Kokkos::deep_copy(cellCoords, cellCoordsHost);

          // points 
          const ValueType exactCubPoints[] = {
            0.5, 0.375, 0.25, 1.5, 0.375, 0.25,
            1.5, 1.125, 0.25, 0.5, 1.125, 0.25,
            0.5, 0.375, 0.75, 1.5, 0.375, 0.75,
            1.5, 1.125, 0.75, 0.5, 1.125, 0.75
          };

          // weights
          const ValueType exactCubWeights[] = {
            0.375, 0.375, 0.375, 0.375, 
            0.375, 0.375, 0.375, 0.375
          };

          // // side points 
          // const ValueType exactSideCubPoints[] = {
          //   1.0, 0.375, 0.25, 1.5, 0.750, 0.25,
          //   1.0, 1.125, 0.25, 0.5, 0.750, 0.25, 
          //   1.0, 0.375, 0.75, 1.5, 0.750, 0.75, 
          //   1.0, 1.125, 0.75, 0.5, 0.750, 0.75,
          //   0.5, 0.375, 0.50, 1.5, 0.375, 0.50,
          //   1.5, 1.125, 0.50, 0.5, 1.125, 0.50
          // };

          // // side weights (these are weighted normals!)
          // const ValueType exactSideCubWeights[] = {
          //   0.375, 0.00, 0.00, 0.00, 0.50, 0.00,
          //   -0.375, 0.00, 0.00, 0.00,-0.50, 0.00,
          //   0.375, 0.00, 0.00, 0.00, 0.50, 0.00,
          //   -0.375, 0.00, 0.00, 0.00,-0.50, 0.00,
          //   0.000, 0.00, 0.75, 0.00, 0.00, 0.75,
          //   0.000, 0.00, 0.75, 0.00, 0.00, 0.75,
          // };

          // // boundary points 
          // const ValueType exactBCCubPoints[] = {
          //   0.5, 0.00, 0.25,
          //   1.5, 0.00, 0.25,
          //   1.5, 0.00, 0.75,
          //   0.5, 0.00, 0.75,
          // };

          // // boundary weights 
          // const ValueType exactBCCubWeights[] = {
          //   0.5, 0.5, 0.5, 0.5
          // };

          // hexahedral primary cells
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<> >());
          errorFlag += Integration_Test25_CubatureControlVolume<ValueType,DeviceSpaceType>(cellTopo,
                                                                                           cellCoords,
                                                                                           exactCubPoints,
                                                                                           exactCubWeights,
                                                                                           outStream);

        }
      } catch (std::exception err) {
        *outStream << err.what() << "\n";
        errorFlag = -10000;
      }  

      // *outStream \
      // << "===============================================================================\n"\
      // << "| TEST 2: comparison of sub-control volume weights and primary cell volume    |\n"\
      // << "===============================================================================\n";

      // // quadrilateral primary cells
      // try {

      //   // set cell topology
      //   cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<> >()));

      //   // define control volume cubature rule
      //   controlVolCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
      //   int numPoints = CubatureControlVolume.getNumPoints();

      //   // define quad coordinates
      //   Intrepid2::FieldContainer<double> cellCoords(1,4,2);
      //   cellCoords(0,0,0) = 0.0; cellCoords(0,0,1) = 0.0;
      //   cellCoords(0,1,0) = 2.4; cellCoords(0,1,1) = 0.0;
      //   cellCoords(0,2,0) = 2.4; cellCoords(0,2,1) = 3.1;
      //   cellCoords(0,3,0) = 0.0; cellCoords(0,3,1) = 3.1;

      //   double exact_area = 2.4*3.1;

      //   // get cubature points and weights 
      //   int numCells = 1; int spaceDim = 2;
      //   Intrepid2::FieldContainer<double> cubPoints(numCells,numPoints,spaceDim);
      //   Intrepid2::FieldContainer<double> cubWeights(numCells,numPoints);
      //   CubatureControlVolume.getCubature(cubPoints,cubWeights,cellCoords);

      //   // loop over number of points (equals number of sub-control volumes) and check total volume
      //   double total_area = 0.0;
      //   for (int i=0; i<numPoints; i++) {
      //      total_area += cubWeights(0,i);
      //   }
      //   if (std::abs(total_area - exact_area) > Intrepid2::INTREPID_TOL) {
      //       errorFlag++;
      //       *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      //       *outStream << " Sum of sub-control volume areas: ";
      //       *outStream << total_area;
      //       *outStream << " does not equal quad primary cell area: " << exact_area << "\n";
      //   }

      // }
      // catch (std::exception err) {
      //   *outStream << err.what() << "\n";
      //   errorFlag = -1;
      // }  

      // // triangle primary cells
      // try {

      //   // set cell topology
      //   cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Triangle<> >()));

      //   // define control volume cubature rule
      //   controlVolCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
      //   int numPoints = CubatureControlVolume.getNumPoints();

      //   // define quad coordinates
      //   Intrepid2::FieldContainer<double> cellCoords(1,3,2);
      //   cellCoords(0,0,0) = 0.0; cellCoords(0,0,1) = 0.0;
      //   cellCoords(0,1,0) = 3.6; cellCoords(0,1,1) = 0.0;
      //   cellCoords(0,2,0) = 0.0; cellCoords(0,2,1) = 2.8;

      //   double exact_area = 0.5*3.6*2.8;

      //   // get cubature points and weights 
      //   int numCells = 1; int spaceDim = 2;
      //   Intrepid2::FieldContainer<double> cubPoints(numCells,numPoints,spaceDim);
      //   Intrepid2::FieldContainer<double> cubWeights(numCells,numPoints);
      //   CubatureControlVolume.getCubature(cubPoints,cubWeights,cellCoords);

      //   // loop over number of points (equals number of sub-control volumes) and check total volume
      //   double total_area = 0.0;
      //   for (int i=0; i<numPoints; i++) {
      //      total_area += cubWeights(0,i);
      //   }
      //   if (std::abs(total_area - exact_area) > Intrepid2::INTREPID_TOL) {
      //       errorFlag++;
      //       *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      //       *outStream << " Sum of sub-control volume areas: ";
      //       *outStream << total_area;
      //       *outStream << " does not equal triangle primary cell area: " << exact_area << "\n";
      //   }

      // }
      // catch (std::exception err) {
      //   *outStream << err.what() << "\n";
      //   errorFlag = -1;
      // }  

      // // hexahedral primary cells
      // try {

      //   // set cell topology
      //   cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<> >()));

      //   // define control volume cubature rule
      //   controlVolCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
      //   int numPoints = CubatureControlVolume.getNumPoints();

      //   // define hexahedron coordinates
      //   Intrepid2::FieldContainer<double> cellCoords(1,8,3);
      //   cellCoords(0,0,0) = 0.0; cellCoords(0,0,1) = 0.0; cellCoords(0,0,2) = 0.0;
      //   cellCoords(0,1,0) = 2.4; cellCoords(0,1,1) = 0.0; cellCoords(0,1,2) = 0.0;
      //   cellCoords(0,2,0) = 2.4; cellCoords(0,2,1) = 3.1; cellCoords(0,2,2) = 0.0;
      //   cellCoords(0,3,0) = 0.0; cellCoords(0,3,1) = 3.1; cellCoords(0,3,2) = 0.0;
      //   cellCoords(0,4,0) = 0.0; cellCoords(0,4,1) = 0.0; cellCoords(0,4,2) = 1.7;
      //   cellCoords(0,5,0) = 2.4; cellCoords(0,5,1) = 0.0; cellCoords(0,5,2) = 1.7;
      //   cellCoords(0,6,0) = 2.4; cellCoords(0,6,1) = 3.1; cellCoords(0,6,2) = 1.7;
      //   cellCoords(0,7,0) = 0.0; cellCoords(0,7,1) = 3.1; cellCoords(0,7,2) = 1.7;

      //   double exact_vol = 2.4*3.1*1.7;

      //   // get cubature points and weights 
      //   int numCells = 1; int spaceDim = 3;
      //   Intrepid2::FieldContainer<double> cubPoints(numCells,numPoints,spaceDim);
      //   Intrepid2::FieldContainer<double> cubWeights(numCells,numPoints);
      //   CubatureControlVolume.getCubature(cubPoints,cubWeights,cellCoords);

      //   // loop over number of points (equals number of sub-control volumes) and check total volume
      //   double total_vol = 0.0;
      //   for (int i=0; i<numPoints; i++) {
      //      total_vol += cubWeights(0,i);
      //   }
      //   if (std::abs(total_vol - exact_vol) > Intrepid2::INTREPID_TOL) {
      //       errorFlag++;
      //       *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      //       *outStream << " Sum of sub-control volumes: ";
      //       *outStream << total_vol;
      //       *outStream << " does not equal hexahedron primary cell volume: " << exact_vol << "\n";
      //   }

      // }
      // catch (std::exception err) {
      //   *outStream << err.what() << "\n";
      //   errorFlag = -1;
      // }  

      // // tetrahedral primary cells
      // try {

      //   // set cell topology
      //   cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<> >()));

      //   // define control volume cubature rule
      //   controlVolCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
      //   int numPoints = CubatureControlVolume.getNumPoints();

      //   // define tetrahedron coordinates
      //   Intrepid2::FieldContainer<double> cellCoords(1,4,3);
      //   cellCoords(0,0,0) = 0.0; cellCoords(0,0,1) = 0.0; cellCoords(0,0,2) = 0.0;
      //   cellCoords(0,1,0) = 3.6; cellCoords(0,1,1) = 0.0; cellCoords(0,1,2) = 0.0;
      //   cellCoords(0,2,0) = 0.0; cellCoords(0,2,1) = 2.8; cellCoords(0,2,2) = 0.0;
      //   cellCoords(0,3,0) = 0.0; cellCoords(0,3,1) = 2.8; cellCoords(0,3,2) = 1.7;

      //   double exact_vol = (0.5*3.6*2.8)*1.7/3.0;

      //   // get cubature points and weights 
      //   int numCells = 1; int spaceDim = 3;
      //   Intrepid2::FieldContainer<double> cubPoints(numCells,numPoints,spaceDim);
      //   Intrepid2::FieldContainer<double> cubWeights(numCells,numPoints);
      //   CubatureControlVolume.getCubature(cubPoints,cubWeights,cellCoords);

      //   // loop over number of points (equals number of sub-control volumes) and check total volume
      //   double total_vol = 0.0;
      //   for (int i=0; i<numPoints; i++) {
      //      total_vol += cubWeights(0,i);
      //   }
      //   if (std::abs(total_vol - exact_vol) > Intrepid2::INTREPID_TOL) {
      //       errorFlag++;
      //       *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      //       *outStream << " Sum of sub-control volumes: ";
      //       *outStream << total_vol;
      //       *outStream << " does not equal tetrahedron primary cell volume: " << exact_vol << "\n";
      //   }

      // }
      // catch (std::exception err) {
      //   *outStream << err.what() << "\n";
      //   errorFlag = -1;
      // }  

      // *outStream \
      // << "===============================================================================\n"\
      // << "| TEST 3: control volume integration                                          |\n"\
      // << "===============================================================================\n";

      // // quadrilateral primary cells
      // try {

      //   // set cell topology
      //   cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<> >()));

      //   // define control volume cubature rule
      //   controlVolCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
      //   int numPoints = CubatureControlVolume.getNumPoints();

      //   // define control volume cubature rule
      //   controlVolSideCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolumeSide<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
      //   int numSidePoints = controlVolSideCub->getNumPoints();

      //   // define quad coordinates - four cells that define a complete control volume around center node
      //   Intrepid2::FieldContainer<double> cellCoords(4,4,2);
      //   cellCoords(0,0,0) = 0.0; cellCoords(0,0,1) = 0.0;
      //   cellCoords(0,1,0) = 0.5; cellCoords(0,1,1) = 0.0;
      //   cellCoords(0,2,0) = 0.41; cellCoords(0,2,1) = 0.58;
      //   cellCoords(0,3,0) = 0.0; cellCoords(0,3,1) = 0.5;
      //   cellCoords(1,0,0) = 0.5; cellCoords(1,0,1) = 0.0;
      //   cellCoords(1,1,0) = 1.0; cellCoords(1,1,1) = 0.0;
      //   cellCoords(1,2,0) = 1.0; cellCoords(1,2,1) = 0.5;
      //   cellCoords(1,3,0) = 0.41; cellCoords(1,3,1) = 0.58;
      //   cellCoords(2,0,0) = 0.0; cellCoords(2,0,1) = 0.5;
      //   cellCoords(2,1,0) = 0.41; cellCoords(2,1,1) = 0.58;
      //   cellCoords(2,2,0) = 0.5; cellCoords(2,2,1) = 1.0;
      //   cellCoords(2,3,0) = 0.0; cellCoords(2,3,1) = 1.0;
      //   cellCoords(3,0,0) = 0.41; cellCoords(3,0,1) = 0.58;
      //   cellCoords(3,1,0) = 1.0; cellCoords(3,1,1) = 0.5;
      //   cellCoords(3,2,0) = 1.0; cellCoords(3,2,1) = 1.0;
      //   cellCoords(3,3,0) = 0.5; cellCoords(3,3,1) = 1.0;

      //   // get cubature points and weights 
      //   int numCells = 4; int spaceDim = 2;
      //   Intrepid2::FieldContainer<double> cubPoints(numCells,numPoints,spaceDim);
      //   Intrepid2::FieldContainer<double> cubWeights(numCells,numPoints);
      //   CubatureControlVolume.getCubature(cubPoints,cubWeights,cellCoords);

      //   Intrepid2::FieldContainer<double> sideCubPoints(numCells,numSidePoints,spaceDim);
      //   Intrepid2::FieldContainer<double> sideCubWeights(numCells,numSidePoints,spaceDim);
      //   controlVolSideCub->getCubature(sideCubPoints,sideCubWeights,cellCoords);

      //   // test cubature rule by checking the equality of  \int_C \nabla \cdot F dV = \int_dC F \cdot n dS
      //   // using F = (a x,b y), (\nabla \cdot F) = a + b.

      //   // first evaluate F at all control volume side points
      //   double a = 2.1; double b = 1.4;
      //   Intrepid2::FieldContainer<double> F(numCells,numSidePoints,spaceDim);
      //   for (int i = 0; i < numCells; i++) {
      //      for (int j = 0; j < numSidePoints; j++) {
      //         F(i,j,0) = a*sideCubPoints(i,j,0);
      //         F(i,j,1) = b*sideCubPoints(i,j,1);
      //      }
      //   }

      //   // gather the correct contributions to the surface integral
      //   double surfaceInt = 0.0;
    
      //   // contributions from first cell 
      //    surfaceInt += - F(0,1,0)*sideCubWeights(0,1,0) - F(0,1,1)*sideCubWeights(0,1,1)
      //                  + F(0,2,0)*sideCubWeights(0,2,0) + F(0,2,1)*sideCubWeights(0,2,1);
    
      //   // contributions from second cell 
      //    surfaceInt += - F(1,2,0)*sideCubWeights(1,2,0) - F(1,2,1)*sideCubWeights(1,2,1)
      //                  + F(1,3,0)*sideCubWeights(1,3,0) + F(1,3,1)*sideCubWeights(1,3,1);

      //   // contributions from third cell 
      //    surfaceInt += - F(2,0,0)*sideCubWeights(2,0,0) - F(2,0,1)*sideCubWeights(2,0,1)
      //                  + F(2,1,0)*sideCubWeights(2,1,0) + F(2,1,1)*sideCubWeights(2,1,1);

      //   // contributions from fourth cell 
      //    surfaceInt += - F(3,3,0)*sideCubWeights(3,3,0) - F(3,3,1)*sideCubWeights(3,3,1)
      //                  + F(3,0,0)*sideCubWeights(3,0,0) + F(3,0,1)*sideCubWeights(3,0,1);

      //   // gather the correct contributions to the volume integral
      //   double volumeInt = 0.0;

      //   // contributions from first cell 
      //     volumeInt += (a+b)*cubWeights(0,2);

      //   // contributions from second cell 
      //     volumeInt += (a+b)*cubWeights(1,3);

      //   // contributions from third cell 
      //     volumeInt += (a+b)*cubWeights(2,1);

      //   // contributions from fourth cell 
      //     volumeInt += (a+b)*cubWeights(3,0);

      //   if (std::abs(surfaceInt - volumeInt) > Intrepid2::INTREPID_TOL) {
      //       errorFlag++;
      //       *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      //       *outStream << " Integral of (F cdot n) over surface : ";
      //       *outStream << surfaceInt;
      //       *outStream << " does not equal integral of div F over volume: " << volumeInt << "\n";
      //   }

      // }
      // catch (std::exception err) {
      //   *outStream << err.what() << "\n";
      //   errorFlag = -1;
      // }  

      // // triangle primary cells
      // try {

      //   // set cell topology
      //   cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Triangle<> >()));

      //   // define control volume cubature rule
      //   controlVolCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
      //   int numPoints = CubatureControlVolume.getNumPoints();

      //   // define control volume cubature rule
      //   controlVolSideCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolumeSide<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
      //   int numSidePoints = controlVolSideCub->getNumPoints();

      //   // define triangle coordinates - four cells that define a complete control volume around center node
      //   Intrepid2::FieldContainer<double> cellCoords(4,3,2);
      //   cellCoords(0,0,0) = 0.0; cellCoords(0,0,1) = 0.0;
      //   cellCoords(0,1,0) = 1.0; cellCoords(0,1,1) = 0.0;
      //   cellCoords(0,2,0) = 0.41; cellCoords(0,2,1) = 0.58;
      //   cellCoords(1,0,0) = 1.0; cellCoords(1,0,1) = 0.0;
      //   cellCoords(1,1,0) = 1.0; cellCoords(1,1,1) = 1.0;
      //   cellCoords(1,2,0) = 0.41; cellCoords(1,2,1) = 0.58;
      //   cellCoords(2,0,0) = 1.0; cellCoords(2,0,1) = 1.0;
      //   cellCoords(2,1,0) = 0.0; cellCoords(2,1,1) = 1.0;
      //   cellCoords(2,2,0) = 0.41; cellCoords(2,2,1) = 0.58;
      //   cellCoords(3,0,0) = 0.0; cellCoords(3,0,1) = 1.0;
      //   cellCoords(3,1,0) = 0.0; cellCoords(3,1,1) = 0.0;
      //   cellCoords(3,2,0) = 0.41; cellCoords(3,2,1) = 0.58;

      //   // get cubature points and weights 
      //   int numCells = 4; int spaceDim = 2;
      //   Intrepid2::FieldContainer<double> cubPoints(numCells,numPoints,spaceDim);
      //   Intrepid2::FieldContainer<double> cubWeights(numCells,numPoints);
      //   CubatureControlVolume.getCubature(cubPoints,cubWeights,cellCoords);

      //   Intrepid2::FieldContainer<double> sideCubPoints(numCells,numSidePoints,spaceDim);
      //   Intrepid2::FieldContainer<double> sideCubWeights(numCells,numSidePoints,spaceDim);
      //   controlVolSideCub->getCubature(sideCubPoints,sideCubWeights,cellCoords);

      //   // test cubature rule by checking the equality of  \int_C \nabla \cdot F dV = \int_dC F \cdot n dS
      //   // using F = (a x,b y), (\nabla \cdot F) = a + b.

      //   // first evaluate F at all control volume side points
      //   double a = 2.1; double b = 1.4;
      //   Intrepid2::FieldContainer<double> F(numCells,numSidePoints,spaceDim);
      //   for (int i = 0; i < numCells; i++) {
      //      for (int j = 0; j < numSidePoints; j++) {
      //         F(i,j,0) = a*sideCubPoints(i,j,0);
      //         F(i,j,1) = b*sideCubPoints(i,j,1);
      //      }
      //   }

      //   // gather the correct contributions to the surface integral
      //   double surfaceInt = 0.0;
    
      //   // contributions from first cell 
      //    surfaceInt += - F(0,1,0)*sideCubWeights(0,1,0) - F(0,1,1)*sideCubWeights(0,1,1)
      //                  + F(0,2,0)*sideCubWeights(0,2,0) + F(0,2,1)*sideCubWeights(0,2,1);
    
      //   // contributions from second cell 
      //    surfaceInt += - F(1,1,0)*sideCubWeights(1,1,0) - F(1,1,1)*sideCubWeights(1,1,1)
      //                  + F(1,2,0)*sideCubWeights(1,2,0) + F(1,2,1)*sideCubWeights(1,2,1);

      //   // contributions from third cell 
      //    surfaceInt += - F(2,1,0)*sideCubWeights(2,1,0) - F(2,1,1)*sideCubWeights(2,1,1)
      //                  + F(2,2,0)*sideCubWeights(2,2,0) + F(2,2,1)*sideCubWeights(2,2,1);

      //   // contributions from fourth cell 
      //    surfaceInt += - F(3,1,0)*sideCubWeights(3,1,0) - F(3,1,1)*sideCubWeights(3,1,1)
      //                  + F(3,2,0)*sideCubWeights(3,2,0) + F(3,2,1)*sideCubWeights(3,2,1);

      //   // gather the correct contributions to the volume integral
      //   double volumeInt = 0.0;

      //   // contributions from first cell 
      //     volumeInt += (a+b)*cubWeights(0,2);

      //   // contributions from second cell 
      //     volumeInt += (a+b)*cubWeights(1,2);

      //   // contributions from third cell 
      //     volumeInt += (a+b)*cubWeights(2,2);

      //   // contributions from fourth cell 
      //     volumeInt += (a+b)*cubWeights(3,2);

      //   if (std::abs(surfaceInt - volumeInt) > Intrepid2::INTREPID_TOL) {
      //       errorFlag++;
      //       *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      //       *outStream << " Integral of (F cdot n) over surface : ";
      //       *outStream << surfaceInt;
      //       *outStream << " does not equal integral of div F over volume: " << volumeInt << "\n";
      //   }

      // }
      // catch (std::exception err) {
      //   *outStream << err.what() << "\n";
      //   errorFlag = -1;
      // }  


      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
      else
        std::cout << "End Result: TEST PASSED\n";

      // reset format state of std::cout
      std::cout.copyfmt(oldFormatState);
      return errorFlag;
    }
  }
}

  
