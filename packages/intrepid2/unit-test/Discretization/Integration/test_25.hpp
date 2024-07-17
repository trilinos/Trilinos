// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "Intrepid2_CubatureControlVolumeSide.hpp"
#include "Intrepid2_CubatureControlVolumeBoundary.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {

  namespace Test {
#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    } catch (std::logic_error &err) {                                    \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    };      

#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)
    
    template<typename ValueType,
             typename DeviceSpaceType,
             typename cubatureType,
             typename cellCoordViewType,
             typename exactValueType,
             typename outStreamType>
    int Integration_Test25_CompareExactValues(const cubatureType cubature,
                                              const cellCoordViewType cellCoords,
                                              const ordinal_type   exactNumPoints,
                                              const exactValueType *exactCubPoints,
                                              const exactValueType *exactCubWeights,
                                              outStreamType outStream) {

      typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;
      
      const auto tol = tolerence();

      int errorFlag = 0;
      try {
        const ordinal_type numCells = cellCoords.extent(0);
        //const auto numVerts = cellCoords.extent(1);
        const ordinal_type spaceDim = cellCoords.extent(2);
        
        const bool isCubatureSide = std::string(cubature.getName()) == std::string("CubatureControlVolumeSide");

        *outStream << " --- Testing cubature values and weights: " << cubature.getName() << "\n";
        
        const ordinal_type numPoints = cubature.getNumPoints();
        if (exactNumPoints > 0) {
          if (numPoints != exactNumPoints) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "Number of cubature points (" << numPoints 
                       << ") does not match to the correct # of points (4)\n";
          }
        }
        
        // get cubature points and weights for volume integration over control volume
        DynRankView cubPoints, cubWeights;
        cubPoints  = DynRankView("Integration::Test25::cubPoints",  numCells, numPoints, spaceDim);

        if (isCubatureSide) 
          cubWeights = DynRankView("Integration::Test25::cubWeights", numCells, numPoints, spaceDim);
        else 
          cubWeights = DynRankView("Integration::Test25::cubWeights", numCells, numPoints);
          
        cubature.getCubature(cubPoints, cubWeights, cellCoords);

        auto cubPointsHost  = Kokkos::create_mirror_view(cubPoints);
        auto cubWeightsHost = Kokkos::create_mirror_view(cubWeights);

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

              // check values of cubature weight for side cubature
              if (isCubatureSide) {
                if (std::isnan(cubWeightsHost(i,j,k)) || 
                    std::abs(cubWeightsHost(i,j,k) - exactCubWeights[p]) > tol) {
                  errorFlag++;
                  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                  *outStream << " At multi-index { ";
                  *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                  *outStream << "}  computed weight: " << cubWeightsHost(i,j,k)
                             << " but reference value: " << exactCubWeights[p] << "\n";
                }
              }
            }
            
            // check values of cubature weights
            if (!isCubatureSide) {
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
          }
      } catch (std::exception &err) {
        *outStream << err.what() << "\n";
        errorFlag = -10000;
      }  
      return errorFlag;
    }

    
    template<typename ValueType,
             typename DeviceSpaceType,
             typename cellCoordViewType,
             typename exactValueType,
             typename outStreamType>
    int Integration_Test25_CompareVolume(const shards::CellTopology cellTopo,
                                         const cellCoordViewType cellCoords,
                                         const exactValueType exactVolume,
                                         outStreamType outStream) {

      typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;
      
      typedef ValueType pointValueType;
      typedef ValueType weightValueType;
      typedef CubatureControlVolume<DeviceSpaceType,pointValueType,weightValueType> CubatureControlVolumeType;
      
      const auto tol = tolerence();

      int errorFlag = 0;
      try {
        const ordinal_type numCells = cellCoords.extent(0);
        //const ordinal_type numVerts = cellCoords.extent(1);
        const ordinal_type spaceDim = cellCoords.extent(2);
        
        CubatureControlVolumeType CubatureControlVolume(cellTopo);
        *outStream << " --- Testing sum of subcontrol volumes: " << CubatureControlVolume.getName() << "\n";
        
        const ordinal_type numPoints = CubatureControlVolume.getNumPoints();

        // get cubature points and weights for volume integration over control volume
        DynRankView cubPoints, cubWeights;
        cubPoints  = DynRankView("Integration::Test25::cubPoints",  numCells, numPoints, spaceDim);
        cubWeights = DynRankView("Integration::Test25::cubWeights", numCells, numPoints);
          
        CubatureControlVolume.getCubature(cubPoints, cubWeights, cellCoords);
        
        auto cubPointsHost  = Kokkos::create_mirror_view(cubPoints);
        auto cubWeightsHost = Kokkos::create_mirror_view(cubWeights);

        Kokkos::deep_copy(cubPointsHost,  cubPoints);
        Kokkos::deep_copy(cubWeightsHost, cubWeights);
        
        // validation
        ValueType volume = 0;
        for (ordinal_type i=0;i<numPoints;++i) {
          volume += cubWeightsHost(0,i);
        }
        if (std::abs(volume - exactVolume) > tol) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Sum of sub-control volume: ";
          *outStream << volume;
          *outStream << " does not equal primary cell volume: " << exactVolume << "\n";
        }

      } catch (std::exception &err) {
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
      typedef Kokkos::DynRankView<ValueType,Kokkos::HostSpace>   DynRankViewHost;

      int errorFlag = 0;

      try {
        *outStream << " - Quadrilateral cell testing \n";
        {
          // quadrilateral primary cells
          const ordinal_type numCells = 2, numVerts = 4, spaceDim = 2, numEdges = 4;
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
          
          // side points 
          const ValueType exactSideCubPoints[] = {
            0.25, 0.25, 0.375, 0.5,
            0.25, 0.75, 0.125, 0.5,
            0.75, 0.25, 0.875, 0.5,
            0.75, 0.75, 0.625, 0.5
          };
          
          // side weights (these are weighted normals!)
          const ValueType exactSideCubWeights[] = {
            0.5, 0.0, 0.0, 0.25,
            -0.5, 0.0, 0.0,-0.25,
            0.5, 0.0, 0.0, 0.25,
            -0.5, 0.0, 0.0,-0.25
          };
          
          // boundary points 
          const ValueType exactBCCubPoints[] = {
            0.375, 1.0, 0.125, 1.0,
            0.875, 1.0, 0.625, 1.0
          };
          
          // boundary weights 
          const ValueType exactBCCubWeights[] = {
            0.25, 0.25, 0.25, 0.25
          };

          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<> >());
          {
            CubatureControlVolume<DeviceSpaceType,ValueType,ValueType> cubature(cellTopo);
            errorFlag += Integration_Test25_CompareExactValues
              <ValueType,DeviceSpaceType>(cubature,
                                          cellCoords,
                                          numVerts,
                                          exactCubPoints,
                                          exactCubWeights,
                                          outStream);
          }
          {
            CubatureControlVolumeSide<DeviceSpaceType,ValueType,ValueType> cubature(cellTopo);
            errorFlag += Integration_Test25_CompareExactValues
              <ValueType,DeviceSpaceType>(cubature,
                                          cellCoords,
                                          numEdges,
                                          exactSideCubPoints,
                                          exactSideCubWeights,
                                          outStream);
          }
          {
            const ordinal_type side = 2;
            const ordinal_type numBoundary = 2;
            CubatureControlVolumeBoundary<DeviceSpaceType,ValueType,ValueType> cubature(cellTopo, side);
            errorFlag += Integration_Test25_CompareExactValues
              <ValueType,DeviceSpaceType>(cubature,
                                          cellCoords,
                                          numBoundary,
                                          exactBCCubPoints,
                                          exactBCCubWeights,
                                          outStream);
          }
        }
        
        *outStream << " - Triangle cell testing \n";
        {
          // triangle primary cells      
          const ordinal_type numCells = 2, numVerts = 3, spaceDim = 2, numEdges = 3; 
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
          
          // side points 
          const ValueType exactSideCubPoints[] = {
            0.2083333333333333, 0.0833333333333333, 0.2083333333333333, 0.2083333333333333,
            0.0833333333333333, 0.2083333333333333, 0.4166666666666667, 0.2916666666666667,
            0.2916666666666667, 0.4166666666666667, 0.2916666666666667, 0.2916666666666667
          };
          
           // side weights (these are weighted normals!)
          const ValueType exactSideCubWeights[] = {
            0.1666666666666667, 0.0833333333333333,-0.0833333333333333, 0.0833333333333333,
           -0.0833333333333333,-0.1666666666666667, 0.0833333333333333, 0.1666666666666667,
           -0.1666666666666667,-0.0833333333333333, 0.0833333333333333,-0.0833333333333333
          };
          
          // boundary points 
          const ValueType exactBCCubPoints[] = {
            0.375, 0.125, 0.125, 0.375,
            0.375, 0.5, 0.125, 0.5
          };
          
          // boundary weights 
          const ValueType exactBCCubWeights[] = {
            0.353553390593274, 0.353553390593274, 0.25, 0.25
          };
          
          // triangle primary cells      
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData< shards::Triangle<> >());
          {
            CubatureControlVolume<DeviceSpaceType,ValueType,ValueType> cubature(cellTopo);
            errorFlag += Integration_Test25_CompareExactValues
              <ValueType,DeviceSpaceType>(cubature,
                                          cellCoords,
                                          numVerts,
                                          exactCubPoints,
                                          exactCubWeights,
                                          outStream);
          }
          {
            CubatureControlVolumeSide<DeviceSpaceType,ValueType,ValueType> cubature(cellTopo);
            errorFlag += Integration_Test25_CompareExactValues
              <ValueType,DeviceSpaceType>(cubature,
                                          cellCoords,
                                          numEdges,
                                          exactSideCubPoints,
                                          exactSideCubWeights,
                                          outStream);
          }
          {
            const ordinal_type side = 1;
            const ordinal_type numBoundary = 2;
            CubatureControlVolumeBoundary<DeviceSpaceType,ValueType,ValueType> cubature(cellTopo, side);
            errorFlag += Integration_Test25_CompareExactValues
              <ValueType,DeviceSpaceType>(cubature,
                                          cellCoords,
                                          numBoundary,
                                          exactBCCubPoints,
                                          exactBCCubWeights,
                                          outStream);
          }
        }
        *outStream << " - Tetrahedron cell testing \n";
        {
          // tetrahedron primary cells      
          const ordinal_type numCells = 1, numVerts = 4, spaceDim = 3, numEdges = 6;
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

          // side points 
          const ValueType exactSideCubPoints[] = {
            0.3541666666666667, 0.1458333333333333, 0.1458333333333333,
            0.3541666666666667, 0.3541666666666667, 0.1458333333333333,
            0.1458333333333333, 0.3541666666666667, 0.1458333333333333,
            0.1458333333333333, 0.1458333333333333, 0.3541666666666667,
            0.3541666666666667, 0.1458333333333333, 0.3541666666666667,
            0.1458333333333333, 0.3541666666666667, 0.3541666666666667
          };

          // side weights (these are weighted normals!)
          const ValueType exactSideCubWeights[] = {
            0.0833333333333333, 0.0416666666666667, 0.041666666666667,
            -0.0416666666666667, 0.0416666666666667, 0.000000000000000,
            -0.0416666666666667,-0.0833333333333333,-0.041666666666667,
            0.0416666666666667, 0.0416666666666667, 0.083333333333333,
            -0.0416666666666667, 0.0000000000000000, 0.041666666666667,
            0.0000000000000000,-0.0416666666666667, 0.041666666666667
          };

          // boundary points 
          const ValueType exactBCCubPoints[] = {
            0.208333333333333, 0.00, 0.208333333333333,
            0.583333333333333, 0.00, 0.208333333333333, 
            0.208333333333333, 0.00, 0.583333333333333, 
          };

          // boundary weights 
          const ValueType exactBCCubWeights[] = {
            0.166666666666667, 0.166666666666667, 0.166666666666667
          };

          // tetrahedral primary cells
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<> >());
          {
            CubatureControlVolume<DeviceSpaceType,ValueType,ValueType> cubature(cellTopo);
            errorFlag += Integration_Test25_CompareExactValues
              <ValueType,DeviceSpaceType>(cubature,
                                          cellCoords,
                                          numVerts,
                                          exactCubPoints,
                                          exactCubWeights,
                                          outStream);
          }
          {
            CubatureControlVolumeSide<DeviceSpaceType,ValueType,ValueType> cubature(cellTopo);
            errorFlag += Integration_Test25_CompareExactValues
              <ValueType,DeviceSpaceType>(cubature,
                                          cellCoords,
                                          numEdges,
                                          exactSideCubPoints,
                                          exactSideCubWeights,
                                          outStream);
          }
          {
            const ordinal_type side = 0;
            const ordinal_type numBoundary = 3;
            CubatureControlVolumeBoundary<DeviceSpaceType,ValueType,ValueType> cubature(cellTopo, side);
            errorFlag += Integration_Test25_CompareExactValues
              <ValueType,DeviceSpaceType>(cubature,
                                          cellCoords,
                                          numBoundary,
                                          exactBCCubPoints,
                                          exactBCCubWeights,
                                          outStream);
          }
         }
        *outStream << " - Hexahedron cell testing \n";
        {
          const ordinal_type numCells = 1, numVerts = 8, spaceDim = 3, numEdges = 12; 
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

          // side points 
          const ValueType exactSideCubPoints[] = {
            1.0, 0.375, 0.25, 1.5, 0.750, 0.25,
            1.0, 1.125, 0.25, 0.5, 0.750, 0.25, 
            1.0, 0.375, 0.75, 1.5, 0.750, 0.75, 
            1.0, 1.125, 0.75, 0.5, 0.750, 0.75,
            0.5, 0.375, 0.50, 1.5, 0.375, 0.50,
            1.5, 1.125, 0.50, 0.5, 1.125, 0.50
          };

          // side weights (these are weighted normals!)
          const ValueType exactSideCubWeights[] = {
            0.375, 0.00, 0.00, 0.00, 0.50, 0.00,
            -0.375, 0.00, 0.00, 0.00,-0.50, 0.00,
            0.375, 0.00, 0.00, 0.00, 0.50, 0.00,
            -0.375, 0.00, 0.00, 0.00,-0.50, 0.00,
            0.000, 0.00, 0.75, 0.00, 0.00, 0.75,
            0.000, 0.00, 0.75, 0.00, 0.00, 0.75,
          };

          // boundary points 
          const ValueType exactBCCubPoints[] = {
            0.5, 0.00, 0.25,
            1.5, 0.00, 0.25,
            1.5, 0.00, 0.75,
            0.5, 0.00, 0.75,
          };

          // boundary weights 
          const ValueType exactBCCubWeights[] = {
            0.5, 0.5, 0.5, 0.5
          };

          // hexahedral primary cells
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<> >());
          {
            CubatureControlVolume<DeviceSpaceType,ValueType,ValueType> cubature(cellTopo);
            errorFlag += Integration_Test25_CompareExactValues
              <ValueType,DeviceSpaceType>(cubature,
                                          cellCoords,
                                          numVerts,
                                          exactCubPoints,
                                          exactCubWeights,
                                          outStream);
          }
          {
            CubatureControlVolumeSide<DeviceSpaceType,ValueType,ValueType> cubature(cellTopo);
            errorFlag += Integration_Test25_CompareExactValues
              <ValueType,DeviceSpaceType>(cubature,
                                          cellCoords,
                                          numEdges,
                                          exactSideCubPoints,
                                          exactSideCubWeights,
                                          outStream);
          }
          {
            const ordinal_type side = 0;
            const ordinal_type numBoundary = 4;
            CubatureControlVolumeBoundary<DeviceSpaceType,ValueType,ValueType> cubature(cellTopo, side);
            errorFlag += Integration_Test25_CompareExactValues
              <ValueType,DeviceSpaceType>(cubature,
                                          cellCoords,
                                          numBoundary,
                                          exactBCCubPoints,
                                          exactBCCubWeights,
                                          outStream);
          }
        }
      } catch (std::exception &err) {
        *outStream << err.what() << "\n";
        errorFlag = -10000;
      }  

      *outStream 
        << "===============================================================================\n" 
        << "| TEST 2: comparison of sub-control volume weights and primary cell volume    |\n"
        << "===============================================================================\n";
      
      // quadrilateral primary cells
      try {
        *outStream << " - Quadrilateral cell testing \n";
        {
          // set cell topology
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<> >());
                    
          // quadrilateral primary cells
          const ordinal_type numCells = 1, numVerts = 4, spaceDim = 2;
          DynRankViewHost ConstructWithLabel(cellCoordsHost, numCells, numVerts, spaceDim);

          cellCoordsHost(0,0,0) = 0.0; cellCoordsHost(0,0,1) = 0.0;
          cellCoordsHost(0,1,0) = 2.4; cellCoordsHost(0,1,1) = 0.0;
          cellCoordsHost(0,2,0) = 2.4; cellCoordsHost(0,2,1) = 3.1;
          cellCoordsHost(0,3,0) = 0.0; cellCoordsHost(0,3,1) = 3.1;
          
          auto cellCoords = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), cellCoordsHost);
          Kokkos::deep_copy(cellCoords, cellCoordsHost);
          
          const ValueType exactVolume = 2.4*3.1;
          errorFlag += Integration_Test25_CompareVolume
            <ValueType,DeviceSpaceType>(cellTopo,
                                        cellCoords,
                                        exactVolume,
                                        outStream);
        }
        *outStream << " - Triangle cell testing \n";
        {
          // set cell topology
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData< shards::Triangle<> >());

          // quadrilateral primary cells
          const ordinal_type numCells = 1, numVerts = 3, spaceDim = 2;
          DynRankViewHost ConstructWithLabel(cellCoordsHost, numCells, numVerts, spaceDim);

          cellCoordsHost(0,0,0) = 0.0; cellCoordsHost(0,0,1) = 0.0;
          cellCoordsHost(0,1,0) = 3.6; cellCoordsHost(0,1,1) = 0.0;
          cellCoordsHost(0,2,0) = 0.0; cellCoordsHost(0,2,1) = 2.8;
          
          auto cellCoords = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), cellCoordsHost);
          Kokkos::deep_copy(cellCoords, cellCoordsHost);
          
          const ValueType exactVolume = 0.5*3.6*2.8;
          errorFlag += Integration_Test25_CompareVolume
            <ValueType,DeviceSpaceType>(cellTopo,
                                        cellCoords,
                                        exactVolume,
                                        outStream);
        }
        *outStream << " - Tetrahedron cell testing \n";
        {
          // set cell topology
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<> >());

          const ordinal_type numCells = 1, numVerts = 4, spaceDim = 3;
          DynRankViewHost ConstructWithLabel(cellCoordsHost, numCells, numVerts, spaceDim);

          cellCoordsHost(0,0,0) = 0.0; cellCoordsHost(0,0,1) = 0.0; cellCoordsHost(0,0,2) = 0.0;
          cellCoordsHost(0,1,0) = 3.6; cellCoordsHost(0,1,1) = 0.0; cellCoordsHost(0,1,2) = 0.0;
          cellCoordsHost(0,2,0) = 0.0; cellCoordsHost(0,2,1) = 2.8; cellCoordsHost(0,2,2) = 0.0;
          cellCoordsHost(0,3,0) = 0.0; cellCoordsHost(0,3,1) = 2.8; cellCoordsHost(0,3,2) = 1.7;
          
          auto cellCoords = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), cellCoordsHost);
          Kokkos::deep_copy(cellCoords, cellCoordsHost);

          const ValueType exactVolume = (0.5*3.6*2.8)*1.7/3.0;
          errorFlag += Integration_Test25_CompareVolume
            <ValueType,DeviceSpaceType>(cellTopo,
                                        cellCoords,
                                        exactVolume,
                                        outStream);
        }
        *outStream << " - Hexahedron cell testing \n";        
        {
          // set cell topology
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<> >());
          
          const ordinal_type numCells = 1, numVerts = 8, spaceDim = 3;
          DynRankViewHost ConstructWithLabel(cellCoordsHost, numCells, numVerts, spaceDim);
          
          cellCoordsHost(0,0,0) = 0.0; cellCoordsHost(0,0,1) = 0.0; cellCoordsHost(0,0,2) = 0.0;
          cellCoordsHost(0,1,0) = 2.4; cellCoordsHost(0,1,1) = 0.0; cellCoordsHost(0,1,2) = 0.0;
          cellCoordsHost(0,2,0) = 2.4; cellCoordsHost(0,2,1) = 3.1; cellCoordsHost(0,2,2) = 0.0;
          cellCoordsHost(0,3,0) = 0.0; cellCoordsHost(0,3,1) = 3.1; cellCoordsHost(0,3,2) = 0.0;
          cellCoordsHost(0,4,0) = 0.0; cellCoordsHost(0,4,1) = 0.0; cellCoordsHost(0,4,2) = 1.7;
          cellCoordsHost(0,5,0) = 2.4; cellCoordsHost(0,5,1) = 0.0; cellCoordsHost(0,5,2) = 1.7;
          cellCoordsHost(0,6,0) = 2.4; cellCoordsHost(0,6,1) = 3.1; cellCoordsHost(0,6,2) = 1.7;
          cellCoordsHost(0,7,0) = 0.0; cellCoordsHost(0,7,1) = 3.1; cellCoordsHost(0,7,2) = 1.7;
          
          auto cellCoords = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), cellCoordsHost);
          Kokkos::deep_copy(cellCoords, cellCoordsHost);
          
          const ValueType exactVolume = 2.4*3.1*1.7;
          errorFlag += Integration_Test25_CompareVolume
            <ValueType,DeviceSpaceType>(cellTopo,
                                        cellCoords,
                                        exactVolume,
                                        outStream);
        }
      } catch (std::exception &err) {
        *outStream << err.what() << "\n";
        errorFlag = -1;
      }  
      
      *outStream                                                     
        << "===============================================================================\n"
        << "| TEST 3: control volume integration                                          |\n"
        << "===============================================================================\n";
      
      // quadrilateral primary cells
      try {
        typedef ValueType pointValueType;
        typedef ValueType weightValueType;
        typedef CubatureControlVolume<DeviceSpaceType,pointValueType,weightValueType> CubatureControlVolumeType;
        typedef CubatureControlVolumeSide<DeviceSpaceType,pointValueType,weightValueType> CubatureControlVolumeSideType;
        
        const auto tol = tolerence();

        *outStream << " - Quadrilateral cell testing \n";
        {
          // set cell topology
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<> >());
                    
          // quadrilateral primary cells
          const ordinal_type numCells = 4, numVerts = 4, spaceDim = 2;
          DynRankViewHost ConstructWithLabel(cellCoordsHost, numCells, numVerts, spaceDim);

          // define quad coordinates - four cells that define a complete control volume around center node
          cellCoordsHost(0,0,0) = 0.0; cellCoordsHost(0,0,1) = 0.0;
          cellCoordsHost(0,1,0) = 0.5; cellCoordsHost(0,1,1) = 0.0;
          cellCoordsHost(0,2,0) = 0.41; cellCoordsHost(0,2,1) = 0.58;
          cellCoordsHost(0,3,0) = 0.0; cellCoordsHost(0,3,1) = 0.5;
          cellCoordsHost(1,0,0) = 0.5; cellCoordsHost(1,0,1) = 0.0;
          cellCoordsHost(1,1,0) = 1.0; cellCoordsHost(1,1,1) = 0.0;
          cellCoordsHost(1,2,0) = 1.0; cellCoordsHost(1,2,1) = 0.5;
          cellCoordsHost(1,3,0) = 0.41; cellCoordsHost(1,3,1) = 0.58;
          cellCoordsHost(2,0,0) = 0.0; cellCoordsHost(2,0,1) = 0.5;
          cellCoordsHost(2,1,0) = 0.41; cellCoordsHost(2,1,1) = 0.58;
          cellCoordsHost(2,2,0) = 0.5; cellCoordsHost(2,2,1) = 1.0;
          cellCoordsHost(2,3,0) = 0.0; cellCoordsHost(2,3,1) = 1.0;
          cellCoordsHost(3,0,0) = 0.41; cellCoordsHost(3,0,1) = 0.58;
          cellCoordsHost(3,1,0) = 1.0; cellCoordsHost(3,1,1) = 0.5;
          cellCoordsHost(3,2,0) = 1.0; cellCoordsHost(3,2,1) = 1.0;
          cellCoordsHost(3,3,0) = 0.5; cellCoordsHost(3,3,1) = 1.0;
          
          auto cellCoords = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), cellCoordsHost);
          Kokkos::deep_copy(cellCoords, cellCoordsHost);
          
          CubatureControlVolumeType     CubatureControlVolume(cellTopo);
          CubatureControlVolumeSideType CubatureControlVolumeSide(cellTopo);
          *outStream << " --- Testing subcontrol volume and side integration: " 
                     << CubatureControlVolume.getName() << ", " 
                     << CubatureControlVolumeSide.getName() << "\n";
          
          const ordinal_type numPointsVolume = CubatureControlVolume.getNumPoints();
          const ordinal_type numPointsSide   = CubatureControlVolumeSide.getNumPoints();
          
          // get cubature points and weights for volume integration over control volume
          DynRankView cubPointsVolume, cubWeightsVolume, cubPointsSide, cubWeightsSide;
          cubPointsVolume  = DynRankView("Integration::Test25::cubPointsVolume",  numCells, numPointsVolume, spaceDim);
          cubWeightsVolume = DynRankView("Integration::Test25::cubWeightsVolume", numCells, numPointsVolume);
          cubPointsSide    = DynRankView("Integration::Test25::cubPointsSide",    numCells, numPointsSide,   spaceDim);
          cubWeightsSide   = DynRankView("Integration::Test25::cubWeights",       numCells, numPointsSide,   spaceDim);
          
          CubatureControlVolume.getCubature    (cubPointsVolume, cubWeightsVolume, cellCoords);
          CubatureControlVolumeSide.getCubature(cubPointsSide,   cubWeightsSide,   cellCoords);
          
          auto cubPointsVolumeHost  = Kokkos::create_mirror_view(cubPointsVolume);
          auto cubWeightsVolumeHost = Kokkos::create_mirror_view(cubWeightsVolume);
          auto cubPointsSideHost    = Kokkos::create_mirror_view(cubPointsSide);
          auto cubWeightsSideHost   = Kokkos::create_mirror_view(cubWeightsSide);
          
          Kokkos::deep_copy(cubPointsVolumeHost,  cubPointsVolume);
          Kokkos::deep_copy(cubWeightsVolumeHost, cubWeightsVolume);
          Kokkos::deep_copy(cubPointsSideHost,    cubPointsSide);
          Kokkos::deep_copy(cubWeightsSideHost,   cubWeightsSide);
          
          // test cubature rule by checking the equality of  \int_C \nabla \cdot F dV = \int_dC F \cdot n dS
          // using F = (a x,b y), (\nabla \cdot F) = a + b.
          
          // first evaluate F at all control volume side points
          const ValueType a = 2.1, b = 1.4;
          DynRankViewHost F("Integration::Test25::F", numCells, numPointsSide, spaceDim);
          for (ordinal_type i=0;i<numCells;++i)
            for (ordinal_type j=0;j<numPointsSide;++j) {
              F(i,j,0) = a*cubPointsSide(i,j,0);
              F(i,j,1) = b*cubPointsSide(i,j,1);
            }
          
          // gather the correct contributions to the surface integral
          ValueType surfaceInt = 0.0;
          
          // contributions from first cell
          surfaceInt += ( - F(0,1,0)*cubWeightsSide(0,1,0) - F(0,1,1)*cubWeightsSide(0,1,1)
                          + F(0,2,0)*cubWeightsSide(0,2,0) + F(0,2,1)*cubWeightsSide(0,2,1) );

          // contributions from second cell
          surfaceInt += ( - F(1,2,0)*cubWeightsSide(1,2,0) - F(1,2,1)*cubWeightsSide(1,2,1)
                          + F(1,3,0)*cubWeightsSide(1,3,0) + F(1,3,1)*cubWeightsSide(1,3,1) );

          // contributions from third cell
          surfaceInt += ( - F(2,0,0)*cubWeightsSide(2,0,0) - F(2,0,1)*cubWeightsSide(2,0,1)
                          + F(2,1,0)*cubWeightsSide(2,1,0) + F(2,1,1)*cubWeightsSide(2,1,1) );

          // contributions from fourth cell
          surfaceInt += ( - F(3,3,0)*cubWeightsSide(3,3,0) - F(3,3,1)*cubWeightsSide(3,3,1)
                          + F(3,0,0)*cubWeightsSide(3,0,0) + F(3,0,1)*cubWeightsSide(3,0,1) );

          // gather the correct contributions to the volume integral
          ValueType volumeInt = 0.0;

          // contributions from first cell
          volumeInt += (a+b)*cubWeightsVolume(0,2);

          // contributions from second cell
          volumeInt += (a+b)*cubWeightsVolume(1,3);

          // contributions from third cell
          volumeInt += (a+b)*cubWeightsVolume(2,1);

          // contributions from fourth cell
          volumeInt += (a+b)*cubWeightsVolume(3,0);

          if (std::abs(surfaceInt - volumeInt) > tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Integral of (F cdot n) over surface : ";
            *outStream << surfaceInt;
            *outStream << " does not equal integral of div F over volume: " << volumeInt << "\n";
          }
        }
        *outStream << " - Triangle cell testing \n";
        {
          // set cell topology
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData< shards::Triangle<> >());
                    
          // quadrilateral primary cells
          const ordinal_type numCells = 4, numVerts = 3, spaceDim = 2;
          DynRankViewHost ConstructWithLabel(cellCoordsHost, numCells, numVerts, spaceDim);

          cellCoordsHost(0,0,0) = 0.0; cellCoordsHost(0,0,1) = 0.0;
          cellCoordsHost(0,1,0) = 1.0; cellCoordsHost(0,1,1) = 0.0;
          cellCoordsHost(0,2,0) = 0.41; cellCoordsHost(0,2,1) = 0.58;
          cellCoordsHost(1,0,0) = 1.0; cellCoordsHost(1,0,1) = 0.0;
          cellCoordsHost(1,1,0) = 1.0; cellCoordsHost(1,1,1) = 1.0;
          cellCoordsHost(1,2,0) = 0.41; cellCoordsHost(1,2,1) = 0.58;
          cellCoordsHost(2,0,0) = 1.0; cellCoordsHost(2,0,1) = 1.0;
          cellCoordsHost(2,1,0) = 0.0; cellCoordsHost(2,1,1) = 1.0;
          cellCoordsHost(2,2,0) = 0.41; cellCoordsHost(2,2,1) = 0.58;
          cellCoordsHost(3,0,0) = 0.0; cellCoordsHost(3,0,1) = 1.0;
          cellCoordsHost(3,1,0) = 0.0; cellCoordsHost(3,1,1) = 0.0;
          cellCoordsHost(3,2,0) = 0.41; cellCoordsHost(3,2,1) = 0.58;
          
          auto cellCoords = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), cellCoordsHost);
          Kokkos::deep_copy(cellCoords, cellCoordsHost);
          
          CubatureControlVolumeType     CubatureControlVolume(cellTopo);
          CubatureControlVolumeSideType CubatureControlVolumeSide(cellTopo);
          *outStream << " --- Testing subcontrol volume and side integration: " 
                     << CubatureControlVolume.getName() << ", " 
                     << CubatureControlVolumeSide.getName() << "\n";
          
          const ordinal_type numPointsVolume = CubatureControlVolume.getNumPoints();
          const ordinal_type numPointsSide   = CubatureControlVolumeSide.getNumPoints();
          
          // get cubature points and weights for volume integration over control volume
          DynRankView cubPointsVolume, cubWeightsVolume, cubPointsSide, cubWeightsSide;
          cubPointsVolume  = DynRankView("Integration::Test25::cubPointsVolume",  numCells, numPointsVolume, spaceDim);
          cubWeightsVolume = DynRankView("Integration::Test25::cubWeightsVolume", numCells, numPointsVolume);
          cubPointsSide    = DynRankView("Integration::Test25::cubPointsSide",    numCells, numPointsSide,   spaceDim);
          cubWeightsSide   = DynRankView("Integration::Test25::cubWeights",       numCells, numPointsSide,   spaceDim);
          
          CubatureControlVolume.getCubature    (cubPointsVolume, cubWeightsVolume, cellCoords);
          CubatureControlVolumeSide.getCubature(cubPointsSide,   cubWeightsSide,   cellCoords);
          
          auto cubPointsVolumeHost  = Kokkos::create_mirror_view(cubPointsVolume);
          auto cubWeightsVolumeHost = Kokkos::create_mirror_view(cubWeightsVolume);
          auto cubPointsSideHost    = Kokkos::create_mirror_view(cubPointsSide);
          auto cubWeightsSideHost   = Kokkos::create_mirror_view(cubWeightsSide);
          
          Kokkos::deep_copy(cubPointsVolumeHost,  cubPointsVolume);
          Kokkos::deep_copy(cubWeightsVolumeHost, cubWeightsVolume);
          Kokkos::deep_copy(cubPointsSideHost,    cubPointsSide);
          Kokkos::deep_copy(cubWeightsSideHost,   cubWeightsSide);
          
          // test cubature rule by checking the equality of  \int_C \nabla \cdot F dV = \int_dC F \cdot n dS
          // using F = (a x,b y), (\nabla \cdot F) = a + b.
          
          // first evaluate F at all control volume side points
          const ValueType a = 2.1, b = 1.4;
          DynRankViewHost F("Integration::Test25::F", numCells, numPointsSide, spaceDim);
          for (ordinal_type i=0;i<numCells;++i)
            for (ordinal_type j=0;j<numPointsSide;++j) {
              F(i,j,0) = a*cubPointsSide(i,j,0);
              F(i,j,1) = b*cubPointsSide(i,j,1);
            }
          
          // gather the correct contributions to the surface integral
          ValueType surfaceInt = 0.0;
          
          // contributions from first cell 
          surfaceInt += ( - F(0,1,0)*cubWeightsSide(0,1,0) - F(0,1,1)*cubWeightsSide(0,1,1)
                          + F(0,2,0)*cubWeightsSide(0,2,0) + F(0,2,1)*cubWeightsSide(0,2,1) );
          
          // contributions from second cell 
          surfaceInt += ( - F(1,1,0)*cubWeightsSide(1,1,0) - F(1,1,1)*cubWeightsSide(1,1,1)
                          + F(1,2,0)*cubWeightsSide(1,2,0) + F(1,2,1)*cubWeightsSide(1,2,1) );
          
          // contributions from third cell 
          surfaceInt += ( - F(2,1,0)*cubWeightsSide(2,1,0) - F(2,1,1)*cubWeightsSide(2,1,1)
                          + F(2,2,0)*cubWeightsSide(2,2,0) + F(2,2,1)*cubWeightsSide(2,2,1) );
          
          // contributions from fourth cell 
          surfaceInt += ( - F(3,1,0)*cubWeightsSide(3,1,0) - F(3,1,1)*cubWeightsSide(3,1,1)
                          + F(3,2,0)*cubWeightsSide(3,2,0) + F(3,2,1)*cubWeightsSide(3,2,1) );
          
          // gather the correct contributions to the volume integral
          ValueType volumeInt = 0.0;
          
          // contributions from first cell 
          volumeInt += (a+b)*cubWeightsVolume(0,2);
          
          // contributions from second cell 
          volumeInt += (a+b)*cubWeightsVolume(1,2);
          
          // contributions from third cell 
          volumeInt += (a+b)*cubWeightsVolume(2,2);
          
          // contributions from fourth cell 
          volumeInt += (a+b)*cubWeightsVolume(3,2);
          
          if (std::abs(surfaceInt - volumeInt) > tol) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Integral of (F cdot n) over surface : ";
            *outStream << surfaceInt;
            *outStream << " does not equal integral of div F over volume: " << volumeInt << "\n";
          }
        }
      } catch (std::exception &err) {
        *outStream << err.what() << "\n";
        errorFlag = -1;
      }          

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
