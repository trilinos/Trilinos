// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Test of the CellTools class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson and Kyungjoo Kim
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_CellTools.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {
  
  namespace Test {
#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      S ;                                                               \
    }                                                                   \
    catch (std::logic_error &err) {                                      \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    };                                                                  
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

    /** \brief  Maps the vertices of the subcell parametrization domain to that subcell. 
        
        Parametrization tests check if the vertices of the parametrization domain are properly 
        mapped to vertices of the resepective reference subcell. Because the parametrization map 
        is a polynomial whose degree depends on the number of vertices, if all vertices are 
        mapped correctly this is sufficient to assert that the parametrization map is correct.
        
        To test reference cells with two different kinds of faces, there are two argument slots
        to pass vertices of parametrization domains "A" and "B". When testing edges always pass
        the same argument because edges have the same parametrization domain [-1,1] (1-cube)
        
        \param  errorFlag       [out] - counts number of errors
        \param  parentCell      [in]  - topology of the reference cell whose 1 and 2-subcells are parametrized
        \param  subcParamVert_A [in]  - vertex coordinates of parametrization domain "A" (2-simplex for faces)
        \param  subcParamVert_A [in]  - vertex coordinates of parametrization domain "B" (2-cube for faces)
        \param  subcDim         [in]  - dimension of the subcells whose parametrizations are tested
        \param  outStream       [in]  - output stream to write
    */
    template<typename ValueType,
             typename DeviceType,
             typename subcParamVertAType,
             typename subcParamVertBType,
             typename tolType,
             typename outStreamPtrType>
    void testSubcellParametrizations( int                        &errorFlag,
                                      const shards::CellTopology  parentCell,
                                      const subcParamVertAType    subcParamVert_A,
                                      const subcParamVertBType    subcParamVert_B,
                                      const int                   subcDim,
                                      const tolType               tol,
                                      const outStreamPtrType      outStreamPtr ) {
      using ct = CellTools<DeviceType>;
      using DynRankView = Kokkos::DynRankView<ValueType,DeviceType>;

      // Get cell dimension and subcell count
      const auto cellDim   = parentCell.getDimension();
      const auto subcCount = parentCell.getSubcellCount(subcDim);

      // Storage for correct reference subcell vertices and for the images of the parametrization domain points
      // set the space large enough .... subcell itself can be a cell
      const ordinal_type maxNodeCount = 256;
      DynRankView ConstructWithLabel(refSubcellNodesMax,  maxNodeCount, cellDim);
      DynRankView ConstructWithLabel(mappedParamVerticesMax, maxNodeCount, cellDim);
      
      // Loop over subcells of the specified dimension
      for(size_type subcOrd=0;subcOrd<subcCount;++subcOrd) {
        const auto subcVertexCount = parentCell.getVertexCount(subcDim, subcOrd);
        const auto subcNodeCount = parentCell.getNodeCount(subcDim, subcOrd);

        const Kokkos::pair<ordinal_type,ordinal_type> nodeRange(0, subcNodeCount);
        const Kokkos::pair<ordinal_type,ordinal_type> vertexRange(0, subcVertexCount);
        
        auto refSubcellNodes  = Kokkos::subdynrankview( refSubcellNodesMax,  nodeRange, Kokkos::ALL() );
        auto mappedParamVertices = Kokkos::subdynrankview( mappedParamVerticesMax, vertexRange, Kokkos::ALL() );
        
        // Retrieve correct reference subcell vertices
        ct::getReferenceSubcellVertices(refSubcellNodes, subcDim, subcOrd, parentCell);

        // Map vertices of the parametrization domain to 1 or 2-subcell with ordinal subcOrd
        switch (subcDim) {
        case 1: {
          // For edges parametrization domain is always 1-cube passed as "subcParamVert_A"
          ct::mapToReferenceSubcell(mappedParamVertices,
                                    subcParamVert_A,
                                    subcDim,
                                    subcOrd,
                                    parentCell);
          break;
        }
        case 2: {
          // For faces need to treat Triangle and Quadrilateral faces separately          
          if (subcVertexCount == 3) // domain "subcParamVert_A" is the standard 2-simplex  
            ct::mapToReferenceSubcell(mappedParamVertices,
                                      subcParamVert_A,
                                      subcDim,
                                      subcOrd,
                                      parentCell);

          else if (subcVertexCount == 4) // Domain "subcParamVert_B" is the standard 2-cube
            ct::mapToReferenceSubcell(mappedParamVertices,
                                      subcParamVert_B,
                                      subcDim,
                                      subcOrd,
                                      parentCell);
          else {
            errorFlag = 1000;
            INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error,
                                          ">>> ERROR (Intrepid2::CellTools::Test01::testSubcellParametrizations): subcell topology is not tri nor quad.");
          }
          break;
        }
        default: {
          errorFlag = 1000;          
          INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error,
                                        ">>> ERROR (Intrepid2::CellTools::Test01::testSubcellParametrizations): subcell dim is not 1 nor 2.");
        }          
        }
        
        //strided subviews from subviews cannot be directly copied to host into a non-strided views
        DynRankView ConstructWithLabel(nonStridedParamVertices, mappedParamVertices.extent(0),mappedParamVertices.extent(1));
        Kokkos::deep_copy(nonStridedParamVertices, mappedParamVertices);
        auto hMappedParamVertices = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), nonStridedParamVertices);

        DynRankView ConstructWithLabel(nonStridedRefNodes, refSubcellNodes.extent(0),refSubcellNodes.extent(1));
        Kokkos::deep_copy(nonStridedRefNodes, refSubcellNodes);
        auto hRefSubcellNodes = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), nonStridedRefNodes);

        // Compare the images of the parametrization domain vertices with the true vertices (test provide vertices only).
        for (size_type subcVertOrd=0;subcVertOrd<subcVertexCount;++subcVertOrd) 
          for (size_type i=0;i<cellDim;++i)
            if (std::abs(hMappedParamVertices(subcVertOrd, i) - hRefSubcellNodes(subcVertOrd, i)) > tol) {
              ++errorFlag; 
              *outStreamPtr 
                << std::setw(70) << "^^^^----FAILURE!" << "\n"
                << " Cell Topology = " << parentCell.getName() << "\n"
                << " Mapped vertex = " << hMappedParamVertices(subcVertOrd, i)
                << " Reference subcell vertex = " << hRefSubcellNodes(subcVertOrd, i) << "\n"
                << " Parametrization of subcell " << subcOrd << " which is "
                << parentCell.getName(subcDim,subcOrd) << " failed for vertex " << subcVertOrd << ":\n"
                << " parametrization map fails to map correctly coordinate " << i << " of that vertex\n\n";
            }
      }
    }
      
      template<typename ValueType, typename DeviceType>
      int CellTools_Test01(const bool verbose) {
        using value_type = ValueType;

        Teuchos::RCP<std::ostream> outStream;
        Teuchos::oblackholestream bhs; // outputs nothing

        if (verbose)
          outStream = Teuchos::rcp(&std::cout, false);
        else
          outStream = Teuchos::rcp(&bhs,       false);

        Teuchos::oblackholestream oldFormatState;
        oldFormatState.copyfmt(std::cout);

        *outStream
          << "===============================================================================\n"
          << "|                                                                             |\n"
          << "|                              Unit Test CellTools                            |\n"
          << "|                                                                             |\n"
          << "|     1) Edge parametrizations                                                |\n"
          << "|     2) Face parametrizations                                                |\n"
          << "|     3) Edge tangents                                                        |\n"
          << "|     4) Face tangents and normals                                            |\n"
          << "|                                                                             |\n"
          << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov)                      |\n"
          << "|                      Denis Ridzal (dridzal@sandia.gov), or                  |\n"
          << "|                      Kara Peterson(kjpeter@sandia.gov), or                  |\n"
          << "|                      Kyungjoo Kim (kyukim@sandia.gov)                       |\n"
          << "|                                                                             |\n"
          << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
          << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
          << "|                                                                             |\n"
          << "===============================================================================\n";
  
        typedef CellTools<DeviceType> ct;
        typedef Kokkos::DynRankView<value_type,DeviceType> DynRankView;

        const value_type tol = tolerence<value_type>()*100.0;

        int errorFlag = 0;
      
        // Vertices of the parametrization domain for 1-subcells: standard 1-cube [-1,1]
        DynRankView ConstructWithLabel(cube_1, 2, 1);
        auto hCube_1 = Kokkos::create_mirror_view(cube_1);
        hCube_1(0,0) = -1.0;
        hCube_1(1,0) = 1.0;
        Kokkos::deep_copy(cube_1, hCube_1);
  
        // Vertices of the parametrization domain for triangular faces: the standard 2-simplex
        DynRankView ConstructWithLabel(simplex_2, 3, 2);
        auto hSimplex_2 = Kokkos::create_mirror_view(simplex_2);
        hSimplex_2(0, 0) = 0.0;   hSimplex_2(0, 1) = 0.0;
        hSimplex_2(1, 0) = 1.0;   hSimplex_2(1, 1) = 0.0;
        hSimplex_2(2, 0) = 0.0;   hSimplex_2(2, 1) = 1.0;
        Kokkos::deep_copy(simplex_2, hSimplex_2);
      
        // Vertices of the parametrization domain for quadrilateral faces: the standard 2-cube
        DynRankView ConstructWithLabel(cube_2, 4, 2);
        auto hCube_2 = Kokkos::create_mirror_view(cube_2);
        hCube_2(0, 0) =  -1.0;    hCube_2(0, 1) =  -1.0;
        hCube_2(1, 0) =   1.0;    hCube_2(1, 1) =  -1.0;
        hCube_2(2, 0) =   1.0;    hCube_2(2, 1) =   1.0;
        hCube_2(3, 0) =  -1.0;    hCube_2(3, 1) =   1.0;
        Kokkos::deep_copy(cube_2, hCube_2);
      
        try {
          // Pull all available topologies from Shards
          std::vector<shards::CellTopology> allTopologies;
          shards::getTopologies(allTopologies);

          const auto topoSize = allTopologies.size();

          *outStream
            << "\n"
            << "===============================================================================\n" 
            << "| Test 1: edge parametrizations:                                              |\n" 
            << "===============================================================================\n\n";
          {
            const auto subcDim = 1;
          
            // Loop over the cell topologies
            for (size_type topoOrd=0;topoOrd<topoSize;++topoOrd) 
              // Test only 2D and 3D topologies that have reference cells, e.g., exclude Line, Pentagon, etc.
              if ( allTopologies[topoOrd].getDimension() > 1 && 
                   ct::hasReferenceCell(allTopologies[topoOrd]) ) {
                *outStream << " Testing edge parametrization for " <<  allTopologies[topoOrd].getName() <<"\n";
                testSubcellParametrizations<value_type,DeviceType>( errorFlag,
                                                                         allTopologies[topoOrd],
                                                                         cube_1,
                                                                         cube_1,
                                                                         subcDim,
                                                                         tol,
                                                                         outStream );
              }
          }
    
          *outStream
            << "\n"
            << "===============================================================================\n" 
            << "| Test 2: face parametrizations:                                              |\n" 
            << "===============================================================================\n\n";
    
          {
            const auto subcDim = 2;
          
            // Loop over the cell topologies
            for (size_type topoOrd=0;topoOrd<topoSize;++topoOrd) 
              // Test only 3D topologies that have reference cells
              if ( allTopologies[topoOrd].getDimension() > 2 && 
                   ct::hasReferenceCell(allTopologies[topoOrd]) ) {
                *outStream << " Testing face parametrization for cell topology " <<  allTopologies[topoOrd].getName() <<"\n";
                testSubcellParametrizations<value_type,DeviceType>( errorFlag,
                                                                         allTopologies[topoOrd],
                                                                         simplex_2,
                                                                         cube_2,
                                                                         subcDim,
                                                                         tol,
                                                                         outStream );
              }
          }          
        } catch (std::logic_error &err) {
          //============================================================================================//
          // Wrap up test: check if the test broke down unexpectedly due to an exception                //
          //============================================================================================//
          *outStream << err.what() << "\n";
          errorFlag = -1000;
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
