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
    \brief  Test of the CellTools class.
    \author Created by P. Bochev, D. Ridzal, K. Peterson and Kyungjoo Kim
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_DefaultCubatureFactory.hpp"
#include "Intrepid2_CellTools.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {
  
  namespace Test {
#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      S ;                                                               \
    }                                                                   \
    catch (std::logic_error err) {                                      \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    };                                                                  


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
    template<typename CellToolType,
             typename subcParamVertAType,
             typename subcParamVertAType,
             typename outStreamPtrType,
    void testSubcellParametrizations( int&                       errorFlag,
                                      const shards::CellTopology parentCell,
                                      const subcParamVertAType   subcParamVert_A,
                                      const subcParamVertBType   subcParamVert_B,
                                      const int                  subcDim,
                                      const outStreamPtrType     outStreamPtr ) {
      
      // Get cell dimension and subcell count
      int cellDim      = parentCell.getDimension();
      int subcCount    = parentCell.getSubcellCount(subcDim);
      
      
      // Loop over subcells of the specified dimension
      for(int subcOrd = 0; subcOrd < subcCount; subcOrd++){
        int subcVertexCount = parentCell.getVertexCount(subcDim, subcOrd);
        
        
        // Storage for correct reference subcell vertices and for the images of the parametrization domain points
        DynRankView ConstructWithLabel refSubcellVertices(subcVertexCount, cellDim);
        DynRankView ConstructWithLabel mappedParamVertices(subcVertexCount, cellDim);
        
        
        // Retrieve correct reference subcell vertices
        CellTools<double>::getReferenceSubcellVertices(refSubcellVertices, subcDim, subcOrd, parentCell);
        
        
        // Map vertices of the parametrization domain to 1 or 2-subcell with ordinal subcOrd
        // For edges parametrization domain is always 1-cube passed as "subcParamVert_A"
        if(subcDim == 1) {
          CellTools<double>::mapToReferenceSubcell(mappedParamVertices,
                                                   subcParamVert_A,
                                                   subcDim,
                                                   subcOrd,
                                                   parentCell);
        }
        // For faces need to treat Triangle and Quadrilateral faces separately
        else if(subcDim == 2) {
          
          // domain "subcParamVert_A" is the standard 2-simplex  
          if(subcVertexCount == 3){
            CellTools<double>::mapToReferenceSubcell(mappedParamVertices,
                                                     subcParamVert_A,
                                                     subcDim,
                                                     subcOrd,
                                                     parentCell);
          }
          // Domain "subcParamVert_B" is the standard 2-cube
          else if(subcVertexCount == 4){
            CellTools<double>::mapToReferenceSubcell(mappedParamVertices,
                                                     subcParamVert_B,
                                                     subcDim,
                                                     subcOrd,
                                                     parentCell);
          }
        }
        
        // Compare the images of the parametrization domain vertices with the true vertices.
        for(int subcVertOrd = 0; subcVertOrd < subcVertexCount; subcVertOrd++){
          for(int dim = 0; dim <  cellDim; dim++){
            
            if(mappedParamVertices(subcVertOrd, dim) != refSubcellVertices(subcVertOrd, dim) ) {
              errorFlag++; 
              *outStream 
                << std::setw(70) << "^^^^----FAILURE!" << "\n"
                << " Cell Topology = " << parentCell.getName() << "\n"
                << " Parametrization of subcell " << subcOrd << " which is "
                << parentCell.getName(subcDim,subcOrd) << " failed for vertex " << subcVertOrd << ":\n"
                << " parametrization map fails to map correctly coordinate " << dim << " of that vertex\n\n";
              
            }//if
          }// for dim 
        }// for subcVertOrd      
      }// for subcOrd
    }
      
    template<typename ValueType, typename DeviceSpaceType>
    int CellTools_Test01(const bool verbose) {
      typedef ValueType value_type;

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

      *outStream << "DeviceSpace::  "; DeviceSpaceType::print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType::print_configuration(*outStream, false);
      
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
  
      typedef CellTools<DeviceSpaceType> ct;
      typedef Kokkos::DynRankView<value_type,DeviceSpaceType> DynRankView;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

      int errorFlag = 0;

      
      // Vertices of the parametrization domain for 1-subcells: standard 1-cube [-1,1]
      DynRankView ConstructWithLabel(cube_1, 2, 1);
      cube_1(0,0) = -1.0; 
      cube_1(1,0) = 1.0;
      
  
      // Vertices of the parametrization domain for triangular faces: the standard 2-simplex
      DynRankView ConstructWithLabel(simplex_2, 3, 2);
      simplex_2(0, 0) = 0.0;   simplex_2(0, 1) = 0.0;
      simplex_2(1, 0) = 1.0;   simplex_2(1, 1) = 0.0;
      simplex_2(2, 0) = 0.0;   simplex_2(2, 1) = 1.0;
      
      
      // Vertices of the parametrization domain for quadrilateral faces: the standard 2-cube
      DynRankView ConstructWithLabel(cube_2, 4, 2);
      cube_2(0, 0) =  -1.0;    cube_2(0, 1) =  -1.0;
      cube_2(1, 0) =   1.0;    cube_2(1, 1) =  -1.0;
      cube_2(2, 0) =   1.0;    cube_2(2, 1) =   1.0;
      cube_2(3, 0) =  -1.0;    cube_2(3, 1) =   1.0;
      
  
      
      try{
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
          for (auto topoOrd=0;topoOrd<topoSize;++topoOrd) {
            
            // Test only 2D and 3D topologies that have reference cells, e.g., exclude Line, Pentagon, etc.
            if ( allTopologies[topoOrd].getDimension() > 1 && 
                 ct::hasReferenceCell(allTopologies[topoOrd]) ) {
              *outStream << " Testing edge parametrization for " <<  allTopologies[topoOrd].getName() <<"\n";
              testSubcellParametrizations( errorFlag,
                                           allTopologies[topoOrd],
                                           cube_1,
                                           cube_1,
                                           subcDim,
                                           outStream );
            }
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
          for (auto topoOrd=0;topoOrd<topoSize;++topoOrd) {

            // Test only 3D topologies that have reference cells
            if ( allTopologies[topoOrd].getDimension() > 2 && 
                 CellTools::hasReferenceCell(allTopologies[topoOrd]) ) {
              *outStream << " Testing face parametrization for cell topology " <<  allTopologies[topoOrd].getName() <<"\n";
              testSubcellParametrizations( errorFlag,
                                           allTopologies[topoOrd],
                                           simplex_2,
                                           cube_2,
                                           subcDim,
                                           outStream );
            }
          }
        }          
      } catch (std::logic_error err) {
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
