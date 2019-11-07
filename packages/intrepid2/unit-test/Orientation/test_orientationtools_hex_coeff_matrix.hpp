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
    \brief  Unit test of experimental high order assembly
    \author Created by Kyungjoo Kim
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Orientation.hpp"
#include "Intrepid2_OrientationTools.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {

  namespace Test {

#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    }                                                                   \
    catch (std::exception &err) {                                        \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }
    
    template<typename DeviceSpaceType>
    int OrientationToolsHexCoeffMatrix(const bool verbose) {
      
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
      *outStream << "\n";

      
      *outStream
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|                 Unit Test (OrientationTools, initialize)                    |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      int errorFlag = 0;
      //const double tol = tolerence();
      constexpr ordinal_type maxOrder = Parameters::MaxOrder;

      typedef OrientationTools<DeviceSpaceType> ots;
      try {

        const ordinal_type testOrderBegin = 1, testOrderEnd = std::min(4, maxOrder);
        for (ordinal_type testOrder=testOrderBegin;testOrder<=testOrderEnd;++testOrder) {
          *outStream << "\n -- Testing order " << testOrder << "\n" 
                     << "===============================================================================\n"
                     << "\n";
          {
            // I do not know how I can provide reference permutation for hexa or any 3D elements
            const auto cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
            const ordinal_type order = testOrder;
          
            {
              *outStream << "\n -- Testing Hexahedral HGRAD \n\n";
              Basis_HGRAD_HEX_Cn_FEM<DeviceSpaceType> cellBasis(order);
              if (cellBasis.requireOrientation()) {
                const auto matData = ots::createCoeffMatrix(&cellBasis);
                
                auto matDataHost = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), matData);
                Kokkos::deep_copy(matDataHost, matData);
                
                // check face only
                const ordinal_type faceDim = 2, numEdge = 12, numFace = 6, numOrt = 8;
                for (auto faceId=0;faceId<numFace;++faceId) {
                  if (cellBasis.getDofCount(faceDim, faceId)) {
                    const auto ordFace = cellBasis.getDofOrdinal(faceDim, faceId, 0);
                    const auto ndofFace = cellBasis.getDofTag(ordFace)(3);
                    
                    const Kokkos::pair<ordinal_type,ordinal_type> range(0, ndofFace);            
                    for (auto faceOrt=0;faceOrt<numOrt;++faceOrt) {
                      *outStream << "\n faceId = " << faceId << " faceOrt = " << faceOrt << "\n";
                      const auto mat = Kokkos::subview(matDataHost, numEdge + faceId, faceOrt, range, range);
                      for (auto i=0;i<ndofFace;++i) {
                        for (auto j=0;j<ndofFace;++j)
                          *outStream << std::setw(5) << std::fixed << std::setprecision(1) << mat(i,j);
                        *outStream << "\n";
                      }
                    }
                  } else {
                    *outStream << "  " << cellBasis.getName() << "  face " << faceId <<" does not have DOFs\n";
                  }
                }
              } else {
                *outStream << "  " << cellBasis.getName() << " does not require orientations\n";                
              }
            }
            ots::clearCoeffMatrix();
          }

          {
            // I do not know how I can provide reference permutation for hexa or any 3D elements
            const auto cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
            const ordinal_type order = testOrder;
          
            {
              *outStream << "\n -- Testing Hexahedral HCURL \n\n";
              Basis_HCURL_HEX_In_FEM<DeviceSpaceType> cellBasis(order);
              if (cellBasis.requireOrientation()) {
                const auto matData = ots::createCoeffMatrix(&cellBasis);
                
                auto matDataHost = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), matData);
                Kokkos::deep_copy(matDataHost, matData);
                
                // check face only
                const ordinal_type faceDim = 2, numEdge = 12, numFace = 6, numOrt = 8;
                for (auto faceId=0;faceId<numFace;++faceId) {
                  if (cellBasis.getDofCount(faceDim, faceId)) {
                    const auto ordFace = cellBasis.getDofOrdinal(faceDim, faceId, 0);
                    const auto ndofFace = cellBasis.getDofTag(ordFace)(3);
                    
                    const Kokkos::pair<ordinal_type,ordinal_type> range(0, ndofFace);            
                    for (auto faceOrt=0;faceOrt<numOrt;++faceOrt) {
                      *outStream << "\n faceId = " << faceId << " faceOrt = " << faceOrt << "\n";
                      const auto mat = Kokkos::subview(matDataHost, numEdge + faceId, faceOrt, range, range);
                      for (auto i=0;i<ndofFace;++i) {
                        for (auto j=0;j<ndofFace;++j)
                          *outStream << std::setw(5) << std::fixed << std::setprecision(1) << mat(i,j);
                        *outStream << "\n";
                      }
                    } 
                  } else {
                    *outStream << "  " << cellBasis.getName() << "  face " << faceId <<" does not have DOFs\n";                    
                  }
                }
              } else {
                *outStream << "  " << cellBasis.getName() << " does not require orientations\n";
              }
            }
            ots::clearCoeffMatrix();
          }

          {
            // I do not know how I can provide reference permutation for hexa or any 3D elements
            const auto cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
            const ordinal_type order = testOrder;
          
            {
              *outStream << "\n -- Testing Hexahedral HDIV \n\n";
              Basis_HDIV_HEX_In_FEM<DeviceSpaceType> cellBasis(order);
              if (cellBasis.requireOrientation()) {
                const auto matData = ots::createCoeffMatrix(&cellBasis);

                auto matDataHost = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), matData);
                Kokkos::deep_copy(matDataHost, matData);
                
                // check face only
                const ordinal_type faceDim = 2, numFace = 6, numOrt = 8;
                for (auto faceId=0;faceId<numFace;++faceId) {
                  if (cellBasis.getDofCount(faceDim, faceId)) {
                    const auto ordFace = cellBasis.getDofOrdinal(faceDim, faceId, 0);
                    const auto ndofFace = cellBasis.getDofTag(ordFace)(3);
                    
                    const Kokkos::pair<ordinal_type,ordinal_type> range(0, ndofFace);            
                    for (auto faceOrt=0;faceOrt<numOrt;++faceOrt) {
                      *outStream << "\n faceId = " << faceId << " faceOrt = " << faceOrt << "\n";
                      const auto mat = Kokkos::subview(matDataHost, faceId, faceOrt, range, range);
                      for (auto i=0;i<ndofFace;++i) {
                        for (auto j=0;j<ndofFace;++j)
                          *outStream << std::setw(5) << std::fixed << std::setprecision(1) << mat(i,j);
                        *outStream << "\n";
                      }
                    } 
                  } else {
                    *outStream << "  " << cellBasis.getName() << "  face " << faceId <<" does not have DOFs\n";
                  }
                }
              } else {
                *outStream << "  " << cellBasis.getName() << " does not require orientations\n";
              }
            }
            ots::clearCoeffMatrix();
          }


        }

      } catch (std::exception &err) {
        std::cout << " Exeption\n";
        *outStream << err.what() << "\n\n";
        errorFlag = -1000;
      }

      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED = " << errorFlag << "\n";
      else
        std::cout << "End Result: TEST PASSED\n";
      
      // reset format state of std::cout
      std::cout.copyfmt(oldFormatState);
      
      return errorFlag;
    }
  }
}
    
