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

#define INTREPID2_TEST_ERROR_EXPECTED( S )              \
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
    
    template<typename QuadBasisType, typename DeviceSpaceType>
    int OrientationToolsModifyBasis_QUAD_HCURL(const bool verbose) {

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
        << "|            Unit Test (OrientationTools, getModifiedHcurl/Hdiv_In_Basis)     |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      int errorFlag = 0;
      //const double tol = tolerence();
      constexpr ordinal_type maxOrder = Parameters::MaxOrder ;

      typedef OrientationTools<DeviceSpaceType> ots;
      try {
        constexpr ordinal_type order = 3;
        if(order <= maxOrder)
        {
          *outStream << "\n -- Testing Quadrilateral \n\n";


          QuadBasisType cellBasis(order);

          *outStream << "\n --- " << cellBasis.getName() << "\n";

          const auto cellTopo = cellBasis.getBaseCellTopology();
          const ordinal_type ndofBasis = cellBasis.getCardinality();
          
          // 
          // 9 12 13 16
          // 4  3 11 15
          // 5  2  8 14
          // 1  6  7 10
          ordinal_type refMesh[9][4] = { { 1, 6, 2, 5 },
                                         { 6, 7, 8, 2 },
                                         { 7,10,14, 8 },
                                         { 5, 2, 3, 4 },
                                         { 2, 8,11, 3 },
                                         { 8,14,15,11 },
                                         { 4, 3,12, 9 },
                                         { 3,11,13,12 },
                                         {11,15,16,13 } };

          const ordinal_type numCells = 9, numVerts = 4, numEdges = 4;

          // view to import refMesh from host          
          Kokkos::DynRankView<ordinal_type,Kokkos::LayoutRight,HostSpaceType> 
            elemNodesHost(&refMesh[0][0], numCells, numVerts);
          auto elemNodes = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), elemNodesHost);
          Kokkos::deep_copy(elemNodes, elemNodesHost);
          
          // compute orientations for cells (one time computation)
          Kokkos::DynRankView<Orientation,DeviceSpaceType> elemOrts("elemOrts", numCells);
          ots::getOrientation(elemOrts, elemNodes, cellTopo);

          auto elemOrtsHost = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), elemOrts);
          Kokkos::deep_copy(elemOrtsHost, elemOrts);
          
          // cell specific modified basis 
          Kokkos::DynRankView<double,DeviceSpaceType> outValues("outValues", numCells, ndofBasis);
          Kokkos::DynRankView<double,DeviceSpaceType> refValues("refValues", numCells, ndofBasis);
          
          auto refValuesHost = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), refValues);
          for (ordinal_type cell=0;cell<numCells;++cell)           
            for (ordinal_type bf=0;bf<ndofBasis;++bf) 
              refValuesHost(cell, bf) = (bf+1);
          Kokkos::deep_copy(refValues, refValuesHost);

          // modify refValues accounting for orientations
          ots::modifyBasisByOrientation(outValues,
                                        refValues,
                                        elemOrts,
                                        &cellBasis);

          auto outValuesHost = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), outValues);
          Kokkos::deep_copy(outValuesHost, outValues);

          for (ordinal_type cell=0;cell<numCells;++cell) {
            int flag = 0 ;
            std::stringstream s1, s2;
            s1 << std::fixed << std::setprecision(0);
            s2 << std::fixed << std::setprecision(0);

            ordinal_type orts[numEdges];
            elemOrtsHost(cell).getEdgeOrientation(orts, numEdges);
            
            // const double ortVal[4][2] = { {  1.0 , -1.0 },
            //                               {  1.0 , -1.0 },
            //                               { -1.0 ,  1.0 },
            //                               { -1.0 ,  1.0 } };
            
            s1 << " :: edge(0000) = " ;
            s2 << " :: edge(" << orts[0] << orts[1] << orts[2] << orts[3] << ") = ";
            for (ordinal_type edgeId=0;edgeId<numEdges;++edgeId) {
              //const ordinal_type ort = orts[edgeId];
              const ordinal_type ndof = cellBasis.getDofTag(cellBasis.getDofOrdinal(1, edgeId, 0))(3);
              for (ordinal_type i=0;i<ndof;++i) {
                //const double signVal = ortVal[edgeId][ort];
                //const ordinal_type ii = (signVal < 0 ? (ndof - i - 1) : i);

                //const ordinal_type refOrd = cellBasis.getDofOrdinal(1, edgeId, ii);
                const ordinal_type outOrd = cellBasis.getDofOrdinal(1, edgeId, i);
                s1 << std::setw(4) << refValuesHost(cell, outOrd);              
                s2 << std::setw(4) << outValuesHost(cell, outOrd);              


                //flag += (std::abs(signVal*outValuesHost(cell, outOrd) - refValuesHost(cell, refOrd)) > tol);
              }

              s1 << " // ";
              s2 << " // ";
            }

            *outStream << "\n cell = " << cell << "\n"
                       << " - refValues = " << s1.str() << "\n"
                       << " - outValues = " << s2.str() << "\n";
            if (flag) {
              *outStream << "                      ^^^^^^^^^^^^ FAILURE\n";
              errorFlag += flag;
            }
          }
          ots::clearCoeffMatrix();
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
    
