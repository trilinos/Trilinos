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

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {

  namespace Test {

#define INTREPID2_TEST_ERROR_EXPECTED( S )              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    }                                                                   \
    catch (std::exception err) {                                        \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }
    
    template<typename DeviceSpaceType>
    int Orientation_Test01(const bool verbose) {

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
        << "|                 Unit Test (Orientation)                                     |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      int errorFlag = 0;

      try {
        ordinal_type nthrow = 0, ncatch = 0;
        {
          Orientation ort;
          INTREPID2_TEST_ERROR_EXPECTED( if (ort.isAlignedToReference()) \
                                           throw std::logic_error("Default Orientation is not zero"); );
        }

        if (nthrow != ncatch) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << nthrow << ")\n";
        }

        {
          *outStream << "\n -- Testing Triangle \n\n";

          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >() );
          const ordinal_type elemNodes[6][3] = { { 1, 2, 3 },
                                                 { 1, 3, 2 },
                                                 { 2, 1, 3 },
                                                 { 2, 3, 1 },
                                                 { 3, 1, 2 },
                                                 { 3, 2, 1 } };

          const ordinal_type refEdgeOrts[6][3] = { { 0, 0, 1 },
                                                   { 0, 1, 1 },
                                                   { 1, 0, 1 },
                                                   { 0, 1, 0 },
                                                   { 1, 0, 0 },
                                                   { 1, 1, 0 } };
          
          for (auto i=0;i<6;++i) {
            // find orientation
            const auto nodes = Kokkos::View<const ordinal_type[3],HostSpaceType>(elemNodes[i]);
            const auto ort = Orientation::getOrientation(cellTopo, nodes);
            
            // decode orientation
            ordinal_type edgeOrt[3] = {};
            for (auto edgeId=0;edgeId<3;++edgeId) 
              ort.getEdgeOrientation(edgeOrt, 3);
            
            *outStream << " elemNodes = " 
                       << elemNodes[i][0] << " " 
                       << elemNodes[i][1] << " " 
                       << elemNodes[i][2] << " :: "
                       << " computed edgeOrts = " 
                       << edgeOrt[0] << " "
                       << edgeOrt[1] << " "
                       << edgeOrt[2] << " :: "
                       << " reference edgeOrts = " 
                       << refEdgeOrts[i][0] << " "
                       << refEdgeOrts[i][1] << " "
                       << refEdgeOrts[i][2] << " ::\n";
            
            if (edgeOrt[0] != refEdgeOrts[i][0] || 
                edgeOrt[1] != refEdgeOrts[i][1] || 
                edgeOrt[2] != refEdgeOrts[i][2])  {
              *outStream << "        ^^^^^^^^^^^^^^^^ FAILURE\n";
              ++errorFlag;
            }
          }
        }

        {
          *outStream << "\n -- Testing Quadrilateral \n\n";
          
          const auto cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
          const ordinal_type elemNodes[24][4] = { { 1 ,  2 ,  3 ,  4 },
                                                  { 2 ,  1 ,  3 ,  4 },
                                                  { 1 ,  3 ,  2 ,  4 },
                                                  { 2 ,  3 ,  1 ,  4 },
                                                  { 3 ,  1 ,  2 ,  4 },
                                                  { 3 ,  2 ,  1 ,  4 },
                                                  { 1 ,  2 ,  4 ,  3 },
                                                  { 2 ,  1 ,  4 ,  3 },
                                                  { 1 ,  3 ,  4 ,  2 },
                                                  { 2 ,  3 ,  4 ,  1 },
                                                  { 3 ,  1 ,  4 ,  2 },
                                                  { 3 ,  2 ,  4 ,  1 },
                                                  { 1 ,  4 ,  2 ,  3 },
                                                  { 2 ,  4 ,  1 ,  3 },
                                                  { 1 ,  4 ,  3 ,  2 },
                                                  { 2 ,  4 ,  3 ,  1 },
                                                  { 3 ,  4 ,  1 ,  2 },
                                                  { 3 ,  4 ,  2 ,  1 },
                                                  { 4 ,  1 ,  2 ,  3 },
                                                  { 4 ,  2 ,  1 ,  3 },
                                                  { 4 ,  1 ,  3 ,  2 },
                                                  { 4 ,  2 ,  3 ,  1 },
                                                  { 4 ,  3 ,  1 ,  2 },
                                                  { 4 ,  3 ,  2 ,  1 } };
          
          const ordinal_type refEdgeOrts[24][4] = { { 0, 0, 0, 1 },
                                                    { 1, 0, 0, 1 },
                                                    { 0, 1, 0, 1 },
                                                    { 0, 1, 0, 1 },
                                                    { 1, 0, 0, 1 },
                                                    { 1, 1, 0, 1 },
                                                    { 0, 0, 1, 1 },
                                                    { 1, 0, 1, 1 },
                                                    { 0, 0, 1, 1 },
                                                    { 0, 0, 1, 0 },
                                                    { 1, 0, 1, 0 },
                                                    { 1, 0, 1, 0 },
                                                    { 0, 1, 0, 1 },
                                                    { 0, 1, 0, 1 },
                                                    { 0, 1, 1, 1 },
                                                    { 0, 1, 1, 0 },
                                                    { 0, 1, 0, 0 },
                                                    { 0, 1, 1, 0 },
                                                    { 1, 0, 0, 0 },
                                                    { 1, 1, 0, 0 },
                                                    { 1, 0, 1, 0 },
                                                    { 1, 0, 1, 0 },
                                                    { 1, 1, 0, 0 },
                                                    { 1, 1, 1, 0 } };

          for (auto i=0;i<24;++i) {
            // find orientation
            const auto nodes = Kokkos::View<const ordinal_type[4],HostSpaceType>(elemNodes[i]);
            const auto ort = Orientation::getOrientation(cellTopo, nodes);
            
            // decode orientation
            ordinal_type edgeOrt[4] = {};
            for (auto edgeId=0;edgeId<4;++edgeId) 
              ort.getEdgeOrientation(edgeOrt, 4);
            
            *outStream << " elemNodes = " 
                       << elemNodes[i][0] << " " 
                       << elemNodes[i][1] << " " 
                       << elemNodes[i][2] << " " 
                       << elemNodes[i][3] << " :: "
                       << " computed edgeOrts = " 
                       << edgeOrt[0] << " "
                       << edgeOrt[1] << " "
                       << edgeOrt[2] << " "
                       << edgeOrt[3] << " :: "
                       << " reference edgeOrts = " 
                       << refEdgeOrts[i][0] << " "
                       << refEdgeOrts[i][1] << " "
                       << refEdgeOrts[i][2] << " "
                       << refEdgeOrts[i][3] << " ::\n";
            
            if (edgeOrt[0] != refEdgeOrts[i][0] || 
                edgeOrt[1] != refEdgeOrts[i][1] || 
                edgeOrt[2] != refEdgeOrts[i][2] ||
                edgeOrt[3] != refEdgeOrts[i][3])  {
              *outStream << "        ^^^^^^^^^^^^^^^^ FAILURE\n";
              ++errorFlag;
            }
          }
        }

      } catch (std::exception err) {
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
      
      Kokkos::finalize();
      
      return errorFlag;
    }
  }
}
    
