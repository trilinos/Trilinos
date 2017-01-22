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

/** \file test_01.hpp
    \brief  Unit test for the PointTools class.
    \author Created by R. Kirby and Kyungjoo Kim.
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_PointTools.hpp"
#include "Teuchos_oblackholestream.hpp"

namespace Intrepid2 {
  
  namespace Test {

#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    {                                                                   \
      try {                                                             \
        S ;                                                             \
      }                                                                 \
      catch (std::logic_error err) {                                    \
        *outStream << "Expected Error ----------------------------------------------------------------\n"; \
        *outStream << err.what() << '\n';                               \
        *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
      };                                                                \
    }
    
    template<typename ValueType, typename DeviceSpaceType>
    int PointTools_Test01(const bool verbose) {
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
        << "|                       Unit Test (PointTools)                                |\n"
        << "|                                                                             |\n"
        << "|     1) Construction of equispaced and warped lattices on simplices          |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n"
        << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n"
        << "|                      Robert Kirby (robert.c.kirby@ttu.edu) or               |\n"
        << "|                      Kyungjoo Kim (kyukim@sandia.gov)                       |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";
      typedef PointTools pts;
      typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)
      
      int errorFlag = 0;     
      
      *outStream 
        << "\n"
        << "===============================================================================\n" 
        << "| TEST 1: size of lattices                                                    |\n" 
        << "===============================================================================\n";

      try {
        const shards::CellTopology line( shards::getCellTopologyData< shards::Line<2> >() );      
        const ordinal_type order[2] = { 4, 3 }, offset[2] = { 0, 1 }, val[2] = { 5, 2}; 
        for (auto i=0;i<2;++i) {
          const auto lsize = pts::getLatticeSize( line, order[i], offset[i] ); 
          if (lsize != val[i]) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << "  lattice size of line (order,offset,exact): " 
                       << order[i] << ", " << offset[i] << ", " << val[i] << ":: computed val = " 
                       << lsize << "\n";
          }
        }

      } catch (std::exception &err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };
      
      *outStream
        << "\n"
        << "===============================================================================\n" 
        << "| TEST 2: check for unsupported cell types                                     \n" 
        << "===============================================================================\n";
#ifdef HAVE_INTREPID2_DEBUG
      try {
        ordinal_type nthrow = 0, ncatch = 0;
        INTREPID2_TEST_ERROR_EXPECTED((pts::getLatticeSize(shards::getCellTopologyData<shards::Quadrilateral<4> >(), 3, 0)));
        INTREPID2_TEST_ERROR_EXPECTED((pts::getLatticeSize(shards::getCellTopologyData<shards::Hexahedron<8> >(), 3, 0 )));

        if (nthrow != ncatch) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << ncatch << ")\n";
        }
      } catch (std::exception &err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };
#endif

      *outStream                                
        << "\n"
        << "===============================================================================\n" 
        << "| TEST 2: malformed point arrays                                               \n" 
        << "===============================================================================\n";
#ifdef HAVE_INTREPID2_DEBUG
      try {
        const shards::CellTopology line( shards::getCellTopologyData< shards::Line<2> >() );      
        ordinal_type nthrow = 0, ncatch = 0;
        {
          DynRankView ConstructWithLabel(points, 3, 1);
          INTREPID2_TEST_ERROR_EXPECTED( PointTools::getLattice( points, line , 5 , 0 ) );
        }
        {
          DynRankView ConstructWithLabel(points, 6, 2);
          INTREPID2_TEST_ERROR_EXPECTED( PointTools::getLattice( points, line , 5 , 0 ) );
        }
        if (nthrow != ncatch) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << ncatch << ")\n";
        }
      } catch (std::exception &err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };       
#endif
      
      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
      else
        std::cout << "End Result: TEST PASSED\n";
      
      return errorFlag;
    }
  }
}
