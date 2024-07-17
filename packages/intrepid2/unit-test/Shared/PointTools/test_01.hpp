// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
      catch (std::logic_error &err) {                                    \
        *outStream << "Expected Error ----------------------------------------------------------------\n"; \
        *outStream << err.what() << '\n';                               \
        *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
      };                                                                \
    }
    
    template<typename ValueType, typename DeviceType>
    int PointTools_Test01(const bool verbose) {
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
      typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;
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
