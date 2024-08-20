// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_Utils_ExtDataDef.hpp
    \brief  Definition file for utility functions for handling external data in tests
    \author Created by P. Bochev and D. Ridzal and Kyungjoo Kim.
*/

#ifndef __INTREPID2_UTILS_EXTDATA_DEF_HPP__
#define __INTREPID2_UTILS_EXTDATA_DEF_HPP__

namespace Intrepid2 {

  /***************************************************************************************************
   *                                                                                                 *
   *               Utility functions for handling external data in tests                             *
   *                                                                                                 *
   ***************************************************************************************************/

  template<typename ValueType,
           class ...testMatProperties>
  ordinal_type compareToAnalytic( std::ifstream &inputFile,
                                  const Kokkos::DynRankView<ValueType,testMatProperties...> testMat,
                                  const ValueType reltol,
                                  const ordinal_type iprint,
                                  const TypeOfExactData analyticDataType ) {
    INTREPID2_TEST_FOR_EXCEPTION( testMat.rank() != 2, std::invalid_argument,
                                  ">>> ERROR (compareToAnalytic): testMat must have rank 2");

    Teuchos::RCP<std::ostream> outStream;
    Teuchos::oblackholestream  outNothing;
    if (iprint > 0)
      outStream = Teuchos::rcp(&std::cout, false);
    else
      outStream = Teuchos::rcp(&outNothing, false);

    // Save the format state of the original std::cout.
    Teuchos::oblackholestream oldFormatState;
    oldFormatState.copyfmt(std::cout);

    std::string line;
    ValueType testentry;
    ValueType abstol;
    ValueType absdiff;
    ordinal_type i = 0, j = 0;
    ordinal_type err = 0;

    while (! inputFile.eof() && i < static_cast<ordinal_type>(testMat.extent(0)) ) {
      std::getline(inputFile,line);
      std::istringstream linestream(line);
      std::string chunk;
      j = 0;
      while( linestream >> chunk ) {
        ordinal_type num1;
        ordinal_type num2;
        std::string::size_type loc = chunk.find( "/", 0);
        if( loc != std::string::npos ) {
          chunk.replace( loc, 1, " ");
          std::istringstream chunkstream(chunk);
          chunkstream >> num1;
          chunkstream >> num2;
          testentry = (ValueType)(num1)/(ValueType)(num2);
          abstol = ( std::fabs(testentry) < reltol ? reltol : std::fabs(reltol*testentry) );
          absdiff = std::fabs(testentry - testMat(i, j));
          if (absdiff > abstol) {
            ++err;
            *outStream << "FAILURE --> ";
          }
          *outStream << "entry[" << i << "," << j << "]:" << "   "
                     << testMat(i, j) << "   " << num1 << "/" << num2 << "   "
                     << absdiff << "   " << "<?" << "   " << abstol << "\n";
        }
        else {
          std::istringstream chunkstream(chunk);
          if (analyticDataType == INTREPID2_UTILS_FRACTION) {
            chunkstream >> num1;
            testentry = (ValueType)(num1);
          }
          else if (analyticDataType == INTREPID2_UTILS_SCALAR)
            chunkstream >> testentry;
          abstol = ( std::fabs(testentry) < reltol ?reltol : std::fabs(reltol*testentry) );
          absdiff = std::fabs(testentry - testMat(i, j));
          if (absdiff > abstol) {
            ++err;
            *outStream << "FAILURE --> ";
          }
          *outStream << "entry[" << i << "," << j << "]:" << "   "
                     << testMat(i, j) << "   " << testentry << "   "
                     << absdiff << "   " << "<?" << "   " << abstol << "\n";
        }
        ++j;
      }
      ++i;
    }
    
    // reset format state of std::cout
    std::cout.copyfmt(oldFormatState);
    
    return err;
  } // end compareToAnalytic
  
  template<typename ValueType,
           class ...testMatProperties>
  void getAnalytic( Kokkos::DynRankView<ValueType,testMatProperties...> testMat,
                    std::ifstream &inputFile,
                    const TypeOfExactData analyticDataType ) {
    INTREPID2_TEST_FOR_EXCEPTION( testMat.rank() != 2, std::invalid_argument,
                                  ">>> ERROR (getToAnalytic): testMat must have rank 2");

    Teuchos::oblackholestream oldFormatState;
    oldFormatState.copyfmt(std::cout);
    
    std::string line;
    ValueType testentry;
    ordinal_type i = 0, j = 0;

    while (! inputFile.eof() && i < static_cast<ordinal_type>(testMat.extent(0)) ) {
      std::getline (inputFile,line);
      std::istringstream linestream(line);
      std::string chunk;
      j = 0;
      while( linestream >> chunk ) {
        ordinal_type num1;
        ordinal_type num2;
        std::string::size_type loc = chunk.find( "/", 0);
        if( loc != std::string::npos ) {
          chunk.replace( loc, 1, " ");
          std::istringstream chunkstream(chunk);
          chunkstream >> num1;
          chunkstream >> num2;
          testentry = (ValueType)(num1)/(ValueType)(num2);
          testMat(i, j) = testentry;
        }
        else {
          std::istringstream chunkstream(chunk);
          if (analyticDataType == INTREPID2_UTILS_FRACTION) {
            chunkstream >> num1;
            testentry = (ValueType)(num1);
          }
          else if (analyticDataType == INTREPID2_UTILS_SCALAR)
            chunkstream >> testentry;
          testMat(i, j) = testentry;
        }
        ++j;
      }
      ++i;
    }
    
    // reset format state of std::cout
    std::cout.copyfmt(oldFormatState);
  } // end getAnalytic

} // end namespace Intrepid2

#endif
