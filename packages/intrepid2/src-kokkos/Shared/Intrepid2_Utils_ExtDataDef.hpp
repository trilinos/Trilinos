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

/** \file   Intrepid_Utils.cpp
    \brief  Implementation file for Intrepid_Utils.hpp.
    \author Created by P. Bochev and D. Ridzal Kyungjoo Kim.
*/

#ifndef __INTREPID2_UTILS_EXTDATA_DEF_HPP__
#define __INTREPID2_UTILS_EXTDATA_DEF_HPP__

#include "Intrepid2_Utils_ExtData.hpp"

namespace Intrepid2 {

  /***************************************************************************************************
   *                                                                                                 *
   *               Utility functions for handling external data in tests                             *
   *                                                                                                 *
   ***************************************************************************************************/

  template<class Scalar>
  ordinal_type compareToAnalytic(const Teuchos::Array< Teuchos::Array<Scalar> > testMat,
                                 std::ifstream & inputFile,
                                 Scalar reltol,
                                 ordinal_type iprint,
                                 TypeOfExactData analyticDataType ) {

    // This little trick lets us prordinal_type to std::cout only if
    // iprordinal_type > 0.
    Teuchos::RCP<std::ostream> outStream;
    Teuchos::oblackholestream bhs; // outputs nothing
    if (iprordinal_type > 0)
      outStream = Teuchos::rcp(&std::cout, false);
    else
      outStream = Teuchos::rcp(&bhs, false);

    // Save the format state of the original std::cout.
    Teuchos::oblackholestream oldFormatState;
    oldFormatState.copyfmt(std::cout);

    std::string line;
    Scalar testentry;
    Scalar abstol;
    Scalar absdiff;
    ordinal_type i=0, j=0;
    ordinal_type err = 0;

    while (! inputFile.eof() )
      {
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
            testentry = (Scalar)(num1)/(Scalar)(num2);
            abstol = ( std::fabs(testentry) < reltol ? reltol : std::fabs(reltol*testentry) );
            absdiff = std::fabs(testentry - testMat[i][j]);
            if (absdiff > abstol) {
              err++;
              *outStream << "FAILURE --> ";
            }
            *outStream << "entry[" << i << "," << j << "]:" << "   "
                       << testMat[i][j] << "   " << num1 << "/" << num2 << "   "
                       << absdiff << "   " << "<?" << "   " << abstol << "\n";
          }
          else {
            std::istringstream chunkstream(chunk);
            if (analyticDataType == INTREPID_UTILS_FRACTION) {
              chunkstream >> num1;
              testentry = (Scalar)(num1);
            }
            else if (analyticDataType == INTREPID2_UTILS_SCALAR)
              chunkstream >> testentry;
            abstol = ( std::fabs(testentry) < reltol ?reltol : std::fabs(reltol*testentry) );
            absdiff = std::fabs(testentry - testMat[i][j]);
            if (absdiff > abstol) {
              err++;
              *outStream << "FAILURE --> ";
            }
            *outStream << "entry[" << i << "," << j << "]:" << "   "
                       << testMat[i][j] << "   " << testentry << "   "
                       << absdiff << "   " << "<?" << "   " << abstol << "\n";
          }
          j++;
        }
        i++;
      }

    // reset format state of std::cout
    std::cout.copyfmt(oldFormatState);

    return err;
  } // end compareToAnalytic

  template<class Scalar>
  ordinal_type compareToAnalytic(const Scalar * testMat,
                                 std::ifstream & inputFile,
                                 Scalar reltol,
                                 ordinal_type iprint,
                                 TypeOfExactData analyticDataType ) {

    // This little trick lets us prordinal_type to std::cout only if
    // iprordinal_type > 0.
    Teuchos::RCP<std::ostream> outStream;
    Teuchos::oblackholestream bhs; // outputs nothing
    if (iprordinal_type > 0)
      outStream = Teuchos::rcp(&std::cout, false);
    else
      outStream = Teuchos::rcp(&bhs, false);

    // Save the format state of the original std::cout.
    Teuchos::oblackholestream oldFormatState;
    oldFormatState.copyfmt(std::cout);

    std::string line;
    Scalar testentry;
    Scalar abstol;
    Scalar absdiff;
    ordinal_type i=0, j=0, offset=0;
    ordinal_type err = 0;

    while (! inputFile.eof() )
      {
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
            testentry = (Scalar)(num1)/(Scalar)(num2);
            abstol = ( std::fabs(testentry) < reltol ? reltol : std::fabs(reltol*testentry) );
            absdiff = std::fabs(testentry - testMat[i*offset+j]);
            if (absdiff > abstol) {
              err++;
              *outStream << "FAILURE --> ";
            }
            *outStream << "entry[" << i << "," << j << "]:" << "   "
                       << testMat[i*offset+j] << "   " << num1 << "/" << num2 << "   "
                       << absdiff << "   " << "<?" << "   " << abstol << "\n";
          }
          else {
            std::istringstream chunkstream(chunk);
            if (analyticDataType == INTREPID_UTILS_FRACTION) {
              chunkstream >> num1;
              testentry = (Scalar)(num1);
            }
            else if (analyticDataType == INTREPID2_UTILS_SCALAR)
              chunkstream >> testentry;
            abstol = ( std::fabs(testentry) < reltol ?reltol : std::fabs(reltol*testentry) );
            absdiff = std::fabs(testentry - testMat[i*offset+j]);
            if (absdiff > abstol) {
              err++;
              *outStream << "FAILURE --> ";
            }
            *outStream << "entry[" << i << "," << j << "]:" << "   "
                       << testMat[i*offset+j] << "   " << testentry << "   "
                       << absdiff << "   " << "<?" << "   " << abstol << "\n";
          }
          j++;
        }
        i++;
        offset = j;
      }

    // reset format state of std::cout
    std::cout.copyfmt(oldFormatState);

    return err;
  } // end compareToAnalytic

  template<class Scalar>
  void getAnalytic(Teuchos::Array< Teuchos::Array<Scalar> > & testMat,
                   std::ifstream & inputFile,
                   TypeOfExactData analyticDataType ) {

    // Save the format state of the original std::cout.
    Teuchos::oblackholestream oldFormatState;
    oldFormatState.copyfmt(std::cout);

    std::string line;
    Scalar testentry;
    ordinal_type i=0, j=0;

    while (! inputFile.eof() )
      {
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
            testentry = (Scalar)(num1)/(Scalar)(num2);
            testMat[i][j] = testentry;
          }
          else {
            std::istringstream chunkstream(chunk);
            if (analyticDataType == INTREPID_UTILS_FRACTION) {
              chunkstream >> num1;
              testentry = (Scalar)(num1);
            }
            else if (analyticDataType == INTREPID2_UTILS_SCALAR)
              chunkstream >> testentry;
            testMat[i][j] = testentry;
          }
          j++;
        }
        i++;
      }

    // reset format state of std::cout
    std::cout.copyfmt(oldFormatState);
  } // end getAnalytic

  template<class Scalar>
  void getAnalytic(Scalar * testMat,
                   std::ifstream & inputFile,
                   TypeOfExactData analyticDataType) {

    // Save the format state of the original std::cout.
    Teuchos::oblackholestream oldFormatState;
    oldFormatState.copyfmt(std::cout);

    std::string line;
    Scalar testentry;
    ordinal_type i=0, j=0, offset=0;

    while (! inputFile.eof() )
      {
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
            testentry = (Scalar)(num1)/(Scalar)(num2);
            testMat[i*offset+j] = testentry;
          }
          else {
            std::istringstream chunkstream(chunk);
            if (analyticDataType == INTREPID_UTILS_FRACTION) {
              chunkstream >> num1;
              testentry = (Scalar)(num1);
            }
            else if (analyticDataType == INTREPID2_UTILS_SCALAR)
              chunkstream >> testentry;
            testMat[i*offset+j] = testentry;
          }
          j++;
        }
        i++;
        offset = j;
      }

    // reset format state of std::cout
    std::cout.copyfmt(oldFormatState);
  } // end getAnalytic

} // end namespace Intrepid2

#endif
