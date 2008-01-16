// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file
\brief Intrepid utilities.
\author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_UTILS_HPP
#define INTREPID_UTILS_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid {


enum TypeOfExactData{
  INTREPID_UTILS_FRACTION=0,
  INTREPID_UTILS_SCALAR
};


template<class Scalar>
int compareToAnalytic(const Teuchos::Array< Teuchos::Array<Scalar> > testMat,
                      std::ifstream & inputFile,
                      Scalar reltol,
                      int iprint,
                      TypeOfExactData analyticDataType = INTREPID_UTILS_FRACTION) {

  // This little trick lets us print to std::cout only if
  // iprint > 0.
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  std::string line;
  std::string chunk;
  Scalar testentry;
  Scalar abstol;
  Scalar absdiff;
  int i=0, j=0;
  int err = 0;

  while (! inputFile.eof() )
  {
    std::getline (inputFile,line);
    std::istringstream linestream(line);
    std::string chunk;
    j = 0;
    while( linestream >> chunk ) {
      int num1;
      int num2;
      std::string::size_type loc = chunk.find( "/", 0);
      if( loc != std::string::npos ) {
        chunk.replace( loc, 1, " ");
        std::istringstream chunkstream(chunk);
        chunkstream >> num1;
        chunkstream >> num2;
        testentry = (Scalar)(num1)/(Scalar)(num2);
        abstol = ( testentry == 0.0 ? reltol : std::fabs(reltol*testentry) );
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
        else if (analyticDataType == INTREPID_UTILS_SCALAR)
          chunkstream >> testentry;
        abstol = ( testentry == 0.0 ? reltol : std::fabs(reltol*testentry) );
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
void getAnalytic(Teuchos::Array< Teuchos::Array<Scalar> > & testMat,
                 std::ifstream & inputFile,
                 TypeOfExactData analyticDataType = INTREPID_UTILS_FRACTION) {

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  std::string line;
  std::string chunk;
  Scalar testentry;
  int i=0, j=0;

  while (! inputFile.eof() )
  {
    std::getline (inputFile,line);
    std::istringstream linestream(line);
    std::string chunk;
    j = 0;
    while( linestream >> chunk ) {
      int num1;
      int num2;
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
        else if (analyticDataType == INTREPID_UTILS_SCALAR)
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

} // end namespace Intrepid

#endif
