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
        else if (analyticDataType == INTREPID_UTILS_SCALAR)
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


/** \brief A helper function that fills enumeration to DoF tag and DoF tag to enumeration lookup
tables for the concrete basis classes in Intrepid.

\param tagToEnum  [out]  - Lookup table for the DoF's local enumeration (DoF Id) by its local DoF tag
\param enumToTag  [out]  - Lookup table for the DoF's local DoF tag by its local enumeration (DoF Id)
\param tags       [in]   - a set of basis-dependent local DoF tags in flat array format.
\param numBf      [in]   - number of basis functions in the basis set
\param tagSize    [in]   - number of fields in the local DoF tag
\param posScDim   [in]   - position in the tag, counting from 0, of the subcell dim 
\param posScId    [in]   - position in the tag, counting from 0, of the subcell id
\param posBfId    [in]   - position in the tag, counting from 0, of DoF Id relative to the subcell
*/
void setEnumTagData(Teuchos::Array<Teuchos::Array<Teuchos::Array<int> > > & tagToEnum,
                    Teuchos::Array<LocalDofTag> & enumToTag,
                    const int *tags,
                    const int numBf,
                    const int tagSize,
                    const int posScDim,
                    const int posScId,
                    const int posBfId) {
  
  
  // Build one-dimensional enumToTag array: its dimension equals the number of basis functions
  enumToTag.resize(numBf);
  for (int i = 0; i < numBf; i++) {
    for (int j = 0; j < tagSize; j++) {
      enumToTag[i].tag_[j] = tags[i*tagSize+j];
    }
  }
  
  // Build the three-dimensional tagToEnum array:   
  // The leading dimension of tagToEnum is the max value of the 1st column + 1 (max subcell dim in the tag)
  int maxScDim = 0; 
  for (int i = 0; i < numBf; i++) { 
    if (maxScDim < tags[i*tagSize + posScDim]) {
      maxScDim = tags[i*tagSize + posScDim];
    }
  }
  maxScDim += 1;
  
  // The 2nd dimension of tagToEnum is the max value of the 2nd column + 1 (max subcell id in the tag) 
  int maxScId = 0;
  for (int i = 0; i < numBf; i++) { 
    if (maxScId < tags[i*tagSize + posScId]) {
      maxScId = tags[i*tagSize + posScId];
    }
  }
  maxScId += 1;
  
  // The 3rd dimension of tagToEnum is the max value of the 3rd column +1 (max subcell BfId in the tag) 
  int maxBfId = 0;
  for (int i = 0; i < numBf; i++) { 
    if (maxBfId < tags[i*tagSize + posBfId]) {
      maxBfId = tags[i*tagSize + posBfId];
    }
  }
  maxBfId += 1;
  
  // Create 1 x maxBfId array (dimensioned by the 3rd dimension of tagToEnum) filled with -1
  Teuchos::Array<int> level1(maxBfId, -1);
  
  // Create maxScId x maxBfId array (dimensioned by the 2nd and 3rd dimensions of tagToEnum)
  Teuchos::Array<Teuchos::Array<int> > level2(maxScId, level1);
  
  // Create the maxScDim x maxScId x maxBfId tagToEnum array
  tagToEnum.assign(maxScDim, level2);
  
  // Overwrite elements of the array corresponding to tags with local DoF Id's, leave all other = -1
  for (int i = 0; i < numBf; i++) {
    tagToEnum[tags[i*tagSize]][tags[i*tagSize+1]][tags[i*tagSize+2]] = i;
  }
}

} // end namespace Intrepid

#endif
