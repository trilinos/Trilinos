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
\brief  Unit test of LexContainer class
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_LexContainer.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

using namespace Intrepid;

int main(int argc, char *argv[]) {
  
  // This little trick lets us print to cout  only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original cout .
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream  \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|                           Unit Test LexContainer                            |\n" \
    << "|                                                                             |\n" \
    << "|     1) Shaping LexContainer based on fieldType and operatorType             |\n" \
    << "|     2) Testing exception handling if intrepid-debug is enabled              |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";
  
  // Test initializations
  int errorFlag  = 0;
  
  // Upper ranges for numer of fields and number of evaluation points
  int maxNumPoints = 10;
  int maxNumFields = 10;
  
  // These tests checks if a LexContainer has the correct size for different field & operator sets
  LexContainer<double> myContainer;
  int correctSize;
  int spaceDim;
  
  *outStream  \
    << "===============================================================================\n"\
    << "| TEST 1: LexContainers for scalar fields                                     |\n"\
    << "===============================================================================\n\n";
  // Alternate fieldType between FORM_0 and FORM_3
  
  // - to shape LexContainer to store VALUE(FIELD_FORM_0) in 2D (dimension not actually used in this case):
  *outStream << " Shape LexContainer to store VALUE(FIELD_FORM_0) in 2D: \n";
  spaceDim  = 2;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_0,
                                 OPERATOR_VALUE,
                                 spaceDim);
      
      // Size should equal numPoints*numFields
      correctSize = numPoints*numFields;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (1)" << "\n";     
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store GRAD(FIELD_FORM_3) in 2D (dimension is used in this case!):
  *outStream << " Shape  LexContainer to store GRAD(FIELD_FORM_3) in 2D: \n";
  spaceDim  = 2;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_3,
                                 OPERATOR_GRAD,
                                 spaceDim);
      
      // Size should equal numPoints*numFields*spaceDim
      correctSize =  numPoints*numFields*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (2)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store D2(FIELD_FORM_0) in 3D (dimension is used in this case!):
  *outStream << " Shape  LexContainer  to store D2(FIELD_FORM_0) in 3D: \n";
  spaceDim = 3;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_0,
                                 OPERATOR_D2,
                                 spaceDim);
      
      // Size should equal numPoints*numFields*spaceDim*spaceDim
      correctSize = numPoints*numFields*spaceDim*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (3)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store D3(FIELD_FORM_3) in 3D (dimension is used in this case!):
  *outStream << " Shape  LexContainer  to store D3(FIELD_FORM_3) in 3D: \n";
  spaceDim = 3;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_3,
                                 OPERATOR_D3,
                                 spaceDim);
      
      // Size should equal numPoints*numFields*spaceDim*spaceDim*spaceDim
      correctSize = numPoints*numFields*spaceDim*spaceDim*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (4)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store D4(FIELD_FORM_0) in 3D (dimension is used in this case!):
  *outStream << " Shape  LexContainer  to store D4(FIELD_FORM_0) in 3D: \n";
  spaceDim = 3;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_0,
                                 OPERATOR_D4,
                                 spaceDim);
      
      // Size should equal numPoints*numFields*spaceDim*spaceDim*spaceDim*spaceDim
      correctSize = numPoints*numFields*spaceDim*spaceDim*spaceDim*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (5)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store D5(FIELD_FORM_3) in 3D (dimension is used in this case!):
  *outStream << " Shape  LexContainer  to store D5(FIELD_FORM_3) in 3D: \n";
  spaceDim = 3;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_3,
                                 OPERATOR_D5,
                                 spaceDim);
      
      // Size should equal numPoints*numFields*spaceDim*spaceDim*spaceDim*spaceDim*spaceDim
      correctSize = numPoints*numFields*spaceDim*spaceDim*spaceDim*spaceDim*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (6)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store CURL(FIELD_FORM_0) in 2D (dimension is used in this case!):
  *outStream << " Shape  LexContainer  to store CURL(FIELD_FORM_0) in 2D: \n";  
  spaceDim = 2;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_0,
                                 OPERATOR_CURL,
                                 spaceDim);      
      
      // Size should equal numPoints*numFields*spaceDim
      correctSize = numPoints*numFields*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (7)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  *outStream  << "\n" \
    << "===============================================================================\n"\
    << "| TEST 2: LexContainers for vector fields                                     |\n"\
    << "===============================================================================\n\n";
  
  // - to shape LexContainer to store VALUE(FIELD_FORM_1) in 2D:
  *outStream << " Shape  LexContainer  to store VALUE(FIELD_FORM_1) in 2D: \n";
  spaceDim = 2;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_1,
                                 OPERATOR_VALUE,
                                 spaceDim);
      
      // Size should equal numPoints*numFields*spaceDim
      correctSize = numPoints*numFields*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (8)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store CURL(FIELD_FORM_2) values in 2D:
  *outStream << " Shape  LexContainer  to store CURL(FIELD_FORM_2) in 2D: \n";
  spaceDim = 2;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_2,
                                 OPERATOR_CURL,
                                 spaceDim);
      
      // Size should equal numPoints*numFields
      correctSize = numPoints*numFields;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (9)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store CURL(FIELD_VECTOR) values in 3D:
  *outStream << " Shape  LexContainer  to store CURL(FIELD_VECTOR) in 3D: \n";
  spaceDim = 3;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_VECTOR,
                                 OPERATOR_CURL,
                                 spaceDim);
      
      // Size should equal numPoints*numFields*spaceDim
      correctSize = numPoints*numFields*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (10)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store D1(FIELD_FORM_1) values in 3D:
  *outStream << " Shape  LexContainer  to store D1(FIELD_FORM_1) in 3D: \n";
  spaceDim = 3;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_1,
                                 OPERATOR_D1,
                                 spaceDim);
      
      // Size should equal numPoints*numFields*spaceDim*spaceDim
      correctSize = numPoints*numFields*spaceDim*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (11)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store D2(FIELD_FORM_2) values in 3D:
  *outStream << " Shape  LexContainer  to store D2(FIELD_FORM_2) in 3D: \n";
  spaceDim = 3;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_2,
                                 OPERATOR_D2,
                                 spaceDim);
      
      // Size should equal numPoints*numFields*spaceDim*spaceDim*spaceDim
      correctSize = numPoints*numFields*spaceDim*spaceDim*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (12)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store D3(FIELD_VECTOR) values in 3D:
  *outStream << " Shape  LexContainer  to store D3(FIELD_VECTOR) in 3D: \n";
  spaceDim = 3;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_VECTOR,
                                 OPERATOR_D3,
                                 spaceDim);
      
      // Size should equal numPoints*numFields*spaceDim*spaceDim*spaceDim*spaceDim
      correctSize = numPoints*numFields*spaceDim*spaceDim*spaceDim*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (13)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store D4(FIELD_FORM_1) values in 3D:
  *outStream << " Shape  LexContainer  to store D4(FIELD_FORM_1) in 3D: \n";
  spaceDim = 3;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_1,
                                 OPERATOR_D4,
                                 spaceDim);
      
      // Size should equal numPoints*numFields*spaceDim*spaceDim*spaceDim*spaceDim*spaceDim
      correctSize = numPoints*numFields*spaceDim*spaceDim*spaceDim*spaceDim*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (14)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  // - to shape LexContainer to store D5(FIELD_FORM_2) values in 2D:
  *outStream << " Shape  LexContainer  to store D5(FIELD_FORM_2) in 2D: \n";
  spaceDim = 2;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_2,
                                 OPERATOR_D5,
                                 spaceDim);
      
      // Size should equal numPoints*numFields*spaceDim*spaceDim*spaceDim*spaceDim*spaceDim*spaceDim
      correctSize = numPoints*numFields*spaceDim*spaceDim*spaceDim*spaceDim*spaceDim*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (14)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store DIV(FIELD_FORM_1) values in 3D: (checking on fieldType = FORM_1)
  *outStream << " Shape  LexContainer  to store DIV(FIELD_FORM_1) in 3D: \n";
  spaceDim = 3;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_1,
                                 OPERATOR_DIV,
                                 spaceDim);
      
      // Size should equal numPoints*numFields
      correctSize = numPoints*numFields;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (15)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store DIV(FIELD_FORM_2) values in 2D: (checking on fieldType = FORM_2)
  *outStream << " Shape  LexContainer  to store DIV(FIELD_FORM_2) in 2D: \n";
  spaceDim = 2;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_2,
                                 OPERATOR_DIV,
                                 spaceDim);
      
      // Size should equal numPoints*numFields
      correctSize = numPoints*numFields;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (16)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store DIV(FIELD_VECTOR) values in 2D: (checking on fieldType = VECTOR)
  *outStream << " Shape  LexContainer  to store DIV(FIELD_VECTOR) in 2D: \n";
  spaceDim = 2;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_VECTOR,
                                 OPERATOR_DIV,
                                 spaceDim);
      
      // Size should equal numPoints*numFields
      correctSize = numPoints*numFields;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (17)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  *outStream  << "\n" \
    << "===============================================================================\n"\
    << "| TEST 3: LexContainers for tensor fields                                     |\n"\
    << "===============================================================================\n\n";
  
  // - to shape LexContainer to store VALUE(TENSOR) values in 2D:
  *outStream << " Shape  LexContainer  to store VALUE(TENSOR) values in 2D: \n";
  spaceDim = 3;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_TENSOR,
                                 OPERATOR_VALUE,
                                 spaceDim);
      
      // Size should equal numPoints*numFields*spaceDim*spaceDim
      correctSize = numPoints*numFields*spaceDim*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (18)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store DIV(TENSOR) values in 2D:
  *outStream << " Shape  LexContainer  to store DIV(TENSOR) values in 2D: \n";
  spaceDim = 2;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_TENSOR,
                                 OPERATOR_DIV,
                                 spaceDim);
      
      
      // Size should equal numPoints*numFields*spaceDim
      correctSize = numPoints*numFields*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (19)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store CURL(TENSOR) values in 3D:
  *outStream << " Shape  LexContainer  to store CURL(TENSOR) values in 3D: \n";
  spaceDim = 3;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_TENSOR,
                                 OPERATOR_CURL,
                                 spaceDim);
      
      // Size should equal numPoints*numFields*spaceDim*spaceDim
      correctSize = numPoints*numFields*spaceDim*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (20)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store D1(TENSOR) values in 3D:
  *outStream << " Shape  LexContainer  to store D1(TENSOR) values in 3D: \n";
  spaceDim = 3;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_TENSOR,
                                 OPERATOR_D1,
                                 spaceDim);
      
      
      // Size should equal numPoints*numFields*spaceDim*spaceDim*spaceDim
      correctSize = numPoints*numFields*spaceDim*spaceDim*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (21)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store D2(TENSOR) values in 3D:
  *outStream << " Shape  LexContainer  to store D2(TENSOR) values in 3D: \n";
  spaceDim = 3;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_TENSOR,
                                 OPERATOR_D2,
                                 spaceDim);
      
      // Size should equal numPoints*numFields*spaceDim*spaceDim*spaceDim*spaceDim
      correctSize = numPoints*numFields*spaceDim*spaceDim*spaceDim*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (22)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // - to shape LexContainer to store D3(TENSOR) values in 3D:
  *outStream << " Shape  LexContainer  to store D3(TENSOR) values in 3D: \n";
  spaceDim = 3;
  for(int numPoints = 0; numPoints < maxNumPoints; numPoints++) {
    for (int numFields = 0; numFields <maxNumFields; numFields++) {  
      myContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_TENSOR,
                                 OPERATOR_D3,
                                 spaceDim);
      
      // Size should equal numPoints*numFields*spaceDim*spaceDim*spaceDim*spaceDim*spaceDim
      correctSize = numPoints*numFields*spaceDim*spaceDim*spaceDim*spaceDim*spaceDim;
      if( myContainer.getSize() != correctSize ){
        errorFlag++;
        *outStream  << std::setw(70) << "^^^^----FAILURE! (23)" << "\n";      
        *outStream << " Container size = " << myContainer.getSize() << "\n";
        *outStream << "   Correct size = " << correctSize << "\n";
      }
    }
  }
  
  
  // These tests should only run if intrepid was configured with --enable-intrepid-debug option
  // Each test is designed to cause an exception. The total number of all caught exceptions should
  // be the same as the number of tests in this section.
#ifdef HAVE_INTREPID_DEBUG
  
  *outStream << "\n" \
    << "===============================================================================\n"\
    << "| TEST 4: Catching exceptions                                                 |\n"\
    << "===============================================================================\n\n";
  
  int numOfPoints = 1;
  int numOfFields = 1;
  
  int numTestException = 8;
  int beginThrowNumber = TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + numTestException;
  
  try{ // Outer try block contains all tests for exception
    
    try{ // catch exception (1)
      
      // Trying to shape container for DIV(FIELD_VECTOR) in 1D is not allowed:
      *outStream << "\n" \
      << "===============================================================================\n"\
      << " Trying to shape container for DIV(FIELD_VECTOR) in 1D is not allowed \n";
      spaceDim = 1;
      myContainer.shapeContainer(numOfPoints,
                                 numOfFields,
                                 FIELD_VECTOR,
                                 OPERATOR_DIV,
                                 spaceDim);
    }
    catch (std::logic_error err) {
      *outStream  << err.what() << "\n";
    };
    
    
    try{ // catch exception (2)
      
      // Trying to shape container for DIV(FIELD_FORM_3) in > 1D is not allowed:
      *outStream << "\n" \
      << "===============================================================================\n"\
      << " Trying to shape container for DIV(FIELD_FORM_3) in > 1D is not allowed \n";
      spaceDim = 3;
      myContainer.shapeContainer(numOfPoints,
                                 numOfFields,
                                 FIELD_FORM_3,
                                 OPERATOR_DIV,
                                 spaceDim);
    }
    catch (std::logic_error err) {
      *outStream  << err.what() << "\n";
    };
    
 
    try{  // catch exception (3)
      
      // Trying to shape container for CURL(FIELD_VECTOR) in  1D is not allowed:
      *outStream << "\n" \
      << "===============================================================================\n"\
      <<  " Trying to shape container for CURL(FIELD_VECTOR) in 1D is not allowed \n";
      spaceDim = 1;
      myContainer.shapeContainer(numOfPoints,
                                 numOfFields,
                                 FIELD_VECTOR,
                                 OPERATOR_CURL,
                                 spaceDim);
    }
    catch (std::logic_error err) {
      *outStream  << err.what() << "\n";
    };
    
    
    try{  // catch exception (4)
      
      // Trying to shape container for CURL(FIELD_FORM_3) in  3D is not allowed:
      *outStream << "\n" \
      << "===============================================================================\n"\
      <<  " Trying to shape container for CURL(FIELD_FORM_3) in 3D is not allowed \n";
      spaceDim = 3;
      myContainer.shapeContainer(numOfPoints,
                                 numOfFields,
                                 FIELD_FORM_3,
                                 OPERATOR_CURL,
                                 spaceDim);
    }
    catch (std::logic_error err) {
      *outStream  << err.what() << "\n";
    };
    
    
    try{ // catch exception (5)
      
      // Trying to shape container for CURL(FIELD_TENSOR) in  2D is not allowed:
      *outStream << "\n" \
      << "===============================================================================\n"\
      << " Trying to shape container for CURL(FIELD_TENSOR) in 2D is not allowed \n";
      spaceDim = 2;
      myContainer.shapeContainer(numOfPoints,
                                 numOfFields,
                                 FIELD_TENSOR,
                                 OPERATOR_CURL,
                                 spaceDim);
    }
    catch (std::logic_error err) {
      *outStream  << err.what() << "\n";
    }
    
    
    
    try{ // catch exception (6)
      *outStream << "\n" \
      << "===============================================================================\n"\
      << " Trying to shape container using negative number of points \n";
      numOfPoints = -1;
      numOfFields =  1;
      spaceDim    =  2;
      myContainer.shapeContainer(numOfPoints,
                                 numOfFields,
                                 FIELD_TENSOR,
                                 OPERATOR_CURL,
                                 spaceDim);
    }
    catch(std::logic_error err) {
     *outStream << err.what() << "\n";
    }
    
    
    try{ // catch exception (7)
      *outStream << "\n" \
      << "===============================================================================\n"\
      << " Trying to shape container using negative number of fields \n";
      numOfPoints =  1;
      numOfFields = -1;
      spaceDim    =  2;
      myContainer.shapeContainer(numOfPoints,
                                 numOfFields,
                                 FIELD_TENSOR,
                                 OPERATOR_CURL,
                                 spaceDim);
    }
    catch(std::logic_error err) {
      *outStream << err.what() << "\n"; 
    }

    
    try{ // catch exception (8)
      *outStream << "\n" \
      << "===============================================================================\n"\
      << " Trying to shape container using invalid dimension \n";
      numOfPoints =  1;
      numOfFields =  1;
      spaceDim    =  0;
      myContainer.shapeContainer(numOfPoints,
                                 numOfFields,
                                 FIELD_TENSOR,
                                 OPERATOR_CURL,
                                 spaceDim);
    }
    catch(std::logic_error err) {
      *outStream << err.what() << "\n"; 
    }
    
    
    // Check if number of caught exceptions matches the expected number
    if (TestForException_getThrowNumber() != endThrowNumber) {
      errorFlag++;
    }
  } // outer try block
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream  << err.what() << "\n";
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  }
#endif
  
  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";
  
  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  
  return errorFlag;
}
