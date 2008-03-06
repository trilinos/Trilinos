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
\brief  Illustrates use of the LexContainer class.
\author Created by P. Bochev and D. Ridzal
*/

#include "Intrepid_LexContainer.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {
  cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                  Example use of the LexContainer class                      |\n" \
  << "|                                                                             |\n" \
  << "|    1) Shaping LexContainers to store basis function evaluations             |\n" \
  << "|    3) using LexContainer in DEBUG mode                                      |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n\n";
  
  // Define variables to create and use LexContainers
  Teuchos::Array<int> indexRange;
  Teuchos::Array<int> multiIndex;
  
  cout \
    << "===============================================================================\n"\
    << "| EXAMPLE 1: LexContainers for scalar fields                               |\n"\
    << "===============================================================================\n\n";
  
  LexContainer<double> scalarContainer;

  // Suppose that we have numFields scalar functions and numPoints evaluation points
  int numPoints = 3;
  int numFields = 0;
  int spaceDim  = 2;
  
  // - to shape LexContainer to store function values in 2D (dimension not actually used in this case):
  scalarContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_0,
                                 OPERATOR_VALUE,
                                 spaceDim);
  std::cout 
    << "Stores " << numFields << " VALUE(SCALAR) evaluated at " 
    << numPoints << " points, in " << spaceDim << " dimensions \n";
  std::cout << scalarContainer;
  
  // - to shape LexContainer to store values of gradient of the scalar in 2D (dimension is used in this case!):
  scalarContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_0,
                                 OPERATOR_GRAD,
                                 spaceDim);
  std::cout 
    << "Stores " << numFields << " GRAD(SCALAR) evaluated at " 
    << numPoints << " points, in " << spaceDim << " dimensions \n";
  std::cout << scalarContainer;
  
  // - to shape LexContainer to store values of D2 of the scalar in 2D (dimension is used in this case!):
  scalarContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_0,
                                 OPERATOR_D2,
                                 spaceDim);
  std::cout 
    << "Stores " << numFields << " D2(SCALAR) evaluated at " 
    << numPoints << " points, in " << spaceDim << " dimensions \n";
  std::cout << scalarContainer;
  
  
  // - to shape LexContainer to store values of CURL of the scalar in 2D (dimension is used in this case!):
  scalarContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_0,
                                 OPERATOR_CURL,
                                 spaceDim);
  std::cout 
    << "Stores " << numFields << " CURL(SCALAR) evaluated at " 
    << numPoints << " points, in " << spaceDim << " dimensions \n";
  std::cout << scalarContainer;
  

  cout << "\n" \
  << "===============================================================================\n"\
  << "| EXAMPLE 2: LexContainers for vector fields                                  |\n"\
  << "===============================================================================\n\n";
  
  // - to shape LexContainer to store vector field values in 2D:
  scalarContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_1,
                                 OPERATOR_VALUE,
                                 spaceDim);
  std::cout 
    << "Stores " << numFields << " VALUE(VECTOR) evaluated at " 
    << numPoints << " points, in " << spaceDim << " dimensions \n";
  std::cout << scalarContainer;

  // - to shape LexContainer to store curl(vector) values in 2D:
  scalarContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_1,
                                 OPERATOR_CURL,
                                 spaceDim);
  std::cout 
    << "Stores " << numFields << " CURL(VECTOR) evaluated at " 
    << numPoints << " points, in " << spaceDim << " dimensions \n";
  std::cout << scalarContainer;
  
  // - to shape LexContainer to store curl(vector) values in 3D:
  spaceDim = 3;
  scalarContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_1,
                                 OPERATOR_CURL,
                                 spaceDim);
  std::cout 
    << "Stores " << numFields << " CURL(VECTOR) evaluated at " 
    << numPoints << " points, in " << spaceDim << " dimensions \n";
  std::cout << scalarContainer;
  
  // - to shape LexContainer to store D1(vector) values in 3D:
  spaceDim = 3;
  scalarContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_1,
                                 OPERATOR_D1,
                                 spaceDim);
  std::cout 
    << "Stores " << numFields << " D1(VECTOR) evaluated at " 
    << numPoints << " points, in " << spaceDim << " dimensions \n";
  std::cout << scalarContainer;
  
  // - to shape LexContainer to store DIV(vector) values in 3D:
  spaceDim = 3;
  scalarContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_FORM_1,
                                 OPERATOR_DIV,
                                 spaceDim);
  std::cout 
    << "Stores " << numFields << " DIV(VECTOR) evaluated at " 
    << numPoints << " points, in " << spaceDim << " dimensions \n";
  std::cout << scalarContainer;
  
  
  
  cout << "\n" \
  << "===============================================================================\n"\
  << "| EXAMPLE 3: LexContainers for tensor fields                                  |\n"\
  << "===============================================================================\n\n";
  
  // - to shape LexContainer to store VALUE(TENSOR) values in 2D:
  spaceDim = 3;
  scalarContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_TENSOR,
                                 OPERATOR_VALUE,
                                 spaceDim);
  std::cout 
    << "Stores " << numFields << " VALUE(TENSOR) evaluated at " 
    << numPoints << " points, in " << spaceDim << " dimensions \n";
  std::cout << scalarContainer;
  

  // - to shape LexContainer to store DIV(TENSOR) values in 2D:
  spaceDim = 2;
  scalarContainer.shapeContainer(numPoints,
                                 numFields,
                                 FIELD_TENSOR,
                                 OPERATOR_DIV,
                                 spaceDim);
  std::cout 
    << "Stores " << numFields << " DIV(TENSOR) evaluated at " 
    << numPoints << " points, in " << spaceDim << " dimensions \n";
  std::cout << scalarContainer;
  
  
  cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 4: Catching exceptions                                              |\n"\
    << "===============================================================================\n\n";
  
  try{
    
    // Trying to shape container for DIV(VECTOR) in 1D is not allowed:
    cout << "\n" \
    << "===============================================================================\n"\
    << " Trying to shape container for DIV(VECTOR) in 1D is not allowed \n";
    spaceDim = 1;
    scalarContainer.shapeContainer(numPoints,
                                   numFields,
                                   FIELD_VECTOR,
                                   OPERATOR_DIV,
                                   spaceDim);
  }
  catch (std::logic_error err) {
    std::cout << err.what() << "\n";
  };
  
  
  try{
    
    // Trying to shape container for DIV(SCALAR) in > 1D is not allowed:
    cout << "\n" \
    << "===============================================================================\n"\
    << " Trying to shape container for DIV(SCALAR) in > 1D is not allowed \n";
    spaceDim = 3;
    scalarContainer.shapeContainer(numPoints,
                                   numFields,
                                   FIELD_FORM_3,
                                   OPERATOR_DIV,
                                   spaceDim);
  }
  catch (std::logic_error err) {
    std::cout << err.what() << "\n";
  };
  
  
  try{
    
    // Trying to shape container for CURL(SCALAR) in  3D is not allowed:
    cout << "\n" \
    << "===============================================================================\n"\
    <<  " Trying to shape container for CURL(SCALAR) in 3D is not allowed \n";
    spaceDim = 3;
    scalarContainer.shapeContainer(numPoints,
                                   numFields,
                                   FIELD_FORM_3,
                                   OPERATOR_CURL,
                                   spaceDim);
  }
  catch (std::logic_error err) {
    std::cout << err.what() << "\n";
  };
 
  
  try{
    
    // Trying to shape container for CURL(TENSOR) in  2D is not allowed:
    cout << "\n" \
    << "===============================================================================\n"\
    << " Trying to shape container for CURL(TENSOR) in 2D is not allowed \n";
    spaceDim = 2;
    scalarContainer.shapeContainer(numPoints,
                                   numFields,
                                   FIELD_TENSOR,
                                   OPERATOR_CURL,
                                   spaceDim);
  }
  catch (std::logic_error err) {
    std::cout << err.what() << "\n";
  };
   return 0;
}
