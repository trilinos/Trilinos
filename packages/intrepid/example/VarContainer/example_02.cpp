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

#include "Intrepid_VarContainer.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {
  cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                  Example use of the VarContainer class                      |\n" \
  << "|                                                                             |\n" \
  << "|    1) Creating and resizing VarContainer objects                            |\n" \
  << "|    2) using VarContainer in DEBUG mode                                      |\n" \
  << "|       requires intrepid to be configured with --enable-intrepid-debug       |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n\n";
  
  // Define variables to make and use VarContainers
  VarContainer<double> myContainer;
  Teuchos::Array<int>  indexRange;
  Teuchos::Array<int>  multiIndex;
  Teuchos::Array<int>  partialMult;
  int numPoints;
  int numFields;
  int spaceDim;
  
  cout \
    << "===============================================================================\n"\
    << "| EXAMPLE 1: VarContainers for scalar fields                                  |\n"\
    << "===============================================================================\n\n";
  
  // - reset VarContainer for function VALUE in 2D:
  numPoints = 3;
  numFields = 0;
  spaceDim  = 3;
  myContainer.reset(numPoints,
                    numFields,
                    FIELD_FORM_0,
                    OPERATOR_VALUE,
                    spaceDim);
  std::cout << myContainer;
  
  // - reset VarContainer for values of GRAD of the scalar in 2D:
  numPoints = 3;
  numFields = 2;
  spaceDim  = 2;
  myContainer.reset(numPoints,
                    numFields,
                    FIELD_FORM_3,
                    OPERATOR_GRAD,
                    spaceDim);
  std::cout << myContainer;
  
  // - reset VarContainer for values of D4 of the scalar in 3D:
  numPoints = 1;
  numFields = 1;
  spaceDim  = 3;
  myContainer.reset(numPoints,
                    numFields,
                    FIELD_FORM_0,
                    OPERATOR_D4,
                    spaceDim);
  std::cout << myContainer;
  
  // - reset VarContainer for values of CURL of the scalar in 2D:
  numPoints = 3;
  numFields = 2;
  spaceDim  = 2;
  myContainer.reset(numPoints,
                    numFields,
                    FIELD_FORM_3,
                    OPERATOR_CURL,
                    spaceDim);
  std::cout << myContainer;
  
  cout \
    << "===============================================================================\n"\
    << "| EXAMPLE 2: VarContainers for vector fields                                  |\n"\
    << "===============================================================================\n\n";
  
  // - reset VarContainer for VECTOR values in 2D:
  numPoints = 4;
  numFields = 1;
  spaceDim  = 2;
  myContainer.reset(numPoints,
                    numFields,
                    FIELD_FORM_1,
                    OPERATOR_VALUE,
                    spaceDim);
  std::cout << myContainer;
  
  // - reset VarContainer for CURL OF VECTOR values in 3D:
  numPoints = 1;
  numFields = 2;
  spaceDim  = 3;
  myContainer.reset(numPoints,
                    numFields,
                    FIELD_FORM_2,
                    OPERATOR_CURL,
                    spaceDim);
  std::cout << myContainer;
  
  // - reset VarContainer for CURL OF VECTOR values in 2D:
  numPoints = 1;
  numFields = 2;
  spaceDim  = 2;
  myContainer.reset(numPoints,
                    numFields,
                    FIELD_FORM_2,
                    OPERATOR_CURL,
                    spaceDim);
  std::cout << myContainer;
  
  // - reset VarContainer for DIV of VECTOR values in 3D:
  numPoints = 5;
  numFields = 3;
  spaceDim  = 3;
  myContainer.reset(numPoints,
                    numFields,
                    FIELD_VECTOR,
                    OPERATOR_DIV,
                    spaceDim);
  std::cout << myContainer;
  
  // - reset VarContainer for D2 of VECTOR values in 2D:
  numPoints = 2;
  numFields = 1;
  spaceDim  = 2;
  myContainer.reset(numPoints,
                    numFields,
                    FIELD_VECTOR,
                    OPERATOR_D2,
                    spaceDim);
  std::cout << myContainer;

  
  cout \
    << "===============================================================================\n"\
    << "| EXAMPLE 3: VarContainers for tensor fields                                  |\n"\
    << "===============================================================================\n\n";
  
  // - reset VarContainer for TENSOR VALUE in 2D:
  numPoints = 2;
  numFields = 3;
  spaceDim  = 2;
  myContainer.reset(numPoints,
                    numFields,
                    FIELD_TENSOR,
                    OPERATOR_VALUE,
                    spaceDim);
  std::cout << myContainer;
  
  // - reset VarContainer for DIV of TENSOR values in 3D:
  numPoints = 1;
  numFields = 2;
  spaceDim  = 3;
  myContainer.reset(numPoints,
                    numFields,
                    FIELD_TENSOR,
                    OPERATOR_DIV,
                    spaceDim);
  std::cout << myContainer;
  
  // - reset VarContainer for D4 of TENSOR values in 3D:
  numPoints = 1;
  numFields = 2;
  spaceDim  = 3;
  myContainer.reset(numPoints,
                    numFields,
                    FIELD_TENSOR,
                    OPERATOR_D4,
                    spaceDim);
  std::cout << myContainer;
  
  
  cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 4: Catching exceptions                                              |\n"\
    << "===============================================================================\n\n";
  try{
    
    // Trying to reset container for DIV(VECTOR) in 1D is not allowed:
    cout << "\n" \
    << "===============================================================================\n"\
    << " Trying to reset container for DIV(VECTOR) in 1D is not allowed \n";
    spaceDim = 1;
    myContainer.reset(numPoints,
                      numFields,
                      FIELD_VECTOR,
                      OPERATOR_DIV,
                      spaceDim);
  }
  catch (std::logic_error err) {
    std::cout << err.what() << "\n";
  };
  
  
  try{
    
    // Trying to reset container for DIV(SCALAR) in > 1D is not allowed:
    cout << "\n" \
    << "===============================================================================\n"\
    << " Trying to reset container for DIV(SCALAR) in > 1D is not allowed \n";
    spaceDim = 3;
    myContainer.reset(numPoints,
                      numFields,
                      FIELD_FORM_3,
                      OPERATOR_DIV,
                      spaceDim);
  }
  catch (std::logic_error err) {
    std::cout << err.what() << "\n";
  };
  
  
  try{
    
    // Trying to reset container for CURL(SCALAR) in  3D is not allowed:
    cout << "\n" \
    << "===============================================================================\n"\
    <<  " Trying to reset container for CURL(SCALAR) in 3D is not allowed \n";
    spaceDim = 3;
    myContainer.reset(numPoints,
                      numFields,
                      FIELD_FORM_3,
                      OPERATOR_CURL,
                      spaceDim);
  }
  catch (std::logic_error err) {
    std::cout << err.what() << "\n";
  };
  return 0;
}


