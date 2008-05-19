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

/** \file test_01.cpp
\brief  Unit tests for the Intrepid::Basis_F0_QUAD_I1_FEM_DEFAULT class.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_DefaultBasisFactory.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {
  
  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|                 Unit Test (Basis_F0_QUAD_I1_FEM_DEFAULT)                     |\n" \
    << "|                                                                             |\n" \
    << "|     1) Basis creation, conversion of Dof tags into enumeration and back     |\n" \
    << "|     2) Basis values for VALUE, GRAD, CURL, and Dk operators                 |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n"\
    << "| TEST 1: Basis creation, exception testing                                   |\n"\
    << "===============================================================================\n";
  
  int errorFlag = 0;
  int beginThrowNumber = TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + 1;
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  endThrowNumber += 5;
#endif
  
  // Reference element points are the 4 vertices of the QUAD cell and the center point (0,0):
  Teuchos::Array< Point<double> > quadNodes;
  Point<double> tempPt(2, FRAME_REFERENCE);
  quadNodes.assign(5, tempPt);  
  quadNodes[0] = Point<double>( -1.0, -1.0, FRAME_REFERENCE);
  quadNodes[1] = Point<double>(  1.0, -1.0, FRAME_REFERENCE);
  quadNodes[2] = Point<double>(  1.0,  1.0, FRAME_REFERENCE);
  quadNodes[3] = Point<double>( -1.0,  1.0, FRAME_REFERENCE);
  quadNodes[4] = Point<double>(  0.0,  0.0, FRAME_REFERENCE);
  
  try{
    FieldContainer<double> vals;
    DefaultBasisFactory<double> BFactory;
    Teuchos::RCP<Basis<double> > quadBasis = BFactory.create(FIELD_FORM_0, 
                                                             CELL_QUAD, 
                                                             RECONSTRUCTION_SPACE_INCOMPLETE, 
                                                             1, 
                                                             BASIS_FEM_DEFAULT, 
                                                             COORDINATES_CARTESIAN);
    // Exception 1: DIV cannot be applied to scalar functions
    try {
      quadBasis->getValues(vals, quadNodes, OPERATOR_DIV);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    
    
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    // Exceptions 2-6: all tags/bd Ids below are wrong and should cause getLocalDofEnumeration() and 
    // getLocalDofTag() to access invalid array elements thereby causing Teuchos bounds check exception
    try {      
      LocalDofTag myTag = {{3,0,0,1}};
      quadBasis -> getLocalDofEnumeration(myTag);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    
    try {
      LocalDofTag myTag = {{1,1,1,0}};
      quadBasis -> getLocalDofEnumeration(myTag);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    
    try {
      LocalDofTag myTag = {{0,4,0,0}};
      quadBasis -> getLocalDofEnumeration(myTag);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    
    try {
      int bfId = 5;
      quadBasis -> getLocalDofTag(bfId);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    
    try {
      int bfId = -1;
      quadBasis -> getLocalDofTag(bfId);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif
  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };
  
  // Check if number of thrown excpetions matches the one we expect (1 + 5)
  if (TestForException_getThrowNumber() != endThrowNumber) {
    errorFlag++;
    *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
  }
  
  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 2: correctness of tag to enum and enum to tag lookups                  |\n"\
    << "===============================================================================\n";
  
  try{
    DefaultBasisFactory<double> BFactory;
    Teuchos::RCP<Basis<double> > quadBasis = BFactory.create(FIELD_FORM_0, 
                                                             CELL_QUAD, 
                                                             RECONSTRUCTION_SPACE_INCOMPLETE, 
                                                             1, 
                                                             BASIS_FEM_DEFAULT, 
                                                             COORDINATES_CARTESIAN);
    Teuchos::Array<LocalDofTag> allTags;
    allTags = quadBasis -> getAllLocalDofTags();
    
    // Loop over all tags, lookup the associated dof enumeration and then lookup the tag again
    for (unsigned i = 0; i < allTags.size(); i++) {
      int        bfId  = quadBasis -> getLocalDofEnumeration(allTags[i]);
      LocalDofTag myTag = quadBasis -> getLocalDofTag(bfId);
       if( !( (myTag.tag_[0] == allTags[i].tag_[0]) &&
             (myTag.tag_[1] == allTags[i].tag_[1]) &&
             (myTag.tag_[2] == allTags[i].tag_[2]) &&
             (myTag.tag_[3] == allTags[i].tag_[3]) ) ) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " getLocalDofEnumeration( {" 
          << allTags[i].tag_[0] << ", " 
          << allTags[i].tag_[1] << ", " 
          << allTags[i].tag_[2] << ", " 
          << allTags[i].tag_[3] << "}) = " << bfId <<" but \n";   
        *outStream << " getLocalDofTag(" << bfId << ") = { "
          << myTag.tag_[0] << ", " 
          << myTag.tag_[1] << ", " 
          << myTag.tag_[2] << ", " 
          << myTag.tag_[3] << "}\n";        
      }
    }
    
    // Now do the same but loop over basis functions
    for( int bfId = 0; bfId < quadBasis -> getNumLocalDof(); bfId++) {
      LocalDofTag myTag  = quadBasis -> getLocalDofTag(bfId);
      int myBfId = quadBasis -> getLocalDofEnumeration(myTag);
      if( bfId != myBfId) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " getLocalDofTag(" << bfId << ") = { "
          << myTag.tag_[0] << ", " 
          << myTag.tag_[1] << ", " 
          << myTag.tag_[2] << ", " 
          << myTag.tag_[3] << "} but getLocalDofEnumeration({" 
          << myTag.tag_[0] << ", " 
          << myTag.tag_[1] << ", " 
          << myTag.tag_[2] << ", " 
          << myTag.tag_[3] << "} ) = " << myBfId << "\n";
      }
    }
  }
  catch (std::logic_error err){
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };
  
  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 3: correctness of basis function values                                |\n"\
    << "===============================================================================\n";
  
  outStream -> precision(20);
  
  // VALUE: Each row gives the 4 correct basis set values at an evaluation point
  double basisValues[] = {
    1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 1.0,
    0.25,0.25,0.25,0.25
  };
  
  // GRAD and D1: each row gives the 8 correct values of the gradients of the 4 basis functions
  double basisGrads[] = {
    -0.5, -0.5,    0.5,  0.0,    0.0,  0.0,    0.0,  0.5,
    -0.5,  0.0,    0.5, -0.5,    0.0,  0.5,    0.0,  0.0,
     0.0,  0.0,    0.0, -0.5,    0.5,  0.5,   -0.5,  0.0,
     0.0, -0.5,    0.0,  0.0,    0.5,  0.0,   -0.5,  0.5,
    -0.25,-0.25,   0.25,-0.25,   0.25, 0.25,  -0.25, 0.25 
  };
  
  // CURL: each row gives the 8 correct values of the curls of the 4 basis functions
  double basisCurls[] = {
    -0.5,  0.5,    0.0, -0.5,    0.0,  0.0,    0.5,  0.0,
     0.0,  0.5,   -0.5, -0.5,    0.5,  0.0,    0.0,  0.0,
     0.0,  0.0,   -0.5,  0.0,    0.5, -0.5,    0.0,  0.5,
    -0.5,  0.0,    0.0,  0.0,    0.0, -0.5,    0.5,  0.5,    
    -0.25, 0.25,  -0.25,-0.25,   0.25,-0.25,   0.25, 0.25 
  };
  
  //D2: each row gives the 12 correct values of all 2nd derivatives of the 4 basis functions
  double basisD2[] = {
    0.0, 0.25, 0.0,   0.0,-0.25, 0.0,   0.0, 0.25, 0.0,   0.0,-0.25, 0.0,
    0.0, 0.25, 0.0,   0.0,-0.25, 0.0,   0.0, 0.25, 0.0,   0.0,-0.25, 0.0,
    0.0, 0.25, 0.0,   0.0,-0.25, 0.0,   0.0, 0.25, 0.0,   0.0,-0.25, 0.0,
    0.0, 0.25, 0.0,   0.0,-0.25, 0.0,   0.0, 0.25, 0.0,   0.0,-0.25, 0.0,
    0.0, 0.25, 0.0,   0.0,-0.25, 0.0,   0.0, 0.25, 0.0,   0.0,-0.25, 0.0,
  };
  
  try{
    FieldContainer<double> vals;
    DefaultBasisFactory<double> BFactory;
    Teuchos::RCP<Basis<double> > quadBasis = BFactory.create(FIELD_FORM_0, 
                                                             CELL_QUAD, 
                                                             RECONSTRUCTION_SPACE_INCOMPLETE, 
                                                             1, 
                                                             BASIS_FEM_DEFAULT, 
                                                             COORDINATES_CARTESIAN);
    
    // Check VALUE of basis functions
    quadBasis -> getValues(vals, quadNodes, OPERATOR_VALUE);
    for (int i=0; i < vals.getSize(); i++) {
      if (std::abs(vals[i] - basisValues[i]) > INTREPID_TOL) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        
        // Get the multi-index of the value where the error is:
        Teuchos::Array<int> myIndex;
        vals.getMultiIndex(myIndex,i);
        *outStream << " At multi-index { ";
        for(int j = 0; j < vals.getRank(); j++) {
          *outStream << myIndex[j] << " ";
        }
        *outStream << "}  computed value: " << vals[i] 
          << " but reference value: " << basisValues[i] << "\n";
      }
    }
    
    // Check GRAD of basis function
    quadBasis -> getValues(vals, quadNodes, OPERATOR_GRAD);
    for (int i=0; i < vals.getSize(); i++) {
      if (std::abs(vals[i] - basisGrads[i]) > INTREPID_TOL) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        
        // Get the multi-index of the value where the error is:
        Teuchos::Array<int> myIndex;
        vals.getMultiIndex(myIndex,i);
        *outStream << " At multi-index { ";
        for(int j = 0; j < vals.getRank(); j++) {
          *outStream << myIndex[j] << " ";
        }
        *outStream << "}  computed grad component: " << vals[i] 
          << " but reference grad component: " << basisGrads[i] << "\n";
      }
    }
    
    // Check D1 of basis function
    quadBasis -> getValues(vals, quadNodes, OPERATOR_D1);
    for (int i=0; i < vals.getSize(); i++) {
      if (std::abs(vals[i] - basisGrads[i]) > INTREPID_TOL) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        
        // Get the multi-index of the value where the error is:
        Teuchos::Array<int> myIndex;
        vals.getMultiIndex(myIndex,i);
        *outStream << " At multi-index { ";
        for(int j = 0; j < vals.getRank(); j++) {
          *outStream << myIndex[j] << " ";
        }
        *outStream << "}  computed D1 component: " << vals[i] 
          << " but reference D1 component: " << basisGrads[i] << "\n";
      }
    }
    
    
    // Check CURL of basis function
    quadBasis -> getValues(vals, quadNodes, OPERATOR_CURL);
    for (int i=0; i < vals.getSize(); i++) {
      if (std::abs(vals[i] - basisCurls[i]) > INTREPID_TOL) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        
        // Get the multi-index of the value where the error is:
        Teuchos::Array<int> myIndex;
        vals.getMultiIndex(myIndex,i);
        *outStream << " At multi-index { ";
        for(int j = 0; j < vals.getRank(); j++) {
          *outStream << myIndex[j] << " ";
        }
        *outStream << "}  computed curl component: " << vals[i] 
          << " but reference curl component: " << basisCurls[i] << "\n";
      }
    }
    
    // Check D2 of basis function
    quadBasis -> getValues(vals, quadNodes, OPERATOR_D2);
    for (int i=0; i < vals.getSize(); i++) {
      if (std::abs(vals[i] - basisD2[i]) > INTREPID_TOL) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        
        // Get the multi-index of the value where the error is:
        Teuchos::Array<int> myIndex;
        vals.getMultiIndex(myIndex,i);
        *outStream << " At multi-index { ";
        for(int j = 0; j < vals.getRank(); j++) {
          *outStream << myIndex[j] << " ";
        }
        *outStream << "}  computed D2 component: " << vals[i] 
          << " but reference D2 component: " << basisD2[i] << "\n";
      }
    }
    
    // Check all higher derivatives - must be zero
    for(EOperator op = OPERATOR_D3; op < OPERATOR_MAX; op++) {
      quadBasis -> getValues(vals, quadNodes, op);
      for (int i=0; i < vals.getSize(); i++) {
        if (std::abs(vals[i]) > INTREPID_TOL) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          
          // Get the multi-index of the value where the error is and the operator order
          Teuchos::Array<int> myIndex;
          vals.getMultiIndex(myIndex,i);
          int ord = Intrepid::getOperatorOrder(op);
          *outStream << " At multi-index { ";
          for(int j = 0; j < vals.getRank(); j++) {
            *outStream << myIndex[j] << " ";
          }
          *outStream << "}  computed D"<< ord <<" component: " << vals[i] 
            << " but reference D" << ord << " component:  0 \n";
        }
      }
    }    
  }
  
  // Catch unexpected errors
  catch (std::logic_error err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };
  
  
  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";
  
  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  
  return errorFlag;
}
