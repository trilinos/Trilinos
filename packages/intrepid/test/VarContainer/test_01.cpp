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

#include "Intrepid_VarContainer.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

using namespace Intrepid;

int main(int argc, char *argv[]) {
  
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
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
    << "|                           Unit Test VarContainer                            |\n" \
    << "|                                                                             |\n" \
    << "|     1) testing related global functions for field and operator rank         |\n" \
    << "|     2) testing global functions for cardinality and enumeration of Dk ops   |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";
  
  int errorFlag  = 0;
  
  try{

    *outStream \
      << "\n"
    << "===============================================================================\n"\
    << "| Test 1: Global function for rank                                           |\n"\
    << "===============================================================================\n\n";

    for(EField fieldType = FIELD_FORM_0; fieldType < FIELD_MAX; fieldType++) {
      int fieldRank = getFieldRank(fieldType);
      int trueRank = -1;
      switch(fieldType) {
        case FIELD_FORM_0:
        case FIELD_FORM_3:
          trueRank = 0;
          break;
        case FIELD_FORM_1:
        case FIELD_FORM_2:
        case FIELD_VECTOR:
          trueRank = 1;
          break;
        case FIELD_TENSOR:
          trueRank = 2;
          break;
        default:
          errorFlag++;
          *outStream << std::setw(50) << "vvvv----FAILURE!" << "\n";
          *outStream << "     Field type = " << EFieldToString(fieldType) << " should not be reached\n";
      }
      if( trueRank != fieldRank ) {
        errorFlag++;
        *outStream << std::setw(50) << "vvvv----FAILURE!" << "\n";
        *outStream << "     True rank of " << EFieldToString(fieldType) << " is " << trueRank << "\n";
        *outStream << " rank by getFieldRank = " << fieldRank << "\n";
      }
    }
    
    *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 2: Global function for operator rank                                   |\n"\
      << "===============================================================================\n";

    // counts number of undefined cases for operator rank
    int undefined=0;
    for(int spDim = 1; spDim <=3; spDim++){
      for(int fRank = 0; fRank <=2; fRank++) {
        for(EOperator opType = OPERATOR_VALUE; opType < OPERATOR_MAX; opType++){
          try{
            int opRank=getOperatorRank(opType,fRank,spDim);
            
            // OPERATOR_VALUE should always have rank 0
            if( (opType == OPERATOR_VALUE) && !(opRank == 0) ) {
              errorFlag++;
              *outStream << std::setw(50) << "vvvv----FAILURE!" << "\n";
              *outStream << " Operator rank for OPERATOR_VALUE is 0 for all field ranks \n";
              *outStream << " Operator rank evaluated to " << opRank 
                << " when applied to field with rank " << fRank<< " in dimension " << spDim << "\n";
            }
            
            // For all these operators rank has to be 1
            if( ( (opType == OPERATOR_GRAD) ||
                  (opType == OPERATOR_D1)   ||
                  (opType == OPERATOR_D2)   ||
                  (opType == OPERATOR_D3)   ||
                  (opType == OPERATOR_D4)   ||
                  (opType == OPERATOR_D5)   ||
                  (opType == OPERATOR_D6)   ||
                  (opType == OPERATOR_D7)   ||
                  (opType == OPERATOR_D8)   ||
                  (opType == OPERATOR_D9)   ||
                  (opType == OPERATOR_D10) ) && !(opRank == 1) ) {
              errorFlag++;
              *outStream << std::setw(50) << "vvvv----FAILURE!" << "\n";
              //*outStream << " Operator rank for " << OperatorNames[opType] << " is 1 for all field ranks \n";
              *outStream << " Operator rank for " << EOperatorToString(opType) << " is 1 for all field ranks \n";
              *outStream << " Operator rank evaluated to " << opRank 
                << " when applied to field with rank " << fRank<< " in dimension " << spDim << "\n";
            }
            
            // Rank is -1 for DIV in all dimensions and CURL in 2 dimensions applied to fields of rank 1,2
            if( ( ( (opType == OPERATOR_DIV) || (opType == OPERATOR_CURL) ) && 
                  (spDim == 2) && 
                  (fRank == 1 || fRank == 2)  && 
                  !(opRank == -1) )  ) {
              errorFlag++;
              *outStream << std::setw(50) << "vvvv----FAILURE!" << "\n";
              //*outStream << " Operator rank for " << OperatorNames[opType] 
              *outStream << " Operator rank for " << EOperatorToString(opType) 
                << " when applied to field with rank " << fRank<< " in dimension " << spDim 
                << " has to be -1 \n";
              *outStream << " Operator rank evaluated to " << opRank << "\n";
            }
                
                // Rank is 1 for CURL in 2D when applied to field with rank 0
                if( (opType == OPERATOR_CURL && fRank == 0) && !(opRank == 1) ) {
                  errorFlag++;
                  *outStream << std::setw(50) << "vvvv----FAILURE!" << "\n";
                  //*outStream << " Operator rank for " << OperatorNames[opType] 
                  *outStream << " Operator rank for " << EOperatorToString(opType) 
                    << " when applied to field with rank " << fRank<< " in dimension " << spDim 
                    << " has to be -1 \n";
                  *outStream << " Operator rank evaluated to " << opRank << "\n";
                }
          }
          catch(std::logic_error) {
            undefined++;
          }
        }
      }
    }
    
    // undefined must equal 31
    if( undefined != 31 ) {
      errorFlag++;
      *outStream << std::setw(50) << "vvvv----FAILURE!" << "\n";
      *outStream << " Number of undefined operator rank cases is  31 \n";
      *outStream << " Number generated by test is = " << undefined << "\n";
    }
    
  } //try 
  
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
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
