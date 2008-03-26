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
\brief  Illustrates use of the VarContainer class.
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
  << "|    1) Testing global functions for rank, order, etc                         |\n" \
  << "|    2) Global functions in debug mode                                        |\n" \
  << "|       requires intrepid to be configured with --enable-intrepid-debug       |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n\n";
  
  // Define variables 
  Teuchos::Array<int> partialMult;
  int spaceDim  = 3;
  int derivativeOrder;
  int derivativeCardinality;
  
  cout \
    << "===============================================================================\n"\
    << "| EXAMPLE 1: Using global functions for rank, order, multiplicities, etc.     |\n"\
    << "===============================================================================\n\n";
  
  // Showing field ranks:
  std::cout 
    << "===============================================================================\n"
    << " *** Using getFieldRank : \n"
    << " Field type        Rank \n"
    <<  "-------------------------------------------------------------------------------\n";
  for(EField fieldType = FIELD_FORM_0; fieldType < FIELD_MAX; fieldType++) {
    std::cout << "\t\t" << std::setw(10) << std::left << EFieldToString(fieldType) 
    << "\t\t\t\t\t" <<  getFieldRank(fieldType) << "\n";
  }
  std::cout << "\n\n";
  
  
  // Showing operator orders:
  std::cout << "===============================================================================\n"
    << " *** Using getOperatorOrder : \n"
    << " Operator type     Order \n"
    << "-------------------------------------------------------------------------------\n"
    << std::left;
  for(EOperator opType = OPERATOR_VALUE; opType < OPERATOR_MAX; opType++){
    std::cout
    << "\t\t" << std::setw(10) << EOperatorToString(opType) 
    << "\t\t\t\t\t" << getOperatorOrder(opType) << "\n";
  }
  std::cout << "\n\n";
  
  
  // Showing operator rank
  std::cout 
    << "===============================================================================\n"
    << " *** Using getOperatorRank : \n"
    << " Operator type    Field rank    Space dim     Operator rank \n";
  
  for(int spDim = 1; spDim <=3; spDim++){
    for(int fRank = 0; fRank <=2; fRank++) {
      std::cout << "-------------------------------------------------------------------------------\n";
      for(EOperator opType = OPERATOR_VALUE; opType < OPERATOR_MAX; opType++){
        try{
          std::cout 
          << "\t\t" 
          << std::setw(20) << EOperatorToString(opType) 
          << std::setw(12) << fRank 
          << std::setw(15) << spDim
          << getOperatorRank(opType,fRank,spDim) << "\n";
        }
        catch(logic_error) {
          std::cout << std::left << " undefined \n";
        }
      }
    }
  }
  std::cout << "\n\n";
  
  
  // Showing cardinality of the set of partial derivatives of order k in 1D, 2D and 3D
  std::cout 
    << "===============================================================================\n"
    << " *** Using getDkCardinality : \n"
    <<  "-------------------------------------------------------------------------------\n"
    << "  order Dk  |   1D         2D        3D        <-    cardinality | Dk |\n"
    << "===============================================================================\n";
  for(int ordDk = 1; ordDk <= INTREPID_MAX_DERIVATIVE; ordDk++) {
    std::cout << "\t\t" << std::setw(14) <<  ordDk;
    for(int spDim = 1; spDim <= 3; spDim++){
      std::cout << std::setw(10) << getDkCardinality(ordDk,spDim);
    }
    std::cout <<"\n";
  }
  std::cout << "\n\n";
  
  
  
  // Showing multiplicities of dx, dy and dz for partial derivatives of order 6 in 3D 
  derivativeOrder = getOperatorOrder(OPERATOR_D6);
  derivativeCardinality = getDkCardinality(derivativeOrder,spaceDim);
  
  std::cout 
    << "===============================================================================\n"
    << " *** Using getDkMultiplicities and getDkEnumeration : \n"
    << " Space dimension        = " << spaceDim << "\n"
    << " Derivative order       = " << derivativeOrder << "\n"
    << " Derivative cardinality = " << derivativeCardinality << "\n"
    <<  "-------------------------------------------------------------------------------\n"
    << "  Enumeration    ->   multiplicity   ->   enumeration    \n"
    << "===============================================================================\n";
  
  for(int derivativeEnum = 0; derivativeEnum < derivativeCardinality; derivativeEnum++){
    getDkMultiplicities(partialMult,
                        derivativeEnum,
                        derivativeOrder,
                        spaceDim);
    std::cout 
    << "\t\t" << std::setw(20) << derivativeEnum << " {" 
    << partialMult[0] << "," << partialMult[1] <<","<< partialMult[2] << std::setw(13) << "}"
    << "\t\t" << getDkEnumeration(partialMult[0],partialMult[1],partialMult[2]) << "\n";
  }
  std::cout << "\n\n";
  
  
  
  return 0;
}
