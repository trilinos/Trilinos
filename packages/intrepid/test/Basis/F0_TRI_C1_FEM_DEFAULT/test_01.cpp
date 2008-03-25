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
\brief  Unit tests for the Point class.
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
  << "|                 Unit Test (Basis_F0_TRI_C1_FEM_DEFAULT)                     |\n" \
  << "|                                                                             |\n" \
  << "|     1) Basis creation, computation of basis function values                 |\n" \
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
  double basisvals[] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  
  VarContainer<double> vals;
  DefaultBasisFactory<double> BFactory;
  Teuchos::RCP<Basis<double> > tribasis =  BFactory.create(
    FIELD_FORM_0, CELL_TRI, RECONSTRUCTION_SPACE_COMPLETE, 1, BASIS_FEM_DEFAULT, COORDINATES_CARTESIAN);
  Teuchos::Array< Point<double> > points01;
  Point<double> pt(2, FRAME_REFERENCE);
  points01.assign(3, pt);
  double ptarray[2];
  ptarray[0] = 0.0; ptarray[1] = 0.0;
  (points01[0]).setCoordinates(ptarray, 2);
  ptarray[0] = 1.0; ptarray[1] = 0.0;
  (points01[1]).setCoordinates(ptarray, 2);
  ptarray[0] = 0.0; ptarray[1] = 1.0;
  (points01[2]).setCoordinates(ptarray, 2);

  try{

    try {
      tribasis->getValues(vals, points01, OPERATOR_DIV);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };

    tribasis->getValues(vals, points01, OPERATOR_VALUE);
    *outStream << vals << "\n";
    tribasis->getValues(vals, points01, OPERATOR_GRAD);
    *outStream << vals << "\n";
    tribasis->getValues(vals, points01, OPERATOR_CURL);
    *outStream << vals << "\n";

    for (EOperator op = OPERATOR_D1; op < OPERATOR_MAX; op++) {
      tribasis->getValues(vals, points01, op);
      *outStream << vals << "\n";
    }

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    try {
      LocalDofTag myTag = {{3,0,0,1}};
      tribasis->getLocalDofEnumeration(myTag);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    try {
      LocalDofTag myTag = {{1,1,1,0}};
      tribasis->getLocalDofEnumeration(myTag);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    try {
      LocalDofTag myTag = {{0,3,0,0}};
      tribasis->getLocalDofEnumeration(myTag);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    try {
      int bfId = 4;
      tribasis->getLocalDofTag(bfId);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    try {
      int bfId = -1;
      tribasis->getLocalDofTag(bfId);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    Teuchos::Array<LocalDofTag> allTags;
    tribasis->getAllLocalDofTags(allTags);
    for (unsigned i=0; i<allTags.size(); i++) {
      LocalDofTag myTag = allTags[i];
      int bfId = tribasis->getLocalDofEnumeration(myTag);
      myTag = tribasis->getLocalDofTag(bfId);
      *outStream << "Local Bf Id = " << bfId << "  -->  " << "DoF Tag = ";
      *outStream << "{" << myTag.tag_[0] << ", " << myTag.tag_[1] << ", " << myTag.tag_[2] << ", " << myTag.tag_[3] << "}\n";
    }

  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

  if (TestForException_getThrowNumber() != endThrowNumber)
    errorFlag++;

  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 2: correctness of basis function values                                |\n"\
  << "===============================================================================\n";

  outStream->precision(20);

  try{
    tribasis->getValues(vals, points01, OPERATOR_VALUE);
    for (int i=0; i<vals.getSize(); i++) {
      *outStream << "  Computed value: " << vals[i] << "   Reference value: " << basisvals[i] << "\n";
      if (std::abs(vals[i]-basisvals[i]) > INTREPID_TOL)
        errorFlag++;
    }
  }
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
