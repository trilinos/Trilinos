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
\brief  Unit test for the ArrayTools class.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

using namespace std;
using namespace Intrepid;

#define INTREPID_TEST_COMMAND( S )                                                                                  \
{                                                                                                                   \
  try {                                                                                                             \
    S ;                                                                                                             \
  }                                                                                                                 \
  catch (std::logic_error err) {                                                                                    \
      *outStream << "Expected Error ----------------------------------------------------------------\n";            \
      *outStream << err.what() << '\n';                                                                             \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";    \
  };                                                                                                                \
}


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
  << "|                       Unit Test (ArrayTools)                                |\n" \
  << "|                                                                             |\n" \
  << "|     1) Array operations: multiplication, contractions                       |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n";


  int errorFlag = 0;
#ifdef HAVE_INTREPID_DEBUG
  int beginThrowNumber = TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + 162;
#endif

  typedef ArrayTools art; 
  typedef RealSpaceTools<double> rst; 
#ifdef HAVE_INTREPID_DEBUG
  ArrayTools atools;
#endif

  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 1: exceptions                                                          |\n"\
  << "===============================================================================\n";

  try{

#ifdef HAVE_INTREPID_DEBUG
    FieldContainer<double> a_2_2(2, 2);
    FieldContainer<double> a_10_2(10, 2);
    FieldContainer<double> a_10_3(10, 3);
    FieldContainer<double> a_10_2_2(10, 2, 2);
    FieldContainer<double> a_10_2_3(10, 2, 3);
    FieldContainer<double> a_10_3_2(10, 3, 2);
    FieldContainer<double> a_9_2_2(9, 2, 2);

    *outStream << "-> contractScalar:\n";
    INTREPID_TEST_COMMAND( atools.contractScalar<double>(a_2_2, a_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractScalar<double>(a_2_2, a_10_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractScalar<double>(a_2_2, a_10_2_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractScalar<double>(a_10_2_2, a_9_2_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractScalar<double>(a_10_2_2, a_10_2_2, a_10_2_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractScalar<double>(a_9_2_2, a_10_3_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractScalar<double>(a_10_3_2, a_10_2_2, a_10_3_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractScalar<double>(a_10_2_2, a_10_2_2, a_10_3_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractScalar<double>(a_10_2_3, a_10_2_2, a_10_3_2, COMP_ENGINE_MAX) );
    INTREPID_TEST_COMMAND( atools.contractScalar<double>(a_10_2_3, a_10_2_2, a_10_3_2, COMP_CPP) );

    FieldContainer<double> a_10_2_2_2(10, 2, 2, 2);
    FieldContainer<double> a_9_2_2_2(9, 2, 2, 2);
    FieldContainer<double> a_10_3_2_2(10, 3, 2, 2);
    FieldContainer<double> a_10_2_3_2(10, 2, 3, 2);
    FieldContainer<double> a_10_2_2_3(10, 2, 2, 3);

    *outStream << "-> contractVector:\n";
    INTREPID_TEST_COMMAND( atools.contractVector<double>(a_2_2, a_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVector<double>(a_2_2, a_10_2_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVector<double>(a_2_2, a_10_2_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVector<double>(a_10_2_2, a_9_2_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVector<double>(a_10_2_2, a_10_2_2_2, a_10_2_3_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVector<double>(a_10_2_2, a_10_2_2_2, a_10_2_2_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVector<double>(a_9_2_2, a_10_2_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVector<double>(a_10_3_2, a_10_2_2_2, a_10_3_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVector<double>(a_10_2_2, a_10_2_2_2, a_10_3_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVector<double>(a_10_2_3, a_10_2_2_2, a_10_3_2_2, COMP_ENGINE_MAX) );
    INTREPID_TEST_COMMAND( atools.contractVector<double>(a_10_2_3, a_10_2_2_2, a_10_3_2_2, COMP_CPP) );

    FieldContainer<double> a_10_2_2_2_2(10, 2, 2, 2, 2);
    FieldContainer<double> a_9_2_2_2_2(9, 2, 2, 2, 2);
    FieldContainer<double> a_10_3_2_2_2(10, 3, 2, 2, 2);
    FieldContainer<double> a_10_2_3_2_2(10, 2, 3, 2, 2);
    FieldContainer<double> a_10_2_2_3_2(10, 2, 2, 3, 2);
    FieldContainer<double> a_10_2_2_2_3(10, 2, 2, 2, 3);

    *outStream << "-> contractTensor:\n";
    INTREPID_TEST_COMMAND( atools.contractTensor<double>(a_2_2, a_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensor<double>(a_2_2, a_10_2_2_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensor<double>(a_2_2, a_10_2_2_2_2, a_10_2_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensor<double>(a_10_2_2, a_9_2_2_2_2, a_10_2_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensor<double>(a_10_2_2, a_10_2_2_2_2, a_10_2_3_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensor<double>(a_10_2_2, a_10_2_2_2_2, a_10_2_2_3_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensor<double>(a_10_2_2, a_10_2_2_2_2, a_10_2_2_2_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensor<double>(a_9_2_2, a_10_2_2_2_2, a_10_2_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensor<double>(a_10_3_2, a_10_2_2_2_2, a_10_2_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensor<double>(a_10_2_2, a_10_2_2_2_2, a_10_3_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensor<double>(a_10_2_3, a_10_2_2_2_2, a_10_3_2_2_2, COMP_ENGINE_MAX) );
    INTREPID_TEST_COMMAND( atools.contractTensor<double>(a_10_2_3, a_10_2_2_2_2, a_10_3_2_2_2, COMP_CPP) );

    FieldContainer<double> a_9_2(9, 2);
    FieldContainer<double> a_10_1(10, 1);

    *outStream << "-> contractScalarData:\n";
    INTREPID_TEST_COMMAND( atools.contractScalarData<double>(a_2_2, a_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractScalarData<double>(a_2_2, a_10_2_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractScalarData<double>(a_10_2_2, a_10_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractScalarData<double>(a_2_2, a_9_2_2, a_10_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractScalarData<double>(a_2_2, a_10_2_2, a_10_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractScalarData<double>(a_9_2, a_10_2_2, a_10_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractScalarData<double>(a_10_2, a_10_3_2, a_10_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractScalarData<double>(a_10_2, a_10_2_2, a_10_2, COMP_ENGINE_MAX) );
    INTREPID_TEST_COMMAND( atools.contractScalarData<double>(a_10_2, a_10_2_2, a_10_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractScalarData<double>(a_10_2, a_10_2_2, a_10_1, COMP_CPP) );

    FieldContainer<double> a_10_1_2(10, 1, 2);
    FieldContainer<double> a_10_1_3(10, 1, 3);

    *outStream << "-> contractVectorData:\n";
    INTREPID_TEST_COMMAND( atools.contractVectorData<double>(a_2_2, a_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVectorData<double>(a_2_2, a_10_2_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVectorData<double>(a_10_2_2, a_10_2_2_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVectorData<double>(a_2_2, a_9_2_2_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVectorData<double>(a_2_2, a_10_2_3_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVectorData<double>(a_2_2, a_10_2_2_3, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVectorData<double>(a_9_2, a_10_2_2_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVectorData<double>(a_10_2, a_10_3_2_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVectorData<double>(a_10_2, a_10_2_2_2, a_10_2_2, COMP_ENGINE_MAX) );
    INTREPID_TEST_COMMAND( atools.contractVectorData<double>(a_10_2, a_10_2_2_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractVectorData<double>(a_10_2, a_10_2_2_2, a_10_1_2, COMP_CPP) );

    FieldContainer<double> a_10_1_2_2(10, 1, 2, 2);

    *outStream << "-> contractTensorData:\n";
    INTREPID_TEST_COMMAND( atools.contractTensorData<double>(a_2_2, a_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensorData<double>(a_2_2, a_10_2_2_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensorData<double>(a_10_2_2, a_10_2_2_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensorData<double>(a_2_2, a_9_2_2_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensorData<double>(a_2_2, a_10_2_3_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensorData<double>(a_2_2, a_10_2_2_3_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensorData<double>(a_2_2, a_10_2_2_2_3, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensorData<double>(a_9_2, a_10_2_2_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensorData<double>(a_10_2, a_10_3_2_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensorData<double>(a_10_2, a_10_2_2_2_2, a_10_2_2_2, COMP_ENGINE_MAX) );
    INTREPID_TEST_COMMAND( atools.contractTensorData<double>(a_10_2, a_10_2_2_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractTensorData<double>(a_10_2, a_10_2_2_2_2, a_10_1_2_2, COMP_CPP) );

    FieldContainer<double> a_2_3_2_2(2, 3, 2, 2);
    FieldContainer<double> a_2_2_2_2(2, 2, 2, 2);
    FieldContainer<double> a_2_10(2, 10);
    FieldContainer<double> a_2(2);

    *outStream << "-> multiplyScalarData:\n";
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_2_2, a_10_2_2, a_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_2_2, a_2_2, a_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_2_2, a_2_2, a_10_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_2_2, a_2_2, a_10_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_2_2, a_10_3, a_10_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_9_2_2_2_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_3_2_2_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_2_3_2_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_2_2_3_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_2_2_2_3, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_2_2_2_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_2_2_2_2, a_10_1, a_10_2_2_2_2) );
    //
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_2_2, a_10_2_2, a_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_2_2_2, a_2_2, a_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_2_2_2, a_2_2, a_10_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_2_2, a_10_2, a_2_10) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_2_2_2_2, a_9_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_3_2_2_2, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_2_3_2_2, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_2_2_3_2, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_2_2_2_3, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_2_2_2_2, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyScalarData<double>(a_10_2_2_2_2, a_10_1, a_2_2_2_2) );

    *outStream << "-> divideByScalarData:\n";
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_2_2, a_10_2_2, a_2_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_2_2, a_2_2, a_2_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_2_2, a_2_2, a_10_2_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_2_2, a_2_2, a_10_2_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_2_2, a_10_3, a_10_2_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_9_2_2_2_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_3_2_2_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_2_3_2_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_2_2_3_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_2_2_2_3, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_2_2_2_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_2_2_2_2, a_10_1, a_10_2_2_2_2) );
    //
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_2_2, a_10_2_2, a_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_2_2_2, a_2_2, a_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_2_2_2, a_2_2, a_10_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_2_2, a_10_2, a_2_10) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_2_2_2_2, a_9_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_3_2_2_2, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_2_3_2_2, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_2_2_3_2, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_2_2_2_3, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_2_2_2_2, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.divideByScalarData<double>(a_10_2_2_2_2, a_10_1, a_2_2_2_2) );

    FieldContainer<double> a_10_10_2(10, 10, 2);
    FieldContainer<double> a_10_10_2_2(10, 10, 2, 2);
    FieldContainer<double> a_9_10_2(9, 10, 2);
    FieldContainer<double> a_10_10_3_2(10, 10, 3, 2);
    FieldContainer<double> a_10_2_10_2(10, 2, 10, 2);
    FieldContainer<double> a_10_2_10_3(10, 2, 10, 3);
    FieldContainer<double> a_10_3_2_10(10, 3, 2, 10);
    FieldContainer<double> a_2_10_2_2(2, 10, 2, 2);
    FieldContainer<double> a_2_10_3_2(2, 10, 3, 2);
    FieldContainer<double> a_2_10_2_3(2, 10, 2, 3);

    *outStream << "-> multiplyVectorData:\n";
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_2_2, a_2_2, a_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_2_2, a_10_2_2, a_10_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_2, a_10_2_2, a_10_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_2_2_2, a_9_2_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_2_2_2, a_10_3_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_2_2_2, a_10_2_3, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_9_2_2_2, a_10_2_3, a_10_2_2_3_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_3_2_2, a_10_2_3, a_10_2_2_3_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_2_3_2, a_10_2_3, a_10_2_2_3_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_2_2_3, a_10_2_3, a_10_2_2_3_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_2_2_2, a_10_2_3, a_10_2_2_3_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_2_2_2, a_10_1_3, a_10_2_2_3_2) );
    //
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_2_2, a_2_2, a_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_2_2, a_10_2_2, a_10_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_2_2_2_2, a_10_2_2, a_10_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_2_2_2, a_10_2_2, a_2_10_2_3) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_2_2_2, a_10_10_2, a_2_10_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_2_10_2, a_10_10_2, a_2_10_3_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_2_10_2, a_9_10_2, a_2_10_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_2_10_2, a_10_10_2, a_10_10_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_2_10_3, a_10_10_2, a_2_10_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_2_10_2, a_10_10_2, a_2_10_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyVectorData<double>(a_10_2_10_2, a_10_1_2, a_2_10_2_2) );

    FieldContainer<double> a_10_10_2_2_2(10, 10, 2, 2, 2);
    FieldContainer<double> a_10_10_3_2_2(10, 10, 3, 2, 2);
    FieldContainer<double> a_10_10_2_3_2(10, 10, 2, 3, 2);
    FieldContainer<double> a_10_10_2_3(10, 10, 2, 3);
    FieldContainer<double> a_2_10_2(2, 10, 2);
    FieldContainer<double> a_2_10_3(2, 10, 3);
    FieldContainer<double> a_10_2_10_2_2(10, 2, 10, 2, 2);
    FieldContainer<double> a_10_2_10_3_2(10, 2, 10, 3, 2);
    FieldContainer<double> a_10_2_10_2_3(10, 2, 10, 2, 3);

    *outStream << "-> multiplyTensorData:\n";
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_2_2, a_10_2_2_2_2, a_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_2_2, a_2_2, a_10_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_2_2, a_2_2, a_10_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_2_2, a_10_2, a_10_2_2_2, 'r') );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_2_2, a_9_2_2, a_10_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_2_2, a_10_3_2, a_10_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_2_2, a_10_2_3, a_10_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_2_2, a_10_2_2_3, a_10_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_2_2_2, a_10_2_2_2, a_10_2_2_2_3) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_9_2_2_2_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_3_2_2_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_3_2_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_2_3_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_2_2_3, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_2_2_2, a_10_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_2_2_2, a_10_2_2_2, a_10_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_2_2_2, a_10_1_2_2, a_10_2_2_2_2) );
    //
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_2, a_10_2_2_2_2, a_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_2, a_10_2_2_2, a_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_2_2_2, a_10_2_2_2, a_2_10_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_2_2, a_10_2_2_2, a_2_10_2, 'f') );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_9_2_2_2, a_10_2_2_2, a_2_10_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_2_2, a_10_2_2_2, a_2_10_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_10_2, a_10_10_2_2, a_2_10_3) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_10_2, a_10_10_2_3, a_2_10_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_10_2_2, a_10_10_2_2, a_2_10_2_3) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_3_2_2, a_10_10_2_2, a_2_10_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_2_2, a_10_10_2_2, a_2_10_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_10_3_2, a_10_10_2_2, a_2_10_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_10_2_3, a_10_10_2_2, a_2_10_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_10_2_2, a_10_10_2_2, a_2_10_2_2) );
    INTREPID_TEST_COMMAND( atools.multiplyTensorData<double>(a_10_2_10_2_2, a_10_1_2_2, a_2_10_2_2) );

    *outStream << "-> cloneValues:\n";
    INTREPID_TEST_COMMAND( atools.cloneValues<double>(a_10_2_2_2, a_2) );
    INTREPID_TEST_COMMAND( atools.cloneValues<double>(a_10_2_2_2, a_10_2) );
    INTREPID_TEST_COMMAND( atools.cloneValues<double>(a_10_3_2, a_2_2) );
    INTREPID_TEST_COMMAND( atools.cloneValues<double>(a_10_2_2_2_2, a_2_3_2_2) );
    INTREPID_TEST_COMMAND( atools.cloneValues<double>(a_10_2_2_3_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.cloneValues<double>(a_10_2_2_2_3, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.cloneValues<double>(a_10_2_2, a_2_2) );
    INTREPID_TEST_COMMAND( atools.cloneValues<double>(a_10_2_2_2_2, a_2_2_2_2) );

    *outStream << "-> cloneScaleValues:\n";
    INTREPID_TEST_COMMAND( atools.cloneScaleValues<double>(a_10_2_2_2, a_2, a_2) );
    INTREPID_TEST_COMMAND( atools.cloneScaleValues<double>(a_10_2_2_2, a_10_2, a_2) );
    INTREPID_TEST_COMMAND( atools.cloneScaleValues<double>(a_10_2_2_2, a_10_2, a_10_2) );
    INTREPID_TEST_COMMAND( atools.cloneScaleValues<double>(a_10_2_2, a_9_2, a_10_2) );
    INTREPID_TEST_COMMAND( atools.cloneScaleValues<double>(a_10_2_2, a_10_3, a_10_2) );
    INTREPID_TEST_COMMAND( atools.cloneScaleValues<double>(a_10_3_2_2_2, a_10_3, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.cloneScaleValues<double>(a_10_2_3_2_2, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.cloneScaleValues<double>(a_10_2_2_3_2, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.cloneScaleValues<double>(a_10_2_2_2_3, a_10_2, a_2_2_2_2) );
    INTREPID_TEST_COMMAND( atools.cloneScaleValues<double>(a_10_2_2, a_10_2, a_2_2) );
    INTREPID_TEST_COMMAND( atools.cloneScaleValues<double>(a_10_2_2_2_2, a_10_2, a_2_2_2_2) );

    *outStream << "-> scaleValues:\n";
    INTREPID_TEST_COMMAND( atools.scaleValues<double>(a_10_2_2_2, a_2) );
    INTREPID_TEST_COMMAND( atools.scaleValues<double>(a_10_2, a_2_2) );
    INTREPID_TEST_COMMAND( atools.scaleValues<double>(a_10_2_2, a_2_2) );
    INTREPID_TEST_COMMAND( atools.scaleValues<double>(a_10_3_2, a_10_2) );
    INTREPID_TEST_COMMAND( atools.scaleValues<double>(a_10_3_2_2, a_10_2) );
    INTREPID_TEST_COMMAND( atools.scaleValues<double>(a_10_3_2_2_2, a_10_2) );
    INTREPID_TEST_COMMAND( atools.scaleValues<double>(a_10_2_2_2_2, a_10_2) );
#endif

  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

#ifdef HAVE_INTREPID_DEBUG
  if (TestForException_getThrowNumber() != endThrowNumber)
    errorFlag++;
#endif

  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 2: correctness of math operations                                      |\n"\
  << "===============================================================================\n";

  outStream->precision(20);

  try {
      { // start scope
      *outStream << "\n************ Checking contractScalar ************\n";

      int c=5, p=9, l=3, r=7;

      FieldContainer<double> in_c_l_p(c, l, p);
      FieldContainer<double> in_c_r_p(c, r, p);
      FieldContainer<double> out1_c_l_r(c, l, r);
      FieldContainer<double> out2_c_l_r(c, l, r);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
      for (int i=0; i<in_c_l_p.size(); i++) {
        in_c_l_p[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_c_r_p.size(); i++) {
        in_c_r_p[i] = Teuchos::ScalarTraits<double>::random();
      }

      art::contractScalar<double>(out1_c_l_r, in_c_l_p, in_c_r_p, COMP_CPP);
      art::contractScalar<double>(out2_c_l_r, in_c_l_p, in_c_r_p, COMP_BLAS);
      rst::subtract(&out1_c_l_r[0], &out2_c_l_r[0], out2_c_l_r.size());
      if (rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractScalar (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking contractVector ************\n";

      int c=5, p=9, l=3, r=7, d=13;

      FieldContainer<double> in_c_l_p_d(c, l, p, d);
      FieldContainer<double> in_c_r_p_d(c, r, p, d);
      FieldContainer<double> out1_c_l_r(c, l, r);
      FieldContainer<double> out2_c_l_r(c, l, r);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
      for (int i=0; i<in_c_l_p_d.size(); i++) {
        in_c_l_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_c_r_p_d.size(); i++) {
        in_c_r_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }

      art::contractVector<double>(out1_c_l_r, in_c_l_p_d, in_c_r_p_d, COMP_CPP);
      art::contractVector<double>(out2_c_l_r, in_c_l_p_d, in_c_r_p_d, COMP_BLAS);

      rst::subtract(&out1_c_l_r[0], &out2_c_l_r[0], out2_c_l_r.size());
      if (rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractVector (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking contractTensor ************\n";

      int c=5, p=9, l=3, r=7, d1=13, d2=5;

      FieldContainer<double> in_c_l_p_d_d(c, l, p, d1, d2);
      FieldContainer<double> in_c_r_p_d_d(c, r, p, d1, d2);
      FieldContainer<double> out1_c_l_r(c, l, r);
      FieldContainer<double> out2_c_l_r(c, l, r);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
      for (int i=0; i<in_c_l_p_d_d.size(); i++) {
        in_c_l_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_c_r_p_d_d.size(); i++) {
        in_c_r_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }

      art::contractTensor<double>(out1_c_l_r, in_c_l_p_d_d, in_c_r_p_d_d, COMP_CPP);
      art::contractTensor<double>(out2_c_l_r, in_c_l_p_d_d, in_c_r_p_d_d, COMP_BLAS);

      rst::subtract(&out1_c_l_r[0], &out2_c_l_r[0], out2_c_l_r.size());
      if (rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractTensor (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking contractScalarData ************\n";

      int c=5, p=9, l=7;

      FieldContainer<double> in_c_l_p(c, l, p);
      FieldContainer<double> data_c_p(c, p);
      FieldContainer<double> data_c_1(c, 1);
      FieldContainer<double> out1_c_l(c, l);
      FieldContainer<double> out2_c_l(c, l);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
      for (int i=0; i<in_c_l_p.size(); i++) {
        in_c_l_p[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<data_c_p.size(); i++) {
        data_c_p[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<data_c_1.size(); i++) {
        data_c_1[i] = Teuchos::ScalarTraits<double>::random();
      }

      // nonconstant data
      art::contractScalarData<double>(out1_c_l, in_c_l_p, data_c_p, COMP_CPP);
      art::contractScalarData<double>(out2_c_l, in_c_l_p, data_c_p, COMP_BLAS);
      rst::subtract(&out1_c_l[0], &out2_c_l[0], out2_c_l.size());
      if (rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractScalarData (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) << "\n\n";
        errorFlag = -1000;
      }
      // constant data
      art::contractScalarData<double>(out1_c_l, in_c_l_p, data_c_1, COMP_CPP);
      art::contractScalarData<double>(out2_c_l, in_c_l_p, data_c_1, COMP_BLAS);
      rst::subtract(&out1_c_l[0], &out2_c_l[0], out2_c_l.size());
      if (rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractScalarData (2): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking contractVectorData ************\n";

      int c=5, p=9, l=7, d=3;

      FieldContainer<double> in_c_l_p_d(c, l, p, d);
      FieldContainer<double> data_c_p_d(c, p, d);
      FieldContainer<double> data_c_1_d(c, 1, d);
      FieldContainer<double> out1_c_l(c, l);
      FieldContainer<double> out2_c_l(c, l);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
      for (int i=0; i<in_c_l_p_d.size(); i++) {
        in_c_l_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<data_c_p_d.size(); i++) {
        data_c_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<data_c_1_d.size(); i++) {
        data_c_1_d[i] = Teuchos::ScalarTraits<double>::random();
      }

      // nonconstant data
      art::contractVectorData<double>(out1_c_l, in_c_l_p_d, data_c_p_d, COMP_CPP);
      art::contractVectorData<double>(out2_c_l, in_c_l_p_d, data_c_p_d, COMP_BLAS);
      rst::subtract(&out1_c_l[0], &out2_c_l[0], out2_c_l.size());
      if (rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractVectorData (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) << "\n\n";
        errorFlag = -1000;
      }
      // constant data
      art::contractVectorData<double>(out1_c_l, in_c_l_p_d, data_c_1_d, COMP_CPP);
      art::contractVectorData<double>(out2_c_l, in_c_l_p_d, data_c_1_d, COMP_BLAS);
      rst::subtract(&out1_c_l[0], &out2_c_l[0], out2_c_l.size());
      if (rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractVectorData (2): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking contractTensorData ************\n";

      int c=5, p=9, l=7, d1=3, d2=13;

      FieldContainer<double> in_c_l_p_d_d(c, l, p, d1, d2);
      FieldContainer<double> data_c_p_d_d(c, p, d1, d2);
      FieldContainer<double> data_c_1_d_d(c, 1, d1, d2);
      FieldContainer<double> out1_c_l(c, l);
      FieldContainer<double> out2_c_l(c, l);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
      for (int i=0; i<in_c_l_p_d_d.size(); i++) {
        in_c_l_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<data_c_p_d_d.size(); i++) {
        data_c_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<data_c_1_d_d.size(); i++) {
        data_c_1_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }

      // nonconstant data
      art::contractTensorData<double>(out1_c_l, in_c_l_p_d_d, data_c_p_d_d, COMP_CPP);
      art::contractTensorData<double>(out2_c_l, in_c_l_p_d_d, data_c_p_d_d, COMP_BLAS);
      rst::subtract(&out1_c_l[0], &out2_c_l[0], out2_c_l.size());
      if (rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractTensorData (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) << "\n\n";
        errorFlag = -1000;
      }
      // constant data
      art::contractTensorData<double>(out1_c_l, in_c_l_p_d_d, data_c_1_d_d, COMP_CPP);
      art::contractTensorData<double>(out2_c_l, in_c_l_p_d_d, data_c_1_d_d, COMP_BLAS);
      rst::subtract(&out1_c_l[0], &out2_c_l[0], out2_c_l.size());
      if (rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractTensorData (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking multiplyScalarData (branch 1) ************\n";

      int c=5, p=9, f=7, d1=7, d2=13;

      FieldContainer<double> in_c_f_p(c, f, p);
      FieldContainer<double> in_c_f_p_d(c, f, p, d1);
      FieldContainer<double> in_c_f_p_d_d(c, f, p, d1, d2);
      FieldContainer<double> data_c_p(c, p);
      FieldContainer<double> datainv_c_p(c, p);
      FieldContainer<double> data_c_1(c, 1);
      FieldContainer<double> datainv_c_1(c, 1);
      FieldContainer<double> out_c_f_p(c, f, p);
      FieldContainer<double> outi_c_f_p(c, f, p);
      FieldContainer<double> out_c_f_p_d(c, f, p, d1);
      FieldContainer<double> outi_c_f_p_d(c, f, p, d1);
      FieldContainer<double> out_c_f_p_d_d(c, f, p, d1, d2);
      FieldContainer<double> outi_c_f_p_d_d(c, f, p, d1, d2);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
      for (int i=0; i<in_c_f_p.size(); i++) {
        in_c_f_p[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_c_f_p_d.size(); i++) {
        in_c_f_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_c_f_p_d_d.size(); i++) {
        in_c_f_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<data_c_p.size(); i++) {
        data_c_p[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_p[i] = 1.0 / data_c_p[i];
      }
      for (int i=0; i<data_c_1.size(); i++) {
        data_c_1[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_1[i] = 1.0 / data_c_1[i];
      }

      art::multiplyScalarData<double>(out_c_f_p, data_c_p, in_c_f_p);
      art::multiplyScalarData<double>(outi_c_f_p, datainv_c_p, out_c_f_p);
      rst::subtract(&outi_c_f_p[0], &in_c_f_p[0], outi_c_f_p.size());
      if (rst::vectorNorm(&outi_c_f_p[0], outi_c_f_p.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyScalarData (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyScalarData<double>(out_c_f_p_d, data_c_p, in_c_f_p_d);
      art::multiplyScalarData<double>(outi_c_f_p_d, datainv_c_p, out_c_f_p_d);
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyScalarData (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyScalarData<double>(out_c_f_p_d_d, data_c_p, in_c_f_p_d_d);
      art::multiplyScalarData<double>(outi_c_f_p_d_d, datainv_c_p, out_c_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyScalarData (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyScalarData<double>(out_c_f_p_d_d, data_c_1, in_c_f_p_d_d);
      art::multiplyScalarData<double>(outi_c_f_p_d_d, datainv_c_1, out_c_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyScalarData (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
      for (int i=0; i<in_c_f_p_d_d.size(); i++) {
        in_c_f_p_d_d[i] = 1.0;
      }
      for (int i=0; i<data_c_p.size(); i++) {
        data_c_p[i] = 5.0;
      }
      for (int i=0; i<data_c_1.size(); i++) {
        data_c_1[i] = 5.0;
      }
      art::multiplyScalarData<double>(out_c_f_p_d_d, data_c_p, in_c_f_p_d_d);
      if (std::abs(rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) - data_c_p[0]*in_c_f_p_d_d[0]*c*p*f*d1*d2) > zero) {
        *outStream << "\n\nINCORRECT multiplyScalarData (5): check result: "
                   << rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) << " != "
                   << data_c_p[0]*in_c_f_p_d_d[0]*c*f*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::multiplyScalarData<double>(out_c_f_p_d_d, data_c_1, in_c_f_p_d_d);
      if (std::abs(rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) - data_c_1[0]*in_c_f_p_d_d[0]*c*p*f*d1*d2) > zero) {
        *outStream << "\n\nINCORRECT multiplyScalarData (6): check result: "
                   << rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) << " != "
                   << data_c_p[0]*in_c_f_p_d_d[0]*c*f*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking multiplyScalarData (branch 2) ************\n";

      int c=5, p=9, f=7, d1=7, d2=13;

      FieldContainer<double> in_f_p(f, p);
      FieldContainer<double> in_f_p_d(f, p, d1);
      FieldContainer<double> in_f_p_d_d(f, p, d1, d2);
      FieldContainer<double> in_c_f_p(c, f, p);
      FieldContainer<double> in_c_f_p_d(c, f, p, d1);
      FieldContainer<double> in_c_f_p_d_d(c, f, p, d1, d2);
      FieldContainer<double> data_c_p(c, p);
      FieldContainer<double> datainv_c_p(c, p);
      FieldContainer<double> data_c_1(c, 1);
      FieldContainer<double> datainv_c_1(c, 1);
      FieldContainer<double> data_c_p_one(c, p);
      FieldContainer<double> data_c_1_one(c, 1);
      FieldContainer<double> out_c_f_p(c, f, p);
      FieldContainer<double> outi_c_f_p(c, f, p);
      FieldContainer<double> out_c_f_p_d(c, f, p, d1);
      FieldContainer<double> outi_c_f_p_d(c, f, p, d1);
      FieldContainer<double> out_c_f_p_d_d(c, f, p, d1, d2);
      FieldContainer<double> outi_c_f_p_d_d(c, f, p, d1, d2);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
      for (int i=0; i<in_f_p.size(); i++) {
        in_f_p[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_f_p_d.size(); i++) {
        in_f_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_f_p_d_d.size(); i++) {
        in_f_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<data_c_p.size(); i++) {
        data_c_p[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_p[i] = 1.0 / data_c_p[i];
        data_c_p_one[i] = 1.0;
      }
      for (int i=0; i<data_c_1.size(); i++) {
        data_c_1[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_1[i] = 1.0 / data_c_1[i];
        data_c_1_one[i] = 1.0;
      }

      art::multiplyScalarData<double>(out_c_f_p, data_c_p, in_f_p);
      art::multiplyScalarData<double>(outi_c_f_p, datainv_c_p, out_c_f_p);
      art::multiplyScalarData<double>(in_c_f_p, data_c_p_one, in_f_p);
      rst::subtract(&outi_c_f_p[0], &in_c_f_p[0], outi_c_f_p.size());
      if (rst::vectorNorm(&outi_c_f_p[0], outi_c_f_p.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyScalarData (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyScalarData<double>(out_c_f_p_d, data_c_p, in_f_p_d);
      art::multiplyScalarData<double>(outi_c_f_p_d, datainv_c_p, out_c_f_p_d);
      art::multiplyScalarData<double>(in_c_f_p_d, data_c_p_one, in_f_p_d);
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyScalarData (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyScalarData<double>(out_c_f_p_d_d, data_c_p, in_f_p_d_d);
      art::multiplyScalarData<double>(outi_c_f_p_d_d, datainv_c_p, out_c_f_p_d_d);
      art::multiplyScalarData<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyScalarData (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyScalarData<double>(out_c_f_p_d_d, data_c_1, in_f_p_d_d);
      art::multiplyScalarData<double>(outi_c_f_p_d_d, datainv_c_1, out_c_f_p_d_d);
      art::multiplyScalarData<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyScalarData (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
      for (int i=0; i<in_f_p_d_d.size(); i++) {
        in_f_p_d_d[i] = 1.0;
      }
      for (int i=0; i<data_c_p.size(); i++) {
        data_c_p[i] = 5.0;
      }
      for (int i=0; i<data_c_1.size(); i++) {
        data_c_1[i] = 5.0;
      }
      art::multiplyScalarData<double>(out_c_f_p_d_d, data_c_p, in_f_p_d_d);
      if (std::abs(rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) - data_c_p[0]*in_f_p_d_d[0]*c*p*f*d1*d2) > zero) {
        *outStream << "\n\nINCORRECT multiplyScalarData (5): check result: "
                   << rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) << " != "
                   << data_c_p[0]*in_f_p_d_d[0]*c*f*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::multiplyScalarData<double>(out_c_f_p_d_d, data_c_1, in_f_p_d_d);
      if (std::abs(rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) - data_c_1[0]*in_f_p_d_d[0]*c*p*f*d1*d2) > zero) {
        *outStream << "\n\nINCORRECT multiplyScalarData (6): check result: "
                   << rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) << " != "
                   << data_c_1[0]*in_f_p_d_d[0]*c*f*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking divideByScalarData (branch 1) ************\n";

      int c=5, p=9, f=7, d1=7, d2=13;

      FieldContainer<double> in_c_f_p(c, f, p);
      FieldContainer<double> in_c_f_p_d(c, f, p, d1);
      FieldContainer<double> in_c_f_p_d_d(c, f, p, d1, d2);
      FieldContainer<double> data_c_p(c, p);
      FieldContainer<double> datainv_c_p(c, p);
      FieldContainer<double> data_c_1(c, 1);
      FieldContainer<double> datainv_c_1(c, 1);
      FieldContainer<double> out_c_f_p(c, f, p);
      FieldContainer<double> outi_c_f_p(c, f, p);
      FieldContainer<double> out_c_f_p_d(c, f, p, d1);
      FieldContainer<double> outi_c_f_p_d(c, f, p, d1);
      FieldContainer<double> out_c_f_p_d_d(c, f, p, d1, d2);
      FieldContainer<double> outi_c_f_p_d_d(c, f, p, d1, d2);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
      for (int i=0; i<in_c_f_p.size(); i++) {
        in_c_f_p[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_c_f_p_d.size(); i++) {
        in_c_f_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_c_f_p_d_d.size(); i++) {
        in_c_f_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<data_c_p.size(); i++) {
        data_c_p[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_p[i] = 1.0 / data_c_p[i];
      }
      for (int i=0; i<data_c_1.size(); i++) {
        data_c_1[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_1[i] = 1.0 / data_c_1[i];
      }

      art::divideByScalarData<double>(out_c_f_p, data_c_p, in_c_f_p);
      art::divideByScalarData<double>(outi_c_f_p, datainv_c_p, out_c_f_p);
      rst::subtract(&outi_c_f_p[0], &in_c_f_p[0], outi_c_f_p.size());
      if (rst::vectorNorm(&outi_c_f_p[0], outi_c_f_p.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT divideByScalarData (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::divideByScalarData<double>(out_c_f_p_d, data_c_p, in_c_f_p_d);
      art::divideByScalarData<double>(outi_c_f_p_d, datainv_c_p, out_c_f_p_d);
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT divideByScalarData (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::divideByScalarData<double>(out_c_f_p_d_d, data_c_p, in_c_f_p_d_d);
      art::divideByScalarData<double>(outi_c_f_p_d_d, datainv_c_p, out_c_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT divideByScalarData (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::divideByScalarData<double>(out_c_f_p_d_d, data_c_1, in_c_f_p_d_d);
      art::divideByScalarData<double>(outi_c_f_p_d_d, datainv_c_1, out_c_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT divideByScalarData (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
      for (int i=0; i<in_c_f_p_d_d.size(); i++) {
        in_c_f_p_d_d[i] = 1.0;
      }
      for (int i=0; i<data_c_p.size(); i++) {
        data_c_p[i] = 5.0;
      }
      for (int i=0; i<data_c_1.size(); i++) {
        data_c_1[i] = 5.0;
      }
      art::divideByScalarData<double>(out_c_f_p_d_d, data_c_p, in_c_f_p_d_d);
      if (std::abs(rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) - (1.0/data_c_p[0])*in_c_f_p_d_d[0]*c*p*f*d1*d2)/rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT divideByScalarData (5): check result: "
                   << rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) << " != "
                   << (1.0/data_c_p[0])*in_c_f_p_d_d[0]*c*p*f*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::divideByScalarData<double>(out_c_f_p_d_d, data_c_1, in_c_f_p_d_d);
      if (std::abs(rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) - (1.0/data_c_1[0])*in_c_f_p_d_d[0]*c*p*f*d1*d2)/rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT divideByScalarData (6): check result: "
                   << rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) << " != "
                   << (1.0/data_c_p[0])*in_c_f_p_d_d[0]*c*p*f*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking divideByScalarData (branch 2) ************\n";

      int c=5, p=9, f=7, d1=7, d2=13;

      FieldContainer<double> in_f_p(f, p);
      FieldContainer<double> in_f_p_d(f, p, d1);
      FieldContainer<double> in_f_p_d_d(f, p, d1, d2);
      FieldContainer<double> in_c_f_p(c, f, p);
      FieldContainer<double> in_c_f_p_d(c, f, p, d1);
      FieldContainer<double> in_c_f_p_d_d(c, f, p, d1, d2);
      FieldContainer<double> data_c_p(c, p);
      FieldContainer<double> datainv_c_p(c, p);
      FieldContainer<double> data_c_1(c, 1);
      FieldContainer<double> datainv_c_1(c, 1);
      FieldContainer<double> data_c_p_one(c, p);
      FieldContainer<double> data_c_1_one(c, 1);
      FieldContainer<double> out_c_f_p(c, f, p);
      FieldContainer<double> outi_c_f_p(c, f, p);
      FieldContainer<double> out_c_f_p_d(c, f, p, d1);
      FieldContainer<double> outi_c_f_p_d(c, f, p, d1);
      FieldContainer<double> out_c_f_p_d_d(c, f, p, d1, d2);
      FieldContainer<double> outi_c_f_p_d_d(c, f, p, d1, d2);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
      for (int i=0; i<in_f_p.size(); i++) {
        in_f_p[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_f_p_d.size(); i++) {
        in_f_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_f_p_d_d.size(); i++) {
        in_f_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<data_c_p.size(); i++) {
        data_c_p[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_p[i] = 1.0 / data_c_p[i];
        data_c_p_one[i] = 1.0;
      }
      for (int i=0; i<data_c_1.size(); i++) {
        data_c_1[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_1[i] = 1.0 / data_c_1[i];
        data_c_1_one[i] = 1.0;
      }

      art::divideByScalarData<double>(out_c_f_p, data_c_p, in_f_p);
      art::divideByScalarData<double>(outi_c_f_p, datainv_c_p, out_c_f_p);
      art::multiplyScalarData<double>(in_c_f_p, data_c_p_one, in_f_p);
      rst::subtract(&outi_c_f_p[0], &in_c_f_p[0], outi_c_f_p.size());
      if (rst::vectorNorm(&outi_c_f_p[0], outi_c_f_p.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT divideByScalarData (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::divideByScalarData<double>(out_c_f_p_d, data_c_p, in_f_p_d);
      art::divideByScalarData<double>(outi_c_f_p_d, datainv_c_p, out_c_f_p_d);
      art::multiplyScalarData<double>(in_c_f_p_d, data_c_p_one, in_f_p_d);
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT divideByScalarData (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::divideByScalarData<double>(out_c_f_p_d_d, data_c_p, in_f_p_d_d);
      art::divideByScalarData<double>(outi_c_f_p_d_d, datainv_c_p, out_c_f_p_d_d);
      art::multiplyScalarData<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT divideByScalarData (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::divideByScalarData<double>(out_c_f_p_d_d, data_c_1, in_f_p_d_d);
      art::divideByScalarData<double>(outi_c_f_p_d_d, datainv_c_1, out_c_f_p_d_d);
      art::multiplyScalarData<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT divideByScalarData (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
      for (int i=0; i<in_f_p_d_d.size(); i++) {
        in_f_p_d_d[i] = 1.0;
      }
      for (int i=0; i<data_c_p.size(); i++) {
        data_c_p[i] = 5.0;
      }
      for (int i=0; i<data_c_1.size(); i++) {
        data_c_1[i] = 5.0;
      }
      art::divideByScalarData<double>(out_c_f_p_d_d, data_c_p, in_f_p_d_d);
      if (std::abs(rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) - (1.0/data_c_p[0])*in_f_p_d_d[0]*c*p*f*d1*d2)/rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT divideByScalarData (5): check result: "
                   << rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) << " != "
                   << (1.0/data_c_p[0])*in_f_p_d_d[0]*c*p*f*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::divideByScalarData<double>(out_c_f_p_d_d, data_c_1, in_f_p_d_d);
      if (std::abs(rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) - (1.0/data_c_1[0])*in_f_p_d_d[0]*c*p*f*d1*d2)/rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT divideByScalarData (6): check result: "
                   << rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) << " != "
                   << (1.0/data_c_1[0])*in_f_p_d_d[0]*c*p*f*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking multiplyVectorData (branch 1) ************\n";

      int c=5, p=9, f=7, d1=6, d2=13;

      FieldContainer<double> in_c_f_p_d(c, f, p, d1);
      FieldContainer<double> in_c_f_p_d_d(c, f, p, d1, d2);
      FieldContainer<double> data_c_p_d(c, p, d1);
      FieldContainer<double> data_c_1_d(c, 1, d1);
      FieldContainer<double> out_c_f_p(c, f, p);
      FieldContainer<double> out_c_f_p_d(c, f, p, d2);
      double zero = INTREPID_TOL*10000.0;

      // fill with alternating 1's and -1's
      for (int i=0; i<in_c_f_p_d.size(); i++) {
        in_c_f_p_d[i] = std::pow((double)-1.0, (int)i);
      }
      for (int i=0; i<in_c_f_p_d_d.size(); i++) {
        in_c_f_p_d_d[i] = std::pow((double)-1.0, (int)i);
      }
      // fill with 1's
      for (int i=0; i<data_c_p_d.size(); i++) {
        data_c_p_d[i] = 1.0;
      }
      // fill with 1's
      for (int i=0; i<data_c_1_d.size(); i++) {
        data_c_1_d[i] = 1.0;
      }

      art::multiplyVectorData<double>(out_c_f_p, data_c_p_d, in_c_f_p_d);
      art::multiplyVectorData<double>(out_c_f_p_d, data_c_p_d, in_c_f_p_d_d);
      if (rst::vectorNorm(&out_c_f_p[0], out_c_f_p.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyVectorData (1): check contraction of orthogonal vector components\n\n";
        errorFlag = -1000;
      }
      if (rst::vectorNorm(&out_c_f_p_d[0], out_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyVectorData (2): check contraction of orthogonal vector components\n\n";
        errorFlag = -1000;
      }
      art::multiplyVectorData<double>(out_c_f_p, data_c_1_d, in_c_f_p_d);
      art::multiplyVectorData<double>(out_c_f_p_d, data_c_1_d, in_c_f_p_d_d);
      if (rst::vectorNorm(&out_c_f_p[0], out_c_f_p.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyVectorData (3): check contraction of orthogonal vector components\n\n";
        errorFlag = -1000;
      }
      if (rst::vectorNorm(&out_c_f_p_d[0], out_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyVectorData (4): check contraction of orthogonal vector components\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking multiplyVectorData (branch 2) ************\n";

      int c=5, p=9, f=7, d1=6, d2=13;

      FieldContainer<double> in_f_p_d(f, p, d1);
      FieldContainer<double> in_f_p_d_d(f, p, d1, d2);
      FieldContainer<double> data_c_p_d(c, p, d1);
      FieldContainer<double> data_c_1_d(c, 1, d1);
      FieldContainer<double> out_c_f_p(c, f, p);
      FieldContainer<double> out_c_f_p_d(c, f, p, d2);
      double zero = INTREPID_TOL*10000.0;

      // fill with alternating 1's and -1's
      for (int i=0; i<in_f_p_d.size(); i++) {
        in_f_p_d[i] = std::pow((double)-1.0, (int)i);
      }
      for (int i=0; i<in_f_p_d_d.size(); i++) {
        in_f_p_d_d[i] = std::pow((double)-1.0, (int)i);
      }
      // fill with 1's
      for (int i=0; i<data_c_p_d.size(); i++) {
        data_c_p_d[i] = 1.0;
      }
      // fill with 1's
      for (int i=0; i<data_c_1_d.size(); i++) {
        data_c_1_d[i] = 1.0;
      }

      art::multiplyVectorData<double>(out_c_f_p, data_c_p_d, in_f_p_d);
      art::multiplyVectorData<double>(out_c_f_p_d, data_c_p_d, in_f_p_d_d);
      if (rst::vectorNorm(&out_c_f_p[0], out_c_f_p.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyVectorData (1): check contraction of orthogonal vector components\n\n";
        errorFlag = -1000;
      }
      if (rst::vectorNorm(&out_c_f_p_d[0], out_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyVectorData (2): check contraction of orthogonal vector components\n\n";
        errorFlag = -1000;
      }
      art::multiplyVectorData<double>(out_c_f_p, data_c_1_d, in_f_p_d);
      art::multiplyVectorData<double>(out_c_f_p_d, data_c_1_d, in_f_p_d_d);
      if (rst::vectorNorm(&out_c_f_p[0], out_c_f_p.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyVectorData (3): check contraction of orthogonal vector components\n\n";
        errorFlag = -1000;
      }
      if (rst::vectorNorm(&out_c_f_p_d[0], out_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyVectorData (4): check contraction of orthogonal vector components\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking multiplyTensorData (branch 1) ************\n";

      // d1 refers to spatial dimension and should be 1, 2 or 3
      // if d1>3, the RealSpaceTools function 'inverse' will fail
      int c=5, p=9, f=7, d1=3;

      FieldContainer<double> in_c_f_p_d(c, f, p, d1);
      FieldContainer<double> in_c_f_p_d_d(c, f, p, d1, d1);
      FieldContainer<double> data_c_p(c, p);
      FieldContainer<double> datainv_c_p(c, p);
      FieldContainer<double> data_c_1(c, 1);
      FieldContainer<double> datainv_c_1(c, 1);
      FieldContainer<double> data_c_p_d(c, p, d1);
      FieldContainer<double> datainv_c_p_d(c, p, d1);
      FieldContainer<double> data_c_1_d(c, 1, d1);
      FieldContainer<double> datainv_c_1_d(c, 1, d1);
      FieldContainer<double> data_c_p_d_d(c, p, d1, d1);
      FieldContainer<double> datainv_c_p_d_d(c, p, d1, d1);
      FieldContainer<double> datainvtrn_c_p_d_d(c, p, d1, d1);
      FieldContainer<double> data_c_1_d_d(c, 1, d1, d1);
      FieldContainer<double> datainv_c_1_d_d(c, 1, d1, d1);
      FieldContainer<double> datainvtrn_c_1_d_d(c, 1, d1, d1);
      FieldContainer<double> out_c_f_p_d(c, f, p, d1);
      FieldContainer<double> outi_c_f_p_d(c, f, p, d1);
      FieldContainer<double> out_c_f_p_d_d(c, f, p, d1, d1);
      FieldContainer<double> outi_c_f_p_d_d(c, f, p, d1, d1);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
      for (int i=0; i<in_c_f_p_d.size(); i++) {
        in_c_f_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_c_f_p_d_d.size(); i++) {
        in_c_f_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<data_c_p.size(); i++) {
        data_c_p[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_p[i] = 1.0 / data_c_p[i];
      }
      for (int i=0; i<data_c_1.size(); i++) {
        data_c_1[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_1[i] = 1.0 / data_c_1[i];
      }

      art::multiplyTensorData<double>(out_c_f_p_d, data_c_p, in_c_f_p_d);
      art::multiplyTensorData<double>(out_c_f_p_d, datainv_c_p, out_c_f_p_d);
      rst::subtract(&out_c_f_p_d[0], &in_c_f_p_d[0], out_c_f_p_d.size());
      if (rst::vectorNorm(&out_c_f_p_d[0], out_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_p, in_c_f_p_d_d);
      art::multiplyTensorData<double>(out_c_f_p_d_d, datainv_c_p, out_c_f_p_d_d);
      rst::subtract(&out_c_f_p_d_d[0], &in_c_f_p_d_d[0], out_c_f_p_d_d.size());
      if (rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(out_c_f_p_d, data_c_1, in_c_f_p_d);
      art::multiplyTensorData<double>(out_c_f_p_d, datainv_c_1, out_c_f_p_d);
      rst::subtract(&out_c_f_p_d[0], &in_c_f_p_d[0], out_c_f_p_d.size());
      if (rst::vectorNorm(&out_c_f_p_d[0], out_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_1, in_c_f_p_d_d);
      art::multiplyTensorData<double>(out_c_f_p_d_d, datainv_c_1, out_c_f_p_d_d);
      rst::subtract(&out_c_f_p_d_d[0], &in_c_f_p_d_d[0], out_c_f_p_d_d.size());
      if (rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with random numbers
      for (int i=0; i<in_c_f_p_d.size(); i++) {
        in_c_f_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_c_f_p_d_d.size(); i++) {
        in_c_f_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<data_c_p_d.size(); i++) {
        data_c_p_d[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_p_d[i] = 1.0 / data_c_p_d[i];
      }
      for (int i=0; i<data_c_1_d.size(); i++) {
        data_c_1_d[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_1_d[i] = 1.0 / data_c_1_d[i];
      }

      art::multiplyTensorData<double>(out_c_f_p_d, data_c_p_d, in_c_f_p_d);
      art::multiplyTensorData<double>(out_c_f_p_d, datainv_c_p_d, out_c_f_p_d);
      rst::subtract(&out_c_f_p_d[0], &in_c_f_p_d[0], out_c_f_p_d.size());
      if (rst::vectorNorm(&out_c_f_p_d[0], out_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (5): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_p_d, in_c_f_p_d_d);
      art::multiplyTensorData<double>(out_c_f_p_d_d, datainv_c_p_d, out_c_f_p_d_d);
      rst::subtract(&out_c_f_p_d_d[0], &in_c_f_p_d_d[0], out_c_f_p_d_d.size());
      if (rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (6): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(out_c_f_p_d, data_c_1_d, in_c_f_p_d);
      art::multiplyTensorData<double>(out_c_f_p_d, datainv_c_1_d, out_c_f_p_d);
      rst::subtract(&out_c_f_p_d[0], &in_c_f_p_d[0], out_c_f_p_d.size());
      if (rst::vectorNorm(&out_c_f_p_d[0], out_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (7): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_1_d, in_c_f_p_d_d);
      art::multiplyTensorData<double>(out_c_f_p_d_d, datainv_c_1_d, out_c_f_p_d_d);
      rst::subtract(&out_c_f_p_d_d[0], &in_c_f_p_d_d[0], out_c_f_p_d_d.size());
      if (rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (8): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with random numbers
      for (int i=0; i<in_c_f_p_d.size(); i++) {
        in_c_f_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_c_f_p_d_d.size(); i++) {
        in_c_f_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int ic=0; ic < c; ic++) {
        for (int ip=0; ip < p; ip++) {
          for (int i=0; i<d1*d1; i++) {
            (&data_c_p_d_d(ic, ip, 0, 0))[i] = Teuchos::ScalarTraits<double>::random();
          }
          RealSpaceTools<double>::inverse(&datainv_c_p_d_d(ic, ip, 0, 0), &data_c_p_d_d(ic, ip, 0, 0), d1);
        }
      }
      for (int ic=0; ic < c; ic++) {
        for (int ip=0; ip < 1; ip++) {
          for (int i=0; i<d1*d1; i++) {
            (&data_c_1_d_d(ic, 0, 0, 0))[i] = Teuchos::ScalarTraits<double>::random();
          }
          RealSpaceTools<double>::inverse(&datainv_c_1_d_d(ic, 0, 0, 0), &data_c_1_d_d(ic, 0, 0, 0), d1);
        }
      }

      art::multiplyTensorData<double>(out_c_f_p_d, data_c_p_d_d, in_c_f_p_d);
      art::multiplyTensorData<double>(outi_c_f_p_d, datainv_c_p_d_d, out_c_f_p_d);
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (9): check matrix inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_p_d_d, in_c_f_p_d_d);
      art::multiplyTensorData<double>(outi_c_f_p_d_d, datainv_c_p_d_d, out_c_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (10): check matrix inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(out_c_f_p_d, data_c_1_d_d, in_c_f_p_d);
      art::multiplyTensorData<double>(outi_c_f_p_d, datainv_c_1_d_d, out_c_f_p_d);
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (11): check matrix inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_1_d_d, in_c_f_p_d_d);
      art::multiplyTensorData<double>(outi_c_f_p_d_d, datainv_c_1_d_d, out_c_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (12): check matrix inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(out_c_f_p_d, data_c_p_d_d, in_c_f_p_d, 't');
      art::multiplyTensorData<double>(outi_c_f_p_d, datainv_c_p_d_d, out_c_f_p_d, 't');
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (13): check matrix inverse property, w/ double transpose\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_p_d_d, in_c_f_p_d_d, 't');
      art::multiplyTensorData<double>(outi_c_f_p_d_d, datainv_c_p_d_d, out_c_f_p_d_d, 't');
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (14): check matrix inverse property, w/ double transpose\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(out_c_f_p_d, data_c_1_d_d, in_c_f_p_d, 't');
      art::multiplyTensorData<double>(outi_c_f_p_d, datainv_c_1_d_d, out_c_f_p_d, 't');
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (15): check matrix inverse property, w/ double transpose\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_1_d_d, in_c_f_p_d_d, 't');
      art::multiplyTensorData<double>(outi_c_f_p_d_d, datainv_c_1_d_d, out_c_f_p_d_d, 't');
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (16): check matrix inverse property, w/ double transpose\n\n";
        errorFlag = -1000;
      }

      // fill with random numbers
      for (int i=0; i<in_c_f_p_d.size(); i++) {
        in_c_f_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_c_f_p_d_d.size(); i++) {
        in_c_f_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int ic=0; ic < c; ic++) {
        for (int ip=0; ip < p; ip++) {
          for (int i=0; i<d1*d1; i++) {
            (&data_c_p_d_d(ic, ip, 0, 0))[i] = Teuchos::ScalarTraits<double>::random();
          }
          RealSpaceTools<double>::inverse(&datainv_c_p_d_d(ic, ip, 0, 0), &data_c_p_d_d(ic, ip, 0, 0), d1);
          RealSpaceTools<double>::transpose(&datainvtrn_c_p_d_d(ic, ip, 0, 0), &datainv_c_p_d_d(ic, ip, 0, 0), d1);
        }
      }
      for (int ic=0; ic < c; ic++) {
        for (int ip=0; ip < 1; ip++) {
          for (int i=0; i<d1*d1; i++) {
            (&data_c_1_d_d(ic, 0, 0, 0))[i] = Teuchos::ScalarTraits<double>::random();
          }
          RealSpaceTools<double>::inverse(&datainv_c_1_d_d(ic, 0, 0, 0), &data_c_1_d_d(ic, 0, 0, 0), d1);
          RealSpaceTools<double>::transpose(&datainvtrn_c_1_d_d(ic, 0, 0, 0), &datainv_c_1_d_d(ic, 0, 0, 0), d1);
        }
      }

      art::multiplyTensorData<double>(out_c_f_p_d, data_c_p_d_d, in_c_f_p_d);
      art::multiplyTensorData<double>(outi_c_f_p_d, datainvtrn_c_p_d_d, out_c_f_p_d, 't');
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (17): check matrix inverse transpose property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_p_d_d, in_c_f_p_d_d);
      art::multiplyTensorData<double>(outi_c_f_p_d_d, datainvtrn_c_p_d_d, out_c_f_p_d_d, 't');
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (18): check matrix inverse transpose property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(out_c_f_p_d, data_c_1_d_d, in_c_f_p_d);
      art::multiplyTensorData<double>(outi_c_f_p_d, datainvtrn_c_1_d_d, out_c_f_p_d, 't');
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (19): check matrix inverse transpose property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_1_d_d, in_c_f_p_d_d);
      art::multiplyTensorData<double>(outi_c_f_p_d_d, datainvtrn_c_1_d_d, out_c_f_p_d_d, 't');
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (20): check matrix inverse transpose property\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking multiplyTensorData (branch 2) ************\n";

      // d1 refers to spatial dimension and should be 1, 2 or 3
      // if d1>3, the RealSpaceTools function 'inverse' will fail
      int c=5, p=9, f=7, d1=3;

      FieldContainer<double> in_f_p_d(f, p, d1);
      FieldContainer<double> in_f_p_d_d(f, p, d1, d1);
      FieldContainer<double> in_c_f_p_d(c, f, p, d1);
      FieldContainer<double> in_c_f_p_d_d(c, f, p, d1, d1);
      FieldContainer<double> data_c_p(c, p);
      FieldContainer<double> datainv_c_p(c, p);
      FieldContainer<double> data_c_1(c, 1);
      FieldContainer<double> datainv_c_1(c, 1);
      FieldContainer<double> data_c_p_d(c, p, d1);
      FieldContainer<double> datainv_c_p_d(c, p, d1);
      FieldContainer<double> data_c_1_d(c, 1, d1);
      FieldContainer<double> datainv_c_1_d(c, 1, d1);
      FieldContainer<double> data_c_p_d_d(c, p, d1, d1);
      FieldContainer<double> datainv_c_p_d_d(c, p, d1, d1);
      FieldContainer<double> datainvtrn_c_p_d_d(c, p, d1, d1);
      FieldContainer<double> data_c_1_d_d(c, 1, d1, d1);
      FieldContainer<double> datainv_c_1_d_d(c, 1, d1, d1);
      FieldContainer<double> datainvtrn_c_1_d_d(c, 1, d1, d1);
      FieldContainer<double> data_c_p_one(c, p);
      FieldContainer<double> data_c_1_one(c, 1);
      FieldContainer<double> out_c_f_p_d(c, f, p, d1);
      FieldContainer<double> outi_c_f_p_d(c, f, p, d1);
      FieldContainer<double> out_c_f_p_d_d(c, f, p, d1, d1);
      FieldContainer<double> outi_c_f_p_d_d(c, f, p, d1, d1);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
      for (int i=0; i<in_f_p_d.size(); i++) {
        in_f_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_f_p_d_d.size(); i++) {
        in_f_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<data_c_p.size(); i++) {
        data_c_p[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_p[i] = 1.0 / data_c_p[i];
        data_c_p_one[i] = 1.0;
      }
      for (int i=0; i<data_c_1.size(); i++) {
        data_c_1[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_1[i] = 1.0 / data_c_1[i];
      }

      art::multiplyTensorData<double>(in_c_f_p_d, data_c_p_one, in_f_p_d);
      art::multiplyTensorData<double>(out_c_f_p_d, data_c_p, in_f_p_d);
      art::multiplyTensorData<double>(outi_c_f_p_d, datainv_c_p, out_c_f_p_d);
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_p, in_f_p_d_d);
      art::multiplyTensorData<double>(outi_c_f_p_d_d, datainv_c_p, out_c_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(in_c_f_p_d, data_c_p_one, in_f_p_d);
      art::multiplyTensorData<double>(out_c_f_p_d, data_c_1, in_f_p_d);
      art::multiplyTensorData<double>(outi_c_f_p_d, datainv_c_1, out_c_f_p_d);
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_1, in_f_p_d_d);
      art::multiplyTensorData<double>(outi_c_f_p_d_d, datainv_c_1, out_c_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with random numbers
      for (int i=0; i<in_f_p_d.size(); i++) {
        in_f_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_f_p_d_d.size(); i++) {
        in_f_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<data_c_p_d.size(); i++) {
        data_c_p_d[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_p_d[i] = 1.0 / data_c_p_d[i];
      }
      for (int i=0; i<data_c_1_d.size(); i++) {
        data_c_1_d[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_1_d[i] = 1.0 / data_c_1_d[i];
      }

      art::multiplyTensorData<double>(in_c_f_p_d, data_c_p_one, in_f_p_d);
      art::multiplyTensorData<double>(out_c_f_p_d, data_c_p_d, in_f_p_d);
      art::multiplyTensorData<double>(outi_c_f_p_d, datainv_c_p_d, out_c_f_p_d);
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (5): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_p_d, in_f_p_d_d);
      art::multiplyTensorData<double>(outi_c_f_p_d_d, datainv_c_p_d, out_c_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (6): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(in_c_f_p_d, data_c_p_one, in_f_p_d);
      art::multiplyTensorData<double>(out_c_f_p_d, data_c_1_d, in_f_p_d);
      art::multiplyTensorData<double>(outi_c_f_p_d, datainv_c_1_d, out_c_f_p_d);
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (7): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_1_d, in_f_p_d_d);
      art::multiplyTensorData<double>(outi_c_f_p_d_d, datainv_c_1_d, out_c_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (8): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with random numbers
      for (int i=0; i<in_f_p_d.size(); i++) {
        in_f_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_f_p_d_d.size(); i++) {
        in_f_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int ic=0; ic < c; ic++) {
        for (int ip=0; ip < p; ip++) {
          for (int i=0; i<d1*d1; i++) {
            (&data_c_p_d_d(ic, ip, 0, 0))[i] = Teuchos::ScalarTraits<double>::random();
          }
          RealSpaceTools<double>::inverse(&datainv_c_p_d_d(ic, ip, 0, 0), &data_c_p_d_d(ic, ip, 0, 0), d1);
        }
      }
      for (int ic=0; ic < c; ic++) {
        for (int ip=0; ip < 1; ip++) {
          for (int i=0; i<d1*d1; i++) {
            (&data_c_1_d_d(ic, 0, 0, 0))[i] = Teuchos::ScalarTraits<double>::random();
          }
          RealSpaceTools<double>::inverse(&datainv_c_1_d_d(ic, 0, 0, 0), &data_c_1_d_d(ic, 0, 0, 0), d1);
        }
      }

      art::multiplyTensorData<double>(in_c_f_p_d, data_c_p_one, in_f_p_d);
      art::multiplyTensorData<double>(out_c_f_p_d, data_c_p_d_d, in_f_p_d);
      art::multiplyTensorData<double>(outi_c_f_p_d, datainv_c_p_d_d, out_c_f_p_d);
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (9): check matrix inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_p_d_d, in_f_p_d_d);
      art::multiplyTensorData<double>(outi_c_f_p_d_d, datainv_c_p_d_d, out_c_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (10): check matrix inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(in_c_f_p_d, data_c_p_one, in_f_p_d);
      art::multiplyTensorData<double>(out_c_f_p_d, data_c_1_d_d, in_f_p_d);
      art::multiplyTensorData<double>(outi_c_f_p_d, datainv_c_1_d_d, out_c_f_p_d);
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (11): check matrix inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_1_d_d, in_f_p_d_d);
      art::multiplyTensorData<double>(outi_c_f_p_d_d, datainv_c_1_d_d, out_c_f_p_d_d);
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (12): check matrix inverse property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(in_c_f_p_d, data_c_p_one, in_f_p_d);
      art::multiplyTensorData<double>(out_c_f_p_d, data_c_p_d_d, in_f_p_d, 't');
      art::multiplyTensorData<double>(outi_c_f_p_d, datainv_c_p_d_d, out_c_f_p_d, 't');
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (13): check matrix inverse property, w/ double transpose\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_p_d_d, in_f_p_d_d, 't');
      art::multiplyTensorData<double>(outi_c_f_p_d_d, datainv_c_p_d_d, out_c_f_p_d_d, 't');
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (14): check matrix inverse property, w/ double transpose\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(in_c_f_p_d, data_c_p_one, in_f_p_d);
      art::multiplyTensorData<double>(out_c_f_p_d, data_c_1_d_d, in_f_p_d, 't');
      art::multiplyTensorData<double>(outi_c_f_p_d, datainv_c_1_d_d, out_c_f_p_d, 't');
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (15): check matrix inverse property, w/ double transpose\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_1_d_d, in_f_p_d_d, 't');
      art::multiplyTensorData<double>(outi_c_f_p_d_d, datainv_c_1_d_d, out_c_f_p_d_d, 't');
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (16): check matrix inverse property, w/ double transpose\n\n";
        errorFlag = -1000;
      }

      // fill with random numbers
      for (int i=0; i<in_f_p_d.size(); i++) {
        in_f_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_f_p_d_d.size(); i++) {
        in_f_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int ic=0; ic < c; ic++) {
        for (int ip=0; ip < p; ip++) {
          for (int i=0; i<d1*d1; i++) {
            (&data_c_p_d_d(ic, ip, 0, 0))[i] = Teuchos::ScalarTraits<double>::random();
          }
          RealSpaceTools<double>::inverse(&datainv_c_p_d_d(ic, ip, 0, 0), &data_c_p_d_d(ic, ip, 0, 0), d1);
          RealSpaceTools<double>::transpose(&datainvtrn_c_p_d_d(ic, ip, 0, 0), &datainv_c_p_d_d(ic, ip, 0, 0), d1);
        }
      }
      for (int ic=0; ic < c; ic++) {
        for (int ip=0; ip < 1; ip++) {
          for (int i=0; i<d1*d1; i++) {
            (&data_c_1_d_d(ic, 0, 0, 0))[i] = Teuchos::ScalarTraits<double>::random();
          }
          RealSpaceTools<double>::inverse(&datainv_c_1_d_d(ic, 0, 0, 0), &data_c_1_d_d(ic, 0, 0, 0), d1);
          RealSpaceTools<double>::transpose(&datainvtrn_c_1_d_d(ic, 0, 0, 0), &datainv_c_1_d_d(ic, 0, 0, 0), d1);
        }
      }

      art::multiplyTensorData<double>(in_c_f_p_d, data_c_p_one, in_f_p_d);
      art::multiplyTensorData<double>(out_c_f_p_d, data_c_p_d_d, in_f_p_d);
      art::multiplyTensorData<double>(outi_c_f_p_d, datainvtrn_c_p_d_d, out_c_f_p_d, 't');
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (17): check matrix inverse transpose property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_p_d_d, in_f_p_d_d);
      art::multiplyTensorData<double>(outi_c_f_p_d_d, datainvtrn_c_p_d_d, out_c_f_p_d_d, 't');
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (18): check matrix inverse transpose property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(in_c_f_p_d, data_c_p_one, in_f_p_d);
      art::multiplyTensorData<double>(out_c_f_p_d, data_c_1_d_d, in_f_p_d);
      art::multiplyTensorData<double>(outi_c_f_p_d, datainvtrn_c_1_d_d, out_c_f_p_d, 't');
      rst::subtract(&outi_c_f_p_d[0], &in_c_f_p_d[0], outi_c_f_p_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d[0], outi_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (19): check matrix inverse transpose property\n\n";
        errorFlag = -1000;
      }
      art::multiplyTensorData<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      art::multiplyTensorData<double>(out_c_f_p_d_d, data_c_1_d_d, in_f_p_d_d);
      art::multiplyTensorData<double>(outi_c_f_p_d_d, datainvtrn_c_1_d_d, out_c_f_p_d_d, 't');
      rst::subtract(&outi_c_f_p_d_d[0], &in_c_f_p_d_d[0], outi_c_f_p_d_d.size());
      if (rst::vectorNorm(&outi_c_f_p_d_d[0], outi_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT multiplyTensorData (20): check matrix inverse transpose property\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking cloneValues ************\n";

      int c=5, p=9, f=7, d1=7, d2=13;

      FieldContainer<double> in_f_p(f, p);
      FieldContainer<double> in_f_p_d(f, p, d1);
      FieldContainer<double> in_f_p_d_d(f, p, d1, d2);
      FieldContainer<double> in_c_f_p(c, f, p);
      FieldContainer<double> in_c_f_p_d(c, f, p, d1);
      FieldContainer<double> in_c_f_p_d_d(c, f, p, d1, d2);
      FieldContainer<double> data_c_p_one(c, p);
      FieldContainer<double> out_c_f_p(c, f, p);
      FieldContainer<double> out_c_f_p_d(c, f, p, d1);
      FieldContainer<double> out_c_f_p_d_d(c, f, p, d1, d2);
      double zero = INTREPID_TOL*100.0;

      // fill with random numbers
      for (int i=0; i<in_f_p.size(); i++) {
        in_f_p[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_f_p_d.size(); i++) {
        in_f_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<in_f_p_d_d.size(); i++) {
        in_f_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<data_c_p_one.size(); i++) {
        data_c_p_one[i] = 1.0;
      }

      art::cloneValues<double>(out_c_f_p, in_f_p);
      art::multiplyScalarData<double>(in_c_f_p, data_c_p_one, in_f_p);
      rst::subtract(&out_c_f_p[0], &in_c_f_p[0], out_c_f_p.size());
      if (rst::vectorNorm(&out_c_f_p[0], out_c_f_p.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT cloneValues (1): check multiplyScalarData vs. cloneValues\n\n";
        errorFlag = -1000;
      }
      art::cloneValues<double>(out_c_f_p_d, in_f_p_d);
      art::multiplyScalarData<double>(in_c_f_p_d, data_c_p_one, in_f_p_d);
      rst::subtract(&out_c_f_p_d[0], &in_c_f_p_d[0], out_c_f_p_d.size());
      if (rst::vectorNorm(&out_c_f_p_d[0], out_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT cloneValues (2): check multiplyScalarData vs. cloneValues\n\n";
        errorFlag = -1000;
      }
      art::cloneValues<double>(out_c_f_p_d_d, in_f_p_d_d);
      art::multiplyScalarData<double>(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      rst::subtract(&out_c_f_p_d_d[0], &in_c_f_p_d_d[0], out_c_f_p_d_d.size());
      if (rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT cloneValues (3): check multiplyScalarData vs. cloneValues\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking cloneScaleValues ************\n";
      int c=5, p=9, f=7, d1=7, d2=13;

      FieldContainer<double> in_f_p(f, p);
      FieldContainer<double> in_f_p_d(f, p, d1);
      FieldContainer<double> in_f_p_d_d(f, p, d1, d2);
      FieldContainer<double> data_c_f(c, f);
      FieldContainer<double> datainv_c_f(c, f);
      FieldContainer<double> c_f_p_one(c, f, p);
      FieldContainer<double> c_f_p_d_one(c, f, p, d1);
      FieldContainer<double> c_f_p_d_d_one(c, f, p, d1, d2);
      FieldContainer<double> out_c_f_p(c, f, p);
      FieldContainer<double> outi_c_f_p(c, f, p);
      FieldContainer<double> out_c_f_p_d(c, f, p, d1);
      FieldContainer<double> outi_c_f_p_d(c, f, p, d1);
      FieldContainer<double> out_c_f_p_d_d(c, f, p, d1, d2);
      FieldContainer<double> outi_c_f_p_d_d(c, f, p, d1, d2);
      double zero = INTREPID_TOL*100.0;

      // fill with 1's
      for (int i=0; i<in_f_p.size(); i++) {
        in_f_p[i] = 1.0;
      }
      for (int i=0; i<in_f_p_d.size(); i++) {
        in_f_p_d[i] = 1.0;
      }
      for (int i=0; i<in_f_p_d_d.size(); i++) {
        in_f_p_d_d[i] = 1.0;
      }
      for (int i=0; i<c_f_p_one.size(); i++) {
        c_f_p_one[i] = 1.0;
      }
      for (int i=0; i<c_f_p_d_one.size(); i++) {
        c_f_p_d_one[i] = 1.0;
      }
      for (int i=0; i<c_f_p_d_d_one.size(); i++) {
        c_f_p_d_d_one[i] = 1.0;
      }
      // fill with random numbers
      for (int i=0; i<data_c_f.size(); i++) {
        data_c_f[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_f[i] = 1.0 / data_c_f[i];
      }

      art::cloneScaleValues<double>(out_c_f_p, data_c_f, in_f_p);
      art::cloneScaleValues<double>(outi_c_f_p, datainv_c_f, in_f_p);
      for (int i=0; i<out_c_f_p.size(); i++) {
        out_c_f_p[i] *= outi_c_f_p[i];
      }
      rst::subtract(&out_c_f_p[0], &c_f_p_one[0], out_c_f_p.size());
      if (rst::vectorNorm(&out_c_f_p[0], out_c_f_p.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT cloneScaleValue (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      art::cloneScaleValues<double>(out_c_f_p_d, data_c_f, in_f_p_d);
      art::cloneScaleValues<double>(outi_c_f_p_d, datainv_c_f, in_f_p_d);
      for (int i=0; i<out_c_f_p_d.size(); i++) {
        out_c_f_p_d[i] *= outi_c_f_p_d[i];
      }
      rst::subtract(&out_c_f_p_d[0], &c_f_p_d_one[0], out_c_f_p_d.size());
      if (rst::vectorNorm(&out_c_f_p_d[0], out_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT cloneScaleValue (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      art::cloneScaleValues<double>(out_c_f_p_d_d, data_c_f, in_f_p_d_d);
      art::cloneScaleValues<double>(outi_c_f_p_d_d, datainv_c_f, in_f_p_d_d);
      for (int i=0; i<out_c_f_p_d_d.size(); i++) {
        out_c_f_p_d_d[i] *= outi_c_f_p_d_d[i];
      }
      rst::subtract(&out_c_f_p_d_d[0], &c_f_p_d_d_one[0], out_c_f_p_d_d.size());
      if (rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT cloneScaleValue (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking scaleValues ************\n";
      int c=5, p=9, f=7, d1=7, d2=13;

      FieldContainer<double> data_c_f(c, f);
      FieldContainer<double> datainv_c_f(c, f);
      FieldContainer<double> out_c_f_p(c, f, p);
      FieldContainer<double> outi_c_f_p(c, f, p);
      FieldContainer<double> out_c_f_p_d(c, f, p, d1);
      FieldContainer<double> outi_c_f_p_d(c, f, p, d1);
      FieldContainer<double> out_c_f_p_d_d(c, f, p, d1, d2);
      FieldContainer<double> outi_c_f_p_d_d(c, f, p, d1, d2);
      double zero = INTREPID_TOL*100.0;

      // fill with random numbers
      for (int i=0; i<out_c_f_p.size(); i++) {
        out_c_f_p[i] = Teuchos::ScalarTraits<double>::random();
        outi_c_f_p[i] = out_c_f_p[i];
      }
      for (int i=0; i<out_c_f_p_d.size(); i++) {
        out_c_f_p_d[i] = Teuchos::ScalarTraits<double>::random();
        outi_c_f_p_d[i] = out_c_f_p_d[i];
      }
      for (int i=0; i<out_c_f_p_d_d.size(); i++) {
        out_c_f_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
        outi_c_f_p_d_d[i] = out_c_f_p_d_d[i];
      }
      for (int i=0; i<data_c_f.size(); i++) {
        data_c_f[i] = Teuchos::ScalarTraits<double>::random();
        datainv_c_f[i] = 1.0 / data_c_f[i];
      }

      art::scaleValues<double>(out_c_f_p, data_c_f);
      art::scaleValues<double>(out_c_f_p, datainv_c_f);
      rst::subtract(&out_c_f_p[0], &outi_c_f_p[0], out_c_f_p.size());
      if (rst::vectorNorm(&out_c_f_p[0], out_c_f_p.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scaleValue (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      art::scaleValues<double>(out_c_f_p_d, data_c_f);
      art::scaleValues<double>(out_c_f_p_d, datainv_c_f);
      rst::subtract(&out_c_f_p_d[0], &outi_c_f_p_d[0], out_c_f_p_d.size());
      if (rst::vectorNorm(&out_c_f_p_d[0], out_c_f_p_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT scaleValue (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      art::scaleValues<double>(out_c_f_p_d_d, data_c_f);
      art::scaleValues<double>(out_c_f_p_d_d, datainv_c_f);
      rst::subtract(&out_c_f_p_d_d[0], &outi_c_f_p_d_d[0], out_c_f_p_d_d.size());
      if (rst::vectorNorm(&out_c_f_p_d_d[0], out_c_f_p_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT cloneScaleValue (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      } // end scope

      /******************************************/
      *outStream << "\n";
  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
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
