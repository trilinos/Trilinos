// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
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
#include "Teuchos_GlobalMPISession.hpp"

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

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

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
  << "|     1) Array operations: contractions                                       |\n" \
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
  int endThrowNumber = beginThrowNumber + 81;
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

    *outStream << "-> contractFieldFieldScalar:\n";
    INTREPID_TEST_COMMAND( atools.contractFieldFieldScalar<double>(a_2_2, a_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldScalar<double>(a_2_2, a_10_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldScalar<double>(a_2_2, a_10_2_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldScalar<double>(a_10_2_2, a_9_2_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldScalar<double>(a_10_2_2, a_10_2_2, a_10_2_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldScalar<double>(a_9_2_2, a_10_3_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldScalar<double>(a_10_3_2, a_10_2_2, a_10_3_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldScalar<double>(a_10_2_2, a_10_2_2, a_10_3_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldScalar<double>(a_10_2_3, a_10_2_2, a_10_3_2, COMP_ENGINE_MAX) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldScalar<double>(a_10_2_3, a_10_2_2, a_10_3_2, COMP_CPP) );

    FieldContainer<double> a_10_2_2_2(10, 2, 2, 2);
    FieldContainer<double> a_9_2_2_2(9, 2, 2, 2);
    FieldContainer<double> a_10_3_2_2(10, 3, 2, 2);
    FieldContainer<double> a_10_2_3_2(10, 2, 3, 2);
    FieldContainer<double> a_10_2_2_3(10, 2, 2, 3);

    *outStream << "-> contractFieldFieldVector:\n";
    INTREPID_TEST_COMMAND( atools.contractFieldFieldVector<double>(a_2_2, a_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldVector<double>(a_2_2, a_10_2_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldVector<double>(a_2_2, a_10_2_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldVector<double>(a_10_2_2, a_9_2_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldVector<double>(a_10_2_2, a_10_2_2_2, a_10_2_3_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldVector<double>(a_10_2_2, a_10_2_2_2, a_10_2_2_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldVector<double>(a_9_2_2, a_10_2_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldVector<double>(a_10_3_2, a_10_2_2_2, a_10_3_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldVector<double>(a_10_2_2, a_10_2_2_2, a_10_3_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldVector<double>(a_10_2_3, a_10_2_2_2, a_10_3_2_2, COMP_ENGINE_MAX) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldVector<double>(a_10_2_3, a_10_2_2_2, a_10_3_2_2, COMP_CPP) );

    FieldContainer<double> a_10_2_2_2_2(10, 2, 2, 2, 2);
    FieldContainer<double> a_9_2_2_2_2(9, 2, 2, 2, 2);
    FieldContainer<double> a_10_3_2_2_2(10, 3, 2, 2, 2);
    FieldContainer<double> a_10_2_3_2_2(10, 2, 3, 2, 2);
    FieldContainer<double> a_10_2_2_3_2(10, 2, 2, 3, 2);
    FieldContainer<double> a_10_2_2_2_3(10, 2, 2, 2, 3);

    *outStream << "-> contractFieldFieldTensor:\n";
    INTREPID_TEST_COMMAND( atools.contractFieldFieldTensor<double>(a_2_2, a_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldTensor<double>(a_2_2, a_10_2_2_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldTensor<double>(a_2_2, a_10_2_2_2_2, a_10_2_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldTensor<double>(a_10_2_2, a_9_2_2_2_2, a_10_2_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldTensor<double>(a_10_2_2, a_10_2_2_2_2, a_10_2_3_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldTensor<double>(a_10_2_2, a_10_2_2_2_2, a_10_2_2_3_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldTensor<double>(a_10_2_2, a_10_2_2_2_2, a_10_2_2_2_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldTensor<double>(a_9_2_2, a_10_2_2_2_2, a_10_2_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldTensor<double>(a_10_3_2, a_10_2_2_2_2, a_10_2_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldTensor<double>(a_10_2_2, a_10_2_2_2_2, a_10_3_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldTensor<double>(a_10_2_3, a_10_2_2_2_2, a_10_3_2_2_2, COMP_ENGINE_MAX) );
    INTREPID_TEST_COMMAND( atools.contractFieldFieldTensor<double>(a_10_2_3, a_10_2_2_2_2, a_10_3_2_2_2, COMP_CPP) );

    FieldContainer<double> a_9_2(9, 2);
    FieldContainer<double> a_10_1(10, 1);

    *outStream << "-> contractDataFieldScalar:\n";
    INTREPID_TEST_COMMAND( atools.contractDataFieldScalar<double>(a_2_2, a_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldScalar<double>(a_2_2, a_10_2_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldScalar<double>(a_10_2_2, a_2_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldScalar<double>(a_2_2, a_10_2, a_9_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldScalar<double>(a_2_2, a_10_3, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldScalar<double>(a_9_2, a_10_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldScalar<double>(a_10_2, a_10_2, a_10_3_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldScalar<double>(a_10_2, a_10_2, a_10_2_2, COMP_ENGINE_MAX) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldScalar<double>(a_10_2, a_10_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldScalar<double>(a_10_2, a_10_1, a_10_2_2, COMP_CPP) );

    FieldContainer<double> a_10_1_2(10, 1, 2);
    FieldContainer<double> a_10_1_3(10, 1, 3);

    *outStream << "-> contractDataFieldVector:\n";
    INTREPID_TEST_COMMAND( atools.contractDataFieldVector<double>(a_2_2, a_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldVector<double>(a_2_2, a_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldVector<double>(a_10_2_2, a_10_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldVector<double>(a_2_2, a_10_2_2, a_9_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldVector<double>(a_2_2, a_10_2_2, a_10_2_3_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldVector<double>(a_2_2, a_10_2_2, a_10_2_2_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldVector<double>(a_9_2, a_10_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldVector<double>(a_10_2, a_10_2_2, a_10_3_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldVector<double>(a_10_2, a_10_2_2, a_10_2_2_2, COMP_ENGINE_MAX) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldVector<double>(a_10_2, a_10_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldVector<double>(a_10_2, a_10_1_2, a_10_2_2_2, COMP_CPP) );

    FieldContainer<double> a_10_1_2_2(10, 1, 2, 2);

    *outStream << "-> contractDataFieldTensor:\n";
    INTREPID_TEST_COMMAND( atools.contractDataFieldTensor<double>(a_2_2, a_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldTensor<double>(a_2_2, a_2_2, a_10_2_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldTensor<double>(a_10_2_2, a_10_2_2_2, a_10_2_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldTensor<double>(a_2_2, a_10_2_2_2, a_9_2_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldTensor<double>(a_2_2, a_10_2_2_2, a_10_2_3_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldTensor<double>(a_2_2, a_10_2_2_2, a_10_2_2_3_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldTensor<double>(a_2_2, a_10_2_2_2, a_10_2_2_2_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldTensor<double>(a_9_2, a_10_2_2_2, a_10_2_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldTensor<double>(a_10_2, a_10_2_2_2, a_10_3_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldTensor<double>(a_10_2, a_10_2_2_2, a_10_2_2_2_2, COMP_ENGINE_MAX) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldTensor<double>(a_10_2, a_10_2_2_2, a_10_2_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataFieldTensor<double>(a_10_2, a_10_1_2_2, a_10_2_2_2_2, COMP_CPP) );

    FieldContainer<double> a_2(2);
    FieldContainer<double> a_10(10);

    *outStream << "-> contractDataDataScalar:\n";
    INTREPID_TEST_COMMAND( atools.contractDataDataScalar<double>(a_2_2, a_10_2_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataScalar<double>(a_2_2, a_10_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataScalar<double>(a_2_2, a_10_2, a_10_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataScalar<double>(a_2, a_9_2, a_10_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataScalar<double>(a_2, a_10_2, a_10_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataScalar<double>(a_2, a_10_2, a_10_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataScalar<double>(a_10, a_10_2, a_10_2, COMP_ENGINE_MAX) );
    INTREPID_TEST_COMMAND( atools.contractDataDataScalar<double>(a_10, a_10_2, a_10_2, COMP_CPP) );

    *outStream << "-> contractDataDataVector:\n";
    INTREPID_TEST_COMMAND( atools.contractDataDataVector<double>(a_2_2, a_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataVector<double>(a_2_2, a_10_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataVector<double>(a_2_2, a_10_2_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataVector<double>(a_10, a_9_2_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataVector<double>(a_10, a_10_3_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataVector<double>(a_10, a_10_2_3, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataVector<double>(a_2, a_10_2_2, a_10_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataVector<double>(a_10, a_10_2_2, a_10_2_2, COMP_ENGINE_MAX) );
    INTREPID_TEST_COMMAND( atools.contractDataDataVector<double>(a_10, a_10_2_2, a_10_2_2, COMP_CPP) );

    *outStream << "-> contractDataDataTensor:\n";
    INTREPID_TEST_COMMAND( atools.contractDataDataTensor<double>(a_2_2, a_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataTensor<double>(a_2_2, a_10_2_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataTensor<double>(a_2_2, a_10_2_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataTensor<double>(a_10, a_9_2_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataTensor<double>(a_10, a_10_2_2_2, a_10_3_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataTensor<double>(a_10, a_10_2_2_2, a_10_2_3_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataTensor<double>(a_10, a_10_2_2_2, a_10_2_2_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataTensor<double>(a_2, a_10_2_2_2, a_10_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( atools.contractDataDataTensor<double>(a_10, a_10_2_2_2, a_10_2_2_2, COMP_ENGINE_MAX) );
    INTREPID_TEST_COMMAND( atools.contractDataDataTensor<double>(a_10, a_10_2_2_2, a_10_2_2_2, COMP_CPP) );

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
      *outStream << "\n************ Checking contractFieldFieldScalar ************\n";

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

      art::contractFieldFieldScalar<double>(out1_c_l_r, in_c_l_p, in_c_r_p, COMP_CPP);
      art::contractFieldFieldScalar<double>(out2_c_l_r, in_c_l_p, in_c_r_p, COMP_BLAS);
      rst::subtract(&out1_c_l_r[0], &out2_c_l_r[0], out2_c_l_r.size());
      if (rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractFieldFieldScalar (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }
      // with sumInto:
      out1_c_l_r.initialize(2.0); out2_c_l_r.initialize(2.0);
      art::contractFieldFieldScalar<double>(out1_c_l_r, in_c_l_p, in_c_r_p, COMP_CPP, true);
      art::contractFieldFieldScalar<double>(out2_c_l_r, in_c_l_p, in_c_r_p, COMP_BLAS, true);
      rst::subtract(&out1_c_l_r[0], &out2_c_l_r[0], out2_c_l_r.size());
      if (rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractFieldFieldScalar (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking contractFieldFieldVector ************\n";

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

      art::contractFieldFieldVector<double>(out1_c_l_r, in_c_l_p_d, in_c_r_p_d, COMP_CPP);
      art::contractFieldFieldVector<double>(out2_c_l_r, in_c_l_p_d, in_c_r_p_d, COMP_BLAS);

      rst::subtract(&out1_c_l_r[0], &out2_c_l_r[0], out2_c_l_r.size());
      if (rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractFieldFieldVector (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }

      // with sumInto:
      out1_c_l_r.initialize(2.0); out2_c_l_r.initialize(2.0);

      art::contractFieldFieldVector<double>(out1_c_l_r, in_c_l_p_d, in_c_r_p_d, COMP_CPP, true);
      art::contractFieldFieldVector<double>(out2_c_l_r, in_c_l_p_d, in_c_r_p_d, COMP_BLAS, true);

      rst::subtract(&out1_c_l_r[0], &out2_c_l_r[0], out2_c_l_r.size());
      if (rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractFieldFieldVector (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking contractFieldFieldTensor ************\n";

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

      art::contractFieldFieldTensor<double>(out1_c_l_r, in_c_l_p_d_d, in_c_r_p_d_d, COMP_CPP);
      art::contractFieldFieldTensor<double>(out2_c_l_r, in_c_l_p_d_d, in_c_r_p_d_d, COMP_BLAS);

      rst::subtract(&out1_c_l_r[0], &out2_c_l_r[0], out2_c_l_r.size());
      if (rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractFieldFieldTensor (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }

      // with sumInto:
      out1_c_l_r.initialize(2.0); out2_c_l_r.initialize(2.0);

      art::contractFieldFieldTensor<double>(out1_c_l_r, in_c_l_p_d_d, in_c_r_p_d_d, COMP_CPP, true);
      art::contractFieldFieldTensor<double>(out2_c_l_r, in_c_l_p_d_d, in_c_r_p_d_d, COMP_BLAS, true);

      rst::subtract(&out1_c_l_r[0], &out2_c_l_r[0], out2_c_l_r.size());
      if (rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractFieldFieldTensor (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l_r[0], out1_c_l_r.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking contractDataFieldScalar ************\n";

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
      art::contractDataFieldScalar<double>(out1_c_l, data_c_p, in_c_l_p, COMP_CPP);
      art::contractDataFieldScalar<double>(out2_c_l, data_c_p, in_c_l_p, COMP_BLAS);
      rst::subtract(&out1_c_l[0], &out2_c_l[0], out2_c_l.size());
      if (rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractDataFieldScalar (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }
      // constant data
      art::contractDataFieldScalar<double>(out1_c_l, data_c_1, in_c_l_p, COMP_CPP);
      art::contractDataFieldScalar<double>(out2_c_l, data_c_1, in_c_l_p, COMP_BLAS);
      rst::subtract(&out1_c_l[0], &out2_c_l[0], out2_c_l.size());
      if (rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractDataFieldScalar (2): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }
      // nonconstant data with sumInto
      out1_c_l.initialize(2.0); out2_c_l.initialize(2.0);
      art::contractDataFieldScalar<double>(out1_c_l, data_c_p, in_c_l_p, COMP_CPP, true);
      art::contractDataFieldScalar<double>(out2_c_l, data_c_p, in_c_l_p, COMP_BLAS, true);
      rst::subtract(&out1_c_l[0], &out2_c_l[0], out2_c_l.size());
      if (rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractDataFieldScalar (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking contractDataFieldVector ************\n";

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
      art::contractDataFieldVector<double>(out1_c_l, data_c_p_d, in_c_l_p_d, COMP_CPP);
      art::contractDataFieldVector<double>(out2_c_l, data_c_p_d, in_c_l_p_d, COMP_BLAS);
      rst::subtract(&out1_c_l[0], &out2_c_l[0], out2_c_l.size());
      if (rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractDataFieldVector (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }
      // constant data
      art::contractDataFieldVector<double>(out1_c_l, data_c_1_d, in_c_l_p_d, COMP_CPP);
      art::contractDataFieldVector<double>(out2_c_l, data_c_1_d, in_c_l_p_d, COMP_BLAS);
      rst::subtract(&out1_c_l[0], &out2_c_l[0], out2_c_l.size());
      if (rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractDataFieldVector (2): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }
      // nonconstant data with sumInto
      out1_c_l.initialize(2.0); out2_c_l.initialize(2.0);
      art::contractDataFieldVector<double>(out1_c_l, data_c_p_d, in_c_l_p_d, COMP_CPP, true);
      art::contractDataFieldVector<double>(out2_c_l, data_c_p_d, in_c_l_p_d, COMP_BLAS, true);
      rst::subtract(&out1_c_l[0], &out2_c_l[0], out2_c_l.size());
      if (rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractDataFieldVector (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking contractDataFieldTensor ************\n";

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
      art::contractDataFieldTensor<double>(out1_c_l, data_c_p_d_d, in_c_l_p_d_d, COMP_CPP);
      art::contractDataFieldTensor<double>(out2_c_l, data_c_p_d_d, in_c_l_p_d_d, COMP_BLAS);
      rst::subtract(&out1_c_l[0], &out2_c_l[0], out2_c_l.size());
      if (rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractDataFieldTensor (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }
      // constant data
      art::contractDataFieldTensor<double>(out1_c_l, data_c_1_d_d, in_c_l_p_d_d, COMP_CPP);
      art::contractDataFieldTensor<double>(out2_c_l, data_c_1_d_d, in_c_l_p_d_d, COMP_BLAS);
      rst::subtract(&out1_c_l[0], &out2_c_l[0], out2_c_l.size());
      if (rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractDataFieldTensor (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }
      // nonconstant data with sumInto
      out1_c_l.initialize(2.0); out2_c_l.initialize(2.0);
      art::contractDataFieldTensor<double>(out1_c_l, data_c_p_d_d, in_c_l_p_d_d, COMP_CPP, true);
      art::contractDataFieldTensor<double>(out2_c_l, data_c_p_d_d, in_c_l_p_d_d, COMP_BLAS, true);
      rst::subtract(&out1_c_l[0], &out2_c_l[0], out2_c_l.size());
      if (rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractDataFieldTensor (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c_l[0], out1_c_l.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking contractDataDataScalar ************\n";

      int c=5, p=9;

      FieldContainer<double> inl_c_p(c, p);
      FieldContainer<double> inr_c_p(c, p);
      FieldContainer<double> out1_c(c);
      FieldContainer<double> out2_c(c);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
      for (int i=0; i<inl_c_p.size(); i++) {
        inl_c_p[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<inr_c_p.size(); i++) {
        inr_c_p[i] = Teuchos::ScalarTraits<double>::random();
      }

      art::contractDataDataScalar<double>(out1_c, inl_c_p, inr_c_p, COMP_CPP);
      art::contractDataDataScalar<double>(out2_c, inl_c_p, inr_c_p, COMP_BLAS);
      rst::subtract(&out1_c[0], &out2_c[0], out2_c.size());
      if (rst::vectorNorm(&out1_c[0], out1_c.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractDataDataScalar (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c[0], out1_c.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }
      // with sumInto:
      out1_c.initialize(2.0); out2_c.initialize(2.0);
      art::contractDataDataScalar<double>(out1_c, inl_c_p, inr_c_p, COMP_CPP, true);
      art::contractDataDataScalar<double>(out2_c, inl_c_p, inr_c_p, COMP_BLAS, true);
      rst::subtract(&out1_c[0], &out2_c[0], out2_c.size());
      if (rst::vectorNorm(&out1_c[0], out1_c.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractDataDataScalar (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c[0], out1_c.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking contractDataDataVector ************\n";

      int c=5, p=9, d=13;

      FieldContainer<double> inl_c_p_d(c, p, d);
      FieldContainer<double> inr_c_p_d(c, p, d);
      FieldContainer<double> out1_c(c);
      FieldContainer<double> out2_c(c);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
      for (int i=0; i<inl_c_p_d.size(); i++) {
        inl_c_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<inr_c_p_d.size(); i++) {
        inr_c_p_d[i] = Teuchos::ScalarTraits<double>::random();
      }

      art::contractDataDataVector<double>(out1_c, inl_c_p_d, inr_c_p_d, COMP_CPP);
      art::contractDataDataVector<double>(out2_c, inl_c_p_d, inr_c_p_d, COMP_BLAS);

      rst::subtract(&out1_c[0], &out2_c[0], out2_c.size());
      if (rst::vectorNorm(&out1_c[0], out1_c.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractDataDataVector (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c[0], out1_c.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }

      // with sumInto:
      out1_c.initialize(2.0); out2_c.initialize(2.0);

      art::contractDataDataVector<double>(out1_c, inl_c_p_d, inr_c_p_d, COMP_CPP, true);
      art::contractDataDataVector<double>(out2_c, inl_c_p_d, inr_c_p_d, COMP_BLAS, true);

      rst::subtract(&out1_c[0], &out2_c[0], out2_c.size());
      if (rst::vectorNorm(&out1_c[0], out1_c.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractDataDataVector (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c[0], out1_c.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking contractDataDataTensor ************\n";

      int c=5, p=9, d1=13, d2=5;

      FieldContainer<double> inl_c_p_d_d(c, p, d1, d2);
      FieldContainer<double> inr_c_p_d_d(c, p, d1, d2);
      FieldContainer<double> out1_c(c);
      FieldContainer<double> out2_c(c);
      double zero = INTREPID_TOL*10000.0;

      // fill with random numbers
      for (int i=0; i<inl_c_p_d_d.size(); i++) {
        inl_c_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<inr_c_p_d_d.size(); i++) {
        inr_c_p_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }

      art::contractDataDataTensor<double>(out1_c, inl_c_p_d_d, inr_c_p_d_d, COMP_CPP);
      art::contractDataDataTensor<double>(out2_c, inl_c_p_d_d, inr_c_p_d_d, COMP_BLAS);

      rst::subtract(&out1_c[0], &out2_c[0], out2_c.size());
      if (rst::vectorNorm(&out1_c[0], out1_c.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractDataDataTensor (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c[0], out1_c.size(), NORM_ONE) << "\n\n";
        errorFlag++;
      }

      // with sumInto:
      out1_c.initialize(2.0); out2_c.initialize(2.0);

      art::contractDataDataTensor<double>(out1_c, inl_c_p_d_d, inr_c_p_d_d, COMP_CPP, true);
      art::contractDataDataTensor<double>(out2_c, inl_c_p_d_d, inr_c_p_d_d, COMP_BLAS, true);

      rst::subtract(&out1_c[0], &out2_c[0], out2_c.size());
      if (rst::vectorNorm(&out1_c[0], out1_c.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT contractDataDataTensor (1): check COMP_CPP vs. COMP_BLAS; "
                   << " diff-1norm = " << rst::vectorNorm(&out1_c[0], out1_c.size(), NORM_ONE) << "\n\n";
        errorFlag++;
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
