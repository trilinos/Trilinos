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
\brief  Unit test for the FunctionSpaceTools class, testing H-grad.
\author Created by D. Ridzal, P. Bochev, and K. Peterson.
*/

#include "Intrepid_TensorProductSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_HGRAD_QUAD_Cn_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
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
  << "|                      Unit Test (TensorProductSpaceTools)                    |\n" \
  << "|                                                                             |\n" \
  << "|     1) exception testing                                                    |\n" \
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
  const int numTotalExceptions = 7;
  int beginThrowNumber = Teuchos::TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + numTotalExceptions;
#endif

  typedef TensorProductSpaceTools tpst; 

  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 1: exceptions                                                          |\n"\
  << "===============================================================================\n";

  try{
#ifdef HAVE_INTREPID_DEBUG
    FieldContainer<double> a_2_2(2,2);
    FieldContainer<double> a_2(2);
    FieldContainer<double> a_4(4);
    FieldContainer<double> a_4_4(4,4);
    FieldContainer<double> a_4_1_4(4,1,4);
    FieldContainer<double> a_4_2(4,2);
    Array<RCP<FieldContainer<double> > > bases_2(2);
    Array<RCP<FieldContainer<double> > > bases_1(1);
    bases_2[0] = Teuchos::rcp( &a_2_2 , false );
    bases_2[1] = Teuchos::rcp( &a_2_2 , false );
    bases_1[0] = Teuchos::rcp( &a_2_2 , false );

    // work with evaluate

    // first arg wrong shape
    INTREPID_TEST_COMMAND( tpst::evaluate<double>( a_4 , a_4_1_4 , bases_2 ) );
    // second arg wrong shape
    INTREPID_TEST_COMMAND( tpst::evaluate<double>( a_4_4 , a_4_4 , bases_2 ) );
    // bases wrong length
    INTREPID_TEST_COMMAND( tpst::evaluate<double>( a_4_4 , a_4_4, bases_1 ) );
    // internal bases wrong shape
    bases_2[1] = Teuchos::rcp( & a_2 , false );
    INTREPID_TEST_COMMAND( tpst::evaluate<double>( a_4_4 , a_4_4, bases_2 ) );
    // different number of cells in first two arguments
    bases_2[1] = Teuchos::rcp( &a_2_2 , false );
    INTREPID_TEST_COMMAND( tpst::evaluate<double>( a_2_2 , a_4_4, bases_2 ) );
    // points and bases incompatible
    INTREPID_TEST_COMMAND( tpst::evaluate<double>( a_4_2 , a_4_4 , bases_2 ) );
    // coeffs and bases incompatible
    INTREPID_TEST_COMMAND( tpst::evaluate<double>( a_4_4 , a_4_2 , bases_2 ) );

    // need to test gradients and other methods, too
    

#endif
  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

#ifdef HAVE_INTREPID_DEBUG
  if (Teuchos::TestForException_getThrowNumber() != endThrowNumber)
    errorFlag++;
#endif

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return errorFlag;
}
