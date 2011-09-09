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
\brief  Unit tests for the Intrepid::Basis_HGRAD_LINE_Cn_FEM_JACOBI class.
\author Created by P. Bochev, D. Ridzal, K. Peterson.
*/

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HGRAD_LINE_Cn_FEM_JACOBI.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace std;
using namespace Intrepid;

#define INTREPID_TEST_COMMAND( S , throwCounter, nException )                                                              \
{                                                                                                                          \
  ++nException;                                                                                                            \
  try {                                                                                                                    \
    S ;                                                                                                                    \
  }                                                                                                                        \
  catch (std::logic_error err) {                                                                                           \
      ++throwCounter;                                                                                                      \
      *outStream << "Expected Error " << nException << " -------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                                                                    \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";           \
  };                                                                                                                       \
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
    << "|               Unit Test (Basis_HGRAD_LINE_Cn_FEM_JACOBI)                    |\n" \
    << "|                                                                             |\n" \
    << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n" \
    << "|     2) Basis values for VALUE, GRAD, CURL, and Dk operators                 |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n"\
    << "| TEST 1: Basis creation, exception testing                                   |\n"\
    << "===============================================================================\n";

  
  // Define basis and error flag
  double alpha = 0.0, beta = 0.0;
  Basis_HGRAD_LINE_Cn_FEM_JACOBI<double, FieldContainer<double> > lineBasis(5, alpha, beta);
  int errorFlag = 0;

  // Initialize throw counter for exception testing
  int nException     = 0;
  int throwCounter   = 0;
  
  // Define array containing vertices of the reference Line and a few other points   
  int numIntervals = 100;
  FieldContainer<double> lineNodes(numIntervals+1, 1);
  for (int i=0; i<numIntervals+1; i++) {
    lineNodes(i,0) = -1.0+(2.0*(double)i)/(double)numIntervals;
  }

  // Generic array for the output values; needs to be properly resized depending on the operator type
  FieldContainer<double> vals;

  try{
    // Exceptions 1-5: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
    // getDofTag() to access invalid array elements thereby causing bounds check exception
    // exception #1
    INTREPID_TEST_COMMAND( lineBasis.getDofOrdinal(2,0,0), throwCounter, nException );
    // exception #2
    INTREPID_TEST_COMMAND( lineBasis.getDofOrdinal(1,1,1), throwCounter, nException );
    // exception #3
    INTREPID_TEST_COMMAND( lineBasis.getDofOrdinal(1,0,7), throwCounter, nException );
    // not an exception
    INTREPID_TEST_COMMAND( lineBasis.getDofOrdinal(1,0,5), throwCounter, nException ); --nException;
    // exception #4
    INTREPID_TEST_COMMAND( lineBasis.getDofTag(6), throwCounter, nException );
    // exception #5
    INTREPID_TEST_COMMAND( lineBasis.getDofTag(-1), throwCounter, nException );
    // not an exception
    INTREPID_TEST_COMMAND( lineBasis.getDofTag(5), throwCounter, nException ); --nException;
#ifdef HAVE_INTREPID_DEBUG
    // Exceptions 6-16 test exception handling with incorrectly dimensioned input/output arrays
    // exception #6: input points array must be of rank-2
    FieldContainer<double> badPoints1(4, 5, 3);
    INTREPID_TEST_COMMAND( lineBasis.getValues(vals, badPoints1, OPERATOR_VALUE), throwCounter, nException );

    // exception #7: dimension 1 in the input point array must equal space dimension of the cell
    FieldContainer<double> badPoints2(4, 3);
    INTREPID_TEST_COMMAND( lineBasis.getValues(vals, badPoints2, OPERATOR_VALUE), throwCounter, nException );

    // exception #8: output values must be of rank-2 for OPERATOR_VALUE
    FieldContainer<double> badVals1(4, 3, 1);
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals1, lineNodes, OPERATOR_VALUE), throwCounter, nException );

    // exception #9: output values must be of rank-3 for OPERATOR_GRAD
    FieldContainer<double> badVals2(4, 3);
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals2, lineNodes, OPERATOR_GRAD), throwCounter, nException );

    // exception #10: output values must be of rank-3 for OPERATOR_CURL
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals2, lineNodes, OPERATOR_CURL), throwCounter, nException );

    // exception #11: output values must be of rank-2 for OPERATOR_DIV
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals2, lineNodes, OPERATOR_DIV), throwCounter, nException );

    // exception #12: output values must be of rank-2 for OPERATOR_D1
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals2, lineNodes, OPERATOR_D1), throwCounter, nException );

    // exception #13: incorrect 0th dimension of output array (must equal number of basis functions)
    FieldContainer<double> badVals3(lineBasis.getCardinality() + 1, lineNodes.dimension(0));
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals3, lineNodes, OPERATOR_VALUE), throwCounter, nException );

    // exception #14: incorrect 1st dimension of output array (must equal number of points)
    FieldContainer<double> badVals4(lineBasis.getCardinality(), lineNodes.dimension(0) + 1);
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals4, lineNodes, OPERATOR_VALUE), throwCounter, nException );

    // exception #15: incorrect 2nd dimension of output array (must equal spatial dimension)
    FieldContainer<double> badVals5(lineBasis.getCardinality(), lineNodes.dimension(0), 2);
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals5, lineNodes, OPERATOR_GRAD), throwCounter, nException );

    // not an exception
    FieldContainer<double> goodVals2(lineBasis.getCardinality(), lineNodes.dimension(0));
    INTREPID_TEST_COMMAND( lineBasis.getValues(goodVals2, lineNodes, OPERATOR_VALUE), throwCounter, nException ); --nException;
#endif

  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

  // Check if number of thrown exceptions matches the one we expect
  if (throwCounter != nException) {
    errorFlag++;
    *outStream << std::setw(70) << "FAILURE! Incorrect number of exceptions." << "\n";
  }


  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 3: orthogonality of basis functions                                    |\n"\
    << "===============================================================================\n";

  outStream -> precision(20);

  try {

    // Check orthogonality property for Legendre polynomials.
    int maxorder = 10;

    DefaultCubatureFactory<double>  cubFactory;                                   // create factory
    shards::CellTopology line(shards::getCellTopologyData< shards::Line<> >());   // create cell topology
    for (int ordi=0; ordi < maxorder; ordi++) {
      //create left basis
      Teuchos::RCP<Basis<double,FieldContainer<double> > > lineBasisLeft =
        Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM_JACOBI<double,FieldContainer<double> >(ordi) );

      for (int ordj=0; ordj < maxorder; ordj++) {

        //create right basis
        Teuchos::RCP<Basis<double,FieldContainer<double> > > lineBasisRight =
          Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM_JACOBI<double,FieldContainer<double> >(ordj) );

        // get cubature points and weights
        Teuchos::RCP<Cubature<double> > lineCub = cubFactory.create(line, ordi+ordj);
        int numPoints = lineCub->getNumPoints();
        FieldContainer<double> cubPoints (numPoints, lineCub->getDimension());
        FieldContainer<double> cubWeights(numPoints);
        FieldContainer<double> cubWeightsC(1, numPoints);
        lineCub->getCubature(cubPoints, cubWeights);
        // "reshape" weights
        for (int i=0; i<numPoints; i++) { cubWeightsC(0,i) = cubWeights(i); }
        

        // get basis values
        int numFieldsLeft  = lineBasisLeft ->getCardinality();
        int numFieldsRight = lineBasisRight->getCardinality();
        FieldContainer<double> valsLeft(numFieldsLeft,numPoints),
                               valsRight(numFieldsRight,numPoints);
        lineBasisLeft ->getValues(valsLeft,  cubPoints, OPERATOR_VALUE);
        lineBasisRight->getValues(valsRight, cubPoints, OPERATOR_VALUE);

        // reshape by cloning and integrate
        FieldContainer<double> valsLeftC(1, numFieldsLeft,numPoints),
                               valsRightC(1, numFieldsRight,numPoints),
                               massMatrix(1, numFieldsLeft, numFieldsRight);
        ArrayTools::cloneFields<double>(valsLeftC, valsLeft);
        ArrayTools::cloneFields<double>(valsRightC, valsRight);
        ArrayTools::scalarMultiplyDataField<double>(valsRightC, cubWeightsC, valsRightC);
        FunctionSpaceTools::integrate<double>(massMatrix, valsLeftC, valsRightC, COMP_CPP);

        // check orthogonality property
        for (int i=0; i<numFieldsLeft; i++) {
          for (int j=0; j<numFieldsRight; j++) {

            if (i==j) {
              if ( std::abs(massMatrix(0,i,j)-(double)(2.0/(2.0*j+1.0))) > INTREPID_TOL ) {
                *outStream << "Incorrect ii (\"diagonal\") value for i=" << i << ", j=" << j << ": "
                           << massMatrix(0,i,j) << " != " << "2/(2*" << j << "+1)\n\n";
                errorFlag++;
              }
            }
            else {
              if ( std::abs(massMatrix(0,i,j)) > INTREPID_TOL ) {
                *outStream << "Incorrect ij (\"off-diagonal\") value for i=" << i << ", j=" << j << ": "
                           << massMatrix(0,i,j) << " != " << "0\n\n";
                errorFlag++;
              }
            }
          }
        }

      }
    }

  }
  // Catch unexpected errors
  catch (std::logic_error err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };

  *outStream \
    << "\n"
    << "===============================================================================\n"\
    << "| TEST 4: correctness of basis function derivatives                           |\n"\
    << "===============================================================================\n";

  outStream -> precision(20);

  // function values stored by bf, then pt
  double basisValues[] = {
    1.000000000000000, 1.000000000000000, 1.000000000000000,	\
    1.000000000000000, -1.000000000000000, -0.3333333333333333, \
    0.3333333333333333, 1.000000000000000, 1.000000000000000,	\
    -0.3333333333333333, -0.3333333333333333, 1.000000000000000,	\
    -1.000000000000000, 0.4074074074074074, -0.4074074074074074,	\
    1.000000000000000};

  double basisD1Values[] = 
    {0, 0, 0, 0, 1.000000000000000, 1.000000000000000, 1.000000000000000, \
     1.000000000000000, -3.000000000000000, -1.000000000000000,		\
     1.000000000000000, 3.000000000000000, 6.000000000000000,		\
     -0.6666666666666667, -0.6666666666666667, 6.000000000000000};

  double basisD2Values[] = 
    {0, 0, 0, 0, 0, 0, 0, 0, 3.000000000000000, 3.000000000000000,	\
     3.000000000000000, 3.000000000000000, -15.00000000000000,		\
     -5.000000000000000, 5.000000000000000, 15.00000000000000};

  double basisD3Values[] = 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 15.00000000000000,		\
     15.00000000000000, 15.00000000000000, 15.00000000000000};
  


  try {
    Basis_HGRAD_LINE_Cn_FEM_JACOBI<double, FieldContainer<double> > lineBasis3(3, alpha, beta);
    int numIntervals = 3;
    FieldContainer<double> lineNodes3(numIntervals+1, 1);
    FieldContainer<double> vals;
    for (int i=0; i<numIntervals+1; i++) {
      lineNodes3(i,0) = -1.0+(2.0*(double)i)/(double)numIntervals;
    }
    int numFields = lineBasis3.getCardinality();
    int numPoints = lineNodes3.dimension(0);

    // test basis values
    vals.resize(numFields, numPoints);
    lineBasis3.getValues(vals,lineNodes3,OPERATOR_VALUE);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
        
        // Compute offset for (F,P) container
        int l =  j + i * numPoints;
	if (std::abs(vals(i,j) - basisValues[l]) > INTREPID_TOL) {
	  errorFlag++;
	  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
	  
	  // Output the multi-index of the value where the error is:
	  *outStream << " At multi-index { ";
	  *outStream << i << " ";*outStream << j << " ";
	  *outStream << "}  computed value: " << vals(i,j)
		     << " but reference value: " << basisValues[l] << "\n";
	}
      }
    }

    // test basis derivatives
    vals.resize(numFields, numPoints,1);
    lineBasis3.getValues(vals,lineNodes3,OPERATOR_D1);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
        
        // Compute offset for (F,P) container
        int l =  j + i * numPoints;
	if (std::abs(vals(i,j,0) - basisD1Values[l]) > INTREPID_TOL) {
	  errorFlag++;
	  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
	  
	  // Output the multi-index of the value where the error is:
	  *outStream << " At multi-index { ";
	  *outStream << i << " ";*outStream << j << " ";
	  *outStream << "}  computed value: " << vals(i,j,0)
		     << " but reference value: " << basisD1Values[l] << "\n";
	}
      }
    }

    vals.resize(numFields, numPoints,1);
    lineBasis3.getValues(vals,lineNodes3,OPERATOR_D2);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
        
        // Compute offset for (F,P) container
        int l =  j + i * numPoints;
	if (std::abs(vals(i,j,0) - basisD2Values[l]) > INTREPID_TOL) {
	  errorFlag++;
	  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
	  
	  // Output the multi-index of the value where the error is:
	  *outStream << " At multi-index { ";
	  *outStream << i << " ";*outStream << j << " ";
	  *outStream << "}  computed value: " << vals(i,j,0)
		     << " but reference value: " << basisD2Values[l] << "\n";
	}
      }
    }

    vals.resize(numFields, numPoints,1);
    lineBasis3.getValues(vals,lineNodes3,OPERATOR_D3);
    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
        
        // Compute offset for (F,P) container
        int l =  j + i * numPoints;
	if (std::abs(vals(i,j,0) - basisD3Values[l]) > INTREPID_TOL) {
	  errorFlag++;
	  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
	  
	  // Output the multi-index of the value where the error is:
	  *outStream << " At multi-index { ";
	  *outStream << i << " ";*outStream << j << " ";
	  *outStream << "}  computed value: " << vals(i,j,0)
		     << " but reference value: " << basisD3Values[l] << "\n";
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
