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
\brief  Unit tests for the Intrepid::Basis_HGRAD_LINE_Hermite_FEM class.
\author Created by G. von Winckel.
*/

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HGRAD_LINE_Hermite_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_PointTools.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iomanip>

using namespace std;
using namespace Intrepid;

#define INTREPID_TEST_COMMAND( S , throwCounter, nException )                                                              \
{                                                                                                                          \
  ++nException;                                                                                                            \
  try {                                                                                                                    \
    S ;                                                                                                                    \
  }                                                                                                                        \
  catch (const std::logic_error & err) {                                                                                           \
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
    << "|               Unit Test (Basis_HGRAD_LINE_Hermite_FEM)                      |\n" \
    << "|                                                                             |\n" \
    << "|     1) Conversion of Dof tags into Dof ordinals and back                    |\n" \
    << "|     2) Basis values for VALUE, GRAD, CURL, and Dk operators                 |\n" \
    << "|     3) Numerical differentiation with Herite Interpolants                   |\n" \
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
  shards::CellTopology line(shards::getCellTopologyData< shards::Line<> >());   // create cell topology
  const int deg = 5;

  FieldContainer<double> pts(PointTools::getLatticeSize(line,deg),1);
  PointTools::getLattice<double,FieldContainer<double> >(pts,line,deg);
  Basis_HGRAD_LINE_Hermite_FEM<double, FieldContainer<double> > lineBasis(pts);

  *outStream << std::endl;

  lineBasis.printTags( *outStream );

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

  *outStream << "lineBasis.getCardinality()      = " << lineBasis.getCardinality()      << std::endl;
  *outStream << "lineBasis.getDegree()           = " << lineBasis.getDegree()           << std::endl;
  *outStream << "lineBasis.getBaseCellTopology() = " << lineBasis.getBaseCellTopology() << std::endl;
  *outStream << "lineBasis.getBasisType()        = " << lineBasis.getBasisType()        << std::endl;
  *outStream << "lineBasis.getCoordinateSystem() = " << lineBasis.getCoordinateSystem() << std::endl;
  *outStream << std::endl;

  try {

#ifdef HAVE_INTREPID_DEBUG
    // exception #1: DIV cannot be applied to scalar functions
    // resize vals to rank-2 container with dimensions (num. points, num. basis functions)
    vals.resize(lineBasis.getCardinality(), lineNodes.dimension(0) );
    INTREPID_TEST_COMMAND( lineBasis.getValues(vals, lineNodes, OPERATOR_DIV), throwCounter, nException );
#endif // HAVE_INTREPID_DEBUG

    // Exceptions #2-8: all bf tags/bf Ids below are wrong and should cause getDofOrdinal() and
    // getDofTag() to access invalid array elements thereby causing bounds check exception

    // exception #2 - There are no two-dimensional subcells
    INTREPID_TEST_COMMAND( lineBasis.getDofOrdinal(2,0,0), throwCounter, nException );
 
    // exception #3 - There are at most two zero-dimensional subcells 
    INTREPID_TEST_COMMAND( lineBasis.getDofOrdinal(0,2,0), throwCounter, nException );

    // exception #4 - Zero-dimensional subcells have at most 2 DoF
    INTREPID_TEST_COMMAND( lineBasis.getDofOrdinal(0,0,2), throwCounter, nException );

    // exception #5 - There is at most one one-dimensional subcell
    INTREPID_TEST_COMMAND( lineBasis.getDofOrdinal(1,1,0), throwCounter, nException );

    // exception #6 - One-dimensional subcell cannot have DoF ordinal larger than 
    //                cardinality-1
    int C = lineBasis.getCardinality();
    INTREPID_TEST_COMMAND( lineBasis.getDofOrdinal(1,0,C), throwCounter, nException );
    
    // exception #7 - DoF cannot exceed cardinality
    INTREPID_TEST_COMMAND( lineBasis.getDofTag(C),  throwCounter, nException );

    // exception #8 - No negative indices
    INTREPID_TEST_COMMAND( lineBasis.getDofTag(-1), throwCounter, nException );

    // not an exception
    INTREPID_TEST_COMMAND( lineBasis.getDofTag(5),  throwCounter, nException ); --nException;

#ifdef HAVE_INTREPID_DEBUG
    // Exceptions 9-16 test exception handling with incorrectly dimensioned input/output arrays
    // exception #9: input points array must be of rank-2
    FieldContainer<double> badPoints1(4, 5, 3);
    INTREPID_TEST_COMMAND( lineBasis.getValues(vals, badPoints1, OPERATOR_VALUE), throwCounter, nException );

    // exception #10: dimension 1 in the input point array must equal space dimension of the cell
    FieldContainer<double> badPoints2(4, 3);
    INTREPID_TEST_COMMAND( lineBasis.getValues(vals, badPoints2, OPERATOR_VALUE), throwCounter, nException );

    // exception #11: output values must be of rank-2 for OPERATOR_VALUE
    FieldContainer<double> badVals1(4, 3, 1);
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals1, lineNodes, OPERATOR_VALUE), throwCounter, nException );

    // exception #12: output values must be of rank-3 for OPERATOR_GRAD
    FieldContainer<double> badVals2(4, 3);
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals2, lineNodes, OPERATOR_GRAD), throwCounter, nException );

    // exception #13: output values must be of rank-2 for OPERATOR_D1
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals2, lineNodes, OPERATOR_D1), throwCounter, nException );

    // exception #14: incorrect 0th dimension of output array (must equal number of basis functions)
    FieldContainer<double> badVals3(lineBasis.getCardinality() + 1, lineNodes.dimension(0));
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals3, lineNodes, OPERATOR_VALUE), throwCounter, nException );
    
    // exception #15: incorrect 1st dimension of output array (must equal number of points)
    FieldContainer<double> badVals4(lineBasis.getCardinality(), lineNodes.dimension(0) + 1);
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals4, lineNodes, OPERATOR_VALUE), throwCounter, nException );

    // exception #16: incorrect 2nd dimension of output array (must equal spatial dimension)
    FieldContainer<double> badVals5(lineBasis.getCardinality(), lineNodes.dimension(0), 2);
    INTREPID_TEST_COMMAND( lineBasis.getValues(badVals5, lineNodes, OPERATOR_GRAD), throwCounter, nException );
    
    // not an exception
    FieldContainer<double> goodVals2(lineBasis.getCardinality(), lineNodes.dimension(0));
    INTREPID_TEST_COMMAND( lineBasis.getValues(goodVals2, lineNodes, OPERATOR_VALUE), throwCounter, nException ); --nException;
#endif // HAVE_INTREPID_DEBUG
  
  }

  
  catch (const std::logic_error & err) {
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
     << "===============================================================================\n" \
     << "| TEST 2: Correctness of basis function values                                |\n" \
     << "===============================================================================\n";
   outStream -> precision(20);

   try {

     int npts=5;
     shards::CellTopology line(shards::getCellTopologyData< shards::Line<> >());   // create cell topology
     FieldContainer<double> pts(PointTools::getLatticeSize(line,npts),1);
     PointTools::getLattice<double,FieldContainer<double> >(pts,line,npts);
     Basis_HGRAD_LINE_Hermite_FEM<double, FieldContainer<double> > lineBasis(pts);

     FieldContainer<double> vals(lineBasis.getCardinality(),npts+1);
     FieldContainer<double> ders(lineBasis.getCardinality(),npts+1,1);

     lineBasis.getValues(vals,pts,OPERATOR_VALUE);
     lineBasis.getValues(ders,pts,OPERATOR_D1);

     int C = lineBasis.getCardinality();
     int n = C/2;

     // Loop over basis functions
     for (int i=0; i<n; i++) {

      // Loop over interpolation points
      for (int j=0; j<n; j++) {

        if ( i == j ) {
  
          if( std::abs( vals(2*i,j) - 1.0 ) > INTREPID_TOL ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Basis function " << 2*i << " does not have unit value at its node\n";
          }

          if( std::abs( ders(2*i,j,0) ) > INTREPID_TOL ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Basis function " << 2*i+1 << " does not have zero first derivative at its node\n";
          }

          if( std::abs( vals(2*i+1,j) ) > INTREPID_TOL ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Basis function " << 2*i+1 << " does not have zero value at its node\n";
          }
   
          if( std::abs( ders(2*i+1,j,0) - 1.0 ) > INTREPID_TOL ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Basis function " << 2*i+1 << " does not have unit first derivative at its node\n";
          }
        }
        else { // i != j

          if( std::abs( vals(2*i,j) ) > INTREPID_TOL ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Basis function " << 2*i << " does not vanish at node " << j << "\n";
          }
           
          if( std::abs( ders(2*i,j,0) ) > INTREPID_TOL ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Basis function " << 2*i << " does not have zero first derivative at node " << j << "\n";
          }

          if( std::abs( vals(2*i+1,j) ) > INTREPID_TOL ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Basis function " << 2*i+1 << " does not vanish at node " << j << "\n";
          }
           
          if( std::abs( ders(2*i+1,j,0) ) > INTREPID_TOL ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            *outStream << " Basis function " << 2*i+1 << " does not have zero first derivative at node " << j << "\n";
          }

        } // end else i != j

      } // end loop over interpolation points

    }  // end loop over basis functions

 
    *outStream << std::setprecision(4);

    *outStream << "\n\nBasis function values on nodes\n" << std::endl;

    // Loop over basis functions
    for (int i=0; i<C; i++) {

      // Loop over interpolation points
      for (int j=0; j<n; j++) {
 
        *outStream << std::setw(12) << vals(i,j);
      }
 
      *outStream << std::endl;
    }


    *outStream << "\n\nBasis function derivatives on nodes\n" << std::endl;

    // Loop over basis functions
    for (int i=0; i<C; i++) {

      // Loop over interpolation points
      for (int j=0; j<n; j++) {
 
        *outStream << std::setw(12) << ders(i,j,0);
      }
 
      *outStream << std::endl;
    }

  } // end try

  // Catch unexpected errors
  catch (const std::logic_error & err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  };

   *outStream \
     << "\n"
     << "===============================================================================\n" \
     << "| TEST 3: Correctness of basis function higher derivatives via numerical      |\n" \
     << "|         differentiation.                                                    |\n" \
     << "|                                                                             |\n" \
     << "| Let f(x_i) = sin(x_i), f'(x_i) = cos(x_i)                                   |\n" \
     << "|                                                                             |\n" \
     << "| and compare the second and third derivatives obtained by differentiating    |\n" \
     << "| the Hermite interpolating polynomial with analytical values                 |\n" \
     << "|                                                                             |\n" \
     << "===============================================================================\n";
   outStream -> precision(20);


  try {

     int npts = 6;
     shards::CellTopology line(shards::getCellTopologyData< shards::Line<> >());   // create cell topology
     FieldContainer<double> pts(PointTools::getLatticeSize(line,npts),1);
     PointTools::getGaussPoints<double,FieldContainer<double> >(pts,npts);
     Basis_HGRAD_LINE_Hermite_FEM<double, FieldContainer<double> > lineBasis(pts);

     int C = lineBasis.getCardinality();
     int n = C/2;

     FieldContainer<double> f0(n);
     FieldContainer<double> f1(n);
     FieldContainer<double> f2(n);
     FieldContainer<double> f3(n);

     FieldContainer<double> der2(C,n,1);
     lineBasis.getValues(der2,pts,OPERATOR_D2);

     FieldContainer<double> der3(C,n,1);
     lineBasis.getValues(der3,pts,OPERATOR_D3);

     // Loop over interpolation points
     for( int j=0; j<n; ++j ) {
       f0(j) = std::sin(pts(j,0));
       f1(j) = std::cos(pts(j,0));
     }

     double error2 = 0;
     double error3 = 0;

     for( int j=0; j<n; ++j ) {
        for( int i=0; i<n; ++i ) {
           f2(j) += f0(i)*der2(2*i,j,0);
           f2(j) += f1(i)*der2(2*i+1,j,0);
 
           f3(j) += f0(i)*der3(2*i,j,0);
           f3(j) += f1(i)*der3(2*i+1,j,0);
       }
       
       error2 += std::pow(f0(j)+f2(j),2);
       error3 += std::pow(f1(j)+f3(j),2);
     }
    
     error2 = std::sqrt(error2);
     error3 = std::sqrt(error3);

     *outStream << std::setprecision(16);

     int width = 24;
     std::string bar(20,'-');


     *outStream << "\n\n"
                << std::setw(width) << "x_i" 
                << std::setw(width) << "exact f(x_i)" 
                << std::setw(width) << "exact f'(x_i)" 
                << std::setw(width) << "computed f\"(x_i)" 
                << std::setw(width) << "computed f'\"(x_i)" << std::endl;

     *outStream << std::setw(width) << bar
                << std::setw(width) << bar
                << std::setw(width) << bar
                << std::setw(width) << bar
                << std::setw(width) << bar << std::endl;


     // Loop over interpolation points
     for (int j=0; j<n; j++) {
 
       *outStream << std::setw(width) << pts(j,0)
                  << std::setw(width) << f0(j)
                  << std::setw(width) << f1(j)
                  << std::setw(width) << f2(j)
                  << std::setw(width) << f3(j) << std::endl;
     }
 
     double errtol = 1e-9;


     *outStream << std::endl;
     *outStream << "|f+f\"| = " << error2 << std::endl;
     *outStream << "|f'+f'\"| = " << error3 << std::endl;

     if( error2 > errtol ) {
       errorFlag++;
       *outStream << std::setw(70) << "FAILURE! Second derivative not within tolerance " << errtol << "\n";
     }
     if( error3 > errtol ) {
       errorFlag++;
       *outStream << std::setw(70) << "FAILURE! Third derivative not within tolerance " << errtol << "\n";
     }

  }

  // Catch unexpected errors
  catch (const std::logic_error & err) {
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
