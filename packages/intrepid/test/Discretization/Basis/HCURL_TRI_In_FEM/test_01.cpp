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
//                    Denis Ridzal (dridzal@sandia.gov) or
//                    Robert Kirby (robert.c.kirby@ttu.edu)
//
// ************************************************************************
// @HEADER


/** \file
\brief  Unit test of Raviart-Thomas class for triangles.
\author Created by R. Kirby
*/

#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Intrepid_PointTools.hpp"
#include "Intrepid_HDIV_TRI_In_FEM.hpp"
#include "Intrepid_HCURL_TRI_In_FEM.hpp"
#include "Shards_CellTopology.hpp"

#include <iostream>
using namespace Intrepid;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  
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
    << "|                           Unit Test HDIV_TRI_In_FEM                         |\n" \
    << "|                                                                             |\n" \
    << "|     1) Tests triangular Raviart-Thomas basis                                |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" \
    << "|                      Robert Kirby (robert.c.kirby@ttu.edu)                  |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";
  
  int errorFlag  = 0;

  // test for basis values, compare against H(div) basis, for they should be related by rotation
  // in the lowest order case
  try {
    const int deg = 1;
    Basis_HDIV_TRI_In_FEM<double,FieldContainer<double> >  myBasisDIV( deg , POINTTYPE_EQUISPACED );
    Basis_HCURL_TRI_In_FEM<double,FieldContainer<double> > myBasisCURL( deg , POINTTYPE_EQUISPACED );

    // Get a lattice
    const int np_lattice = PointTools::getLatticeSize( myBasisDIV.getBaseCellTopology() , deg , 0 );
    FieldContainer<double> lattice( np_lattice , 2 );

    FieldContainer<double> myBasisValuesDIV( myBasisDIV.getCardinality() , np_lattice , 2 );
    FieldContainer<double> myBasisValuesCURL( myBasisDIV.getCardinality() , np_lattice , 2 );
    PointTools::getLattice<double,FieldContainer<double> >( lattice , 
							    myBasisDIV.getBaseCellTopology() , 
							    deg , 
							    0 , 
							    POINTTYPE_EQUISPACED );    

    myBasisDIV.getValues( myBasisValuesDIV , lattice , OPERATOR_VALUE );
    myBasisCURL.getValues( myBasisValuesCURL, lattice , OPERATOR_VALUE );

    int cur=0;
    for (int i=0;i<myBasisValuesDIV.dimension(0);i++) {
      for (int j=0;j<myBasisValuesDIV.dimension(1);j++) {
	if (std::abs( myBasisValuesDIV(i,j,1) + myBasisValuesCURL(i,j,0) ) > INTREPID_TOL ) {
	  errorFlag++;
	  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";	  
	  // Output the multi-index of the value where the error is:
	  *outStream << " At multi-index { ";
	  *outStream << i << " " << j << " and component 0";
	  *outStream << "}  computed value: " << myBasisValuesCURL(i,j,0)
		     << " but correct value: " << -myBasisValuesDIV(i,j,1) << "\n";
          *outStream << "Difference: " << std::abs( myBasisValuesCURL(i,j,0) + myBasisValuesDIV(i,j,1) ) << "\n";
	}
	if (std::abs( myBasisValuesDIV(i,j,0) - myBasisValuesCURL(i,j,1) ) > INTREPID_TOL ) {
	  errorFlag++;
	  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";	  
	  // Output the multi-index of the value where the error is:
	  *outStream << " At multi-index { ";
	  *outStream << i << " " << j << " and component 1";
	  *outStream << "}  computed value: " << myBasisValuesCURL(i,j,1)
		     << " but correct value: " << myBasisValuesDIV(i,j,0) << "\n";
          *outStream << "Difference: " << std::abs( myBasisValuesCURL(i,j,1) - myBasisValuesDIV(i,j,1) ) << "\n";
	}
      }
    }
  }
  catch (std::exception err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }

  // next test: code-code comparison with FIAT 
  try {
    const int deg = 2;
    Basis_HCURL_TRI_In_FEM<double,FieldContainer<double> > myBasis( deg , POINTTYPE_EQUISPACED );
    const int np_lattice = PointTools::getLatticeSize( myBasis.getBaseCellTopology() , deg , 0 );
    FieldContainer<double> lattice( np_lattice , 2 );
    PointTools::getLattice<double,FieldContainer<double> >( lattice , 
							    myBasis.getBaseCellTopology() , 
							    deg , 
							    0 , 
							    POINTTYPE_EQUISPACED );        
    FieldContainer<double> myBasisValues( myBasis.getCardinality() , np_lattice , 2 );

    
    myBasis.getValues( myBasisValues, lattice , OPERATOR_VALUE );

    const double fiat_vals[] = {
      2.000000000000003e+00, 5.194805194804781e-02,
      4.999999999999997e-01, 2.272727272726738e-02,
      -9.999999999999971e-01, -1.233766233766239e+00,
      1.493506493506511e-01, -6.493506493506069e-03,
      -3.279220779220769e-01, -2.987012987012970e-01,
      2.857142857142894e-01, 5.194805194804817e-02,
      -1.000000000000003e+00, 2.337662337662381e-01,
      5.000000000000004e-01, 4.772727272727326e-01,
      1.999999999999998e+00, 1.948051948051953e+00,
      -3.279220779220798e-01, -2.922077922077961e-02,
      1.493506493506487e-01, 1.558441558441543e-01,
      2.857142857142822e-01, 2.337662337662372e-01,
      -3.136380044566067e-15, 2.337662337662380e-01,
      2.498001805406602e-16, 4.772727272727328e-01,
      -2.831068712794149e-15, 1.948051948051952e+00,
      6.720779220779202e-01, -2.922077922077976e-02,
      -3.506493506493512e-01, 1.558441558441542e-01,
      1.285714285714283e+00, 2.337662337662378e-01,
      1.110223024625157e-16, -2.857142857142858e-01,
      0.000000000000000e+00, -4.999999999999999e-01,
      -1.110223024625157e-16, -7.142857142857142e-01,
      -3.214285714285714e-01, 3.571428571428573e-02,
      -3.214285714285715e-01, 1.428571428571428e-01,
      -2.571428571428573e+00, -2.857142857142857e-01,
      4.440892098500626e-16, 7.142857142857150e-01,
      5.551115123125783e-17, 5.000000000000003e-01,
      0.000000000000000e+00, 2.857142857142864e-01,
      -3.214285714285709e-01, -4.642857142857139e-01,
      -3.214285714285718e-01, -3.571428571428574e-01,
      -2.571428571428572e+00, -2.285714285714286e+00,
      2.525757381022231e-15, -1.948051948051953e+00,
      -3.608224830031759e-16, -4.772727272727331e-01,
      2.414735078559715e-15, -2.337662337662392e-01,
      -3.506493506493492e-01, -5.064935064935062e-01,
      6.720779220779231e-01, 7.012987012987030e-01,
      1.285714285714289e+00, 1.051948051948047e+00,
      -2.220446049250313e-15, -3.376623376623336e-01,
      4.440892098500626e-16, -5.227272727272672e-01,
      -3.330669073875470e-15, 5.194805194805240e-01,
      1.529220779220778e+00, 4.220779220779172e-02,
      5.064935064935057e-01, 4.415584415584397e-01,
      -8.571428571428609e-01, -3.376623376623332e-01,
      5.689893001203927e-15, -1.818181818181906e-01,
      -7.179069302488778e-16, 1.045454545454534e+00,
      5.467848396278896e-15, -1.818181818181913e-01,
      -1.022727272727269e+00, 2.272727272727360e-02,
      1.022727272727274e+00, 1.045454545454549e+00,
      6.541777753180100e-15, -1.818181818181904e-01
    };

    int cur=0;
    for (int i=0;i<myBasisValues.dimension(0);i++) {
      for (int j=0;j<myBasisValues.dimension(1);j++) {
        for (int k=0;k<myBasisValues.dimension(2);k++) {
          if (std::abs( myBasisValues(i,j,k) - fiat_vals[cur] ) > INTREPID_TOL ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            
            // Output the multi-index of the value where the error is:
            *outStream << " At multi-index { ";
            *outStream << i << " " << j << " " << k;
            *outStream << "}  computed value: " << myBasisValues(i,j,k)
                       << " but correct value: " << fiat_vals[cur] << "\n";
          *outStream << "Difference: " << std::abs( myBasisValues(i,j,k) - fiat_vals[cur] ) << "\n";
          }
          cur++;
        }
      }
    }
  }
  catch (std::exception err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }

  try {
    const int deg = 2;
    Basis_HCURL_TRI_In_FEM<double,FieldContainer<double> > myBasis( deg , POINTTYPE_EQUISPACED );
    const int np_lattice = PointTools::getLatticeSize( myBasis.getBaseCellTopology() , deg , 0 );
    FieldContainer<double> lattice( np_lattice , 2 );
    PointTools::getLattice<double,FieldContainer<double> >( lattice , 
							    myBasis.getBaseCellTopology() , 
							    deg , 
							    0 , 
							    POINTTYPE_EQUISPACED );        
    FieldContainer<double> myBasisValues( myBasis.getCardinality() , np_lattice );

    
    myBasis.getValues( myBasisValues, lattice , OPERATOR_CURL );

    const double fiat_curls[] = {
      -6.857142857142859e+00,
      -2.357142857142857e+00,
      2.142857142857144e+00,
      -2.357142857142857e+00,
      2.142857142857144e+00,
      2.142857142857144e+00,
      2.142857142857145e+00,
      -2.357142857142857e+00,
      -6.857142857142859e+00,
      2.142857142857143e+00,
      -2.357142857142859e+00,
      2.142857142857141e+00,
      2.142857142857142e+00,
      -2.357142857142858e+00,
      -6.857142857142857e+00,
      2.142857142857143e+00,
      -2.357142857142856e+00,
      2.142857142857145e+00,
      1.714285714285715e+00,
      1.714285714285715e+00,
      1.714285714285715e+00,
      -2.785714285714287e+00,
      -2.785714285714286e+00,
      -7.285714285714288e+00,
      1.714285714285717e+00,
      1.714285714285714e+00,
      1.714285714285713e+00,
      -2.785714285714285e+00,
      -2.785714285714287e+00,
      -7.285714285714286e+00,
      -6.857142857142859e+00,
      -2.357142857142857e+00,
      2.142857142857144e+00,
      -2.357142857142858e+00,
      2.142857142857144e+00,
      2.142857142857143e+00,
      8.571428571428569e+00,
      4.071428571428571e+00,
      -4.285714285714279e-01,
      -4.285714285714297e-01,
      -4.928571428571429e+00,
      -9.428571428571431e+00,
      -8.999999999999996e+00,
      8.881784197001252e-16,
      9.000000000000000e+00,
      -4.499999999999999e+00,
      4.500000000000000e+00,
      -8.881784197001252e-16
    };

    int cur=0;
    for (int i=0;i<myBasisValues.dimension(0);i++) {
      for (int j=0;j<myBasisValues.dimension(1);j++) {
	if (std::abs( myBasisValues(i,j) - fiat_curls[cur] ) > INTREPID_TOL ) {
	  errorFlag++;
	  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
	  
	  // Output the multi-index of the value where the error is:
	  *outStream << " At multi-index { ";
	  *outStream << i << " " << j;
	  *outStream << "}  computed value: " << myBasisValues(i,j)
		     << " but correct value: " << fiat_curls[cur] << "\n";
          *outStream << "Difference: " << std::abs( myBasisValues(i,j) - fiat_curls[cur] ) << "\n";
	}
	cur++;
      }
    }
  }
  catch (std::exception err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }  

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";
  
  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  
  return errorFlag;
}
