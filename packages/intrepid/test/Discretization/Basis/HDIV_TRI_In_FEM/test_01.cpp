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
#include "Shards_CellTopology.hpp"

#include <iostream>
using namespace Intrepid;

/** \brief Performs a code-code comparison to FIAT for Raviart-Thomas bases on triangles (values and divs)
    \param argc [in] - number of command-line arguments
    \param argv [in] - command-line arguments
 */
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

  // test for basis values, compare against fiat
  try {
    const int deg = 2;
    Basis_HDIV_TRI_In_FEM<double,FieldContainer<double> >  myBasis( deg , POINTTYPE_EQUISPACED );

    // Get a lattice
    const int np_lattice = PointTools::getLatticeSize( myBasis.getBaseCellTopology() , deg , 0 );
    FieldContainer<double> lattice( np_lattice , 2 );
    FieldContainer<double> myBasisValues( myBasis.getCardinality() , np_lattice , 2 );
    PointTools::getLattice<double,FieldContainer<double> >( lattice , 
							    myBasis.getBaseCellTopology() , 
							    deg , 
							    0 , 
							    POINTTYPE_EQUISPACED );    

    myBasis.getValues( myBasisValues , lattice , OPERATOR_VALUE );

    const double fiat_vals[] = {
      5.194805194805292e-02, -2.000000000000000e+00,
      2.272727272727657e-02, -5.000000000000001e-01,
      -1.233766233766231e+00, 1.000000000000000e+00,
      -6.493506493506596e-03, -1.493506493506507e-01,
      -2.987012987013012e-01, 3.279220779220801e-01,
      5.194805194805287e-02, -2.857142857142830e-01,
      2.337662337662333e-01, 9.999999999999997e-01,
      4.772727272727243e-01, -5.000000000000000e-01,
      1.948051948051946e+00, -2.000000000000000e+00,
      -2.922077922077916e-02, 3.279220779220787e-01,
      1.558441558441582e-01, -1.493506493506513e-01,
      2.337662337662332e-01, -2.857142857142884e-01,
      2.337662337662333e-01, -1.110223024625157e-16,
      4.772727272727246e-01, 2.775557561562891e-17,
      1.948051948051946e+00, -5.551115123125783e-17,
      -2.922077922077895e-02, -6.720779220779212e-01,
      1.558441558441579e-01, 3.506493506493489e-01,
      2.337662337662326e-01, -1.285714285714288e+00,
      -2.857142857142860e-01, -5.551115123125783e-17,
      -5.000000000000008e-01, -5.551115123125783e-17,
      -7.142857142857150e-01, -5.551115123125783e-17,
      3.571428571428574e-02, 3.214285714285722e-01,
      1.428571428571434e-01, 3.214285714285712e-01,
      -2.857142857142860e-01, 2.571428571428572e+00,
      7.142857142857146e-01, -1.665334536937735e-16,
      4.999999999999997e-01, -1.110223024625157e-16,
      2.857142857142857e-01, -5.551115123125783e-17,
      -4.642857142857141e-01, 3.214285714285715e-01,
      -3.571428571428567e-01, 3.214285714285710e-01,
      -2.285714285714285e+00, 2.571428571428571e+00,
      -1.948051948051947e+00, 0.000000000000000e+00,
      -4.772727272727238e-01, -1.110223024625157e-16,
      -2.337662337662305e-01, -1.110223024625157e-16,
      -5.064935064935067e-01, 3.506493506493494e-01,
      7.012987012986992e-01, -6.720779220779199e-01,
      1.051948051948052e+00, -1.285714285714282e+00,
      -1.818181818181798e-01, 2.220446049250313e-16,
      1.045454545454552e+00, 5.648254401576814e-17,
      -1.818181818181759e-01, 2.220446049250313e-16,
      2.272727272727237e-02, 1.022727272727270e+00,
      1.045454545454541e+00, -1.022727272727269e+00,
      -1.818181818181795e-01, 2.953271850160422e-15,
      3.376623376623393e-01, 2.220446049250313e-16,
      5.227272727272771e-01, 1.110223024625157e-16,
      -5.194805194805161e-01, 2.220446049250313e-16,
      -4.220779220779241e-02, 1.529220779220777e+00,
      -4.415584415584443e-01, 5.064935064935090e-01,
      3.376623376623393e-01, -8.571428571428555e-01
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
    Basis_HDIV_TRI_In_FEM<double,FieldContainer<double> >  myBasis( deg , POINTTYPE_EQUISPACED );

    // Get a lattice
    const int np_lattice = PointTools::getLatticeSize( myBasis.getBaseCellTopology() , deg , 0 );
    FieldContainer<double> lattice( np_lattice , 2 );
    FieldContainer<double> myBasisDivs( myBasis.getCardinality() , np_lattice );
    PointTools::getLattice<double,FieldContainer<double> >( lattice , 
							    myBasis.getBaseCellTopology() , 
							    deg , 
							    0 , 
							    POINTTYPE_EQUISPACED );    

    myBasis.getValues( myBasisDivs , lattice , OPERATOR_DIV );

    const double fiat_divs[] = {
      6.857142857142859e+00,
      2.357142857142859e+00,
      -2.142857142857144e+00,
      2.357142857142858e+00,
      -2.142857142857142e+00,
      -2.142857142857140e+00,
      -2.142857142857147e+00,
      2.357142857142856e+00,
      6.857142857142859e+00,
      -2.142857142857144e+00,
      2.357142857142856e+00,
      -2.142857142857143e+00,
      -2.142857142857145e+00,
      2.357142857142855e+00,
      6.857142857142858e+00,
      -2.142857142857144e+00,
      2.357142857142855e+00,
      -2.142857142857142e+00,
      -1.714285714285716e+00,
      -1.714285714285716e+00,
      -1.714285714285718e+00,
      2.785714285714288e+00,
      2.785714285714286e+00,
      7.285714285714281e+00,
      -1.714285714285718e+00,
      -1.714285714285716e+00,
      -1.714285714285716e+00,
      2.785714285714286e+00,
      2.785714285714286e+00,
      7.285714285714278e+00,
      6.857142857142857e+00,
      2.357142857142859e+00,
      -2.142857142857141e+00,
      2.357142857142858e+00,
      -2.142857142857140e+00,
      -2.142857142857136e+00,
      9.000000000000004e+00,
      3.996802888650562e-15,
      -9.000000000000000e+00,
      4.500000000000000e+00,
      -4.499999999999996e+00,
      8.881784197001252e-16,
      8.571428571428573e+00,
      4.071428571428576e+00,
      -4.285714285714253e-01,
      -4.285714285714295e-01,
      -4.928571428571429e+00,
      -9.428571428571420e+00
    };

    int cur=0;
    for (int i=0;i<myBasisDivs.dimension(0);i++) {
      for (int j=0;j<myBasisDivs.dimension(1);j++) {
	if (std::abs( myBasisDivs(i,j) - fiat_divs[cur] ) > INTREPID_TOL ) {
	  errorFlag++;
	  *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
	  
	  // Output the multi-index of the value where the error is:
	  *outStream << " At multi-index { ";
	  *outStream << i << " " << j;
	  *outStream << "}  computed value: " << myBasisDivs(i,j)
		     << " but correct value: " << fiat_divs[cur] << "\n";
          *outStream << "Difference: " << std::abs( myBasisDivs(i,j) - fiat_divs[cur] ) << "\n";
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
