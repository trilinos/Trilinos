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
\brief  Unit test of Nedelec class for triangles.
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

/** \brief Performs a code-code comparison to FIAT for Nedelec bases on triangles (values and curls)
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
    << "|                           Unit Test HCURL_TRI_In_FEM                        |\n" \
    << "|                                                                             |\n" \
    << "|     1) Tests triangular Nedelec-Thomas basis                                |\n" \
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
2.000000000000000e+00, 0.000000000000000e+00,
5.000000000000000e-01, 2.500000000000000e-01,
-1.000000000000000e+00, -1.000000000000000e+00,
2.500000000000000e-01, 0.000000000000000e+00,
-5.000000000000000e-01, -5.000000000000000e-01,
0.000000000000000e+00, 0.000000000000000e+00,
-1.000000000000000e+00, 0.000000000000000e+00,
5.000000000000000e-01, 2.500000000000000e-01,
2.000000000000000e+00, 2.000000000000000e+00,
-5.000000000000000e-01, 0.000000000000000e+00,
2.500000000000000e-01, 2.500000000000000e-01,
0.000000000000000e+00, 0.000000000000000e+00,
0.000000000000000e+00, 0.000000000000000e+00,
0.000000000000000e+00, 2.500000000000000e-01,
0.000000000000000e+00, 2.000000000000000e+00,
5.000000000000000e-01, 0.000000000000000e+00,
-2.500000000000000e-01, 2.500000000000000e-01,
1.000000000000000e+00, 0.000000000000000e+00,
0.000000000000000e+00, 0.000000000000000e+00,
0.000000000000000e+00, -5.000000000000000e-01,
0.000000000000000e+00, -1.000000000000000e+00,
-2.500000000000000e-01, 0.000000000000000e+00,
-2.500000000000000e-01, 2.500000000000000e-01,
-2.000000000000000e+00, 0.000000000000000e+00,
0.000000000000000e+00, 1.000000000000000e+00,
0.000000000000000e+00, 5.000000000000000e-01,
0.000000000000000e+00, 0.000000000000000e+00,
-2.500000000000000e-01, -5.000000000000000e-01,
-2.500000000000000e-01, -2.500000000000000e-01,
-2.000000000000000e+00, -2.000000000000000e+00,
0.000000000000000e+00, -2.000000000000000e+00,
0.000000000000000e+00, -2.500000000000000e-01,
0.000000000000000e+00, 0.000000000000000e+00,
-2.500000000000000e-01, -5.000000000000000e-01,
5.000000000000000e-01, 5.000000000000000e-01,
1.000000000000000e+00, 1.000000000000000e+00,
0.000000000000000e+00, 0.000000000000000e+00,
0.000000000000000e+00, -7.500000000000000e-01,
0.000000000000000e+00, 0.000000000000000e+00,
1.500000000000000e+00, 0.000000000000000e+00,
7.500000000000000e-01, 7.500000000000000e-01,
0.000000000000000e+00, 0.000000000000000e+00,
0.000000000000000e+00, 0.000000000000000e+00,
0.000000000000000e+00, 1.500000000000000e+00,
0.000000000000000e+00, 0.000000000000000e+00,
-7.500000000000000e-01, 0.000000000000000e+00,
7.500000000000000e-01, 7.500000000000000e-01,
0.000000000000000e+00, 0.000000000000000e+00
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
-7.000000000000000e+00,
-2.500000000000000e+00,
2.000000000000000e+00,
-2.500000000000000e+00,
2.000000000000000e+00,
2.000000000000000e+00,
2.000000000000000e+00,
-2.500000000000000e+00,
-7.000000000000000e+00,
2.000000000000000e+00,
-2.500000000000000e+00,
2.000000000000000e+00,
2.000000000000000e+00,
-2.500000000000000e+00,
-7.000000000000000e+00,
2.000000000000000e+00,
-2.500000000000000e+00,
2.000000000000000e+00,
2.000000000000000e+00,
2.000000000000000e+00,
2.000000000000000e+00,
-2.500000000000000e+00,
-2.500000000000000e+00,
-7.000000000000000e+00,
2.000000000000000e+00,
2.000000000000000e+00,
2.000000000000000e+00,
-2.500000000000000e+00,
-2.500000000000000e+00,
-7.000000000000000e+00,
-7.000000000000000e+00,
-2.500000000000000e+00,
2.000000000000000e+00,
-2.500000000000000e+00,
2.000000000000000e+00,
2.000000000000000e+00,
9.000000000000000e+00,
4.500000000000000e+00,
0.000000000000000e+00,
0.000000000000000e+00,
-4.500000000000000e+00,
-9.000000000000000e+00,
-9.000000000000000e+00,
0.000000000000000e+00,
9.000000000000000e+00,
-4.500000000000000e+00,
4.500000000000000e+00,
0.000000000000000e+00
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
