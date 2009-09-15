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
    \brief  Unit test of Nedelec class for tetrahedra.
    \author Created by R. Kirby
*/

#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Intrepid_PointTools.hpp"
#include "Intrepid_HCURL_TET_In_FEM.hpp"
#include "Shards_CellTopology.hpp"

#include <iostream>
using namespace Intrepid;

/** \brief Performs a code-code comparison to FIAT for Nedelec bases on tets (values and curls)
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
    << "|                           Unit Test HCURL_TET_In_FEM                        |\n" \
    << "|                                                                             |\n" \
    << "|     1) Tests tetrahedral Nedelec basis                                      |\n" \
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

  // code-code comparison with FIAT 
  try {
    const int deg = 1;
    Basis_HCURL_TET_In_FEM<double,FieldContainer<double> > myBasis( deg , POINTTYPE_EQUISPACED );

    const int np_lattice = PointTools::getLatticeSize( myBasis.getBaseCellTopology() , deg , 0 );
    FieldContainer<double> lattice( np_lattice , 3 );
    FieldContainer<double> myBasisValues( myBasis.getCardinality() , np_lattice , 3 );
    PointTools::getLattice<double,FieldContainer<double> >( lattice , 
                                                            myBasis.getBaseCellTopology() , 
                                                            deg , 
                                                            0 , 
                                                            POINTTYPE_EQUISPACED );    

    myBasis.getValues( myBasisValues , lattice , OPERATOR_VALUE );

    const double fiat_vals[] = {
      1.000000000000000e+00, 2.499999999999997e-01, 3.214285714285716e-01,
      9.999999999999998e-01, 7.500000000000002e-01, 6.785714285714289e-01,
      -4.999999999999997e-01, -2.500000000000003e-01, -2.500000000000002e-01,
      -6.428571428571425e-01, -3.214285714285715e-01, -3.214285714285716e-01,
      3.053113317719180e-16, 2.499999999999998e-01, 1.984126984126992e-02,
      -2.775557561562891e-17, 7.500000000000002e-01, 9.126984126984138e-02,
      -1.500000000000000e+00, -2.500000000000000e-01, -3.611111111111112e-01,
      7.142857142857156e-02, -1.309523809523809e-01, -1.984126984126987e-02,
      2.220446049250313e-16, -7.500000000000001e-01, -9.126984126984120e-02,
      -1.110223024625157e-16, -2.499999999999995e-01, -1.984126984126991e-02,
      -1.499999999999999e+00, -1.250000000000000e+00, -1.138888888888889e+00,
      7.142857142857162e-02, 2.023809523809523e-01, 9.126984126984133e-02,
      -1.665334536937735e-16, -5.551115123125783e-17, 5.873015873015872e-01,
      5.551115123125783e-17, 5.551115123125783e-17, 3.015873015873011e-01,
      -1.665334536937735e-16, -5.551115123125783e-17, 1.111111111111109e-01,
      1.714285714285714e+00, 1.523809523809524e+00, 1.412698412698413e+00,
      2.775557561562891e-16, -1.249000902703301e-16, 3.015873015873016e-01,
      5.551115123125783e-17, -5.551115123125783e-17, 5.873015873015877e-01,
      -5.551115123125783e-17, -1.942890293094024e-16, 1.111111111111111e-01,
      -1.714285714285714e+00, -1.904761904761906e-01, -3.015873015873017e-01,
      5.994282761129235e-17, 5.551115123125783e-17, 1.111111111111111e-01,
      5.994282761129230e-17, 5.551115123125783e-17, 1.111111111111112e-01,
      -1.065906260824812e-16, -5.551115123125783e-17, 7.777777777777778e-01,
      -4.414058094731610e-17, -1.333333333333333e+00, -1.111111111111112e-01
    };

    int cur=0;
    for (int i=0;i<myBasisValues.dimension(0);i++) {
      for (int j=0;j<myBasisValues.dimension(1);j++) {
        for (int k=0;k<myBasisValues.dimension(2);k++) {
          if (std::abs( myBasisValues(i,j,k) - fiat_vals[cur] ) > 10.0*INTREPID_TOL ) {
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
    const int deg = 1;
    Basis_HCURL_TET_In_FEM<double,FieldContainer<double> > myBasis( deg , POINTTYPE_EQUISPACED );

    const int np_lattice = PointTools::getLatticeSize( myBasis.getBaseCellTopology() , deg , 0 );
    FieldContainer<double> lattice( np_lattice , 3 );
    FieldContainer<double> myBasisValues( myBasis.getCardinality() , np_lattice , 3 );
    PointTools::getLattice<double,FieldContainer<double> >( lattice , 
                                                            myBasis.getBaseCellTopology() , 
                                                            deg , 
                                                            0 , 
                                                            POINTTYPE_EQUISPACED );    

    myBasis.getValues( myBasisValues , lattice , OPERATOR_CURL );

    const double fiat_curls[] = {
      -5.551115123125783e-16, -2.000000000000000e+00, 2.000000000000000e+00,
      -5.551115123125783e-16, -2.000000000000000e+00, 2.000000000000000e+00,
      -5.551115123125783e-16, -2.000000000000000e+00, 2.000000000000000e+00,
      -5.551115123125783e-16, -2.000000000000000e+00, 2.000000000000000e+00,
      -4.440892098500626e-16, -2.775557561562891e-16, 2.000000000000000e+00,
      -4.440892098500626e-16, -2.775557561562891e-16, 2.000000000000000e+00,
      -4.440892098500626e-16, -2.775557561562891e-16, 2.000000000000000e+00,
      -4.440892098500626e-16, -2.775557561562891e-16, 2.000000000000000e+00,
      -2.000000000000000e+00, 5.551115123125783e-17, 2.000000000000000e+00,
      -2.000000000000000e+00, 5.551115123125783e-17, 2.000000000000000e+00,
      -2.000000000000000e+00, 5.551115123125783e-17, 2.000000000000000e+00,
      -2.000000000000000e+00, 5.551115123125783e-17, 2.000000000000000e+00,
      -2.000000000000000e+00, 2.000000000000000e+00, 9.861075762086680e-17,
      -2.000000000000000e+00, 2.000000000000000e+00, 9.861075762086680e-17,
      -2.000000000000000e+00, 2.000000000000000e+00, 9.861075762086680e-17,
      -2.000000000000000e+00, 2.000000000000000e+00, 9.861075762086680e-17,
      -2.775557561562891e-17, -2.000000000000000e+00, 4.287451790760826e-16,
      -2.775557561562891e-17, -2.000000000000000e+00, 4.287451790760826e-16,
      -2.775557561562891e-17, -2.000000000000000e+00, 4.287451790760826e-16,
      -2.775557561562891e-17, -2.000000000000000e+00, 4.287451790760826e-16,
      2.000000000000000e+00, -2.185751579730777e-16, 1.526556658859590e-16,
      2.000000000000000e+00, -2.185751579730777e-16, 1.526556658859590e-16,
      2.000000000000000e+00, -2.185751579730777e-16, 1.526556658859590e-16,
      2.000000000000000e+00, -2.185751579730777e-16, 1.526556658859590e-16
    };

    int cur=0;
    for (int i=0;i<myBasisValues.dimension(0);i++) {
      for (int j=0;j<myBasisValues.dimension(1);j++) {
        for (int k=0;k<myBasisValues.dimension(2);k++) {
          if (std::abs( myBasisValues(i,j,k) - fiat_curls[cur] ) > 10.0*INTREPID_TOL ) {
            errorFlag++;
            *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
            
            // Output the multi-index of the value where the error is:
            *outStream << " At multi-index { ";
            *outStream << i << " " << j << " " << k;
            *outStream << "}  computed value: " << myBasisValues(i,j,k)
                       << " but correct value: " << fiat_curls[cur] << "\n";
	    *outStream << "Difference: " << std::abs( myBasisValues(i,j,k) - fiat_curls[cur] ) << "\n";
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


  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";
  
  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  
  return errorFlag;
}
