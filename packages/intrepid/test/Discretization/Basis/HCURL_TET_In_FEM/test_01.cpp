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
1.000000000000001e+00, -2.498001805406602e-16, -1.665334536937735e-16,
9.999999999999998e-01, 1.000000000000000e+00, 1.000000000000000e+00,
5.828670879282072e-16, 1.110223024625157e-16, 2.498001805406602e-16,
7.771561172376096e-16, 8.326672684688674e-17, 1.110223024625157e-16,
2.081668171172169e-16, -2.914335439641036e-16, 1.280865063236792e-16,
-3.191891195797325e-16, 1.000000000000000e+00, -4.293998586504916e-17,
-9.999999999999994e-01, 2.081668171172169e-16, 2.400576428367544e-16,
2.220446049250313e-16, -5.551115123125783e-17, 1.084013877651281e-16,
3.469446951953614e-16, -1.000000000000000e+00, 1.387778780781446e-16,
-1.804112415015879e-16, 1.942890293094024e-16, -1.387778780781446e-16,
-9.999999999999993e-01, -9.999999999999996e-01, -9.999999999999998e-01,
5.551115123125783e-17, -2.220446049250313e-16, -8.326672684688674e-17,
-2.220446049250313e-16, -5.551115123125783e-17, 9.999999999999999e-01,
1.665334536937735e-16, 1.110223024625157e-16, -6.383782391594650e-16,
1.110223024625157e-16, 1.110223024625157e-16, -1.110223024625157e-16,
9.999999999999990e-01, 9.999999999999994e-01, 9.999999999999996e-01,
1.387778780781446e-16, -2.496931404305374e-17, -1.665334536937735e-16,
-2.498001805406602e-16, -2.149987498083074e-16, 1.000000000000000e+00,
8.326672684688674e-17, -3.769887250591415e-17, 8.326672684688674e-17,
-9.999999999999994e-01, 1.556977698723022e-16, 2.220446049250313e-16,
-9.422703950001342e-18, 1.665334536937735e-16, -2.359223927328458e-16,
-9.422703950001268e-18, -8.326672684688674e-17, 1.387778780781446e-17,
-7.525083148581445e-17, 2.775557561562891e-17, 1.000000000000000e+00,
2.789513560273035e-16, -9.999999999999998e-01, -5.551115123125783e-17
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
