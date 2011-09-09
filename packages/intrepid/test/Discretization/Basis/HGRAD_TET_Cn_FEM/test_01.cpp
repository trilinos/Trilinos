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
\brief  Unit test of Lagrange basis class
\author Created by R. Kirby
*/

#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Intrepid_PointTools.hpp"
#include "Intrepid_HGRAD_TET_Cn_FEM.hpp"
#include "Shards_CellTopology.hpp"

#include <iostream>
using namespace Intrepid;


/** \brief Tests for Lagrange basis on tets.  Tests Kronecker property of basis and basic execution
          of differentiation and dof-tab lookup
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
    << "|                           Unit Test HGRAD_TET_Cn_FEM                        |\n" \
    << "|                                                                             |\n" \
    << "|     1) Tests triangular Lagrange basis                                      |\n" \
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

  // Let's instantiate a basis
  try {
    const int deg = 3;
    Basis_HGRAD_TET_Cn_FEM<double,FieldContainer<double> >  myBasis( deg , POINTTYPE_WARPBLEND );

    // Get a lattice
    const int np_lattice = PointTools::getLatticeSize( myBasis.getBaseCellTopology() , deg , 0 );
    const int nbf = myBasis.getCardinality();
    FieldContainer<double> lattice( np_lattice , 3 );
    PointTools::getLattice<double,FieldContainer<double> >( lattice , 
                                                            myBasis.getBaseCellTopology() , 
                                                            deg , 
                                                            0 , 
                                                            POINTTYPE_WARPBLEND );         
    FieldContainer<double> vals( nbf , np_lattice );

    myBasis.getValues( vals , lattice , OPERATOR_VALUE );

    // test for Kronecker property
    for (int i=0;i<nbf;i++) {
      for (int j=0;j<np_lattice;j++) {
        if ( i==j && std::abs( vals(i,j) - 1.0 ) > 100.0 * INTREPID_TOL ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Basis function " << i << " does not have unit value at its node\n";
        }
        if ( i!=j && std::abs( vals(i,j) ) > 100.0 * INTREPID_TOL ) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << " Basis function " << i << " does not vanish at node " << j << "\n";
          *outStream << " Basis function value is " << vals(i,j) << "\n";
        }
      }
    }
  }
  catch (std::exception err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }

  try {
    const int deg = 4;
    Basis_HGRAD_TET_Cn_FEM<double,FieldContainer<double> >  myBasis( deg , POINTTYPE_WARPBLEND );
    const std::vector<std::vector<std::vector<int> > >& dofdata = myBasis.getDofOrdinalData();
    for (unsigned d=0;d<dofdata.size();d++) {
      std::cout << "Dimension " << d << "\n";
      for (unsigned f=0;f<dofdata[d].size();f++) {
        std::cout << "\tFacet number " << f << "\n";
        std::cout << "\t\tDOFS:\n";
        for (unsigned n=0;n<dofdata[d][f].size();n++) {
          if ( dofdata[d][f][n] >= 0 ) {
            std::cout << "\t\t\t" << dofdata[d][f][n] << "\n";
          }
        }
      }
    }
  }
  catch (std::exception err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
  }

  try {
    const int deg = 3;
    Basis_HGRAD_TET_Cn_FEM<double,FieldContainer<double> >  myBasis( deg , POINTTYPE_EQUISPACED );

    // Get a lattice
    const int np_lattice = PointTools::getLatticeSize( myBasis.getBaseCellTopology() , deg , 0 );
    const int nbf = myBasis.getCardinality();
    FieldContainer<double> lattice( np_lattice , 3 );
    PointTools::getLattice<double,FieldContainer<double> >( lattice , 
                                                            myBasis.getBaseCellTopology() , 
                                                            deg , 
                                                            0 , 
                                                            POINTTYPE_EQUISPACED );         
    FieldContainer<double> vals( nbf , np_lattice , 3 );

    myBasis.getValues( vals , lattice , OPERATOR_GRAD );

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
