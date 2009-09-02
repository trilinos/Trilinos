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
