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
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file test_01.cpp
\brief  Unit test for the PointTools class.
\author Created by R. Kirby.
*/

#include "Intrepid_PointTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Shards_CellTopology.hpp"

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
  << "|                       Unit Test (PointTools)                                |\n" \
  << "|                                                                             |\n" \
  << "|     1) Construction of equispaced and warped lattices on simplices          |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" \
  << "|                      Robert Kirby (robert.c.kirby@ttu.edu)                  |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n";



  int errorFlag = 0;
#ifdef HAVE_INTREPID_DEBUG
  int beginThrowNumber = TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + 50;
#endif

  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 1: size of lattices                                                   |\n"\
  << "===============================================================================\n";

  try{
    // first try the lattices with offset = 0.  This is a spot-check

    shards::CellTopology myLine_2( shards::getCellTopologyData< shards::Line<2> >() );
    shards::CellTopology myTri_3( shards::getCellTopologyData< shards::Triangle<3> >() );
    shards::CellTopology myTet_4( shards::getCellTopologyData< shards::Tetrahedron<4> >() );

    if (PointTools::getLatticeSize( myLine_2 , 4 , 0 ) != 5) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " size of 4th order lattice on a line with no offset: " << PointTools::getLatticeSize( myLine_2 , 4 , 0 )  << "\n";
      *outStream << " should be 5\n";
    }


    if (PointTools::getLatticeSize( myTri_3 , 3 , 0 ) != 10) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " size of 3rd order lattice on a line with no offset: " << PointTools::getLatticeSize( myTri_3 , 3 , 0 )  << "\n";
      *outStream << " should be 10\n";    
    }


    if (PointTools::getLatticeSize( myTet_4 , 3 , 0 ) != 20) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " size of 3rd order lattice on a tet with no offset: " << PointTools::getLatticeSize( myTet_4 , 3 , 0 )  << "\n";
      *outStream << " should be 20\n"; 
    }

			
    // check with the offset = 1
    if (PointTools::getLatticeSize( myLine_2 , 3 , 1 ) != 2) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " size of 3rd order lattice on a line with 1 offset: " << PointTools::getLatticeSize( myLine_2 , 3 , 1 )  << "\n";
      *outStream << " should be 2\n";       
    }

    if (PointTools::getLatticeSize( myTri_3 , 4 , 1 ) != 3) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " size of 4th order lattice on a triangle with 1 offset: " << PointTools::getLatticeSize( myTri_3 , 4 , 1 )  << "\n";
      *outStream << " should be 3\n";           
    }

    if (PointTools::getLatticeSize( myTet_4 , 5 , 1 ) != 4) {
      errorFlag++;
      *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      *outStream << " size of 5th order lattice on a tetrahedron with 1 offset: " << PointTools::getLatticeSize( myTet_4 , 5 , 1 )  << "\n";
      *outStream << " should be 4\n";   
    }

  }
  catch (std::exception &err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

  // Now verify that we throw an exception on some of the non-supported cell types.

  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 2: check for unsupported cell types                                     \n"\
  << "===============================================================================\n";
  try{
    try {
      PointTools::getLatticeSize( shards::getCellTopologyData< shards::Quadrilateral<4> >() , 3 , 0 );
    }
    catch (std::invalid_argument err) {
      *outStream << err.what() << "\n";
    }

    try {
      PointTools::getLatticeSize( shards::getCellTopologyData< shards::Hexahedron<8> >() , 3 , 0 );
    }
    catch (std::invalid_argument err) {
      *outStream << err.what() << "\n";
    }

  }
  catch (std::exception &err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
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
