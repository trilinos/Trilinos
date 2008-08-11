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


/** \file
\brief  Unit test (CubatureDirect,CubatureTensor): correctness of volume
        computations for reference cells.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_CubatureDirect.hpp"
#include "Intrepid_CubatureTensor.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

using namespace Intrepid;

/*
  Computes volume of a given reference cell.
*/
double computeRefVolume(ECell cellType, int cubDegree) {

  Teuchos::RCP< Cubature<double> > myCub;
  double vol = 0.0;

  switch (cellType) {

    case CELL_EDGE:
    case CELL_TRI:
    case CELL_TET:
        myCub = Teuchos::rcp(new CubatureDirect<double>(cellType, cubDegree));
      break;

    case CELL_QUAD:
    case CELL_HEX:
    case CELL_TRIPRISM:
        myCub = Teuchos::rcp(new CubatureTensor<double>(cellType, cubDegree));
      break;

    default:
      TEST_FOR_EXCEPTION((cellType != CELL_EDGE) && (cellType != CELL_TRI) && (cellType != CELL_TET) &&
                         (cellType != CELL_QUAD) && (cellType != CELL_HEX) && (cellType != CELL_TRIPRISM),
                          std::invalid_argument,
                          ">>> ERROR (Unit Test -- Cubature -- Volume): Invalid cell type.");
  } // end switch

  int numCubPoints = 0;

  Teuchos::Array< Point<double> > cubPoints;
  Teuchos::Array<double> cubWeights;

  myCub->getCubature(numCubPoints, cubPoints, cubWeights);

  for (int i=0; i<numCubPoints; i++)
    vol += cubWeights[i];

  return vol;
}


int main(int argc, char *argv[]) {

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
  << "|                  Unit Test (CubatureDirect,CubatureTensor)                  |\n" \
  << "|                                                                             |\n" \
  << "|     1) Computing volumes of reference cells                                 |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| TEST 1: volume computations                                                 |\n"\
  << "===============================================================================\n";

  int errorFlag  = 0;
  double testVol = 0.0;
  double tol     = 100.0 * INTREPID_TOL;

  // list of analytic volume values, listed in the enumerated reference cell order up to CELL_HEXAPRISM
  double volumeList[] = {0.0, 2.0, 1.0/2.0, 4.0, 1.0/6.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0};

  *outStream << "\nReference cell volumes:\n";

  try {
    for (int deg=0; deg<INTREPID_MAX_CUBATURE_DEGREE_EDGE; deg++) {
      testVol = computeRefVolume(CELL_EDGE, deg);
      *outStream << std::setw(24) << "EDGE volume --> " << std::setw(10) << std::scientific << testVol <<
                    std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - volumeList[CELL_EDGE]) << "\n";
      if (std::abs(testVol - volumeList[CELL_EDGE]) > tol) {
        errorFlag = 1;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      }
    }
    *outStream << std::endl;
    for (int deg=0; deg<INTREPID_MAX_CUBATURE_DEGREE_EDGE; deg++) {
      testVol = computeRefVolume(CELL_QUAD, deg);
      *outStream << std::setw(24) << "QUAD volume --> " << std::setw(10) << std::scientific << testVol <<
                    std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - volumeList[CELL_QUAD]) << "\n";
      if (std::abs(testVol - volumeList[CELL_QUAD]) > tol) {
        errorFlag = 1;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      }
    }
    *outStream << std::endl;
    for (int deg=0; deg<INTREPID_MAX_CUBATURE_DEGREE_TRI; deg++) {
      testVol = computeRefVolume(CELL_TRI, deg);
      *outStream << std::setw(24) << "TRI volume --> " << std::setw(10) << std::scientific << testVol <<
                    std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - volumeList[CELL_TRI]) << "\n";
      if (std::abs(testVol - volumeList[CELL_TRI]) > tol) {
        errorFlag = 1;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      }
    }
    *outStream << std::endl;
    for (int deg=0; deg<INTREPID_MAX_CUBATURE_DEGREE_TET; deg++) {
      testVol = computeRefVolume(CELL_TET, deg);
      *outStream << std::setw(24) << "TET volume --> " << std::setw(10) << std::scientific << testVol <<
                    std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - volumeList[CELL_TET]) << "\n";
      if (std::abs(testVol - volumeList[CELL_TET]) > tol) {
        errorFlag = 1;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      }
    }
    *outStream << std::endl;
    for (int deg=0; deg<INTREPID_MAX_CUBATURE_DEGREE_EDGE; deg++) {
      testVol = computeRefVolume(CELL_HEX, deg);
      *outStream << std::setw(24) << "HEX volume --> " << std::setw(10) << std::scientific << testVol <<
                    std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - volumeList[CELL_HEX]) << "\n";
      if (std::abs(testVol - volumeList[CELL_HEX]) > tol) {
        errorFlag = 1;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      }
    }
    *outStream << std::endl;
    for (int deg=0; deg<std::min(INTREPID_MAX_CUBATURE_DEGREE_EDGE,INTREPID_MAX_CUBATURE_DEGREE_TRI); deg++) {
      testVol = computeRefVolume(CELL_TRIPRISM, deg);
      *outStream << std::setw(24) << "TRIPRISM volume --> " << std::setw(10) << std::scientific << testVol <<
                    std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - volumeList[CELL_TRIPRISM]) << "\n";
      if (std::abs(testVol - volumeList[CELL_TRIPRISM]) > tol) {
        errorFlag = 1;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      }
    }
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1;
  };


  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return errorFlag;
}
