// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file test_07.cpp
\brief  Unit test (CubatureSparse): correctness of
        integration of monomials for 2D reference cells.
\author Created by P. Bochev, D. Ridzal, and Matthew Keegan.
*/

#include "Intrepid2_CubatureSparse.hpp"
#include "Intrepid2_Utils.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace Intrepid2;

//Teuchos::RCP<std::ostream> outStream;

/*
  Monomial evaluation.
    in 1D, for point p(x)    : x^xDeg
    in 2D, for point p(x,y)  : x^xDeg * y^yDeg
    in 3D, for point p(x,y,z): x^xDeg * y^yDeg * z^zDeg
*/
double computeMonomial(FieldContainer<double> & p, int xDeg, int yDeg=0, int zDeg=0) {
  double val = 1.0;
  int polydeg[3];
  polydeg[0] = xDeg; polydeg[1] = yDeg; polydeg[2] = zDeg;
  for (int i=0; i<p.extent(0); i++) {
    val *= std::pow(p(i),polydeg[i]);
  }
  return val;
}


/*
  Computes integrals of monomials over a given reference cell.
*/
double computeIntegral(int cubDegree, int xDeg, int yDeg) {

  CubatureSparse<double,2> myCub(cubDegree);

  double val       = 0.0;
  int cubDim       = myCub.getDimension();
  int numCubPoints = myCub.getNumPoints();

  FieldContainer<double> point(cubDim);
  FieldContainer<double> cubPoints(numCubPoints, cubDim);
  FieldContainer<double> cubWeights(numCubPoints);

  myCub.getCubature(cubPoints, cubWeights);

  for (int i=0; i<numCubPoints; i++) {
    for (int j=0; j<cubDim; j++) {
      point(j) = cubPoints(i,j);
    }
    val += computeMonomial(point, xDeg, yDeg)*cubWeights(i);
  }

  return val;
}


int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
Kokkos::initialize();
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
  << "|                         Unit Test (CubatureSparse)                          |\n" \
  << "|                                                                             |\n" \
  << "|     1) Computing integrals of monomials on reference cells in 2D            |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov),                     |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" \
  << "|                      Matthew Keegan (mskeega@sandia.gov)                    |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| TEST 1: integrals of monomials in 2D for Sparse Grid Construction           |\n"\
  << "===============================================================================\n";

  // internal variables:
  int                                      errorFlag = 0;
  int                                      polyCt = 0;
  int                                      offset = 0;
  Teuchos::Array< Teuchos::Array<double> > testInt;
  Teuchos::Array< Teuchos::Array<double> > analyticInt;
  Teuchos::Array<double>                   tmparray(1);
  double                                   reltol = 1.0e+03 * INTREPID_TOL;
  int maxDeg                             = 30; // can be as large as INTREPID2_CUBATURE_SPARSE2D_GAUSS_MAX, but runtime is excessive
  int maxOffset                          = INTREPID2_CUBATURE_LINE_GAUSS_MAX;
  int numPoly                            = (maxDeg+1)*(maxDeg+2)/2;
  int numAnalytic                        = (maxOffset+1)*(maxOffset+2)/2;
  testInt.assign(numPoly, tmparray);
  analyticInt.assign(numAnalytic, tmparray);

  // get names of files with analytic values
  std::string basedir = "./data";
  std::stringstream namestream;
  std::string filename;
  namestream << basedir << "/QUAD_integrals" << ".dat";
  namestream >> filename;

  // compute and compare integrals
  try {
      *outStream << "\nIntegrals of monomials:\n";

      std::ifstream filecompare(&filename[0]);
      // compute integrals
      for (int cubDeg=0; cubDeg <= maxDeg; cubDeg++) {
        polyCt = 0;
        testInt[cubDeg].resize((cubDeg+1)*(cubDeg+2)/2);
        for (int xDeg=0; xDeg <= cubDeg; xDeg++) {
          for (int yDeg=0; yDeg <= cubDeg-xDeg; yDeg++) {
            testInt[cubDeg][polyCt] = computeIntegral(cubDeg, xDeg, yDeg);
            polyCt++; 
          }
        }
      }

      // get analytic values
      if (filecompare.is_open()) {
        getAnalytic(analyticInt, filecompare);
        // close file
        filecompare.close();
      }
      // perform comparison
      for (int cubDeg=0; cubDeg <= maxDeg; cubDeg++) {
        polyCt = 0;
        offset = 0;
        for (int xDeg=0; xDeg <= cubDeg; xDeg++) {
          for (int yDeg=0; yDeg <= cubDeg-xDeg; yDeg++) {
            double abstol = ( analyticInt[polyCt+offset][0] == 0.0 ? reltol : std::fabs(reltol*analyticInt[polyCt+offset][0]) );
            double absdiff = std::fabs(analyticInt[polyCt+offset][0] - testInt[cubDeg][polyCt]);
            *outStream << "Cubature order " << std::setw(2) << std::left << cubDeg << " integrating "
                       << "x^" << std::setw(2) << std::left << xDeg << " * y^" << std::setw(2) << yDeg << ":" << "   "
                       << std::scientific << std::setprecision(16) << testInt[cubDeg][polyCt] << "   " << analyticInt[polyCt+offset][0] << "   "
                       << std::setprecision(4) << absdiff << "   " << "<?" << "   " << abstol << "\n";
            if (absdiff > abstol) {
              errorFlag++;
              *outStream << std::right << std::setw(111) << "^^^^---FAILURE!\n";
            }
            polyCt++;
          }
          offset = offset + maxOffset - cubDeg;
        }
        *outStream << "\n";
      }
      *outStream << "\n";
  }
  catch (std::logic_error &err) {
    *outStream << err.what() << "\n";
    errorFlag = -1;
  };


  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
Kokkos::finalize();
  return errorFlag;
}
