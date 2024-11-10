// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
\brief  Unit test (CubatureGenSparse): correctness of
        integration of monomials for 3D reference cells.
\author Created by P. Bochev, D. Ridzal and M. Keegan
*/

#include "Intrepid2_CubatureGenSparse.hpp"
#include "Intrepid2_Utils.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace Intrepid2;


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
void computeIntegral(Teuchos::Array<double>& testIntFixDeg, int cubDegree) {

  CubatureGenSparse<double,3> myCub(cubDegree);

  int cubDim       = myCub.getDimension();
  int numCubPoints = myCub.getNumPoints();
  int numPolys     = (cubDegree+1)*(cubDegree+2)*(cubDegree+3)/6;

  FieldContainer<double> point(cubDim);
  FieldContainer<double> cubPoints(numCubPoints, cubDim);
  FieldContainer<double> cubWeights(numCubPoints);
  FieldContainer<double> functValues(numCubPoints, numPolys);

  myCub.getCubature(cubPoints, cubWeights);

  int polyCt = 0;
  for (int xDeg=0; xDeg <= cubDegree; xDeg++) {
    for (int yDeg=0; yDeg <= cubDegree-xDeg; yDeg++) {
      for (int zDeg=0; zDeg <= cubDegree-xDeg-yDeg; zDeg++) {
        for (int i=0; i<numCubPoints; i++) {
          for (int j=0; j<cubDim; j++) {
            point(j) = cubPoints(i,j);
          }
          functValues(i,polyCt) = computeMonomial(point, xDeg, yDeg, zDeg);
        }
        polyCt++;
      }
    }
  }

  Teuchos::BLAS<int, double> myblas;
  int inc = 1;
  double alpha = 1.0;
  double beta  = 0.0;
  myblas.GEMV(Teuchos::NO_TRANS, numPolys, numCubPoints, alpha, &functValues(0,0), numPolys,
              &cubWeights(0), inc, beta, &testIntFixDeg[0], inc);
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
  << "|                      Unit Test (CubatureGenSparse)                          |\n" \
  << "|                                                                             |\n" \
  << "|     1) Computing integrals of monomials on reference cells in 3D            |\n" \
  << "|                         - using Level 2 BLAS -                              |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov),                     |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" \
  << "|                      Matthew Keegan (mskeega@sandia.gov).                   |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| TEST 1: integrals of monomials in 3D (Level 2 BLAS version) using           |\n"\
  << "|           Generalized Sparse Grid Construction                              |\n"\
  << "===============================================================================\n";

  // internal variables:
  int                                      errorFlag = 0;
  int                                      polyCt = 0;
  int                                      offset = 0;
  Teuchos::Array< Teuchos::Array<double> > testInt;
  Teuchos::Array< Teuchos::Array<double> > analyticInt;
  Teuchos::Array<double>                   tmparray(1);
  double                                   reltol = 1.0e+04 * INTREPID_TOL;
  int maxDeg                             = 20; // can be as large as INTREPID2_CUBATURE_SPARSE3D_GAUSS_MAX, but runtime is excessive
  int maxOffset                          = INTREPID2_CUBATURE_LINE_GAUSS_MAX;
  int numPoly                            = (maxDeg+1)*(maxDeg+2)*(maxDeg+3)/6;
  int numAnalytic                        = (maxOffset+1)*(maxOffset+2)*(maxOffset+3)/6;
  testInt.assign(numPoly, tmparray);
  analyticInt.assign(numAnalytic, tmparray);

  // get names of files with analytic values
  std::string basedir = "./data";
  std::stringstream namestream;
  std::string filename;
  namestream << basedir << "/HEX_integrals" << ".dat";
  namestream >> filename;

  // format of data files with analytic values
  TypeOfExactData dataFormat = INTREPID_UTILS_FRACTION;

  // compute and compare integrals
  try {
      *outStream << "\nIntegrals of monomials:\n";
      std::ifstream filecompare(&filename[0]);
      // compute integrals
      for (int cubDeg=0; cubDeg <= maxDeg; cubDeg++) {
        int numMonomials = (cubDeg+1)*(cubDeg+2)*(cubDeg+3)/6; 
        testInt[cubDeg].resize(numMonomials);
        computeIntegral(testInt[cubDeg], cubDeg);
      }
      // get analytic values
      if (filecompare.is_open()) {
        getAnalytic(analyticInt, filecompare, dataFormat);
        // close file
        filecompare.close();
      }
      // perform comparison
      for (int cubDeg=0; cubDeg <= maxDeg; cubDeg++) {
        polyCt = 0;
        offset = 0;
        int oldErrorFlag = errorFlag;
        for (int xDeg=0; xDeg <= cubDeg; xDeg++) {
          for (int yDeg=0; yDeg <= cubDeg-xDeg; yDeg++) {
            for (int zDeg=0; zDeg <= cubDeg-xDeg-yDeg; zDeg++) {
              double abstol = ( analyticInt[polyCt+offset][0] == 0.0 ? reltol : std::fabs(reltol*analyticInt[polyCt+offset][0]) );
              double absdiff = std::fabs(analyticInt[polyCt+offset][0] - testInt[cubDeg][polyCt]);
              if (absdiff > abstol) {
                *outStream << "Cubature order " << std::setw(2) << std::left << cubDeg << " integrating "
                           << "x^" << std::setw(2) << std::left << xDeg << " * y^" << std::setw(2) << yDeg
                           << " * z^" << std::setw(2) << zDeg << ":" << "   "
                           << std::scientific << std::setprecision(16)
                           << testInt[cubDeg][polyCt] << "   " << analyticInt[polyCt+offset][0] << "   "
                           << std::setprecision(4) << absdiff << "   " << "<?" << "   " << abstol << "\n";
                errorFlag++;
                *outStream << std::right << std::setw(118) << "^^^^---FAILURE!\n";
              }
              polyCt++;
            }
            offset = offset + maxOffset - cubDeg;
          }
          offset = offset + (maxOffset - cubDeg)*(maxOffset - cubDeg + 1)/2;
        }
        *outStream << "Cubature order " << std::setw(2) << std::left << cubDeg;
        if (errorFlag == oldErrorFlag)
         *outStream << " passed.\n";
        else
         *outStream << " failed.\n";
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
