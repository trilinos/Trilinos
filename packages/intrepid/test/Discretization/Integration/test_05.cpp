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
\brief  Unit test (CubatureDirect,CubatureTensor,DefaultCubatureFactory): correctness of
        integration of monomials for 3D reference cells.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace Intrepid;

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
  for (int i=0; i<p.dimension(0); i++) {
    val *= std::pow(p(i),polydeg[i]);
  }
  return val;
}


/*
  Computes integrals of monomials over a given reference cell.
*/
double computeIntegral(shards::CellTopology & cellTopology, int cubDegree, int xDeg, int yDeg, int zDeg) {

  DefaultCubatureFactory<double>  cubFactory;                                         // create factory
  Teuchos::RCP<Cubature<double> > myCub = cubFactory.create(cellTopology, cubDegree); // create default cubature

  double val       = 0.0;
  int cubDim       = myCub->getDimension();
  int numCubPoints = myCub->getNumPoints();

  FieldContainer<double> point(cubDim);
  FieldContainer<double> cubPoints(numCubPoints, cubDim);
  FieldContainer<double> cubWeights(numCubPoints);
  FieldContainer<double> functValues(numCubPoints);

  myCub->getCubature(cubPoints, cubWeights);

  for (int i=0; i<numCubPoints; i++) {
    for (int j=0; j<cubDim; j++) {
      point(j) = cubPoints(i,j);
    }
    functValues(i) = computeMonomial(point, xDeg, yDeg, zDeg);
  }

  Teuchos::BLAS<int, double> myblas;
  int inc = 1;
  val = myblas.DOT(numCubPoints, &functValues[0], inc, &cubWeights[0], inc);

  return val;
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
  << "|      Unit Test (CubatureDirect,CubatureTensor,DefaultCubatureFactory)       |\n" \
  << "|                                                                             |\n" \
  << "|     1) Computing integrals of monomials on reference cells in 3D            |\n" \
  << "|                         - using Level 1 BLAS -                              |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| TEST 1: integrals of monomials in 3D (Level 1 BLAS version)                 |\n"\
  << "===============================================================================\n";

  // internal variables:
  int                                      errorFlag = 0;
  int                                      polyCt = 0;
  int                                      offset = 0;
  Teuchos::Array< Teuchos::Array<double> > testInt;
  Teuchos::Array< Teuchos::Array<double> > analyticInt;
  Teuchos::Array<double>                   tmparray(1);
  double                                   reltol = 1.0e+04 * INTREPID_TOL;
  int                                      maxDeg[3];
  int                                      maxOffset[3];
  int                                      numPoly[3];
  int                                      numAnalytic[3];
  // max polynomial degree tested, per cell type:
  maxDeg[0]                              = INTREPID_CUBATURE_TET_DEFAULT_MAX;
  maxDeg[1]                              = 20; // can be as large as INTREPID_CUBATURE_LINE_GAUSS_MAX, but runtime is excessive
  maxDeg[2]                              = INTREPID_CUBATURE_TRI_DEFAULT_MAX;
  // max polynomial degree recorded in analytic comparison files, per cell type:
  maxOffset[0]                           = INTREPID_CUBATURE_TET_DEFAULT_MAX;
  maxOffset[1]                           = INTREPID_CUBATURE_LINE_GAUSS_MAX;
  maxOffset[2]                           = INTREPID_CUBATURE_TRI_DEFAULT_MAX;
  for (int i=0; i<3; i++) {
    numPoly[i] = (maxDeg[i]+1)*(maxDeg[i]+2)*(maxDeg[i]+3)/6;
  }
  for (int i=0; i<3; i++) {
    numAnalytic[i] = (maxOffset[i]+1)*(maxOffset[i]+2)*(maxOffset[i]+3)/6;
  }

  // get names of files with analytic values
  std::string basedir = "./data";
  std::stringstream namestream[3];
  std::string filename[3];
  namestream[0] << basedir << "/TET_integrals" << ".dat";
  namestream[0] >> filename[0];
  namestream[1] << basedir << "/HEX_integrals" << ".dat";
  namestream[1] >> filename[1];
  namestream[2] << basedir << "/TRIPRISM_integrals" << ".dat";
  namestream[2] >> filename[2];

  // reference cells tested
  shards::CellTopology cellType[] = {shards::getCellTopologyData< shards::Tetrahedron<> >(),
                                     shards::getCellTopologyData< shards::Hexahedron<> >(),
                                     shards::getCellTopologyData< shards::Wedge<> >()};
  // format of data files with analytic values
  TypeOfExactData dataFormat[] = {INTREPID_UTILS_SCALAR, INTREPID_UTILS_FRACTION, INTREPID_UTILS_FRACTION};

  // compute and compare integrals
  try {
    for (int cellCt=0; cellCt < 3; cellCt++) {
      testInt.assign(numPoly[cellCt], tmparray);
      analyticInt.assign(numAnalytic[cellCt], tmparray);

      *outStream << "\nIntegrals of monomials on a reference " << cellType[cellCt].getBaseCellTopologyData()->name << ":\n";
      std::ifstream filecompare(&filename[cellCt][0]);
      // compute integrals
      for (int cubDeg=0; cubDeg <= maxDeg[cellCt]; cubDeg++) {
        polyCt = 0;
        testInt[cubDeg].resize((cubDeg+1)*(cubDeg+2)*(cubDeg+3)/6);
        for (int xDeg=0; xDeg <= cubDeg; xDeg++) {
          for (int yDeg=0; yDeg <= cubDeg-xDeg; yDeg++) {
            for (int zDeg=0; zDeg <= cubDeg-xDeg-yDeg; zDeg++) {
              testInt[cubDeg][polyCt] = computeIntegral(cellType[cellCt], cubDeg, xDeg, yDeg, zDeg);
              polyCt++; 
            }
          }
        }
      }
      // get analytic values
      if (filecompare.is_open()) {
        getAnalytic(analyticInt, filecompare, dataFormat[cellCt]);
        // close file
        filecompare.close();
      }
      // perform comparison
      for (int cubDeg=0; cubDeg <= maxDeg[cellCt]; cubDeg++) {
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
            offset = offset + maxOffset[cellCt] - cubDeg;
          }
          offset = offset + (maxOffset[cellCt] - cubDeg)*(maxOffset[cellCt] - cubDeg + 1)/2;
        }
        *outStream << "Cubature order " << std::setw(2) << std::left << cubDeg;
        if (errorFlag == oldErrorFlag)
         *outStream << " passed.\n";
        else
         *outStream << " failed.\n";
      }
      *outStream << "\n";
    }  // end for cellCt
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
