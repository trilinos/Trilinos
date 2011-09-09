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
\brief  Unit test (CubatureDirect,CubatureTensor): correctness of volume
        computations for reference cells.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_CubatureDirectLineGauss.hpp"
#include "Intrepid_CubatureDirectTriDefault.hpp"
#include "Intrepid_CubatureDirectTetDefault.hpp"
#include "Intrepid_CubatureTensor.hpp"
#include "Shards_CellTopology.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace Intrepid;

#define INTREPID_TEST_COMMAND( S )                                                                                  \
{                                                                                                                   \
  try {                                                                                                             \
    S ;                                                                                                             \
  }                                                                                                                 \
  catch (std::logic_error err) {                                                                                    \
      *outStream << "Expected Error ----------------------------------------------------------------\n";            \
      *outStream << err.what() << '\n';                                                                             \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";    \
  };                                                                                                                \
}

/*
  Computes volume of a given reference cell.
*/
double computeRefVolume(shards::CellTopology & cellTopology, int cubDegree) {
  Teuchos::RCP< Cubature<double> > myCub;
  double vol = 0.0;

  switch (cellTopology.getBaseCellTopologyData()->key) {

    case shards::Line<>::key:
        myCub = Teuchos::rcp(new CubatureDirectLineGauss<double>(cubDegree));
      break;
    case shards::Triangle<>::key:
        myCub = Teuchos::rcp(new CubatureDirectTriDefault<double>(cubDegree));
      break;
    case shards::Tetrahedron<>::key:
        myCub = Teuchos::rcp(new CubatureDirectTetDefault<double>(cubDegree));
      break;
    case shards::Quadrilateral<>::key: { 
        std::vector< Teuchos::RCP< Cubature<double> > > lineCubs(2);
        lineCubs[0] = Teuchos::rcp(new CubatureDirectLineGauss<double>(cubDegree));
        lineCubs[1] = Teuchos::rcp(new CubatureDirectLineGauss<double>((cubDegree+7) % INTREPID_CUBATURE_LINE_GAUSS_MAX));
        myCub = Teuchos::rcp(new CubatureTensor<double>(lineCubs));
        }
      break;
    case shards::Hexahedron<>::key: { 
        std::vector< Teuchos::RCP< Cubature<double> > > lineCubs(2);
        std::vector< Teuchos::RCP< Cubature<double> > > miscCubs(2);
        lineCubs[0] = Teuchos::rcp(new CubatureDirectLineGauss<double>(cubDegree));
        lineCubs[1] = Teuchos::rcp(new CubatureDirectLineGauss<double>((cubDegree+7) % INTREPID_CUBATURE_LINE_GAUSS_MAX));
        miscCubs[0] = Teuchos::rcp(new CubatureTensor<double>(lineCubs));
        miscCubs[1] = Teuchos::rcp(new CubatureDirectLineGauss<double>((cubDegree+25) % INTREPID_CUBATURE_LINE_GAUSS_MAX));
        myCub = Teuchos::rcp(new CubatureTensor<double>(miscCubs));
        }
      break;
    case shards::Wedge<>::key: {
        Teuchos::RCP<CubatureDirect<double> > triCub  = Teuchos::rcp(new CubatureDirectTriDefault<double>(cubDegree));
        Teuchos::RCP<CubatureDirect<double> > lineCub = Teuchos::rcp(new CubatureDirectLineGauss<double>(cubDegree));
        myCub = Teuchos::rcp(new CubatureTensor<double>(triCub,lineCub));
        }
      break;

    default:
      TEST_FOR_EXCEPTION( ( (cellTopology.getBaseCellTopologyData()->key != shards::Line<>::key),
                            (cellTopology.getBaseCellTopologyData()->key != shards::Triangle<>::key),
                            (cellTopology.getBaseCellTopologyData()->key != shards::Tetrahedron<>::key),
                            (cellTopology.getBaseCellTopologyData()->key != shards::Quadrilateral<>::key),
                            (cellTopology.getBaseCellTopologyData()->key != shards::Hexahedron<>::key),
                            (cellTopology.getBaseCellTopologyData()->key != shards::Wedge<>::key) ),
                          std::invalid_argument,
                          ">>> ERROR (Unit Test -- Cubature -- Volume): Invalid cell type.");
  } // end switch

  int numCubPoints = myCub->getNumPoints();

  FieldContainer<double> cubPoints( numCubPoints, myCub->getDimension() );
  FieldContainer<double> cubWeights( numCubPoints );
  myCub->getCubature(cubPoints, cubWeights);

  for (int i=0; i<numCubPoints; i++)
    vol += cubWeights[i];

  return vol;
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
  << "| TEST 1: construction and basic functionality                                |\n"\
  << "===============================================================================\n";

  int errorFlag  = 0;

  int beginThrowNumber = TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + 7;  

  try {
    /* Line cubature. */
    INTREPID_TEST_COMMAND( CubatureDirectLineGauss<double> lineCub(-1) );
    INTREPID_TEST_COMMAND( CubatureDirectLineGauss<double> lineCub(INTREPID_CUBATURE_LINE_GAUSS_MAX+1) );
    INTREPID_TEST_COMMAND( CubatureDirectLineGauss<double> lineCub;
                           std::string testName    = "INTREPID_CUBATURE_LINE_GAUSS";
                           std::string lineCubName = lineCub.getName();
                           *outStream << "\nComparing strings: " << testName << " and " << lineCubName << "\n\n";
                           TEST_FOR_EXCEPTION( (testName != lineCubName), std::logic_error, "Name mismatch!" ) );
    INTREPID_TEST_COMMAND( CubatureDirectLineGauss<double> lineCub;
                           std::vector<int> accuracy;
                           lineCub.getAccuracy(accuracy);
                           TEST_FOR_EXCEPTION( (accuracy[0] != 0), std::logic_error, "Check member getAccuracy!" ) );
    INTREPID_TEST_COMMAND( CubatureDirectLineGauss<double> lineCub(55);
                           TEST_FOR_EXCEPTION( (lineCub.getNumPoints() != 28), std::logic_error, "Check member getNumPoints!" ) );
    INTREPID_TEST_COMMAND( CubatureDirectLineGauss<double> lineCub;
                           TEST_FOR_EXCEPTION( (lineCub.getDimension() != 1),
                                               std::logic_error,
                                               "Check member dimension!" ) );
    /* Triangle cubature. */
    INTREPID_TEST_COMMAND( CubatureDirectTriDefault<double> triCub(-1) );
    INTREPID_TEST_COMMAND( CubatureDirectTriDefault<double> triCub(INTREPID_CUBATURE_TRI_DEFAULT_MAX+1) );
    INTREPID_TEST_COMMAND( CubatureDirectTriDefault<double> triCub;
                           std::string testName    = "INTREPID_CUBATURE_TRI_DEFAULT";
                           std::string triCubName = triCub.getName();
                           *outStream << "\nComparing strings: " << testName << " and " << triCubName << "\n\n";
                           TEST_FOR_EXCEPTION( (testName != triCubName), std::logic_error, "Name mismatch!" ) );
    INTREPID_TEST_COMMAND( CubatureDirectTriDefault<double> triCub;
                           std::vector<int> accuracy;
                           triCub.getAccuracy(accuracy);
                           TEST_FOR_EXCEPTION( (accuracy[0] != 0), std::logic_error, "Check member getAccuracy!" ) );
    INTREPID_TEST_COMMAND( CubatureDirectTriDefault<double> triCub(17);
                           TEST_FOR_EXCEPTION( (triCub.getNumPoints() != 61), std::logic_error, "Check member getNumPoints!" ) );
    INTREPID_TEST_COMMAND( CubatureDirectTriDefault<double> triCub;
                           TEST_FOR_EXCEPTION( (triCub.getDimension() != 2),
                                               std::logic_error,
                                               "Check member dimension!" ) );
    /* Tetrahedron cubature. */
    INTREPID_TEST_COMMAND( CubatureDirectTetDefault<double> tetCub(-1) );
    INTREPID_TEST_COMMAND( CubatureDirectTetDefault<double> tetCub(INTREPID_CUBATURE_TET_DEFAULT_MAX+1) );
    INTREPID_TEST_COMMAND( CubatureDirectTetDefault<double> tetCub;
                           std::string testName    = "INTREPID_CUBATURE_TET_DEFAULT";
                           std::string tetCubName = tetCub.getName();
                           *outStream << "\nComparing strings: " << testName << " and " << tetCubName << "\n\n";
        std::vector< Teuchos::RCP< Cubature<double> > > lineCubs(2);
                           TEST_FOR_EXCEPTION( (testName != tetCubName), std::logic_error, "Name mismatch!" ) );
    INTREPID_TEST_COMMAND( CubatureDirectTetDefault<double> tetCub;
                           std::vector<int> accuracy;
                           tetCub.getAccuracy(accuracy);
                           TEST_FOR_EXCEPTION( (accuracy[0] != 0), std::logic_error, "Check member getAccuracy!" ) );
    INTREPID_TEST_COMMAND( CubatureDirectTetDefault<double> tetCub(17);
                           TEST_FOR_EXCEPTION( (tetCub.getNumPoints() != 495), std::logic_error, "Check member getNumPoints!" ) );
    INTREPID_TEST_COMMAND( CubatureDirectTetDefault<double> tetCub;
                           TEST_FOR_EXCEPTION( (tetCub.getDimension() != 3),
                                               std::logic_error,
                                               "Check member getCellTopology!" ) );
    /* Tensor cubature. */
    INTREPID_TEST_COMMAND( std::vector< Teuchos::RCP< Cubature<double> > > lineCubs(0);
                           CubatureTensor<double> quadCub(lineCubs) );
    INTREPID_TEST_COMMAND( std::vector< Teuchos::RCP< Cubature<double> > > miscCubs(3);
                           std::vector< Teuchos::RCP< Cubature<double> > > lineCubs(2);
                           lineCubs[0] = Teuchos::rcp(new CubatureDirectLineGauss<double>(3));
                           lineCubs[1] = Teuchos::rcp(new CubatureDirectLineGauss<double>(16));
                           miscCubs[0] = Teuchos::rcp(new CubatureTensor<double>(lineCubs));
                           miscCubs[1] = Teuchos::rcp(new CubatureDirectTriDefault<double>);
                           miscCubs[2] = Teuchos::rcp(new CubatureDirectTetDefault<double>(19));
                           CubatureTensor<double> tensorCub(miscCubs);
                           std::vector<int> a(4); a[0]=3; a[1]=16; a[2]=0; a[3]=19;
                           std::vector<int> atest(4);
                           tensorCub.getAccuracy(atest);
                           TEST_FOR_EXCEPTION( (a != atest), std::logic_error, "Check member getAccuracy!" ) );
    INTREPID_TEST_COMMAND( std::vector< Teuchos::RCP< Cubature<double> > > lineCubs(2);
                           lineCubs[0] = Teuchos::rcp(new CubatureDirectLineGauss<double>(15));
                           lineCubs[1] = Teuchos::rcp(new CubatureDirectLineGauss<double>(11));
                           CubatureTensor<double> tensorCub(lineCubs);
                           TEST_FOR_EXCEPTION( (tensorCub.getNumPoints() != 48), std::logic_error, "Check member getNumPoints!" ) );
    INTREPID_TEST_COMMAND( std::vector< Teuchos::RCP< Cubature<double> > > miscCubs(3);
                           miscCubs[0] = Teuchos::rcp(new CubatureDirectLineGauss<double>);
                           miscCubs[1] = Teuchos::rcp(new CubatureDirectTriDefault<double>);
                           miscCubs[2] = Teuchos::rcp(new CubatureDirectTetDefault<double>);
                           CubatureTensor<double> tensorCub(miscCubs);
                           TEST_FOR_EXCEPTION( (tensorCub.getDimension() != 6), std::logic_error, "Check member dimension!" ) );
    INTREPID_TEST_COMMAND( std::vector< Teuchos::RCP< Cubature<double> > > miscCubs(3);
                           miscCubs[0] = Teuchos::rcp(new CubatureDirectLineGauss<double>(3));
                           miscCubs[1] = Teuchos::rcp(new CubatureDirectLineGauss<double>(7));
                           miscCubs[2] = Teuchos::rcp(new CubatureDirectLineGauss<double>(5));
                           CubatureTensor<double> tensorCub(miscCubs);
                           FieldContainer<double> points(tensorCub.getNumPoints(), tensorCub.getDimension());
                           FieldContainer<double> weights(tensorCub.getNumPoints());
                           tensorCub.getCubature(points, weights)
                         )
    INTREPID_TEST_COMMAND( Teuchos::RCP< CubatureDirect<double> > lineCub = Teuchos::rcp(new CubatureDirectLineGauss<double>(15));
                           Teuchos::RCP< CubatureDirect<double> > triCub = Teuchos::rcp(new CubatureDirectTriDefault<double>(12));
                           CubatureTensor<double> tensorCub(lineCub, triCub);
                           std::vector<int> a(2); a[0] = 15; a[1] = 12;
                           std::vector<int> atest(2);
                           tensorCub.getAccuracy(atest);
                           TEST_FOR_EXCEPTION( (tensorCub.getDimension() != 3) || (a != atest),
                                               std::logic_error,
                                               "Check constructormembers dimension and getAccuracy!" ) );
    INTREPID_TEST_COMMAND( Teuchos::RCP< CubatureDirect<double> > lineCub = Teuchos::rcp(new CubatureDirectLineGauss<double>(15));
                           Teuchos::RCP< CubatureDirect<double> > triCub = Teuchos::rcp(new CubatureDirectTriDefault<double>(12));
                           CubatureTensor<double> tensorCub(triCub, lineCub, triCub);
                           std::vector<int> a(3); a[0] = 12; a[1] = 15; a[2] = 12;
                           std::vector<int> atest(3);
                           tensorCub.getAccuracy(atest);
                           TEST_FOR_EXCEPTION( (tensorCub.getDimension() != 5) || (a != atest),
                                               std::logic_error,
                                               "Check constructor and members dimension and getAccuracy!" ) );
    INTREPID_TEST_COMMAND( Teuchos::RCP< CubatureDirect<double> > triCub = Teuchos::rcp(new CubatureDirectTriDefault<double>(12));
                           CubatureTensor<double> tensorCub(triCub, 5);
                           std::vector<int> a(5); a[0] = 12; a[1] = 12; a[2] = 12; a[3] = 12; a[4] = 12;
                           std::vector<int> atest(5);
                           tensorCub.getAccuracy(atest);
                           TEST_FOR_EXCEPTION( (tensorCub.getDimension() != 10) || (a != atest),
                                               std::logic_error,
                                               "Check constructor and members dimension and getAccuracy!" ) );
    if (TestForException_getThrowNumber() != endThrowNumber) {
      errorFlag = -1000;
    }
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  };
 
  *outStream \
  << "===============================================================================\n"\
  << "| TEST 2: volume computations                                                 |\n"\
  << "===============================================================================\n";

  double testVol = 0.0;
  double tol     = 100.0 * INTREPID_TOL;

  // list of analytic volume values, listed in the enumerated reference cell order up to CELL_HEXAPRISM
  double volumeList[] = {0.0, 2.0, 1.0/2.0, 4.0, 1.0/6.0, 8.0, 1.0, 32.0};

  *outStream << "\nReference cell volumes:\n\n";

  try {
    shards::CellTopology line(shards::getCellTopologyData< shards::Line<> >());
    for (int deg=0; deg<=INTREPID_CUBATURE_LINE_GAUSS_MAX; deg++) {
      testVol = computeRefVolume(line, deg);
      *outStream << std::setw(30) << "Line volume --> " << std::setw(10) << std::scientific << testVol <<
                    std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - volumeList[1]) << "\n";
      if (std::abs(testVol - volumeList[1]) > tol) {
        errorFlag = 1;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      }
    }
    *outStream << "\n\n";
    shards::CellTopology tri(shards::getCellTopologyData< shards::Triangle<> >());
    for (int deg=0; deg<=INTREPID_CUBATURE_TRI_DEFAULT_MAX; deg++) {
      testVol = computeRefVolume(tri, deg);
      *outStream << std::setw(30) << "Triangle volume --> " << std::setw(10) << std::scientific << testVol <<
                    std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - volumeList[2]) << "\n";
      if (std::abs(testVol - volumeList[2]) > tol) {
        errorFlag = 1;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      }
    }
    *outStream << "\n\n";
    shards::CellTopology quad(shards::getCellTopologyData< shards::Quadrilateral<> >());
    for (int deg=0; deg<=INTREPID_CUBATURE_LINE_GAUSS_MAX; deg++) {
      testVol = computeRefVolume(quad, deg);
      *outStream << std::setw(30) << "Quadrilateral volume --> " << std::setw(10) << std::scientific << testVol <<
                    std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - volumeList[3]) << "\n";
      if (std::abs(testVol - volumeList[3]) > tol) {
        errorFlag = 1;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      }
    }
    *outStream << "\n\n";
    shards::CellTopology tet(shards::getCellTopologyData< shards::Tetrahedron<> >());
    for (int deg=0; deg<=INTREPID_CUBATURE_TET_DEFAULT_MAX; deg++) {
      testVol = computeRefVolume(tet, deg);
      *outStream << std::setw(30) << "Tetrahedron volume --> " << std::setw(10) << std::scientific << testVol <<
                    std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - volumeList[4]) << "\n";
      if (std::abs(testVol - volumeList[4]) > tol) {
        errorFlag = 1;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      }
    }
    *outStream << "\n\n";
    shards::CellTopology hex(shards::getCellTopologyData< shards::Hexahedron<> >());
    for (int deg=0; deg<=INTREPID_CUBATURE_LINE_GAUSS_MAX; deg++) {
      testVol = computeRefVolume(hex, deg);
      *outStream << std::setw(30) << "Hexahedron volume --> " << std::setw(10) << std::scientific << testVol <<
                    std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - volumeList[5]) << "\n";
      if (std::abs(testVol - volumeList[5]) > tol) {
        errorFlag = 1;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      }
    }
    *outStream << "\n\n";
    shards::CellTopology wedge(shards::getCellTopologyData< shards::Wedge<> >());
    for (int deg=0; deg<=std::min(INTREPID_CUBATURE_LINE_GAUSS_MAX,INTREPID_CUBATURE_TRI_DEFAULT_MAX); deg++) {
      testVol = computeRefVolume(wedge, deg);
      *outStream << std::setw(30) << "Wedge volume --> " << std::setw(10) << std::scientific << testVol <<
                    std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - volumeList[6]) << "\n";
      if (std::abs(testVol - volumeList[6]) > tol) {
        errorFlag = 1;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
      }
    }
    *outStream << "\n\n";
    for (int deg=0; deg<=20; deg++) {
      Teuchos::RCP<CubatureDirectLineGauss<double> > lineCub = Teuchos::rcp(new CubatureDirectLineGauss<double>(deg));
      CubatureTensor<double> hypercubeCub(lineCub, 5);
      int numCubPoints = hypercubeCub.getNumPoints();
      FieldContainer<double> cubPoints( numCubPoints, hypercubeCub.getDimension() );
      FieldContainer<double> cubWeights( numCubPoints );
      hypercubeCub.getCubature(cubPoints, cubWeights);
      testVol = 0;
      for (int i=0; i<numCubPoints; i++)
        testVol += cubWeights[i];
      *outStream << std::setw(30) << "5-D Hypercube volume --> " << std::setw(10) << std::scientific << testVol <<
                    std::setw(10) << "diff = " << std::setw(10) << std::scientific << std::abs(testVol - volumeList[7]) << "\n";
      if (std::abs(testVol - volumeList[7])/std::abs(testVol) > tol) {
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
