// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file
\brief  Unit test (CubatureControlVolume, CubatureControlVolume Side): 
        correctness of values.
\author Created by K. Peterson, P. Bochev and D. Ridzal.
*/

#include "Intrepid2_CubatureControlVolume.hpp"
#include "Intrepid2_CubatureControlVolumeSide.hpp"
#include "Intrepid2_Utils.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"


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
  << "|                         Unit Test (CubatureControlVolume)                   |\n" \
  << "|                                   (CubatureControlVolumeSide)               |\n" \
  << "|                                                                             |\n" \
  << "|     1) Correctness of cubature points and weights                           |\n"\
  << "|     2) Comparison of sub-control volume weights and primary cell volume     |\n"\
  << "|     3) Control volume integration                                           |\n"\
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" \
  << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| TEST 1: correctness of cubature points and weights                          |\n"\
  << "===============================================================================\n";

  int errorFlag = 0;

  Teuchos::RCP<shards::CellTopology> cellTopo;
  Teuchos::RCP<Intrepid2::Cubature<double,Intrepid2::FieldContainer<double>  > > controlVolCub;
  Teuchos::RCP<Intrepid2::Cubature<double,Intrepid2::FieldContainer<double>  > > controlVolSideCub;


  // quadrilateral primary cells
  try {

    // set cell topology
    cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<> >()));

    // define control volume cubature rule
    controlVolCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
    int numPoints = controlVolCub->getNumPoints();

    // define control volume side cubature rule
    controlVolSideCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolumeSide<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
    int numSidePoints = controlVolSideCub->getNumPoints();

    // define quad coordinates
    Intrepid2::FieldContainer<double> cellCoords(2,4,2);
    cellCoords(0,0,0) = 0.0; cellCoords(0,0,1) = 0.0;
    cellCoords(0,1,0) = 0.5; cellCoords(0,1,1) = 0.0;
    cellCoords(0,2,0) = 0.5; cellCoords(0,2,1) = 1.0;
    cellCoords(0,3,0) = 0.0; cellCoords(0,3,1) = 1.0;
    cellCoords(1,0,0) = 0.5; cellCoords(1,0,1) = 0.0;
    cellCoords(1,1,0) = 1.0; cellCoords(1,1,1) = 0.0;
    cellCoords(1,2,0) = 1.0; cellCoords(1,2,1) = 1.0;
    cellCoords(1,3,0) = 0.5; cellCoords(1,3,1) = 1.0;

     // points 
    double exactCubPoints[] = {
      0.125, 0.25, 0.375, 0.25,
      0.375, 0.75, 0.125, 0.75,
      0.625, 0.25, 0.875, 0.25,
      0.875, 0.75, 0.625, 0.75
    };

     // weights
    double exactCubWeights[] = {
      0.125, 0.125, 0.125, 0.125,
      0.125, 0.125, 0.125, 0.125
    };

     // side points 
    double exactSideCubPoints[] = {
      0.25, 0.25, 0.375, 0.5,
      0.25, 0.75, 0.125, 0.5,
      0.75, 0.25, 0.875, 0.5,
      0.75, 0.75, 0.625, 0.5
    };

     // side weights (these are weighted normals!)
    double exactSideCubWeights[] = {
      0.5, 0.0, 0.0, 0.25,
     -0.5, 0.0, 0.0,-0.25,
      0.5, 0.0, 0.0, 0.25,
     -0.5, 0.0, 0.0,-0.25
    };

    int exactPoints = 4;
    if (numPoints != exactPoints) {
       errorFlag++;
       *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
       *outStream << "Number of cubature points: " << numPoints;
       *outStream << "   does not equal correct number: " << exactPoints << "\n";
    }

    // get cubature points and weights for volume integration over control volume
    int numCells = 2; int spaceDim = 2;
    Intrepid2::FieldContainer<double> cubPoints(numCells,numPoints,spaceDim);
    Intrepid2::FieldContainer<double> cubWeights(numCells,numPoints);
    controlVolCub->getCubature(cubPoints,cubWeights,cellCoords);

    // get cubature points and weights for surface integration over control volume
    // (Note: the weights here are really weighted normals)
    Intrepid2::FieldContainer<double> sideCubPoints(numCells,numSidePoints,spaceDim);
    Intrepid2::FieldContainer<double> sideCubWeights(numCells,numSidePoints,spaceDim);
    controlVolSideCub->getCubature(sideCubPoints,sideCubWeights,cellCoords);

    int countp = 0;
    int countw = 0;
    for (int i=0; i<numCells; i++) {
       for (int j=0; j<numPoints; j++) {
          for (int k=0; k<spaceDim; k++) {

             // check values of cubature points
             if (std::abs(cubPoints(i,j,k) - exactCubPoints[countp]) > Intrepid2::INTREPID_TOL) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed point: " << cubPoints(i,j,k)
                 << " but reference value: " << exactCubPoints[countp] << "\n";
             }

             // check values of side cubature points
             if (std::abs(sideCubPoints(i,j,k) - exactSideCubPoints[countp]) > Intrepid2::INTREPID_TOL) {
                errorFlag++;
                 *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed point: " << sideCubPoints(i,j,k)
                 << " but reference value: " << exactSideCubPoints[countp] << "\n";
             }

             // check values of side cubature weights
             if (std::abs(sideCubWeights(i,j,k) - exactSideCubWeights[countp]) > Intrepid2::INTREPID_TOL) {
                errorFlag++;
                 *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed weight: " << sideCubWeights(i,j,k)
                 << " but reference value: " << exactSideCubWeights[countp] << "\n";
             }

            countp++;

          }

          // check values of cubature weights
          if (std::abs(cubWeights(i,j) - exactCubWeights[countw]) > Intrepid2::INTREPID_TOL) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";
             *outStream << "}  computed weight: " << cubWeights(i,j)
               << " but reference value: " << exactCubWeights[countw] << "\n";
          }
          countw++;
       }
    }

  }
  catch (std::exception err) {
    *outStream << err.what() << "\n";
    errorFlag = -1;
  }  

  // triangle primary cells
  try {

    // set cell topology
    cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Triangle<> >()));

    // define control volume cubature rule
    controlVolCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
    int numPoints = controlVolCub->getNumPoints();

    // define control volume side cubature rule
    controlVolSideCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolumeSide<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
    int numSidePoints = controlVolSideCub->getNumPoints();

    // define triangle coordinates
    Intrepid2::FieldContainer<double> cellCoords(2,3,2);
    cellCoords(0,0,0) = 0.0; cellCoords(0,0,1) = 0.0;
    cellCoords(0,1,0) = 0.5; cellCoords(0,1,1) = 0.0;
    cellCoords(0,2,0) = 0.0; cellCoords(0,2,1) = 0.5;
    cellCoords(1,0,0) = 0.5; cellCoords(1,0,1) = 0.0;
    cellCoords(1,1,0) = 0.5; cellCoords(1,1,1) = 0.5;
    cellCoords(1,2,0) = 0.0; cellCoords(1,2,1) = 0.5;

     // points 
    double exactCubPoints[] = {
      0.1041666666666667, 0.1041666666666667, 0.2916666666666667, 0.1041666666666667,
      0.1041666666666667, 0.2916666666666667, 0.3958333333333333, 0.2083333333333333,
      0.3958333333333333, 0.3958333333333333, 0.2083333333333333, 0.3958333333333333
    };

     // weights
    double exactCubWeights[] = {
      0.0416666666666667, 0.0416666666666667, 0.0416666666666667, 
      0.0416666666666667, 0.0416666666666667, 0.0416666666666667 
    };

     // side points 
    double exactSideCubPoints[] = {
      0.2083333333333333, 0.0833333333333333, 0.2083333333333333, 0.2083333333333333,
      0.0833333333333333, 0.2083333333333333, 0.4166666666666667, 0.2916666666666667,
      0.2916666666666667, 0.4166666666666667, 0.2916666666666667, 0.2916666666666667
    };

     // side weights (these are weighted normals!)
    double exactSideCubWeights[] = {
      0.1666666666666667, 0.0833333333333333,-0.0833333333333333, 0.0833333333333333,
     -0.0833333333333333,-0.1666666666666667, 0.0833333333333333, 0.1666666666666667,
     -0.1666666666666667,-0.0833333333333333, 0.0833333333333333,-0.0833333333333333
    };

    // same number of points for volume and side in 2-D
    int exactPoints = 3;
    if (numPoints != exactPoints) {
       errorFlag++;
       *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
       *outStream << "Number of cubature points: " << numPoints;
       *outStream << "   does not equal correct number: " << exactPoints << "\n";
    }
    if (numSidePoints != exactPoints) {
       errorFlag++;
       *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
       *outStream << "Number of side cubature points: " << numSidePoints;
       *outStream << "   does not equal correct number: " << exactPoints << "\n";
    }

    // get cubature points and weights for volume integration over control volume
    int numCells = 2; int spaceDim = 2;
    Intrepid2::FieldContainer<double> cubPoints(numCells,numPoints,spaceDim);
    Intrepid2::FieldContainer<double> cubWeights(numCells,numPoints);
    controlVolCub->getCubature(cubPoints,cubWeights,cellCoords);

    // get cubature points and weights for surface integration over control volume
    // (Note: the weights here are really weighted normals)
    Intrepid2::FieldContainer<double> sideCubPoints(numCells,numSidePoints,spaceDim);
    Intrepid2::FieldContainer<double> sideCubWeights(numCells,numSidePoints,spaceDim);
    controlVolSideCub->getCubature(sideCubPoints,sideCubWeights,cellCoords);

    int countp = 0;
    int countw = 0;
    for (int i=0; i<numCells; i++) {
       for (int j=0; j<numPoints; j++) {
          for (int k=0; k<spaceDim; k++) {

             // check values of cubature points
             if (std::abs(cubPoints(i,j,k) - exactCubPoints[countp]) > Intrepid2::INTREPID_TOL) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed point: " << cubPoints(i,j,k)
                 << " but reference value: " << exactCubPoints[countp] << "\n";
             }

             // check values of side cubature points
             if (std::abs(sideCubPoints(i,j,k) - exactSideCubPoints[countp]) > Intrepid2::INTREPID_TOL) {
                errorFlag++;
                 *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed point: " << sideCubPoints(i,j,k)
                 << " but reference value: " << exactSideCubPoints[countp] << "\n";
             }

             // check values of side cubature weights
             if (std::abs(sideCubWeights(i,j,k) - exactSideCubWeights[countp]) > Intrepid2::INTREPID_TOL) {
                errorFlag++;
                 *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed weight: " << sideCubWeights(i,j,k)
                 << " but reference value: " << exactSideCubWeights[countp] << "\n";
             }

            countp++;

          }

          // check values of cubature weights
          if (std::abs(cubWeights(i,j) - exactCubWeights[countw]) > Intrepid2::INTREPID_TOL) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";
             *outStream << "}  computed weight: " << cubWeights(i,j)
               << " but reference value: " << exactCubWeights[countw] << "\n";
          }
          countw++;
       }
    }

  }
  catch (std::exception err) {
    *outStream << err.what() << "\n";
    errorFlag = -1;
  }  

  // tetrahedral primary cells
  try {

    // set cell topology
    cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<> >()));

    // define control volume cubature rule
    controlVolCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
    int numPoints = controlVolCub->getNumPoints();

    // define control volume side cubature rule
    controlVolSideCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolumeSide<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
    int numSidePoints = controlVolSideCub->getNumPoints();

    // define tetrahedron coordinates
    Intrepid2::FieldContainer<double> cellCoords(1,4,3);
    cellCoords(0,0,0) = 0.0; cellCoords(0,0,1) = 0.0; cellCoords(0,0,2) = 0.0;
    cellCoords(0,1,0) = 1.0; cellCoords(0,1,1) = 0.0; cellCoords(0,1,2) = 0.0;
    cellCoords(0,2,0) = 0.0; cellCoords(0,2,1) = 1.0; cellCoords(0,2,2) = 0.0;
    cellCoords(0,3,0) = 0.0; cellCoords(0,3,1) = 0.0; cellCoords(0,3,2) = 1.0;

     // points 
    double exactCubPoints[] = {
      0.17708333333333333,  0.17708333333333333,  0.17708333333333333, 
      0.46875000000000000,  0.17708333333333333,  0.17708333333333333, 
      0.17708333333333333,  0.46875000000000000,  0.17708333333333333, 
      0.17708333333333333,  0.17708333333333333,  0.46875000000000000
    };

     // weights
    double exactCubWeights[] = {
      0.0416666666666667, 0.0416666666666667, 
      0.0416666666666667, 0.0416666666666667 
    };

     // side points 
    double exactSideCubPoints[] = {
      0.3541666666666667, 0.1458333333333333, 0.1458333333333333,
      0.3541666666666667, 0.3541666666666667, 0.1458333333333333,
      0.1458333333333333, 0.3541666666666667, 0.1458333333333333,
      0.1458333333333333, 0.1458333333333333, 0.3541666666666667,
      0.3541666666666667, 0.1458333333333333, 0.3541666666666667,
      0.1458333333333333, 0.3541666666666667, 0.3541666666666667
    };

     // side weights (these are weighted normals!)
    double exactSideCubWeights[] = {
      0.0833333333333333, 0.0416666666666667, 0.041666666666667,
     -0.0416666666666667, 0.0416666666666667, 0.000000000000000,
     -0.0416666666666667,-0.0833333333333333,-0.041666666666667,
      0.0416666666666667, 0.0416666666666667, 0.083333333333333,
     -0.0416666666666667, 0.0000000000000000, 0.041666666666667,
      0.0000000000000000,-0.0416666666666667, 0.041666666666667
    };

    // volume points
    int exactPoints = 4;
    if (numPoints != exactPoints) {
       errorFlag++;
       *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
       *outStream << "Number of cubature points: " << numPoints;
       *outStream << "   does not equal correct number: " << exactPoints << "\n";
    }
    // side points
    exactPoints = 6;
    if (numSidePoints != exactPoints) {
       errorFlag++;
       *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
       *outStream << "Number of side cubature points: " << numSidePoints;
       *outStream << "   does not equal correct number: " << exactPoints << "\n";
    }

    // get cubature points and weights for volume integration over control volume
    int numCells = 1; int spaceDim = 3;
    Intrepid2::FieldContainer<double> cubPoints(numCells,numPoints,spaceDim);
    Intrepid2::FieldContainer<double> cubWeights(numCells,numPoints);
    controlVolCub->getCubature(cubPoints,cubWeights,cellCoords);

    // get cubature points and weights for surface integration over control volume
    // (Note: the weights here are really weighted normals)
    Intrepid2::FieldContainer<double> sideCubPoints(numCells,numSidePoints,spaceDim);
    Intrepid2::FieldContainer<double> sideCubWeights(numCells,numSidePoints,spaceDim);
    controlVolSideCub->getCubature(sideCubPoints,sideCubWeights,cellCoords);


    int countp = 0;
    int countw = 0;
    for (int i=0; i<numCells; i++) {
       for (int j=0; j<numPoints; j++) {
          for (int k=0; k<spaceDim; k++) {

             // check values of cubature points
             if (std::abs(cubPoints(i,j,k) - exactCubPoints[countp]) > Intrepid2::INTREPID_TOL) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed point: " << cubPoints(i,j,k)
                 << " but reference value: " << exactCubPoints[countp] << "\n";
             }

             // check values of side cubature points
             if (std::abs(sideCubPoints(i,j,k) - exactSideCubPoints[countp]) > Intrepid2::INTREPID_TOL) {
                errorFlag++;
                 *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed point: " << sideCubPoints(i,j,k)
                 << " but reference value: " << exactSideCubPoints[countp] << "\n";
             }

             // check values of side cubature weights
             if (std::abs(sideCubWeights(i,j,k) - exactSideCubWeights[countp]) > Intrepid2::INTREPID_TOL) {
                errorFlag++;
                 *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed weight: " << sideCubWeights(i,j,k)
                 << " but reference value: " << exactSideCubWeights[countp] << "\n";
             }

            countp++;

          }

          // check values of cubature weights
          if (std::abs(cubWeights(i,j) - exactCubWeights[countw]) > Intrepid2::INTREPID_TOL) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";
             *outStream << "}  computed weight: " << cubWeights(i,j)
               << " but reference value: " << exactCubWeights[countw] << "\n";
                 *outStream << cubWeights(i,j)-exactCubWeights[countp] << "\n";
          }
          countw++;
       }
    }


  }
  catch (std::exception err) {
    *outStream << err.what() << "\n";
    errorFlag = -1;
  }  

  // hexahedral primary cells
  try {

    // set cell topology
    cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<> >()));

    // define control volume cubature rule
    controlVolCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
    int numPoints = controlVolCub->getNumPoints();

    // define control volume side cubature rule
    controlVolSideCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolumeSide<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
    int numSidePoints = controlVolSideCub->getNumPoints();

    // define hexahedron coordinates
    Intrepid2::FieldContainer<double> cellCoords(1,8,3);
    cellCoords(0,0,0) = 0.0; cellCoords(0,0,1) = 0.0; cellCoords(0,0,2) = 0.0;
    cellCoords(0,1,0) = 2.0; cellCoords(0,1,1) = 0.0; cellCoords(0,1,2) = 0.0;
    cellCoords(0,2,0) = 2.0; cellCoords(0,2,1) = 1.5; cellCoords(0,2,2) = 0.0;
    cellCoords(0,3,0) = 0.0; cellCoords(0,3,1) = 1.5; cellCoords(0,3,2) = 0.0;
    cellCoords(0,4,0) = 0.0; cellCoords(0,4,1) = 0.0; cellCoords(0,4,2) = 1.0;
    cellCoords(0,5,0) = 2.0; cellCoords(0,5,1) = 0.0; cellCoords(0,5,2) = 1.0;
    cellCoords(0,6,0) = 2.0; cellCoords(0,6,1) = 1.5; cellCoords(0,6,2) = 1.0;
    cellCoords(0,7,0) = 0.0; cellCoords(0,7,1) = 1.5; cellCoords(0,7,2) = 1.0;

     // points 
    double exactCubPoints[] = {
     0.5, 0.375, 0.25, 1.5, 0.375, 0.25,
     1.5, 1.125, 0.25, 0.5, 1.125, 0.25,
     0.5, 0.375, 0.75, 1.5, 0.375, 0.75,
     1.5, 1.125, 0.75, 0.5, 1.125, 0.75
    };

     // weights
    double exactCubWeights[] = {
     0.375, 0.375, 0.375, 0.375, 
     0.375, 0.375, 0.375, 0.375
    };

     // side points 
    double exactSideCubPoints[] = {
     1.0, 0.375, 0.25, 1.5, 0.750, 0.25,
     1.0, 1.125, 0.25, 0.5, 0.750, 0.25, 
     1.0, 0.375, 0.75, 1.5, 0.750, 0.75, 
     1.0, 1.125, 0.75, 0.5, 0.750, 0.75,
     0.5, 0.375, 0.50, 1.5, 0.375, 0.50,
     1.5, 1.125, 0.50, 0.5, 1.125, 0.50
    };

     // side weights (these are weighted normals!)
    double exactSideCubWeights[] = {
     0.375, 0.00, 0.00, 0.00, 0.50, 0.00,
    -0.375, 0.00, 0.00, 0.00,-0.50, 0.00,
     0.375, 0.00, 0.00, 0.00, 0.50, 0.00,
    -0.375, 0.00, 0.00, 0.00,-0.50, 0.00,
     0.000, 0.00, 0.75, 0.00, 0.00, 0.75,
     0.000, 0.00, 0.75, 0.00, 0.00, 0.75,
    };

    // volume points
    int exactPoints = 8;
    if (numPoints != exactPoints) {
       errorFlag++;
       *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
       *outStream << "Number of cubature points: " << numPoints;
       *outStream << "   does not equal correct number: " << exactPoints << "\n";
    }
    // side points
    exactPoints = 12;
    if (numSidePoints != exactPoints) {
       errorFlag++;
       *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
       *outStream << "Number of side cubature points: " << numSidePoints;
       *outStream << "   does not equal correct number: " << exactPoints << "\n";
    }

    // get cubature points and weights for volume integration over control volume
    int numCells = 1; int spaceDim = 3;
    Intrepid2::FieldContainer<double> cubPoints(numCells,numPoints,spaceDim);
    Intrepid2::FieldContainer<double> cubWeights(numCells,numPoints);
    controlVolCub->getCubature(cubPoints,cubWeights,cellCoords);

    // get cubature points and weights for surface integration over control volume
    // (Note: the weights here are really weighted normals)
    Intrepid2::FieldContainer<double> sideCubPoints(numCells,numSidePoints,spaceDim);
    Intrepid2::FieldContainer<double> sideCubWeights(numCells,numSidePoints,spaceDim);
    controlVolSideCub->getCubature(sideCubPoints,sideCubWeights,cellCoords);


    int countp = 0;
    int countw = 0;
    for (int i=0; i<numCells; i++) {
       for (int j=0; j<numPoints; j++) {
          for (int k=0; k<spaceDim; k++) {

             // check values of cubature points
             if (std::abs(cubPoints(i,j,k) - exactCubPoints[countp]) > Intrepid2::INTREPID_TOL) {
                errorFlag++;
                *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed point: " << cubPoints(i,j,k)
                 << " but reference value: " << exactCubPoints[countp] << "\n";
             }

             // check values of side cubature points
             if (std::abs(sideCubPoints(i,j,k) - exactSideCubPoints[countp]) > Intrepid2::INTREPID_TOL) {
                errorFlag++;
                 *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed point: " << sideCubPoints(i,j,k)
                 << " but reference value: " << exactSideCubPoints[countp] << "\n";
             }

             // check values of side cubature weights
             if (std::abs(sideCubWeights(i,j,k) - exactSideCubWeights[countp]) > Intrepid2::INTREPID_TOL) {
                errorFlag++;
                 *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

                // Output the multi-index of the value where the error is:
                *outStream << " At multi-index { ";
                *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
                *outStream << "}  computed weight: " << sideCubWeights(i,j,k)
                 << " but reference value: " << exactSideCubWeights[countp] << "\n";
             }

            countp++;

          }

          // check values of cubature weights
          if (std::abs(cubWeights(i,j) - exactCubWeights[countw]) > Intrepid2::INTREPID_TOL) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";
             *outStream << "}  computed weight: " << cubWeights(i,j)
               << " but reference value: " << exactCubWeights[countw] << "\n";
                 *outStream << cubWeights(i,j)-exactCubWeights[countp] << "\n";
          }
          countw++;
       }
    }


  }
  catch (std::exception err) {
    *outStream << err.what() << "\n";
    errorFlag = -1;
  }  

  *outStream \
  << "===============================================================================\n"\
  << "| TEST 2: comparison of sub-control volume weights and primary cell volume    |\n"\
  << "===============================================================================\n";

  // quadrilateral primary cells
  try {

    // set cell topology
    cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<> >()));

    // define control volume cubature rule
    controlVolCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
    int numPoints = controlVolCub->getNumPoints();

    // define quad coordinates
    Intrepid2::FieldContainer<double> cellCoords(1,4,2);
    cellCoords(0,0,0) = 0.0; cellCoords(0,0,1) = 0.0;
    cellCoords(0,1,0) = 2.4; cellCoords(0,1,1) = 0.0;
    cellCoords(0,2,0) = 2.4; cellCoords(0,2,1) = 3.1;
    cellCoords(0,3,0) = 0.0; cellCoords(0,3,1) = 3.1;

    double exact_area = 2.4*3.1;

    // get cubature points and weights 
    int numCells = 1; int spaceDim = 2;
    Intrepid2::FieldContainer<double> cubPoints(numCells,numPoints,spaceDim);
    Intrepid2::FieldContainer<double> cubWeights(numCells,numPoints);
    controlVolCub->getCubature(cubPoints,cubWeights,cellCoords);

    // loop over number of points (equals number of sub-control volumes) and check total volume
    double total_area = 0.0;
    for (int i=0; i<numPoints; i++) {
       total_area += cubWeights(0,i);
    }
    if (std::abs(total_area - exact_area) > Intrepid2::INTREPID_TOL) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Sum of sub-control volume areas: ";
        *outStream << total_area;
        *outStream << " does not equal quad primary cell area: " << exact_area << "\n";
    }

  }
  catch (std::exception err) {
    *outStream << err.what() << "\n";
    errorFlag = -1;
  }  

  // triangle primary cells
  try {

    // set cell topology
    cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Triangle<> >()));

    // define control volume cubature rule
    controlVolCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
    int numPoints = controlVolCub->getNumPoints();

    // define quad coordinates
    Intrepid2::FieldContainer<double> cellCoords(1,3,2);
    cellCoords(0,0,0) = 0.0; cellCoords(0,0,1) = 0.0;
    cellCoords(0,1,0) = 3.6; cellCoords(0,1,1) = 0.0;
    cellCoords(0,2,0) = 0.0; cellCoords(0,2,1) = 2.8;

    double exact_area = 0.5*3.6*2.8;

    // get cubature points and weights 
    int numCells = 1; int spaceDim = 2;
    Intrepid2::FieldContainer<double> cubPoints(numCells,numPoints,spaceDim);
    Intrepid2::FieldContainer<double> cubWeights(numCells,numPoints);
    controlVolCub->getCubature(cubPoints,cubWeights,cellCoords);

    // loop over number of points (equals number of sub-control volumes) and check total volume
    double total_area = 0.0;
    for (int i=0; i<numPoints; i++) {
       total_area += cubWeights(0,i);
    }
    if (std::abs(total_area - exact_area) > Intrepid2::INTREPID_TOL) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Sum of sub-control volume areas: ";
        *outStream << total_area;
        *outStream << " does not equal triangle primary cell area: " << exact_area << "\n";
    }

  }
  catch (std::exception err) {
    *outStream << err.what() << "\n";
    errorFlag = -1;
  }  

  // hexahedral primary cells
  try {

    // set cell topology
    cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Hexahedron<> >()));

    // define control volume cubature rule
    controlVolCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
    int numPoints = controlVolCub->getNumPoints();

    // define hexahedron coordinates
    Intrepid2::FieldContainer<double> cellCoords(1,8,3);
    cellCoords(0,0,0) = 0.0; cellCoords(0,0,1) = 0.0; cellCoords(0,0,2) = 0.0;
    cellCoords(0,1,0) = 2.4; cellCoords(0,1,1) = 0.0; cellCoords(0,1,2) = 0.0;
    cellCoords(0,2,0) = 2.4; cellCoords(0,2,1) = 3.1; cellCoords(0,2,2) = 0.0;
    cellCoords(0,3,0) = 0.0; cellCoords(0,3,1) = 3.1; cellCoords(0,3,2) = 0.0;
    cellCoords(0,4,0) = 0.0; cellCoords(0,4,1) = 0.0; cellCoords(0,4,2) = 1.7;
    cellCoords(0,5,0) = 2.4; cellCoords(0,5,1) = 0.0; cellCoords(0,5,2) = 1.7;
    cellCoords(0,6,0) = 2.4; cellCoords(0,6,1) = 3.1; cellCoords(0,6,2) = 1.7;
    cellCoords(0,7,0) = 0.0; cellCoords(0,7,1) = 3.1; cellCoords(0,7,2) = 1.7;

    double exact_vol = 2.4*3.1*1.7;

    // get cubature points and weights 
    int numCells = 1; int spaceDim = 3;
    Intrepid2::FieldContainer<double> cubPoints(numCells,numPoints,spaceDim);
    Intrepid2::FieldContainer<double> cubWeights(numCells,numPoints);
    controlVolCub->getCubature(cubPoints,cubWeights,cellCoords);

    // loop over number of points (equals number of sub-control volumes) and check total volume
    double total_vol = 0.0;
    for (int i=0; i<numPoints; i++) {
       total_vol += cubWeights(0,i);
    }
    if (std::abs(total_vol - exact_vol) > Intrepid2::INTREPID_TOL) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Sum of sub-control volumes: ";
        *outStream << total_vol;
        *outStream << " does not equal hexahedron primary cell volume: " << exact_vol << "\n";
    }

  }
  catch (std::exception err) {
    *outStream << err.what() << "\n";
    errorFlag = -1;
  }  

  // tetrahedral primary cells
  try {

    // set cell topology
    cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Tetrahedron<> >()));

    // define control volume cubature rule
    controlVolCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
    int numPoints = controlVolCub->getNumPoints();

    // define tetrahedron coordinates
    Intrepid2::FieldContainer<double> cellCoords(1,4,3);
    cellCoords(0,0,0) = 0.0; cellCoords(0,0,1) = 0.0; cellCoords(0,0,2) = 0.0;
    cellCoords(0,1,0) = 3.6; cellCoords(0,1,1) = 0.0; cellCoords(0,1,2) = 0.0;
    cellCoords(0,2,0) = 0.0; cellCoords(0,2,1) = 2.8; cellCoords(0,2,2) = 0.0;
    cellCoords(0,3,0) = 0.0; cellCoords(0,3,1) = 2.8; cellCoords(0,3,2) = 1.7;

    double exact_vol = (0.5*3.6*2.8)*1.7/3.0;

    // get cubature points and weights 
    int numCells = 1; int spaceDim = 3;
    Intrepid2::FieldContainer<double> cubPoints(numCells,numPoints,spaceDim);
    Intrepid2::FieldContainer<double> cubWeights(numCells,numPoints);
    controlVolCub->getCubature(cubPoints,cubWeights,cellCoords);

    // loop over number of points (equals number of sub-control volumes) and check total volume
    double total_vol = 0.0;
    for (int i=0; i<numPoints; i++) {
       total_vol += cubWeights(0,i);
    }
    if (std::abs(total_vol - exact_vol) > Intrepid2::INTREPID_TOL) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Sum of sub-control volumes: ";
        *outStream << total_vol;
        *outStream << " does not equal tetrahedron primary cell volume: " << exact_vol << "\n";
    }

  }
  catch (std::exception err) {
    *outStream << err.what() << "\n";
    errorFlag = -1;
  }  

  *outStream \
  << "===============================================================================\n"\
  << "| TEST 3: control volume integration                                          |\n"\
  << "===============================================================================\n";

  // quadrilateral primary cells
  try {

    // set cell topology
    cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Quadrilateral<> >()));

    // define control volume cubature rule
    controlVolCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
    int numPoints = controlVolCub->getNumPoints();

    // define control volume cubature rule
    controlVolSideCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolumeSide<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
    int numSidePoints = controlVolSideCub->getNumPoints();

    // define quad coordinates - four cells that define a complete control volume around center node
    Intrepid2::FieldContainer<double> cellCoords(4,4,2);
    cellCoords(0,0,0) = 0.0; cellCoords(0,0,1) = 0.0;
    cellCoords(0,1,0) = 0.5; cellCoords(0,1,1) = 0.0;
    cellCoords(0,2,0) = 0.41; cellCoords(0,2,1) = 0.58;
    cellCoords(0,3,0) = 0.0; cellCoords(0,3,1) = 0.5;
    cellCoords(1,0,0) = 0.5; cellCoords(1,0,1) = 0.0;
    cellCoords(1,1,0) = 1.0; cellCoords(1,1,1) = 0.0;
    cellCoords(1,2,0) = 1.0; cellCoords(1,2,1) = 0.5;
    cellCoords(1,3,0) = 0.41; cellCoords(1,3,1) = 0.58;
    cellCoords(2,0,0) = 0.0; cellCoords(2,0,1) = 0.5;
    cellCoords(2,1,0) = 0.41; cellCoords(2,1,1) = 0.58;
    cellCoords(2,2,0) = 0.5; cellCoords(2,2,1) = 1.0;
    cellCoords(2,3,0) = 0.0; cellCoords(2,3,1) = 1.0;
    cellCoords(3,0,0) = 0.41; cellCoords(3,0,1) = 0.58;
    cellCoords(3,1,0) = 1.0; cellCoords(3,1,1) = 0.5;
    cellCoords(3,2,0) = 1.0; cellCoords(3,2,1) = 1.0;
    cellCoords(3,3,0) = 0.5; cellCoords(3,3,1) = 1.0;

    // get cubature points and weights 
    int numCells = 4; int spaceDim = 2;
    Intrepid2::FieldContainer<double> cubPoints(numCells,numPoints,spaceDim);
    Intrepid2::FieldContainer<double> cubWeights(numCells,numPoints);
    controlVolCub->getCubature(cubPoints,cubWeights,cellCoords);

    Intrepid2::FieldContainer<double> sideCubPoints(numCells,numSidePoints,spaceDim);
    Intrepid2::FieldContainer<double> sideCubWeights(numCells,numSidePoints,spaceDim);
    controlVolSideCub->getCubature(sideCubPoints,sideCubWeights,cellCoords);

    // test cubature rule by checking the equality of  \int_C \nabla \cdot F dV = \int_dC F \cdot n dS
    // using F = (a x,b y), (\nabla \cdot F) = a + b.

    // first evaluate F at all control volume side points
    double a = 2.1; double b = 1.4;
    Intrepid2::FieldContainer<double> F(numCells,numSidePoints,spaceDim);
    for (int i = 0; i < numCells; i++) {
       for (int j = 0; j < numSidePoints; j++) {
          F(i,j,0) = a*sideCubPoints(i,j,0);
          F(i,j,1) = b*sideCubPoints(i,j,1);
       }
    }

    // gather the correct contributions to the surface integral
    double surfaceInt = 0.0;
    
    // contributions from first cell 
     surfaceInt += - F(0,1,0)*sideCubWeights(0,1,0) - F(0,1,1)*sideCubWeights(0,1,1)
                   + F(0,2,0)*sideCubWeights(0,2,0) + F(0,2,1)*sideCubWeights(0,2,1);
    
    // contributions from second cell 
     surfaceInt += - F(1,2,0)*sideCubWeights(1,2,0) - F(1,2,1)*sideCubWeights(1,2,1)
                   + F(1,3,0)*sideCubWeights(1,3,0) + F(1,3,1)*sideCubWeights(1,3,1);

    // contributions from third cell 
     surfaceInt += - F(2,0,0)*sideCubWeights(2,0,0) - F(2,0,1)*sideCubWeights(2,0,1)
                   + F(2,1,0)*sideCubWeights(2,1,0) + F(2,1,1)*sideCubWeights(2,1,1);

    // contributions from fourth cell 
     surfaceInt += - F(3,3,0)*sideCubWeights(3,3,0) - F(3,3,1)*sideCubWeights(3,3,1)
                   + F(3,0,0)*sideCubWeights(3,0,0) + F(3,0,1)*sideCubWeights(3,0,1);

    // gather the correct contributions to the volume integral
    double volumeInt = 0.0;

    // contributions from first cell 
      volumeInt += (a+b)*cubWeights(0,2);

    // contributions from second cell 
      volumeInt += (a+b)*cubWeights(1,3);

    // contributions from third cell 
      volumeInt += (a+b)*cubWeights(2,1);

    // contributions from fourth cell 
      volumeInt += (a+b)*cubWeights(3,0);

    if (std::abs(surfaceInt - volumeInt) > Intrepid2::INTREPID_TOL) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Integral of (F cdot n) over surface : ";
        *outStream << surfaceInt;
        *outStream << " does not equal integral of div F over volume: " << volumeInt << "\n";
    }

  }
  catch (std::exception err) {
    *outStream << err.what() << "\n";
    errorFlag = -1;
  }  

  // triangle primary cells
  try {

    // set cell topology
    cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData< shards::Triangle<> >()));

    // define control volume cubature rule
    controlVolCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
    int numPoints = controlVolCub->getNumPoints();

    // define control volume cubature rule
    controlVolSideCub  = Teuchos::rcp(new Intrepid2::CubatureControlVolumeSide<double,Intrepid2::FieldContainer<double>,Intrepid2::FieldContainer<double> >(cellTopo));
    int numSidePoints = controlVolSideCub->getNumPoints();

    // define triangle coordinates - four cells that define a complete control volume around center node
    Intrepid2::FieldContainer<double> cellCoords(4,3,2);
    cellCoords(0,0,0) = 0.0; cellCoords(0,0,1) = 0.0;
    cellCoords(0,1,0) = 1.0; cellCoords(0,1,1) = 0.0;
    cellCoords(0,2,0) = 0.41; cellCoords(0,2,1) = 0.58;
    cellCoords(1,0,0) = 1.0; cellCoords(1,0,1) = 0.0;
    cellCoords(1,1,0) = 1.0; cellCoords(1,1,1) = 1.0;
    cellCoords(1,2,0) = 0.41; cellCoords(1,2,1) = 0.58;
    cellCoords(2,0,0) = 1.0; cellCoords(2,0,1) = 1.0;
    cellCoords(2,1,0) = 0.0; cellCoords(2,1,1) = 1.0;
    cellCoords(2,2,0) = 0.41; cellCoords(2,2,1) = 0.58;
    cellCoords(3,0,0) = 0.0; cellCoords(3,0,1) = 1.0;
    cellCoords(3,1,0) = 0.0; cellCoords(3,1,1) = 0.0;
    cellCoords(3,2,0) = 0.41; cellCoords(3,2,1) = 0.58;

    // get cubature points and weights 
    int numCells = 4; int spaceDim = 2;
    Intrepid2::FieldContainer<double> cubPoints(numCells,numPoints,spaceDim);
    Intrepid2::FieldContainer<double> cubWeights(numCells,numPoints);
    controlVolCub->getCubature(cubPoints,cubWeights,cellCoords);

    Intrepid2::FieldContainer<double> sideCubPoints(numCells,numSidePoints,spaceDim);
    Intrepid2::FieldContainer<double> sideCubWeights(numCells,numSidePoints,spaceDim);
    controlVolSideCub->getCubature(sideCubPoints,sideCubWeights,cellCoords);

    // test cubature rule by checking the equality of  \int_C \nabla \cdot F dV = \int_dC F \cdot n dS
    // using F = (a x,b y), (\nabla \cdot F) = a + b.

    // first evaluate F at all control volume side points
    double a = 2.1; double b = 1.4;
    Intrepid2::FieldContainer<double> F(numCells,numSidePoints,spaceDim);
    for (int i = 0; i < numCells; i++) {
       for (int j = 0; j < numSidePoints; j++) {
          F(i,j,0) = a*sideCubPoints(i,j,0);
          F(i,j,1) = b*sideCubPoints(i,j,1);
       }
    }

    // gather the correct contributions to the surface integral
    double surfaceInt = 0.0;
    
    // contributions from first cell 
     surfaceInt += - F(0,1,0)*sideCubWeights(0,1,0) - F(0,1,1)*sideCubWeights(0,1,1)
                   + F(0,2,0)*sideCubWeights(0,2,0) + F(0,2,1)*sideCubWeights(0,2,1);
    
    // contributions from second cell 
     surfaceInt += - F(1,1,0)*sideCubWeights(1,1,0) - F(1,1,1)*sideCubWeights(1,1,1)
                   + F(1,2,0)*sideCubWeights(1,2,0) + F(1,2,1)*sideCubWeights(1,2,1);

    // contributions from third cell 
     surfaceInt += - F(2,1,0)*sideCubWeights(2,1,0) - F(2,1,1)*sideCubWeights(2,1,1)
                   + F(2,2,0)*sideCubWeights(2,2,0) + F(2,2,1)*sideCubWeights(2,2,1);

    // contributions from fourth cell 
     surfaceInt += - F(3,1,0)*sideCubWeights(3,1,0) - F(3,1,1)*sideCubWeights(3,1,1)
                   + F(3,2,0)*sideCubWeights(3,2,0) + F(3,2,1)*sideCubWeights(3,2,1);

    // gather the correct contributions to the volume integral
    double volumeInt = 0.0;

    // contributions from first cell 
      volumeInt += (a+b)*cubWeights(0,2);

    // contributions from second cell 
      volumeInt += (a+b)*cubWeights(1,2);

    // contributions from third cell 
      volumeInt += (a+b)*cubWeights(2,2);

    // contributions from fourth cell 
      volumeInt += (a+b)*cubWeights(3,2);

    if (std::abs(surfaceInt - volumeInt) > Intrepid2::INTREPID_TOL) {
        errorFlag++;
        *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
        *outStream << " Integral of (F cdot n) over surface : ";
        *outStream << surfaceInt;
        *outStream << " does not equal integral of div F over volume: " << volumeInt << "\n";
    }

  }
  catch (std::exception err) {
    *outStream << err.what() << "\n";
    errorFlag = -1;
  }  


  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
Kokkos::finalize();
  return errorFlag;
}
