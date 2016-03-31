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

/** \file test_01.cpp
\brief  Unit tests for the Intrepid2::G_TET_COMP12_FEM class.
\author Created by P. Bochev, D. Ridzal, K. Peterson, and J. Ostien.
*/
#include "Intrepid2_FieldContainer.hpp"
#include "Intrepid2_HGRAD_TET_COMP12_FEM.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace std;
using namespace Intrepid2;

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
    << "|             Unit Test (Basis_HGRAD_TET_COMP12_FEM)                          |\n" \
    << "|                                                                             |\n" \
    << "|     1) Evaluation of Basis Function Values                                  |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov),                    |\n" \
    << "|                      Jake Ostien   (jtostie@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n"\
    << "| TEST 1: Basis creation, values                                              |\n"\
    << "===============================================================================\n";

  // Define basis and error flag
  Basis_HGRAD_TET_COMP12_FEM<double, FieldContainer<double> > tetBasis;
  int errorFlag = 0;

  // Initialize throw counter for exception testing
  //int nException     = 0;
  //int throwCounter   = 0;

  // Define array containing the 10 vertices of the reference TET
  FieldContainer<double> tetNodes(10, 3);
  tetNodes(0,0) = 0.0;  tetNodes(0,1) = 0.0;  tetNodes(0,2) = 0.0;
  tetNodes(1,0) = 1.0;  tetNodes(1,1) = 0.0;  tetNodes(1,2) = 0.0;
  tetNodes(2,0) = 0.0;  tetNodes(2,1) = 1.0;  tetNodes(2,2) = 0.0;
  tetNodes(3,0) = 0.0;  tetNodes(3,1) = 0.0;  tetNodes(3,2) = 1.0;
  tetNodes(4,0) = 0.5;  tetNodes(4,1) = 0.0;  tetNodes(4,2) = 0.0;
  tetNodes(5,0) = 0.5;  tetNodes(5,1) = 0.5;  tetNodes(5,2) = 0.0;
  tetNodes(6,0) = 0.0;  tetNodes(6,1) = 0.5;  tetNodes(6,2) = 0.0;
  tetNodes(7,0) = 0.0;  tetNodes(7,1) = 0.0;  tetNodes(7,2) = 0.5;
  tetNodes(8,0) = 0.5;  tetNodes(8,1) = 0.0;  tetNodes(8,2) = 0.5;
  tetNodes(9,0) = 0.0;  tetNodes(9,1) = 0.5;  tetNodes(9,2) = 0.5;
  // Define array containing 5 integration points
  FieldContainer<double> tetPoints(9, 3);
  // from the 5 point integration
  tetPoints(0,0) = 0.25;     tetPoints(0,1) = 0.25;     tetPoints(0,2) = 0.25;
  tetPoints(1,0) = 0.5;      tetPoints(1,1) = (1./6.);  tetPoints(1,2) = (1./6.);
  tetPoints(2,0) = (1./6.);  tetPoints(2,1) = 0.5;      tetPoints(2,2) = (1./6.);
  tetPoints(3,0) = (1./6.);  tetPoints(3,1) = (1./6.);  tetPoints(3,2) = 0.5;
  tetPoints(4,0) = (1./6.);  tetPoints(4,1) = (1./6.);  tetPoints(4,2) = (1./6.);
  // from the 4 point integration
  tetPoints(5,0) = 0.1381966011250105151795413165634361882280;
  tetPoints(5,1) = 0.1381966011250105151795413165634361882280;
  tetPoints(5,2) = 0.1381966011250105151795413165634361882280;

  tetPoints(6,0) = 0.5854101966249684544613760503096914353161;
  tetPoints(6,1) = 0.1381966011250105151795413165634361882280;
  tetPoints(6,2) = 0.1381966011250105151795413165634361882280;

  tetPoints(7,0) = 0.1381966011250105151795413165634361882280;
  tetPoints(7,1) = 0.5854101966249684544613760503096914353161;
  tetPoints(7,2) = 0.1381966011250105151795413165634361882280;

  tetPoints(8,0) = 0.1381966011250105151795413165634361882280;
  tetPoints(8,1) = 0.1381966011250105151795413165634361882280;
  tetPoints(8,2) = 0.5854101966249684544613760503096914353161;

  // output precision
  outStream -> precision(20);

  // VALUE: Each row gives the 10 correct basis set values at an evaluation point
  double nodalBasisValues[] = {
    // first 4 vertices
    1.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    // second 6 vertices
    0.0, 0.0, 0.0, 0.0,  1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,  0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 1.0
  };
  double pointBasisValues[] = {
    // pt 0 {0.25, 0.25, 0.25}
    0.0, 0.0, 0.0, 0.0,  1./6., 1./6., 1./6., 1./6., 1./6., 1./6.,
    // pt 1 {0.5, 1/6, 1/6}
    0.0, 0.0, 0.0, 0.0,  1./3., 1./3., 0.0, 0.0, 1./3., 0.0,
    // pt 2 {1/6, 0.5, 0.1/6}
    0.0, 0.0, 0.0, 0.0,  0.0, 1./3., 1./3., 0.0, 0.0, 1./3.,
    // pt 3 {1/6, 1/6, 0.5}
    0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 1./3., 1./3., 1./3.,
    // pt 4 {1/6, 1/6, 1/6}
    0.0, 0.0, 0.0, 0.0,  1./3., 0.0, 1./3., 1./3., 0.0, 0.0,
    // pt 5
    0.170820393249936908922752100619382870632, 0.0, 0.0, 0.0, 0.276393202250021030359082633126872376456, 0.0, 0.276393202250021030359082633126872376456, 0.276393202250021030359082633126872376456, 0.0, 0.0,
    // pt 6
    0.0, 0.170820393249936908922752100619382870632, 0.0, 0.0, 0.276393202250021030359082633126872376456, 0.276393202250021030359082633126872376456, 0.0, 0.0, 0.276393202250021030359082633126872376456, 0.0,
    // pt 7
    0.0, 0.0, 0.170820393249936908922752100619382870632, 0.0, 0.0, 0.276393202250021030359082633126872376456, 0.276393202250021030359082633126872376456, 0.0, 0.0, 0.276393202250021030359082633126872376456,
    // pt 8
    0.0, 0.0, 0.0, 0.170820393249936908922752100619382870632, 0.0, 0.0, 0.0, 0.276393202250021030359082633126872376456, 0.276393202250021030359082633126872376456, 0.276393202250021030359082633126872376456,
  };

  // GRAD and D1: each row gives the 3x10 correct values of the gradients of the 10 basis functions
  double pointBasisGrads[] = {
    // point 0
      -1./4.,   -1./4.,   -1./4., \
       1./4.,      0.0,      0.0, \
         0.0,    1./4.,      0.0, \
         0.0,      0.0,    1./4., \
         0.0,   -3./4.,   -3./4., \
       3./4.,    3./4.,      0.0, \
      -3./4.,      0.0,   -3./4., \
      -3./4.,   -3./4.,      0.0, \
       3./4.,      0.0,    3./4., \
         0.0,    3./4.,    3./4., \

    // point 1
     -1./24.,  -1./24.,  -1./24., \
       7./8.,      0.0,      0.0, \
         0.0,   1./24.,      0.0, \
         0.0,      0.0,   1./24., \
    -35./36., -19./12., -19./12., \
     11./18.,  19./12.,      0.0, \
    -17./36.,      0.0,   -1./3., \
    -17./36.,   -1./3.,      0.0, \
     11./18.,      0.0,  19./12., \
     -5./36.,    1./3.,    1./3., \

    // point 2
     -1./24.,  -1./24.,  -1./24., \
      1./24.,      0.0,      0.0, \
         0.0,    7./8.,      0.0, \
         0.0,      0.0,   1./24., \
         0.0, -17./36.,   -1./3., \
     19./12.,  11./18.,      0.0, \
    -19./12., -35./36., -19./12., \
      -1./3., -17./36.,      0.0, \
       1./3.,  -5./36.,    1./3., \
         0.0,  11./18.,  19./12., \

    // point 3
     -1./24.,  -1./24.,  -1./24., \
      1./24.,      0.0,      0.0, \
         0.0,   1./24.,      0.0, \
         0.0,      0.0,    7./8., \
         0.0,   -1./3., -17./36., \
       1./3.,    1./3.,  -5./36., \
      -1./3.,      0.0, -17./36., \
    -19./12., -19./12., -35./36., \
     19./12.,      0.0,  11./18., \
         0.0,  19./12.,  11./18., \

    // point 4
      -7./8.,   -7./8.,   -7./8., \
      1./24.,      0.0,      0.0, \
         0.0,   1./24.,      0.0, \
         0.0,      0.0,   1./24., \
     35./36., -11./18., -11./18., \
     17./36.,  17./36.,   5./36., \
    -11./18.,  35./36., -11./18., \
    -11./18., -11./18.,  35./36., \
     17./36.,   5./36.,  17./36., \
      5./36.,  17./36.,  17./36., \

    // point 5
      -1.088525491562421136153440125774228588290, -1.088525491562421136153440125774228588290, -1.088525491562421136153440125774228588290, \
      -0.029508497187473712051146708591409529430, 0.0, 0.0, \
      0.0, -0.029508497187473712051146708591409529430, 0.0, \
      0.0, 0.0, -0.029508497187473712051146708591409529430, \
      1.30437298687487732290535130675991113734, -0.563661001875017525299235527605726980380, -0.563661001875017525299235527605726980380, \
      0.377322003750035050598471055211453960760, 0.377322003750035050598471055211453960760, 0.186338998124982474700764472394273019620, \
      -0.563661001875017525299235527605726980380, 1.30437298687487732290535130675991113734, -0.563661001875017525299235527605726980380, \
      -0.563661001875017525299235527605726980380, -0.563661001875017525299235527605726980380, 1.30437298687487732290535130675991113734, \
      0.377322003750035050598471055211453960760, 0.186338998124982474700764472394273019620, 0.377322003750035050598471055211453960760, \
      0.186338998124982474700764472394273019620, 0.377322003750035050598471055211453960760, 0.377322003750035050598471055211453960760, \

    // point 6
      0.029508497187473712051146708591409529430, 0.029508497187473712051146708591409529430, 0.029508497187473712051146708591409529430, \
      1.088525491562421136153440125774228588290, 0.0, 0.0, \
      0.0, -0.029508497187473712051146708591409529430, 0.0, \
      0.0, 0.0, -0.029508497187473712051146708591409529430, \
      -1.30437298687487732290535130675991113734, -1.868033988749894848204586834365638117720, -1.868033988749894848204586834365638117720, \
      0.563661001875017525299235527605726980380, 1.868033988749894848204586834365638117720, 0.0, \
      -0.377322003750035050598471055211453960760, 0.0, -0.190983005625052575897706582817180941140, \
      -0.377322003750035050598471055211453960760, -0.190983005625052575897706582817180941140, 0.0, \
      0.563661001875017525299235527605726980380, 0.0, 1.868033988749894848204586834365638117720, \
      -0.186338998124982474700764472394273019620, 0.190983005625052575897706582817180941140, 0.19098300562505257589770658281718094114, \

    // point 7
      0.029508497187473712051146708591409529430, 0.029508497187473712051146708591409529430, 0.029508497187473712051146708591409529430, \
      -0.029508497187473712051146708591409529430, 0.0, 0.0, \
      0.0, 1.088525491562421136153440125774228588290, 0.0, \
      0.0, 0.0, -0.029508497187473712051146708591409529430, \
      0.0, -0.377322003750035050598471055211453960760, -0.190983005625052575897706582817180941140, \
      1.868033988749894848204586834365638117720, 0.563661001875017525299235527605726980380, 0.0, \
      -1.868033988749894848204586834365638117720, -1.30437298687487732290535130675991113734, -1.868033988749894848204586834365638117720, \
      -0.190983005625052575897706582817180941140, -0.377322003750035050598471055211453960760, 0.0, \
      0.190983005625052575897706582817180941140, -0.186338998124982474700764472394273019620, 0.190983005625052575897706582817180941140, \
      0.0, 0.563661001875017525299235527605726980380, 1.868033988749894848204586834365638117720, \

    // point 8
      0.029508497187473712051146708591409529430, 0.029508497187473712051146708591409529430, 0.029508497187473712051146708591409529430, \
      -0.029508497187473712051146708591409529430, 0.0, 0.0, \
      0.0, -0.029508497187473712051146708591409529430, 0.0, \
      0.0, 0.0, 1.088525491562421136153440125774228588290, \
      0.0, -0.190983005625052575897706582817180941140, -0.377322003750035050598471055211453960760, \
      0.190983005625052575897706582817180941140, 0.190983005625052575897706582817180941140, -0.186338998124982474700764472394273019620, \
      -0.190983005625052575897706582817180941140, 0.0, -0.377322003750035050598471055211453960760, \
      -1.868033988749894848204586834365638117720, -1.868033988749894848204586834365638117720, -1.30437298687487732290535130675991113734,
      1.868033988749894848204586834365638117720, 0.0, 0.563661001875017525299235527605726980380, \
      0.0, 1.868033988749894848204586834365638117720, 0.563661001875017525299235527605726980380, \
  };

  try{

    // Dimensions for the output arrays:
    int numFields = tetBasis.getCardinality();
    int numNodes  = tetNodes.dimension(0);
    int spaceDim  = tetBasis.getBaseCellTopology().getDimension();

    // Generic array for values, grads, curls, etc. that will be properly sized before each call
    FieldContainer<double> vals;

    // Check VALUE of basis functions at nodes: resize vals to rank-2 container:\n";
    *outStream << " check VALUE of basis functions at nodes\n";
    vals.resize(numFields, numNodes);
    tetBasis.getValues(vals, tetNodes, OPERATOR_VALUE);

    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numNodes; j++) {
          int l =  i + j * numFields;
          if (std::abs(vals(i,j) - nodalBasisValues[l]) > INTREPID_TOL) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

             // Output the multi-index of the value where the error is:
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";
             *outStream << "}  computed value: " << vals(i,j)
               << " but reference value: " << nodalBasisValues[l] << "\n";
         }
      }
    }

    // Check VALUE of basis functions at points: resize vals to rank-2 container:\n";
    *outStream << " check VALUE of basis functions at points\n";
    int numPoints = tetPoints.dimension(0);
    vals.resize(numFields, numPoints);
    vals.initialize(0.0);
    tetBasis.getValues(vals, tetPoints, OPERATOR_VALUE);

    for (int i = 0; i < numFields; i++) {
      for (int j = 0; j < numPoints; j++) {
          int l =  i + j * numFields;
          if (std::abs(vals(i,j) - pointBasisValues[l]) > INTREPID_TOL) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

             // Output the multi-index of the value where the error is:
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";
             *outStream << "}  computed value: " << vals(i,j)
               << " but reference value: " << pointBasisValues[l] << "\n";
         }
      }
    }

    // Check GRAD of basis functions at points: resize vals to rank-3 container:\n";
    numPoints = tetPoints.dimension(0);
    vals.resize(numFields, numPoints, spaceDim);
    vals.initialize(0.0);
    tetBasis.getValues(vals, tetPoints, OPERATOR_GRAD);

     for (int i = 0; i < numFields; i++) {
       for (int j = 0; j < numPoints; j++) {
         for (int k = 0; k < spaceDim; k++) {
           int l = k + i * spaceDim + j * spaceDim * numFields;
           if (std::abs(vals(i,j,k) - pointBasisGrads[l]) > INTREPID_TOL) {
             errorFlag++;
             *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";

             // Output the multi-index of the value where the error is:
             *outStream << " At multi-index { ";
             *outStream << i << " ";*outStream << j << " ";*outStream << k << " ";
             *outStream << "}  computed grad component: " << vals(i,j,k)
                        << " but reference grad component: " << pointBasisGrads[l] << "\n";
           }
         }
       }
     }
  }
  // Catch unexpected errors
  catch (std::logic_error err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1000;
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
