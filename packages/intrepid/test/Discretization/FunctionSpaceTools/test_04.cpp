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

/** \file test_04.cpp
\brief  Unit test for the FunctionSpaceTools class, testing dataIntegral
        via integrate.
\author Created by D. Ridzal, P. Bochev, and K. Peterson.
*/

#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace std;
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
  << "|                      Unit Test (FunctionSpaceTools)                         |\n" \
  << "|                                                                             |\n" \
  << "|     1) volume integration on tetrahedra, testing dataIntegral               |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n";


  int errorFlag = 0;

  typedef FunctionSpaceTools fst; 

  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 1: correctness of cell volumes                                         |\n"\
  << "===============================================================================\n";

  outStream->precision(20);

  try {
      shards::CellTopology cellType = shards::getCellTopologyData< shards::Tetrahedron<> >();   // cell type: tet

      /* Related to cubature. */
      DefaultCubatureFactory<double> cubFactory;                                                // create cubature factory
      int cubDegree = 0;                                                                        // cubature degree
      Teuchos::RCP<Cubature<double> > myCub = cubFactory.create(cellType, cubDegree);           // create default cubature
      int spaceDim = myCub->getDimension();                                                     // get spatial dimension 
      int numCubPoints = myCub->getNumPoints();                                                 // get number of cubature points

      /* Cell geometries. */
      int numCells    = 4;
      int numNodes    = 4;
      int numCellData = numCells*numNodes*spaceDim;
      double tetnodes[] = {
        // tet 0
        0.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        // tet 1
        4.0, 5.0, 1.0,
       -6.0, 2.0, 0.0,
        4.0, -3.0, -1.0,
        0.0, 2.0, 5.0,
        // tet 2
        -6.0, -3.0, 1.0,
        9.0, 2.0, 1.0,
        8.9, 2.1, 0.9,
        8.9, 2.1, 1.1,
        // tet 3
        -6.0, -3.0, 1.0,
        12.0, 3.0, 1.0,
        2.9, 0.1, 0.9,
        2.9, 0.1, 1.1
      };

      /* Analytic volumes. */
      double tetvols[] = {1.0/6.0, 194.0/3.0, 1.0/15.0, 2.0/25.0};

      /* Computational arrays. */
      FieldContainer<double> cub_points(numCubPoints, spaceDim);
      FieldContainer<double> cub_weights(numCubPoints);
      FieldContainer<double> cell_nodes(numCells, numNodes, spaceDim);
      FieldContainer<double> jacobian(numCells, numCubPoints, spaceDim, spaceDim);
      FieldContainer<double> jacobian_det(numCells, numCubPoints);
      FieldContainer<double> weighted_measure(numCells, numCubPoints);
      FieldContainer<double> data_one(numCells, numCubPoints);
      FieldContainer<double> volumes(numCells);

      /******************* START COMPUTATION ***********************/

      // get cubature points and weights
      myCub->getCubature(cub_points, cub_weights);

      // fill cell vertex array
      cell_nodes.setValues(tetnodes, numCellData);

      // compute geometric cell information
      CellTools<double>::setJacobian(jacobian, cub_points, cell_nodes, cellType);
      CellTools<double>::setJacobianDet(jacobian_det, jacobian);

      // compute weighted measure
      fst::computeCellMeasure<double>(weighted_measure, jacobian_det, cub_weights);

      // set data to 1.0
      for (int cell=0; cell<data_one.dimension(0); cell++) { 
        for (int qp=0; qp<data_one.dimension(1); qp++) {
          data_one(cell,qp) = 1.0;
        }
      } 

      // compute volumes
      fst::integrate<double>(volumes, data_one, weighted_measure, COMP_CPP);

      /******************* STOP COMPUTATION ***********************/


      /******************* START COMPARISON ***********************/
      for (int cell_id = 0; cell_id < numCells; cell_id++) {
        *outStream << "Volume of cell " << cell_id << " = " << volumes(cell_id) << "    vs.    Analytic value =  " << tetvols[cell_id] << "\n";
        if (std::fabs(volumes(cell_id)-tetvols[cell_id]) > INTREPID_TOL) {
          errorFlag++;
        }
      }
      /******************* STOP COMPARISON ***********************/

      *outStream << "\n";
  }
  catch (std::logic_error err) {
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
