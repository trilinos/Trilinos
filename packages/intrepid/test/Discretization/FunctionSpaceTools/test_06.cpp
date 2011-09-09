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

/** \file test_06.cpp
\brief  Unit test for the FunctionSpaceTools class, testing H-grad, with
        varying types of arrays. Similar to test_01. 
\author Created by D. Ridzal, P. Bochev, and K. Peterson.
*/

#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_HGRAD_TET_C1_FEM.hpp"
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
  << "|     1) basic operator transformations and integration in HGRAD              |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n";


  int errorFlag = 0;
#ifdef HAVE_INTREPID_DEBUG
  int beginThrowNumber = TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + 28;
#endif

  typedef FunctionSpaceTools fst; 

  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 1: exceptions                                                          |\n"\
  << "===============================================================================\n";

  try{
#ifdef HAVE_INTREPID_DEBUG
    FieldContainer<double,  1> a_2(2);
    FieldContainer<double,  2> a_2_2(2, 2);
    FieldContainer<double,  3> a_2_3(2, 3);
    FieldContainer<double,  4> a_3_2(3, 2);
    FieldContainer<double,  5> a_2_2_3(2, 2, 3);
    FieldContainer<double,  6> a_2_2_3_3(2, 2, 3, 3);
    FieldContainer<double,  7> a_2_2_2(2, 2, 2);
    FieldContainer<double,  8> a_2_2_2_3_3(2, 2, 2, 3, 3);
    FieldContainer<double,  9> a_2_2_2_2_2(2, 2, 2, 2, 2);
    FieldContainer<double, 10> a_2_2_2_2(2, 2, 2, 2);
    FieldContainer<double, 11> a_3_2_2_2(3, 2, 2, 2);
    FieldContainer<double, 12> a_2_3_2_2(2, 3, 2, 2);
    FieldContainer<double, 13> a_2_2_3_2(2, 2, 3, 2);
    FieldContainer<double, 14> a_2_2_2_3(2, 2, 2, 3);

    *outStream << "\n >>>>> TESTING computeCellMeasure:\n";
    INTREPID_TEST_COMMAND( fst::computeCellMeasure<double>(a_2_2, a_2, a_2) );
    INTREPID_TEST_COMMAND( fst::computeCellMeasure<double>(a_2_2, a_2_2, a_2) );

    *outStream << "\n >>>>> TESTING computeFaceMeasure:\n";
    INTREPID_TEST_COMMAND( fst::computeFaceMeasure<double>(a_2_2, a_2, a_2, 0, shards::getCellTopologyData< shards::Tetrahedron<> >()) );
    INTREPID_TEST_COMMAND( fst::computeFaceMeasure<double>(a_2_2, a_2_2_3_3, a_2, 0, shards::getCellTopologyData< shards::Tetrahedron<> >()) );

    *outStream << "\n >>>>> TESTING computeEdgeMeasure:\n";
    INTREPID_TEST_COMMAND( fst::computeEdgeMeasure<double>(a_2_2, a_2, a_2, 0, shards::getCellTopologyData< shards::Triangle<> >()) );
    INTREPID_TEST_COMMAND( fst::computeEdgeMeasure<double>(a_2_2, a_2_2_2_2, a_2, 0, shards::getCellTopologyData< shards::Triangle<> >()) );

    *outStream << "\n >>>>> TESTING integrate:\n";
    INTREPID_TEST_COMMAND( fst::integrate<double>(a_2_2_2_2, a_2_2_2, a_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::integrate<double>(a_2, a_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::integrate<double>(a_2, a_2_2_3, a_2_2_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::integrate<double>(a_2, a_2_2_3_3, a_2_2_3_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::integrate<double>(a_2_2, a_2_2, a_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::integrate<double>(a_2_2, a_2_2_3, a_2_2_2_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::integrate<double>(a_2_2, a_2_2_3_3, a_2_2_2_3_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::integrate<double>(a_2_2_2, a_2_2_2, a_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::integrate<double>(a_2_2_2, a_2_2_2_3, a_2_2_2_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::integrate<double>(a_2_2_2, a_2_2_2_3_3, a_2_2_2_3_3, COMP_CPP) );

    *outStream << "\n >>>>> TESTING operatorIntegral:\n";
    INTREPID_TEST_COMMAND( fst::operatorIntegral<double>(a_2_2_2, a_2_2, a_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::operatorIntegral<double>(a_2_2_2, a_2_2_2, a_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::operatorIntegral<double>(a_2_2_2, a_2_2_2_3, a_2_2_2_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::operatorIntegral<double>(a_2_2_2, a_2_2_2_3_3, a_2_2_2_3_3, COMP_CPP) );

    *outStream << "\n >>>>> TESTING functionalIntegral:\n";
    INTREPID_TEST_COMMAND( fst::functionalIntegral<double>(a_2_2, a_2_2_2_3_3, a_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::functionalIntegral<double>(a_2_2, a_2_2, a_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::functionalIntegral<double>(a_2_2, a_2_2_3, a_2_2_2_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::functionalIntegral<double>(a_2_2, a_2_2_3_3, a_2_2_2_3_3, COMP_CPP) );

    *outStream << "\n >>>>> TESTING dataIntegral:\n";
    INTREPID_TEST_COMMAND( fst::dataIntegral<double>(a_2, a_2, a_2_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::dataIntegral<double>(a_2, a_2_2, a_2_2, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::dataIntegral<double>(a_2, a_2_2_3, a_2_2_3, COMP_CPP) );
    INTREPID_TEST_COMMAND( fst::dataIntegral<double>(a_2, a_2_2_3_3, a_2_2_3_3, COMP_CPP) );

    *outStream << "\n >>>>> TESTING applyLeftFieldSigns:\n";
    INTREPID_TEST_COMMAND( fst::applyLeftFieldSigns<double>(a_2, a_2) );
    INTREPID_TEST_COMMAND( fst::applyLeftFieldSigns<double>(a_2_2_2, a_2) );
    INTREPID_TEST_COMMAND( fst::applyLeftFieldSigns<double>(a_2_2_2, a_3_2) );
    INTREPID_TEST_COMMAND( fst::applyLeftFieldSigns<double>(a_2_2_2, a_2_3) );
    INTREPID_TEST_COMMAND( fst::applyLeftFieldSigns<double>(a_2_2_2, a_2_2) );

    *outStream << "\n >>>>> TESTING applyRightFieldSigns:\n";
    INTREPID_TEST_COMMAND( fst::applyRightFieldSigns<double>(a_2, a_2) );
    INTREPID_TEST_COMMAND( fst::applyRightFieldSigns<double>(a_2_2_2, a_2) );
    INTREPID_TEST_COMMAND( fst::applyRightFieldSigns<double>(a_2_2_2, a_3_2) );
    INTREPID_TEST_COMMAND( fst::applyRightFieldSigns<double>(a_2_2_2, a_2_3) );
    INTREPID_TEST_COMMAND( fst::applyRightFieldSigns<double>(a_2_2_2, a_2_2) );

    *outStream << "\n >>>>> TESTING applyFieldSigns:\n";
    INTREPID_TEST_COMMAND( fst::applyFieldSigns<double>(a_2, a_2) );
    INTREPID_TEST_COMMAND( fst::applyFieldSigns<double>(a_2_2, a_2) );
    INTREPID_TEST_COMMAND( fst::applyFieldSigns<double>(a_2_2, a_3_2) );
    INTREPID_TEST_COMMAND( fst::applyFieldSigns<double>(a_2_2, a_2_3) );
    INTREPID_TEST_COMMAND( fst::applyFieldSigns<double>(a_2_2_2_3_3, a_2_2) );

    *outStream << "\n >>>>> TESTING evaluate:\n";
    INTREPID_TEST_COMMAND( fst::evaluate<double>(a_2, a_2, a_2_2) );
    INTREPID_TEST_COMMAND( fst::evaluate<double>(a_2, a_2, a_2_2_2_3_3) );
    INTREPID_TEST_COMMAND( fst::evaluate<double>(a_2, a_2_2, a_2_2_2_3_3) );
    INTREPID_TEST_COMMAND( fst::evaluate<double>(a_2_2_3_3, a_3_2, a_2_2_2_3_3) );
    INTREPID_TEST_COMMAND( fst::evaluate<double>(a_2_2_3_3, a_2_3, a_2_2_2_3_3) );
    INTREPID_TEST_COMMAND( fst::evaluate<double>(a_3_2_2_2, a_2_2, a_2_2_2_2_2) );
    INTREPID_TEST_COMMAND( fst::evaluate<double>(a_2_3_2_2, a_2_2, a_2_2_2_2_2) );
    INTREPID_TEST_COMMAND( fst::evaluate<double>(a_2_2_3_2, a_2_2, a_2_2_2_2_2) );
    INTREPID_TEST_COMMAND( fst::evaluate<double>(a_2_2_2_3, a_2_2, a_2_2_2_2_2) );
    INTREPID_TEST_COMMAND( fst::evaluate<double>(a_2_2_2_2, a_2_2, a_2_2_2_2_2) );
#endif
  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

#ifdef HAVE_INTREPID_DEBUG
  if (TestForException_getThrowNumber() != endThrowNumber)
    errorFlag++;
#endif

  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 2: correctness of math operations                                      |\n"\
  << "===============================================================================\n";

  outStream->precision(20);

  try {
      // cell type: tet
      shards::CellTopology cellType = shards::getCellTopologyData< shards::Tetrahedron<> >();

      /* Related to cubature. */
      // create cubature factory
      DefaultCubatureFactory<double, FieldContainer<double, 1>, FieldContainer<double, 2> > cubFactory;
      // cubature degree
      int cubDegree = 2;
      // create default cubature
      Teuchos::RCP<Cubature<double, FieldContainer<double, 1>, FieldContainer<double, 2> > > myCub = cubFactory.create(cellType, cubDegree);
      // get spatial dimension 
      int spaceDim = myCub->getDimension();
      // get number of cubature points
      int numCubPoints = myCub->getNumPoints();

      /* Related to basis. */
      // create tet basis
      Basis_HGRAD_TET_C1_FEM<double, FieldContainer<double, 1> > tetBasis;
      // get basis cardinality
      int numFields = tetBasis.getCardinality();
 
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

      /* Computational arrays. */
      FieldContainer<double,  1> cub_points(numCubPoints, spaceDim);
      FieldContainer<double,  2> cub_weights(numCubPoints);
      FieldContainer<double,  3> cell_nodes(numCells, numNodes, spaceDim);
      FieldContainer<double,  4> jacobian(numCells, numCubPoints, spaceDim, spaceDim);
      FieldContainer<double,  5> jacobian_inv(numCells, numCubPoints, spaceDim, spaceDim);
      FieldContainer<double,  6> jacobian_det(numCells, numCubPoints);
      FieldContainer<double,  7> weighted_measure(numCells, numCubPoints);

      FieldContainer<double,  1> grad_of_basis_at_cub_points(numFields, numCubPoints, spaceDim);
      FieldContainer<double,  9> transformed_grad_of_basis_at_cub_points(numCells, numFields, numCubPoints, spaceDim);
      FieldContainer<double, 10> weighted_transformed_grad_of_basis_at_cub_points(numCells, numFields, numCubPoints, spaceDim);
      FieldContainer<double, 11> stiffness_matrices(numCells, numFields, numFields);

      FieldContainer<double,  1> value_of_basis_at_cub_points(numFields, numCubPoints);
      FieldContainer<double, 13> transformed_value_of_basis_at_cub_points(numCells, numFields, numCubPoints);
      FieldContainer<double, 14> weighted_transformed_value_of_basis_at_cub_points(numCells, numFields, numCubPoints);
      FieldContainer<double, 15> mass_matrices(numCells, numFields, numFields);

      /******************* START COMPUTATION ***********************/

      // get cubature points and weights
      myCub->getCubature(cub_points, cub_weights);

      // fill cell vertex array
      cell_nodes.setValues(tetnodes, numCellData);

      // compute geometric cell information
      CellTools<double>::setJacobian(jacobian, cub_points, cell_nodes, cellType);
      CellTools<double>::setJacobianInv(jacobian_inv, jacobian);
      CellTools<double>::setJacobianDet(jacobian_det, jacobian);

      // compute weighted measure
      fst::computeCellMeasure<double>(weighted_measure, jacobian_det, cub_weights);

      // Computing stiffness matrices:
      // tabulate gradients of basis functions at (reference) cubature points
      tetBasis.getValues(grad_of_basis_at_cub_points, cub_points, OPERATOR_GRAD);

      // transform gradients of basis functions
      fst::HGRADtransformGRAD<double>(transformed_grad_of_basis_at_cub_points,
                                      jacobian_inv,
                                      grad_of_basis_at_cub_points);

      // multiply with weighted measure
      fst::multiplyMeasure<double>(weighted_transformed_grad_of_basis_at_cub_points,
                                   weighted_measure,
                                   transformed_grad_of_basis_at_cub_points);

      // compute stiffness matrices
      fst::integrate<double>(stiffness_matrices,
                             transformed_grad_of_basis_at_cub_points,
                             weighted_transformed_grad_of_basis_at_cub_points,
                             COMP_CPP);


      // Computing mass matrices:
      // tabulate values of basis functions at (reference) cubature points
      tetBasis.getValues(value_of_basis_at_cub_points, cub_points, OPERATOR_VALUE);

      // transform values of basis functions
      fst::HGRADtransformVALUE<double>(transformed_value_of_basis_at_cub_points,
                                       value_of_basis_at_cub_points);

      // multiply with weighted measure
      fst::multiplyMeasure<double>(weighted_transformed_value_of_basis_at_cub_points,
                                   weighted_measure,
                                   transformed_value_of_basis_at_cub_points);

      // compute mass matrices
      fst::integrate<double>(mass_matrices,
                             transformed_value_of_basis_at_cub_points,
                             weighted_transformed_value_of_basis_at_cub_points,
                             COMP_CPP);

      /*******************  STOP COMPUTATION ***********************/


      /******************* START COMPARISON ***********************/
      string basedir = "./testdata";
      for (int cell_id = 0; cell_id < numCells; cell_id++) {

        stringstream namestream;
        string filename;
        namestream <<  basedir << "/mass_TET_FEM_P1" << "_" << "0" << cell_id+1 << ".dat";
        namestream >> filename;

        ifstream massfile(&filename[0]);
        if (massfile.is_open()) {
          if (compareToAnalytic<double>(&mass_matrices(cell_id, 0, 0), massfile, 1e-10, 0) > 0)
            errorFlag++;
          massfile.close();
        }
        else {
          errorFlag = -1;
          std::cout << "End Result: TEST FAILED\n";
          return errorFlag;
        }

        namestream.clear();
        namestream << basedir << "/stiff_TET_FEM_P1" << "_" << "0" << cell_id+1 << ".dat";
        namestream >> filename;

        ifstream stifffile(&filename[0]);
        if (stifffile.is_open())
        {
          if (compareToAnalytic<double>(&stiffness_matrices(cell_id, 0, 0), stifffile, 1e-10, 0) > 0)
            errorFlag++;
          stifffile.close();
        }
        else {
          errorFlag = -1;
          std::cout << "End Result: TEST FAILED\n";
          return errorFlag;
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
