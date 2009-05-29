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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov), 
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file test_02.cpp
\brief  Patch test for the Intrepid::Basis_HGRAD_LINE_Cn_FEM_JACOBI class.
\author Created by P. Bochev, D. Ridzal, K. Peterson.
*/

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HGRAD_LINE_Cn_FEM_JACOBI.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace std;
using namespace Intrepid;

void rhsFunc(FieldContainer<double>, const FieldContainer<double>, int);
void u_exact(FieldContainer<double>, const FieldContainer<double>, int);

void rhsFunc(FieldContainer<double> result, const FieldContainer<double> points, int degree) {
  if (degree == 0) {
    for (int cell=0; cell<result.dimension(0); cell++) {
      for (int pt=0; pt<result.dimension(1); pt++) {
        result(cell,pt) = 1.0;
      }
    }
  }
  else if (degree == 1) {
    for (int cell=0; cell<result.dimension(0); cell++) {
      for (int pt=0; pt<result.dimension(1); pt++) {
        result(cell,pt) = points(cell,pt,0);
      }
    }
  }
  else {
    for (int cell=0; cell<result.dimension(0); cell++) {
      for (int pt=0; pt<result.dimension(1); pt++) {
        result(cell,pt) = pow(points(cell,pt,0), degree) + degree*(degree-1)*pow(points(cell,pt,0), degree-2);
      }
    }
  }
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
    << "|               Unit Test (Basis_HGRAD_LINE_Cn_FEM_JACOBI)                    |\n" \
    << "|                                                                             |\n" \
    << "|     1) Patch test involving mass and stiffness matrices,                    |\n" \
    << "|        for the Neumann problem                                              |\n" \
    << "|                                                                             |\n" \
    << "|            - u'' + u = f  in (-1,1),  u' = g at -1,1                        |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n"\
    << "| TEST 1: Patch test                                                          |\n"\
    << "===============================================================================\n";

  
  int errorFlag = 0;
  outStream -> precision(20);


  try {

    int maxorder = 10;

    // Define computational cells
    int numCells = 1;
    int numNodes = 2;

    // Define array containing points at which the solution is evaluated
    int numIntervals = 100;
    FieldContainer<double> evalPoints(numIntervals+1, 1);
    for (int i=0; i<numIntervals+1; i++) {
      evalPoints(i,0) = -1.0+(2.0*(double)i)/(double)numIntervals;
    }
    
    DefaultCubatureFactory<double>  cubFactory;                                   // create factory
    shards::CellTopology line(shards::getCellTopologyData< shards::Line<> >());   // create cell topology

    for (int ordi=0; ordi < maxorder; ordi++) {
      //create basis
      Teuchos::RCP<Basis<double,FieldContainer<double> > > lineBasis =
        Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM_JACOBI<double,FieldContainer<double> >(ordi) );
      int numFields = lineBasis->getCardinality();

      // create cubature
      Teuchos::RCP<Cubature<double> > lineCub = cubFactory.create(line, 2*ordi);
      int numCubPoints = lineCub->getNumPoints();
      int spaceDim = lineCub->getDimension();

      /* Computational arrays. */
      FieldContainer<double> cub_points(numCubPoints, spaceDim);
      FieldContainer<double> cub_weights(numCubPoints);
      FieldContainer<double> cell_nodes(numCells, numNodes, spaceDim);
      FieldContainer<double> jacobian(numCells, numCubPoints, spaceDim, spaceDim);
      FieldContainer<double> jacobian_inv(numCells, numCubPoints, spaceDim, spaceDim);
      FieldContainer<double> jacobian_det(numCells, numCubPoints);
      FieldContainer<double> weighted_measure(numCells, numCubPoints);

      FieldContainer<double> value_of_basis_at_cub_points(numFields, numCubPoints);
      FieldContainer<double> transformed_value_of_basis_at_cub_points(numCells, numFields, numCubPoints);
      FieldContainer<double> weighted_transformed_value_of_basis_at_cub_points(numCells, numFields, numCubPoints);
      FieldContainer<double> mass_matrices(numCells, numFields, numFields);

      FieldContainer<double> grad_of_basis_at_cub_points(numFields, numCubPoints, spaceDim);
      FieldContainer<double> transformed_grad_of_basis_at_cub_points(numCells, numFields, numCubPoints, spaceDim);
      FieldContainer<double> weighted_transformed_grad_of_basis_at_cub_points(numCells, numFields, numCubPoints, spaceDim);
      FieldContainer<double> stiffness_matrices(numCells, numFields, numFields);


      /******************* START COMPUTATION ***********************/

      // get cubature points and weights
      lineCub->getCubature(cub_points, cub_weights);

      // fill cell vertex array
      for (int i=0; i<numCells; i++) {
        for (int j=0; j<numNodes; j++) {
          cell_nodes(i,j,0) = -1.0+(2.0*(double)(i+j))/(double)numCells;
        }
      }

      // compute geometric cell information
      CellTools<double>::setJacobian(jacobian, cub_points, cell_nodes, line);
      CellTools<double>::setJacobianInv(jacobian_inv, jacobian);
      CellTools<double>::setJacobianDet(jacobian_det, jacobian);

      // compute weighted measure
      FunctionSpaceTools::computeMeasure<double>(weighted_measure, jacobian_det, cub_weights);

      // Computing mass matrices:
      // tabulate values of basis functions at (reference) cubature points
      lineBasis->getValues(value_of_basis_at_cub_points, cub_points, OPERATOR_VALUE);

      // transform gradients of basis functions into physical space
      FunctionSpaceTools::HGRADtransformVALUE<double>(transformed_value_of_basis_at_cub_points,
                                                      value_of_basis_at_cub_points);

      // multiply with weighted measure
      FunctionSpaceTools::multiplyMeasure<double>(weighted_transformed_value_of_basis_at_cub_points,
                                                  weighted_measure,
                                                  transformed_value_of_basis_at_cub_points);

      // compute mass matrices
      FunctionSpaceTools::integrate<double>(mass_matrices,
                                            transformed_value_of_basis_at_cub_points,
                                            weighted_transformed_value_of_basis_at_cub_points,
                                            COMP_CPP);

      // Computing stiffness matrices:
      // tabulate gradients of basis functions at (reference) cubature points
      lineBasis->getValues(grad_of_basis_at_cub_points, cub_points, OPERATOR_GRAD);

      // transform gradients of basis functions into physical space
      FunctionSpaceTools::HGRADtransformGRAD<double>(transformed_grad_of_basis_at_cub_points,
                                                     jacobian_inv,
                                                     grad_of_basis_at_cub_points);

      // multiply with weighted measure
      FunctionSpaceTools::multiplyMeasure<double>(weighted_transformed_grad_of_basis_at_cub_points,
                                                  weighted_measure,
                                                  transformed_grad_of_basis_at_cub_points);

      // compute stiffness matrices
      FunctionSpaceTools::integrate<double>(stiffness_matrices,
                                            transformed_grad_of_basis_at_cub_points,
                                            weighted_transformed_grad_of_basis_at_cub_points,
                                            COMP_CPP);
    } // end for
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

  return errorFlag;
}
