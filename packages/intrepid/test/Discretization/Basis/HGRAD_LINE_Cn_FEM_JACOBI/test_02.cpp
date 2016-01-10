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

/** \file test_02.cpp
\brief  Patch test for the Intrepid::Basis_HGRAD_LINE_Cn_FEM_JACOBI class.
\author Created by P. Bochev, D. Ridzal, K. Peterson.
*/

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HGRAD_LINE_Cn_FEM_JACOBI.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"

using namespace std;
using namespace Intrepid;

void rhsFunc(FieldContainer<double> &, const FieldContainer<double> &, int);
void neumann(FieldContainer<double> &, const FieldContainer<double> &, const FieldContainer<double> &, int);
void u_exact(FieldContainer<double> &, const FieldContainer<double> &, int);

/// right-hand side function
void rhsFunc(FieldContainer<double> & result, const FieldContainer<double> & points, int degree) {
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
        result(cell,pt) = pow(points(cell,pt,0), degree) - degree*(degree-1)*pow(points(cell,pt,0), degree-2);
      }
    }
  }
}

/// neumann boundary conditions
void neumann(FieldContainer<double> & g_phi, const FieldContainer<double> & phi1, const FieldContainer<double> & phi2, int degree) {
  double g_at_one, g_at_minus_one;
  int num_fields;

  if (degree == 0) {
    g_at_one = 0.0;
    g_at_minus_one = 0.0;
  }
  else {
    g_at_one = degree*pow(1.0, degree-1);
    g_at_minus_one = degree*pow(-1.0, degree-1);
  }

  num_fields = phi1.dimension(0);

  for (int i=0; i<num_fields; i++) {
    g_phi(0,i) = g_at_minus_one*phi1(i,0);
    g_phi(1,i) = g_at_one*phi2(i,0);
  }
}

/// exact solution
void u_exact(FieldContainer<double> & result, const FieldContainer<double> & points, int degree) {
  for (int cell=0; cell<result.dimension(0); cell++) {
    for (int pt=0; pt<result.dimension(1); pt++) {
      result(cell,pt) = pow(points(pt,0), degree);
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
    << "|        for the Neumann problem on a REFERENCE line:                         |\n" \
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
  double zero = 100*INTREPID_TOL;
  outStream -> precision(20);


  try {

    int max_order = 10;  // max total order of polynomial solution

    // Define array containing points at which the solution is evaluated
    int numIntervals = 100;
    int numInterpPoints = numIntervals + 1;
    FieldContainer<double> interp_points(numInterpPoints, 1);
    for (int i=0; i<numInterpPoints; i++) {
      interp_points(i,0) = -1.0+(2.0*(double)i)/(double)numIntervals;
    }
    
    DefaultCubatureFactory<double>  cubFactory;                                   // create factory
    shards::CellTopology line(shards::getCellTopologyData< shards::Line<> >());   // create cell topology

    for (int soln_order=1; soln_order <= max_order; soln_order++) {

      // evaluate exact solution
      FieldContainer<double> exact_solution(1, numInterpPoints);
      u_exact(exact_solution, interp_points, soln_order);

      for (int basis_order=soln_order; basis_order <= max_order; basis_order++) {

        //create basis
        Teuchos::RCP<Basis<double,FieldContainer<double> > > lineBasis =
          Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM_JACOBI<double,FieldContainer<double> >(basis_order) );
        int numFields = lineBasis->getCardinality();

        // create cubature
        Teuchos::RCP<Cubature<double> > lineCub = cubFactory.create(line, 2*basis_order-2);
        int numCubPoints = lineCub->getNumPoints();
        int spaceDim = lineCub->getDimension();

        /* Computational arrays. */
        FieldContainer<double> cub_points(numCubPoints, spaceDim);
        FieldContainer<double> cub_points_physical(1, numCubPoints, spaceDim);
        FieldContainer<double> cub_weights(numCubPoints);
        FieldContainer<double> cell_nodes(1, 2, spaceDim);
        FieldContainer<double> jacobian(1, numCubPoints, spaceDim, spaceDim);
        FieldContainer<double> jacobian_inv(1, numCubPoints, spaceDim, spaceDim);
        FieldContainer<double> jacobian_det(1, numCubPoints);
        FieldContainer<double> weighted_measure(1, numCubPoints);

        FieldContainer<double> value_of_basis_at_cub_points(numFields, numCubPoints);
        FieldContainer<double> transformed_value_of_basis_at_cub_points(1, numFields, numCubPoints);
        FieldContainer<double> weighted_transformed_value_of_basis_at_cub_points(1, numFields, numCubPoints);
        FieldContainer<double> grad_of_basis_at_cub_points(numFields, numCubPoints, spaceDim);
        FieldContainer<double> transformed_grad_of_basis_at_cub_points(1, numFields, numCubPoints, spaceDim);
        FieldContainer<double> weighted_transformed_grad_of_basis_at_cub_points(1, numFields, numCubPoints, spaceDim);
        FieldContainer<double> fe_matrix(1, numFields, numFields);

        FieldContainer<double> rhs_at_cub_points_physical(1, numCubPoints);
        FieldContainer<double> rhs_and_soln_vector(1, numFields);

        FieldContainer<double> one_point(1, 1);
        FieldContainer<double> value_of_basis_at_one(numFields, 1);
        FieldContainer<double> value_of_basis_at_minusone(numFields, 1);
        FieldContainer<double> bc_neumann(2, numFields);

        FieldContainer<double> value_of_basis_at_interp_points(numFields, numInterpPoints);
        FieldContainer<double> transformed_value_of_basis_at_interp_points(1, numFields, numInterpPoints);
        FieldContainer<double> interpolant(1, numInterpPoints);

        FieldContainer<int> ipiv(numFields);

        /******************* START COMPUTATION ***********************/

        // get cubature points and weights
        lineCub->getCubature(cub_points, cub_weights);

        // fill cell vertex array
        cell_nodes(0, 0, 0) = -1.0;
        cell_nodes(0, 1, 0) = 1.0;

        // compute geometric cell information
        CellTools<double>::setJacobian(jacobian, cub_points, cell_nodes, line);
        CellTools<double>::setJacobianInv(jacobian_inv, jacobian);
        CellTools<double>::setJacobianDet(jacobian_det, jacobian);

        // compute weighted measure
        FunctionSpaceTools::computeCellMeasure<double>(weighted_measure, jacobian_det, cub_weights);

        ///////////////////////////
        // Computing mass matrices:
        // tabulate values of basis functions at (reference) cubature points
        lineBasis->getValues(value_of_basis_at_cub_points, cub_points, OPERATOR_VALUE);

        // transform values of basis functions
        FunctionSpaceTools::HGRADtransformVALUE<double>(transformed_value_of_basis_at_cub_points,
                                                        value_of_basis_at_cub_points);

        // multiply with weighted measure
        FunctionSpaceTools::multiplyMeasure<double>(weighted_transformed_value_of_basis_at_cub_points,
                                                    weighted_measure,
                                                    transformed_value_of_basis_at_cub_points);

        // compute mass matrices
        FunctionSpaceTools::integrate<double>(fe_matrix,
                                              transformed_value_of_basis_at_cub_points,
                                              weighted_transformed_value_of_basis_at_cub_points,
                                              COMP_CPP);
        ///////////////////////////

        ////////////////////////////////
        // Computing stiffness matrices:
        // tabulate gradients of basis functions at (reference) cubature points
        lineBasis->getValues(grad_of_basis_at_cub_points, cub_points, OPERATOR_GRAD);

        // transform gradients of basis functions
        FunctionSpaceTools::HGRADtransformGRAD<double>(transformed_grad_of_basis_at_cub_points,
                                                       jacobian_inv,
                                                       grad_of_basis_at_cub_points);

        // multiply with weighted measure
        FunctionSpaceTools::multiplyMeasure<double>(weighted_transformed_grad_of_basis_at_cub_points,
                                                    weighted_measure,
                                                    transformed_grad_of_basis_at_cub_points);

        // compute stiffness matrices and sum into fe_matrix
        FunctionSpaceTools::integrate<double>(fe_matrix,
                                              transformed_grad_of_basis_at_cub_points,
                                              weighted_transformed_grad_of_basis_at_cub_points,
                                              COMP_CPP,
                                              true);
        ////////////////////////////////

        ///////////////////////////////
        // Computing RHS contributions:
        // map (reference) cubature points to physical space
        CellTools<double>::mapToPhysicalFrame(cub_points_physical, cub_points, cell_nodes, line);

        // evaluate rhs function
        rhsFunc(rhs_at_cub_points_physical, cub_points_physical, soln_order);

        // compute rhs
        FunctionSpaceTools::integrate<double>(rhs_and_soln_vector,
                                              rhs_at_cub_points_physical,
                                              weighted_transformed_value_of_basis_at_cub_points,
                                              COMP_CPP);

        // compute neumann b.c. contributions and adjust rhs
        one_point(0,0) = 1.0;   lineBasis->getValues(value_of_basis_at_one, one_point, OPERATOR_VALUE);
        one_point(0,0) = -1.0;  lineBasis->getValues(value_of_basis_at_minusone, one_point, OPERATOR_VALUE);
        neumann(bc_neumann, value_of_basis_at_minusone, value_of_basis_at_one, soln_order);
        for (int i=0; i<numFields; i++) {
          rhs_and_soln_vector(0, i) -= bc_neumann(0, i);
          rhs_and_soln_vector(0, i) += bc_neumann(1, i);
        }
        ///////////////////////////////

        /////////////////////////////
        // Solution of linear system:
        int info = 0;
        Teuchos::LAPACK<int, double> solver;
        //solver.GESV(numRows, 1, &fe_mat(0,0), numRows, &ipiv(0), &fe_vec(0), numRows, &info);
        solver.GESV(numFields, 1, &fe_matrix[0], numFields, &ipiv(0), &rhs_and_soln_vector[0], numFields, &info);
        /////////////////////////////

        ////////////////////////
        // Building interpolant:
        // evaluate basis at interpolation points
        lineBasis->getValues(value_of_basis_at_interp_points, interp_points, OPERATOR_VALUE);
        // transform values of basis functions
        FunctionSpaceTools::HGRADtransformVALUE<double>(transformed_value_of_basis_at_interp_points,
                                                        value_of_basis_at_interp_points);
        FunctionSpaceTools::evaluate<double>(interpolant, rhs_and_soln_vector, transformed_value_of_basis_at_interp_points);
        ////////////////////////

        /******************* END COMPUTATION ***********************/
      
        RealSpaceTools<double>::subtract(interpolant, exact_solution);

        *outStream << "\nNorm-2 difference between exact solution polynomial of order "
                   << soln_order << " and finite element interpolant of order " << basis_order << ": "
                   << RealSpaceTools<double>::vectorNorm(&interpolant[0], interpolant.dimension(1), NORM_TWO) << "\n";

        if (RealSpaceTools<double>::vectorNorm(&interpolant[0], interpolant.dimension(1), NORM_TWO) > zero) {
          *outStream << "\n\nPatch test failed for solution polynomial order "
                     << soln_order << " and basis order " << basis_order << "\n\n";
          errorFlag++;
        }

      } // end for basis_order

   } // end for soln_order

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
