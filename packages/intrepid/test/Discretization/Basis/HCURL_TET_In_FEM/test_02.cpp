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

/** \file   test_02.cpp
    \brief  Mass matrix patch test for the Intrepid::Basis_HCURL_TET_In_FEM class.
    \author Created by P. Bochev, R. Kirby, D. Ridzal, K. Peterson.
*/

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HCURL_TET_In_FEM.hpp"
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

void rhsFunc( FieldContainer<double> &, const FieldContainer<double> &, int, int, int, int );
void u_exact( FieldContainer<double> &, const FieldContainer<double> &, int , int, int, int );

void u_exact( FieldContainer<double> &result, 
              const FieldContainer<double> &points,
	      int comp, 
              int xd,
              int yd,
              int zd)
{
  for (int cell=0;cell<result.dimension(0);cell++){
    for (int pt=0;pt<result.dimension(1);pt++) {
      result(cell,pt,comp) = std::pow(points(cell,pt,0),xd)*std::pow(points(cell,pt,1),yd)
        *std::pow(points(cell,pt,2),zd);
    }
  }
  return;
}

void rhsFunc( FieldContainer<double> & result , 
	      const FieldContainer<double> &points ,
	      int comp,
	      int xd,
	      int yd,
	      int zd )
{
  u_exact( result , points , comp , xd , yd , zd );
}

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  
  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint = argc - 1;
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
    << "|                  Unit Test (Basis_HCURL_TET_In_FEM)                         |\n" \
    << "|                                                                             |\n" \
    << "| 1) Patch test involving H(curl) matrices                                    |\n" \
    << "|                                                                             |\n" \
    << "|   Questions? Contact Pavel Bochev (pbboche@sandia.gov),                     |\n" \
    << "|                      Robert Kirby (robert.c.kirby@ttu.edu),                 |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov),                     |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n" \
    << "| TEST 2: Patch test for mass matrices                                        |\n" \
    << "===============================================================================\n";
  
  
  int errorFlag = 0;
  
  outStream -> precision(16);
  
  try {
    DefaultCubatureFactory<double> cubFactory;                                          // create cubature factory
    shards::CellTopology cell(shards::getCellTopologyData< shards::Tetrahedron<> >());  // create parent cell topology
    
    int cellDim = cell.getDimension();
    
    int min_order = 1;
    int max_order = 5;
    
    int numIntervals = max_order;
    int numInterpPoints = ((numIntervals + 1)*(numIntervals + 2)*(numIntervals+3))/6;
    FieldContainer<double> interp_points_ref(numInterpPoints, cellDim);
    int counter = 0;
    for (int j=0; j<=numIntervals; j++) {
      for (int i=0; i<=numIntervals-j; i++) {
        for (int k=0;k<numIntervals-j-i;k++) {
          interp_points_ref(counter,0) = i*(1.0/numIntervals);
          interp_points_ref(counter,1) = j*(1.0/numIntervals);
          interp_points_ref(counter,2) = k*(1.0/numIntervals);
          counter++;
        }
      }
    }
    
    for (int basis_order=min_order;basis_order<=max_order;basis_order++) {
      // create basis
      Teuchos::RCP<Basis<double,FieldContainer<double> > > basis =
        Teuchos::rcp(new Basis_HCURL_TET_In_FEM<double,FieldContainer<double> >(basis_order,POINTTYPE_EQUISPACED) );
      
      int numFields = basis->getCardinality();
      
      // create cubatures
      Teuchos::RCP<Cubature<double> > cellCub = cubFactory.create(cell, 2*(basis_order+1));
      
      int numCubPointsCell = cellCub->getNumPoints();
      
      // hold cubature information
      FieldContainer<double> cub_points_cell(numCubPointsCell, cellDim);
      FieldContainer<double> cub_weights_cell(numCubPointsCell);
      
      // hold basis function information on refcell
      FieldContainer<double> value_of_basis_at_cub_points_cell(numFields, numCubPointsCell, cellDim );
      FieldContainer<double> w_value_of_basis_at_cub_points_cell(1, numFields, numCubPointsCell, cellDim);

      // holds rhs data
      FieldContainer<double> rhs_at_cub_points_cell(1,numCubPointsCell,cellDim);
      
      // FEM mass matrix
      FieldContainer<double> fe_matrix_bak(1,numFields,numFields);
      FieldContainer<double> fe_matrix(1,numFields,numFields);
      FieldContainer<double> rhs_and_soln_vec(1,numFields);
      
      FieldContainer<int> ipiv(numFields);
      FieldContainer<double> value_of_basis_at_interp_points( numFields , numInterpPoints , cellDim);
      FieldContainer<double> interpolant( 1, numInterpPoints , cellDim );

      int info = 0;
      Teuchos::LAPACK<int, double> solver;
      
      // set test tolerance
      double zero = (basis_order+1)*(basis_order+1)*1000*INTREPID_TOL;
      
      // build matrices outside the loop, and then just do the rhs
      // for each iteration
      cellCub->getCubature(cub_points_cell, cub_weights_cell);
      
      // need the vector basis
      basis->getValues(value_of_basis_at_cub_points_cell,
		       cub_points_cell,
		       OPERATOR_VALUE);
      basis->getValues( value_of_basis_at_interp_points ,
			interp_points_ref ,
			OPERATOR_VALUE );
			
			


      // construct mass matrix
      cub_weights_cell.resize(1,numCubPointsCell);
      FunctionSpaceTools::multiplyMeasure<double>(w_value_of_basis_at_cub_points_cell ,
                                                  cub_weights_cell ,
                                                  value_of_basis_at_cub_points_cell ); 
      cub_weights_cell.resize(numCubPointsCell);
      
      
      value_of_basis_at_cub_points_cell.resize( 1 , numFields , numCubPointsCell , cellDim );
      FunctionSpaceTools::integrate<double>(fe_matrix_bak,
                                            w_value_of_basis_at_cub_points_cell ,
                                            value_of_basis_at_cub_points_cell ,
                                            COMP_BLAS );
      value_of_basis_at_cub_points_cell.resize( numFields , numCubPointsCell , cellDim );
      

      //std::cout << fe_matrix_bak << std::endl;

      for (int x_order=0;x_order<basis_order;x_order++) {
        for (int y_order=0;y_order<basis_order-x_order;y_order++) {
          for (int z_order=0;z_order<basis_order-x_order-y_order;z_order++) {
	    for (int comp=0;comp<cellDim;comp++) {
	      fe_matrix.initialize();
	      // copy mass matrix 
	      for (int i=0;i<numFields;i++) {
		for (int j=0;j<numFields;j++) {
		  fe_matrix(0,i,j) = fe_matrix_bak(0,i,j);
		}
	      }
	      
	      // clear old vector data
	      rhs_and_soln_vec.initialize();
	      
	      // now get rhs vector
	      
	      cub_points_cell.resize(1,numCubPointsCell,cellDim);
	     
	      rhs_at_cub_points_cell.initialize();
	      rhsFunc(rhs_at_cub_points_cell,
		      cub_points_cell,
		      comp, 
		      x_order,
		      y_order,
		      z_order);
            
	      cub_points_cell.resize(numCubPointsCell,cellDim);

	      cub_weights_cell.resize(numCubPointsCell);

	      FunctionSpaceTools::integrate<double>(rhs_and_soln_vec,
						    rhs_at_cub_points_cell,
						    w_value_of_basis_at_cub_points_cell,
						    COMP_BLAS);
	      
	      // solve linear system

// 	      solver.GESV(numFields, 1, &fe_matrix[0], numFields, &ipiv(0), &rhs_and_soln_vec[0], 
// 			  numFields, &info);
	      solver.POTRF('L',numFields,&fe_matrix[0],numFields,&info);
	      solver.POTRS('L',numFields,1,&fe_matrix[0],numFields,&rhs_and_soln_vec[0],numFields,&info);
	      
	      interp_points_ref.resize(1,numInterpPoints,cellDim);
	      // get exact solution for comparison
	      FieldContainer<double> exact_solution(1,numInterpPoints,cellDim);
	      exact_solution.initialize();
	      u_exact( exact_solution , interp_points_ref , comp , x_order, y_order, z_order);
	      interp_points_ref.resize(numInterpPoints,cellDim);

	      // compute interpolant
	      // first evaluate basis at interpolation points
	      value_of_basis_at_interp_points.resize(1,numFields,numInterpPoints,cellDim);
	      FunctionSpaceTools::evaluate<double>( interpolant , 
						    rhs_and_soln_vec ,
						    value_of_basis_at_interp_points );
	      value_of_basis_at_interp_points.resize(numFields,numInterpPoints,cellDim);
	      
	      RealSpaceTools<double>::subtract(interpolant,exact_solution);
	      
	      double nrm= RealSpaceTools<double>::vectorNorm(&interpolant[0],interpolant.dimension(1), NORM_TWO);
	      
	      *outStream << "\nNorm-2 error between scalar components of exact solution of order ("
			 << x_order << ", " << y_order << ", " << z_order
			 << ") in component " << comp 
			 << " and finite element interpolant of order " << basis_order << ": "
			 << nrm << "\n";
	      
	      if (nrm > zero) {
		*outStream << "\n\nPatch test failed for solution polynomial order ("
			   << x_order << ", " << y_order << ", " << z_order << ") and basis order (scalar, vector)  ("
			   << basis_order << ", " << basis_order+1 << ")\n\n";
		errorFlag++;
	      }
            }
          }
        }
      }
    }
    
  }
  
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
