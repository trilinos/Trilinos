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
    \brief  Patch test for the Intrepid::Basis_HDIV_TRI_In_FEM_ORTH class.
    \author Created by P. Bochev, R. Kirby, D. Ridzal, K. Peterson.
*/

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HGRAD_TRI_Cn_FEM_ORTH.hpp"
#include "Intrepid_HDIV_TRI_In_FEM.hpp"
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

void rhsFunc( FieldContainer<double> &, const FieldContainer<double> &, int, int );
void u_exact( FieldContainer<double> &, const FieldContainer<double> &, int, int );

// This is the rhs for (div tau,w) = (f,w),
// which makes f the negative Laplacian of scalar solution
void rhsFunc( FieldContainer<double> &result, 
              const FieldContainer<double> &points,
              int xd,
              int yd )
{
  for (int cell=0;cell<result.dimension(0);cell++) {
    for (int pt=0;pt<result.dimension(1);pt++) {
      result(cell,pt) = 0.0;
      if (xd >=2) {
        result(cell,pt) += xd*(xd-1)*pow(points(cell,pt,0),xd-2)*pow(points(cell,pt,1),yd);
      }
      if (yd >=2) {
        result(cell,pt) += yd*(yd-1)*pow(points(cell,pt,0),xd)*pow(points(cell,pt,1),yd-2);
      }
    }
  }
}

void u_exact( FieldContainer<double> &result, 
              const FieldContainer<double> &points,
              int xd,
              int yd)
{
  for (int cell=0;cell<result.dimension(0);cell++){
    for (int pt=0;pt<result.dimension(1);pt++) {
      result(cell,pt) = std::pow(points(cell,pt,0),xd)*std::pow(points(cell,pt,1),yd);
    }
  }
  return;
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
    << "|                Unit Test (Basis_HGRAD_TRI_Cn_FEM_ORTH)                      |\n" \
    << "|                                                                             |\n" \
    << "|     1) Patch test involving H(div) matrices                                 |\n" \
    << "|        for the Dirichlet problem on a triangular patch                      |\n" \
    << "|        Omega with boundary Gamma.                                           |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Robert Kirby  (robert.c.kirby@ttu.edu),                |\n" \
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

  outStream -> precision(16);

  try {
    DefaultCubatureFactory<double>  cubFactory;                                      // create cubature factory
    shards::CellTopology cell(shards::getCellTopologyData< shards::Triangle<> >());  // create parent cell topology
    shards::CellTopology side(shards::getCellTopologyData< shards::Line<> >());      // create relevant subcell (side) topology
    
    int cellDim = cell.getDimension();
    int sideDim = side.getDimension();

    int min_order = 0;
    int max_order = 9;

    int numIntervals = max_order;
    int numInterpPoints = ((numIntervals + 1)*(numIntervals + 2))/2;
    FieldContainer<double> interp_points_ref(numInterpPoints, 2);
    int counter = 0;
    for (int j=0; j<=numIntervals; j++) {
      for (int i=0; i<=numIntervals; i++) {
        if (i <= numIntervals-j) {
          interp_points_ref(counter,0) = i*(1.0/numIntervals);
          interp_points_ref(counter,1) = j*(1.0/numIntervals);
          counter++;
        }
      }
    }

    interp_points_ref.resize(numInterpPoints,2);

    for (int basis_order=min_order;basis_order<=max_order;basis_order++) {
      // create bases
      Teuchos::RCP<Basis<double,FieldContainer<double> > > vectorBasis =
        Teuchos::rcp(new Basis_HDIV_TRI_In_FEM<double,FieldContainer<double> >(basis_order+1,POINTTYPE_EQUISPACED) );
      Teuchos::RCP<Basis<double,FieldContainer<double> > > scalarBasis =
        Teuchos::rcp(new Basis_HGRAD_TRI_Cn_FEM_ORTH<double,FieldContainer<double> >(basis_order) );

      int numVectorFields = vectorBasis->getCardinality();
      int numScalarFields = scalarBasis->getCardinality();
      int numTotalFields = numVectorFields + numScalarFields;

      // create cubatures
      Teuchos::RCP<Cubature<double> > cellCub = cubFactory.create(cell, 2*(basis_order+1));
      Teuchos::RCP<Cubature<double> > sideCub = cubFactory.create(side, 2*(basis_order+1));

      int numCubPointsCell = cellCub->getNumPoints();
      int numCubPointsSide = sideCub->getNumPoints();

      // hold cubature information
      FieldContainer<double> cub_points_cell(numCubPointsCell, cellDim);
      FieldContainer<double> cub_weights_cell(numCubPointsCell);
      FieldContainer<double> cub_points_side( numCubPointsSide, sideDim );
      FieldContainer<double> cub_weights_side( numCubPointsSide );
      FieldContainer<double> cub_points_side_refcell( numCubPointsSide , cellDim );

      // hold basis function information on refcell
      FieldContainer<double> value_of_v_basis_at_cub_points_cell(numVectorFields, numCubPointsCell, cellDim );
      FieldContainer<double> w_value_of_v_basis_at_cub_points_cell(1, numVectorFields, numCubPointsCell, cellDim);
      FieldContainer<double> div_of_v_basis_at_cub_points_cell( numVectorFields, numCubPointsCell );
      FieldContainer<double> w_div_of_v_basis_at_cub_points_cell( 1, numVectorFields , numCubPointsCell );
      FieldContainer<double> value_of_s_basis_at_cub_points_cell(numScalarFields,numCubPointsCell);
      FieldContainer<double> w_value_of_s_basis_at_cub_points_cell(1,numScalarFields,numCubPointsCell);

      // containers for side integration:
      // I just need the normal component of the vector basis
      // and the exact solution at the cub points
      FieldContainer<double> value_of_v_basis_at_cub_points_side(numVectorFields,numCubPointsSide,cellDim);
      FieldContainer<double> n_of_v_basis_at_cub_points_side(numVectorFields,numCubPointsSide);
      FieldContainer<double> w_n_of_v_basis_at_cub_points_side(1,numVectorFields,numCubPointsSide);
      FieldContainer<double> diri_data_at_cub_points_side(1,numCubPointsSide);
      FieldContainer<double> side_normal(cellDim);

      // holds rhs data
      FieldContainer<double> rhs_at_cub_points_cell(1,numCubPointsCell);

      // FEM matrices and vectors
      FieldContainer<double> fe_matrix_M(1,numVectorFields,numVectorFields);
      FieldContainer<double> fe_matrix_B(1,numVectorFields,numScalarFields);
      FieldContainer<double> fe_matrix(1,numTotalFields,numTotalFields);

      FieldContainer<double> rhs_vector_vec(1,numVectorFields);
      FieldContainer<double> rhs_vector_scal(1,numScalarFields);
      FieldContainer<double> rhs_and_soln_vec(1,numTotalFields);
                        
      FieldContainer<int> ipiv(numTotalFields);
      FieldContainer<double> value_of_s_basis_at_interp_points( numScalarFields , numInterpPoints);
      FieldContainer<double> interpolant( 1 , numInterpPoints );

      // set test tolerance
      double zero = (basis_order+1)*(basis_order+1)*100*INTREPID_TOL;

      // build matrices outside the loop, and then just do the rhs
      // for each iteration

      cellCub->getCubature(cub_points_cell, cub_weights_cell);
      sideCub->getCubature(cub_points_side, cub_weights_side);
      
      // need the vector basis & its divergences
      vectorBasis->getValues(value_of_v_basis_at_cub_points_cell,
                            cub_points_cell,
                            OPERATOR_VALUE);
      vectorBasis->getValues(div_of_v_basis_at_cub_points_cell,
                            cub_points_cell,
                            OPERATOR_DIV);
      
      // need the scalar basis as well
      scalarBasis->getValues(value_of_s_basis_at_cub_points_cell,
                            cub_points_cell,
                            OPERATOR_VALUE);
                                        
      // construct mass matrix      
      cub_weights_cell.resize(1,numCubPointsCell);
      FunctionSpaceTools::multiplyMeasure<double>(w_value_of_v_basis_at_cub_points_cell ,
                                                  cub_weights_cell ,
                                                  value_of_v_basis_at_cub_points_cell ); 
      cub_weights_cell.resize(numCubPointsCell);


      value_of_v_basis_at_cub_points_cell.resize( 1 , numVectorFields , numCubPointsCell , cellDim );
      FunctionSpaceTools::integrate<double>(fe_matrix_M,
                                            w_value_of_v_basis_at_cub_points_cell ,
                                            value_of_v_basis_at_cub_points_cell ,
                                            COMP_BLAS );
      value_of_v_basis_at_cub_points_cell.resize( numVectorFields , numCubPointsCell , cellDim );

      // div matrix
      cub_weights_cell.resize(1,numCubPointsCell);
      FunctionSpaceTools::multiplyMeasure<double>(w_div_of_v_basis_at_cub_points_cell,
                                                  cub_weights_cell,
                                                  div_of_v_basis_at_cub_points_cell);
      cub_weights_cell.resize(numCubPointsCell);
      
      value_of_s_basis_at_cub_points_cell.resize(1,numScalarFields,numCubPointsCell);
      FunctionSpaceTools::integrate<double>(fe_matrix_B,
                                            w_div_of_v_basis_at_cub_points_cell ,
                                            value_of_s_basis_at_cub_points_cell ,
                                            COMP_BLAS );
      value_of_s_basis_at_cub_points_cell.resize(numScalarFields,numCubPointsCell);

      
      // construct div matrix
      
      for (int x_order=0;x_order<=basis_order;x_order++) {
        for (int y_order=0;y_order<=basis_order-x_order;y_order++) {

          // reset global matrix since I destroyed it in LU factorization.
          fe_matrix.initialize();
          // insert mass matrix into global matrix
          for (int i=0;i<numVectorFields;i++) {
            for (int j=0;j<numVectorFields;j++) {
              fe_matrix(0,i,j) = fe_matrix_M(0,i,j);
            }
          }

          // insert div matrix into global matrix
          for (int i=0;i<numVectorFields;i++) {
            for (int j=0;j<numScalarFields;j++) {
              fe_matrix(0,i,numVectorFields+j)=-fe_matrix_B(0,i,j);
              fe_matrix(0,j+numVectorFields,i)=fe_matrix_B(0,i,j);
            }
          }

          // clear old vector data
          rhs_vector_vec.initialize();
          rhs_vector_scal.initialize();
          rhs_and_soln_vec.initialize();
                                        
          // now get rhs vector
          // rhs_vector_scal is just (rhs,w) for w in the scalar basis
          // I already have the scalar basis tabulated.
          cub_points_cell.resize(1,numCubPointsCell,cellDim);
          rhsFunc(rhs_at_cub_points_cell,
                  cub_points_cell,
                  x_order,
                  y_order);

          cub_points_cell.resize(numCubPointsCell,cellDim);

          cub_weights_cell.resize(1,numCubPointsCell);
          FunctionSpaceTools::multiplyMeasure<double>(w_value_of_s_basis_at_cub_points_cell,
                                                      cub_weights_cell,
                                                      value_of_s_basis_at_cub_points_cell);
          cub_weights_cell.resize(numCubPointsCell);
          FunctionSpaceTools::integrate<double>(rhs_vector_scal,
                                                rhs_at_cub_points_cell,
                                                w_value_of_s_basis_at_cub_points_cell,
                                                COMP_BLAS);

          for (int i=0;i<numScalarFields;i++) {
            rhs_and_soln_vec(0,numVectorFields+i) = rhs_vector_scal(0,i);
          }


          // now get <u,v.n> on boundary
          for (unsigned side_cur=0;side_cur<3;side_cur++) {
            // map side cubature to current side
            CellTools<double>::mapToReferenceSubcell( cub_points_side_refcell ,
                                                      cub_points_side ,
                                                      sideDim ,
                                                      (int)side_cur ,
                                                      cell );
                                                                                                
            // Evaluate dirichlet data
            cub_points_side_refcell.resize(1,numCubPointsSide,cellDim);
            u_exact(diri_data_at_cub_points_side,
                    cub_points_side_refcell,x_order,y_order);
            cub_points_side_refcell.resize(numCubPointsSide,cellDim);
                                                                                                                              
            // get normal direction, this has the edge weight factored into it already
            CellTools<double>::getReferenceSideNormal(side_normal , 
                                                      (int)side_cur,cell );

            // v.n at cub points on side
            vectorBasis->getValues(value_of_v_basis_at_cub_points_side ,
                                  cub_points_side_refcell ,
                                  OPERATOR_VALUE );
                                                                                      
                                                        
            for (int i=0;i<numVectorFields;i++) {
              for (int j=0;j<numCubPointsSide;j++) {
                n_of_v_basis_at_cub_points_side(i,j) = 0.0;
                for (int k=0;k<cellDim;k++) {
                  n_of_v_basis_at_cub_points_side(i,j) += side_normal(k) * 
                    value_of_v_basis_at_cub_points_side(i,j,k);
                }
              }
            }
                                                
            cub_weights_side.resize(1,numCubPointsSide);
            FunctionSpaceTools::multiplyMeasure<double>(w_n_of_v_basis_at_cub_points_side,
                                                        cub_weights_side,
                                                        n_of_v_basis_at_cub_points_side);
            cub_weights_side.resize(numCubPointsSide);
                                                
            FunctionSpaceTools::integrate<double>(rhs_vector_vec,
                                                  diri_data_at_cub_points_side,
                                                  w_n_of_v_basis_at_cub_points_side,
                                                  COMP_BLAS,
                                                  false);
            for (int i=0;i<numVectorFields;i++) {
              rhs_and_soln_vec(0,i) -= rhs_vector_vec(0,i);
            }
                                                
          }

          // solve linear system
          int info = 0;
          Teuchos::LAPACK<int, double> solver;
          solver.GESV(numTotalFields, 1, &fe_matrix[0], numTotalFields, &ipiv(0), &rhs_and_soln_vec[0], 
                      numTotalFields, &info);

          // compute interpolant; the scalar entries are last
          scalarBasis->getValues(value_of_s_basis_at_interp_points,
                                interp_points_ref,
                                OPERATOR_VALUE);
          for (int pt=0;pt<numInterpPoints;pt++) {
            interpolant(0,pt)=0.0;
            for (int i=0;i<numScalarFields;i++) {
              interpolant(0,pt) += rhs_and_soln_vec(0,numVectorFields+i)
                * value_of_s_basis_at_interp_points(i,pt);
            }
          }

          interp_points_ref.resize(1,numInterpPoints,cellDim);
          // get exact solution for comparison
          FieldContainer<double> exact_solution(1,numInterpPoints);
          u_exact( exact_solution , interp_points_ref , x_order, y_order);
          interp_points_ref.resize(numInterpPoints,cellDim);

          RealSpaceTools<double>::add(interpolant,exact_solution);

          double nrm= RealSpaceTools<double>::vectorNorm(&interpolant[0],interpolant.dimension(1), NORM_TWO);

          *outStream << "\nNorm-2 error between scalar components of exact solution polynomial of order ("
                     << x_order << ", " << y_order << ") and finite element interpolant of order " << basis_order << ": "
                     << nrm << "\n";

          if (nrm > zero) {
            *outStream << "\n\nPatch test failed for solution polynomial order ("
                       << x_order << ", " << y_order << ") and basis order (scalar / vector)  ("
                       << basis_order << ", " << basis_order + 1 << ")\n\n";
            errorFlag++;
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
