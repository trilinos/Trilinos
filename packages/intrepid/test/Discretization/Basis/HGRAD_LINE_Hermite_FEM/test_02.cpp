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
\brief  Unit tests for the Intrepid::Basis_HGRAD_LINE_Hermite_FEM class.
\author Created by G. von Winckel.
*/

// Intrepid Includes
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_HGRAD_LINE_Hermite_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_PointTools.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"

// Teuchos Includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_SerialDenseVector.hpp"

#include <iomanip>

using namespace std;
using namespace Intrepid;


int main(int argc, char *argv[]) {

  using FC  = FieldContainer<double>;
  using CT  = CellTools<double>;
  using FST = FunctionSpaceTools;
  using AT  = ArrayTools;

  using Teuchos::RCP;
  using Teuchos::rcpFromRef;
  using Vector = Teuchos::SerialDenseVector<int,double>;
  using Matrix = Teuchos::SerialDenseMatrix<int,double>;
  using Solver = Teuchos::SerialDenseSolver<int,double>;

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
    << "=============================================================================\n" \
    << "|                                                                           |\n" \
    << "|                Unit Test (Basis_HGRAD_LINE_Hermite_FEM)                   |\n" \
    << "|                                                                           |\n" \
    << "|   Solve the cantilevered Euler-Bernoulli static beam equation with unit   |\n" \
    << "|   second moment of area and verify using a manufactured solution.         |\n" \
    << "|                                                                           |\n" \
    << "|                          D^2[E(x) D^2 w(x) = q(x)                         |\n" \
    << "|                                                                           |\n" \
    << "|   with clamped boundary conditions    w(0) = 0,       w'(0) = 0           |\n" \
    << "|   stress boundary condition           E(1)w\"(1)=-6                        |\n" \
    << "|   and shear force boundary condition  [Ew\"]'(1)=-6                        |\n" \
    << "|                                                                           |\n" \
    << "|   The exact deflection is w(x) = 3x^2-2*x^3                               |\n" \
    << "|   The elastic modulus is  E(x) = 2-x                                      |\n" \
    << "|   The forcing term is     q(x) = 24                                       |\n" \
    << "|                                                                           |\n" \
    << "|   Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                 |\n" \
    << "|                       Denis Ridzal  (dridzal@sandia.gov),                 |\n" \
    << "|                       Kara Peterson (kjpeter@sandia.gov).                 |\n" \
    << "|                                                                           |\n" \
    << "|   Intrepid's website: http://trilinos.sandia.gov/packages/intrepid        |\n" \
    << "|   Trilinos website:   http://trilinos.sandia.gov                          |\n" \
    << "|                                                                           |\n" \
    << "=============================================================================\n";

  int errorFlag = 0;


  try {
  
    shards::CellTopology line(shards::getCellTopologyData<shards::Line<2>>());

    DefaultCubatureFactory<double> cubFactory;

    int numCells  = 10;        // Number of cells
    int numVert   = 2;         // Number of vertices per cell
    int cubOrder  = 8;         // Highest order of polynomial to integrate exactly
    int numPts    = 3;         // Number of interpolation points per cell
    int numFields = 2*numPts;  // Number of basis functions per cell
    int spaceDim  = 1;         // Number of spatial dimensions

    double length = 1.0;                   // Computatonal domain length
    double cellLength = length/numCells;

    *outStream << "\n\nDiscretization Details"      << std::endl;
    *outStream << "-------------------------------" << std::endl; 
    *outStream << "Number of cells           = "    << numCells   << std::endl;
    *outStream << "Cubature order            = "    << cubOrder   << std::endl;
    *outStream << "Number of basis functions = "    << numFields  << std::endl;
    *outStream << "Physical cell length      = "    << cellLength << std::endl;

    // Total degrees of freedom 
    // Exclude 2 for the clamped boundary condition at x=0
    // Exclude 2 per cell for value and derivative node condensation
    int numDof = numCells*(numFields-2);

    *outStream << "Total degrees of freedom  = " << numDof << std::endl; 

    FC cellVert(numCells,numVert,spaceDim);  // Elemental end points
 
    // Set cell vertices
    for(int cell=0; cell<numCells; ++cell ) {
      cellVert(cell,0,0) = cell*cellLength;
      cellVert(cell,1,0) = (cell+1)*cellLength;   
    }   
 
   
    /*****************************************/
    /*   CREATE CELL AND PHYSICAL CUBATURE   */
    /*****************************************/

    RCP<Cubature<double>> cellCub = cubFactory.create(line,cubOrder);

    FC cubPts(cubOrder, spaceDim);               // Reference cell cubature points
    FC cubWts(cubOrder);                         // Reference cell cubature weights
    cellCub->getCubature(cubPts,cubWts);
 
    // Determine how many points are used and resize accordingly
    int numCubPts = cellCub->getNumPoints();

    *outStream << "Number of cubature points = " << numCubPts << std::endl; 

    cubPts.resize(numCubPts,spaceDim);
    cubWts.resize(numCubPts);

    FC physCubPts(numCells,numCubPts, spaceDim); // Physical cubature points
    FC wtdMeasure(numCells,numCubPts);          
    
    CellTools<double>::mapToPhysicalFrame(physCubPts,cubPts,cellVert,line);

    *outStream << std::setprecision(5) << std::endl;
    *outStream << "Cell Vertices:" << std::endl;;
    for(int cell=0; cell<numCells; ++cell) {
      *outStream << std::setw(5) << cellVert(cell,0,0); 
    }
    *outStream << std::setw(5) << cellVert(numCells-1,1,0) << std::endl; 
 
    *outStream << "\nReference cell cubature points:" << std::endl;
    for(int pt=0; pt<numCubPts; ++pt) {
      *outStream << std::setw(10) << cubPts(pt,0);
    }
    *outStream << std::endl;
   
    *outStream << "\nReference cell cubature weights:" << std::endl;
    for( int pt=0; pt<numCubPts; ++pt) {
      *outStream << std::setw(10) << cubWts(pt);
    }
    *outStream << std::endl;

    *outStream << "\nPhysical cubature points:\n" << std::endl;
    *outStream << std::setw(7) << "Cell\\Pt | ";
    for(int pt=0; pt<numCubPts; ++pt) {
      *outStream << std::setw(10) << pt;
    }  
    *outStream << std::endl;
    
    *outStream << std::string(10*(1+numCubPts),'-') << std::endl;
    for(int cell=0; cell<numCells; ++cell){
      *outStream << std::setw(7) << cell << " | ";
      for(int pt=0; pt<numCubPts; ++pt) {
        *outStream << std::setw(10) << physCubPts(cell,pt,0);   
      } 
      *outStream << std::endl;
    }


    /********************************************/
    /*   ELASTIC MODULUS AND FORCING FUNCTION   */
    /********************************************/

    FC elasmod(numCells,numCubPts); 
    FC qforce(numCells,numCubPts);

    for(int cell=0; cell<numCells; ++cell) {

      for(int pt=0; pt<numCubPts; ++pt) {
        double x = physCubPts(cell,pt,0);
        elasmod(cell,pt) = 2.0-x; //std::exp(-x);
        qforce(cell,pt)  = 24.0; // 4.0-3.0*x; //6*x;  // (x-2.0)*std::exp(-x);
      }
    }   

    /****************************************/
    /*   CREATE HERMITE INTERPOLANT BASIS   */
    /****************************************/


    FC pts(PointTools::getLatticeSize(line,numPts-1),1);
    PointTools::getLattice<double,FC>(pts,line,numPts-1);

    *outStream << "\nReference cell interpolation points:" << std::endl;
    for(int pt=0; pt<numPts; ++pt) {
      *outStream << std::setw(10) << pts(pt,0);
    }
    *outStream << std::endl;
    

    FC physPts(numCells, numPts, spaceDim);  // Physical interpolation points
    CellTools<double>::mapToPhysicalFrame(physPts,pts,cellVert,line);

    *outStream << "\nPhysical interpolation points:\n" << std::endl;
    *outStream << std::setw(7) << "Cell\\Pt | ";
    for(int pt=0; pt<numPts; ++pt) {
      *outStream << std::setw(10) << pt;
    }  
    *outStream << std::endl;
    
    *outStream << std::string(10*(1+numPts),'-') << std::endl;
    for(int cell=0; cell<numCells; ++cell){
      *outStream << std::setw(7) << cell << " | ";
      for(int pt=0; pt<numPts; ++pt) {
        *outStream << std::setw(10) << physPts(cell,pt,0);   
      } 
      *outStream << std::endl;
    }

    Basis_HGRAD_LINE_Hermite_FEM<double,FC> hermiteBasis(pts);

    FC valsCubPts(numFields,numCubPts);
    FC der2CubPts(numFields,numCubPts,spaceDim);
    
    hermiteBasis.getValues(valsCubPts,cubPts,OPERATOR_VALUE);
    hermiteBasis.getValues(der2CubPts,cubPts,OPERATOR_D2);

    FC jacobian(numCells,numCubPts,spaceDim,spaceDim);
    FC jacInv(numCells,numCubPts,spaceDim,spaceDim);
    FC jacDet(numCells,numCubPts);

    FC tranValsCubPts(numCells,numFields,numCubPts);
    FC tranDer2CubPts(numCells,numFields,numCubPts,spaceDim);
    FC wtdTranValsCubPts(numCells,numFields,numCubPts);
    FC wtdTranDer2CubPts(numCells,numFields,numCubPts,spaceDim);    

    CT::setJacobian(jacobian,cubPts,cellVert,line);
    CT::setJacobianInv(jacInv,jacobian);
    CT::setJacobianDet(jacDet,jacobian);

    FST::computeCellMeasure<double>(wtdMeasure,jacDet,cubWts);
    FST::HGRADtransformVALUE<double>(tranValsCubPts,valsCubPts); 

    // There is no predefined transform for second derivatives
    // Note that applying the Jacobian inverse twice is only valid because of the
    // affine mapping between reference and physical cells. For general nonlinear
    // mappings, second order terms would be needed. 

    // Apply once
    AT::matvecProductDataField<double>(tranDer2CubPts,jacInv,der2CubPts);
    FC temp_Der2CubPts(tranDer2CubPts);

    // Apply twice
    AT::matvecProductDataField<double>(tranDer2CubPts,jacInv,temp_Der2CubPts);


    // Scale derivative interpolants by cell length

    for( int cell=0; cell<numCells; ++cell ) {
      double scale = (cellVert(cell,1,0)-cellVert(cell,0,0))/2.0;
      for( int field=0; field<numFields/2; ++field ) {
        for( int pt=0; pt<numCubPts; ++pt ) { 
          tranValsCubPts(cell,2*field+1,pt)   *= scale;
          tranDer2CubPts(cell,2*field+1,pt,0) *= scale;
        }
      }
    }
 
    /********************************************/
    /*   EVALUATE FORCING AND STIFFNESS TERMS   */
    /********************************************/

    FST::multiplyMeasure<double>(wtdTranValsCubPts,wtdMeasure,tranValsCubPts);

    FST::multiplyMeasure<double>(wtdTranDer2CubPts,wtdMeasure,tranDer2CubPts);
    FC temp_wtdTranDer2CubPts(wtdTranDer2CubPts);
    FST::multiplyMeasure<double>(wtdTranDer2CubPts,elasmod,temp_wtdTranDer2CubPts);

    FC loadVectors(numCells,numFields);
    FC stiffnessMatrices(numCells,numFields,numFields);
    
    FST::integrate(loadVectors, qforce, wtdTranValsCubPts, COMP_CPP);
    FST::integrate(stiffnessMatrices, tranDer2CubPts, wtdTranDer2CubPts, COMP_CPP); 

    /***********************************************************/
    /*   ASSEMBLY OF GLOBAL STIFFNESS MATRIX AND LOAD VECTOR   */
    /***********************************************************/

    Vector q(numDof);           // Global Load Vector
    Vector w(numDof);           // Global Displacement Vector (solution)
    Matrix K(numDof,numDof);    // Global Stiffness Matrix
    
    // For the first cell, we exclude the first two fields to enforce the clamped
    // boundary condition at x=0

    for( int row=0; row<numFields-2; ++row ) {
      q(row) = loadVectors(0,row+2);
      for(int col=0; col<numFields-2; ++col ) {
        K(row,col) = stiffnessMatrices(0,row+2,col+2);
      }  
    }

    for( int cell=1; cell<numCells; ++cell ) {
      for( int rf=0; rf<numFields; ++rf ) {
        int row = rf + (numFields-2)*cell-2;
        q(row) += loadVectors(cell,rf);

        for( int cf=0; cf<numFields; ++cf ) {
          int col = cf + (numFields-2)*cell-2;
          K(row,col) += stiffnessMatrices(cell,rf,cf); 
        } 
      }
    } 

    // Boundary conditions
    q(numDof-2) += 6.0;  // Stress boundary condition
    q(numDof-1) -= 6.0;  // Shear force boundary condition

    Solver solver;
    solver.setMatrix(rcpFromRef(K));
    solver.factorWithEquilibration(true);
    solver.factor();
    solver.setVectors(rcpFromRef(w),rcpFromRef(q));
    solver.solve();
  
    int dim = 1+numDof/2;
    Vector w0( dim );
    Vector w1( dim );

    // Separate solution into value and derivative
    for(int i=1; i<dim; ++i) {
       w0(i) = w(2*i-2);   // Extract deflection values
       w1(i) = w(2*i-1);   // Extract deflection derivatives
    }

    // Create exact solution and its derivative
    Vector w0_exact( dim );
    Vector w1_exact( dim );

    int row=0;
    for( int cell=0; cell<numCells; ++cell ) {
      for( int pt=0; pt<numPts-1; ++pt ) {
        double x = physPts(cell,pt,0);
        w0_exact(row) = (3.0-2*x)*x*x;
        w1_exact(row) = 6.0*x*(1.0-x);
        row++;
      }
    }

    w0_exact(dim-1) = 1.0;

    double error0 = 0;
    double error1 = 0;

    for( int i=0; i<dim; ++i ) {
      error0 += std::pow(w0(i)-w0_exact(i),2);
      error1 += std::pow(w1(i)-w1_exact(i),2);
    }
 
    error0 = std::sqrt(error0);
    error1 = std::sqrt(error1);

    *outStream << "\n\n";
    *outStream << "|w-w_exact|   = " << error0 << std::endl;
    *outStream << "|w'-w'_exact| = " << error1 << std::endl;

    double tolerance = 2e-10;

    if( error0 > tolerance ) {
      *outStream << "Solution failed to converge within tolerance" << std::endl;
      errorFlag++;
    }

    if( error1 > tolerance ) {
      *outStream << "Derivative of solution failed to converge within tolerance" << std::endl;
      errorFlag++;
    }

  }

  // Catch unexpected errors
  catch (const std::logic_error & err) {
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
