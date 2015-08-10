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

/** \file   example_03.cpp
    \brief  Example building stiffness matrix and right hand side for a Poisson equation 
            using nodal (Hgrad) elements.

    \verbatim
             div grad u = f in Omega
                      u = 0 on Gamma 

     Discrete linear system for nodal coefficients(x):
        
                 Kx = b

            K - HGrad stiffness matrix
            b - right hand side vector 
                
    \endverbatim

    \author Created by P. Bochev, D. Ridzal and K. Peterson.

    
     \remark Usage
     \verbatim

     ./Intrepid_example_Drivers_Example_03.exe NX NY NZ verbose

        int NX              - num intervals in x direction (assumed box domain, 0,1)
        int NY              - num intervals in y direction (assumed box domain, 0,1)
        int NZ              - num intervals in z direction (assumed box domain, 0,1)
        verbose (optional)  - any character, indicates verbose output

     \endverbatim

    \remark Sample command line
    \code   ./Intrepid_example_Drivers_Example_03.exe 10 10 10 \endcode
*/

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

// Epetra includes
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_SerialComm.h"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

using namespace std;
using namespace Intrepid;

// Functions to evaluate exact solution and derivatives
double evalu(double & x, double & y, double & z);
int evalGradu(double & x, double & y, double & z, double & gradu1, double & gradu2, double & gradu3);
double evalDivGradu(double & x, double & y, double & z);

int main(int argc, char *argv[]) {

  //Check number of arguments
   if (argc < 4) {
      std::cout <<"\n>>> ERROR: Invalid number of arguments.\n\n";
      std::cout <<"Usage:\n\n";
      std::cout <<"  ./Intrepid_example_Drivers_Example_03.exe NX NY NZ verbose\n\n";
      std::cout <<" where \n";
      std::cout <<"   int NX              - num intervals in x direction (assumed box domain, 0,1) \n";
      std::cout <<"   int NY              - num intervals in y direction (assumed box domain, 0,1) \n";
      std::cout <<"   int NZ              - num intervals in z direction (assumed box domain, 0,1) \n";
      std::cout <<"   verbose (optional)  - any character, indicates verbose output \n\n";
      exit(1);
   }
  
  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 3)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|  Example: Generate Stiffness Matrix and Right Hand Side Vector for          |\n" \
    << "|                   Poisson Equation on Hexahedral Mesh                       |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
    << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
    << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";


// ************************************ GET INPUTS **************************************

    int NX            = atoi(argv[1]);  // num intervals in x direction (assumed box domain, 0,1)
    int NY            = atoi(argv[2]);  // num intervals in y direction (assumed box domain, 0,1)
    int NZ            = atoi(argv[3]);  // num intervals in z direction (assumed box domain, 0,1)

// *********************************** CELL TOPOLOGY **********************************

   // Get cell topology for base hexahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology hex_8(shards::getCellTopologyData<shards::Hexahedron<8> >() );

   // Get dimensions 
    int numNodesPerElem = hex_8.getNodeCount();
    int spaceDim = hex_8.getDimension();

// *********************************** GENERATE MESH ************************************

    *outStream << "Generating mesh ... \n\n";

    *outStream << "   NX" << "   NY" << "   NZ\n";
    *outStream << std::setw(5) << NX <<
                 std::setw(5) << NY <<
                 std::setw(5) << NZ << "\n\n";

   // Print mesh information
    int numElems = NX*NY*NZ;
    int numNodes = (NX+1)*(NY+1)*(NZ+1);
    *outStream << " Number of Elements: " << numElems << " \n";
    *outStream << "    Number of Nodes: " << numNodes << " \n\n";

   // Cube
    double leftX = 0.0, rightX = 1.0;
    double leftY = 0.0, rightY = 1.0;
    double leftZ = 0.0, rightZ = 1.0;

   // Mesh spacing
    double hx = (rightX-leftX)/((double)NX);
    double hy = (rightY-leftY)/((double)NY);
    double hz = (rightZ-leftZ)/((double)NZ);

   // Get nodal coordinates
    FieldContainer<double> nodeCoord(numNodes, spaceDim);
    FieldContainer<int> nodeOnBoundary(numNodes);
    int inode = 0;
    for (int k=0; k<NZ+1; k++) {
      for (int j=0; j<NY+1; j++) {
        for (int i=0; i<NX+1; i++) {
          nodeCoord(inode,0) = leftX + (double)i*hx;
          nodeCoord(inode,1) = leftY + (double)j*hy;
          nodeCoord(inode,2) = leftZ + (double)k*hz;
          if (k==0 || j==0 || i==0 || k==NZ || j==NY || i==NX){
             nodeOnBoundary(inode)=1;
          }
          else {
             nodeOnBoundary(inode)=0;
          }
          inode++;
        }
      }
    }
#define DUMP_DATA
#ifdef DUMP_DATA
   // Print nodal coords
    ofstream fcoordout("coords.dat");
    for (int i=0; i<numNodes; i++) {
       fcoordout << nodeCoord(i,0) <<" ";
       fcoordout << nodeCoord(i,1) <<" ";
       fcoordout << nodeCoord(i,2) <<"\n";
    }
    fcoordout.close();
#endif


  // Element to Node map
    FieldContainer<int> elemToNode(numElems, numNodesPerElem);
    int ielem = 0;
    for (int k=0; k<NZ; k++) {
      for (int j=0; j<NY; j++) {
        for (int i=0; i<NX; i++) {
          elemToNode(ielem,0) = (NY + 1)*(NX + 1)*k + (NX + 1)*j + i;
          elemToNode(ielem,1) = (NY + 1)*(NX + 1)*k + (NX + 1)*j + i + 1;
          elemToNode(ielem,2) = (NY + 1)*(NX + 1)*k + (NX + 1)*(j + 1) + i + 1;
          elemToNode(ielem,3) = (NY + 1)*(NX + 1)*k + (NX + 1)*(j + 1) + i;
          elemToNode(ielem,4) = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*j + i;
          elemToNode(ielem,5) = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*j + i + 1;
          elemToNode(ielem,6) = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*(j + 1) + i + 1;
          elemToNode(ielem,7) = (NY + 1)*(NX + 1)*(k + 1) + (NX + 1)*(j + 1) + i;
          ielem++;
        }
      }
    }
#ifdef DUMP_DATA
   // Output connectivity
    ofstream fe2nout("elem2node.dat");
    for (int k=0; k<NZ; k++) {
      for (int j=0; j<NY; j++) {
        for (int i=0; i<NX; i++) {
          int ielem = i + j * NX + k * NX * NY;
          for (int m=0; m<numNodesPerElem; m++){
              fe2nout << elemToNode(ielem,m) <<"  ";
           }
          fe2nout <<"\n";
        }
      }
    }
    fe2nout.close();
#endif


// ************************************ CUBATURE **************************************

    *outStream << "Getting cubature ... \n\n";

   // Get numerical integration points and weights
    DefaultCubatureFactory<double>  cubFactory;                                   
    int cubDegree = 2;
    Teuchos::RCP<Cubature<double> > hexCub = cubFactory.create(hex_8, cubDegree); 

    int cubDim       = hexCub->getDimension();
    int numCubPoints = hexCub->getNumPoints();

    FieldContainer<double> cubPoints(numCubPoints, cubDim);
    FieldContainer<double> cubWeights(numCubPoints);

    hexCub->getCubature(cubPoints, cubWeights);


// ************************************** BASIS ***************************************

     *outStream << "Getting basis ... \n\n";

   // Define basis 
     Basis_HGRAD_HEX_C1_FEM<double, FieldContainer<double> > hexHGradBasis;
     int numFieldsG = hexHGradBasis.getCardinality();
     FieldContainer<double> hexGVals(numFieldsG, numCubPoints); 
     FieldContainer<double> hexGrads(numFieldsG, numCubPoints, spaceDim); 

  // Evaluate basis values and gradients at cubature points
     hexHGradBasis.getValues(hexGVals, cubPoints, OPERATOR_VALUE);
     hexHGradBasis.getValues(hexGrads, cubPoints, OPERATOR_GRAD);


// ******** LOOP OVER ELEMENTS TO CREATE LOCAL STIFFNESS MATRIX *************

    *outStream << "Building stiffness matrix and right hand side ... \n\n";

 // Settings and data structures for mass and stiffness matrices
    typedef CellTools<double>  CellTools;
    typedef FunctionSpaceTools fst;
    int numCells = 1; 

   // Container for nodes
    FieldContainer<double> hexNodes(numCells, numNodesPerElem, spaceDim);
   // Containers for Jacobian
    FieldContainer<double> hexJacobian(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> hexJacobInv(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> hexJacobDet(numCells, numCubPoints);
   // Containers for element HGRAD stiffness matrix
    FieldContainer<double> localStiffMatrix(numCells, numFieldsG, numFieldsG);
    FieldContainer<double> weightedMeasure(numCells, numCubPoints);
    FieldContainer<double> hexGradsTransformed(numCells, numFieldsG, numCubPoints, spaceDim);
    FieldContainer<double> hexGradsTransformedWeighted(numCells, numFieldsG, numCubPoints, spaceDim);
   // Containers for right hand side vectors
    FieldContainer<double> rhsData(numCells, numCubPoints);
    FieldContainer<double> localRHS(numCells, numFieldsG);
    FieldContainer<double> hexGValsTransformed(numCells, numFieldsG, numCubPoints);
    FieldContainer<double> hexGValsTransformedWeighted(numCells, numFieldsG, numCubPoints);
   // Container for cubature points in physical space
    FieldContainer<double> physCubPoints(numCells, numCubPoints, cubDim);

   // Global arrays in Epetra format 
    Epetra_SerialComm Comm;
    Epetra_Map globalMapG(numNodes, 0, Comm);
    Epetra_FECrsMatrix StiffMatrix(Copy, globalMapG, numFieldsG);
    Epetra_FEVector rhs(globalMapG);

 // *** Element loop ***
    for (int k=0; k<numElems; k++) {

     // Physical cell coordinates
      for (int i=0; i<numNodesPerElem; i++) {
         hexNodes(0,i,0) = nodeCoord(elemToNode(k,i),0);
         hexNodes(0,i,1) = nodeCoord(elemToNode(k,i),1);
         hexNodes(0,i,2) = nodeCoord(elemToNode(k,i),2);
      }

    // Compute cell Jacobians, their inverses and their determinants
       CellTools::setJacobian(hexJacobian, cubPoints, hexNodes, hex_8);
       CellTools::setJacobianInv(hexJacobInv, hexJacobian );
       CellTools::setJacobianDet(hexJacobDet, hexJacobian );

// ************************** Compute element HGrad stiffness matrices *******************************
  
     // transform to physical coordinates 
      fst::HGRADtransformGRAD<double>(hexGradsTransformed, hexJacobInv, hexGrads);
      
     // compute weighted measure
      fst::computeCellMeasure<double>(weightedMeasure, hexJacobDet, cubWeights);

     // multiply values with weighted measure
      fst::multiplyMeasure<double>(hexGradsTransformedWeighted,
                                   weightedMeasure, hexGradsTransformed);

     // integrate to compute element stiffness matrix
      fst::integrate<double>(localStiffMatrix,
                             hexGradsTransformed, hexGradsTransformedWeighted, COMP_BLAS);

      // assemble into global matrix
      for (int row = 0; row < numFieldsG; row++){
        for (int col = 0; col < numFieldsG; col++){
            int rowIndex = elemToNode(k,row);
            int colIndex = elemToNode(k,col);
            double val = localStiffMatrix(0,row,col);
            StiffMatrix.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
         }
      }

// ******************************* Build right hand side ************************************

      // transform integration points to physical points
       CellTools::mapToPhysicalFrame(physCubPoints, cubPoints, hexNodes, hex_8);

      // evaluate right hand side function at physical points
       for (int nPt = 0; nPt < numCubPoints; nPt++){

          double x = physCubPoints(0,nPt,0);
          double y = physCubPoints(0,nPt,1);
          double z = physCubPoints(0,nPt,2);

          rhsData(0,nPt) = evalDivGradu(x, y, z);
       }

     // transform basis values to physical coordinates 
      fst::HGRADtransformVALUE<double>(hexGValsTransformed, hexGVals);

     // multiply values with weighted measure
      fst::multiplyMeasure<double>(hexGValsTransformedWeighted,
                                   weightedMeasure, hexGValsTransformed);

     // integrate rhs term
      fst::integrate<double>(localRHS, rhsData, hexGValsTransformedWeighted, 
                             COMP_BLAS);

    // assemble into global vector
     for (int row = 0; row < numFieldsG; row++){
           int rowIndex = elemToNode(k,row);
           double val = -localRHS(0,row);
           rhs.SumIntoGlobalValues(1, &rowIndex, &val);
      }
     
 } // *** end element loop ***


  // Assemble global matrices
   StiffMatrix.GlobalAssemble(); StiffMatrix.FillComplete();
   rhs.GlobalAssemble();

 
  // Adjust stiffness matrix and rhs based on boundary conditions
   for (int row = 0; row<numNodes; row++){
       if (nodeOnBoundary(row)) {
          int rowindex = row;
          for (int col=0; col<numNodes; col++){
              double val = 0.0;
              int colindex = col;
              StiffMatrix.ReplaceGlobalValues(1, &rowindex, 1, &colindex, &val);
          }
          double val = 1.0;
          StiffMatrix.ReplaceGlobalValues(1, &rowindex, 1, &rowindex, &val);
          val = 0.0;
          rhs.ReplaceGlobalValues(1, &rowindex, &val);
       }
    }

#ifdef DUMP_DATA
   // Dump matrices to disk
     EpetraExt::RowMatrixToMatlabFile("stiff_matrix.dat",StiffMatrix);
     EpetraExt::MultiVectorToMatrixMarketFile("rhs_vector.dat",rhs,0,0,false);
#endif

   std::cout << "End Result: TEST PASSED\n";
   
   // reset format state of std::cout
   std::cout.copyfmt(oldFormatState);
   
   return 0;
}


// Calculates value of exact solution u
 double evalu(double & x, double & y, double & z)
 {
 /*
   // function1
    double exactu = sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
 */

   // function2
   double exactu = sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z);

   return exactu;
 }

// Calculates gradient of exact solution u
 int evalGradu(double & x, double & y, double & z, double & gradu1, double & gradu2, double & gradu3)
 {
 /*
   // function 1
       gradu1 = M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
       gradu2 = M_PI*sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*z);
       gradu3 = M_PI*sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*z);
 */

   // function2
       gradu1 = (M_PI*cos(M_PI*x)+sin(M_PI*x))
                  *sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z);
       gradu2 = (M_PI*cos(M_PI*y)+sin(M_PI*y))
                  *sin(M_PI*x)*sin(M_PI*z)*exp(x+y+z);
       gradu3 = (M_PI*cos(M_PI*z)+sin(M_PI*z))
                  *sin(M_PI*x)*sin(M_PI*y)*exp(x+y+z);
  
   return 0;
 }

// Calculates Laplacian of exact solution u
 double evalDivGradu(double & x, double & y, double & z)
 {
 /*
   // function 1
    double divGradu = -3.0*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
 */

   // function 2
   double divGradu = -3.0*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z)
                    + 2.0*M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z)
                    + 2.0*M_PI*cos(M_PI*y)*sin(M_PI*x)*sin(M_PI*z)*exp(x+y+z)
                    + 2.0*M_PI*cos(M_PI*z)*sin(M_PI*x)*sin(M_PI*y)*exp(x+y+z)
                    + 3.0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z);
   
   return divGradu;
 }

