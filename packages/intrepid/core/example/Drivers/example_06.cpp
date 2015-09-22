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

/** \file   example_06.cpp
    \brief  Matrix-free application of the Laplace stiffness matrix
            for polynomials of degree d on an NX x NY mesh.
	    We are using a reference element stiffness matrix
	    and level 3 BLAS for the application, but not using any
	    tensor-product decomposition.

    \verbatim
             div grad u = f in Omega
                      u = 0 on Gamma 

     Discrete linear system for nodal coefficients(x):
        
                 Kx = b

            K - HGrad stiffness matrix
            b - right hand side vector 
                
    \endverbatim

    \author Created by P. Bochev, R. Kirby, D. Ridzal and K. Peterson.

    
     \remark Usage
     \verbatim

     ./Intrepid_example_Drivers_Example_06.exe N verbose
        int degree          - polynomial degree
        int NX              - num intervals in x direction (assumed box domain, 0,1)
        int NY              - num intervals in x direction (assumed box domain, 0,1)
        verbose (optional)  - any character, indicates verbose output

     \endverbatim

    \remark Sample command line
    \code   ./Intrepid_example_Drivers_Example_06.exe 2 10 10 \endcode
*/

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_HGRAD_QUAD_Cn_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

// Epetra includes
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_FEVector.h"
#include "Epetra_SerialComm.h"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_BLAS_types.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_MultiVectorOut.h"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {

  //Check number of arguments
   if (argc < 4) {
      std::cout <<"\n>>> ERROR: Invalid number of arguments.\n\n";
      std::cout <<"Usage:\n\n";
      std::cout <<"  ./Intrepid_example_Drivers_Example_06.exe deg NX NY verbose\n\n";
      std::cout <<" where \n";
      std::cout <<"   int deg             - polynomial degree to be used (assumed > 1) \n";
      std::cout <<"   int NX              - num intervals in x direction (assumed box domain, 0,1) \n";
      std::cout <<"   int NY              - num intervals in y direction (assumed box domain, 0,1) \n";
      std::cout <<"   verbose (optional)  - any character, indicates verbose output \n\n";
      exit(1);
   }
  
  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 2)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|  Example: Apply Stiffness Matrix for                                        |\n" \
    << "|                   Poisson Equation on Quadrilateral Mesh                    |\n" \
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
  
  int deg          = atoi(argv[1]);  // polynomial degree to use
  int NX            = atoi(argv[2]);  // num intervals in x direction (assumed box domain, 0,1)
  int NY            = atoi(argv[3]);  // num intervals in y direction (assumed box domain, 0,1)
  

  // *********************************** CELL TOPOLOGY **********************************
  
  // Get cell topology for base hexahedron
  typedef shards::CellTopology    CellTopology;
  CellTopology quad_4(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
  
  // Get dimensions 
  int numNodesPerElem = quad_4.getNodeCount();
  int spaceDim = quad_4.getDimension();
  
  // *********************************** GENERATE MESH ************************************
  
  *outStream << "Generating mesh ... \n\n";
  
  *outStream << "   NX" << "   NY\n";
  *outStream << std::setw(5) << NX <<
    std::setw(5) << NY << "\n\n";
  
  // Print mesh information
  int numElems = NX*NY;
  int numNodes = (NX+1)*(NY+1);
  *outStream << " Number of Elements: " << numElems << " \n";
  *outStream << "    Number of Nodes: " << numNodes << " \n\n";
  
  // Square
  double leftX = 0.0, rightX = 1.0;
  double leftY = 0.0, rightY = 1.0;

  // Mesh spacing
  double hx = (rightX-leftX)/((double)NX);
  double hy = (rightY-leftY)/((double)NY);

  // Get nodal coordinates
  FieldContainer<double> nodeCoord(numNodes, spaceDim);
  FieldContainer<int> nodeOnBoundary(numNodes);
  int inode = 0;
  for (int j=0; j<NY+1; j++) {
    for (int i=0; i<NX+1; i++) {
      nodeCoord(inode,0) = leftX + (double)i*hx;
      nodeCoord(inode,1) = leftY + (double)j*hy;
      if (j==0 || i==0 || j==NY || i==NX){
	nodeOnBoundary(inode)=1;
      }
      else {
	nodeOnBoundary(inode)=0;
      }
      inode++;
    }
  }
#define DUMP_DATA
#ifdef DUMP_DATA
  // Print nodal coords
  ofstream fcoordout("coords.dat");
  for (int i=0; i<numNodes; i++) {
    fcoordout << nodeCoord(i,0) <<" ";
    fcoordout << nodeCoord(i,1) <<"\n";
  }
  fcoordout.close();
#endif
  
  
  // Element to Node map
  // We'll keep it around, but this is only the DOFMap if you are in the lowest order case.
  FieldContainer<int> elemToNode(numElems, numNodesPerElem);
  int ielem = 0;
  for (int j=0; j<NY; j++) {
    for (int i=0; i<NX; i++) {
      elemToNode(ielem,0) = (NX + 1)*j + i;
      elemToNode(ielem,1) = (NX + 1)*j + i + 1;
      elemToNode(ielem,2) = (NX + 1)*(j + 1) + i + 1;
      elemToNode(ielem,3) = (NX + 1)*(j + 1) + i;
      ielem++;
    }
  }
#ifdef DUMP_DATA
  // Output connectivity
  ofstream fe2nout("elem2node.dat");
  for (int j=0; j<NY; j++) {
    for (int i=0; i<NX; i++) {
      int ielem = i + j * NX;
      for (int m=0; m<numNodesPerElem; m++){
	fe2nout << elemToNode(ielem,m) <<"  ";
      }
      fe2nout <<"\n";
    }
  }
  fe2nout.close();
#endif
  
  // ************************************ CUBATURE ************************************** 
  *outStream << "Getting cubature ... \n\n";
  
  // Get numerical integration points and weights
  DefaultCubatureFactory<double>  cubFactory;                                   
  int cubDegree = 2*deg;
  Teuchos::RCP<Cubature<double> > quadCub = cubFactory.create(quad_4, cubDegree); 
  
  int cubDim       = quadCub->getDimension();
  int numCubPoints = quadCub->getNumPoints();
  
  FieldContainer<double> cubPoints(numCubPoints, cubDim);
  FieldContainer<double> cubWeights(numCubPoints);
  
  quadCub->getCubature(cubPoints, cubWeights);
  

  // ************************************** BASIS ***************************************
  
  *outStream << "Getting basis ... \n\n";
  
  // Define basis 
  Basis_HGRAD_QUAD_Cn_FEM<double, FieldContainer<double> > quadHGradBasis(deg,POINTTYPE_SPECTRAL);
  int numFieldsG = quadHGradBasis.getCardinality();
  FieldContainer<double> quadGVals(numFieldsG, numCubPoints); 
  FieldContainer<double> quadGrads(numFieldsG, numCubPoints, spaceDim); 
  
  // Evaluate basis values and gradients at cubature points
  quadHGradBasis.getValues(quadGVals, cubPoints, OPERATOR_VALUE);
  quadHGradBasis.getValues(quadGrads, cubPoints, OPERATOR_GRAD);

  // create the local-global mapping for higher order elements
  FieldContainer<int> ltgMapping(numElems,numFieldsG);
  const int numDOF = (NX*deg+1)*(NY*deg+1);
  ielem=0;
  for (int j=0;j<NY;j++) {
    for (int i=0;i<NX;i++) {
      const int start = deg * j * ( NX * deg + 1 ) + i * deg;
      // loop over local dof on this cell
      int local_dof_cur=0;
      for (int vertical=0;vertical<=deg;vertical++) {
	for (int horizontal=0;horizontal<=deg;horizontal++) {
	  ltgMapping(ielem,local_dof_cur) = start + vertical*(NX*deg+1)+horizontal;
	  local_dof_cur++;
	}
      }
      ielem++;
    }
  }
#ifdef DUMP_DATA
  // Output ltg mapping
//   ofstream ltgout("ltg.dat");
//   for (int j=0; j<NY; j++) {
//     for (int i=0; i<NX; i++) {
//       int ielem = i + j * NX;
//       for (int m=0; m<numFieldsG; m++){
// 	ltgout << ltgMapping(ielem,m) <<"  ";
//       }
//       ltgout <<"\n";
//     }
//   }
//   ltgout.close();
#endif
  
  // ******** CREATE A SINGLE STIFFNESS MATRIX, WHICH IS REPLICATED ON ALL ELEMENTS *********
  *outStream << "Applying stiffness matrix and right hand side ... \n\n";

  // Settings and data structures for mass and stiffness matrices
  typedef CellTools<double>  CellTools;
  typedef FunctionSpaceTools fst;
  int numCells = 1; 

  // Container for nodes
  FieldContainer<double> refQuadNodes(numCells, numNodesPerElem, spaceDim);
  // Containers for Jacobian
  FieldContainer<double> refQuadJacobian(numCells, numCubPoints, spaceDim, spaceDim);
  FieldContainer<double> refQuadJacobInv(numCells, numCubPoints, spaceDim, spaceDim);
  FieldContainer<double> refQuadJacobDet(numCells, numCubPoints);
  // Containers for element HGRAD stiffness matrix
  FieldContainer<double> localStiffMatrix(numCells, numFieldsG, numFieldsG);
  FieldContainer<double> weightedMeasure(numCells, numCubPoints);
  FieldContainer<double> quadGradsTransformed(numCells, numFieldsG, numCubPoints, spaceDim);
  FieldContainer<double> quadGradsTransformedWeighted(numCells, numFieldsG, numCubPoints, spaceDim);
  // Containers for right hand side vectors
  FieldContainer<double> rhsData(numCells, numCubPoints);
  FieldContainer<double> localRHS(numCells, numFieldsG);
  FieldContainer<double> quadGValsTransformed(numCells, numFieldsG, numCubPoints);
  FieldContainer<double> quadGValsTransformedWeighted(numCells, numFieldsG, numCubPoints);
  // Container for cubature points in physical space
  FieldContainer<double> physCubPoints(numCells, numCubPoints, cubDim);
  
  // Global arrays in Epetra format 
  Epetra_SerialComm Comm;
  Epetra_Map globalMapG(numDOF, 0, Comm);
  Epetra_FEVector u(globalMapG);
  Epetra_FEVector Ku(globalMapG);
  u.Random();

  std::cout << "About to start ref element matrix\n";

  // ************************** Compute element HGrad stiffness matrices *******************************  
  refQuadNodes(0,0,0) = 0.0;
  refQuadNodes(0,0,1) = 0.0;
  refQuadNodes(0,1,0) = hx;
  refQuadNodes(0,1,1) = 0.0;
  refQuadNodes(0,2,0) = hx;
  refQuadNodes(0,2,1) = hy;
  refQuadNodes(0,3,0) = 0.0;
  refQuadNodes(0,3,1) = hy;

  // Compute cell Jacobians, their inverses and their determinants
  CellTools::setJacobian(refQuadJacobian, cubPoints, refQuadNodes, quad_4);
  CellTools::setJacobianInv(refQuadJacobInv, refQuadJacobian );
  CellTools::setJacobianDet(refQuadJacobDet, refQuadJacobian );
  
  // transform from [-1,1]^2 to [0,hx]x[0,hy]
  fst::HGRADtransformGRAD<double>(quadGradsTransformed, refQuadJacobInv, quadGrads);
      
  // compute weighted measure
  fst::computeCellMeasure<double>(weightedMeasure, refQuadJacobDet, cubWeights);

  // multiply values with weighted measure
  fst::multiplyMeasure<double>(quadGradsTransformedWeighted,
			       weightedMeasure, quadGradsTransformed);

  // integrate to compute element stiffness matrix
  fst::integrate<double>(localStiffMatrix,
			 quadGradsTransformed, quadGradsTransformedWeighted, COMP_BLAS);

  std::cout << "Finished with reference element matrix\n";

  
  // now we will scatter global degrees of freedom, apply the local stiffness matrix 
  // with BLAS, and then gather the results
  FieldContainer<double> uScattered(numElems,numFieldsG);
  FieldContainer<double> KuScattered(numElems,numFieldsG);

  // to extract info from u

  u.GlobalAssemble();

  Epetra_Time multTimer(Comm);

  Ku.PutScalar(0.0);
  Ku.GlobalAssemble();

  double *uVals = u[0];
  double *KuVals = Ku[0];

  Teuchos::BLAS<int,double> blas;
  Epetra_Time scatterTime(Comm);
  std::cout << "Scattering\n";
  // Scatter
  for (int k=0; k<numElems; k++) 
    {
      for (int i=0;i<numFieldsG;i++) 
	{
	  uScattered(k,i) = uVals[ltgMapping(k,i)];
	}
    }
  const double scatTime = scatterTime.ElapsedTime();
  std::cout << "Scattered in time " << scatTime << "\n";

  Epetra_Time blasTimer(Comm);
  blas.GEMM(Teuchos::NO_TRANS , Teuchos::NO_TRANS , 
	    numFieldsG , numElems, numFieldsG  , 
	    1.0 , 
	    &localStiffMatrix(0,0,0) , 
	    numFieldsG ,
	    &uScattered(0,0) , 
	    numFieldsG , 
	    0.0 , 
	     &KuScattered(0,0) , 
	    numFieldsG );
  const double blasTime = blasTimer.ElapsedTime();
  std::cout << "Element matrices applied in " << blasTime << "\n";

  Epetra_Time gatherTimer(Comm);
  // Gather
  for (int k=0;k<numElems;k++)
    {
      for (int i=0;i<numFieldsG;i++)
	{
	  KuVals[ltgMapping(k,i)] += KuScattered(k,i);
	}
    }

  const double gatherTime = gatherTimer.ElapsedTime();
  std::cout << "Gathered in " << gatherTime << "\n";
  

  const double applyTime = gatherTime + blasTime + scatTime;
  std::cout << "Time to do matrix-free product: " << applyTime << std::endl;
  
  
  std::cout << "End Result: TEST PASSED\n";
  
  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  
  return 0;
}

