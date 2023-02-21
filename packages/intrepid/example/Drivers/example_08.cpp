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

/** \file   example_08.cpp
    \brief  Example building stiffness matrix and right hand side for a Poisson equation 
            using nodal (Hgrad) elements on squares.
	    This code transforms the basis function gradients to each cell and performs
	    quadrature. 

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

     ./Intrepid_example_Drivers_Example_05.exe N verbose
        int deg             - polynomial degree
        int NX              - num intervals in x direction (assumed box domain, 0,1)
        int NY              - num intervals in x direction (assumed box domain, 0,1)
        verbose (optional)  - any character, indicates verbose output

     \endverbatim

    \remark Sample command line
    \code   ./Intrepid_example_Drivers_Example_05.exe 2 10 10 \endcode
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
      std::cout <<"  ./Intrepid_example_Drivers_Example_05.exe deg NX NY verbose\n\n";
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
    << "|  Example: Generate Stiffness Matrix and Right Hand Side Vector for          |\n" \
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
  int ielem=0;
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
      ielem = i + j * NX;
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
  ofstream ltgout("ltg.dat");
  for (int j=0; j<NY; j++) {
    for (int i=0; i<NX; i++) {
      ielem = i + j * NX;
      for (int m=0; m<numFieldsG; m++){
	ltgout << ltgMapping(ielem,m) <<"  ";
      }
      ltgout <<"\n";
    }
  }
  ltgout.close();
#endif
  
  // ******** CREATE ALL LOCAL STIFFNESS MATRICES *********
  *outStream << "Building stiffness matrix and right hand side ... \n\n";

  // Settings and data structures for mass and stiffness matrices
  typedef CellTools<double>  CellTools;
  typedef FunctionSpaceTools fst;
  int numCells = numElems; 

  // Container for nodes
  FieldContainer<double> cellVertices(numCells,numNodesPerElem,spaceDim);

  // Containers for Jacobian
  FieldContainer<double> cellJacobian(numCells, numCubPoints, spaceDim, spaceDim);
  FieldContainer<double> cellJacobInv(numCells, numCubPoints, spaceDim, spaceDim);
  FieldContainer<double> cellJacobDet(numCells, numCubPoints);

  // Containers for element HGRAD stiffness matrices
  FieldContainer<double> localStiffMatrices(numCells, numFieldsG, numFieldsG);
  FieldContainer<double> transformedBasisGradients(numCells,numFieldsG,numCubPoints,spaceDim);
  FieldContainer<double> weightedTransformedBasisGradients(numCells,numFieldsG,numCubPoints,spaceDim);
  FieldContainer<double> weightedMeasure(numCells, numCubPoints);

  
  // Global arrays in Epetra format 
  Epetra_SerialComm Comm;
  Epetra_Map globalMapG(numDOF, 0, Comm);

  Epetra_Time graphTimer(Comm);
  Epetra_CrsGraph grph( Copy , globalMapG , 4 * numFieldsG );
  for (int k=0;k<numElems;k++) 
    {
      for (int i=0;i<numFieldsG;i++)
        {
          grph.InsertGlobalIndices(ltgMapping(k,i),numFieldsG,&ltgMapping(k,0));
        }
    }
  grph.FillComplete();

  const double graphTime = graphTimer.ElapsedTime();
  std::cout << "Graph computed in " << graphTime << "\n";

  Epetra_Time instantiateTimer( Comm );
  Epetra_FECrsMatrix StiffMatrix( Copy , grph );
  const double instantiateTime = instantiateTimer.ElapsedTime(  );
  std::cout << "Matrix instantiated in " << instantiateTime << "\n";

  Epetra_FEVector u(globalMapG);
  Epetra_FEVector Ku(globalMapG);

  u.Random();
    
  // ************************** Compute element HGrad stiffness matrices *******************************  
  // Get vertices of all the cells

  for (int i=0;i<numElems;i++)
    {
      for (int j=0;j<4;j++)
	{
	  const int nodeCur = elemToNode(i,j);
	  for (int k=0;k<spaceDim;k++) 
	    {
	      cellVertices(i,j,k) = nodeCoord(nodeCur,k);
	    }
	}
    }

  Epetra_Time localConstructTimer(Comm);

  // Compute all cell Jacobians, their inverses and their determinants
  CellTools::setJacobian(cellJacobian, cubPoints, cellVertices, quad_4);
  CellTools::setJacobianInv(cellJacobInv, cellJacobian );
  CellTools::setJacobianDet(cellJacobDet, cellJacobian );
  
  // transform reference element gradients to each cell
  fst::HGRADtransformGRAD<double>(transformedBasisGradients, cellJacobInv, quadGrads);
      
  // compute weighted measure
  fst::computeCellMeasure<double>(weightedMeasure, cellJacobDet, cubWeights);

  // multiply values with weighted measure
  fst::multiplyMeasure<double>(weightedTransformedBasisGradients,
			       weightedMeasure, transformedBasisGradients);

  // integrate to compute element stiffness matrix
  fst::integrate<double>(localStiffMatrices,
			 transformedBasisGradients, weightedTransformedBasisGradients , COMP_BLAS);

  const double localConstructTime = localConstructTimer.ElapsedTime();
  std::cout << "Time to build local matrices (including Jacobian computation): "<< localConstructTime << "\n";

  Epetra_Time insertionTimer(Comm);

  // *** Element loop ***
   for (int k=0; k<numElems; k++) 
     {
       // assemble into global matrix
       StiffMatrix.InsertGlobalValues(numFieldsG,&ltgMapping(k,0),numFieldsG,&ltgMapping(k,0),&localStiffMatrices(k,0,0));

     }
   StiffMatrix.GlobalAssemble(); StiffMatrix.FillComplete();
   const double insertionTime = insertionTimer.ElapsedTime( );
   
   std::cout << "Time to assemble global matrix from local matrices: " << insertionTime << "\n";


#ifdef DUMP_DATA
   // Dump matrices to disk
//    EpetraExt::RowMatrixToMatlabFile("stiff_matrix.dat",StiffMatrix);
//    EpetraExt::MultiVectorToMatrixMarketFile("rhs_vector.dat",rhs,0,0,false);
#endif

   
   std::cout << "End Result: TEST PASSED\n";

   // reset format state of std::cout
   std::cout.copyfmt(oldFormatState);
   
   return 0;
}

