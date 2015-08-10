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
    \brief  Example of applying stiffness matrix for a Poisson equation 
            using nodal (Hgrad) elements on squares.  We are using a
	    tensor-product decomposition of the element Poisson operator.


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

     ./Intrepid_example_Drivers_Example_09.exe N verbose
        int deg             - polynomial degree
        int NX              - num intervals in x direction (assumed box domain, 0,1)
        int NY              - num intervals in x direction (assumed box domain, 0,1)
        verbose (optional)  - any character, indicates verbose output

     \endverbatim

    \remark Sample command line
    \code   ./Intrepid_example_Drivers_Example_09.exe 2 10 10 \endcode
*/

// Intrepid includes
#include "Intrepid_CubaturePolylib.hpp"
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

int main(int argc, char *argv[]) {

  //Check number of arguments
   if (argc < 4) {
      std::cout <<"\n>>> ERROR: Invalid number of arguments.\n\n";
      std::cout <<"Usage:\n\n";
      std::cout <<"  ./Intrepid_example_Drivers_Example_09.exe deg NX NY verbose\n\n";
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
  
  int deg           = atoi(argv[1]);  // polynomial degree to use
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
  // I only need this on the line since I'm doing tensor products 
  DefaultCubatureFactory<double>  cubFactory;                                   

  Teuchos::RCP<Cubature<double,FieldContainer<double>,FieldContainer<double> > > glcub
    = Teuchos::rcp(new CubaturePolylib<double,FieldContainer<double>,FieldContainer<double> >(2*deg-1,PL_GAUSS_LOBATTO) );
      
  const int numCubPoints = glcub->getNumPoints();

  FieldContainer<double> cubPoints1D(numCubPoints, 1);
  FieldContainer<double> cubWeights1D(numCubPoints);
  
  glcub->getCubature(cubPoints1D,cubWeights1D);

  
  // ************************************** BASIS ***************************************
  *outStream << "Getting basis ... \n\n";
  
  // Define basis: I only need this on the line also
  Basis_HGRAD_LINE_Cn_FEM<double, FieldContainer<double> > lineHGradBasis(deg,POINTTYPE_SPECTRAL);
  int numLineFieldsG = lineHGradBasis.getCardinality();
  FieldContainer<double> lineGrads(numLineFieldsG, numCubPoints, 1); 
  
  // Evaluate basis values and gradients at cubature points
  lineHGradBasis.getValues(lineGrads, cubPoints1D, OPERATOR_GRAD);

  // ************************************** LTG mapping **********************************
  FieldContainer<int> ltgMapping(numElems,numLineFieldsG*numLineFieldsG);
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
      int ielem = i + j * NX;
      for (int m=0; m<numLineFieldsG; m++){
	ltgout << ltgMapping(ielem,m) <<"  ";
      }
      ltgout <<"\n";
    }
  }
  ltgout.close();
#endif
  

  // Global arrays in Epetra format 
  Epetra_SerialComm Comm;
  Epetra_Map globalMapG(numDOF, 0, Comm);

  Epetra_FEVector u(globalMapG);
  Epetra_FEVector Ku(globalMapG);

  u.Random();

    
  // ************************** Compute element HGrad stiffness matrices *******************************  
//   // Get vertices of all the cells
//   for (int i=0;i<numElems;i++)
//     {
//       for (int j=0;j<4;j++)
// 	{
// 	  const int nodeCur = elemToNode(i,j);
// 	  for (int k=0;k<spaceDim;k++) 
// 	    {
// 	      cellVertices(i,j,k) = nodeCoord(nodeCur,k);
// 	    }
// 	}
//     }

  FieldContainer<double> uScattered(numElems,numLineFieldsG*numLineFieldsG);
  FieldContainer<double> KuScattered(numElems,numLineFieldsG*numLineFieldsG);

  // need storage for derivatives of u on each cell
  // the number of line dof should be the same as the
  // number of cub points.
  // This is indexed by Du(q2,q1)
  FieldContainer<double> Du(numCubPoints,numCubPoints);



  double *uVals = u[0];
  double *KuVals = Ku[0];
  Epetra_Time scatterTime(Comm);
  *outStream << "Scattering\n";
  // Scatter
  for (int k=0; k<numElems; k++) 
    {
      for (int i=0;i<numLineFieldsG*numLineFieldsG;i++) 
        {
          uScattered(k,i) = uVals[ltgMapping(k,i)];
        }
    }
  const double scatTime = scatterTime.ElapsedTime();
  *outStream << "Scattered in time " << scatTime << "\n";
 
  Epetra_Time applyTimer(Comm);
  
  uScattered.resize(numElems,numLineFieldsG,numLineFieldsG);

  for (int k=0;k<numElems;k++)
    {
      // local operation: x-derivative term of stiffness matrix
      //    evaluate x derivative of u on cell k.
      for (int j=0;j<numLineFieldsG;j++)
	{
	  for (int i=0;i<numLineFieldsG;i++)
	    {
	      Du(j,i) = 0.0;
	      for (int q=0;q<numLineFieldsG;q++)
		{
		  Du(j,i) += uScattered(k,j,i) * lineGrads(i,q,0);
		}
	    }
	}

      // initialize Ku
      for (int i=0;i<numLineFieldsG*numLineFieldsG;i++)
	{
	  KuScattered(k,i) = 0.0;
	}

      // loop over basis functions for x term
      int cur = 0;
      for (int j=0;j<numLineFieldsG;j++)
	{
	  for (int i=0;i<numLineFieldsG;i++)
	    {
	      // do the quadrature 
	      for (int q1=0;q1<numCubPoints;q1++)
		{
		  KuScattered(k,cur) += cubWeights1D(j) * cubWeights1D(q1) * Du(j,q1) * lineGrads(i,q1,0);
		}
	      cur ++;
	    }
	}

      // local operation: y-derivative term of stiffness matrix, store in Du
      for (int j=0;j<numLineFieldsG;j++)
	{
	  for (int i=0;i<numLineFieldsG;i++)
	    {
	      Du(j,i) = 0.0;
	      for (int q=0;q<numLineFieldsG;q++)
		{
		  Du(j,i) += uScattered(k,j,i) * lineGrads(j,q,0);
		}
	    }
	}


      // evaluate y-derivatives of u
      cur = 0;
      for (int j=0;j<numLineFieldsG;j++)
	{
	  for (int i=0;i<numLineFieldsG;i++)
	    {
	      for (int q2=0;q2<numCubPoints;q2++)
		{
		  KuScattered(k,cur) += cubWeights1D(q2) * Du(q2,i) * lineGrads(j,q2,0) * cubWeights1D(i);
		}
	    }
	}
    }
  
  uScattered.resize(numElems,numLineFieldsG*numLineFieldsG);
  
  const double applyTime = applyTimer.ElapsedTime();

  *outStream << "Local application: " << applyTime << "\n";
      
  // gather
  Epetra_Time gatherTimer(Comm);
  // Gather
  for (int k=0;k<numElems;k++)
    {
      for (int i=0;i<numLineFieldsG*numLineFieldsG;i++)
        {
          KuVals[ltgMapping(k,i)] += KuScattered(k,i);
        }
    }

  const double gatherTime = gatherTimer.ElapsedTime();
  *outStream << "Gathered in " << gatherTime << "\n";



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

