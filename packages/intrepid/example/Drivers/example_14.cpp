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

/** \file   example_14.cpp
    \brief  Application of Laplace operator on a uniform
    hexahedral mesh using arbitrary-degree elements by
    using tensor product structure and Gauss-Lobatto
    quadrature

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

    ./Intrepid_example_Drivers_Example_14.exe N verbose
    int degree          - polynomial degree
    int NX              - num intervals in x direction (assumed box domain, 0,1)
    int NY              - num intervals in x direction (assumed box domain, 0,1)
    int NZ              - num intervals in x direction (assumed box domain, 0,1)
    verbose (optional)  - any character, indicates verbose output

    \endverbatim

    \remark Sample command line
    \code   ./Intrepid_example_Drivers_Example_14.exe 2 10 10 10 \endcode
*/

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_CubaturePolylib.hpp"
//#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_HGRAD_LINE_Cn_FEM.hpp"
//#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_Utils.hpp"

// Epetra includes
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_FEVector.h"
#include "Epetra_SerialComm.h"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
//#include "Teuchos_BLAS.hpp"
//#include "Teuchos_BLAS_types.hpp"

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
    std::cout <<"  ./Intrepid_example_Drivers_Example_14.exe deg NX NY NZ verbose\n\n";
    std::cout <<" where \n";
    std::cout <<"   int deg             - polynomial degree to be used (assumed >= 1) \n";
    std::cout <<"   int NX              - num intervals in x direction (assumed box domain, 0,1) \n";
    std::cout <<"   int NY              - num intervals in y direction (assumed box domain, 0,1) \n";
    std::cout <<"   int NZ              - num intervals in y direction (assumed box domain, 0,1) \n";
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
  
  *outStream								\
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|  Example: Apply Stiffness Matrix for                                        |\n" \
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
  
  int deg          = atoi(argv[1]);  // polynomial degree to use
  int NX           = atoi(argv[2]);  // num intervals in x direction (assumed box domain, 0,1)
  int NY           = atoi(argv[3]);  // num intervals in y direction (assumed box domain, 0,1)
  int NZ           = atoi(argv[4]);  // num intervals in y direction (assumed box domain, 0,1)
  

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
    std::setw(5) << NY << std::setw(5) << NZ << "\n\n";
  
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
  for (int k=0; k<NZ+1; k++) 
    {
      for (int j=0; j<NY+1; j++) 
	{
	  for (int i=0; i<NX+1; i++) 
	    {
	      nodeCoord(inode,0) = leftX + (double)i*hx;
	      nodeCoord(inode,1) = leftY + (double)j*hy;
	      nodeCoord(inode,2) = leftZ + (double)k*hz;
	      if (k==0 || k==NZ || j==0 || i==0 || j==NY || i==NX)
		{
		  nodeOnBoundary(inode)=1;
		}
	      else 
		{
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
  // We'll keep it around, but this is only the DOFMap if you are in the lowest order case.
  FieldContainer<int> elemToNode(numElems, numNodesPerElem);
  int ielem = 0;
  for (int k=0; k<NZ; k++) 
    {
      for (int j=0; j<NY; j++) 
	{
	  for (int i=0; i<NX; i++) 
	    {
 	      elemToNode(ielem,0) = k * ( NX + 1 ) * ( NY + 1 ) + j * ( NX + 1 ) + i;
 	      elemToNode(ielem,1) = k * ( NX + 1 ) * ( NY + 1 ) + j * ( NX + 1 ) + i + 1;
 	      elemToNode(ielem,2) = k * ( NX + 1 ) * ( NY + 1 ) + ( j + 1 ) * ( NX + 1 ) + i + 1;
 	      elemToNode(ielem,3) = k * ( NX + 1 ) * ( NY + 1 ) + ( j + 1 ) * ( NX + 1 ) + i;
 	      elemToNode(ielem,4) = ( k + 1 ) * ( NX + 1 ) * ( NY + 1 ) + j * ( NX + 1 ) + i;
 	      elemToNode(ielem,5) = ( k + 1 ) * ( NX + 1 ) * ( NY + 1 ) + j * ( NX + 1 ) + i + 1;
 	      elemToNode(ielem,6) = ( k + 1 ) * ( NX + 1 ) * ( NY + 1 ) + ( j + 1 ) * ( NX + 1 ) + i + 1;
 	      elemToNode(ielem,7) = ( k + 1 ) * ( NX + 1 ) * ( NY + 1 ) + ( j + 1 ) * ( NX + 1 ) + i;
	      ielem++;
	    }
	}
    }
#ifdef DUMP_DATA
  // Output connectivity
  ofstream fe2nout("elem2node.dat");
  for (int k=0;k<NZ;k++)
    {
      for (int j=0; j<NY; j++) 
	{
	  for (int i=0; i<NX; i++) 
	    {
	      ielem = i + j * NX + k * NY * NY;
	      for (int m=0; m<numNodesPerElem; m++)
		{
		  fe2nout << elemToNode(ielem,m) <<"  ";
		}
	      fe2nout <<"\n";
	    }
	}
    }
  fe2nout.close();
#endif
  
  // ********************************* 1-D CUBATURE AND BASIS *********************************** 
  *outStream << "Getting cubature and basis ... \n\n";
  
  // Get numerical integration points and weights
  // I only need this on the line since I'm doing tensor products 
  Teuchos::RCP<Cubature<double,FieldContainer<double>,FieldContainer<double> > > glcub
    = Teuchos::rcp(new CubaturePolylib<double,FieldContainer<double>,FieldContainer<double> >(2*deg-1,PL_GAUSS_LOBATTO) );
      
  const int numCubPoints = glcub->getNumPoints();

  FieldContainer<double> cubPoints1D(numCubPoints, 1);
  FieldContainer<double> cubWeights1D(numCubPoints);
  
  glcub->getCubature(cubPoints1D,cubWeights1D);
  // Define basis: I only need this on the line also
  Basis_HGRAD_LINE_Cn_FEM<double, FieldContainer<double> > lineHGradBasis(deg,POINTTYPE_SPECTRAL);
  int numLineFieldsG = lineHGradBasis.getCardinality();
  FieldContainer<double> lineGrads(numLineFieldsG, numCubPoints, 1); 
  
  // Evaluate basis values and gradients at cubature points
  lineHGradBasis.getValues(lineGrads, cubPoints1D, OPERATOR_GRAD);


  // ********************************* 3-D LOCAL-TO-GLOBAL MAPPING *******************************
  FieldContainer<int> ltgMapping(numElems,numLineFieldsG*numLineFieldsG*numLineFieldsG);
  const int numDOF = (NX*deg+1)*(NY*deg+1)*(NZ*deg+1);
  ielem=0;
  for (int k=0;k<NZ;k++) 
    {
      for (int j=0;j<NY;j++) 
	{
	  for (int i=0;i<NX;i++) 
	    {
	      const int start = k * ( NY * deg + 1 ) * ( NX * deg + 1 ) + j * ( NX * deg + 1 ) + i * deg;
	      // loop over local dof on this cell
	      int local_dof_cur=0;
	      for (int kloc=0;kloc<=deg;kloc++) 
		{
		  for (int jloc=0;jloc<=deg;jloc++) 
		    {
		      for (int iloc=0;iloc<=deg;iloc++)
			{
			  ltgMapping(ielem,local_dof_cur) = start 
			    + kloc * ( NX * deg + 1 ) * ( NY * deg + 1 )
			    + jloc * ( NX * deg + 1 )
			    + iloc;
			  local_dof_cur++;
			}
		    }
		}
	      ielem++;
	    }
	}
    }
#ifdef DUMP_DATA
  // Output ltg mapping 
  ofstream ltgout("ltg.dat");
  for (int k=0;k<NZ;k++)  
    {
      for (int j=0; j<NY; j++) 
	{
	  for (int i=0; i<NX; i++) 
	    {
	      ielem = i + j * NX + k * NX * NY;
	      for (int m=0; m<numLineFieldsG*numLineFieldsG*numLineFieldsG; m++)
		{
		  ltgout << ltgMapping(ielem,m) <<"  ";
		}
	      ltgout <<"\n";
	    }
	}
    }
  ltgout.close();
#endif

  // ********** DECLARE GLOBAL OBJECTS *************
  Epetra_SerialComm Comm;
  Epetra_Map globalMapG(numDOF, 0, Comm);
  Epetra_FEVector u(globalMapG);  u.Random();
  Epetra_FEVector Ku(globalMapG);


  // ************* MATRIX-FREE APPLICATION 
  FieldContainer<double> uScattered(numElems,numLineFieldsG*numLineFieldsG*numLineFieldsG);
  FieldContainer<double> KuScattered(numElems,numLineFieldsG*numLineFieldsG*numLineFieldsG);

  u.GlobalAssemble();

  Epetra_Time multTimer(Comm);
  Teuchos::BLAS<int,double> blas;
  Ku.PutScalar(0.0);
  Ku.GlobalAssemble();

  double *uVals = u[0];
  double *KuVals = Ku[0];

  Epetra_Time scatterTimer(Comm);
  std::cout << "Scattering\n";
  // Scatter
  for (int k=0; k<numElems; k++) 
    {
      for (int i=0;i<numLineFieldsG*numLineFieldsG*numLineFieldsG;i++) 
        {
          uScattered(k,i) = uVals[ltgMapping(k,i)];
        }
    }


  const double scatterTime = scatterTimer.ElapsedTime();



  FieldContainer<double> Du(numLineFieldsG,numLineFieldsG,numLineFieldsG);

  Epetra_Time localAppTimer(Comm);

  uScattered.resize(numElems,numLineFieldsG,numLineFieldsG,numLineFieldsG);


  int cur;
  double hcur;

  for (ielem=0;ielem<numElems;ielem++)
    {
      // X-COMPONENT OF ELEMENT LAPLACIAN

      // zero out derivative
      for (int k=0;k<numLineFieldsG;k++)
	{
	  for (int j=0;j<numLineFieldsG;j++)
	    {
	      for (int i=0;i<numLineFieldsG;i++)
		{
		  Du(k,j,i) = 0.0;
		}
	    }
	}


      // compute x derivative
      // this could be a very simple dgemm call
      for (int k=0;k<numLineFieldsG;k++)
	{
	  for (int j=0;j<numLineFieldsG;j++)
	    {
	      for (int i=0;i<numLineFieldsG;i++)
		{
		  for (int q=0;q<numLineFieldsG;q++)
		    {
		      Du(k,j,i) += uScattered(ielem,k,j,i) * lineGrads(i,q,0);
		    }
		}
	    }
	}

      // integration loop for x derivative term
      cur = 0;
      hcur = hy * hz / hx;
      for (int k=0;k<numLineFieldsG;k++)
	{
	  const double wt1 = hcur * cubWeights1D(k);
	  for (int j=0;j<numLineFieldsG;j++)
	    {
	      const double wtstuff = wt1 * cubWeights1D(j);
	      for (int i=0;i<numLineFieldsG;i++)
		{
		  for (int q=0;q<numLineFieldsG;q++)
		    {
		      KuScattered(ielem,cur) += wtstuff
			 * cubWeights1D(q) * Du(k,j,q) * lineGrads(i,q,0);
		    }
		  cur++;
		}
	    }
	}


      // Y-COMPONENT OF ELEMENT LAPLACIAN

      // zero out derivative
      for (int k=0;k<numLineFieldsG;k++)
	{
	  for (int j=0;j<numLineFieldsG;j++)
	    {
	      for (int i=0;i<numLineFieldsG;i++)
		{
		  Du(k,j,i) = 0.0;
		}
	    }
	}

      // compute y derivative
      for (int k=0;k<numLineFieldsG;k++)
	{
	  for (int j=0;j<numLineFieldsG;j++)
	    {
	      for (int i=0;i<numLineFieldsG;i++)
		{
		  for (int q=0;q<numLineFieldsG;q++)
		    {
		      Du(k,j,i) += uScattered(ielem,k,j,i) * lineGrads(j,q,0);
		    }
		}
	    }
	}

      // integration loop for y derivative term
      cur = 0;
      hcur = hx * hz / hy;
      for (int k=0;k<numLineFieldsG;k++)
	{
	  const double wt1 = hcur * cubWeights1D(k);
	  for (int j=0;j<numLineFieldsG;j++)
	    {
	      for (int i=0;i<numLineFieldsG;i++)
		{
		  const double wtstuff = cubWeights1D(i) * wt1;
		  for (int q=0;q<numLineFieldsG;q++)
		    {
		      KuScattered(ielem,cur) += wtstuff * cubWeights1D(q) * Du(k,q,i) * lineGrads(j,q,0);
		    }
		  cur++;
		}
	    }
	}


      // Z-COMPONENT OF ELEMENT LAPLACIAN

      // zero out derivative
      for (int k=0;k<numLineFieldsG;k++)
	{
	  for (int j=0;j<numLineFieldsG;j++)
	    {
	      for (int i=0;i<numLineFieldsG;i++)
		{
		  Du(k,j,i) = 0.0;
		}
	    }
	}

      // compute z derivative
      for (int k=0;k<numLineFieldsG;k++)
	{
	  for (int j=0;j<numLineFieldsG;j++)
	    {
	      for (int i=0;i<numLineFieldsG;i++)
		{
		  for (int q=0;q<numLineFieldsG;q++)
		    {
		      Du(k,j,i) += uScattered(ielem,k,j,i) * lineGrads(k,q,0);
		    }
		}
	    }
	}
      
      // integration loop for z derivative term.
      cur = 0;
      hcur = hx * hy / hz;
      for (int k=0;k<numLineFieldsG;k++)
	{
	  for (int j=0;j<numLineFieldsG;j++)
	    {
	      const double wt1 = hcur * cubWeights1D(j);
	      for (int i=0;i<numLineFieldsG;i++)
		{
		  const double wtstuff = cubWeights1D(i) * wt1;
		  for (int q=0;q<numLineFieldsG;q++)
		    {
		      KuScattered(ielem,cur) += wtstuff * cubWeights1D(q) * Du(q,j,i) * lineGrads(k,q,0);
		    }
		  cur++;
		}
	    }
	}

    }

  const double localAppTime = localAppTimer.ElapsedTime();

  Epetra_Time gatherTimer(Comm);
  // Gather
  for (int k=0;k<numElems;k++)
    {
      for (int i=0;i<numLineFieldsG*numLineFieldsG*numLineFieldsG;i++)
        {
          KuVals[ltgMapping(k,i)] += KuScattered(k,i);
        }
    }
  
  const double gatherTime = gatherTimer.ElapsedTime();
  

  *outStream << "Time to scatter " << scatterTime << "\n";
  *outStream << "Time for local application " << localAppTime << "\n";
  *outStream << "Time to gather " << gatherTime << "\n";
  *outStream << "Total matrix-free time " << scatterTime + localAppTime + gatherTime << "\n";
 

  *outStream << "End Result: TEST PASSED\n";
  
  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  
  return 0;
}

