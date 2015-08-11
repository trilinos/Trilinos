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

/** \file   example_16.cpp
    \brief  Application of Laplace operator on a 
    hexahedral mesh using arbitrary-degree elements by
    using TensorProductSpaceTools

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
#include "Intrepid_FunctionSpaceToolsInPlace.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_CubaturePolylib.hpp"
#include "Intrepid_CubatureTensor.hpp"
#include "Intrepid_HGRAD_HEX_Cn_FEM.hpp"
#include "Intrepid_TensorProductSpaceTools.hpp"
#include "Intrepid_Utils.hpp"

// Epetra includes
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_FEVector.h"
#include "Epetra_SerialComm.h"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Time.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_MultiVectorOut.h"

using namespace std;
using namespace Intrepid;
using Teuchos::rcp;

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
//#define DUMP_DATA
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
  
  
  
  // ********************************* CUBATURE AND BASIS *********************************** 
  *outStream << "Getting cubature and basis ... \n\n";
  
  // Get numerical integration points and weights
  // I only need this on the line since I'm doing tensor products 
  Teuchos::RCP<Cubature<double,FieldContainer<double>,FieldContainer<double> > > glcub
    = Teuchos::rcp(new CubaturePolylib<double,FieldContainer<double>,FieldContainer<double> >(2*deg-1,PL_GAUSS_LOBATTO) );
      
  const int numCubPoints1D = glcub->getNumPoints();

  FieldContainer<double> cubPoints1D(numCubPoints1D, 1);
  FieldContainer<double> cubWeights1D(numCubPoints1D);
  
  glcub->getCubature(cubPoints1D,cubWeights1D);

  std::vector<Teuchos::RCP<Cubature<double,FieldContainer<double>,FieldContainer<double> > > >
    cub_to_tensor;  
  cub_to_tensor.push_back( glcub );
  cub_to_tensor.push_back( glcub );
  cub_to_tensor.push_back( glcub );

  Array<RCP<FieldContainer<double> > > wts_by_dim(3);
  wts_by_dim[0] = rcp( &cubWeights1D , false ); wts_by_dim[1] = wts_by_dim[0]; wts_by_dim[2] = wts_by_dim[1];

  CubatureTensor<double,FieldContainer<double>,FieldContainer<double> > cubhex( cub_to_tensor );

  Basis_HGRAD_HEX_Cn_FEM<double, FieldContainer<double> > hexBasis( deg , POINTTYPE_SPECTRAL );

  Array< Array< RCP< Basis< double , FieldContainer<double> > > > > &bases = hexBasis.getBases();

  // get the bases tabulated at the quadrature points, dimension-by-dimension

  Array< RCP< FieldContainer<double> > > basisVals( 3 );
  FieldContainer<double> bvals1D( bases[0][0]->getCardinality() , numCubPoints1D );
  bases[0][0]->getValues( bvals1D , cubPoints1D , OPERATOR_VALUE );
  basisVals[0] = rcp( &bvals1D , false ); basisVals[1] = basisVals[0]; basisVals[2] = basisVals[0];
  
  Array< RCP< FieldContainer<double> > > basisDVals( 3 );
  FieldContainer<double> bdvals1D( bases[0][0]->getCardinality() , numCubPoints1D , 1);
  bases[0][0]->getValues( bdvals1D , cubPoints1D , OPERATOR_D1 );
  basisDVals[0] = rcp( &bdvals1D , false ); basisDVals[1] = basisDVals[0]; basisDVals[2] = basisDVals[0];


  const int numCubPoints = cubhex.getNumPoints();
  FieldContainer<double> cubPoints3D( numCubPoints , 3 );
  FieldContainer<double> cubWeights3D( numCubPoints );
  cubhex.getCubature( cubPoints3D , cubWeights3D );
  

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
              int ielem = i + j * NX + k * NY * NY;
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


  // ********************************* 3-D LOCAL-TO-GLOBAL MAPPING *******************************
  FieldContainer<int> ltgMapping(numElems,hexBasis.getCardinality());
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
  ielem = 0;
  ofstream ltgout("ltg.dat");
  for (int k=0;k<NZ;k++)  
    {
      for (int j=0; j<NY; j++) 
	{
	  for (int i=0; i<NX; i++) 
	    {
	      int ielem = i + j * NX + k * NX * NY;
	      for (int m=0; m<hexBasis.getCardinality(); m++)
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

  // ************* For Jacobians **********************
  FieldContainer<double> cellVertices(numElems,numNodesPerElem,spaceDim);
  FieldContainer<double> cellJacobian(numElems,numCubPoints,spaceDim,spaceDim);
  FieldContainer<double> cellJacobInv(numElems,numCubPoints,spaceDim,spaceDim);
  FieldContainer<double> cellJacobDet(numElems,numCubPoints);


  // get vertices of cells (for computing Jacobians)
  for (int i=0;i<numElems;i++)
    {
      for (int j=0;j<numNodesPerElem;j++)
        {
          const int nodeCur = elemToNode(i,j);
          for (int k=0;k<spaceDim;k++) 
            {
              cellVertices(i,j,k) = nodeCoord(nodeCur,k);
            }
        }
    }

  // jacobian evaluation 
  CellTools<double>::setJacobian(cellJacobian,cubPoints3D,cellVertices,hex_8);
  CellTools<double>::setJacobianInv(cellJacobInv, cellJacobian );
  CellTools<double>::setJacobianDet(cellJacobDet, cellJacobian );


  // ************* MATRIX-FREE APPLICATION 
  FieldContainer<double> uScattered(numElems,1,hexBasis.getCardinality());
  FieldContainer<double> KuScattered(numElems,1,hexBasis.getCardinality());
  FieldContainer<double> gradU(numElems,1,hexBasis.getCardinality(),3);

  u.GlobalAssemble();



  Ku.PutScalar(0.0);
  Ku.GlobalAssemble();

  double *uVals = u[0];
  double *KuVals = Ku[0];

  Teuchos::Time full_timer( "Time to apply operator matrix-free:" );
  Teuchos::Time scatter_timer( "Time to scatter dof:" );
  Teuchos::Time elementwise_timer( "Time to do elementwise computation:" ); 
  Teuchos::Time grad_timer( "Time to compute gradients:" );
  Teuchos::Time pointwise_timer( "Time to do pointwise transformations:" );
  Teuchos::Time moment_timer( "Time to compute moments:" );
  Teuchos::Time gather_timer( "Time to gather dof:" );

  full_timer.start();

  scatter_timer.start();
  for (int k=0; k<numElems; k++) 
    {
      for (int i=0;i<hexBasis.getCardinality();i++) 
        {
          uScattered(k,0,i) = uVals[ltgMapping(k,i)];
        }
    }
  scatter_timer.stop();

  elementwise_timer.start();

  grad_timer.start();
  Intrepid::TensorProductSpaceTools::evaluateGradient<double>( gradU , uScattered ,basisVals , basisDVals );
  grad_timer.stop();
  pointwise_timer.start();
  Intrepid::FunctionSpaceToolsInPlace::HGRADtransformGRAD<double>( gradU , cellJacobian );
  Intrepid::FunctionSpaceToolsInPlace::HGRADtransformGRADDual<double>( gradU , cellJacobian );
  Intrepid::FunctionSpaceToolsInPlace::multiplyMeasure<double>( gradU , cellJacobDet );
  pointwise_timer.stop();
  moment_timer.start();
  Intrepid::TensorProductSpaceTools::momentsGrad<double>( KuScattered , gradU , basisVals , basisDVals , wts_by_dim );
  moment_timer.stop();
  elementwise_timer.stop();
  gather_timer.start();
  for (int k=0;k<numElems;k++)
    {
      for (int i=0;i<hexBasis.getCardinality();i++)
        {
          KuVals[ltgMapping(k,i)] += KuScattered(k,0,i);
        }
    }
  gather_timer.stop();
  full_timer.stop();

  *outStream << full_timer.name() << " " << full_timer.totalElapsedTime() << " sec\n";
  *outStream << "\t" << scatter_timer.name() << " " << scatter_timer.totalElapsedTime() << " sec\n";
  *outStream << "\t" << elementwise_timer.name() << " " << elementwise_timer.totalElapsedTime() << " sec\n";
  *outStream << "\t\t" << grad_timer.name() << " " << grad_timer.totalElapsedTime() << " sec\n";
  *outStream << "\t\t" << pointwise_timer.name() << " " << pointwise_timer.totalElapsedTime() << " sec\n";
  *outStream << "\t\t" << moment_timer.name() << " " << moment_timer.totalElapsedTime() << " sec\n";
  *outStream << "\t" << gather_timer.name() << " " << gather_timer.totalElapsedTime() << " sec\n";


  *outStream << "End Result: TEST PASSED\n";
  
  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  
  return 0;
}

