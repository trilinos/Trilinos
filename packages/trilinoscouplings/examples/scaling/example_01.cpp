// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov),
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   example_01.cpp
    \brief  Example creation of mass and stiffness matrices for div-curl system on a hexadedral mesh using curl-conforming elements.
    \author Created by P. Bochev, D. Ridzal and K. Peterson.
 
    \remark Sample command line
    \code   ./example_01.exe 10 10 10 false 1.0 10.0 0.0 1.0 -1.0 1.0 -1.0 1.0 \endcode
*/


// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_HCURL_HEX_I1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

// Epetra includes
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_LinearProblem.h"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

// AztecOO includes
#include "AztecOO.h"

// ML Includes
#include "ml_MultiLevelPreconditioner.h"
#include "ml_RefMaxwell_11_Operator.h"
#include "ml_RefMaxwell.h"
#include "ml_EdgeMatrixFreePreconditioner.h"
#include "ml_epetra_utils.h"

// Pamgen includes
#include "create_inline_mesh.h"
#include "../mesh_spec_lt/im_exodusII.h"
#include "../mesh_spec_lt/im_ne_nemesisI.h"

#define ABS(x) ((x)>0?(x):-(x))

void TestMultiLevelPreconditioner_CurlLSFEM(char ProblemType[],
                                           Teuchos::ParameterList   & MLList,
                                           Epetra_CrsMatrix   & CurlCurl,
                                           Epetra_CrsMatrix   & D0clean,
                                           Epetra_CrsMatrix   & M0inv,
                                           Epetra_CrsMatrix   & M1,
                                           Epetra_MultiVector & xh,
                                           Epetra_MultiVector & b,
                                           double & TotalErrorResidual,
                                             double & TotalErrorExactSol);

using namespace std;
using namespace Intrepid;


// Functions to evaluate exact solution and derivatives
int evalu(double & uExact0, double & uExact1, double & uExact2, double & x, double & y, double & z);
double evalDivu(double & x, double & y, double & z);
int evalCurlu(double & curlu0, double & curlu1, double & curlu2, double & x, double & y, double & z);
int evalGradDivu(double & gradDivu0, double & gradDivu1, double & gradDivu2, double & x, double & y, double & z);

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
#endif

   //Check number of arguments
    TEST_FOR_EXCEPTION( ( argc < 13 ),
                      std::invalid_argument,
                      ">>> ERROR (example_01): Invalid number of arguments. See code listing for requirements.");
  
  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 12)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|          Example: Div-Curl System on Hexahedral Mesh                        |\n" \
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

  /* In the implementation for discontinuous material properties only the boundaries for
     region 1, associated with mu1, are input. The remainder of the grid is assumed to use mu2.
     Note that the material properties are assigned using the undeformed grid. */

    int NX            = atoi(argv[1]);  // num intervals in x direction (assumed box domain, -1,1)
    int NY            = atoi(argv[2]);  // num intervals in y direction (assumed box domain, -1,1)
    int NZ            = atoi(argv[3]);  // num intervals in z direction (assumed box domain, -1,1)
    int randomMesh    = atoi(argv[4]);  // 1 if mesh randomizer is to be used 0 if not
    double mu1        = atof(argv[5]);  // material property value for region 1
    double mu2        = atof(argv[6]);  // material property value for region 2
    double mu1LeftX   = atof(argv[7]);  // left X boundary for region 1
    double mu1RightX  = atof(argv[8]);  // right X boundary for region 1
    double mu1LeftY   = atof(argv[9]);  // left Y boundary for region 1
    double mu1RightY  = atof(argv[10]); // right Y boundary for region 1
    double mu1LeftZ   = atof(argv[11]); // left Z boundary for region 1
    double mu1RightZ  = atof(argv[12]); // right Z boundary for region 1

// *********************************** CELL TOPOLOGY **********************************

   // Get cell topology for base hexahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology hex_8(shards::getCellTopologyData<shards::Hexahedron<8> >() );

   // Get dimensions 
    int numNodesPerElem = hex_8.getNodeCount();
    int numEdgesPerElem = hex_8.getEdgeCount();
    int numFacesPerElem = hex_8.getSideCount();
    int numNodesPerEdge = 2;
    int numNodesPerFace = 4;
    int numEdgesPerFace = 4;
    int spaceDim = hex_8.getDimension();

   // Build reference element edge to node map
    FieldContainer<int> refEdgeToNode(numEdgesPerElem,numNodesPerEdge);
    for (int i=0; i<numEdgesPerElem; i++){
        refEdgeToNode(i,0)=hex_8.getNodeMap(1, i, 0);
        refEdgeToNode(i,1)=hex_8.getNodeMap(1, i, 1);
    }

// *********************************** GENERATE MESH ************************************

    std::cout << "Generating mesh ... \n\n";

    std::cout << "    NX" << "   NY" << "   NZ\n";
    std::cout << std::setw(5) << NX <<
                 std::setw(5) << NY <<
                 std::setw(5) << NZ << "\n\n";

   // Cube
    double leftX = -1.0, rightX = 1.0;
    double leftY = -1.0, rightY = 1.0;
    double leftZ = -1.0, rightZ = 1.0;

   // Mesh spacing
    double hx = (rightX-leftX)/((double)NX);
    double hy = (rightY-leftY)/((double)NY);
    double hz = (rightZ-leftZ)/((double)NZ);

   // Create Pamgen input file
    stringstream ss;
    ss.clear();
    ss << "mesh \n";
    ss << "  rectilinear \n"; 
    ss << "     nx = " << NX << "\n";
    ss << "     ny = " << NY << "\n"; 
    ss << "     nz = " << NZ << "\n"; 
    ss << "     bx = 1\n";
    ss << "     by = 1\n"; 
    ss << "     bz = 1\n"; 
    ss << "     gmin = " << leftX << " " << leftY << " " << leftZ << "\n";
    ss << "     gmax = " << rightX << " " << rightY << " " << rightZ << "\n";
    ss << "  end \n";
    ss << "  set assign \n";
    ss << "     sideset, ilo, 1\n"; 
    ss << "     sideset, jlo, 2\n"; 
    ss << "     sideset, klo, 3\n"; 
    ss << "     sideset, ihi, 4\n"; 
    ss << "     sideset, jhi, 5\n"; 
    ss << "     sideset, khi, 6\n"; 
    ss << "  end \n";
    ss << "end \n";

    string meshInput = ss.str();
    std::cout << meshInput <<"\n";

   // Generate mesh with Pamgen
    int dim=3;
    int rank=0;
    int numProcs=1;
    long int maxInt = 100000000;
    Create_Pamgen_Mesh(meshInput.c_str(), dim, rank, numProcs, maxInt);
    
   // Get mesh size info
    char title[100];
    int numDim;
    int numNodes;
    int numElems;
    int numElemBlk;
    int numNodeSets;
    int numSideSets;
    int id = 0;

    im_ex_get_init(id, title, &numDim, &numNodes, 
                                &numElems, &numElemBlk, &numNodeSets,
                                &numSideSets);

   // Read node coordinates and place in field container
    FieldContainer<double> nodeCoord(numNodes,dim);
    double * nodeCoordx = new double [numNodes];
    double * nodeCoordy = new double [numNodes];
    double * nodeCoordz = new double [numNodes];
    im_ex_get_coord(id,nodeCoordx,nodeCoordy,nodeCoordz);
    for (int i=0; i<numNodes; i++) {          
      nodeCoord(i,0)=nodeCoordx[i];
      nodeCoord(i,1)=nodeCoordy[i];
      nodeCoord(i,2)=nodeCoordz[i];
    }

   // Get node-element connectivity
    FieldContainer<int> elemToNode(numElems,numNodesPerElem);
    int elemBlkId = 1;
    int * connect = new int [numNodesPerElem*numElems];
    im_ex_get_elem_conn(id,elemBlkId,connect);
    for (int i=0; i<numElems; i++) {
       for (int j=0; j<numNodesPerElem; j++) {
           elemToNode(i,j)=connect[i*numNodesPerElem + j] - 1;
       }
    }
    delete [] connect;


   // Print mesh information  
    int numEdges = (NX)*(NY + 1)*(NZ + 1) + (NX + 1)*(NY)*(NZ + 1) + (NX + 1)*(NY + 1)*(NZ);
    int numFaces = (NX)*(NY)*(NZ + 1) + (NX)*(NY + 1)*(NZ) + (NX + 1)*(NY)*(NZ);
    std::cout << " Number of Elements: " << numElems << " \n";
    std::cout << "    Number of Nodes: " << numNodes << " \n";
    std::cout << "    Number of Edges: " << numEdges << " \n";
    std::cout << "    Number of Faces: " << numFaces << " \n\n";
  
  // Get edge connectivity
    FieldContainer<int> edgeToNode(numEdges, numNodesPerEdge);
    FieldContainer<int> elemToEdge(numElems, numEdgesPerElem);
    int ielem;
    int iedge = 0;
    int inode = 0;
    for (int k=0; k<NZ+1; k++) {
      for (int j=0; j<NY+1; j++) {
        for (int i=0; i<NX+1; i++) {
           if (i < NX){
               edgeToNode(iedge,0) = inode;
               edgeToNode(iedge,1) = inode + 1;
               if (j < NY && k < NZ){
                  ielem=i+j*NX+k*NX*NY;
                  elemToEdge(ielem,0) = iedge;
                  if (j > 0)
                     elemToEdge(ielem-NX,2) = iedge; 
                  if (k > 0)
                     elemToEdge(ielem-NX*NY,4) = iedge; 
                  if (j > 0 && k > 0)
                     elemToEdge(ielem-NX*NY-NX,6) = iedge; 
                }
               else if (j == NY && k == NZ){
                  ielem=i+(NY-1)*NX+(NZ-1)*NX*NY;
                  elemToEdge(ielem,6) = iedge;
                }
               else if (k == NZ && j < NY){
                  ielem=i+j*NX+(NZ-1)*NX*NY;
                  elemToEdge(ielem,4) = iedge;
                  if (j > 0)
                    elemToEdge(ielem-NX,6) = iedge;
                }
               else if (k < NZ && j == NY){
                  ielem=i+(NY-1)*NX+k*NX*NY;
                  elemToEdge(ielem,2) = iedge;
                  if (k > 0)
                     elemToEdge(ielem-NX*NY,6) = iedge;
                }
               iedge++;
            }
           if (j < NY){
               edgeToNode(iedge,0) = inode;
               edgeToNode(iedge,1) = inode + NX+1;
               if (i < NX && k < NZ){
                  ielem=i+j*NX+k*NX*NY;
                  elemToEdge(ielem,3) = iedge;
                  if (i > 0)
                     elemToEdge(ielem-1,1) = iedge; 
                  if (k > 0)
                     elemToEdge(ielem-NX*NY,7) = iedge; 
                  if (i > 0 && k > 0)
                     elemToEdge(ielem-NX*NY-1,5) = iedge; 
                }
               else if (i == NX && k == NZ){
                  ielem=NX-1+j*NX+(NZ-1)*NX*NY;
                  elemToEdge(ielem,5) = iedge;
                }
               else if (k == NZ && i < NX){
                  ielem=i+j*NX+(NZ-1)*NX*NY;
                  elemToEdge(ielem,7) = iedge;
                  if (i > 0)
                    elemToEdge(ielem-1,5) = iedge;
                }
               else if (k < NZ && i == NX){
                  ielem=NX-1+j*NX+k*NX*NY;
                  elemToEdge(ielem,1) = iedge;
                  if (k > 0)
                     elemToEdge(ielem-NX*NY,5) = iedge;
                }
               iedge++;
            }
           if (k < NZ){
               edgeToNode(iedge,0) = inode;
               edgeToNode(iedge,1) = inode + (NX+1)*(NY+1);
               if (i < NX && j < NY){
                  ielem=i+j*NX+k*NX*NY;
                  elemToEdge(ielem,8) = iedge;
                  if (i > 0)
                     elemToEdge(ielem-1,9) = iedge; 
                  if (j > 0)
                     elemToEdge(ielem-NX,11) = iedge; 
                  if (i > 0 && j > 0)
                     elemToEdge(ielem-NX-1,10) = iedge; 
                }
               else if (i == NX && j == NY){
                  ielem=NX-1+(NY-1)*NX+k*NX*NY;
                  elemToEdge(ielem,10) = iedge;
                }
               else if (j == NY && i < NX){
                  ielem=i+(NY-1)*NX+k*NX*NY;
                  elemToEdge(ielem,11) = iedge;
                  if (i > 0)
                    elemToEdge(ielem-1,10) = iedge;
                }
               else if (j < NY && i == NX){
                  ielem=NX-1+j*NX+k*NX*NY;
                  elemToEdge(ielem,9) = iedge;
                  if (j > 0)
                     elemToEdge(ielem-NX,10) = iedge;
                }
               iedge++;
            }
            inode++;
         }
      }
   }
 
   // Output element to edge connectivity
#define DUMP_DATA
#ifdef DUMP_DATA
    ofstream fout("elem2edge.dat");
    ofstream fout2("elem2node.dat");
    for (int k=0; k<NZ; k++) {
      for (int j=0; j<NY; j++) {
        for (int i=0; i<NX; i++) {
          int ielem = i + j * NX + k * NX * NY;
          for (int m=0; m<numNodesPerElem; m++){
              fout2 << elemToNode(ielem,m) <<"  ";
           }
          fout2 <<"\n";
          for (int l=0; l<numEdgesPerElem; l++) {
             fout << elemToEdge(ielem,l) << "  ";
          } 
          fout << "\n";
        }
      }
    }
    fout.close();
    fout2.close();

    ofstream fNodeOut("edge2node.dat");
    for (int i=0; i<numEdges; i++) {
       fNodeOut << edgeToNode(i,0) <<" ";
       fNodeOut << edgeToNode(i,1) <<"\n";
    }
    fNodeOut.close();
#endif

   // Get face connectivity
    FieldContainer<int> faceToNode(numFaces, numNodesPerFace);
    FieldContainer<int> elemToFace(numElems, numFacesPerElem);
    FieldContainer<int> faceToEdge(numFaces, numEdgesPerFace);
    int iface = 0;
    inode = 0;
    for (int k=0; k<NZ+1; k++) {
      for (int j=0; j<NY+1; j++) {
        for (int i=0; i<NX+1; i++) {
           if (i < NX && k < NZ) {
              faceToNode(iface,0)=inode;
              faceToNode(iface,1)=inode + 1;
              faceToNode(iface,2)=inode + (NX+1)*(NY+1)+1;
              faceToNode(iface,3)=inode + (NX+1)*(NY+1);
              if (j < NY) {
                 ielem=i+j*NX+k*NX*NY;
                 faceToEdge(iface,0)=elemToEdge(ielem,0);
                 faceToEdge(iface,1)=elemToEdge(ielem,9);
                 faceToEdge(iface,2)=elemToEdge(ielem,4);
                 faceToEdge(iface,3)=elemToEdge(ielem,8);
                 elemToFace(ielem,0)=iface;
                 if (j > 0) {
                    elemToFace(ielem-NX,2)=iface;
                 }
              }
              else if (j == NY) {
                 ielem=i+(NY-1)*NX+k*NX*NY;
                 faceToEdge(iface,0)=elemToEdge(ielem,2);
                 faceToEdge(iface,1)=elemToEdge(ielem,10);
                 faceToEdge(iface,2)=elemToEdge(ielem,6);
                 faceToEdge(iface,3)=elemToEdge(ielem,11);
                 elemToFace(ielem,2)=iface;
              }
              iface++;
           }
           if (j < NY && k < NZ) {
              faceToNode(iface,0)=inode;
              faceToNode(iface,1)=inode + NX+1;
              faceToNode(iface,2)=inode + (NX+1)*(NY+1) + NX+1;
              faceToNode(iface,3)=inode + (NX+1)*(NY+1);
              if (i < NX) {
                 ielem=i+j*NX+k*NX*NY;
                 faceToEdge(iface,0)=elemToEdge(ielem,3);
                 faceToEdge(iface,1)=elemToEdge(ielem,11);
                 faceToEdge(iface,2)=elemToEdge(ielem,7);
                 faceToEdge(iface,3)=elemToEdge(ielem,8);
                 elemToFace(ielem,3)=iface;
                 if (i > 0) {
                    elemToFace(ielem-1,1)=iface;
                 }
              }
              else if (i == NX) {
                 ielem=NX-1+j*NX+k*NX*NY;
                 faceToEdge(iface,0)=elemToEdge(ielem,1);
                 faceToEdge(iface,1)=elemToEdge(ielem,10);
                 faceToEdge(iface,2)=elemToEdge(ielem,5);
                 faceToEdge(iface,3)=elemToEdge(ielem,9);
                 elemToFace(ielem,1)=iface;
              }
              iface++;
           }
           if (i < NX && j < NY) {
              faceToNode(iface,0)=inode;
              faceToNode(iface,1)=inode + 1;
              faceToNode(iface,2)=inode + NX+2;
              faceToNode(iface,3)=inode + NX+1;
              if (k < NZ) {
                 ielem=i+j*NX+k*NX*NY;
                 faceToEdge(iface,0)=elemToEdge(ielem,0);
                 faceToEdge(iface,1)=elemToEdge(ielem,1);
                 faceToEdge(iface,2)=elemToEdge(ielem,2);
                 faceToEdge(iface,3)=elemToEdge(ielem,3);
                 elemToFace(ielem,4)=iface;
                 if (k > 0) {
                    elemToFace(ielem-NX*NY,5)=iface;
                 }
              }
              else if (k == NZ) {
                 ielem=i+j*NX+(NZ-1)*NX*NY;
                 faceToEdge(iface,0)=elemToEdge(ielem,4);
                 faceToEdge(iface,1)=elemToEdge(ielem,5);
                 faceToEdge(iface,2)=elemToEdge(ielem,6);
                 faceToEdge(iface,3)=elemToEdge(ielem,7);
                 elemToFace(ielem,5)=iface;
              }
              iface++;
           }
          inode++;
         }
      }
   }

  // Get boundary (side set) information
   // Side set 1 - left,  Side set 2 - front, Side set 3 - bottom,
   // Side set 4 - right, Side set 5 - back,  Side set 6 - top
    int * sideSetIds = new int [numSideSets];
    FieldContainer<int> numElemsOnBoundary(numSideSets);
    int numSidesinSet;
    int numDFinSet;
    im_ex_get_side_set_ids(id,sideSetIds);
    for (int i=0; i<numSideSets; i++) {
        im_ex_get_side_set_param(id,sideSetIds[i],&numSidesinSet,&numDFinSet);
        numElemsOnBoundary(i)=numSidesinSet;
     }

   // Container indicating whether an edge or face is on the boundary (1-yes 0-no)
    FieldContainer<int> edgeOnBoundary(numEdges);
    FieldContainer<int> faceOnBoundary(numFaces);

   // Side set 1: left
    if (numElemsOnBoundary(0) > 0){
     int * sideSetElemList1 = new int [numElemsOnBoundary(0)];
     int * sideSetSideList1 = new int [numElemsOnBoundary(0)];
     im_ex_get_side_set(id,sideSetIds[0],sideSetElemList1,sideSetSideList1);
     for (int i=0; i<numElemsOnBoundary(0); i++) {
          edgeOnBoundary(elemToEdge(sideSetElemList1[i]-1,3))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList1[i]-1,8))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList1[i]-1,7))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList1[i]-1,11))=1;
          faceOnBoundary(elemToFace(sideSetElemList1[i]-1,3))=1;
     }
     delete [] sideSetElemList1;
     delete [] sideSetSideList1;
   }

  // Side set 2: front
    if (numElemsOnBoundary(1) > 0){
     int * sideSetElemList2 = new int [numElemsOnBoundary(1)];
     int * sideSetSideList2 = new int [numElemsOnBoundary(1)];
     im_ex_get_side_set(id,sideSetIds[1],sideSetElemList2,sideSetSideList2);
     for (int i=0; i<numElemsOnBoundary(1); i++) {
          edgeOnBoundary(elemToEdge(sideSetElemList2[i]-1,0))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList2[i]-1,9))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList2[i]-1,4))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList2[i]-1,8))=1;
          faceOnBoundary(elemToFace(sideSetElemList2[i]-1,0))=1;
     }
     delete [] sideSetElemList2;
     delete [] sideSetSideList2;
    }

   // Side set 3: bottom
    if (numElemsOnBoundary(2) > 0){
     int * sideSetElemList3 = new int [numElemsOnBoundary(2)];
     int * sideSetSideList3 = new int [numElemsOnBoundary(2)];
     im_ex_get_side_set(id,sideSetIds[2],sideSetElemList3,sideSetSideList3);
     for (int i=0; i<numElemsOnBoundary(2); i++) {
          edgeOnBoundary(elemToEdge(sideSetElemList3[i]-1,0))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList3[i]-1,1))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList3[i]-1,2))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList3[i]-1,3))=1;
          faceOnBoundary(elemToFace(sideSetElemList3[i]-1,4))=1;
     }
     delete [] sideSetElemList3;
     delete [] sideSetSideList3;
    }

   // Side set 4: right
    if (numElemsOnBoundary(3) > 0){
     int * sideSetElemList4 = new int [numElemsOnBoundary(3)];
     int * sideSetSideList4 = new int [numElemsOnBoundary(3)];
     im_ex_get_side_set(id,sideSetIds[3],sideSetElemList4,sideSetSideList4);
     for (int i=0; i<numElemsOnBoundary(3); i++) {
          edgeOnBoundary(elemToEdge(sideSetElemList4[i]-1,1))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList4[i]-1,9))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList4[i]-1,5))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList4[i]-1,10))=1;
          faceOnBoundary(elemToFace(sideSetElemList4[i]-1,1))=1;
     }
     delete [] sideSetElemList4;
     delete [] sideSetSideList4;
    }
  // Side set 5: back
    if (numElemsOnBoundary(4) > 0){
     int * sideSetElemList5 = new int [numElemsOnBoundary(4)];
     int * sideSetSideList5 = new int [numElemsOnBoundary(4)];
     im_ex_get_side_set(id,sideSetIds[4],sideSetElemList5,sideSetSideList5);
     for (int i=0; i<numElemsOnBoundary(4); i++) {
          edgeOnBoundary(elemToEdge(sideSetElemList5[i]-1,2))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList5[i]-1,10))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList5[i]-1,6))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList5[i]-1,11))=1;
          faceOnBoundary(elemToFace(sideSetElemList5[i]-1,2))=1;
     }
     delete [] sideSetElemList5;
     delete [] sideSetSideList5;
    }

   // Side set 6: top
    if (numElemsOnBoundary(5) > 0){
     int * sideSetElemList6 = new int [numElemsOnBoundary(5)];
     int * sideSetSideList6 = new int [numElemsOnBoundary(5)];
     im_ex_get_side_set(id,sideSetIds[5],sideSetElemList6,sideSetSideList6);
     for (int i=0; i<numElemsOnBoundary(5); i++) {
          edgeOnBoundary(elemToEdge(sideSetElemList6[i]-1,4))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList6[i]-1,5))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList6[i]-1,6))=1;
          edgeOnBoundary(elemToEdge(sideSetElemList6[i]-1,7))=1;
          faceOnBoundary(elemToFace(sideSetElemList6[i]-1,5))=1;
     }
     delete [] sideSetElemList6;
     delete [] sideSetSideList6;
    }

    delete [] sideSetIds;

   //TEMP
    ofstream fEdgeout("edgeOnBndy.dat");
    for (int i=0; i<numEdges; i++){
       fEdgeout << edgeOnBoundary(i) <<"\n";
    }
    fEdgeout.close();

   // Container indicating whether a node is on the boundary (1-yes 0-no)
    FieldContainer<int> nodeOnBoundary(numNodes);
     int numEdgeOnBndy=0;
     for (int i=0; i<numEdges; i++) {
        if (edgeOnBoundary(i)){
           nodeOnBoundary(edgeToNode(i,0))=1;
           nodeOnBoundary(edgeToNode(i,1))=1;
           numEdgeOnBndy++;
        }
     }   

    
   // Set material properties using undeformed grid assuming each element has only one value of mu
    FieldContainer<double> muVal(numElems);
    for (int k=0; k<NZ; k++) {
      for (int j=0; j<NY; j++) {
        for (int i=0; i<NX; i++) {
          int ielem = i + j * NX + k * NX * NY;
          double midElemX = nodeCoord(elemToNode(ielem,0),0) + hx/2.0;
          double midElemY = nodeCoord(elemToNode(ielem,0),1) + hy/2.0;
          double midElemZ = nodeCoord(elemToNode(ielem,0),2) + hz/2.0;
          if ( (midElemX > mu1LeftX) && (midElemY > mu1LeftY) && (midElemZ > mu1LeftZ) &&
               (midElemX <= mu1RightX) && (midElemY <= mu1RightY) && (midElemZ <= mu1RightZ) ){
             muVal(ielem) = mu1;
          }
           else {
             muVal(ielem) = mu2;
          }
        }
      }
    }

   // Perturb mesh coordinates (only interior nodes)
    if (randomMesh){
      for (int k=1; k<NZ; k++) {
        for (int j=1; j<NY; j++) {
          for (int i=1; i<NX; i++) {
            int inode = i + j * (NX + 1) + k * (NX + 1) * (NY + 1);
           // random numbers between -1.0 and 1.0
            double rx = 2.0 * (double)rand()/RAND_MAX - 1.0;
            double ry = 2.0 * (double)rand()/RAND_MAX - 1.0;
            double rz = 2.0 * (double)rand()/RAND_MAX - 1.0; 
           // limit variation to 1/4 edge length
            nodeCoord(inode,0) = nodeCoord(inode,0) + 0.125 * hx * rx;
            nodeCoord(inode,1) = nodeCoord(inode,1) + 0.125 * hy * ry;
            nodeCoord(inode,2) = nodeCoord(inode,2) + 0.125 * hz * rz;
          }
        }
      }
    }

#ifdef DUMP_DATA
    // Print coords
    FILE *f=fopen("coords.dat","w");
    inode=0;
    for (int k=0; k<NZ+1; k++) {
      for (int j=0; j<NY+1; j++) {
        for (int i=0; i<NX+1; i++) {          
          fprintf(f,"%22.16e %22.16e %22.16e\n",nodeCoord(inode,0),nodeCoord(inode,1),nodeCoord(inode,2));
          inode++;
        }
      }
    }
    fclose(f);
#endif

// **************************** INCIDENCE MATRIX **************************************

   // Node to edge incidence matrix
    std::cout << "Building incidence matrix ... \n\n";

#ifdef HAVE_MPI
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif

    Epetra_Map globalMapC(numEdges, 0, Comm);
    Epetra_Map globalMapG(numNodes, 0, Comm);
    Epetra_FECrsMatrix DGrad(Copy, globalMapC, globalMapG, 2);

    double vals[2];
    vals[0]=-0.5; vals[1]=0.5;
    for (int j=0; j<numEdges; j++){
        int rowNum = j;
        int colNum[2];
        colNum[0] = edgeToNode(j,0);
        colNum[1] = edgeToNode(j,1);
        DGrad.InsertGlobalValues(1, &rowNum, 2, colNum, vals);
    }


// ************************************ CUBATURE **************************************

   // Get numerical integration points and weights for cell
    std::cout << "Getting cubature ... \n\n";

    DefaultCubatureFactory<double>  cubFactory;                                   
    int cubDegree = 2;
    Teuchos::RCP<Cubature<double> > hexCub = cubFactory.create(hex_8, cubDegree); 

    int cubDim       = hexCub->getDimension();
    int numCubPoints = hexCub->getNumPoints();

    FieldContainer<double> cubPoints(numCubPoints, cubDim);
    FieldContainer<double> cubWeights(numCubPoints);

    hexCub->getCubature(cubPoints, cubWeights);

  
   // Get numerical integration points and weights for hexahedron face
    //             (needed for rhs boundary term)

    // Define topology of the face parametrization domain as [-1,1]x[-1,1]
    CellTopology paramQuadFace(shards::getCellTopologyData<shards::Quadrilateral<4> >() );

    // Define cubature 
    DefaultCubatureFactory<double>  cubFactoryFace;
    Teuchos::RCP<Cubature<double> > hexFaceCubature = cubFactoryFace.create(paramQuadFace, 3);
    int cubFaceDim    = hexFaceCubature -> getDimension();
    int numFacePoints = hexFaceCubature -> getNumPoints();

    // Define storage for cubature points and weights on [-1,1]x[-1,1]
    FieldContainer<double> paramGaussWeights(numFacePoints);
    FieldContainer<double> paramGaussPoints(numFacePoints,cubFaceDim);

    // Define storage for cubature points on workset faces
    hexFaceCubature -> getCubature(paramGaussPoints, paramGaussWeights);


// ************************************** BASIS ***************************************

   // Define basis 
    std::cout << "Getting basis ... \n\n";
    Basis_HCURL_HEX_I1_FEM<double, FieldContainer<double> > hexHCurlBasis;
    Basis_HGRAD_HEX_C1_FEM<double, FieldContainer<double> > hexHGradBasis;

    int numFieldsC = hexHCurlBasis.getCardinality();
    int numFieldsG = hexHGradBasis.getCardinality();

  // Evaluate basis at cubature points
     FieldContainer<double> hexGVals(numFieldsG, numCubPoints); 
     FieldContainer<double> hexCVals(numFieldsC, numCubPoints, spaceDim); 
     FieldContainer<double> hexCurls(numFieldsC, numCubPoints, spaceDim); 
     FieldContainer<double> worksetCVals(numFieldsC, numFacePoints, spaceDim); 

     hexHCurlBasis.getValues(hexCVals, cubPoints, OPERATOR_VALUE);
     hexHCurlBasis.getValues(hexCurls, cubPoints, OPERATOR_CURL);
     hexHGradBasis.getValues(hexGVals, cubPoints, OPERATOR_VALUE);


// ******** LOOP OVER ELEMENTS TO CREATE LOCAL MASS and STIFFNESS MATRICES *************


    std::cout << "Building mass and stiffness matrices ... \n\n";

 // Settings and data structures for mass and stiffness matrices
    typedef CellTools<double>  CellTools;
    typedef FunctionSpaceTools fst;
    typedef ArrayTools art;
    int numCells = 1; 

   // Containers for nodes and edge signs 
    FieldContainer<double> hexNodes(numCells, numNodesPerElem, spaceDim);
    FieldContainer<double> hexEdgeSigns(numCells, numFieldsC);
   // Containers for Jacobian
    FieldContainer<double> hexJacobian(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> hexJacobInv(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> hexJacobDet(numCells, numCubPoints);
   // Containers for element HGRAD mass matrix
    FieldContainer<double> massMatrixG(numCells, numFieldsG, numFieldsG);
    FieldContainer<double> weightedMeasure(numCells, numCubPoints);
    FieldContainer<double> weightedMeasureMuInv(numCells, numCubPoints);
    FieldContainer<double> hexGValsTransformed(numCells, numFieldsG, numCubPoints);
    FieldContainer<double> hexGValsTransformedWeighted(numCells, numFieldsG, numCubPoints);
   // Containers for element HCURL mass matrix
    FieldContainer<double> massMatrixC(numCells, numFieldsC, numFieldsC);
    FieldContainer<double> hexCValsTransformed(numCells, numFieldsC, numCubPoints, spaceDim);
    FieldContainer<double> hexCValsTransformedWeighted(numCells, numFieldsC, numCubPoints, spaceDim);
   // Containers for element HCURL stiffness matrix
    FieldContainer<double> stiffMatrixC(numCells, numFieldsC, numFieldsC);
    FieldContainer<double> weightedMeasureMu(numCells, numCubPoints);    
    FieldContainer<double> hexCurlsTransformed(numCells, numFieldsC, numCubPoints, spaceDim);
    FieldContainer<double> hexCurlsTransformedWeighted(numCells, numFieldsC, numCubPoints, spaceDim);
   // Containers for right hand side vectors
    FieldContainer<double> rhsDatag(numCells, numCubPoints, cubDim);
    FieldContainer<double> rhsDatah(numCells, numCubPoints, cubDim);
    FieldContainer<double> gC(numCells, numFieldsC);
    FieldContainer<double> hC(numCells, numFieldsC);
    FieldContainer<double> hCBoundary(numCells, numFieldsC);
    FieldContainer<double> refGaussPoints(numFacePoints,spaceDim);
    FieldContainer<double> worksetGaussPoints(numCells,numFacePoints,spaceDim);
    FieldContainer<double> worksetJacobians(numCells, numFacePoints, spaceDim, spaceDim);
    FieldContainer<double> worksetJacobInv(numCells, numFacePoints, spaceDim, spaceDim);
    FieldContainer<double> worksetFaceTu(numCells, numFacePoints, spaceDim);
    FieldContainer<double> worksetFaceTv(numCells, numFacePoints, spaceDim);
    FieldContainer<double> worksetFaceN(numCells, numFacePoints, spaceDim);
    FieldContainer<double> worksetVFieldVals(numCells, numFacePoints, spaceDim);
    FieldContainer<double> worksetCValsTransformed(numCells, numFieldsC, numFacePoints, spaceDim);
    FieldContainer<double> divuFace(numCells, numFacePoints);
    FieldContainer<double> worksetFieldDotNormal(numCells, numFieldsC, numFacePoints);

   // Container for cubature points in physical space
    FieldContainer<double> physCubPoints(numCells,numCubPoints, cubDim);

    
   // Global arrays in Epetra format
    Epetra_FECrsMatrix MassG(Copy, globalMapG, numFieldsG);
    Epetra_FECrsMatrix MassC(Copy, globalMapC, numFieldsC);
    Epetra_FECrsMatrix StiffC(Copy, globalMapC, numFieldsC);
    Epetra_FEVector rhsC(globalMapC);

#ifdef DUMP_DATA
    ofstream fSignsout("edgeSigns.dat");
#endif
 // *** Element loop ***
    for (int k=0; k<numElems; k++) {

     // Physical cell coordinates
      for (int i=0; i<numNodesPerElem; i++) {
         hexNodes(0,i,0) = nodeCoord(elemToNode(k,i),0);
         hexNodes(0,i,1) = nodeCoord(elemToNode(k,i),1);
         hexNodes(0,i,2) = nodeCoord(elemToNode(k,i),2);
      }

     // Edge signs
      for (int j=0; j<numEdgesPerElem; j++) {
          if (elemToNode(k,refEdgeToNode(j,0))==edgeToNode(elemToEdge(k,j),0) &&
              elemToNode(k,refEdgeToNode(j,1))==edgeToNode(elemToEdge(k,j),1))
              hexEdgeSigns(0,j) = 1.0;
          else 
              hexEdgeSigns(0,j) = -1.0;
#ifdef DUMP_DATA
          fSignsout << hexEdgeSigns(0,j) << "  ";
#endif
      } 
#ifdef DUMP_DATA     
       fSignsout << "\n";
#endif

       // Compute cell Jacobians, their inverses and their determinants
       CellTools::setJacobian(hexJacobian, cubPoints, hexNodes, hex_8);
       CellTools::setJacobianInv(hexJacobInv, hexJacobian );
       CellTools::setJacobianDet(hexJacobDet, hexJacobian );

// ************************** Compute element HGrad mass matrices *******************************
  
     // transform to physical coordinates 
      fst::HGRADtransformVALUE<double>(hexGValsTransformed, hexGVals);
      
     // compute weighted measure
      fst::computeMeasure<double>(weightedMeasure, hexJacobDet, cubWeights);

      // combine mu value with weighted measure
      for (int nC = 0; nC < numCells; nC++){
        for (int nPt = 0; nPt < numCubPoints; nPt++){
          weightedMeasureMuInv(nC,nPt) = weightedMeasure(nC,nPt) / muVal(k);
        }
      }
      
     // multiply values with weighted measure
      fst::multiplyMeasure<double>(hexGValsTransformedWeighted,
                                   weightedMeasureMuInv, hexGValsTransformed);

     // integrate to compute element mass matrix
      fst::integrate<double>(massMatrixG,
                             hexGValsTransformed, hexGValsTransformedWeighted, COMP_CPP);

      // assemble into global matrix
      int err = 0;
      for (int row = 0; row < numFieldsG; row++){
        for (int col = 0; col < numFieldsG; col++){
            int rowIndex = elemToNode(k,row);
            int colIndex = elemToNode(k,col);
            double val = massMatrixG(0,row,col);
            MassG.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
         }
      }

// ************************** Compute element HCurl mass matrices *******************************

     // transform to physical coordinates 
      fst::HCURLtransformVALUE<double>(hexCValsTransformed, hexJacobInv, 
                                   hexCVals);

     // multiply by weighted measure
      fst::multiplyMeasure<double>(hexCValsTransformedWeighted,
                                   weightedMeasure, hexCValsTransformed);

     // integrate to compute element mass matrix
      fst::integrate<double>(massMatrixC,
                             hexCValsTransformed, hexCValsTransformedWeighted,
                             COMP_CPP);

     // apply edge signs
      fst::applyLeftFieldSigns<double>(massMatrixC, hexEdgeSigns);
      fst::applyRightFieldSigns<double>(massMatrixC, hexEdgeSigns);


     // assemble into global matrix
      err = 0;
      for (int row = 0; row < numFieldsC; row++){
        for (int col = 0; col < numFieldsC; col++){
            int rowIndex = elemToEdge(k,row);
            int colIndex = elemToEdge(k,col);
            double val = massMatrixC(0,row,col);
            MassC.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
         }
      }

// ************************ Compute element HCurl stiffness matrices *****************************

      // transform to physical coordinates 
      fst::HCURLtransformCURL<double>(hexCurlsTransformed, hexJacobian, hexJacobDet, 
                                      hexCurls);

      // combine mu value with weighted measure
      for (int nC = 0; nC < numCells; nC++){
        for (int nPt = 0; nPt < numCubPoints; nPt++){
          weightedMeasureMu(nC,nPt) = weightedMeasure(nC,nPt) / muVal(k);
         }
      }

     // multiply by weighted measure
      fst::multiplyMeasure<double>(hexCurlsTransformedWeighted,
                                   weightedMeasureMu, hexCurlsTransformed);

     // integrate to compute element stiffness matrix
      fst::integrate<double>(stiffMatrixC,
                             hexCurlsTransformed, hexCurlsTransformedWeighted,
                             COMP_CPP);

     // apply edge signs
      fst::applyLeftFieldSigns<double>(stiffMatrixC, hexEdgeSigns);
      fst::applyRightFieldSigns<double>(stiffMatrixC, hexEdgeSigns);

     // assemble into global matrix
      err = 0;
      for (int row = 0; row < numFieldsC; row++){
        for (int col = 0; col < numFieldsC; col++){
            int rowIndex = elemToEdge(k,row);
            int colIndex = elemToEdge(k,col);
            double val = stiffMatrixC(0,row,col);
            StiffC.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
         }
      }

// ******************************* Build right hand side ************************************

      // transform integration points to physical points
       CellTools::mapToPhysicalFrame(physCubPoints, cubPoints, hexNodes, hex_8);

      // evaluate right hand side functions at physical points
       for (int nPt = 0; nPt < numCubPoints; nPt++){

          double x = physCubPoints(0,nPt,0);
          double y = physCubPoints(0,nPt,1);
          double z = physCubPoints(0,nPt,2);
          double du1, du2, du3;

          evalCurlu(du1, du2, du3, x, y, z);
          rhsDatag(0,nPt,0) = du1;
          rhsDatag(0,nPt,1) = du2;
          rhsDatag(0,nPt,2) = du3;
         
          evalGradDivu(du1, du2, du3,  x, y, z);
          rhsDatah(0,nPt,0) = du1;
          rhsDatah(0,nPt,1) = du2;
          rhsDatah(0,nPt,2) = du3;
       }

     // integrate (g,curl w) term
      fst::integrate<double>(gC, rhsDatag, hexCurlsTransformedWeighted,
                             COMP_CPP);

     // integrate (h,div w) term
      fst::integrate<double>(hC, rhsDatah, hexCValsTransformedWeighted,
                             COMP_CPP);

     // apply signs
      fst::applyFieldSigns<double>(gC, hexEdgeSigns);
      fst::applyFieldSigns<double>(hC, hexEdgeSigns);

    // evaluate boundary term
      for (int i = 0; i < numFacesPerElem; i++){
        if (faceOnBoundary(elemToFace(k,i))){
           
         // map Gauss points on quad to reference face: paramGaussPoints -> refGaussPoints
            CellTools::mapToReferenceSubcell(refGaussPoints,
                                   paramGaussPoints,
                                   2, i, hex_8);

         // get basis values at points on reference cell
           hexHCurlBasis.getValues(worksetCVals, refGaussPoints, OPERATOR_VALUE);

         // compute Jacobians at Gauss pts. on reference face for all parent cells
           CellTools::setJacobian(worksetJacobians,
                         refGaussPoints,
                         hexNodes, hex_8);
           CellTools::setJacobianInv(worksetJacobInv, worksetJacobians );

         // transform to physical coordinates 
            fst::HCURLtransformVALUE<double>(worksetCValsTransformed, worksetJacobInv, 
                                   worksetCVals);

         // map Gauss points on quad from ref. face to face workset: refGaussPoints -> worksetGaussPoints
            CellTools::mapToPhysicalFrame(worksetGaussPoints,
                                refGaussPoints,
                                hexNodes, hex_8);

         // compute face tangents
            CellTools::getPhysicalFaceTangents(worksetFaceTu,
                                     worksetFaceTv,
                                     paramGaussPoints,
                                     worksetJacobians,
                                     i, hex_8);

         // face outer normals (relative to parent cell) are uTan x vTan
            RealSpaceTools<double>::vecprod(worksetFaceN, worksetFaceTu, worksetFaceTv);


         // evaluate div u at face points
           for(int nPt = 0; nPt < numFacePoints; nPt++){

             double x = worksetGaussPoints(0, nPt, 0);
             double y = worksetGaussPoints(0, nPt, 1);
             double z = worksetGaussPoints(0, nPt, 2);

             divuFace(0,nPt)=evalDivu(x, y, z);
           }

          // compute the dot product and multiply by Gauss weights
           for (int nF = 0; nF < numFieldsC; nF++){
              for(int nPt = 0; nPt < numFacePoints; nPt++){
                 worksetFieldDotNormal(0,nF,nPt)=0.0;
                  for (int dim = 0; dim < spaceDim; dim++){
                      worksetFieldDotNormal(0,nF,nPt) += worksetCValsTransformed(0,nF,nPt,dim)
                                              * worksetFaceN(0,nPt,dim) * paramGaussWeights(nPt);
                  } //dim
              } //nPt
           } //nF

          // integrate 
          fst::integrate<double>(hCBoundary, divuFace, worksetFieldDotNormal,
                             COMP_CPP);

          // apply signs
           fst::applyFieldSigns<double>(hCBoundary, hexEdgeSigns);

          // add into hC term
            for (int nF = 0; nF < numFieldsC; nF++){
                hC(0,nF) = hC(0,nF) - hCBoundary(0,nF);
            }
          
        } // if faceOnBoundary
      } // numFaces


    // assemble into global vector
     for (int row = 0; row < numFieldsC; row++){
           int rowIndex = elemToEdge(k,row);
           double val = gC(0,row)-hC(0,row);
           rhsC.SumIntoGlobalValues(1, &rowIndex, &val);
     }
 
     
    } // *** end element loop ***

    // Assemble over multiple processors, if necessary
    MassG.GlobalAssemble();  MassG.FillComplete();
    MassC.GlobalAssemble();  MassC.FillComplete();
    StiffC.GlobalAssemble(); StiffC.FillComplete();
    rhsC.GlobalAssemble();

    DGrad.GlobalAssemble(); DGrad.FillComplete(MassG.RowMap(),MassC.RowMap());    

    EpetraExt::RowMatrixToMatlabFile("k1_0.dat",StiffC);
    EpetraExt::MultiVectorToMatrixMarketFile("rhsC0.dat",rhsC,0,0,false);
  
     int numBCEdges=0;
     for (int i=0; i<numEdges; i++){
         if (edgeOnBoundary(i)){
             numBCEdges++;
         }
      }
      int * BCEdges = new int [numBCEdges];
      int indbc=0;
      for (int i=0; i<numEdges; i++){
         if (edgeOnBoundary(i)){
            BCEdges[indbc]=i;
            indbc++;
            rhsC[0][i]=0;
         }
      }
      ML_Epetra::Apply_OAZToMatrix(BCEdges, numBCEdges, StiffC);

      delete [] BCEdges;
  
            
    EpetraExt::RowMatrixToMatlabFile("m1_0.dat",MassC);
    EpetraExt::RowMatrixToMatlabFile("m0.dat",MassG);


    
   // Build the inverse diagonal for MassG
   Epetra_Vector DiagG(MassG.RowMap());
   DiagG.PutScalar(1.0);
   MassG.Multiply(false,DiagG,DiagG);
   for(int i=0;i<DiagG.MyLength();i++) DiagG[i]=1.0/DiagG[i];
   Epetra_CrsMatrix MassGinv(Copy,MassG.RowMap(),MassG.RowMap(),1);
   //   MassGinv.ReplaceDiagonalValues(DiagG);
   for(int i=0;i<DiagG.MyLength();i++) {
     int CID=MassG.GCID(i);
     MassGinv.InsertGlobalValues(MassG.GRID(i),1,&(DiagG[i]),&CID);
   }
   MassGinv.FillComplete();
   
   for(int i=0;i<numNodes;i++) {
     if (nodeOnBoundary(i)){
      double val=0.0;
      MassGinv.ReplaceGlobalValues(i,1,&val,&i);
     }
   }
  

   // Solve!
   Teuchos::ParameterList MLList,dummy;
   double TotalErrorResidual=0, TotalErrorExactSol=0;   
   ML_Epetra::SetDefaultsRefMaxwell(MLList);
   Teuchos::ParameterList MLList2=MLList.get("refmaxwell: 11list",MLList);
   MLList2.set("x-coordinates",nodeCoordx);
   MLList2.set("y-coordinates",nodeCoordy);
   MLList2.set("z-coordinates",nodeCoordz);   
   MLList2.set("ML output",10);
   MLList2.set("smoother: sweeps (level 0)",2);
   MLList2.set("smoother: sweeps",2);
   MLList2.set("smoother: type","Chebyshev");
   //   MLList2.set("eigen-analysis: max iters", 20);
   MLList2.set("eigen-analysis: type", "power-method");
   MLList2.get("edge matrix free: coarse",dummy);
   dummy.set("PDE equations",3);
   dummy.set("ML output",10);
   dummy.set("smoother: sweeps",2);
   dummy.set("smoother: type","symmetric Gauss-Seidel");
   dummy.set("max levels",1);
   dummy.set("coarse: type","Amesos-KLU");
   MLList2.set("edge matrix free: coarse",dummy);


   //   MLList2.set("refmaxwell: disable addon",false);
   //   MLList2.set("print hierarchy",true);
   
   cout<<MLList2<<endl;
   
   //Epetra_FEVector xexact(rhsC);
   //xexact.PutScalar(0.0);//haq

   Epetra_FEVector xh(rhsC);

   MassC.SetLabel("M1");
   StiffC.SetLabel("K1");
   DGrad.SetLabel("D0");
   MassGinv.SetLabel("M0^{-1}");
   

   TestMultiLevelPreconditioner_CurlLSFEM("curl-lsfem",MLList2,StiffC,
                                          DGrad,MassGinv,MassC,
                                          xh,rhsC,
                                          TotalErrorResidual, TotalErrorExactSol);

    // ********  Calculate Error in Solution ***************
     double L2err = 0.0;
     double HCurlerr = 0.0;
     double Linferr = 0.0;

   // Import solution onto current processor
    // Epetra_Map  solnMap(numEdgesGlobal, numEdgesGlobal, 0, Comm);
    // Epetra_Import  solnImporter(solnMap, globalMapG);
    // Epetra_Vector  uCoeff(solnMap);
    // uCoeff.Import(xh, solnImporter, Insert);

   // Get cubature points and weights for error calc (may be different from previous)
     DefaultCubatureFactory<double>  cubFactoryErr;
     int cubDegErr = 3;
     Teuchos::RCP<Cubature<double> > hexCubErr = cubFactoryErr.create(hex_8, cubDegErr);
     int cubDimErr       = hexCubErr->getDimension();
     int numCubPointsErr = hexCubErr->getNumPoints();
     FieldContainer<double> cubPointsErr(numCubPointsErr, cubDimErr);
     FieldContainer<double> cubWeightsErr(numCubPointsErr);
     hexCubErr->getCubature(cubPointsErr, cubWeightsErr);

   // Containers for Jacobian
     FieldContainer<double> hexJacobianE(numCells, numCubPointsErr, spaceDim, spaceDim);
     FieldContainer<double> hexJacobInvE(numCells, numCubPointsErr, spaceDim, spaceDim);
     FieldContainer<double> hexJacobDetE(numCells, numCubPointsErr);
     FieldContainer<double> weightedMeasureE(numCells, numCubPointsErr);

 // Evaluate basis values and curls at cubature points
     FieldContainer<double> uhCVals(numFieldsC, numCubPointsErr, spaceDim);
     FieldContainer<double> uhCValsTrans(numCells,numFieldsC, numCubPointsErr, spaceDim);
     FieldContainer<double> uhCurls(numFieldsC, numCubPointsErr, spaceDim);
     FieldContainer<double> uhCurlsTrans(numCells, numFieldsC, numCubPointsErr, spaceDim);
     hexHCurlBasis.getValues(uhCVals, cubPointsErr, OPERATOR_VALUE);
     hexHCurlBasis.getValues(uhCurls, cubPointsErr, OPERATOR_CURL);


   // Loop over elements
    for (int k=0; k<numElems; k++){

      double L2errElem = 0.0;
      double HCurlerrElem = 0.0;
      double uExact1, uExact2, uExact3;
      double curluExact1, curluExact2, curluExact3;

     // physical cell coordinates
      for (int i=0; i<numNodesPerElem; i++) {
         hexNodes(0,i,0) = nodeCoord(elemToNode(k,i),0);
         hexNodes(0,i,1) = nodeCoord(elemToNode(k,i),1);
         hexNodes(0,i,2) = nodeCoord(elemToNode(k,i),2);
      }
     // Edge signs
      for (int j=0; j<numEdgesPerElem; j++) {
          if (elemToNode(k,refEdgeToNode(j,0))==edgeToNode(elemToEdge(k,j),0) &&
              elemToNode(k,refEdgeToNode(j,1))==edgeToNode(elemToEdge(k,j),1))
              hexEdgeSigns(0,j) = 1.0;
          else 
              hexEdgeSigns(0,j) = -1.0;
      } 

    // compute cell Jacobians, their inverses and their determinants
       CellTools::setJacobian(hexJacobianE, cubPointsErr, hexNodes, hex_8);
       CellTools::setJacobianInv(hexJacobInvE, hexJacobianE );
       CellTools::setJacobianDet(hexJacobDetE, hexJacobianE );

      // transform integration points to physical points
       CellTools::mapToPhysicalFrame(physCubPoints, cubPointsErr, hexNodes, hex_8);

      // transform basis values to physical coordinates
       fst::HCURLtransformVALUE<double>(uhCValsTrans, hexJacobInvE, uhCVals);
       fst::HCURLtransformCURL<double>(uhCurlsTrans, hexJacobianE, hexJacobDetE, uhCurls);

      // compute weighted measure
       fst::computeMeasure<double>(weightedMeasureE, hexJacobDetE, cubWeightsErr);

     // loop over cubature points
       for (int nPt = 0; nPt < numCubPoints; nPt++){

         // get exact solution and curls
          double x = physCubPoints(0,nPt,0);
          double y = physCubPoints(0,nPt,1);
          double z = physCubPoints(0,nPt,2);
          evalu(uExact1, uExact2, uExact3, x, y, z);
          evalCurlu(curluExact1, curluExact2, curluExact3, x, y, z);

         // calculate approximate solution and curls
          double uApprox1 = 0.0;
          double uApprox2 = 0.0;
          double uApprox3 = 0.0;
          double curluApprox1 = 0.0;
          double curluApprox2= 0.0;
          double curluApprox3 = 0.0;
          for (int i = 0; i < numFieldsC; i++){
             int rowIndex = elemToEdge(k,i);
             double uh1 = xh.Values()[rowIndex];
             uApprox1 += uh1*uhCValsTrans(0,i,nPt,0)*hexEdgeSigns(0,i);
             uApprox2 += uh1*uhCValsTrans(0,i,nPt,1)*hexEdgeSigns(0,i);
             uApprox3 += uh1*uhCValsTrans(0,i,nPt,2)*hexEdgeSigns(0,i);
             curluApprox1 += uh1*uhCurlsTrans(0,i,nPt,0)*hexEdgeSigns(0,i);
             curluApprox2 += uh1*uhCurlsTrans(0,i,nPt,1)*hexEdgeSigns(0,i);
             curluApprox3 += uh1*uhCurlsTrans(0,i,nPt,2)*hexEdgeSigns(0,i);
          }

         // evaluate the error at cubature points
          Linferr = max(Linferr, abs(uExact1 - uApprox1));
          Linferr = max(Linferr, abs(uExact2 - uApprox2));
          Linferr = max(Linferr, abs(uExact3 - uApprox3));
          L2errElem+=(uExact1 - uApprox1)*(uExact1 - uApprox1)*weightedMeasureE(0,nPt);
          L2errElem+=(uExact2 - uApprox2)*(uExact2 - uApprox2)*weightedMeasureE(0,nPt);
          L2errElem+=(uExact3 - uApprox3)*(uExact3 - uApprox3)*weightedMeasureE(0,nPt);
          HCurlerrElem+=((curluExact1 - curluApprox1)*(curluExact1 - curluApprox1))
                     *weightedMeasureE(0,nPt);
          HCurlerrElem+=((curluExact2 - curluApprox2)*(curluExact2 - curluApprox2))
                     *weightedMeasureE(0,nPt);
          HCurlerrElem+=((curluExact3 - curluApprox3)*(curluExact3 - curluApprox3))
                     *weightedMeasureE(0,nPt);
        }

       L2err+=L2errElem;
       HCurlerr+=HCurlerrElem;
     }

   // sum over all processors
    //Comm.SumAll(&L2err,&L2errTot,1);
    //Comm.SumAll(&H1err,&H1errTot,1);
    //Comm.MaxAll(&Linferr,&LinferrTot,1);

    std::cout << "\n" << "L2 Error:  " << sqrt(L2err) <<"\n";
    std::cout << "HCurl Error:  " << sqrt(L2err+HCurlerr) <<"\n";
    std::cout << "LInf Error:  " << Linferr <<"\n\n";


    delete [] nodeCoordx;
    delete [] nodeCoordy;
    delete [] nodeCoordz;

#ifdef DUMP_DATA
   fSignsout.close();
#endif
 // delete mesh
 Delete_Pamgen_Mesh();

 // reset format state of std::cout
 std::cout.copyfmt(oldFormatState);
 
 return 0;
}




void solution_test(string msg, const Epetra_Operator &A,const Epetra_MultiVector &lhs,const Epetra_MultiVector &rhs,const Epetra_MultiVector &xexact,Epetra_Time & Time, double & TotalErrorExactSol, double& TotalErrorResidual){
  // ==================================================== //  
  // compute difference between exact solution and ML one //
  // ==================================================== //  
  double d = 0.0, d_tot = 0.0;  
  for( int i=0 ; i<lhs.Map().NumMyElements() ; ++i )
    d += (lhs[0][i] - xexact[0][i]) * (lhs[0][i] - xexact[0][i]);
  
  A.Comm().SumAll(&d,&d_tot,1);
  
  // ================== //
  // compute ||Ax - b|| //
  // ================== //
  double Norm;
  Epetra_Vector Ax(rhs.Map());
  A.Apply(lhs, Ax);
  Ax.Update(1.0, rhs, -1.0);
  Ax.Norm2(&Norm);
  
  if (A.Comm().MyPID() == 0) {
    cout << msg << "......Using " << A.Comm().NumProc() << " processes" << endl;
    cout << msg << "......||A x - b||_2 = " << Norm << endl;
    cout << msg << "......||x_exact - x||_2 = " << sqrt(d_tot) << endl;
    cout << msg << "......Total Time = " << Time.ElapsedTime() << endl;
  }
  
  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm;
}


void TestMultiLevelPreconditioner_CurlLSFEM(char ProblemType[],
                                           Teuchos::ParameterList   & MLList,
                                           Epetra_CrsMatrix   & CurlCurl,
                                           Epetra_CrsMatrix   & D0clean,
                                           Epetra_CrsMatrix   & M0inv,
                                           Epetra_CrsMatrix   & M1,
                                           Epetra_MultiVector & xh,
                                           Epetra_MultiVector & b,
                                           double & TotalErrorResidual,
                                           double & TotalErrorExactSol){
  
  /* Nuke M1 for D0, OAZ*/
  Epetra_CrsMatrix D0(D0clean);                                             
  ML_Epetra::Apply_BCsToGradient(CurlCurl,D0);  
  
  /* Get the BC edges*/
  int numBCedges;  
  int* BCedges=ML_Epetra::FindLocalDiricheltRowsFromOnesAndZeros(CurlCurl,numBCedges);  
  ML_Epetra::Apply_OAZToMatrix(BCedges,numBCedges,M1);

  /* Build the (1,1) Block Operator */
  ML_Epetra::ML_RefMaxwell_11_Operator Operator11(CurlCurl,D0,M0inv,M1);
  
  /* Build the AztecOO stuff */
  Epetra_MultiVector x(xh);
  x.PutScalar(0.0);
  
  Epetra_LinearProblem Problem(&Operator11,&x,&b); 
  Epetra_MultiVector* lhs = Problem.GetLHS();
  Epetra_MultiVector* rhs = Problem.GetRHS();
  
  Epetra_Time Time(CurlCurl.Comm());
    
  /* Build the aggregation guide matrix */
  Epetra_CrsMatrix *TMT_Agg_Matrix;
  ML_Epetra::ML_Epetra_PtAP(M1,D0clean,TMT_Agg_Matrix,false);
  
  /* Approximate the diagonal for EMFP: Double the Diagonal of CurlCurl */
  Epetra_Vector Diagonal(CurlCurl.DomainMap(),false);
  CurlCurl.ExtractDiagonalCopy(Diagonal);
  for(int i=0;i<CurlCurl.NumMyRows();i++) Diagonal[i]*=2;

  /* Build the EMFP Preconditioner */  
  ML_Epetra::EdgeMatrixFreePreconditioner EMFP(Operator11,Diagonal,D0,D0clean,*TMT_Agg_Matrix,BCedges,numBCedges,MLList);

  /* Solve! */
  AztecOO solver(Problem);  
  solver.SetPrecOperator(&EMFP);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 32);
  solver.Iterate(200, 1e-10);
  //  solver.Iterate(1, 1e-10);

  Epetra_MultiVector xexact(xh);
  xexact.PutScalar(0.0);
  
  // accuracy check
  string msg = ProblemType;
  solution_test(msg,Operator11,*lhs,*rhs,xexact,Time,TotalErrorExactSol,TotalErrorResidual);

  #define DUMP_OUT_MATRICES
#ifdef DUMP_OUT_MATRICES
   EpetraExt::RowMatrixToMatlabFile("mag_m0inv_matrix.dat",M0inv);
   EpetraExt::RowMatrixToMatlabFile("mag_m1_matrix.dat",M1);
   EpetraExt::RowMatrixToMatlabFile("mag_k1_matrix.dat",CurlCurl);
   EpetraExt::RowMatrixToMatlabFile("mag_t_matrix.dat",D0);
   EpetraExt::RowMatrixToMatlabFile("mag_t_clean_matrix.dat",D0clean);
   EpetraExt::MultiVectorToMatrixMarketFile("mag_rhs.dat",*rhs,0,0,false);
   EpetraExt::MultiVectorToMatrixMarketFile("mag_lhs.dat",*lhs,0,0,false);
#endif

     xh = *lhs;

}

// Calculates value of exact solution u
 int evalu(double & uExact0, double & uExact1, double & uExact2, double & x, double & y, double & z)
 {
 /*
    uExact0 = exp(x+y+z)*(y+1.0)*(y-1.0)*(z+1.0)*(z-1.0);
    uExact1 = exp(x+y+z)*(x+1.0)*(x-1.0)*(z+1.0)*(z-1.0);
    uExact2 = exp(x+y+z)*(x+1.0)*(x-1.0)*(y+1.0)*(y-1.0);
    
    uExact0 = cos(M_PI*x)*exp(y*z)*(y+1.0)*(y-1.0)*(z+1.0)*(z-1.0);
    uExact1 = cos(M_PI*y)*exp(x*z)*(x+1.0)*(x-1.0)*(z+1.0)*(z-1.0);
    uExact2 = cos(M_PI*z)*exp(x*y)*(x+1.0)*(x-1.0)*(y+1.0)*(y-1.0);
   
    uExact0 = cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
    uExact1 = sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*z);
    uExact2 = sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*z);
 */
 
 
    uExact0 = x*(y*y - 1.0)*(z*z-1.0);
    uExact1 = y*(x*x - 1.0)*(z*z-1.0);
    uExact2 = z*(x*x - 1.0)*(y*y-1.0);

 /*
    uExact0 = (y*y - 1.0)*(z*z-1.0);
    uExact1 = (x*x - 1.0)*(z*z-1.0);
    uExact2 = (x*x - 1.0)*(y*y-1.0);
 
 */
  
   return 0;
 }

// Calculates divergence of exact solution u
 double evalDivu(double & x, double & y, double & z)
 {
  /*
   double divu = exp(x+y+z)*(y+1.0)*(y-1.0)*(z+1.0)*(z-1.0)
                 + exp(x+y+z)*(x+1.0)*(x-1.0)*(z+1.0)*(z-1.0)
                 + exp(x+y+z)*(x+1.0)*(x-1.0)*(y+1.0)*(y-1.0);

   double divu = -M_PI*sin(M_PI*x)*exp(y*z)*(y+1.0)*(y-1.0)*(z+1.0)*(z-1.0)
                 -M_PI*sin(M_PI*y)*exp(x*z)*(x+1.0)*(x-1.0)*(z+1.0)*(z-1.0)
                 -M_PI*sin(M_PI*z)*exp(x*y)*(x+1.0)*(x-1.0)*(y+1.0)*(y-1.0);

   double divu = -3.0*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);

  */
   double divu = (y+1.0)*(y-1.0)*(z+1.0)*(z-1.0)
                 + (x+1.0)*(x-1.0)*(z+1.0)*(z-1.0)
                 + (x+1.0)*(x-1.0)*(y+1.0)*(y-1.0);

  // double divu = 0.0;

   return divu;
 }


// Calculates curl of exact solution u
 int evalCurlu(double & curlu0, double & curlu1, double & curlu2, double & x, double & y, double & z)
 {
  
 /*
   double duxdy = exp(x+y+z)*(z*z-1.0)*(y*y+2.0*y-1.0);
   double duxdz = exp(x+y+z)*(y*y-1.0)*(z*z+2.0*z-1.0);
   double duydx = exp(x+y+z)*(z*z-1.0)*(x*x+2.0*x-1.0);
   double duydz = exp(x+y+z)*(x*x-1.0)*(z*z+2.0*z-1.0);
   double duzdx = exp(x+y+z)*(y*y-1.0)*(x*x+2.0*x-1.0);
   double duzdy = exp(x+y+z)*(x*x-1.0)*(y*y+2.0*y-1.0);

 
   double duxdy = cos(M_PI*x)*exp(y*z)*(z+1.0)*(z-1.0)*(z*(y+1.0)*(y-1.0) + 2.0*y);
   double duxdz = cos(M_PI*x)*exp(y*z)*(y+1.0)*(y-1.0)*(y*(z+1.0)*(z-1.0) + 2.0*z);
   double duydx = cos(M_PI*y)*exp(x*z)*(z+1.0)*(z-1.0)*(z*(x+1.0)*(x-1.0) + 2.0*x);
   double duydz = cos(M_PI*y)*exp(x*z)*(x+1.0)*(x-1.0)*(x*(z+1.0)*(z-1.0) + 2.0*z);
   double duzdx = cos(M_PI*z)*exp(x*y)*(y+1.0)*(y-1.0)*(y*(x+1.0)*(x-1.0) + 2.0*x);
   double duzdy = cos(M_PI*z)*exp(x*y)*(x+1.0)*(x-1.0)*(x*(y+1.0)*(y-1.0) + 2.0*y);

 
   double duxdy = M_PI*cos(M_PI*x)*cos(M_PI*y)*sin(M_PI*z);
   double duxdz = M_PI*cos(M_PI*x)*sin(M_PI*y)*cos(M_PI*z);
   double duydx = M_PI*cos(M_PI*x)*cos(M_PI*y)*sin(M_PI*z);
   double duydz = M_PI*sin(M_PI*x)*cos(M_PI*y)*cos(M_PI*z);
   double duzdx = M_PI*cos(M_PI*x)*sin(M_PI*y)*cos(M_PI*z);
   double duzdy = M_PI*sin(M_PI*x)*cos(M_PI*y)*cos(M_PI*z);
 */
 

   double duxdy = 2.0*x*y*(z*z-1);
   double duxdz = 2.0*x*z*(y*y-1);
   double duydx = 2.0*y*x*(z*z-1);
   double duydz = 2.0*y*z*(x*x-1);
   double duzdx = 2.0*z*x*(y*y-1);
   double duzdy = 2.0*z*y*(x*x-1);
 /*
  
   double duxdy = 2.0*y*(z*z-1);
   double duxdz = 2.0*z*(y*y-1);
   double duydx = 2.0*x*(z*z-1);
   double duydz = 2.0*z*(x*x-1);
   double duzdx = 2.0*x*(y*y-1);
   double duzdy = 2.0*y*(x*x-1);
   
*/

   curlu0 = duzdy - duydz;
   curlu1 = duxdz - duzdx;
   curlu2 = duydx - duxdy;

   return 0;
 }

// Calculates gradient of the divergence of exact solution u
 int evalGradDivu(double & gradDivu0, double & gradDivu1, double & gradDivu2, double & x, double & y, double & z)
{
   
 
  /*
   gradDivu0 = exp(x+y+z)*((y*y-1.0)*(z*z-1.0)+(x*x+2.0*x-1.0)*(z*z-1.0)+(x*x+2.0*x-1.0)*(y*y-1.0));
   gradDivu1 = exp(x+y+z)*((y*y+2.0*y-1.0)*(z*z-1.0)+(x*x-1.0)*(z*z-1.0)+(x*x-1.0)*(y*y+2.0*y-1.0));
   gradDivu2 = exp(x+y+z)*((y*y-1.0)*(z*z+2.0*z-1.0)+(x*x-1.0)*(z*z+2.0*z-1.0)+(x*x-1.0)*(y*y-1.0));
 
   
    gradDivu0 = -M_PI*M_PI*cos(M_PI*x)*exp(y*z)*(y+1.0)*(y-1.0)*(z+1.0)*(z-1.0)
                  -M_PI*sin(M_PI*y)*exp(x*z)*(z+1.0)*(z-1.0)*(z*(x+1.0)*(x-1.0)+2.0*x)
                  -M_PI*sin(M_PI*z)*exp(x*y)*(y+1.0)*(y-1.0)*(y*(x+1.0)*(x-1.0)+2.0*x);
    gradDivu1 = -M_PI*sin(M_PI*x)*exp(y*z)*(z+1.0)*(z-1.0)*(z*(y+1.0)*(y-1.0)+2.0*y)
                  -M_PI*M_PI*cos(M_PI*y)*exp(x*z)*(x+1.0)*(x-1.0)*(z+1.0)*(z-1.0)
                  -M_PI*sin(M_PI*z)*exp(x*y)*(x+1.0)*(x-1.0)*(x*(y+1.0)*(y-1.0)+2.0*y);
    gradDivu2 = -M_PI*sin(M_PI*x)*exp(y*z)*(y+1.0)*(y-1.0)*(y*(z+1.0)*(z-1.0)+2.0*z)
                  -M_PI*sin(M_PI*y)*exp(x*z)*(x+1.0)*(x-1.0)*(x*(z+1.0)*(z-1.0)+2.0*z)
                  -M_PI*M_PI*cos(M_PI*z)*exp(x*y)*(x+1.0)*(x-1.0)*(y+1.0)*(y-1.0);
  

   gradDivu0 = -3.0*M_PI*M_PI*cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
   gradDivu1 = -3.0*M_PI*M_PI*sin(M_PI*x)*cos(M_PI*y)*sin(M_PI*z);
   gradDivu2 = -3.0*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*cos(M_PI*z);

 */

   gradDivu0 = 2.0*x*((z*z-1.0)+(y*y-1.0));
   gradDivu1 = 2.0*y*((z*z-1.0)+(x*x-1.0));
   gradDivu2 = 2.0*z*((x*x-1.0)+(y*y-1.0));

 
 /*
   gradDivu0 = 0;
   gradDivu1 = 0;
   gradDivu2 = 0;
  */
   

   return 0;
}

