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
    \brief  Example creation of mass and stiffness matrices for div-curl system on a hexadedral mesh using div-conforming elements.
    \author Created by P. Bochev, D. Ridzal and K. Peterson.
 
    \remark Sample command line
    \code   ./example_DivLSFEM.exe 10 10 10 \endcode
*/

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_HCURL_HEX_I1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HDIV_HEX_I1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

// Epetra includes
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_SerialComm.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

// Pamgen includes
#include "create_inline_mesh.h"
#include "../mesh_spec_lt/im_exodusII.h"
#include "../mesh_spec_lt/im_ne_nemesisI.h"

using namespace std;
using namespace Intrepid;
using namespace shards;

// Functions to evaluate exact solution and its derivatives
int evalu(double & uExact0, double & uExact1, double & uExact2, double & x, double & y, double & z);
double evalDivu(double & x, double & y, double & z);
int evalCurlu(double & curlu0, double & curlu1, double & curlu2, double & x, double & y, double & z);
int evalCurlCurlu(double & curlCurlu0, double & curlCurlu1, double & curlCurlu2, double & x, double & y, double & z);

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

   //Check number of arguments
    TEST_FOR_EXCEPTION( ( argc < 4 ),
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

// *********************************** CELL TOPOLOGY **********************************

   // Get cell topology for base hexahedron
    typedef shards::CellTopology    CellTopology;
    CellTopology hex_8(shards::getCellTopologyData<Hexahedron<8> >() );

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

   // Build reference element face to node map
    FieldContainer<int> refFaceToNode(numFacesPerElem,numNodesPerFace);
    for (int i=0; i<numFacesPerElem; i++){
        refFaceToNode(i,0)=hex_8.getNodeMap(2, i, 0);
        refFaceToNode(i,1)=hex_8.getNodeMap(2, i, 1);
        refFaceToNode(i,2)=hex_8.getNodeMap(2, i, 2);
        refFaceToNode(i,3)=hex_8.getNodeMap(2, i, 3);
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
  
    char inputArray[1000];
    int ntmp = meshInput.size();
    for (int i=0; i<ntmp; i++) {          
      inputArray[i]=meshInput[i];
    }

   // Generate mesh with Pamgen
    int dim=3;
    int rank=0;
    int numProcs=1;
    long int maxInt = 100000000;
    Create_Pamgen_Mesh(inputArray, dim, rank, numProcs, maxInt);
    
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
    delete [] nodeCoordx;
    delete [] nodeCoordy;
    delete [] nodeCoordz;

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
                 faceToEdge(iface,1)=elemToEdge(ielem,8);
                 faceToEdge(iface,2)=elemToEdge(ielem,7);
                 faceToEdge(iface,3)=elemToEdge(ielem,11);
                 elemToFace(ielem,3)=iface;
                 if (i > 0) {
                    elemToFace(ielem-1,1)=iface;
                 }
              }
              else if (i == NX) {
                 ielem=NX-1+j*NX+k*NX*NY;
                 faceToEdge(iface,0)=elemToEdge(ielem,1);
                 faceToEdge(iface,1)=elemToEdge(ielem,9);
                 faceToEdge(iface,2)=elemToEdge(ielem,5);
                 faceToEdge(iface,3)=elemToEdge(ielem,10);
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
 
   // Output element to face connectivity
    ofstream el2fout("elem2face.dat");
    for (int k=0; k<NZ; k++) {
      for (int j=0; j<NY; j++) {
        for (int i=0; i<NX; i++) {
          int ielem = i + j * NX + k * NX * NY;
          for (int l=0; l<numFacesPerElem; l++) {
             el2fout << elemToFace(ielem,l) << "  ";
          } 
          el2fout << "\n";
        }
      }
    }
    el2fout.close();

   // Output face to edge and face to node connectivity
    ofstream f2edout("face2edge.dat");
    ofstream f2nout("face2node.dat");
    for (int k=0; k<numFaces; k++) {
       for (int i=0; i<numEdgesPerFace; i++) {
           f2edout << faceToEdge(k,i) << "  ";
       } 
       for (int j=0; j<numNodesPerFace; j++) {
           f2nout << faceToNode(k,j) << "  ";
       } 
       f2edout << "\n";
       f2nout << "\n";
    }
    f2edout.close();
    f2nout.close();

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

   // Container indicating whether a face is on the boundary (1-yes 0-no)
    FieldContainer<int> faceOnBoundary(numEdges);

   // Side set 1: left
    if (numElemsOnBoundary(0) > 0){
     int * sideSetElemList1 = new int [numElemsOnBoundary(0)];
     int * sideSetSideList1 = new int [numElemsOnBoundary(0)];
     im_ex_get_side_set(id,sideSetIds[0],sideSetElemList1,sideSetSideList1);
     for (int i=0; i<numElemsOnBoundary(0); i++) {
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
          faceOnBoundary(elemToFace(sideSetElemList6[i]-1,5))=1;
     }
     delete [] sideSetElemList6;
     delete [] sideSetSideList6;
    }

    delete [] sideSetIds;

   //TEMP
    ofstream fFaceout("faceOnBndy.dat");
    for (int i=0; i<numFaces; i++){
       fFaceout << faceOnBoundary(i) <<"  ";
    }
    fFaceout.close();




/*
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
*/

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


// **************************** INCIDENCE MATRIX **************************************

   // Edge to face incidence matrix
    std::cout << "Building incidence matrix ... \n\n";

    Epetra_SerialComm Comm;
    Epetra_Map globalMapD(numFaces, 0, Comm);
    Epetra_Map globalMapC(numEdges, 0, Comm);
    Epetra_FECrsMatrix DCurl(Copy, globalMapD, globalMapC, 2);

    double vals[4];
    vals[0]=1; vals[1]=1; vals[2]=-1; vals[3]=-1;
    for (int j=0; j<numFaces; j++){
        int rowNum = j;
        int colNum[4];
        colNum[0] = faceToEdge(j,0);
        colNum[1] = faceToEdge(j,1);
        colNum[2] = faceToEdge(j,2);
        colNum[3] = faceToEdge(j,3);
        DCurl.InsertGlobalValues(1, &rowNum, 4, colNum, vals);
    }


// ************************************ CUBATURE **************************************

   // Get numerical integration points and weights
    std::cout << "Getting cubature ... \n\n";

    DefaultCubatureFactory<double>  cubFactory;                                   
    int cubDegree = 2;
    Teuchos::RCP<Cubature<double> > hexCub = cubFactory.create(hex_8, cubDegree); 

    int cubDim       = hexCub->getDimension();
    int numCubPoints = hexCub->getNumPoints();

    FieldContainer<double> cubPoints(numCubPoints, cubDim);
    FieldContainer<double> cubWeights(numCubPoints);

    hexCub->getCubature(cubPoints, cubWeights);


// ************************************** BASIS ***************************************

   // Define basis 
    std::cout << "Getting basis ... \n\n";
    Basis_HCURL_HEX_I1_FEM<double, FieldContainer<double> > hexHCurlBasis;
    Basis_HDIV_HEX_I1_FEM<double, FieldContainer<double> > hexHDivBasis;

    int numFieldsC = hexHCurlBasis.getCardinality();
    int numFieldsD = hexHDivBasis.getCardinality();

  // Evaluate basis at cubature points
     FieldContainer<double> hexCVals(numFieldsC, numCubPoints, spaceDim); 
     FieldContainer<double> hexDVals(numFieldsD, numCubPoints, spaceDim); 
     FieldContainer<double> hexDivs(numFieldsD, numCubPoints); 

     hexHCurlBasis.getValues(hexCVals, cubPoints, OPERATOR_VALUE);
     hexHDivBasis.getValues(hexDVals, cubPoints, OPERATOR_VALUE);
     hexHDivBasis.getValues(hexDivs, cubPoints, OPERATOR_DIV);


// ******** LOOP OVER ELEMENTS TO CREATE LOCAL MASS and STIFFNESS MATRICES *************


    std::cout << "Building mass and stiffness matrices ... \n\n";

 // Settings and data structures for mass and stiffness matrices
    typedef CellTools<double>  CellTools;
    typedef FunctionSpaceTools fst;
    int numCells = 1; 

   // Containers for nodes, edge and face signs 
    FieldContainer<double> hexNodes(numCells, numNodesPerElem, spaceDim);
    FieldContainer<double> hexEdgeSigns(numCells, numFieldsC);
    FieldContainer<double> hexFaceSigns(numCells, numFieldsD);
   // Containers for Jacobian
    FieldContainer<double> hexJacobian(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> hexJacobInv(numCells, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> hexJacobDet(numCells, numCubPoints);
   // Containers for element HCURL mass matrix
    FieldContainer<double> massMatrixC(numCells, numFieldsC, numFieldsC);
    FieldContainer<double> weightedMeasure(numCells, numCubPoints);
    FieldContainer<double> weightedMeasureMuInv(numCells, numCubPoints);
    FieldContainer<double> hexCValsTransformed(numCells, numFieldsC, numCubPoints, spaceDim);
    FieldContainer<double> hexCValsTransformedWeighted(numCells, numFieldsC, numCubPoints, spaceDim);
   // Containers for element HDIV mass matrix
    FieldContainer<double> massMatrixD(numCells, numFieldsD, numFieldsD);
    FieldContainer<double> weightedMeasureMu(numCells, numCubPoints);    
    FieldContainer<double> hexDValsTransformed(numCells, numFieldsD, numCubPoints, spaceDim);
    FieldContainer<double> hexDValsTransformedWeighted(numCells, numFieldsD, numCubPoints, spaceDim);
   // Containers for element HDIV stiffness matrix
    FieldContainer<double> stiffMatrixD(numCells, numFieldsD, numFieldsD);
    FieldContainer<double> hexDivsTransformed(numCells, numFieldsD, numCubPoints);
    FieldContainer<double> hexDivsTransformedWeighted(numCells, numFieldsD, numCubPoints);
   // Containers for right hand side vectors
    FieldContainer<double> rhsDatag(numCells, numCubPoints, cubDim);
    FieldContainer<double> rhsDatah(numCells, numCubPoints);
    FieldContainer<double> gD(numCells, numFieldsD);
    FieldContainer<double> hD(numCells, numFieldsD);
   // Container for cubature points in physical space
    FieldContainer<double> physCubPoints(numCells,numCubPoints, cubDim);

    
   // Global arrays in Epetra format
    Epetra_FECrsMatrix MassC(Copy, globalMapC, numFieldsC);
    Epetra_FECrsMatrix MassD(Copy, globalMapD, numFieldsD);
    Epetra_FECrsMatrix StiffD(Copy, globalMapD, numFieldsD);
    Epetra_FEVector rhsD(globalMapD);

    ofstream fSignsout("faceSigns.dat");

 // *** Element loop ***
    for (int k=0; k<numElems; k++) {

     // Physical cell coordinates
      for (int i=0; i<numNodesPerElem; i++) {
         hexNodes(0,i,0) = nodeCoord(elemToNode(k,i),0);
         hexNodes(0,i,1) = nodeCoord(elemToNode(k,i),1);
         hexNodes(0,i,2) = nodeCoord(elemToNode(k,i),2);
      }

     // Face signs
      for (int j=0; j<numFacesPerElem; j++) {
         hexFaceSigns(0,j) = -1.0;
         for (int i=0; i<numNodesPerFace; i++) {
           int indf=i+1;
           if (indf > numNodesPerFace) indf=0;
           if (elemToNode(k,refFaceToNode(j,0))==faceToNode(elemToFace(k,j),i) &&
               elemToNode(k,refFaceToNode(j,1))==faceToNode(elemToFace(k,j),indf))
                hexFaceSigns(0,j) = 1.0;
          }
         fSignsout << hexFaceSigns(0,j) << "  ";
       }
       fSignsout << "\n";

     // Edge signs
      for (int j=0; j<numEdgesPerElem; j++) {
          if (elemToNode(k,refEdgeToNode(j,0))==edgeToNode(elemToEdge(k,j),0) &&
              elemToNode(k,refEdgeToNode(j,1))==edgeToNode(elemToEdge(k,j),1))
              hexEdgeSigns(0,j) = 1.0;
          else 
              hexEdgeSigns(0,j) = -1.0;
       }

    // Compute cell Jacobians, their inverses and their determinants
       CellTools::setJacobian(hexJacobian, cubPoints, hexNodes, hex_8);
       CellTools::setJacobianInv(hexJacobInv, hexJacobian );
       CellTools::setJacobianDet(hexJacobDet, hexJacobian );


// ************************** Compute element HCurl mass matrices *******************************

     // transform to physical coordinates 
      fst::HCURLtransformVALUE<double>(hexCValsTransformed, hexJacobInv, 
                                   hexCVals);

    // compute weighted measure
      fst::computeMeasure<double>(weightedMeasure, hexJacobDet, cubWeights);

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
      int err = 0;
      for (int row = 0; row < numFieldsC; row++){
        for (int col = 0; col < numFieldsC; col++){
            int rowIndex = elemToEdge(k,row);
            int colIndex = elemToEdge(k,col);
            double val = massMatrixC(0,row,col);
            err = MassC.SumIntoGlobalValues(1, &rowIndex, 1, &colIndex, &val);
            if (err > 0) {
                MassC.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
            }
         }
      }

// ************************** Compute element HDiv mass matrices *******************************

     // transform to physical coordinates 
      fst::HDIVtransformVALUE<double>(hexDValsTransformed, hexJacobian, hexJacobDet,
                                   hexDVals);

     // multiply by weighted measure
      fst::multiplyMeasure<double>(hexDValsTransformedWeighted,
                                   weightedMeasure, hexDValsTransformed);

     // integrate to compute element mass matrix
      fst::integrate<double>(massMatrixD,
                             hexDValsTransformed, hexDValsTransformedWeighted,
                             COMP_CPP);
     // apply face signs
      fst::applyLeftFieldSigns<double>(massMatrixD, hexFaceSigns);
      fst::applyRightFieldSigns<double>(massMatrixD, hexFaceSigns);

     // assemble into global matrix
      err = 0;
      for (int row = 0; row < numFieldsD; row++){
        for (int col = 0; col < numFieldsD; col++){
            int rowIndex = elemToFace(k,row);
            int colIndex = elemToFace(k,col);
            double val = massMatrixD(0,row,col);
            err = MassD.SumIntoGlobalValues(1, &rowIndex, 1, &colIndex, &val);
            if (err > 0) {
                MassD.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
            }
         }
      }

// ************************ Compute element HDiv stiffness matrices *****************************

      // transform to physical coordinates 
      fst::HDIVtransformDIV<double>(hexDivsTransformed, hexJacobDet,
                                    hexDivs);

     // multiply by weighted measure
      fst::multiplyMeasure<double>(hexDivsTransformedWeighted,
                                   weightedMeasure, hexDivsTransformed);

     // integrate to compute element stiffness matrix
      fst::integrate<double>(stiffMatrixD,
                             hexDivsTransformed, hexDivsTransformedWeighted,
                             COMP_CPP);

     // apply face signs
      fst::applyLeftFieldSigns<double>(stiffMatrixD, hexFaceSigns);
      fst::applyRightFieldSigns<double>(stiffMatrixD, hexFaceSigns);

     // assemble into global matrix
      err = 0;
      for (int row = 0; row < numFieldsD; row++){
        for (int col = 0; col < numFieldsD; col++){
            int rowIndex = elemToFace(k,row);
            int colIndex = elemToFace(k,col);
            double val = stiffMatrixD(0,row,col);
            err = StiffD.SumIntoGlobalValues(1, &rowIndex, 1, &colIndex, &val);
            if (err > 0) {
                StiffD.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
            }
         }
      }

// ******************************* Build right hand side ************************************

      // transform integration points to physical points
       FieldContainer<double> physCubPoints(numCells,numCubPoints, cubDim);
       CellTools::mapToPhysicalFrame(physCubPoints, cubPoints, hexNodes, hex_8);

      // evaluate right hand side functions at physical points
       FieldContainer<double> rhsDatag(numCells, numCubPoints, cubDim);
       FieldContainer<double> rhsDatah(numCells, numCubPoints);
       for (int nPt = 0; nPt < numCubPoints; nPt++){

          double x = physCubPoints(0,nPt,0);
          double y = physCubPoints(0,nPt,1);
          double z = physCubPoints(0,nPt,2);
          double du1, du2, du3;

          evalCurlCurlu(du1, du2, du3, x, y, z);
          rhsDatag(0,nPt,0) = du1;
          rhsDatag(0,nPt,1) = du2;
          rhsDatag(0,nPt,2) = du3;
         
          rhsDatah(0,nPt) = evalDivu(x, y, z);
       }

     // integrate (g,curl w) term
      fst::integrate<double>(gD, rhsDatag, hexDValsTransformedWeighted,
                             COMP_CPP);

     // integrate (h,div w) term
      fst::integrate<double>(hD, rhsDatah, hexDivsTransformedWeighted,
                             COMP_CPP);

    // assemble into global vector
     for (int row = 0; row < numFieldsD; row++){
           int rowIndex = elemToEdge(k,row);
           double val = gD(0,row)-hD(0,row);
           err = rhsD.SumIntoGlobalValues(1, &rowIndex, &val);
     }
 
     
 } // *** end element loop ***

  // Assemble over multiple processors, if necessary
   DCurl.GlobalAssemble(); DCurl.FillComplete(MassC.RowMap(),MassD.RowMap()); 
   MassC.GlobalAssemble();  MassC.FillComplete();
   MassD.GlobalAssemble();  MassD.FillComplete();
   StiffD.GlobalAssemble(); StiffD.FillComplete();
   rhsD.GlobalAssemble();
   
  // Dump matrices to disk
   EpetraExt::RowMatrixToMatlabFile("mag_m1_matrix.dat",MassC);
   EpetraExt::RowMatrixToMatlabFile("mag_m2_matrix.dat",MassD);
   EpetraExt::RowMatrixToMatlabFile("mag_k2_matrix.dat",StiffD);
   EpetraExt::RowMatrixToMatlabFile("mag_t1_matrix.dat",DCurl);
   EpetraExt::MultiVectorToMatlabFile("rhs2_vector.dat",rhsD);

   fSignsout.close();

 // delete mesh
 Delete_Pamgen_Mesh();

 // reset format state of std::cout
 std::cout.copyfmt(oldFormatState);
 
 return 0;
}

// Calculates value of exact solution u
 int evalu(double & uExact0, double & uExact1, double & uExact2, double & x, double & y, double & z)
 {
    uExact0 = cos(M_PI*y)*cos(M_PI*z)*(x+1.0)*(x-1.0);
    uExact1 = cos(M_PI*x)*cos(M_PI*z)*(y+1.0)*(y-1.0);
    uExact2 = cos(M_PI*x)*cos(M_PI*y)*(z+1.0)*(z-1.0);

   return 0;
 }

// Calculates divergence of exact solution u
 double evalDivu(double & x, double & y, double & z)
 {
   double divu = 2.0*x*cos(M_PI*y)*cos(M_PI*z) + 2.0*y*cos(M_PI*x)*cos(M_PI*z)
                  + 2.0*z*cos(M_PI*x)*cos(M_PI*y);
   return divu;
 }


// Calculates curl of exact solution u
 int evalCurlu(double & curlu0, double & curlu1, double & curlu2, double & x, double & y, double & z)
 {
   double duxdy = -M_PI*sin(M_PI*y)*cos(M_PI*z)*(x+1.0)*(x-1.0);
   double duxdz = -M_PI*sin(M_PI*z)*cos(M_PI*y)*(x+1.0)*(x-1.0);
   double duydx = -M_PI*sin(M_PI*x)*cos(M_PI*z)*(y+1.0)*(y-1.0);
   double duydz = -M_PI*sin(M_PI*z)*cos(M_PI*x)*(y+1.0)*(y-1.0);
   double duzdx = -M_PI*sin(M_PI*x)*cos(M_PI*y)*(z+1.0)*(z-1.0);
   double duzdy = -M_PI*sin(M_PI*y)*cos(M_PI*x)*(z+1.0)*(z-1.0);

   curlu0 = duzdy - duydz;
   curlu1 = duxdz - duzdx;
   curlu2 = duydx - duxdy;

   return 0;
 }

// Calculates curl of the curl of exact solution u
 int evalCurlCurlu(double & curlCurlu0, double & curlCurlu1, double & curlCurlu2, double & x, double & y, double & z)
{
    double dcurlu0dy = -M_PI*M_PI*cos(M_PI*y)*cos(M_PI*x)*(z+1.0)*(z-1.0)
                           + 2.0*y*M_PI*sin(M_PI*z)*cos(M_PI*x);
    double dcurlu0dz = -2.0*z*M_PI*sin(M_PI*y)*cos(M_PI*x)
                          + M_PI*M_PI*cos(M_PI*z)*cos(M_PI*x)*(y+1.0)*(y-1.0);
    double dcurlu1dx = -2.0*x*M_PI*sin(M_PI*z)*cos(M_PI*y)
                          + M_PI*M_PI*cos(M_PI*x)*cos(M_PI*y)*(z+1.0)*(z-1.0);
    double dcurlu1dz = -M_PI*M_PI*cos(M_PI*z)*cos(M_PI*y)*(x+1.0)*(x-1.0)
                           + 2.0*z*M_PI*sin(M_PI*x)*cos(M_PI*y);
    double dcurlu2dx = -M_PI*M_PI*cos(M_PI*x)*cos(M_PI*z)*(y+1.0)*(y-1.0)
                           + 2.0*x*M_PI*sin(M_PI*y)*cos(M_PI*z);
    double dcurlu2dy = -2.0*y*M_PI*sin(M_PI*x)*cos(M_PI*z)
                          + M_PI*M_PI*cos(M_PI*y)*cos(M_PI*z)*(x+1.0)*(x-1.0);
                       
    curlCurlu0 = dcurlu2dy - dcurlu1dz;
    curlCurlu1 = dcurlu0dz - dcurlu2dx;
    curlCurlu2 = dcurlu1dx - dcurlu0dy;

    return 0;
}
