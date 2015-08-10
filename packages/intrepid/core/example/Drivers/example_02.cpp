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

/** \file   example_02.cpp
    \brief  Example building mass and stiffness matrices and right hand side
            for a div-curl system on a hexahedral mesh using div-conforming
            (face) elements.

    \verbatim
                       curl u = g  in Omega
                        div u = h  in Omega
                          u.n = 0  on Gamma

            Discrete linear system for face element coeficients (x):

                      (Kd + Md*Dc*McInv*Dc'*Md)x = b

                      Kd    - Hdiv stiffness matrix
                      Md    - Hdiv mass matrix
                      Dc    - Edge to Face incidence matrix
                      McInv - Hcurl mass matrix inverse
                      b     - right hand side vector

    \endverbatim

    \author Created by P. Bochev, D. Ridzal and K. Peterson.
 
     \remark Usage
     \verbatim

     ./Intrepid_example_Drivers_Example_02.exe NX NY NZ randomMesh mu1 mu2 mu1LX mu1RX mu1LY mu1RY mu1LZ mu1RZ verbose

        int NX              - num intervals in x direction (assumed box domain, -1,1) 
        int NY              - num intervals in y direction (assumed box domain, -1,1)
        int NZ              - num intervals in z direction (assumed box domain, -1,1)
        int randomMesh      - 1 if mesh randomizer is to be used 0 if not 
        double mu1          - material property value for region 1 
        double mu2          - material property value for region 2 
        double mu1LX        - left X boundary for region 1 
        double mu1RX        - right X boundary for region 1 
        double mu1LY        - left Y boundary for region 1 
        double mu1RY        - right Y boundary for region 1 
        double mu1LZ        - bottom Z boundary for region 1 
        double mu1RZ        - top Z boundary for region 1 
        verbose (optional)  - any character, indicates verbose output 

     \endverbatim

    \remark Sample command line
    \code   ./Intrepid_example_Drivers_Example_02.exe 10 10 10 0 1.0 1.0 -1.0 1.0 -1.0 1.0 -1.0 1.0 \endcode
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
#include "Epetra_Vector.h"

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

// Functions to evaluate exact solution and its derivatives
int evalu(double & uExact0, double & uExact1, double & uExact2, double & x, double & y, double & z);
double evalDivu(double & x, double & y, double & z);
int evalCurlu(double & curlu0, double & curlu1, double & curlu2, double & x, double & y, double & z);
int evalCurlCurlu(double & curlCurlu0, double & curlCurlu1, double & curlCurlu2, double & x, double & y, double & z);

int main(int argc, char *argv[]) {

   //Check number of arguments
   if (argc < 13) {
      std::cout <<"\n>>> ERROR: Invalid number of arguments.\n\n";
      std::cout <<"Usage:\n\n";
      std::cout <<"  ./Intrepid_example_Drivers_Example_02.exe NX NY NZ randomMesh mu1 mu2 mu1LX mu1RX mu1LY mu1RY mu1LZ mu1RZ verbose\n\n";
      std::cout <<" where \n";
      std::cout <<"   int NX              - num intervals in x direction (assumed box domain, -1,1) \n";
      std::cout <<"   int NY              - num intervals in y direction (assumed box domain, -1,1) \n";
      std::cout <<"   int NZ              - num intervals in z direction (assumed box domain, -1,1) \n";
      std::cout <<"   int randomMesh      - 1 if mesh randomizer is to be used 0 if not \n";
      std::cout <<"   double mu1          - material property value for region 1 \n";
      std::cout <<"   double mu2          - material property value for region 2 \n";
      std::cout <<"   double mu1LX        - left X boundary for region 1 \n";
      std::cout <<"   double mu1RX        - right X boundary for region 1 \n";
      std::cout <<"   double mu1LY        - left Y boundary for region 1 \n";
      std::cout <<"   double mu1RY        - right Y boundary for region 1 \n";
      std::cout <<"   double mu1LZ        - bottom Z boundary for region 1 \n";
      std::cout <<"   double mu1RZ        - top Z boundary for region 1 \n";
      std::cout <<"   verbose (optional)  - any character, indicates verbose output \n\n";
      exit(1);
   }

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
    << "|   Example: Generate Mass and Stiffness Matrices and Right-Hand Side Vector  |\n"
    << "|     for Div-Curl System on Hexahedral Mesh with Div-Conforming Elements     |\n" \
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

   // Build reference element face to node map
    FieldContainer<int> refFaceToNode(numFacesPerElem,numNodesPerFace);
    for (int i=0; i<numFacesPerElem; i++){
        refFaceToNode(i,0)=hex_8.getNodeMap(2, i, 0);
        refFaceToNode(i,1)=hex_8.getNodeMap(2, i, 1);
        refFaceToNode(i,2)=hex_8.getNodeMap(2, i, 2);
        refFaceToNode(i,3)=hex_8.getNodeMap(2, i, 3);
    }

// *********************************** GENERATE MESH ************************************

    *outStream << "Generating mesh ... \n\n";

    *outStream << "    NX" << "   NY" << "   NZ\n";
    *outStream << std::setw(5) << NX <<
                 std::setw(5) << NY <<
                 std::setw(5) << NZ << "\n\n";

   // Print mesh information
    int numElems = NX*NY*NZ;
    int numNodes = (NX+1)*(NY+1)*(NZ+1);
    int numEdges = (NX)*(NY + 1)*(NZ + 1) + (NX + 1)*(NY)*(NZ + 1) + (NX + 1)*(NY + 1)*(NZ);
    int numFaces = (NX)*(NY)*(NZ + 1) + (NX)*(NY + 1)*(NZ) + (NX + 1)*(NY)*(NZ);
    *outStream << " Number of Elements: " << numElems << " \n";
    *outStream << "    Number of Nodes: " << numNodes << " \n";
    *outStream << "    Number of Edges: " << numEdges << " \n";
    *outStream << "    Number of Faces: " << numFaces << " \n\n";

   // Cube
    double leftX = -1.0, rightX = 1.0;
    double leftY = -1.0, rightY = 1.0;
    double leftZ = -1.0, rightZ = 1.0;

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
          inode++;
        }
      }
    }

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

  // Get edge connectivity
    FieldContainer<int> edgeToNode(numEdges, numNodesPerEdge);
    FieldContainer<int> elemToEdge(numElems, numEdgesPerElem);
    int iedge = 0;
    inode = 0;
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

   // Find boundary edges
    FieldContainer<int> edgeOnBoundary(numEdges);
    for (int i=0; i<numEdges; i++){
       if (nodeOnBoundary(edgeToNode(i,0)) && nodeOnBoundary(edgeToNode(i,1))){
           edgeOnBoundary(i)=1;
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

   // Find boundary faces
    FieldContainer<int> faceOnBoundary(numFaces);
    for (int i=0; i<numFaces; i++){
       if (nodeOnBoundary(faceToNode(i,0)) && nodeOnBoundary(faceToNode(i,1))
          && nodeOnBoundary(faceToNode(i,2)) && nodeOnBoundary(faceToNode(i,3))){
           faceOnBoundary(i)=1;
       }
    }

#define DUMP_DATA
#ifdef DUMP_DATA
   // Output connectivity
    ofstream fe2nout("elem2node.dat");
    ofstream fe2fout("elem2face.dat");
    for (int k=0; k<NZ; k++) {
      for (int j=0; j<NY; j++) {
        for (int i=0; i<NX; i++) {
          int ielem = i + j * NX + k * NX * NY;
          for (int m=0; m<numNodesPerElem; m++){
              fe2nout << elemToNode(ielem,m) <<"  ";
           }
          fe2nout <<"\n";
          for (int n=0; n<numFacesPerElem; n++) {
             fe2fout << elemToFace(ielem,n) << "  ";
          }
          fe2fout << "\n";
        }
      }
    }
    fe2nout.close();
    fe2fout.close();
#endif

#ifdef DUMP_DATA_EXTRA
    ofstream ff2nout("face2node.dat");
    ofstream ff2eout("face2edge.dat");
    for (int i=0; i<numFaces; i++) {
       for (int j=0; j<numNodesPerFace; j++) {
           ff2nout << faceToNode(i,j) << "  ";
       }
       for (int k=0; k<numEdgesPerFace; k++) {
           ff2eout << faceToEdge(i,k) << "  ";
       }
       ff2nout << "\n";
       ff2eout << "\n";
    }
    ff2nout.close();
    ff2eout.close();

    ofstream fBnodeout("nodeOnBndy.dat");
    ofstream fBfaceout("faceOnBndy.dat");
    for (int i=0; i<numNodes; i++) {
        fBnodeout << nodeOnBoundary(i) <<"\n";
    }
    for (int i=0; i<numFaces; i++) {
        fBfaceout << faceOnBoundary(i) <<"\n";
    }
    fBnodeout.close();
    fBfaceout.close();
#endif

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
   // Print nodal coords
    ofstream fcoordout("coords.dat");
    for (int i=0; i<numNodes; i++) {
       fcoordout << nodeCoord(i,0) <<" ";
       fcoordout << nodeCoord(i,1) <<" ";
       fcoordout << nodeCoord(i,2) <<"\n";
    }
    fcoordout.close();
#endif


// **************************** INCIDENCE MATRIX **************************************

   // Edge to face incidence matrix
    *outStream << "Building incidence matrix ... \n\n";

    Epetra_SerialComm Comm;
    Epetra_Map globalMapD(numFaces, 0, Comm);
    Epetra_Map globalMapC(numEdges, 0, Comm);
    Epetra_FECrsMatrix DCurl(Copy, globalMapD, globalMapC, 2);

    double vals[4];
    vals[0]=1.0; vals[1]=1.0; vals[2]=-1.0; vals[3]=-1.0;
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
    *outStream << "Getting cubature ... \n\n";

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
    *outStream << "Getting basis ... \n\n";
    Basis_HCURL_HEX_I1_FEM<double, FieldContainer<double> > hexHCurlBasis;
    Basis_HDIV_HEX_I1_FEM<double, FieldContainer<double> > hexHDivBasis;

    int numFieldsC = hexHCurlBasis.getCardinality();
    int numFieldsD = hexHDivBasis.getCardinality();

  // Evaluate basis at cubature points
     FieldContainer<double> hexCVals(numFieldsC, numCubPoints, spaceDim); 
     FieldContainer<double> hexDVals(numFieldsD, numCubPoints, spaceDim); 
     FieldContainer<double> hexDivs(numFieldsD, numCubPoints); 
     FieldContainer<double> worksetDVals(numFieldsD, numFacePoints, spaceDim);

     hexHCurlBasis.getValues(hexCVals, cubPoints, OPERATOR_VALUE);
     hexHDivBasis.getValues(hexDVals, cubPoints, OPERATOR_VALUE);
     hexHDivBasis.getValues(hexDivs, cubPoints, OPERATOR_DIV);


// ******** LOOP OVER ELEMENTS TO CREATE LOCAL MASS and STIFFNESS MATRICES *************


    *outStream << "Building mass and stiffness matrices ... \n\n";

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
    FieldContainer<double> gDBoundary(numCells, numFieldsD);
    FieldContainer<double> refGaussPoints(numFacePoints,spaceDim);
    FieldContainer<double> worksetGaussPoints(numCells,numFacePoints,spaceDim);
    FieldContainer<double> worksetJacobians(numCells, numFacePoints, spaceDim, spaceDim);
    FieldContainer<double> worksetJacobDet(numCells, numFacePoints);
    FieldContainer<double> worksetFaceN(numCells, numFacePoints, spaceDim);
    FieldContainer<double> worksetVFieldVals(numCells, numFacePoints, spaceDim);
    FieldContainer<double> worksetDValsTransformed(numCells, numFieldsD, numFacePoints, spaceDim);
    FieldContainer<double> curluFace(numCells, numFacePoints, spaceDim);
    FieldContainer<double> worksetDataCrossField(numCells, numFieldsD, numFacePoints, spaceDim);
   // Container for cubature points in physical space
    FieldContainer<double> physCubPoints(numCells,numCubPoints, cubDim);

    
   // Global arrays in Epetra format
    Epetra_FECrsMatrix MassC(Copy, globalMapC, numFieldsC);
    Epetra_FECrsMatrix MassD(Copy, globalMapD, numFieldsD);
    Epetra_FECrsMatrix StiffD(Copy, globalMapD, numFieldsD);
    Epetra_FEVector rhsD(globalMapD);

#ifdef DUMP_DATA
    ofstream fSignsout("faceSigns.dat");
#endif

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
#ifdef DUMP_DATA
         fSignsout << hexFaceSigns(0,j) << "  ";
#endif
       }
#ifdef DUMP_DATA
       fSignsout << "\n";
#endif

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
      fst::computeCellMeasure<double>(weightedMeasure, hexJacobDet, cubWeights);

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
      for (int row = 0; row < numFieldsC; row++){
        for (int col = 0; col < numFieldsC; col++){
            int rowIndex = elemToEdge(k,row);
            int colIndex = elemToEdge(k,col);
            double val = massMatrixC(0,row,col);
            MassC.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
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
      for (int row = 0; row < numFieldsD; row++){
        for (int col = 0; col < numFieldsD; col++){
            int rowIndex = elemToFace(k,row);
            int colIndex = elemToFace(k,col);
            double val = massMatrixD(0,row,col);
            MassD.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
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
      for (int row = 0; row < numFieldsD; row++){
        for (int col = 0; col < numFieldsD; col++){
            int rowIndex = elemToFace(k,row);
            int colIndex = elemToFace(k,col);
            double val = stiffMatrixD(0,row,col);
            StiffD.InsertGlobalValues(1, &rowIndex, 1, &colIndex, &val);
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

          evalCurlCurlu(du1, du2, du3, x, y, z);
          rhsDatag(0,nPt,0) = du1;
          rhsDatag(0,nPt,1) = du2;
          rhsDatag(0,nPt,2) = du3;
         
          rhsDatah(0,nPt) = evalDivu(x, y, z);
       }

     // integrate (curl g, w) term
      fst::integrate<double>(gD, rhsDatag, hexDValsTransformedWeighted,
                             COMP_CPP);

     // integrate (h,div w) term
      fst::integrate<double>(hD, rhsDatah, hexDivsTransformedWeighted,
                             COMP_CPP);

     // apply signs
      fst::applyFieldSigns<double>(gD, hexFaceSigns);
      fst::applyFieldSigns<double>(hD, hexFaceSigns);

     // calculate boundary term, (g x w, n)_{\Gamma}
      for (int i = 0; i < numFacesPerElem; i++){
        if (faceOnBoundary(elemToFace(k,i))){

         // map Gauss points on quad to reference face: paramGaussPoints -> refGaussPoints
            CellTools::mapToReferenceSubcell(refGaussPoints,
                                   paramGaussPoints,
                                   2, i, hex_8);

         // get basis values at points on reference cell
           hexHDivBasis.getValues(worksetDVals, refGaussPoints, OPERATOR_VALUE);

         // compute Jacobians at Gauss pts. on reference face for all parent cells
           CellTools::setJacobian(worksetJacobians, refGaussPoints,
                         hexNodes, hex_8);
           CellTools::setJacobianDet(worksetJacobDet, worksetJacobians);

         // transform to physical coordinates
            fst::HDIVtransformVALUE<double>(worksetDValsTransformed, worksetJacobians,
                                   worksetJacobDet, worksetDVals);

         // map Gauss points on quad from ref. face to face workset: refGaussPoints -> worksetGaussPoints
            CellTools::mapToPhysicalFrame(worksetGaussPoints,
                                refGaussPoints,
                                hexNodes, hex_8);

        // compute face normals
            CellTools::getPhysicalFaceNormals(worksetFaceN,
                                              worksetJacobians,
                                              i, hex_8);

         // evaluate curl u at face points
           for(int nPt = 0; nPt < numFacePoints; nPt++){

             double x = worksetGaussPoints(0, nPt, 0);
             double y = worksetGaussPoints(0, nPt, 1);
             double z = worksetGaussPoints(0, nPt, 2);

             evalCurlu(curluFace(0,nPt,0), curluFace(0,nPt,1), curluFace(0,nPt,2), x, y, z);
           }

          // compute the cross product of curluFace with basis and multiply by weights
           for (int nF = 0; nF < numFieldsD; nF++){
              for(int nPt = 0; nPt < numFacePoints; nPt++){
                  worksetDataCrossField(0,nF,nPt,0) = (curluFace(0,nPt,1)*worksetDValsTransformed(0,nF,nPt,2)
                                 - curluFace(0,nPt,2)*worksetDValsTransformed(0,nF,nPt,1))
                                  * paramGaussWeights(nPt);
                  worksetDataCrossField(0,nF,nPt,1) = (curluFace(0,nPt,2)*worksetDValsTransformed(0,nF,nPt,0)
                                 - curluFace(0,nPt,0)*worksetDValsTransformed(0,nF,nPt,2))
                                  * paramGaussWeights(nPt);
                  worksetDataCrossField(0,nF,nPt,2) = (curluFace(0,nPt,0)*worksetDValsTransformed(0,nF,nPt,1)
                                 - curluFace(0,nPt,1)*worksetDValsTransformed(0,nF,nPt,0))
                                  *paramGaussWeights(nPt);
              } //nPt
           } //nF

          // Integrate
           fst::integrate<double>(gDBoundary, worksetFaceN, worksetDataCrossField,
                             COMP_CPP);

          // apply signs
           fst::applyFieldSigns<double>(gDBoundary, hexFaceSigns);

          // add into gD term
            for (int nF = 0; nF < numFieldsD; nF++){
                gD(0,nF) = gD(0,nF) - gDBoundary(0,nF);
            }

        } // if faceOnBoundary
      } // numFaces


    // assemble into global vector
     for (int row = 0; row < numFieldsD; row++){
           int rowIndex = elemToFace(k,row);
           double val = hD(0,row)+gD(0,row);
           rhsD.SumIntoGlobalValues(1, &rowIndex, &val);
     }
 
     
 } // *** end element loop ***

  // Assemble over multiple processors, if necessary
   DCurl.GlobalAssemble(); DCurl.FillComplete(MassC.RowMap(),MassD.RowMap());  
   MassC.GlobalAssemble();  MassC.FillComplete();
   MassD.GlobalAssemble();  MassD.FillComplete();
   StiffD.GlobalAssemble(); StiffD.FillComplete();
   rhsD.GlobalAssemble();
   
#ifdef DUMP_DATA
   fSignsout.close();
#endif

  // Build the inverse diagonal for MassC
   Epetra_CrsMatrix MassCinv(Copy,MassC.RowMap(),MassC.RowMap(),1);
   Epetra_Vector DiagC(MassC.RowMap());

   DiagC.PutScalar(1.0);
   MassC.Multiply(false,DiagC,DiagC);
   for(int i=0; i<DiagC.MyLength(); i++) {
     DiagC[i]=1.0/DiagC[i];
   }
   for(int i=0; i<DiagC.MyLength(); i++) {
     int CID=MassC.GCID(i);
     MassCinv.InsertGlobalValues(MassC.GRID(i),1,&(DiagC[i]),&CID);
   }
   MassCinv.FillComplete();

  // Set value to zero on diagonal that cooresponds to boundary edge
   for(int i=0;i<numEdges;i++) {
     if (edgeOnBoundary(i)){
      double val=0.0;
      MassCinv.ReplaceGlobalValues(i,1,&val,&i);
     }
   }

    int numEntries;
    double *values;
    int *cols;

  // Adjust matrices and rhs due to boundary conditions
   for (int row = 0; row<numFaces; row++){
      MassD.ExtractMyRowView(row,numEntries,values,cols);
        for (int i=0; i<numEntries; i++){
           if (faceOnBoundary(cols[i])) {
             values[i]=0;
          }
       }
      StiffD.ExtractMyRowView(row,numEntries,values,cols);
        for (int i=0; i<numEntries; i++){
           if (faceOnBoundary(cols[i])) {
             values[i]=0;
          }
       }
    }
   for (int row = 0; row<numFaces; row++){
       if (faceOnBoundary(row)) {
          int rowindex = row;
          StiffD.ExtractMyRowView(row,numEntries,values,cols);
          for (int i=0; i<numEntries; i++){
             values[i]=0;
          }
          MassD.ExtractMyRowView(row,numEntries,values,cols);
          for (int i=0; i<numEntries; i++){
             values[i]=0;
          }
         rhsD[0][row]=0;
         double val = 1.0;
         StiffD.ReplaceGlobalValues(1, &rowindex, 1, &rowindex, &val);
       }
    }

#ifdef DUMP_DATA
  // Dump matrices to disk
   EpetraExt::RowMatrixToMatlabFile("mag_m1inv_matrix.dat",MassCinv);
   EpetraExt::RowMatrixToMatlabFile("mag_m2_matrix.dat",MassD);
   EpetraExt::RowMatrixToMatlabFile("mag_k2_matrix.dat",StiffD);
   EpetraExt::RowMatrixToMatlabFile("mag_t1_matrix.dat",DCurl);
   EpetraExt::MultiVectorToMatrixMarketFile("mag_rhs2_vector.dat",rhsD,0,0,false);
#endif

   std::cout << "End Result: TEST PASSED\n";

 // reset format state of std::cout
 std::cout.copyfmt(oldFormatState);
 
 return 0;
}
// Calculates value of exact solution u
 int evalu(double & uExact0, double & uExact1, double & uExact2, double & x, double & y, double & z)
 {

   // function 1
    uExact0 = exp(y+z)*(x+1.0)*(x-1.0);
    uExact1 = exp(x+z)*(y+1.0)*(y-1.0);
    uExact2 = exp(x+y)*(z+1.0)*(z-1.0);

 /*
   // function 2
    uExact0 = cos(M_PI*y)*cos(M_PI*z)*(x+1.0)*(x-1.0);
    uExact1 = cos(M_PI*x)*cos(M_PI*z)*(y+1.0)*(y-1.0);
    uExact2 = cos(M_PI*x)*cos(M_PI*y)*(z+1.0)*(z-1.0);

 */
 /*
   // function 3
    uExact0 = x*x-1.0;
    uExact1 = y*y-1.0;
    uExact2 = z*z-1.0;

   // function 4
    uExact0 = sin(M_PI*x);
    uExact1 = sin(M_PI*y);
    uExact2 = sin(M_PI*z);
 */

   return 0;
 }
// Calculates divergence of exact solution u
 double evalDivu(double & x, double & y, double & z)
 {

   // function 1
    double divu = 2.0*x*exp(y+z)+2.0*y*exp(x+z)+2.0*z*exp(x+y);

   // function 2
  //  double divu = 2.0*x*cos(M_PI*y)*cos(M_PI*z) + 2.0*y*cos(M_PI*x)*cos(M_PI*z)
  //                 + 2.0*z*cos(M_PI*x)*cos(M_PI*y);

   // function 3
   // double divu = 2.0*(x + y + z);

   // function 4
   // double divu = M_PI*(cos(M_PI*x)+cos(M_PI*y)+cos(M_PI*z));

   return divu;
 }
// Calculates curl of exact solution u
 int evalCurlu(double & curlu0, double & curlu1, double & curlu2, double & x, double & y, double & z)
 {

   // function 1
    double duxdy = exp(y+z)*(x+1.0)*(x-1.0);
    double duxdz = exp(y+z)*(x+1.0)*(x-1.0);
    double duydx = exp(x+z)*(y+1.0)*(y-1.0);
    double duydz = exp(x+z)*(y+1.0)*(y-1.0);
    double duzdx = exp(x+y)*(z+1.0)*(z-1.0);
    double duzdy = exp(x+y)*(z+1.0)*(z-1.0);

  /*

   // function 2
    double duxdy = -M_PI*sin(M_PI*y)*cos(M_PI*z)*(x+1.0)*(x-1.0);
    double duxdz = -M_PI*sin(M_PI*z)*cos(M_PI*y)*(x+1.0)*(x-1.0);
    double duydx = -M_PI*sin(M_PI*x)*cos(M_PI*z)*(y+1.0)*(y-1.0);
    double duydz = -M_PI*sin(M_PI*z)*cos(M_PI*x)*(y+1.0)*(y-1.0);
    double duzdx = -M_PI*sin(M_PI*x)*cos(M_PI*y)*(z+1.0)*(z-1.0);
    double duzdy = -M_PI*sin(M_PI*y)*cos(M_PI*x)*(z+1.0)*(z-1.0);
  */

    curlu0 = duzdy - duydz;
    curlu1 = duxdz - duzdx;
    curlu2 = duydx - duxdy;

  /*
   // function 3 and 4
    curlu0 = 0;
    curlu1 = 0;
    curlu2 = 0;
  */

   return 0;
 }

// Calculates curl of the curl of exact solution u
 int evalCurlCurlu(double & curlCurlu0, double & curlCurlu1, double & curlCurlu2, double & x, double & y, double & z)
{

   // function 1
    double dcurlu0dy = exp(x+y)*(z+1.0)*(z-1.0) - 2.0*y*exp(x+z);
    double dcurlu0dz = 2.0*z*exp(x+y) - exp(x+z)*(y+1.0)*(y-1.0);
    double dcurlu1dx = 2.0*x*exp(y+z) - exp(x+y)*(z+1.0)*(z-1.0);
    double dcurlu1dz = exp(y+z)*(x+1.0)*(x-1.0) - 2.0*z*exp(x+y);
    double dcurlu2dx = exp(x+z)*(y+1.0)*(y-1.0) - 2.0*x*exp(y+z);
    double dcurlu2dy = 2.0*y*exp(x+z) - exp(y+z)*(x+1.0)*(x-1.0);

 /*
   // function 2
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
 */

    curlCurlu0 = dcurlu2dy - dcurlu1dz;
    curlCurlu1 = dcurlu0dz - dcurlu2dx;
    curlCurlu2 = dcurlu1dx - dcurlu0dy;

 /*
   // function 3 and 4
    curlCurlu0 = 0.0;
    curlCurlu1 = 0.0;
    curlCurlu2 = 0.0;
 */

    return 0;
}


