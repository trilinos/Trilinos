//**************************************************************************
//
//                                 NOTICE
//
// This software is a result of the research described in the report
//
// " A comparison of algorithms for modal analysis in the absence 
//   of a sparse direct method", P. ArbenZ, R. Lehoucq, and U. Hetmaniuk,
//  Sandia National Laboratories, Technical report SAND2003-1028J.
//
// It is based on the Epetra, AztecOO, and ML packages defined in the Trilinos
// framework ( http://software.sandia.gov/trilinos/ ).
//
// The distribution of this software follows also the rules defined in Trilinos.
// This notice shall be marked on anY reproduction of this software, in whole or
// in part.
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// Code Authors: U. Hetmaniuk (ulhetma@sandia.gov), R. Lehoucq (rblehou@sandia.gov)
//
//**************************************************************************

#include "ModeLaplace3DQ2.h"


const int ModeLaplace3DQ2::dofEle = 27;
const int ModeLaplace3DQ2::maxConnect = 125;
#ifndef M_PI
const double ModeLaplace3DQ2::M_PI = 3.14159265358979323846;
#endif


ModeLaplace3DQ2::ModeLaplace3DQ2(const Epetra_Comm &_Comm, double _Lx, int _nX,
                                 double _Ly, int _nY, double _Lz, int _nZ)
                : myVerify(_Comm),
                  MyComm(_Comm),
                  mySort(),
                  Map(0),
                  K(0),
                  M(0),
                  Lx(_Lx),
                  nX(_nX),
                  Ly(_Ly),
                  nY(_nY),
                  Lz(_Lz),
                  nZ(_nZ),
                  x(0),
                  y(0),
                  z(0)
                {

  preProcess();

}


ModeLaplace3DQ2::~ModeLaplace3DQ2() {

  if (Map)
    delete Map;
  Map = 0;

  if (K)
    delete K;
  K = 0;

  if (M)
    delete M;
  M = 0;

  if (x)
    delete[] x;
  x = 0;

  if (y)
    delete[] y;
  y = 0;

  if (z)
    delete[] z;
  z = 0;

}


void ModeLaplace3DQ2::preProcess() {

  // Create the distribution of equations across processors
  makeMap();

  // Count the number of elements touched by this processor
  bool *isTouched = new bool[nX*nY*nZ];
  int i, j, k;
  for (k = 0; k < nZ; ++k)
    for (j = 0; j < nY; ++j)
      for (i=0; i<nX; ++i)
        isTouched[i + j*nX + k*nX*nY] = false;

  int numEle = countElements(isTouched);

  // Create the mesh
  int *elemTopo = new int[dofEle*numEle];
  makeMyElementsTopology(elemTopo, isTouched);

  delete[] isTouched;

  // Get the number of nonZeros per row
  int localSize = Map->NumMyElements();
  int *numNz = new int[localSize];
  int *connectivity = new int[localSize*maxConnect];
  makeMyConnectivity(elemTopo, numEle, connectivity, numNz);

  // Make the stiffness matrix
  makeStiffness(elemTopo, numEle, connectivity, numNz);

  // Assemble the mass matrix
  makeMass(elemTopo, numEle, connectivity, numNz);

  // Free some memory
  delete[] elemTopo;
  delete[] numNz;
  delete[] connectivity;

  // Get the geometrical coordinates of the managed nodes
  double hx = Lx/nX;
  double hy = Ly/nY;
  double hz = Lz/nZ;

  x = new double[localSize];
  y = new double[localSize];
  z = new double[localSize];

  for (k=0; k<2*nZ-1; ++k) {
    for (j=0; j<2*nY-1; ++j) {
      for (i=0; i<2*nX-1; ++i) {
        int node = i + j*(2*nX-1) + k*(2*nX-1)*(2*nY-1);
        if (Map->LID(node) > -1) {
          x[Map->LID(node)] = (i+1)*hx*0.5;
          y[Map->LID(node)] = (j+1)*hy*0.5;
          z[Map->LID(node)] = (k+1)*hz*0.5;
        }
      }
    }
  }

}


void ModeLaplace3DQ2::makeMap() {

  int numProc = MyComm.NumProc();
  int globalSize = (2*nX - 1)*(2*nY - 1)*(2*nZ - 1);
  assert(globalSize > numProc);

#ifdef _USE_CHACO
  // Use the partitioner Chaco to distribute the unknowns
  int *start = new int[globalSize+1];
  memset(start, 0, (globalSize+1)*sizeof(int));

  int i, j, k;
  for (k=0; k<2*nZ-1; ++k) {
    for (j=0; j<2*nY-1; ++j) {
      for (i=0; i<2*nX-1; ++i) {
        int node = i + j*(2*nX-1) + k*(2*nX-1)*(2*nY-1);
        int connect = 1;
        if (i%2 == 1) {
          connect *= ((i == 1) || (i == 2*nX-3)) ? 4 : 5;
        }
        else {
          connect *= ((i == 0) || (i == 2*nX-2)) ? 2 : 3;
        }
        if (j%2 == 1) {
          connect *= ((j == 1) || (j == 2*nY-3)) ? 4 : 5;
        }
        else {
          connect *= ((j == 0) || (j == 2*nY-2)) ? 2 : 3;
        }
        if (k%2 == 1) {
          connect *= ((k == 1) || (k == 2*nZ-3)) ? 4 : 5;
        }
        else {
          connect *= ((k == 0) || (k == 2*nZ-2)) ? 2 : 3;
        }
        // Don't count the node itself for Chaco
        start[node+1] = connect - 1;
      }
    }
  }

  for (i = 0; i < globalSize; ++i)
    start[i+1] += start[i];

  int *adjacency = new int[start[globalSize]];
  memset(adjacency, 0, start[globalSize]*sizeof(int));

  int *elemTopo = new int[dofEle];
  for (k=0; k<nZ; ++k) {
    for (j=0; j<nY; ++j) {
      for (i=0; i<nX; ++i) {
        int middleMiddle = 2*i + 2*j*(2*nX - 1) + 2*k*(2*nX-1)*(2*nY-1);
        elemTopo[26] = middleMiddle;
        elemTopo[25] = (i == nX-1) ? -1 : middleMiddle + 1;
        elemTopo[24] = (i == 0) ?    -1 : middleMiddle - 1;
        elemTopo[23] = (j == nY-1) ? -1 : middleMiddle + 2*nX - 1;
        elemTopo[22] = (j == 0) ?    -1 : middleMiddle - 2*nX + 1;
        elemTopo[21] = (k == nZ-1) ? -1 : middleMiddle + (2*nX-1)*(2*nY-1) ;
        elemTopo[20] = (k == 0) ?    -1 : middleMiddle - (2*nX-1)*(2*nY-1);
        elemTopo[19] = ((i == 0) || (j == nY-1)) ? -1 :
                                     elemTopo[23] - 1;
        elemTopo[18] = ((i == nX-1) || (j == nY-1)) ? -1 :
                                     elemTopo[23] + 1;
        elemTopo[17] = ((i == nX-1) || (j == 0)) ? -1 :
                                     elemTopo[22] + 1;
        elemTopo[16] = ((i == 0) || (j == 0)) ? -1 :
                                     elemTopo[22] - 1;
        elemTopo[15] = ((i == 0) || (k == nZ-1)) ? -1 :
                                     elemTopo[21] - 1;
        elemTopo[14] = ((j == nY-1) || (k == nZ-1)) ? -1 :
                                     elemTopo[21] + 2*nX - 1;
        elemTopo[13] = ((i == nX-1) || (k == nZ-1)) ? -1 :
                                     elemTopo[21] + 1;
        elemTopo[12] = ((j == 0) || (k == nZ-1)) ? -1 :
                                     elemTopo[21] - 2*nX + 1;
        elemTopo[11] = ((i == 0) || (k == 0)) ? -1 :
                                     elemTopo[20] - 1;
        elemTopo[10] = ((j == nY-1) || (k == 0)) ? -1 :
                                     elemTopo[20] + 2*nX - 1;
        elemTopo[ 9] = ((i == nX-1) || (k == 0)) ? -1 :
                                     elemTopo[20] + 1;
        elemTopo[ 8] = ((j == 0) || (k == 0)) ? -1 :
                                     elemTopo[20] - 2*nX + 1;
        elemTopo[ 7] = ((i == 0) || (j == nY-1) || (k == nZ-1)) ? -1 :
                                     elemTopo[21] + 2*nX - 2;
        elemTopo[ 6] = ((i == nX-1) || (j == nY-1) || (k == nZ-1)) ? -1 :
                                     elemTopo[21] + 2*nX;
        elemTopo[ 5] = ((i == nX-1) || (j == 0) || (k == nZ-1)) ? -1 :
                                     elemTopo[21] - 2*nX + 2;
        elemTopo[ 4] = ((i == 0) || (j == 0) || (k == nZ-1)) ? -1 :
                                     elemTopo[21] - 2*nX;
        elemTopo[ 3] = ((i == 0) || (j == nY-1) || (k == 0)) ? -1 :
                                     elemTopo[20] + 2*nX - 2;
        elemTopo[ 2] = ((i == nX-1) || (j == nY-1) || (k == 0)) ? -1 :
                                     elemTopo[20] + 2*nX;
        elemTopo[ 1] = ((i == nX-1) || (j == 0) || (k == 0)) ? -1 :
                                     elemTopo[20] - 2*nX + 2;
        elemTopo[ 0] = ((i == 0) || (j == 0) || (k == 0)) ? -1 :
                                     elemTopo[20] - 2*nX;
        for (int iD = 0; iD < dofEle; ++iD) {
          if (elemTopo[iD] == -1)
            continue;
          for (int jD = 0; jD < dofEle; ++jD) {
            int neighbor = elemTopo[jD];
            // Don't count the node itself for Chaco
            if ((neighbor == -1) || (iD == jD))
              continue;
            // Check if this neighbor is already stored
            for (int l = start[elemTopo[iD]]; l < start[elemTopo[iD]+1]; ++l) {
              // Note that Chaco uses a Fortran numbering
              if (adjacency[l] == neighbor + 1)
                break;
              if (adjacency[l] == 0) {
                // Store the node ID
                // Note that Chaco uses a Fortran numbering
                adjacency[l] = elemTopo[jD] + 1;
                break;
              }
            } // for (int l = start[elemTopo[iD]]; l < start[elemTopo[iD]+1]; ++l)
          } // for (int jD = 0; jD < dofEle; ++jD)
        } // for (int iD = 0; iD < dofEle; ++iD)
      }
    }
  }
  delete[] elemTopo;

  int nDir[3];
  nDir[0] = numProc;
  nDir[1] = 0;
  nDir[2] = 0;
  short int *partition = new short int[globalSize];
  // Call Chaco to partition the matrix
  interface(globalSize, start, adjacency, 0, 0, 0, 0, 0, 0, 0,
            partition, 1, 1, nDir, 0, 1, 1, 1, 250, 1, 0.001, 7654321L);
  // Define the Epetra_Map
  int localSize = 0;
  int myPid = MyComm.MyPID();
  for (i = 0; i < globalSize; ++i) {
    if (partition[i] == myPid)
      localSize += 1;
  }
  int *myRows = new int[localSize];
  localSize = 0;
  for (i = 0; i < globalSize; ++i) { 
    if (partition[i] == myPid) {
      myRows[localSize] = i;
      localSize +=1 ;
    }
  }
  delete[] partition;
  delete[] adjacency;
  delete[] start;
  Map = new Epetra_Map(globalSize, localSize, myRows, 0, MyComm);
  delete[] myRows;
#else
  // Create a uniform distribution of the unknowns across processors
  Map = new Epetra_Map(globalSize, 0, MyComm);
#endif

}


int ModeLaplace3DQ2::countElements(bool *isTouched) {

  // This routine counts and flags the elements that contain the nodes
  // on this processor.

  int i, j, k;
  int numEle = 0;

  for (k=0; k<nZ; ++k) {
    for (j=0; j<nY; ++j) {
      for (i=0; i<nX; ++i) {
        isTouched[i + j*nX + k*nX*nY] = false;
        int middleMiddle = 2*i + 2*j*(2*nX - 1) + 2*k*(2*nX-1)*(2*nY-1);
        int node;
        node = middleMiddle;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = (i == nX-1) ? -1 : middleMiddle + 1;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = (i == 0) ?    -1 : middleMiddle - 1;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = (j == nY-1) ? -1 : middleMiddle + 2*nX - 1;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = (j == 0) ?    -1 : middleMiddle - 2*nX + 1;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = (k == nZ-1) ? -1 : middleMiddle + (2*nX-1)*(2*nY-1) ;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = (k == 0) ?    -1 : middleMiddle - (2*nX-1)*(2*nY-1);
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((i == 0) || (j == nY-1)) ? -1 :
                                     middleMiddle + 2*nX - 1 - 1;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((i == nX-1) || (j == nY-1)) ? -1 :
                                     middleMiddle + 2*nX - 1 + 1;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((i == nX-1) || (j == 0)) ? -1 : middleMiddle - 2*nX + 1 + 1;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((i == 0) || (j == 0)) ? -1 : middleMiddle - 2*nX + 1 - 1;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((i == 0) || (k == nZ-1)) ? -1 : middleMiddle + (2*nX-1)*(2*nY-1) - 1;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((j == nY-1) || (k == nZ-1)) ? -1 : middleMiddle + (2*nX-1)*(2*nY-1) + 2*nX-1;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((i == nX-1) || (k == nZ-1)) ? -1 : middleMiddle + (2*nX-1)*(2*nY-1) + 1;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((j == 0) || (k == nZ-1)) ? -1 : middleMiddle + (2*nX-1)*(2*nY-1) - 2*nX+1;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((i == 0) || (k == 0)) ? -1 : middleMiddle - (2*nX-1)*(2*nY-1) - 1;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((j == nY-1) || (k == 0)) ? -1 : middleMiddle - (2*nX-1)*(2*nY-1) + 2*nX-1;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((i == nX-1) || (k == 0)) ? -1 : middleMiddle - (2*nX-1)*(2*nY-1) + 1;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((j == 0) || (k == 0)) ? -1 : middleMiddle - (2*nX-1)*(2*nY-1) - 2*nX+1;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((i == 0) || (j == nY-1) || (k == nZ-1)) ? -1 : 
                                     middleMiddle + (2*nX-1)*(2*nY-1) + 2*nX-2;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((i == nX-1) || (j == nY-1) || (k == nZ-1)) ? -1 :
                                     middleMiddle + (2*nX-1)*(2*nY-1) + 2*nX;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((i == nX-1) || (j == 0) || (k == nZ-1)) ? -1 :
                                     middleMiddle + (2*nX-1)*(2*nY-1) - 2*nX + 2;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((i == 0) || (j == 0) || (k == nZ-1)) ? -1 :
                                     middleMiddle + (2*nX-1)*(2*nY-1) - 2*nX;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((i == 0) || (j == nY-1) || (k == 0)) ? -1 :
                                     middleMiddle - (2*nX-1)*(2*nY-1) + 2*nX - 2;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((i == nX-1) || (j == nY-1) || (k == 0)) ? -1 :
                                     middleMiddle - (2*nX-1)*(2*nY-1) + 2*nX;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((i == nX-1) || (j == 0) || (k == 0)) ? -1 :
                                     middleMiddle - (2*nX-1)*(2*nY-1) - 2*nX + 2;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
        node = ((i == 0) || (j == 0) || (k == 0)) ? -1 :
                                     middleMiddle - (2*nX-1)*(2*nY-1) - 2*nX;
        if ((node > -1) && (Map->LID(node) > -1)) {
          isTouched[i + j*nX + k*nX*nY] = true;
          numEle += 1;
          continue;
        }
      }
    }
  }

  return numEle;

}


void ModeLaplace3DQ2::makeMyElementsTopology(int *elemTopo, bool *isTouched) {

  // Create the element topology for the elements containing nodes for this processor
  // Note: Put the flag -1 when the node has a Dirichlet boundary condition

  int i, j, k;
  int numEle = 0;

  for (k=0; k<nZ; ++k) {
    for (j=0; j<nY; ++j) {
      for (i=0; i<nX; ++i) {
        if (isTouched[i + j*nX + k*nX*nY] == false)
          continue;
        int middleMiddle = 2*i + 2*j*(2*nX - 1) + 2*k*(2*nX-1)*(2*nY-1);
        elemTopo[dofEle*numEle+26] = middleMiddle;
        elemTopo[dofEle*numEle+25] = (i == nX-1) ? -1 : middleMiddle + 1;
        elemTopo[dofEle*numEle+24] = (i == 0) ?    -1 : middleMiddle - 1;
        elemTopo[dofEle*numEle+23] = (j == nY-1) ? -1 : middleMiddle + 2*nX - 1;
        elemTopo[dofEle*numEle+22] = (j == 0) ?    -1 : middleMiddle - 2*nX + 1;
        elemTopo[dofEle*numEle+21] = (k == nZ-1) ? -1 : middleMiddle + (2*nX-1)*(2*nY-1) ;
        elemTopo[dofEle*numEle+20] = (k == 0) ?    -1 : middleMiddle - (2*nX-1)*(2*nY-1);
        elemTopo[dofEle*numEle+19] = ((i == 0) || (j == nY-1)) ? -1 :
                                     elemTopo[dofEle*numEle+23] - 1;
        elemTopo[dofEle*numEle+18] = ((i == nX-1) || (j == nY-1)) ? -1 :
                                     elemTopo[dofEle*numEle+23] + 1;
        elemTopo[dofEle*numEle+17] = ((i == nX-1) || (j == 0)) ? -1 :
                                     elemTopo[dofEle*numEle+22] + 1;
        elemTopo[dofEle*numEle+16] = ((i == 0) || (j == 0)) ? -1 :
                                     elemTopo[dofEle*numEle+22] - 1;
        elemTopo[dofEle*numEle+15] = ((i == 0) || (k == nZ-1)) ? -1 :
                                     elemTopo[dofEle*numEle+21] - 1;
        elemTopo[dofEle*numEle+14] = ((j == nY-1) || (k == nZ-1)) ? -1 :
                                     elemTopo[dofEle*numEle+21] + 2*nX - 1;
        elemTopo[dofEle*numEle+13] = ((i == nX-1) || (k == nZ-1)) ? -1 :
                                     elemTopo[dofEle*numEle+21] + 1;
        elemTopo[dofEle*numEle+12] = ((j == 0) || (k == nZ-1)) ? -1 :
                                     elemTopo[dofEle*numEle+21] - 2*nX + 1;
        elemTopo[dofEle*numEle+11] = ((i == 0) || (k == 0)) ? -1 :
                                     elemTopo[dofEle*numEle+20] - 1;
        elemTopo[dofEle*numEle+10] = ((j == nY-1) || (k == 0)) ? -1 :
                                     elemTopo[dofEle*numEle+20] + 2*nX - 1;
        elemTopo[dofEle*numEle+ 9] = ((i == nX-1) || (k == 0)) ? -1 :
                                     elemTopo[dofEle*numEle+20] + 1;
        elemTopo[dofEle*numEle+ 8] = ((j == 0) || (k == 0)) ? -1 :
                                     elemTopo[dofEle*numEle+20] - 2*nX + 1;
        elemTopo[dofEle*numEle+ 7] = ((i == 0) || (j == nY-1) || (k == nZ-1)) ? -1 :
                                     elemTopo[dofEle*numEle+21] + 2*nX - 2;
        elemTopo[dofEle*numEle+ 6] = ((i == nX-1) || (j == nY-1) || (k == nZ-1)) ? -1 :
                                     elemTopo[dofEle*numEle+21] + 2*nX;
        elemTopo[dofEle*numEle+ 5] = ((i == nX-1) || (j == 0) || (k == nZ-1)) ? -1 :
                                     elemTopo[dofEle*numEle+21] - 2*nX + 2;
        elemTopo[dofEle*numEle+ 4] = ((i == 0) || (j == 0) || (k == nZ-1)) ? -1 :
                                     elemTopo[dofEle*numEle+21] - 2*nX;
        elemTopo[dofEle*numEle+ 3] = ((i == 0) || (j == nY-1) || (k == 0)) ? -1 :
                                     elemTopo[dofEle*numEle+20] + 2*nX - 2;
        elemTopo[dofEle*numEle+ 2] = ((i == nX-1) || (j == nY-1) || (k == 0)) ? -1 :
                                     elemTopo[dofEle*numEle+20] + 2*nX;
        elemTopo[dofEle*numEle+ 1] = ((i == nX-1) || (j == 0) || (k == 0)) ? -1 :
                                     elemTopo[dofEle*numEle+20] - 2*nX + 2;
        elemTopo[dofEle*numEle   ] = ((i == 0) || (j == 0) || (k == 0)) ? -1 :
                                     elemTopo[dofEle*numEle+20] - 2*nX;
        numEle += 1;
      }
    }
  }

}


void ModeLaplace3DQ2::makeMyConnectivity(int *elemTopo, int numEle, int *connectivity,
                                         int *numNz) {

  // This routine creates the connectivity of each managed node
  // from the element topology.

  int i, j;
  int localSize = Map->NumMyElements();

  for (i=0; i<localSize; ++i) {
    numNz[i] = 0;
    for (j=0; j<maxConnect; ++j) {
      connectivity[i*maxConnect + j] = -1; 
    }
  }

  for (i=0; i<numEle; ++i) {
    for (j=0; j<dofEle; ++j) {
      if (elemTopo[dofEle*i+j] == -1)
        continue;
      int node = Map->LID(elemTopo[dofEle*i+j]);
      if (node > -1) {
        int k;
        for (k=0; k<dofEle; ++k) {
          int neighbor = elemTopo[dofEle*i+k];
          if (neighbor == -1)
            continue;
          // Check if this neighbor is already stored
          int l;
          for (l=0; l<maxConnect; ++l) {
            if (neighbor == connectivity[node*maxConnect + l])
              break;
            if (connectivity[node*maxConnect + l] == -1) {
              connectivity[node*maxConnect + l] = neighbor;
              numNz[node] += 1;
              break;
            }
          } // for (l = 0; l < maxConnect; ++l)
        } // for (k = 0; k < dofEle; ++k)
      } // if (node > -1)
    } // for (j = 0; j < dofEle; ++j)
  } // for (i = 0; i < numEle; ++i)

}


void ModeLaplace3DQ2::makeStiffness(int *elemTopo, int numEle, int *connectivity,
                                    int *numNz) {

  // Create Epetra_Matrix for stiffness
  Epetra_CrsMatrix *KK = new Epetra_CrsMatrix(Copy, *Map, numNz);

  int i;
  int localSize = Map->NumMyElements();

  double *values = new double[maxConnect];
  for (i=0; i<maxConnect; ++i) 
    values[i] = 0.0;

  for (i=0; i<localSize; ++i) {
    assert(KK->InsertGlobalValues(Map->GID(i), numNz[i], values, 
           connectivity+maxConnect*i) == 0);
  }

  // Define the elementary matrix
  double *kel = new double[dofEle*dofEle];
  makeElementaryStiffness(kel);

  // Assemble the matrix
  int *indices = new int[dofEle];
  int numEntries;
  int j;
  for (i=0; i<numEle; ++i) {
    for (j=0; j<dofEle; ++j) {
      if (elemTopo[dofEle*i + j] == -1)
        continue;
      if (Map->LID(elemTopo[dofEle*i+j]) == -1)
        continue;
      numEntries = 0;
      int k;
      for (k=0; k<dofEle; ++k) {
        if (elemTopo[dofEle*i+k] == -1)
          continue;
        indices[numEntries] = elemTopo[dofEle*i+k];
        values[numEntries] = kel[dofEle*j + k];
        numEntries += 1;
      }
      assert(KK->SumIntoGlobalValues(elemTopo[dofEle*i+j], numEntries, values, indices) == 0);
    }
  }

  delete[] kel;
  delete[] values;
  delete[] indices;

  assert(KK->FillComplete()== 0);
  assert(KK->OptimizeStorage() == 0);

  K = KK;

}


void ModeLaplace3DQ2::makeElementaryStiffness(double *kel) const {

  int i, j, k;

  double hx = Lx/nX;
  double hy = Ly/nY;
  double hz = Lz/nZ;

  for (i=0; i<dofEle; ++i) {
    for (j=0; j<dofEle; ++j) {
      kel[j + dofEle*i] = 0.0;
    }
  }

  double gaussP[3], gaussW[3];
  gaussP[0] = - sqrt(3.0/5.0); gaussP[1] = 0.0; gaussP[2] = - gaussP[0];
  gaussW[0] = 5.0/9.0; gaussW[1] = 8.0/9.0; gaussW[2] = 5.0/9.0;
  double jac = hx*hy*hz/8.0;
  double *qx = new double[dofEle];
  double *qy = new double[dofEle];
  double *qz = new double[dofEle];
  for (i=0; i<3; ++i) {
    double xi = gaussP[i];
    double wi = gaussW[i];
    for (j=0; j<3; ++j) {
      double eta = gaussP[j];
      double wj = gaussW[j];
      for (k=0; k<3; ++k) {
        double zeta = gaussP[k];
        double wk = gaussW[k];

        // Get the shape functions

        qx[ 0] = 2.0/hx*(xi-0.5)*0.5*eta*(eta-1.0)*0.5*zeta*(zeta-1.0);
        qy[ 0] = 0.5*xi*(xi-1.0)*2.0/hy*(eta-0.5)*0.5*zeta*(zeta-1.0);
        qz[ 0] = 0.5*xi*(xi-1.0)*0.5*eta*(eta-1.0)*2.0/hz*(zeta-0.5);

        qx[ 1] = 2.0/hx*(xi+0.5)*0.5*eta*(eta-1.0)*0.5*zeta*(zeta-1.0);
        qy[ 1] = 0.5*xi*(xi+1.0)*2.0/hy*(eta-0.5)*0.5*zeta*(zeta-1.0);
        qz[ 1] = 0.5*xi*(xi+1.0)*0.5*eta*(eta-1.0)*2.0/hz*(zeta-0.5);

        qx[ 2] = 2.0/hx*(xi+0.5)*0.5*eta*(eta+1.0)*0.5*zeta*(zeta-1.0);
        qy[ 2] = 0.5*xi*(xi+1.0)*2.0/hy*(eta+0.5)*0.5*zeta*(zeta-1.0);
        qz[ 2] = 0.5*xi*(xi+1.0)*0.5*eta*(eta+1.0)*2.0/hz*(zeta-0.5);

        qx[ 3] = 2.0/hx*(xi-0.5)*0.5*eta*(eta+1.0)*0.5*zeta*(zeta-1.0);
        qy[ 3] = 0.5*xi*(xi-1.0)*2.0/hy*(eta+0.5)*0.5*zeta*(zeta-1.0);
        qz[ 3] = 0.5*xi*(xi-1.0)*0.5*eta*(eta+1.0)*2.0/hz*(zeta-0.5);

        qx[ 4] = 2.0/hx*(xi-0.5)*0.5*eta*(eta-1.0)*0.5*zeta*(zeta+1.0);
        qy[ 4] = 0.5*xi*(xi-1.0)*2.0/hy*(eta-0.5)*0.5*zeta*(zeta+1.0);
        qz[ 4] = 0.5*xi*(xi-1.0)*0.5*eta*(eta-1.0)*2.0/hz*(zeta+0.5);

        qx[ 5] = 2.0/hx*(xi+0.5)*0.5*eta*(eta-1.0)*0.5*zeta*(zeta+1.0);
        qy[ 5] = 0.5*xi*(xi+1.0)*2.0/hy*(eta-0.5)*0.5*zeta*(zeta+1.0);
        qz[ 5] = 0.5*xi*(xi+1.0)*0.5*eta*(eta-1.0)*2.0/hz*(zeta+0.5);

        qx[ 6] = 2.0/hx*(xi+0.5)*0.5*eta*(eta+1.0)*0.5*zeta*(zeta+1.0);
        qy[ 6] = 0.5*xi*(xi+1.0)*2.0/hy*(eta+0.5)*0.5*zeta*(zeta+1.0);
        qz[ 6] = 0.5*xi*(xi+1.0)*0.5*eta*(eta+1.0)*2.0/hz*(zeta+0.5);

        qx[ 7] = 2.0/hx*(xi-0.5)*0.5*eta*(eta+1.0)*0.5*zeta*(zeta+1.0);
        qy[ 7] = 0.5*xi*(xi-1.0)*2.0/hy*(eta+0.5)*0.5*zeta*(zeta+1.0);
        qz[ 7] = 0.5*xi*(xi-1.0)*0.5*eta*(eta+1.0)*2.0/hz*(zeta+0.5);

        qx[ 8] = 2.0/hx*(-2.0*xi)*0.5*eta*(eta-1.0)*0.5*zeta*(zeta-1.0);
        qy[ 8] = (xi+1.0)*(1.0-xi)*2.0/hy*(eta-0.5)*0.5*zeta*(zeta-1.0);
        qz[ 8] = (xi+1.0)*(1.0-xi)*0.5*eta*(eta-1.0)*2.0/hz*(zeta-0.5);

        qx[ 9] = 2.0/hx*(xi+0.5)*(eta+1.0)*(1.0-eta)*0.5*zeta*(zeta-1.0);
        qy[ 9] = 0.5*xi*(xi+1.0)*2.0/hy*(-2.0*eta)*0.5*zeta*(zeta-1.0);
        qz[ 9] = 0.5*xi*(xi+1.0)*(eta+1.0)*(1.0-eta)*2.0/hz*(zeta-0.5);

        qx[10] = 2.0/hx*(-2.0*xi)*0.5*eta*(eta+1.0)*0.5*zeta*(zeta-1.0);
        qy[10] = (xi+1.0)*(1.0-xi)*2.0/hy*(eta+0.5)*0.5*zeta*(zeta-1.0);
        qz[10] = (xi+1.0)*(1.0-xi)*0.5*eta*(eta+1.0)*2.0/hz*(zeta-0.5);

        qx[11] = 2.0/hx*(xi-0.5)*(eta+1.0)*(1.0-eta)*0.5*zeta*(zeta-1.0);
        qy[11] = 0.5*xi*(xi-1.0)*2.0/hy*(-2.0*eta)*0.5*zeta*(zeta-1.0);
        qz[11] = 0.5*xi*(xi-1.0)*(eta+1.0)*(1.0-eta)*2.0/hz*(zeta-0.5);

        qx[12] = 2.0/hx*(-2.0*xi)*0.5*eta*(eta-1.0)*0.5*zeta*(zeta+1.0);
        qy[12] = (xi+1.0)*(1.0-xi)*2.0/hy*(eta-0.5)*0.5*zeta*(zeta+1.0);
        qz[12] = (xi+1.0)*(1.0-xi)*0.5*eta*(eta-1.0)*2.0/hz*(zeta+0.5);

        qx[13] = 2.0/hx*(xi+0.5)*(eta+1.0)*(1.0-eta)*0.5*zeta*(zeta+1.0);
        qy[13] = 0.5*xi*(xi+1.0)*2.0/hy*(-2.0*eta)*0.5*zeta*(zeta+1.0);
        qz[13] = 0.5*xi*(xi+1.0)*(eta+1.0)*(1.0-eta)*2.0/hz*(zeta+0.5);

        qx[14] = 2.0/hx*(-2.0*xi)*0.5*eta*(eta+1.0)*0.5*zeta*(zeta+1.0);
        qy[14] = (xi+1.0)*(1.0-xi)*2.0/hy*(eta+0.5)*0.5*zeta*(zeta+1.0);
        qz[14] = (xi+1.0)*(1.0-xi)*0.5*eta*(eta+1.0)*2.0/hz*(zeta+0.5);

        qx[15] = 2.0/hx*(xi-0.5)*(eta+1.0)*(1.0-eta)*0.5*zeta*(zeta+1.0);
        qy[15] = 0.5*xi*(xi-1.0)*2.0/hy*(-2.0*eta)*0.5*zeta*(zeta+1.0);
        qz[15] = 0.5*xi*(xi-1.0)*(eta+1.0)*(1.0-eta)*2.0/hz*(zeta+0.5);

        qx[16] = 2.0/hx*(xi-0.5)*0.5*eta*(eta-1.0)*(zeta+1.0)*(1.0-zeta);
        qy[16] = 0.5*xi*(xi-1.0)*2.0/hy*(eta-0.5)*(zeta+1.0)*(1.0-zeta);
        qz[16] = 0.5*xi*(xi-1.0)*0.5*eta*(eta-1.0)*2.0/hz*(-2.0*zeta);

        qx[17] = 2.0/hx*(xi+0.5)*0.5*eta*(eta-1.0)*(zeta+1.0)*(1.0-zeta);
        qy[17] = 0.5*xi*(xi+1.0)*2.0/hy*(eta-0.5)*(zeta+1.0)*(1.0-zeta);
        qz[17] = 0.5*xi*(xi+1.0)*0.5*eta*(eta-1.0)*2.0/hz*(-2.0*zeta);

        qx[18] = 2.0/hx*(xi+0.5)*0.5*eta*(eta+1.0)*(zeta+1.0)*(1.0-zeta);
        qy[18] = 0.5*xi*(xi+1.0)*2.0/hy*(eta+0.5)*(zeta+1.0)*(1.0-zeta);
        qz[18] = 0.5*xi*(xi+1.0)*0.5*eta*(eta+1.0)*2.0/hz*(-2.0*zeta);

        qx[19] = 2.0/hx*(xi-0.5)*0.5*eta*(eta+1.0)*(zeta+1.0)*(1.0-zeta);
        qy[19] = 0.5*xi*(xi-1.0)*2.0/hy*(eta+0.5)*(zeta+1.0)*(1.0-zeta);
        qz[19] = 0.5*xi*(xi-1.0)*0.5*eta*(eta+1.0)*2.0/hz*(-2.0*zeta);

        qx[20] = 2.0/hx*(-2.0*xi)*(eta+1.0)*(1.0-eta)*0.5*zeta*(zeta-1.0);
        qy[20] = (xi+1.0)*(1.0-xi)*2.0/hy*(-2.0*eta)*0.5*zeta*(zeta-1.0);
        qz[20] = (xi+1.0)*(1.0-xi)*(eta+1.0)*(1.0-eta)*2.0/hz*(zeta-0.5);

        qx[21] = 2.0/hx*(-2.0*xi)*(eta+1.0)*(1.0-eta)*0.5*zeta*(zeta+1.0);
        qy[21] = (xi+1.0)*(1.0-xi)*2.0/hy*(-2.0*eta)*0.5*zeta*(zeta+1.0);
        qz[21] = (xi+1.0)*(1.0-xi)*(eta+1.0)*(1.0-eta)*2.0/hz*(zeta+0.5);

        qx[22] = 2.0/hx*(-2.0*xi)*0.5*eta*(eta-1.0)*(zeta+1.0)*(1.0-zeta);
        qy[22] = (xi+1.0)*(1.0-xi)*2.0/hy*(eta-0.5)*(zeta+1.0)*(1.0-zeta);
        qz[22] = (xi+1.0)*(1.0-xi)*0.5*eta*(eta-1.0)*2.0/hz*(-2.0*zeta);

        qx[23] = 2.0/hx*(-2.0*xi)*0.5*eta*(eta+1.0)*(zeta+1.0)*(1.0-zeta);
        qy[23] = (xi+1.0)*(1.0-xi)*2.0/hy*(eta+0.5)*(zeta+1.0)*(1.0-zeta);
        qz[23] = (xi+1.0)*(1.0-xi)*0.5*eta*(eta+1.0)*2.0/hz*(-2.0*zeta);

        qx[24] = 2.0/hx*(xi-0.5)*(eta+1.0)*(1.0-eta)*(zeta+1.0)*(1.0-zeta);
        qy[24] = 0.5*xi*(xi-1.0)*2.0/hy*(-2.0*eta)*(zeta+1.0)*(1.0-zeta);
        qz[24] = 0.5*xi*(xi-1.0)*(eta+1.0)*(1.0-eta)*2.0/hz*(-2.0*zeta);

        qx[25] = 2.0/hx*(xi+0.5)*(eta+1.0)*(1.0-eta)*(zeta+1.0)*(1.0-zeta);
        qy[25] = 0.5*xi*(xi+1.0)*2.0/hy*(-2.0*eta)*(zeta+1.0)*(1.0-zeta);
        qz[25] = 0.5*xi*(xi+1.0)*(eta+1.0)*(1.0-eta)*2.0/hz*(-2.0*zeta);

        qx[26] = 2.0/hx*(-2.0*xi)*(eta+1.0)*(1.0-eta)*(zeta+1.0)*(1.0-zeta);
        qy[26] = (xi+1.0)*(1.0-xi)*2.0/hy*(-2.0*eta)*(zeta+1.0)*(1.0-zeta);
        qz[26] = (xi+1.0)*(1.0-xi)*(eta+1.0)*(1.0-eta)*2.0/hz*(-2.0*zeta);

        // Add in the elementary matrix
        int ii, jj;
        for (ii=0; ii<dofEle; ++ii) {
          for (jj=ii; jj<dofEle; ++jj) {
            kel[dofEle*ii + jj] += wi*wj*wk*jac*(qx[ii]*qx[jj] + qy[ii]*qy[jj] + qz[ii]*qz[jj]);
            kel[dofEle*jj + ii] = kel[dofEle*ii + jj];
          }
        }
      }
    }
  }
  delete[] qx;
  delete[] qy;
  delete[] qz;

}


void ModeLaplace3DQ2::makeMass(int *elemTopo, int numEle, int *connectivity,
                               int *numNz) {

  // Create Epetra_Matrix for mass
  Epetra_CrsMatrix *MM = new Epetra_CrsMatrix(Copy, *Map, numNz);

  int i;
  int localSize = Map->NumMyElements();

  double *values = new double[maxConnect];
  for (i=0; i<maxConnect; ++i) 
    values[i] = 0.0;
  for (i=0; i<localSize; ++i) 
    assert(MM->InsertGlobalValues(Map->GID(i), numNz[i], values,
                                   connectivity + maxConnect*i) == 0); 

  // Define the elementary matrix
  double *mel = new double[dofEle*dofEle];
  makeElementaryMass(mel);

  // Assemble the matrix
  int *indices = new int[dofEle];
  int numEntries;
  int j;
  for (i=0; i<numEle; ++i) {
    for (j=0; j<dofEle; ++j) {
      if (elemTopo[dofEle*i + j] == -1)
        continue;
      if (Map->LID(elemTopo[dofEle*i+j]) == -1)
        continue;
      numEntries = 0;
      int k;
      for (k=0; k<dofEle; ++k) {
        if (elemTopo[dofEle*i+k] == -1)
          continue;
        indices[numEntries] = elemTopo[dofEle*i+k];
        values[numEntries] = mel[dofEle*j + k];
        numEntries += 1;
      }
      assert(MM->SumIntoGlobalValues(elemTopo[dofEle*i+j], numEntries, values, indices) == 0);
    }
  }

  delete[] mel;
  delete[] values;
  delete[] indices;

  assert(MM->FillComplete()== 0);
  assert(MM->OptimizeStorage() == 0);

  M = MM;

}


void ModeLaplace3DQ2::makeElementaryMass(double *mel) const {

  int i, j, k;

  double hx = Lx/nX;
  double hy = Ly/nY;
  double hz = Lz/nZ;

  for (i=0; i<dofEle; ++i) {
    for (j=0; j<dofEle; ++j) {
      mel[j + dofEle*i] = 0.0;
    }
  }

  double gaussP[3], gaussW[3];
  gaussP[0] = - sqrt(3.0/5.0); gaussP[1] = 0.0; gaussP[2] = - gaussP[0];
  gaussW[0] = 5.0/9.0; gaussW[1] = 8.0/9.0; gaussW[2] = 5.0/9.0;
  double jac = hx*hy*hz/8.0;
  double *q = new double[dofEle];
  for (i=0; i<3; ++i) {
    double xi = gaussP[i];
    double wi = gaussW[i];
    for (j=0; j<3; ++j) {
      double eta = gaussP[j];
      double wj = gaussW[j];
      for (k=0; k<3; ++k) {
        double zeta = gaussP[k];
        double wk = gaussW[k];
        // Get the shape functions
        q[ 0] = 0.5*xi*(xi-1.0)*0.5*eta*(eta-1.0)*0.5*zeta*(zeta-1.0);
        q[ 1] = 0.5*xi*(xi+1.0)*0.5*eta*(eta-1.0)*0.5*zeta*(zeta-1.0);
        q[ 2] = 0.5*xi*(xi+1.0)*0.5*eta*(eta+1.0)*0.5*zeta*(zeta-1.0);
        q[ 3] = 0.5*xi*(xi-1.0)*0.5*eta*(eta+1.0)*0.5*zeta*(zeta-1.0);
        q[ 4] = 0.5*xi*(xi-1.0)*0.5*eta*(eta-1.0)*0.5*zeta*(zeta+1.0);
        q[ 5] = 0.5*xi*(xi+1.0)*0.5*eta*(eta-1.0)*0.5*zeta*(zeta+1.0);
        q[ 6] = 0.5*xi*(xi+1.0)*0.5*eta*(eta+1.0)*0.5*zeta*(zeta+1.0);
        q[ 7] = 0.5*xi*(xi-1.0)*0.5*eta*(eta+1.0)*0.5*zeta*(zeta+1.0);
        q[ 8] = (xi+1.0)*(1.0-xi)*0.5*eta*(eta-1.0)*0.5*zeta*(zeta-1.0);
        q[ 9] = 0.5*xi*(xi+1.0)*(eta+1.0)*(1.0-eta)*0.5*zeta*(zeta-1.0);
        q[10] = (xi+1.0)*(1.0-xi)*0.5*eta*(eta+1.0)*0.5*zeta*(zeta-1.0);
        q[11] = 0.5*xi*(xi-1.0)*(eta+1.0)*(1.0-eta)*0.5*zeta*(zeta-1.0);
        q[12] = (xi+1.0)*(1.0-xi)*0.5*eta*(eta-1.0)*0.5*zeta*(zeta+1.0);
        q[13] = 0.5*xi*(xi+1.0)*(eta+1.0)*(1.0-eta)*0.5*zeta*(zeta+1.0);
        q[14] = (xi+1.0)*(1.0-xi)*0.5*eta*(eta+1.0)*0.5*zeta*(zeta+1.0);
        q[15] = 0.5*xi*(xi-1.0)*(eta+1.0)*(1.0-eta)*0.5*zeta*(zeta+1.0);
        q[16] = 0.5*xi*(xi-1.0)*0.5*eta*(eta-1.0)*(zeta+1.0)*(1.0-zeta);
        q[17] = 0.5*xi*(xi+1.0)*0.5*eta*(eta-1.0)*(zeta+1.0)*(1.0-zeta);
        q[18] = 0.5*xi*(xi+1.0)*0.5*eta*(eta+1.0)*(zeta+1.0)*(1.0-zeta);
        q[19] = 0.5*xi*(xi-1.0)*0.5*eta*(eta+1.0)*(zeta+1.0)*(1.0-zeta);
        q[20] = (xi+1.0)*(1.0-xi)*(eta+1.0)*(1.0-eta)*0.5*zeta*(zeta-1.0);
        q[21] = (xi+1.0)*(1.0-xi)*(eta+1.0)*(1.0-eta)*0.5*zeta*(zeta+1.0);
        q[22] = (xi+1.0)*(1.0-xi)*0.5*eta*(eta-1.0)*(zeta+1.0)*(1.0-zeta);
        q[23] = (xi+1.0)*(1.0-xi)*0.5*eta*(eta+1.0)*(zeta+1.0)*(1.0-zeta);
        q[24] = 0.5*xi*(xi-1.0)*(eta+1.0)*(1.0-eta)*(zeta+1.0)*(1.0-zeta);
        q[25] = 0.5*xi*(xi+1.0)*(eta+1.0)*(1.0-eta)*(zeta+1.0)*(1.0-zeta);
        q[26] = (xi+1.0)*(1.0-xi)*(eta+1.0)*(1.0-eta)*(zeta+1.0)*(1.0-zeta);
        // Add in the elementary matrix
        int ii, jj;
        for (ii=0; ii<dofEle; ++ii) {
          for (jj=ii; jj<dofEle; ++jj) {
            mel[dofEle*ii + jj] += wi*wj*wk*jac*q[ii]*q[jj];
            mel[dofEle*jj + ii] = mel[dofEle*ii + jj];
          }
        }
      }
    }
  }
  delete[] q;

}


double ModeLaplace3DQ2::getFirstMassEigenValue() const {

  double hx = Lx/nX;
  double hy = Ly/nY;
  double hz = Lz/nZ;

  // Compute the coefficient alphaz
  double cosk = cos(M_PI*hz/2/Lz);
  double a = 2.0*cosk;
  double b = 4.0 + cos(M_PI*hz/Lz);
  double c = -2.0*cosk;
  double delta = sqrt(b*b - 4*a*c);
  double alphaz = (-b - delta)*0.5/a;

  // Compute the coefficient alphay
  double cosj = cos(M_PI*hy/2/Ly);
  a = 2.0*cosj;
  b = 4.0 + cos(M_PI*hy/Ly);
  c = -2.0*cosj;
  delta = sqrt(b*b - 4*a*c);
  double alphay = (-b - delta)*0.5/a;

  // Compute the coefficient alphax
  double cosi = cos(M_PI*hx/2/Lx);
  a = 2.0*cosi;
  b = 4.0 + cos(M_PI*hx/Lx);
  c = -2.0*cosi;
  delta = sqrt(b*b - 4*a*c);
  double alphax = (-b - delta)*0.5/a;

  double discrete = hx/15.0*(8.0+2*alphax*cosi);
  discrete *= hy/15.0*(8.0+2*alphay*cosj);
  discrete *= hz/15.0*(8.0+2*alphaz*cosk);

  return discrete;

}


int ModeLaplace3DQ2::eigenCheck(const Epetra_MultiVector &Q, double *lambda, 
                                double *normWeight) const { 

  int info = 0;
  int qc = Q.NumVectors();
  int myPid = MyComm.MyPID();

  cout.precision(2);
  cout.setf(ios::scientific, ios::floatfield);

  // Check orthonormality of eigenvectors
  double tmp = myVerify.errorOrthonormality(&Q, M);
  if (myPid == 0)
    cout << " Maximum coefficient in matrix Q^T M Q - I = " << tmp << endl;

  // Print out norm of residuals
  myVerify.errorEigenResiduals(Q, lambda, K, M, normWeight);

  // Check the eigenvalues
  int numX = (int) ceil(sqrt(Lx*Lx*lambda[qc-1]/M_PI/M_PI));
  numX = (numX > 2*nX) ? 2*nX : numX;
  int numY = (int) ceil(sqrt(Ly*Ly*lambda[qc-1]/M_PI/M_PI));
  numY = (numY > 2*nY) ? 2*nY : numY;
  int numZ = (int) ceil(sqrt(Lz*Lz*lambda[qc-1]/M_PI/M_PI));
  numZ = (numZ > 2*nZ) ? 2*nZ : numZ;
  int newSize = (numX-1)*(numY-1)*(numZ-1);
  double *discrete = new (nothrow) double[2*newSize];
  if (discrete == 0) {
    return -1;
  }
  double *continuous = discrete + newSize;

  double hx = Lx/nX;
  double hy = Ly/nY;
  double hz = Lz/nZ;

  int i, j, k;
  for (k = 1; k < numZ; ++k) {
    // Compute the coefficient alphaz
    double cosk = cos(k*M_PI*hz/2/Lz);
    double a = cosk*(92.0 - 12.0*cos(k*M_PI*hz/Lz));
    double b = 48.0 + 32.0*cos(k*M_PI*hz/Lz);
    double c = -160.0*cosk;
    double delta = sqrt(b*b - 4*a*c);
    double alphaz = ((-b - delta)*0.5/a < 0.0) ? (-b + delta)*0.5/a : (-b - delta)*0.5/a;
    for (j = 1; j < numY; ++j) {
      // Compute the coefficient alphay
      double cosj = cos(j*M_PI*hy/2/Ly);
      a = cosj*(92.0 - 12.0*cos(j*M_PI*hy/Ly));
      b = 48.0 + 32.0*cos(j*M_PI*hy/Ly);
      c = -160.0*cosj;
      delta = sqrt(b*b - 4*a*c);
      double alphay = ((-b - delta)*0.5/a < 0.0) ? (-b + delta)*0.5/a : (-b - delta)*0.5/a;
      for (i = 1; i < numX; ++i) {
        // Compute the coefficient alphax
        double cosi = cos(i*M_PI*hx/2/Lx);
        a = cosi*(92.0 - 12.0*cos(i*M_PI*hx/Lx));
        b = 48.0 + 32.0*cos(i*M_PI*hx/Lx);
        c = -160.0*cosi;
        delta = sqrt(b*b - 4*a*c);
        double alphax = ((-b - delta)*0.5/a < 0.0) ? (-b + delta)*0.5/a : (-b - delta)*0.5/a;
        // Compute the continuous eigenvalue
        int pos = i-1 + (j-1)*(numX-1) + (k-1)*(numX-1)*(numY-1);
        continuous[pos] = M_PI*M_PI*(i*i/(Lx*Lx) + j*j/(Ly*Ly) + k*k/(Lz*Lz));
        // Compute the discrete eigenvalue
        discrete[pos] = 240.0*(1.0-alphax*cosi)/((8.0+2*alphax*cosi)*(3.0*hx*hx));
        discrete[pos] += 240.0*(1.0-alphay*cosj)/((8.0+2*alphay*cosj)*(3.0*hy*hy));
        discrete[pos] += 240.0*(1.0-alphaz*cosk)/((8.0+2*alphaz*cosk)*(3.0*hz*hz));
      }
    }
  }

  // Sort the eigenvalues in ascending order
  mySort.sortScalars(newSize, continuous);

  int *used = new (nothrow) int[newSize];
  if (used == 0) {
    delete[] discrete;
    return -1;
  }

  mySort.sortScalars(newSize, discrete, used);

  int *index = new (nothrow) int[newSize];
  if (index == 0) {
    delete[] discrete;
    delete[] used;
    return -1;
  }

  for (i=0; i<newSize; ++i) {
    index[used[i]] = i;
  }
  delete[] used;

  int nMax = myVerify.errorLambda(continuous, discrete, newSize, lambda, qc);

  // Define the exact discrete eigenvectors
  int localSize = Map->NumMyElements();
  double *vQ = new (nothrow) double[(nMax+1)*localSize + nMax];
  if (vQ == 0) {
    delete[] discrete;
    delete[] index;
    info = -1;
    return info;
  }

  double *normL2 = vQ + (nMax+1)*localSize;
  Epetra_MultiVector Qex(View, *Map, vQ, localSize, nMax);

  if ((myPid == 0) && (nMax > 0)) {
    cout << endl;
    cout << " --- Relative discretization errors for exact eigenvectors ---" << endl;
    cout << endl;
    cout << "       Cont. Values   Disc. Values     Error      H^1 norm   L^2 norm\n";
  }

  for (k=1; k < numZ; ++k) {
    for (j=1; j < numY; ++j) {
      for (i=1; i < numX; ++i) {
        int pos = i-1 + (j-1)*(numX-1) + (k-1)*(numX-1)*(numY-1);
        if (index[pos] < nMax) {
          int ii;
          for (ii=0; ii<localSize; ++ii) {
             // Compute the coefficient alphaz
            double cosk = cos(k*M_PI*hz/2/Lz);
            double a = cosk*(92.0 - 12.0*cos(k*M_PI*hz/Lz));
            double b = 48.0 + 32.0*cos(k*M_PI*hz/Lz);
            double c = -160.0*cosk;
            double delta = sqrt(b*b - 4*a*c);
            double alphaz = ((-b - delta)*0.5/a < 0.0) ? (-b + delta)*0.5/a : (-b - delta)*0.5/a;
            // Compute the coefficient alphay
            double cosj = cos(j*M_PI*hy/2/Ly);
            a = cosj*(92.0 - 12.0*cos(j*M_PI*hy/Ly));
            b = 48.0 + 32.0*cos(j*M_PI*hy/Ly);
            c = -160.0*cosj;
            delta = sqrt(b*b - 4*a*c);
            double alphay = ((-b - delta)*0.5/a < 0.0) ? (-b + delta)*0.5/a : (-b - delta)*0.5/a;
            // Compute the coefficient alphax
            double cosi = cos(i*M_PI*hx/2/Lx);
            a = cosi*(92.0 - 12.0*cos(i*M_PI*hx/Lx));
            b = 48.0 + 32.0*cos(i*M_PI*hx/Lx);
            c = -160.0*cosi;
            delta = sqrt(b*b - 4*a*c);
            double alphax = ((-b - delta)*0.5/a < 0.0) ? (-b + delta)*0.5/a : (-b - delta)*0.5/a;
            // Get the value for this eigenvector
            double coeff = sin(i*(M_PI/Lx)*x[ii])*sin(j*(M_PI/Ly)*y[ii])*sin(k*(M_PI/Lz)*z[ii]);
            if (fabs(x[ii] - floor(x[ii]/hx+0.5)*hx) < 0.25*hx)
              coeff *= alphax;
            if (fabs(y[ii] - floor(y[ii]/hy+0.5)*hy) < 0.25*hy)
              coeff *= alphay;
            if (fabs(z[ii] - floor(z[ii]/hz+0.5)*hz) < 0.25*hz)
              coeff *= alphaz;
            Qex.ReplaceMyValue(ii, index[pos], coeff);
          }
          // Normalize Qex against the mass matrix
          Epetra_MultiVector MQex(View, *Map, vQ + nMax*localSize, localSize, 1);
          Epetra_MultiVector Qi(View, Qex, index[pos], 1);
          M->Apply(Qi, MQex);
          double mnorm = 0.0;
          Qi.Dot(MQex, &mnorm); 
          Qi.Scale(1.0/sqrt(mnorm));
          // Compute the L2 norm
          Epetra_MultiVector shapeInt(View, *Map, vQ + nMax*localSize, localSize, 1);
          for (ii=0; ii<localSize; ++ii) {
            double iX, iY, iZ;
            if (fabs(x[ii] - floor(x[ii]/hx+0.5)*hx) < 0.25*hx)
              iX = 2.0*sin(i*(M_PI/Lx)*x[ii])/(hx*hx*i*(M_PI/Lx)*i*(M_PI/Lx)*i*(M_PI/Lx))*
                   sqrt(2.0/Lx)*( 3*hx*i*(M_PI/Lx) - 4*sin(i*(M_PI/Lx)*hx) +
                                  cos(i*(M_PI/Lx)*hx)*hx*i*(M_PI/Lx) );
            else
              iX = 8.0*sin(i*(M_PI/Lx)*x[ii])/(hx*hx*i*(M_PI/Lx)*i*(M_PI/Lx)*i*(M_PI/Lx))*
                   sqrt(2.0/Lx)*( 2*sin(i*(M_PI/Lx)*0.5*hx) - 
                                  cos(i*(M_PI/Lx)*0.5*hx)*hx*i*(M_PI/Lx));
            if (fabs(y[ii] - floor(y[ii]/hy+0.5)*hy) < 0.25*hy)
              iY = 2.0*sin(j*(M_PI/Ly)*y[ii])/(hy*hy*j*(M_PI/Ly)*j*(M_PI/Ly)*j*(M_PI/Ly))*
                   sqrt(2.0/Ly)*( 3*hy*j*(M_PI/Ly) - 4*sin(j*(M_PI/Ly)*hy) +
                                  cos(j*(M_PI/Ly)*hy)*hy*j*(M_PI/Ly) );
            else
              iY = 8.0*sin(j*(M_PI/Ly)*y[ii])/(hy*hy*j*(M_PI/Ly)*j*(M_PI/Ly)*j*(M_PI/Ly))*
                   sqrt(2.0/Ly)*( 2*sin(j*(M_PI/Ly)*0.5*hy) - 
                                  cos(j*(M_PI/Ly)*0.5*hy)*hy*j*(M_PI/Ly));
            if (fabs(z[ii] - floor(z[ii]/hz+0.5)*hz) < 0.25*hz)
              iZ = 2.0*sin(k*(M_PI/Lz)*z[ii])/(hz*hz*k*(M_PI/Lz)*k*(M_PI/Lz)*k*(M_PI/Lz))*
                   sqrt(2.0/Lz)*( 3*hz*k*(M_PI/Lz) - 4*sin(k*(M_PI/Lz)*hz) +
                                  cos(k*(M_PI/Lz)*hz)*hz*k*(M_PI/Lz) );
            else
              iZ = 8.0*sin(k*(M_PI/Lz)*z[ii])/(hz*hz*k*(M_PI/Lz)*k*(M_PI/Lz)*k*(M_PI/Lz))*
                   sqrt(2.0/Lz)*( 2*sin(k*(M_PI/Lz)*0.5*hz) - 
                                  cos(k*(M_PI/Lz)*0.5*hz)*hz*k*(M_PI/Lz));
            shapeInt.ReplaceMyValue(ii, 0, iX*iY*iZ);
          }
          Qi.Dot(shapeInt, normL2 + index[pos]);
        } // if index[pos] < nMax)
      } // for (i=1; i < numX; ++i)
    } // for (j=1; j < numY; ++j)
  } // for (k=1; k < numZ; ++k)

  if (myPid == 0) {
    for (i = 0; i < nMax; ++i) {
      double normH1 = continuous[i]*(1.0 - 2.0*normL2[i]) + discrete[i];
      normL2[i] = 2.0 - 2.0*normL2[i];
      normH1+= normL2[i];
      // Print out the result
      if (myPid == 0) {
        cout << " ";
        cout.width(4);
        cout << i+1 << ". ";
        cout.setf(ios::scientific, ios::floatfield);
        cout.precision(8);
        cout << continuous[i] << " " << discrete[i] << "  ";
        cout.precision(3);
        cout << fabs(discrete[i] - continuous[i])/continuous[i] << "  ";
        cout << sqrt(fabs(normH1)/(continuous[i]+1.0)) << "  ";
        cout << sqrt(fabs(normL2[i])) << endl;
      }
    } // for (i = 0; i < nMax; ++i)
  } // if (myPid == 0)

  delete[] discrete;
  delete[] index;

  // Check the angles between exact discrete eigenvectors and computed

  myVerify.errorSubspaces(Q, Qex, M);

  delete[] vQ;

  return info;

}


void ModeLaplace3DQ2::memoryInfo() const {

  int myPid = MyComm.MyPID();

  Epetra_RowMatrix *Mat = dynamic_cast<Epetra_RowMatrix *>(M);
  if ((myPid == 0) && (Mat)) {
    cout << " Total number of nonzero entries in mass matrix      = ";
    cout.width(15);
    cout << Mat->NumGlobalNonzeros() << endl;
    double memSize = Mat->NumGlobalNonzeros()*(sizeof(double) + sizeof(int));
    memSize += 2*Mat->NumGlobalRows()*sizeof(int);
    cout << " Memory requested for mass matrix per processor      = (EST) ";
    cout.precision(2);
    cout.width(6);
    cout.setf(ios::fixed, ios::floatfield);
    cout << memSize/1024.0/1024.0/MyComm.NumProc() << " MB " << endl;
    cout << endl;
  }

  Mat = dynamic_cast<Epetra_RowMatrix *>(K);
  if ((myPid == 0) && (Mat)) {
    cout << " Total number of nonzero entries in stiffness matrix = ";
    cout.width(15);
    cout << Mat->NumGlobalNonzeros() << endl;
    double memSize = Mat->NumGlobalNonzeros()*(sizeof(double) + sizeof(int));
    memSize += 2*Mat->NumGlobalRows()*sizeof(int);
    cout << " Memory requested for stiffness matrix per processor = (EST) ";
    cout.precision(2);
    cout.width(6);
    cout.setf(ios::fixed, ios::floatfield);
    cout << memSize/1024.0/1024.0/MyComm.NumProc() << " MB " << endl;
    cout << endl;
  }

}


void ModeLaplace3DQ2::problemInfo() const { 

  int myPid = MyComm.MyPID();

  if (myPid == 0) {
    cout.precision(2);
    cout.setf(ios::fixed, ios::floatfield);
    cout << " --- Problem definition ---\n\n";
    cout << " >> Laplace equation in 3D with homogeneous Dirichlet condition\n";
    cout << " >> Domain = [0, " << Lx << "] x [0, " << Ly << "] x [0, " << Lz << "]\n";
    cout << " >> Orthogonal mesh uniform per direction with Q2 elements (27 nodes)\n";
    cout << endl;
    cout << " Global size = " << Map->NumGlobalElements() << endl;
    cout << endl;
    cout << " Number of elements in [0, " << Lx << "] (X-direction): " << nX << endl;
    cout << " Number of elements in [0, " << Ly << "] (Y-direction): " << nY << endl;
    cout << " Number of elements in [0, " << Lz << "] (Z-direction): " << nZ << endl;
    cout << endl;
    cout << " Number of interior nodes in the X-direction: " << 2*nX-1 << endl;
    cout << " Number of interior nodes in the Y-direction: " << 2*nY-1 << endl;
    cout << " Number of interior nodes in the Z-direction: " << 2*nZ-1 << endl;
    cout << endl;
  }

}


