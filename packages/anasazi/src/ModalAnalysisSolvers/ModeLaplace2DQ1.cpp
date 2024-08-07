// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// Code Authors: U. Hetmaniuk (ulhetma@sandia.gov), R. Lehoucq (rblehou@sandia.gov)
//
//**************************************************************************

#include "ModeLaplace2DQ1.h"
#include "Teuchos_Assert.hpp"


const int ModeLaplace2DQ1::dofEle = 4;
const int ModeLaplace2DQ1::maxConnect = 9;
#ifndef M_PI
const double ModeLaplace2DQ1::M_PI = 3.14159265358979323846;
#endif


ModeLaplace2DQ1::ModeLaplace2DQ1(const Epetra_Comm &_Comm, double _Lx, int _nX,
                                 double _Ly, int _nY)
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
                  x(0),
                  y(0)
                {

  preProcess();

}


ModeLaplace2DQ1::~ModeLaplace2DQ1() {

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

}


void ModeLaplace2DQ1::preProcess() {

  // Create the distribution of equations across processors
  makeMap();

  // Count the number of elements touched by this processor
  bool *isTouched = new bool[nX*nY];
  int i, j;
  for (j = 0; j < nY; ++j)
    for (i=0; i<nX; ++i)
      isTouched[i + j*nX] = false;

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

  x = new double[localSize];
  y = new double[localSize];

  for (j=0; j<nY-1; ++j) {
    for (i=0; i<nX-1; ++i) {
      int node = i + j*(nX-1);
      if (Map->LID(node) > -1) {
        x[Map->LID(node)] = (i+1)*hx;
        y[Map->LID(node)] = (j+1)*hy;
      }
    }
  }

}


void ModeLaplace2DQ1::makeMap() {

  int numProc = MyComm.NumProc();
  int globalSize = (nX - 1)*(nY - 1);
  assert(globalSize > numProc);

#ifdef _USE_CHACO
  // Use the partitioner Chaco to distribute the unknowns
  int *start = new int[globalSize+1];
  memset(start, 0, (globalSize+1)*sizeof(int));

  int i, j;
  for (j=0; j<nY-1; ++j) {
    for (i=0; i<nX-1; ++i) {
      int node = i + j*(nX-1);
      int connect = 1;
      connect *= ((i == 0) || (i == nX-2)) ? 2 : 3;
      connect *= ((j == 0) || (j == nY-2)) ? 2 : 3;
      // Don't count the node itself for Chaco
      start[node+1] = connect - 1;
    }
  }

  for (i = 0; i < globalSize; ++i)
    start[i+1] += start[i];

  int *adjacency = new int[start[globalSize]];
  memset(adjacency, 0, start[globalSize]*sizeof(int));

  int *elemTopo = new int[dofEle];
  for (j=0; j<nY; ++j) {
    for (i=0; i<nX; ++i) {
      int bottomLeft = i-1 + (j-1)*(nX-1);
      elemTopo[0] = ((i==0) || (j==0))   ? -1 : bottomLeft;
      elemTopo[1] = ((i==nX-1) || (j==0))  ? -1 : bottomLeft + 1;
      elemTopo[2] = ((i==nX-1) || (j==nY-1)) ? -1 : bottomLeft + nX;
      elemTopo[3] = ((i==0) || (j==nY-1))  ? -1 : bottomLeft + nX - 1;
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


int ModeLaplace2DQ1::countElements(bool *isTouched) {

  // This routine counts and flags the elements that contain the nodes
  // on this processor.

  int i, j;
  int numEle = 0;

  for (j=0; j<nY; ++j) {
    for (i=0; i<nX; ++i) {
      isTouched[i+j*nX] = false;
      int bottomLeft = i-1 + (j-1)*(nX-1);
      int node;
      node = ((i==0) || (j==0))   ? -1 : bottomLeft;
      if ((node > -1) && (Map->LID(node) > -1)) {
        isTouched[i+j*nX] = true;
        numEle += 1;
        continue;
      }
      node = ((i==nX-1) || (j==0))  ? -1 : bottomLeft + 1;
      if ((node > -1) && (Map->LID(node) > -1)) {
        isTouched[i+j*nX] = true;
        numEle += 1;
        continue;
      }
      node = ((i==nX-1) || (j==nY-1)) ? -1 : bottomLeft + nX;
      if ((node > -1) && (Map->LID(node) > -1)) {
        isTouched[i+j*nX] = true;
        numEle += 1;
        continue;
      }
      node = ((i==0) || (j==nY-1))  ? -1 : bottomLeft + nX - 1;
      if ((node > -1) && (Map->LID(node) > -1)) {
        isTouched[i+j*nX] = true;
        numEle += 1;
        continue;
      }
    }
  }

  return numEle;

}


void ModeLaplace2DQ1::makeMyElementsTopology(int *elemTopo, bool *isTouched) {

  // Create the element topology for the elements containing nodes for this processor
  // Note: Put the flag -1 when the node has a Dirichlet boundary condition

  int i, j;
  int numEle = 0;

  for (j=0; j<nY; ++j) {
    for (i=0; i<nX; ++i) {
      if (isTouched[i + j*nX] == false)
        continue;
      int bottomLeft = i-1 + (j-1)*(nX-1);
      elemTopo[dofEle*numEle]   = ((i==0) || (j==0))   ? -1 : bottomLeft;
      elemTopo[dofEle*numEle+1] = ((i==nX-1) || (j==0))  ? -1 : bottomLeft + 1;
      elemTopo[dofEle*numEle+2] = ((i==nX-1) || (j==nY-1)) ? -1 : bottomLeft + nX;
      elemTopo[dofEle*numEle+3] = ((i==0) || (j==nY-1))  ? -1 : bottomLeft + nX - 1;
      numEle += 1;
    }
  }

}


void ModeLaplace2DQ1::makeMyConnectivity(int *elemTopo, int numEle, int *connectivity,
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


void ModeLaplace2DQ1::makeStiffness(int *elemTopo, int numEle, int *connectivity,
                                    int *numNz) {

  // Create Epetra_Matrix for stiffness
  K = new Epetra_CrsMatrix(Copy, *Map, numNz);

  int i;
  int localSize = Map->NumMyElements();

  double *values = new double[maxConnect];
  for (i=0; i<maxConnect; ++i) 
    values[i] = 0.0;

  for (i=0; i<localSize; ++i) {
    int info = K->InsertGlobalValues(Map->GID(i), numNz[i], values, connectivity+maxConnect*i);
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
        "ModeLaplace2DQ1::makeMass(): InsertGlobalValues() returned error code " << info);
  }

  // Define the elementary matrix
  double hx = Lx/nX;
  double hy = Ly/nY;

  double *kel = new double[dofEle*dofEle];
  kel[0] = (hx/hy + hy/hx)/3.0; kel[1] = (hx/hy - 2.0*hy/hx)/6.0;
  kel[2] = -(hx/hy + hy/hx)/6.0; kel[3] = (hy/hx - 2.0*hx/hy)/6.0;
  kel[ 4] = kel[1]; kel[ 5] = kel[0]; kel[ 6] = kel[3]; kel[ 7] = kel[2]; 
  kel[ 8] = kel[2]; kel[ 9] = kel[3]; kel[10] = kel[0]; kel[11] = kel[1]; 
  kel[12] = kel[3]; kel[13] = kel[2]; kel[14] = kel[1]; kel[15] = kel[0]; 

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
      int info = K->SumIntoGlobalValues(elemTopo[dofEle*i+j], numEntries, values, indices);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
          "ModeLaplace2DQ1::makeStiffness(): SumIntoGlobalValues() returned error code " << info);
    }
  }

  delete[] kel;
  delete[] values;
  delete[] indices;

  int info = K->FillComplete();
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
      "ModeLaplace2DQ1::makeStiffness(): FillComplete() returned error code " << info);
  info = K->OptimizeStorage();
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
      "ModeLaplace2DQ1::makeStiffness(): OptimizeStorage() returned error code " << info);

}


void ModeLaplace2DQ1::makeMass(int *elemTopo, int numEle, int *connectivity,
                               int *numNz) {

  // Create Epetra_Matrix for mass
  M = new Epetra_CrsMatrix(Copy, *Map, numNz);

  int i;
  int localSize = Map->NumMyElements();

  double *values = new double[maxConnect];
  for (i=0; i<maxConnect; ++i) 
    values[i] = 0.0;
  for (i=0; i<localSize; ++i) {
    int info = M->InsertGlobalValues(Map->GID(i), numNz[i], values, connectivity + maxConnect*i);
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
        "ModeLaplace2DQ1::makeMass(): InsertGlobalValues() returned error code " << info);
  }

  // Define the elementary matrix
  double hx = Lx/nX;
  double hy = Ly/nY;

  double *mel = new double[dofEle*dofEle];
  mel[0] = hx*hy/9.0; mel[1] = hx*hy/18.0; mel[2] = hx*hy/36.0; mel[3] = hx*hy/18.0;
  mel[ 4] = mel[1]; mel[ 5] = mel[0]; mel[ 6] = mel[3]; mel[ 7] = mel[2]; 
  mel[ 8] = mel[2]; mel[ 9] = mel[3]; mel[10] = mel[0]; mel[11] = mel[1]; 
  mel[12] = mel[3]; mel[13] = mel[2]; mel[14] = mel[1]; mel[15] = mel[0]; 

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
      int info = M->SumIntoGlobalValues(elemTopo[dofEle*i+j], numEntries, values, indices);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
          "ModeLaplace2DQ1::makeMass(): SumIntoGlobalValues() returned error code " << info);
    }
  }

  delete[] mel;
  delete[] values;
  delete[] indices;

  int info = M->FillComplete();
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
      "ModeLaplace2DQ1::makeMass(): FillComplete() returned error code " << info);
  info = M->OptimizeStorage();
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
      "ModeLaplace2DQ1::makeMass(): OptimizeStorage() returned error code " << info);

}


double ModeLaplace2DQ1::getFirstMassEigenValue() const {

  return Lx/(3.0*nX)*(2.0-cos(M_PI/nX))*Ly/(3.0*nY)*(2.0-cos(M_PI/nY));

}


int ModeLaplace2DQ1::eigenCheck(const Epetra_MultiVector &Q, double *lambda, 
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
  numX = (numX > nX) ? nX : numX;
  int numY = (int) ceil(sqrt(Ly*Ly*lambda[qc-1]/M_PI/M_PI));
  numY = (numY > nY) ? nY : numY;
  int newSize = (numX-1)*(numY-1);
  double *discrete = new (nothrow) double[2*newSize];
  if (discrete == 0) {
    return -1;
  }
  double *continuous = discrete + newSize;

  double hx = Lx/nX;
  double hy = Ly/nY;

  int i, j;
  for (j = 1; j < numY; ++j) {
    for (i = 1; i < numX; ++i) {
      int pos = i-1 + (j-1)*(numX-1);
      continuous[pos] = M_PI*M_PI*(i*i/(Lx*Lx) + j*j/(Ly*Ly));
      discrete[pos] = 6.0*(1.0-cos(i*(M_PI/Lx)*hx))/(2.0+cos(i*(M_PI/Lx)*hx))/hx/hx;
      discrete[pos] += 6.0*(1.0-cos(j*(M_PI/Ly)*hy))/(2.0+cos(j*(M_PI/Ly)*hy))/hy/hy;
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

  for (j=1; j < numY; ++j) {
    for (i=1; i < numX; ++i) {
      int pos = i-1 + (j-1)*(numX-1);
      if (index[pos] < nMax) {
        double coeff = (2.0 + cos(i*M_PI/Lx*hx))*Lx/6.0;
        coeff *= (2.0 + cos(j*M_PI/Ly*hy))*Ly/6.0;
        coeff = 1.0/sqrt(coeff);
        int ii;
        for (ii=0; ii<localSize; ++ii) {
          Qex.ReplaceMyValue(ii, index[pos], coeff*sin(i*(M_PI/Lx)*x[ii])
                                  *sin(j*(M_PI/Ly)*y[ii]) );
        }
        // Compute the L2 norm
        Epetra_MultiVector shapeInt(View, *Map, vQ + nMax*localSize, localSize, 1);
        Epetra_MultiVector Qi(View, Qex, index[pos], 1);
        for (ii=0; ii<localSize; ++ii) {
          double iX = 4.0*sqrt(2.0/Lx)*sin(i*(M_PI/Lx)*x[ii])/hx*
                      pow(sin(i*(M_PI/Lx)*0.5*hx)/(i*M_PI/Lx), 2.0);
          double iY = 4.0*sqrt(2.0/Ly)*sin(j*(M_PI/Ly)*y[ii])/hy*
                      pow(sin(j*(M_PI/Ly)*0.5*hy)/(j*M_PI/Ly), 2.0);
          shapeInt.ReplaceMyValue(ii, 0, iX*iY);
        }
        normL2[index[pos]] = 0.0;
        Qi.Dot(shapeInt, normL2 + index[pos]);
      } // if index[pos] < nMax)
    } // for (i=1; i < numX; ++i)
  } // for (j=1; j < numY; ++j)

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


void ModeLaplace2DQ1::memoryInfo() const {

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


void ModeLaplace2DQ1::problemInfo() const { 

  int myPid = MyComm.MyPID();

  if (myPid == 0) {
    cout.precision(2);
    cout.setf(ios::fixed, ios::floatfield);
    cout << " --- Problem definition ---\n\n";
    cout << " >> Laplace equation in 2D with homogeneous Dirichlet condition\n";
    cout << " >> Domain = [0, " << Lx << "] x [0, " << Ly << "]\n";
    cout << " >> Orthogonal mesh uniform per direction with Q1 elements\n";
    cout << endl;
    cout << " Global size = " << Map->NumGlobalElements() << endl;
    cout << endl;
    cout << " Number of elements in [0, " << Lx << "] (X-direction): " << nX << endl;
    cout << " Number of elements in [0, " << Ly << "] (Y-direction): " << nY << endl;
    cout << endl;
    cout << " Number of interior nodes in the X-direction: " << nX-1 << endl;
    cout << " Number of interior nodes in the Y-direction: " << nY-1 << endl;
    cout << endl;
  }

}


