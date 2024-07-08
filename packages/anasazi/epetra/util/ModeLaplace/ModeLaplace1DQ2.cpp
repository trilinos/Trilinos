// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This software is a result of the research described in the report
//
//     "A comparison of algorithms for modal analysis in the absence
//     of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//     Sandia National Laboratories, Technical report SAND2003-1028J.
//
// It is based on the Epetra, AztecOO, and ML packages defined in the Trilinos
// framework ( http://trilinos.org/ ).

#include "ModeLaplace1DQ2.h"
#include "Teuchos_Assert.hpp"


const int ModeLaplace1DQ2::dofEle = 3;
const int ModeLaplace1DQ2::maxConnect = 5;
#ifndef M_PI
const double ModeLaplace1DQ2::M_PI = 3.14159265358979323846;
#endif


ModeLaplace1DQ2::ModeLaplace1DQ2(const Epetra_Comm &_Comm, double _Lx, int _nX)
                : myVerify(_Comm),
                  MyComm(_Comm),
                  mySort(),
                  Map(0),
                  K(0),
                  M(0),
                  Lx(_Lx),
                  nX(_nX),
                  x(0)
                {

  preProcess();

}


ModeLaplace1DQ2::~ModeLaplace1DQ2() {

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

}


void ModeLaplace1DQ2::preProcess() {

  // Create the distribution of equations across processors
  makeMap();

  // Count the number of elements touched by this processor
  bool *isTouched = new bool[nX];
  int i;
  for (i=0; i<nX; ++i)
    isTouched[i] = false;

  int numEle = countElements(isTouched);

  // Create the mesh
  int *elemTopo = new int[dofEle*numEle];
  makeMyElementsTopology(elemTopo, isTouched);

  delete[] isTouched;

  // Get the number of nonzeros per row
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

  x = new double[localSize];

  for (i=0; i<2*nX-1; ++i) {
    int node = i;
    if (Map->LID(node) > -1) {
      x[Map->LID(node)] = (i+1)*hx*0.5;
    }
  }

}


void ModeLaplace1DQ2::makeMap() {

  int globalSize = (2*nX - 1);
  TEUCHOS_TEST_FOR_EXCEPTION(globalSize <= MyComm.NumProc(),std::logic_error,"Parameter error in ModeLaplace1DQ2.");

  // Create a uniform distribution of the unknowns across processors
  Map = new Epetra_Map(globalSize, 0, MyComm);

}


int ModeLaplace1DQ2::countElements(bool *isTouched) {

  // This routine counts and flags the elements that contain the nodes
  // on this processor.

  int i;
  int numEle = 0;

  for (i=0; i<nX; ++i) {
    isTouched[i] = false;
    int node;
    node = (i==0)  ? -1 : 2*i-1;
    if ((node > -1) && (Map->LID(node) > -1)) {
      isTouched[i] = true;
      numEle += 1;
      continue;
    }
    node = (i==nX-1) ? -1 : 2*i+1;
    if ((node > -1) && (Map->LID(node) > -1)) {
      isTouched[i] = true;
      numEle += 1;
      continue;
    }
    node = 2*i;
    if ((node > -1) && (Map->LID(node) > -1)) {
      isTouched[i] = true;
      numEle += 1;
      continue;
    }
  }

  return numEle;

}


void ModeLaplace1DQ2::makeMyElementsTopology(int *elemTopo, bool *isTouched) {

  // Create the element topology for the elements containing nodes for this processor
  // Note: Put the flag -1 when the node has a Dirichlet boundary condition

  int i;
  int numEle = 0;

  for (i=0; i<nX; ++i) {
    if (isTouched[i] == false)
      continue;
    elemTopo[dofEle*numEle]   = (i==0)  ? -1 : 2*i-1;
    elemTopo[dofEle*numEle+1] = (i==nX-1) ? -1 : 2*i+1;
    elemTopo[dofEle*numEle+2] = 2*i;
    numEle += 1;
  }

}


void ModeLaplace1DQ2::makeMyConnectivity(int *elemTopo, int numEle, int *connectivity,
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


void ModeLaplace1DQ2::makeStiffness(int *elemTopo, int numEle, int *connectivity,
                                    int *numNz) {

  // Create Epetra_Matrix for stiffness
  K = new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, numNz);

  int i;
  int localSize = Map->NumMyElements();

  double *values = new double[maxConnect];
  for (i=0; i<maxConnect; ++i) 
    values[i] = 0.0;

  for (i=0; i<localSize; ++i) {
    int info = K->InsertGlobalValues(Map->GID(i), numNz[i], values, connectivity+maxConnect*i);
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
        "ModeLaplace1DQ2::makeStiffness(): InsertGlobalValues() returned error code " << info);
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
      int info = K->SumIntoGlobalValues(elemTopo[dofEle*i+j], numEntries, values, indices);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
          "ModeLaplace1DQ2::makeStiffness(): SumIntoGlobalValues() returned error code " << info);
    }
  }

  delete[] kel;
  delete[] values;
  delete[] indices;

  int info = K->FillComplete();
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
      "ModeLaplace1DQ2::makeStiffness(): FillComplete() returned error code " << info);
  info = K->OptimizeStorage();
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
      "ModeLaplace1DQ2::makeStiffness(): OptimizeStorage() returned error code " << info);

}


void ModeLaplace1DQ2::makeElementaryStiffness(double *kel) const {

  double hx = Lx/nX;

  kel[0] = 7.0/(3.0*hx); kel[1] = 1.0/(3.0*hx); kel[2] = - 8.0/(3.0*hx);
  kel[3] = 1.0/(3.0*hx); kel[4] = 7.0/(3.0*hx); kel[5] = - 8.0/(3.0*hx);
  kel[6] = - 8.0/(3.0*hx); kel[7] = - 8.0/(3.0*hx); kel[8] = 16.0/(3.0*hx);

}


void ModeLaplace1DQ2::makeMass(int *elemTopo, int numEle, int *connectivity,
                               int *numNz) {

  // Create Epetra_Matrix for mass
  M = new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *Map, numNz);

  int i;
  int localSize = Map->NumMyElements();

  double *values = new double[maxConnect];
  for (i=0; i<maxConnect; ++i) 
    values[i] = 0.0;
  for (i=0; i<localSize; ++i) {
    int info = M->InsertGlobalValues(Map->GID(i), numNz[i], values, connectivity + maxConnect*i);
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
        "ModeLaplace1DQ2::makeMass(): InsertGlobalValues() returned error code " << info);
  }

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
      int info = M->SumIntoGlobalValues(elemTopo[dofEle*i+j], numEntries, values, indices);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
          "ModeLaplace1DQ2::makeMass(): SumIntoGlobalValues() returned error code " << info);
    }
  }

  delete[] mel;
  delete[] values;
  delete[] indices;

  int info = M->FillComplete();
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
      "ModeLaplace1DQ2::makeMass(): FillComplete() returned error code " << info);
  info = M->OptimizeStorage();
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
      "ModeLaplace1DQ2::makeMass(): OptimizeStorage() returned error code " << info);
}


void ModeLaplace1DQ2::makeElementaryMass(double *mel) const {

  double hx = Lx/nX;

  mel[0] = 2.0*hx/15.0; mel[1] = -hx/30.0; mel[2] = hx/15.0;
  mel[3] = -hx/30.0; mel[4] = 2.0*hx/15.0; mel[5] = hx/15.0;
  mel[6] = hx/15.0; mel[7] = hx/15.0; mel[8] = 8.0*hx/15.0;

}


double ModeLaplace1DQ2::getFirstMassEigenValue() const {

  double hx = Lx/nX;

  // Compute the coefficient alphax
  double cosi = cos(M_PI*hx/2/Lx);
  double a = 2.0*cosi;
  double b = 4.0 + cos(M_PI*hx/Lx);
  double c = -2.0*cosi;
  double delta = sqrt(b*b - 4*a*c);
  double alphax = (-b - delta)*0.5/a;

  double discrete = hx/15.0*(8.0+2*alphax*cosi);

  return discrete;

}


int ModeLaplace1DQ2::eigenCheck(const Epetra_MultiVector &Q, double *lambda, 
                                double *normWeight, bool smallest) const {

  using std::cout;
  using std::ios;

  int info = 0;
  int qc = Q.NumVectors();
  int myPid = MyComm.MyPID();

  std::cout.precision(2);
  std::cout.setf(ios::scientific, ios::floatfield);

  // Check orthonormality of eigenvectors
  double tmp = myVerify.errorOrthonormality(&Q, M);
  if (myPid == 0)
    std::cout << " Maximum coefficient in matrix Q^T M Q - I = " << tmp << std::endl;

  // Print out norm of residuals
  myVerify.errorEigenResiduals(Q, lambda, K, M, normWeight);

  // Check the eigenvalues
  int numX = (int) ceil(sqrt(Lx*Lx*lambda[qc-1]/M_PI/M_PI));
  numX = (numX > 2*nX) ? 2*nX : numX;
  int newSize = (numX-1);
  double *discrete = new (std::nothrow) double[2*newSize];
  if (discrete == 0) {
    return -1;
  }
  double *continuous = discrete + newSize;

  double hx = Lx/nX;

  int i;
  for (i = 1; i < numX; ++i) {
    continuous[i-1] = (M_PI/Lx)*(M_PI/Lx)*i*i;
    // Compute the coefficient alphaX
    double cosi = cos(i*M_PI*hx/2.0/Lx);
    double a = cosi*(92.0 - 12.0*cos(i*M_PI*hx/Lx));
    double b = 48.0 + 32.0*cos(i*M_PI*hx/Lx);
    double c = -160.0*cosi;
    double delta = sqrt(b*b - 4*a*c);
    double alphaX = (-b - delta)*0.5/a;
    alphaX = (alphaX < 0.0) ? (-b + delta)*0.5/a : alphaX;
    // Compute the discrete eigenvalue
    discrete[i-1] = 240.0*(1.0 - alphaX*cosi)/( (8.0 + 2*alphaX*cosi)*(3.0*hx*hx) );
  }

  // Sort the eigenvalues in ascending order
  mySort.sortScalars(newSize, continuous);

  int *used = new (std::nothrow) int[newSize];
  if (used == 0) {
    delete[] discrete;
    return -1;
  }

  mySort.sortScalars(newSize, discrete, used);

  int *index = new (std::nothrow) int[newSize];
  if (index == 0) {
    delete[] discrete;
    delete[] used;
    return -1;
  }

  for (i=0; i<newSize; ++i) {
    index[used[i]] = i;
  }
  delete[] used;

  // sort the eigenvalues/vectors in ascending order
  double *lambdasorted = new (std::nothrow) double[qc];
  if (lambdasorted == 0) {
    delete[] discrete;
    return -1;
  }
  for (i=0; i<qc; i++) {
    lambdasorted[i] = lambda[i];
  }
  mySort.sortScalars(qc, lambdasorted);
  // compare the discrete eigenvalues with the user eigenvalues
  int nMax = myVerify.errorLambda(continuous, discrete, newSize, lambdasorted, qc, smallest);
  // 0 <= nMax < newSize
  // if smallest == true, nMax is the rightmost value in continuous that we matched against
  // if smallest == false, nMax is the leftmost value in continuous that we matched against

  // Define the exact discrete eigenvectors
  int localSize = Map->NumMyElements();
  double *vQ = new (std::nothrow) double[(nMax+2)*localSize + nMax];
  if (vQ == 0) {
    delete[] discrete;
    delete[] index;
    info = -1;
    return info;
  }

  double *normL2 = vQ + (nMax+1)*localSize;
  int Qexdim;
  if (smallest) {
    Qexdim = nMax+1;
  }
  else {
    Qexdim = numX-nMax-1;
  }

  Epetra_MultiVector Qex(Epetra_DataAccess::View, *Map, vQ, localSize, Qexdim);

  if ((myPid == 0) && (Qexdim > 0)) {
    std::cout << std::endl;
    std::cout << " --- Relative discretization errors for exact eigenvectors ---" << std::endl;
    std::cout << std::endl;
    std::cout << "       Cont. Values   Disc. Values     Error      H^1 norm   L^2 norm\n";
  }

  if (smallest) {
    for (i=1; i < numX; ++i) {
      if (index[i-1] <= nMax) {
        // Compute the coefficient alphaX
        double cosi = cos(i*M_PI*hx/2/Lx);
        double a = cosi*(92.0 - 12.0*cos(i*M_PI*hx/Lx));
        double b = 48.0 + 32.0*cos(i*M_PI*hx/Lx);
        double c = -160.0*cosi;
        double delta = sqrt(b*b - 4*a*c);
        double alphaX = (-b - delta)*0.5/a;
        alphaX = (alphaX < 0.0) ? (-b + delta)*0.5/a : alphaX;
        int ii;
        for (ii=0; ii<localSize; ++ii) {
          if (fabs(x[ii] - floor(x[ii]/hx+0.5)*hx) < 0.25*hx) {
            Qex.ReplaceMyValue(ii, index[i-1], alphaX*sin(i*(M_PI/Lx)*x[ii]));
          }
          else {
            Qex.ReplaceMyValue(ii, index[i-1], sin(i*(M_PI/Lx)*x[ii]));
          }
        }
        // Normalize Qex against the mass matrix
        Epetra_MultiVector MQex(Epetra_DataAccess::View, *Map, vQ + (nMax+1)*localSize, localSize, 1);
        Epetra_MultiVector Qi(Epetra_DataAccess::View, Qex, index[i-1], 1);
        M->Apply(Qi, MQex);
        double mnorm = 0.0;
        Qi.Dot(MQex, &mnorm); 
        Qi.Scale(1.0/sqrt(mnorm));
        // Compute the L2 norm
        Epetra_MultiVector shapeInt(Epetra_DataAccess::View, *Map, vQ + (nMax+1)*localSize, localSize, 1);
        for (ii=0; ii<localSize; ++ii) {
          double iX;
          if (fabs(x[ii] - floor(x[ii]/hx+0.5)*hx) < 0.25*hx)
            iX = 2.0*sin(i*(M_PI/Lx)*x[ii])/(hx*hx*i*(M_PI/Lx)*i*(M_PI/Lx)*i*(M_PI/Lx))*
                 sqrt(2.0/Lx)*( 3*hx*i*(M_PI/Lx) - 4*sin(i*(M_PI/Lx)*hx) +
                                cos(i*(M_PI/Lx)*hx)*hx*i*(M_PI/Lx) );
          else
            iX = 8.0*sin(i*(M_PI/Lx)*x[ii])/(hx*hx*i*(M_PI/Lx)*i*(M_PI/Lx)*i*(M_PI/Lx))*
                 sqrt(2.0/Lx)*( 2*sin(i*(M_PI/Lx)*0.5*hx) - 
                                cos(i*(M_PI/Lx)*0.5*hx)*hx*i*(M_PI/Lx));
          shapeInt.ReplaceMyValue(ii, 0, iX);
        }
        Qi.Dot(shapeInt, normL2+index[i-1]);
      } // if index[i-1] <= nMax)
    } // for (i=1; i < numX; ++i)
    if (myPid == 0) {
      for (i = 0; i <= nMax; ++i) {
        double normH1 = continuous[i]*(1.0 - 2.0*normL2[i]) + discrete[i];
        normL2[i] = 2.0 - 2.0*normL2[i];
        normH1+= normL2[i];
        // Print out the result
        if (myPid == 0) {
          std::cout << " ";
          std::cout.width(4);
          std::cout << i+1 << ". ";
          std::cout.setf(ios::scientific, ios::floatfield);
          std::cout.precision(8);
          std::cout << continuous[i] << " " << discrete[i] << "  ";
          std::cout.precision(3);
          std::cout << fabs(discrete[i] - continuous[i])/continuous[i] << "  ";
          std::cout << sqrt(fabs(normH1)/(continuous[i]+1.0)) << "  ";
          std::cout << sqrt(fabs(normL2[i])) << std::endl;
        }
      } // for (i = 0; i <= nMax; ++i)
    } // if (myPid == 0)
  }
  else {
    for (i=1; i < numX; ++i) {
      if (index[i-1] >= nMax) {
        // Compute the coefficient alphaX
        double cosi = cos(i*M_PI*hx/2/Lx);
        double a = cosi*(92.0 - 12.0*cos(i*M_PI*hx/Lx));
        double b = 48.0 + 32.0*cos(i*M_PI*hx/Lx);
        double c = -160.0*cosi;
        double delta = sqrt(b*b - 4*a*c);
        double alphaX = (-b - delta)*0.5/a;
        alphaX = (alphaX < 0.0) ? (-b + delta)*0.5/a : alphaX;
        int ii;
        for (ii=0; ii<localSize; ++ii) {
          if (fabs(x[ii] - floor(x[ii]/hx+0.5)*hx) < 0.25*hx) {
            Qex.ReplaceMyValue(ii, index[i-1]-nMax, alphaX*sin(i*(M_PI/Lx)*x[ii]));
          }
          else {
            Qex.ReplaceMyValue(ii, index[i-1]-nMax, sin(i*(M_PI/Lx)*x[ii]));
          }
        }
        // Normalize Qex against the mass matrix
        Epetra_MultiVector MQex(Epetra_DataAccess::View, *Map, vQ + (numX-nMax-1)*localSize, localSize, 1);
        Epetra_MultiVector Qi(Epetra_DataAccess::View, Qex, index[i-1]-nMax, 1);
        M->Apply(Qi, MQex);
        double mnorm = 0.0;
        Qi.Dot(MQex, &mnorm); 
        Qi.Scale(1.0/sqrt(mnorm));
        // Compute the L2 norm
        Epetra_MultiVector shapeInt(Epetra_DataAccess::View, *Map, vQ + (numX-nMax-1)*localSize, localSize, 1);
        for (ii=0; ii<localSize; ++ii) {
          double iX;
          if (fabs(x[ii] - floor(x[ii]/hx+0.5)*hx) < 0.25*hx)
            iX = 2.0*sin(i*(M_PI/Lx)*x[ii])/(hx*hx*i*(M_PI/Lx)*i*(M_PI/Lx)*i*(M_PI/Lx))*
                 sqrt(2.0/Lx)*( 3*hx*i*(M_PI/Lx) - 4*sin(i*(M_PI/Lx)*hx) +
                                cos(i*(M_PI/Lx)*hx)*hx*i*(M_PI/Lx) );
          else
            iX = 8.0*sin(i*(M_PI/Lx)*x[ii])/(hx*hx*i*(M_PI/Lx)*i*(M_PI/Lx)*i*(M_PI/Lx))*
                 sqrt(2.0/Lx)*( 2*sin(i*(M_PI/Lx)*0.5*hx) - 
                                cos(i*(M_PI/Lx)*0.5*hx)*hx*i*(M_PI/Lx));
          shapeInt.ReplaceMyValue(ii, 0, iX);
        }
        Qi.Dot(shapeInt, normL2+index[i-1]);
      } // if index[i-1] <= nMax)
    } // for (i=1; i < numX; ++i)
    if (myPid == 0) {
      for (i = nMax; i < numX-1; ++i) {
        double normH1 = continuous[i]*(1.0 - 2.0*normL2[i]) + discrete[i];
        normL2[i] = 2.0 - 2.0*normL2[i];
        normH1+= normL2[i];
        // Print out the result
        if (myPid == 0) {
          std::cout << " ";
          std::cout.width(4);
          std::cout << i+1 << ". ";
          std::cout.setf(ios::scientific, ios::floatfield);
          std::cout.precision(8);
          std::cout << continuous[i] << " " << discrete[i] << "  ";
          std::cout.precision(3);
          std::cout << fabs(discrete[i] - continuous[i])/continuous[i] << "  ";
          std::cout << sqrt(fabs(normH1)/(continuous[i]+1.0)) << "  ";
          std::cout << sqrt(fabs(normL2[i])) << std::endl;
        }
      } // for (i = nMax; i < numX-1; ++i)
    } // if (myPid == 0)
  }


  // Check the angles between exact discrete eigenvectors and computed
  myVerify.errorSubspaces(Q, Qex, M);

  delete[] vQ;

  delete[] discrete;
  delete[] index;

  return info;
}


void ModeLaplace1DQ2::memoryInfo() const {
  
  using std::cout;
  using std::ios;

  int myPid = MyComm.MyPID();

  Epetra_RowMatrix *Mat = dynamic_cast<Epetra_RowMatrix *>(M);
  if ((myPid == 0) && (Mat)) {
    std::cout << " Total number of nonzero entries in mass matrix      = ";
    std::cout.width(15);
    std::cout << Mat->NumGlobalNonzeros() << std::endl;
    double memSize = Mat->NumGlobalNonzeros()*(sizeof(double) + sizeof(int));
    memSize += 2*Mat->NumGlobalRows()*sizeof(int);
    std::cout << " Memory requested for mass matrix per processor      = (EST) ";
    std::cout.precision(2);
    std::cout.width(6);
    std::cout.setf(ios::fixed, ios::floatfield);
    std::cout << memSize/1024.0/1024.0/MyComm.NumProc() << " MB " << std::endl;
    std::cout << std::endl;
  }

  Mat = dynamic_cast<Epetra_RowMatrix *>(K);
  if ((myPid == 0) && (Mat)) {
    std::cout << " Total number of nonzero entries in stiffness matrix = ";
    std::cout.width(15);
    std::cout << Mat->NumGlobalNonzeros() << std::endl;
    double memSize = Mat->NumGlobalNonzeros()*(sizeof(double) + sizeof(int));
    memSize += 2*Mat->NumGlobalRows()*sizeof(int);
    std::cout << " Memory requested for stiffness matrix per processor = (EST) ";
    std::cout.precision(2);
    std::cout.width(6);
    std::cout.setf(ios::fixed, ios::floatfield);
    std::cout << memSize/1024.0/1024.0/MyComm.NumProc() << " MB " << std::endl;
    std::cout << std::endl;
  }
}


void ModeLaplace1DQ2::problemInfo() const {

  using std::cout;
  using std::ios;

  int myPid = MyComm.MyPID();

  if (myPid == 0) {
    std::cout.precision(2);
    std::cout.setf(ios::fixed, ios::floatfield);
    std::cout << " --- Problem definition ---\n\n";
    std::cout << " >> Laplace equation in 1D with homogeneous Dirichlet condition\n";
    std::cout << " >> Domain = [0, " << Lx << "]\n";
    std::cout << " >> Orthogonal mesh uniform per direction with Q2 elements (3 nodes)\n";
    std::cout << std::endl;
    std::cout << " Global size = " << Map->NumGlobalElements() << std::endl;
    std::cout << std::endl;
    std::cout << " Number of elements in [0, " << Lx << "] (X-direction): " << nX << std::endl;
    std::cout << std::endl;
    std::cout << " Number of interior nodes in the X-direction: " << 2*nX-1 << std::endl;
    std::cout << std::endl;
  }
}
