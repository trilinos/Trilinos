//**************************************************************************
//
//                                 NOTICE
//
// This software is a result of the research described in the report
//
// " A comparison of algorithms for modal analysis in the absence 
//   of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//  Sandia National Laboratories, Technical report SAND2003-1028J.
//
// It is based on the Epetra, AztecOO, and ML packages defined in the Trilinos
// framework ( http://software.sandia.gov/trilinos/ ).
//
// The distribution of this software follows also the rules defined in Trilinos.
// This notice shall be marked on any reproduction of this software, in whole or
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

#include "ModeFileMatrices.h"


ModeFileMatrices::ModeFileMatrices(const Epetra_Comm &_Comm, const char *_stif,
                                   const char *_mass)
                 : MyComm(_Comm),
                 Map(0),
                 K(0),
                 M(0),
                 stiffnessFileName(_stif),
                 massFileName(_mass),
                 kScale(1.0),
                 mScale(1.0),
                 shift(0.0) {

  preProcess();

}


ModeFileMatrices::~ModeFileMatrices() {

  if (Map)
    delete Map;
  Map = 0;

  if (K)
    delete K;
  K = 0;

  if (M)
    delete M;
  M = 0;

}


void ModeFileMatrices::preProcess() {

  int numProc = MyComm.NumProc();
  int myPid = MyComm.MyPID();

  const int defaultSize = 16;

  int *massNNZ = 0;
  int *start = 0;
  int *adjacency = 0;

  int i = 0, j = 0;
  int globalSize = 0;
  double v = 0.0;

  // Read the information from the binary file
  int iProc;
  for (iProc = 0; iProc < numProc; ++iProc) {

    MyComm.Barrier();  

    if (myPid == iProc) {
      // Open the mass input file
      ifstream fileM(massFileName);
      if (!fileM) {
        cerr << "\n\n !! In 'ModeFileMatrices', mass input file opening failed !! \n\n\n";
        exit(1);
      }

      // Get the global size of the problem
      fileM.read((char *) &globalSize, sizeof(int));
      massNNZ = new int[globalSize];
      memset(massNNZ, 0, globalSize*sizeof(int));

      int **tmpAdj = new int*[globalSize];
      for (i = 0; i < globalSize; ++i) {
        tmpAdj[i] = new int[defaultSize];
        memset(tmpAdj[i], 0, defaultSize*sizeof(int)); 
      }

      // Read the adjacency of the mass matrix
      // Assumption: the mass matrix has more zeros than the stiffness
      // Warning: the input matrix files have Fortran-like indices
      i = 0;
      j = 0;
      while (! ((i == globalSize) && (j == globalSize)) ) {
        fileM.read((char *) &i, sizeof(int));
        fileM.read((char *) &j, sizeof(int));
        fileM.read((char *) &v, sizeof(double));
        if (i != j) {
          if (massNNZ[i-1] < defaultSize) {
            tmpAdj[i-1][massNNZ[i-1]] = j;
          }
          else {
            int *tmp = new int[massNNZ[i-1]+1];
            memcpy(tmp, tmpAdj[i-1], massNNZ[i-1]*sizeof(int));
            tmp[massNNZ[i-1]] = j;
            delete[] tmpAdj[i-1];
            tmpAdj[i-1] = tmp;
          }
          massNNZ[i-1] += 1;
          if (massNNZ[j-1] < defaultSize) {
            tmpAdj[j-1][massNNZ[j-1]] = i;
          }
          else {
            int *tmp = new int[massNNZ[j-1]+1];
            memcpy(tmp, tmpAdj[j-1], massNNZ[j-1]*sizeof(int));
            tmp[massNNZ[j-1]] = i;
            delete[] tmpAdj[j-1];
            tmpAdj[j-1] = tmp;
          }
          massNNZ[j-1] += 1;
        } // if (i!= j)
      } // while (! ((i == globalSize) && (j == globalSize)) )
      fileM.close();

      // Open the stiffness file
      ifstream fileK(stiffnessFileName);
      if (!fileK) {
        cerr << "\n\n !! Stiffness input file opening failed !! \n\n\n";
        exit(1);
      }

      fileK.read((char *) &globalSize, sizeof(int));

      start = new int[globalSize+1];
      start[0] = 0;
      memcpy(start + 1, massNNZ, globalSize*sizeof(int));

      // Read the adjacency for the stiffness matrix
      // Warning: the input matrix files have Fortran-like indices
      i = 0;
      j = 0;
      while (! ((i == globalSize) && (j == globalSize)) ) {
        fileK.read((char *) &i, sizeof(int));
        fileK.read((char *) &j, sizeof(int));
        fileK.read((char *) &v, sizeof(double));
        if (i != j) {
          // Check if this entry is already in tmpAdj
          int k;
          for (k = 0; k < massNNZ[i-1]; ++k) {
            if (tmpAdj[i-1][k] == j)
              break;
          }
          if (k == massNNZ[i-1]) {
            // Add a new entry to the global adjacency
            start[i] += 1;
            start[j] += 1;
          }
        } // if (i != j)
      } // while (! ((i == globalSize) && (j == globalSize)) )
      fileK.close();

      // Count the total number of nonzeros
      for (i = 1; i <= globalSize; ++i) 
        start[i] += start[i-1];

      // Put the adjacency of M in the global adjacency
      adjacency = new int[start[globalSize]];
      for (i = 0; i < globalSize; ++i) {
        memcpy(adjacency + start[i], tmpAdj[i], massNNZ[i]*sizeof(int));
        delete[] tmpAdj[i];
      }
      delete[] tmpAdj;

      // Add the adjacency of K to the global adjacency
      fileK.clear();
      fileK.open(stiffnessFileName);

      fileK.read((char *) &globalSize, sizeof(int));

      int *count = new int[globalSize];
      memcpy(count, massNNZ, globalSize*sizeof(int));

      i = 0;
      j = 0;
      while (! ((i == globalSize) && (j == globalSize)) ) {
        fileK.read((char *) &i, sizeof(int));
        fileK.read((char *) &j, sizeof(int));
        fileK.read((char *) &v, sizeof(double));
        if (i != j) {
          // Check if this entry is already stored
          int k;
          int *pointer = adjacency + start[i-1];
          for (k = 0; k < count[i-1]; ++k) {
            if (pointer[k] == j)
              break;
          }
          if (k == count[i-1]) {
            // Add a new entry to the global adjacency
            adjacency[start[i-1]+count[i-1]] = j;
            count[i-1] += 1;
            adjacency[start[j-1]+count[j-1]] = i;
            count[j-1] += 1;
          }
        } // if (i != j)
      } // while (! ((i == globalSize) && (j == globalSize)) )
      fileK.close();

      delete[] count;

    }  // if (myPid == iProc)

    MyComm.Barrier();

  } // for (iProc = 0; iProc < numProc; ++iProc)

  int localSize;
#ifdef _USE_CHACO
  int nDir[3];
  nDir[0] = numProc;
  nDir[1] = 0;
  nDir[2] = 0;
  short int *partition = new short int[globalSize];
  // Call Chaco to partition the matrix
  interface(globalSize, start, adjacency, 0, 0, 0, 0, 0, 0, 0,
            partition, 1, 1, nDir, 0, 1, 1, 1, 250, 1, 0.001, 7654321L);
  // Define the Epetra_Map
  localSize = 0;
  for (i = 0; i < globalSize; ++i) 
    if (partition[i] == myPid)
      localSize += 1;
  int *myRows = new int[localSize];
  localSize = 0;
  for (i = 0; i < globalSize; ++i) { 
    if (partition[i] == myPid) {
      myRows[localSize] = i;
      localSize += 1;
    }
  }
  delete[] partition;
  Map = new Epetra_Map(globalSize, localSize, myRows, 0, MyComm);
  delete[] myRows;
#else
  // Define a uniform distribution of rows
  Map = new Epetra_Map(globalSize, 0, MyComm);
#endif

  localSize = Map->NumMyElements();

  // Convert the adjacency to C numbering
  for (i = 0; i < start[globalSize]; ++i)
    adjacency[i] -= 1;

  // Get the number of non-zero entries for the mass
  int *myMassNNZ = new int[localSize];
  for (i = 0; i < localSize; ++i) {
    int globalID = Map->GID(i);
    // Add one for the diagonal term
    myMassNNZ[i] = massNNZ[globalID] + 1;
  }

  // Define the graph for M as an Epetra_CrsGraph object
  Epetra_CrsMatrix *Mat;
  Mat = new Epetra_CrsMatrix(Copy, *Map, myMassNNZ);
  delete[] myMassNNZ;

  for (i = 0; i < localSize; ++i) {
    int globalID = Map->GID(i);
    // Get the number of non-zeros (add one for the diagonal term)
    int numNonZeros = massNNZ[globalID] + 1;
    int *tmpRow = new int[numNonZeros];
    tmpRow[0] = globalID;
    memcpy(tmpRow + 1, adjacency + start[globalID], (numNonZeros-1)*sizeof(int));
    double *values = new double[numNonZeros];
    memset(values, 0, numNonZeros*sizeof(double));
    assert(Mat->InsertGlobalValues(globalID, numNonZeros, values, tmpRow) == 0);
    delete[] tmpRow;
    delete[] values;
  }
  M = Mat;

  delete[] massNNZ;

  // Get the number of non-zero entries for the stiffness
  // Note that we use the global adjacency as we shift the stiffness
  int *myStifNNZ = new int[localSize];
  for (i = 0; i < localSize; ++i) {
    int globalID = Map->GID(i);
    // Add one for the diagonal term
    myStifNNZ[i] = start[globalID+1] - start[globalID] + 1;
  }

  // Define the graph for K as an Epetra_CrsGraph object
  Mat = new Epetra_CrsMatrix(Copy, *Map, myStifNNZ);
  delete[] myStifNNZ;

  for (i = 0; i < localSize; ++i) {
    int globalID = Map->GID(i);
    // Get the number of non-zeros (add one for the diagonal term)
    int numNonZeros = start[globalID+1] - start[globalID] + 1;
    int *tmpRow = new int[numNonZeros];
    tmpRow[0] = globalID;
    memcpy(tmpRow + 1, adjacency + start[globalID], (numNonZeros-1)*sizeof(int));
    double *values = new double[numNonZeros];
    memset(values, 0, numNonZeros*sizeof(double));
    assert(Mat->InsertGlobalValues(globalID, numNonZeros, values, tmpRow) == 0);
    delete[] tmpRow;
    delete[] values;
  }
  K = Mat;

  delete[] start;
  delete[] adjacency;

  // Fill the entries of the mass and stiffness matrices
  fillMatrices();

}


void ModeFileMatrices::fillMatrices() {

  int globalSize = Map->NumGlobalElements();
  int numProc = MyComm.NumProc();
  int myPid = MyComm.MyPID();

  int i = 0, j = 0;
  int iProc;
  double v = 0.0;

  // Fill the mass matrix
  // Read the information from the mass binary file
  double maxDiagMass = 0.0;
  Epetra_CrsMatrix *MM = dynamic_cast<Epetra_CrsMatrix*>(M);
  for (iProc = 0; iProc < numProc; ++iProc) {
    MyComm.Barrier();  
    if (myPid == iProc) {

      ifstream fileM(massFileName);
      fileM.read((char *) &globalSize, sizeof(int));

      i = 0;
      j = 0;
      while (! ((i == globalSize) && (j == globalSize)) ) {
        fileM.read((char *) &i, sizeof(int));
        fileM.read((char *) &j, sizeof(int));
        fileM.read((char *) &v, sizeof(double));
        if (i == j) {
          maxDiagMass = (fabs(v) > maxDiagMass) ? fabs(v) : maxDiagMass;
        }
        int localI = Map->LID(i-1);
        int localJ = Map->LID(j-1);
        if ((localI < 0) && (localJ < 0))
          continue;
        if (localI > -1) {
          int posJ = j - 1;
          assert(MM->SumIntoGlobalValues(i-1, 1, &v, &posJ) == 0);
        }
        if ((i != j) && (localJ > -1)) {
          int posI = i - 1;
          assert(MM->SumIntoGlobalValues(j-1, 1, &v, &posI) == 0);
        }
      } // while (! ((i == globalSize) && (j == globalSize)) )

    } // if (myPid == iProc)
    MyComm.Barrier();
  } // for (iProc = 0; iProc < numProc; ++iProc)

  // Read the information from the stiffness binary file
  Epetra_CrsMatrix *KK = dynamic_cast<Epetra_CrsMatrix*>(K);
  double maxDiagStif = 0.0;
  for (iProc = 0; iProc < numProc; ++iProc) {
    MyComm.Barrier();
    if (myPid == iProc) {

      ifstream fileK(stiffnessFileName);
      fileK.read((char *) &globalSize, sizeof(int));

      i = 0;
      j = 0;
      while (! ((i == globalSize) && (j == globalSize)) ) {
        fileK.read((char *) &i, sizeof(int));
        fileK.read((char *) &j, sizeof(int));
        fileK.read((char *) &v, sizeof(double));
        if (i == j) {
          maxDiagStif = (fabs(v) > maxDiagStif) ? fabs(v) : maxDiagStif;
        }
        if ((Map->LID(i-1) < 0) && (Map->LID(j-1) < 0))
          continue;
        if (Map->LID(i-1) > -1) {
          int posJ = j - 1;
          assert(KK->SumIntoGlobalValues(i-1, 1, &v, &posJ) == 0);
        }
        if ((i != j) && (Map->LID(j-1) > -1)) {
          int posI = i - 1;
          assert(KK->SumIntoGlobalValues(j-1, 1, &v, &posI) == 0);
        }
      } // while (! ((i == globalSize) && (j == globalSize)) )

    } // if (myPid == iProc)
    MyComm.Barrier();
  } // for (iProc = 0; iProc < numProc; ++iProc)

  // Shift the stiffness matrix with the scaled mass
  //shift = (maxDiagMass == 0.0) ? 0.0 : sqrt(maxDiagStif/maxDiagMass);
  shift = 0.0;
  if (shift != 0.0) {
    int localSize = Map->NumMyElements();
    for (i = 0; i < localSize; ++i) {
      int globalID = Map->GID(i);
      int length = MM->NumGlobalEntries(globalID);
      int numEntries = 0;
      double *values = new double[length];
      int *indices = new int[length];
      assert(MM->ExtractGlobalRowCopy(globalID, length, numEntries, values, indices)==0);
      for (j = 0; j < numEntries; ++j)
        values[j] *= shift;
      assert(KK->SumIntoGlobalValues(globalID, numEntries, values, indices)==0);
      delete[] values;
      delete[] indices;
    }
  }

  // Finish the assembly
  assert(MM->FillComplete()== 0);
  assert(MM->OptimizeStorage() == 0);

  assert(KK->FillComplete()== 0);
  assert(KK->OptimizeStorage() == 0);

  // Scale the stiffness matrix
  kScale = 1.0;
//  kScale = KK->NormInf();
  KK->Scale(kScale);

  // Scale the mass matrix
  mScale = 1.0;
//  mScale = MM->NormInf();
  MM->Scale(mScale);

}


int ModeFileMatrices::eigenCheck(const Epetra_MultiVector &Q, double *lambda,
                                 double *normWeight) const {

  int myPid = MyComm.MyPID();

  // Check orthonormality of eigenvectors
  double maxDot = 0.0;
  int qc = Q.NumVectors();

  int i, j;
  for (i = 0; i < qc; ++i) {
    Epetra_Vector MQi(Copy, Q, i);
    M->Apply(*(Q(i)), MQi);
    double dot = 0.0;
    for (j = 0; j < qc; ++j) {
      Q(j)->Dot(MQi, &dot);
      dot = (i == j) ? fabs(dot - 1.0) : fabs(dot);
      maxDot = (dot > maxDot) ? dot : maxDot;
    }
  }

  if (myPid == 0) {
    cout.precision(2);
    cout.setf(ios::scientific, ios::floatfield);
    cout << " Maximum coefficient in matrix Q^T M Q - I = " << maxDot << endl;
  }

  // Print out norm of residuals
  if (myPid == 0) {
    cout << endl;
    cout << " --- Norms of residuals for computed eigenmodes ---\n";
    cout << endl;
    cout << "        Eigenvalue";
    if (normWeight)
      cout << "     User Norm     Scaled User N.";
    cout << "     2-Norm     Scaled 2-Nor.\n";
  }

  double maxUserNorm = 0.0;
  double minUserNorm = 1.0e+16;
  double maxL2Norm = 0.0;
  double minL2Norm = 1.0e+16;

  // Compute the residuals and norms
  // Check how to incorporate the possible scaling and shift
  for (j=0; j<qc; ++j) {
    Epetra_Vector MQj(Copy, Q, j);
    M->Apply(*(Q(j)), MQj);
    Epetra_Vector KQj(Copy, Q, j);
    K->Apply(*(Q(j)), KQj);

    KQj.Update(-lambda[j], MQj, 1.0);

    double residualL2 = 0.0;
    KQj.Norm2(&residualL2);
    double residualUser = 0.0;
    if (normWeight) {
      Epetra_Vector vectWeight(View, Q.Map(), normWeight);
      KQj.NormWeighted(vectWeight, &residualUser);
    }

    if (myPid == 0) {
      cout << " ";
      cout.width(4);
      cout.precision(8);
      cout.setf(ios::scientific, ios::floatfield);
      cout << j+1 << ". " << lambda[j] << " ";
      if (normWeight)
        cout << residualUser << " " << residualUser/lambda[j] << " ";
      cout << residualL2 << " " << residualL2/lambda[j] << " ";
      cout << endl;
    }
    maxL2Norm = (residualL2/lambda[j] > maxL2Norm) ? residualL2/lambda[j] : maxL2Norm;
    minL2Norm = (residualL2/lambda[j] < minL2Norm) ? residualL2/lambda[j] : minL2Norm;
    if (normWeight) {
      maxUserNorm = (residualUser/lambda[j] > maxUserNorm) ? residualUser/lambda[j]
                                                           : maxUserNorm;
      minUserNorm = (residualUser/lambda[j] < minUserNorm) ? residualUser/lambda[j]
                                                           : minUserNorm;
    }
  } // for j=0; j<qc; ++j) 

  if (myPid == 0) {
    cout << endl;
    if (normWeight) {
      cout << " >> Minimum scaled user-norm of residuals = " << minUserNorm << endl;
      cout << " >> Maximum scaled user-norm of residuals = " << maxUserNorm << endl;
      cout << endl;
    }
    cout << " >> Minimum scaled L2 norm of residuals   = " << minL2Norm << endl;
    cout << " >> Maximum scaled L2 norm of residuals   = " << maxL2Norm << endl;
    cout << endl;
  }

  return 0;

}


void ModeFileMatrices::memoryInfo() const {

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


void ModeFileMatrices::problemInfo() const { 

  int myPid = MyComm.MyPID();

  if (myPid == 0) {
    cout.precision(2);
    cout.setf(ios::fixed, ios::floatfield);
    cout << " --- Problem definition ---\n\n";
    cout << " >> Stiffness matrice from file: " << stiffnessFileName << endl;
    cout << " >> Mass matrix from file      : " << massFileName << endl;
    cout << endl;
    cout << " Global size = " << Map->NumGlobalElements() << endl;
    cout << endl;
  }

}


