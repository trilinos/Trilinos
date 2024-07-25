// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Utils.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_BlockMap.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"

void Ifpack2_PrintLine()
{
  cout << "================================================================================" << endl;
}

//============================================================================
void Ifpack2_BreakForDebugger(Tpetra_Comm& Comm)
{
  char hostname[80];
  char buf[80];
  if (Comm.MyPID()  == 0) cout << "Host and Process Ids for tasks" << endl;
  for (int i = 0; i <Comm.NumProc() ; i++) {
    if (i == Comm.MyPID() ) {
#if defined(TFLOP) || defined(JANUS_STLPORT)
      sprintf(buf, "Host: %s   PID: %d", "janus", getpid());
#elif defined(__INTEL_COMPILER) && defined(_WIN32)
      sprintf(buf,"Intel compiler, unknown hostname and PID!");
#else
      gethostname(hostname, sizeof(hostname));
      sprintf(buf, "Host: %s\tComm.MyPID(): %d\tPID: %d",
              hostname, Comm.MyPID(), getpid());
#endif
      printf("%s\n",buf);
      fflush(stdout);
#if !( defined(__INTEL_COMPILER) && defined(_WIN32) )
      sleep(1);
#endif
    }
  }
  if(Comm.MyPID() == 0) {
    printf("\n");
    printf("** Pausing to attach debugger...\n");
    printf("** You may now attach debugger to the processes listed above.\n");
    printf( "**\n");
    printf( "** Enter a character to continue > "); fflush(stdout);
    char go;
    scanf("%c",&go);
  }

  Comm.Barrier();

}

//============================================================================
Tpetra_CrsMatrix* Ifpack2_CreateOverlappingCrsMatrix(const Tpetra_RowMatrix* Matrix,
                                                    const int OverlappingLevel)
{

  if (OverlappingLevel == 0) 
    return(0); // All done
  if (Matrix->Comm().NumProc() == 1) 
    return(0); // All done

  Tpetra_CrsMatrix* OverlappingMatrix;
  OverlappingMatrix = 0;
  Tpetra_Map* OverlappingMap;
  OverlappingMap = (Tpetra_Map*)&(Matrix->RowMatrixRowMap());

  const Tpetra_RowMatrix* OldMatrix;
  const Tpetra_Map* DomainMap = &(Matrix->OperatorDomainMap());
  const Tpetra_Map* RangeMap = &(Matrix->OperatorRangeMap());

  for (int level = 1; level <= OverlappingLevel ; ++level) {

    if (OverlappingMatrix)
      OldMatrix = OverlappingMatrix;
    else
      OldMatrix = Matrix;

    Tpetra_Import* OverlappingImporter;
    OverlappingImporter = (Tpetra_Import*)OldMatrix->RowMatrixImporter();
    int NumMyElements = OverlappingImporter->TargetMap().NumMyElements();
    int* MyGlobalElements = OverlappingImporter->TargetMap().MyGlobalElements();

    // need to build an Tpetra_Map in this way because Tpetra_CrsMatrix
    // requires Tpetra_Map and not Tpetra_BlockMap

    OverlappingMap = new Tpetra_Map(-1,NumMyElements,MyGlobalElements,
                                    0, Matrix->Comm());

    if (level < OverlappingLevel)
      OverlappingMatrix = new Tpetra_CrsMatrix(Copy, *OverlappingMap, 0);
    else
      // On last iteration, we want to filter out all columns except 
      // those that correspond
      // to rows in the graph.  This assures that our matrix is square
      OverlappingMatrix = new Tpetra_CrsMatrix(Copy, *OverlappingMap, 
                                               *OverlappingMap, 0);

    OverlappingMatrix->Import(*OldMatrix, *OverlappingImporter, Insert);
    if (level < OverlappingLevel) {
      OverlappingMatrix->FillComplete(*DomainMap, *RangeMap);
    }
    else {
      OverlappingMatrix->FillComplete(*DomainMap, *RangeMap);
    }

    delete OverlappingMap;

    if (level > 1) {
      delete OldMatrix;
    }
    OverlappingMatrix->FillComplete();

  }

  return(OverlappingMatrix);
}

//============================================================================
Tpetra_CrsGraph* Ifpack2_CreateOverlappingCrsMatrix(const Tpetra_CrsGraph* Graph,
                                                   const int OverlappingLevel)
{

  if (OverlappingLevel == 0) 
    return(0); // All done
  if (Graph->Comm().NumProc() == 1) 
    return(0); // All done

  Tpetra_CrsGraph* OverlappingGraph;
  Tpetra_BlockMap* OverlappingMap;
  OverlappingGraph = const_cast<Tpetra_CrsGraph*>(Graph);
  OverlappingMap = const_cast<Tpetra_BlockMap*>(&(Graph->RowMap()));

  Tpetra_CrsGraph* OldGraph;
  Tpetra_BlockMap* OldMap;
  const Tpetra_BlockMap* DomainMap = &(Graph->DomainMap());
  const Tpetra_BlockMap* RangeMap = &(Graph->RangeMap());

  for (int level = 1; level <= OverlappingLevel ; ++level) {

    OldGraph = OverlappingGraph;
    OldMap = OverlappingMap;

    Tpetra_Import* OverlappingImporter;
    OverlappingImporter = const_cast<Tpetra_Import*>(OldGraph->Importer());
    OverlappingMap = new Tpetra_BlockMap(OverlappingImporter->TargetMap());

    if (level < OverlappingLevel)
      OverlappingGraph = new Tpetra_CrsGraph(Copy, *OverlappingMap, 0);
    else
      // On last iteration, we want to filter out all columns except 
      // those that correspond
      // to rows in the graph.  This assures that our matrix is square
      OverlappingGraph = new Tpetra_CrsGraph(Copy, *OverlappingMap, 
                                          *OverlappingMap, 0);

    OverlappingGraph->Import(*OldGraph, *OverlappingImporter, Insert);
    if (level < OverlappingLevel) 
      OverlappingGraph->FillComplete(*DomainMap, *RangeMap);
    else {
      // Copy last OverlapImporter because we will use it later
      OverlappingImporter = new Tpetra_Import(*OverlappingMap, *DomainMap);
      OverlappingGraph->FillComplete(*DomainMap, *RangeMap);
    }

    if (level > 1) {
      delete OldGraph;
      delete OldMap;
    }

    delete OverlappingMap;
    OverlappingGraph->FillComplete();

  }

  return(OverlappingGraph);
}

//============================================================================
string Ifpack2_toString(const int& x)
{
  char s[100];
  sprintf(s, "%d", x);
  return string(s);
}

//============================================================================
string Ifpack2_toString(const double& x)
{
  char s[100];
  sprintf(s, "%g", x);
  return string(s);
}

//============================================================================
int Ifpack2_PrintResidual(char* Label, const Tpetra_RowMatrix& A,
                         const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&Y)
{
  if (X.Comm().MyPID() == 0) {
    cout << "***** " << Label << endl;
  }
  Ifpack2_PrintResidual(0,A,X,Y);

  return(0);
}

//============================================================================
int Ifpack2_PrintResidual(const int iter, const Tpetra_RowMatrix& A,
                         const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&Y)
{
  Tpetra_MultiVector RHS(X);
  std::vector<double> Norm2;
  Norm2.resize(X.NumVectors());

  IFPACK2_CHK_ERR(A.Multiply(false,X,RHS));
  RHS.Update(1.0, Y, -1.0);

  RHS.Norm2(&Norm2[0]);

  if (X.Comm().MyPID() == 0) {
    cout << "***** iter: " << iter << ":  ||Ax - b||_2 = " 
         << Norm2[0] << endl;
  }

  return(0);
}

//============================================================================
// very simple debugging function that prints on screen the structure
// of the input matrix.
void Ifpack2_PrintSparsity_Simple(const Tpetra_RowMatrix& A)
{
  int MaxEntries = A.MaxNumEntries();
  vector<int> Indices(MaxEntries);
  vector<double> Values(MaxEntries);
  vector<bool> FullRow(A.NumMyRows());

  cout << "+-";
  for (int j = 0 ; j < A.NumMyRows() ; ++j)
    cout << '-';
  cout << "-+" << endl;

  for (int i = 0 ; i < A.NumMyRows() ; ++i) {

    int Length;
    A.ExtractMyRowCopy(i,MaxEntries,Length,
                       &Values[0], &Indices[0]);

    for (int j = 0 ; j < A.NumMyRows() ; ++j)
      FullRow[j] = false;

    for (int j = 0 ; j < Length ; ++j) {
      FullRow[Indices[j]] = true;
    }

    cout << "| ";
    for (int j = 0 ; j < A.NumMyRows() ; ++j) {
      if (FullRow[j])
        cout << '*';
      else
        cout << ' ';
    }
    cout << " |" << endl;
  }

  cout << "+-";
  for (int j = 0 ; j < A.NumMyRows() ; ++j)
    cout << '-';
  cout << "-+" << endl << endl;

}

//============================================================================

double Ifpack2_FrobeniusNorm(const Tpetra_RowMatrix& A)
{
  double MyNorm = 0.0, GlobalNorm;

  vector<int> colInd(A.MaxNumEntries());
  vector<double> colVal(A.MaxNumEntries());

  for (int i = 0 ; i < A.NumMyRows() ; ++i) {

    int Nnz;
    IFPACK2_CHK_ERR(A.ExtractMyRowCopy(i,A.MaxNumEntries(),Nnz,
                                      &colVal[0],&colInd[0]));

    for (int j = 0 ; j < Nnz ; ++j) {
      MyNorm += colVal[j] * colVal[j];
    }
  }

  A.Comm().SumAll(&MyNorm,&GlobalNorm,1);

  return(sqrt(GlobalNorm));
}

static void print()
{
  printf("\n");
}

#include <iomanip>
template<class T>
static void print(char* str, T val)
{
  cout.width(30); cout.setf(ios::left);
  cout << str;
  cout << " = " << val << endl;
}

template<class T>
static void print(char* str, T val, double percentage)
{
  cout.width(30); cout.setf(ios::left);
  cout << str;
  cout << " = ";
  cout.width(20); cout.setf(ios::left);
  cout << val;
  cout << " ( " << percentage << " %)" << endl;
}
template<class T>
static void print(char* str, T one, T two, T three, bool equal = true)
{
  cout.width(30); cout.setf(ios::left);
  cout << str;
  if (equal) 
    cout << " = ";
  else
    cout << "   ";
  cout.width(15); cout.setf(ios::left);
  cout << one;
  cout.width(15); cout.setf(ios::left);
  cout << two;
  cout.width(15); cout.setf(ios::left);
  cout << three;
  cout << endl;
}

//============================================================================
#include "limits.hpp"
#include "float.hpp"
#include "Tpetra_FECrsMatrix.hpp"

int Ifpack2_Analyze(const Tpetra_RowMatrix& A, const bool Cheap,
                   const int NumPDEEqns)
{

  int NumMyRows = A.NumMyRows();
  int NumGlobalRows = A.NumGlobalRows();
  int NumGlobalCols = A.NumGlobalCols();
  int MyBandwidth = 0, GlobalBandwidth;
  int MyLowerNonzeros = 0, MyUpperNonzeros = 0;
  int GlobalLowerNonzeros, GlobalUpperNonzeros;
  int MyDiagonallyDominant = 0, GlobalDiagonallyDominant;
  int MyWeaklyDiagonallyDominant = 0, GlobalWeaklyDiagonallyDominant;
  double MyMin, MyAvg, MyMax;
  double GlobalMin, GlobalAvg, GlobalMax;
  int GlobalStorage;

  bool verbose = (A.Comm().MyPID() == 0);

  GlobalStorage = sizeof(int*) * NumGlobalRows + 
    sizeof(int) * A.NumGlobalNonzeros() + 
    sizeof(double) * A.NumGlobalNonzeros();

  if (verbose) {
    print();
    Ifpack2_PrintLine();
    print<const char*>("Label", A.Label());
    print<int>("Global rows", NumGlobalRows);
    print<int>("Global columns", NumGlobalCols);
    print<int>("Stored nonzeros", A.NumGlobalNonzeros());
    print<int>("Nonzeros / row", A.NumGlobalNonzeros() / NumGlobalRows);
    print<double>("Estimated storage (Mbytes)", 1.0e-6 * GlobalStorage);
  }

  int NumMyActualNonzeros = 0, NumGlobalActualNonzeros;
  int NumMyEmptyRows = 0, NumGlobalEmptyRows;
  int NumMyDirichletRows = 0, NumGlobalDirichletRows;

  vector<int> colInd(A.MaxNumEntries());
  vector<double> colVal(A.MaxNumEntries());

  Tpetra_Vector Diag(A.RowMatrixRowMap());
  Tpetra_Vector RowSum(A.RowMatrixRowMap());
  Diag.PutScalar(0.0);
  RowSum.PutScalar(0.0);

  for (int i = 0 ; i < NumMyRows ; ++i) {

    int GRID = A.RowMatrixRowMap().GID(i);
    int Nnz;
    IFPACK2_CHK_ERR(A.ExtractMyRowCopy(i,A.MaxNumEntries(),Nnz,
                                      &colVal[0],&colInd[0]));

    if (Nnz == 0)
      NumMyEmptyRows++;

    if (Nnz == 1)
      NumMyDirichletRows++;

    for (int j = 0 ; j < Nnz ; ++j) {

      double v = colVal[j];
      if (v < 0) v = -v;
      if (colVal[j] != 0.0)
        NumMyActualNonzeros++;

      int GCID = A.RowMatrixColMap().GID(colInd[j]);

      if (GCID != GRID)
        RowSum[i] += v;
      else
        Diag[i] = v;

      if (GCID < GRID) 
        MyLowerNonzeros++;
      else if (GCID > GRID) 
        MyUpperNonzeros++;
      int b = GCID - GRID;
      if (b < 0) b = -b;
      if (b > MyBandwidth)
        MyBandwidth = b;
    }

    if (Diag[i] > RowSum[i])
      MyDiagonallyDominant++;

    if (Diag[i] >= RowSum[i])
      MyWeaklyDiagonallyDominant++;

    RowSum[i] += Diag[i];
  }

  // ======================== //
  // summing up global values //
  // ======================== //
 
  A.Comm().SumAll(&MyDiagonallyDominant,&GlobalDiagonallyDominant,1);
  A.Comm().SumAll(&MyWeaklyDiagonallyDominant,&GlobalWeaklyDiagonallyDominant,1);
  A.Comm().SumAll(&NumMyActualNonzeros, &NumGlobalActualNonzeros, 1);
  A.Comm().SumAll(&NumMyEmptyRows, &NumGlobalEmptyRows, 1);
  A.Comm().SumAll(&NumMyDirichletRows, &NumGlobalDirichletRows, 1);
  A.Comm().SumAll(&MyBandwidth, &GlobalBandwidth, 1);
  A.Comm().SumAll(&MyLowerNonzeros, &GlobalLowerNonzeros, 1);
  A.Comm().SumAll(&MyUpperNonzeros, &GlobalUpperNonzeros, 1);
  A.Comm().SumAll(&MyDiagonallyDominant, &GlobalDiagonallyDominant, 1);
  A.Comm().SumAll(&MyWeaklyDiagonallyDominant, &GlobalWeaklyDiagonallyDominant, 1);
 
  double NormOne = A.NormOne();
  double NormInf = A.NormInf();
  double NormF   = Ifpack2_FrobeniusNorm(A);

  if (verbose) {
    print();
    print<int>("Actual nonzeros", NumGlobalActualNonzeros);
    print<int>("Nonzeros in strict lower part", GlobalLowerNonzeros);
    print<int>("Nonzeros in strict upper part", GlobalUpperNonzeros);
    print();
    print<int>("Empty rows", NumGlobalEmptyRows,
               100.0 * NumGlobalEmptyRows / NumGlobalRows);
    print<int>("Dirichlet rows", NumGlobalDirichletRows,
               100.0 * NumGlobalDirichletRows / NumGlobalRows);
    print<int>("Diagonally dominant rows", GlobalDiagonallyDominant,
               100.0 * GlobalDiagonallyDominant / NumGlobalRows);
    print<int>("Weakly diag. dominant rows", 
               GlobalWeaklyDiagonallyDominant,
               100.0 * GlobalWeaklyDiagonallyDominant / NumGlobalRows);
    print();
    print<int>("Maximum bandwidth", GlobalBandwidth);

    print();
    print("", "one-norm", "inf-norm", "Frobenius", false);
    print("", "========", "========", "=========", false);
    print();

    print<double>("A", NormOne, NormInf, NormF);
  }

  if (Cheap == false) {

    // create A + A^T and A - A^T

    Tpetra_FECrsMatrix AplusAT(Copy, A.RowMatrixRowMap(), 0);
    Tpetra_FECrsMatrix AminusAT(Copy, A.RowMatrixRowMap(), 0);

    for (int i = 0 ; i < NumMyRows ; ++i) {

      int GRID = A.RowMatrixRowMap().GID(i);
      assert (GRID != -1);

      int Nnz;
      IFPACK2_CHK_ERR(A.ExtractMyRowCopy(i,A.MaxNumEntries(),Nnz,
                                        &colVal[0],&colInd[0]));

      for (int j = 0 ; j < Nnz ; ++j) {

        int GCID         = A.RowMatrixColMap().GID(colInd[j]);
        assert (GCID != -1);

        double plus_val  = colVal[j];
        double minus_val = -colVal[j];

        if (AplusAT.SumIntoGlobalValues(1,&GRID,1,&GCID,&plus_val) != 0) {
          IFPACK2_CHK_ERR(AplusAT.InsertGlobalValues(1,&GRID,1,&GCID,&plus_val));
        }

        if (AplusAT.SumIntoGlobalValues(1,&GCID,1,&GRID,&plus_val) != 0) {
          IFPACK2_CHK_ERR(AplusAT.InsertGlobalValues(1,&GCID,1,&GRID,&plus_val));
        }

        if (AminusAT.SumIntoGlobalValues(1,&GRID,1,&GCID,&plus_val) != 0) {
          IFPACK2_CHK_ERR(AminusAT.InsertGlobalValues(1,&GRID,1,&GCID,&plus_val));
        }

        if (AminusAT.SumIntoGlobalValues(1,&GCID,1,&GRID,&minus_val) != 0) {
          IFPACK2_CHK_ERR(AminusAT.InsertGlobalValues(1,&GCID,1,&GRID,&minus_val));
        }

      }
    }

    AplusAT.FillComplete();
    AminusAT.FillComplete();

    AplusAT.Scale(0.5);
    AminusAT.Scale(0.5);

    NormOne = AplusAT.NormOne();
    NormInf = AplusAT.NormInf();
    NormF   = Ifpack2_FrobeniusNorm(AplusAT);

    if (verbose) {
      print<double>("A + A^T", NormOne, NormInf, NormF);
    }

    NormOne = AminusAT.NormOne();
    NormInf = AminusAT.NormInf();
    NormF   = Ifpack2_FrobeniusNorm(AminusAT);

    if (verbose) {
      print<double>("A - A^T", NormOne, NormInf, NormF);
    }
  }

  if (verbose) {
    print();
    print<char*>("", "min", "avg", "max", false);
    print<char*>("", "===", "===", "===", false);
  }

  MyMax = -DBL_MAX;
  MyMin = DBL_MAX;
  MyAvg = 0.0;

  for (int i = 0 ; i < NumMyRows ; ++i) {

    int Nnz;
    IFPACK2_CHK_ERR(A.ExtractMyRowCopy(i,A.MaxNumEntries(),Nnz,
                                      &colVal[0],&colInd[0]));

    for (int j = 0 ; j < Nnz ; ++j) {
      MyAvg += colVal[j];
      if (colVal[j] > MyMax) MyMax = colVal[j];
      if (colVal[j] < MyMin) MyMin = colVal[j];
    }
  }

  A.Comm().MaxAll(&MyMax, &GlobalMax, 1);
  A.Comm().MinAll(&MyMin, &GlobalMin, 1);
  A.Comm().SumAll(&MyAvg, &GlobalAvg, 1);
  GlobalAvg /= A.NumGlobalNonzeros();

  if (verbose) {
    print();
    print<double>(" A(i,j)", GlobalMin, GlobalAvg, GlobalMax);
  }

  MyMax = 0.0;
  MyMin = DBL_MAX;
  MyAvg = 0.0;

  for (int i = 0 ; i < NumMyRows ; ++i) {

    int Nnz;
    IFPACK2_CHK_ERR(A.ExtractMyRowCopy(i,A.MaxNumEntries(),Nnz,
                                      &colVal[0],&colInd[0]));

    for (int j = 0 ; j < Nnz ; ++j) {
      double v = colVal[j];
      if (v < 0) v = -v;
      MyAvg += v;
      if (colVal[j] > MyMax) MyMax = v;
      if (colVal[j] < MyMin) MyMin = v;
    }
  }

  A.Comm().MaxAll(&MyMax, &GlobalMax, 1);
  A.Comm().MinAll(&MyMin, &GlobalMin, 1);
  A.Comm().SumAll(&MyAvg, &GlobalAvg, 1);
  GlobalAvg /= A.NumGlobalNonzeros();

  if (verbose) {
    print<double>("|A(i,j)|", GlobalMin, GlobalAvg, GlobalMax);
  }

  // ================= //
  // diagonal elements //
  // ================= //

  Diag.MinValue(&GlobalMin);
  Diag.MaxValue(&GlobalMax);
  Diag.MeanValue(&GlobalAvg);

  if (verbose) {
    print();
    print<double>(" A(k,k)", GlobalMin, GlobalAvg, GlobalMax);
  }

  Diag.Abs(Diag);
  Diag.MinValue(&GlobalMin);
  Diag.MaxValue(&GlobalMax);
  Diag.MeanValue(&GlobalAvg);
  if (verbose) {
    print<double>("|A(k,k)|", GlobalMin, GlobalAvg, GlobalMax);
  }
  
  // ============================================== //
  // cycle over all equations for diagonal elements //
  // ============================================== //

  if (NumPDEEqns > 1 ) {

    if (verbose) print();

    for (int ie = 0 ; ie < NumPDEEqns ; ie++) {

      MyMin = DBL_MAX;
      MyMax = -DBL_MAX;
      MyAvg = 0.0;

      for (int i = ie ; i < Diag.MyLength() ; i += NumPDEEqns) {
        double d = Diag[i];
        MyAvg += d;
        if (d < MyMin)
          MyMin = d;
        if (d > MyMax)
          MyMax = d;
      }
      A.Comm().MinAll(&MyMin, &GlobalMin, 1);
      A.Comm().MaxAll(&MyMax, &GlobalMax, 1);
      A.Comm().SumAll(&MyAvg, &GlobalAvg, 1);
      // does not really work fine if the number of global
      // elements is not a multiple of NumPDEEqns
      GlobalAvg /= (Diag.GlobalLength() / NumPDEEqns);

      if (verbose) {
        char str[80];
        sprintf(str, " A(k,k), eq %d", ie);
        print<double>(str, GlobalMin, GlobalAvg, GlobalMax);
      }
    }
  }

  // ======== //
  // row sums //
  // ======== //
  
  RowSum.MinValue(&GlobalMin);
  RowSum.MaxValue(&GlobalMax);
  RowSum.MeanValue(&GlobalAvg);

  if (verbose) {
    print();
    print<double>(" sum_j A(k,j)", GlobalMin, GlobalAvg, GlobalMax);
  }

  // ===================================== //
  // cycle over all equations for row sums //
  // ===================================== //

  if (NumPDEEqns > 1 ) {

    if (verbose) print();

    for (int ie = 0 ; ie < NumPDEEqns ; ie++) {

      MyMin = DBL_MAX;
      MyMax = -DBL_MAX;
      MyAvg = 0.0;

      for (int i = ie ; i < Diag.MyLength() ; i += NumPDEEqns) {
        double d = RowSum[i];
        MyAvg += d;
        if (d < MyMin)
          MyMin = d;
        if (d > MyMax)
          MyMax = d;
      }
      A.Comm().MinAll(&MyMin, &GlobalMin, 1);
      A.Comm().MaxAll(&MyMax, &GlobalMax, 1);
      A.Comm().SumAll(&MyAvg, &GlobalAvg, 1);
      // does not really work fine if the number of global
      // elements is not a multiple of NumPDEEqns
      GlobalAvg /= (Diag.GlobalLength() / NumPDEEqns);

      if (verbose) {
        char str[80];
        sprintf(str, " sum_j A(k,j), eq %d", ie);
        print<double>(str, GlobalMin, GlobalAvg, GlobalMax);
      }
    }
  }

  if (verbose)
    Ifpack2_PrintLine();

  return(0);
}

int Ifpack2_AnalyzeVectorElements(const Tpetra_Vector& Diagonal,
                                 const bool abs, const int steps)
{

  bool verbose = (Diagonal.Comm().MyPID() == 0);
  double min_val =  DBL_MAX;
  double max_val = -DBL_MAX;

  for (int i = 0 ; i < Diagonal.MyLength() ; ++i) {
    double v = Diagonal[i];
    if (abs)
      if (v < 0) v = -v;
    if (v > max_val)
      max_val = v;
    if (v < min_val)
      min_val = v;
  }

  if (verbose) {
    cout << endl;
    Ifpack2_PrintLine();
    cout << "Vector label = " << Diagonal.Label() << endl;
    cout << endl;
  }

  double delta = (max_val - min_val) / steps;
  for (int k = 0 ; k < steps ; ++k) {

    double below = delta * k + min_val;
    double above = below + delta;
    int MyBelow = 0, GlobalBelow;

    for (int i = 0 ; i < Diagonal.MyLength() ; ++i) {
      double v = Diagonal[i];
      if (v < 0) v = -v;
      if (v >= below && v < above) MyBelow++;
    }

    Diagonal.Comm().SumAll(&MyBelow, &GlobalBelow, 1);

    if (verbose) {
      printf("Elements in [%+7e, %+7e) = %10d ( = %5.2f %%)\n",
             below, above, GlobalBelow,
             100.0 * GlobalBelow / Diagonal.GlobalLength());
    }
  }
  
  if (verbose) {
    Ifpack2_PrintLine();
    cout << endl;
  }

  return(0);
}

// ====================================================================== 

int Ifpack2_AnalyzeMatrixElements(const Tpetra_RowMatrix& A,
                                 const bool abs, const int steps)
{

  bool verbose = (A.Comm().MyPID() == 0);
  double min_val =  DBL_MAX;
  double max_val = -DBL_MAX;

  vector<int>    colInd(A.MaxNumEntries());
  vector<double> colVal(A.MaxNumEntries());
  
  for (int i = 0 ; i < A.NumMyRows() ; ++i) {

    int Nnz;
    IFPACK2_CHK_ERR(A.ExtractMyRowCopy(i,A.MaxNumEntries(),Nnz,
                                      &colVal[0],&colInd[0]));

    for (int j = 0 ; j < Nnz ; ++j) {
      double v = colVal[j];
      if (abs)
        if (v < 0) v = -v;
      if (v < min_val)
        min_val = v;
      if (v > max_val)
        max_val = v;
    }
  }

  if (verbose) {
    cout << endl;
    Ifpack2_PrintLine();
    cout << "Label of matrix = " << A.Label() << endl;
    cout << endl;
  }

  double delta = (max_val - min_val) / steps;
  for (int k = 0 ; k < steps ; ++k) {

    double below = delta * k + min_val;
    double above = below + delta;
    int MyBelow = 0, GlobalBelow;

    for (int i = 0 ; i < A.NumMyRows() ; ++i) {

      int Nnz;
      IFPACK2_CHK_ERR(A.ExtractMyRowCopy(i,A.MaxNumEntries(),Nnz,
                                        &colVal[0],&colInd[0]));

      for (int j = 0 ; j < Nnz ; ++j) {
        double v = colVal[j];
        if (abs)
          if (v < 0) v = -v;
        if (v >= below && v < above) MyBelow++;
      }

    }
    A.Comm().SumAll(&MyBelow, &GlobalBelow, 1);
    if (verbose) {
      printf("Elements in [%+7e, %+7e) = %10d ( = %5.2f %%)\n",
             below, above, GlobalBelow,
             100.0 * GlobalBelow / A.NumGlobalNonzeros());
    }
  }

  if (verbose) {
    Ifpack2_PrintLine();
    cout << endl;
  }

  return(0);
}

// ====================================================================== 
int Ifpack2_PrintSparsity(const Tpetra_RowMatrix& A, const char* InputFileName, 
                         const int NumPDEEqns)
{

  int m,nc,nr,maxdim,ltit;
  double lrmrgn,botmrgn,xtit,ytit,ytitof,fnstit,siz = 0.0;
  double xl,xr, yb,yt, scfct,u2dot,frlw,delt,paperx;
  bool square = false;
  /*change square to .true. if you prefer a square frame around
    a rectangular matrix */
  double conv = 2.54;
  char munt = 'E'; /* put 'E' for centimeters, 'U' for inches */
  int ptitle = 0; /* position of the title, 0 under the drawing,
                     else above */
  FILE* fp = NULL;
  int NumMyRows;
  //int NumMyCols;
  int NumGlobalRows;
  int NumGlobalCols;
  int MyPID;
  int NumProc;
  char FileName[1024];
  char title[1024];
  
  const Tpetra_Comm& Comm = A.Comm();

  /* --------------------- execution begins ---------------------- */

  if (strlen(A.Label()) != 0)
    strcpy(title, A.Label());
  else
    sprintf(title, "%s", "matrix");

  if (InputFileName == 0)
    sprintf(FileName, "%s.ps", title);
  else
    strcpy(FileName, InputFileName);

  MyPID = Comm.MyPID();
  NumProc = Comm.NumProc();

  NumMyRows = A.NumMyRows();
  //NumMyCols = A.NumMyCols();

  NumGlobalRows = A.NumGlobalRows();
  NumGlobalCols = A.NumGlobalCols();

  if (NumGlobalRows != NumGlobalCols)
    IFPACK2_CHK_ERR(-1); // never tested

  /* to be changed for rect matrices */
  maxdim = (NumGlobalRows>NumGlobalCols)?NumGlobalRows:NumGlobalCols;
  maxdim /= NumPDEEqns;

  m = 1 + maxdim;
  nr = NumGlobalRows / NumPDEEqns + 1;
  nc = NumGlobalCols / NumPDEEqns + 1;

  if (munt == 'E') {
    u2dot = 72.0/conv;
    paperx = 21.0;
    siz = 10.0;
  }
  else {
    u2dot = 72.0;
    paperx = 8.5*conv;
    siz = siz*conv;
  }

  /* left and right margins (drawing is centered) */

  lrmrgn = (paperx-siz)/2.0;

  /* bottom margin : 2 cm */

  botmrgn = 2.0;
  /* c scaling factor */
  scfct = siz*u2dot/m;
  /* matrix frame line witdh */
  frlw = 0.25;
  /* font size for title (cm) */
  fnstit = 0.5;
  if (title) ltit = strlen(title);
  else       ltit = 0;
  /* position of title : centered horizontally */
  /*                     at 1.0 cm vertically over the drawing */
  ytitof = 1.0;
  xtit = paperx/2.0;
  ytit = botmrgn+siz*nr/m + ytitof;
  /* almost exact bounding box */
  xl = lrmrgn*u2dot - scfct*frlw/2;
  xr = (lrmrgn+siz)*u2dot + scfct*frlw/2;
  yb = botmrgn*u2dot - scfct*frlw/2;
  yt = (botmrgn+siz*nr/m)*u2dot + scfct*frlw/2;
  if (ltit == 0) {
    yt = yt + (ytitof+fnstit*0.70)*u2dot;
  } 
  /* add some room to bounding box */
  delt = 10.0;
  xl = xl-delt;
  xr = xr+delt;
  yb = yb-delt;
  yt = yt+delt;

  /* correction for title under the drawing */
  if ((ptitle == 0) && (ltit == 0)) {
    ytit = botmrgn + fnstit*0.3;
    botmrgn = botmrgn + ytitof + fnstit*0.7;
  }

  /* begin of output */

  if (MyPID == 0) {

    fp = fopen(FileName,"w");

    fprintf(fp,"%%%%!PS-Adobe-2.0\n");
    fprintf(fp,"%%%%Creator: TIFPACK\n");
    fprintf(fp,"%%%%BoundingBox: %f %f %f %f\n",
            xl,yb,xr,yt);
    fprintf(fp,"%%%%EndComments\n");
    fprintf(fp,"/cm {72 mul 2.54 div} def\n");
    fprintf(fp,"/mc {72 div 2.54 mul} def\n");
    fprintf(fp,"/pnum { 72 div 2.54 mul 20 string ");
    fprintf(fp,"cvs print ( ) print} def\n");
    fprintf(fp,"/Cshow {dup stringwidth pop -2 div 0 rmoveto show} def\n");

    /* we leave margins etc. in cm so it is easy to modify them if
       needed by editing the output file */
    fprintf(fp,"gsave\n");
    if (ltit != 0) {
      fprintf(fp,"/Helvetica findfont %e cm scalefont setfont\n",
              fnstit);
      fprintf(fp,"%f cm %f cm moveto\n",
              xtit,ytit);
      fprintf(fp,"(%s) Cshow\n", title);
      fprintf(fp,"%f cm %f cm translate\n",
              lrmrgn,botmrgn);
    }
    fprintf(fp,"%f cm %d div dup scale \n",
            siz,m);
    /* draw a frame around the matrix */

    fprintf(fp,"%f setlinewidth\n",
            frlw);
    fprintf(fp,"newpath\n");
    fprintf(fp,"0 0 moveto ");
    if (square) {
      printf("------------------- %d\n",m);
      fprintf(fp,"%d %d lineto\n",
              m,0);
      fprintf(fp,"%d %d lineto\n",
              m, m);
      fprintf(fp,"%d %d lineto\n",
              0, m);
    } 
    else {
      fprintf(fp,"%d %d lineto\n",
              nc, 0);
      fprintf(fp,"%d %d lineto\n",
              nc, nr);
      fprintf(fp,"%d %d lineto\n",
              0, nr);
    }
    fprintf(fp,"closepath stroke\n");

    /* plotting loop */

    fprintf(fp,"1 1 translate\n");
    fprintf(fp,"0.8 setlinewidth\n");
    fprintf(fp,"/p {moveto 0 -.40 rmoveto \n");
    fprintf(fp,"           0  .80 rlineto stroke} def\n");

    fclose(fp);
  }
  
  int MaxEntries = A.MaxNumEntries();
  vector<int> Indices(MaxEntries);
  vector<double> Values(MaxEntries);

  for (int pid = 0 ; pid < NumProc ; ++pid) {

    if (pid == MyPID) {

      fp = fopen(FileName,"a");
      if( fp == NULL ) {
        fprintf(stderr,"ERROR\n");
        exit(EXIT_FAILURE);
      }

      for (int i = 0 ; i < NumMyRows ; ++i) {

        if (i % NumPDEEqns) continue;

        int Nnz;
        A.ExtractMyRowCopy(i,MaxEntries,Nnz,&Values[0],&Indices[0]);

        int grow = A.RowMatrixRowMap().GID(i);

        for (int j = 0 ; j < Nnz ; ++j) {
          int col = Indices[j];
          if (col % NumPDEEqns == 0) {
            int gcol = A.RowMatrixColMap().GID(Indices[j]);
            grow /= NumPDEEqns;
            gcol /= NumPDEEqns;
            fprintf(fp,"%d %d p\n",
                    gcol, NumGlobalRows - grow - 1); 
          }
        }
      }

      fprintf(fp,"%%end of data for this process\n");

      if( pid == NumProc - 1 )
        fprintf(fp,"showpage\n");

      fclose(fp);
    }
    Comm.Barrier();
  }

  return(0);
}
