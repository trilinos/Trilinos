#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Import.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

//============================================================================
void Ifpack_BreakForDebugger(Epetra_Comm& Comm)
{
  char hostname[80];
  char buf[80];
  if (Comm.MyPID()  == 0) cout << "Host and Process Ids for tasks" << endl;
  for (int i = 0; i <Comm.NumProc() ; i++) {
    if (i == Comm.MyPID() ) {
      gethostname(hostname, sizeof(hostname));
      sprintf(buf, "Host: %s\tComm.MyPID(): %d\tPID: %d",
              hostname, Comm.MyPID(), getpid());
      printf("%s\n",buf);
      fflush(stdout);
      sleep(1);
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
Epetra_CrsMatrix* Ifpack_CreateOverlappingCrsMatrix(const Epetra_RowMatrix* Matrix,
                                                    const int OverlappingLevel)
{

  if (OverlappingLevel == 0) 
    return(0); // All done
  if (Matrix->Comm().NumProc() == 1) 
    return(0); // All done

  Epetra_CrsMatrix* OverlappingMatrix;
  OverlappingMatrix = 0;
  Epetra_Map* OverlappingMap;
  OverlappingMap = (Epetra_Map*)&(Matrix->RowMatrixRowMap());

  const Epetra_RowMatrix* OldMatrix;
  const Epetra_Map* DomainMap = &(Matrix->OperatorDomainMap());
  const Epetra_Map* RangeMap = &(Matrix->OperatorRangeMap());

  for (int level = 1; level <= OverlappingLevel ; ++level) {

    if (OverlappingMatrix)
      OldMatrix = OverlappingMatrix;
    else
      OldMatrix = Matrix;

    Epetra_Import* OverlappingImporter;
    OverlappingImporter = (Epetra_Import*)OldMatrix->RowMatrixImporter();
    int NumMyElements = OverlappingImporter->TargetMap().NumMyElements();
    int* MyGlobalElements = OverlappingImporter->TargetMap().MyGlobalElements();

    // need to build an Epetra_Map in this way because Epetra_CrsMatrix
    // requires Epetra_Map and not Epetra_BlockMap

    OverlappingMap = new Epetra_Map(-1,NumMyElements,MyGlobalElements,
                                    0, Matrix->Comm());

    if (level < OverlappingLevel)
      OverlappingMatrix = new Epetra_CrsMatrix(Copy, *OverlappingMap, 0);
    else
      // On last iteration, we want to filter out all columns except 
      // those that correspond
      // to rows in the graph.  This assures that our matrix is square
      OverlappingMatrix = new Epetra_CrsMatrix(Copy, *OverlappingMap, 
                                               *OverlappingMap, 0);

    OverlappingMatrix->Import(*OldMatrix, *OverlappingImporter, Insert);
    if (level < OverlappingLevel) {
      OverlappingMatrix->TransformToLocal(DomainMap, RangeMap);
    }
    else {
      OverlappingMatrix->TransformToLocal(DomainMap, RangeMap);
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
Epetra_CrsGraph* Ifpack_CreateOverlappingCrsMatrix(const Epetra_CrsGraph* Graph,
                                                   const int OverlappingLevel)
{

  if (OverlappingLevel == 0) 
    return(0); // All done
  if (Graph->Comm().NumProc() == 1) 
    return(0); // All done

  Epetra_CrsGraph* OverlappingGraph;
  Epetra_BlockMap* OverlappingMap;
  OverlappingGraph = const_cast<Epetra_CrsGraph*>(Graph);
  OverlappingMap = const_cast<Epetra_BlockMap*>(&(Graph->RowMap()));

  Epetra_CrsGraph* OldGraph;
  Epetra_BlockMap* OldMap;
  const Epetra_BlockMap* DomainMap = &(Graph->DomainMap());
  const Epetra_BlockMap* RangeMap = &(Graph->RangeMap());

  for (int level = 1; level <= OverlappingLevel ; ++level) {

    OldGraph = OverlappingGraph;
    OldMap = OverlappingMap;

    Epetra_Import* OverlappingImporter;
    OverlappingImporter = const_cast<Epetra_Import*>(OldGraph->Importer());
    OverlappingMap = new Epetra_BlockMap(OverlappingImporter->TargetMap());

    if (level < OverlappingLevel)
      OverlappingGraph = new Epetra_CrsGraph(Copy, *OverlappingMap, 0);
    else
      // On last iteration, we want to filter out all columns except 
      // those that correspond
      // to rows in the graph.  This assures that our matrix is square
      OverlappingGraph = new Epetra_CrsGraph(Copy, *OverlappingMap, 
                                          *OverlappingMap, 0);

    OverlappingGraph->Import(*OldGraph, *OverlappingImporter, Insert);
    if (level < OverlappingLevel) 
      OverlappingGraph->TransformToLocal(DomainMap, RangeMap);
    else {
      // Copy last OverlapImporter because we will use it later
      OverlappingImporter = new Epetra_Import(*OverlappingMap, *DomainMap);
      OverlappingGraph->TransformToLocal(DomainMap, RangeMap);
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
string Ifpack_toString(const int& x)
{
  char s[100];
  sprintf(s, "%d", x);
  return string(s);
}

//============================================================================
string Ifpack_toString(const double& x)
{
  char s[100];
  sprintf(s, "%g", x);
  return string(s);
}

//============================================================================
int Ifpack_PrintResidual(char* Label, const Epetra_RowMatrix& A,
                         const Epetra_MultiVector& X, const Epetra_MultiVector&Y)
{
  if (X.Comm().MyPID() == 0) {
    cout << "***** " << Label << endl;
  }
  Ifpack_PrintResidual(0,A,X,Y);

  return(0);
}

//============================================================================
int Ifpack_PrintResidual(const int iter, const Epetra_RowMatrix& A,
                         const Epetra_MultiVector& X, const Epetra_MultiVector&Y)
{
  Epetra_MultiVector RHS(X);
  std::vector<double> Norm2;
  Norm2.resize(X.NumVectors());

  IFPACK_CHK_ERR(A.Multiply(false,X,RHS));
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
void Ifpack_PrintSparsity_Simple(const Epetra_RowMatrix& A)
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
// First written for ML_operators, this explains the C style
#include "limits.h"
#include "float.h"
int Ifpack_Analyze(Epetra_RowMatrix& A, const int NumEquations)
{

  int i,j;
  double MyFrobeniusNorm = 0.0; 
  double FrobeniusNorm = 0.0;
  double MyMinElement = DBL_MAX; 
  double MinElement = DBL_MAX;
  double MyMaxElement = DBL_MIN; 
  double MaxElement = DBL_MIN;
  double MyMinAbsElement = DBL_MAX; 
  double MinAbsElement = DBL_MAX;
  double MyMaxAbsElement = 0.0; 
  double MaxAbsElement = 0.0;
  int Nonzeros; /* nonzero elements in a row */
  int MyMinNonzeros = INT_MAX;
  int MyMaxNonzeros = 0;
  int MinNonzeros = INT_MAX;
  int MaxNonzeros = 0;
  int NumMyRows = 0, NumGlobalRows = 0;
  
  int MyPID;
  int NumProc;

  double Element, AbsElement; /* generic nonzero element and its abs value */
  int DiagonallyDominant = 0;
  int MyDiagonallyDominant = 0;
  int WeaklyDiagonallyDominant = 0;
  int MyWeaklyDiagonallyDominant = 0;
  int MyDirichletRows = 0;
  int DirichletRows = 0;
  int MyLowerNonzeros = 0;
  int MyUpperNonzeros = 0;
  int LowerNonzeros = 0;
  int UpperNonzeros = 0;
  double Min;
  double MyMin;
  double Max;
  double MyMax;
  int MyBandwidth = 0;
  int Bandwidth = 0;
  int grow;
  int gcol;
  int MyActualNonzeros = 0;
  int ActualNonzeros = 0;
  const Epetra_Comm& Comm = A.Comm();

  /* ---------------------- execution begins ------------------------------ */

  // only square matrices. The code should work for
  // rectangular matrices as well, but it has never been tested.
  if( A.NumGlobalRows() != A.NumGlobalCols() ) 
    IFPACK_CHK_ERR(-1);

  /* set up data and allocate memory */
  MyPID = A.Comm().MyPID();
  NumProc = A.Comm().NumProc();
  NumMyRows = A.NumMyRows();
  NumGlobalRows = A.NumGlobalRows();

  vector<int> colInd(A.MaxNumEntries());
  vector<double> colVal(A.MaxNumEntries());

  vector<double> Diagonal(NumMyRows);
  vector<double> SumOffDiagonal(NumMyRows);
  
  for (i = 0 ; i < NumMyRows ; ++i) {
    Diagonal[i] = 0.0;
    SumOffDiagonal[i] = 0.0;
  }
 
  // zero-out values
  MyActualNonzeros = 0;
  
  /* cycle over all matrix rows */

  for (int i = 0 ; i < NumMyRows ; ++i) {

    grow = A.RowMatrixRowMap().GID(i);
    assert (grow != -1);

    int Nnz;
    IFPACK_CHK_ERR(A.ExtractMyRowCopy(i,A.MaxNumEntries(),Nnz,
                                      &colVal[0],&colInd[0]));
                                      
    /* compute the real number of nonzero elements */
    
    int count = 0;
    for (j = 0 ; j < Nnz ; ++j) {
      /* compute the real number of nonzeros only */
      if (colVal[j] != 0.0) {
        ++count;
        gcol = A.RowMatrixColMap().GID(colInd[j]);
        assert (gcol != -1);
        if (gcol < grow) 
          MyLowerNonzeros++;
        else if (gcol > grow) 
          MyUpperNonzeros++;
        /* compute bandwidth */
        if (IFPACK_ABS((gcol - grow)) > MyBandwidth)
          MyBandwidth = IFPACK_ABS((gcol - grow));
      }
    }
    MyActualNonzeros += count;

    if (count == 1) MyDirichletRows++;
    
    /* start looking for min/max element (with and without abs()) */
    /* I consider the following:
     * - min element;
     * - min abs element;
     * - max element;
     * - max abs element;
     * - diagonal element;
     * - sum off abs off-diagonal elements;
     * Here I compute the local ones, with prefix `My'. The global
     * ones (without prefix) will be computed later
     */

    for (j = 0 ; j < Nnz ; ++j) {
      Element = colVal[j];
      if (Element != 0.0) {
        AbsElement = abs(Element);
        if ((Element < MyMinElement) && (Element != 0.0))
          MyMinElement = Element;
        if (Element > MyMaxElement) 
          MyMaxElement = Element;
        if ((AbsElement < MyMinAbsElement) && (AbsElement != 0.0))
          MyMinAbsElement = AbsElement;
        if (AbsElement > MyMaxAbsElement) 
          MyMaxAbsElement = AbsElement;
        if (colInd[j] == i)
          Diagonal[i] = AbsElement;
        else
          SumOffDiagonal[i] += abs(Element);
        MyFrobeniusNorm += Element*Element;
      }
    }
  } /* for over all matrix rows */

  /* compute the min/max of important quantities over all processes */

  Comm.MinAll(&MyMinElement,&MinElement,1);
  Comm.MaxAll(&MyMaxElement,&MaxElement,1);
  Comm.MinAll(&MyMinAbsElement,&MinAbsElement,1);
  Comm.MaxAll(&MyMaxAbsElement,&MaxAbsElement,1);
  
  Comm.SumAll(&MyLowerNonzeros,&LowerNonzeros,1);
  Comm.SumAll(&MyUpperNonzeros,&UpperNonzeros,1);

  Comm.SumAll(&MyFrobeniusNorm,&FrobeniusNorm,1);
  Comm.SumAll(&MyDirichletRows,&DirichletRows,1);
  Comm.SumAll(&MyActualNonzeros,&ActualNonzeros,1);

  Comm.MaxAll(&MyBandwidth,&Bandwidth,1);

  /* a test to see if matrix is diagonally-dominant */

  MyDiagonallyDominant = 0;
  MyWeaklyDiagonallyDominant = 0;

  for (i = 0 ; i < NumMyRows ; ++i) {
    if (Diagonal[i] > SumOffDiagonal[i]) 
      ++MyDiagonallyDominant;
    else if (Diagonal[i] == SumOffDiagonal[i]) 
      ++MyWeaklyDiagonallyDominant;
    /* else nothing to track */
  }

  Comm.SumAll(&MyDiagonallyDominant,&DiagonallyDominant,1);
  Comm.SumAll(&MyWeaklyDiagonallyDominant,&WeaklyDiagonallyDominant,1);

  /* simply no output for MyPID>0, only proc 0 write on os */
  if (MyPID == 0) {

    cout << "\n\n\t*** IFPACK Analysis of Epetra_RowMatrix `" << A.Label() << "'" << endl;
    printf("\t%-50s = %d\n", 
           "Number of global rows", NumGlobalRows);
    printf("\t%-50s = %d\n",
           "Number of stored elements", A.NumGlobalNonzeros());
    printf("\t%-50s = %d\n",
           "Number of nonzero elements", ActualNonzeros);
    printf("\t%-50s = %f\n",
           "Average number of nonzero elements/rows", 
           1.0 * A.NumGlobalNonzeros() / NumGlobalRows);
    printf("\t%-50s = %d\n",
           "Nonzero elements in strict lower part", LowerNonzeros);
    printf("\t%-50s = %d\n",
           "Nonzero elements in strict upper part", UpperNonzeros);
    printf("\t%-50s = %d\n",
           "Max |i-j|, a(i,j) != 0",Bandwidth);
    printf("\t%-50s = %d (= %5.2f%%)\n",
           "Number of diagonally dominant rows",
           DiagonallyDominant,
           100.0*DiagonallyDominant/NumGlobalRows); 
    printf("\t%-50s = %d (= %5.2f%%)\n",
           "Number of weakly diagonally dominant rows",
           WeaklyDiagonallyDominant,
           100.0*WeaklyDiagonallyDominant/NumGlobalRows);
    printf("\t%-50s = %d (= %5.2f%%)\n",
           "Number of Dirichlet rows",
           DirichletRows,
           100.0*DirichletRows/NumGlobalRows);
    printf("\t%-50s = %f\n",
           "||A||_F",sqrt(FrobeniusNorm));
    printf("\t%-50s = %f\n",
           "Min_{i,j} ( a(i,j) )", MinElement);
    printf("\t%-50s = %f\n",
           "Max_{i,j} ( a(i,j) )", MaxElement);
    printf("\t%-50s = %f\n",
           "Min_{i,j} ( abs(a(i,j)) )", MinAbsElement);
    printf("\t%-50s = %f\n",
           "Max_{i,j} ( abs(a(i,j)) )", MaxAbsElement);

  }

  /* Analyze elements on diagonal for the entire matrix */

  MyMin = DBL_MAX, MyMax = 0.0;
  for( i=0 ; i<NumMyRows ; ++i ) {
    if (Diagonal[i] < MyMin ) {
      if (Diagonal[i] != 0.0) 
        MyMin = Diagonal[i];
    }
    if (Diagonal[i] > MyMax) 
      MyMax = Diagonal[i];
  }
  Comm.MinAll(&MyMin,&Min,1);
  Comm.MaxAll(&MyMax,&Max,1);

  if (MyPID == 0) {
    printf("\t%-50s = %f\n",
           "Min_i ( abs(a(i,i)) )", Min);
    printf("\t%-50s = %f\n",
           "Max_i ( abs(a(i,i)) )", Max);
  }

  /* Analyze elements off diagonal for the entire matrix */

  MyMin = DBL_MAX, MyMax = 0.0;
  for (i = 0 ; i < NumMyRows ; ++i) {
    if (SumOffDiagonal[i] < MyMin ) 
      if (SumOffDiagonal[i] != 0.0) 
        MyMin = SumOffDiagonal[i];
    if (SumOffDiagonal[i] > MyMax) 
      MyMax = SumOffDiagonal[i];
  }
  Comm.MinAll(&MyMin,&Min,1);
  Comm.MaxAll(&MyMax,&Max,1);

  if (MyPID == 0) {
    printf("\t%-50s = %f\n",
           "Min_i ( \\sum_{j!=i} abs(a(i,j)) )", Min);
    printf("\t%-50s = %f\n",
           "Max_i ( \\sum_{j!=i} abs(a(i,j)) )", Max);
  }

  /* cycle over all equations and analyze diagonal elements. 
   * This may show that the matrix is badly scaled */
  
  if (NumEquations > 1) {

    for (int Equation = 0 ; Equation < NumEquations ; ++Equation) {

      /* Analyze elements on diagonal */

      MyMin = DBL_MAX, MyMax = 0.0;
      for (i = Equation ; i < NumMyRows ; i += NumEquations) {
        if (Diagonal[i] < MyMin) {
          if (Diagonal[i] != 0.0)
            MyMin = Diagonal[i];
        }
        if (Diagonal[i] > MyMax) 
          MyMax = Diagonal[i];
      }
      Comm.MinAll(&MyMin,&Min,1);
      Comm.MaxAll(&MyMax,&Max,1);

      if (MyPID == 0) {
        printf("\t(Eq %2d) %-42s = %f\n",
               Equation,
               "Min_i ( abs(a(i,i)) )", 
               Min);
        printf("\t(Eq %2d) %-42s = %f\n",
               Equation,
               "Max_i ( abs(a(i,i)) )", 
               Max);
      }
    }
  }

  if (MyPID == 0)
    cout << endl;

  return 0;

}

// ====================================================================== 
int Ifpack_PrintSparsity(const Epetra_RowMatrix& A, char* title,
                         char* FileName,
                         int NumPDEEqns)
{

  int m,nc,nr,maxdim,ltit;
  double lrmrgn,botmrgn,xtit,ytit,ytitof,fnstit,siz;
  double xl,xr, yb,yt, scfct,u2dot,frlw,delt,paperx,xx,yy;
  bool square = false;
  /*change square to .true. if you prefer a square frame around
    a rectangular matrix */
  double haf = 0.5, zero = 0.0, conv = 2.54;
  char munt = 'E'; /* put 'E' for centimeters, 'U' for inches */
  int ptitle = 0; /* position of the title, 0 under the drawing,
                     else above */
  FILE* fp = NULL;
  int NumMyRows;
  int NumMyCols;
  int NumGlobalRows;
  int NumGlobalCols;
  int MyPID;
  int NumProc;
  
  const Epetra_Comm& Comm = A.Comm();

  /* --------------------- execution begins ---------------------- */

  MyPID = Comm.MyPID();
  NumProc = Comm.NumProc();

  NumMyRows = A.NumMyRows();
  NumMyCols = A.NumMyCols();

  NumGlobalRows = A.NumGlobalRows();
  NumGlobalCols = A.NumGlobalCols();

  if (NumGlobalRows != NumGlobalCols)
    IFPACK_CHK_ERR(-1); // never tested

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
    fprintf(fp,"%%%%Creator: IFPACK\n");
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
