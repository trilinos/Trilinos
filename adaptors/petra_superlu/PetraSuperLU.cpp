#include "PetraSuperLU.h"

int PetraSuperLU(Epetra_LinearProblem * Problem) {

  // Build Epetra_Comm object

  Epetra_CrsMatrix * Epetra_A = dynamic_cast<Epetra_CrsMatrix *> (Problem->GetOperator());
  Epetra_MultiVector * Epetra_X = Problem->GetLHS();
  Epetra_MultiVector * Epetra_B = Problem->GetRHS();

  Epetra_CrsGraph AG = Epetra_A->Graph(); // Get matrix graph

  int NumIndices;
  int GlobalMaxNumIndices = AG.GlobalMaxNumIndices();
  int * Indices = new int[GlobalMaxNumIndices];
  double * Values = new double[GlobalMaxNumIndices];
  int NumMyCols = AG.NumMyCols();
  int NumMyEquations = Epetra_A->NumMyRows();
  int * TransNumNz = new int[NumMyCols];
  for (int i=0;i<NumMyCols; i++) TransNumNz[i] = 0;
  for (int i=0; i<NumMyEquations; i++) {
    assert(AG.ExtractMyRowCopy(i, GlobalMaxNumIndices, NumIndices, Indices)==0); // Get copy of ith row
    for (int j=0; j<NumIndices; j++) ++TransNumNz[Indices[j]];
  }

  int ** TransIndices = new int*[NumMyCols];
  double ** TransValues = new double*[NumMyCols];

  // Count total number of nonzeros
  int TotalNumNz = 0;
  for(int i=0; i<NumMyCols; i++) TotalNumNz += TransNumNz[i];

  // Allocate space for matrix values and indices in arrays that SuperLU will use
  int * RowIndices = new int[TotalNumNz];
  double * RowValues = new double[TotalNumNz];
  int * ColPointers = new int[NumMyCols+1];

  // Set pointers into the RowIndices and Values arrays, define ColPointers
  NumIndices = 0;
  for(int i=0; i<NumMyCols; i++) {
    ColPointers[i] = NumIndices;
    TransIndices[i] = RowIndices + NumIndices;
    TransValues[i] = RowValues + NumIndices;
    NumIndices += TransNumNz[i];
  }
  ColPointers[NumMyCols] = NumIndices;

  // Now copy values and global indices into newly create transpose storage

  for (int i=0;i<NumMyCols; i++) TransNumNz[i] = 0; // Reset transpose NumNz counter
  for (int i=0; i<NumMyEquations; i++) {
    assert(Epetra_A->ExtractMyRowCopy(i, GlobalMaxNumIndices, NumIndices, Values, Indices)==0);
    int ii = Epetra_A->GRID(i);
    for (int j=0; j<NumIndices; j++) {
      int TransRow = Indices[j];
      int loc = TransNumNz[TransRow];
      TransIndices[TransRow][loc] = ii;
      TransValues[TransRow][loc] = Values[j];
      ++TransNumNz[TransRow]; // increment counter into current transpose row
    }
  }

  delete [] Indices;
  delete [] Values;
  delete [] TransIndices;
  delete [] TransValues;


  SuperMatrix A, L, U, B;

  int m = Epetra_A->NumGlobalRows();
  int n = Epetra_A->NumGlobalCols();
  int nnz = TotalNumNz;
  double * a = RowValues;
  int * asub = RowIndices;
  int * xa = ColPointers;
  /* Create matrix A in the format expected by SuperLU. */
  dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, NC, _D, GE);
  
  /* Create right-hand side matrix B. */
  int nrhs = 1;
  double * rhs;
  int LDA;
  Epetra_X->Update(1.0, *Epetra_B, 0.0); // Copy B into X
  Epetra_X->ExtractView(&rhs, &LDA);
  dCreate_Dense_Matrix(&B, m, nrhs, rhs, LDA, DN, _D, GE);
  if (m<200) dPrint_Dense_Matrix("B", &B);
  
  int * perm_r = new int[m];
  int * perm_c = new int[n];
  
  /*
   * Get column permutation vector perm_c[], according to permc_spec:
   *   permc_spec = 0: use the natural ordering
   *   permc_spec = 1: use minimum degree ordering on structure of A'*A
   *   permc_spec = 2: use minimum degree ordering on structure of A'+A
   */
  Epetra_Time timer(Epetra_A->Comm());
  int permc_spec = 2;
  get_perm_c(permc_spec, &A, perm_c);
  double permTime = timer.ElapsedTime();
  
  int info = 0;
  double t0 = timer.ElapsedTime();
  dgssv(&A, perm_c, perm_r, &L, &U, &B, &info);
  double solveTime = timer.ElapsedTime() - t0;

  cout << "Permutation time = " << permTime << endl
       << "Solve time       = " << solveTime << endl;
  
  if (m<200) dPrint_CompCol_Matrix("A", &A);
  if (m<200) dPrint_CompCol_Matrix("U", &U);
  if (m<200) dPrint_SuperNode_Matrix("L", &L);
  if (m<200) PrintInt10("\nperm_r", m, perm_r);
  if (m<200) dPrint_Dense_Matrix("X", &B);
  
  delete [] RowValues;
  delete [] RowIndices;
  delete [] ColPointers;
  delete [] perm_r;
  delete [] perm_c;
  /* De-allocate storage */
  // Destroy_CompCol_Matrix(&A);  This get rid of RowValues, etc.  Already done
  Destroy_SuperMatrix_Store(&A);
  Destroy_SuperMatrix_Store(&B);
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);
  return 0;
}
