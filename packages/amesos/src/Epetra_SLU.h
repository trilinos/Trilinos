#ifndef _EPETRA_SLU_H_
#define _EPETRA_SLU_H_

class Epetra_LinearProblem;
class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_CrsGraph;

class SLUData;

class Epetra_SLU
{

 public:

  Epetra_SLU( Epetra_LinearProblem * Problem,
              int fill_fac = -1,
              int panel_size = -1,
              int relax = -1 );

  ~Epetra_SLU();

  int Solve( bool Verbose = false,
             bool Equil = true,
             bool Factor = true,
             int perm_type = 2,
             double pivot_thresh = -1,
             bool Refact = true,
             bool Trans = false );

 private:

  void Copy();

  Epetra_CrsMatrix * A_;
  Epetra_MultiVector * X_;
  Epetra_MultiVector * B_;

  double * R_;
  double * C_;

  Epetra_CrsGraph * AG_;

  SLUData * data_;

  int * Indices_;
  double * Values_;

  int * TransNumNZ_;
  int TotNumNZ_;

  int ** TransIndices_;
  double ** TransValues_;

  int * RowIndices_;
  double * RowValues_;
  int * ColPointers_;

  int * perm_r_;
  int * perm_c_;

  int count_;

  int * etree_;

  double * ferr_;
  double * berr_;
};

#endif
