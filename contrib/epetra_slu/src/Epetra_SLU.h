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
              const int fill_fac = -1,
              const int panel_size = -1,
              const int relax = -1 );

  ~Epetra_SLU();

  int Solve( const bool Verbose = false,
             const bool Equil = true,
             const bool Factor = true,
             const int perm_type = 2,
             const double pivot_thresh = -1 );

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
