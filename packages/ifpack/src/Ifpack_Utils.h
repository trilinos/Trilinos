#ifndef IFPACK_UTILS_H
#define IFPACK_UTILS_H

#include "Ifpack_ConfigDefs.h"
#include "Epetra_Comm.h"
#include "unistd.h"
class Epetra_RowMatrix;
class Epetra_CrsMatrix;
class Epetra_CrsGraph;
class Epetra_RowMatrix;
class Epetra_MultiVector;

//! Stops the execution of code, so that a debugger can be attached.
void Ifpack_BreakForDebugger(Epetra_Comm& Comm);

//! Creates an overlapping Epetra_CrsMatrix. Returns 0 if OverlappingLevel is 0.
Epetra_CrsMatrix* Ifpack_CreateOverlappingCrsMatrix(const Epetra_RowMatrix* Matrix,
						    const int OverlappingLevel);

//! Creates an overlapping Epetra_CrsGraph. Returns 0 if OverlappingLevel is 0.
Epetra_CrsGraph* Ifpack_CreateOverlappingCrsMatrix(const Epetra_CrsGraph* Graph,
						   const int OverlappingLevel);

//! Convertes an integer to string.
string Ifpack_toString(const int& x);

//! Convertes a double to string.
string Ifpack_toString(const double& x);

//! Prints on cout the true residual.
int Ifpack_PrintResidual(char* Label,  const Epetra_RowMatrix& A,
                         const Epetra_MultiVector& X, const Epetra_MultiVector&Y);

int Ifpack_PrintResidual(const int iter, const Epetra_RowMatrix& A,
                         const Epetra_MultiVector& X, const Epetra_MultiVector&Y);

void Ifpack_PrintSparsity_Simple(const Epetra_RowMatrix& A);

int Ifpack_Analyze(const Epetra_RowMatrix& A, const int NumEquations = 1);
int Ifpack_PrintSparsity(const Epetra_RowMatrix& A, char* title,
                       char* FileName,
                       int NumPDEEqns = 1);

//==============================================================================
class Ifpack_Element {

public:
  Ifpack_Element() {};

  Ifpack_Element(const Ifpack_Element& rhs) {
    i_ = rhs.Index();
    val_ = rhs.Value();
    aval_ = rhs.AbsValue();
  }

  inline int Index() const {
    return(i_);
  }

  inline double Value() const {
    return(val_);
  }

  inline double AbsValue() const {
    return(aval_);
  }

  inline void SetIndex(const int i)
  {
    i_ = i;
  }

  inline void SetValue(const double val)
  {
    val_ = val;
    aval_ = IFPACK_ABS(val_);
  }

  inline bool operator <(const Ifpack_Element& rhs) const 
  {
    if (rhs.AbsValue() > AbsValue())
      return(false);
    else if (rhs.AbsValue() < AbsValue())
      return(true);
    else if (rhs.Index() < Index())
        return(true);
  }

private:
  int i_;
  double val_;
  double aval_;

};

#endif // IFPACK_UTILS_H
