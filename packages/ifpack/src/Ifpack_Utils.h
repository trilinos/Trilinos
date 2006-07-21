#ifndef IFPACK_UTILS_H
#define IFPACK_UTILS_H

#include "Ifpack_ConfigDefs.h"
#include "Epetra_Comm.h"
#if !( defined(__INTEL_COMPILER) && defined(_WIN32) )
#  include "unistd.h" // Not a standard header file!
#endif
class Epetra_RowMatrix;
class Epetra_CrsMatrix;
class Epetra_CrsGraph;
class Epetra_RowMatrix;
class Epetra_MultiVector;
class Epetra_Vector;

/*! \file Ifpack_Utils.h
 */

//! Prints a line of `=' on cout
void Ifpack_PrintLine();

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

//! Analyzes the basic properties of the input matrix A; see \ref ifp_analyze.
int Ifpack_Analyze(const Epetra_RowMatrix& A, const bool Cheap = false,
                   const int NumPDEEqns = 1);

//! Analyzes the distribution of values of the input matrix A.
/*!
 \param A - (In) matrix to be analyzed.
 \param abs - (In) if \c true, the function will analyze matrix
              B, whose elements are defined as \f$ B_{i,i} = | A_{i,i}| \f$.
 \param steps - (In) number of intervals for the analysis.

 An example of output is reported \ref ifp_matrix.
 */
int Ifpack_AnalyzeMatrixElements(const Epetra_RowMatrix& A,
                                 const bool abs = false, 
                                 const int steps = 10);

//! Analyzes the distribution of values of the input vector Diagonal.
/*!
 \param Diagonal - (In) Vector to be analyzed.
 \param abs - (In) if \c true, the function will analyze vector
              B, whose elements are defined as \f$ B_{i} = | D_{i}| \f$.
 \param steps - (In) number of intervals for the analysis.

 An example of output is reported \ref ifp_vector.
 */
int Ifpack_AnalyzeVectorElements(const Epetra_Vector& Diagonal,
                                 const bool abs = false, 
                                 const int steps = 10);

//! Plots the sparsity pattern of an Epetra_RowMatrix into a PS file.
/*!
 \param A (In) - Epetra_RowMatrix whose sparsity pattern will be plotted.

 \param FileName (In) - char string containing the filename.
                        If 0, then the matrix label is used as file name,
                        after appending .ps.

 \param NumPDEEqns (In) - number of PDE equations. The function will plot
               the block structure of the matrix if NumPDEEqns > 1

 \name Largely inspired from Yousef Saad's SPARSKIT plot function. 
 */
int Ifpack_PrintSparsity(const Epetra_RowMatrix& A, const char* FileName = 0, 
                         const int NumPDEEqns = 1);

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
    return(false);
  }

private:
  int i_;
  double val_;
  double aval_;

};

#endif // IFPACK_UTILS_H
