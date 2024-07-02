// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_UTILS_HPP
#define IFPACK2_UTILS_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include <cmath>
#include "Tpetra_Comm.hpp"
#if !( defined(__INTEL_COMPILER) && defined(_WIN32) )
#  include "unistd.hpp" // Not a standard header file!
#endif
class Tpetra_RowMatrix;
class Tpetra_CrsMatrix;
class Tpetra_CrsGraph;
class Tpetra_RowMatrix;
class Tpetra_MultiVector;
class Tpetra_Vector;

/*! \file Ifpack2_Utils.hpp
 */

//! Prints a line of `=' on cout
void Ifpack2_PrintLine();

//! Stops the execution of code, so that a debugger can be attached.
void Ifpack2_BreakForDebugger(Tpetra_Comm& Comm);

//! Creates an overlapping Tpetra_CrsMatrix. Returns 0 if OverlappingLevel is 0.
Tpetra_CrsMatrix* Ifpack2_CreateOverlappingCrsMatrix(const Tpetra_RowMatrix* Matrix,
						    const int OverlappingLevel);

//! Creates an overlapping Tpetra_CrsGraph. Returns 0 if OverlappingLevel is 0.
Tpetra_CrsGraph* Ifpack2_CreateOverlappingCrsMatrix(const Tpetra_CrsGraph* Graph,
						   const int OverlappingLevel);

//! Convertes an integer to string.
string Ifpack2_toString(const int& x);

//! Convertes a double to string.
string Ifpack2_toString(const double& x);

//! Prints on cout the true residual.
int Ifpack2_PrintResidual(char* Label,  const Tpetra_RowMatrix& A,
                         const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&Y);

int Ifpack2_PrintResidual(const int iter, const Tpetra_RowMatrix& A,
                         const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, const Tpetra_MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&Y);

void Ifpack2_PrintSparsity_Simple(const Tpetra_RowMatrix& A);

//! Analyzes the basic properties of the input matrix A; see \ref ifp_analyze.
int Ifpack2_Analyze(const Tpetra_RowMatrix& A, const bool Cheap = false,
                   const int NumPDEEqns = 1);

//! Analyzes the distribution of values of the input matrix A.
/*!
 \param A - (In) matrix to be analyzed.
 \param abs - (In) if \c true, the function will analyze matrix
              B, whose elements are defined as \f$ B_{i,i} = | A_{i,i}| \f$.
 \param steps - (In) number of intervals for the analysis.

 An example of output is reported \ref ifp_matrix.
 */
int Ifpack2_AnalyzeMatrixElements(const Tpetra_RowMatrix& A,
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
int Ifpack2_AnalyzeVectorElements(const Tpetra_Vector& Diagonal,
                                 const bool abs = false, 
                                 const int steps = 10);

//! Plots the sparsity pattern of an Tpetra_RowMatrix into a PS file.
/*!
 \param A (In) - Tpetra_RowMatrix whose sparsity pattern will be plotted.

 \param FileName (In) - char string containing the filename.
                        If 0, then the matrix label is used as file name,
                        after appending .ps.

 \param NumPDEEqns (In) - number of PDE equations. The function will plot
               the block structure of the matrix if NumPDEEqns > 1

 \name Largely inspired from Yousef Saad's SPARSKIT plot function. 
 */
int Ifpack2_PrintSparsity(const Tpetra_RowMatrix& A, const char* FileName = 0, 
                         const int NumPDEEqns = 1);

//==============================================================================
class Ifpack2_Element {

public:
  Ifpack2_Element() {};

  Ifpack2_Element(const Ifpack2_Element& rhs) {
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
    aval_ = std::abs(val_);
  }

  inline bool operator <(const Ifpack2_Element& rhs) const 
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

#endif // IFPACK2_UTILS_HPP
