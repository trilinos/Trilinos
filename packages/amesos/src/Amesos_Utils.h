#ifndef AMESOS_UTILS_H
#define AMESOS_UTILS_H

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"

/*!
 \class Amesos_Utils
 
 \brief Amesos_Utils: Collections of basic utilities.

 \author Marzio Sala, SNL 9214

 \date Last updated on 24-May-05 (Champions' League Final day)
*/

class Amesos_Utils
{
public:
  //! Default constructor.
  Amesos_Utils() {}
  
  //! Default destructor.
  ~Amesos_Utils() {}

  //! Computes the true residual, B - Matrix * X, and prints the results.
  void ComputeTrueResidual(const Epetra_RowMatrix& Matrix,
                           const Epetra_MultiVector& X,
                           const Epetra_MultiVector& B,
                           const bool UseTranspose,
                           const std::string prefix) const
  {
    double Norm;
    Epetra_Vector Ax(B.Map());
    int NumVectors = X.NumVectors();

    for (int i = 0 ; i < NumVectors ; ++i) 
    {
      Matrix.Multiply(UseTranspose, *X(i), Ax);
      Ax.Update(1.0, *B(i), -1.0);
      Ax.Norm2(&Norm);

      if (Matrix.Comm().MyPID() == 0) 
        std::cout << prefix << " : vector " << i << ", ||Ax - b|| = " 
          << Norm << std::endl;
    }
  }

  //! Computes the norms of X and B and print the results.
  void ComputeVectorNorms(const Epetra_MultiVector& X,
                          const Epetra_MultiVector& B,
                          const std::string prefix) const
  {
    double NormLHS;
    double NormRHS;
    int NumVectors = X.NumVectors();

    for (int i = 0 ; i < NumVectors ; ++i) 
    {
      X(i)->Norm2(&NormLHS);
      B(i)->Norm2(&NormRHS);
      if (X.Comm().MyPID() == 0) 
        std::cout << prefix << " : vector " << i << ", ||x|| = " << NormLHS
          << ", ||b|| = " << NormRHS << std::endl;
    }
  }

  //! Prints line on std::cout.
  void PrintLine() const
  {
    std::cout << "--------------------------------------------";
    std::cout << "--------------------------------" << std::endl;
  }

  void SetMaxProcesses(int& MaxProcesses, const Epetra_RowMatrix& A)
  {
    int MaxProcs = A.Comm().NumProc();

    switch(MaxProcesses) {
    case -3:
      MaxProcesses = MaxProcs;
      break;
    case -2:
      MaxProcesses = (int) sqrt(1.0 * MaxProcs);
      break;
    case -1:
      MaxProcesses = 1 + static_cast<int>(EPETRA_MAX(A.NumGlobalRows64() / 10000, 
                                    A.NumGlobalNonzeros64() / 1000000));
      break;
    }

    if (MaxProcesses <= 0) MaxProcesses = 1;
    if (MaxProcesses > MaxProcs) MaxProcesses = MaxProcs;

    return;
  }

};
#endif
