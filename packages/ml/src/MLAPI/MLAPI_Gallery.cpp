#include "ml_common.h"
#ifdef HAVE_ML_MLAPI
#include "ml_include.h"
#include <iostream>
#include "MLAPI_Error.h"
#include "MLAPI_Space.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Gallery.h"
#include "MLAPI_Operator.h"
#include "MLAPI_DistributedMatrix.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace std;

namespace MLAPI {

// ====================================================================== 
Operator Gallery(const string ProblemType,
                 const Space& MySpace)
{
#if defined(HAVE_ML_TRIUTILS)
  int NumGlobalElements = MySpace.GetNumGlobalElements();
  Trilinos_Util::CrsMatrixGallery Gallery(ProblemType.c_str(), GetEpetra_Comm());
  Gallery.Set("problem_size", NumGlobalElements);
  Epetra_CrsMatrix* EpetraA = new Epetra_CrsMatrix(*(Gallery.GetMatrix()));
  
  EpetraA->FillComplete();

  Operator A(MySpace,MySpace,EpetraA);
  return (A);
#else
  ML_THROW("Configure with --enable-triutils", -1);
#endif

}

// ====================================================================== 
Operator GetShiftedLaplacian1D(const int NumGlobalElements,
                               const double Factor)
{

  double Pi = atan(1.0) * 4;
  double LambdaMin = 4.0 * pow((sin(1.0 / (NumGlobalElements + 1) * Pi / 2)),
                               2.0) * Factor;

  Space FineSpace(NumGlobalElements);

  DistributedMatrix* MatA = new DistributedMatrix(FineSpace, FineSpace);

  if (GetMyPID() == 0) {
    for (int i = 0 ; i < NumGlobalElements ; ++i) {
      MatA->SetElement(i, i, 2.0 - LambdaMin);
      if (i)
        MatA->SetElement(i, i - 1, - 1.0);
      if (i != NumGlobalElements - 1)
        MatA->SetElement(i, i + 1, - 1.0);
    }
  }

  MatA->FillComplete();

  Operator A(FineSpace, FineSpace, MatA, true);

  return(A);
}

// ====================================================================== 
Operator GetShiftedLaplacian2D(const int NX, const int NY, 
                               const double Factor,
                               const bool RandomScale)
{

  int NumGlobalElements = NX * NY;

  double Pi = atan(1.0) * 4;
  double LambdaMin = (4.0 * pow((sin(1.0 / (NX + 1) * Pi / 2)), 2.0) + 4.0 * pow((sin(1.0 / (NY + 1) * Pi / 2)), 2.0)) * Factor;

  Space FineSpace(NumGlobalElements);
  int n = 0;
  if (GetMyPID() == 0)
    n = NumGlobalElements;
  Space LocalizedSpace(-1, n);

  MultiVector Scale(LocalizedSpace);

  DistributedMatrix* MatA = new DistributedMatrix(FineSpace, FineSpace);

  // assemble the matrix on processor 0 only
  if (GetMyPID() == 0) {

    if (RandomScale) {

      Scale.Random();
      for (int i = 0 ; i < Scale.GetMyLength() ; ++i)
        Scale(i) = 1.0001 + Scale(i);

      double sr, sc;
      for (int i = 0 ; i < NX ; ++i) {
        for (int j = 0 ; j < NY ; ++j) {
          int row = i + j * NX;
          sr = Scale(row);
          MatA->SetElement(row, row, sr * sr * (4.0 - LambdaMin));
          if (i > 0) {
            sc = Scale(row - 1);
            MatA->SetElement(row, row - 1, sr * sc * (-1.0));
          }
          if (i < NX - 1) {
            sc = Scale(row + 1);
            MatA->SetElement(row, row + 1, sr * sc * (-1.0));
          }
          if (j > 0) {
            sc = Scale(row - NX);
            MatA->SetElement(row, row - NX, sr * sc * (-1.0));
          }
          if (j < NY - 1) {
            sc = Scale(row + NX);
            MatA->SetElement(row, row + NX, sr * sc * (-1.0));
          }
        }
      }
    }
    else {

      for (int i = 0 ; i < NX ; ++i) {
        for (int j = 0 ; j < NY ; ++j) {
          int row = i + j * NX;
          MatA->SetElement(row, row, 4.0 - LambdaMin);
          if (i > 0) 
            MatA->SetElement(row, row - 1, -1.0);
          if (i < NX - 1)
            MatA->SetElement(row, row + 1, -1.0);
          if (j > 0) 
            MatA->SetElement(row, row - NX, -1.0);
          if (j < NY - 1)
            MatA->SetElement(row, row + NX, -1.0);
        }
      }
    }

  } // mypid == 0

  MatA->FillComplete();

  Operator A(FineSpace, FineSpace, MatA, true);

  return (A);
}

// ====================================================================== 
Operator ReadMatrix(const char* FileName)
{

  int NumGlobalElements;
  int NumGlobalNonzeros;
  std::ifstream File;
  File.open(FileName);

  if (!File.good() )
    ML_THROW("Error opening file", -1);

  File >> NumGlobalElements;
  File >> NumGlobalNonzeros;

  if (NumGlobalElements <= 0)
    ML_THROW("Invalid number of global elements (" + GetString(NumGlobalElements) + ")", -1);

  if (NumGlobalNonzeros <= 0)
    ML_THROW("Invalid number of global nonzeros" + GetString(NumGlobalNonzeros) + ")", -1);

  if (GetMyPID()) File.close();

  Space FineSpace(NumGlobalElements);

  DistributedMatrix* MatA = new DistributedMatrix(FineSpace, FineSpace);

  if (GetMyPID() == 0) {
    int row, col;
    double val;
    for (int i = 0 ; i < NumGlobalNonzeros ; ++i) {
      File >> row;
      File >> col;
      File >> val;

      if (row < 0 || row >= NumGlobalElements)
        ML_THROW("Invalid row number (" + GetString(row) + ")", -1);
      if (col < 0 || col >= NumGlobalElements)
        ML_THROW("Invalid col number (" + GetString(col) + ")", -1);

      MatA->SetElement(row, col, val);
    }
  }

  MatA->FillComplete();

  Operator A(FineSpace, FineSpace, MatA, true);

  return(A);
}

// ====================================================================== 
} // namespace MLAPI

#endif // HAVE_ML_MLAPI
