#ifndef MLAPI_GALLERY_H
#define MLAPI_GALLERY_H
#include "ml_include.h"
#include <iostream>
#include "MLAPI_Space.h"
#include "MLAPI_DoubleVector.h"
#include "MLAPI_Workspace.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace std;

namespace MLAPI {

Operator Gallery(const string ProblemType,
                 const Space& MySpace)
{
  int NumGlobalElements = MySpace.NumGlobalElements();
  Trilinos_Util::CrsMatrixGallery Gallery(ProblemType.c_str(), GetEpetraComm());
  Gallery.Set("problem_size", NumGlobalElements);
  Epetra_CrsMatrix* EpetraA = new Epetra_CrsMatrix(*(Gallery.GetMatrix()));
  
  EpetraA->FillComplete();

  Operator A(MySpace,MySpace,EpetraA);
  return (A);

}
} // namespace MLAPI
#endif // MLAPI_GALLERY_H
