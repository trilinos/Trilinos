#ifndef MLAPI_GALLERY_H
#define MLAPI_GALLERY_H

#include "ml_common.h"
#ifdef HAVE_ML_MLAPI

#include "ml_include.h"
#include <iostream>
#include "MLAPI_Error.h"
#include "MLAPI_Space.h"
#include "MLAPI_Workspace.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace std;

namespace MLAPI {

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
} // namespace MLAPI

#endif // HAVE_ML_MLAPI

#endif // MLAPI_GALLERY_H
