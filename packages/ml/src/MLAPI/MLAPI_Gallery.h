#ifndef MLAPI_GALLERY_H
#define MLAPI_GALLERY_H

#include "ml_common.h"
#ifdef HAVE_ML_MLAPI

#include "ml_include.h"
#include <iostream>

namespace MLAPI {

  class Space;
  class Operator;

  Operator Gallery(const string ProblemType, const Space& MySpace);
}

#endif // HAVE_ML_MLAPI

#endif // MLAPI_GALLERY_H
