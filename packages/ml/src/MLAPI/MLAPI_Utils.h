#ifndef MLAPI_UTILS_H
#define MLAPI_UTILS_H

#include "ml_include.h"
#include "MLAPI_Operator.h"

namespace MLAPI {

static ML_Comm* MLAPI_Comm_ = 0;

void Init() {
  ML_Comm_Create(&MLAPI_Comm_);
}

void Finalize() {
  ML_Comm_Destroy(&MLAPI_Comm_);
}

ML_Comm* GetMLComm() {
  return(MLAPI_Comm_);
}

} // namespace MLAPI

#endif
