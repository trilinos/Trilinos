#ifndef ML_THROW
#define ML_THROW(str,val) { \
  std::cerr << "ERROR: In function/method " << __func__ << "()" << endl; \
  std::cerr << "ERROR: File " << __FILE__ << ", line " << __LINE__ << endl; \
  std::cerr << "ERROR: " << str << endl; \
  throw(val); \
  }
#endif

#include "MLAPI_Workspace.h"
#include "MLAPI_BaseObject.h"
#include "MLAPI_Space.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_MultiVector_Utils.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Operator_Utils.h"
#include "MLAPI_InverseOperator.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_BaseOperator.h"
#include "MLAPI_EpetraBaseOperator.h"
#include "MLAPI_MultiLevelSA.h"
#include "MLAPI_MATLABStream.h"
#include "MLAPI_Gallery.h"
#include "MLAPI_Krylov.h"
#include "MLAPI_Aggregation.h"
#include "MLAPI_Eig.h"
#include "MLAPI_Matrix.h"
