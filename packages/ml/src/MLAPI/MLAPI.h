/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#include "MLAPI_Error.h"
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
#include "MLAPI_MultiLevelAdaptiveSA.h"
#include "MLAPI_MATLABStream.h"
#include "MLAPI_Gallery.h"
#include "MLAPI_Krylov.h"
#include "MLAPI_Aggregation.h"
#include "MLAPI_Eig.h"
#include "MLAPI_SerialMatrix.h"
#include "MLAPI_DistributedMatrix.h"
#include "MLAPI_BaseLinearCombination.h"
#include "MLAPI_LinearCombinations.h"
#include "MLAPI_Defaults.h"

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif
