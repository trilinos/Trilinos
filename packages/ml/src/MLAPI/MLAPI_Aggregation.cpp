#include "ml_common.h"
#ifdef HAVE_ML_MLAPI

#include "ml_include.h"
#include <iostream>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "MLAPI_Error.h"
#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Aggregation.h"

using namespace std;

namespace MLAPI {

#include "ml_aggregate.h"
#include "ml_agg_METIS.h"

// ====================================================================== 
void GetPtent(const Operator& A, Teuchos::ParameterList& List,
              const MultiVector& ThisNS, 
              Operator& Ptent, MultiVector& NextNS)
{

  string CoarsenType     = List.get("aggregation: type", "Uncoupled");
  int    NodesPerAggr    = List.get("aggregation: per aggregate", 64);
  double Threshold       = List.get("aggregation:", 0.0);
  int    NumPDEEquations = List.get("PDE equations", 1);

  ML_Aggregate* agg_object;
  ML_Aggregate_Create(&agg_object);
  ML_Aggregate_Set_MaxLevels(agg_object,2);
  ML_Aggregate_Set_StartLevel(agg_object,0);
  ML_Aggregate_Set_Threshold(agg_object,Threshold);
  //agg_object->curr_threshold = 0.0;
  
  ML_Operator* ML_Ptent = 0;
  ML_Ptent = ML_Operator_Create(GetML_Comm());

  if (ThisNS.GetNumVectors() == 0)
    ML_THROW("zero-dimension null space", -1);
             
  int size = ThisNS.GetMyLength();

  double* null_vect = 0;
  ML_memory_alloc((void **)&null_vect, sizeof(double) * size * ThisNS.GetNumVectors(), "ns");

  int incr = 1;
  for (int v = 0 ; v < ThisNS.GetNumVectors() ; ++v)
    DCOPY_F77(&size, (double*)ThisNS.GetValues(v), &incr,
              null_vect + v * ThisNS.GetMyLength(), &incr);

  ML_Aggregate_Set_NullSpace(agg_object, NumPDEEquations,
                             ThisNS.GetNumVectors(), null_vect, 
                             ThisNS.GetMyLength());

  //DELETE agg_object->nullspace_dim = ThisNS.GetNumVectors();
  //agg_object->num_PDE_eqns = NumPDEEquations;

  int NextSize;
  
  if (CoarsenType == "Uncoupled") {
    NextSize = ML_Aggregate_CoarsenUncoupled(agg_object, A.GetML_Operator(),
                                             &ML_Ptent, GetML_Comm());
  }
  else if (CoarsenType == "MIS") {
    NextSize = ML_Aggregate_CoarsenMIS(agg_object, A.GetML_Operator(),
                                       &ML_Ptent, GetML_Comm());
  }
  else if (CoarsenType == "METIS") {
    ML ml_object;
    ml_object.ML_num_levels = 1; // crap for line below
    ML_Aggregate_Set_NodesPerAggr(&ml_object,agg_object,0,NodesPerAggr);
    NextSize = ML_Aggregate_CoarsenMETIS(agg_object, A.GetML_Operator(),
                                         &ML_Ptent, GetML_Comm());
  }
  else {
    ML_THROW("Requested aggregation scheme (" + CoarsenType +
             ") not recognized", -1);
  }

  int NumMyElements = NextSize;
  Space CoarseSpace(-1,NumMyElements);
  Ptent.Reshape(CoarseSpace,A.GetRangeSpace(),ML_Ptent,true);

  assert (NextSize * ThisNS.GetNumVectors() != 0);

  NextNS.Reshape(CoarseSpace, ThisNS.GetNumVectors());

  size = NextNS.GetMyLength();
  for (int v = 0 ; v < NextNS.GetNumVectors() ; ++v)
    DCOPY_F77(&size, agg_object->nullspace_vect + v * size, &incr,
              NextNS.GetValues(v), &incr);

  ML_Aggregate_Destroy(&agg_object);
  ML_memory_free(&null_vect);

}

// ====================================================================== 
void GetPtent(const Operator& A, Teuchos::ParameterList& List, Operator& Ptent)
{
  MultiVector FineNS(A.GetDomainSpace(),1);
  FineNS = 1.0;
  MultiVector CoarseNS;

  GetPtent(A, List, FineNS, Ptent, CoarseNS);
  
}

} // namespace MLAPI

#endif // HAVE_ML_MLAPI
