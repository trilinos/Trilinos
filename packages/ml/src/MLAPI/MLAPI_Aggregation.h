#ifndef MLAPI_AGGREGATION_H
#define MLAPI_AGGREGATION_H
#include "ml_include.h"
#include <iostream>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_DataBase.h"

/*!
\file MLAPI_Aggregation.h

\brief Functions to create aggregation-based prolongator operators.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/

using namespace std;

namespace Teuchos {
  class ParameterList;
}

namespace MLAPI {

#include "ml_aggregate.h"
#include "ml_agg_METIS.h"

//! Builds the tentative prolongator using aggregation.
void BuildPtent(const Operator& A, const AggregationDataBase& Data,
                const NullSpace& ThisNS, 
                Operator& Ptent, NullSpace& NextNS)
{

  string CoarsenType  = Data.GetType();
  int    NodesPerAggr = Data.GetNodesPerAggregate();
  double Threshold    = Data.GetThreshold();

  ML_Aggregate* agg_object;
  ML_Aggregate_Create(&agg_object);
  ML_Aggregate_Set_MaxLevels(agg_object,2);
  ML_Aggregate_Set_StartLevel(agg_object,0);
  ML_Aggregate_Set_Threshold(agg_object,Threshold);
  //agg_object->curr_threshold = 0.0;
  
  ML_Operator* ML_Ptent = 0;
  ML_Ptent = ML_Operator_Create(GetML_Comm());

  if (ThisNS.GetDimension() == 0)
    ML_THROW("zero-dimension null space", -1);
             
  int size = A.GetDomainSpace().GetNumMyElements() * ThisNS.GetDimension();
  // HOW TO FREE THIS MEMORY??
  ML_memory_alloc((void **)&(agg_object->nullspace_vect), 
                  sizeof(double) * size, "ns");
  cout << size << endl;
  for (int i = 0 ; i < size ; ++i)
    agg_object->nullspace_vect[i] = 1.0 ; //////ThisNS.Values()[i];

  agg_object->nullspace_dim = ThisNS.GetDimension();

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

  // FIXME: this is broken
  assert (NextSize * ThisNS.GetDimension() != 0);

  double* NextNullSpace = new double[NextSize * ThisNS.GetDimension()];

  for (int i = 0 ; i < NextSize * ThisNS.GetDimension() ; ++i)
    NextNullSpace[i] = agg_object->nullspace_vect[i];

  ML_Aggregate_Destroy(&agg_object);

  // FIXME: what is the new dimension of the null space??
  NextNS.Reshape(CoarseSpace, ThisNS.GetDimension(), NextNullSpace,
                 true);
}

//! Builds the tentative prolongator using aggregation and default null space
void BuildPtent(const Operator& A, const AggregationDataBase& Data,
                Operator& Ptent)
{
  NullSpace FineNS(A.GetDomainSpace(),1);
  NullSpace CoarseNS;

  BuildPtent(A, Data, FineNS, Ptent, CoarseNS);
  
}

} // namespace MLAPI

#endif // MLAPI_AGGREGATION_H
