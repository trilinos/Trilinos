#ifndef MLAPI_AGGREGATION_H
#define MLAPI_AGGREGATION_H

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

// ====================================================================== 
//! Builds the tentative prolongator using aggregation.
// ====================================================================== 

void GetPtent(const Operator& A, Teuchos::ParameterList& List,
              const MultiVector& ThisNS, 
              Operator& Ptent, MultiVector& NextNS)
{

  // FIXME-RST
  // Ray, I would appreciate if you check and fix (or tell me how to fix)
  // the following:
  // - number of PDEs
  // - dimension of the null space
  // - settings of all the parameters to call ML_Aggregate_Coarsen()
  // - set/get of null space, before and after the call to ML_Aggregate_Coarsen()
  // - for output, how to set the current level?
  // - something else has to be fixed? Memory leaks?
  // Thanks
  
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
             
  int size = A.GetDomainSpace().GetNumMyElements() * ThisNS.GetNumVectors();

  // FIXME-RST HOW TO FREE THIS MEMORY??
  ML_memory_alloc((void **)&(agg_object->nullspace_vect), 
                  sizeof(double) * size, "ns");
  for (int i = 0 ; i < size ; ++i)
    agg_object->nullspace_vect[i] = ThisNS.GetValues()[i];

  agg_object->nullspace_dim = ThisNS.GetNumVectors();
  agg_object->num_PDE_eqns = NumPDEEquations;

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
  assert (NextSize * ThisNS.GetNumVectors() != 0);

  NextNS.Reshape(CoarseSpace, ThisNS.GetNumVectors());

  // FIXME: is it correct?
  for (int i = 0 ; i < NextNS.GetMyTotalLength() ; ++i)
    NextNS.GetValues()[i] = agg_object->nullspace_vect[i];

  ML_Aggregate_Destroy(&agg_object);

}

// ====================================================================== 
//! Builds the tentative prolongator with default null space.
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

#endif // MLAPI_AGGREGATION_H
