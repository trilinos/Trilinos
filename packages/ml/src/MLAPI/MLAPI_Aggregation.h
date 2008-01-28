#ifndef MLAPI_AGGREGATION_H
#define MLAPI_AGGREGATION_H

#include "ml_common.h"

/*!
\file MLAPI_Aggregation.h

\brief Functions to create aggregation-based prolongator operators.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/

namespace Teuchos {
  class ParameterList;
}

namespace MLAPI {

class Operator;
class MultiVector;

//! Builds the tentative prolongator using aggregation.
void GetPtent(const Operator& A, Teuchos::ParameterList& List,
              const MultiVector& ThisNS, 
              Operator& Ptent, MultiVector& NextNS);

//! Builds the tentative prolongator with default null space.
void GetPtent(const Operator& A, Teuchos::ParameterList& List, Operator& Ptent);

/*!
\brief Call ML aggregation on A according to parameters supplied in List. Return
       aggregates in aggrinfo
       
       On input, map of aggrinfo has to map row map of A. On output, aggrinfo[i]
       contains number of aggregate the row belongs to, where aggregates are
       numbered starting from 0. 
       Return value is the processor-local number of aggregates build.
       If aggrinfo[i] >= return-value, then i is a processor local row
       of a row that ML has detected to be on a Dirichlet BC.
       
\param A (in): Matrix to be aggregated on
\param List (in): ParameterList containing ML options
\param ThisNS (in): nullspace
\param aggrinfo(out): vector containing aggregation information

\note Map of aggrinfo has to match rowmap of A on input.

\return returns processor-local number of aggregates

\author Michael Gee (gee@lnm.mw.tum.de)
*/
int GetAggregates(const Operator& A, Teuchos::ParameterList& List,
                   const MultiVector& ThisNS, Epetra_IntVector& aggrinfo);

} // namespace MLAPI

#endif // MLAPI_AGGREGATION_H
