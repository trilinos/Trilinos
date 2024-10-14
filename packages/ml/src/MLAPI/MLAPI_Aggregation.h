#ifndef MLAPI_AGGREGATION_H
#define MLAPI_AGGREGATION_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_common.h"

/*!
\file MLAPI_Aggregation.h

\brief Functions to create aggregation-based prolongator operators.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

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
\brief Builds the tentative prolongator using aggregation.

Build Ptent and NextNS as usual but as Epetra objects.

\param A (in): Matrix to be aggregated on
\param List (in): ParameterList containing ML options
\param thisns (in): nullspace in format ML accepts
\param Ptent(out): Matrix containing tentative prolongator
\param NextNS (out): MultiVector containing next level nullspace.
\param domainoffset (in,optional): give an offset such that the domainmap
                                   of Ptent starts global numbering from domainoffset
                                   instead from zero. This is useful to
                                   create block matrices.

\author Michael Gee (gee\@lnm.mw.tum.de)
*/
void GetPtent(const Epetra_RowMatrix& A, Teuchos::ParameterList& List,
              double* thisns, Teuchos::RCP<Epetra_CrsMatrix>& Ptent,
              Teuchos::RCP<Epetra_MultiVector>& NextNS, const int domainoffset = 0);

/*!
\brief Builds the tentative prolongator using aggregation.

Build Ptent and NextNS as usual but as Epetra objects.

\param A (in): Matrix to be aggregated on
\param List (in): ParameterList containing ML options
\param thisns (in): nullspace in format ML accepts
\param Ptent(out): Matrix containing tentative prolongator
\param domainoffset (in,optional): give an offset such that the domainmap
                                   of Ptent starts global numbering from domainoffset
                                   instead from zero. This is useful to
                                   create block matrices.

\author Michael Gee (gee\@lnm.mw.tum.de)
*/
void GetPtent(const Epetra_RowMatrix& A, Teuchos::ParameterList& List,
              double* thisns, Teuchos::RCP<Epetra_CrsMatrix>& Ptent,
              const int domainoffset = 0);

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
\param thisns (in): nullspace
\param aggrinfo(out): vector containing aggregation information

\note Map of aggrinfo has to match rowmap of A on input.

\return returns processor-local number of aggregates

\author Michael Gee (gee\@lnm.mw.tum.de)
*/
int GetAggregates(Epetra_RowMatrix& A, Teuchos::ParameterList& List,
                   double* thisns, Epetra_IntVector& aggrinfo);

/*!
\brief Call ML aggregation on A according to parameters supplied in List. Return
       aggregates in aggrinfo

       On input, map of aggrinfo has to map row map of A. On output, aggrinfo[i]
       contains number of global aggregate the row belongs to, where aggregates are
       numbered starting from 0 globally.
       Return value is the processor-local number of aggregates build.
       If aggrinfo[i] < 0, then i is a processor local row
       that ML has detected to be on a Dirichlet BC.
       if aggrinfo[i] >= 0, then i is a processor local row and aggrinfo[i] is
       a global aggregate id.

\param A (in): Matrix to be aggregated on
\param List (in): ParameterList containing ML options
\param thisns (in): nullspace
\param aggrinfo(out): vector containing aggregation information in global numbering

\note Map of aggrinfo has to match rowmap of A on input.

\return returns processor-local number of aggregates

\author Michael Gee (gee\@lnm.mw.tum.de)
*/
int GetGlobalAggregates(Epetra_RowMatrix& A, Teuchos::ParameterList& List,
                        double* thisns, Epetra_IntVector& aggrinfo);
} // namespace MLAPI

#endif // MLAPI_AGGREGATION_H
