/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_common.h"
#ifdef HAVE_ML_MLAPI

#include "ml_include.h"
#include <iostream>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_LAPACK.h"
#include "Epetra_CrsMatrix.h"
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
  /* old version
  int    NodesPerAggr    = List.get("aggregation: per aggregate", 64);
  */
  double Threshold       = List.get("aggregation: threshold", 0.0);
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
  ML_memory_alloc((void **)(&null_vect), sizeof(double) * size * ThisNS.GetNumVectors(), "ns");

  int incr = 1;
  for (int v = 0 ; v < ThisNS.GetNumVectors() ; ++v)
    DCOPY_F77(&size, (double*)ThisNS.GetValues(v), &incr,
              null_vect + v * ThisNS.GetMyLength(), &incr);


  ML_Aggregate_Set_NullSpace(agg_object, NumPDEEquations,
                             ThisNS.GetNumVectors(), null_vect, 
                             ThisNS.GetMyLength());

  if (CoarsenType == "Uncoupled") 
    agg_object->coarsen_scheme = ML_AGGR_UNCOUPLED;
  else if (CoarsenType == "Uncoupled-MIS")
    agg_object->coarsen_scheme = ML_AGGR_HYBRIDUM;
  else if (CoarsenType == "MIS") {
   /* needed for MIS, otherwise it sets the number of equations to
    * the null space dimension */
    agg_object->max_levels  = -7;
    agg_object->coarsen_scheme = ML_AGGR_MIS;
  }
  else if (CoarsenType == "METIS")
    agg_object->coarsen_scheme = ML_AGGR_METIS;
  else {
    ML_THROW("Requested aggregation scheme (" + CoarsenType +
             ") not recognized", -1);
  }

  int NextSize = ML_Aggregate_Coarsen(agg_object, A.GetML_Operator(), 
                                      &ML_Ptent, GetML_Comm());

  /* This is the old version
  int NextSize;
  
  if (CoarsenType == "Uncoupled") {
    NextSize = ML_Aggregate_CoarsenUncoupled(agg_object, A.GetML_Operator(),
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
  */

  ML_Operator_ChangeToSinglePrecision(ML_Ptent);

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
  ML_memory_free((void**)(&null_vect));
}

// ====================================================================== 
void GetPtent(const Operator& A, Teuchos::ParameterList& List, Operator& Ptent)
{
  MultiVector FineNS(A.GetDomainSpace(),1);
  FineNS = 1.0;
  MultiVector CoarseNS;

  GetPtent(A, List, FineNS, Ptent, CoarseNS);
}

// ====================================================================== 
int GetAggregates(Epetra_RowMatrix& A, Teuchos::ParameterList& List,
                  double* thisns, Epetra_IntVector& aggrinfo)
{
  if (!A.RowMatrixRowMap().SameAs(aggrinfo.Map()))
    ML_THROW("map of aggrinfo must match row map of operator", -1);
  
  string CoarsenType     = List.get("aggregation: type", "Uncoupled");
  double Threshold       = List.get("aggregation: threshold", 0.0);
  int    NumPDEEquations = List.get("PDE equations", 1);
  int    nsdim           = List.get("null space: dimension",-1);
  if (nsdim==-1) ML_THROW("dimension of nullspace not set", -1);
  int size = A.RowMatrixRowMap().NumMyElements();

  ML_Aggregate* agg_object;
  ML_Aggregate_Create(&agg_object);
  ML_Aggregate_KeepInfo(agg_object,1);
  ML_Aggregate_Set_MaxLevels(agg_object,2);
  ML_Aggregate_Set_StartLevel(agg_object,0);
  ML_Aggregate_Set_Threshold(agg_object,Threshold);
  //agg_object->curr_threshold = 0.0;
  
  ML_Operator* ML_Ptent = 0;
  ML_Ptent = ML_Operator_Create(GetML_Comm());
  
  if (!thisns)
    ML_THROW("nullspace is NULL", -1);
             
  ML_Aggregate_Set_NullSpace(agg_object, NumPDEEquations,
                             nsdim, thisns,size);

  if (CoarsenType == "Uncoupled") 
    agg_object->coarsen_scheme = ML_AGGR_UNCOUPLED;
  else if (CoarsenType == "Uncoupled-MIS")
    agg_object->coarsen_scheme = ML_AGGR_HYBRIDUM;
  else if (CoarsenType == "MIS") {
   /* needed for MIS, otherwise it sets the number of equations to
    * the null space dimension */
    agg_object->max_levels  = -7;
    agg_object->coarsen_scheme = ML_AGGR_MIS;
  }
  else if (CoarsenType == "METIS")
    agg_object->coarsen_scheme = ML_AGGR_METIS;
  else {
    ML_THROW("Requested aggregation scheme (" + CoarsenType +
             ") not recognized", -1);
  }

  ML_Operator* ML_A = ML_Operator_Create(GetML_Comm());
  ML_Operator_WrapEpetraMatrix(&A,ML_A);

  int NextSize = ML_Aggregate_Coarsen(agg_object, ML_A, 
                                      &ML_Ptent, GetML_Comm());

  int* aggrmap = NULL;
  ML_Aggregate_Get_AggrMap(agg_object,0,&aggrmap);
  if (!aggrmap) ML_THROW("aggr_info not available", -1);

#if 0 // debugging  
  fflush(stdout);
  for (int proc=0; proc<A.GetRowMatrix()->Comm().NumProc(); ++proc)
  {
    if (A.GetRowMatrix()->Comm().MyPID()==proc)
    {
      cout << "Proc " << proc << ":" << endl;
      cout << "aggrcount " << aggrcount << endl;
      cout << "NextSize " << NextSize << endl;
      for (int i=0; i<size; ++i)
        cout << "aggrmap[" << i << "] = " << aggrmap[i] << endl;
      fflush(stdout);
    }
    A.GetRowMatrix()->Comm().Barrier();
  }
#endif
  
  assert (NextSize * nsdim != 0);
  for (int i=0; i<size; ++i) aggrinfo[i] = aggrmap[i];

  ML_Aggregate_Destroy(&agg_object);

  return (NextSize/nsdim);
}

// ====================================================================== 
int GetGlobalAggregates(Epetra_RowMatrix& A, Teuchos::ParameterList& List,
                        double* thisns, Epetra_IntVector& aggrinfo)
{
  int naggregates = GetAggregates(A,List,thisns,aggrinfo);
  const Epetra_Comm& comm = A.Comm();
  vector<int> local(comm.NumProc());
  vector<int> global(comm.NumProc());
  for (int i=0; i<comm.NumProc(); ++i) local[i] = 0;
  local[comm.MyPID()] = naggregates;
  comm.SumAll(&local[0],&global[0],comm.NumProc());
  int offset = 0;
  for (int i=0; i<comm.MyPID(); ++i) offset += global[i];
  for (int i=0; i<aggrinfo.MyLength(); ++i)
    if (aggrinfo[i]<naggregates) aggrinfo[i] += offset;
    else                         aggrinfo[i] = -1;
  return naggregates;
}

// ====================================================================== 
void GetPtent(const Epetra_RowMatrix& A, Teuchos::ParameterList& List,
              double* thisns, Teuchos::RCP<Epetra_CrsMatrix>& Ptent, 
              Teuchos::RCP<Epetra_MultiVector>& NextNS, const int domainoffset)
{
  const int nsdim = List.get<int>("null space: dimension",-1);
  if (nsdim<=0) ML_THROW("null space dimension not given",-1);
  const Epetra_Map& rowmap = A.RowMatrixRowMap();
  const int mylength = rowmap.NumMyElements();
  
  // wrap nullspace into something more handy
  Epetra_MultiVector ns(View,rowmap,thisns,mylength,nsdim);
  
  // vector to hold aggregation information
  Epetra_IntVector aggs(rowmap,false);
  
  // aggregation with global aggregate numbering
  int naggregates = GetGlobalAggregates(const_cast<Epetra_RowMatrix&>(A),List,thisns,aggs);
  
  // build a domain map for Ptent
  // find first aggregate on proc
  int firstagg = -1;
  int offset = -1;
  for (int i=0; i<mylength; ++i)
    if (aggs[i]>=0)
    {
      offset = firstagg = aggs[i];
      break;
    }
  offset *= nsdim;
  if (offset<0) ML_THROW("could not find any aggregate on proc",-2);
  vector<int> coarsegids(naggregates*nsdim);
  for (int i=0; i<naggregates; ++i)
    for (int j=0; j<nsdim; ++j)
    {
      coarsegids[i*nsdim+j] = offset + domainoffset;
      ++offset;
    }
  Epetra_Map pdomainmap(-1,naggregates*nsdim,&coarsegids[0],0,A.Comm());
  
  // loop aggregates and build ids for dofs
  map<int,vector<int> > aggdofs;
  map<int,vector<int> >::iterator fool;
  for (int i=0; i<naggregates; ++i)
  {
    vector<int> gids(0);
    aggdofs.insert(pair<int,vector<int> >(firstagg+i,gids));
  }
  for (int i=0; i<mylength; ++i)
  {
    if (aggs[i]<0) continue;
    vector<int>& gids = aggdofs[aggs[i]];
    gids.push_back(aggs.Map().GID(i));
  }
  
#if 0 // debugging output    
  for (int proc=0; proc<A.Comm().NumProc(); ++proc)
  {
    if (proc==A.Comm().MyPID())
    {
      for (fool=aggdofs.begin(); fool!=aggdofs.end(); ++fool)
      {
        cout << "Proc " << proc << " Aggregate " << fool->first << " Dofs ";
        vector<int>& gids = fool->second;
        for (int i=0; i<(int)gids.size(); ++i) cout << gids[i] << " ";
        cout << endl;
      }
    }
    fflush(stdout);
    A.Comm().Barrier();
  }
#endif
  
  // coarse level nullspace to be filled
  NextNS = Teuchos::rcp(new Epetra_MultiVector(pdomainmap,nsdim,true));
  Epetra_MultiVector& nextns = *NextNS;
  
  // matrix Ptent
  Ptent = Teuchos::rcp(new Epetra_CrsMatrix(Copy,rowmap,nsdim));
  
  // loop aggregates and extract the appropriate slices of the nullspace.
  // do QR and assemble Q and R to Ptent and NextNS
  for (fool=aggdofs.begin(); fool!=aggdofs.end(); ++fool)
  {
    // extract aggregate-local junk of nullspace
    const int aggsize = (int)fool->second.size();
    Epetra_SerialDenseMatrix Bagg(aggsize,nsdim);
    for (int i=0; i<aggsize; ++i)
      for (int j=0; j<nsdim; ++j)
        Bagg(i,j) = (*ns(j))[ns.Map().LID(fool->second[i])];
    
    // Bagg = Q*R
    int m = Bagg.M();
    int n = Bagg.N();
    int lwork = n*10;
    int info = 0;
    int k = min(m,n);
    if (k!=n) ML_THROW("Aggregate too small, fatal!",-1);
    vector<double> work(lwork);
    vector<double> tau(k);
    Epetra_LAPACK lapack;
    lapack.GEQRF(m,n,Bagg.A(),m,&tau[0],&work[0],lwork,&info);
    if (info) ML_THROW("Lapack dgeqrf returned nonzero",info);
    if (work[0]>lwork)
    {
      lwork = (int)work[0];
      work.resize(lwork);
    }
    
    // get R (stored on Bagg) and assemble it into nextns
    int agg_cgid = fool->first*nsdim;
    if (!nextns.Map().MyGID(agg_cgid+domainoffset))
      ML_THROW("Missing coarse column id on this proc",-1);
    for (int i=0; i<n; ++i)
      for (int j=i; j<n; ++j)
        (*nextns(j))[nextns.Map().LID(domainoffset+agg_cgid+i)] = Bagg(i,j);
    
    // get Q and assemble it into Ptent
    lapack.ORGQR(m,n,k,Bagg.A(),m,&tau[0],&work[0],lwork,&info);
    if (info) ML_THROW("Lapack dorgqr returned nonzero",info);
    for (int i=0; i<aggsize; ++i)
    {
      const int actgrow = fool->second[i];
      for (int j=0; j<nsdim; ++j)
      {
        int actgcol = fool->first*nsdim+j+domainoffset;
        int errone = Ptent->SumIntoGlobalValues(actgrow,1,&Bagg(i,j),&actgcol);
        if (errone>0)
        {
          int errtwo = Ptent->InsertGlobalValues(actgrow,1,&Bagg(i,j),&actgcol);
          if (errtwo<0) ML_THROW("Epetra_CrsMatrix::InsertGlobalValues returned negative nonzero",errtwo);
        }
        else if (errone) ML_THROW("Epetra_CrsMatrix::SumIntoGlobalValues returned negative nonzero",errone);
      }
    } // for (int i=0; i<aggsize; ++i)
  } // for (fool=aggdofs.begin(); fool!=aggdofs.end(); ++fool)
  int err = Ptent->FillComplete(pdomainmap,rowmap);
  if (err) ML_THROW("Epetra_CrsMatrix::FillComplete returned nonzero",err);
  err = Ptent->OptimizeStorage();
  if (err) ML_THROW("Epetra_CrsMatrix::OptimizeStorage returned nonzero",err);
  
  return;
}              

// ====================================================================== 
void GetPtent(const Epetra_RowMatrix& A, Teuchos::ParameterList& List,
              double* thisns, Teuchos::RCP<Epetra_CrsMatrix>& Ptent, 
              const int domainoffset)
{
  Teuchos::RCP<Epetra_MultiVector> tmp;
  MLAPI::GetPtent(A,List,thisns,Ptent,tmp,domainoffset);
  return;
}              
              


} // namespace MLAPI

#endif // HAVE_ML_MLAPI
