/*!
 *  \file ml_MultiLevelPreconditioner.cpp
 *
 *  \brief ML black-box preconditioner for Epetra_RowMatrix derived classes.
 *
 *  \author Marzio Sala, SNL, 9214
 *
 *  \date Last update to Doxygen: 07-Jul-2006
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* Notes on memory analyzer (MS):
 * - if mallinfo() is available, ML_MALLINFO is automatically defined
 *   by configure. This should be reasonably cheap and accurate.
 * - define ML_MEM_CHECK to an expensive analysis of memory corruption
 *   and memory leaks, for all blocks allocated using ML's allocating
 *   functions. Users have to set this variable explicitly (not handled by
 *   configure).
 * - define ML_MALLOC to estimate the maximum free memory by calling
 *   malloc(). This may be severely slow on some machine. Not defined by
 *   default. Users have to set this variable explicitly (not handled by
 *   configure). This option gives sometimes "strange" results, to be
 *   used only as last resort!
 */

/*#############################################################################
# CVS File Information
#    Current revision: $Revision$
#    Branch:           $Branch$
#    Last modified:    $Date$
#    Modified by:      $Author$
#############################################################################*/

#include "ml_common.h"
#include "ml_include.h"
#include "ml_RowMatrix.h"

#ifdef IFPACK_NODE_AWARE_CODE
extern int ML_NODE_ID; //FIXME
#endif

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#include "ml_memory.h"
#include "ml_DD_prec.h"
#include <iostream>
#include <iomanip>

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Time.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "ml_amesos_wrap.h"
#include "ml_agg_METIS.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_agg_ParMETIS.h"
#include "ml_anasazi.h"
#include "ml_FilterType.h"
#include "ml_ValidateParameters.h"

#ifdef HAVE_ML_EPETRAEXT
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_OperatorOut.h"
#endif

#include "ml_ifpack_wrap.h"
#include "ml_viz_stats.h"
using namespace Teuchos;
using namespace std;

// ================================================ ====== ==== ==== == =
void ML_Epetra::MultiLevelPreconditioner::PrintMem(char *fmt, int min, int sum, int max)
{

  if (Comm().MyPID() == 0)
    printf(fmt,min,sum / Comm().NumProc(),max,sum);
  puts(" (Mb)");

  return;
}

// ================================================ ====== ==== ==== == =
int ML_Epetra::MultiLevelPreconditioner::DestroyPreconditioner()
{

  {
    int NumDestroy = OutputList_.get("number of destruction phases", 0);
    OutputList_.set("number of destruction phases", ++NumDestroy);
  }

  if (verbose_ && mlpLabel_ != "not-set")
    std::cout << "***" << std::endl
         << "*** destroying ML_Epetra::MultiLevelPreconditioner ["
         << mlpLabel_ << "]" << std::endl
         << "***" << std::endl;

  if (SubMatMLPrec_ != NULL) {
    for (int i = 0; i < NBlocks_ ; i++) delete SubMatMLPrec_[i];
    ML_free(SubMatMLPrec_);
  }

  ML_Aggregate_VizAndStats_Clean(ml_);
  if (ml_nodes_ != 0) ML_Aggregate_VizAndStats_Clean(ml_nodes_);

  // destroy main objects
  if (agg_ != 0) {
    // destroy aggregate information
    if ((agg_)->aggr_info != NULL) {
      for (int i = 0 ; i < NumLevels_ ; ++i) {
        if ((agg_)->aggr_info[i] != NULL)
          ML_memory_free((void **)&((agg_)->aggr_info[i]));
      }
    }
    ML_Aggregate_Destroy(&agg_); agg_ = 0;
  }

  if (TMatrixML_ != 0) {
    ML_Operator_Destroy(&TMatrixML_);
    TMatrixML_ = 0;
  }

  if (ML_Kn_ != 0) {
    ML_Operator_Destroy(&ML_Kn_);
    ML_Kn_ = 0;
  }

  if (TtATMatrixML_ != 0) {
    ML_Operator_Destroy(&TtATMatrixML_);
    TtATMatrixML_ = 0;
  }

  if (TMatrixTransposeML_ != 0) {
    ML_Operator_Destroy(&TMatrixTransposeML_);
    TMatrixTransposeML_ = 0;
  }

  if (Tmat_array != 0) {
    ML_MGHierarchy_ReitzingerDestroy(LevelID_[0]-1, &Tmat_array,
                     &Tmat_trans_array);
    Tmat_array = 0;
    Tmat_trans_array = 0;
  }

  if (nodal_args_ != 0) {
    ML_Smoother_Arglist_Delete(&nodal_args_);
    nodal_args_ = 0;
  }

  if (edge_args_  != 0) {
    ML_Smoother_Arglist_Delete(&edge_args_);
    edge_args_ = 0;
  }

  if (MassMatrix_array != 0) {
    for (int i=0; i<MaxLevels_; i++)
      ML_Operator_Destroy(MassMatrix_array+LevelID_[i]);
    delete [] MassMatrix_array;
    MassMatrix_array = 0;
  }

  if (CurlCurlMatrix_array != 0) {
    for (int i=0; i<MaxLevels_; i++)
      ML_Operator_Destroy(CurlCurlMatrix_array+LevelID_[i]);
    delete [] CurlCurlMatrix_array;
    CurlCurlMatrix_array = 0;
  }

  // These should be destroyed after all the individual operators, as the
  // communicator may be used in profiling some of the individual ML_Operators
  // in the hierarchy at when they are destroyed.
  if (ml_ != 0) { ML_Destroy(&ml_); ml_ = 0; }
  if (ml_nodes_ != 0) { ML_Destroy(&ml_nodes_); ml_nodes_ = 0; }

  if (Label_) {
    delete [] Label_;
    Label_ = 0;
  }

  if (CreatedTMatrix_ == true) {
    delete TMatrix_;
    TMatrix_ = NULL;
  }

  if (CreatedML_Kn_ == true) {
    delete ML_Kn_;
    ML_Kn_ = NULL;
  }

  if (CreatedNodeMatrix_ == true) {
    delete NodeMatrix_;
    NodeMatrix_ = NULL;
  }

  if (CreatedEdgeMatrix_) {
    delete RowMatrix_;
    RowMatrix_ = NULL;
  }

  if (ml_comm_ != 0) { ML_Comm_Destroy(&ml_comm_); ml_comm_ = 0; }

  int min[ML_MEM_SIZE], max[ML_MEM_SIZE], sum[ML_MEM_SIZE];
  for( int i=0 ; i<ML_MEM_SIZE ; ++i ) sum[i] = 0;

  if (AnalyzeMemory_) {
    memory_[ML_MEM_TOT2] = memory_[ML_MEM_TOT1] + memory_[ML_MEM_PREC_FIRST];
#ifdef ML_MALLOC
    memory_[ML_MEM_TOT2_MALLOC] = memory_[ML_MEM_TOT1_MALLOC] + memory_[ML_MEM_PREC_FIRST_MALLOC];
#endif
    Comm().MinAll(memory_,min,ML_MEM_SIZE);
    Comm().MaxAll(memory_,max,ML_MEM_SIZE);
    Comm().SumAll(memory_,sum,ML_MEM_SIZE);
  }

  if (verbose_ && NumApplications_)
    ReportTime();

  if (verbose_ && AnalyzeMemory_) {

    // print memory usage

    std::cout << std::endl;
    std::cout << "   ML memory information:                 min     avg     max       tot"
         << std::endl << std::endl;
#ifdef ML_MALLINFO
    std::cout << "   1- estimated ML memory usage, using mallinfo()" << std::endl;
    PrintMem("      for the hierarchy              = %5d   %5d   %5d   %7d",
         min[ML_MEM_HIERARCHY],sum[ML_MEM_HIERARCHY],max[ML_MEM_HIERARCHY]);
    PrintMem("      for the smoother(s)            = %5d   %5d   %5d   %7d",
         min[ML_MEM_SMOOTHER],sum[ML_MEM_SMOOTHER],max[ML_MEM_SMOOTHER]);
    PrintMem("      for the coarse solver          = %5d   %5d   %5d   %7d",
         min[ML_MEM_COARSE],sum[ML_MEM_COARSE],max[ML_MEM_COARSE]);
    PrintMem("      preconditioning                = %5d   %5d   %5d   %7d",
         min[ML_MEM_PREC_FIRST],sum[ML_MEM_PREC_FIRST],max[ML_MEM_PREC_FIRST]);
    PrintMem("      total (w/o other prec data)    = %5d   %5d   %5d   %7d",
         min[ML_MEM_TOT1],sum[ML_MEM_TOT1],max[ML_MEM_TOT1]);
    PrintMem("      total (w/  other prec data)    = %5d   %5d   %5d   %7d",
         min[ML_MEM_TOT2],sum[ML_MEM_TOT2],max[ML_MEM_TOT2]);
    std::cout << std::endl;
#endif
#ifdef ML_MALLOC
    std::cout << "   3- estimated ML memory usage, using malloc()" << std::endl;
    PrintMem("      for the hierarchy              = %5d   %5d   %5d   %7d",
         min[ML_MEM_HIERARCHY_MALLOC],sum[ML_MEM_HIERARCHY_MALLOC],max[ML_MEM_HIERARCHY_MALLOC]);
    PrintMem("      for the smoother(s)            = %5d   %5d   %5d   %7d",
         min[ML_MEM_SMOOTHER_MALLOC],sum[ML_MEM_SMOOTHER_MALLOC],max[ML_MEM_SMOOTHER_MALLOC]);
    PrintMem("      for the coarse solver          = %5d   %5d   %5d   %7d",
         min[ML_MEM_COARSE_MALLOC],sum[ML_MEM_COARSE_MALLOC],max[ML_MEM_COARSE_MALLOC]);
    PrintMem("      preconditioning                = %5d   %5d   %5d   %7d",
         min[ML_MEM_PREC_FIRST_MALLOC],sum[ML_MEM_PREC_FIRST_MALLOC],max[ML_MEM_PREC_FIRST_MALLOC]);
    PrintMem("      total (w/o other prec data)    = %5d   %5d   %5d   %7d",
         min[ML_MEM_TOT1_MALLOC],sum[ML_MEM_TOT1_MALLOC],max[ML_MEM_TOT1_MALLOC]);
    PrintMem("      total (w/  other prec data)    = %5d   %5d   %5d   %7d",
         min[ML_MEM_TOT2_MALLOC],sum[ML_MEM_TOT2_MALLOC],max[ML_MEM_TOT2_MALLOC]);
    std::cout << std::endl;
#endif
    std::cout << "      (memory for Aztec smoothers reported in `preconditioning')" << std::endl;
    std::cout << "      (warning: data may be incorrect for small runs, or" << std::endl
     << "       if more processes share the same physical memory)" << std::endl
     << std::endl;
  }

  if (verbose_ && NumApplications_)
    ML_print_line("-",78);

  if (NullSpaceToFree_ != 0) {
    delete [] NullSpaceToFree_;
    NullSpaceToFree_ = 0;
  }

  if (RowMatrixAllocated_) {
    delete RowMatrixAllocated_;
    RowMatrixAllocated_ = 0;
  }
  if (AllocatedRowMatrix_) {
    delete RowMatrix_;
    RowMatrixAllocated_ = 0;
  }

  // filtering stuff

  if (flt_ml_) {
    ML_Destroy(&flt_ml_);
    flt_ml_ = 0;
  }

  if (flt_agg_) {
    ML_Aggregate_Destroy(&flt_agg_);
    flt_agg_ = 0;
  }

  IsComputePreconditionerOK_ = false;

  return 0;
}

// ================================================ ====== ==== ==== == =

ML_Epetra::MultiLevelPreconditioner::
MultiLevelPreconditioner(const Epetra_RowMatrix & inRowMatrix,
             const bool ComputePrec) :
  RowMatrix_(&inRowMatrix),
  RowMatrixAllocated_(0),
  AllocatedRowMatrix_(false)
{

  ParameterList NewList;
  List_ = NewList;
  SetDefaults("SA",List_,(int *)0, (double *)0);

  ML_CHK_ERRV(Initialize());

  // construct hierarchy
  if (ComputePrec == true)
    ML_CHK_ERRV(ComputePreconditioner());
}

// ================================================ ====== ==== ==== == =

ML_Epetra::MultiLevelPreconditioner::
MultiLevelPreconditioner( const Epetra_RowMatrix & inRowMatrix,
             const ParameterList & List, const bool ComputePrec) :
  RowMatrix_(&inRowMatrix),
  RowMatrixAllocated_(0),
  AllocatedRowMatrix_(false)
{

  List_ = List;

  ML_CHK_ERRV(Initialize());

  // construct hierarchy
  if (ComputePrec == true)
    ML_CHK_ERRV(ComputePreconditioner());
}

// ================================================ ====== ==== ==== == =
ML_Epetra::MultiLevelPreconditioner::
MultiLevelPreconditioner(ML_Operator * Operator,
             const ParameterList & List, Epetra_RowMatrix **DiagOperators,
             Teuchos::ParameterList *DiagLists, int NBlocks,
             const bool ComputePrec) :
  DiagOperators_(DiagOperators),
  DiagLists_(DiagLists),
  RowMatrixAllocated_(0),
  AllocatedRowMatrix_(true)
{

  // need to wrap an Epetra_RowMatrix around Operator.
  // This is quite not the best approach for small matrices

  if (NBlocks != 0) Comm_ = &(DiagOperators[0]->Comm());
  else Comm_ = NULL;

  // matrix must be freed by destructor.

  RowMatrix_ = new ML_Epetra::RowMatrix(Operator,Comm_);
  NBlocks_ = NBlocks;
  List_ = List;

  ML_CHK_ERRV(Initialize());
  AMGSolver_ = ML_COMPOSITE;
  AfineML_ = Operator;

  // construct hierarchy
  if (ComputePrec == true)
    ML_CHK_ERRV(ComputePreconditioner());
}

// ================================================ ====== ==== ==== == =

/*! The constructor for the Maxwell equations requires three Epetra_RowMatrices.
 * Two conditions are required on their maps:
 * - TMatrix.OperatorDomainMap() == NodeMatrix.OperatorRangeMap()
 * - TMatrix.OperatorRangeMap()  == EdgeMatrix.OperatorDomainMap()
 */
ML_Epetra::MultiLevelPreconditioner::
MultiLevelPreconditioner(const Epetra_RowMatrix & EdgeMatrix,
             const Epetra_RowMatrix & TMatrix,
             const Epetra_RowMatrix & NodeMatrix,
             const ParameterList & List,
             const bool ComputePrec,
             const bool UseNodeMatrixForSmoother):
  RowMatrix_(&EdgeMatrix),
  RowMatrixAllocated_(0),
  AllocatedRowMatrix_(false)
{

  // some sanity checks
  if (! TMatrix.OperatorDomainMap().SameAs(NodeMatrix.OperatorRangeMap()) ) {
    std::cerr << ErrorMsg_ << "discrete grad DomainMap != node RangeMap..." << std::endl;
    ML_CHK_ERRV(-1); // error on discrete grad
  }

  if (! TMatrix.OperatorRangeMap().SameAs(EdgeMatrix.OperatorDomainMap()) ) {
    std::cerr << ErrorMsg_ << "discrete grad RangeMap != edge DomainMap..." <<std::endl;
    ML_CHK_ERRV(-2); // error on discrete grad
  }

  List_ = List;

  ML_CHK_ERRV(Initialize());

  // set Maxwell here.
  // NOTE: RowMatrix_ and EdgeMatrix_ pointer to the same Epetra_RowMatrix
  AMGSolver_ = ML_MAXWELL;
  NodeMatrix_ = & NodeMatrix;
  TMatrix_ = & TMatrix;
  EdgeMatrix_ = & EdgeMatrix;
  UseNodeMatrixForSmoother_=UseNodeMatrixForSmoother;

  // construct hierarchy
  if (ComputePrec == true)
    ML_CHK_ERRV(ComputePreconditioner());

}

// ================================================ ====== ==== ==== == =

/*! Constructor for the Maxwell equations.  This version takes the stiffness
 * and edge mass matrices separately.
 * Two conditions are required on their maps:
 * - TMatrix.OperatorDomainMap() == NodeMatrix.OperatorRangeMap()
 * - TMatrix.OperatorRangeMap()  == EdgeMatrix.OperatorDomainMap()
 */
ML_Epetra::MultiLevelPreconditioner::
MultiLevelPreconditioner(const Epetra_RowMatrix & CurlCurlMatrix,
                         const Epetra_RowMatrix & MassMatrix,
                         const Epetra_RowMatrix & TMatrix,
                         const Epetra_RowMatrix & NodeMatrix,
                         const ParameterList & List,
                         const bool ComputePrec) :
  RowMatrixAllocated_(0),
  AllocatedRowMatrix_(false)
{

  // Check compatibility of CurlCurl and  TMatrix.
  if (! TMatrix.OperatorDomainMap().SameAs(NodeMatrix.OperatorRangeMap()) ) {
    std::cerr << ErrorMsg_ << "discrete grad DomainMap != node RangeMap..." << std::endl;
    ML_CHK_ERRV(-1); // error on discrete grad
  }

  if (! TMatrix.OperatorRangeMap().SameAs(CurlCurlMatrix.OperatorDomainMap())) {
    std::cerr << ErrorMsg_ << "discrete grad RangeMap != edge DomainMap..." <<std::endl;
    ML_CHK_ERRV(-2); // error on discrete grad
  }

  List_ = List;

  // later on, this will be reset to the sum of the curl-curl and mass
  RowMatrix_ = & CurlCurlMatrix;

  ML_CHK_ERRV(Initialize());

  AMGSolver_ = ML_MAXWELL;
  EdgeMatrix_ = RowMatrix_;
  CurlCurlMatrix_ = & CurlCurlMatrix;
  MassMatrix_ = & MassMatrix;
  NodeMatrix_ = & NodeMatrix;
  TMatrix_ = & TMatrix;

  // construct hierarchy
  if (ComputePrec == true)
    ML_CHK_ERRV(ComputePreconditioner());

}
#ifdef NewStuff
// ================================================ ====== ==== ==== == =
/*! Constructor for multiphysics problems with variable dofs per node.  This version uses a
 *  discrete laplacian and padding on coarse grids to make a hierarchy. It also allows for
 *  the removal of column nonzeros associated with Dirichlet points. To use this option
 *  the rhs and initial guess must be provided.
 */
#ifdef outputter
extern int *update, *update_index, *external, *extern_index;
#endif
ML_Epetra::MultiLevelPreconditioner::
MultiLevelPreconditioner(Epetra_RowMatrix & RowMatrix,
                         const Teuchos::ParameterList & List,
                         const int & nNodes,
                         const int & maxDofPerNode,
                         bool * NotMLVecDofPresent,   // how do I get a const back ..
                         Epetra_Vector & Lhs,
                         Epetra_Vector & Rhs,
                         const bool  rhsAndsolProvided,
                         const bool ComputePrec) :
  RowMatrixAllocated_(0),
  AllocatedRowMatrix_(false)
{
  Epetra_CrsMatrix *Acrs =
              dynamic_cast<Epetra_CrsMatrix *>(&RowMatrix);

  TEUCHOS_TEST_FOR_EXCEPT_MSG(Acrs == NULL,
             ErrorMsg_ << "matrix must be crs to use this constructor\n");

  int nDofs = Acrs->NumMyRows();
  int *dofGlobals;

  Acrs->RowMap().MyGlobalElementsPtr(dofGlobals);

  List_ = List;

  int  nGlobalNodes;
  int itemp = nNodes;

  MLVec<bool> dofPresent(NotMLVecDofPresent,NotMLVecDofPresent+nNodes*maxDofPerNode);

  RowMatrix.Comm().SumAll(&itemp,&nGlobalNodes,1);

  double *XXX = NULL, *YYY = NULL, *ZZZ = NULL;

  XXX = List_.get("x-coordinates",(double *) 0);
  YYY = List_.get("y-coordinates",(double *) 0);
  ZZZ = List_.get("z-coordinates",(double *) 0);

  TEUCHOS_TEST_FOR_EXCEPT_MSG( (XXX == NULL) && (nDofs != 0),
             ErrorMsg_ << "Must supply coordinates to use multiphysics variable dof constructor\n");

  struct wrappedCommStruct epetraFramework;

  epetraFramework.whichOne = epetraType;
  epetraFramework.data     = (void *) Acrs;
  epetraFramework.nProcs   = RowMatrix.Comm().NumProc();
  epetraFramework.myPid    = RowMatrix.Comm().MyPID();
  int nGhost = MLnMyGhost(epetraFramework);

  // allocate vectors

  MLVec<bool>   dirOrNot(nDofs + nGhost);
  MLVec<double> theDiag(nDofs + nGhost);
  MLVec<int>    map(nDofs);
  MLVec<int>    amalgRowPtr, amalgCols;
  MLVec<int>    myLocalNodeIds(nDofs+nGhost);
  MLVec<int>    amalgRowMap;
  MLVec<int>    amalgColMap;

  // Make MLVec versions of Crs arrays by wrapping them.

  double *vtemp;  int *ctemp,*rtemp;
  Acrs->ExtractCrsDataPointers(rtemp, ctemp, vtemp );

  MLVec<int>    rowPtr(rtemp,rtemp + nDofs + 1);
  MLVec<int>    cols(ctemp,ctemp   + rtemp[nDofs]);
  MLVec<double> vals(vtemp,vtemp   + rtemp[nDofs]);

  // grab diagonal and detect Dirichlets

  double dirDropTol = List_.get("ML advanced Dirichlet: threshold",1.0e-5);

  MLextractDiag(rowPtr,   cols, vals, theDiag, epetraFramework);
  MLfindDirichlets(rowPtr,cols, vals, theDiag, dirDropTol, dirOrNot, epetraFramework);

  // if the rhs and solution are provided, remove column entries
  // associated with Dirichlet rows. Here, we overwrite the original
  // matrix arrays to reflect removal.

  if (rhsAndsolProvided) {

     double *ptrRHS, *ptrLHS;

     Lhs.ExtractView(&ptrLHS);
     Rhs.ExtractView(&ptrRHS);

     MLVec<double> myvecRHS(ptrRHS, ptrRHS + nDofs );
     MLVec<double> myvecLHS(ptrLHS, ptrLHS + nDofs );

     MLrmDirichletCols(rowPtr, cols, vals, theDiag, true,
                     myvecLHS, myvecRHS, dirOrNot, epetraFramework);
   }


   MLbuildMap(dofPresent, map, nDofs);

   for (int i = 0; i < nDofs; i++)
      myLocalNodeIds[i] = (int) floor( map[i]/maxDofPerNode);

   int nLocalNodes, nLocalPlusGhostNodes;

   MLassignGhostLocalNodeIds( myLocalNodeIds, nDofs, nDofs+nGhost,
        epetraFramework, nLocalNodes, nLocalPlusGhostNodes);

   MLfillNodalMaps(amalgRowMap, amalgColMap, myLocalNodeIds,
        nDofs,  epetraFramework,  nLocalNodes, nLocalPlusGhostNodes);

   double amalgDropTol = List_.get("variable DOF  amalgamation: threshold",1.8e-9);
   MLvariableDofAmalg(nDofs+nGhost, rowPtr, cols, vals, nNodes,maxDofPerNode,
                    map, theDiag, amalgDropTol, amalgRowPtr, amalgCols,
                    epetraFramework, myLocalNodeIds);

   theDiag.resize(0);   // free space

   // remove (i,j) connnections in amalgamated matrix associated with nodes
   // where the number of dofs per node is different

   MLrmDifferentDofsCrossings(dofPresent,maxDofPerNode,amalgRowPtr,amalgCols,nDofs+nGhost, epetraFramework, myLocalNodeIds);

   int iiii = amalgColMap.size();

   MLVec<double> ghostedXXX(iiii);
   if (YYY == NULL) iiii = 0;
   MLVec<double> ghostedYYY(iiii);
   if (ZZZ == NULL) iiii = 0;
   MLVec<double> ghostedZZZ(iiii);

   for (int i = 0; i < nNodes; i++) ghostedXXX[i] = XXX[i];
   nodalComm(ghostedXXX, myLocalNodeIds, epetraFramework);

   int localDim=1, globalDim;

   if (YYY != NULL) localDim = 2;
   if (ZZZ != NULL) localDim = 3;
   RowMatrix.Comm().MaxAll(&localDim,&globalDim,1);

   if (globalDim > 1) {
      for (int i = 0; i < nNodes; i++) ghostedYYY[i] = YYY[i];
      nodalComm(ghostedYYY, myLocalNodeIds, epetraFramework);
   }
   if (globalDim > 2) {
      for (int i = 0; i < nNodes; i++) ghostedZZZ[i] = ZZZ[i];
      nodalComm(ghostedZZZ, myLocalNodeIds, epetraFramework);
   }

   MLVec<double> lapVals(amalgRowPtr[nNodes]);

   MLbuildLaplacian(amalgRowPtr, amalgCols, lapVals, ghostedXXX,ghostedYYY,ghostedZZZ);

   MLsortCols(amalgRowPtr, amalgCols, lapVals);

   MLVec<char> status(nNodes*maxDofPerNode);
   MLfineStatus(dofPresent, map, dirOrNot, status);

   dirOrNot.resize(0);      // free space

   map.resize(0);           // free space
   Epetra_Map* epetAmalgRowMap=new Epetra_Map(-1,nNodes, amalgRowMap.getptr(), 0, RowMatrix.Comm());
   Epetra_Map* epetAmalgColMap=new Epetra_Map(-1,amalgColMap.size(),amalgColMap.getptr(),0,RowMatrix.Comm());

   // Need to take CRS vectors and create an epetra matrix

   // If we don't want to supply a column map, things should work by
   // uncommenting all toggle lines below and commenting out the
   // corresponding colmap lines.

   Epetra_CrsMatrix* LapMat = new Epetra_CrsMatrix(Copy,*epetAmalgRowMap,*epetAmalgColMap,0);
// Epetra_CrsMatrix* LapMat = new Epetra_CrsMatrix(Copy,*epetAmalgRowMap,0); // toggle

   for (int i = 0 ; i < nNodes ; ++i) {
     LapMat->InsertMyValues(i,amalgRowPtr[i+1]-amalgRowPtr[i],
                     &(lapVals[amalgRowPtr[i]]), &(amalgCols[amalgRowPtr[i]]));
     /*    toggle
     for (int j = amalgRowPtr[i]; j < amalgRowPtr[i+1]; j++) {
       LapMat->InsertGlobalValues(amalgRowMap[i],1, &(lapVals[j]),&(amalgColMap[amalgCols[j]]));
     } */
   }
   amalgRowPtr.resize(0); amalgCols.resize(0); lapVals.resize(0); // free space

   LapMat->FillComplete(*epetAmalgRowMap,*epetAmalgRowMap);
// LapMat->FillComplete();      // toggle
// LapMat->OptimizeStorage();   // toggle

#ifdef outputter
// this stuff is just for forcing ML do do the same thing on one and many procs
// the code is hidden behind an ifdef in ML_agg_MIS.c. You can use run in parallel
// turning on ML_AGGR_OUTAGGR. This will output a bunch of files (one per proc per
// level). You need to cat together files associated with each level (e.g.,
// cat agg*_0 > agg_0 ; cat agg*_1 > agg_1 ; ....). Then, recompile with
// ML_AGGR_OUTAGGR undefined but ML_AGGR_INAGGR defined. Not run in serial.
// If the smoother is jacobi, and the max eigenvalue is fixed independent of the
// number of procs (e.g., hacking in return 2.0; inside ML_Krylov_Get_MaxEigenvalue()),
// then the output should be the same as the parallel run

// if you take play with this, don't forget to ifdef out the extern int above
   LapMat->ColMap().MyGlobalElementsPtr(update);
update_index = (int *) ML_allocate(sizeof(int)*(nDofs+nGhost));
extern_index = (int *) ML_allocate(sizeof(int)*(nDofs+nGhost));
for (int i = 0; i < nNodes; i++) update_index[i] = i;
for (int i = 0; i < nGhost; i++) extern_index[i] = i+nNodes;
external = &(update[nNodes]);
#endif

  RowMatrix_ = LapMat;
  AllocatedRowMatrix_ = false;

  ML_CHK_ERRV(Initialize());

  DontSetSmoothers_ = true;   // don't make smoothers based on Laplacian
                              // smoothers on the real matrix will be done below
  // construct hierarchy
  if (ComputePrec == true)
    ML_CHK_ERRV(ComputePreconditioner());

  DontSetSmoothers_ = false;

  ML *altml;
  altml = ml_;

Comm_ = &(RowMatrix.Comm());
DomainMap_ = &(RowMatrix.OperatorDomainMap());
RangeMap_ = &(RowMatrix.OperatorRangeMap());


  int levelIndex, levelIncr = -1;

  levelIndex = ml_->ML_finest_level;
  if (levelIndex == 0) levelIncr = 1;

  // shove the real matrix into ML

  int numMyRows = Acrs->NumMyRows();
  int N_ghost   = Acrs->NumMyCols() - numMyRows;

  ML_Operator_Clean(&(altml->Amat[levelIndex]));
  ML_Operator_Init(&(altml->Amat[levelIndex]), altml->comm);
  delete LapMat;
  delete epetAmalgColMap;
  delete epetAmalgRowMap;
  ML_Init_Amatrix(altml,levelIndex,numMyRows, numMyRows, (void *) Acrs);
  altml->Amat[levelIndex].type = ML_TYPE_CRS_MATRIX;
  ML_Set_Amatrix_Getrow(altml, levelIndex, ML_Epetra_CrsMatrix_getrow,
                        ML_Epetra_CrsMatrix_comm_wrapper,
                        numMyRows+N_ghost);

  ML_Set_Amatrix_Matvec(altml, levelIndex, ML_Epetra_CrsMatrix_matvec);

  int nNodesCoarse;
  int nDofsCoarse;
  int nNodesFine = nNodes;

  for (int ii = 0; ii < ml_->ML_num_actual_levels - 1 ; ii++) {
    struct ML_CSR_MSRdata* ptr= (struct ML_CSR_MSRdata *) ml_->Pmat[levelIndex+levelIncr].data;

     // wrap crs stuff into std vectors

    nNodesCoarse = ml_->Pmat[levelIndex+levelIncr].invec_leng;

    MLVec<int>     PAmalgRowPtr(ptr->rowptr,ptr->rowptr+nNodesFine+1);
    MLVec<int>     PAmalgCols(ptr->columns,ptr->columns+(ptr->rowptr)[nNodesFine] );
    MLVec<double>  PAmalgVals(ptr->values,ptr->values + (ptr->rowptr)[nNodesFine] );
    MLVec<int>     newPRowPtr(nDofs+1);

    MLVec<int>     newPCols(PAmalgRowPtr[nNodesFine]*maxDofPerNode);
    MLVec<double>  newPVals(PAmalgRowPtr[nNodesFine]*maxDofPerNode);

    MLunamalgP(PAmalgRowPtr, PAmalgCols, PAmalgVals,
             maxDofPerNode, status, (ii != 0), newPRowPtr, newPCols, newPVals);

    int nGhostNew = 0;

    if ( (ml_->Pmat[levelIndex+levelIncr].getrow != NULL) &&
         (ml_->Pmat[levelIndex+levelIncr].getrow->pre_comm != NULL) )
      nGhostNew = ML_CommInfoOP_Compute_TotalRcvLength(ml_->Pmat[levelIndex+levelIncr].getrow->pre_comm);


    struct wrappedCommStruct mlFramework;
    mlFramework.whichOne = mlType;
    mlFramework.data     = (void *) &(ml_->Pmat[levelIndex+levelIncr]);
    mlFramework.nProcs   = RowMatrix.Comm().NumProc();
    mlFramework.myPid    = RowMatrix.Comm().MyPID();
    mlFramework.maxDofPerNode = maxDofPerNode;
    mlFramework.vecSize    = maxDofPerNode*(nGhostNew+ml_->Pmat[levelIndex+levelIncr].invec_leng);

    MLShove(&(altml->Pmat[levelIndex+levelIncr]),newPRowPtr,newPCols,newPVals,
            ml_->Pmat[levelIndex+levelIncr].invec_leng*maxDofPerNode,
            dofCommUsingMlNodalMatrix, mlFramework, maxDofPerNode*nGhostNew);

    ML_Operator_Set_1Levels(&(altml->Pmat[levelIndex+levelIncr]),
                           &(altml->SingleLevel[levelIndex+levelIncr]),
                           &(altml->SingleLevel[levelIndex]));

    // take the transpose and push in Rmat

    ML_Operator_Clean(&(altml->Rmat[levelIndex]));
    ML_Operator_Init(&(altml->Rmat[levelIndex]), altml->comm);
    ML_Gen_Restrictor_TransP(altml, levelIndex, levelIndex+levelIncr, NULL);
    ML_Operator_Clean(&(altml->Amat[levelIndex+levelIncr]));
    ML_Operator_Init(&(altml->Amat[levelIndex+levelIncr]), altml->comm);
    ML_Gen_AmatrixRAP(altml, levelIndex, levelIndex+levelIncr);
    nDofsCoarse = altml->Amat[levelIndex+levelIncr].outvec_leng;

    struct ML_CSR_MSRdata* newptr= (struct ML_CSR_MSRdata *) altml->Amat[levelIndex+levelIncr].data;

    int nCoarse = altml->Amat[levelIndex+levelIncr].outvec_leng;

    // the coarse level matrix might have some empty rows due to the
    // padding intended to force all coarse matrices to have the same
    // number of dofs per node. For these empty rows, we want to
    // stick in some Dirichlet points. It turns out that this would
    // be real easy with ML's MSR matrices. We would just looks for
    // zeros on the diagonal and replace them by ones. However, I'd
    // like to prepare the code for MueLu and CRS matrices. So ...
    // I instead copy the MSR data to a buncy of CRS arrays. I then
    // invoke MLfindEmptyRows() & MLreplaceEmptyByDir() that are intended
    // for CRS matrices. Finally, I need to copy these back to the
    // MSR that ML wants.

    // For a code like MueLu that works with CRS matrices, we could
    // perhaps just wrap the CRS vectors into MLVec's instead of copy as is
    // done for ML. We'd have to be careful that there is actually
    // enough storage for the new 1's entries.
    //
    //  MLVec<int>    AcoarseRowPtr(newptr->rowptr, newptr->rowptr +  nCoarse+1);
    //  MLVec<int>    AcoarseCols(  newptr->columns,newptr->columns + (newptr->rowptr)[nCoarse] );
    //  MLVec<double> AcoarseVals(  newptr->values, newptr->values  + (newptr->rowptr)[nCoarse] );
    //
    //

    // count the number of nonzeros (for ML)

    int count = 0;
    for (int i = 0; i < (newptr->columns)[nCoarse]; i++)
      if (newptr->values[i] != 0.0) count++;

    MLVec<int>      AcoarseRowPtr( nCoarse+1);
    MLVec<int>     AcoarseCols(  count );
    MLVec<double>  AcoarseVals(  count );

    // grab MSR stuff and put it into CRS style arrays

    count = 0;
    for (int i = 0; i < nCoarse; i++) {
       AcoarseRowPtr[i] = count;
       if (newptr->values[i] != 0.0) {
          AcoarseCols[count  ] = i;
          AcoarseVals[count++] = newptr->values[i];
       }
       for (int j = (newptr->columns)[i]; j < (newptr->columns)[i+1]; j++) {
          if (newptr->values[j] != 0.0) {
             AcoarseCols[count  ] = newptr->columns[j];
             AcoarseVals[count++] = newptr->values[j];
          }
       }
    }
    AcoarseRowPtr[nCoarse] = count;

    MLVec<bool> rowEmpty;
    MLfindEmptyRows(AcoarseRowPtr,AcoarseCols, AcoarseVals, rowEmpty);
    MLreplaceEmptyByDirichlet(AcoarseRowPtr,AcoarseCols,AcoarseVals, rowEmpty);

    // grabs CRS arrays and put results back into MSR format
    count = (newptr->columns)[0];
    for (int i = 0; i < nCoarse; i++) {
       for (int j = AcoarseRowPtr[i]; j < AcoarseRowPtr[i+1]; j++) {
          if ( AcoarseCols[j] == i)   (newptr->values)[i] = AcoarseVals[j];
          else {
             (newptr->columns)[count ] =  AcoarseCols[j];
             (newptr->values)[count++] =  AcoarseVals[j];
          }
       }
       (newptr->columns)[i+1] =  count;
    }
    AcoarseRowPtr.resize(0); AcoarseCols.resize(0); AcoarseVals.resize(0);

    // compute the status array for the next level indicating whether dofs
    // are padded or real.

    MLVec<bool>   empty;
    if (ii != ml_->ML_num_actual_levels - 2)
      MLcoarseStatus(rowEmpty, empty , status);

    rowEmpty.resize(0);   // MLVec.resize(0) actually frees space

    levelIndex = levelIndex + levelIncr;
    nNodesFine = nNodesCoarse;
    nDofs      = nDofsCoarse;

  }

  SetSmoothers();
  CreateLabel();



}
#endif
#ifdef NewStuff
// ================================================ ====== ==== ==== == =
/*! Constructor for scalar PDE problems based on applying AMG to the distance
 *  Laplacian operator when constructing grid transfers. The main unique
 *  feature is that there may be some dofs that correspond to the same node
 *  location. These shared dofs fall into two categories. If these dofs are
 *  strongly connected to each other (as determined by tol), they are
 *  explicitly elminated from the Laplacian (merged into a supernode). Once
 *  a P is obtained, this P is then expanded to account for shared nodes
 *  by simply duplicating the supernodes row of P for each of the individual
 *  vertices that contribute to the supernode. If share dofs are weakly
 *  connected (or not connected at all), nothing special is done (other than
 *  the ususal ignoring of weak connections). One last trick is employed,
 *  connections between supernodes and non-supernodes (i.e., regular nodes)
 *  are always assumed to be weak. Shared nodes are often used to capture
 *  interfaces or other features. By breaking these connections, the AMG
 *  can better maintain these features throughout the hierarchy. Finally, the
 *  constructor also allows for *  the removal of column nonzeros associated
 *  with Dirichlet points. To use this option the rhs and initial guess must be
 *  provided. Modification of the matrix, rhs, and initial guess must be
 *  allowable to use this option.
 */
ML_Epetra::MultiLevelPreconditioner::
MultiLevelPreconditioner(Epetra_RowMatrix & RowMatrix,
    const Teuchos::ParameterList & List,
    const double distTol, // two points are at the same location when
                          // || (x_1,y_1,z_1) -  (x_2,y_2,z_2)||_2 < distTol
    const double tol,     // ignore values when
                          //       A(i,j)^2 < A(i,i)*A(j,j)*tol^2
    Epetra_Vector & Lhs,
    Epetra_Vector & Rhs,
    const bool  rhsAndsolProvided,
    const bool ComputePrec) :

  RowMatrixAllocated_(0),
  AllocatedRowMatrix_(false)
{


  Epetra_CrsMatrix *Acrs =
              dynamic_cast<Epetra_CrsMatrix *>(&RowMatrix);

  TEUCHOS_TEST_FOR_EXCEPT_MSG(Acrs == NULL,
             ErrorMsg_ << "matrix must be crs to use this constructor\n");

  int nDofs = Acrs->NumMyRows();
  int *dofGlobals;

  Acrs->RowMap().MyGlobalElementsPtr(dofGlobals);

  List_ = List;

  double *XXX = NULL, *YYY = NULL, *ZZZ = NULL;

  XXX = List_.get("x-coordinates",(double *) 0);
  YYY = List_.get("y-coordinates",(double *) 0);
  ZZZ = List_.get("z-coordinates",(double *) 0);

  TEUCHOS_TEST_FOR_EXCEPT_MSG(XXX == NULL,
             ErrorMsg_ << "Must supply coordinates to use shaerd dof constructor\n");

  struct wrappedCommStruct epetraFramework;

  epetraFramework.whichOne = epetraType;
  epetraFramework.data     = (void *) Acrs;
  epetraFramework.nProcs   = RowMatrix.Comm().NumProc();
  epetraFramework.myPid    = RowMatrix.Comm().MyPID();
  int nGhost = MLnMyGhost(epetraFramework);

  // allocate vectors

  MLVec<bool>   dirOrNot(nDofs + nGhost);
  MLVec<double> theDiag(nDofs + nGhost);
  MLVec<int>    map(nDofs);
  MLVec<int>    amalgRowPtr, amalgCols;
  MLVec<int>    amalgRowMap;
  MLVec<int>    amalgColMap;

  // Make MLVec versions of Crs arrays by wrapping them.

  double *vtemp;  int *ctemp,*rtemp;
  Acrs->ExtractCrsDataPointers(rtemp, ctemp, vtemp );

  MLVec<int>    rowPtr(rtemp,rtemp + nDofs + 1);
  MLVec<int>    cols(ctemp,ctemp   + rtemp[nDofs]);
  MLVec<double> vals(vtemp,vtemp   + rtemp[nDofs]);

  // grab diagonal and detect Dirichlets
  //
  double dirDropTol = List_.get("ML advanced Dirichlet: threshold",1.0e-5);
  MLextractDiag(rowPtr,   cols, vals, theDiag, epetraFramework);
  MLfindDirichlets(rowPtr,cols, vals, theDiag, dirDropTol, dirOrNot, epetraFramework);

  // if the rhs and solution are provided, remove column entries
  // associated with Dirichlet rows. Here, we overwrite the original
  // matrix arrays to reflect removal.

  if (rhsAndsolProvided) {

     double *ptrRHS, *ptrLHS;

     Lhs.ExtractView(&ptrLHS);
     Rhs.ExtractView(&ptrRHS);

     MLVec<double> myvecRHS(ptrRHS, ptrRHS + nDofs );
     MLVec<double> myvecLHS(ptrLHS, ptrLHS + nDofs );

     MLrmDirichletCols(rowPtr, cols, vals, theDiag, true,
                     myvecLHS, myvecRHS, dirOrNot, epetraFramework);
   }

    MLVec<int>    rowZeros(nDofs+nGhost);
    MLVec<int>    colZeros(nDofs+nGhost);

    int *rowmapIds, *colmapIds;
    Acrs->RowMap().MyGlobalElementsPtr(rowmapIds);
    Acrs->ColMap().MyGlobalElementsPtr(colmapIds);
    MLVec<int>    myLocalNodeIds(nDofs+nGhost);

    for (int i = 0; i < nDofs; i++) myLocalNodeIds[i] = i;

    int nLocalNodes, nLocalPlusGhostNodes;

    MLassignGhostLocalNodeIds( myLocalNodeIds, nDofs, nDofs+nGhost,
        epetraFramework, nLocalNodes, nLocalPlusGhostNodes);

    int iiii = Acrs->ColMap().NumMyElements();

    MLVec<double> xx(iiii);   if (YYY == NULL) iiii = 0;
    MLVec<double> yy(iiii);   if (ZZZ == NULL) iiii = 0;
    MLVec<double> zz(iiii);

    for (int i = 0; i < nDofs; i++) xx[i] = XXX[i];
    nodalComm(xx, myLocalNodeIds, epetraFramework);

    if (YYY != NULL) {
      for (int i = 0; i < nDofs; i++) yy[i] = YYY[i];
      nodalComm(yy, myLocalNodeIds, epetraFramework);
    }
    if (ZZZ != NULL) {
      for (int i = 0; i < nDofs; i++) zz[i] = ZZZ[i];
      nodalComm(zz, myLocalNodeIds, epetraFramework);
    }

    MLVec<int> newGlobals(iiii);

    ZeroDist(xx,yy,zz,rowPtr,cols,vals,theDiag,tol,rowZeros,colZeros,distTol);

    MLVec<int> groupHead(nDofs+nGhost);
    MLVec<int> groupNext(nDofs+nGhost);

    MergeShared(cols, rowZeros, colZeros, groupHead, groupNext);
    for (int i = 0; i < nDofs; i++) newGlobals[i] = rowmapIds[i];
    // shared dofs will have some global ids
    for (int i = 0; i < nDofs; i++) {
       if (groupHead[i] > -1) newGlobals[i] = newGlobals[groupHead[i]];
    }
    nodalComm(newGlobals, myLocalNodeIds, epetraFramework);

    MLVec<int> newMap(nDofs);
    int        nNonShared;

    nNonShared = BuildNonSharedMap(newMap, groupHead, groupNext);

    // fill the ghost part of newMap
    newMap.resize(nDofs+nGhost);
    nodalComm(groupHead, myLocalNodeIds, epetraFramework);

    int Ncols = nNonShared;
    for (int i = nDofs; i < newMap.size(); i++) {
       if  (groupHead[i] == -1) {
          newMap[i] = Ncols;
          Ncols++;
       }
       else if (groupHead[i] >= 0) {// first time group encountered
          newMap[i] = Ncols;
          for (int j = i+1; j < newMap.size(); j++) {
             if ( newGlobals[j] == newGlobals[i]) {
                newMap[j] = Ncols;
                groupHead[j] = -2;
             }
          }
          Ncols++;
       }
    }

    MLVec<int> newRowPtr(nNonShared+1);
    MLVec<int> newCols(cols.size());

    buildCompressedA(rowPtr, cols, vals, theDiag, groupHead, groupNext,
                     newRowPtr, newCols, newMap, nNonShared, tol);

    MLVec<double> newVals(newCols.size());

    MLVec<double> xxCompressed(Ncols);
    MLVec<double> yyCompressed(Ncols);
    MLVec<double> zzCompressed(Ncols);
    MLVec<int>    newGlobsCompressed(Ncols);

    for (int i = 0; i < newMap.size(); i++)  {
       xxCompressed[newMap[i]] = xx[i];
       yyCompressed[newMap[i]] = yy[i];
       zzCompressed[newMap[i]] = zz[i];
       newGlobsCompressed[newMap[i]] = newGlobals[i];
    }

    MLbuildLaplacian(newRowPtr, newCols, newVals, xxCompressed, yyCompressed, zzCompressed);

   MLsortCols(newRowPtr, newCols, newVals);


   dirOrNot.resize(0);      // free space

   Epetra_Map* epetNewRowMap=new Epetra_Map(-1,nNonShared,newGlobsCompressed.getptr(),0,RowMatrix.Comm());
   Epetra_Map* epetNewColMap=new Epetra_Map(-1,Ncols,newGlobsCompressed.getptr(),0,RowMatrix.Comm());

   // Need to take CRS vectors and create an epetra matrix

   Epetra_CrsMatrix* LapMat = new Epetra_CrsMatrix(Copy,*epetNewRowMap,*epetNewColMap,0);
   for (int i = 0 ; i < nNonShared; ++i) {
     LapMat->InsertMyValues(i,newRowPtr[i+1]-newRowPtr[i],
                     &(newVals[newRowPtr[i]]), &(newCols[newRowPtr[i]]));
   }
   newRowPtr.resize(0); newCols.resize(0); newVals.resize(0); // free space

   LapMat->FillComplete(*epetNewRowMap,*epetNewRowMap);

#ifdef outputter
   LapMat->ColMap().MyGlobalElementsPtr(update);
update_index = (int *) ML_allocate(sizeof(int)*(nNonShared));
extern_index = (int *) ML_allocate(sizeof(int)*(nNonShared));
for (int i = 0; i < nNonShared; i++) update_index[i] = i;
for (int i = 0; i < nGhost; i++) extern_index[i] = i+nNonShared;
external = &(update[nNonShared]);
#endif

  RowMatrix_ = LapMat;
  AllocatedRowMatrix_ = false;

  ML_CHK_ERRV(Initialize());

  DontSetSmoothers_ = true;   // don't make smoothers based on Laplacian
                              // smoothers on the real matrix will be done below
  // construct hierarchy
  if (ComputePrec == true)
    ML_CHK_ERRV(ComputePreconditioner());

  DontSetSmoothers_ = false;

  ML *altml;
  altml = ml_;

Comm_ = &(RowMatrix.Comm());
DomainMap_ = &(RowMatrix.OperatorDomainMap());
RangeMap_ = &(RowMatrix.OperatorRangeMap());


  int levelIndex, levelIncr = -1;

  levelIndex = ml_->ML_finest_level;
  if (levelIndex == 0) levelIncr = 1;

  // shove actual discretization matrix into ML

  int numMyRows = Acrs->NumMyRows();
  int N_ghost   = Acrs->NumMyCols() - numMyRows;

  ML_Operator_Clean(&(altml->Amat[levelIndex]));
  ML_Operator_Init(&(altml->Amat[levelIndex]), altml->comm);
  delete LapMat;
  delete epetNewColMap;
  delete epetNewRowMap;
  ML_Init_Amatrix(altml,levelIndex,numMyRows, numMyRows, (void *) Acrs);
  altml->Amat[levelIndex].type = ML_TYPE_CRS_MATRIX;
  ML_Set_Amatrix_Getrow(altml, levelIndex, ML_Epetra_CrsMatrix_getrow,
                        ML_Epetra_CrsMatrix_comm_wrapper,
                        numMyRows+N_ghost);
  ML_Set_Amatrix_Matvec(altml, levelIndex, ML_Epetra_CrsMatrix_matvec);

  /* uncompress nonshared prolongator to reflect shared dofs */

  struct ML_CSR_MSRdata* ptr= (struct ML_CSR_MSRdata *) ml_->Pmat[levelIndex+levelIncr].data;
  double               *Pvals  = ptr->values;
  int                  *Pcols  = ptr->columns;
  int                  *Prowptr= ptr->rowptr;

  // count the number of nonzeros in uncompressed P

  int nnz = 0;
  for (int i = 0; i < nDofs; i++) {
     int j = newMap[i];
     nnz += ( Prowptr[j+1] - Prowptr[j]);
  }

  // allocated new CSR arrays and fill with uncompressed data

  MLVec<int> newPRowPtr(nDofs+1);
  MLVec<int>     newPCols(nnz);
  MLVec<double>  newPVals(nnz);

  nnz = 0;
  newPRowPtr[0] = 0;
  for (int i = 0; i < nDofs; i++) {
    for (int j= Prowptr[newMap[i]]; j < Prowptr[newMap[i]+1]; j++){
      newPCols[nnz] = Pcols[j];
      newPVals[nnz] = Pvals[j];
      nnz++;
    }
    newPRowPtr[i+1] = nnz;
  }

  // stick new CSR arrays into the hierarchy

  struct ML_CSR_MSRdata* Newptr = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
  if (Newptr == NULL) pr_error("no space for NewPtr\n");

  Newptr->rowptr  = newPRowPtr.getptr();     // getptr() + relinquishData
  Newptr->columns = newPCols.getptr();       // effectively pulls data
  Newptr->values  = newPVals.getptr();       // out of a MLVec and empties
  newPRowPtr.relinquishData();               // MLVec's contents
  newPCols.relinquishData();
  newPVals.relinquishData();


  // old communication widget should still be good as uncompressed matrix has
  // the same ghost information. So, we'll take it out of the hierarchy, then
  // set it to null within the hierarchy so that ML_Operator_Clean() does not
  // free it. Then later will manuall assign the save point to the new
  // uncompressed matrix within the hierarchy

  int insize  = ml_->Pmat[levelIndex+levelIncr].invec_leng;
  ML_CommInfoOP *newPrecComm = altml->Pmat[levelIndex+levelIncr].getrow->pre_comm;
  altml->Pmat[levelIndex+levelIncr].getrow->pre_comm = NULL;
  ML_Operator_Clean(&(altml->Pmat[levelIndex+levelIncr]));
  ML_Operator_Init(&(altml->Pmat[levelIndex+levelIncr]), altml->comm);
  ML_Operator_Set_ApplyFuncData(&(altml->Pmat[levelIndex+levelIncr]), insize,
                                nDofs, Newptr, nDofs, NULL,0);
  ML_Operator_Set_Getrow(&(altml->Pmat[levelIndex+levelIncr]),nDofs,CSR_getrow);
  ML_Operator_Set_ApplyFunc(&(altml->Pmat[levelIndex+levelIncr]), CSR_matvec);
  altml->Pmat[levelIndex+levelIncr].data_destroy     = ML_CSR_MSRdata_Destroy;
  altml->Pmat[levelIndex+levelIncr].N_nonzeros       = nnz;
  altml->Pmat[levelIndex+levelIncr].getrow->pre_comm = newPrecComm;

  ML_Operator_Set_1Levels(&(altml->Pmat[levelIndex+levelIncr]),
                          &(altml->SingleLevel[levelIndex+levelIncr]),
                          &(altml->SingleLevel[levelIndex]));

  // take the transpose and push in Rmat

  ML_Operator_Clean(&(altml->Rmat[levelIndex]));
  ML_Operator_Init(&(altml->Rmat[levelIndex]), altml->comm);
  ML_Gen_Restrictor_TransP(altml, levelIndex, levelIndex+levelIncr, NULL);

  // re-do all RAPs with the real A (as opposesd to Laplacian)
  for (int ii = 0; ii < ml_->ML_num_actual_levels - 1 ; ii++) {

    ML_Operator_Clean(&(altml->Amat[levelIndex+levelIncr]));
    ML_Operator_Init(&(altml->Amat[levelIndex+levelIncr]), altml->comm);
    ML_Gen_AmatrixRAP(altml, levelIndex, levelIndex+levelIncr);

    levelIndex = levelIndex + levelIncr;
  }

  SetSmoothers();
  CreateLabel();

}
#endif

// ================================================ ====== ==== ==== == =

#ifdef HAVE_ML_AZTECOO
/*! Another constructor for Maxwell equations that takes an Epetra_MsrMatrix,
 * an ML_Operator, and an AZ_MATRIX type.  The Epetra_MsrMatrix type is defined
 * in aztecoo.
 * Two conditions are required on their maps:
 * - TMatrix.OperatorDomainMap() == NodeMatrix.OperatorRangeMap()
 * - TMatrix.OperatorRangeMap()  == EdgeMatrix.OperatorDomainMap()
 */
ML_Epetra::MultiLevelPreconditioner::
MultiLevelPreconditioner(const Epetra_MsrMatrix & EdgeMatrix,
                         ML_Operator * ML_TMatrix,
                         AZ_MATRIX * AZ_NodeMatrix,
                         int       * proc_config,
                         const ParameterList & List,
                         const bool ComputePrec):
  RowMatrix_(&EdgeMatrix),
  RowMatrixAllocated_(0),
  AllocatedRowMatrix_(false)
{
  ML_Comm *ml_comm;
  ML_Operator *loc_ML_Kn;

  int MaxNumNonzeros;
  double CPUTime;
  Epetra_CrsMatrix *NodeMatrix, *TMatrix;

  //next-to-last argument must be false because T may have zero rows due
  //to Dirichlet b.c.'s.
  ML_Operator2EpetraCrsMatrix(ML_TMatrix,TMatrix,MaxNumNonzeros,
                  false,CPUTime);

  ML_Comm_Create(&ml_comm);
#ifdef ML_MPI
  const Epetra_MpiComm *epcomm = dynamic_cast<const Epetra_MpiComm*>(&(EdgeMatrix.Comm()));
  // Get the MPI communicator, as it may not be MPI_COMM_W0RLD, and update the ML comm object
  if (epcomm) ML_Comm_Set_UsrComm(ml_comm,epcomm->Comm());
#endif
  loc_ML_Kn = ML_Operator_Create(ml_comm);
  AZ_convert_aztec_matrix_2ml_matrix(AZ_NodeMatrix,loc_ML_Kn,proc_config);
  ML_Operator2EpetraCrsMatrix(loc_ML_Kn, NodeMatrix,MaxNumNonzeros,
                  false,CPUTime);

  // some sanity checks
  if (! TMatrix->OperatorDomainMap().SameAs(NodeMatrix->OperatorRangeMap()) ) {
    std::cerr << ErrorMsg_ << "discrete grad DomainMap != node RangeMap..." << std::endl;
    ML_CHK_ERRV(-1); // error on discrete grad
  }

  if (! TMatrix->OperatorRangeMap().SameAs(EdgeMatrix.OperatorDomainMap()) ) {
    std::cerr << ErrorMsg_ << "discrete grad RangeMap != edge DomainMap..." <<std::endl;
    ML_CHK_ERRV(-2); // error on discrete grad
  }

  List_ = List;

  ML_CHK_ERRV(Initialize());

  // set Maxwell here.
  // NOTE: RowMatrix_ and EdgeMatrix_ pointer to the same Epetra_RowMatrix
  AMGSolver_ = ML_MAXWELL;
  NodeMatrix_ = NodeMatrix;
  CreatedNodeMatrix_ = true;
  TMatrix_ = TMatrix;
  CreatedTMatrix_ = true;
  EdgeMatrix_ = & EdgeMatrix;
  ML_Kn_ = loc_ML_Kn;
  CreatedML_Kn_ = true;

  // construct hierarchy
  if (ComputePrec == true)
    ML_CHK_ERRV(ComputePreconditioner());

  ML_Comm_Destroy(&ml_comm);
}
#endif

// ================================================ ====== ==== ==== == =
// FIXME: should I be deleted??
ML_Epetra::MultiLevelPreconditioner::
MultiLevelPreconditioner(ML_Operator * Operator,
             const ParameterList & List, const bool ComputePrec):
  AllocatedRowMatrix_(false)
{

  // need to wrap an Epetra_RowMatrix around Operator.
  // This is quite not the best approach for small matrices

  int MaxNumNonzeros;
  double CPUTime;

  // FIXME: substitute with a better version
  ML_Operator2EpetraCrsMatrix(Operator,RowMatrixAllocated_,MaxNumNonzeros,
                  true,CPUTime);

  // this matrix must be freed by dtor. Keep trace of it in this pointer
  RowMatrix_ = RowMatrixAllocated_;

  // from now on as for the other constructors
  List_ = List;

  ML_CHK_ERRV(Initialize());

  // construct hierarchy
  if (ComputePrec == true)
    ML_CHK_ERRV(ComputePreconditioner());
}

// ================================================ ====== ==== ==== == =
/*! - set to 0 all allocatable pointers
 *  - put default values in Aztec vectors.
 *  - zero-out timing
 */
int ML_Epetra::MultiLevelPreconditioner::Initialize()
{

  Comm_ = &(RowMatrix_->Comm());
  DomainMap_ = &(RowMatrix_->OperatorDomainMap());
  RangeMap_ = &(RowMatrix_->OperatorRangeMap());
  AfineML_ = 0;
  SubMatMLPrec_ = 0;

#ifdef NewStuff
  DontSetSmoothers_ = false;
#endif
  verbose_ = false;
  MaxLevels_ = 20;
  IsComputePreconditionerOK_ = false;

  NumPDEEqns_ = 1;

  NullSpaceToFree_ = 0;

  Label_ = 0;

  ml_ = 0;
  agg_ = 0;

  sprintf(ErrorMsg_,"%s","*ML*ERR* : ");
  PrintMsg_ = "";

#ifdef HAVE_ML_AZTECOO
  SmootherOptions_ = Teuchos::rcp(new std::vector<int>(AZ_OPTIONS_SIZE));
  SmootherParams_  = Teuchos::rcp(new std::vector<double>(AZ_PARAMS_SIZE));
  AZ_defaults(&(*SmootherOptions_)[0],&(*SmootherParams_)[0]);
  (*SmootherOptions_)[AZ_precond] = AZ_dom_decomp;
  (*SmootherOptions_)[AZ_subdomain_solve] = AZ_ilu;
  (*SmootherOptions_)[AZ_overlap] = 0;
#endif

  // By default we assume that a smoothed aggregation-like method is being used.
  AMGSolver_ = ML_SA_FAMILY;
  // Maxwell stuff is off by default
  NodeMatrix_ = 0;
  CreatedNodeMatrix_ = false;
  ML_Kn_ = 0;
  CreatedML_Kn_ = false;
  EdgeMatrix_ = 0;
  CurlCurlMatrix_ = 0;
  CreatedEdgeMatrix_ = false;
  MassMatrix_ = 0;
  TMatrix_ = 0;
  TtATMatrixML_=0;
  UseNodeMatrixForSmoother_=false;
  ml_nodes_ = 0;
  TMatrixML_ = 0;
  CreatedTMatrix_ = false;
  TMatrixTransposeML_ = 0;
  Tmat_array = 0;
  Tmat_trans_array = 0;
  nodal_args_ = 0;
  edge_args_ = 0;
  MassMatrix_array = 0;
  CurlCurlMatrix_array = 0;

  // timing
  NumApplications_ = 0;
  ApplicationTime_ = 0.0;
  FirstApplication_ = true;
  FirstApplicationTime_ = 0.0;
  NumConstructions_ = 0;
  ConstructionTime_ = 0.0;

  // some tracking here
  int NumInitializations = OutputList_.get("number of initialization phases", 0);
  OutputList_.set("number of initialization phases", ++NumInitializations);

  // memory
  AnalyzeMemory_ = false;

  for (int i = 0 ; i < ML_MEM_SIZE ; ++i)
    memory_[i] = 0;

  // filtering vectors
  flt_ml_ = 0;
  flt_agg_ = 0;

  // CheckPreconditioner stuff
  RateOfConvergence_ = -1.0;

  return 0;
}
int ML_Epetra::MultiLevelPreconditioner::SetEigenScheme()
{
  // ================================================================== //
  // Define scheme for computing the spectral radius of the matrix.     //
  // ================================================================== //

  std::string str = List_.get("eigen-analysis: type","power-method");

  if( verbose_ ) std::cout << PrintMsg_ << "Using `" << str << "' for eigen-computations" << std::endl;

  if( str == "cg" )                ML_Set_SpectralNormScheme_Calc(ml_);
  else if( str == "Anorm" )        ML_Set_SpectralNormScheme_Anorm(ml_);
  else if( str == "Anasazi" )      ML_Set_SpectralNormScheme_Anasazi(ml_);
  else if( str == "power-method" ) ML_Set_SpectralNormScheme_PowerMethod(ml_);
  else {
    if( Comm().MyPID() == 0 ) {
      std::cerr << ErrorMsg_ << "parameter `eigen-analysis: type' has an incorrect value"
       << "(" << str << ")" << std::endl;
      std::cerr << ErrorMsg_ << "It should be: " << std::endl
       << ErrorMsg_ << "<cg> / <Anorm> / <Anasazi> / <power-method>" << std::endl;
    }
    ML_EXIT(-10); // wrong input parameter
  }
  int NumEigenIts = List_.get("eigen-analysis: iterations",10);
  ML_Set_SpectralNorm_Iterations(ml_, NumEigenIts);

   return 0;
}
int ML_Epetra::MultiLevelPreconditioner::SetLevelIds(int Direction)
{
  // level traversal (increasing of descreasing), default:  ML_INCREASING.

  int FinestLevel;
  ML_Set_LevelID(ml_,Direction);  // These are internal to ml_ . They are
                                  // also reset (to the same numbers)
                                  // when using smoothed aggregation.

  LevelID_.resize(NumLevels_);    // These owned by MultiLevelPreconditioner()
  if (Direction == ML_INCREASING) {
    FinestLevel = 0;
    for( int i=0 ; i<NumLevels_ ; ++i ) LevelID_[i] = FinestLevel+i;
  } else {
    FinestLevel = NumLevels_-1;
    for( int i=0 ; i<NumLevels_ ; ++i ) LevelID_[i] = FinestLevel-i;
  }

  // check no levels are negative
  for (int i = 0 ; i < NumLevels_ ; ++i)
    if (LevelID_[i] < 0) {
      std::cerr << ErrorMsg_ << "Level " << i << " has a negative ID" << std::endl;
      ML_EXIT(EXIT_FAILURE);
    }

   return 0;
}

void ML_Epetra::MultiLevelPreconditioner::Complexities(double &complexity, double &fineNnz)
{
  if ( agg_ != 0 && agg_->fine_complexity != 0.0 ) {
    complexity = agg_->operator_complexity / agg_->fine_complexity;
    fineNnz = agg_->fine_complexity;
  } else {
    complexity = -1.0;
    fineNnz = -1.0;
  }
}

int ML_Epetra::MultiLevelPreconditioner::ComputeAndPrintComplexities()
{
  // ====================================================================== //
  // Recompute complexities                                                 //
  // ====================================================================== //

  if (AMGSolver_ != ML_MAXWELL) {

    double RowComplexity = 0.0, RowZero = 0.0;
    double NnzComplexity = 0.0, NnzZero = 0.0;

    for (int i = 0 ; i < NumLevels_ ; ++i)
    {
      double local[2];
      double global[2];
      local[0] = ml_->Amat[LevelID_[i]].invec_leng;
      local[1] = ml_->Amat[LevelID_[i]].N_nonzeros;
      Comm().SumAll(local,global,2);
      if (i == 0) {
        RowZero = global[0];
        NnzZero = global[1];
      }
      RowComplexity += global[0];
      NnzComplexity += global[1];
    }
    if (verbose_)
    {
      std::cout << PrintMsg_ << std::endl;
      std::cout << PrintMsg_ << "sum n_i   / n_finest   = " << RowComplexity / RowZero << std::endl;
      std::cout << PrintMsg_ << "sum nnz_i / nnz_finest = " << NnzComplexity / NnzZero << std::endl;
    }
  }
  return 0;
}
int ML_Epetra::MultiLevelPreconditioner::MatrixDumper()
{

#ifdef HAVE_ML_EPETRAEXT
  // ============================================================== //
  // dump matrix in matrix market format using EpetraExt            //
  // ============================================================== //

  if (List_.get("dump matrix: enable", false)) {
    std::string FileName = List_.get("dump matrix: file name", "ml-matrix.mm");
    char FileName2[80];
    static int count = 0;
    char comment[80];
    sprintf(comment,"printed by MultiLevelPreconditioner, %d processors",
            Comm().NumProc());

    if (AMGSolver_ == ML_MAXWELL) {
      sprintf(FileName2,"ml-matrix-%s.mm","edge");
      EpetraExt::RowMatrixToMatrixMarketFile(FileName2, *RowMatrix_,
                              "Maxwell -- edge finite element matrix",
                              comment);
      sprintf(FileName2,"ml-matrix-%s.mm","curlcurl");
      if (CurlCurlMatrix_)
        EpetraExt::RowMatrixToMatrixMarketFile(FileName2, *CurlCurlMatrix_,
                                "Maxwell -- (curl,curl) matrix",
                                comment);
      sprintf(FileName2,"ml-matrix-%s.mm","mass");
      if (MassMatrix_)
        EpetraExt::RowMatrixToMatrixMarketFile(FileName2, *MassMatrix_,
                                "Maxwell -- mass matrix",
                                comment);
      sprintf(FileName2,"ml-matrix-%s.mm","gradient");
      EpetraExt::RowMatrixToMatrixMarketFile(FileName2, *TMatrix_,
                                "Maxwell -- gradient matrix",
                                comment);
      sprintf(FileName2,"ml-matrix-%s.mm","node");
      EpetraExt::RowMatrixToMatrixMarketFile(FileName2, *NodeMatrix_,
                                "Maxwell -- auxiliary nodal FE matrix",
                                comment);
    }
    else {
      if (FileName == "auto-count")
        sprintf(FileName2,"ml-matrix-%d.mm",count++);
      else
        sprintf(FileName2,"%s",FileName.c_str());

      if (verbose_)
        std::cout << PrintMsg_ << "Print matrix on file `" << FileName2
             << "'..." << std::endl;
      EpetraExt::RowMatrixToMatrixMarketFile(FileName2, *RowMatrix_,
                     "Epetra_RowMatrix", comment);
    }
  }
#endif
  // ====================================================================== //
  // One may decide to print out the entire hierarchy (to be analyzed in    //
  // MATLAB, for instance).                                                 //
  // ====================================================================== //
  int whichLevels = List_.get("print hierarchy", -2);
  if ( whichLevels > -2) Print(whichLevels);

  return 0;
}
int ML_Epetra::MultiLevelPreconditioner::SetFinestLevelMatrix()
{
  // avoid possible integer overflow in Epetra's global nnz count
  double localNnz = RowMatrix_->NumMyNonzeros();
  double globalNnz=0;
  Comm().SumAll(&localNnz,&globalNnz,1);

  if (verbose_) {
    std::cout << PrintMsg_ << "*** " << std::endl;
    std::cout << PrintMsg_ << "*** ML_Epetra::MultiLevelPreconditioner";
    if (mlpLabel_ != "not-set")
      std::cout << " [" << mlpLabel_ << "]";
    std::cout << std::endl << PrintMsg_ << "***" << std::endl;
    std::cout << PrintMsg_ << "Matrix has " << RowMatrix_->NumGlobalRows()
     << " rows and " << globalNnz
         << " nonzeros, distributed over " << Comm().NumProc() << " process(es)" << std::endl;
    {
      const Epetra_CrsMatrix * dummy = dynamic_cast<const Epetra_CrsMatrix *>(RowMatrix_);
      if( dummy ) std::cout << "The linear system matrix is an Epetra_CrsMatrix" << std::endl;
    }
    {
      const Epetra_VbrMatrix * dummy = dynamic_cast<const Epetra_VbrMatrix *>(RowMatrix_);
      if( dummy ) std::cout << "The linear system matrix is an Epetra_VbrMatrix" << std::endl;
    }
  }

  if (AMGSolver_ != ML_MAXWELL) {
#ifdef HAVE_ML_EPETRAEXT
    RowMatrix_ = ModifyEpetraMatrixColMap(*RowMatrix_,RowMatrixColMapTrans_,
                                          "Main linear system");
#endif
  }
  else {
    // if curl-curl and mass matrix were given separately, add them
    if (CurlCurlMatrix_) {
#ifdef HAVE_ML_EPETRAEXT
      CurlCurlMatrix_ = ModifyEpetraMatrixColMap(*CurlCurlMatrix_,
                                    CurlCurlMatrixColMapTrans_,"(curl,curl)");
      MassMatrix_ = ModifyEpetraMatrixColMap(*MassMatrix_,
                                             MassMatrixColMapTrans_,"mass");
#endif
      Epetra_RowMatrix *cc = const_cast<Epetra_RowMatrix*>(CurlCurlMatrix_);
      Epetra_RowMatrix *mm = const_cast<Epetra_RowMatrix*>(MassMatrix_);
      RowMatrix_ = Epetra_MatrixAdd(cc,mm,1.0);
      CreatedEdgeMatrix_ = true;
    }
#ifdef HAVE_ML_EPETRAEXT
    EdgeMatrix_ = ModifyEpetraMatrixColMap(*RowMatrix_,
                       RowMatrixColMapTrans_,"edge element matrix");
#endif
    RowMatrix_ = EdgeMatrix_;
  }
    int numMyRows;
    numMyRows = RowMatrix_->NumMyRows();
    int N_ghost = RowMatrix_->NumMyCols() - numMyRows;

    if (N_ghost < 0) N_ghost = 0;  // A->NumMyCols() = 0 for an empty matrix

    const Epetra_VbrMatrix *Avbr =
              dynamic_cast<const Epetra_VbrMatrix *>(RowMatrix_);
    const Epetra_CrsMatrix *Acrs =
              dynamic_cast<const Epetra_CrsMatrix *>(RowMatrix_);
    if (AMGSolver_ == ML_COMPOSITE) {
       struct MLBlkMat *widget = NULL;
       widget = (struct MLBlkMat *) AfineML_->data;

       ML_CommInfoOP_Clone(&(ml_->Amat[LevelID_[0]].getrow->pre_comm),
                           AfineML_->getrow->pre_comm);
       ML_Operator_Set_ApplyFuncData(&(ml_->Amat[LevelID_[0]]), widget->invec,
                                     widget->outvec, widget, widget->outvec,
                                     ML_Operator_BlkMatMatvec,0);
       ML_Operator_Set_Getrow(&(ml_->Amat[LevelID_[0]]), widget->outvec,
                              ML_Operator_BlkMatGetrow);
      ml_->Amat[LevelID_[0]].type = -666;
    }
    else if (Avbr) {
      ML_Init_Amatrix(ml_,LevelID_[0],numMyRows, numMyRows, (void *) Avbr);
      ml_->Amat[LevelID_[0]].type = ML_TYPE_VBR_MATRIX;
    }
    else if (Acrs) {
      ML_Init_Amatrix(ml_,LevelID_[0],numMyRows, numMyRows, (void *) Acrs);
      ml_->Amat[LevelID_[0]].type = ML_TYPE_CRS_MATRIX;
    }
    else {
      ML_Init_Amatrix(ml_,LevelID_[0],numMyRows, numMyRows,(void *) RowMatrix_);
      ml_->Amat[LevelID_[0]].type = ML_TYPE_ROW_MATRIX;
    }
    ml_->Amat[LevelID_[0]].N_nonzeros = RowMatrix_->NumMyNonzeros();

    ML_Epetra::FilterType FT;
    FT = List_.get("filter: type", ML_Epetra::ML_NO_FILTER);
    double* Mask = 0;
    Mask = List_.get("filter: mask", Mask);
    if (FT != ML_Epetra::ML_NO_FILTER) {
      List_.set("filter: equations", NumPDEEqns_);
      ML_Set_Filter(List_);
      ML_Set_Amatrix_Getrow(ml_, LevelID_[0], ML_Epetra_getrow_Filter,
                            ML_Epetra_comm_wrapper, numMyRows+N_ghost);
      ML_Set_Amatrix_Matvec(ml_, LevelID_[0], ML_Epetra_matvec_Filter);
    }
    else {
      if (AMGSolver_ == ML_COMPOSITE) {
        ML_Set_Amatrix_Matvec(ml_, LevelID_[0], ML_Operator_BlkMatMatvec);
      }
      else if (Avbr) {
        // check that the number of PDEs is constant
        if( Avbr->NumMyRows() % Avbr->NumMyBlockRows() != 0 ){
          std::cerr << "Error : NumPDEEqns does not seem to be constant" << std::endl;
          ML_CHK_ERR(-1);
        }
        ML_Set_Amatrix_Getrow(ml_, LevelID_[0], ML_Epetra_VbrMatrix_getrow,
                              ML_Epetra_VbrMatrix_comm_wrapper,
                              numMyRows+N_ghost);

        ML_Set_Amatrix_Matvec(ml_, LevelID_[0], ML_Epetra_VbrMatrix_matvec);
      }
      else if (Acrs) {
        ML_Set_Amatrix_Getrow(ml_, LevelID_[0], ML_Epetra_CrsMatrix_getrow,
                              ML_Epetra_CrsMatrix_comm_wrapper,
                              numMyRows+N_ghost);

        ML_Set_Amatrix_Matvec(ml_, LevelID_[0], ML_Epetra_CrsMatrix_matvec);
      }
      else {
        ML_Set_Amatrix_Getrow(ml_, LevelID_[0], ML_Epetra_RowMatrix_getrow,
                              ML_Epetra_comm_wrapper, numMyRows+N_ghost);

        ML_Set_Amatrix_Matvec(ml_, LevelID_[0], ML_Epetra_matvec);
      }
    }

  return 0;
}
int ML_Epetra::MultiLevelPreconditioner::
ConditionallyDestroyPreconditioner(const bool CheckPreconditioner)
{
  // ==============================================================  //
  // check whether the old filtering is still ok for the new matrix. //
  // This corresponds to ml_MultiLevelPreconditioner_Filtering.cpp   //
  // where there is code to determine if an old preconditioner is    //
  // still performing acceptably.
  // ==============================================================  //

  if (CheckPreconditioner == true && RateOfConvergence_ != -1.0) {

    // If the previous preconditioner was computed with option
    // "reuse: enable" == true, we know the rate of convergence
    // with the previous matrix (and this preconditioner. Now,
    // we recompute this ratio, and compare it with the previous one.
    // This requires an AztecOO object to be defined (ML must have been
    // configured with --enable-aztecoo)

    if (CheckPreconditionerKrylov() == false) {
      ML_CHK_ERR(DestroyPreconditioner());
    }
    else
      return 0;

  } else if (CheckPreconditioner == false &&
             IsComputePreconditionerOK_ == true) {

    // get rid of what done before
    ML_CHK_ERR(DestroyPreconditioner());

  } // nothing else is left
  else if ((verbose_) && (CheckPreconditioner == true) && (RateOfConvergence_ == -1.0) &&
           (Comm().MyPID() == 0)){
     printf("ML Warning: it appears as if 'reuse: enable' has not been set, but\n");
     printf("             ComputePreconditioner() has been invoked with a 'true'\n");
     printf("             argument indicating that the preconditioner is to be\n");
     printf("             checked. This might cause a memory leak!\n");
  }

  return(1);
}

// ================================================ ====== ==== ==== == =

int ML_Epetra::MultiLevelPreconditioner::
ComputePreconditioner(const bool CheckPreconditioner)
{

 try {

  if (ConditionallyDestroyPreconditioner(CheckPreconditioner) == 0) return 0;

  Epetra_Time Time(Comm());
  Epetra_Time InitialTime(Comm());
  Epetra_Time TotalTime(Comm());
  {
    int NumCompute = OutputList_.get("number of construction phases", 0);
    OutputList_.set("number of construction phases", ++NumCompute);
  }


  // =======================================================================//
  //              MPI Stuff                                                 //
  // =======================================================================//
#ifdef HAVE_ML_AZTECOO
#ifdef HAVE_MPI
  const Epetra_MpiComm & MpiComm =
    dynamic_cast<const Epetra_MpiComm&>(Comm());
  AZ_set_proc_config(ProcConfig_,MpiComm.Comm());
#else
  AZ_set_proc_config(ProcConfig_,AZ_NOT_MPI);
#endif
#endif
  // =======================================================================//
  //              End MPI Stuff                                             //
  // =======================================================================//

  if (List_.name() == "ANONYMOUS") List_.setName("ML preconditioner");

  bool init_seed = List_.get("initialize random seed",false);
  if (init_seed) ML_random_init();

  if (List_.get("read XML", false)) {
    std::string XMLFileName = List_.get("XML input file","ml_ParameterList.xml");
    ReadXML(XMLFileName,List_,Comm());
  }

  if (List_.get("ML debug mode", false)) ML_BreakForDebugger(*Comm_);
  int ProcID;
  if ((ProcID = List_.get("ML print initial list",-2)) > -2)
    if ((Comm().MyPID() == ProcID || ProcID == -1)) PrintList();

  // ========================================================================//
  //               Setting Basic Options                                     //
  // ========================================================================//
  int OutputLevel = List_.get("ML output", 0);                               //
  ML_Set_PrintLevel(OutputLevel);                                            //
                                                                             //
  verbose_             =((5 < ML_Get_PrintLevel()) && (Comm().MyPID() == 0));//
  FirstApplication_    =true;                                                //
  PrintMsg_            =List_.get("output prefix",PrintMsg_);                //
  NumLevels_           =List_.get("max levels",10);                          //
  MaxLevels_           =NumLevels_;                                          //
  AnalyzeMemory_       =List_.get("analyze memory", false);                  //
  CycleApplications_   =List_.get("cycle applications", 1);                  //
  ZeroStartingSolution_=List_.get("zero starting solution", true);           //
  mlpLabel_            = List_.get("ML label","not-set");                    //
  std::string IsIncreasing  =List_.get("increasing or decreasing", "increasing"); //


  if (AMGSolver_ == ML_MAXWELL) IsIncreasing = "decreasing";
          // JJH increasing not supported for Maxwell. Among other problems
          // Zoltan partitioning does not work with increasing+Maxwell
  int Direction;
  if (IsIncreasing == "increasing") Direction = ML_INCREASING;
  else                              Direction = ML_DECREASING;

  if (List_.get("low memory usage", false)) ML_Enable_LowMemory();
  else ML_Disable_LowMemory();

  if (Label_) delete[] Label_;
  Label_ = new char[80];

  const Epetra_VbrMatrix* VbrMatrix;
  VbrMatrix = dynamic_cast<const Epetra_VbrMatrix *>(RowMatrix_);
  if (VbrMatrix == 0) {
    NumPDEEqns_ = List_.get("PDE equations", 1);
  }
  else {
    int NumBlockRows = VbrMatrix->RowMap().NumGlobalElements();
    int NumRows = VbrMatrix->RowMap().NumGlobalPoints();
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        NumRows % NumBlockRows != 0,
        "# rows must be a multiple of # block rows ("
        << NumRows << "," << NumBlockRows << ")"
        );
    NumPDEEqns_ = NumRows/NumBlockRows;
  }
  // ========================================================================//
  //               End Setting Basic Options                                 //
  // ========================================================================//

  if( verbose_ ) ML_print_line("-",78);

  int call1 = 0, call2 = 0;
#ifdef ML_MALLOC
  int call1_malloc = 0, call2_malloc = 0;
#endif
  if( AnalyzeMemory_ ) {
    memory_[ML_MEM_INITIAL] = ML_MaxMemorySize();
    call1 = memory_[ML_MEM_INITIAL];
#ifdef ML_MALLOC
    call1_malloc = ML_MaxAllocatableSize();
    memory_[ML_MEM_INITIAL_MALLOC] = call1_malloc;
    if (verbose_) std::cout << "Memory : max allocatable block = " << call1_malloc << " Mbytes" << std::endl;
#endif
  }

  // Creates new list with level-specific smoother, level-specific
  // aggregation, and coarse options now in sublists.
  ParameterList newList;
  ML_CreateSublists(List_,newList);
  List_ = newList;
  // Validate Parameter List
  const int depth = List_.get("ML validate depth",5);
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      List_.get("ML validate parameter list",true)
      && !ValidateMLPParameters(List_,depth),
      "ERROR: ML's Teuchos::ParameterList contains incorrect parameter!\n\n"
      << "** IMPORTANT **\n\n"
      << "ML copies your parameter list and modifies copy:\n"
      << "   1) Level-specific smoother and aggregation options are placed in sublists\n"
      << "      called \"smoother: list (level XX)\" and \"aggregation: list (level XX)\",\n"
      << "      respectively.\n"
      << "   2) Coarse options are placed in a sublist called \"coarse: list\".\n"
      << "   3) In \"coarse: list\", any option that started with \"coarse:\" now starts\n"
      << "      with \"smoother:\".\n\n"
      );

  int MaxCreationLevels = NumLevels_;

  ML_Comm_Create(&ml_comm_);
#ifdef ML_MPI
  const Epetra_MpiComm *epcomm = dynamic_cast<const Epetra_MpiComm*>(&(RowMatrix_->Comm()));
  // Get the MPI communicator, as it may not be MPI_COMM_W0RLD, and update the ML comm object
  if (epcomm) ML_Comm_Set_UsrComm(ml_comm_,epcomm->Comm());
#endif
  ML_Create(&ml_,MaxCreationLevels);
  ml_->output_level = OutputLevel;
  ml_->comm->ML_nprocs = Comm().NumProc();
  ml_->comm->ML_mypid = Comm().MyPID();
#ifdef ML_MPI
  const Epetra_MpiComm* C2 = dynamic_cast<const Epetra_MpiComm*>(&Comm());
  ml_->comm->USR_comm = C2->Comm();
#endif
  profileIterations_ = List_.get("profile: operator iterations", 0);
  ML_Operator_Profile_SetIterations(profileIterations_);

  // CHECK IF x-coordinates WAS SET TO NULL AND REPARTITONTING IS USED. This could cause ML to hang
  // during repartitioning. The basic issue is that the coordinate pointers are checked to see if
  // they are null or not. If they are null, it is assumed that the user did not supply coordinates
  // and so the logic is a bit different. However, an empty processor might decided to set these
  // to null because it has no data. In this case, some proessors are viewed as supplying coordinates
  // while others are not. This messes up ML ... so we require that either all processors supply
  // non-null coordinate ptrs or none of them supply coordinates. For empty procs, one could simple
  // allocate length 1 vectors for the coordinates.

  int NumNullCoord = 0;
  if ( List_.isParameter("x-coordinates") ) {
    if ( (List_.get("x-coordinates",(double *) 0) == NULL) && (List_.get("repartition: enable",0)))
      NumNullCoord = 1;
  }
  NumNullCoord = ML_gsum_int(NumNullCoord, ml_comm_);
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
      (NumNullCoord != Comm().NumProc()) && (NumNullCoord != 0),
      "ERROR: ML's Teuchos::ParameterList should not have x-coordinates set to NULL even \n"
      <<"ERROR: if a processor has no matrix rows. Either a nonzero pointer must be given \n"
      <<"ERROR: (even for empty processors) or if the user does not wish to supply coordinates, \n"
      <<"ERROR: then ALL processors should simply not invoke list.set(\"x-coordinates\",...).\n"
      <<"ERROR: Otherwise ML is confused on some nodes as to whether or not a user supplied coordinates\n"
      );

  SetEigenScheme();

  SetLevelIds(Direction);

  SetFinestLevelMatrix();

  /* Set the RAP storage type */
  if(List_.get("use crs matrix storage",false)) ml_->RAP_storage_type=ML_CSR_MATRIX;
  else ml_->RAP_storage_type=ML_MSR_MATRIX;

  if (verbose_) {
    if( List_.isParameter("default values") ){
      std::cout << PrintMsg_ << "Default values for `" << List_.get("default values", "DD") << "'" << std::endl;
    }
    std::cout << PrintMsg_ << "Maximum number of levels = " << NumLevels_ << std::endl;
    if( IsIncreasing == "increasing" ) std::cout << PrintMsg_ << "Using increasing levels. ";
    else                               std::cout << PrintMsg_ << "Using decreasing levels. ";
    std::cout << "Finest level  = " << LevelID_[0];
    std::cout << ", coarsest level = " << LevelID_[NumLevels_-1] << std::endl;
    std::cout << "Number of applications of the ML cycle = " << CycleApplications_ << std::endl;
    std::cout << "Number of PDE equations = " << NumPDEEqns_ << std::endl;
  }

  if (AMGSolver_ == ML_MAXWELL) {

    // ====================================================================== //
    // create Maxwell nodal hierarchy.
    // ====================================================================== //

    ML_Create(&ml_nodes_,MaxCreationLevels);
    ml_nodes_->comm->ML_nprocs = Comm().NumProc();
    ml_nodes_->comm->ML_mypid  = Comm().MyPID();
    ml_nodes_->output_level = OutputLevel;
#ifdef ML_MPI
    ml_nodes_->comm->USR_comm = dynamic_cast<const Epetra_MpiComm*>(&Comm())->Comm();
#endif
    ML_Set_Label(ml_nodes_, const_cast<char*>("nodes"));
#ifdef HAVE_ML_EPETRAEXT
    NodeMatrix_ = ModifyEpetraMatrixColMap(*NodeMatrix_,
                                           NodeMatrixColMapTrans_,"Node");
#endif
    int numMyRows = NodeMatrix_->NumMyRows();
    int N_ghost   = NodeMatrix_->NumMyCols() - numMyRows;
    if (N_ghost < 0) N_ghost = 0;  // A->NumMyCols() = 0 for an empty matrix

    const Epetra_CrsMatrix *Acrs=dynamic_cast<const Epetra_CrsMatrix*>(NodeMatrix_);
    if (Acrs != 0) {
      ML_Init_Amatrix(ml_nodes_,LevelID_[0],numMyRows,numMyRows,
                      (void *) (const_cast<Epetra_CrsMatrix*>(Acrs)));
      ML_Set_Amatrix_Getrow(ml_nodes_, LevelID_[0], ML_Epetra_CrsMatrix_getrow,
                ML_Epetra_CrsMatrix_comm_wrapper, numMyRows+N_ghost);
      ML_Set_Amatrix_Matvec(ml_nodes_, LevelID_[0], ML_Epetra_CrsMatrix_matvec);
      ml_nodes_->Amat[LevelID_[0]].type = ML_TYPE_CRS_MATRIX;
    }
    else {
      ML_Init_Amatrix(ml_nodes_,LevelID_[0],numMyRows, numMyRows,
                      (void *) NodeMatrix_);
      ML_Set_Amatrix_Getrow(ml_nodes_, LevelID_[0], ML_Epetra_getrow,
                ML_Epetra_comm_wrapper, numMyRows+N_ghost);
      ML_Set_Amatrix_Matvec(ml_nodes_, LevelID_[0], ML_Epetra_matvec);
      ml_nodes_->Amat[LevelID_[0]].type = ML_TYPE_ROW_MATRIX;
    }
    ml_nodes_->Amat[LevelID_[0]].N_nonzeros = NodeMatrix_->NumMyNonzeros();

  } // if (AMGSolver_ == ML_MAXWELL)

  int ReqAggrePerProc = 128;

  if ((AMGSolver_ == ML_SA_FAMILY) || (AMGSolver_ == ML_MAXWELL)) {
     ML_Aggregate_Create(&agg_);
     ML_Aggregate_Set_MaxLevels(agg_, MaxCreationLevels);
     ML_Aggregate_VizAndStats_Setup(ml_);

     SetAggregation();
     double AggressiveFactor = .5;
     AggressiveFactor = List_.get("aggregation: aggressive", AggressiveFactor);
     ML_Aggregate_Set_Phase3AggregateCreationAggressiveness(agg_,AggressiveFactor);

     double Threshold = 0.0;
     Threshold = List_.get("aggregation: threshold", Threshold);
     ML_Aggregate_Set_Threshold(agg_,Threshold);

     double RowSumThreshold = -1.0;
     RowSumThreshold = List_.get("aggregation: rowsum threshold", RowSumThreshold);
     ML_Aggregate_Set_RowSum_Threshold(agg_,RowSumThreshold);

     int MaxCoarseSize = 128;
     if (List_.isSublist("coarse: list")) {
       ParameterList &coarseList = List_.sublist("coarse: list");
       MaxCoarseSize = coarseList.get("smoother: max size", MaxCoarseSize);
     }
     ML_Aggregate_Set_MaxCoarseSize(agg_, MaxCoarseSize );

     if( List_.isParameter("aggregation: req aggregates per process") )
       ReqAggrePerProc=List_.get("aggregation: req aggregates per process",
                                 ReqAggrePerProc);
     else {
       ReqAggrePerProc=List_.get("aggregation: next-level aggregates per process",
                                 ReqAggrePerProc);
     }

     ML_Aggregate_Set_ReqLocalCoarseSize(ml_->ML_num_levels,agg_,-1,
                                         ReqAggrePerProc);

     if( verbose_ ) {
        std::cout << PrintMsg_ << "Aggregation threshold = " << Threshold << std::endl;
        std::cout << PrintMsg_ << "Max coarse size = " << MaxCoarseSize << std::endl;

     }
  }
  if (AMGSolver_ == ML_SA_FAMILY) {
    ML_Aggregate_Set_Reuse(agg_);   // keep tentative prolongator in aux matrix
    agg_->keep_agg_information = 1; // or ReComputePreconditioner() case
    bool UseSymmetrize = List_.get("aggregation: symmetrize", false);
    if (UseSymmetrize == true) ML_Set_Symmetrize(ml_, ML_YES);
    else                       ML_Set_Symmetrize(ml_, ML_NO);

    if (!List_.get("energy minimization: enable", false))
      ML_CHK_ERR(SetSmoothingDamping());

    if (List_.get("aggregation: aux: enable", false)) {
      double Threshold   = List_.get("aggregation: aux: threshold", 0.0);
      int MaxAuxLevels   = List_.get("aggregation: aux: max levels", 10);
      ml_->Amat[LevelID_[0]].aux_data->threshold = Threshold;
      ml_->Amat[LevelID_[0]].aux_data->enable = 1;
      ml_->Amat[LevelID_[0]].aux_data->max_level = MaxAuxLevels;
      ml_->Amat[LevelID_[0]].num_PDEs = NumPDEEqns_;
    }

    if (List_.get("aggregation: block scaling", false) && NumPDEEqns_ != 1)
      agg_->block_scaled_SA = 1; // Unadvertised

    if (List_.get("aggregation: use tentative restriction", false))
      agg_->minimizing_energy = -1;

    // structured grid semicoarsening in z direction
    agg_->semicoarsen_levels = List_.get("semicoarsen: number of levels", -1);
    if (agg_->semicoarsen_levels != -1) {
      agg_->coarsen_rate      =  List_.get("semicoarsen: coarsen rate", -1);
      ml_->Amat[LevelID_[0]].num_PDEs = NumPDEEqns_;
      ml_->Amat[LevelID_[0]].NumZDir= List_.get("semicoarsen: line direction nodes", -1);
      std::string orient =List_.get("semicoarsen: line orientation", "use coordinates"); //
      ml_->Amat[LevelID_[0]].Zorientation= -1;
      if (orient == "vertical")    ml_->Amat[LevelID_[0]].Zorientation= 1;
      if (orient == "horizontal")  ml_->Amat[LevelID_[0]].Zorientation= 2;
            /* -1: not specified */ /*  1: vertical      */ /*  2: horizontal    */
      if (agg_->coarsen_rate == -1) {
         if (Comm().MyPID() == 0)
           std::cerr << ErrorMsg_ << "Must specify 'semicoarsen: coarsen rate' when using semicoarsening" << std::endl;
         ML_EXIT(-1);
      }
      if ((ml_->Amat[LevelID_[0]].NumZDir == -1) && (ml_->Amat[LevelID_[0]].Zorientation != -1)) {
         if (Comm().MyPID() == 0)
           std::cerr << ErrorMsg_ << "Must specify 'semicoarsen: line direction nodes' when using semicoarsening and orientation is given (i.e., not deduced from coordinates)" << std::endl;
         ML_EXIT(-1);
      }
      std::string coord_dir = List_.get("semicoarsen: coarsen coordinate","z");
      if (coord_dir == "x")       { agg_->semicoarsen_coordinate = 'x'; ml_->Amat[LevelID_[0]].coarsencoord= 'x'; }
      else if (coord_dir == "y")  { agg_->semicoarsen_coordinate = 'y'; ml_->Amat[LevelID_[0]].coarsencoord= 'y'; }
      else if (coord_dir == "z")  { agg_->semicoarsen_coordinate = 'z'; ml_->Amat[LevelID_[0]].coarsencoord= 'z'; }
      else {
           std::cerr << ErrorMsg_ << "Warning: Unrecognized coordinate direction for semicoarsening" << std::endl;
           std::cerr << ErrorMsg_ << coord_dir << std::endl;
      }
    }


    // energy minimization
    if (List_.get("energy minimization: enable", false)) {
      if (verbose_) {
        std::cout << std::endl;
        std::cout << "Warning: `energy minimization' safer with"  << std::endl;
        std::cout << "Warning: Uncoupled aggregation. Other "     << std::endl;
        std::cout << "Warning: schemes may crash the code."       << std::endl;
    //  std::cout << "Warning: Remember to use block scaling for" << std::endl;
    //  std::cout << "Warning: std::vector problem by setting"    << std::endl;
    //  std::cout << "Warning: `aggregation: block scaling'= true"<< std::endl;
        std::cout << std::endl;
      }
      agg_->minimizing_energy = List_.get("energy minimization: type", 2);
      agg_->minimizing_energy_droptol=List_.get("energy minimization: droptol",0.0);
      if (List_.get("energy minimization: cheap", false))
         agg_->cheap_minimizing_energy = 1;
    }
    ML_CHK_ERR(SetNullSpace());
  }

  if (AMGSolver_ == ML_MAXWELL) {
    ML_Aggregate_VizAndStats_Setup(ml_nodes_);
    ML_Aggregate_Set_ReqLocalCoarseSize(ml_nodes_->ML_num_levels,agg_,
                                       -1,ReqAggrePerProc);
    ML_Aggregate_Set_DampingFactor( agg_, 0.0);
  }

  ML_CHK_ERR(SetupCoordinates());

  if (List_.get("RAP: sort columns",0))                                     //
    ml_->sortColumnsAfterRAP = 1;
  // ========================================================================//
  //               Setting Repartitioning                                    //
  // ========================================================================//
                                                                             //
  if (List_.get("repartition: enable",0)) {                                  //
    ML_Repartition_Activate(ml_);                                            //
                                                                             //
    double minmax = List_.get("repartition: max min ratio", 1.3);            //
    ML_Repartition_Set_LargestMinMaxRatio(ml_,minmax);                       //
    int minperproc = List_.get("repartition: min per proc", 512);            //
    ML_Repartition_Set_MinPerProc(ml_,minperproc);                           //
    ML_Repartition_Set_PutOnSingleProc(ml_,                                  //
                    List_.get("repartition: put on single proc", 5000));     //
    int startLevel = List_.get("repartition: start level", 1);               //
    ML_Repartition_Set_StartLevel(ml_,startLevel);                           //
                                                                             //
    std::string Repartitioner=List_.get("repartition: partitioner","Zoltan");//
    if (Repartitioner == "Zoltan") {                                         //
       ML_Repartition_Set_Partitioner(ml_,ML_USEZOLTAN);                     //
       int NumDimensions = List_.get("repartition: Zoltan dimensions", 0);   //
       if (agg_ != NULL) ML_Aggregate_Set_Dimensions(agg_, NumDimensions);   //
       int zoltan_est_its=List_.get("repartition: estimated iterations",40); //
       bool zoltan_timers=List_.get("repartition: output timings",false);    //
       std::string zoltan_type=List_.get("repartition: Zoltan type","RCB");       //
       int smoother_steps = List_.get("smoother: sweeps", 2);                //
       if (List_.get("smoother: pre or post","post") == "pre or post")       //
         smoother_steps*=2;                                                  //
       int int_zoltan_type=ML_ZOLTAN_TYPE_RCB;                               //
       if(zoltan_type=="RCB")                  int_zoltan_type = ML_ZOLTAN_TYPE_RCB;
       else if(zoltan_type=="hypergraph")      int_zoltan_type = ML_ZOLTAN_TYPE_HYPERGRAPH;
       else if(zoltan_type=="fast hypergraph") int_zoltan_type = ML_ZOLTAN_TYPE_FAST_HYPERGRAPH;
       for(int i=0;i<ml_->ML_num_levels;i++){                                //
         ML_Aggregate_Viz_Stats * grid_info =(ML_Aggregate_Viz_Stats *)ml_->Grid[i].Grid;
         grid_info->zoltan_type          = int_zoltan_type;                  //
         grid_info->zoltan_estimated_its = zoltan_est_its;                   //
         grid_info->zoltan_timers        = (int)zoltan_timers;               //
         grid_info->smoothing_steps      = smoother_steps;                   //
       }                                                                     //
    }                                                                        //
    else if (Repartitioner == "ParMETIS")                                    //
      ML_Repartition_Set_Partitioner(ml_,ML_USEPARMETIS);                    //
    else {                                                                   //
      if (Comm().MyPID() == 0) {                                             //
        std::cerr << ErrorMsg_ << "Unrecognized partitioner `"
                  << Repartitioner << "'" << std::endl;
        std::cerr << ErrorMsg_ << "It should be: " << std::endl;
        std::cerr << ErrorMsg_ << "<Zoltan> / <ParMETIS>" << std::endl;
      }                                                                      //
      ML_EXIT(-1);                                                           //
    }                                                                        //
                                                                             //
    if (AMGSolver_ == ML_MAXWELL) {                                          //
      ML_Repartition_Activate(ml_nodes_);                                    //
      minmax = List_.get("repartition: node max min ratio", 1.3);            //
      ML_Repartition_Set_LargestMinMaxRatio(ml_nodes_,minmax);               //
      minperproc = List_.get("repartition: node min per proc", 170);         //
      ML_Repartition_Set_MinPerProc(ml_nodes_,minperproc);                   //
      if (Repartitioner == "Zoltan") {                                       //
        ML_Repartition_Set_Partitioner(ml_nodes_,ML_USEZOLTAN);              //
      }                                                                      //
    }                                                                        //
  }                                                                          //
  // ========================================================================//
  //               End of Repartition Settings                               //
  // ========================================================================//


  OutputList_.set("time: initial phase", InitialTime.ElapsedTime()
                  + OutputList_.get("time: initial phase", 0.0));
  InitialTime.ResetStartTime();
  Time.ResetStartTime();


  if (AMGSolver_ == ML_COMPOSITE) {

     // Create MG hierarchies for individual operators

    SubMatMLPrec_=(ML_Epetra::MultiLevelPreconditioner **) ML_allocate(NBlocks_*
                                 sizeof(ML_Epetra::MultiLevelPreconditioner *));

    for (int i = 0; i < NBlocks_ ; i++) {
      SubMatMLPrec_[i] = new ML_Epetra::MultiLevelPreconditioner(
                                          *DiagOperators_[i], DiagLists_[i]);
     }

     // Get the ml pointers, the finest level grids and compute the
     // maximum actual number of levels from the ml structures

     int ActualLevels = 0;
     void *tptr;
     ML *mlptr;
     int *CurrentLevel = NULL;
     CurrentLevel = (int *) ML_allocate(sizeof(int)*NBlocks_);
     for (int i = 0; i < NBlocks_; i++) {
        tptr = (void *) SubMatMLPrec_[i]->GetML();
        mlptr  = (ML *) tptr;
        CurrentLevel[i] =  mlptr->ML_finest_level;
        if ( mlptr->ML_num_actual_levels > ActualLevels)
          ActualLevels = mlptr->ML_num_actual_levels;
     }

     /********************************************************************/
     /* Build composite Rmat and Pmat operators by pulling out grid      */
     /* transfers from individual hierarchies and inserting them within  */
     /* the block matrices that define the composite AMG grid transfers  */
     /********************************************************************/

     ML_Operator *Rmat,*Pmat;
     ml_->ML_finest_level = LevelID_[0];
     for (int i = 0; i < ActualLevels-1; i++) {
         // NOTE: ML by default builds weird transpose operators in parallel
         //       for the restriction operator. ML's matmat mult should work
         //       we these ... but no one else will. My guess is that we
         //       should probably make an option for ML to compute standard
         //       matrices for restriction when we are doing this composite
         //       thing so that everyone can play with them.

         ML_Operator_BlkMatInit(&(ml_->Rmat[LevelID_[i]]), ml_comm_, NBlocks_,
                                     NBlocks_, ML_DESTROY_SHALLOW);
         ML_Operator_BlkMatInit(&(ml_->Pmat[LevelID_[i+1]]), ml_comm_, NBlocks_,
                                     NBlocks_, ML_DESTROY_SHALLOW);

         for (int j=0; j < NBlocks_; j++) {
            tptr = (void *) SubMatMLPrec_[j]->GetML();
            mlptr = (ML *) tptr;
            if (i < mlptr->ML_num_levels) {
               Rmat = &(mlptr->Rmat[CurrentLevel[j]]);
               if (Rmat->to != NULL) {
                  CurrentLevel[j] = Rmat->to->levelnum;
                  Pmat = &(mlptr->Pmat[CurrentLevel[j]]);
               }
               else Pmat = NULL;
            }
            else { Rmat = NULL; Pmat = NULL; }
            if (Rmat != NULL)
              ML_Operator_BlkMatInsert(&(ml_->Rmat[LevelID_[i]]),Rmat,j,j);
            if (Pmat != NULL)
              ML_Operator_BlkMatInsert(&(ml_->Pmat[LevelID_[i+1]]),Pmat,j,j);
         }
         ML_Operator_BlkMatFinalize(&(ml_->Rmat[LevelID_[i]]));
         ML_Operator_BlkMatFinalize(&(ml_->Pmat[LevelID_[i+1]]));
         ML_Operator_Set_1Levels(&(ml_->Pmat[LevelID_[i+1]]),
                        &(ml_->SingleLevel[LevelID_[i+1]]),
                        &(ml_->SingleLevel[LevelID_[i]]));
         ML_Operator_Set_1Levels(&(ml_->Rmat[LevelID_[i]]),
                        &(ml_->SingleLevel[LevelID_[i]]),
                        &(ml_->SingleLevel[LevelID_[i+1]]));
     }
     // Now build coarse discretization operators. Right now this is not
     // super efficient in that sometimes the diagonal blocks have already
     // been computed in the submatrix hiearchies ... and now we are going
     // to recompute  them.
     for (int i = 1; i < ActualLevels; i++) {
        ML_Blkrap(&(ml_->Rmat[LevelID_[i-1]]),&(ml_->Amat[LevelID_[i-1]]),
                  &(ml_->Pmat[LevelID_[i]]),&(ml_->Amat[LevelID_[i]]),
                  ML_CSR_MATRIX);
     }
     NumLevels_ = ActualLevels;
     if (CurrentLevel != NULL) ML_free(CurrentLevel);
  }
  if (AMGSolver_ == ML_SA_FAMILY) {
    // ==================================================================== //
    // Build smoothed aggregation hierarchy                                 //
    // ==================================================================== //
    NumLevels_ = ML_Gen_MultiLevelHierarchy_UsingAggregation(ml_,LevelID_[0],
                                                             Direction,agg_);

    if (verbose_)
      std::cout << PrintMsg_ << "Time to build the hierarchy = "
           << Time.ElapsedTime() << " (s)" << std::endl;

  }
  if (AMGSolver_ == ML_MAXWELL) {
    // ==================================================================== //
    // Build Maxwell hierarchy                                              //
    // ==================================================================== //

#ifdef HAVE_ML_EPETRAEXT
    TMatrix_ = ModifyEpetraMatrixColMap(*TMatrix_,
                                        TMatrixColMapTrans_,"Gradient");
#endif
    Apply_BCsToGradient( *EdgeMatrix_, *TMatrix_);

    if (CurlCurlMatrix_ != 0) CheckNullSpace();

    if( TMatrixML_ == 0 ) {
      // convert TMatrix to ML_Operator
      TMatrixML_ = ML_Operator_Create(ml_->comm);
      ML_Operator_WrapEpetraMatrix(const_cast<Epetra_RowMatrix*>(TMatrix_),
                                   TMatrixML_);
    }

    TMatrixTransposeML_ = ML_Operator_Create(ml_->comm);
    ML_Operator_Transpose_byrow(TMatrixML_,TMatrixTransposeML_);

    double EdgeDampingFactor =
              List_.get("aggregation: damping factor",ML_DDEFAULT);
    double enrichBeta = List_.get("aggregation: enrich beta",0.0);

    double PeDropThreshold = List_.get("aggregation: edge prolongator drop threshold",ML_DDEFAULT);

    // Set up arrays for holding the curl-curl and mass matrices separately.
    typedef ML_Operator* pML_Operator;  // to avoid compiler warning
    CurlCurlMatrix_array = new pML_Operator[NumLevels_];
    for (int i=0; i< NumLevels_; i++) CurlCurlMatrix_array[i] = NULL;
    if (CurlCurlMatrix_) {
      CurlCurlMatrix_array[LevelID_[0]] = ML_Operator_Create(ml_->comm);
      ML_Operator_WrapEpetraMatrix(
                            const_cast<Epetra_RowMatrix*>(CurlCurlMatrix_),
                            CurlCurlMatrix_array[LevelID_[0]] );
    }
    MassMatrix_array = new pML_Operator[NumLevels_];
    for (int i=0; i< NumLevels_; i++) MassMatrix_array[i] = NULL;
    if (MassMatrix_) {
      MassMatrix_array[LevelID_[0]] = ML_Operator_Create(ml_->comm);
      ML_Operator_WrapEpetraMatrix( const_cast<Epetra_RowMatrix*>(MassMatrix_),
                                  MassMatrix_array[LevelID_[0]] );
    }

    NumLevels_ = ML_Gen_MGHierarchy_UsingReitzinger(ml_,&ml_nodes_,
                            LevelID_[0], Direction,agg_,
                            CurlCurlMatrix_array,
                            MassMatrix_array,
                            TMatrixML_,TMatrixTransposeML_,
                            &Tmat_array,&Tmat_trans_array,
                            EdgeDampingFactor,
                            enrichBeta, PeDropThreshold);

  }
  {
    int NL2 = OutputList_.get("max number of levels", 0); //what the hell???
    OutputList_.set("max number of levels", NL2+NumLevels_);
  }

  if( AnalyzeMemory_ ) {
    call2 = ML_MaxMemorySize();
    memory_[ML_MEM_HIERARCHY] = call2-call1;
    call1 = call2;
#ifdef ML_MALLOC
    call2_malloc = ML_MaxAllocatableSize();
    if (verbose_) std::cout << "Memory : max allocatable block = " << call2_malloc << " Mbytes" << std::endl;
    memory_[ML_MEM_HIERARCHY_MALLOC] = call1_malloc-call2_malloc;
    call1_malloc = call2_malloc;
#endif
  }

  if( verbose_ ) std::cout << PrintMsg_ << "Number of actual levels : " << NumLevels_ << std::endl;

  OutputList_.set("time: hierarchy", InitialTime.ElapsedTime()
                  + OutputList_.get("time: hierarchy", 0.0));
  InitialTime.ResetStartTime();

  // ====================================================================== //
  // Generate all smoothers and coarse grid solver.                         //
  // ====================================================================== //

#ifdef NewStuff
  if (!DontSetSmoothers_)
#endif
     ML_CHK_ERR(SetSmoothers());

  InitialTime.ResetStartTime();

  if (AnalyzeMemory_) {
    call2 = ML_MaxMemorySize();
    memory_[ML_MEM_SMOOTHER] = call2 - call1;
    call1 = call2;
#ifdef ML_MALLOC
    call2_malloc = ML_MaxAllocatableSize();
    if (verbose_) std::cout << "Memory : max allocatable block = " << call2_malloc << " Mbytes" << std::endl;
    memory_[ML_MEM_SMOOTHER_MALLOC] = call1_malloc - call2_malloc;
    call1_malloc = call2_malloc;
#endif
  }

  ownership_ = false;

  if (AnalyzeMemory_) {
    call2 = ML_MaxMemorySize();
    memory_[ML_MEM_COARSE] = call2 - call1;
    call1 = call2;
#ifdef ML_MALLOC
    call2_malloc = ML_MaxAllocatableSize();
    if (verbose_) std::cout << "Memory : max allocatable block = " << call2_malloc << " Mbytes" << std::endl;
    memory_[ML_MEM_COARSE_MALLOC] = call1_malloc - call2_malloc;
    call1_malloc = call2_malloc;
#endif
  }

  ML_Gen_Solver(ml_, ML_MGV, LevelID_[0], LevelID_[NumLevels_-1]);

  // ====================================================================== //
  // Use of filtering functions, here called `filtering: enable'            //
  // If this option is true, the code detects the non-converging modes of   //
  // I - ML^{-1}A (where ML is the preconditioner we have just built), and  //
  // creates a new V cycle to be added to the preconditioner. This part is  //
  // equivalent to the GGB files of Haim Waisman (files ml_struct.c and     //
  // ml_ggb.c, in the Main subdirectory).                                   //
  // ====================================================================== //

  ML_CHK_ERR(SetFiltering());

  MatrixDumper();

  CreateLabel();

  SetPreconditioner();

  IsComputePreconditionerOK_ = true;

  // ====================================================================== //
  // Compute the rate of convergence (for reuse preconditioners)            //
  // ====================================================================== //

  if( List_.get("reuse: enable", false) == true )
    CheckPreconditionerKrylov();

  if (AnalyzeMemory_) {
    memory_[ML_MEM_FINAL] = ML_MaxMemorySize();
    memory_[ML_MEM_TOT1] = memory_[ML_MEM_FINAL] - memory_[ML_MEM_INITIAL];
#ifdef ML_MALLOC
    memory_[ML_MEM_FINAL_MALLOC] = ML_MaxAllocatableSize();
    memory_[ML_MEM_TOT1_MALLOC] = memory_[ML_MEM_INITIAL_MALLOC] - memory_[ML_MEM_FINAL_MALLOC];
#endif
  }

  // print unused parameters
  if (List_.isParameter("print unused")) {
    ProcID = List_.get("print unused",-2);
    if ((Comm().MyPID() == ProcID || ProcID == -1) && verbose_) PrintUnused();
  }
  if ((ProcID = List_.get("ML print final list",-2)) > -2) {
    if ((Comm().MyPID() == ProcID || ProcID == -1)) PrintList();
  }

  // ===================================================================== //
  // compute the coordinates for each level (that is, the center of        //
  // gravity, as no mesh is really available for coarser levels)           //
  // ===================================================================== //

  OutputList_.set("time: final setup", InitialTime.ElapsedTime()
                  + OutputList_.get("time: final setup", 0.0));
  InitialTime.ResetStartTime();

  if (ML_Get_PrintLevel() >= 10 && Comm().MyPID() == 0) {
    std::cout << std::endl;
    std::cout << "Cumulative timing for construction so far: " << std::endl;
    std::cout << PrintMsg_ << "- for initial setup   = "
         << OutputList_.get("time: initial phase", 0.0) << " (s)" << std::endl;
    std::cout << PrintMsg_ << "- for hierarchy setup = "
         << OutputList_.get("time: hierarchy", 0.0) << " (s)" << std::endl;
    std::cout << PrintMsg_ << "- for smoothers setup = "
         << OutputList_.get("time: smoothers setup", 0.0) << " (s)" << std::endl;
    std::cout << PrintMsg_ << "- for coarse setup    = "
         << OutputList_.get("time: coarse solver setup", 0.0) << " (s)" << std::endl;
    std::cout << PrintMsg_ << "- for final setup     = "
         << OutputList_.get("time: final setup", 0.0) << " (s)" << std::endl;
    std::cout << PrintMsg_ << "Total for this setup  = "
         << TotalTime.ElapsedTime() << " (s)" << std::endl;
  }


  if( verbose_ ) ML_print_line("-",78);

  ConstructionTime_ += Time.ElapsedTime();
  OutputList_.set("time: construction", ConstructionTime_);
  ++NumConstructions_;

  VisualizeAggregates();
 }
 catch(const std::exception& e)
 {
   if (Comm().MyPID() == 0) {
     const char * what = e.what();

     fprintf(stderr,"%s","\n*********************************************************\n");
     fprintf(stderr,"%s","ML failed to compute the multigrid preconditioner. The\n");
     fprintf(stderr,"%s","most common problem is an incorrect  data type in ML's\n");
     fprintf(stderr,"%s","parameter list (e.g. 'int' instead of 'bool').\n\n");
     fprintf(stderr,"%s","Note: List.set(\"ML print initial list\",X) might help\nfigure out the bad one on pid X.\n");
     fprintf(stderr,"%s","\nThe following exception was thrown:\n");
     fprintf(stderr,"%s",what);
     fprintf(stderr,"%s","*********************************************************\n\n");
   }
   ML_CHK_ERR(-1);
 }
 catch(...)
 {
   if (Comm().MyPID() == 0) {
     fprintf(stderr,"%s","\n*********************************************************\n");
     fprintf(stderr,"%s","ML failed to compute the multigrid preconditioner. The\n");
     fprintf(stderr,"%s","most common problem is an incorrect  data type in ML's\n");
     fprintf(stderr,"%s","parameter list (e.g. 'int' instead of 'bool').\n\n");
     fprintf(stderr,"%s","Note: List.set(\"ML print initial list\",X) might help\nfigure out the bad one on pid X.\n");
     fprintf(stderr,"%s","*********************************************************\n\n");
   }

   ML_CHK_ERR(-1);
 }

 return(0);
}

// ================================================ ====== ==== ==== == =

int ML_Epetra::MultiLevelPreconditioner::
ReComputePreconditioner(bool keepFineLevelSmoother)
{

 try{

  if (AMGSolver_ == ML_MAXWELL)
    ML_CHK_ERR(-1);

  if (IsPreconditionerComputed() == false)
    ML_CHK_ERR(-2);

  IsComputePreconditionerOK_ = false;

  // =========================== //
  // re-build the preconditioner //
  // =========================== //

  int OutputLevel = List_.get("ML output", -47);
  if (OutputLevel == -47) OutputLevel = List_.get("output", 0);
  ML_Set_PrintLevel(OutputLevel);

  Epetra_Time Time(Comm());
  Epetra_Time InitialTime(Comm());
  InitialTime.ResetStartTime();
  Time.ResetStartTime();

  {
    int NumCompute = OutputList_.get("number of construction phases", 0);
    OutputList_.set("number of construction phases", ++NumCompute);
  }

  for (int i = 0; i < ml_->ML_num_levels; i++)
  {
    if (i!=0 || keepFineLevelSmoother==false) {
      ML_Smoother_Clean(&(ml_->pre_smoother[i]));
      ML_Smoother_Init(&(ml_->pre_smoother[i]), &(ml_->SingleLevel[i]));
      ML_Smoother_Clean(&(ml_->post_smoother[i]));
      ML_Smoother_Init(&(ml_->post_smoother[i]), &(ml_->SingleLevel[i]));
      ML_CSolve_Clean(&(ml_->csolve[i]));
      ML_CSolve_Init(&(ml_->csolve[i]));
    }
  }

  if( verbose_ )
    ML_print_line("-",78);

  FirstApplication_ = true;

  if (verbose_) {
    std::cout << PrintMsg_ << "Re-computing the preconditioner..." << std::endl;
  }

  profileIterations_ = List_.get("profile: operator iterations", 0);
  ML_Operator_Profile_SetIterations(profileIterations_);
  ML_CHK_ERR(SetNullSpace());
  NumPDEEqns_ = List_.get("PDE equations", 1);
  ml_->Amat[LevelID_[0]].num_PDEs = NumPDEEqns_;
  ML_Gen_MultiLevelHierarchy_UsingSmoothedAggr_ReuseExistingAgg(ml_, agg_);

  if (verbose_)
    std::cout << PrintMsg_ << "Time to re-build the hierarchy = "
         << Time.ElapsedTime() << " (s)" << std::endl;
  Time.ResetStartTime();

  if (verbose_)
    std::cout << PrintMsg_ << "Number of actual levels : " << NumLevels_ << std::endl;

  OutputList_.set("time: hierarchy", Time.ElapsedTime()
                  + OutputList_.get("time: hierarchy", 0.0));
  Time.ResetStartTime();

  // ====================================================================== //
  // Generate all smoothers and coarse grid solver.                         //
  // ====================================================================== //

  try{SetSmoothers(keepFineLevelSmoother);}
  //try{SetSmoothers();}
  catch(...) {
    if (Comm().MyPID() == 0) {
      fprintf(stderr,"%s","\n**************************\n");
      fprintf(stderr,"%s","Problem setting smoothers.\n");
      fprintf(stderr,"%s","****************************\n\n");
      ML_CHK_ERR(-1);
    }
  }
  Time.ResetStartTime();

  ML_Gen_Solver(ml_, ML_MGV, LevelID_[0], LevelID_[NumLevels_-1]);

  IsComputePreconditionerOK_ = true;

  OutputList_.set("time: final setup", Time.ElapsedTime()
                  + OutputList_.get("time: final setup", 0.0));
  Time.ResetStartTime();

  if (ML_Get_PrintLevel() == 10 && Comm().MyPID() == 0) {
    std::cout << std::endl;
    std::cout << "Cumulative timing for construction so far: " << std::endl;
    std::cout << PrintMsg_ << "- for initial setup   = "
      << OutputList_.get("time: initial phase", 0.0) << " (s)" << std::endl;
    std::cout << PrintMsg_ << "- for hierarchy setup = "
      << OutputList_.get("time: hierarchy", 0.0) << " (s)" << std::endl;
    std::cout << PrintMsg_ << "- for smoothers setup = "
      << OutputList_.get("time: smoothers setup", 0.0) << " (s)" << std::endl;
    std::cout << PrintMsg_ << "- for coarse setup    = "
      << OutputList_.get("time: coarse solver setup", 0.0) << " (s)" << std::endl;
    std::cout << PrintMsg_ << "- for final setup     = "
      << OutputList_.get("time: final setup", 0.0) << " (s)" << std::endl;
  }

  /* ------------------- that's all folks --------------------------------- */

  if (verbose_)
    ML_print_line("-",78);

  ConstructionTime_ += InitialTime.ElapsedTime();
  OutputList_.set("time: construction", ConstructionTime_);
  ++NumConstructions_;
  OutputList_.set("number of constructions", NumConstructions_);

 }
 catch(...)
 {
   ML_CHK_ERR(-1);
 }

 return(0);

}

// ================================================ ====== ==== ==== == =

/*! Print the individual operators in the multigrid hierarchy.

 \param level (In) By default, this method prints the entire
                   multigrid hierarchy. If you desire only the operators
                   associated with a particular level, pass in the level number.
                   The fine level is 0, coarser levels are positive integers.
*/

void ML_Epetra::MultiLevelPreconditioner::
Print(int level)
{
  // ====================================================================== //
  // One may decide to print out the entire hierarchy (to be analyzed in    //
  // MATLAB, for instance).                                                 //
  // ====================================================================== //

  if ( (level < -1) || (level > NumLevels_-1) )
    return;

  if(Comm().MyPID()==0) {
    if (level == -1) {
      std::cout << std::endl;
      std::cout << PrintMsg_ << "You are printing the entire multigrid hierarchy," << std::endl
                << PrintMsg_ << "from finest level (" << LevelID_[0]
                << ") to coarsest (" << LevelID_[NumLevels_-1] << ")." << std::endl
                << PrintMsg_ << "MATLAB can be used to load the matrices, using spconvert()" << std::endl;
      std::cout << std::endl;
    } else {
      std::cout << PrintMsg_ << "You are printing the multigrid operators on level " << LevelID_[level]
                             << "." << std::endl
                << PrintMsg_ << "(The finest level is " << LevelID_[0] << ".)" << std::endl
                << PrintMsg_ << "MATLAB can be used to load the matrices, using spconvert()" << std::endl;
      std::cout << std::endl;
    }
  }

  ML* mlptr = ml_;
  char name[80];

  if (level > -1) {

    sprintf(name,"Amat_%d", LevelID_[level]);
    ML_Operator_Print_UsingGlobalOrdering(mlptr->Amat+LevelID_[level], name, NULL,NULL);
    if (level > 0) {
      sprintf(name,"Pmat_%d", LevelID_[level]);
      ML_Operator_Print_UsingGlobalOrdering(mlptr->Pmat+LevelID_[level+1], name, NULL,NULL);
    }
    if (level < NumLevels_-1) {
      sprintf(name,"Rmat_%d", LevelID_[level]);
      ML_Operator_Print_UsingGlobalOrdering(mlptr->Rmat+LevelID_[level], name, NULL,NULL);
    }

  } else {

    // Amat (one for each level)
    for( int i=0 ; i<NumLevels_ ; ++i ) {
      sprintf(name,"Amat_%d", LevelID_[i]);
      ML_Operator_Print_UsingGlobalOrdering(mlptr->Amat+LevelID_[i], name, NULL,NULL);
    }

    // Pmat (one for each level, except last)
    for( int i=1 ; i<NumLevels_ ; ++i ) {
      sprintf(name,"Pmat_%d", LevelID_[i]);
      ML_Operator_Print_UsingGlobalOrdering(mlptr->Pmat+LevelID_[i], name, NULL,NULL);
    }

    // Rmat (one for each level, except first)
    for( int i=0 ; i<NumLevels_-1 ; ++i ) {
      sprintf(name,"Rmat_%d", LevelID_[i]);
      ML_Operator_Print_UsingGlobalOrdering(mlptr->Rmat+LevelID_[i], name, NULL,NULL);
    }

  } //if-then-else

} //Print() method

// ============================================================================

void ML_Epetra::MultiLevelPreconditioner::PrintUnused(const int MyPID) const
{
  if( Comm().MyPID() == MyPID ) {
    ML_print_line("-",78);
    std::cout << PrintMsg_ << "Unused parameters:" << std::endl;
    PrintUnused();
    ML_print_line("-",78);
  }
}

// ============================================================================

void ML_Epetra::MultiLevelPreconditioner::PrintList()
{
  ML_print_line("+",78);
  std::cout << "+++ Printing ML parameter list \"" << List_.name()
       << "\" on pid " << Comm().MyPID() << std::endl;
  ML_print_line("+",78);
  List_.print(std::cout);
  ML_print_line("-",49);
  std::cout << "----------- end of ML parameter list ------------" << std::endl;
  ML_print_line("-",49);
  std::cout << std::endl;
}

// ============================================================================

int ML_Epetra::MultiLevelPreconditioner::
SetParameterList(const ParameterList & List)
{
  if( IsComputePreconditionerOK_ == true ) DestroyPreconditioner();
  List_ = List;
  return 0;

}

// ============================================================================

int ML_Epetra::MultiLevelPreconditioner::CreateLabel()
{

  char finest[80];
  char coarsest[80];
  finest[0] = '\0';
  coarsest[0] = '\0';
  char * label;

  int i = ml_->ML_finest_level;

  if (ml_->pre_smoother[i].smoother->func_ptr != NULL) {
    label = ml_->pre_smoother[i].label;
    if( strncmp(label,"PreS_",4) == 0 ) sprintf(finest, "%s", "~");
    else                                sprintf(finest, "%s", label);
  } else                                sprintf(finest, "%s", "~");

  if (ml_->post_smoother[i].smoother->func_ptr != NULL) {
    label = ml_->post_smoother[i].label;
    if( strncmp(label,"PostS_", 5) == 0 ) sprintf(finest, "%s/~", finest);
    else                                  sprintf(finest, "%s/%s", finest, label);
  } else                                  sprintf(finest, "%s/~", finest);

  if (i != ml_->ML_coarsest_level) {
    i = ml_->ML_coarsest_level;
    if ( ML_CSolve_Check( &(ml_->csolve[i]) ) == 1 ) {
    sprintf(coarsest, "%s", ml_->csolve[i].label);
    }

    else {
      if (ml_->pre_smoother[i].smoother->func_ptr != NULL) {
    label = ml_->pre_smoother[i].label;
    if( strncmp(label,"PreS_",4) == 0 ) sprintf(coarsest, "%s", "~");
    else                                sprintf(coarsest, "%s", label);
      } else                                sprintf(coarsest, "%s", "~");
      if (ml_->post_smoother[i].smoother->func_ptr != NULL) {
    label = ml_->post_smoother[i].label;
    if( strncmp(label,"PostS_", 5) == 0 ) sprintf(coarsest, "%s/~", coarsest);
    else                                  sprintf(coarsest, "%s/%s",coarsest, label);
      } else                                  sprintf(coarsest, "%s/~", coarsest);
    }
  }

  if (AMGSolver_ != ML_MAXWELL)
    sprintf(Label_, "ML (L=%d, %s, %s)",
     ml_->ML_num_actual_levels, finest, coarsest);
  else
    sprintf(Label_, "ML (Maxwell, L=%d, %s, %s)",
     ml_->ML_num_actual_levels, finest, coarsest);

  return 0;

} //MLP::CreateLabel()

// ============================================================================

int ML_Epetra::MultiLevelPreconditioner::
ApplyInverse(const Epetra_MultiVector& X,
         Epetra_MultiVector& Y) const
{

  int before = 0, after = 0;
#ifdef ML_MALLOC
  int before_malloc = 0, after_malloc = 0;
#endif

  if (AnalyzeMemory_) {
#ifdef ML_MALLOC
    before_malloc = ML_MaxAllocatableSize();
#endif
    before = ML_MaxMemorySize();
  }

  Epetra_Time Time(Comm());

  // ========================================================= //
  // I do not perfom checks on maps, like the following,       //
  // if (!X.Map().SameAs(OperatorDomainMap())) ML_CHK_ERR(-1)  //
  // if (!Y.Map().SameAs(OperatorRangeMap())) ML_CHK_ERR(-1)   //
  // because they can be expensive                             //
  //                                                           //
  // We allocate xtmp, needed in case X is scaled in solver or //
  // X = Y. We cannot pre-allocate xtmp, since the user might  //
  // give in input different number of vectors.                //
  //                                                           //
  // Finally, since does not handle multivectors, we need to   //
  // extract and iterate one at a time.                        //
  //                                                           //
  // FIXME: the Y.PutScalar(0.0) can probably be skipped.      //
  // ========================================================= //

  if (Y.NumVectors() != X.NumVectors())
    ML_CHK_ERR(-3);
  if( !IsPreconditionerComputed() )
    ML_CHK_ERR(-10);

  Epetra_MultiVector xtmp(X);

  if (ZeroStartingSolution_) Y.PutScalar(0.0);

  double** xvectors;
  double** yvectors;
  ML_CHK_ERR(xtmp.ExtractView(&xvectors));
  ML_CHK_ERR(Y.ExtractView(&yvectors));

  for (int i = 0; i < X.NumVectors(); ++i) {

    // ================================================== //
    // Added support for multiple cycles on 08-Mar-05     //
    // Tested for ML_Cycle_MG only                        //
    //                                                    //
    // note that all the DD stuff is "delicate"           //
    // in the sense that the coarse solver must be Amesos //
    //                                                    //
    // I am not sure if ML_NONZERO applies to something   //
    // else than ML_Cycle_MG().                           //
    //                                                    //
    // The flt_ml_ is the filtering, which requires a     //
    // suitable setup phase. Note that in the current     //
    // implementation the                                 //
    // resulting preconditioner is always non-symmetric   //
    // ================================================== //

    for (int ia = 0 ; ia < CycleApplications_ ; ++ia) {

      int StartingSolution;
      if (ia || !ZeroStartingSolution_)
        StartingSolution = ML_NONZERO;
      else
        StartingSolution = ML_ZERO;

      switch(ml_->ML_scheme) {
      case(ML_MGFULLV):
        ML_Solve_MGFull(ml_, xvectors[i], yvectors[i]);
        break;
      case(ML_SAAMG): //Marian Brezina's solver
        ML_Solve_AMGV(ml_, xvectors[i], yvectors[i]);
        break;
      case(ML_PAMGV): //V-cycle, where modes are projected out before & after
        ML_Solve_ProjectedAMGV(ml_, xvectors[i], yvectors[i]);
        break;
      case(ML_ONE_LEVEL_DD):
        ML_DD_OneLevel(&(ml_->SingleLevel[ml_->ML_finest_level]),
                       yvectors[i], xvectors[i],
                       ML_ZERO, ml_->comm, ML_NO_RES_NORM, ml_);
        break;
      case(ML_TWO_LEVEL_DD_ADD):
        ML_DD_Additive(&(ml_->SingleLevel[ml_->ML_finest_level]),
                       yvectors[i], xvectors[i],
                       ML_ZERO, ml_->comm, ML_NO_RES_NORM, ml_);
        break;
      case(ML_TWO_LEVEL_DD_HYBRID):
        ML_DD_Hybrid(&(ml_->SingleLevel[ml_->ML_finest_level]),
                     yvectors[i], xvectors[i],
                     ML_ZERO, ml_->comm, ML_NO_RES_NORM, ml_);
        break;
    case(ML_TWO_LEVEL_DD_HYBRID_2):
      ML_DD_Hybrid_2(&(ml_->SingleLevel[ml_->ML_finest_level]),
             yvectors[i], xvectors[i],
             ML_ZERO, ml_->comm, ML_NO_RES_NORM, ml_);
      break;
    default:
#ifdef ReverseOrder
       ML *ml_ggb = (ML *) ml_->void_options;
       throw "not done\n";
       /*some code would need to go here to change the order of the ML Cycles */
       /*there would probably also need to be some additional changes elsewhere*/
#else
      ML_Cycle_MG(&(ml_->SingleLevel[ml_->ML_finest_level]),
                  yvectors[i], xvectors[i], StartingSolution,
                  ml_->comm, ML_NO_RES_NORM, ml_);
#endif
      }
    }

    if (flt_ml_) {
      ML_Cycle_MG(&(flt_ml_->SingleLevel[flt_ml_->ML_finest_level]),
          yvectors[i], xvectors[i],
          ML_NONZERO, flt_ml_->comm, ML_NO_RES_NORM, flt_ml_);
    }
  }
  // ====== //
  // timing //
  // ====== //

  MultiLevelPreconditioner * This = const_cast<MultiLevelPreconditioner *>(this);

  if (AnalyzeMemory_) {
#ifdef ML_MALLOC
    after_malloc = ML_MaxAllocatableSize();
#endif
    after = ML_MaxMemorySize();
  }

  double t = Time.ElapsedTime();
  if (FirstApplication_) {
    This->FirstApplication_ = false;
    This->FirstApplicationTime_ += t;
    This->OutputList_.set("time: first application", This->FirstApplicationTime_);
    This->memory_[ML_MEM_PREC_FIRST] = after - before;
#ifdef ML_MALLOC
    This->memory_[ML_MEM_PREC_FIRST_MALLOC] = before_malloc - after_malloc;
#endif
  }
  else {
    This->memory_[ML_MEM_PREC_OTHER] = after - before;
#ifdef ML_MALLOC
    This->memory_[ML_MEM_PREC_OTHER_MALLOC] = before_malloc - after_malloc;
#endif
    This->ApplicationTime_ += t;
    This->OutputList_.set("time: total apply", This->FirstApplicationTime_+This->ApplicationTime_);
  }

  ++(This->NumApplications_);
  This->OutputList_.set("number of applications", This->NumApplications_);

  return 0;
} //MLP::ApplyInverse

// ============================================================================
/*! Values for \c "coarse: type"
 * - \c Jacobi
 * - \c Gauss-Seidel
 * - \c symmetric Gauss-Seidel
 * - \c MLS
 * - \c Hiptmair (Maxwell only)
 * - \c SuperLU (deprecated)
 * - \c Amesos-KLU
 * - \c Amesos-UMFPACK
 * - \c Amesos-Superludist
 * - \c Amesos-Superlu
 * - \c Amesos-MUMPS
 * - \c Amesos-ScALAPACK (under development in Amesos)
 * - \c do-nothing
 */
int ML_Epetra::MultiLevelPreconditioner::SetCoarse()
{

  std::string CoarseSolution = List_.get("coarse: type", "Amesos-KLU");
  int NumSmootherSteps = List_.get("coarse: sweeps", 2);
  double Omega = List_.get("coarse: damping factor", 1.0);
  double AddToDiag = List_.get("coarse: add to diag", 1e-12);
  std::string PreOrPostSmoother = List_.get("coarse: pre or post","post");
  int pre_or_post;

  if( PreOrPostSmoother == "pre" ) pre_or_post = ML_PRESMOOTHER;
  else if( PreOrPostSmoother == "both" ) pre_or_post = ML_BOTH;
  else pre_or_post = ML_POSTSMOOTHER;

  // rst: Changing polynomial interface:
  //    1) polynomial degree is set from "coarse: sweeps"
  //    2) "coarse: Chebyshev" also calls ML's Chebyshev
  //    3) "coarse: Chebshev alpha" also sets alpha
  //    4) "coarse: node sweeps" and "coarse: edge sweeps" now set polynomial
  //                degree within Hiptmair
  //
  // For backward compatiblity, we still take the degree from
  // "coarse: MLS polynomial order" if set, still recognize MLS

  double ChebyshevAlpha = List_.get("coarse: MLS alpha",-2.0);
  if ( ChebyshevAlpha == -2.) ChebyshevAlpha = List_.get("coarse: Chebyshev alpha", 30.);

  int ChebyshevPolyOrder = List_.get("coarse: MLS polynomial order",-7);
  if (ChebyshevPolyOrder == -7) ChebyshevPolyOrder = NumSmootherSteps;

  char msg[80];
  sprintf(msg, "Coarse solve (level %d) : ", LevelID_[NumLevels_-1]);

  int MaxProcs = List_.get("coarse: max processes", -1);

  if( CoarseSolution == "Jacobi" ) {
    if( verbose_ ) std::cout << msg << "Jacobi (sweeps="
                        << NumSmootherSteps << ",omega=" << Omega
                        << "," << PreOrPostSmoother << ")" << std::endl;
    ML_Gen_Smoother_Jacobi(ml_, LevelID_[NumLevels_-1], pre_or_post,
               NumSmootherSteps, Omega);

  } else if( CoarseSolution == "Gauss-Seidel" ) {
    if( verbose_ ) std::cout << msg << "Gauss-Seidel (sweeps="
                        << NumSmootherSteps << ",omega=" << Omega
                        << "," << PreOrPostSmoother << ")" << std::endl;
#ifdef HAVE_ML_IFPACK
    if (ml_->Amat[LevelID_[NumLevels_-1]].type == ML_TYPE_CRS_MATRIX) {
      if (verbose_)
        std::cout << msg << "Epetra_CrsMatrix detected, using "
             << "Ifpack implementation" << std::endl;
      std::string IfpackType = "point relaxation stand-alone";
      ParameterList& IfpackList = List_.sublist("smoother: ifpack list");
      IfpackList.set("relaxation: type", "Gauss-Seidel");
      IfpackList.set("relaxation: sweeps", NumSmootherSteps);
      IfpackList.set("relaxation: damping factor", Omega);
      ML_Gen_Smoother_Ifpack(ml_, IfpackType.c_str(),
                             0, LevelID_[NumLevels_-1], pre_or_post,
                             (void*)(&IfpackList),(void *) Comm_);
    }
    else
#endif
      ML_Gen_Smoother_GaussSeidel(ml_, LevelID_[NumLevels_-1], pre_or_post,
               NumSmootherSteps, Omega);

  } else if( CoarseSolution == "ML Gauss-Seidel" ) {
    if( verbose_ ) std::cout << msg << "Gauss-Seidel (sweeps="
                        << NumSmootherSteps << ",omega=" << Omega
                        << "," << PreOrPostSmoother << ")" << std::endl;
    ML_Gen_Smoother_GaussSeidel(ml_, LevelID_[NumLevels_-1], pre_or_post,
                                NumSmootherSteps, Omega);

  } else if( CoarseSolution == "symmetric Gauss-Seidel" ) {
    if( verbose_ ) std::cout << msg << "symmetric Gauss-Seidel (sweeps="
                        << NumSmootherSteps << ",omega=" << Omega
                        << "," << PreOrPostSmoother << ")" << std::endl;
#ifdef HAVE_ML_IFPACK
    if (ml_->Amat[LevelID_[NumLevels_-1]].type == ML_TYPE_CRS_MATRIX) {
      if (verbose_)
        std::cout << msg << "Epetra_CrsMatrix detected, using "
             << "Ifpack implementation" << std::endl;
      std::string IfpackType = "point relaxation stand-alone";
      ParameterList& IfpackList = List_.sublist("smoother: ifpack list");
      IfpackList.set("relaxation: type", "symmetric Gauss-Seidel");
      IfpackList.set("relaxation: sweeps", NumSmootherSteps);
      IfpackList.set("relaxation: damping factor", Omega);
      ML_Gen_Smoother_Ifpack(ml_, IfpackType.c_str(),
                             0, LevelID_[NumLevels_-1], pre_or_post,
                             //IfpackList,*Comm_);
                             (void*)&IfpackList,(void*)Comm_);
    }
    else
#endif
      ML_Gen_Smoother_SymGaussSeidel(ml_, LevelID_[NumLevels_-1],
                     pre_or_post, NumSmootherSteps, Omega);

  } else if( CoarseSolution == "ML symmetric Gauss-Seidel" ) {
    if( verbose_ ) std::cout << msg << "symmetric Gauss-Seidel (sweeps="
                        << NumSmootherSteps << ",omega=" << Omega
                        << "," << PreOrPostSmoother << ")" << std::endl;
    ML_Gen_Smoother_SymGaussSeidel(ml_, LevelID_[NumLevels_-1],
                                   pre_or_post, NumSmootherSteps, Omega);

  } else if(( CoarseSolution == "MLS" ) || (CoarseSolution == "Chebyshev")) {
  if( verbose_ ) std::cout << msg << "Chebyshev (degree="
                        << ChebyshevPolyOrder << ",alpha=" << ChebyshevAlpha
                        << "," << PreOrPostSmoother << ")" << std::endl;
     ML_Gen_Smoother_Cheby(ml_, LevelID_[NumLevels_-1], pre_or_post,
                         ChebyshevAlpha, ChebyshevPolyOrder);

  } else if( CoarseSolution == "Hiptmair" ) {

      if (AMGSolver_ != ML_MAXWELL) {
        if (Comm().MyPID() == 0) {
          std::cerr << ErrorMsg_ << "Hiptmair smoothing is only supported" << std::endl;
          std::cerr << ErrorMsg_ << "for solving eddy current equations." << std::endl;
          std::cerr << ErrorMsg_ << "Choose another smoother." << std::endl;
        }
        ML_EXIT(EXIT_FAILURE);
      }

      std::string SubSmootherType = List_.get("coarse: subsmoother type","MLS");
      int nodal_its = List_.get("coarse: node sweeps", 2);
      int edge_its = List_.get("coarse: edge sweeps", 2);

      int logical_level = LevelID_[NumLevels_-1];
      void *edge_smoother = 0, *nodal_smoother = 0;

#if defined(__GNUC__) && defined(__GNUC_MINOR__) && defined(__GNUC_PATCHLEVEL__)
#define GCC_VERSION __GNUC__*100+__GNUC_MINOR__*10+__GNUC_PATCHLEVEL__
#endif

#if defined(GCC_VERSION) && GCC_VERSION >= 460
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-pedantic"
#endif
      double omega;
      char subsmDetails[80];
      // only MLS and SGS are currently supported
      if ( (SubSmootherType == "MLS") || (SubSmootherType == "Chebyshev"))
      {
        nodal_smoother=(void *) ML_Gen_Smoother_Cheby;
        nodal_args_ = ML_Smoother_Arglist_Create(2);
        edge_args_ = ML_Smoother_Arglist_Create(2);
        // set polynomial degree
        ML_Smoother_Arglist_Set(nodal_args_, 0, &nodal_its);
        edge_smoother=(void *) ML_Gen_Smoother_Cheby;
        ML_Smoother_Arglist_Set(edge_args_, 0, &edge_its);

        // set alpha (lower bound of smoothing interval (alpha,beta) )

        ML_Smoother_Arglist_Set(edge_args_, 1, &ChebyshevAlpha);
        ML_Smoother_Arglist_Set(nodal_args_, 1, &ChebyshevAlpha);
        sprintf(subsmDetails,"edge sweeps=%d, node sweeps=%d,alpha=%3.1e",
                edge_its,nodal_its,ChebyshevAlpha);
      }
      else if (SubSmootherType == "symmetric Gauss-Seidel") {
        omega = List_.get("coarse: damping factor",1.0);
        nodal_smoother=(void *) ML_Gen_Smoother_SymGaussSeidel;
        nodal_args_ = ML_Smoother_Arglist_Create(2);
        edge_args_ = ML_Smoother_Arglist_Create(2);
        ML_Smoother_Arglist_Set(nodal_args_, 0, &nodal_its);
        ML_Smoother_Arglist_Set(nodal_args_, 1, &omega);
        edge_smoother=(void *) ML_Gen_Smoother_SymGaussSeidel;
        ML_Smoother_Arglist_Set(edge_args_, 0, &edge_its);
        ML_Smoother_Arglist_Set(edge_args_, 1, &omega);
        sprintf(subsmDetails,"edge sweeps=%d,node sweeps=%d,omega=%3.1e",
                edge_its,nodal_its,omega);
      }
      else if (Comm().MyPID() == 0)
        std::cerr << ErrorMsg_ << "Only Chebyshev(or MLS) and SGS are supported as "
             << "Hiptmair subsmoothers ... not " << SubSmootherType << std::endl;

      int hiptmair_type = (int)
             List_.get("smoother: Hiptmair efficient symmetric", true);

     if( verbose_ ) std::cout << msg << "Hiptmair " << "(outer sweeps="
                          << NumSmootherSteps << ","
                          << PreOrPostSmoother << "," << SubSmootherType
                          << "," << subsmDetails << ")" << std::endl;
     ML_Gen_Smoother_Hiptmair2(ml_, logical_level,pre_or_post,
                                NumSmootherSteps, Tmat_array, Tmat_trans_array, NULL,
                                MassMatrix_array,TtATMatrixML_,
                                edge_smoother, edge_args_, nodal_smoother, nodal_args_,
                                hiptmair_type);

      ML_Smoother_Arglist_Delete(&nodal_args_);
      ML_Smoother_Arglist_Delete(&edge_args_);

#if defined(GCC_VERSION) && GCC_VERSION >= 460
#pragma GCC diagnostic pop
#endif

  } else if( CoarseSolution == "SuperLU" )
    ML_Gen_CoarseSolverSuperLU( ml_, LevelID_[NumLevels_-1]);
  else if( CoarseSolution == "Amesos-LAPACK" ) {
    ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1],
                           ML_AMESOS_LAPACK, MaxProcs, AddToDiag);
  }
  else if( CoarseSolution == "Amesos-KLU" ) {
    ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1],
                           ML_AMESOS_KLU, MaxProcs, AddToDiag);
  } else if( CoarseSolution == "Amesos-UMFPACK" )
    ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1],
                           ML_AMESOS_UMFPACK, MaxProcs, AddToDiag);
  else if(  CoarseSolution == "Amesos-Superludist" )
    ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1],
                           ML_AMESOS_SUPERLUDIST, MaxProcs, AddToDiag);
  else if(  CoarseSolution == "Amesos-Superlu" )
    ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1],
                           ML_AMESOS_SUPERLU, MaxProcs, AddToDiag);
  else if( CoarseSolution == "Amesos-MUMPS" )
    ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1],
                           ML_AMESOS_MUMPS, MaxProcs, AddToDiag);
  else if( CoarseSolution == "Amesos-ScaLAPACK" )
    ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1],
                           ML_AMESOS_SCALAPACK, MaxProcs, AddToDiag);

 else if( CoarseSolution == "block Gauss-Seidel" ) {
    if( verbose_ ) std::cout << msg << "block Gauss-Seidel (sweeps="
                        << NumSmootherSteps << ",omega=" << Omega
                        << "," << PreOrPostSmoother << "numPDEs="
                        << NumPDEEqns_ << ")" << std::endl;
    ML_Gen_Smoother_BlockGaussSeidel(ml_, LevelID_[NumLevels_-1], pre_or_post,
                                        NumSmootherSteps, Omega, NumPDEEqns_);

  } else if( CoarseSolution == "symmetric block Gauss-Seidel" ) {
    if( verbose_ ) std::cout << msg << "symmetric block Gauss-Seidel (sweeps="
                        << NumSmootherSteps << ",omega=" << Omega
                        << "," << PreOrPostSmoother << "numPDEs="
                        << NumPDEEqns_ << ")" << std::endl;
    ML_Gen_Smoother_SymBlockGaussSeidel(ml_, LevelID_[NumLevels_-1],pre_or_post,
                                        NumSmootherSteps, Omega, NumPDEEqns_);

  } else if(CoarseSolution == "user-defined"||CoarseSolution == "user defined"){
    // ============ //
    // user-defined //
    // ============ //

    int (*userSmootherPtr)(ML_Smoother *, int, double *, int, double *);
    userSmootherPtr = NULL;
    userSmootherPtr = List_.get("coarse: user-defined function",
                                userSmootherPtr);
    std::string userSmootherName;
    userSmootherName = List_.get("coarse: user-defined name", "User-defined");

    if( verbose_ ) std::cout << msg << userSmootherName << " (sweeps="
			 << NumSmootherSteps << ", " << PreOrPostSmoother << ")" << std::endl;

    if (userSmootherPtr == NULL) {
      if (Comm().MyPID() == 0)
        std::cerr << ErrorMsg_
             << "No pointer to user-defined smoother function found." << std::endl;
      ML_EXIT(EXIT_FAILURE);
    }
    ML_Operator *data;
    ML_Get_Amatrix(ml_, LevelID_[NumLevels_-1], &data);
    ML_Set_Smoother(ml_, LevelID_[NumLevels_-1], pre_or_post, data,
                    userSmootherPtr,
                    const_cast<char *>(userSmootherName.c_str()));

  } else if( CoarseSolution == "do-nothing" ) {
    if( verbose_ ) std::cout << msg << "No Coarse Solve" << std::endl;

    // do nothing, ML will not use any coarse solver
  } else {
    if( Comm().MyPID() == 0 ) {
      std::cout << ErrorMsg_ << "specified options for coarse solver ("
           << CoarseSolution << ") not valid. Should be:" << std::endl;
      std::cout << ErrorMsg_ << "<Jacobi> / <Gauss-Seidel> / <symmetric Gauss-Seidel /" << std::endl;
      std::cout << ErrorMsg_ << "<MLS> / <Hiptmair> / <Amesos-LAPACK> / <Amesos-KLU> /" << std::endl;
      std::cout << ErrorMsg_ << "<Amesos-UMFPACK> / <Amesos-Superludist> / <Amesos-Superlu> /" << std::endl;
      std::cout << ErrorMsg_ << "<Amesos-MUMPS> / <block Gauss-Seidel> /" << std::endl;
      std::cout << ErrorMsg_ << "<symmetric block Gauss-Seidel> / <user-defined>"
           << std::endl;
    }

    ML_EXIT(-1);
  }

  return 0;
}

// ============================================================================
/*! Values for \c "aggregation: type"
 * - \c Uncoupled-MIS
 * - \c METIS
 * - \c ParMETIS
 * - \c Uncoupled
 * - \c Coupled (deprecated)
 * - \c MIS
 * - \c user
 * - \c greedy
 */
#include "ml_agg_user.h"
int ML_Epetra::MultiLevelPreconditioner::SetAggregation()
{
  char aggListName[80];

  int value = -777; /* pagina 777 di televideo */
  std::string CoarsenScheme = List_.get("aggregation: type","Uncoupled");

  if (CoarsenScheme == "VBMETIS")   // Not advertised in manual
  {
     ML_Aggregate_Set_CoarsenScheme_VBMETIS(agg_);
     int  nblocks   = List_.get("aggregation: nblocks",0);
     int  blockdim  = List_.get("aggregation: length blocks",0);
     int* blocks    = List_.get("aggregation: blocks", (int*)NULL);
     int* block_pde = List_.get("aggregation: block_pde", (int*)NULL);
     ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS(agg_,0,NumLevels_,nblocks,
                                                    blocks,block_pde,blockdim);
     int nodesperagg = List_.get("aggregation: nodes per aggregate",0);
     if (!nodesperagg) {
       std::cout << "*ML*WRN* aggregation: nodes per aggregate no set, using 9\n";
       nodesperagg = 9;
     }
     for (int i=0; i<NumLevels_; ++i)
       ML_Aggregate_Set_NodesPerAggr(ml_,agg_,i,nodesperagg);
  }
  else {
     for( int level=0 ; level<NumLevels_-1 ; ++level ) {
       sprintf(aggListName,"aggregation: list (level %d)",LevelID_[level]);
       ParameterList &aggList = List_.sublist(aggListName);

       std::string MyCoarsenScheme = aggList.get("aggregation: type",CoarsenScheme);

       if (MyCoarsenScheme == "METIS")
         ML_Aggregate_Set_CoarsenSchemeLevel_METIS(level,NumLevels_,agg_);
       else if (MyCoarsenScheme == "ParMETIS")
         ML_Aggregate_Set_CoarsenSchemeLevel_ParMETIS(level,NumLevels_,agg_);
       else if (MyCoarsenScheme == "Zoltan")
         ML_Aggregate_Set_CoarsenSchemeLevel_Zoltan(level,NumLevels_,agg_);
       else if (MyCoarsenScheme == "MIS")
         ML_Aggregate_Set_CoarsenSchemeLevel_MIS(level,NumLevels_,agg_);
       else if (MyCoarsenScheme == "Uncoupled")
         ML_Aggregate_Set_CoarsenSchemeLevel_Uncoupled(level,NumLevels_,agg_);
       else if (MyCoarsenScheme == "Uncoupled-MIS")
         ML_Aggregate_Set_CoarsenSchemeLevel_UncoupledMIS(level,NumLevels_,agg_);
       else if (MyCoarsenScheme == "Coupled")
         ML_Aggregate_Set_CoarsenSchemeLevel_Coupled(level,NumLevels_,agg_);
       else if (MyCoarsenScheme == "user")
         ML_Aggregate_Set_CoarsenSchemeLevel_User(level,NumLevels_,agg_);
#ifdef MARZIO
       else if (MyCoarsenScheme == "greedy") {
         // MS // this is delicate as it burns the "user" aggregation,
         // MS // should we fix it?
         // MS // 05-Dec-04
         ML_SetUserLabel(ML_Aggregate_GreedyLabel);
         ML_SetUserPartitions(ML_Aggregate_Greedy);
         ML_Aggregate_Set_CoarsenSchemeLevel_User(level,NumLevels_,agg_);
       }
#endif
       else {
         TEUCHOS_TEST_FOR_EXCEPT_MSG(
             true,
             ErrorMsg_ << "specified options ("
             << MyCoarsenScheme << ") not valid. Should be:" << std::endl
             << ErrorMsg_ << "<METIS> / <ParMETIS> / <Zoltan> /<MIS> / <Uncoupled> / <Coupled> / <user>" << std::endl
             );
       }

       /* FIXME then DELETEME: check that user still works,
        * then delete this part
       if (CoarsenScheme == "Zoltan" || CoarsenScheme == "user") {
         // This copies the coordinates if the aggregation scheme
         // of at least one level is Zoltan. Coordinates will be
         // projected for ALL levels independently of the
         // aggregation scheme.

         double * coord = List_.get("aggregation: coordinates", (double *)0);
         int NumDimensions = List_.get("aggregation: dimensions", 0);

         ML_Aggregate_Set_NodalCoordinates(ml_, agg_, coord);
         ML_Aggregate_Set_Dimensions(agg_, NumDimensions);

       }
       */

       if( MyCoarsenScheme == "METIS" || MyCoarsenScheme == "ParMETIS" ||
           MyCoarsenScheme == "Zoltan" ) {

         bool isSet = false;

         // first look for parameters without any level specification

         if( List_.isParameter("aggregation: global aggregates") ){
           value = -777; // simply means not set
           value = List_.get("aggregation: global aggregates",value);
           if( value != -777 ) {
             ML_Aggregate_Set_GlobalNumber(ml_,agg_,LevelID_[level],value );
             isSet = true;
           }
         }

         if( List_.isParameter("aggregation: local aggregates") ){
           value = -777;
           value = List_.get("aggregation: local aggregates",value);
           if( value != -777 ) {
             ML_Aggregate_Set_LocalNumber(ml_,agg_,LevelID_[level],value );
             isSet = true;
           }
             }

         if( List_.isParameter("aggregation: nodes per aggregate") ){
           value = -777;
           value = List_.get("aggregation: nodes per aggregate",value);
           if( value != -777 ) {
             ML_Aggregate_Set_NodesPerAggr(ml_,agg_,LevelID_[level],value );
             isSet = true;
           }
         }

         // now for level-specific data

         if( aggList.isParameter("aggregation: global aggregates") ){
           value = -777; // simply means not set
           value = aggList.get("aggregation: global aggregates",value);
           if( value != -777 ) {
             ML_Aggregate_Set_GlobalNumber(ml_,agg_,LevelID_[level],value );
             isSet = true;
           }
         }

         if( aggList.isParameter("aggregation: local aggregates") ){
           value = -777;
           value = aggList.get("aggregation: local aggregates",value);
           if( value != -777 ) {
             ML_Aggregate_Set_LocalNumber(ml_,agg_,LevelID_[level],value );
             isSet = true;
           }
         }

         if( aggList.isParameter("aggregation: nodes per aggregate") ){
           value = -777;
           value = aggList.get("aggregation: nodes per aggregate",value);
           if( value != -777 ) {
             ML_Aggregate_Set_NodesPerAggr(ml_,agg_,LevelID_[level],value );
             isSet = true;
           }
         }

         if( isSet == false ) {
           // put default values
           value = aggList.get("aggregation: local aggregates",1);
           ML_Aggregate_Set_LocalNumber(ml_,agg_,LevelID_[level],value);
         }

       } // if( CoarsenScheme == "METIS" || CoarsenScheme == "ParMETIS" )

     } /* for */
  } /* else */

  return 0;
}

// ================================================ ====== ==== ==== == =

int ML_Epetra::MultiLevelPreconditioner::SetPreconditioner()
{
  std::string str = List_.get("prec type","MGV");

  if( str == "one-level-postsmoothing" ) {

    sprintf(Label_,  "%s","1-level postsmoothing only");
    ml_->ML_scheme = ML_ONE_LEVEL_DD;

  } else if( str == "two-level-additive" ) {

    sprintf(Label_,  "%s","two-level additive DD");
    ml_->ML_scheme = ML_TWO_LEVEL_DD_ADD;
    if( NumLevels_ != 2 ) {
      if( Comm().MyPID() == 0 ) {
    std::cerr << ErrorMsg_ << "You asked for `two-level additive DD' but you don't have" << std::endl
         << ErrorMsg_ << "exacty two levels. Now continue, but check you input..." << std::endl;
      }
    }

  } else if( str == "two-level-hybrid") {

    sprintf(Label_,  "%s","two-level hybrid DD");
    ml_->ML_scheme = ML_TWO_LEVEL_DD_HYBRID;
    if( NumLevels_ != 2 ) {
      if( Comm().MyPID() == 0 ) {
    std::cerr << ErrorMsg_ << "You asked for `two-level hybrid DD' but you don't have" << std::endl
         << ErrorMsg_ << "exacty two levels. Now continue, but check you input..." << std::endl;
      }
    }

  } else if( str == "two-level-hybrid2") {

    sprintf(Label_,  "%s","two-level hybrid DD (2)");
    ml_->ML_scheme = ML_TWO_LEVEL_DD_HYBRID_2;
    if( NumLevels_ != 2 ) {
      if( Comm().MyPID() == 0 ) {
    std::cerr << ErrorMsg_ << "You asked for `two-level hybrid DD (2)' but you don't have" << std::endl
         << ErrorMsg_ << "exacty two levels. Now continue, but check you input..." << std::endl;
      }
    }

  } else if( str == "full-MGV" ) {
    ml_->ML_scheme = ML_MGFULLV;

  } else if( str == "MGW" ) {
    ml_->ML_scheme = ML_MGW;

  } else if( str == "projected MGV" ) {
    sprintf(Label_,  "%s","projected MGV");
    ml_->ML_scheme = ML_PAMGV;

    int numModes = List_.get("number of projected modes",0);
    double **periodicModes = List_.get("projected modes", (double **) 0);

    // Check that the number of modes is 1-3, and that the modes are there.
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        numModes < 1 || numModes > 3,
        ErrorMsg_ <<
        "You have chosen `projected MGV', but `number of projected modes'\n"
        << ErrorMsg_ <<
        " has an incorrect value.  It should be 1, 2 or 3."
        );
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        periodicModes == 0,
        ErrorMsg_ << "You have chosen `projected MGV', but `projected modes' is NULL."
        );
    for (int i=0; i<numModes; i++) {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
          periodicModes[i] == 0,
          ErrorMsg_ << "You have chosen `projected MGV', but mode " << i+1 << " is NULL.\n"
          );
    }

    //JJH 7-22-05  I think NumMyRows is the right value...
    ML_Operator_SetSubspace(ml_, periodicModes, numModes,
                            RowMatrix_->NumMyRows());

  } else if( str == "MGV" ) {
    // it is the default
    ml_->ML_scheme = ML_MGV;
  } else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        true,
        ErrorMsg_ << "`prec type' has an incorrect value of '"
        << str << "'. It should be\n"
        << ErrorMsg_ << "<one-level-postsmoothing> / <two-level-additive>\n"
        << ErrorMsg_ << "<two-level-hybrid> / <two-level-hybrid2>"
        );
  }

  return 0;
}

// ================================================ ====== ==== ==== == =

int ML_Epetra::MultiLevelPreconditioner::SetSmoothingDamping()
{

  Epetra_Time Time(Comm());

  // Check for additional smoothing of prolongator.  This should be used in
  // conjunction with more aggressive coarsening.
  int PSmSweeps = List_.get("aggregation: smoothing sweeps", 1);
  char aggListName[80];
  for (int i=0; i<MaxLevels_; i++) {
    sprintf(aggListName,"aggregation: list (level %d)",LevelID_[i]);
    ParameterList &aggList = List_.sublist(aggListName);
    int MyPSmSweeps = aggList.get("aggregation: smoothing sweeps",PSmSweeps);
    ML_Aggregate_Set_DampingSweeps(agg_,MyPSmSweeps,LevelID_[i]);
  }
  //ML_Aggregate_Set_BlockDiagScaling(agg_);

  /* ********************************************************************** */
  /* Strategies to determine the field-of-values.                           */
  /* almost everything here is experimental ;)                              */
  /* ********************************************************************** */

  std::string RandPSmoothing = List_.get("R and P smoothing: type", "classic");

  /* start looping over different options */

  if( RandPSmoothing == "classic" ) {

    /* ********************************************************************** */
    /* For "classical" approach to determine lambda_max only.                 */
    /* ********************************************************************** */

    ML_CHK_ERR(SetSmoothingDampingClassic());

  } else if( RandPSmoothing == "advanced" ) {

    /* ********************************************************************** */
    /* This is the new way, based on Anasazi to compute eigen-widgets         */
    /* ********************************************************************** */

    //    if( verbose_ )
    //      std::cout << PrintMsg_ << "Use A to smooth restriction operator" << std::endl;
    agg_->Restriction_smoothagg_transpose = ML_TRUE;

    // fix default values in List_
    List_.set("eigen-analysis: use symmetric algorithm",false);
    List_.set("eigen-analysis: tolerance", 1e-2);
    List_.set("eigen-analysis: use diagonal scaling", true);
    List_.set("eigen-analysis: restart", 100);
    List_.set("eigen-analysis: length", 20);
    List_.set("field-of-values: tolerance", 1e-2);
    List_.set("field-of-values: use diagonal scaling", true);
    List_.set("field-of-values: restart", 100);
    List_.set("field-of-values: ", 20);
    List_.set("field-of-values: print current status", false);
    List_.set("ML output",0);

    struct ML_Field_Of_Values* field_of_values;

    // stick default values (undefined)
    field_of_values = (struct ML_Field_Of_Values *) ML_allocate( sizeof(struct ML_Field_Of_Values) );
    field_of_values->eta     = 0.0;
    field_of_values->real_max= -1.0;
    field_of_values->imag_max= -1.0;
    field_of_values->poly_order = 0;

    if( List_.get("aggregation: compute field of values",true) )
      field_of_values->compute_field_of_values = ML_YES;
    else
      field_of_values->compute_field_of_values = ML_NO;

    if( List_.get("aggreation: compute field of values for non-scaled",false) )
      field_of_values->compute_field_of_values_non_scaled = ML_YES;
    else
      field_of_values->compute_field_of_values_non_scaled = ML_NO;

    // following values for choice:
    // -1 : undefined
    //  0 : do nothing, put eta = 0 (used in non-smoothed aggregation)
    //  1 : compute the box of the field of values (eta = imag_max/real_max)
    //  2 : compute ||lambda_max|| (eta = std::sqrt(imag_max^2 + real_max^2))
    field_of_values->choice     = -1;
    // and this is a pointer for the object's interal ParameterList
    field_of_values->EigenList = (void *) &List_;

    // still to set up polynomial coeffiecients
    std::string DampingType =  List_.get("R and P smoothing: damping", "default");

    if( DampingType == "non-smoothed" ) {

      if( verbose_ )
    std::cout << PrintMsg_ << "R and P smoothing : non-smoothed aggregation" << std::endl;

      field_of_values->choice     =  0;
      field_of_values->poly_order =  0;
      // I don't really need them, smoothing will be set to zero
      field_of_values->R_coeff[0] =  0.0;
      field_of_values->R_coeff[1] =  0.0;
      field_of_values->R_coeff[2] =  0.0;

      field_of_values->P_coeff[0] =  0.0;
      field_of_values->P_coeff[1] =  0.0;
      field_of_values->P_coeff[2] =  0.0;

      ML_Aggregate_Set_DampingFactor(agg_,0.0);

      agg_->Restriction_smoothagg_transpose = ML_FALSE;

    } else if( DampingType == "almost-non-smoothed" ) {

      if( verbose_ )
    std::cout << PrintMsg_ << "R and P smoothing : almost non-smoothed aggregation" << std::endl;

      field_of_values->choice     =  1;
      field_of_values->poly_order =  0;
      // I don't really need them, smoothing will be set to zero
      field_of_values->R_coeff[0] =  0.000000001;
      field_of_values->R_coeff[1] =  0.0;
      field_of_values->R_coeff[2] =  0.0;

      field_of_values->P_coeff[0] =  0.000000001;
      field_of_values->P_coeff[1] =  0.0;
      field_of_values->P_coeff[2] =  0.0;

      ML_Aggregate_Set_DampingFactor(agg_,0.0);

    } else if( DampingType == "fov-1" ) {

      // those are the coefficients proposed by Ray
      if( verbose_ )
    std::cout << PrintMsg_ << "R and P smoothing : Using `fov-1' values" << std::endl;

      field_of_values->choice     =  1;
      field_of_values->poly_order =  2;

      field_of_values->R_coeff[0] =  1.107;
      field_of_values->R_coeff[1] =  0.285;
      field_of_values->R_coeff[2] =  0.718;

      field_of_values->P_coeff[0] =  1.878;
      field_of_values->P_coeff[1] = -2.515;
      field_of_values->P_coeff[2] =  0.942;

      ML_Aggregate_Set_DampingFactor(agg_,0.0);

    } else if( DampingType == "fov-2" ) {

      if( verbose_ )
    std::cout << PrintMsg_ << "R and P smoothing : Using `fov-2' values" << std::endl;

      field_of_values->choice     =  1;
      field_of_values->poly_order =  2;

      field_of_values->R_coeff[0] =  1.138;
      field_of_values->R_coeff[1] =  1.162;
      field_of_values->R_coeff[2] = -2.384;

      field_of_values->P_coeff[0] =  2.143;
      field_of_values->P_coeff[1] = -2.179;
      field_of_values->P_coeff[2] =  0.101;

      ML_Aggregate_Set_DampingFactor(agg_,0.0);

    } else if( DampingType == "fov-5" ) {

      if( verbose_ )
    std::cout << PrintMsg_ << "R and P smoothing : Using `fov-5' values" << std::endl;

      field_of_values->choice     =  1;
      field_of_values->poly_order =  2;

      field_of_values->R_coeff[0] = 1.631;
      field_of_values->R_coeff[1] = -3.9015;
      field_of_values->R_coeff[2] = 2.5957;

      field_of_values->P_coeff[0] = 1.00145;
      field_of_values->P_coeff[1] = -1.4252;
      field_of_values->P_coeff[2] = 0.6627;

      ML_Aggregate_Set_DampingFactor(agg_,0.0);

    } else if( DampingType == "fov-10" ) {

      // those are Marzio's best values ;^)
      if( verbose_ )
    std::cout << PrintMsg_ << "R and P smoothing : Using `fov-10' values" << std::endl;

      field_of_values->choice     =  1;
      field_of_values->poly_order =  2;


      field_of_values->R_coeff[0] = 1.768909e+00;
      field_of_values->R_coeff[1] = -4.132227e+00;
      field_of_values->R_coeff[2] = 2.669318e+00;
      field_of_values->P_coeff[0] = 1.619455e+00;
      field_of_values->P_coeff[1] = -2.347773e+00;
      field_of_values->P_coeff[2] = 8.652273e-01;

      ML_Aggregate_Set_DampingFactor(agg_,0.0);

    } else if( DampingType == "random" ) {

      // only to play with
      if( verbose_ )
        std::cout << PrintMsg_ << "R and P smoothing : Using `random' values" <<std::endl;

      field_of_values->choice     =  1;
      field_of_values->poly_order =  2;

      // initialize seed
      unsigned int s = (int)(Time.ElapsedTime()*10000);
      srand(s);

      // put random values
      field_of_values->R_coeff[0] =  (double)rand()/RAND_MAX;
      field_of_values->R_coeff[1] =  (double)rand()/RAND_MAX;
      field_of_values->R_coeff[2] =  (double)rand()/RAND_MAX;

      field_of_values->P_coeff[0] =  (double)rand()/RAND_MAX;
      field_of_values->P_coeff[1] =  (double)rand()/RAND_MAX;
      field_of_values->P_coeff[2] =  (double)rand()/RAND_MAX;

      if( verbose_ ) {
        std::cout << PrintMsg_ << "Random coefficients for R and P:" << std::endl
             << PrintMsg_ << field_of_values->R_coeff[0] << "   "
             << field_of_values->R_coeff[1] << "   "
             << field_of_values->R_coeff[2] << std::endl
             << PrintMsg_ << field_of_values->P_coeff[0] << "   "
             << field_of_values->P_coeff[1] << "   "
             << field_of_values->P_coeff[2] << "   (seed = "
             << s << ")" << std::endl;
      }

      ML_Aggregate_Set_DampingFactor(agg_,0.0);

    } else if( DampingType == "classic" ) {

      // This must be as with classic ML approach.
      // It can be used to estimate the field of value, tough.

      if( verbose_ )
        std::cout << PrintMsg_ << "R and P smoothing : Using `classic'" << std::endl;

      // First set damping as usual
      SetSmoothingDampingClassic();

      agg_->Restriction_smoothagg_transpose = ML_FALSE;

    }  else if( DampingType == "classic-use-A" ) {

      if( verbose_ )
        std::cout << PrintMsg_ << "R and P smoothing : Using `classic-use-A'" <<std::endl;

      field_of_values->choice     =  2;
      field_of_values->poly_order =  2;

      field_of_values->R_coeff[0] =  1.333;
      field_of_values->R_coeff[1] =  0.0;
      field_of_values->R_coeff[2] =  0.0;

      field_of_values->P_coeff[0] =  1.333;
      field_of_values->P_coeff[1] =  0.0;
      field_of_values->P_coeff[2] =  0.0;

      ML_Aggregate_Set_DampingFactor(agg_,0.0);

    } else if( DampingType == "user-defined") {

      if( verbose_ )
    std::cout << PrintMsg_ << "R and P smoothing : Using `user-defined'" << std::endl;

      // user may specify each coefficient. Default values as for fov-1
      field_of_values->choice     =  1;
      field_of_values->poly_order =  2;

      // get them from parameters' list.
      field_of_values->R_coeff[0] =  List_.get("R and P smoothing: c_0",  1.107);
      field_of_values->R_coeff[1] =  List_.get("R and P smoothing: c_1",  0.285);
      field_of_values->R_coeff[2] =  List_.get("R and P smoothing: c_2",  0.718);

      field_of_values->P_coeff[0] =  List_.get("R and P smoothing: g_0",  1.878);
      field_of_values->P_coeff[1] =  List_.get("R and P smoothing: g_1", -2.515);
      field_of_values->P_coeff[2] =  List_.get("R and P smoothing: g_2",  0.942);

      ML_Aggregate_Set_DampingFactor(agg_,0.0);

    } else {

      if( Comm().MyPID() == 0 ) {
    std::cerr << std::endl;
    std::cerr << ErrorMsg_ << "Parameter for `R and P smoothing : damping' not recognized" << std::endl
         << ErrorMsg_ << "It is: `" << DampingType << "'. It should be one of:" << std::endl
         << ErrorMsg_ << "<fov-1> / <fov-2> / <fov-5> / <fov-10> / <user-defined>" << std::endl;
      }

      ML_EXIT(-10); // wrong input parameter
    }

    agg_->field_of_values = (void*) field_of_values;

  } else {

    if( Comm().MyPID() == 0 ) {
      std::cerr << std::endl;
      std::cerr << ErrorMsg_ << "Parameter for `R and P smoothing : type' not recognized" << std::endl
       << ErrorMsg_ << "It is: `" << RandPSmoothing << "'. It should be one of:" << std::endl
       << ErrorMsg_ << "<classic> / <advanced>" << std::endl;
    }

    ML_EXIT(-11); // wrong input parameter

  }

  return 0;
}

// ============================================================================

int ML_Epetra::MultiLevelPreconditioner::SetSmoothingDampingClassic()
{

  double DampingFactor = 1.333;
  if (AMGSolver_ == ML_MAXWELL) DampingFactor = 0.0;

  DampingFactor = List_.get("aggregation: damping factor",
                DampingFactor);
  ML_Aggregate_Set_DampingFactor( agg_, DampingFactor );

  if( verbose_ ) {
    std::cout << PrintMsg_ << "R and P smoothing : P = (I-\\omega A) P_t, R = P^T" << std::endl;
    std::cout << PrintMsg_ << "R and P smoothing : \\omega = " << DampingFactor << "/lambda_max" <<std::endl;
  }

#ifdef no_longer_done_here
  /*********************************************************************/
  /* Changed the code so that the eigen-analysis type is no longer     */
  /* connected to aggregates but is instead connected to matrices. This*/
  /* means that the eigen-analysis stuff is parsed and set elsewhere.  */
  /*********************************************************************/

  std::string str = List_.get("eigen-analysis: type","cg");

  if( verbose_ ) std::cout << PrintMsg_ << "Using `" << str << "' scheme for eigen-computations" << std::endl;

  if( str == "cg" )                ML_Aggregate_Set_SpectralNormScheme_Calc(agg_);
  else if( str == "Anorm" )        ML_Aggregate_Set_SpectralNormScheme_Anorm(agg_);
  else if( str == "Anasazi" )      ML_Aggregate_Set_SpectralNormScheme_Anasazi(agg_);
  else if( str == "power-method" ) ML_Aggregate_Set_SpectralNormScheme_PowerMethod(agg_);
  else {
    if( Comm().MyPID() == 0 ) {
      std::cerr << ErrorMsg_ << "parameter `eigen-analysis: type' has an incorrect value"
       << "(" << str << ")" << std::endl;
      std::cerr << ErrorMsg_ << "It should be: " << std::endl
       << ErrorMsg_ << "<cg> / <Anorm> / <Anasazi> / <power-method>" << std::endl;
    }
    ML_EXIT(-10); // wrong input parameter
  }
#endif

  return 0;
}

// ============================================================================

int ML_Epetra::MultiLevelPreconditioner::
PrintStencil2D(const int nx, const int ny,
           int NodeID,
           const int EquationID)
{

  if (nx <= 0)
    ML_CHK_ERR(-1); // need nodes along the X-axis

  if (ny <= 0)
    ML_CHK_ERR(-2); // need nodes along Y-axis


  if (RowMatrix_ == 0)
    ML_CHK_ERR(-3); // matrix still not set

  // automatically compute NodeID, somewhere in the middle of the grid
  if (NodeID == -1) {
    if (ny == 1)
      NodeID = (int)(nx/2);
    else
      NodeID = (int)(ny*(nx/2) + nx/2);
  }

  // need to convert from NodeID (BlockRowID) to PointRowID
  int GID = NodeID * NumPDEEqns_;

  int LID = RowMatrix_->RowMatrixRowMap().LID(GID);

  // only processor having this node will go on
  if (LID == -1)
    return(0);

  int MaxPerRow = RowMatrix_->MaxNumEntries();
  int NumEntriesRow;   // local entries on each row
  std::vector<double> Values;  Values.resize(MaxPerRow);
  std::vector<int>    Indices; Indices.resize(MaxPerRow);

  int ierr = RowMatrix_->ExtractMyRowCopy(LID, MaxPerRow, NumEntriesRow,
                      &Values[0], &Indices[0]);

  if (ierr)
    ML_CHK_ERR(-4);

  // cycle over nonzero elements, look for elements in positions that we
  // can understand

  Epetra_IntSerialDenseMatrix StencilInd(3,3);
  Epetra_SerialDenseMatrix StencilVal(3,3);
  for( int i=0 ; i<3 ; ++i )
    for( int j=0 ; j<3 ; ++j ) {
      StencilVal(i,j) = 0.0;
    }

  // look for the following positions
  StencilInd(0,0) = RowMatrix_->RowMatrixColMap().LID(NodeID-1-nx);
  StencilInd(1,0) = RowMatrix_->RowMatrixColMap().LID(NodeID-nx);
  StencilInd(2,0) = RowMatrix_->RowMatrixColMap().LID(NodeID+1-nx);
  StencilInd(0,1) = RowMatrix_->RowMatrixColMap().LID(NodeID-1);
  StencilInd(1,1) = RowMatrix_->RowMatrixColMap().LID(NodeID);
  StencilInd(2,1) = RowMatrix_->RowMatrixColMap().LID(NodeID+1);
  StencilInd(0,2) = RowMatrix_->RowMatrixColMap().LID(NodeID-1+nx);
  StencilInd(1,2) = RowMatrix_->RowMatrixColMap().LID(NodeID+nx);
  StencilInd(2,2) = RowMatrix_->RowMatrixColMap().LID(NodeID+1+nx);

  for( int i=0 ; i<NumEntriesRow ; ++i ) {
    // get only the required equation
    if( Indices[i]%NumPDEEqns_ ) continue;
    // convert into block row
    int LocalColID = Indices[i]/NumPDEEqns_;
    // look for known positions
    for( int ix=0 ; ix<3 ; ++ix ) {
      for( int iy=0 ; iy<3 ; ++iy ) {
    if( StencilInd(ix,iy) == LocalColID ) {
      StencilVal(ix,iy) = Values[i];
    }
      }
    }
  }

  std::cout << "2D computational stencil for equation " << EquationID << " at node " << NodeID
       << " (grid is " << nx << " x " << ny << ")" << std::endl;
  std::cout << std::endl;
  for( int iy=0 ; iy<3 ; ++iy ) {
    std::cout << "\t";
    for( int ix=0 ; ix<3 ; ++ix ) {
      std::cout << " " << std::setw(15) << StencilVal(ix,iy);
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  return 0;
}

// ============================================================================

void ML_Epetra::MultiLevelPreconditioner::Apply_BCsToGradient(
             const Epetra_RowMatrix & iEdgeMatrix,
             const Epetra_RowMatrix & iGrad)
{
  /* This function zeros out *rows* of T that correspond to Dirichlet rows in
     the curl-curl matrix.  It mimics what was done previously in
     ML_Tmat_applyDirichletBC().
     Input:
         EdgeMatrix         curl-curl matrix
         Grad               gradient matrix
     Output:
         Grad               gradient matrix with *rows* zeroed out

     Comments: The graph of Grad is unchanged.
  */

  const Epetra_CrsMatrix *EdgeMatrix = dynamic_cast<const Epetra_CrsMatrix*>
                                           (&iEdgeMatrix );
  const Epetra_CrsMatrix *Grad = dynamic_cast<const Epetra_CrsMatrix*>(&iGrad );

  if (EdgeMatrix == 0 || Grad == 0) {
    if (verbose_)
      std::cout << "Not applying Dirichlet boundary conditions to gradient "
           << "because cast failed." << std::endl;
    return;
  }

  // locate Dirichlet edges
  int *dirichletEdges = new int[EdgeMatrix->NumMyRows()];
  int numBCEdges = 0;
  for (int i=0; i<EdgeMatrix->NumMyRows(); i++) {
    int numEntries, *cols;
    double *vals;
    int ierr = EdgeMatrix->ExtractMyRowView(i,numEntries,vals,cols);
    if (ierr == 0) {
      int nz=0;
      for (int j=0; j<numEntries; j++) if (vals[j] != 0.0) nz++;
      if (nz == 1) {
        dirichletEdges[numBCEdges++] = i;
      }
    }
  }
  // -------------------------
  // now zero out the rows
  // -------------------------
  for (int i=0; i < numBCEdges; i++) {
    int numEntries;
    double *vals;
    int *cols;
    Grad->ExtractMyRowView(dirichletEdges[i],numEntries,vals,cols);
    for (int j=0; j < numEntries; j++)
      vals[j] = 0.0;
  }
  delete [] dirichletEdges;
} //Apply_BCsToGradient

#ifdef ML_VERSION_THAT_ZEROS_COLUMNS
// ============================================================================

void ML_Epetra::MultiLevelPreconditioner::Apply_BCsToGradient(
             const Epetra_RowMatrix & iEdgeMatrix,
             const Epetra_RowMatrix & iGrad)
{
  /* This function zeros out columns of T that are no longer in the null space
     of the curl-curl matrix because of Dirichlet boundary conditions.
     Input:
         EdgeMatrix         curl-curl matrix
         Grad               gradient matrix
     Output:
         Grad               gradient matrix with columns zeroed out

     Comments: The graph of Grad is unchanged.
  */

  const Epetra_CrsMatrix *EdgeMatrix = dynamic_cast<const Epetra_CrsMatrix*>
                                           (&iEdgeMatrix );
  const Epetra_CrsMatrix *Grad = dynamic_cast<const Epetra_CrsMatrix*>(&iGrad );

  if (EdgeMatrix == 0 || Grad == 0) {
    if (verbose_)
      std::cout << "Not applying Dirichlet boundary conditions to gradient "
           << "because cast failed." << std::endl;
    return;
  }

  // locate Dirichlet edges
  int *dirichletEdges = new int[EdgeMatrix->NumMyRows()];
  int numBCEdges = 0;
  for (int i=0; i<EdgeMatrix->NumMyRows(); i++) {
    int numEntries, *cols;
    double *vals;
    int ierr = EdgeMatrix->ExtractMyRowView(i,numEntries,vals,cols);
    if (ierr == 0) {
      int nz=0;
      for (int j=0; j<numEntries; j++) if (vals[j] != 0.0) nz++;
      if (nz == 1) {
        dirichletEdges[numBCEdges++] = i;
      }
    }
  }

  const Epetra_Map & ColMap = Grad->ColMap();
  int indexBase = ColMap.IndexBase();
  Epetra_Map globalMap(Grad->NumGlobalCols(),indexBase,EdgeMatrix->Comm());

  // create the exporter from this proc's column map to global 1-1 column map
  Epetra_Export Exporter(ColMap,globalMap);

  // create a std::vector of global column indices that we will export to
  Epetra_Vector globColsToZero(globalMap);
  // create a std::vector of local column indices that we will export from
  Epetra_Vector myColsToZero(ColMap);
  myColsToZero.PutScalar(0);

  // for each local column j in a local dirichlet row, set myColsToZero[j]=1
  for (int i=0; i < numBCEdges; i++) {
    int numEntries;
    double *vals;
    int *cols;
    Grad->ExtractMyRowView(dirichletEdges[i],numEntries,vals,cols);
    for (int j=0; j < numEntries; j++)
      myColsToZero[ cols[j] ] = 1;
  }

  // export to the global column map
  globColsToZero.Export(myColsToZero,Exporter,Add);
  // now import from the global column map to the local column map
  myColsToZero.Import(globColsToZero,Exporter,Insert);

  // -------------------------
  // now zero out the columns
  // -------------------------

  for (int i=0; i < Grad->NumMyRows(); i++) {
    int numEntries;
    double *vals;
    int *cols;
    Grad->ExtractMyRowView(i,numEntries,vals,cols);
    for (int j=0; j < numEntries; j++) {
      if (myColsToZero[ cols[j] ] > 0) {
        vals[j] = 0.0;
      }
    }
  }
  delete [] dirichletEdges;
} //Apply_BCsToGradient
#endif //ifdef ML_VERSION_THAT_ZEROS_COLUMNS

// ============================================================================

void ML_Epetra::MultiLevelPreconditioner::
CheckNullSpace()
{
  // Check that the null space is correct.
  Epetra_Vector XXX(TMatrix_->OperatorDomainMap());
  Epetra_Vector YYY(TMatrix_->OperatorRangeMap());
  Epetra_Vector ZZZ(TMatrix_->OperatorRangeMap());
  XXX.Random();
  double norm;
  XXX.Norm2(&norm);
  XXX.Scale(1.0/norm);
  TMatrix_->Multiply(false,XXX,YYY);
  CurlCurlMatrix_->Multiply(false,YYY,ZZZ);
  ZZZ.Norm2(&norm);
  int mypid = CurlCurlMatrix_->Comm().MyPID();
  if (mypid ==0 && ML_Get_PrintLevel() > 0)
    std::cout << "Checking curl/gradient relationship" << std::endl;
  double norminf = CurlCurlMatrix_->NormInf();
  if (mypid ==0 && ML_Get_PrintLevel() > 0) {
    if (norm > (1e-12 * norminf))
      std::cout << std::endl
           << "**WARNING** ||curlcurl * grad * vrand||_2 = " << norm << std::endl
           << "**WARNING** Either the curl-curl or the null space may be wrong."
           << std::endl << std::endl;
    else {
      std::cout << "||curlcurl||_inf              = " << norminf << std::endl;
      std::cout << "||curlcurl * grad * vrand||_2 = " << norm << std::endl;
    }
  }
} //CheckNullSpace()

// ============================================================================

#ifdef HAVE_ML_EPETRAEXT
Epetra_RowMatrix* ML_Epetra::MultiLevelPreconditioner::
ModifyEpetraMatrixColMap(const Epetra_RowMatrix &A,
                         EpetraExt::CrsMatrix_SolverMap &transform,
                         const char *matrixLabel)
{
    Epetra_RowMatrix *B;
    Epetra_CrsMatrix *Acrs;

    const Epetra_CrsMatrix *Atmp = dynamic_cast<const Epetra_CrsMatrix*>(&A);
    if (Atmp != 0) {
      Acrs = const_cast<Epetra_CrsMatrix*>(Atmp);
      B = &(transform(*Acrs));
    }
    else {
      B = const_cast<Epetra_RowMatrix *>(&A);
      //Some classes may not implement the Importer method, e.g., Ifpack_LocalFilter
      if (A.RowMatrixImporter() != 0) {
        //Even if the matrix cannot be cast to an Epetra_CrsMatrix, we
        //still test for missing column indices and warn if any are found.
        const Epetra_BlockMap & DomainMap = A.RowMatrixImporter()->SourceMap();
        const Epetra_Map & ColMap = A.RowMatrixColMap();
        int NumCols = DomainMap.NumMyElements();
        int Match = 0;
        for( int i = 0; i < NumCols; ++i )
          if( DomainMap.GID(i) != ColMap.GID(i) ) {
            Match = 1;
            break;
          }
        int MatchAll = 0;
        A.Comm().SumAll( &Match, &MatchAll, 1 );
        if( MatchAll ) {
          if (verbose_ && !A.Comm().MyPID())
            printf("WARNING: Matrix %s is missing local column indices, but is not an Epetra_CrsMatrix, so doing nothing.\n",matrixLabel);
        }
        else if (verbose_ && !A.Comm().MyPID())
          printf("** Matrix is not an Epetra_CrsMatrix; column map is not missing any entries.\n");
      } //if (A.RowMatrixImporter() != 0)
    }

    if (verbose_ && !A.Comm().MyPID()) {
      if (B != &A)
        printf("** Transforming column map of %s matrix\n", matrixLabel);
      else
        printf("** Leaving column map of %s matrix unchanged\n",matrixLabel);
    }

    return B;
} //ModifyEpetraMatrixColMap()
#endif

// ============================================================================

void ML_Epetra::MultiLevelPreconditioner::
ReportTime()
{
  ML_print_line("-",78);
  double TotalTime = FirstApplicationTime_ + ApplicationTime_;
  std::cout << PrintMsg_ << "   ML time information (seconds)          total          avg" << std::endl << std::endl
   << PrintMsg_ << "   1- Construction                  = "
       << setw(10) << ConstructionTime_ << "  "
       << setw(10) << ConstructionTime_ / NumConstructions_ << std::endl;
  std::cout << PrintMsg_ << "   2- Preconditioner apply          = "
       << setw(10) << TotalTime << std::endl;
  std::cout << PrintMsg_ << "     a- first application(s) only   = "
       << setw(10) << FirstApplicationTime_ << "  "
       << setw(10) << FirstApplicationTime_ / NumConstructions_
       << std::endl;
  std::cout << PrintMsg_ << "     b- remaining applications      = "
       << setw(10) << ApplicationTime_ << "  "
       << setw(10) << ApplicationTime_ / NumApplications_ << std::endl;
  std::cout << std::endl;
  std::cout << PrintMsg_ << "   3- Total time required by ML so far is "
       << ConstructionTime_ + TotalTime << " seconds" << std::endl;
  std::cout << PrintMsg_ << "      (constr + all applications)" << std::endl;
} //MLP::ReportTime()

// ============================================================================

/*
void ML_Epetra::MultiLevelPreconditioner::
ResetTime()
{
  FirstApplicationTime_ = 0.0;
  ApplicationTime_      = 0.0;
  NumApplications_      = 0;
} //MLP::ResetTime()
*/

#ifdef NewStuff
// All of these functions are used by the multiphysics variable dof per node
// precondiitoner. They should be made as private functions of this class.

int MLcolGlobalIds(struct wrappedCommStruct& framework, MLVec<int>& myGids)
{
   /***************************************************************************
    Gets global column ids from an epetra matrx.
   ***************************************************************************/

   int nLocal;

   if (framework.whichOne == mlType) {
       std::cout << "MLcolGlobalIds not implemented for mlType " << std::endl;
   }
   else {
      Epetra_CrsMatrix *Amat = (Epetra_CrsMatrix *) framework.data;
      nLocal = Amat->ColMap().NumMyPoints();

      myGids.resize(nLocal);
      int *iptr = myGids.getptr();
      Amat->ColMap().MyGlobalElements( iptr );
      return(nLocal);
   }
   return 0;
}



int dofCommUsingMlNodalMatrix(double *data, void *widget)
{
   /***************************************************************************
    Performs communication to update an unamalgamated/dof vector using
    an ML style amalgamated/nodal matrix to define the
    communication/import pattern. It assumes that all nodes have
    maxDofPerNode dofs and repeatedly employs the nodal matrix's
    communication function to update each the 1st dof within all
    nodes, the 2nd dof within all nodes, ...
   ***************************************************************************/

    struct wrappedCommStruct *framework = (struct wrappedCommStruct *) widget;
    int n = framework->vecSize;

    MLVec<double> vector(data,data+n);

    int maxDofPerNode = framework->maxDofPerNode;
    ML_Operator *Amat = (ML_Operator *) framework->data;
    if (Amat == NULL) return 0;
    if (Amat->getrow == NULL) return 0;
    if (Amat->getrow->pre_comm == NULL) return 0;

    MLVec<double> temp(vector.size()/maxDofPerNode);

    for (int j = 0; j < maxDofPerNode; j++) {
       for (int i = 0; i < vector.size()/maxDofPerNode; i++)
          temp[i] = vector[i*maxDofPerNode+j];

       ML_exchange_bdry(temp.getptr(),Amat->getrow->pre_comm,
                 Amat->invec_leng,Amat->comm,ML_OVERWRITE,NULL);

       for (int i = 0; i < vector.size()/maxDofPerNode; i++)
          vector[i*maxDofPerNode+j] = temp[i];
    }



    return 0;
}

int MLnMyGhost(struct wrappedCommStruct& framework)
{
   /***************************************************************************
    Determine the number of ghost nodes based on either an ML matrix or an
    epetra matrix.
   ***************************************************************************/

   if (framework.whichOne == mlType) {
      ML_Operator *Amat = (ML_Operator *) framework.data;
      if (Amat == NULL) return 0;
      if (Amat->getrow == NULL) return 0;
      if (Amat->getrow->pre_comm == NULL) return 0;
      return( ML_CommInfoOP_Compute_TotalRcvLength(Amat->getrow->pre_comm));
   }
   else {
      Epetra_CrsMatrix *Amat = (Epetra_CrsMatrix *) framework.data;
      return( Amat->RowMatrixColMap().NumMyElements() -
              Amat->OperatorDomainMap().NumMyElements());
   }
}

/*!
  @brief Take Crs std::vec's, make arrays and shove them into an ML Operator

  This includes allocating double arrays and a widget to point to them and
  calling all the appropriate ML functions to clean out what used to be in the
  operator before inserting the new matrix.

  @param[in] rowPtr,cols,vals    Crs matrix std:vectors
  @param     Mat                 On input, a populated ML Operator. On output, old ML Operator is cleaned and a new Crs matrix is filled into the operator
*/
int MLShove(ML_Operator *Mat, MLVec<int>& rowPtr, MLVec<int>& cols, MLVec<double>& vals, int invec_leng, int (*commfunc  )(double *vec, void *data), struct wrappedCommStruct& framework, int nGhost)
{
    int outvec_leng = rowPtr.size() - 1;

    struct ML_CSR_MSRdata* ptr = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
    if (ptr == NULL) pr_error("ML_Operator_Add: no space for temp\n");

    ptr->rowptr = rowPtr.getptr();         // getptr() + relinquishData
    rowPtr.relinquishData();               // effectively pulls data
    ptr->columns = cols.getptr();          // out of a MLVec and empties
    cols.relinquishData();                 // MLVec's contents
    ptr->values = vals.getptr();
    vals.relinquishData();

    ML_Comm * temp_comm = Mat->comm;

    // Need to create new communication widget before cleaning out Mat
    // as Mat may be used within commfunc.

    ML_CommInfoOP *newPrecComm = NULL;
    if (commfunc != NULL) {
        ML_CommInfoOP_Generate( &newPrecComm, commfunc, (void *) &framework,
                              Mat->comm, invec_leng, nGhost);
    }
    ML_Operator_Clean(Mat);
    ML_Operator_Init(Mat, temp_comm);
    ML_Operator_Set_ApplyFuncData(Mat, invec_leng, outvec_leng,
                                  ptr, outvec_leng, NULL,0);
    ML_Operator_Set_Getrow(Mat, outvec_leng, CSR_getrow);
    ML_Operator_Set_ApplyFunc(Mat, CSR_matvec);
    Mat->data_destroy = ML_CSR_MSRdata_Destroy;
    Mat->N_nonzeros     = cols.size();

    Mat->getrow->pre_comm = newPrecComm;

  return(0);
}
/******************************************************************************
 *****************************************************************************/
int MLextractDiag(const MLVec<int>& rowPtr, const MLVec<int>& cols, const MLVec<double>& vals, MLVec<double>& diagonal, struct wrappedCommStruct &framework)
{
/******************************************************************************
 * extract matrix diagonal and store in diagonal.
 *****************************************************************************/

   for (int i = 0; i < (int) rowPtr.size()-1; i++) {
      for (int j = rowPtr[i]; j < rowPtr[i+1]; j++) {
         if (i == cols[j])  diagonal[i] = vals[j];
      }
   }
   dofComm(diagonal, framework);
   return(0);
}
/******************************************************************************
 *  *****************************************************************************/
int MLfindDirichlets(const MLVec<int>& rowPtr, const MLVec<int>& cols, const MLVec<double>& vals, const MLVec<double>& diagonal, double tol, MLVec<bool>& dirOrNot, struct wrappedCommStruct &framework)
{
/******************************************************************************
 * Look at each matrix row and mark it as Dirichlet if there is only one
 * "not small" nonzero on the diagonal. In determining whether a nonzero is
 * "not small" use
 *           abs(A(i,j))/sqrt(abs(diag[i]*diag[j])) > tol
 *
 * On output,  dirOrNot is set to to true for Dirichlet rows, false for non-
 * Dirichlet rows
 *****************************************************************************/
  int zeroOnDiag = 0,hasDiag, count;

  for (int i = 0; i < (int) rowPtr.size()-1; i++) {
     dirOrNot[i] = false;
     count   = 0;
     hasDiag = 0;
     for (int j = rowPtr[i]; j < rowPtr[i+1]; j++) {
        if (i != cols[j]) {
           if (fabs(vals[j]/sqrt(fabs(diagonal[i]*diagonal[cols[j]]))) >tol)
              count++;
        }
        else hasDiag = 1;
      }
      if (hasDiag == 0)    { zeroOnDiag = 1; }
      else if (count == 0) dirOrNot[i] = true;
   }
   dofComm(dirOrNot, framework);

   return(zeroOnDiag);
}
/******************************************************************************
 *  *****************************************************************************/
int MLrmDirichletCols(MLVec<int>& rowPtr, MLVec<int>& cols, MLVec<double>& vals, const MLVec<double>& diagonal, bool squeeze, MLVec<double>& solution, MLVec<double>& rhs, const MLVec<bool>& dirOrNot, struct wrappedCommStruct &framework)
{
/******************************************************************************
 *  Remove any matrix col entries associated Dirichlet boundary conditions.
 *  In doing this the right hand side must be properly adjusted to take into
 *  account that the Dirichlet bcs are no longer in the matrix. Additionally,
 *  the solution at Dirichlet points is also set if a solution is provided.
 *  When removing matrix entries, we can either set a zero in vals[k] or
 *  if squeeze == 1 we can actually really remove the entries and change
 *  rowPtr to reflect the removal.
 *
 *  On output:  rowPtr, cols, vals, rhs, and possible solution
 *    are updated to reflect the removed entries.
 *
 *  ******************************************************************************/

  if (!solution.empty()) {
    for (int i = 0; i < (int) rowPtr.size()-1; i++)
       if (dirOrNot[i]) solution[i] = rhs[i]/diagonal[i];
  }
  MLVec<double> rhsCopy(dirOrNot.size());
  for (int i = 0; i < (int) rowPtr.size()-1; i++) rhsCopy[i] = rhs[i];
  dofComm(rhsCopy, framework);

  for (int i = 0; i < (int) rowPtr.size()-1; i++) {
     for (int j = rowPtr[i]; j < rowPtr[i+1]; j++) {
        if (i != cols[j]) {
            if (dirOrNot[cols[j]]) {
               rhs[i] = rhs[i] - vals[j]*rhsCopy[cols[j]]/diagonal[cols[j]];
               vals[j] = 0.0;
            }
        }
     }
  }
  /* squeeze out zeros */

  MLVec<bool> empty;
  if (squeeze) MLsqueezeOutNnzs(rowPtr, cols, vals, empty);

  return(0);
}
/******************************************************************************
 *  *****************************************************************************/
int MLsqueezeOutNnzs(MLVec<int>& rowPtr, MLVec<int>& cols, MLVec<double>& vals, const MLVec<bool>& keep)
{
/******************************************************************************
 * Get rid of nonzero entries that have 0's in them and properly change
 * rowPtr to reflect this removal (either vals == NULL & vals != NULL
 * or the contrary)
 *
 *****************************************************************************/
  int count, newStart;

  count = 0;
  int nRows = rowPtr.size()-1;

  if (!vals.empty()) {
     for (int i = 0; i < nRows; i++) {
        newStart = count;
        for (int j = rowPtr[i]; j < rowPtr[i+1]; j++) {
           if (vals[j] != 0.0) {
              cols[count  ] = cols[j];
              vals[count++] = vals[j];
           }
        }
        rowPtr[i] = newStart;
     }
   }
   else {
     for (int i = 0; i < nRows; i++) {
        newStart = count;
        for (int j = rowPtr[i]; j < rowPtr[i+1]; j++) {
           if (keep[j]) {
              cols[count++] = cols[j];
           }
        }
        rowPtr[i] = newStart;
     }
   }

   rowPtr[nRows] = count;

  return(0);
}
/******************************************************************************
 *****************************************************************************/
int MLbuildMap(const MLVec<bool>& dofPresent, MLVec<int>& map, int nDofs)
{
/******************************************************************************
 *  Compute map that maps dofs in a variable dof matrix to dofs in a padded
 *  matrix. Specifically, on input:
 *     dofPresent[ i*maxDofPerNode+j] indicates whether or not the jth dof
 *                                    at the ith node is present in the
 *                                    variable dof matrix (e.g., the ith node
 *                                    has an air pressure dof). true means the
 *                                    dof is present while false means it is not
 *
 *
 * On output:
 *    map[k]                          indicates that the kth dof in the variable
 *                                    dof matrix would correspond to the
 *                                    map[k]th dof in a padded system. This
 *                                    padded system is not ever created but
 *                                    would be the associated matrix if every
 *                                    node had maxDofPerNode dofs.
 *
 *    nDofs                           length of map that should equal the
 *                                    dimension of the variable dof matrix.
 *                                    It is returned here to double check this.
 *  NOTE: in parallel we would want to do something so that map[j] also
 *  has padded dofs for j's associated with ghost information.
 *****************************************************************************/
   int count;

   count = 0;
   for (int i = 0; i < (int) dofPresent.size(); i++) {
     if (dofPresent[i]) map[count++] = i;
   }

   TEUCHOS_TEST_FOR_EXCEPTION(nDofs != count,std::logic_error,
      "MLbuildMap error, # dofs in dofPresent does not match expected value:"
      << nDofs << " vs. " << count << "\n" );

   return(0);
}


int MLassignGhostLocalNodeIds(MLVec<int> &myLocalNodeIds, int nLocalDofs, int nLocalPlusGhostDofs, struct wrappedCommStruct &framework, int &nLocalNodes, int &nLocalPlusGhostNodes)
{
/******************************************************************************
    Assign the ghost portion of myLocalNodeIds[], which will be the local
    amalgamated node id associated with the kth local non-amalgamated dof.
    It is assumed that the non-ghost portion is properly set.

    For ghosts, we assign a "unique" number from nLocalNodes to
    nLocalNodes+nGhostNodes-1 (which is actually nLocalPlusGhostNodes-1
    in the code). Here, nLocalNodes refers to the number of mesh nodes owned
    by the processor. nGhostNodes is the total number of ghost mesh nodes
    associated with all of its ghost dofs. When I say "unique", I mean that
    I get to pick any local numbering of the ghost nodes so long as the
    following is true:
       a) no two ghost dofs that are associated with different
          mesh nodes have the same nodal id.
       b) all ghost nodes that are associated with a processor
          have consecutive numbering

    This is done is as follows:

       a) communicate so that for each ghost dof, a processor knows the local
          node id on the owning proc and knows the processor that owns it

       b) for each owning proc ...
            i) copy the owning proc local node id into a new array, and record
               as well where this information came from. That is,

                  tempId[  j] = myLocalNodeIds[i];
                  location[j] = i;

            ii) sort tempId to make it easier to determine whether two ghosts
                belong to the same node. Also shuffle location[] in the same
                way as tempId[] so tempId[j] corresonds to myLocalNodeIds[ location[j] ].

       c) go through tempId[] and assign a unique ghost node id to tempId[k] if
          it has a different owning proc or a different local node id on the
          owning proc. That is, increment the number of ghost node ids when
          tempId[k] is a different node from tempId[k-1]
*****************************************************************************/

   MLVec<int> myProc(nLocalPlusGhostDofs);
   for (int i = 0; i < nLocalDofs; i++) myProc[i] = framework.myPid;

   dofComm(myLocalNodeIds, framework);
   dofComm(myProc, framework);

   // At this point, the ghost part of myLocalNodeIds corresponds to the
   // local ids associated with the owning processor. We want to convert
   // these to local ids associated with the processor on which these are
   // ghosts. Thus, we have to re-number them. In doing this re-numbering
   // we must make sure that we find all ghosts with the same id & proc
   // and assign a unique local id to this group. To do this find, we sort
   // sort all ghost myLocalNodeIds that are owned by the same processor.
   // Then we can look for duplicates (i.e. several ghost entries
   // corresponding to dofs with the same node id) easily and make sure these
   // are all assigned the same local id. To do the sorting
   // we'll make a temporary copy of the ghosts via tempId and tempProc and
   // sort this multiple times for each group owned by the same processor

   int* location= (int *) ML_allocate(sizeof(int)*(nLocalPlusGhostDofs - nLocalDofs + 1));
   int* tempId  = (int *) ML_allocate(sizeof(int)*(nLocalPlusGhostDofs - nLocalDofs + 1));
   int* tempProc= (int *) ML_allocate(sizeof(int)*(nLocalPlusGhostDofs - nLocalDofs + 1));

   TEUCHOS_TEST_FOR_EXCEPTION(tempProc == NULL,std::logic_error,
   "poo error, not enough space for vector of size " << nLocalPlusGhostDofs-nLocalDofs << "\n");

   int notProcessed = nLocalDofs;
   int tempIndex = 0;
   int first, neighbor;

   while ( notProcessed < nLocalPlusGhostDofs ) {

     neighbor = myProc[notProcessed];
     first = tempIndex;
     location[tempIndex  ] = notProcessed;
     tempId[  tempIndex++] = myLocalNodeIds[notProcessed];
     myProc[notProcessed]  = -1 - neighbor; // temporarily mark by negative # to
                                           // indicate this guy has been visited

     for (int i = notProcessed+1; i < nLocalPlusGhostDofs; i++) {
        if ( myProc[i] == neighbor ) {
           location[tempIndex  ] = i;
           tempId[  tempIndex++] = myLocalNodeIds[i];
           myProc[     i       ]= -1; // mark as visited
        }
     }
     ML_az_sort(&(tempId[first]), tempIndex-first, &(location[first]), NULL);

     for (int i = first; i < tempIndex; i++) tempProc[i] = neighbor;

     // find next notProcessed corresponding to first non-visited guy

     notProcessed++;
     while ( (notProcessed < nLocalPlusGhostDofs) && (myProc[notProcessed] < 0) )
       notProcessed++;
   }
   TEUCHOS_TEST_FOR_EXCEPTION(tempIndex != nLocalPlusGhostDofs-nLocalDofs,std::logic_error,
   "poo error, number of nonzero ghosts is not consistent\n");


   // now assign ids to all ghost nodes (giving the same id to those with
   // the same myProc[] and the same local id (on the proc
   // that actually owns the variable associated with the ghost)

   int lagged = -1;
   nLocalNodes = 0;

   if (nLocalDofs > 0) nLocalNodes = myLocalNodeIds[nLocalDofs-1] + 1;
   nLocalPlusGhostNodes = nLocalNodes;

   if (nLocalDofs < nLocalPlusGhostDofs) nLocalPlusGhostNodes++; // 1st ghost node is unique
                                               // (not already accounted for)

   // check if two adjacent ghost dofs correspond to different nodes. To do
   // this, check if they are from different processors or whether they have
   // different local node ids

   for (int i = nLocalDofs+1; i < nLocalPlusGhostDofs; i++) {

      lagged = nLocalPlusGhostNodes-1;

      // i is a new unique ghost node (not already accounted for)
      if ((tempId[i-nLocalDofs] != tempId[i-1-nLocalDofs]) ||
          (tempProc[i-nLocalDofs]         != tempProc[i-1-nLocalDofs]))
         nLocalPlusGhostNodes++;

      tempId[i-1-nLocalDofs] = lagged;
   }
   if (nLocalPlusGhostDofs > nLocalDofs)
      tempId[nLocalPlusGhostDofs-1-nLocalDofs] = nLocalPlusGhostNodes-1;


   // now copy the nodal ids to the proper spot in myLocalNodeIds[]

   for (int i = nLocalDofs; i < nLocalPlusGhostDofs; i++)
      myLocalNodeIds[location[i-nLocalDofs]] = tempId[i-nLocalDofs];

   ML_free(location);
   ML_free(tempProc);
   ML_free(tempId);

   return 0;
}

int MLfillNodalMaps(MLVec<int> &amalgRowMap, MLVec<int> &amalgColMap,
        MLVec<int> &myLocalNodeIds, int nLocalDofs,
        struct wrappedCommStruct &framework,
        int nLocalNodes, int nLocalPlusGhostNodes)
{
/*
    Fill amalgRowMap and amalgColMap.

    amalgRowMap is easy. It corresponds to the global id of the 1st dof of each node
    (owned by this processor). It is easy because it doesn't involve communication.

    amalgColMap is essentially the colmap version of amalgRowMap. We basically
    assign the local portion of amalgColMap in the same way as amalgRowMap. To
    get the ghost part of amalgColMap, we use nodalComm(). This routine uses
    myLocalNodeIds[] in conjunction with a function to perform amalgamated dof
    communication.
*/

   MLVec<int> myGids;

   MLcolGlobalIds(framework, myGids);

   amalgRowMap.resize(nLocalNodes);
   amalgColMap.resize(nLocalPlusGhostNodes);

   int count = 0;
   if (nLocalDofs > 0) {
      amalgRowMap[count] = myGids[0];
      amalgColMap[count] = myGids[0];
      count++;
   }
   // dofs belonging to same node must be consecutive only for the local
   // part of myLocalNodeIds[]

   for (int i = 1; i < nLocalDofs; i++) {
      if (myLocalNodeIds[i] != myLocalNodeIds[i-1]) {
         amalgRowMap[count] = myGids[i];
         amalgColMap[count] = myGids[i];
         count++;
      }
   }
   nodalComm(amalgColMap,  myLocalNodeIds, framework);

   return 0;
}

/******************************************************************************
 *****************************************************************************/
int MLvariableDofAmalg(int nCols, const MLVec<int>& rowPtr,
                     const MLVec<int>& cols, const MLVec<double>& vals,
                     int nNodes, int maxDofPerNode, const MLVec<int>& map,
                     const MLVec<double>& diag, double tol,
                     MLVec<int>& amalgRowPtr, MLVec<int>& amalgCols,
                     struct wrappedCommStruct &/* framework */,
                     MLVec<int>& myLocalNodeIds)
{
/******************************************************************************
 *  Amalgmate crs matrix (rowPtr,cols,vals) and store result in
 *  amalgRowPtr, amalgCols.  Optionally, small vals in the non-amalgamated
 *  matrix can be dropped when performing this amalgamation  when
 *
 *           abs(A(i,j))/sqrt(abs(diag[i]*diag[j])) < tol
 *
 *  Here, the local non-amalgamated matrix is of size nRows x nCols and the
 *  resulting amalgamated matrix should have nNodes local rows.  On input,
 *  map(i) indicates that the ith dof in the non-amalgamated matrix corresponds
 *  to the map(i)th dof in a padded version of the non-amalgamated matrix.
 *  A padded version of the matrix would have maxDofPerNode dofs at every grid
 *  node.
 *****************************************************************************/
   int  blockRow = 0, blockColumn, newNzs, oldBlockRow;
   int  nNonZeros, doNotDrop;
   int  nLocal = rowPtr.size()-1;



   amalgRowPtr.resize(nNodes+1);
   amalgCols.resize(rowPtr[nLocal]);
   nNonZeros   = 0;
   MLVec<bool> isNonZero(nCols);
   MLVec<int> nonZeroList(nCols);
   for (int i = 0; i < nCols; i++) isNonZero[i] = false;

   oldBlockRow = 0;
   newNzs = 0;
   amalgRowPtr[0] = newNzs;

   doNotDrop = 0;
   if  (tol ==     0.0) doNotDrop = 1;
   if  (vals.empty())  doNotDrop = 1;

   for (int i = 0; i < nLocal; i++ ) {
      blockRow = (int) floor( map[i]/maxDofPerNode);
      if (blockRow != oldBlockRow) {
         /* zero out info recording nonzeros in oldBlockRow */

         for (int j = 0; j < nNonZeros; j++) isNonZero[nonZeroList[j]] = false;
         nNonZeros = 0;

        amalgRowPtr[blockRow] = newNzs;  /* record start of next row */
      }

      for (int j = rowPtr[i]; j < rowPtr[i+1]; j++) {
        if (doNotDrop ||
            ( fabs(vals[j]/sqrt(fabs(diag[i]*diag[cols[j]]))) >= tol)) {
           blockColumn = myLocalNodeIds[cols[j]];
           if  (isNonZero[blockColumn] == false) {
              isNonZero[blockColumn]    = true;
              nonZeroList[nNonZeros++ ] = blockColumn;
              amalgCols[newNzs++] = blockColumn;
           }
        }
      }
      oldBlockRow = blockRow;
   }
   amalgRowPtr[blockRow+1] = newNzs;

   TEUCHOS_TEST_FOR_EXCEPTION((blockRow+1 != nNodes)&&(nNodes != 0),std::logic_error,
      "MLvariableDofAmalg error, computed # block rows (" << blockRow+1 <<
      ") != nNodes (" << nNodes << ")\n");

   amalgCols.resize(amalgRowPtr[nNodes]);

   return(0);
}
/******************************************************************************
 *****************************************************************************/
int MLrmDifferentDofsCrossings(const MLVec<bool>& dofPresent, int maxDofPerNode, MLVec<int>& rowPtr, MLVec<int>& cols, int nCols, struct wrappedCommStruct& framework, MLVec<int> &myLocalNodeIds)
{
/******************************************************************************
 *  Remove matrix entries (i,j) where the ith node and the jth node
 *  have different dofs that are 'present'.
 *
 *  Specifically, on input:
 *    dofPresent[ i*maxDofPerNode+k] indicates whether or not the kth dof
 *                                   at the ith node is present in the
 *                                   variable dof matrix (e.g., the ith node
 *                                   has an air pressure dof). true means the
 *                                   dof is present while false means it is not
 *
 *  We create a unique id for the ith node (i.e. uniqueId[i]) via
 *      sum_{k=0 to maxDofPerNode-1} dofPresent[i*maxDofPerNode+k]*2^k
 *  and use this unique idea to remove entries (i,j) when
 *  uniqueId[i] != uniqueId[j]
 *
 *****************************************************************************/
   int temp, ii;
   int nRows = rowPtr.size()-1;


   MLVec<int> uniqueId(nCols);
   MLVec<bool> keep(rowPtr[nRows]);

   for (int i = 0; i < rowPtr[nRows]; i++) keep[i] = true;
   ii = 0;
   for (int i = 0; i < nRows; i++) {
       temp = 1;
       uniqueId[i] = 0;
       for (int j = 0; j < maxDofPerNode; j++) {
          if (dofPresent[ii++]) uniqueId[i] += temp;
          temp = temp*2;
       }
   }
   nodalComm(uniqueId, myLocalNodeIds, framework);

   for (int i = 0; i < nRows; i++) {
      for (int j = rowPtr[i]; j < rowPtr[i+1]; j++) {
         if  (uniqueId[i] !=  uniqueId[cols[j]])  keep[j] = false;
      }
   }

   MLVec<double>   empty;
   MLsqueezeOutNnzs(rowPtr, cols, empty, keep);

   return(0);
}


/******************************************************************************
 *****************************************************************************/
int MLbuildLaplacian(const MLVec<int>& rowPtr, const MLVec<int>& cols, MLVec<double>& vals, const MLVec<double>& x, const MLVec<double>& y, const MLVec<double>& z)
{
/******************************************************************************
 *  Given an inputMatrix (rowPtr,cols) Build a Laplacian matrix defined by
 *
 *  L(i,j)  = 0                     if  inputMatrix(i,j) = 0
 *  L(i,j)  = 1/nodalDistance(i,j)  otherwise
 *
 *  RST: MIGHT WANT TO DO SOMETHING WITH ALL DIRICHLETS SO THAT ROWSUMS THAT
 *  WERE PREVIOUSLY NONZERO REMAIN THAT WAY?
 *****************************************************************************/
   int    diag;
   double sum;

   if ( z.empty() ) {

      for (int i= 0 ; i < (int) rowPtr.size()-1; i++) {
         sum = 0.0;
         diag = -1;
         for (int j= rowPtr[i] ; j < rowPtr[i+1]; j++) {
            if ( cols[j] != i ) {
               vals[j] = sqrt( (x[i]-x[cols[j]])*(x[i]-x[cols[j]]) +
                                 (y[i]-y[cols[j]])*(y[i]-y[cols[j]]));

               TEUCHOS_TEST_FOR_EXCEPTION(vals[j] == 0.0,std::logic_error,
               "MLbuildLaplacian error, " << i << " and " << cols[j] <<
               " have same coordinates: " << x[i] << ", " << y[i] << "\n");

               vals[j] = -1./vals[j];
               sum = sum - vals[j];
            }
            else diag = j;
         }
         if (sum == 0.0) sum = 1.0;

         TEUCHOS_TEST_FOR_EXCEPTION(diag == -1,std::logic_error,
               "MLbuildLaplacian error, row " << i << " has zero diagonal\n");

         vals[diag] = sum;
      }
   }
   else {
      for (int i= 0 ; i < (int) rowPtr.size()-1; i++) {
         sum = 0.0;
         diag = -1;
         for (int j= rowPtr[i] ; j < rowPtr[i+1]; j++) {
            if ( cols[j] != i ) {
               vals[j] = sqrt( (x[i]-x[cols[j]])*(x[i]-x[cols[j]]) +
                                 (y[i]-y[cols[j]])*(y[i]-y[cols[j]]) +
                                 (z[i]-z[cols[j]])*(z[i]-z[cols[j]]));
               TEUCHOS_TEST_FOR_EXCEPTION(vals[j] == 0.0,std::logic_error,
               "MLbuildLaplacian error, " << i << " and " << cols[j] <<
               " have same coordinates\n");
               vals[j] = -1./vals[j];
               sum = sum - vals[j];
            }
            else diag = j;
         }
         if (sum == 0.0) sum = 1.0;
         TEUCHOS_TEST_FOR_EXCEPTION(diag == -1,std::logic_error,
               "MLbuildLaplacian error, row " << i << " has zero diagonal\n");
         vals[diag] = sum;
      }
   }
   return(0);
}


/*! @brief Unamalgamate prolongator so that it is suitable for PDE systems
 *

Unamalgamate Pmat (represented as a CRS matrix via amalgRowPtr, amalgCols, amalgVals
for two different situations:
  Case 1:  Fine level matrix is padded
           In this case, we basically replicate the unamalgamated
           operator, except that padded dofs and Dirichlet points use
           injection. Thus, PUnAmalg is equivalent to something like
                       (  Pmat      0       0  )
                       (    0     Pmat      0  )
                       (    0       0     Pmat )
           where entries associated with fine level padded dofs or Dir. BCs
           are replaced by rows with a single nonzero (whose
           value equals one) and the entire matrix is instead ordered
           so that dofs associated with nodes are consecutive.

  Case 2:  Fine level matrix is not padded.
           Thus, PUnAmalg is equivalent to
                       (  Pmat      0       0  )
                       (    0     Pmat      0  )
                       (    0       0     Pmat )
           where rows associated with dofs not present on finest level (e.g.,
           an air pressure dof within water region) are removed and rows
           associated with Dir. BCs are replaced by rows with a single
           nonzero (whose value equals one) and the entire matrix is
           instead ordered so dofs associated with nodes are consecutive.

In both of the above cases, the coarse level discretization matrix is assumed to be padded.

@params amalgNRows local number of rows for amalgamted prolongator
@params amalgRowPtr  CRS local row pointer for amalgamted prolongator
@params amalgCols  CRS local cols for amalgamted prolongator
@params amalgVals  CRS vals for amalgamted prolongator
@params maxDofPerNode maximum number of degrees-of-freedom at any mesh node
@params status    status[i*maxDofPerNode+j] refers to the jth dof at the ith node.  status()==s ==> standard element: present in fine operator and not a Dirichlet BC. status()==d ==> element corresponds to Dirichlet BC.  status()==p ==> element not present or is padded dof.
@params fineIsPadded Indicates whether fine grid matrix includes padded dofs for those not really present in the PDE system
@params rowPtr   CRS local row pointer for resulting unamalgamated prolongator
@params cols   CRS local cols for resulting unamalgamated prolongator
@params vals   CRS vals for resulting unamalgamated prolongator
*/

int MLunamalgP(const MLVec<int>& amalgRowPtr, const MLVec<int>& amalgCols, const MLVec<double>& amalgVals,
 int maxDofPerNode, const MLVec<char>& status, bool fineIsPadded, MLVec<int>& rowPtr, MLVec<int>& cols, MLVec<double>& vals)
{

   int paddedNrows;
   int Cnt, rowLength, rowCount = 0;

   paddedNrows = (amalgRowPtr.size()-1)*maxDofPerNode;

   if (fineIsPadded) {
      /* Build Pmat for padded fine level matrices. Note: padded fine level */
      /* dofs are transfered by injection. That is, these interpolation     */
      /* stencils do not take averages of coarse level variables. Further,  */
      /* fine level Dirichlet points also use injection.                    */

      Cnt  = 0;
      for (int i=0; i < (int) amalgRowPtr.size()-1; i++) {
       rowLength = amalgRowPtr[i+1] -  amalgRowPtr[i];
       for (int j = 0; j < maxDofPerNode; j++) {
          rowPtr[i*maxDofPerNode+j] = Cnt;
          if (status[i*maxDofPerNode+j] == 's') {
             for (int k = 0; k < rowLength; k++) {
                cols[Cnt  ] = amalgCols[k+amalgRowPtr[i]]*maxDofPerNode+j;
                vals[Cnt++] = amalgVals[k+amalgRowPtr[i]];
             }
          }
       }

      }
      rowPtr[paddedNrows] = Cnt;
      rowCount = paddedNrows;

   }
   else {
      /* Build Pmat for non-padded fine level matrices.  Need to map from   */
      /* non-padded dofs to padded dofs. For this, look at status array and */
      /* skip padded dofs.                                                  */

     Cnt = 0;
     for (int i=0; i < (int) amalgRowPtr.size()-1; i++) {
        rowLength = amalgRowPtr[i+1] -  amalgRowPtr[i];

        for (int j = 0; j < maxDofPerNode; j++) {
           /* no interpolation for padded fine nodes (as these don't exist) */

          if (status[i*maxDofPerNode+j] == 's') {
             rowPtr[rowCount++] = Cnt;
             for (int k = 0; k < rowLength; k++) {
                cols[Cnt  ] = amalgCols[k+amalgRowPtr[i]]*maxDofPerNode+j;
                vals[Cnt++] = amalgVals[k+amalgRowPtr[i]];
             }
          }
          if (status[i*maxDofPerNode+j] == 'd')
             rowPtr[rowCount++] = Cnt;
        }
     }
     rowPtr[rowCount] = Cnt;

   }


   return(rowCount);
}

int MLfindEmptyRows(const MLVec<int>& rowPtr, const MLVec<int>& /* cols */, const MLVec<double>& vals, MLVec<bool>& rowEmpty)
{
/******************************************************************************
  Find rows that have no nonzero entries.

  Padded matrices might give rise to empty rows (e.g., after RAP) that we want
  to replace with Dirichlet rows (where the Dirichlet point has no connections
  to the rest of the matrix).

 *****************************************************************************/

   rowEmpty.resize(rowPtr.size()-1);

   for (int i = 0; i < (int) rowPtr.size()-1; i++) {
      rowEmpty[i] = true;
      for (int j = rowPtr[i]; j < rowPtr[i+1]; j++) {
         if (vals[j] != 0.0) {
            rowEmpty[i] = false;
            break;
         }
      }
   }
   return(0);
}
int MLreplaceEmptyByDirichlet(MLVec<int>& rowPtr, MLVec<int>& cols, MLVec<double>& vals, const MLVec<bool>& rowEmpty)
{
/*****************************************************************************

  Take any empty rows and stuff a Dirichlet condition into the row.

  Padded matrices might give rise to empty rows (e.g., after RAP) that we want
  to replace with Dirichlet rows (where the Dirichlet point has no connections
  to the rest of the matrix).

 *****************************************************************************/
   int nEmpties, rowStart, rowEnd, lastMoved;
   int nRows = rowPtr.size()-1;


   /* count empties */

   nEmpties = 0;

   for (int i = 0; i < nRows; i++)
      if (rowEmpty[i]) nEmpties++;

   /* allocate extra space for empties */

   cols.resize(rowPtr[nRows]+nEmpties);
   vals.resize(rowPtr[nRows]+nEmpties);

   lastMoved   = rowPtr[nRows]+nEmpties;
   rowEnd = rowPtr[nRows]-1;
   rowPtr[nRows] = lastMoved;
   for (int i = nRows ; i > 0; i--) {
      rowStart = rowPtr[i-1];
      /* first copy row */
      for (int j = rowEnd; j >= rowStart; j--) {
           lastMoved--;
           cols[lastMoved] = cols[j];
           vals[lastMoved] = vals[j];
      }
      /* now see if we need to add someone */
      if (rowEmpty[i-1]) {
         lastMoved--;
         cols[lastMoved] = i-1;
         vals[ lastMoved] = 1.0;
      }
      rowPtr[i-1] = lastMoved;
      rowEnd    = rowStart-1;
   }
   TEUCHOS_TEST_FOR_EXCEPTION(lastMoved != 0,std::logic_error,
               "MLreplaceEmptyByDirichlet error, number of empties"
               " seems to be wrong\n");
   return(0);
}

int MLfineStatus(const MLVec<bool>& dofPresent, const MLVec<int>& map, const MLVec<bool>& dirOrNot, MLVec<char>& status)
{
/*****************************************************************************
  Fill the status array on the finest level based on the information in both
  dofPresent and dirOrNot.  status[i*maxDofPerNode+j] refers to the jth dof
  at the ith node.  status()==s ==> standard element: present in fine operator
                                    and not a Dirichlet BC.
                    status()==d ==> element corresponds to Dirichlet BC.
                    status()==p ==> element not present (or on coarse levels
                                                         is padded)
 ****************************************************************************/

   for (int i = 0; i < (int) status.size(); i++) status[i] = 's';
   for (int i = 0; i < (int) status.size(); i++) {
      if (dofPresent[i] == false) status[i] = 'p';
   }
   if (!dirOrNot.empty()) {
      for (int i = 0; i < (int) map.size(); i++) {
         if (dirOrNot[i]) status[map[i]] = 'd';
      }
   }
   return(0);
}
int MLcoarseStatus(const MLVec<bool>& rowEmpty, const MLVec<bool>& dirOrNot, MLVec<char>& status)
{
/*****************************************************************************
  Fill the status array on a coarse level based on the information in dirOrNot.
  status[i*maxDofPerNode+j] refers to the jth dof at the ith node.

                    status()==s ==> standard element: present in fine operator
                                    and not a Dirichlet BC.
                    status()==p ==> element is padded dof.

  The main difference with MLfineStatus() is that it is assumed that the
  dofsPerNode at all nodes is equal to maxDofsPerNode (that is padding
  has already been used so that this is true) and that no Dirichlet
  points remain.
 ****************************************************************************/

   status.resize(rowEmpty.size());
   for (int i = 0; i < (int) rowEmpty.size(); i++) {
      status[i] = 's';
      if (rowEmpty[i]) status[i] = 'p';
   }
   if (!dirOrNot.empty()) {
      for (int i = 0; i < (int) rowEmpty.size(); i++) {
         if (dirOrNot[i]) status[i] = 'p';
      }
   }

   return(0);
}

int MLsortCols(const MLVec<int>& ARowPtr, MLVec<int>& ACols, MLVec<double>& AVals)
{
   int i,j;

   for (i = 0; i < (int) ARowPtr.size()-1; i++) {
      j = ARowPtr[i];
      ML_az_sort(&(ACols[j]), ARowPtr[i+1]-j, NULL, &(AVals[j]));
   }
   return(0);
}

/****************************************************************************/
/* Given arrays denoting pairwise vertices that are shared, make a set of
 * linked lists that groups all shared vertices together. For example,
 * consider
 *     rowZeros[0] = i; colZeros[0] = j;  // indicates i and j are shared
 *     rowZeros[1] = j; colZeros[1] = k;  // indicates j and k are shared
 *     rowZeros[2] = k; colZeros[2] = m;  // indicates k and m are shared
 *
 * Create a linked list indicated that i,j,k, and m are all shared. More
 * generally, we might have several linked lists indicating different
 * sets of shared vertices. On output groupHead[] and groupNext[] encode
 * all of the linked lists. groupHead[] and groupNext[] are described in
 * the comments to buildCompressedA().
 */
/****************************************************************************/
int MergeShared(MLVec<int>& /* cols */, MLVec<int>& rowZeros, MLVec<int>& colZeros,
                MLVec<int>& groupHead, MLVec<int>& groupNext)
{
   int vertOne, vertTwo, vertOneHead;
   int secondGuy, currentListTwo;
   int prior;

   for (int i = 0; i < groupHead.size(); i++) groupHead[i] = -1;
   for (int i = 0; i < groupNext.size(); i++) groupNext[i] = -1;


   for (int i = 0; i < rowZeros.size(); i++) {
      vertOne = rowZeros[i];
      vertTwo = colZeros[i];

      if (groupHead[vertOne] == -1) {
         if (groupHead[vertTwo] != -1) {
           // group already exists, so add vertOne to the list
           // of vertices already shared with vertTwo
            groupHead[vertOne] = groupHead[vertTwo];
            groupNext[vertOne] = groupNext[groupHead[vertTwo]];
            groupNext[groupHead[vertTwo]] = vertOne;
         }
         else { // new group found
            groupHead[vertOne] = vertOne;
            groupHead[vertTwo] = vertOne;
            groupNext[vertOne] = vertTwo;
         }
      }
      else {
         if (groupHead[vertTwo] == -1) {
            // group already exists, so we need to add vertTwo to the list
            // of vertices already shared with vertOne
            groupHead[vertTwo] = groupHead[vertOne];
            groupNext[vertTwo] = groupNext[groupHead[vertOne]];
            groupNext[groupHead[vertOne]] = vertTwo;
         }
         else {
            // both vertOne and vertTwo belong to a shared group.
            // if these are the same shared group, then there is nothing
            // else to do. If they are different groups, we must merge
            // the two lists together into a single group.

            if (groupHead[vertOne] != groupHead[vertTwo]) {
               // insert all of vertTwo's list just after the head of
               // vertOne's list

               vertOneHead = groupHead[vertOne];
               secondGuy   = groupNext[vertOneHead];
               currentListTwo = groupHead[vertTwo];
               groupNext[vertOneHead] = currentListTwo;

               while (currentListTwo != -1) {
                  groupHead[currentListTwo] = vertOneHead;
                  prior          = currentListTwo;
                  currentListTwo = groupNext[currentListTwo];
               }
               groupNext[prior] = secondGuy;
            } // if groupHead(vertOne) ~= groupHead(vertTwo)
         }
      }
   }
   return 0;
}
/****************************************************************************/
/* Look at the graph described by (rowPtr,cols,vals) and calculate the
 * distance between any two adjacent vertices. If any of these distances
 * are less than a tol, then record these two vertices in
 * rowZeros and colZeros.
 *
 * Note: adjacent vertices that have small values are ignored. Specifically,
 * dropped entries correspond to
 *           a(i,j) <= a(i,i)*a(j,j)*tol^2
 */
/****************************************************************************/

int ZeroDist(MLVec<double>& xxx, MLVec<double>& yyy, MLVec<double>& zzz,
             MLVec<int>& rowPtr, MLVec<int>& cols, MLVec<double>& vals,
             MLVec<double>& diagonal, double tol, MLVec<int>& rowZeros,  MLVec<int>& colZeros,
             double disttol)
{
   double tol2 = tol*tol, disttol2 = disttol*disttol, normsq, dtemp;
   int    vertTwo;
   int count = 0;

   for (int vertOne = 0; vertOne < rowPtr.size()-1; vertOne++) {
     for (int j = rowPtr[vertOne]; j < rowPtr[vertOne+1]; j++) {
       vertTwo = cols[j];
       if ((vertOne != vertTwo) && (vals[j]*vals[j] > diagonal[vertOne]*diagonal[vertTwo]*tol2)) {
         dtemp = xxx[vertOne] - xxx[vertTwo];
         normsq = dtemp*dtemp;
         if (yyy.empty() == false) { dtemp = yyy[vertOne] - yyy[vertTwo]; normsq += (dtemp*dtemp);}
         if (zzz.empty() == false) { dtemp = zzz[vertOne] - zzz[vertTwo]; normsq += (dtemp*dtemp);}

         if ( normsq <= disttol2) {
            if (count >= rowZeros.size()) {
               int itemp = (int) ceil(count*1.4);
               rowZeros.resize(itemp);
               colZeros.resize(itemp);
            }
            rowZeros[count  ] = vertOne; colZeros[count++] = vertTwo;
         }
       }
     }
   }
   rowZeros.resize(count);
   colZeros.resize(count);

   return 0;
}

/****************************************************************************/
/* Given a list of shared dofs (encoded within groupHead and groupNext),
 * compute an array that maps original rows to a compressed set of rows
 * where all shared dofs correspond to the same row in the compressed set.
 * The list of shared dofs is given by a strange linked list that is described
 * in the comments for buildCompressedA().
 */
/****************************************************************************/
int BuildNonSharedMap(MLVec<int>& newMap, MLVec<int>& groupHead, MLVec<int>& groupNext) {

    for (int i = 0; i < newMap.size(); i++) newMap[i] = 0;

    int Ncols = 0, current;

    for (int i = 0; i < newMap.size(); i++) {
       if  (groupHead[i] == -1) {
          newMap[i] = Ncols;
          Ncols++;
       }
       else {
          if  (groupHead[i] >= 0) {  // first time a new group is encountered
             current = groupHead[i];
             while (current != -1) {
                newMap[current] = Ncols;
                groupHead[current] = -2 - groupHead[current]; // encode as visited
                current = groupNext[current];
             }
             Ncols++;
          }
       }
    }
    // restore groupHead

    for (int i = 0; i < groupHead.size(); i++) {
       if (groupHead[i] < -1) groupHead[i] = -2 - groupHead[i];
    }
    return Ncols;
}

/****************************************************************************/
/*
Given a set of shared dofs, build a compressed graph of the matrix. This
compressed version is equivalent to taking any rows or columns that are
shared and merging them together into a single row or a single column
within a compressed matrix. However, small entries in the original matrix
are dropped when building the compressed graph. Specifically, dropped entries
correspond to
             a(i,j) <= a(i,i)*a(j,j)*tol^2

The list of shared dofs is given by a strange linked list. In particular,
if groupHead[i] == -1, then i is not a shared dof. If groupHead[i] > -1,
then i is shared. All the other dofs that are shared with i are given
by the linked list whose head is groupHead[i]. For example, we might have
the following when four nodes are shared (i,j,k,m):

    ______________      ______________     ______________     ______________
   | groupHead[i] |    | groupHead[j] |   | groupHead[k] |   | groupHead[m] |
   | groupNext[i]-|--->| groupNext[j]-|-->| groupNext[k]-|-->| groupNext[m] |
    --------------      --------------     --------------      --------------
where
    groupHead[i] =    groupHead[j] =    groupHead[k] = groupHead[m]
    groupNext[i] = j; groupNext[j] = k; groupNext[k] = m;
    groupNext[m] = -1;

Additionally, it is assumed that map[] has already been computed such that
map[ii] gives the row/col number in the compressed matrix corresponding to ii.
Thus, in the above example  map[i]=map[j]=map[k]=map[m] as i,j,k,m all map
to the same row. Such a map can be computed by invoking int BuildNonSharedMap()
given groupHead[] and groupNext[].  Finally, newN is also given on input and
it corresponds to the total number of rows in the compressed matrix.
*/
/****************************************************************************/
int buildCompressedA(MLVec<int>& inputRowPtr, MLVec<int>& inputCols,
                     MLVec<double>& inputVals, MLVec<double>& diagonal,
                     MLVec<int>& groupHead, MLVec<int>& groupNext,
                     MLVec<int>& outputRowPtr, MLVec<int>& outputCols,
                     MLVec<int>& map, int newN, double tol)
{
   MLVec<int>   result(newN);
   MLVec<bool>  filled(newN);

   int    count = 0, nFilled, current = 0, newRow = -1, newCol;
   double tol2  = tol*tol;
   int oldRow = 0;


   for (int i=0; i < filled.size() ; i++) filled[i] = false;

   outputRowPtr[0] = 0;
   for (int i=0; i < inputRowPtr.size()-1 ; i++) {
      nFilled = 0;
      current = i;
      if (groupHead[i] > -1) current = groupHead[i];
      if (groupHead[i] < -1) current = -1;

      while (current != -1) {
         newRow = map[current];
         for (int j=inputRowPtr[current]; j < inputRowPtr[current+1]; j++) {
            newCol = map[inputCols[j]];
            // remove connections between shared dofs and nonshared dofs

           if  ( ((groupHead[current]==-1)&&(groupHead[inputCols[j]] == -1)) ||
                 ((groupHead[current]!=-1)&&(groupHead[inputCols[j]] != -1))) {

              if (inputVals[j]*inputVals[j] > diagonal[current]*diagonal[inputCols[j]]*tol2) {
                 if (filled[newCol] == false) {
                    filled[newCol] = true;
                    result[nFilled] = newCol;
                    nFilled++;
                    if (newRow < oldRow) printf("bad\n");
                    for (int k = oldRow; k < newRow; k++) {outputRowPtr[k+1] = count; oldRow = newRow;}
                    outputCols[count] = newCol;
                    count++;
                 } // if filled(newCol) == 0,
              }
           }
         } // for (int j=inputRowPtr[current]; ...

         if (groupHead[current] > -1) groupHead[current]= -2 - groupHead[current]; // visited
         current = groupNext[current];
      } // while (current != -1)
      for (int k=0; k < nFilled; k++)  filled[result[k]] = 0;
      nFilled = 0;
   }
   outputRowPtr[outputRowPtr.size()-1] = count;

   outputCols.resize(count);
   // restore groupHead

   for (int i = 0; i < groupHead.size(); i++) {
      if (groupHead[i] < -1) groupHead[i] = -2 - groupHead[i];
   }
   return 0;
}
#endif
#endif /*ifdef HAVE_ML_EPETRA && HAVE_ML_TEUCHOS */
