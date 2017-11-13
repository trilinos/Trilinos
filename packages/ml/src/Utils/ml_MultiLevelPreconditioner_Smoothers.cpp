/*!
 *  \file ml_MultiLevelPreconditioner_Smoothers.cpp
 *
 *  \brief ML black-box preconditioner for Epetra_RowMatrix derived classes.
 *
 *  \authors Marzio Sala, Ray Tuminaro, Jonathan Hu, Michael Gee, Chris Siefert
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_common.h"
#include "ml_include.h"
#include <algorithm>
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "ml_ifpack_wrap.h"
#include "ml_self_wrap.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "ml_epetra.h"
#include "ml_epetra_utils.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_viz_stats.h"
#include "ml_utils.h"

#ifdef HAVE_ML_IFPACK
#include "Ifpack_Preconditioner.h"
#include "Ifpack_Chebyshev.h"
#endif
#include "ml_petsc.h"

#ifdef HAVE_ML_TekoSmoothers
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace Teko {
class InverseLibrary;
}

extern "C"
// int ML_Gen_Smoother_Teko(ML *ml, int level, int pre_or_post, int ntimes, const std::string & filename, const std::string & inverse,bool isBlocked);
int ML_Gen_Smoother_Teko(ML *ml, int level, int pre_or_post, int ntimes, const Teuchos::RCP<const Teuchos::ParameterList> & tekoPList,
                         const Teuchos::RCP<const Teko::InverseLibrary> & invLib,const std::string & inverse,bool isBlocked);
#endif

extern "C" {
extern int ML_Anasazi_Get_SpectralNorm_Anasazi(ML_Operator * Amat,
                                               ML_Smoother* Smoother,
                                               int MaxIters, double Tolerance,
                                               int IsProblemSymmetric,
                                               int UseDiagonalScaling,
                                               double * LambdaMax );
}

double ML_Smoother_ChebyshevAlpha(double, ML*, int, int);

using namespace Teuchos;


// ===========================================================================

void ML_Operator_Blocked_Getrow(ML_Operator *Amat, int block_size, int requested_block_row, int allocated_space, int columns[],double values[], int row_lengths[]) {
  // Get all of the individual rows
  int offset = 0;
  for (int i = 0; i < block_size; i++) {
    int row = requested_block_row*block_size+i;
    int size=0;
    int status = ML_Operator_Getrow(Amat, 1, &row, allocated_space, &(columns[offset]),&(values[offset]), &size );
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        status == 0,
        "ML_Operator_Blocked_Getrow: Insufficient memory allocated"
        );
    offset+=size;
  }


  // Switch columns to block IDs
  for(int i=0; i<offset; i++) {
    columns[i] = (int) (columns[i] / block_size);
    values[i]  = 1.0;
  }

  // Sort & merge
  std::sort(columns, columns+offset);
  row_lengths[0] = std::unique(columns, columns+offset) - columns;
}


// ============================================================================
inline void local_automatic_line_search(ML * ml, int currentLevel, int NumEqns, int * blockIndices, int last, int next, int LineID, double tol, int *itemp, double * dtemp) {
  double *xvals= NULL, *yvals = NULL, *zvals = NULL;
  int N = ml->Amat[currentLevel].outvec_leng;
  ML_Aggregate_Viz_Stats *grid_info = (ML_Aggregate_Viz_Stats *) ml->Grid[currentLevel].Grid;
  if (grid_info != NULL) xvals = grid_info->x;
  if (grid_info != NULL) yvals = grid_info->y;
  if (grid_info != NULL) zvals = grid_info->z;

  ML_Operator * Amat = &(ml->Amat[currentLevel]);
  int allocated_space = NumEqns*(Amat->max_nz_per_row)+2;
  int * cols    = itemp;
  int * indices = &itemp[allocated_space];
  double * vals = dtemp;
  double * dist = &dtemp[allocated_space];

  // Note: "next" is now a block index, not a point one
  //  printf("[%d] next = %d/%d\n",Amat->comm->ML_mypid,next,N/NumEqns);

  while (blockIndices[next*NumEqns] == -1) {
    // Get the next row
    int n=0;
    int neighbors_in_line=0;

    if(NumEqns==1) ML_Operator_Getrow(Amat,1,&next,allocated_space, cols,vals,&n);
    else           ML_Operator_Blocked_Getrow(Amat,NumEqns,next,allocated_space,cols,vals,&n);

    double x0 = (xvals) ? xvals[next] : 0.0;
    double y0 = (yvals) ? yvals[next] : 0.0;
    double z0 = (zvals) ? zvals[next] : 0.0;

    // Calculate neighbor distances & sort
    int neighbor_len=0;
    //    for(int i=0; i<n; i+=NumEqns) {
    for(int i=0; i<n; i++) {
      double mydist = 0.0;

      if(cols[i]*NumEqns > N) continue; // Check for off-proc entries
      int nn=cols[i];
      if(blockIndices[nn*NumEqns]==LineID) neighbors_in_line++;
      if(xvals!=NULL) mydist += (x0 - xvals[nn]) * (x0 - xvals[nn]);
      if(yvals!=NULL) mydist += (y0 - yvals[nn]) * (y0 - yvals[nn]);
      if(zvals!=NULL) mydist += (z0 - zvals[nn]) * (z0 - zvals[nn]);
      dist[neighbor_len] = sqrt(mydist);
      indices[neighbor_len]=cols[i];
      neighbor_len++;
    }
    // If more than one of my neighbors is already in this line.  I
    // can't be because I'd create a cycle
    if(neighbors_in_line > 1) break;

    // Otherwise add me to the line
    for(int k=0; k<NumEqns; k++)
      blockIndices[next*NumEqns + k] = LineID;

    // Try to find the next guy in the line (only check the closest two that aren't element 0 (diagonal)
    ML_az_dsort2(dist,neighbor_len,indices);

    if(neighbor_len > 2 && indices[1] != last && blockIndices[indices[1]*NumEqns] == -1 && dist[1]/dist[neighbor_len-1] < tol) {
      last=next;
      next=indices[1];
    }
    else if(neighbor_len > 3 && indices[2] != last && blockIndices[indices[2]*NumEqns] == -1 && dist[2]/dist[neighbor_len-1] < tol) {
      last=next;
      next=indices[2];
    }
    else {
      // I have no further neighbors in this line
      break;
    }
  }
}

// ============================================================================
int ML_Compute_Blocks_AutoLine(ML * ml, int currentLevel, int NumEqns, double tol,  int * blockIndices) {
  ML_Operator * Amat = &(ml->Amat[currentLevel]);
  int N = ml->Amat[currentLevel].outvec_leng;
  int allocated_space = NumEqns*(Amat->max_nz_per_row)+2;
  double *xvals= NULL, *yvals = NULL, *zvals = NULL;
  ML_Aggregate_Viz_Stats *grid_info = (ML_Aggregate_Viz_Stats *) ml->Grid[currentLevel].Grid;
  if (grid_info != NULL) xvals = grid_info->x;
  if (grid_info != NULL) yvals = grid_info->y;
  if (grid_info != NULL) zvals = grid_info->z;

  int * cols    = (int    *) ML_allocate(2*allocated_space*sizeof(int   ));
  int * indices = &cols[allocated_space];
  double * vals = (double *) ML_allocate(2*allocated_space*sizeof(double));
  double * dist = &vals[allocated_space];

  int * itemp   = (int    *) ML_allocate(2*allocated_space*sizeof(int   ));
  double *dtemp = (double *) ML_allocate(2*allocated_space*sizeof(double));

  int num_lines = 0;

  // Have everyone check their send lists for block correctness, and error out if one is deficient.
  // We want to make sure that if a processor has a one column in a block for it's column map, it needs
  // all columns in said block.
  ML_CommInfoOP * comm_info = Amat->getrow->pre_comm;
  if(comm_info) {
    bool sends_all_cols=true;
    for(int i=0; i<comm_info->N_neighbors; i++) {
      ML_NeighborList *neighbor = &(comm_info->neighbors[i]);
      int * sends = comm_info->neighbors[i].send_list;
      for (int j=0; sends_all_cols && j<neighbor->N_send; j++) {
        if( ! (j % NumEqns ==0 || sends[j] - sends[j-1] == 1))
          sends_all_cols=false;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
        !sends_all_cols,
        "ML_Compute_Blocks_AutoLine: Incomplete ghost block detected.  Ghost blocks must be complete for line detecion to work"
        );
  }


  // Loop over all of the blocks
  for(int i=0; i<N; i+=NumEqns) {
    int nz=0;
    int ii = i / NumEqns;

    // Short circuit if I've already been blocked
    if(blockIndices[i] !=-1) continue;

    // Get neighbors and sort by distance
    if(NumEqns==1) ML_Operator_Getrow(Amat,1,&ii,allocated_space, cols,vals,&nz);
    else           ML_Operator_Blocked_Getrow(Amat,NumEqns,ii,allocated_space,cols,vals,&nz);

    double x0 = (xvals) ? xvals[ii] : 0.0;
    double y0 = (yvals) ? yvals[ii] : 0.0;
    double z0 = (zvals) ? zvals[ii] : 0.0;
    int neighbor_len=0;
    for(int j=0; j<nz; j++) {
      double mydist = 0.0;
      int nn = cols[j];
      if(cols[j]*NumEqns > N) continue; // Check for off-proc entries
      if(xvals!=NULL) mydist += (x0 - xvals[nn]) * (x0 - xvals[nn]);
      if(yvals!=NULL) mydist += (y0 - yvals[nn]) * (y0 - yvals[nn]);
      if(zvals!=NULL) mydist += (z0 - zvals[nn]) * (z0 - zvals[nn]);
      dist[neighbor_len] = sqrt(mydist);
      indices[neighbor_len]=cols[j];
      neighbor_len++;
    }
    ML_az_dsort2(dist,neighbor_len,indices);

    // Number myself
    for(int k=0; k<NumEqns; k++)
      blockIndices[i + k] = num_lines;

    // Fire off a neighbor line search (nearest neighbor)
    if(neighbor_len > 2 && dist[1]/dist[neighbor_len-1] < tol) {
      local_automatic_line_search(ml,currentLevel,NumEqns,blockIndices,ii,indices[1],num_lines,tol,itemp,dtemp);
    }
    // Fire off a neighbor line search (second nearest neighbor)
    if(neighbor_len > 3 && dist[2]/dist[neighbor_len-1] < tol) {
      local_automatic_line_search(ml,currentLevel,NumEqns,blockIndices,ii,indices[2],num_lines,tol,itemp,dtemp);
    }
    num_lines++;
  }

  // Cleanup
  ML_free(cols);
  ML_free(vals);
  ML_free(itemp);
  ML_free(dtemp);

  return num_lines;
}


// ============================================================================
/*! Values for \c "smoother: type"
 * - \c Jacobi
 * - \c Gauss-Seidel
 * - \c symmetric Gauss-Seidel
 * - \c block Gauss-Seidel
 * - \c symmetric block Gauss-Seidel
 * - \c MLS
 * - \c Chebyshev
 * - \c self
 * - \c Aztec
 * - \c IFPACK
 * - \c Hiptmair
 * - \c ParaSails
 * - \c user-defined
 * - \c do-nothing
 */
int ML_Epetra::MultiLevelPreconditioner::SetSmoothers(bool keepFineLevelSmoother)
{
  Epetra_Time Time(Comm());

  double smooTime=0, coarseTime=0, perLevelTime=0, totalTime=0;

  int num_smoother_steps = List_.get("smoother: sweeps", 2);

  double omega = List_.get("smoother: damping factor",1.0);
  ml_->Cheby_eig_boost = List_.get("smoother: Chebyshev eig boost", 1.1);

  int pre_or_post = 0;
  std::string PreOrPostSmoother = List_.get("smoother: pre or post","both");

  std::string Smoother = List_.get("smoother: type","Chebyshev");

  int  minSizeForCoarse = List_.get("smoother: min size for coarse",-1);
  bool CoarsestIsFine = false;
#ifdef HAVE_ML_AZTECOO
  RCP<std::vector<int> > aztecOptions = List_.get("smoother: Aztec options",SmootherOptions_);
  RCP<std::vector<double> > aztecParams = List_.get("smoother: Aztec params",SmootherParams_);

  bool AztecSmootherAsASolver = List_.get("smoother: Aztec as solver",false);
  int aztec_its;
#endif

  // rst: Changing polynomial interface:
  //    1) polynomial degree is set from "smoother: sweeps"
  //    2) "smoother: Chebyshev" also calls ML's Chebyshev
  //    3) "smoother: Chebshev alpha" also sets alpha
  //    4) "smoother: node sweeps" and "smoother: edge sweeps" now set
  //                  polynomial degree within Hiptmair
  //
  // For backward compatiblity, we still take the degree from
  // "smoother: MLS polynomial order" or "smoother: polynomial order"
  // if set. We also still recognize MLS.
  //
  // Note: At this point, ChebyshevPolyOrder & ChebyshevAlpha have bogus values if
  // they are not set. These get fixed when checking level specific options.

  int ChebyshevPolyOrder = List_.get("smoother: MLS polynomial order",-97);
  if (ChebyshevPolyOrder == -97)
     ChebyshevPolyOrder = List_.get("smoother: polynomial order",-97);

  double ChebyshevAlpha = List_.get("smoother: MLS alpha",-2.0);
  if (ChebyshevAlpha == -2.) ChebyshevAlpha = List_.get("smoother: Chebyshev alpha", -2.0);

  int SmootherLevels = NumLevels_;

  int NumVerticalNodes = List_.get("smoother: line direction nodes",-1);
  std::string MeshNumbering = List_.get("smoother: line orientation","use coordinates");
  std::string GroupDofsString= List_.get("smoother: line group dofs","separate");
  std::string GSType   = List_.get("smoother: line GS Type","symmetric");
  double LineDetectionThreshold = List_.get("smoother: line detection threshold",-1.0);
                                  /* 1: group all dofs per node within a line */
                                  /*    into a single block. Current version  */
                                  /*    is not efficient.                     */
                                  /* 0: don't group dofs per node within a    */
                                  /*    line. Instead if we have n dofs per   */
                                  /*    node, we create n tridiagonal solves  */
                                  /*    for each line.                        */
  int ParaSailsN = List_.get("smoother: ParaSails levels",0);

  // this can be:
  // 0) nonsymmetric and/or indefinite
  // 1) SPD
  // (The ParaSails manual seems to allow also 2, but not ML)
  int ParaSailsSym = List_.get("smoother: ParaSails matrix",0);
  double ParaSailsThresh = List_.get("smoother: ParaSails threshold",0.01);
  double ParaSailsFilter = List_.get("smoother: ParaSails filter",0.05);
  double ParaSailsLB = List_.get("smoother: ParaSails load balancing",0.0);
  int ParaSailsFactorized = List_.get("smoother: ParaSails factorized",0);

  // Ifpack-specific
  std::string IfpackType = List_.get("smoother: ifpack type", "Amesos");
  int IfpackOverlap = List_.get("smoother: ifpack overlap",0);
  // note: lof has different meanings for IC and ICT.  For IC and ILU, we
  // will cast it to an integer later.
  double IfpackLOF=0.0;
  if(IfpackType=="ILUT" || IfpackType=="ICT") IfpackLOF=List_.get("smoother: ifpack level-of-fill",1.);
  else IfpackLOF=List_.get("smoother: ifpack level-of-fill",0.);
  double IfpackRelThreshold=List_.get("smoother: ifpack relative threshold",1.);
  double IfpackAbsThreshold=List_.get("smoother: ifpack absolute threshold",0.);

  // Block Chebyshev parameters
  int cheby_nBlocks=List_.get("smoother: Block Chebyshev number of blocks",-1);
  int *cheby_blockIndices=List_.get("smoother: Block Chebyshev block list",(int*)0);
  int *cheby_blockStarts=List_.get("smoother: Block Chebyshev block starts",(int*)0);

  // Chebyshev-NE parameters
  bool cheby_NE=List_.get("smoother: chebyshev solve normal equations",false);

  // Hiptmair-specific declarations
  std::string SubSmType,NodeSubSmType,EdgeSubSmType;
  int NodeSubSmIts = 1, EdgeSubSmIts = 1;
  double EdgeSubSmLOF=0., NodeSubSmLOF=0.;
  int EdgeSubSmOverlap=0, NodeSubSmOverlap=0;
  double EdgeSubSmOmega=0., NodeSubSmOmega=0.;
  double EdgeSubSmRelThreshold=0., NodeSubSmRelThreshold=0.;
  double EdgeSubSmAbsThreshold=0., NodeSubSmAbsThreshold=0.;

  if (AMGSolver_ == ML_MAXWELL) {
    if (Comm().NumProc() == 1) EdgeSubSmOmega = 1.0;
    else                       EdgeSubSmOmega = ML_DDEFAULT;
    EdgeSubSmOmega = List_.get("subsmoother: damping factor",EdgeSubSmOmega);
    EdgeSubSmOmega = List_.get("subsmoother: edge damping factor",EdgeSubSmOmega);
    if (Comm().NumProc() == 1) NodeSubSmOmega = 1.0;
    else                       NodeSubSmOmega = ML_DDEFAULT;
    NodeSubSmOmega = List_.get("subsmoother: damping factor",NodeSubSmOmega);
    NodeSubSmOmega = List_.get("subsmoother: node damping factor",NodeSubSmOmega);
    SubSmType = List_.get("subsmoother: type","MLS");

    // Grab or set subsmoother options that are not level specific.
    EdgeSubSmType    = List_.get("subsmoother: edge type",SubSmType);
    NodeSubSmType    = List_.get("subsmoother: node type",SubSmType);
    NodeSubSmIts     = List_.get("subsmoother: node sweeps", 2);
    EdgeSubSmIts     = List_.get("subsmoother: edge sweeps", 2);
    EdgeSubSmLOF     = List_.get("subsmoother: edge level-of-fill",0.0);
    NodeSubSmLOF     = List_.get("subsmoother: node level-of-fill",0.0);
    EdgeSubSmOverlap = List_.get("subsmoother: edge overlap",0);
    NodeSubSmOverlap = List_.get("subsmoother: node overlap",0);
    /*
       According to the Ifpack manual, a modified matrix B is factored, where
                 B_ij = A_ij for i \neq j
                 B_ii = \alpha * sgn(A_ii) + \rho * A_ii;
       where
            \alpha = absolute threshold
            \rho   = relative threshold.

       The defaults here are chosen so that B = A.
    */
    EdgeSubSmRelThreshold =List_.get("subsmoother: edge relative threshold",1.);
    NodeSubSmRelThreshold =List_.get("subsmoother: node relative threshold",1.);
    EdgeSubSmAbsThreshold =List_.get("subsmoother: edge absolute threshold",0.);
    NodeSubSmAbsThreshold =List_.get("subsmoother: node absolute threshold",0.);
  }

  totalTime = Time.ElapsedTime();

  // ===================== //
  // cycle over all levels //
  // ===================== //

  char smListName[80];
  int startLevel=0;
  int issueSmootherReuseWarning = keepFineLevelSmoother;
  for (int level = startLevel ; level < SmootherLevels ; ++level) {
    CoarsestIsFine = false;

    if (verbose_) std::cout << std::endl;

    Time.ResetStartTime();

    int coarseLevel  = LevelID_[NumLevels_-1];
    int currentLevel = LevelID_[level];

    long int local[2];
    long int global[2];
    local[0] = ml_->Amat[currentLevel].invec_leng;
    local[1] = ml_->Amat[currentLevel].N_nonzeros;

    if (minSizeForCoarse != -1) {
       // check if the coarsest grid is larger than some user specified value
       // (which might happen if ML aborted its coarsening process). In this
       // case, use the fine smoother instead of the coarse solver
       
       Comm().SumAll(local,global,2);
       if (global[0] > minSizeForCoarse)  CoarsestIsFine = true;
    }

    // general parameters for more than one smoother

    if ( (currentLevel != coarseLevel) || CoarsestIsFine)
      sprintf(smListName,"smoother: list (level %d)",currentLevel);
    else
      strcpy(smListName,"coarse: list");
    ParameterList &smList = List_.sublist(smListName);

    int Mynum_smoother_steps = smList.get("smoother: sweeps",num_smoother_steps);
    double Myomega = smList.get("smoother: damping factor",omega);

    std::string MyPreOrPostSmoother = smList.get("smoother: pre or post", PreOrPostSmoother);

    if( MyPreOrPostSmoother      == "post" ) pre_or_post = ML_POSTSMOOTHER;
    else if( MyPreOrPostSmoother == "pre"  ) pre_or_post = ML_PRESMOOTHER;
    else if( MyPreOrPostSmoother == "both" ) pre_or_post = ML_BOTH;
    else
      std::cerr << ErrorMsg_ << "smoother not recognized (" << MyPreOrPostSmoother << ")\n";

    std::string MySmoother = smList.get("smoother: type",Smoother);

    // If we don't have a level-specific ifpack list, copy the global one
    if(!smList.isSublist("smoother: ifpack list") && List_.isSublist("smoother: ifpack list"))
      ML_Epetra::UpdateList(List_.sublist("smoother: ifpack list"),smList.sublist("smoother: ifpack list"),true);

    char msg[80];
    double AddToDiag = smList.get("smoother: add to diag", 1e-12);
    int MaxProcs=-1;


    if ( (currentLevel != coarseLevel) || CoarsestIsFine) 
    {
      // smoothers
      sprintf(msg,"Smoother (level %d) : ", currentLevel);

      if (minSizeForCoarse == -1) Comm().SumAll(local,global,2);

      // minor information about matrix size on each level
      if (verbose_) {
        int i = std::cout.precision(0);
        std::cout.setf(std::ios::fixed);
        std::cout << msg << "# global rows = " << global[0]
             << ", # estim. global nnz = " << global[1];
        std::cout.precision(2);
        std::cout << ", # nnz per row = " << ((double)global[1]) / ((double) global[0]) << std::endl;
        std::cout.precision(i);
        std::cout.unsetf(std::ios::fixed);
      }
    } else {
      // coarse grid solver
      sprintf(msg,"Coarse solve (level %d) : ", currentLevel);
/*
      MyPreOrPostSmoother = smList.get("smoother: pre or post","post");
      if      (MyPreOrPostSmoother == "post") pre_or_post = ML_POSTSMOOTHER;
      else if (MyPreOrPostSmoother == "pre")  pre_or_post = ML_PRESMOOTHER;
      else if (MyPreOrPostSmoother == "both") pre_or_post = ML_BOTH;
      MySmoother = smList.get("smoother: type", "Amesos-KLU");
      AddToDiag = smList.get("smoother: add to diag", 1e-12);
      Mynum_smoother_steps =  List_.get("coarse: sweeps", 2);
      Myomega = List_.get("coarse: damping factor", 1.0);
      AddToDiag = List_.get("coarse: add to diag", 1e-12);

      ChebyshevAlpha = List_.get("coarse: MLS alpha",-2.0);
      if ( ChebyshevAlpha == -2.)
        ChebyshevAlpha = List_.get("coarse: Chebyshev alpha", 30.);
      ChebyshevPolyOrder = List_.get("coarse: MLS polynomial order",-7);
      if (ChebyshevPolyOrder == -7) ChebyshevPolyOrder = Mynum_smoother_steps;
      MaxProcs = List_.get("coarse: max processes", -1);
*/
/*
//TODO Still need to fix user-defined coarse grid solver
      userSmootherPtr = List_.get("coarse: user-defined function",
                                userSmootherPtr);
      std::string userSmootherName;
      userSmootherName = List_.get("coarse: user-defined name", "User-defined");

//TODO coarse: node sweeps
//TODO coarse: edge sweeps
//TODO coarse: damping factor
*/
    }

    if( MySmoother == "Jacobi" ) {

      // ============ //
      // point Jacobi //
      // ============ //

      if( verbose_ ) std::cout << msg << "Jacobi (sweeps="
                          << Mynum_smoother_steps << ",omega=" << Myomega << ","
                          << MyPreOrPostSmoother << ")" << std::endl;
      ML_Gen_Smoother_Jacobi(ml_, currentLevel, pre_or_post,
                             Mynum_smoother_steps, Myomega);

    } else if( MySmoother == "Gauss-Seidel"){
      // ================== //
      // point Gauss-Seidel //
      // ================== //

      bool gs_type = List_.get("smoother: Gauss-Seidel efficient symmetric",false);
      bool use_l1  = List_.get("smoother: use l1 Gauss-Seidel",false);
      bool use_ml_smoother = ml_->Amat[currentLevel].type != ML_TYPE_CRS_MATRIX && ml_->Amat[currentLevel].type != ML_TYPE_ROW_MATRIX
		&& ml_->Amat[currentLevel].type != ML_TYPE_VBR_MATRIX;
#ifndef HAVE_ML_IFPACK
      use_ml_smoother=true;
#endif

      if( verbose_ ) {
	if(use_ml_smoother)
	  std::cout << msg << "ML Gauss-Seidel (sweeps="
		    << Mynum_smoother_steps << ",omega=" << Myomega << ","
		    << MyPreOrPostSmoother
		    << (gs_type ? ",efficient symmetric" : "" )
		    << ")" <<std::endl;
	else
	  std::cout << msg << "Gauss-Seidel (sweeps="
		    << Mynum_smoother_steps << ",omega=" << Myomega << ","
		    << MyPreOrPostSmoother
		    << (gs_type ? ",efficient symmetric" : "" )
		    << (use_l1  ? ",l1 damping" : "" )
		    << ")" <<std::endl;
      }

      if(!use_ml_smoother) {
#ifdef HAVE_ML_IFPACK
	std::string MyIfpackType = "point relaxation stand-alone";
	ParameterList& MyIfpackList = smList.sublist("smoother: ifpack list");
	MyIfpackList.set("relaxation: type", "Gauss-Seidel");
	MyIfpackList.set("relaxation: sweeps", Mynum_smoother_steps);
	MyIfpackList.set("relaxation: damping factor", Myomega);
	MyIfpackList.set("relaxation: use l1",use_l1);
	int smoothing_indices=0;
	if(MyIfpackList.isParameter("relaxation: number of local smoothing indices"))
	  smoothing_indices = MyIfpackList.get("relaxation: number of local smoothing indices",0);

	if (verbose_) {
	  if (ml_->Amat[currentLevel].type == ML_TYPE_CRS_MATRIX)
	    std::cout << msg << "Epetra_CrsMatrix detected, using "
		      << "Ifpack implementation" << std::endl;
	  else
	    std::cout << msg << "Wrapping to use "
		      << "Ifpack implementation" << std::endl;
	  if (smoothing_indices)
	    std::cout << msg << "Local/reordered smoothing with " << smoothing_indices<<" indices" << std::endl;
	}

	if(gs_type){
	  if(pre_or_post==ML_PRESMOOTHER || pre_or_post==ML_BOTH) {
	    ML_Gen_Smoother_Ifpack(ml_, MyIfpackType.c_str(),
				   IfpackOverlap, currentLevel, ML_PRESMOOTHER,
				   (void*)&MyIfpackList,(void*)Comm_);
	  }
	  if(pre_or_post==ML_POSTSMOOTHER || pre_or_post==ML_BOTH) {
	    ParameterList& BackwardSmoothingList_= MyIfpackList;
	    BackwardSmoothingList_.set("relaxation: backward mode",true);
	    ML_Gen_Smoother_Ifpack(ml_, MyIfpackType.c_str(),
				   IfpackOverlap, currentLevel,  ML_POSTSMOOTHER,
				   (void*)&BackwardSmoothingList_,(void*)Comm_);
	  }
	}
	else{
	  ML_Gen_Smoother_Ifpack(ml_, MyIfpackType.c_str(),
				 IfpackOverlap, currentLevel, pre_or_post,
				 //MyIfpackList,*Comm_);
				 (void*)&MyIfpackList,(void*)Comm_);
	}
#endif
      }
      else {
        if(gs_type)
          ML_Gen_Smoother_EffSymGaussSeidel(ml_, currentLevel, pre_or_post,
                                            Mynum_smoother_steps, Myomega);

        else
          ML_Gen_Smoother_GaussSeidel(ml_, currentLevel, pre_or_post,
                                      Mynum_smoother_steps, Myomega);
      }

    } else if( MySmoother == "ML Gauss-Seidel"){
      // ======================= //
      // ML's point Gauss-Seidel //
      // ======================= //

      bool gs_type = List_.get("smoother: Gauss-Seidel efficient symmetric",false);
      if( verbose_ ) std::cout << msg << "Gauss-Seidel (sweeps="
                         << Mynum_smoother_steps << ",omega=" << Myomega << ","
                         << MyPreOrPostSmoother << (gs_type ? ",efficient symmetric)" : ")") << std::endl;

      if(gs_type)
        ML_Gen_Smoother_EffSymGaussSeidel(ml_, currentLevel, pre_or_post,
                                          Mynum_smoother_steps, Myomega);
      else
        ML_Gen_Smoother_GaussSeidel(ml_, currentLevel, pre_or_post,
                                    Mynum_smoother_steps, Myomega);


    } else if( MySmoother == "symmetric Gauss-Seidel" ) {

      // ====================== //
      // symmetric Gauss-Seidel //
      // ====================== //
      bool use_l1  = List_.get("smoother: use l1 Gauss-Seidel",false);
      bool use_ml_smoother = ml_->Amat[currentLevel].type != ML_TYPE_CRS_MATRIX && ml_->Amat[currentLevel].type != ML_TYPE_ROW_MATRIX
	&& ml_->Amat[currentLevel].type != ML_TYPE_VBR_MATRIX;
#ifndef HAVE_ML_IFPACK
      use_ml_smoother=true;
#endif
      if( verbose_ ) {
	if(use_ml_smoother)
	  std::cout << msg << "ML symmetric Gauss-Seidel (sweeps="
		    << Mynum_smoother_steps << ",omega=" << Myomega << ","
		    << MyPreOrPostSmoother
		    << ")" <<std::endl;
	else
	  std::cout << msg << "symmetric Gauss-Seidel (sweeps="
		    << Mynum_smoother_steps << ",omega=" << Myomega << ","
		    << MyPreOrPostSmoother
		    << (use_l1  ? ",l1 damping" : "" )
		    << ")" <<std::endl;
      }

      if(!use_ml_smoother) {
#ifdef HAVE_ML_IFPACK
	std::string MyIfpackType = "point relaxation stand-alone";
	ParameterList& MyIfpackList = smList.sublist("smoother: ifpack list");;
	MyIfpackList.set("relaxation: type", "symmetric Gauss-Seidel");
	MyIfpackList.set("relaxation: sweeps", Mynum_smoother_steps);
	MyIfpackList.set("relaxation: damping factor", Myomega);
	MyIfpackList.set("relaxation: use l1",use_l1);
	int smoothing_indices=0;
	if(MyIfpackList.isParameter("relaxation: number of local smoothing indices"))
	  smoothing_indices = MyIfpackList.get("relaxation: number of local smoothing indices",0);


	if (verbose_) {
	  if (ml_->Amat[currentLevel].type == ML_TYPE_CRS_MATRIX)
	    std::cout << msg << "Epetra_CrsMatrix detected, using "
		      << "Ifpack implementation" << std::endl;
	  else
	    std::cout << msg << "Wrapping to use "
		      << "Ifpack implementation" << std::endl;
	  if (smoothing_indices)
	    std::cout << msg << "Local/reordered smoothing with " << smoothing_indices<<" indices" << std::endl;
	}

	ML_Gen_Smoother_Ifpack(ml_, MyIfpackType.c_str(),
			       IfpackOverlap, currentLevel, pre_or_post,
			       (void*)&MyIfpackList,(void*)Comm_);
#endif
      }
      else
	ML_Gen_Smoother_SymGaussSeidel(ml_, currentLevel, pre_or_post,
				       Mynum_smoother_steps, Myomega);
    } else if( MySmoother == "ML symmetric Gauss-Seidel" ) {

      // =========================== //
      // ML's symmetric Gauss-Seidel //
      // ============================//

      if( verbose_ ) std::cout << msg << "ML symmetric Gauss-Seidel (sweeps="
                          << Mynum_smoother_steps << ",omega=" << Myomega << ","
                          << MyPreOrPostSmoother << ")" << std::endl;
        ML_Gen_Smoother_SymGaussSeidel(ml_, currentLevel, pre_or_post,
                                       Mynum_smoother_steps, Myomega);

    } else if( MySmoother == "block Gauss-Seidel" ) {

      // ================== //
      // block Gauss-Seidel //
      // ================== //

      if( verbose_ ) std::cout << msg << "block Gauss-Seidel (sweeps="
                          << Mynum_smoother_steps << ",omega=" << Myomega << ","
                          << MyPreOrPostSmoother << ")" << std::endl;
      ML_Gen_Smoother_BlockGaussSeidel(ml_, currentLevel, pre_or_post,
                       Mynum_smoother_steps, Myomega, NumPDEEqns_);

    } else if( MySmoother == "symmetric block Gauss-Seidel" ) {

      // ============================ //
      // symmetric block Gauss-Seidel //
      // ============================ //

      if( verbose_ ) std::cout << msg << "symmetric block Gauss-Seidel (sweeps="
                          << Mynum_smoother_steps << ",omega=" << Myomega << ","
                          << MyPreOrPostSmoother << ")" << std::endl;
      ML_Gen_Smoother_SymBlockGaussSeidel(ml_, currentLevel, pre_or_post,
                                          Mynum_smoother_steps, Myomega, NumPDEEqns_);

    } else if (( MySmoother == "line Gauss-Seidel" )||( MySmoother == "line Jacobi" )){

      if( verbose_ ) std::cout << msg << MySmoother << "(sweeps="
                          << Mynum_smoother_steps << ",omega=" << Myomega << ","
                          << MyPreOrPostSmoother << ")" << std::endl;

       int nnn = ml_->Amat[currentLevel].outvec_leng;
       int MyNumVerticalNodes = smList.get("smoother: line direction nodes",NumVerticalNodes);
       double MyLineDetectionThreshold = List_.get("smoother: line detection threshold",LineDetectionThreshold);
       std::string MyMeshNumbering = smList.get("smoother: line orientation",MeshNumbering);
       std::string MyGSType = smList.get("smoother: line GS Type",GSType);
       std::string MyGroupDofsString= List_.get("smoother: line group dofs",GroupDofsString);
       int MyGroupDofsInLine = 0;
       int NumEqnsOnLevel = ml_->Amat[currentLevel].num_PDEs;


       if (GroupDofsString == "separate") MyGroupDofsInLine = 0;
       if (GroupDofsString == "grouped")  MyGroupDofsInLine = 1;
                                  /* 1: group all dofs per node within a line */
                                  /*    into a single block. Current version  */
                                  /*    is not efficient.                     */
                                  /* 0: don't group dofs per node within a    */
                                  /*    line. Instead if we have n dofs per   */
                                  /*    node, we create n tridiagonal solves  */
                                  /*    for each line.                        */
       ML_GS_SWEEP_TYPE GS_type;
       if      (GSType == "standard") GS_type = ML_GS_standard;
       else if (GSType == "symmetric")  GS_type = ML_GS_symmetric;
       else if (GSType == "efficient symmetric") GS_type = ML_GS_efficient_symmetric;
       else {
         TEUCHOS_TEST_FOR_EXCEPT_MSG(
             true,
             ErrorMsg_ << "smoother: line GS Type not recognized ==>" << MyGSType<< "\n"
             );
       }

       /* the number of nodes per line and even the orientation can change level */
       /* by level, so if these are in the P (via the semicoarsening        */
       /* option), then we should use these values. Basically, this code started  */
       /* out as only a fine level smoother ... so I set up things like          */
       /* smoother:line orientation. However, for multilevel work it seems better*/
       /* to grab these from A, which should be properly set up if one sets the  */
       /* proper semi coarsening options.                                        */
       if (ml_->Amat[currentLevel].NumZDir != -1) {
          MyNumVerticalNodes = ml_->Amat[currentLevel].NumZDir;
          if     (ml_->Amat[currentLevel].Zorientation== 1) MyMeshNumbering= "vertical";
          else if(ml_->Amat[currentLevel].Zorientation== 2) MyMeshNumbering= "horizontal";
          else MyMeshNumbering = "use coordinates";
       }
       if (ml_->Pmat[currentLevel].NumZDir != -1) {
          MyNumVerticalNodes = ml_->Pmat[currentLevel].NumZDir;
          if     (ml_->Pmat[currentLevel].Zorientation== 1) MyMeshNumbering= "vertical";
          else if(ml_->Pmat[currentLevel].Zorientation== 2) MyMeshNumbering= "horizontal";
          else MyMeshNumbering = "use coordinates";
       }
       if (ml_->Pmat[currentLevel].NumZDir == -7) {
          MyNumVerticalNodes = 1;
          MyMeshNumbering = "vertical";
       }

       double *xvals= NULL, *yvals = NULL, *zvals = NULL, *ztemp = NULL;
       ML_Aggregate_Viz_Stats *grid_info = NULL;

       if ((MyMeshNumbering != "horizontal") && (MyMeshNumbering != "vertical")) {

          grid_info = (ML_Aggregate_Viz_Stats *) ml_->Grid[currentLevel].Grid;
          if (grid_info != NULL) xvals = grid_info->x;
          if (grid_info != NULL) yvals = grid_info->y;
          if (grid_info != NULL) zvals = grid_info->z;

          ztemp = zvals;
          if (ml_->Amat[currentLevel].coarsencoord == 'x')
             {zvals = xvals; if (ztemp == NULL) xvals=yvals; else xvals= ztemp;}
          if (ml_->Amat[currentLevel].coarsencoord == 'y')
             {zvals = yvals; if (ztemp == NULL) yvals=xvals; else yvals= ztemp;}
          if (ml_->Amat[currentLevel].coarsencoord == 'z') {
            if ( (zvals == NULL) && (xvals != NULL) && (yvals != NULL) ) {
              printf("Cannot coarsen 2D problems in z direction. Must set semicoarsen_coordinate to x or y\n");
            }
          }

          TEUCHOS_TEST_FOR_EXCEPT_MSG(
              (nnn != 0) && ((xvals == NULL) || (yvals == NULL) || (zvals == NULL)),
              ErrorMsg_ << "line smoother: must supply either coordinates or orientation should be either 'horizontal' or 'vertical'\n"
              );
       }
       else {
          TEUCHOS_TEST_FOR_EXCEPT_MSG(
              MyNumVerticalNodes == -1,
              ErrorMsg_ << "must supply 'line direction nodes' unless line orientation is not supplied and deduced from coordinates" << MySmoother << "\n"
              );
          TEUCHOS_TEST_FOR_EXCEPT_MSG(
              nnn % MyNumVerticalNodes != 0,
              "mod(nnn = " << nnn << ", "
              << "MyNumVerticalNodes = "<< MyNumVerticalNodes
              << ") must be zero\n"
              );
       }
       int *blockOffset  = NULL;
       int *blockIndices = (int *) ML_allocate(sizeof(int)*(nnn+1));

       for (int i = 0; i < nnn;  i++) blockIndices[i] = -1;

       int tempi;
       int NumBlocks;

       if(MyLineDetectionThreshold > 0.0) {
	 // Use Mavriplis-inspired line detection
	 NumBlocks = ML_Compute_Blocks_AutoLine(ml_,currentLevel,NumEqnsOnLevel,MyLineDetectionThreshold,blockIndices);
	 int GlobalBlocks=0;
	 Comm().SumAll(&NumBlocks,&GlobalBlocks,1);
	 if( verbose_ ) std::cout << msg << MySmoother << ": using automatic line detection ("<<GlobalBlocks<<" blocks found)"<<std::endl;
       }
       else {
	 // Use Tuminaro's line detection

	 if (MyMeshNumbering == "vertical") { /* vertical numbering for nodes */
	   if (MyGroupDofsInLine == 1) { /*  one line for all dofs */
	     for (int dof = 0; dof < NumEqnsOnLevel; dof++) {
	       for (int iii = dof; iii < nnn; iii+= NumEqnsOnLevel) {
		 tempi = iii/(NumEqnsOnLevel*MyNumVerticalNodes);
                blockIndices[iii] = tempi;
	       }
	     }
	   }
	   else { /*  different lines for each dof */
	     for (int dof = 0; dof < NumEqnsOnLevel; dof++) {
	       for (int iii = dof; iii < nnn; iii+= NumEqnsOnLevel) {
                tempi = iii/(NumEqnsOnLevel*MyNumVerticalNodes);
                blockIndices[iii] = NumEqnsOnLevel*tempi+dof;
	       }
	     }
	   }
	 }
	 else if (MyMeshNumbering == "horizontal") {/* horizontal numbering for nodes */
          tempi = nnn/MyNumVerticalNodes;
          if (MyGroupDofsInLine == 1) {/*  one line for all dofs */
            for (int iii = 0; iii < nnn; iii++)
	      blockIndices[iii] = (int) floor(((double)(iii%tempi))/
					      ((double) NumEqnsOnLevel)+.00001);
          }
          else { /* different lines for each dof */
            for (int iii = 0; iii < nnn; iii++) blockIndices[iii] = (iii%tempi);
          }
	 }
	 else {

	   blockOffset = (int *) ML_allocate(sizeof(int)*(nnn+1));
	   for (int i = 0; i < nnn;  i++) blockOffset[i] = 0;

	   int    NumCoords, index, next, subindex, subnext;
	   double xfirst, yfirst;

	   NumCoords = nnn/NumEqnsOnLevel;

	   /* sort coordinates so that we can order things according to lines */

	   double *xtemp, *ytemp, *ztemp;
	   int    *OrigLoc;

	   OrigLoc = (int    *) ML_allocate(sizeof(int   )*(NumCoords+1));
	   xtemp   = (double *) ML_allocate(sizeof(double)*(NumCoords+1));
	   ytemp   = (double *) ML_allocate(sizeof(double)*(NumCoords+1));
	   ztemp   = (double *) ML_allocate(sizeof(double)*(NumCoords+1));

     TEUCHOS_TEST_FOR_EXCEPT_MSG(
         ztemp == NULL,
         "Not enough memory for line smoothers\n"
         );
	   for (int i = 0; i < NumCoords; i++) xtemp[i]= xvals[i];
	   for (int i = 0; i < NumCoords; i++) OrigLoc[i]= i;

	   ML_az_dsort2(xtemp,NumCoords,OrigLoc);
	   for (int i = 0; i < NumCoords; i++) ytemp[i]= yvals[OrigLoc[i]];

	   index = 0;

	   while ( index < NumCoords ) {
             xfirst = xtemp[index];
             next   = index+1;
             while ( (next != NumCoords) && (xtemp[next] == xfirst))
	       next++;
             ML_az_dsort2(&(ytemp[index]),next-index,&(OrigLoc[index]));
             for (int i = index; i < next; i++) ztemp[i]= zvals[OrigLoc[i]];
             /* One final sort so that the ztemps are in order */
             subindex = index;
             while (subindex != next) {
	       yfirst = ytemp[subindex]; subnext = subindex+1;
	       while ( (subnext != next) && (ytemp[subnext] == yfirst)) subnext++;
	       ML_az_dsort2(&(ztemp[subindex]),subnext-subindex,&(OrigLoc[subindex]));
	       subindex = subnext;
             }
             index = next;
	   }

	   /* go through each vertical line and populate blockIndices so all   */
	   /* dofs within a PDE within a vertical line correspond to one block.*/

	   NumBlocks = 0;
	   index = 0;
	   int  NotGrouped;

	   NotGrouped = 1 - MyGroupDofsInLine;
	   while ( index < NumCoords ) {
	     xfirst = xtemp[index];  yfirst = ytemp[index];
	     next = index+1;
	     while ( (next != NumCoords) && (xtemp[next] == xfirst) &&
		     (ytemp[next] == yfirst))
               next++;
	     if (NumBlocks == 0) MyNumVerticalNodes = next-index;
       TEUCHOS_TEST_FOR_EXCEPT_MSG(
           next-index != MyNumVerticalNodes,
           "Error code only works for constant block size now!!! "
           << "A size of " << next-index << " found instead of " << MyNumVerticalNodes
           );
	     int count = 0;
	     for (int i = 0; i < NumEqnsOnLevel; i++) {
               if (MyGroupDofsInLine != 1) count = 0;
               for (int j= index; j < next; j++) {
		 blockIndices[NumEqnsOnLevel*OrigLoc[j]+i] = NumBlocks;
		 blockOffset[NumEqnsOnLevel*OrigLoc[j]+i] = count++;
               }
               NumBlocks += NotGrouped;
	     }
	     NumBlocks += MyGroupDofsInLine;
	     index = next;
	   }
	   ML_free(ztemp);
	   ML_free(ytemp);
	   ML_free(xtemp);
	   ML_free(OrigLoc);
	 }
       }// end line detection

      /* check that everyone was assigned to one block */
       for (int i = 0; i < nnn;  i++) {
          int BadCount = 0;
          if (blockIndices[i] == -1) {
             if (BadCount<5) printf("Warning: did not assign %d to a block?????\n",i);
             BadCount++;
          }
       }

       int nBlocks;
       if(MyLineDetectionThreshold > 0.0)
	 nBlocks = NumBlocks;
       else {
	 if (MyGroupDofsInLine == 1) nBlocks = nnn/(MyNumVerticalNodes*NumEqnsOnLevel);
	 else                        nBlocks = nnn/MyNumVerticalNodes;
       }

       if (MySmoother == "line Jacobi") {
           if (MyGroupDofsInLine == 0 && MyLineDetectionThreshold < 0.0 )
             ML_Gen_Smoother_LineSmoother(ml_ , currentLevel, pre_or_post,
                             Mynum_smoother_steps,Myomega,nBlocks,blockIndices,
                             blockOffset, ML_Smoother_LineJacobi, GS_type);
           else
             ML_Gen_Smoother_VBlockJacobi( ml_ , currentLevel, pre_or_post,
                             Mynum_smoother_steps,Myomega,nBlocks,blockIndices);

       } else {
           if (MyGroupDofsInLine == 0 && MyLineDetectionThreshold < 0.0 )
             ML_Gen_Smoother_LineSmoother(ml_ , currentLevel, pre_or_post,
                             Mynum_smoother_steps,Myomega,nBlocks,blockIndices,
                             blockOffset, ML_Smoother_LineGS, GS_type);
           else {
             ML_Gen_Smoother_VBlockSymGaussSeidel(ml_,currentLevel,pre_or_post,
                             Mynum_smoother_steps,Myomega,nBlocks,blockIndices);
             // real hack
             ml_->pre_smoother[currentLevel].gs_sweep_type=GS_type;
             ml_->post_smoother[currentLevel].gs_sweep_type=GS_type;
           }

       }

       ML_free(blockIndices);
       if (blockOffset != NULL) ML_free(blockOffset);
    } else if( ( MySmoother == "MLS" ) || ( MySmoother == "Chebyshev" )
               || (MySmoother == "Block Chebyshev") ) {

      // ========= //
      // Chebyshev //
      // ========= //

      int thisLevel = currentLevel;     // current level
      int nextLevel = 0;                // next coarser level
      if (currentLevel != coarseLevel)
        nextLevel = LevelID_[level+1];
      int MyChebyshevPolyOrder = smList.get("smoother: MLS polynomial order",ChebyshevPolyOrder);
      if (MyChebyshevPolyOrder == -97)
         MyChebyshevPolyOrder = smList.get("smoother: polynomial order",MyChebyshevPolyOrder);
      if (MyChebyshevPolyOrder== -97) MyChebyshevPolyOrder=Mynum_smoother_steps;

      double MyChebyshevAlpha = smList.get("smoother: MLS alpha",ChebyshevAlpha);
      if (MyChebyshevAlpha == -2.)
        MyChebyshevAlpha = smList.get("smoother: Chebyshev alpha",MyChebyshevAlpha);
      if (MyChebyshevAlpha == -2.) MyChebyshevAlpha = 20.;
      MyChebyshevAlpha = ML_Smoother_ChebyshevAlpha(MyChebyshevAlpha, ml_,
                                                    thisLevel, nextLevel);
      /* Grab the Block-Cheby stuff, if applicable */
      int MyCheby_nBlocks=smList.get("smoother: Block Chebyshev number of blocks",cheby_nBlocks);
      int* MyCheby_blockIndices=smList.get("smoother: Block Chebyshev block list",cheby_blockIndices);

      if (verbose_) {
        if (MySmoother == "Block Chebyshev" && MyCheby_blockIndices && MyCheby_nBlocks>0)
        {
          std::cout << msg << "MLS/Block Chebyshev, polynomial order = "
               <<  MyChebyshevPolyOrder
               << ", alpha = " << MyChebyshevAlpha << ", "
               << MyPreOrPostSmoother << std::endl;

        }
        else if (MySmoother == "Chebyshev" && MyChebyshevPolyOrder > 0)
        {
          std::cout << msg << "MLS/Chebyshev, polynomial order = "
               <<  MyChebyshevPolyOrder
               << ", alpha = " << MyChebyshevAlpha << ", "
               << MyPreOrPostSmoother << std::endl;
        }
        else
        {
          std::cout << msg << "MLS, polynomial order = " << -MyChebyshevPolyOrder
               << ", alpha = " << MyChebyshevAlpha << ", "
               << MyPreOrPostSmoother << std::endl;
        }
      }

      /* Check the input block indices */
      if(MySmoother == "Block Chebyshev" && MyCheby_nBlocks>0 && MyCheby_blockIndices) {
        ML_Gen_Smoother_BlockDiagScaledCheby(ml_, currentLevel, pre_or_post,
                                             MyChebyshevAlpha, MyChebyshevPolyOrder,
                                             cheby_nBlocks,cheby_blockIndices);

      }
      else
        ML_Gen_Smoother_Cheby(ml_, currentLevel, pre_or_post,
                              MyChebyshevAlpha, MyChebyshevPolyOrder);

      if (verbose_) {
        ML_Operator* this_A = &(ml_->Amat[currentLevel]);
        std::cout << msg << "lambda_min = " << this_A->lambda_min
             << ", lambda_max = " << this_A->lambda_max << std::endl;
      }
    } else if( MySmoother == "Aztec" ) {

#ifdef HAVE_ML_AZTECOO
      // ======= //
      // AztecOO //
      // ======= //

      // These should remain int* and double*, rather than Teuchos::RCP's.
      // The user created the options & params arrays, so he is responsible for
      // freeing them.
      RCP<std::vector<int> > myaztecOptions   = smList.get("smoother: Aztec options",SmootherOptions_);
      RCP<std::vector<double> > myaztecParams = smList.get("smoother: Aztec params",SmootherParams_);
#ifdef HARDWIRED_AZTEC_SOLVER
RCP<std::vector<int> > m_smootherAztecOptions = rcp(new std::vector<int>(AZ_OPTIONS_SIZE));
RCP<std::vector<double> > m_smootherAztecParams = rcp(new std::vector<double>(AZ_PARAMS_SIZE));
AZ_defaults(&(*m_smootherAztecOptions)[0],&(*m_smootherAztecParams)[0]);
(*m_smootherAztecOptions)[AZ_max_iter]         = 100;
(*m_smootherAztecOptions)[AZ_solver]         = AZ_cg;
(*m_smootherAztecOptions)[AZ_precond]         = AZ_dom_decomp;
(*m_smootherAztecOptions)[AZ_subdomain_solve] = AZ_icc;
myaztecOptions= m_smootherAztecOptions;  // output set in ml_aztec_utils.c
myaztecParams = m_smootherAztecParams;
#endif
      int* MySmootherOptionsPtr = &(*myaztecOptions)[0];
      double* MySmootherParamsPtr = &(*myaztecParams)[0];
      bool MyAztecSmootherAsASolver = smList.get("smoother: Aztec as solver",AztecSmootherAsASolver);

      if( MyAztecSmootherAsASolver == false ) aztec_its = AZ_ONLY_PRECONDITIONER;
      else                                  aztec_its = Mynum_smoother_steps;

      if( verbose_ ) {
        std::cout << msg << "Aztec";
        if( MyAztecSmootherAsASolver){
        switch (MySmootherOptionsPtr[AZ_solver]){
          case AZ_cg:
          case AZ_cg_condnum:
            std::cout<<"-CG";
            break;
          case AZ_gmres:
          case AZ_gmres_condnum:
            std::cout<<"-GMRES";
            break;
          case AZ_cgs:
            std::cout<<"-CGS";
            break;
          case AZ_tfqmr:
            std::cout<<"-TSQMR";
            break;
          case AZ_bicgstab:
            std::cout<<"-TSQMR";
            break;
          case AZ_GMRESR:
            std::cout<<"-GMRESR";
            break;
          }
          std::cout<<"("<<aztec_its<<")";
        }

        if( MySmootherOptionsPtr[AZ_precond] == AZ_dom_decomp ) {
          std::cout << " DD, overlap=" << MySmootherOptionsPtr[AZ_overlap] << ", ";
          if( MySmootherOptionsPtr[AZ_reorder] == 1 ) std::cout << "reord, ";
          else std::cout << "no reord, ";
          switch( MySmootherOptionsPtr[AZ_subdomain_solve] ) {
          case AZ_lu: std::cout << " LU"; break;
          case AZ_ilu:
            std::cout << "ILU(fill="  << MySmootherOptionsPtr[AZ_graph_fill] << ")";
            break;
          case AZ_ilut:
            std::cout << "ILUT(fill=" << MySmootherParamsPtr[AZ_ilut_fill] << ",drop="
             << MySmootherParamsPtr[AZ_drop] << ")";
            break;
          case AZ_icc:
            std::cout << "ICC(fill="  << MySmootherOptionsPtr[AZ_graph_fill] << ")";
            break;
          case AZ_bilu:
            std::cout << "BILU(fill="  << MySmootherOptionsPtr[AZ_graph_fill] << ")";
            break;
          case AZ_rilu:
            std::cout << "RILU(fill="  << MySmootherOptionsPtr[AZ_graph_fill] << ",omega="
             << MySmootherParamsPtr[AZ_omega] << ")";
            break;
          }
        } else if( MySmootherOptionsPtr[AZ_precond] == AZ_Jacobi ) {
          std::cout << " Jacobi preconditioner, sweeps = " << MySmootherOptionsPtr[AZ_poly_ord];
        } else if( MySmootherOptionsPtr[AZ_precond] == AZ_Neumann ) {
          std::cout << " Neumann preconditioner, order = " << MySmootherOptionsPtr[AZ_poly_ord];
        } else if( MySmootherOptionsPtr[AZ_precond] == AZ_ls ) {
          std::cout << " LS preconditioner, order = " << MySmootherOptionsPtr[AZ_poly_ord];
        } else if( MySmootherOptionsPtr[AZ_precond] == AZ_sym_GS ) {
          std::cout << " symmetric Gauss-Seidel preconditioner, sweeps = " << MySmootherOptionsPtr[AZ_poly_ord];
        } else if( MySmootherOptionsPtr[AZ_precond] == AZ_none ) {
          std::cout << " No preconditioning";
        }
        std::cout << ", "  << MyPreOrPostSmoother << std::endl;
      }

      ML_Gen_SmootherAztec(ml_, currentLevel, MySmootherOptionsPtr, MySmootherParamsPtr,
                           ProcConfig_, SmootherStatus_,
                           aztec_its, pre_or_post, NULL);

#else
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
          true,
          "Please configure ML with --enable-aztecoo to use \n"
          << "AztecOO smoothers"
          );
#endif
    } else if( MySmoother == "IFPACK" || MySmoother == "ILU" ||
               MySmoother == "IC" || MySmoother == "ILUT"   ||
               MySmoother == "ICT" || MySmoother == "SILU") {

      // ====== //
      // IFPACK //
      // ====== //

#ifdef HAVE_ML_IFPACK
      std::string MyIfpackType;
      if (MySmoother == "IFPACK")
      {
        MyIfpackType = smList.get("smoother: ifpack type", IfpackType);
      }
      else
      {
        // MS // ILU and IC added on 08-Aug-06 for WebTrilinos
        // MS // Just a shortcut because sublists are not supported by
        // MS // the web interface.
        MyIfpackType = MySmoother;
      }

      double MyLOF=smList.get("smoother: ifpack level-of-fill",IfpackLOF);
      int MyIfpackOverlap = smList.get("smoother: ifpack overlap", IfpackOverlap);
      double MyIfpackRT = smList.get("smoother: ifpack relative threshold", IfpackRelThreshold);
      double MyIfpackAT = smList.get("smoother: ifpack absolute threshold", IfpackAbsThreshold);

      Teuchos::ParameterList& MyIfpackList=smList.sublist("smoother: ifpack list");
      issueSmootherReuseWarning = 2;
      if ( MyIfpackList.isParameter("reuse symbolic factorization") || MyIfpackList.isParameter("reuse numeric factorization") )
        issueSmootherReuseWarning = 0;
      int NumAggr = ML_Aggregate_Get_AggrCount(agg_,level);
      int* AggrMap = 0;
      ML_CHK_ERR(ML_Aggregate_Get_AggrMap(agg_,level,&AggrMap));
      MyIfpackList.set("ILU: sweeps", Mynum_smoother_steps);

      // set these in the case the user wants "partitioner: type" = "user"
      // (if not, these values are ignored).
      if (MyIfpackList.isParameter("partitioner: type")) {
        std::string  partitionerType = MyIfpackList.get<std::string>("partitioner: type");
        if (partitionerType == "user")
          MyIfpackList.set("partitioner: local parts", NumAggr);
        if (partitionerType == "linear" && !MyIfpackList.isParameter("partitioner: local parts"))
          MyIfpackList.set("partitioner: local parts", ml_->Amat[currentLevel].outvec_leng / NumPDEEqns_);
      }
      MyIfpackList.set("partitioner: map", AggrMap);

      // Set the fact: LOF options, but only if they're not set already... All this sorcery is because level-of-fill
      // is an int for ILU and a double for ILUT.  Lovely.
      if(MyIfpackType=="ILUT" || MyIfpackType=="ICT"){
        MyIfpackList.set("fact: level-of-fill", MyIfpackList.get("fact: level-of-fill",MyLOF));
        MyIfpackList.set("fact: ilut level-of-fill", MyIfpackList.get("fact: ilut level-of-fill",MyLOF));
        MyIfpackList.set("fact: ict level-of-fill", MyIfpackList.get("fact: ict level-of-fill",MyLOF));
        MyLOF=MyIfpackList.get("fact: level-of-fill",MyLOF);
      }
      else{
        MyIfpackList.set("fact: level-of-fill", (int) MyIfpackList.get("fact: level-of-fill",(int)MyLOF));
        MyLOF=MyIfpackList.get("fact: level-of-fill",(int)MyLOF);
      }

      MyIfpackList.set("fact: relative threshold", MyIfpackRT);
      MyIfpackList.set("fact: absolute threshold", MyIfpackAT);

      if( verbose_ ) {
        // SORa needs special handling
        if(MyIfpackType == "SORa"){
            std::cout << msg << "IFPACK/SORa("<<MyIfpackList.get("sora: alpha",1.5)<<","<<MyIfpackList.get("sora: gamma",1.0)<<")"
            << ", sweeps = " <<MyIfpackList.get("sora: sweeps",1)<<std::endl;
            if(MyIfpackList.get("sora: oaz boundaries",false))
              std::cout << msg << "oaz boundary handling enabled"<<std::endl;
            if(MyIfpackList.get("sora: use interproc damping",false))
              std::cout << msg << "interproc damping enabled"<<std::endl;
            if(MyIfpackList.get("sora: use global damping",false))
              std::cout << msg << "global damping enabled"<<std::endl;
        }
        else{
          std::cout << msg << "IFPACK, type=`" << MyIfpackType << "'," << std::endl
               << msg << MyPreOrPostSmoother
               << ",overlap=" << MyIfpackOverlap << std::endl;
          if (MyIfpackType != "Amesos") {
            if (MyIfpackType == "ILU" || MyIfpackType == "IC") {
              std::cout << msg << "level-of-fill=" << MyLOF;
            }
            else {
              std::cout << msg << "level-of-fill=" << MyLOF;
            }
            std::cout << ",rel. threshold=" << MyIfpackRT
                 << ",abs. threshold=" << MyIfpackAT << std::endl;
          }
        }
      }
      ML_Gen_Smoother_Ifpack(ml_, MyIfpackType.c_str(),
                             MyIfpackOverlap, currentLevel, pre_or_post,
                             (void*)&MyIfpackList,(void*)Comm_);

#else
      std::cerr << ErrorMsg_ << "IFPACK not available." << std::endl
           << ErrorMsg_ << "ML must be configured with --enable-ifpack" << std::endl
           << ErrorMsg_ << "to use IFPACK as a smoother" << std::endl
           << ErrorMsg_ << "NO SMOOTHER SET FOR THIS LEVEL" << std::endl;
#endif

    } else if( MySmoother == "IFPACK-Chebyshev"  || MySmoother == "IFPACK-Block Chebyshev" ) {

#ifdef HAVE_ML_IFPACK
      int nextLevel = 0;                // next coarser level
      if (currentLevel != coarseLevel)
        nextLevel = LevelID_[level+1];


      int MyChebyshevPolyOrder = smList.get("smoother: MLS polynomial order",ChebyshevPolyOrder);
      if (MyChebyshevPolyOrder == -97)
         MyChebyshevPolyOrder = smList.get("smoother: polynomial order",MyChebyshevPolyOrder);
      if (MyChebyshevPolyOrder== -97) MyChebyshevPolyOrder=Mynum_smoother_steps;

      double MyChebyshevAlpha = smList.get("smoother: MLS alpha",ChebyshevAlpha);
      if (MyChebyshevAlpha == -2.)
        MyChebyshevAlpha = smList.get("smoother: Chebyshev alpha",MyChebyshevAlpha);
      if (MyChebyshevAlpha == -2.) MyChebyshevAlpha = 20.;

      MyChebyshevAlpha = ML_Smoother_ChebyshevAlpha(MyChebyshevAlpha, ml_,
                                       currentLevel,nextLevel);

      /* Grab the Block-Cheby stuff, if applicable */
      int MyCheby_nBlocks=smList.get("smoother: Block Chebyshev number of blocks",cheby_nBlocks);
      int* MyCheby_blockIndices=smList.get("smoother: Block Chebyshev block list",cheby_blockIndices);
      int* MyCheby_blockStarts=smList.get("smoother: Block Chebyshev block starts",cheby_blockStarts);
      bool MyCheby_NE=smList.get("smoother: chebyshev solve normal equations",cheby_NE);

      if( verbose_ ) {
        if (MySmoother == "IFPACK-Block Chebyshev" && MyCheby_blockIndices && MyCheby_blockStarts)
          std::cout << msg << "IFPACK Block Chebyshev, order = " << MyChebyshevPolyOrder
               << ", alpha = " << MyChebyshevAlpha << ", " << MyPreOrPostSmoother << std::endl;
        else
          std::cout << msg << "IFPACK Chebyshev, order = " << MyChebyshevPolyOrder
               << ", alpha = " << MyChebyshevAlpha << ", " << MyPreOrPostSmoother << std::endl;
      }


      ML_Operator* this_A = &(ml_->Amat[currentLevel]);

      Teuchos::ParameterList IFPACKList;
      if(MySmoother == "IFPACK-Block Chebyshev" && MyCheby_blockIndices && MyCheby_blockStarts){
        // If we're using Block Chebyshev, it can compute it's own eigenvalue estimate..
        Teuchos::ParameterList PermuteList,BlockList;
        BlockList.set("apply mode","invert");
        PermuteList.set("number of local blocks",MyCheby_nBlocks);
        PermuteList.set("block start index",MyCheby_blockStarts);
        //        if(is_lid) PermuteList.set("block entry lids",Blockids_);
        //NTS: Add LID support
        PermuteList.set("block entry gids",MyCheby_blockIndices);
        PermuteList.set("blockdiagmatrix: list",BlockList);

        IFPACKList.set("chebyshev: use block mode",true);
        IFPACKList.set("chebyshev: block list",PermuteList);
        IFPACKList.set("chebyshev: eigenvalue max iterations",this_A->spectral_radius_max_iters);

        // EXPERIMENTAL: Cheby-NE
       IFPACKList.set("chebyshev: solve normal equations",MyCheby_NE);
      }
      else {
        // Regular Chebyshev needs an eigenvalue estimate
        ML_Gimmie_Eigenvalues(this_A, ML_DIAGSCALE,
                              this_A->spectral_radius_scheme, ml_->symmetrize_matrix);
      }

      IFPACKList.set("chebyshev: ratio eigenvalue", MyChebyshevAlpha);
      IFPACKList.set("chebyshev: min eigenvalue", this_A->lambda_min);
      IFPACKList.set("chebyshev: max eigenvalue", this_A->lambda_max);
      IFPACKList.set("chebyshev: degree", MyChebyshevPolyOrder);

      ML_Gen_Smoother_Ifpack(ml_, "Chebyshev", 0, currentLevel,
                             pre_or_post, (void*)&IFPACKList, (void*)Comm_);

      if( verbose_ ) {
        std::cout << msg << "lambda_min = " << this_A->lambda_min
             << ", lambda_max = " << this_A->lambda_max << std::endl;
      }

#else
      std::cerr << ErrorMsg_ << "IFPACK not available." << std::endl
           << ErrorMsg_ << "ML must be configured with --enable-ifpack" << std::endl
           << ErrorMsg_ << "to use IFPACK as a smoother" << std::endl
           << ErrorMsg_ << "NO SMOOTHER SET FOR THIS LEVEL" << std::endl;
#endif
    } else if( MySmoother == "self" ) {

#ifdef HAVE_ML_IFPACK
      int MyIfpackOverlap;
      if(smList.isParameter("smoother: self overlap"))
        MyIfpackOverlap = smList.get("smoother: self overlap",0);
      else
        MyIfpackOverlap = List_.get("smoother: self overlap",0);

      if( verbose_ ) {
        std::cout << msg << "ML as self-smoother ("
             << "cycles=" << Mynum_smoother_steps
             << ",overlap=" << MyIfpackOverlap << ","
             << MyPreOrPostSmoother << ")" << std::endl;
      }

      Teuchos::ParameterList MyIfpackList;
      Teuchos::ParameterList& SelfList = MyIfpackList.sublist("ML list");
      Teuchos::ParameterList& tmpList = List_.sublist("smoother: self list");
      SelfList.setParameters(tmpList);
      MyIfpackList.set( "ML node id",List_.get("ML node id",-1) );
      char procLabel[30];
      sprintf(procLabel,"node id %d",List_.get("ML node id",-1));
      SelfList.set("ML label",procLabel);
      SelfList.set("zero starting solution", false);
      std::string xxx = SelfList.get("SetDefaults", "not-set");
      if (xxx != "not-set") {
        if (verbose_ && Comm().MyPID() == 0)
          std::cout << msg << "Setting default values to type `" << xxx << "'" << std::endl;
        SetDefaults(xxx, SelfList,0,0,false);
      }

      if (verbose_ && SelfList.get("ML output",0) > 0)
        std::cout << msg << "*** * Start of self-smoother generation * ***" << std::endl;
      int currentPrintLevel = ML_Get_PrintLevel();
      ML_Gen_Smoother_Self(ml_, MyIfpackOverlap, currentLevel, pre_or_post,
                           Mynum_smoother_steps, MyIfpackList,*Comm_);
      ML_Set_PrintLevel(currentPrintLevel);
      if (verbose_ && SelfList.get("ML output",0) > 0)
        std::cout << msg << "*** * End of self-smoother generation * ***" << std::endl;

#else
      std::cerr << ErrorMsg_ << "IFPACK not available." << std::endl
           << ErrorMsg_ << "ML must be configured with --enable-ifpack" << std::endl
           << ErrorMsg_ << "to use ML as a smoother" << std::endl
           << ErrorMsg_ << "NO SMOOTHER SET FOR THIS LEVEL" << std::endl;
#endif

    } else if( MySmoother == "ParaSails" ) {

      // ========= //
      // ParaSails //
      // ========= //

      int MyParaSailsN = smList.get("smoother: ParaSails levels",ParaSailsN);

      int MyParaSailsSym = smList.get("smoother: ParaSails matrix",ParaSailsSym);

      double MyParaSailsThresh = smList.get("smoother: ParaSails threshold",ParaSailsThresh);

      double MyParaSailsFilter = smList.get("smoother: ParaSails filter",ParaSailsFilter);

      double MyParaSailsLB = smList.get("smoother: ParaSails load balancing",ParaSailsLB);

      int MyParaSailsFactorized = smList.get("smoother: ParaSails factorized",ParaSailsFactorized);

      if( verbose_ )
        std::cout << msg << "ParaSails "
             << "(n=" << MyParaSailsN
             << ",sym=" << MyParaSailsSym
             << ",thresh=" << MyParaSailsThresh
             << ",filter=" << MyParaSailsFilter
             << ",lb=" << MyParaSailsLB
             << "fact=" << MyParaSailsFactorized
             << ")" << std::endl;

#ifdef HAVE_ML_PARASAILS
      // I am not sure about the ending `0' and of ML
      ML_Gen_Smoother_ParaSails(ml_, currentLevel,
                                pre_or_post, Mynum_smoother_steps,
                                MyParaSailsSym, MyParaSailsThresh,
                                MyParaSailsN,
                                MyParaSailsFilter, (int) MyParaSailsLB,
                                MyParaSailsFactorized);
#else
      std::cerr << ErrorMsg_ << "ParaSails not available." << std::endl
           << ErrorMsg_ << "ML must be configured with --with-ml_parasails" << std::endl
           << ErrorMsg_ << "to use ParaSails as a smoother" << std::endl
           << ErrorMsg_ << "NO SMOOTHER SET FOR THIS LEVEL" << std::endl;
#endif

    } else if( MySmoother == "Hiptmair" ) {

      // ==================================================== //
      // Hiptmair                                             //
      // supported subsmoothers:                              //
      //   Chebyshev, SGS, Ifpack incomplete factorizations   //
      // ==================================================== //
      if (AMGSolver_ != ML_MAXWELL) {
        if (Comm().MyPID() == 0) {
          std::cerr << ErrorMsg_ << "Hiptmair smoothing is only supported" << std::endl;
          std::cerr << ErrorMsg_ << "for solving eddy current equations." << std::endl;
          std::cerr << ErrorMsg_ << "Choose another smoother." << std::endl;
        }
        ML_EXIT(EXIT_FAILURE);
      }

      char EdgeSmootherInfo[80], NodeSmootherInfo[80];
      char *SmInfo=0;

      int thisLevel = currentLevel;   // current level
      int nextLevel = 0;                // next coarser level
      if (currentLevel != coarseLevel)
        nextLevel = LevelID_[level+1];
      void *edge_smoother = 0, *nodal_smoother = 0;
      double edge_coarsening_rate=0.0, node_coarsening_rate=0.0;


      // The following section allows a user to specify node & edge options
      // independently.  These variables must exist until the Hiptmair
      // smoother is created for this level, as they (the variables) are
      // are passed as pointers to the smoother create function.
      // Hence, they are declared outside the FOR loop over the levels.

      std::string MyEdgeSubSmType    = smList.get("subsmoother: edge type",EdgeSubSmType);
      std::string MyNodeSubSmType    = smList.get("subsmoother: node type",NodeSubSmType);
      int MyNodeSubSmIts     = smList.get("subsmoother: node sweeps", NodeSubSmIts);
      int MyEdgeSubSmIts     = smList.get("subsmoother: edge sweeps", EdgeSubSmIts);
      double MyEdgeSubSmLOF     = smList.get("subsmoother: edge level-of-fill",EdgeSubSmLOF);
      double MyNodeSubSmLOF     = smList.get("subsmoother: node level-of-fill",NodeSubSmLOF);
      int MyEdgeSubSmOverlap = smList.get("subsmoother: edge overlap",EdgeSubSmOverlap);
      int MyNodeSubSmOverlap = smList.get("subsmoother: node overlap",NodeSubSmOverlap);
      double MyEdgeSubSmRelThreshold=smList.get("subsmoother: edge relative threshold",EdgeSubSmRelThreshold);
      double MyNodeSubSmRelThreshold=smList.get("subsmoother: node relative threshold",NodeSubSmRelThreshold);
      double MyEdgeSubSmAbsThreshold=smList.get("subsmoother: edge absolute threshold",EdgeSubSmAbsThreshold);
      double MyNodeSubSmAbsThreshold=smList.get("subsmoother: node absolute threshold",NodeSubSmAbsThreshold);
      double MyEdgeSubSmOmega   = smList.get("subsmoother: edge omega",EdgeSubSmOmega);
      double MyNodeSubSmOmega   = smList.get("subsmoother: node omega",NodeSubSmOmega);
      double MyEdgeSubSmAlpha   = smList.get("subsmoother: edge alpha",-2.0);
      double MyNodeSubSmAlpha   = smList.get("subsmoother: node alpha",-2.0);

      Teuchos::ParameterList& edgeList = smList.sublist("edge list");
      Teuchos::ParameterList& nodeList = smList.sublist("node list");

      // ++++++++++++++++++++++++++++++++++++++++++++++
      // Set the node and edge subsmoothers separately.
      // ----------------------------------------------
      enum nodeOrEdge {NODE, EDGE, DONE};
      for (enum nodeOrEdge ne=NODE; ne!= DONE; ne=nodeOrEdge(ne+1)) {

        std::string *MySubSmType=0;
        int *MySubSmIts=0, *MySubSmOverlap=0;
        double MySubSmLOF=0.,MySubSmRelThreshold=0.,MySubSmAbsThreshold=0.,
               *MySubSmOmega=0, MySubSmAlpha=0.;
        Teuchos::ParameterList *ifpackList=0;
        void **argList=0;

        switch(ne) {
          case NODE:
            ifpackList = &nodeList;
            MySubSmType = &MyNodeSubSmType;
            MySubSmIts = &MyNodeSubSmIts;
            MySubSmLOF = MyNodeSubSmLOF;
            MySubSmOverlap = &MyNodeSubSmOverlap;
            MySubSmRelThreshold = MyNodeSubSmRelThreshold;
            MySubSmAbsThreshold = MyNodeSubSmAbsThreshold;
            MySubSmOmega = &MyNodeSubSmOmega;
            MySubSmAlpha = MyNodeSubSmAlpha;
            SmInfo = NodeSmootherInfo;
            break;
          case EDGE:
            ifpackList = &edgeList;
            MySubSmType = &MyEdgeSubSmType;
            MySubSmIts = &MyEdgeSubSmIts;
            MySubSmLOF = MyEdgeSubSmLOF;
            MySubSmOverlap = &MyEdgeSubSmOverlap;
            MySubSmRelThreshold = MyEdgeSubSmRelThreshold;
            MySubSmAbsThreshold = MyEdgeSubSmAbsThreshold;
            MySubSmOmega = &MyEdgeSubSmOmega;
            MySubSmAlpha = MyEdgeSubSmAlpha;
            SmInfo = EdgeSmootherInfo;
            break;
          case DONE:
            pr_error("Something has gone wrong in Hiptmair smoother setup\n");
            break;
        } //switch(ne)

#if defined(__GNUC__) && defined(__GNUC_MINOR__) && defined(__GNUC_PATCHLEVEL__)
#define GCC_VERSION __GNUC__*100+__GNUC_MINOR__*10+__GNUC_PATCHLEVEL__
#endif

#if defined(GCC_VERSION) && GCC_VERSION >= 460
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#endif
        if ( (*MySubSmType == "MLS") || (*MySubSmType == "Chebyshev"))
        {
          // --------------------------------------
          // Chebyshev subsmoother
          // --------------------------------------
          double *coarsening_rate=0;
          ML  *mlptr=0;
          if (ne == EDGE) {
            edge_smoother=(void *) ML_Gen_Smoother_Cheby;
            edge_args_ = ML_Smoother_Arglist_Create(2);
            argList = edge_args_;
            coarsening_rate = &edge_coarsening_rate;
            mlptr = ml_;
          } else if (ne == NODE) {
            nodal_smoother=(void *) ML_Gen_Smoother_Cheby;
            nodal_args_ = ML_Smoother_Arglist_Create(2);
            argList = nodal_args_;
            coarsening_rate = &node_coarsening_rate;
            mlptr = ml_nodes_;
          }
          // This is for backward compatibility
          int itemp = List_.get("subsmoother: MLS polynomial order",-97);
          if (itemp == -97) itemp=List_.get("subsmoother: polynomial order",-97);
          itemp = smList.get("subsmoother: MLS polynomial order",itemp);
          if (itemp != -97) *MySubSmIts = itemp;

          double SubAlpha = List_.get("subsmoother: MLS alpha",MySubSmAlpha);
          if (SubAlpha == -2.) SubAlpha=List_.get("subsmoother: Chebyshev alpha", -2.);
          if (SubAlpha == -2.) SubAlpha = 20.;
          *coarsening_rate = ML_Smoother_ChebyshevAlpha(SubAlpha, mlptr,
                                                        thisLevel, nextLevel);
          ML_Smoother_Arglist_Set(argList, 0, MySubSmIts);
          ML_Smoother_Arglist_Set(argList, 1, coarsening_rate);

          // FIXME:  T could be NULL
          sprintf(SmInfo,"Chebyshev,degree=%d,alpha=%5.3e",
                  *MySubSmIts,*coarsening_rate);

        } else if (*MySubSmType == "symmetric Gauss-Seidel") {

          // --------------------------------------
          // symmetric Gauss-Seidel subsmoother
          // --------------------------------------
          if (ne == EDGE) {
            edge_smoother=(void *) ML_Gen_Smoother_SymGaussSeidel;
            edge_args_ = ML_Smoother_Arglist_Create(2);
            argList = edge_args_;
          } else if (ne == NODE) {
            nodal_smoother=(void *) ML_Gen_Smoother_SymGaussSeidel;
            nodal_args_ = ML_Smoother_Arglist_Create(2);
            argList = nodal_args_;
          }
          ML_Smoother_Arglist_Set(argList, 0, MySubSmIts);
          ML_Smoother_Arglist_Set(argList, 1, MySubSmOmega);
          sprintf(SmInfo,"symmetric Gauss-Seidel,sweeps=%d,omega=%5.3f",
                  *MySubSmIts,*MySubSmOmega);

        } else if ( *MySubSmType == "ILU" ||
                  *MySubSmType == "IC"  ||
                  *MySubSmType == "ILUT" ||
                  *MySubSmType == "ICT")       {
#ifdef HAVE_ML_IFPACK
          // --------------------------------------
          // incomplete factorization subsmoothers
          // --------------------------------------
          if (ne == EDGE) {
              edge_args_ = ML_Smoother_Arglist_Create(4);
              edge_smoother=(void *) ML_Gen_Smoother_Ifpack;
              argList = edge_args_;
          } else if (ne == NODE) {
            nodal_args_ = ML_Smoother_Arglist_Create(4);
            nodal_smoother=(void *) ML_Gen_Smoother_Ifpack;
            argList = nodal_args_;
          }
          //ML uses the same parameter for all of Ifpack's levels-of-fill.
          //Rather than figure out which Ifpack method is really being
          //used, I just set them all.
          ifpackList->set("fact: level-of-fill", (int)(MySubSmLOF));
          ifpackList->set("fact: ict level-of-fill", MySubSmLOF);
          ifpackList->set("fact: ilut level-of-fill", MySubSmLOF);
          ifpackList->set("fact: relative threshold", MySubSmRelThreshold);
          ifpackList->set("fact: absolute threshold", MySubSmAbsThreshold);
          ML_Smoother_Arglist_Set(argList, 0, const_cast<char*>(MySubSmType->c_str()));
          ML_Smoother_Arglist_Set(argList, 1, ifpackList);
          ML_Smoother_Arglist_Set(argList, 2, MySubSmOverlap);
          ML_Smoother_Arglist_Set(argList, 3,const_cast<Epetra_Comm*>(Comm_));

          if ( *MySubSmType == "ILU" || *MySubSmType == "IC")
            sprintf(SmInfo,"%s,overlap=%d,level-of-fill=%d",
                    MySubSmType->c_str(),*MySubSmOverlap,(int)MySubSmLOF);
          else
            sprintf(SmInfo,"%s,overlap=%d,level-of-fill=%3.2e",
                    MySubSmType->c_str(),*MySubSmOverlap,MySubSmLOF);
#else
          pr_error("%sIFPACK subsmoother unavailable.  Configure with --enable-ifpack.\n",
                    ErrorMsg_);
#endif
        } else if (Comm().MyPID() == 0)
          std::cerr << ErrorMsg_
            <<"Only Chebyshev (or MLS), SGS, ILU, IC, ILUT, and ICT" << std::endl
            << "are supported as Hiptmair subsmoothers ... not "
            << *MySubSmType << std::endl;
#if defined(GCC_VERSION) && GCC_VERSION >= 460
#pragma GCC diagnostic pop
#endif

      } //for (enum nodeOrEdge ne=NODE; ne!=DONE ...


      int hiptmair_type = (int)
                  List_.get("smoother: Hiptmair efficient symmetric", true);

      if( verbose_ ) std::cout << msg << "Hiptmair (outer sweeps="
             << Mynum_smoother_steps << ")" << std::endl
             << msg << "edge: " << EdgeSmootherInfo << std::endl
             << msg << "node: " << NodeSmootherInfo << std::endl;

      ML_Gen_Smoother_Hiptmair2(ml_, thisLevel, ML_BOTH,
                                Mynum_smoother_steps, Tmat_array, Tmat_trans_array, NULL,
                                MassMatrix_array,TtATMatrixML_,
                                edge_smoother, edge_args_, nodal_smoother, nodal_args_,
                                hiptmair_type);

      ML_Smoother_Arglist_Delete(&nodal_args_);
      ML_Smoother_Arglist_Delete(&edge_args_);

      bool indefiniteProblem = List_.get("negative conductivity",false);
      /* This corrects for the case of negative sigma        */
      /* I think it has something to do with the eigenvalues */
      /* coming out negative.                                */
      if (indefiniteProblem && MyNodeSubSmType == "MLS") //JJH check this
      {
        if (verbose_ && Comm().MyPID() == 0)
          std::cout << "ML*WRN* "
             << "Resetting nodal smoother on level " << thisLevel << std::endl
             << "ML*WRN* to account for negative mass matrix." << std::endl;
        //pre-smoother

        ML_Sm_Hiptmair_Data *hiptmairSmData =
           (ML_Sm_Hiptmair_Data *) ml_->pre_smoother[thisLevel].smoother->data;
        ML *ml_subproblem = hiptmairSmData->ml_nodal;

        struct MLSthing *widget =
                       (struct MLSthing *) ml_subproblem->pre_smoother->smoother->data;
        double eig_ratio = widget->eig_ratio;
        int degree = widget->mlsDeg;
        ml_subproblem->pre_smoother->data_destroy(
                                 ml_subproblem->pre_smoother->smoother->data);
        ML_Operator *Amat = ml_subproblem->Amat;
        if ( Amat->lambda_max < 0.) {
           double tmp = Amat->lambda_max;
           Amat->lambda_max = fabs(Amat->lambda_min);
           Amat->lambda_min = fabs(tmp);
           ML_Gen_Smoother_Cheby(ml_subproblem,0,ML_PRESMOOTHER,eig_ratio,degree);

           //post-smoother
           hiptmairSmData = (ML_Sm_Hiptmair_Data *)
                          ml_->post_smoother[thisLevel].smoother->data;
           ml_subproblem = hiptmairSmData->ml_nodal;

           // Note:  this is correct because the pre_smoother is the only one
           // used in the subproblem
           widget = (struct MLSthing *) ml_subproblem->pre_smoother->smoother->data;
           eig_ratio = widget->eig_ratio;
           degree = widget->mlsDeg;
           ml_subproblem->pre_smoother->data_destroy(
               ml_subproblem->pre_smoother->smoother->data);

           ML_Gen_Smoother_Cheby(ml_subproblem,0,ML_PRESMOOTHER,eig_ratio,degree);
        }
      }

    } else if( MySmoother == "petsc" ) {

      // ======================================== //
      // PETSc smoother (for PETSc matrices only) //
      // ======================================== //

      // We assume the PETSc application has set up the KSP object entirely.
      // ML just applies it.

#     ifdef HAVE_ML_PETSC


/*
      void *voidPC = 0;
      ML_PetscPC petscPC = 0;
      petscPC = (ML_PetscPC) smList.get("smoother: petsc pc", voidPC);
*/
      void *voidKSP = 0;
      ML_PetscKSP petscKSP = 0;
      petscKSP = (ML_PetscKSP) smList.get("smoother: petsc ksp", voidKSP);
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
          petscKSP == 0,
          ErrorMsg_
          << "You must provide a fully-constructed KSP context to use a PETSc smoother."
          );
      const char* pcName;
      ML_PetscPC petscPC;
      int ierr = KSPGetPC(petscKSP,&petscPC);CHKERRQ(ierr);
      ierr = PCGetType(petscPC,&pcName);
      if( verbose_ ) std::cout << msg << "PETSc smoother (type="
                          << pcName
                          << ",sweeps=" << Mynum_smoother_steps << ","
                          << MyPreOrPostSmoother << ")" << std::endl;

      ML_Gen_Smoother_Petsc(ml_, currentLevel, pre_or_post, Mynum_smoother_steps, petscKSP);

#     else
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
          true,
          ErrorMsg_ << "You must configure ML with PETSc support enabled."
          );
#     endif /*ifdef HAVE_ML_PETSC*/
    } else if( MySmoother == "teko" ) {
#ifdef HAVE_ML_TekoSmoothers
      // ======================================== //
      // Teko smoother (for block matrices only) //
      // ======================================== //
      Teuchos::RCP<const Teko::InverseLibrary> invLib =
            List_.get<Teuchos::RCP<const Teko::InverseLibrary> >("smoother: teko inverse library",Teuchos::null);

      std::string tekoFilename = List_.get<std::string>("smoother: teko filename","teko_smoother.xml");
      Teuchos::RCP<Teuchos::ParameterList> tekoPList
             = List_.get<Teuchos::RCP<Teuchos::ParameterList> >("smoother: teko parameter list",Teuchos::null);
      std::string tekoInverse = List_.get<std::string>("smoother: teko inverse");
      int isBlocked = List_.get<int>("smoother: teko is blocked");

      tekoFilename = smList.get("smoother: teko filename",tekoFilename);
      tekoPList = smList.get("smoother: teko parameter list",tekoPList);
      tekoInverse = smList.get("smoother: teko inverse",tekoInverse);
      isBlocked = smList.get("smoother: teko is blocked",isBlocked);

      // if no parameter list read one from the specified file
      if(tekoPList==Teuchos::null && invLib==Teuchos::null)
        tekoPList = Teuchos::getParametersFromXmlFile(tekoFilename);

      // ML_Gen_Smoother_Teko(ml_, currentLevel, pre_or_post, Mynum_smoother_steps,
      //                      tekoFilename,tekoInverse,isBlocked);
      ML_Gen_Smoother_Teko(ml_, currentLevel, pre_or_post, Mynum_smoother_steps,
                           tekoPList,invLib,tekoInverse,isBlocked);
#else
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
          true,
          ErrorMsg_
          << "You must configure ML with Teko support enabled. "
          << "Enable flag ENABLE_TekoML and add Teko to the end of the library line"
          );
#endif
    } else if( MySmoother == "user-defined" || MySmoother == "user defined" ) {

      // ============ //
      // user-defined //
      // ============ //

      int (*userSmootherPtr)(ML_Smoother *, int, double *, int, double *);
      userSmootherPtr = NULL;
      userSmootherPtr = List_.get("smoother: user-defined function",
                                  userSmootherPtr);
      std::string userSmootherName;
      userSmootherName = List_.get("smoother: user-defined name",
                                   "User-defined");

      if( verbose_ ) std::cout << msg << userSmootherName << " (sweeps="
                          << Mynum_smoother_steps << ","
                          << MyPreOrPostSmoother << ")" << std::endl;

      if (userSmootherPtr == NULL) {
        if (Comm().MyPID() == 0)
          std::cerr << ErrorMsg_
               << "No pointer to user-defined smoother function found." << std::endl;
        ML_EXIT(EXIT_FAILURE);
      }
      ML_Operator *data;
      ML_Get_Amatrix(ml_, currentLevel, &data);
      ML_Set_Smoother(ml_, currentLevel, pre_or_post , data,
                      userSmootherPtr,
                      const_cast<char *>(userSmootherName.c_str()));

    } else if( MySmoother == "SuperLU" ) {
      ML_Gen_CoarseSolverSuperLU( ml_, LevelID_[NumLevels_-1]);

    } else if( MySmoother == "Amesos-LAPACK" ) {
      ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1],
                             ML_AMESOS_LAPACK, MaxProcs, AddToDiag);
    } else if( MySmoother == "Amesos-KLU" ) {
      ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1],
                             ML_AMESOS_KLU, MaxProcs, AddToDiag);
    } else if( MySmoother == "Amesos-UMFPACK" ) {
      ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1],
                             ML_AMESOS_UMFPACK, MaxProcs, AddToDiag);
    } else if(  MySmoother == "Amesos-Superludist" ) {
      ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1],
                             ML_AMESOS_SUPERLUDIST, MaxProcs, AddToDiag);
    } else if(  MySmoother == "Amesos-Superlu" ) {
      ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1],
                             ML_AMESOS_SUPERLU, MaxProcs, AddToDiag);
    } else if( MySmoother == "Amesos-MUMPS" ) {
      ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1],
                             ML_AMESOS_MUMPS, MaxProcs, AddToDiag);
    } else if( MySmoother == "Amesos-ScaLAPACK" ) {
      ML_Gen_Smoother_Amesos(ml_, LevelID_[NumLevels_-1],
                             ML_AMESOS_SCALAPACK, MaxProcs, AddToDiag);

    } else if( MySmoother == "do-nothing" ) {

      // ======================== //
      // do-nothing (no smoother) //
      // ======================== //

      if( verbose_ ) std::cout << msg << "do-nothing smoother" << std::endl;

    } else {

      // ======================================= //
      // error: no smoother found with this name //
      // ======================================= //

      if (Comm().MyPID() == 0)
         std::cerr << ErrorMsg_
              << "Smoother '" << MySmoother << "' not recognized!" << std::endl
              << ErrorMsg_
              << "(file " << __FILE__ << ",line " << __LINE__ << ")" << std::endl
              << ErrorMsg_
              << "You chose: " << MySmoother << ". It should be: " << std::endl
              << ErrorMsg_
              << "<Jacobi> / <Gauss-Seidel> / <block Gauss-Seidel>" << std::endl
              << ErrorMsg_
              << "<symmetric Gauss-Seidel> / <Aztec> / <IFPACK>" << std::endl
              << ErrorMsg_
              << "<Chebyshev> / <ParaSails> / <Hiptmair>" << std::endl
              << ErrorMsg_ << "<user-defined>" << std::endl;
      ML_EXIT(-99); }

    perLevelTime = Time.ElapsedTime();
    if (currentLevel != coarseLevel) {
      smooTime += perLevelTime;
      if (verbose_)
        std::cout << msg << "Setup time : " << perLevelTime << " (s)" << std::endl;
    }
    else
      coarseTime = perLevelTime;

    if (level > 0) issueSmootherReuseWarning = 0;

    if (verbose_) {
      if (issueSmootherReuseWarning == 1) {
        std::cout << std::endl;
        std::cout << "WARNING **********************************************************************" << std::endl
                  << "You have requested to keep the finest level smoother from a previous" << std::endl
                  << "solve.  This option is fragile and intended for reusing part or all of" << std::endl
                  << "an incomplete factorization smoother only. You are attempting to reuse the" << std::endl
                  << "\"" << MySmoother << "\" smoother, which may result in memory leaks or crashes." << std::endl;
      }
      if (issueSmootherReuseWarning == 2) {
        std::cout << std::endl;
        std::cout << "WARNING **********************************************************************" << std::endl
                  << "You have requested to keep the finest level smoother from a previous" << std::endl
                  << "solve.  In order to reuse part or all of the incomplete factorization" << std::endl
                  << "from smoother \"" << MySmoother << "\", you must also specify the smoother options" << std::endl
                  << "\"reuse numeric factorization\" and/or \"reuse symbolic factorization\"." << std::endl
                  << "Not doing so  may result in memory leaks or crashes." << std::endl;
      }
    } //if (verbose_)

  } /* for (int level = 0 ; level < SmootherLevels ; ++level) */

  totalTime += (smooTime + coarseTime);
  OutputList_.set("time: smoothers setup", smooTime
                  + OutputList_.get("time: smoothers setup", 0.0));
  OutputList_.set("time: coarse solver setup", coarseTime
                  + OutputList_.get("time: coarse solver setup", 0.0));
  if(  verbose_ ) std::cout << std::endl;

  return(0);
}

/*------------------------------------------------------------------------
function ML_Smoother_ChebyshevAlpha()
   alpha  -- default alpha for Chebyshev polynomial
   ml     -- pointer to the ML structure
   here   -- level number of this level
   next   -- level number of next coarser level
------------------------------------------------------------------------*/
double ML_Smoother_ChebyshevAlpha(double alpha, ML* ml,int here, int next)
{
  int itmp, Ncoarse,
      Nfine = ml->Amat[here].outvec_leng,
      coarsest_level = ml->ML_coarsest_level;
  double coarsening_rate=0.0;

  ML_gsum_scalar_int(&Nfine, &itmp, ml->comm);
  if (here != coarsest_level) {
    Ncoarse = ml->Amat[next].outvec_leng;
    ML_gsum_scalar_int(&Ncoarse, &itmp, ml->comm);
    if (Ncoarse != 0.0)
      coarsening_rate =  ((double) Nfine) / ((double) Ncoarse);
      //coarsening_rate =  2.0*((double) Nfine) / ((double) Ncoarse);
  }
  //printf("level %d, before, coarsening_rate = %e, nc = %d, nf = %d\n",
  //        here, coarsening_rate, Ncoarse, Nfine);
  if (coarsening_rate < alpha)
    coarsening_rate =  alpha;
  return coarsening_rate;
} //ML_Smoother_ChebyshevAlpha()

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/

