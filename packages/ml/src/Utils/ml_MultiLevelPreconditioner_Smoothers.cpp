/*!
 *  \file ml_MultiLevelPreconditioner_Smoothers.cpp
 *
 *  \brief ML black-box preconditioner for Epetra_RowMatrix derived classes.
 *
 *  \author Marzio Sala, SNL, 9214
 *
 *  \date Last update do Doxygen: 22-Jul-04
 *
 */

#include "ml_common.h"
#include "ml_include.h"
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
#include "ml_epetra.h"
#include "ml_MultiLevelPreconditioner.h"
#ifdef HAVE_ML_IFPACK
#include "Ifpack_Preconditioner.h"
#include "Ifpack_Chebyshev.h"
#endif
extern "C" {
extern int ML_Anasazi_Get_SpectralNorm_Anasazi(ML_Operator * Amat,
                                               ML_Smoother* Smoother,
                                               int MaxIters, double Tolerance,
                                               int IsProblemSymmetric,
                                               int UseDiagonalScaling,
                                               double * LambdaMax );
}

using namespace Teuchos;

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
int ML_Epetra::MultiLevelPreconditioner::SetSmoothers() 
{

  char parameter[80];
  Epetra_Time Time(Comm());

  int num_smoother_steps = List_.get("smoother: sweeps", 1);

  double omega = List_.get("smoother: damping factor",1.0);

  int pre_or_post = 0;
  string PreOrPostSmoother = List_.get("smoother: pre or post","both");

#ifndef HAVE_ML_AZTECOO
  string Smoother = List_.get("smoother: type","symmetric Gauss-Seidel");
#else
  string Smoother = List_.get("smoother: type","Aztec");

  int* SmootherOptionsPtr = (int*)0;
  SmootherOptionsPtr = List_.get("smoother: Aztec options",SmootherOptionsPtr);

  double* SmootherParamsPtr = (double*)0;
  SmootherParamsPtr = List_.get("smoother: Aztec params",SmootherParamsPtr);

  bool AztecSmootherAsASolver = List_.get("smoother: Aztec as solver",false);
  int aztec_its;
#endif

  /* It does not make sense to use both a polynomial order and */
  /* smoother sweeps. If the user has set the sweeps and not   */
  /* the polynomial order, let us just assume that sweeps was  */
  /* really intended as polynomial order.                      */

  int MLSPolynomialOrder = List_.get("smoother: MLS polynomial order",-7);
  if (MLSPolynomialOrder == -7) {
     if (num_smoother_steps != 1) MLSPolynomialOrder = num_smoother_steps;
     else MLSPolynomialOrder = 3;
  }
  double MLSalpha = List_.get("smoother: MLS alpha",30.0);
  
  // IFPACK version of Chebyshev (no MLS here)
  int PolynomialOrder = List_.get("smoother: polynomial order",-7);
  if (PolynomialOrder == -7) {
     if (num_smoother_steps != 1) PolynomialOrder = num_smoother_steps;
     else PolynomialOrder = 3;
  }
  double alpha = List_.get("smoother: alpha", 30.0);

  int SmootherLevels = (NumLevels_>1)?(NumLevels_-1):1;

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

  string IfpackType = List_.get("smoother: ifpack type", "Amesos");
  int IfpackOverlap = List_.get("smoother: ifpack overlap",0);

  string SubSmootherType;
  int nodal_its, edge_its;
  if (SolvingMaxwell_ == true) {
    SubSmootherType = List_.get("subsmoother: type","MLS");
    nodal_its = List_.get("subsmoother: node sweeps", 1);
    edge_its = List_.get("subsmoother: edge sweeps", 1);
  }

  // ===================== //
  // cycle over all levels //
  // ===================== //
 
  for (int level = 0 ; level < SmootherLevels ; ++level) {

    if (verbose_) cout << endl;

    Time.ResetStartTime();

    // general parameters for more than one smoother

    sprintf(parameter,"smoother: sweeps (level %d)", LevelID_[level] );
    num_smoother_steps = List_.get(parameter,num_smoother_steps);

    sprintf(parameter,"smoother: damping factor (level %d)", LevelID_[level] );
    omega = List_.get(parameter,omega);

    sprintf(parameter,"smoother: pre or post (level %d)", LevelID_[level] );
    PreOrPostSmoother = List_.get(parameter, PreOrPostSmoother);
    
    if( PreOrPostSmoother == "post" ) pre_or_post = ML_POSTSMOOTHER;
    else if( PreOrPostSmoother == "pre" ) pre_or_post = ML_PRESMOOTHER;
    else if( PreOrPostSmoother == "both" ) pre_or_post = ML_BOTH;
    else 
      cerr << ErrorMsg_ << "smoother not recognized (" << PreOrPostSmoother << ")\n";
    
    sprintf(parameter,"smoother: type (level %d)", LevelID_[level]);
    Smoother = List_.get(parameter,Smoother);

    char msg[80];
    sprintf(msg,"Smoother (level %d) : ", LevelID_[level]);

    { // minor information about matrix size on each level
      int local[2];
      int global[2];
      local[0] = ml_->Amat[LevelID_[level]].invec_leng;
      local[1] = ml_->Amat[LevelID_[level]].N_nonzeros;
      Comm().SumAll(local,global,2);
      // kludge because it appears that the nonzeros for ML
      // defined operators are only local, while for Epetra
      // are global.
      if (level == 0)
        global[1] = local[1];
      if (verbose_)
        cout << msg << "# global rows = " << global[0] 
             << ", # estim. global nnz = " << global[1] << endl;
      // FIXME? above `est' is estimated, because there is something
      // wrong with the way nonzeros are stored in operators.
    }

    if( Smoother == "Jacobi" ) {

      // ============ //
      // point Jacobi //
      // ============ //

      if( verbose_ ) cout << msg << "Jacobi (sweeps="
			 << num_smoother_steps << ",omega=" << omega << ","
			 << PreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_Jacobi(ml_, LevelID_[level], pre_or_post,
			     num_smoother_steps, omega);
     
    } else if( Smoother == "Gauss-Seidel" ) {

      // ================== //
      // point Gauss-Seidel //
      // ================== //

      if( verbose_ ) cout << msg << "Gauss-Seidel (sweeps="
			 << num_smoother_steps << ",omega=" << omega << ","
			 << PreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_GaussSeidel(ml_, LevelID_[level], pre_or_post,
				  num_smoother_steps, omega);

    } else if( Smoother == "symmetric Gauss-Seidel" ) {

      // ====================== //
      // symmetric Gauss-Seidel //
      // ====================== //

      if( verbose_ ) cout << msg << "symmetric Gauss-Seidel (sweeps="
			  << num_smoother_steps << ",omega=" << omega << ","
			  << PreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_SymGaussSeidel(ml_, LevelID_[level], pre_or_post,
				     num_smoother_steps, omega);
    } else if( Smoother == "block Gauss-Seidel" ) {

      // ================== //
      // block Gauss-Seidel //
      // ================== //
      
      if( verbose_ ) cout << msg << "block Gauss-Seidel (sweeps="
			  << num_smoother_steps << ",omega=" << omega << ","
			  << PreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_BlockGaussSeidel(ml_, LevelID_[level], pre_or_post,
				       num_smoother_steps, omega, NumPDEEqns_);
    } else if( Smoother == "symmetric block Gauss-Seidel" ) {

      // ================== //
      // block Gauss-Seidel //
      // ================== //
      
      if( verbose_ ) cout << msg << "symmetric block Gauss-Seidel (sweeps="
			  << num_smoother_steps << ",omega=" << omega << ","
			  << PreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_SymBlockGaussSeidel(ml_, LevelID_[level], pre_or_post,
				       num_smoother_steps, omega, NumPDEEqns_);
    } else if( Smoother == "MLS" ) {

      // === //
      // MLS //
      // === //

      int logical_level = LevelID_[level];
      sprintf(parameter,"smoother: MLS polynomial order (level %d)",
              logical_level);
      MLSPolynomialOrder = List_.get(parameter,MLSPolynomialOrder);
      if (verbose_) 
        if (MLSPolynomialOrder > 0)
        {
          cout << msg << "MLS/Chebyshev, polynomial order = " << MLSPolynomialOrder
               << ", alpha = " << alpha << ", " << PreOrPostSmoother << endl;
        }
        else
        {
          cout << msg << "MLS, polynomial order = " << -MLSPolynomialOrder
               << ", alpha = " << alpha << ", " << PreOrPostSmoother << endl;
        }

      sprintf(parameter,"smoother: MLS alpha (level %d)",logical_level);
      MLSalpha = List_.get(parameter,MLSalpha);

      if (SolvingMaxwell_) {
        int Nfine_edge = Tmat_array[logical_level]->outvec_leng;
        int itmp, Ncoarse_edge;
        int coarsest_level = ml_->ML_coarsest_level;
        double edge_coarsening_rate=0.0;
        ML_gsum_scalar_int(&Nfine_edge, &itmp, ml_->comm);
        if (logical_level != coarsest_level) {
          Ncoarse_edge = Tmat_array[logical_level-1]->outvec_leng;
          ML_gsum_scalar_int(&Ncoarse_edge, &itmp, ml_->comm);
          if (Ncoarse_edge != 0.0)
            edge_coarsening_rate =  2.*((double) Nfine_edge)/
                                    ((double) Ncoarse_edge);
        }
        printf("level %d, before, edge_coarsening_rate = %e\n",
                logical_level, edge_coarsening_rate);
        if (edge_coarsening_rate < MLSalpha)
          edge_coarsening_rate =  MLSalpha;

/*
        printf("\n\nlevel %d, MLSalpha = %e, edge_coarsening_rate = %e,
  Nfine_edge = %d, Ncoarse_edge = %d\n\n",
               logical_level, MLSalpha, edge_coarsening_rate, Nfine_edge,
               Ncoarse_edge);
*/

        MLSalpha = edge_coarsening_rate;
      } //if (SolvingMaxwell_)

      ML_Gen_Smoother_MLS(ml_, LevelID_[level], pre_or_post,
                          MLSalpha, MLSPolynomialOrder);

      if (verbose_) {
        ML_Operator* this_A = &(ml_->Amat[LevelID_[level]]);
	cout << msg << "lambda_min = " << this_A->lambda_min
	     << ", lambda_max = " << this_A->lambda_max << endl;
      }
    } else if( Smoother == "Aztec" ) {
      
#ifdef HAVE_ML_AZTECOO
      // ======= //
      // AztecOO //
      // ======= //
      
      sprintf(parameter,"smoother: Aztec options (level %d)", LevelID_[level]);
      SmootherOptionsPtr = List_.get(parameter, SmootherOptionsPtr);
      sprintf(parameter,"smoother: Aztec params (level %d)", LevelID_[level]);
      SmootherParamsPtr = List_.get(parameter, SmootherParamsPtr);
      sprintf(parameter,"smoother: Aztec as solver (level %d)", LevelID_[level]);
      AztecSmootherAsASolver = List_.get(parameter,AztecSmootherAsASolver);
     
      if( AztecSmootherAsASolver == false ) aztec_its = AZ_ONLY_PRECONDITIONER;
      else                                  aztec_its = num_smoother_steps;
      
      if( SmootherOptionsPtr == NULL || SmootherParamsPtr == NULL ) {
      	SmootherOptionsPtr = SmootherOptions_;
      	SmootherParamsPtr = SmootherParams_;
      }

      if( verbose_ ) {
	cout << msg << "Aztec";
	if( SmootherOptionsPtr[AZ_precond] == AZ_dom_decomp ) {
	  cout << " DD, overlap=" << SmootherOptionsPtr[AZ_overlap] << ", ";
	  if( SmootherOptionsPtr[AZ_reorder] == 1 ) cout << "reord, ";
	  else cout << "no reord, ";
	  switch( SmootherOptionsPtr[AZ_subdomain_solve] ) {
	  case AZ_lu: cout << " LU"; break;
	  case AZ_ilu:
	    cout << "ILU(fill="  << SmootherOptionsPtr[AZ_graph_fill] << ")";
	    break;
	  case AZ_ilut:
	    cout << "ILUT(fill=" << SmootherParamsPtr[AZ_ilut_fill] << ",drop="
		 << SmootherParamsPtr[AZ_drop] << ")";
	    break;
	  case AZ_icc:
	    cout << "ICC(fill="  << SmootherOptionsPtr[AZ_graph_fill] << ")";
	    break;
	  case AZ_bilu:
	    cout << "BILU(fill="  << SmootherOptionsPtr[AZ_graph_fill] << ")";
	    break;
	  case AZ_rilu:
	    cout << "RILU(fill="  << SmootherOptionsPtr[AZ_graph_fill] << ",omega="
		 << SmootherParamsPtr[AZ_omega] << ")";
	    break;
	  }
	} else if( SmootherOptionsPtr[AZ_precond] == AZ_Jacobi ) {
	  cout << msg << " Jacobi preconditioner";
	} else if( SmootherOptionsPtr[AZ_precond] == AZ_Neumann ) {
	  cout << msg << " Neumann preconditioner, order = " << SmootherOptionsPtr[AZ_poly_ord];
	} else if( SmootherOptionsPtr[AZ_precond] == AZ_ls ) {
	  cout << msg << " LS preconditioner, order = " << SmootherOptionsPtr[AZ_poly_ord];
	} else if( SmootherOptionsPtr[AZ_precond] == AZ_sym_GS ) {
	  cout << msg << " symmetric Gauss-Seidel preconditioner, sweeps = " << SmootherOptionsPtr[AZ_poly_ord];
	} else if( SmootherOptionsPtr[AZ_precond] == AZ_none ) {
	  cout << msg << " with no preconditioning";
	}
	cout << ", "  << PreOrPostSmoother << endl;
      }
      
      ML_Gen_SmootherAztec(ml_, LevelID_[level], SmootherOptionsPtr, SmootherParamsPtr,
			   ProcConfig_, SmootherStatus_,
			   aztec_its, pre_or_post, NULL);
      
#else
      cerr << "Please configure ML with --enable-aztecoo to use" << endl;
      cerr << "AztecOO smoothers" << endl;
      exit(EXIT_FAILURE);
#endif
    } else if( Smoother == "IFPACK" || Smoother == "ILU" || Smoother == "IC") {

      // ====== //
      // IFPACK //
      // ====== //

#ifdef HAVE_ML_IFPACK
      int lof = -1;
      if (Smoother == "IFPACK")
      {
        sprintf(parameter,"smoother: ifpack type (level %d)", LevelID_[level]);
        IfpackType = List_.get(parameter, IfpackType);
      }
      else 
      {
        // MS // ILU and IC added on 08-Aug-06 for WebTrilinos
        // MS // Just a shortcut because sublists are not supported by
        // MS // the web interface.
        IfpackType = Smoother;
        lof = List_.get("smoother: ifpack level-of-fill", 0);
      }

      sprintf(parameter,"smoother: ifpack overlap (level %d)", LevelID_[level]);
      IfpackOverlap = List_.get(parameter, IfpackOverlap);

      if( verbose_ ) {
	cout << msg << "IFPACK, type = `" << IfpackType << "', " << endl
	     << msg << PreOrPostSmoother
	     << ", Overlap = " << IfpackOverlap << endl;
      }

      Teuchos::ParameterList& IfpackList = List_.sublist("smoother: ifpack list");
      int NumAggr = ML_Aggregate_Get_AggrCount(agg_,level);
      int* AggrMap = 0;
      ML_CHK_ERR(ML_Aggregate_Get_AggrMap(agg_,level,&AggrMap));

      // set these in the case the user wants "partitioner: type" = "user"
      // (if not, these values are ignored).
      if (IfpackList.get("partitioner: type", "user") == "user")
        IfpackList.set("partitioner: local parts", NumAggr);
      IfpackList.set("partitioner: map", AggrMap);

      if (lof != -1)
        IfpackList.set("fact: level-of-fill", lof);
                       
      ML_Gen_Smoother_Ifpack(ml_, IfpackType.c_str(),
                             IfpackOverlap, LevelID_[level], pre_or_post,
                             IfpackList,*Comm_);
      
#else
      cerr << ErrorMsg_ << "IFPACK not available." << endl
	   << ErrorMsg_ << "ML must be configured with --enable-ifpack" << endl
	   << ErrorMsg_ << "to use IFPACK as a smoother" << endl
	   << ErrorMsg_ << "NO SMOOTHER SET FOR THIS LEVEL" << endl;
#endif

    } else if( Smoother == "Chebyshev" ) {

#ifdef HAVE_ML_IFPACK
      if (SolvingMaxwell_)
        ML_CHK_ERR(-1); // not supported at this point

      sprintf(parameter,"smoother: polynomial order (level %d)",
              LevelID_[level]);
      PolynomialOrder = List_.get(parameter, PolynomialOrder);

      if( verbose_ ) {
	cout << msg << "IFPACK Chebyshev, order = " << PolynomialOrder
	     << ", alpha = " << alpha << ", " << PreOrPostSmoother << endl;
      }

      ML_Operator* this_A = &(ml_->Amat[LevelID_[level]]);
      ML_Gimmie_Eigenvalues(this_A, ML_DIAGSCALE, 
                  this_A->spectral_radius_scheme, ml_->symmetrize_matrix);

      Teuchos::ParameterList IFPACKList;
      IFPACKList.set("chebyshev: ratio eigenvalue", alpha);
      IFPACKList.set("chebyshev: min eigenvalue", this_A->lambda_min);
      IFPACKList.set("chebyshev: max eigenvalue", this_A->lambda_max);
      IFPACKList.set("chebyshev: degree", PolynomialOrder);

      if( verbose_ ) {
	cout << msg << "lambda_min = " << this_A->lambda_min
	     << ", lambda_max = " << this_A->lambda_max << endl;
      }

      ML_Gen_Smoother_Ifpack(ml_, "Chebyshev", 0, LevelID_[level], 
                             pre_or_post, IFPACKList, *Comm_);
#else
      cerr << ErrorMsg_ << "IFPACK not available." << endl
	   << ErrorMsg_ << "ML must be configured with --enable-ifpack" << endl
	   << ErrorMsg_ << "to use IFPACK as a smoother" << endl
	   << ErrorMsg_ << "NO SMOOTHER SET FOR THIS LEVEL" << endl;
#endif
    } else if( Smoother == "self" ) {

#ifdef HAVE_ML_IFPACK
      int IfpackOverlap = List_.get("smoother: self overlap",0);
      
      if( verbose_ ) {
	cout << msg << "ML as self-smoother, " << endl
	     << msg << PreOrPostSmoother
	     << ", Overlap = " << IfpackOverlap << endl;
      }

      Teuchos::ParameterList& SelfList = List_.sublist("smoother: self list");

      ML_Gen_Smoother_Self(ml_, IfpackOverlap, LevelID_[level], pre_or_post,
                           SelfList,*Comm_);
      
#else
      cerr << ErrorMsg_ << "IFPACK not available." << endl
	   << ErrorMsg_ << "ML must be configured with --enable-ifpack" << endl
	   << ErrorMsg_ << "to use ML as a smoother" << endl
	   << ErrorMsg_ << "NO SMOOTHER SET FOR THIS LEVEL" << endl;
#endif

    } else if( Smoother == "ParaSails" ) {

      // ========= //
      // ParaSails //
      // ========= //

      sprintf(parameter,"smoother: ParaSails levels (level %d)", level);
      ParaSailsN = List_.get(parameter,ParaSailsN);

      sprintf(parameter,"smoother: ParaSails matrix (level %d)", level);
      ParaSailsSym = List_.get(parameter,ParaSailsSym);

      sprintf(parameter,"smoother: ParaSails threshold (level %d)", level);
      ParaSailsThresh = List_.get(parameter,ParaSailsThresh);

      sprintf(parameter,"smoother: ParaSails filter (level %d)", level);
      ParaSailsFilter = List_.get(parameter,ParaSailsFilter);

      sprintf(parameter,"smoother: ParaSails load balancing (level %d)", level);
      ParaSailsLB = List_.get(parameter,ParaSailsLB);

      sprintf(parameter,"smoother: ParaSails factorized (level %d)", level);
      ParaSailsFactorized = List_.get(parameter,ParaSailsFactorized);

      if( verbose_ ) 
	cout << msg << "ParaSails "
	     << "(n=" << ParaSailsN
	     << ",sym=" << ParaSailsSym 
	     << ",thresh=" << ParaSailsThresh 
	     << ",filter=" << ParaSailsFilter 
	     << ",lb=" << ParaSailsLB
	     << ")" << endl;
      
#ifdef HAVE_ML_PARASAILS
      // I am not sure about the ending `0' and of ML
      ML_Gen_Smoother_ParaSails(ml_, LevelID_[level], 
				pre_or_post, num_smoother_steps,
				ParaSailsSym, ParaSailsThresh,
				ParaSailsN,
				ParaSailsFilter, (int)ParaSailsLB, 
				ParaSailsFactorized);
#else
      cerr << ErrorMsg_ << "ParaSails not available." << endl
	   << ErrorMsg_ << "ML must be configured with --with-ml_parasails" << endl
	   << ErrorMsg_ << "to use ParaSails as a smoother" << endl
	   << ErrorMsg_ << "NO SMOOTHER SET FOR THIS LEVEL" << endl;
#endif

    } else if( Smoother == "Hiptmair" ) {
      // ======== //
      // Hiptmair //
      // ======== //
      if (SolvingMaxwell_ == false) {
        if (Comm().MyPID() == 0) {
          cerr << ErrorMsg_ << "Hiptmair smoothing is only supported" << endl;
          cerr << ErrorMsg_ << "for solving eddy current equations." << endl;
          cerr << ErrorMsg_ << "Choose another smoother." << endl;
        }
        ML_EXIT(EXIT_FAILURE);
      }

      double subsmOmega;
      if (Comm().NumProc() == 1) subsmOmega = 1.0;
      else                       subsmOmega = ML_DDEFAULT;
      subsmOmega = List_.get("subsmoother: damping factor",subsmOmega);
      sprintf(parameter,"subsmoother: type (level %d)", level);
      SubSmootherType = List_.get(parameter,SubSmootherType);

      int logical_level = LevelID_[level];
      void *edge_smoother = 0, *nodal_smoother = 0;
      double node_coarsening_rate=0.0;
      double edge_coarsening_rate=0.0;

      nodal_its = List_.get("subsmoother: node sweeps", nodal_its);
      edge_its = List_.get("subsmoother: edge sweeps", edge_its);

      // only MLS and SGS are currently supported
      if (SubSmootherType == "MLS")
      {
        int subMLSPolynomialOrder =
            List_.get("subsmoother: MLS polynomial order",3);
        MLSalpha = List_.get("subsmoother: MLS alpha",27.0);
        sprintf(parameter,"subsmoother: MLS polynomial order (level %d)",level);
        MLSPolynomialOrder = List_.get(parameter,subMLSPolynomialOrder);

        nodal_smoother=(void *) ML_Gen_Smoother_MLS;
        ML_Smoother_Arglist_Set(nodal_args_, 0, &MLSPolynomialOrder);
        edge_smoother=(void *) ML_Gen_Smoother_MLS;
        ML_Smoother_Arglist_Set(edge_args_, 0, &MLSPolynomialOrder);

        // FIXME:  T could be NULL
        int Nfine_edge = Tmat_array[logical_level]->outvec_leng;
        int itmp, Ncoarse_edge;
        int coarsest_level = ml_->ML_coarsest_level;
        ML_gsum_scalar_int(&Nfine_edge, &itmp, ml_->comm);
        if (logical_level != coarsest_level) {
          Ncoarse_edge = Tmat_array[logical_level-1]->outvec_leng;
          ML_gsum_scalar_int(&Ncoarse_edge, &itmp, ml_->comm);
          if (Ncoarse_edge != 0.0)
            edge_coarsening_rate =  2.*((double) Nfine_edge)/
                                    ((double) Ncoarse_edge);
        }
        else edge_coarsening_rate =  0.0;
        if (edge_coarsening_rate < MLSalpha)
          edge_coarsening_rate =  MLSalpha;

        int Nfine_node = Tmat_array[logical_level]->invec_leng;
        int Ncoarse_node;
        ML_Smoother_Arglist_Set(edge_args_, 1, &edge_coarsening_rate);
        if (logical_level != coarsest_level) {
          Ncoarse_node = Tmat_array[logical_level-1]->invec_leng;
          ML_gsum_scalar_int(&Ncoarse_node, &itmp, ml_->comm);
          if (Ncoarse_node != 0.0)
            node_coarsening_rate =  2.*((double) Nfine_node)/
                                       ((double) Ncoarse_node);
        }
        else node_coarsening_rate = 0.0;
        if (node_coarsening_rate < MLSalpha)
          node_coarsening_rate =  MLSalpha;
                                                                                
        ML_Smoother_Arglist_Set(nodal_args_, 1, &node_coarsening_rate);
      }
      else if (SubSmootherType == "symmetric Gauss-Seidel") {
        sprintf(parameter,"subsmoother: damping factor (level %d)",
                logical_level);
        subsmOmega = List_.get(parameter,subsmOmega);
        nodal_smoother=(void *) ML_Gen_Smoother_SymGaussSeidel;
        ML_Smoother_Arglist_Set(nodal_args_, 0, &nodal_its);
        ML_Smoother_Arglist_Set(nodal_args_, 1, &subsmOmega);
        edge_smoother=(void *) ML_Gen_Smoother_SymGaussSeidel;
        ML_Smoother_Arglist_Set(edge_args_, 0, &edge_its);
        ML_Smoother_Arglist_Set(edge_args_, 1, &subsmOmega);
      }
      else if (Comm().MyPID() == 0)
        cerr << ErrorMsg_ << "Only MLS and SGS are supported as Hiptmair subsmoothers." << endl;
    
      int hiptmair_type = (int)
                  List_.get("smoother: Hiptmair efficient symmetric", true);
        
      ML_Gen_Smoother_Hiptmair(ml_, logical_level, ML_BOTH,
                 num_smoother_steps, Tmat_array, Tmat_trans_array, NULL, 
                 MassMatrix_array,
                 edge_smoother, edge_args_, nodal_smoother, nodal_args_,
                 hiptmair_type);

      bool indefiniteProblem = List_.get("negative conductivity",false);
      /* This corrects for the case of negative sigma */
      if (indefiniteProblem && SubSmootherType == "MLS")
      {
        if (verbose_ && Comm().MyPID() == 0)
          cout << "ML*WRN* "
             << "Resetting nodal smoother on level " << logical_level << endl
             << "ML*WRN* to account for negative mass matrix." << endl;
        //pre-smoother

        ML_Sm_Hiptmair_Data *hiptmairSmData =
           (ML_Sm_Hiptmair_Data *) ml_->pre_smoother[logical_level].smoother->data;
        ML *ml_subproblem = hiptmairSmData->ml_nodal;
                                                                                
        struct MLSthing *widget =
                       (struct MLSthing *) ml_subproblem->pre_smoother->smoother->data;
        double eig_ratio = widget->eig_ratio;
        int degree = widget->mlsDeg;
        ml_subproblem->pre_smoother->data_destroy(
                                 ml_subproblem->pre_smoother->smoother->data);
        ML_Operator *Amat = ml_subproblem->Amat;
        double tmp = Amat->lambda_max;
        Amat->lambda_max = fabs(Amat->lambda_min);
        Amat->lambda_min = fabs(tmp);
                                                                                
        ML_Gen_Smoother_MLS(ml_subproblem,0,ML_PRESMOOTHER,eig_ratio,degree);
                                                                                
        //post-smoother
        hiptmairSmData = (ML_Sm_Hiptmair_Data *)
                          ml_->post_smoother[logical_level].smoother->data;
        ml_subproblem = hiptmairSmData->ml_nodal;
                                                                                
        // Note:  this is correct because the pre_smoother is the only one
        // used in the subproblem
        widget = (struct MLSthing *) ml_subproblem->pre_smoother->smoother->data;
        eig_ratio = widget->eig_ratio;
        degree = widget->mlsDeg;
        ml_subproblem->pre_smoother->data_destroy(
               ml_subproblem->pre_smoother->smoother->data);
        Amat = ml_subproblem->Amat;
                                                                                
        ML_Gen_Smoother_MLS(ml_subproblem,0,ML_PRESMOOTHER,eig_ratio,degree);
      }

    } else if( Smoother == "user-defined" || Smoother == "user defined" ) {

      // ============ //
      // user-defined //
      // ============ //

      int (*userSmootherPtr)(ML_Smoother *, int, double *, int, double *);
      userSmootherPtr = NULL;
      userSmootherPtr = List_.get("smoother: user-defined function",
                                  userSmootherPtr);
      string userSmootherName;
      userSmootherName = List_.get("smoother: user-defined name",
                                   "User-defined");

      if( verbose_ ) cout << msg << userSmootherName << " (sweeps=" 
			 << num_smoother_steps << "," << PreOrPostSmoother << ")" << endl;

      if (userSmootherPtr == NULL) {
        if (Comm().MyPID() == 0)
          cerr << ErrorMsg_
               << "No pointer to user-defined smoother function found." << endl;
        ML_EXIT(EXIT_FAILURE);
      }
      ML_Operator *data;
      ML_Get_Amatrix(ml_, LevelID_[level], &data);
      ML_Set_Smoother(ml_, LevelID_[level], pre_or_post , data,
                      userSmootherPtr,
                      const_cast<char *>(userSmootherName.c_str()));

    } else if( Smoother == "do-nothing" ) {

      // ======================== //
      // do-nothing (no smoother) //
      // ======================== //

      if( verbose_ ) cout << msg << "do-nothing smoother" << endl;

    } else {

      // ======================================= //
      // error: no smoother found with this name //
      // ======================================= //

      if (Comm().MyPID() == 0)
         cerr << ErrorMsg_
              << "Smoother '" << Smoother << "' not recognized!" << endl
              << ErrorMsg_
              << "(file " << __FILE__ << ",line " << __LINE__ << ")" << endl
              << ErrorMsg_
              << "Now is: " << Smoother << ". It should be: " << endl
	          << ErrorMsg_
              << "<Jacobi> / <Gauss-Seidel> / <block Gauss-Seidel>" << endl
	          << ErrorMsg_
              << "<symmetric Gauss-Seidel> / <Aztec> / <IFPACK>" << endl
	          << ErrorMsg_ << "<MLS> / <ParaSails>" << endl;
      ML_EXIT(-99); }
    
    if( verbose_ ) 
      cout << msg << "Setup time : " << Time.ElapsedTime() << " (s)" << endl;
    
  } /* for */

  if(  verbose_ ) cout << endl;

  return(0);
}

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/

