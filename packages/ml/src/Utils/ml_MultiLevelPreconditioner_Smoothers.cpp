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
 * - \c MLS
 * - \c Aztec
 * - \c IFPACK
 * - \c do-nothing
 */
int ML_Epetra::MultiLevelPreconditioner::SetSmoothers() 
{

  char parameter[80];
  Epetra_Time Time(Comm());

  int num_smoother_steps = List_.get("smoother: sweeps", 1);

  double omega = List_.get("smoother: damping factor",0.67);

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

  int MLSPolynomialOrder = List_.get("smoother: MLS polynomial order",3);

  double MLSalpha = List_.get("smoother: MLS alpha",30.0);
  
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

  string SubSmootherType;
  int nodal_its, edge_its;
  if (SolvingMaxwell_ == true) {
    sprintf(parameter,"smoother: Hiptmair subsmoother type");
    SubSmootherType = List_.get(parameter,"MLS");
    sprintf(parameter,"smoother: Hiptmair node sweeps");
    nodal_its = List_.get(parameter, 1);
    sprintf(parameter,"smoother: Hiptmair edge sweeps");
    edge_its = List_.get(parameter, 1);
  }

  ML *ml;
  if (SolvingMaxwell_ == true)
    ml = ml_edges_;
  else
    ml = ml_;

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
      local[0] = ml->Amat[LevelID_[level]].invec_leng;
      local[1] = ml->Amat[LevelID_[level]].N_nonzeros;
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
      ML_Gen_Smoother_Jacobi(ml, LevelID_[level], pre_or_post,
			     num_smoother_steps, omega);
     
    } else if( Smoother == "Gauss-Seidel" ) {

      // ================== //
      // point Gauss-Seidel //
      // ================== //

      if( verbose_ ) cout << msg << "Gauss-Seidel (sweeps="
			 << num_smoother_steps << ",omega=" << omega << ","
			 << PreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_GaussSeidel(ml, LevelID_[level], pre_or_post,
				  num_smoother_steps, omega);

    } else if( Smoother == "symmetric Gauss-Seidel" ) {

      // ====================== //
      // symmetric Gauss-Seidel //
      // ====================== //

      if( verbose_ ) cout << msg << "symmetric Gauss-Seidel (sweeps="
			  << num_smoother_steps << ",omega=" << omega << ","
			  << PreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_SymGaussSeidel(ml, LevelID_[level], pre_or_post,
				     num_smoother_steps, omega);
    } else if( Smoother == "block Gauss-Seidel" ) {

      // ================== //
      // block Gauss-Seidel //
      // ================== //
      
      if( verbose_ ) cout << msg << "block Gauss-Seidel (sweeps="
			  << num_smoother_steps << ",omega=" << omega << ","
			  << PreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_BlockGaussSeidel(ml, LevelID_[level], pre_or_post,
				       num_smoother_steps, omega, NumPDEEqns_);
    } else if( Smoother == "MLS" ) {

      // === //
      // MLS //
      // === //

      sprintf(parameter,"smoother: MLS polynomial order (level %d)", LevelID_[level]);
      if( verbose_ ) cout << msg << "MLS,"
			 << PreOrPostSmoother << endl;
      
      ML_Gen_Smoother_MLS(ml, LevelID_[level], pre_or_post, MLSalpha,
			  MLSPolynomialOrder);
      
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
      
      ML_Gen_SmootherAztec(ml, LevelID_[level], SmootherOptionsPtr, SmootherParamsPtr,
			   ProcConfig_, SmootherStatus_,
			   aztec_its, pre_or_post, NULL);
      
#else
      cerr << "Please configure ML with --enable-aztecoo to use" << endl;
      cerr << "AztecOO smoothers" << endl;
      exit(EXIT_FAILURE);
#endif
    } else if( Smoother == "IFPACK" ) {

      // ====== //
      // IFPACK //
      // ====== //

#ifdef HAVE_ML_IFPACK
      string IfpackType = List_.get("smoother: ifpack type", "Amesos");
      int IfpackOverlap = List_.get("smoother: ifpack overlap",0);
      
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
      IfpackList.set("partitioner: local parts", NumAggr);
      IfpackList.set("partitioner: map", AggrMap);
      double Omega = IfpackList.get("relaxation: damping factor", 1.0);

      ML_Gen_Smoother_Ifpack(ml, IfpackType.c_str(),
                             IfpackOverlap, LevelID_[level], pre_or_post,
                             IfpackList,*Comm_);
      
      // omega = -1.0 means compute it through Anasazi, using
      // 10 iterations and 1e-5 as tolerance
#ifdef HAVE_ML_ANASAZIzzz
      if (Omega == -1.0) {
        // the following is broken
        //
        if (Omega == -1.0)
          IfpackList.set("relaxation: damping factor", 1.0);
        double LambdaMax;
        ML_Anasazi_Get_SpectralNorm_Anasazi(&(ml->Amat[LevelID_[level]]),
                                            &(ml->post_smoother[LevelID_[level]]),
                                            10, 1e-5, false, false, &LambdaMax);

        // compute optimal damping parameter
        Omega = 1.0 / LambdaMax;
        IfpackList.set("relaxation: damping factor", Omega);

        // some crap to re-set the damping parameter
        if (pre_or_post == ML_PRESMOOTHER || pre_or_post == ML_BOTH) {
          Ifpack_Preconditioner* Ifp = 
            (Ifpack_Preconditioner*)(ml->pre_smoother[LevelID_[level]].smoother->data);
          assert (Ifp != 0);
          Ifp->SetParameters(IfpackList);
        }
        if (pre_or_post == ML_POSTSMOOTHER || pre_or_post == ML_BOTH) {
          Ifpack_Preconditioner* Ifp = 
            (Ifpack_Preconditioner*)(ml->post_smoother[LevelID_[level]].smoother->data);
          assert (Ifp != 0);
          Ifp->SetParameters(IfpackList);
        }

        if (verbose_)
          cout << msg << "new damping parameter = " << Omega 
               << " (= 1 / " << LambdaMax << ")" << endl; 

        // set to -1 for the next level
        IfpackList.set("relaxation: damping factor", -1.0);
#else
	cerr << ErrorMsg_ << "Please compile with --enable-anasazi" << endl;
	cerr << ErrorMsg_ << "to use `relaxation: damping factor' == -1.0" << endl;
	exit(EXIT_FAILURE);
#endif
      }

#else
      cerr << ErrorMsg_ << "IFPACK not available." << endl
	   << ErrorMsg_ << "ML must be configure with --enable-ifpack" << endl
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

      ML_Gen_Smoother_Self(ml, IfpackOverlap, LevelID_[level], pre_or_post,
                           SelfList,*Comm_);
      
#else
      cerr << ErrorMsg_ << "IFPACK not available." << endl
	   << ErrorMsg_ << "ML must be configure with --enable-ifpack" << endl
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
      ML_Gen_Smoother_ParaSails(ml, LevelID_[level], 
				pre_or_post, num_smoother_steps,
				ParaSailsSym, ParaSailsThresh,
				ParaSailsN,
				ParaSailsFilter, (int)ParaSailsLB, 
				ParaSailsFactorized);
#else
      cerr << ErrorMsg_ << "ParaSails not available." << endl
	   << ErrorMsg_ << "ML must be configure with --with-ml_parasails" << endl
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

      sprintf(parameter,"smoother: Hiptmair subsmoother type (level %d)", level);
      SubSmootherType = List_.get(parameter,SubSmootherType);

      int logical_level = LevelID_[level];
      void *edge_smoother = 0, *nodal_smoother = 0;

      sprintf(parameter,"smoother: Hiptmair node sweeps");
      nodal_its = List_.get(parameter, nodal_its);
      sprintf(parameter,"smoother: Hiptmair edge sweeps");
      edge_its = List_.get(parameter, edge_its);

      // only MLS and SGS are currently supported
      if (SubSmootherType == "MLS")
      {
        sprintf(parameter,"smoother: Hiptmair MLS polynomial order (level %d)", level);
        MLSPolynomialOrder = List_.get(parameter,3);

        nodal_smoother=(void *) ML_Gen_Smoother_MLS;
        ML_Smoother_Arglist_Set(nodal_args_, 0, &MLSPolynomialOrder);
        edge_smoother=(void *) ML_Gen_Smoother_MLS;
        ML_Smoother_Arglist_Set(edge_args_, 0, &MLSPolynomialOrder);

        // FIXME:  T could be NULL
        int Nfine_edge = Tmat_array[logical_level]->outvec_leng;
        int itmp, Ncoarse_edge;
        int coarsest_level = ml_edges_->ML_coarsest_level;
        double node_coarsening_rate, edge_coarsening_rate;
        ML_gsum_scalar_int(&Nfine_edge, &itmp, ml_edges_->comm);
        if (logical_level != coarsest_level) {
          Ncoarse_edge = Tmat_array[logical_level-1]->outvec_leng;
          ML_gsum_scalar_int(&Ncoarse_edge, &itmp, ml_edges_->comm);
          edge_coarsening_rate =  2.*((double) Nfine_edge)/
                                  ((double) Ncoarse_edge);
        }
        else edge_coarsening_rate =  (double) Nfine_edge;
                                                                                
        int Nfine_node = Tmat_array[logical_level]->invec_leng;
        int Ncoarse_node;
        ML_Smoother_Arglist_Set(edge_args_, 1, &edge_coarsening_rate);
        if (logical_level != coarsest_level) {
          Ncoarse_node = Tmat_array[logical_level-1]->invec_leng;
          ML_gsum_scalar_int(&Ncoarse_node, &itmp, ml_edges_->comm);
          node_coarsening_rate =  2.*((double) Nfine_node)/
                                     ((double) Ncoarse_node);
        }
        else node_coarsening_rate = (double) Nfine_node;
                                                                                
        ML_Smoother_Arglist_Set(nodal_args_, 1, &node_coarsening_rate);
      }
      else if (SubSmootherType == "symmetric Gauss-Seidel") {
        sprintf(parameter,"smoother: Hiptmair SGS damping factor (level %d)",
                level);
        omega = List_.get(parameter,1.0);
        nodal_smoother=(void *) ML_Gen_Smoother_SymGaussSeidel;
        ML_Smoother_Arglist_Set(nodal_args_, 0, &nodal_its);
        ML_Smoother_Arglist_Set(nodal_args_, 1, &omega);
        edge_smoother=(void *) ML_Gen_Smoother_SymGaussSeidel;
        ML_Smoother_Arglist_Set(edge_args_, 0, &edge_its);
        ML_Smoother_Arglist_Set(edge_args_, 1, &omega);
      }
      else if (Comm().MyPID() == 0)
        cerr << ErrorMsg_ << "Only MLS and SGS are supported as Hiptmair subsmoothers." << endl;
    
      sprintf(parameter,"smoother: Hiptmair efficient symmetric");
      int hiptmair_type = (int) List_.get(parameter, true);
      sprintf(parameter,"smoother: Hiptmair sweeps");
      int num_smoother_sweeps = List_.get(parameter, 1);
        
      ML_Gen_Smoother_Hiptmair(ml_edges_, logical_level, ML_BOTH,
                 num_smoother_sweeps, Tmat_array, Tmat_trans_array, NULL, 
                 edge_smoother, edge_args_, nodal_smoother, nodal_args_,
                 hiptmair_type);

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
	cerr << ErrorMsg_ << "Smoother not recognized!" << endl
	     << ErrorMsg_ << "(file " << __FILE__ << ",line " << __LINE__ << ")" << endl
	     << ErrorMsg_ << "Now is: " << Smoother << ". It should be: " << endl
	     << ErrorMsg_ << "<Jacobi> / <Gauss-Seidel> / <block Gauss-Seidel>" << endl
	     << ErrorMsg_ << "<symmetric Gauss-Seidel> / <Aztec> / <IFPACK>" << endl
	     << ErrorMsg_ << "<MLS> / <ParaSails>" << endl;
      ML_EXIT(-99);
    }
    
    if( verbose_ ) 
      cout << msg << "Setup time : " << Time.ElapsedTime() << " (s)" << endl;
    
  } /* for */

  if(  verbose_ ) cout << endl;

  return(0);
}

// ============================================================================

/*! Options for Maxwell smoothing:
 *  - \c "smoother: sweeps" (int)
 *  - \c "smoother: type" (string). Possible values:
 *    - \c MLS
 *    - \c Gauss-Seidel
 *  - \c "smoother: node: sweeps" (int): sweeps for nodal hierarchy
 *  - \c "smoother: node: damping factor" (double)
 *  - \c "smoother: edge: sweeps" (int): sweeps for edge hierarchy
 *  - \c "smoother:  edge: damping factor" (double)
 *  
 * This function is still under development. The following options are
 * hardwired:
 * - HALF_HIPTMAIR
 *   
 */

int ML_Epetra::MultiLevelPreconditioner::SetSmoothersMaxwell()
{

  char parameter[80];
  
  int num_smoother_sweeps = List_.get("smoother: sweeps", 1);
  string SmootherType = List_.get("smoother: type","Hiptmair");

  // get user's defined parameters for nodal smoothing
  int nodal_its = List_.get("smoother: node: sweeps", 1);
  double nodal_omega = List_.get("smoother: node: damping factor",1.0);


  // get user's defined parameters for edge smoothing
  int edge_its = List_.get("smoother: edge: sweeps", 1);
  double edge_omega = List_.get("smoother: edge: damping factor",1.0);

  // arguments for nodal smoother
  nodal_args_ = ML_Smoother_Arglist_Create(2);
  int hiptmair_type = (int) List_.get("smoother: half Hiptmair", true);
  num_smoother_sweeps = List_.get(parameter, 1);
  int SmootherLevels = (NumLevels_>1)?(NumLevels_-1):1;
    
  // set the smoother (only two choices at this time)

  SmootherType = "Gauss-Seidel";

  if( SmootherType == "MLS" ) {
#ifdef LATER    
    int coarsest_level = -1;
    cout << "FIX COARSEST_LEVEL!!!" << endl;
    int Ncoarse_edge, Ncoarse_node;
    double edge_coarsening_rate, node_coarsening_rate;
    int Nfine_edge, itmp;

    for( int level=0 ; level<SmootherLevels ; ++level ) {

      if (level != coarsest_level) {
	Ncoarse_edge = Tmat_array[level-1]->outvec_leng;
	ML_gsum_scalar_int(&Ncoarse_edge, &itmp, ml_edges_->comm);
	edge_coarsening_rate =  2.*((double) Nfine_edge)/
                                   ((double) Ncoarse_edge);
      }
      else edge_coarsening_rate =  (double) Nfine_edge;

      ML_Smoother_Arglist_Set(edge_args_, 1, &edge_coarsening_rate);
      Nfine_edge = Ncoarse_edge;

      if (level != coarsest_level) {
	Ncoarse_node = Tmat_array[level-1]->invec_leng;
	ML_gsum_scalar_int(&Ncoarse_node, &itmp, ml_edges_->comm);
	node_coarsening_rate =  2.*((double) Nfine_node)/ 
                                   ((double) Ncoarse_node);
      }
      else node_coarsening_rate = (double) Nfine_node;

      ML_Smoother_Arglist_Set(nodal_args_, 1, &node_coarsening_rate);
      Nfine_node = Ncoarse_node;
    }

    ML_Gen_Smoother_Hiptmair(ml_edges_, LevelID_[level], ML_BOTH, num_smoother_sweeps,
			     Tmat_array, Tmat_trans_array, NULL, 
			     (void *)ML_Gen_Smoother_MLS, edge_args, 
			     (void *)ML_Gen_Smoother_MLS, nodal_args, hiptmair_type);
#endif
  } else if( SmootherType == "Gauss-Seidel") {

    ML_Smoother_Arglist_Set(nodal_args_, 0, &nodal_its);
    ML_Smoother_Arglist_Set(nodal_args_, 1, &nodal_omega);

    // arguments for edge smoother
    edge_args_ = ML_Smoother_Arglist_Create(2);

    ML_Smoother_Arglist_Set(edge_args_, 0, &edge_its);
    ML_Smoother_Arglist_Set(edge_args_, 1, &edge_omega);

    for( int level=0 ; level<SmootherLevels ; ++level ) {

      ML_Gen_Smoother_Hiptmair(ml_edges_, LevelID_[level], ML_BOTH, num_smoother_sweeps,
			       Tmat_array, Tmat_trans_array, NULL, 
			       (void *)ML_Gen_Smoother_SymGaussSeidel, edge_args_, 
			       (void *)ML_Gen_Smoother_SymGaussSeidel, nodal_args_, hiptmair_type);
    }
  } else {
    cerr << ErrorMsg_ << "`smoother: type' has an incorrect value (" << SmootherType << ")" << endl
         << ErrorMsg_ << "For Maxwell, it should be" << endl
	 << ErrorMsg_ << "<MLS> / <Gauss-Seidel>" << endl;
    exit( EXIT_FAILURE );
  }

  return(0);
}

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/

