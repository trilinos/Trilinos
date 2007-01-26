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

  int num_smoother_steps = List_.get("smoother: sweeps", 2);

  double omega = List_.get("smoother: damping factor",1.0);

  int pre_or_post = 0;
  string PreOrPostSmoother = List_.get("smoother: pre or post","both");

  string Smoother = List_.get("smoother: type","Chebyshev");

#ifdef HAVE_ML_AZTECOO
  int* SmootherOptionsPtr = (int*)0;
  SmootherOptionsPtr = List_.get("smoother: Aztec options",SmootherOptionsPtr);

  double* SmootherParamsPtr = (double*)0;
  SmootherParamsPtr = List_.get("smoother: Aztec params",SmootherParamsPtr);

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
  int nodal_its = 1, edge_its = 1;
  if (SolvingMaxwell_ == true) {
    SubSmootherType = List_.get("subsmoother: type","MLS");
    nodal_its = List_.get("subsmoother: node sweeps", 2);
    edge_its = List_.get("subsmoother: edge sweeps", 2);
  }

  // ===================== //
  // cycle over all levels //
  // ===================== //
 
  for (int level = 0 ; level < SmootherLevels ; ++level) {

    if (verbose_) cout << endl;

    Time.ResetStartTime();

    // general parameters for more than one smoother

    sprintf(parameter,"smoother: sweeps (level %d)", LevelID_[level] );
    int Mynum_smoother_steps = List_.get(parameter,num_smoother_steps);

    sprintf(parameter,"smoother: damping factor (level %d)", LevelID_[level] );
    double Myomega = List_.get(parameter,omega);

    sprintf(parameter,"smoother: pre or post (level %d)", LevelID_[level] );
    string MyPreOrPostSmoother = List_.get(parameter, PreOrPostSmoother);
    
    if( MyPreOrPostSmoother      == "post" ) pre_or_post = ML_POSTSMOOTHER;
    else if( MyPreOrPostSmoother == "pre"  ) pre_or_post = ML_PRESMOOTHER;
    else if( MyPreOrPostSmoother == "both" ) pre_or_post = ML_BOTH;
    else 
      cerr << ErrorMsg_ << "smoother not recognized (" << MyPreOrPostSmoother << ")\n";
    
    sprintf(parameter,"smoother: type (level %d)", LevelID_[level]);
    string MySmoother = List_.get(parameter,Smoother);

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

    if( MySmoother == "Jacobi" ) {

      // ============ //
      // point Jacobi //
      // ============ //

      if( verbose_ ) cout << msg << "Jacobi (sweeps="
			 << Mynum_smoother_steps << ",omega=" << Myomega << ","
			 << MyPreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_Jacobi(ml_, LevelID_[level], pre_or_post,
			     Mynum_smoother_steps, Myomega);
     
    } else if( MySmoother == "Gauss-Seidel" ) {

      // ================== //
      // point Gauss-Seidel //
      // ================== //

      if( verbose_ ) cout << msg << "Gauss-Seidel (sweeps="
			 << Mynum_smoother_steps << ",omega=" << Myomega << ","
			 << MyPreOrPostSmoother << ")" << endl;
#ifdef HAVE_ML_IFPACK
      if (ml_->Amat[LevelID_[level]].type == ML_TYPE_CRS_MATRIX) {
        if (verbose_)
          cout << msg << "Epetra_CrsMatrix detected, using "
               << "Ifpack implementation" << endl;
        string MyIfpackType = "point relaxation stand-alone";
        ParameterList& MyIfpackList = List_.sublist("smoother: ifpack list");;
        MyIfpackList.set("relaxation: type", "Gauss-Seidel");
        MyIfpackList.set("relaxation: sweeps", Mynum_smoother_steps);
        MyIfpackList.set("relaxation: damping factor", Myomega);
        ML_Gen_Smoother_Ifpack(ml_, MyIfpackType.c_str(),
                               IfpackOverlap, LevelID_[level], pre_or_post,
                               MyIfpackList,*Comm_);
      }
      else
#endif
      ML_Gen_Smoother_GaussSeidel(ml_, LevelID_[level], pre_or_post,
				  Mynum_smoother_steps, Myomega);

    } else if( Smoother == "ML Gauss-Seidel" ) {

      // ======================= //
      // ML's point Gauss-Seidel //
      // ======================= //

      if( verbose_ ) cout << msg << "Gauss-Seidel (sweeps="
                         << Mynum_smoother_steps << ",omega=" << Myomega << ","
                         << MyPreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_GaussSeidel(ml_, LevelID_[level], pre_or_post,
                                  Mynum_smoother_steps, Myomega);


    } else if( MySmoother == "symmetric Gauss-Seidel" ) {

      // ====================== //
      // symmetric Gauss-Seidel //
      // ====================== //

      if( verbose_ ) cout << msg << "symmetric Gauss-Seidel (sweeps="
			  << Mynum_smoother_steps << ",omega=" << Myomega << ","
			  << MyPreOrPostSmoother << ")" << endl;
#ifdef HAVE_ML_IFPACK
      if (ml_->Amat[LevelID_[level]].type == ML_TYPE_CRS_MATRIX) {
        if (verbose_)
          cout << msg << "Epetra_CrsMatrix detected, using "
               << "Ifpack implementation" << endl;
        string MyIfpackType = "point relaxation stand-alone";
        ParameterList& MyIfpackList = List_.sublist("smoother: ifpack list");;
        MyIfpackList.set("relaxation: type", "symmetric Gauss-Seidel");
        MyIfpackList.set("relaxation: sweeps", Mynum_smoother_steps);
        MyIfpackList.set("relaxation: damping factor", Myomega);
        ML_Gen_Smoother_Ifpack(ml_, MyIfpackType.c_str(),
                               IfpackOverlap, LevelID_[level], pre_or_post,
                               MyIfpackList,*Comm_);
      }
      else
#endif
      ML_Gen_Smoother_SymGaussSeidel(ml_, LevelID_[level], pre_or_post,
				     Mynum_smoother_steps, Myomega);

    } else if( Smoother == "ML symmetric Gauss-Seidel" ) {

      // =========================== //
      // ML's symmetric Gauss-Seidel //
      // ============================//

      if( verbose_ ) cout << msg << "ML symmetric Gauss-Seidel (sweeps="
                          << Mynum_smoother_steps << ",omega=" << Myomega << ","
                          << MyPreOrPostSmoother << ")" << endl;
        ML_Gen_Smoother_SymGaussSeidel(ml_, LevelID_[level], pre_or_post,
                                       Mynum_smoother_steps, Myomega);

    } else if( MySmoother == "block Gauss-Seidel" ) {

      // ================== //
      // block Gauss-Seidel //
      // ================== //
      
      if( verbose_ ) cout << msg << "block Gauss-Seidel (sweeps="
			  << Mynum_smoother_steps << ",omega=" << Myomega << ","
			  << MyPreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_BlockGaussSeidel(ml_, LevelID_[level], pre_or_post,
				       Mynum_smoother_steps, Myomega, NumPDEEqns_);

    } else if( MySmoother == "symmetric block Gauss-Seidel" ) {

      // ============================ //
      // symmetric block Gauss-Seidel //
      // ============================ //
      
      if( verbose_ ) cout << msg << "symmetric block Gauss-Seidel (sweeps="
			  << Mynum_smoother_steps << ",omega=" << Myomega << ","
			  << MyPreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_SymBlockGaussSeidel(ml_, LevelID_[level], pre_or_post,
				    Mynum_smoother_steps, Myomega, NumPDEEqns_);

    } else if( ( MySmoother == "MLS" ) || ( MySmoother == "Chebyshev" )) {

      // === //
      // MLS //
      // === //

      int logical_level = LevelID_[level];
      sprintf(parameter,"smoother: MLS polynomial order (level %d)",
              logical_level);
      int MyChebyshevPolyOrder = List_.get(parameter,ChebyshevPolyOrder);
      if (MyChebyshevPolyOrder == -97) {
         sprintf(parameter,"smoother: polynomial order (level %d)",
                 logical_level);
         MyChebyshevPolyOrder = List_.get(parameter,MyChebyshevPolyOrder);
      }
      if (MyChebyshevPolyOrder== -97) MyChebyshevPolyOrder=Mynum_smoother_steps;

      sprintf(parameter,"smoother: MLS alpha (level %d)",logical_level);
      double MyChebyshevAlpha = List_.get(parameter,ChebyshevAlpha);
      if (MyChebyshevAlpha == -2.) {
        sprintf(parameter,"smoother: Chebyshev alpha (level %d)",logical_level);
         MyChebyshevAlpha = List_.get(parameter,MyChebyshevAlpha);
      }
      if (MyChebyshevAlpha == -2.) MyChebyshevAlpha = 20.;

      if (verbose_) 
        if (MyChebyshevPolyOrder > 0)
        {
          cout << msg << "MLS/Chebyshev, polynomial order = " << MyChebyshevPolyOrder
               << ", alpha = " << MyChebyshevAlpha << ", " << MyPreOrPostSmoother << endl;
        }
        else
        {
          cout << msg << "MLS, polynomial order = " << -MyChebyshevPolyOrder
               << ", alpha = " << MyChebyshevAlpha << ", " << MyPreOrPostSmoother << endl;
        }


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
        if (edge_coarsening_rate < MyChebyshevAlpha)
          edge_coarsening_rate =  MyChebyshevAlpha;

        MyChebyshevAlpha = edge_coarsening_rate;
      } //if (SolvingMaxwell_)

      ML_Gen_Smoother_MLS(ml_, LevelID_[level], pre_or_post,
                          MyChebyshevAlpha, MyChebyshevPolyOrder);

      if (verbose_) {
        ML_Operator* this_A = &(ml_->Amat[LevelID_[level]]);
	cout << msg << "lambda_min = " << this_A->lambda_min
	     << ", lambda_max = " << this_A->lambda_max << endl;
      }
    } else if( MySmoother == "Aztec" ) {
      
#ifdef HAVE_ML_AZTECOO
      // ======= //
      // AztecOO //
      // ======= //
      
      sprintf(parameter,"smoother: Aztec options (level %d)", LevelID_[level]);
      int* MySmootherOptionsPtr = List_.get(parameter, SmootherOptionsPtr);
      sprintf(parameter,"smoother: Aztec params (level %d)", LevelID_[level]);
      double* MySmootherParamsPtr = List_.get(parameter, SmootherParamsPtr);
      sprintf(parameter,"smoother: Aztec as solver (level %d)", LevelID_[level]);
      bool MyAztecSmootherAsASolver = List_.get(parameter,AztecSmootherAsASolver);
     
      if( MyAztecSmootherAsASolver == false ) aztec_its = AZ_ONLY_PRECONDITIONER;
      else                                  aztec_its = Mynum_smoother_steps;
      
      if( MySmootherOptionsPtr == NULL || MySmootherParamsPtr == NULL ) {
      	MySmootherOptionsPtr = SmootherOptions_;
      	MySmootherParamsPtr = SmootherParams_;
      }

      if( verbose_ ) {
	cout << msg << "Aztec";
	if( MySmootherOptionsPtr[AZ_precond] == AZ_dom_decomp ) {
	  cout << " DD, overlap=" << MySmootherOptionsPtr[AZ_overlap] << ", ";
	  if( MySmootherOptionsPtr[AZ_reorder] == 1 ) cout << "reord, ";
	  else cout << "no reord, ";
	  switch( MySmootherOptionsPtr[AZ_subdomain_solve] ) {
	  case AZ_lu: cout << " LU"; break;
	  case AZ_ilu:
	    cout << "ILU(fill="  << MySmootherOptionsPtr[AZ_graph_fill] << ")";
	    break;
	  case AZ_ilut:
	    cout << "ILUT(fill=" << MySmootherParamsPtr[AZ_ilut_fill] << ",drop="
		 << MySmootherParamsPtr[AZ_drop] << ")";
	    break;
	  case AZ_icc:
	    cout << "ICC(fill="  << MySmootherOptionsPtr[AZ_graph_fill] << ")";
	    break;
	  case AZ_bilu:
	    cout << "BILU(fill="  << MySmootherOptionsPtr[AZ_graph_fill] << ")";
	    break;
	  case AZ_rilu:
	    cout << "RILU(fill="  << MySmootherOptionsPtr[AZ_graph_fill] << ",omega="
		 << MySmootherParamsPtr[AZ_omega] << ")";
	    break;
	  }
	} else if( MySmootherOptionsPtr[AZ_precond] == AZ_Jacobi ) {
	  cout << msg << " Jacobi preconditioner";
	} else if( MySmootherOptionsPtr[AZ_precond] == AZ_Neumann ) {
	  cout << msg << " Neumann preconditioner, order = " << MySmootherOptionsPtr[AZ_poly_ord];
	} else if( MySmootherOptionsPtr[AZ_precond] == AZ_ls ) {
	  cout << msg << " LS preconditioner, order = " << MySmootherOptionsPtr[AZ_poly_ord];
	} else if( MySmootherOptionsPtr[AZ_precond] == AZ_sym_GS ) {
	  cout << msg << " symmetric Gauss-Seidel preconditioner, sweeps = " << MySmootherOptionsPtr[AZ_poly_ord];
	} else if( MySmootherOptionsPtr[AZ_precond] == AZ_none ) {
	  cout << msg << " with no preconditioning";
	}
	cout << ", "  << MyPreOrPostSmoother << endl;
      }
      
      ML_Gen_SmootherAztec(ml_, LevelID_[level], MySmootherOptionsPtr, MySmootherParamsPtr,
			   ProcConfig_, SmootherStatus_,
			   aztec_its, pre_or_post, NULL);
      
#else
      cerr << "Please configure ML with --enable-aztecoo to use" << endl;
      cerr << "AztecOO smoothers" << endl;
      exit(EXIT_FAILURE);
#endif
    } else if( MySmoother == "IFPACK" || MySmoother == "ILU" || MySmoother == "IC") {

      // ====== //
      // IFPACK //
      // ====== //

#ifdef HAVE_ML_IFPACK
      int lof = -1;
      string MyIfpackType;
      if (MySmoother == "IFPACK")
      {
        sprintf(parameter,"smoother: ifpack type (level %d)", LevelID_[level]);
        MyIfpackType = List_.get(parameter, IfpackType);
      }
      else 
      {
        // MS // ILU and IC added on 08-Aug-06 for WebTrilinos
        // MS // Just a shortcut because sublists are not supported by
        // MS // the web interface.
        MyIfpackType = MySmoother;
        lof = List_.get("smoother: ifpack level-of-fill", 0);
      }

      sprintf(parameter,"smoother: ifpack overlap (level %d)", LevelID_[level]);
      int MyIfpackOverlap = List_.get(parameter, IfpackOverlap);

      if( verbose_ ) {
	cout << msg << "IFPACK, type = `" << MyIfpackType << "', " << endl
	     << msg << MyPreOrPostSmoother
	     << ", Overlap = " << MyIfpackOverlap << endl;
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
                       
      ML_Gen_Smoother_Ifpack(ml_, MyIfpackType.c_str(),
                             MyIfpackOverlap, LevelID_[level], pre_or_post,
                             IfpackList,*Comm_);
      
#else
      cerr << ErrorMsg_ << "IFPACK not available." << endl
	   << ErrorMsg_ << "ML must be configured with --enable-ifpack" << endl
	   << ErrorMsg_ << "to use IFPACK as a smoother" << endl
	   << ErrorMsg_ << "NO SMOOTHER SET FOR THIS LEVEL" << endl;
#endif

    } else if( MySmoother == "IFPACK-Chebyshev" ) {

#ifdef HAVE_ML_IFPACK
      int logical_level = LevelID_[level];

      if (SolvingMaxwell_)
        ML_CHK_ERR(-1); // not supported at this point

      sprintf(parameter,"smoother: MLS polynomial order (level %d)",
              logical_level);
      int MyChebyshevPolyOrder = List_.get(parameter,ChebyshevPolyOrder);
      if (MyChebyshevPolyOrder == -97) {
         sprintf(parameter,"smoother: polynomial order (level %d)",
                 logical_level);
         MyChebyshevPolyOrder = List_.get(parameter,MyChebyshevPolyOrder);
      }
      if (MyChebyshevPolyOrder== -97) MyChebyshevPolyOrder=Mynum_smoother_steps;

      sprintf(parameter,"smoother: MLS alpha (level %d)",logical_level);
      double MyChebyshevAlpha = List_.get(parameter,ChebyshevAlpha);
      if (MyChebyshevAlpha == -2.) {
        sprintf(parameter,"smoother: Chebyshev alpha (level %d)",logical_level);
         MyChebyshevAlpha = List_.get(parameter,MyChebyshevAlpha);
      }
      if (MyChebyshevAlpha == -2.) MyChebyshevAlpha = 20.;

      if( verbose_ ) {
	cout << msg << "IFPACK Chebyshev, order = " << MyChebyshevPolyOrder
	     << ", alpha = " << MyChebyshevAlpha << ", " << MyPreOrPostSmoother << endl;
      }

      ML_Operator* this_A = &(ml_->Amat[LevelID_[level]]);
      ML_Gimmie_Eigenvalues(this_A, ML_DIAGSCALE, 
                  this_A->spectral_radius_scheme, ml_->symmetrize_matrix);

      Teuchos::ParameterList IFPACKList;
      IFPACKList.set("chebyshev: ratio eigenvalue", MyChebyshevAlpha);
      IFPACKList.set("chebyshev: min eigenvalue", this_A->lambda_min);
      IFPACKList.set("chebyshev: max eigenvalue", this_A->lambda_max);
      IFPACKList.set("chebyshev: degree", MyChebyshevPolyOrder);

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
    } else if( MySmoother == "self" ) {

#ifdef HAVE_ML_IFPACK
      int MyIfpackOverlap = List_.get("smoother: self overlap",0);
      
      if( verbose_ ) {
	cout << msg << "ML as self-smoother, " << endl
	     << msg << MyPreOrPostSmoother
	     << ", Overlap = " << MyIfpackOverlap << endl;
      }

      Teuchos::ParameterList& SelfList = List_.sublist("smoother: self list");

      ML_Gen_Smoother_Self(ml_, MyIfpackOverlap, LevelID_[level], pre_or_post,
                           SelfList,*Comm_);
      
#else
      cerr << ErrorMsg_ << "IFPACK not available." << endl
	   << ErrorMsg_ << "ML must be configured with --enable-ifpack" << endl
	   << ErrorMsg_ << "to use ML as a smoother" << endl
	   << ErrorMsg_ << "NO SMOOTHER SET FOR THIS LEVEL" << endl;
#endif

    } else if( MySmoother == "ParaSails" ) {

      // ========= //
      // ParaSails //
      // ========= //

      sprintf(parameter,"smoother: ParaSails levels (level %d)", level);
      int MyParaSailsN = List_.get(parameter,ParaSailsN);

      sprintf(parameter,"smoother: ParaSails matrix (level %d)", level);
      int MyParaSailsSym = List_.get(parameter,ParaSailsSym);

      sprintf(parameter,"smoother: ParaSails threshold (level %d)", level);
      double MyParaSailsThresh = List_.get(parameter,ParaSailsThresh);

      sprintf(parameter,"smoother: ParaSails filter (level %d)", level);
      double MyParaSailsFilter = List_.get(parameter,ParaSailsFilter);

      sprintf(parameter,"smoother: ParaSails load balancing (level %d)", level);
      double MyParaSailsLB = List_.get(parameter,ParaSailsLB);

      sprintf(parameter,"smoother: ParaSails factorized (level %d)", level);
      int MyParaSailsFactorized = List_.get(parameter,ParaSailsFactorized);

      if( verbose_ ) 
	cout << msg << "ParaSails "
	     << "(n=" << MyParaSailsN
	     << ",sym=" << MyParaSailsSym 
	     << ",thresh=" << MyParaSailsThresh 
	     << ",filter=" << MyParaSailsFilter 
	     << ",lb=" << MyParaSailsLB
	     << "fact=" << MyParaSailsFactorized
	     << ")" << endl;
      
#ifdef HAVE_ML_PARASAILS
      // I am not sure about the ending `0' and of ML
      ML_Gen_Smoother_ParaSails(ml_, LevelID_[level], 
				pre_or_post, Mynum_smoother_steps,
				MyParaSailsSym, MyParaSailsThresh,
				MyParaSailsN,
				MyParaSailsFilter, (int) MyParaSailsLB, 
				MyParaSailsFactorized);
#else
      cerr << ErrorMsg_ << "ParaSails not available." << endl
	   << ErrorMsg_ << "ML must be configured with --with-ml_parasails" << endl
	   << ErrorMsg_ << "to use ParaSails as a smoother" << endl
	   << ErrorMsg_ << "NO SMOOTHER SET FOR THIS LEVEL" << endl;
#endif

    } else if( MySmoother == "Hiptmair" ) {
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
      string MySubSmootherType = List_.get(parameter,SubSmootherType);

      int logical_level = LevelID_[level];
      void *edge_smoother = 0, *nodal_smoother = 0;
      double node_coarsening_rate=0.0;
      double edge_coarsening_rate=0.0;

      sprintf(parameter,"subsmoother: node sweeps (level %d)", level);
      int Mynodal_its = List_.get(parameter, nodal_its);
      sprintf(parameter,"subsmoother: edge sweeps (level %d)", level);
      int Myedge_its = List_.get(parameter, edge_its);

      // only MLS and SGS are currently supported
      if ( (MySubSmootherType == "MLS") || (MySubSmootherType == "Chebyshev"))
      {
        // This is for backward compatibiltiy 
        int itemp = List_.get("subsmoother: MLS polynomial order",-97);
        if (itemp == -97) itemp=List_.get("subsmoother: polynomial order",-97);
        sprintf(parameter,"subsmoother: MLS polynomial order (level %d)",level);
        itemp = List_.get(parameter,itemp);
        if (itemp != -97) { Mynodal_its = itemp; Myedge_its = itemp; }

        double SubAlpha = List_.get("subsmoother: MLS alpha",-2.0);
        if (SubAlpha == -2.) SubAlpha=List_.get("subsmoother: Chebyshev alpha", -2.);
        if (SubAlpha == -2.) SubAlpha = 20.;


        nodal_smoother=(void *) ML_Gen_Smoother_MLS;
        ML_Smoother_Arglist_Set(nodal_args_, 0, &Mynodal_its);
        edge_smoother=(void *) ML_Gen_Smoother_MLS;
        ML_Smoother_Arglist_Set(edge_args_, 0, &Myedge_its);

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
        if (edge_coarsening_rate < SubAlpha)
          edge_coarsening_rate =  SubAlpha;

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
        if (node_coarsening_rate < SubAlpha)
          node_coarsening_rate =  SubAlpha;
                                                                                
        ML_Smoother_Arglist_Set(nodal_args_, 1, &node_coarsening_rate);
      }
      else if (MySubSmootherType == "symmetric Gauss-Seidel") {
        sprintf(parameter,"subsmoother: damping factor (level %d)",
                logical_level);
        double MysubsmOmega = List_.get(parameter,subsmOmega);
        nodal_smoother=(void *) ML_Gen_Smoother_SymGaussSeidel;
        ML_Smoother_Arglist_Set(nodal_args_, 0, &Mynodal_its);
        ML_Smoother_Arglist_Set(nodal_args_, 1, &MysubsmOmega);
        edge_smoother=(void *) ML_Gen_Smoother_SymGaussSeidel;
        ML_Smoother_Arglist_Set(edge_args_, 0, &Myedge_its);
        ML_Smoother_Arglist_Set(edge_args_, 1, &MysubsmOmega);
      }
      else if (Comm().MyPID() == 0)
        cerr << ErrorMsg_ << "Only Chebyshev (or MLS) and SGS are supported as "
             << "Hiptmair subsmoothers ... not " << MySubSmootherType << endl;

    
      int hiptmair_type = (int)
                  List_.get("smoother: Hiptmair efficient symmetric", true);
        
      ML_Gen_Smoother_Hiptmair(ml_, logical_level, ML_BOTH,
                 Mynum_smoother_steps, Tmat_array, Tmat_trans_array, NULL, 
                 MassMatrix_array,
                 edge_smoother, edge_args_, nodal_smoother, nodal_args_,
                 hiptmair_type);

      bool indefiniteProblem = List_.get("negative conductivity",false);
      /* This corrects for the case of negative sigma        */
      /* I think it has something to do with the eigenvalues */
      /* coming out negative.                                */
      if (indefiniteProblem && MySubSmootherType == "MLS")
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
        if ( Amat->lambda_max < 0.) {
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
                                                                                
           ML_Gen_Smoother_MLS(ml_subproblem,0,ML_PRESMOOTHER,eig_ratio,degree);
        }
      }

    } else if( MySmoother == "user-defined" || MySmoother == "user defined" ) {

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
			 << Mynum_smoother_steps << "," << MyPreOrPostSmoother << ")" << endl;

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

    } else if( MySmoother == "do-nothing" ) {

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
              << "Smoother '" << MySmoother << "' not recognized!" << endl
              << ErrorMsg_
              << "(file " << __FILE__ << ",line " << __LINE__ << ")" << endl
              << ErrorMsg_
              << "Now is: " << MySmoother << ". It should be: " << endl
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

