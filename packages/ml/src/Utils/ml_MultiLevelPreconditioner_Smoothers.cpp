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
#include "ml_MultiLevelPreconditioner.h"
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
  Epetra_Time Time(Comm());

  double smooTime=0, coarseTime=0, perLevelTime=0, totalTime=0;

  int num_smoother_steps = List_.get("smoother: sweeps", 2);

  double omega = List_.get("smoother: damping factor",1.0);

  int pre_or_post = 0;
  string PreOrPostSmoother = List_.get("smoother: pre or post","both");

  string Smoother = List_.get("smoother: type","Chebyshev");
  
#ifdef HAVE_ML_AZTECOO
  RCP<std::vector<int> > aztecOptions = List_.get("smoother: Aztec options",SmootherOptions_);
  RCP<std::vector<double> > aztecParams = List_.get("smoother: Aztec params",SmootherParams_);
  int* SmootherOptionsPtr = &(*aztecOptions)[0];
  double* SmootherParamsPtr = &(*aztecParams)[0];

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
  string IfpackType = List_.get("smoother: ifpack type", "Amesos");
  int IfpackOverlap = List_.get("smoother: ifpack overlap",0);
  // note: lof has different meanings for IC and ICT.  For IC and ILU, we
  // will cast it to an integer later.
  double IfpackLOF=List_.get("smoother: ifpack level-of-fill",0.);
  double IfpackRelThreshold=List_.get("smoother: ifpack relative threshold",1.);
  double IfpackAbsThreshold=List_.get("smoother: ifpack absolute threshold",0.);

  // Block Chebyshev parameters
  int cheby_nBlocks=List_.get("smoother: Block Chebyshev number of blocks",-1);
  int *cheby_blockIndices=List_.get("smoother: Block Chebyshev block list",(int*)0);
  int *cheby_blockStarts=List_.get("smoother: Block Chebyshev block starts",(int*)0);
  
  // Chebyshev-NE parameters
  bool cheby_NE=List_.get("smoother: chebyshev solve normal equations",false);

  // Hiptmair-specific declarations
  string SubSmType,NodeSubSmType,EdgeSubSmType;
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
  for (int level = 0 ; level < SmootherLevels ; ++level) {

    if (verbose_) cout << endl;

    Time.ResetStartTime();

    int coarseLevel  = LevelID_[NumLevels_-1];
    int currentLevel = LevelID_[level];

    // general parameters for more than one smoother

    if (currentLevel != coarseLevel)
      sprintf(smListName,"smoother: list (level %d)",currentLevel);
    else
      strcpy(smListName,"coarse: list");
    ParameterList &smList = List_.sublist(smListName);

    int Mynum_smoother_steps = smList.get("smoother: sweeps",num_smoother_steps);
    double Myomega = smList.get("smoother: damping factor",omega);

    string MyPreOrPostSmoother = smList.get("smoother: pre or post", PreOrPostSmoother);
    
    if( MyPreOrPostSmoother      == "post" ) pre_or_post = ML_POSTSMOOTHER;
    else if( MyPreOrPostSmoother == "pre"  ) pre_or_post = ML_PRESMOOTHER;
    else if( MyPreOrPostSmoother == "both" ) pre_or_post = ML_BOTH;
    else 
      cerr << ErrorMsg_ << "smoother not recognized (" << MyPreOrPostSmoother << ")\n";
    
    string MySmoother = smList.get("smoother: type",Smoother);

    char msg[80];
    double AddToDiag = smList.get("smoother: add to diag", 1e-12);
    int MaxProcs=-1;

    if (currentLevel != coarseLevel)
    {
      // smoothers
      sprintf(msg,"Smoother (level %d) : ", currentLevel);
      // minor information about matrix size on each level
      double local[2];
      double global[2];
      local[0] = ml_->Amat[currentLevel].invec_leng;
      local[1] = ml_->Amat[currentLevel].N_nonzeros;
      Comm().SumAll(local,global,2);
      if (verbose_) {
        int i = cout.precision(0);
        cout.setf(std::ios::fixed);
        cout << msg << "# global rows = " << global[0] 
             << ", # estim. global nnz = " << global[1] << endl;
        cout.precision(i);
        cout.unsetf(std::ios::fixed);
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

      if( verbose_ ) cout << msg << "Jacobi (sweeps="
                          << Mynum_smoother_steps << ",omega=" << Myomega << ","
                          << MyPreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_Jacobi(ml_, currentLevel, pre_or_post,
                             Mynum_smoother_steps, Myomega);
     
    } else if( MySmoother == "Gauss-Seidel" ) {

      // ================== //
      // point Gauss-Seidel //
      // ================== //

      if( verbose_ ) cout << msg << "Gauss-Seidel (sweeps="
                          << Mynum_smoother_steps << ",omega=" << Myomega << ","
                          << MyPreOrPostSmoother << ")" << endl;

      bool gs_type = List_.get("smoother: Gauss-Seidel efficient symmetric",false);
      
#ifdef HAVE_ML_IFPACK
      if (ml_->Amat[currentLevel].type == ML_TYPE_CRS_MATRIX) {
        if (verbose_)
          cout << msg << "Epetra_CrsMatrix detected, using "
               << "Ifpack implementation" << endl;
        string MyIfpackType = "point relaxation stand-alone";
        ParameterList& MyIfpackList = List_.sublist("smoother: ifpack list");;
        MyIfpackList.set("relaxation: type", "Gauss-Seidel");
        MyIfpackList.set("relaxation: sweeps", Mynum_smoother_steps);
        MyIfpackList.set("relaxation: damping factor", Myomega);

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
      }
      else
#endif

        if(gs_type)
          ML_Gen_Smoother_EffSymGaussSeidel(ml_, currentLevel, pre_or_post,
                                            Mynum_smoother_steps, Myomega);

        else
          ML_Gen_Smoother_GaussSeidel(ml_, currentLevel, pre_or_post,
                                      Mynum_smoother_steps, Myomega);

    } else if( MySmoother == "ML Gauss-Seidel" ) {

      // ======================= //
      // ML's point Gauss-Seidel //
      // ======================= //

      bool gs_type = List_.get("smoother: Gauss-Seidel efficient symmetric",false);
      if( verbose_ ) cout << msg << "Gauss-Seidel (sweeps="
                         << Mynum_smoother_steps << ",omega=" << Myomega << ","
                         << MyPreOrPostSmoother << ")" << endl;

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
      if( verbose_ ) cout << msg << "symmetric Gauss-Seidel (sweeps="
			  << Mynum_smoother_steps << ",omega=" << Myomega << ","
			  << MyPreOrPostSmoother << ")" << endl;
#ifdef HAVE_ML_IFPACK
      if (ml_->Amat[currentLevel].type == ML_TYPE_CRS_MATRIX) {
        if (verbose_)
          cout << msg << "Epetra_CrsMatrix detected, using "
               << "Ifpack implementation" << endl;
        string MyIfpackType = "point relaxation stand-alone";
        ParameterList& MyIfpackList = List_.sublist("smoother: ifpack list");;
        MyIfpackList.set("relaxation: type", "symmetric Gauss-Seidel");
        MyIfpackList.set("relaxation: sweeps", Mynum_smoother_steps);
        MyIfpackList.set("relaxation: damping factor", Myomega);
        ML_Gen_Smoother_Ifpack(ml_, MyIfpackType.c_str(),
                               IfpackOverlap, currentLevel, pre_or_post,
                               (void*)&MyIfpackList,(void*)Comm_);
      }
      else
#endif
      ML_Gen_Smoother_SymGaussSeidel(ml_, currentLevel, pre_or_post,
				     Mynum_smoother_steps, Myomega);

    } else if( MySmoother == "ML symmetric Gauss-Seidel" ) {

      // =========================== //
      // ML's symmetric Gauss-Seidel //
      // ============================//

      if( verbose_ ) cout << msg << "ML symmetric Gauss-Seidel (sweeps="
                          << Mynum_smoother_steps << ",omega=" << Myomega << ","
                          << MyPreOrPostSmoother << ")" << endl;
        ML_Gen_Smoother_SymGaussSeidel(ml_, currentLevel, pre_or_post,
                                       Mynum_smoother_steps, Myomega);

    } else if( MySmoother == "block Gauss-Seidel" ) {

      // ================== //
      // block Gauss-Seidel //
      // ================== //
      
      if( verbose_ ) cout << msg << "block Gauss-Seidel (sweeps="
			  << Mynum_smoother_steps << ",omega=" << Myomega << ","
			  << MyPreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_BlockGaussSeidel(ml_, currentLevel, pre_or_post,
				       Mynum_smoother_steps, Myomega, NumPDEEqns_);

    } else if( MySmoother == "symmetric block Gauss-Seidel" ) {

      // ============================ //
      // symmetric block Gauss-Seidel //
      // ============================ //
      
      if( verbose_ ) cout << msg << "symmetric block Gauss-Seidel (sweeps="
			  << Mynum_smoother_steps << ",omega=" << Myomega << ","
			  << MyPreOrPostSmoother << ")" << endl;
      ML_Gen_Smoother_SymBlockGaussSeidel(ml_, currentLevel, pre_or_post,
				    Mynum_smoother_steps, Myomega, NumPDEEqns_);

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
          cout << msg << "MLS/Block Chebyshev, polynomial order = "
               <<  MyChebyshevPolyOrder
               << ", alpha = " << MyChebyshevAlpha << ", "
               << MyPreOrPostSmoother << endl;

        }
        else if (MySmoother == "Chebyshev" && MyChebyshevPolyOrder > 0)
        {
          cout << msg << "MLS/Chebyshev, polynomial order = "
               <<  MyChebyshevPolyOrder
               << ", alpha = " << MyChebyshevAlpha << ", "
               << MyPreOrPostSmoother << endl;
        }
        else
        {
          cout << msg << "MLS, polynomial order = " << -MyChebyshevPolyOrder
               << ", alpha = " << MyChebyshevAlpha << ", "
               << MyPreOrPostSmoother << endl;
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
        cout << msg << "lambda_min = " << this_A->lambda_min
             << ", lambda_max = " << this_A->lambda_max << endl;
      }
    } else if( MySmoother == "Aztec" ) {
      
#ifdef HAVE_ML_AZTECOO
      // ======= //
      // AztecOO //
      // ======= //

      // These should remain int* and double*, rather than Teuchos::RCP's.
      // The user created the options & params arrays, so he is responsible for
      // freeing them.
//      int* MySmootherOptionsPtr         = smList.get("smoother: Aztec options", SmootherOptionsPtr);
RCP<std::vector<int> > myaztecOptions   = smList.get("smoother: Aztec options",SmootherOptions_);
RCP<std::vector<double> > myaztecParams = smList.get("smoother: Aztec params",SmootherParams_);
  int* MySmootherOptionsPtr = &(*myaztecOptions)[0];
  double* MySmootherParamsPtr = &(*myaztecParams)[0];
//      double* MySmootherParamsPtr = smList.get("smoother: Aztec params", SmootherParamsPtr);
      bool MyAztecSmootherAsASolver = smList.get("smoother: Aztec as solver",AztecSmootherAsASolver);
     
      if( MyAztecSmootherAsASolver == false ) aztec_its = AZ_ONLY_PRECONDITIONER;
      else                                  aztec_its = Mynum_smoother_steps;
      
      if( verbose_ ) {
	cout << msg << "Aztec";
	if( MyAztecSmootherAsASolver){
	  switch (MySmootherOptionsPtr[AZ_solver]){
	  case AZ_cg:
	  case AZ_cg_condnum:
	    cout<<"-CG";
	    break;
	  case AZ_gmres:
	  case AZ_gmres_condnum:
	    cout<<"-GMRES";
	    break;
	  case AZ_cgs:
	    cout<<"-CGS";
	    break;
	  case AZ_tfqmr:
	    cout<<"-TSQMR";
	    break;
	  case AZ_bicgstab:
	    cout<<"-TSQMR";
	    break;
	  case AZ_GMRESR:
	    cout<<"-GMRESR";
	    break;	    
	  }
	  cout<<"("<<aztec_its<<")";
	}

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
	  cout << " Jacobi preconditioner, sweeps = " << MySmootherOptionsPtr[AZ_poly_ord];
	} else if( MySmootherOptionsPtr[AZ_precond] == AZ_Neumann ) {
	  cout << " Neumann preconditioner, order = " << MySmootherOptionsPtr[AZ_poly_ord];
	} else if( MySmootherOptionsPtr[AZ_precond] == AZ_ls ) {
	  cout << " LS preconditioner, order = " << MySmootherOptionsPtr[AZ_poly_ord];
	} else if( MySmootherOptionsPtr[AZ_precond] == AZ_sym_GS ) {
	  cout << " symmetric Gauss-Seidel preconditioner, sweeps = " << MySmootherOptionsPtr[AZ_poly_ord];
	} else if( MySmootherOptionsPtr[AZ_precond] == AZ_none ) {
	  cout << " No preconditioning";
	}
	cout << ", "  << MyPreOrPostSmoother << endl;
      }
      
      ML_Gen_SmootherAztec(ml_, currentLevel, MySmootherOptionsPtr, MySmootherParamsPtr,
			   ProcConfig_, SmootherStatus_,
			   aztec_its, pre_or_post, NULL);
      
#else
      cerr << "Please configure ML with --enable-aztecoo to use" << endl;
      cerr << "AztecOO smoothers" << endl;
      exit(EXIT_FAILURE);
#endif
    } else if( MySmoother == "IFPACK" || MySmoother == "ILU" ||
               MySmoother == "IC" || MySmoother == "ILUT"   ||
               MySmoother == "ICT" || MySmoother == "SILU") {

      // ====== //
      // IFPACK //
      // ====== //

#ifdef HAVE_ML_IFPACK
      double lof = -1.0;
      string MyIfpackType;
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
        lof = List_.get("smoother: ifpack level-of-fill", IfpackLOF);
      }

      int MyIfpackOverlap = smList.get("smoother: ifpack overlap", IfpackOverlap);
      double MyIfpackRT = smList.get("smoother: ifpack relative threshold", IfpackRelThreshold);
      double MyIfpackAT = smList.get("smoother: ifpack absolute threshold", IfpackAbsThreshold);

/*
      if( verbose_ ) {
        cout << msg << "IFPACK, type=`" << MyIfpackType << "'," << endl
             << msg << MyPreOrPostSmoother
             << ",overlap=" << MyIfpackOverlap << endl;
        if (MyIfpackType != "Amesos") {
          if (MyIfpackType == "ILU" || MyIfpackType == "IC")
            cout << msg << "level-of-fill=" << (int)lof;
          else
            cout << msg << "level-of-fill=" << lof;
          cout << ",rel. threshold=" << MyIfpackRT
               << ",abs. threshold=" << MyIfpackAT << endl;
        }
      }
*/

      Teuchos::ParameterList& IfpackList=List_.sublist("smoother: ifpack list");
      int NumAggr = ML_Aggregate_Get_AggrCount(agg_,level);
      int* AggrMap = 0;
      ML_CHK_ERR(ML_Aggregate_Get_AggrMap(agg_,level,&AggrMap));
      IfpackList.set("ILU: sweeps", Mynum_smoother_steps);

      // set these in the case the user wants "partitioner: type" = "user"
      // (if not, these values are ignored).
      if (IfpackList.get("partitioner: type", "user") == "user")
        IfpackList.set("partitioner: local parts", NumAggr);
      IfpackList.set("partitioner: map", AggrMap);

      if (lof != -1.0) {
        IfpackList.set("fact: level-of-fill", (int) lof);
        IfpackList.set("fact: ict level-of-fill", lof);
        IfpackList.set("fact: ilut level-of-fill", lof);
      }
      IfpackList.set("fact: relative threshold", MyIfpackRT);
      IfpackList.set("fact: absolute threshold", MyIfpackAT);

      if( verbose_ ) {
	// SORa needs special handling
	if(MyIfpackType == "SORa"){
	  cout << msg << "IFPACK/SORa("<<IfpackList.get("sora: alpha",1.5)<<","<<IfpackList.get("sora: gamma",1.0)<<")"
	       << ", sweeps = " <<IfpackList.get("sora: sweeps",1)<<endl;
	  if(IfpackList.get("sora: oaz boundaries",false))
	    cout << msg << "oaz boundary handling enabled"<<endl;
	  if(IfpackList.get("sora: use interproc damping",false))
	    cout << msg << "interproc damping enabled"<<endl;
	  if(IfpackList.get("sora: use global damping",false))
	    cout << msg << "global damping enabled"<<endl;
	}
	else{
	  cout << msg << "IFPACK, type=`" << MyIfpackType << "'," << endl
	       << msg << MyPreOrPostSmoother
	       << ",overlap=" << MyIfpackOverlap << endl;
	  if (MyIfpackType != "Amesos") {
	    if (MyIfpackType == "ILU" || MyIfpackType == "IC") {
	      int myLof = IfpackList.get("fact: level-of-fill", (int) lof);
	      cout << msg << "level-of-fill=" << myLof;
	    }
	    else {
            double myLof = IfpackList.get("fact: level-of-fill", (int) lof);
            cout << msg << "level-of-fill=" << myLof;
	    }
	    cout << ",rel. threshold=" << MyIfpackRT
		 << ",abs. threshold=" << MyIfpackAT << endl;
	  }
	}
      }
      ML_Gen_Smoother_Ifpack(ml_, MyIfpackType.c_str(),
                             MyIfpackOverlap, currentLevel, pre_or_post,
                             (void*)&IfpackList,(void*)Comm_);
      
#else
      cerr << ErrorMsg_ << "IFPACK not available." << endl
	   << ErrorMsg_ << "ML must be configured with --enable-ifpack" << endl
	   << ErrorMsg_ << "to use IFPACK as a smoother" << endl
	   << ErrorMsg_ << "NO SMOOTHER SET FOR THIS LEVEL" << endl;
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
          cout << msg << "IFPACK Block Chebyshev, order = " << MyChebyshevPolyOrder
               << ", alpha = " << MyChebyshevAlpha << ", " << MyPreOrPostSmoother << endl;
        else
          cout << msg << "IFPACK Chebyshev, order = " << MyChebyshevPolyOrder
               << ", alpha = " << MyChebyshevAlpha << ", " << MyPreOrPostSmoother << endl;
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
	cout << msg << "lambda_min = " << this_A->lambda_min
	     << ", lambda_max = " << this_A->lambda_max << endl;
      }

#else
      cerr << ErrorMsg_ << "IFPACK not available." << endl
	   << ErrorMsg_ << "ML must be configured with --enable-ifpack" << endl
	   << ErrorMsg_ << "to use IFPACK as a smoother" << endl
	   << ErrorMsg_ << "NO SMOOTHER SET FOR THIS LEVEL" << endl;
#endif
    } else if( MySmoother == "self" ) {

#ifdef HAVE_ML_IFPACK
      int MyIfpackOverlap;
      if(smList.isParameter("smoother: self overlap"))
        MyIfpackOverlap = smList.get("smoother: self overlap",0);
      else
        MyIfpackOverlap = List_.get("smoother: self overlap",0);
      
      if( verbose_ ) {
        cout << msg << "ML as self-smoother ("
             << "cycles=" << Mynum_smoother_steps
             << ",overlap=" << MyIfpackOverlap << ","
             << MyPreOrPostSmoother << ")" << endl;
      }

      Teuchos::ParameterList IfpackList;
      Teuchos::ParameterList& SelfList = IfpackList.sublist("ML list");
      Teuchos::ParameterList& tmpList = List_.sublist("smoother: self list");
      SelfList.setParameters(tmpList);
      IfpackList.set( "ML node id",List_.get("ML node id",-1) );
      char procLabel[30];
      sprintf(procLabel,"node id %d",List_.get("ML node id",-1));
      SelfList.set("ML label",procLabel);
      SelfList.set("zero starting solution", false);  
      string xxx = SelfList.get("SetDefaults", "not-set");
      if (xxx != "not-set") {
        if (verbose_ && Comm().MyPID() == 0)
          cout << msg << "Setting default values to type `" << xxx << "'" << endl;
        SetDefaults(xxx, SelfList,0,0,false);
      }

      if (verbose_ && SelfList.get("ML output",0) > 0)
        cout << msg << "*** * Start of self-smoother generation * ***" << endl;
      int currentPrintLevel = ML_Get_PrintLevel();
      ML_Gen_Smoother_Self(ml_, MyIfpackOverlap, currentLevel, pre_or_post,
                           Mynum_smoother_steps, IfpackList,*Comm_);
      ML_Set_PrintLevel(currentPrintLevel);
      if (verbose_ && SelfList.get("ML output",0) > 0)
        cout << msg << "*** * End of self-smoother generation * ***" << endl;
      
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

      int MyParaSailsN = smList.get("smoother: ParaSails levels",ParaSailsN);

      int MyParaSailsSym = smList.get("smoother: ParaSails matrix",ParaSailsSym);

      double MyParaSailsThresh = smList.get("smoother: ParaSails threshold",ParaSailsThresh);

      double MyParaSailsFilter = smList.get("smoother: ParaSails filter",ParaSailsFilter);

      double MyParaSailsLB = smList.get("smoother: ParaSails load balancing",ParaSailsLB);

      int MyParaSailsFactorized = smList.get("smoother: ParaSails factorized",ParaSailsFactorized);

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
      ML_Gen_Smoother_ParaSails(ml_, currentLevel, 
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

      // ==================================================== //
      // Hiptmair                                             //
      // supported subsmoothers:                              //
      //   Chebyshev, SGS, Ifpack incomplete factorizations   //
      // ==================================================== //
      if (AMGSolver_ != ML_MAXWELL) {
        if (Comm().MyPID() == 0) {
          cerr << ErrorMsg_ << "Hiptmair smoothing is only supported" << endl;
          cerr << ErrorMsg_ << "for solving eddy current equations." << endl;
          cerr << ErrorMsg_ << "Choose another smoother." << endl;
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

      string MyEdgeSubSmType    = smList.get("subsmoother: edge type",EdgeSubSmType);
      string MyNodeSubSmType    = smList.get("subsmoother: node type",NodeSubSmType);
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

        string *MySubSmType=0;
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
          cerr << ErrorMsg_
            <<"Only Chebyshev (or MLS), SGS, ILU, IC, ILUT, and ICT" << endl
            << "are supported as Hiptmair subsmoothers ... not "
            << *MySubSmType << endl;
    
      } //for (enum nodeOrEdge ne=NODE; ne!=DONE ...


      int hiptmair_type = (int)
                  List_.get("smoother: Hiptmair efficient symmetric", true);

      if( verbose_ ) cout << msg << "Hiptmair (outer sweeps="
             << Mynum_smoother_steps << ")" << endl
             << msg << "edge: " << EdgeSmootherInfo << endl
             << msg << "node: " << NodeSmootherInfo << endl;
        
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
          cout << "ML*WRN* "
             << "Resetting nodal smoother on level " << thisLevel << endl
             << "ML*WRN* to account for negative mass matrix." << endl;
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

#     ifdef HAVE_PETSC

      
/*
      void *voidPC = 0;
      ML_PetscPC petscPC = 0;
      petscPC = (ML_PetscPC) smList.get("smoother: petsc pc", voidPC);
*/
      void *voidKSP = 0;
      ML_PetscKSP petscKSP = 0;
      petscKSP = (ML_PetscKSP) smList.get("smoother: petsc ksp", voidKSP);
      if (petscKSP == 0) {
        if (Comm().MyPID() == 0)
          cerr << ErrorMsg_
               << "You must provide a fully-constructed KSP context to use a PETSc smoother."
               << endl;
        exit(EXIT_FAILURE);
      }
      const char* pcName;
      ML_PetscPC petscPC;
      int ierr = KSPGetPC(petscKSP,&petscPC);CHKERRQ(ierr);
      ierr = PCGetType(petscPC,&pcName);
      if( verbose_ ) cout << msg << "PETSc smoother (type="
              << pcName
			  << ",sweeps=" << Mynum_smoother_steps << ","
			  << MyPreOrPostSmoother << ")" << endl;

      ML_Gen_Smoother_Petsc(ml_, currentLevel, pre_or_post, Mynum_smoother_steps, petscKSP);

#     else
      if (Comm().MyPID() == 0)
       cerr << ErrorMsg_
            << "You must configure ML with PETSc support enabled."
            << endl;
      exit(EXIT_FAILURE);
#     endif /*ifdef HAVE_PETSC*/
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
      if (Comm().MyPID() == 0)
       cerr << ErrorMsg_
            << "You must configure ML with Teko support enabled. "
            << "Enable flag ENABLE_TekoML and add Teko to the end of the library line"
            << endl;
      exit(EXIT_FAILURE);
#endif
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
                          << Mynum_smoother_steps << ","
                          << MyPreOrPostSmoother << ")" << endl;

      if (userSmootherPtr == NULL) {
        if (Comm().MyPID() == 0)
          cerr << ErrorMsg_
               << "No pointer to user-defined smoother function found." << endl;
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
              << "You chose: " << MySmoother << ". It should be: " << endl
	          << ErrorMsg_
              << "<Jacobi> / <Gauss-Seidel> / <block Gauss-Seidel>" << endl
	          << ErrorMsg_
              << "<symmetric Gauss-Seidel> / <Aztec> / <IFPACK>" << endl
	          << ErrorMsg_
              << "<Chebyshev> / <ParaSails> / <Hiptmair>" << endl
	          << ErrorMsg_ << "<user-defined>" << endl;
      ML_EXIT(-99); }
    
    perLevelTime = Time.ElapsedTime();
    if (currentLevel != coarseLevel) {
      smooTime += perLevelTime;
      if (verbose_)
        cout << msg << "Setup time : " << perLevelTime << " (s)" << endl;
    }
    else
      coarseTime = perLevelTime;
    
  } /* for (int level = 0 ; level < SmootherLevels ; ++level) */

  totalTime += (smooTime + coarseTime);
  OutputList_.set("time: smoothers setup", smooTime
                  + OutputList_.get("time: smoothers setup", 0.0));
  OutputList_.set("time: coarse solver setup", coarseTime
                  + OutputList_.get("time: coarse solver setup", 0.0));
  if(  verbose_ ) cout << endl;

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

