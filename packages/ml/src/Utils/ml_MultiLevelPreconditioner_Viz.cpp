/*!
 * }
 *  \file ml_MultiLevelPreconditioner_Viz.cpp
 *
 *  \brief Visualization utilities for MultiLevelPreconditioner class
 *
 *  \author Marzio Sala, SNL, 9214
 *
 *  \date Last update do Doxygen: 09-Aug-04
 *
 */

#include "ml_common.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

#include "ml_include.h"
#ifndef ML_CPP
#ifdef __cplusplus
extern "C" 
{
#endif
#endif
  
extern int ML_Aggregate_VizAndStats_Compute( ML *ml, ML_Aggregate *ag, int MaxMgLevels,
                                     double *x, double *y, double *z, int Ndimensions,
                                     char *base_filename );
extern int ML_Aggregate_Stats_ComputeCoordinates( ML *ml, ML_Aggregate *ag,
						 double *x, double *y, double *z);
extern int ML_Aggregate_Stats_Analyze( ML *ml, ML_Aggregate *ag);
extern int ML_Aggregate_Viz( ML *ml, ML_Aggregate *ag, int choice,
			    double * vector, char * base_filename, int level);
extern int ML_Aggregate_Stats_CleanUp_Amalgamate( ML *ml, ML_Aggregate *ag);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_FECrsMatrix.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "ml_epetra.h"
#include "ml_MultiLevelPreconditioner.h"

// ============================================================================
// Used in VizMePleaze()
void ML_Epetra::MultiLevelPreconditioner::
RandomAndZero(double * tmp_rhs, double * tmp_sol, int size)
{
  // create random numbers between 0.5 and 1.0
  ML_random_vec(tmp_rhs, size, ml_->comm);
  for( int i=0 ; i<size ; ++i ) tmp_rhs[i] = 0.5+0.25*tmp_rhs[i];
  for( int i=0 ; i<size ; ++i ) tmp_sol[i] = 0.0;
}

// ============================================================================
int ML_Epetra::MultiLevelPreconditioner::
VisualizeAggregates()
{
  if (IsPreconditionerComputed() == false)
    ML_CHK_ERR(-1); // need an already computed preconditioner

  ML_CHK_ERR(Visualize(true, false, false, false, -1, -1, -1));

  return(0);
}

// ============================================================================
int ML_Epetra::MultiLevelPreconditioner::
VisualizeSmoothers(int NumPreCycles, int NumPostCycles)
{

  if (IsPreconditionerComputed() == false)
    ML_CHK_ERR(-1); // need an already computed preconditioner

  bool VizPreSmoother = false;
  bool VizPostSmoother = false;

  if (NumPreCycles != 0) 
    VizPreSmoother = true;
  if (NumPostCycles != 0) 
    VizPostSmoother = true;

  int ierr = Visualize(false, VizPreSmoother, VizPostSmoother,
		       false, NumPreCycles, NumPostCycles, -1);

  ML_CHK_ERR(ierr);

  return(0);
}

// ============================================================================
int ML_Epetra::MultiLevelPreconditioner::
VisualizeCycle(int NumCycles)
{

  if (IsPreconditionerComputed() == false)
    ML_CHK_ERR(-1); // need an already computed preconditioner

  int ierr = Visualize(false, false, false, true,
		       -1, -1, NumCycles);

  ML_CHK_ERR(ierr);

  return(0);
}


// ============================================================================
// Visualize aggregates and (for XYZ format) also plot vectors
// date: Aug-04
int ML_Epetra::MultiLevelPreconditioner::
Visualize(bool VizAggre, bool VizPreSmoother,
	  bool VizPostSmoother, bool VizCycle,
	  int NumApplPreSmoother, int NumApplPostSmoother,
	  int NumCycleSmoother)
{

  char filename[80];

  // does not work with Maxwell
  if( ml_ == 0 ) {
    cerr << ErrorMsg_ << "Visualization does not work (yet) with Maxwell..." << endl;
    ML_CHK_ERR(-1);
  }

  string Prefix = Prefix_;

  int NumDimensions = 0;
  double * x_coord = List_.get(Prefix+"viz: x-coordinates", (double *)0);
  double * y_coord = List_.get(Prefix+"viz: y-coordinates", (double *)0);
  double * z_coord = List_.get(Prefix+"viz: z-coordinates", (double *)0);

  if( x_coord ) NumDimensions++;
  if( y_coord ) NumDimensions++;
  if( z_coord ) NumDimensions++;

  assert( NumDimensions != 0 );

  if (VizAggre == true) {

    // stats about level matrix sizes
    if( verbose_ ) 
      cout << endl << "- number of rows for each level's matrix:" << endl << endl;

    for( int ilevel=0 ; ilevel < NumLevels_ ; ++ilevel ) {
      int imin, iavg, imax;
      int Nrows = ml_->Amat[LevelID_[ilevel]].outvec_leng/NumPDEEqns_;
      Comm().MinAll(&Nrows,&imin,1);
      Comm().MaxAll(&Nrows,&imax,1);
      Comm().SumAll(&Nrows,&iavg,1); iavg /= Comm().NumProc();

      if( verbose_ ) {
	printf( "\t(level %d) rows per process (min) = %d\n", ilevel, imin);
	printf( "\t(level %d) rows per process (avg) = %d\n", ilevel, iavg);
	printf( "\t(level %d) rows per process (max) = %d\n", ilevel, imax);
	cout << endl;
      }
    } 
 
    if( verbose_ ) 
      cout << endl << "- analysis of the computational domain (finest level):" 
	<< endl << endl;
    ML_Aggregate_Stats_Analyze(ml_,agg_);

  }

  // prepare output format. Now it can be:
  // - OpenDX (1D/2D/3D)
  // - XD3D (2D only)

  int Format;
  string FileFormat = List_.get(Prefix+"viz: output format", "xyz");
  // you are a cool guy if you plot with "xyz"
  if( FileFormat == "xyz" ) Format = 1;
  // you are a poor man if you need "dx". God bless you.
  else if( FileFormat == "dx" ) Format = 0;
  else {
    cerr << ErrorMsg_ << "Option `viz: output format' has an incorrect" << endl
      << ErrorMsg_ << "value (" << FileFormat << "). Possible values are" << endl
      << ErrorMsg_ << "<dx> / <xyz>" << endl;
    exit( EXIT_FAILURE );
  }

  int ieqn             = List_.get(Prefix+"viz: equation to plot", -1);
  if( ieqn >= NumPDEEqns_ ) ieqn = 0;
  bool PrintStarting   = List_.get(Prefix+"viz: print starting solution", false);

  ML_Smoother * ptr;
  double * tmp_rhs = new double[NumMyRows()]; 
  double * tmp_sol = new double[NumMyRows()]; 
  double * plot_me = new double[NumMyRows()/NumPDEEqns_];
  
  // Note that this requires the new version of the
  // visualization routines. OpenDX cannot visualize vectors.

  if( ( VizPreSmoother || VizPostSmoother || VizCycle ) && ( Format == 0 ) ) {
    cerr << endl;
    cerr << ErrorMsg_ << "Option `viz: output format' == `dx' cannot be used" << endl
         << ErrorMsg_ << "to visualize the effect of smoothers and cycle." << endl;
    cerr << endl;
    VizPreSmoother = false;
    VizPostSmoother = false;
    VizCycle = false;
  }

  if( verbose_ ) 
    cout << endl << "- visualization:" << endl << endl;

  // =============================================================== //
  // cycle over all levels. Note that almost the same thing          //
  // is done for pre-smoothing, post-smoothing, and the effect       //
  // of the cycle itself. For each of these, I plot on file          //
  // the starting solution (before-), the final solution (after-),   //
  // for each equation, and for each level (smoother only).          //
  // All these junk works with XYZ only, and it should work in       //
  // 3D too (although I never tested in 3D).                         //
  // =============================================================== //

  for( int ilevel=0 ; ilevel<NumLevels_ ; ++ilevel ) {

    // =================== //
    // plot the aggregates //
    // =================== //

    if( VizAggre ) 
      ML_Aggregate_Viz(ml_,agg_,Format,NULL,NULL,LevelID_[ilevel]);

    // ============ //
    // pre-smoother //
    // ============ //

    ptr = ((ml_->SingleLevel[LevelID_[ilevel]]).pre_smoother);

    if( ptr != NULL && VizPreSmoother ) {

      RandomAndZero(tmp_sol, tmp_rhs,ml_->Amat[LevelID_[ilevel]].outvec_leng);

      // visualize starting vector
      if( PrintStarting ) {
	if( ieqn != -1 ) {
	  for( int i=0 ; i<NumMyRows() ; i+=NumPDEEqns_ ) 
	    plot_me[i/NumPDEEqns_] = tmp_sol[i+ieqn];
	  sprintf(filename,"before-presmoother-eq%d", ieqn);
	  ML_Aggregate_Viz(ml_,agg_,Format,plot_me,filename,LevelID_[ilevel]);
	} else { // by default, print out all equations
	  for( int eq=0 ; eq<NumPDEEqns_ ; eq++ ) {
	    for( int i=0 ; i<NumMyRows() ; i+=NumPDEEqns_ ) 
	      plot_me[i/NumPDEEqns_] = tmp_sol[i+eq];
	    sprintf(filename,"before-presmoother-eq%d", eq);
	    ML_Aggregate_Viz(ml_,agg_,Format,plot_me,filename,LevelID_[ilevel]);
	  }
	}
      }

      // increase the number of applications of the smoother
      // and run the smoother
      int old_ntimes = ptr->ntimes;
      ptr->ntimes = NumApplPreSmoother;
      ML_Smoother_Apply(ptr, 
			ml_->Amat[LevelID_[ilevel]].outvec_leng,
			tmp_sol,
			ml_->Amat[LevelID_[ilevel]].outvec_leng,
			tmp_rhs, ML_NONZERO);
      ptr->ntimes = old_ntimes;

      // visualize
      // user may have required one spefic equation only
      if( ieqn != -1 ) {
	for( int i=0 ; i<NumMyRows() ; i+=NumPDEEqns_ ) 
	  plot_me[i/NumPDEEqns_] = tmp_sol[i+ieqn];
	sprintf(filename,"after-presmoother-eq%d", ieqn);
	ML_Aggregate_Viz(ml_,agg_,Format,plot_me,filename,LevelID_[ilevel]);
      } else { // by default, print out all equations
	for( int eq=0 ; eq<NumPDEEqns_ ; eq++ ) {
	  for( int i=0 ; i<NumMyRows() ; i+=NumPDEEqns_ ) 
	    plot_me[i/NumPDEEqns_] = tmp_sol[i+eq];
	  sprintf(filename,"after-presmoother-eq%d", eq);
	  ML_Aggregate_Viz(ml_,agg_,Format,plot_me,filename,LevelID_[ilevel]);
	}
      }
    } // VizPreSmoother

    // ============= //
    // post-smoother //
    // ============= //

    ptr = ((ml_->SingleLevel[LevelID_[ilevel]]).post_smoother);
    if( ptr != NULL && VizPostSmoother ) {

      // random solution and 0 rhs
      RandomAndZero(tmp_sol, tmp_rhs,ml_->Amat[LevelID_[ilevel]].outvec_leng);

      // visualize starting vector
      if( PrintStarting ) {
	if( ieqn != -1 ) {
	  for( int i=0 ; i<NumMyRows() ; i+=NumPDEEqns_ ) 
	    plot_me[i/NumPDEEqns_] = tmp_sol[i+ieqn];
	  sprintf(filename,"before-postsmoother-eq%d", ieqn);
	  ML_Aggregate_Viz(ml_,agg_,Format,plot_me,filename,LevelID_[ilevel]);
	} else { // by default, print out all equations
	  for( int eq=0 ; eq<NumPDEEqns_ ; eq++ ) {
	    for( int i=0 ; i<NumMyRows() ; i+=NumPDEEqns_ ) 
	      plot_me[i/NumPDEEqns_] = tmp_sol[i+eq];
	    sprintf(filename,"before-postsmoother-eq%d", eq);
	    ML_Aggregate_Viz(ml_,agg_,Format,plot_me,filename,LevelID_[ilevel]);
	  }
	}
      }

      // increase the number of applications of the smoother
      // and run the smoother
      int old_ntimes = ptr->ntimes;
      ptr->ntimes = NumApplPostSmoother;
      ML_Smoother_Apply(ptr, 
			ml_->Amat[LevelID_[ilevel]].outvec_leng,
			tmp_sol,
			ml_->Amat[LevelID_[ilevel]].outvec_leng,
			tmp_rhs, ML_ZERO);
      ptr->ntimes = old_ntimes;

      // visualize
      // user may have required one spefic equation only
      if( ieqn != -1 ) {
	for( int i=0 ; i<NumMyRows() ; i+=NumPDEEqns_ ) 
	  plot_me[i/NumPDEEqns_] = tmp_sol[i+ieqn];
	printf(filename,"after-postsmoother-eq%d", ieqn);
	ML_Aggregate_Viz(ml_,agg_,Format,plot_me,filename,LevelID_[ilevel]);
      } else { // by default, print out all equations
	for( int eq=0 ; eq<NumPDEEqns_ ; eq++ ) {
	  for( int i=0 ; i<NumMyRows() ; i+=NumPDEEqns_ ) 
	    plot_me[i/NumPDEEqns_] = tmp_sol[i+eq];
	  sprintf(filename,"after-postsmoother-eq%d", eq);
	  ML_Aggregate_Viz(ml_,agg_,Format,plot_me,filename,LevelID_[ilevel]);
	}
      }
    } // VizPostSmoother 
  } // for( ilevel )

  // =============================== //
  // run ML cycle on a random vector //
  // =============================== //

  if( VizCycle ) {

    // random solution and zero rhs
    RandomAndZero(tmp_sol, tmp_rhs,ml_->Amat[LevelID_[0]].outvec_leng);

    // visualize starting vector
    if( PrintStarting ) {
      if( ieqn != -1 ) {
	for( int i=0 ; i<NumMyRows() ; i+=NumPDEEqns_ ) 
	  plot_me[i/NumPDEEqns_] = tmp_sol[i+ieqn];
	sprintf(filename,"before-cycle-eq%d", ieqn);
	ML_Aggregate_Viz(ml_,agg_,Format,plot_me,filename,LevelID_[0]);
      } else { // by default, print out all equations
	for( int eq=0 ; eq<NumPDEEqns_ ; eq++ ) {
	  for( int i=0 ; i<NumMyRows() ; i+=NumPDEEqns_ ) 
	    plot_me[i/NumPDEEqns_] = tmp_sol[i+eq];
	  sprintf(filename,"before-cycle-eq%d", eq);
	  ML_Aggregate_Viz(ml_,agg_,Format,plot_me,filename,LevelID_[0]);
	}
      }
    }

    // run the cycle
    for( int i=0 ; i<NumCycleSmoother ; ++i ) 
      ML_Cycle_MG(&(ml_->SingleLevel[ml_->ML_finest_level]),
		  tmp_sol, tmp_rhs,
		  ML_NONZERO, ml_->comm, ML_NO_RES_NORM, ml_);

    // visualize
    // user may have required one spefic equation only
    if( ieqn != -1 ) {
      for( int i=0 ; i<NumMyRows() ; i+=NumPDEEqns_ ) 
	plot_me[i/NumPDEEqns_] = tmp_sol[i+ieqn];
      sprintf(filename,"after-cycle-eq%d", ieqn);
      ML_Aggregate_Viz(ml_,agg_,Format,plot_me,filename,LevelID_[0]);
    } else { // by default, print out all equations
      for( int eq=0 ; eq<NumPDEEqns_ ; eq++ ) {
	for( int i=0 ; i<NumMyRows() ; i+=NumPDEEqns_ ) 
	  plot_me[i/NumPDEEqns_] = tmp_sol[i+eq];
	sprintf(filename,"after-cycle-eq%d", eq);
	ML_Aggregate_Viz(ml_,agg_,Format,plot_me,filename,LevelID_[0]);
      }
    }
  } // VizCycle

  // =================== //
  // clean up and return //
  // =================== //

  delete [] tmp_sol;
  delete [] tmp_rhs;
  delete [] plot_me;

  return(0);
}

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
