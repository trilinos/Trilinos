/*!
 *  \file ml_MultiLevelPreconditioner_Viz.cpp
 *
 *  \brief Visualization utilities for MultiLevelPreconditioner class
 *
 *  \author Marzio Sala, SNL, 9214
 *
 *  \date Last update do Doxygen: 06-Aug-04
 *
 */

#include "ml_common.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

#include "ml_include.h"
extern "C" {
  
extern int ML_Aggregate_VizAndStats_Setup(ML_Aggregate * agg_, int NumLevels_);
extern int ML_Aggregate_VizAndStats_Compute( ML *ml, ML_Aggregate *ag, int MaxMgLevels,
                                     double *x, double *y, double *z, int Ndimensions,
                                     char *base_filename );
extern int ML_Aggregate_VizAndStats_Setup(ML_Aggregate *, int);
extern int ML_Aggregate_Stats_ComputeCoordinates( ML *ml, ML_Aggregate *ag,
						 double *x, double *y, double *z);
extern int ML_Aggregate_Stats_Analyze( ML *ml, ML_Aggregate *ag);
extern int ML_Aggregate_Viz( ML *ml, ML_Aggregate *ag, int choice,
			    double * vector, char * base_filename, int level);
extern int ML_Aggregate_Stats_CleanUp_Info( ML *ml, ML_Aggregate *ag);
extern int ML_Aggregate_Stats_CleanUp_Amalgamate( ML *ml, ML_Aggregate *ag);

}
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
#include "ml_epetra_preconditioner.h"

// ============================================================================
void ML_Epetra::MultiLevelPreconditioner::RandomAndZero(double * tmp_rhs, double * tmp_sol, int size)
{
  // does not work with Maxwell
  if( ml_ == 0 ) exit( EXIT_FAILURE );

  // create random numbers between 0.5 and 1.0
  ML_random_vec(tmp_rhs, size, ml_->comm);
  for( int i=0 ; i<size ; ++i ) tmp_rhs[i] = 0.5+0.25*tmp_rhs[i];
  for( int i=0 ; i<size ; ++i ) tmp_sol[i] = 0.0;
}

// ============================================================================
// Visualize aggregates and (for XYZ format) also plot vectors
void ML_Epetra::MultiLevelPreconditioner::VizMePleaze()
{

  string Prefix = Prefix_;

  int NumDimensions = 0;
  double * x_coord = List_.get(Prefix+"viz: x-coordinates", (double *)0);
  double * y_coord = List_.get(Prefix+"viz: y-coordinates", (double *)0);
  double * z_coord = List_.get(Prefix+"viz: z-coordinates", (double *)0);

  if( x_coord ) NumDimensions++;
  if( y_coord ) NumDimensions++;
  if( z_coord ) NumDimensions++;

  assert( NumDimensions != 0 );

  // compute the coordinates for each level (that is, the
  // center of gravity, as no mesh is really available for
  // coarser levels)
  ML_Aggregate_Stats_ComputeCoordinates(ml_, agg_,
					x_coord, y_coord, z_coord);
  // few stats
  ML_Aggregate_Stats_Analyze(ml_,agg_);

  // prepare output format. Now it can be:
  // - OpenDX (1D/2D/3D)
  // - XD3D (2D only)

  int Format;
  string FileFormat = List_.get(Prefix+"viz: output format", "xyz");
  // you are a cool guy if you plot with "xyz"
  if( FileFormat == "xyz" ) Format = 1;
  // you are a poor man if you need "dx". God bless us.
  else if( FileFormat == "dx" ) Format = 0;
  else {
    cerr << ErrorMsg_ << "Option `viz: output format' has an incorrect" << endl
      << ErrorMsg_ << "value (" << FileFormat << "). Possible values are" << endl
      << ErrorMsg_ << "<dx> / <xyz>" << endl;
    exit( EXIT_FAILURE );
  }

  bool VizAggre        = List_.get(Prefix+"viz: aggregates", true);
  bool VizPreSmoother  = List_.get(Prefix+"viz: presmoother", false);
  bool VizPostSmoother = List_.get(Prefix+"viz: postsmoother", false);
  bool VizCycle        = List_.get(Prefix+"viz: cycle", false);
  int NumApplSmoother  = List_.get(Prefix+"viz: smoother application", 10);
  int NumCycleSmoother = List_.get(Prefix+"viz: cycle application", 10);
  int ieqn             = List_.get(Prefix+"viz: equation to plot", 0);
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

  // cycle over all levels

  for( int ilevel=0 ; ilevel<NumLevels_ ; ++ilevel ) {

    // =================== //
    // plot the aggregates //
    // =================== //

    if( VizAggre ) 
      ML_Aggregate_Viz(ml_,agg_,Format,NULL,NULL,LevelID_[ilevel]);

    // ============ //
    // pre-smoother //
    // ============ //

    ptr = ((ml_->SingleLevel[ilevel]).pre_smoother);

    if( ptr != NULL && VizPreSmoother ) {

      RandomAndZero(tmp_sol, tmp_rhs,ml_->Amat[ilevel].outvec_leng);
      if( PrintStarting ) 
	ML_Aggregate_Viz(ml_,agg_,Format,tmp_sol,"starting-sol-presmoother",LevelID_[ilevel]);

      // increase the number of applications of the smoother
      int old_ntimes = ptr->ntimes;
      ptr->ntimes = NumApplSmoother;
      ML_Smoother_Apply(ptr, 
			ml_->Amat[ilevel].outvec_leng,
			tmp_sol,
			ml_->Amat[ilevel].outvec_leng,
			tmp_rhs, ML_NONZERO);
      ptr->ntimes = old_ntimes;

      // visualize
      for( int i=0 ; i<NumMyRows() ; i+=NumPDEEqns_ ) 
	plot_me[i/NumPDEEqns_] = tmp_sol[i+ieqn];
      ML_Aggregate_Viz(ml_,agg_,Format,plot_me,"presmoother",LevelID_[ilevel]);
    }

    // ============= //
    // post-smoother //
    // ============= //

    ptr = ((ml_->SingleLevel[ilevel]).post_smoother);
    if( ptr != NULL && VizPostSmoother ) {

      RandomAndZero(tmp_sol, tmp_rhs,ml_->Amat[ilevel].outvec_leng);
      if( PrintStarting ) 
	ML_Aggregate_Viz(ml_,agg_,Format,tmp_sol,"starting-sol-postsmoother",LevelID_[ilevel]);

      // increase the number of applications of the smoother
      int old_ntimes = ptr->ntimes;
      ptr->ntimes = 10;
      ML_Smoother_Apply(ptr, 
			ml_->Amat[ilevel].outvec_leng,
			tmp_sol,
			ml_->Amat[ilevel].outvec_leng,
			tmp_rhs, ML_ZERO);
      ptr->ntimes = old_ntimes;

      // visualize
      for( int i=0 ; i<NumMyRows() ; i+=NumPDEEqns_ ) 
	plot_me[i/NumPDEEqns_] = tmp_sol[i+ieqn];
      ML_Aggregate_Viz(ml_,agg_,Format,plot_me,"postsmoother",LevelID_[ilevel]);
    }

  } // for( ilevel )

  // =============================== //
  // run ML cycle on a random vector //
  // =============================== //

  if( VizCycle ) {

    RandomAndZero(tmp_sol, tmp_rhs,ml_->Amat[0].outvec_leng);
    if( PrintStarting ) 
      ML_Aggregate_Viz(ml_,agg_,Format,tmp_sol,"starting-sol-cycle",LevelID_[0]);

    for( int i=0 ; i<NumCycleSmoother ; ++i ) 
      ML_Cycle_MG(&(ml_->SingleLevel[ml_->ML_finest_level]),
		  tmp_sol, tmp_rhs,
		  ML_NONZERO, ml_->comm, ML_NO_RES_NORM, ml_);

    // visualize the result, this only for level 0
    for( int i=0 ; i<NumMyRows() ; i+=NumPDEEqns_ ) 
      plot_me[i/NumPDEEqns_] = tmp_sol[i+ieqn];
    ML_Aggregate_Viz(ml_,agg_,Format,plot_me,"cycle",LevelID_[0]);

  }

  // =================== //
  // clean up and return //
  // =================== //

  delete [] tmp_sol;
  delete [] tmp_rhs;
  delete [] plot_me;

  ML_Aggregate_Stats_CleanUp_Info(ml_, agg_);
  ML_Aggregate_VizAndStats_Clean( agg_, NumLevels_);

}

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
