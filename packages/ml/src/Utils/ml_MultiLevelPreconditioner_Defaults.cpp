/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/************************************************************************/
/*          Utilities for Trilinos/ML users                             */
/*----------------------------------------------------------------------*/
/* Authors : Mike Heroux (SNL)                                          */
/*           Jonathan Hu  (SNL)                                         */
/*           Ray Tuminaro (SNL)                                         */
/*           Marzio Sala (SNL)                                          */
/************************************************************************/


#include "ml_common.h"
#include "ml_common.h"
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
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

//#include <cstring>
#include "ml_amesos_wrap.h"
#include "ml_ifpack_wrap.h"
#include "ml_agg_METIS.h"
#include "ml_epetra_utils.h"

#include "ml_epetra_preconditioner.h"
#include "ml_agg_ParMETIS.h"

#include "ml_anasazi.h"

#include "ml_ggb.h"

using namespace Teuchos;
using namespace ML_Epetra;

// ============================================================================

/*! Currently supported defaults:

  - SA: classic smoothed aggregation, using Uncoupled/MIS schemes;

  - maxwell: for Maxwell equations (under development);

  - DD: 2-level domain decomposition preconditioner, using METIS, and
  Aztec smoothers with ILUT(0);

  - DD-ML: 3-level domain decomposition preconditioner, using METIS and
  ParMETIS, and Aztec smoothers with ILUT(0);

  - DD-LU:  2-level domain decomposition preconditioner, using METIS, and
  Aztec smoothers with LU.

*/
int ML_Epetra::SetDefaults(string ProblemType, ParameterList & List, char * Prefix_ )
{
  
  int rv = 0;
  
  if( ProblemType == "SA" )
    return( ML_Epetra::SetDefaultsSA(List, Prefix_) );
  else if( ProblemType == "maxwell" )
    return( ML_Epetra::SetDefaultsMaxwell(List, Prefix_ ) );
  else if( ProblemType == "DD-ML" )
    return( ML_Epetra::SetDefaultsDD_3Levels(List, Prefix_ ) );
  else if( ProblemType == "DD-ML-LU" )
    return( ML_Epetra::SetDefaultsDD_3Levels_LU(List, Prefix_ ) );
  else if( ProblemType == "DD" )
    return( ML_Epetra::SetDefaultsDD(List, Prefix_ ) );
  else if( ProblemType == "DD-LU" )
    return( ML_Epetra::SetDefaultsDD_LU(List, Prefix_ ) );
  else {
    cerr << "ERROR: Wrong input parameter in `SetDefaults' ("
	 << ProblemType << "). Should be: " << endl
	 << "ERROR: <SA> / <DD> / <DD-ML> / <maxwell>" << endl;
    exit( EXIT_FAILURE );
    rv = 1;
  }

  EPETRA_CHK_ERR(rv);

  return rv;
  
  
}

// ============================================================================

int ML_Epetra::SetDefaultsDD(ParameterList & List, char * Prefix) 
{

  char parameter[80];

  sprintf(parameter,"%sdefault values", Prefix);
  List.set(parameter,"DD");

  sprintf(parameter,"%smax levels", Prefix);
  List.set(parameter,2);

  sprintf(parameter,"%soutput", Prefix);
  List.set(parameter,8);
  
  sprintf(parameter,"%sincreasing or decreasing", Prefix);
  List.set(parameter,"increasing");

  sprintf(parameter,"%sPDE equations", Prefix);
  List.set(parameter,1);

  sprintf(parameter,"%saggregation: type",Prefix);
  List.set(parameter,"METIS");

  sprintf(parameter,"%saggregation: local aggregates",Prefix);
  List.set(parameter,1);
  
  sprintf(parameter,"%saggregation: damping factor",Prefix);
  List.set(parameter,1.333);

  sprintf(parameter,"%scoarse: max size",Prefix);
  List.set(parameter,128);

  sprintf(parameter,"%saggregation: threshold",Prefix);
  List.set(parameter,0.0);
  
  sprintf(parameter,"%ssmoother: sweeps",Prefix);
  List.set(parameter,2);

  sprintf(parameter,"%ssmoother: damping factor",Prefix);
  List.set(parameter,0.67);

  sprintf(parameter,"%ssmoother: pre or post",Prefix);
  List.set(parameter,"both");

  // -- aztec section -- //
  
  sprintf(parameter,"%ssmoother: type",Prefix);
  List.set(parameter,"Aztec");

  int * SmootherOptionsList = new int[AZ_OPTIONS_SIZE];
  double * SmootherParamsList = new double[AZ_PARAMS_SIZE];
  
  AZ_defaults(SmootherOptionsList,SmootherParamsList);
  SmootherOptionsList[AZ_precond] = AZ_dom_decomp;
  SmootherOptionsList[AZ_scaling] = AZ_none;
  SmootherOptionsList[AZ_subdomain_solve] = AZ_ilut;
  
  sprintf(parameter,"%ssmoother: Aztec options",Prefix);
  List.set(parameter,SmootherOptionsList);
    
  sprintf(parameter,"%ssmoother: Aztec params",Prefix);
  List.set(parameter,SmootherParamsList);
    
  sprintf(parameter,"%ssmoother: Aztec as solver",Prefix);
  List.set(parameter,false);

  // --- coarse solver --- ///
  
  sprintf(parameter,"%scoarse: type",Prefix);
  List.set(parameter,"Amesos-KLU");

  sprintf(parameter,"%sprec type",Prefix);
  List.set(parameter,"MGV");

  sprintf(parameter,"%sprint unused",Prefix);
  List.set(parameter,-2);

  return 0;

}

// ============================================================================

int ML_Epetra::SetDefaultsDD_LU(ParameterList & List, char * Prefix) 
{

  char parameter[80];

  sprintf(parameter,"%sdefault values", Prefix);
  List.set(parameter,"DD-LU");

  sprintf(parameter,"%smax levels", Prefix);
  List.set(parameter,2);

  sprintf(parameter,"%soutput", Prefix);
  List.set(parameter,8);
  
  sprintf(parameter,"%sincreasing or decreasing", Prefix);
  List.set(parameter,"increasing");

  sprintf(parameter,"%sPDE equations", Prefix);
  List.set(parameter,1);

  sprintf(parameter,"%saggregation: type",Prefix);
  List.set(parameter,"METIS");

  sprintf(parameter,"%saggregation: local aggregates",Prefix);
  List.set(parameter,1);
  
  sprintf(parameter,"%saggregation: damping factor",Prefix);
  List.set(parameter,0.01);

  sprintf(parameter,"%scoarse: max size",Prefix);
  List.set(parameter,128);

  sprintf(parameter,"%saggregation: threshold",Prefix);
  List.set(parameter,0.0);
  
  sprintf(parameter,"%ssmoother: sweeps",Prefix);
  List.set(parameter,2);

  sprintf(parameter,"%ssmoother: damping factor",Prefix);
  List.set(parameter,0.67);

  sprintf(parameter,"%ssmoother: pre or post",Prefix);
  List.set(parameter,"both");

  // -- aztec section -- //
  
  sprintf(parameter,"%ssmoother: type",Prefix);
  List.set(parameter,"Aztec");

  int * SmootherOptionsList = new int[AZ_OPTIONS_SIZE];
  double * SmootherParamsList = new double[AZ_PARAMS_SIZE];
  
  AZ_defaults(SmootherOptionsList,SmootherParamsList);
  SmootherOptionsList[AZ_precond] = AZ_dom_decomp;
  SmootherOptionsList[AZ_scaling] = AZ_none;
  SmootherOptionsList[AZ_subdomain_solve] = AZ_lu;
  
  sprintf(parameter,"%ssmoother: Aztec options",Prefix);
  List.set(parameter,SmootherOptionsList);
    
  sprintf(parameter,"%ssmoother: Aztec params",Prefix);
  List.set(parameter,SmootherParamsList);
    
  sprintf(parameter,"%ssmoother: Aztec as solver",Prefix);
  List.set(parameter,false);

  // --- coarse solver --- ///
  
  sprintf(parameter,"%scoarse: type",Prefix);
  List.set(parameter,"Amesos-KLU");

  sprintf(parameter,"%sprec type",Prefix);
  List.set(parameter,"MGV");

  sprintf(parameter,"%sprint unused",Prefix);
  List.set(parameter,-2);

  return 0;

}

// ============================================================================

int ML_Epetra::SetDefaultsDD_3Levels(ParameterList & List, char * Prefix_) 
{

  char parameter[80];

  sprintf(parameter,"%sdefault values", Prefix_);
  List.set(parameter,"DD-ML");

  sprintf(parameter,"%smax levels", Prefix_);
  List.set(parameter,3);

  sprintf(parameter,"%soutput", Prefix_);
  List.set(parameter,8);
  
  sprintf(parameter,"%sincreasing or decreasing", Prefix_);
  List.set(parameter,"increasing");

  sprintf(parameter,"%sPDE equations", Prefix_);
  List.set(parameter,1);

  sprintf(parameter,"%saggregation: type (level 0)",Prefix_);
  List.set(parameter,"METIS");

  sprintf(parameter,"%saggregation: type (level 1)",Prefix_);
  List.set(parameter,"ParMETIS");

  sprintf(parameter,"%saggregation: nodes per aggregate (level 0)",Prefix_);
  List.set(parameter,512);

  sprintf(parameter,"%saggregation: nodes per aggregate (level 1)",Prefix_);
  List.set(parameter,32);
  
  sprintf(parameter,"%saggregation: damping factor",Prefix_);
  List.set(parameter,4.0/3);

  sprintf(parameter,"%scoarse: max size",Prefix_);
  List.set(parameter,128);

  sprintf(parameter,"%saggregation: threshold",Prefix_);
  List.set(parameter,0.0);
  
  sprintf(parameter,"%ssmoother: sweeps (level 0)",Prefix_);
  List.set(parameter,2);

  sprintf(parameter,"%ssmoother: damping factor (level 0)",Prefix_);
  List.set(parameter,0.67);

  sprintf(parameter,"%ssmoother: pre or post (level 0)",Prefix_);
  List.set(parameter,"both");

  // --- Aztec --- //
  
  sprintf(parameter,"%ssmoother: type",Prefix_);
  List.set(parameter,"Aztec");
  
  int * SmootherOptionsList = new int[AZ_OPTIONS_SIZE];
  double * SmootherParamsList = new double[AZ_PARAMS_SIZE];

  AZ_defaults(SmootherOptionsList,SmootherParamsList);
  SmootherOptionsList[AZ_precond] = AZ_dom_decomp;
  SmootherOptionsList[AZ_subdomain_solve] = AZ_ilut;
  SmootherOptionsList[AZ_overlap] = 0;

  sprintf(parameter,"%ssmoother: Aztec options (level 0)",Prefix_);
  List.set(parameter,SmootherOptionsList);
    
  sprintf(parameter,"%ssmoother: Aztec params (level 0)",Prefix_);
  List.set(parameter,SmootherParamsList);
    
  sprintf(parameter,"%ssmoother: Aztec as solver (level 0)",Prefix_);
  List.set(parameter,false);
  
  // --- coarse --- ///
  
  sprintf(parameter,"%scoarse: type",Prefix_);
  List.set(parameter,"Amesos-KLU");

  sprintf(parameter,"%sprec type",Prefix_);
  List.set(parameter,"MGV");

  sprintf(parameter,"%sprint unused",Prefix_);
  List.set(parameter,-2);
  
  return 0;

}

// ============================================================================

int ML_Epetra::SetDefaultsDD_3Levels_LU(ParameterList & List, char * Prefix_) 
{

  char parameter[80];

  sprintf(parameter,"%sdefault values", Prefix_);
  List.set(parameter,"DD-ML-LU");

  sprintf(parameter,"%smax levels", Prefix_);
  List.set(parameter,3);

  sprintf(parameter,"%soutput", Prefix_);
  List.set(parameter,10);
  
  sprintf(parameter,"%sincreasing or decreasing", Prefix_);
  List.set(parameter,"increasing");

  sprintf(parameter,"%sPDE equations", Prefix_);
  List.set(parameter,1);

  sprintf(parameter,"%saggregation: type (level 0)",Prefix_);
  List.set(parameter,"METIS");

  sprintf(parameter,"%saggregation: type (level 1)",Prefix_);
  List.set(parameter,"ParMETIS");

  sprintf(parameter,"%saggregation: nodes per aggregate (level 0)",Prefix_);
  List.set(parameter,512);

  sprintf(parameter,"%saggregation: nodes per aggregate (level 1)",Prefix_);
  List.set(parameter,512);
  
  sprintf(parameter,"%saggregation: damping factor",Prefix_);
  List.set(parameter,4.0/3);

  sprintf(parameter,"%scoarse: max size",Prefix_);
  List.set(parameter,128);

  sprintf(parameter,"%saggregation: threshold",Prefix_);
  List.set(parameter,0.0);
  
  sprintf(parameter,"%ssmoother: sweeps (level 0)",Prefix_);
  List.set(parameter,2);

  sprintf(parameter,"%ssmoother: damping factor (level 0)",Prefix_);
  List.set(parameter,0.67);

  sprintf(parameter,"%ssmoother: pre or post (level 0)",Prefix_);
  List.set(parameter,"both");

  // --- Aztec --- //
  
  sprintf(parameter,"%ssmoother: type",Prefix_);
  List.set(parameter,"Aztec");
  
  int * SmootherOptionsList = new int[AZ_OPTIONS_SIZE];
  double * SmootherParamsList = new double[AZ_PARAMS_SIZE];

  AZ_defaults(SmootherOptionsList,SmootherParamsList);
  SmootherOptionsList[AZ_precond] = AZ_dom_decomp;
  SmootherOptionsList[AZ_subdomain_solve] = AZ_lu;
  SmootherOptionsList[AZ_overlap] = 0;

  sprintf(parameter,"%ssmoother: Aztec options (level 0)",Prefix_);
  List.set(parameter,SmootherOptionsList);
    
  sprintf(parameter,"%ssmoother: Aztec params (level 0)",Prefix_);
  List.set(parameter,SmootherParamsList);
    
  sprintf(parameter,"%ssmoother: Aztec as solver (level 0)",Prefix_);
  List.set(parameter,false);
  
  // --- coarse --- ///
  
  sprintf(parameter,"%scoarse: type",Prefix_);
  List.set(parameter,"Amesos-KLU");

  sprintf(parameter,"%sprec type",Prefix_);
  List.set(parameter,"MGV");

  sprintf(parameter,"%sprint unused",Prefix_);
  List.set(parameter,-2);
  
  return 0;

}

// ============================================================================

int ML_Epetra::SetDefaultsMaxwell(ParameterList & List, char * Prefix_) 
{

  // FIXME : here default values for Maxwell
  
  char parameter[80];
  int MaxLevels = 10;
  
  sprintf(parameter,"%sdefault values", Prefix_);
  List.set(parameter,"maxwell");

  sprintf(parameter,"%smax levels", Prefix_);
  List.set(parameter,MaxLevels);

  sprintf(parameter,"%soutput", Prefix_);
  List.set(parameter,10);
  
  sprintf(parameter,"%sPDE equations", Prefix_);
  List.set(parameter,1);

  sprintf(parameter,"%sincreasing or decreasing", Prefix_);
  List.set(parameter,"decreasing");

  // aggregation: Uncoupled for first levels, then MIS
  sprintf(parameter,"%saggregation: type",Prefix_);
  List.set(parameter,"Uncoupled-MIS");

  // optimal value for smoothed aggregation
  sprintf(parameter,"%saggregation: damping factor",Prefix_);
  List.set(parameter,1.3333);

  // relative small coarse size
  sprintf(parameter,"%scoarse: max size",Prefix_);
  List.set(parameter,16);

  // don't forget any element
  sprintf(parameter,"%saggregation: threshold",Prefix_);
  List.set(parameter,0.0);

  // gauss-seidel for all levels
  sprintf(parameter,"%ssmoother: sweeps",Prefix_);
  List.set(parameter,2);

  sprintf(parameter,"%ssmoother: damping factor",Prefix_);
  List.set(parameter,0.67);

  sprintf(parameter,"%ssmoother: type",Prefix_);
  List.set(parameter,"MLS");

  sprintf(parameter,"%ssmoother: MLS polynomial order", Prefix_);
  List.set(parameter,3);
  
  sprintf(parameter,"%ssmoother: pre or post",Prefix_);
  List.set(parameter,"both");
  
  // simplest solver on coarse problem
  sprintf(parameter,"%scoarse: type",Prefix_);
  List.set(parameter,"SuperLU");
//  Tim Davis' simple serial LU package.  It's part of Amesos
//  itself.
//  List.set(parameter,"Amesos-KLU");

  sprintf(parameter,"%sprec type",Prefix_);
  List.set(parameter,"MGV");

  // print unused parameters on proc 0
  sprintf(parameter,"%sprint unused",Prefix_);
  List.set(parameter,-2);

  return 0;
  
}

// ============================================================================

int ML_Epetra::SetDefaultsSA(ParameterList & List, char * Prefix_) 
{

  char parameter[80];
  int MaxLevels = 16;
  
  sprintf(parameter,"%sdefault values", Prefix_);
  List.set(parameter,"SA");

  sprintf(parameter,"%smax levels", Prefix_);
  List.set(parameter,MaxLevels);

  sprintf(parameter,"%soutput", Prefix_);
  List.set(parameter,8);
  
  sprintf(parameter,"%sPDE equations", Prefix_);
  List.set(parameter,1);

  sprintf(parameter,"%sincreasing or decreasing", Prefix_);
  List.set(parameter,"increasing");

  // aggregation: Uncoupled for first levels, then MIS
  sprintf(parameter,"%saggregation: type",Prefix_);
  List.set(parameter,"Uncoupled-MIS");
  
  // optimal value for smoothed aggregation
  sprintf(parameter,"%saggregation: damping factor",Prefix_);
  List.set(parameter,1.3333);

  // relative small coarse size
  sprintf(parameter,"%scoarse: max size",Prefix_);
  List.set(parameter,16);

  // don't forget any element
  sprintf(parameter,"%saggregation: threshold",Prefix_);
  List.set(parameter,0.0);

  // guass-seidel for all levels
  sprintf(parameter,"%ssmoother: sweeps",Prefix_);
  List.set(parameter,2);

  sprintf(parameter,"%ssmoother: damping factor",Prefix_);
  List.set(parameter,0.67);

  sprintf(parameter,"%ssmoother: type",Prefix_);
  List.set(parameter,"Gauss-Seidel");
  
  sprintf(parameter,"%ssmoother: pre or post",Prefix_);
  List.set(parameter,"both");
  
  // simplest solver on coarse problem
  sprintf(parameter,"%scoarse: type",Prefix_);
  List.set(parameter,"Amesos-KLU");

  sprintf(parameter,"%sprec type",Prefix_);
  List.set(parameter,"MGV");

  // print unused parameters on proc 0
  sprintf(parameter,"%sprint unused",Prefix_);
  List.set(parameter,-2);
  
  return 0;

}

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
