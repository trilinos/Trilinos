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

#include <vector>

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
int ML_Epetra::SetDefaults(string ProblemType, ParameterList & List, 
			   int * options, double * params,
			   const string Prefix) 
{
  
  // allocate some memory if the user is not passing the vectors.
  // This is cute, but it may cause memory leaks.

  if( options == NULL ) options = new int[AZ_OPTIONS_SIZE];
  if( params  == NULL ) params  = new double[AZ_PARAMS_SIZE];

  if( ProblemType == "SA" )
    return( ML_Epetra::SetDefaultsSA(List, Prefix, options, params) );
  else if( ProblemType == "maxwell" )
    return( ML_Epetra::SetDefaultsMaxwell(List, Prefix, options, params ) );
  else if( ProblemType == "DD-ML" )
    return( ML_Epetra::SetDefaultsDD_3Levels(List, Prefix, options, params ) );
  else if( ProblemType == "DD-ML-LU" )
    return( ML_Epetra::SetDefaultsDD_3Levels_LU(List, Prefix, options, params ) );
  else if( ProblemType == "DD" )
    return( ML_Epetra::SetDefaultsDD(List, Prefix, options, params ) );
  else if( ProblemType == "DD-LU" )
    return( ML_Epetra::SetDefaultsDD_LU(List, Prefix, options, params ) );
  else {
    cerr << "ERROR: Wrong input parameter in `SetDefaults' ("
	 << ProblemType << "). Should be: " << endl
	 << "ERROR: <SA> / <DD> / <DD-ML> / <maxwell>" << endl;
    exit( EXIT_FAILURE );
  }

  return 0;
  
  
}

// ============================================================================

int ML_Epetra::SetDefaultsDD(ParameterList & List, const string Prefix,
			     int * options, double * params) 
{

  List.set(Prefix+"default values","DD");

  List.set(Prefix+"max levels",2);

  List.set(Prefix+"output",8);
  
  List.set(Prefix+"increasing or decreasing","increasing");

  List.set(Prefix+"PDE equations",1);

  List.set(Prefix+"aggregation: type","METIS");

  List.set(Prefix+"aggregation: local aggregates",1);
  
  List.set(Prefix+"aggregation: damping factor",1.333);

  List.set(Prefix+"coarse: max size",128);

  List.set(Prefix+"aggregation: threshold",0.0);
  
  List.set(Prefix+"smoother: sweeps",2);

  List.set(Prefix+"smoother: damping factor",0.67);

  List.set(Prefix+"smoother: pre or post","both");

  // -- aztec section -- //
  
  List.set(Prefix+"smoother: type","Aztec");

  AZ_defaults(options,params);
  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_scaling] = AZ_none;
  options[AZ_subdomain_solve] = AZ_ilut;
  
  List.set(Prefix+"smoother: Aztec options",options);
    
  List.set(Prefix+"smoother: Aztec params",params);
    
  List.set(Prefix+"smoother: Aztec as solver",false);

  // --- coarse solver --- ///
  
  List.set(Prefix+"coarse: type","Amesos-KLU");

  List.set(Prefix+"prec type","MGV");

  List.set(Prefix+"print unused",-2);

  return 0;

}

// ============================================================================

int ML_Epetra::SetDefaultsDD_LU(ParameterList & List, const string Prefix,
				int * options, double * params) 
{

  List.set(Prefix+"default values","DD-LU");

  List.set(Prefix+"max levels",2);

  List.set(Prefix+"output",8);
  
  List.set(Prefix+"increasing or decreasing","increasing");

  List.set(Prefix+"PDE equations",1);

  List.set(Prefix+"aggregation: type","METIS");

  List.set(Prefix+"aggregation: local aggregates",1);
  
  List.set(Prefix+"aggregation: damping factor",0.01);

  List.set(Prefix+"coarse: max size",128);

  List.set(Prefix+"aggregation: threshold",0.0);
  
  List.set(Prefix+"smoother: sweeps",2);

  List.set(Prefix+"smoother: damping factor",0.67);

  List.set(Prefix+"smoother: pre or post","both");

  // -- aztec section -- //
  
  List.set(Prefix+"smoother: type","Aztec");

  AZ_defaults(options,params);
  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_scaling] = AZ_none;
  options[AZ_subdomain_solve] = AZ_lu;
  
  List.set(Prefix+"smoother: Aztec options",options);
    
  List.set(Prefix+"smoother: Aztec params",params);
    
  List.set(Prefix+"smoother: Aztec as solver",false);

  // --- coarse solver --- ///
  
  List.set(Prefix+"coarse: type","Amesos-KLU");

  List.set(Prefix+"prec type","MGV");

  List.set(Prefix+"print unused",-2);

  return 0;

}

// ============================================================================

int ML_Epetra::SetDefaultsDD_3Levels(ParameterList & List, const string Prefix,
				     int * options, double * params)
{

  List.set(Prefix+"default values","DD-ML");

  List.set(Prefix+"max levels",3);

  List.set(Prefix+"output",8);
  
  List.set(Prefix+"increasing or decreasing","increasing");

  List.set(Prefix+"PDE equations",1);

  List.set(Prefix+"aggregation: type (level 0)","METIS");

  List.set(Prefix+"aggregation: type (level 1)","ParMETIS");

  List.set(Prefix+"aggregation: nodes per aggregate (level 0)",512);

  List.set(Prefix+"aggregation: nodes per aggregate (level 1)",32);
  
  List.set(Prefix+"aggregation: damping factor",4.0/3);

  List.set(Prefix+"coarse: max size",128);

  List.set(Prefix+"aggregation: threshold",0.0);
  
  List.set(Prefix+"smoother: sweeps (level 0)",2);

  List.set(Prefix+"smoother: damping factor (level 0)",0.67);

  List.set(Prefix+"smoother: pre or post (level 0)","both");

  // --- Aztec --- //
  
  List.set(Prefix+"smoother: type","Aztec");
  
  AZ_defaults(options,params);
  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_subdomain_solve] = AZ_ilut;
  options[AZ_overlap] = 0;

  List.set(Prefix+"smoother: Aztec options (level 0)",options);
    
  List.set(Prefix+"smoother: Aztec params (level 0)",params);
    
  List.set(Prefix+"smoother: Aztec as solver (level 0)",false);
  
  // --- coarse --- ///
  
  List.set(Prefix+"coarse: type","Amesos-KLU");

  List.set(Prefix+"prec type","MGV");

  List.set(Prefix+"print unused",-2);
  
  return 0;

}

// ============================================================================

int ML_Epetra::SetDefaultsDD_3Levels_LU(ParameterList & List, const string Prefix,
					int * options, double * params)
{

  List.set(Prefix+"default values","DD-ML-LU");

  List.set(Prefix+"max levels",3);

  List.set(Prefix+"output",10);
  
  List.set(Prefix+"increasing or decreasing","increasing");

  List.set(Prefix+"PDE equations",1);

  List.set(Prefix+"aggregation: type (level 0)","METIS");

  List.set(Prefix+"aggregation: type (level 1)","ParMETIS");

  List.set(Prefix+"aggregation: nodes per aggregate (level 0)",512);

  List.set(Prefix+"aggregation: nodes per aggregate (level 1)",512);
  
  List.set(Prefix+"aggregation: damping factor",4.0/3);

  List.set(Prefix+"coarse: max size",128);

  List.set(Prefix+"aggregation: threshold",0.0);
  
  List.set(Prefix+"smoother: sweeps (level 0)",2);

  List.set(Prefix+"smoother: damping factor (level 0)",0.67);

  List.set(Prefix+"smoother: pre or post (level 0)","both");

  // --- Aztec --- //
  
  List.set(Prefix+"smoother: type","Aztec");

  AZ_defaults(options,params);
  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_subdomain_solve] = AZ_lu;
  options[AZ_overlap] = 0;

  List.set(Prefix+"smoother: Aztec options (level 0)",options);
    
  List.set(Prefix+"smoother: Aztec params (level 0)",params);
    
  List.set(Prefix+"smoother: Aztec as solver (level 0)",false);
  
  // --- coarse --- ///
  
  List.set(Prefix+"coarse: type","Amesos-KLU");

  List.set(Prefix+"prec type","MGV");

  List.set(Prefix+"print unused",-2);
  
  return 0;

}

// ============================================================================

int ML_Epetra::SetDefaultsMaxwell(ParameterList & List, const string Prefix,
				  int * options, double * params)
{

  List.set(Prefix+"default values","maxwell");

  List.set(Prefix+"max levels",10);

  List.set(Prefix+"output",10);
  
  List.set(Prefix+"PDE equations",1);

  List.set(Prefix+"increasing or decreasing","decreasing");

  // aggregation: Uncoupled for first levels, then MIS
  List.set(Prefix+"aggregation: type","Uncoupled-MIS");

  // optimal value for smoothed aggregation
  List.set(Prefix+"aggregation: damping factor",1.3333);

  // relative small coarse size
  List.set(Prefix+"coarse: max size",16);

  // don't forget any element
  List.set(Prefix+"aggregation: threshold",0.0);

  // gauss-seidel for all levels
  List.set(Prefix+"smoother: sweeps",2);

  List.set(Prefix+"smoother: damping factor",0.67);

  List.set(Prefix+"smoother: type","MLS");

  List.set(Prefix+"smoother: MLS polynomial order",3);
  
  List.set(Prefix+"smoother: pre or post","both");
  
  // simplest solver on coarse problem
  List.set(Prefix+"coarse: type","SuperLU");
//  Tim Davis' simple serial LU package.  It's part of Amesos
//  itself.
//  List.set(Prefix,"Amesos-KLU");

  List.set(Prefix+"prec type","MGV");

  // print unused Prefixs on proc 0
  List.set(Prefix+"print unused",-2);

  return 0;
  
}

// ============================================================================

int ML_Epetra::SetDefaultsSA(ParameterList & List, const string Prefix,
			     int * options, double * params)
{

  List.set(Prefix+"default values","SA");

  List.set(Prefix+"max levels",10);

  List.set(Prefix+"output",8);
  
  List.set(Prefix+"PDE equations",1);

  List.set(Prefix+"increasing or decreasing","increasing");

  // aggregation: Uncoupled for first levels, then MIS
  List.set(Prefix+"aggregation: type","Uncoupled-MIS");
  
  // optimal value for smoothed aggregation
  List.set(Prefix+"aggregation: damping factor",1.3333);

  // relative small coarse size
  List.set(Prefix+"coarse: max size",16);

  // don't forget any element
  List.set(Prefix+"aggregation: threshold",0.0);

  // guass-seidel for all levels
  List.set(Prefix+"smoother: sweeps",2);

  List.set(Prefix+"smoother: damping factor",0.67);

  List.set(Prefix+"smoother: type","Gauss-Seidel");
  
  List.set(Prefix+"smoother: pre or post","both");
  
  // simplest solver on coarse problem
  List.set(Prefix+"coarse: type","Amesos-KLU");

  List.set(Prefix+"prec type","MGV");

  // print unused Prefixs on proc 0
  List.set(Prefix+"print unused",-2);
  
  return 0;

}

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
