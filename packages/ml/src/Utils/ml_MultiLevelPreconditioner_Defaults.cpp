/*!
 *  \file ml_MultiLevelPreconditioner_Defaults.cpp
 *
 *  \brief ML black-box preconditioner for Epetra_RowMatrix derived classes.
 *
 *  \author Marzio Sala, SNL, 9214
 *
 *  \date Last update do Doxygen: 22-Jul-04
 *
 */

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
/*! Sets the following default values for "DD":
  - "default values" = DD"
  - \c "max levels" = 2
  - \c "output" = 8
  - \c "increasing or decreasing" = increasing"
  - \c "PDE equations" = 1
  - \c "aggregation: type" = METIS"
  - \c "aggregation: local aggregates" = 1
  - \c "aggregation: damping factor" = 1.333
  - \c "coarse: max size" = 128
  - \c "aggregation: threshold" = 0.0
  - \c "smoother: sweeps" = 2
  - \c "smoother: damping factor" = 0.67
  - \c "smoother: pre or post" = both"
  - \c "smoother: type" = Aztec". The Aztec smoother has the following features:
    - \c  options[AZ_precond] = \c AZ_dom_decomp
    - \c options[AZ_scaling] = \c AZ_none
    - \c options[AZ_subdomain_solve] = \c AZ_ilut
  - \c "smoother: Aztec options" = options
  - \c "smoother: Aztec params" = params
  - \c "smoother: Aztec as solver" = false
  - \c "coarse: type" = Amesos-KLU"
  - \c "prec type" = MGV"
  - \c "print unused" = 1
 */
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

  List.set(Prefix+"print unused",-1);

  return 0;

}

// ============================================================================
/*! As for SetDefaultsDD(), but used exact LU decompositions on each subdomains.
 */
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

  List.set(Prefix+"print unused",-1);

  return 0;

}

// ============================================================================
/*! Sets the following default values for "DD"
 - \c "default values" = "DD-ML"
 - \c "max levels" = 3
 - \c "output" = 8
 - \c "increasing or decreasing" = "increasing"
 - \c "PDE equations" = 1
 - \c "aggregation: type (level 0)" = "METIS"
 - \c "aggregation: type (level 1)" = "ParMETIS"
 - \c "aggregation: nodes per aggregate (level 0)" = 512
 - \c "aggregation: nodes per aggregate (level 1)" = 32
 - \c "aggregation: damping factor" = 4.0/3
 - \c "coarse: max size" = 128
 - \c "aggregation: threshold" = 0.0
 - \c "smoother: sweeps (level 0)" = 2
 - \c "smoother: damping factor (level 0)" = 0.67
 - \c "smoother: pre or post (level 0)" = "both"
 - \c "smoother: type" = "Aztec".  The Aztec smoother has the following features:
     - \c options[AZ_precond] = \c AZ_dom_decomp
     - \c options[AZ_scaling] = \c AZ_none
     - \c options[AZ_subdomain_solve] = \c AZ_ilut
     - \c options[AZ_overlap] = 0
 - \c "smoother: Aztec options (level 0)" = options
 - \c "smoother: Aztec params (level 0)" = params
 - \c "smoother: Aztec as solver (level 0)" = false
 - \c "coarse: type" = "Amesos-KLU"
 - \c "prec type" = "MGV"
 - \c "print unused" = -1
 */
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

  List.set(Prefix+"print unused",-1);
  
  return 0;

}

// ============================================================================
/*! As for SetDefaultsDD_3Levels, but with LU factorizations on each subdomain
 * for each level
 */
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

  List.set(Prefix+"print unused",-1);
  
  return 0;

}

// ============================================================================
/*! Set values for Maxwell:
 * - \c "default values" = "maxwell"
 * - \c "max levels" = 10
 * - \c "output" = 10
 * - \c "PDE equations" = 1
 * - \c "increasing or decreasing" = "decreasing"
 * - \c "aggregation: type" = "Uncoupled-MIS"
 * - \c "aggregation: damping factor" = 1.3333
 * - \c "coarse: max size" = 16
 * - \c "aggregation: threshold" = 0.0
 * - \c "smoother: sweeps" = 2
 * - \c "smoother: damping factor" = 0.67
 * - \c "smoother: type" = "MLS"
 * - \c "smoother: MLS polynomial order" = 3
 * - \c "smoother: pre or post" = "both"
 * - \c "coarse: type" = "SuperLU"
 * - \c "prec type" = "MGV"
 * - \c "print unused" = -1
 */
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
  List.set(Prefix+"print unused",-1);

  return 0;
  
}

// ============================================================================
/*! Set default values for classical smoothed aggregation:
 * - \c "default values" = "SA"
 * - \c "max levels" = 10
 * - \c "output" = 8
 * - \c "PDE equations" = 1
 * - \c "increasing or decreasing" = "increasing"
 * - \c "aggregation: type" = "Uncoupled-MIS"
 * - \c "aggregation: damping factor" = 1.3333
 * - \c "coarse: max size" = 16
 * - \c "aggregation: threshold" = 0.0
 * - \c "smoother: sweeps" = 2
 * - \c "smoother: damping factor" = 0.67
 * - \c "smoother: type" = "Gauss-Seidel
 * - \c "smoother: pre or post" = "both"
 * - \c "coarse: type" = "Amesos-KLU"
 * - \c "prec type" = "MGV"
 * - \c "print unused" =-1
 */
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

  // gauss-seidel for all levels
  List.set(Prefix+"smoother: sweeps",2);

  List.set(Prefix+"smoother: damping factor",0.67);

  List.set(Prefix+"smoother: type","Gauss-Seidel");
  
  List.set(Prefix+"smoother: pre or post","both");
  
  // simplest solver on coarse problem
  List.set(Prefix+"coarse: type","Amesos-KLU");

  List.set(Prefix+"prec type","MGV");

  // print unused Prefixs on proc 0
  List.set(Prefix+"print unused",-1);
  
  return 0;

}

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
