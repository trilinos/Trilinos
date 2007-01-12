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

#include "ml_MultiLevelPreconditioner.h"
#include "ml_agg_ParMETIS.h"

#include "ml_epetra.h"
#include "ml_anasazi.h"
#include "Teuchos_RefCountPtr.hpp"

using namespace Teuchos;
using namespace ML_Epetra;

// ============================================================================

int ML_Epetra::SetDefaults(string ProblemType, ParameterList & List, 
			   int * options, double * params, bool OverWrite)
{
  
  // allocate some memory if the user is not passing the vectors.
  // This is cute, but it may cause memory leaks.

#ifdef HAVE_ML_AZTECOO
  bool SetDefaults = false;
  if (options == NULL || params == NULL)
    SetDefaults = true;

  if (options == NULL) options = new int[AZ_OPTIONS_SIZE];
  if (params  == NULL) params  = new double[AZ_PARAMS_SIZE];
  if (SetDefaults)
    AZ_defaults(options,params);
#endif

  if( ProblemType == "SA" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsSA(List, options, params, OverWrite) );
  }
  else if( ProblemType == "NSSA" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsNSSA(List, options, params, OverWrite) );
  }
  else if( ProblemType == "maxwell" || ProblemType == "Maxwell" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsMaxwell(List, options, params,OverWrite));
  }
  else if( ProblemType == "DD-ML" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsDD_3Levels(List, options,
                params,OverWrite) );
  }
  else if( ProblemType == "DD-ML-LU" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsDD_3Levels_LU(List, options, params,
                OverWrite ) );
  }
  else if( ProblemType == "DD" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsDD(List, options, params, OverWrite ) );
  }
  else if( ProblemType == "DD-LU" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsDD_LU(List, options, params, OverWrite ));
  }
  else {
    cerr << "ERROR: Wrong input parameter in `SetDefaults' ("
	 << ProblemType << "). Should be: " << endl
	 << "ERROR: <SA> / <DD> / <DD-ML> / <maxwell>" << endl;
    ML_CHK_ERR(-1); // bad input option
  }

  return(0);
  
  
}

// ============================================================================
/*! Sets the following default values for "DD":
  - "default values" = DD"
  - \c "max levels" = 2
  - \c "output" = 0
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
  - \c "smoother: type" = Aztec". 
  - \c "smoother: Aztec options" = options
  - \c "smoother: Aztec params" = params
  - \c "smoother: Aztec as solver" = false
  - \c "coarse: type" = Amesos-KLU"
  - \c "prec type" = MGV"
  - \c "print unused" = -2
 */
int ML_Epetra::SetDefaultsDD(ParameterList & inList, 
			     int * options, double * params, bool OverWrite) 
{
  ParameterList List;

  List.set("default values","DD");

  List.set("max levels",2);

  List.set("output",0);
  
  List.set("increasing or decreasing","increasing");

  List.set("PDE equations",1);

  List.set("aggregation: type","METIS");

  List.set("aggregation: local aggregates",128);
  
  List.set("aggregation: damping factor",1.333);

  List.set("coarse: max size",128);

  List.set("aggregation: threshold",0.0);
  
  List.set("smoother: sweeps",2);

  List.set("smoother: damping factor",0.67);

  List.set("smoother: pre or post","both");

  // -- aztec section -- //
  
  List.set("smoother: type","Aztec");

#ifdef HAVE_ML_AZTECOO
  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_subdomain_solve] = AZ_ilu;
  List.set("smoother: Aztec options",options);
    
  List.set("smoother: Aztec params",params);
    
  List.set("smoother: Aztec as solver",false);
#endif

  // --- coarse solver --- ///
  
  List.set("coarse: type","Amesos-KLU");

  List.set("prec type","MGV");

  List.set("print unused",-2);

  for (ParameterList::ConstIterator param=List.begin(); param!=List.end(); param++)
    if ( inList.isParameter(inList.name(param)) == false || OverWrite )
      inList.setEntry( inList.name(param) , inList.entry(param) );

  return 0;

} //ML_Epetra::SetDefaultsDD()

// ============================================================================
/*! As for SetDefaultsDD(), but used exact LU decompositions on each subdomains.
 */
int ML_Epetra::SetDefaultsDD_LU(ParameterList & inList, 
				int * options, double * params, bool OverWrite) 
{
  ParameterList List;

  List.set("default values","DD-LU");

  List.set("max levels",2);

  List.set("output",0);
  
  List.set("increasing or decreasing","increasing");

  List.set("PDE equations",1);

  List.set("aggregation: type","METIS");

  List.set("aggregation: local aggregates",1);
  
  List.set("aggregation: damping factor",0.01);

  List.set("coarse: max size",128);

  List.set("aggregation: threshold",0.0);
  
  List.set("smoother: sweeps",2);

  List.set("smoother: damping factor",0.67);

  List.set("smoother: pre or post","both");

  // -- aztec section -- //
  
  List.set("smoother: type","Aztec");

#ifdef HAVE_ML_AZTECOO
  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_subdomain_solve] = AZ_lu;

  List.set("smoother: Aztec options",options);
    
  List.set("smoother: Aztec params",params);
    
  List.set("smoother: Aztec as solver",false);
#endif

  // --- coarse solver --- ///
  
  List.set("coarse: type","Amesos-KLU");

  List.set("prec type","MGV");

  List.set("print unused",-2);

  for(ParameterList::ConstIterator param=List.begin(); param!=List.end(); param++)
    if ( inList.isParameter(inList.name(param)) == false || OverWrite )
      inList.setEntry( inList.name(param) , inList.entry(param) );

  return 0;

} //ML_Epetra::SetDefaultsDD_LU()

// ============================================================================
/*! Sets the following default values for "DD"
 - \c "default values" = "DD-ML"
 - \c "max levels" = 3
 - \c "output" = 0
 - \c "increasing or decreasing" = "increasing"
 - \c "PDE equations" = 1
 - \c "aggregation: type" = "METIS"
 - \c "aggregation: type" = "METIS"
 - \c "aggregation: nodes per aggregate" = 64
 - \c "aggregation: damping factor" = 4.0/3
 - \c "coarse: max size" = 128
 - \c "aggregation: threshold" = 0.0
 - \c "smoother: sweeps" = 2
 - \c "smoother: damping factor" = 0.67
 - \c "smoother: pre or post" = "both"
 - \c "smoother: type" = "Aztec".  
 - \c "smoother: Aztec options" = options
 - \c "smoother: Aztec params" = params
 - \c "smoother: Aztec as solver" = false
 - \c "coarse: type" = "Amesos-KLU"
 - \c "prec type" = "MGV"
 - \c "print unused" = -2
 */
int ML_Epetra::SetDefaultsDD_3Levels(ParameterList & inList, 
				     int * options, double * params, bool OverWrite)
{
  ParameterList List;

  List.set("default values","DD-ML");

  List.set("max levels",3);

  List.set("output",0);
  
  List.set("increasing or decreasing","increasing");

  List.set("PDE equations",1);

  List.set("aggregation: type","METIS");

  List.set("aggregation: type","METIS");

  List.set("aggregation: nodes per aggregate",64);

  List.set("aggregation: damping factor",4.0/3);

  List.set("coarse: max size",128);

  List.set("aggregation: threshold",0.0);
  
  List.set("smoother: sweeps",2);

  List.set("smoother: damping factor",0.67);

  List.set("smoother: pre or post","both");

  // --- Aztec --- //
  
  List.set("smoother: type","Aztec");
  
#ifdef HAVE_ML_AZTECOO
  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_subdomain_solve] = AZ_ilu;

  List.set("smoother: Aztec options",options);
    
  List.set("smoother: Aztec params",params);
    
  List.set("smoother: Aztec as solver",false);
#endif
  
  // --- coarse --- ///
  
  List.set("coarse: type","Amesos-KLU");

  List.set("prec type","MGV");

  List.set("print unused",-2);

  for (ParameterList::ConstIterator param=List.begin(); param!=List.end(); param++)
    if ( inList.isParameter(inList.name(param)) == false || OverWrite )
      inList.setEntry( inList.name(param) , inList.entry(param) );
  
  return 0;

} //ML_Epetra::SetDefaultsDD_3Levels()

// ============================================================================
/*! As for SetDefaultsDD_3Levels, but with LU factorizations on each subdomain
 * for each level
 */
int ML_Epetra::SetDefaultsDD_3Levels_LU(ParameterList & inList, 
					int * options, double * params, bool OverWrite)
{
  ParameterList List;

  List.set("default values","DD-ML-LU");

  List.set("max levels",3);

  List.set("output",0);
  
  List.set("increasing or decreasing","increasing");

  List.set("PDE equations",1);

  List.set("aggregation: type","METIS");

  List.set("aggregation: nodes per aggregate",64);
  
  List.set("aggregation: damping factor",4.0/3);

  List.set("coarse: max size",128);

  List.set("aggregation: threshold",0.0);
  
  List.set("smoother: sweeps",2);

  List.set("smoother: damping factor",0.67);

  List.set("smoother: pre or post","both");

  // --- Aztec --- //
  
  List.set("smoother: type","Aztec");

#ifdef HAVE_ML_AZTECOO
  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_subdomain_solve] = AZ_lu;

  List.set("smoother: Aztec options",options);
    
  List.set("smoother: Aztec params",params);
    
  List.set("smoother: Aztec as solver",false);
#endif
  
  // --- coarse --- ///
  
  List.set("coarse: type","Amesos-KLU");

  List.set("prec type","MGV");

  List.set("print unused",-2);

  for(ParameterList::ConstIterator param=List.begin(); param!=List.end(); param++)
    if ( inList.isParameter(inList.name(param)) == false || OverWrite )
      inList.setEntry( inList.name(param) , inList.entry(param) );
  
  return 0;

} //ML_Epetra::SetDefaultsDD_3Levels_LU()

// ============================================================================
/*! Set values for Maxwell:
 * - \c "default values" \c = \c "maxwell"
   - General
       - \c "max levels" \c = \c 10
       - \c "output" \c = \c 0
       - \c "prec type" \c = \c "MGV"
       - \c "PDE equations" \c = \c 1
       - \c "print unused" \c = \c -2
   - Coarsening
       - \c "increasing or decreasing" \c = \c "decreasing"
       - \c "aggregation: type" \c = \c "Uncoupled-MIS"
       - \c "aggregation: damping factor" \c =   1.3333
       - \c "coarse: max size" \c = \c 75
       - \c "aggregation: threshold" \c = 0.0
       - \c "aggregation: edge prolongator drop threshold" \c = 0.0
   - Smoothing
       - \c "smoother: type" \c = \c "Hiptmair"
       - \c "smoother: sweeps" \c = \c 2
       - \c "smoother: Hiptmair efficient symmetric" \c = \c true
       - \c "smoother: damping factor" \c =  1.0
       - \c "smoother: type" \c = \c "MLS"
       - \c "smoother: MLS polynomial order" \c = \c 2
       - \c "smoother: pre or post" \c = \c "both"
   - Coarse Solution
       - \c "coarse: type" \c = \c "Amesos-KLU"
 */
int ML_Epetra::SetDefaultsMaxwell(ParameterList & inList, 
				  int * options, double * params, bool OverWrite)
{
  ParameterList List;

  List.set("default values","maxwell");
  List.set("max levels",10);
  List.set("output",0);
  
  List.set("PDE equations",1);

  List.set("increasing or decreasing","decreasing");

  // aggregation: Uncoupled-MIS for all levels
  List.set("aggregation: type","Uncoupled-MIS");

  // optimal value for smoothed aggregation
  List.set("aggregation: damping factor",1.3333);

  // relative small coarse size
  List.set("coarse: max size",75);

  // don't forget any element
  List.set("aggregation: threshold",0.0);

  // eigen-analysis
  List.set("eigen-analysis: type","power-method");

  // dropping threshold for small entries in edge prolongator
  List.set("aggregation: edge prolongator drop threshold",0.0);

  // Hiptmair for all levels but coarsest
  List.set("smoother: sweeps",2);

  List.set("smoother: damping factor",1.0);

  List.set("smoother: type","Hiptmair");
  List.set("smoother: Hiptmair efficient symmetric",true);

  //Hiptmair subsmoother options
  List.set("subsmoother: type", "MLS");
  List.set("subsmoother: MLS polynomial order",2);
  List.set("subsmoother: MLS alpha",27.0);

  List.set("smoother: pre or post","both");
  
  // direct solver on coarse problem
//  List.set("coarse: type","SuperLU");
//  Tim Davis' simple serial LU package.  It's part of Amesos
//  itself.
  List.set("coarse: type","Amesos-KLU");

  List.set("prec type","MGV");

  // print unused on proc 0
  List.set("print unused",-2);

  for(ParameterList::ConstIterator param=List.begin(); param!=List.end(); param++)
    if ( inList.isParameter(inList.name(param)) == false || OverWrite )
      inList.setEntry( inList.name(param) , inList.entry(param) );

  return 0;
  
} //ML_Epetra::SetDefaultsMaxwell()

// ============================================================================
/*! Set default values for classical smoothed aggregation:
  - \c "default values" \c = \c "SA"
   - General
        - \c "max levels" \c = \c 10
        - \c "output" \c = \c 0
        - \c "PDE equations" \c = \c 1
        - \c "prec type" \c = \c "MGV"
        - \c "print unused" \c = \c -2
   - Coarsening
        - \c "increasing or decreasing" \c = \c "increasing"
        - \c "aggregation: type" \c = \c "Uncoupled-MIS"
        - \c "aggregation: damping factor" \c = \f$\frac{4}{3}\f$
        - \c "coarse: max size" \c = \c 16
        - \c "aggregation: threshold" \c =  0.0
   - Smoothing
        - \c "smoother: sweeps" \c = \c 2
        - \c "smoother: damping factor" \c = 1.0 
        - \c "smoother: type" \c = \c "symmetric Gauss-Seidel"
        - \c "smoother: pre or post" \c = \c "both"
   - Coarse Solution
        - \c "coarse: type" \c = \c "Amesos-KLU"
 */
int ML_Epetra::SetDefaultsSA(ParameterList & inList, 
			     int * options, double * params, bool OverWrite)
{
  ParameterList List;

  List.set("default values","SA");

  List.set("max levels",10);

  List.set("output",0);
  
  List.set("PDE equations",1);

  List.set("increasing or decreasing","increasing");

  // aggregation: Uncoupled for all levels
  List.set("aggregation: type","Uncoupled-MIS");
  
  // optimal value for smoothed aggregation
  List.set("aggregation: damping factor",1.3333);

  // relative small coarse size
  List.set("coarse: max size",16);

  // don't forget any element
  List.set("aggregation: threshold",0.0);

  // eigen-analysis
  List.set("eigen-analysis: type","power-method");

  // gauss-seidel for all levels
  List.set("smoother: sweeps",2);

  List.set("smoother: damping factor",1.0);

  List.set("smoother: type","symmetric Gauss-Seidel");
  
  List.set("smoother: pre or post","both");
  
  // simplest solver on coarse problem
  List.set("coarse: type","Amesos-KLU");

  List.set("prec type","MGV");

  // print unused on proc 0
  List.set("print unused",-2);

  for(ParameterList::ConstIterator param=List.begin(); param!=List.end(); param++)
    if ( inList.isParameter(inList.name(param)) == false || OverWrite )
      inList.setEntry( inList.name(param) , inList.entry(param) );
  return 0;

} //ML_Epetra::SetDefaultsSA()

// ============================================================================
/*! Set default values for smoothed aggregation for nonsymmetric problems:
  - \c "default values" \c = \c "NSSA"
   - General
        - \c "max levels" \c = \c 10
        - \c "output" \c = \c 0
        - \c "PDE equations" \c = \c 1
        - \c "prec type" \c = \c "MGV"
        - \c "print unused" \c = \c -2
   - Coarsening
        - \c "increasing or decreasing" \c = \c "increasing"
        - \c "aggregation: type" \c = \c "Uncoupled-MIS"
        - \c "aggregation: damping factor" \c = \f$\frac{4}{3}\f$
        - \c "coarse: max size" \c = \c 16
        - \c "aggregation: threshold" \c =  0.0
   - Smoothing
        - \c "smoother: sweeps" \c = \c 2
        - \c "smoother: damping factor" \c = 1.0 
        - \c "smoother: type" \c = \c "Gauss-Seidel"
        - \c "smoother: pre or post" \c = \c "post"
   - Coarse Solution
        - \c "coarse: type" \c = \c "Amesos-KLU"
 */
int ML_Epetra::SetDefaultsNSSA(ParameterList & inList, 
			     int * options, double * params, bool OverWrite)
{
  ParameterList List;

  List.set("default values","NSSA");

  List.set("max levels",10);

  List.set("output",0);
  
  List.set("PDE equations",1);

  List.set("increasing or decreasing","increasing");

  // aggregation: Uncoupled for all levels
  List.set("aggregation: type","Uncoupled-MIS");
  
  // optimal value for smoothed aggregation
  List.set("aggregation: damping factor",1.3333);

  // relative small coarse size
  List.set("coarse: max size",16);

  // don't forget any element
  List.set("aggregation: threshold",0.0);

  // eigen-analysis
  List.set("eigen-analysis: type","power-method");

  // gauss-seidel for all levels
  List.set("smoother: sweeps",2);

  List.set("smoother: damping factor",1.0);

  List.set("smoother: type","Gauss-Seidel");
  
  List.set("smoother: pre or post","post");
  
  // simplest solver on coarse problem
  List.set("coarse: type","Amesos-KLU");

  List.set("prec type","MGV");

  // print unused on proc 0
  List.set("print unused",-2);

  for(ParameterList::ConstIterator param=List.begin(); param!=List.end(); param++)
    if ( inList.isParameter(inList.name(param)) == false || OverWrite )
      inList.setEntry( inList.name(param) , inList.entry(param) );
  return 0;

} //ML_Epetra::SetDefaultsNSSA()

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
