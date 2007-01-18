/*!
 *  \file ml_MultiLevelPreconditioner_Defaults.cpp
 *
 *  \brief ML black-box defaults for Epetra_RowMatrix derived classes.
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
  else if( ProblemType == "DD" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsDD(List, options, params, OverWrite ) );
  }
  else if( ProblemType == "DD-ML" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsDD_3Levels(List, options,
                params,OverWrite) );
  }
  else if( ProblemType == "maxwell" || ProblemType == "Maxwell" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsMaxwell(List, options, params,OverWrite));
  }
  else if( ProblemType == "NSSA" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsNSSA(List, options, params, OverWrite) );
  }
  else if( ProblemType == "DD-ML-LU" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsDD_3Levels_LU(List, options, params,
                OverWrite ) );
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
/*! Set default values for classical smoothed aggregation:
  - \c "default values" \c = \c "SA"
   - General
        - \c "max levels" \c = \c 10
        - \c "prec type" \c = \c "MGV"
        - \c "increasing or decreasing" \c = \c "increasing"
   - Grid Transfer
        - \c "aggregation: type" \c = \c "Uncoupled-MIS"
        - \c "aggregation: damping factor" \c = \f$ 1.333\f$
        - \c "eigen-analysis: type" \c = \c "cg"
        - \c "eigen-analysis: iterations" \c = \c 10
   - Smoothing
        - \c "smoother: sweeps" \c = \c 2
        - \c "smoother: damping factor" \c = \f$ 1.0\f$
        - \c "smoother: pre or post" \c = \c "both"
        - \c "smoother: type" \c = \c "symmetric Gauss-Seidel"
   - Coarse Solution
        - \c "coarse: type" \c = \c "Amesos-KLU"
        - \c "coarse: max size" \c = \c 128
 */
int ML_Epetra::SetDefaultsSA(ParameterList & inList, 
			     int * options, double * params, bool OverWrite)
{
  ParameterList List;

  List.set("default values","SA");
  List.set("max levels",10);
  List.set("prec type","MGV");
  List.set("increasing or decreasing","increasing");

  List.set("aggregation: type","Uncoupled-MIS");
  List.set("aggregation: damping factor",1.333);
  List.set("eigen-analysis: type","cg");
  List.set("eigen-analysis: iterations",10);

  List.set("smoother: sweeps",2);
  List.set("smoother: damping factor",1.0);
  List.set("smoother: pre or post","both");
  List.set("smoother: type","symmetric Gauss-Seidel");
  
  List.set("coarse: type","Amesos-KLU");
  List.set("coarse: max size",128);

  for(ParameterList::ConstIterator param=List.begin(); param!=List.end(); param++)
    if ( inList.isParameter(inList.name(param)) == false || OverWrite )
      inList.setEntry( inList.name(param) , inList.entry(param) );
  return 0;

} //ML_Epetra::SetDefaultsSA()

// ============================================================================
/*! Sets the following default values for "DD":
  - "default values" = DD"
   - General
        - \c "max levels" = 2
        - \c "prec type" = MGV"
        - \c "increasing or decreasing" = increasing"
   - Grid Transfer
        - \c "aggregation: type" = METIS"
        - \c "aggregation: local aggregates" = 1
        - \c "aggregation: damping factor" = \f$ 1.333\f$
        - \c "eigen-analysis: type" \c = \c "power-method"
        - \c "eigen-analysis: iterations" \c = \c 20
   - Smoothing
        - \c "smoother: sweeps" = 1
        - \c "smoother: pre or post" = both"
        - \c "smoother: type" = Aztec". 
        - \c "smoother: Aztec options" = options
        - \c "smoother: Aztec params" = params
        - \c "smoother: Aztec as solver" = false
   - Coarse Solution
        - \c "coarse: type" = Amesos-KLU"
        - \c "coarse: max size" = 128
 */
int ML_Epetra::SetDefaultsDD(ParameterList & inList, 
			     int * options, double * params, bool OverWrite) 
{
  ParameterList List;

  List.set("default values","DD");
  List.set("max levels",2);
  List.set("prec type","MGV");
  List.set("increasing or decreasing","increasing");

  List.set("aggregation: type","METIS");
  List.set("aggregation: local aggregates",1);
  List.set("aggregation: damping factor",1.333);
  List.set("eigen-analysis: type","power-method");
  List.set("eigen-analysis: iterations",20);

  List.set("smoother: sweeps",1);
  List.set("smoother: pre or post","both");
#ifdef HAVE_ML_AZTECOO
  List.set("smoother: type","Aztec");
  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_subdomain_solve] = AZ_ilu;
  List.set("smoother: Aztec options",options);
  List.set("smoother: Aztec params",params);
  List.set("smoother: Aztec as solver",false);
#endif

  List.set("coarse: type","Amesos-KLU");
  List.set("coarse: max size",128);

  for (ParameterList::ConstIterator param=List.begin(); param!=List.end(); param++)
    if ( inList.isParameter(inList.name(param)) == false || OverWrite )
      inList.setEntry( inList.name(param) , inList.entry(param) );

  return 0;

} //ML_Epetra::SetDefaultsDD()

// ============================================================================
/*! Sets the following default values for "DD-ML"
 - \c "default values" = "DD-ML"
   - General
        - \c "max levels" = 3
        - \c "prec type" = "MGV"
        - \c "increasing or decreasing" = "increasing"
   - Grid Transfer
        - \c "aggregation: type" = "METIS"
        - \c "aggregation: nodes per aggregate" = 512
        - \c "aggregation: next-level aggregates per process" = 128
        - \c "aggregation: damping factor" = 1.333
        - \c "eigen-analysis: type" \c = \c "power-method"
        - \c "eigen-analysis: iterations" \c = \c 20
   - Smoothing
        - \c "smoother: sweeps" = 1
        - \c "smoother: pre or post" = "both"
        - \c "smoother: type" = "Aztec".  
        - \c "smoother: Aztec options" = options
        - \c "smoother: Aztec params" = params
               - \c "smoother: Aztec as solver" = false
   - Coarse Solution
        - \c "coarse: type" = "Amesos-KLU"
        - \c "coarse: max size" = 128
 */
int ML_Epetra::SetDefaultsDD_3Levels(ParameterList & inList, 
				     int * options, double * params, bool OverWrite)
{
  ParameterList List;

  List.set("default values","DD-ML");

  List.set("max levels",3);
  List.set("prec type","MGV");
  List.set("increasing or decreasing","increasing");

  List.set("aggregation: type","METIS");
  List.set("aggregation: nodes per aggregate",512);
  List.set("aggregation: next-level aggregates per process",128);
  List.set("aggregation: damping factor",1.333);
  List.set("eigen-analysis: type","power-method");
  List.set("eigen-analysis: iterations",20);

  List.set("smoother: sweeps",1);
  List.set("smoother: pre or post","both");
#ifdef HAVE_ML_AZTECOO
  List.set("smoother: type","Aztec");
  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_subdomain_solve] = AZ_ilu;
  List.set("smoother: Aztec options",options);
  List.set("smoother: Aztec params",params);
  List.set("smoother: Aztec as solver",false);
#endif
  
  List.set("coarse: type","Amesos-KLU");
  List.set("coarse: max size",128);

  for (ParameterList::ConstIterator param=List.begin(); param!=List.end(); param++)
    if ( inList.isParameter(inList.name(param)) == false || OverWrite )
      inList.setEntry( inList.name(param) , inList.entry(param) );
  
  return 0;

} //ML_Epetra::SetDefaultsDD_3Levels()

// ============================================================================
/*! Set values for Maxwell:
 * - \c "default values" \c = \c "maxwell"
   - General
        - \c "max levels" \c = \c 10
        - \c "prec type" \c = \c "MGV"
        - \c "increasing or decreasing" \c = \c "decreasing"
   - Grid Transfer
        - \c "aggregation: type" \c = \c "Uncoupled-MIS"
        - \c "aggregation: damping factor" \c =   1.333
        - \c "eigen-analysis: type" \c = \c "cg"
        - \c "eigen-analysis: iterations" \c = \c 10
        - \c "aggregation: edge prolongator drop threshold" \c = 0.0
   - Smoothing
       - \c "smoother: sweeps" \c = \c 1
       - \c "smoother: damping factor" \c =  1.0
       - \c "smoother: pre or post" \c = \c "both"
       - \c "smoother: type" \c = \c "Hiptmair"
       - \c "smoother: Hiptmair efficient symmetric" \c = \c true
       - \c "subsmoother: type" \c = \c "Chebyshev"
       - \c "subsmoother: Chebyshev alpha" \c = \c 27.0
       - \c "subsmoother: node sweeps" \c = \c 4
       - \c "subsmoother: edge sweeps" \c = \c 4
   - Coarse Solution
       - \c "coarse: type" \c = \c "Amesos-KLU"
       - \c "coarse: max size" \c = \c 128
 */
int ML_Epetra::SetDefaultsMaxwell(ParameterList & inList, 
			  int * options, double * params, bool OverWrite)
{
  ParameterList List;

  List.set("default values","maxwell");
  List.set("max levels",10);
  List.set("prec type","MGV");
  List.set("increasing or decreasing","decreasing");

  List.set("aggregation: type","Uncoupled-MIS");
  List.set("aggregation: damping factor",1.333);
  List.set("eigen-analysis: type","cg");
  List.set("eigen-analysis: iterations",10);
  // dropping threshold for small entries in edge prolongator
  List.set("aggregation: edge prolongator drop threshold",0.0);

  List.set("smoother: sweeps",1);
  List.set("smoother: damping factor",1.0);
  List.set("smoother: pre or post","both");
  List.set("smoother: type","Hiptmair");
  List.set("smoother: Hiptmair efficient symmetric",true);
  List.set("subsmoother: type", "Chebyshev");           //Hiptmair subsmoother options
  List.set("subsmoother: Chebyshev alpha",27.0);
  List.set("subsmoother: node sweeps",4);
  List.set("subsmoother: edge sweeps",4);

  // direct solver on coarse problem
  List.set("coarse: type","Amesos-KLU");
  List.set("coarse: max size",128);

  for(ParameterList::ConstIterator param=List.begin(); param!=List.end(); param++)
    if ( inList.isParameter(inList.name(param)) == false || OverWrite )
      inList.setEntry( inList.name(param) , inList.entry(param) );

  return 0;
  
} //ML_Epetra::SetDefaultsMaxwell()

// ============================================================================
/*! Set default values for smoothed aggregation for nonsymmetric problems:
  - \c "default values" \c = \c "NSSA"
   - General
        - \c "max levels" \c = \c 10
        - \c "prec type" \c = \c "MGW"
        - \c "increasing or decreasing" \c = \c "increasing"
   - Grid Transfer
        - \c "aggregation: type" \c = \c "Uncoupled-MIS"
        - \c "energy minimization: enable" \c = \c true
        - \c "eigen-analysis: type" \c = \c "power-method"
        - \c "eigen-analysis: iterations" \c = \c 20
   - Smoothing
        - \c "smoother: sweeps" \c = \c 4
        - \c "smoother: damping factor" \c = 0.67
        - \c "smoother: pre or post" \c = \c "post"
        - \c "smoother: type" \c = \c "symmetric Gauss-Seidel"
   - Coarse Solution
        - \c "coarse: type" \c = \c "Amesos-KLU"
        - \c "coarse: max size" \c = \c 256
 */
int ML_Epetra::SetDefaultsNSSA(ParameterList & inList, 
			     int * options, double * params, bool OverWrite)
{
  ParameterList List;

  List.set("default values","NSSA");
  List.set("max levels",10);
  List.set("prec type","MGW");
  List.set("increasing or decreasing","increasing");

  List.set("aggregation: type","Uncoupled-MIS");
  List.set("energy minimization: enable",true);
  List.set("eigen-analysis: type","power-method");
  List.set("eigen-analysis: iterations",20);

  List.set("smoother: sweeps",4);
  List.set("smoother: damping factor",.67);
  List.set("smoother: pre or post","post");
  List.set("smoother: type","symmetric Gauss-Seidel");

  List.set("coarse: type","Amesos-KLU");
  List.set("coarse: max size",256);


  for(ParameterList::ConstIterator param=List.begin(); param!=List.end(); param++)
    if ( inList.isParameter(inList.name(param)) == false || OverWrite )
      inList.setEntry( inList.name(param) , inList.entry(param) );
  return 0;

} //ML_Epetra::SetDefaultsNSSA()

// ============================================================================
/*! Same as SetDefaultsDD(), but used exact LU decompositions on subdomains.
  - "default values" = DD-LU"
   - General
        - \c "max levels" = 2
        - \c "prec type" = MGV"
        - \c "increasing or decreasing" = increasing"
   - Grid Transfer
        - \c "aggregation: type" = METIS"
        - \c "aggregation: local aggregates" = 1
        - \c "aggregation: damping factor" = \f$ 1.333\f$
        - \c "eigen-analysis: type" \c = \c "power-method"
        - \c "eigen-analysis: iterations" \c = \c 20
   - Smoothing
        - \c "smoother: sweeps" = 1
        - \c "smoother: pre or post" = both"
        - \c "smoother: type" = Aztec". 
        - \c "smoother: Aztec options" = options
        - \c "smoother: Aztec params" = params
        - \c "smoother: Aztec as solver" = false
   - Coarse Solution
        - \c "coarse: type" = Amesos-KLU"
        - \c "coarse: max size" = 128
 */
int ML_Epetra::SetDefaultsDD_LU(ParameterList & inList, 
				int * options, double * params, bool OverWrite) 
{
  ParameterList List;

  List.set("default values","DD-LU");
  List.set("max levels",2);
  List.set("prec type","MGV");
  List.set("increasing or decreasing","increasing");

  List.set("aggregation: type","METIS");
  List.set("aggregation: local aggregates",1);
  List.set("aggregation: damping factor",1.333);
  List.set("eigen-analysis: type","power-method");
  List.set("eigen-analysis: iterations",20);

  List.set("smoother: sweeps",1);
  List.set("smoother: pre or post","both");

#ifdef HAVE_ML_AZTECOO
  List.set("smoother: type","Aztec");
  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_subdomain_solve] = AZ_lu;
  List.set("smoother: Aztec options",options);
  List.set("smoother: Aztec params",params);
  List.set("smoother: Aztec as solver",false);
#endif

  List.set("coarse: type","Amesos-KLU");
  List.set("coarse: max size",128);


  for(ParameterList::ConstIterator param=List.begin(); param!=List.end(); param++)
    if ( inList.isParameter(inList.name(param)) == false || OverWrite )
      inList.setEntry( inList.name(param) , inList.entry(param) );

  return 0;

} //ML_Epetra::SetDefaultsDD_LU()


// ============================================================================
/*! Same as SetDefaultsDD_3Levels but with LU factorizations subdomains.
 - \c "default values" = "DD-ML-LU"
   - General
        - \c "max levels" = 3
        - \c "prec type" = "MGV"
        - \c "increasing or decreasing" = "increasing"
   - Grid Transfer
        - \c "aggregation: type" = "METIS"
        - \c "aggregation: nodes per aggregate" = 512
        - \c "aggregation: next-level aggregates per process" = 128
        - \c "aggregation: damping factor" = 1.333
        - \c "eigen-analysis: type" \c = \c "power-method"
        - \c "eigen-analysis: iterations" \c = \c 20
   - Smoothing
        - \c "smoother: sweeps" = 1
        - \c "smoother: pre or post" = "both"
        - \c "smoother: type" = "Aztec".  
        - \c "smoother: Aztec options" = options
        - \c "smoother: Aztec params" = params
               - \c "smoother: Aztec as solver" = false
   - Coarse Solution
        - \c "coarse: type" = "Amesos-KLU"
        - \c "coarse: max size" = 128
 */
int ML_Epetra::SetDefaultsDD_3Levels_LU(ParameterList & inList, 
					int * options, double * params, bool OverWrite)
{
  ParameterList List;

  List.set("default values","DD-ML-LU");
  List.set("max levels",3);
  List.set("prec type","MGV");
  List.set("increasing or decreasing","increasing");

  List.set("aggregation: type","METIS");
  List.set("aggregation: nodes per aggregate",512);
  List.set("aggregation: next-level aggregates per process",128);
  List.set("aggregation: damping factor",1.333);
  
  List.set("smoother: sweeps",1);
  List.set("smoother: pre or post","both");
#ifdef HAVE_ML_AZTECOO
  List.set("smoother: type","Aztec");
  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_subdomain_solve] = AZ_lu;
  List.set("smoother: Aztec options",options);
  List.set("smoother: Aztec params",params);
  List.set("smoother: Aztec as solver",false);
#endif
  List.set("coarse: type","Amesos-KLU");
  List.set("coarse: max size",128);

  for(ParameterList::ConstIterator param=List.begin(); param!=List.end(); param++)
    if ( inList.isParameter(inList.name(param)) == false || OverWrite )
      inList.setEntry( inList.name(param) , inList.entry(param) );
  
  return 0;

} //ML_Epetra::SetDefaultsDD_3Levels_LU()



#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
