/*!
 *  \file ml_MultiLevelPreconditioner_Defaults.cpp
 *
 *  \brief ML black-box defaults for Epetra_RowMatrix derived classes.
 *
 *  \authors Marzio Sala, Ray Tuminaro, Jonathan Hu, Chris Siefert, Michael Gee
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include "ml_common.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#include "ml_epetra.h"
#include "ml_MultiLevelPreconditioner.h"  // to define ML_Epetra namespace
#include "Teuchos_RefCountPtr.hpp"

using namespace Teuchos;
using namespace ML_Epetra;

// ============================================================================

int ML_Epetra::SetDefaults(string ProblemType, ParameterList & List, 
			   int * ioptions, double * iparams, const bool OverWrite)
{
  RCP<std::vector<int> >    options;
  RCP<std::vector<double> > params;

#ifdef HAVE_ML_AZTECOO
  bool SetDefaults = false;
  if (ioptions == NULL || iparams == NULL)
    SetDefaults = true;

  if (ioptions == NULL)
    options = rcp(new std::vector<int>(AZ_OPTIONS_SIZE));
  else
    options = rcp(new std::vector<int>(ioptions,ioptions+AZ_OPTIONS_SIZE));
  if (iparams  == NULL)
    params  = rcp(new std::vector<double>(AZ_PARAMS_SIZE));
  else
    params = rcp(new std::vector<double>(iparams,iparams+AZ_PARAMS_SIZE));
  if (SetDefaults)
    AZ_defaults(&(*options)[0],&(*params)[0]);
#endif

  if( ProblemType == "SA" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsSA(List, options, params, OverWrite) );
  }
  else if( ProblemType == "DD" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsDD(List, options, params, OverWrite ) );
  }
  else if( ProblemType == "DD-ML" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsDD_3Levels(List, options, params,OverWrite) );
  }
  else if( ProblemType == "maxwell" || ProblemType == "Maxwell" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsMaxwell(List, options, params,OverWrite));
  }
  else if( ProblemType == "NSSA" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsNSSA(List, options, params, OverWrite) );
  }
  else if( ProblemType == "DD-ML-LU" ) {
    ML_CHK_ERR( ML_Epetra::SetDefaultsDD_3Levels_LU(List, options, params, OverWrite ) );
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

void ML_OverwriteDefaults(ParameterList &inList, ParameterList &List,bool OverWrite);

// ============================================================================
/*! Set default values for classical smoothed aggregation.
<br>
    This method should <b>not</b> be called directly.  Instead,
    ML_Epetra::SetDefaults("SA",...) should be used.
   - General
        - \c "default values" \c = \c "SA"
        - \c "max levels" \c = \c 10
        - \c "prec type" \c = \c "MGV"
        - \c "increasing or decreasing" \c = \c "increasing"
   - Grid Transfer
        - \c "aggregation: type" \c = \c "Uncoupled-MIS"
        - \c "aggregation: damping factor" \c = <tt>1.333</tt>
        - \c "eigen-analysis: type" \c = \c "cg"
        - \c "eigen-analysis: iterations" \c = \c 10
   - Smoothing
        - \c "smoother: sweeps" \c = \c 2
        - \c "smoother: damping factor" \c = <tt>1.0</tt>
        - \c "smoother: pre or post" \c = \c "both"
        - \c "smoother: type" \c = \c "symmetric Gauss-Seidel"
   - Coarse Solution
        - \c "coarse: type" \c = \c "Amesos-KLU"
        - \c "coarse: max size" \c = \c 128
 */
int ML_Epetra::SetDefaultsSA(ParameterList & inList, 
                 Teuchos::RCP<std::vector<int> > &options,
                 Teuchos::RCP<std::vector<double> > &params,
                 bool OverWrite)
{
  ParameterList List;

  inList.setName("SA default values");
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
  List.set("coarse: pre or post","post");
  List.set("coarse: sweeps",1);

  ML_OverwriteDefaults(inList, List, OverWrite);
  return 0;
} //ML_Epetra::SetDefaultsSA()

// ============================================================================
/*! Sets the following default values for "DD".
<br>
    Note: This method should <b>not</b> be called directly.  Instead,
    ML_Epetra::SetDefaults("DD",...) should be used.
   - General
        - \c "default values" \c = \c "DD"
        - \c "max levels" \c = \c 2
        - \c "prec type" \c = \c "MGV"
        - \c "increasing or decreasing" \c = \c "increasing"
   - Grid Transfer
        - \c "aggregation: type" \c = \c "METIS"
        - \c "aggregation: local aggregates" \c = \c 1
        - \c "aggregation: damping factor" \c = <tt>1.333</tt>
        - \c "eigen-analysis: type" \c = \c "power-method"
        - \c "eigen-analysis: iterations" \c = \c 20
   - Smoothing
        - \c "smoother: sweeps" \c = \c 1
        - \c "smoother: pre or post" \c = \c "both"
        - \c "smoother: type" \c = \c Aztec". 
        - \c "smoother: Aztec options" \c = \c options
        - \c "smoother: Aztec params" \c = \c params
        - \c "smoother: Aztec as solver" \c = \c false
   - Coarse Solution
        - \c "coarse: type" \c = \c "Amesos-KLU"
        - \c "coarse: max size" \c = \c 128
 */
int ML_Epetra::SetDefaultsDD(ParameterList & inList, 
                 Teuchos::RCP<std::vector<int> > &options,
                 Teuchos::RCP<std::vector<double> > &params,
                 bool OverWrite)
{
  ParameterList List;

  inList.setName("DD default values");
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
  (*options)[AZ_precond] = AZ_dom_decomp;
  (*options)[AZ_subdomain_solve] = AZ_ilu;
  List.set("smoother: Aztec options",options);
  List.set("smoother: Aztec params",params);
  List.set("smoother: Aztec as solver",false);
#endif

  List.set("coarse: type","Amesos-KLU");
  List.set("coarse: max size",128);
  List.set("coarse: pre or post","post");
  List.set("coarse: sweeps",1);

  ML_OverwriteDefaults(inList, List, OverWrite);
  return 0;
} //ML_Epetra::SetDefaultsDD()

// ============================================================================
/*! Sets the following default values for "DD-ML".
<br>
    This method should <b>not</b> be called directly.  Instead,
    ML_Epetra::SetDefaults("DD-ML",...) should be used.
   - General
        - \c "default values" \c = \c "DD-ML"
        - \c "max levels" \c = \c 3
        - \c "prec type" \c = \c "MGV"
        - \c "increasing or decreasing" \c = \c "increasing"
   - Grid Transfer
        - \c "aggregation: type" \c = \c "METIS"
        - \c "aggregation: nodes per aggregate" \c = \c 512
        - \c "aggregation: next-level aggregates per process" \c = \c 128
        - \c "aggregation: damping factor" \c = <tt>1.333</tt>
        - \c "eigen-analysis: type" \c = \c "power-method"
        - \c "eigen-analysis: iterations" \c = \c 20
   - Smoothing
        - \c "smoother: sweeps" \c = \c 1
        - \c "smoother: pre or post" \c = \c "both"
        - \c "smoother: type" \c = \c "Aztec"
        - \c "smoother: Aztec options" \c = \c options
        - \c "smoother: Aztec params" \c = \c params
               - \c "smoother: Aztec as solver" \c = \c false
   - Coarse Solution
        - \c "coarse: type" \c = \c "Amesos-KLU"
        - \c "coarse: max size" \c = \c 128
 */
int ML_Epetra::SetDefaultsDD_3Levels(ParameterList & inList, 
                 Teuchos::RCP<std::vector<int> > &options,
                 Teuchos::RCP<std::vector<double> > &params,
                 bool OverWrite)
{
  ParameterList List;

  inList.setName("DD-ML default values");
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
  (*options)[AZ_precond] = AZ_dom_decomp;
  (*options)[AZ_subdomain_solve] = AZ_ilu;
  List.set("smoother: Aztec options",options);
  List.set("smoother: Aztec params",params);
  List.set("smoother: Aztec as solver",false);
#endif
  
  List.set("coarse: type","Amesos-KLU");
  List.set("coarse: max size",128);
  List.set("coarse: pre or post","post");
  List.set("coarse: sweeps",1);

  ML_OverwriteDefaults(inList, List, OverWrite);
  return 0;
} //ML_Epetra::SetDefaultsDD_3Levels()

// ============================================================================
/*! Set values for Maxwell:
<br>
    Note: This method should <b>not</b> be called directly.  Instead,
    ML_Epetra::SetDefaults("maxwell",...) should be used.
   - General
        - \c "default values" \c = \c "maxwell"
        - \c "max levels" \c = \c 10
        - \c "prec type" \c = \c "MGV"
        - \c "increasing or decreasing" \c = \c "decreasing"
   - Grid Transfer
        - \c "aggregation: type" \c = \c "Uncoupled-MIS"
        - \c "aggregation: damping factor" \c = <tt>1.333</tt>
        - \c "eigen-analysis: type" \c = \c "cg"
        - \c "eigen-analysis: iterations" \c = \c 10
        - \c "aggregation: edge prolongator drop threshold" \c = <tt>0.0</tt>
   - Smoothing
       - \c "smoother: sweeps" \c = \c 1
       - \c "smoother: damping factor" \c =  <tt>1.0</tt>
       - \c "smoother: pre or post" \c = \c "both"
       - \c "smoother: type" \c = \c "Hiptmair"
       - \c "smoother: Hiptmair efficient symmetric" \c = \c true
       - \c "subsmoother: type" \c = \c "Chebyshev"
       - \c "subsmoother: Chebyshev alpha" \c = \c <tt>20.0</tt>
       - \c "subsmoother: node sweeps" \c = \c 4
       - \c "subsmoother: edge sweeps" \c = \c 4
   - Coarse Solution
       - \c "coarse: type" \c = \c "Amesos-KLU"
       - \c "coarse: max size" \c = \c 128
 */
int ML_Epetra::SetDefaultsMaxwell(ParameterList & inList, 
                 Teuchos::RCP<std::vector<int> > &options,
                 Teuchos::RCP<std::vector<double> > &params,
                 bool OverWrite)
{
  ParameterList List;

  inList.setName("Maxwell default values");
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
  List.set("subsmoother: Chebyshev alpha",20.0);
  List.set("subsmoother: node sweeps",4);
  List.set("subsmoother: edge sweeps",4);

  // direct solver on coarse problem
  List.set("coarse: type","Amesos-KLU");
  List.set("coarse: max size",128);
  List.set("coarse: pre or post","post");
  List.set("coarse: sweeps",1);

  ML_OverwriteDefaults(inList, List, OverWrite);
  return 0;
} //ML_Epetra::SetDefaultsMaxwell()

// ============================================================================
/*! Set default values for smoothed aggregation for nonsymmetric problems:
<br>
    Note: This method should <b>not</b> be called directly.  Instead,
    ML_Epetra::SetDefaults("NSSA",...) should be used.
   - General
        - \c "default values" \c = \c "NSSA"
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
        - \c "smoother: damping factor" \c = <tt>0.67</tt>
        - \c "smoother: pre or post" \c = \c "post"
        - \c "smoother: type" \c = \c "symmetric Gauss-Seidel"
   - Coarse Solution
        - \c "coarse: type" \c = \c "Amesos-KLU"
        - \c "coarse: max size" \c = \c 256
 */
int ML_Epetra::SetDefaultsNSSA(ParameterList & inList, 
                 Teuchos::RCP<std::vector<int> > &options,
                 Teuchos::RCP<std::vector<double> > &params,
                 bool OverWrite)
{
  ParameterList List;

  inList.setName("NSSA default values");
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
  List.set("coarse: pre or post","post");
  List.set("coarse: sweeps",1);

  ML_OverwriteDefaults(inList, List, OverWrite);
  return 0;
} //ML_Epetra::SetDefaultsNSSA()

// ============================================================================
/*! Same as SetDefaultsDD(), but used exact LU decompositions on subdomains.
<br>
    Note: This method should <b>not</b> be called directly.  Instead,
    ML_Epetra::SetDefaults("DD-LU",...) should be used.
   - General
        - "default values" \c = \c "DD-LU"
        - \c "max levels" \c = \c 2
        - \c "prec type" \c = \c "MGV"
        - \c "increasing or decreasing" \c = \c "increasing"
   - Grid Transfer
        - \c "aggregation: type" \c = \c "METIS"
        - \c "aggregation: local aggregates" \c = \c 1
        - \c "aggregation: damping factor" \c = <tt>1.333</tt>
        - \c "eigen-analysis: type" \c = \c "power-method"
        - \c "eigen-analysis: iterations" \c = \c 20
   - Smoothing
        - \c "smoother: sweeps" \c = \c 1
        - \c "smoother: pre or post" \c = \c "both"
        - \c "smoother: type" \c = \c "Aztec" 
        - \c "smoother: Aztec options" \c = \c options
        - \c "smoother: Aztec params" \c = \c params
        - \c "smoother: Aztec as solver" \c = \c false
   - Coarse Solution
        - \c "coarse: type" \c = \c "Amesos-KLU"
        - \c "coarse: max size" \c = \c 128
 */
int ML_Epetra::SetDefaultsDD_LU(ParameterList & inList, 
                 Teuchos::RCP<std::vector<int> > &options,
                 Teuchos::RCP<std::vector<double> > &params,
                 bool OverWrite)
{
  ParameterList List;

  inList.setName("DD-LU default values");
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
  (*options)[AZ_precond] = AZ_dom_decomp;
  (*options)[AZ_subdomain_solve] = AZ_lu;
  List.set("smoother: Aztec options",options);
  List.set("smoother: Aztec params",params);
  List.set("smoother: Aztec as solver",false);
#endif

  List.set("coarse: type","Amesos-KLU");
  List.set("coarse: max size",128);
  List.set("coarse: pre or post","post");
  List.set("coarse: sweeps",1);

  ML_OverwriteDefaults(inList, List, OverWrite);
  return 0;
} //ML_Epetra::SetDefaultsDD_LU()


// ============================================================================
/*! Same as SetDefaultsDD_3Levels but with LU factorizations subdomains.
<br>
    Note: This method should <b>not</b> be called directly.  Instead,
    ML_Epetra::SetDefaults("DD-ML-LU",...) should be used.
   - General
        - \c "default values" \c = \c "DD-ML-LU"
        - \c "max levels" \c = \c 3
        - \c "prec type" \c = \c "MGV"
        - \c "increasing or decreasing" \c = \c "increasing"
   - Grid Transfer
        - \c "aggregation: type" \c = \c "METIS"
        - \c "aggregation: nodes per aggregate" \c = \c 512
        - \c "aggregation: next-level aggregates per process" \c = \c 128
        - \c "aggregation: damping factor" \c = <tt>1.333</tt>
        - \c "eigen-analysis: type" \c = \c "power-method"
        - \c "eigen-analysis: iterations" \c = \c 20
   - Smoothing
        - \c "smoother: sweeps" \c = \c 1
        - \c "smoother: pre or post" \c = \c "both"
        - \c "smoother: type" \c = \c "Aztec"
        - \c "smoother: Aztec options" \c = \c options
        - \c "smoother: Aztec params" \c = \c params
               - \c "smoother: Aztec as solver" \c = \c false
   - Coarse Solution
        - \c "coarse: type" \c = \c "Amesos-KLU"
        - \c "coarse: max size" \c = \c 128
 */
int ML_Epetra::SetDefaultsDD_3Levels_LU(ParameterList & inList, 
			     Teuchos::RCP<std::vector<int> > &options,
                 Teuchos::RCP<std::vector<double> > &params,
                 bool OverWrite)
{
  ParameterList List;

  inList.setName("DD-ML-LU default values");
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
  (*options)[AZ_precond] = AZ_dom_decomp;
  (*options)[AZ_subdomain_solve] = AZ_lu;
  List.set("smoother: Aztec options",options);
  List.set("smoother: Aztec params",params);
  List.set("smoother: Aztec as solver",false);
#endif
  List.set("coarse: type","Amesos-KLU");
  List.set("coarse: max size",128);
  List.set("coarse: pre or post","post");
  List.set("coarse: sweeps",1);

  ML_OverwriteDefaults(inList, List, OverWrite);
  return 0;
} //ML_Epetra::SetDefaultsDD_3Levels_LU()

void ML_OverwriteDefaults(ParameterList &inList, ParameterList &List, bool OverWrite)
{
  ParameterList *coarseList=0;
  //Don't create the coarse list if it doesn't already exist!
  if (inList.isSublist("coarse: list"))
    coarseList = &(inList.sublist("coarse: list"));
  for(ParameterList::ConstIterator param=List.begin(); param!=List.end(); param++) {
    string pname = List.name(param);
    if (coarseList && pname.find("coarse: ",0) != string::npos) {
      if (!coarseList->isParameter(pname) || OverWrite)
        coarseList->setEntry(pname , List.entry(param));
    }
    else if ( !inList.isParameter(pname) || OverWrite ) {
      inList.setEntry( pname , List.entry(param) );
    }
  }
} //ML_OverwriteDefaults()



#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
