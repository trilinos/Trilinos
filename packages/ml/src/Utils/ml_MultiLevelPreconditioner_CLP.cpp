/*!
 *  \file ml_MultiLevelPreconditioner_CLP.cpp
 *
 *  \brief Sets up a parameters list from Teuchos CommandLineParser.
 *
 *  \author Marzio Sala, SNL, 9214
 *
 *  \date Last updated on Nov-04
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
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "ml_MultiLevelPreconditioner.h"

using namespace Teuchos;
using namespace ML_Epetra;

// ================================================ ====== ==== ==== == =

int ML_Epetra::SetParameters(int argc, char* argv[],
                             Teuchos::ParameterList& List)
{
 
  Teuchos::CommandLineProcessor CLP;
  string ml_defaults = "not-set";
  CLP.setOption("ml-defaults",&ml_defaults,"Default values for ML");

  int ml_levels = List.get("max levels",5);
  CLP.setOption("ml-levels",&ml_levels,"Maximum number of levels");
  int ml_output = List.get("output",5);
  CLP.setOption("ml-output",&ml_output,"Output level (from 0 to 10)");
  // damping factor
  double ml_damping = List.get("aggregation: damping factor", 1.333);
  CLP.setOption("ml-damping",&ml_damping,"Damping for prolongator");
  // smoothers type
  string ml_smoother_type = List.get("smoother: type","Gauss-Seidel");
  CLP.setOption("ml-smoother-type",&ml_smoother_type,"Smoother type");
  string ml_smoother_type_0 = List.get("smoother: type (level 0)","Gauss-Seidel");
  CLP.setOption("ml-smoother-type-0",&ml_smoother_type_0,"Smoother type for level 0");
  string ml_smoother_type_1 = List.get("smoother: type (level 1)","Gauss-Seidel");
  CLP.setOption("ml-smoother-type-1",&ml_smoother_type_1,"Smoother type for level 1");
  string ml_smoother_type_2 = List.get("smoother: type (level 2)","Gauss-Seidel");
  CLP.setOption("ml-smoother-type-2",&ml_smoother_type_2,"Smoother type for level 2");
  string ml_smoother_type_3 = List.get("smoother: type (level 3)","Gauss-Seidel");
  CLP.setOption("ml-smoother-type-3",&ml_smoother_type_3,"Smoother type for level 2");
  string ml_smoother_type_4 = List.get("smoother: type (level 4)","Gauss-Seidel");
  CLP.setOption("ml-smoother-type-4",&ml_smoother_type_4,"Smoother type for level 2");
  // smoother sweeps
  int ml_smoother_sweeps = List.get("smoother: sweeps",1);
  CLP.setOption("ml-smoother-sweeps",&ml_smoother_sweeps,"Sweeps");
  // smoother damping
  double ml_smoother_damping = List.get("smoother: damping factor", 0.67);
  CLP.setOption("ml-smoother-damping",&ml_smoother_damping,"damping factor");
  double ml_smoother_damping_0 = List.get("smoother: damping factor (level 0)", 0.67);
  CLP.setOption("ml-smoother-damping-0",&ml_smoother_damping_0,"Damping for level 0");
  double ml_smoother_damping_1 = List.get("smoother: damping factor (level 1)", 0.67);
  CLP.setOption("ml-smoother-damping-1",&ml_smoother_damping_1,"Damping for level 1");
  double ml_smoother_damping_2 = List.get("smoother: damping factor (level 2)", 0.67);
  CLP.setOption("ml-smoother-damping-2",&ml_smoother_damping_2,"Damping for level 2");
  double ml_smoother_damping_3 = List.get("smoother: damping factor (level 3)", 0.67);
  CLP.setOption("ml-smoother-damping-3",&ml_smoother_damping_3,"Damping for level 3");
  double ml_smoother_damping_4 = List.get("smoother: damping factor (level 4)", 0.67);
  CLP.setOption("ml-smoother-damping-4",&ml_smoother_damping_4,"Damping for level 4");
  // pre or post smoother
  string ml_smoother_pre_or_post = List.get("smoother: pre or post","both");
  CLP.setOption("ml-smoother-pre-or-post",&ml_smoother_pre_or_post,"pre, post smoother or both");
  // aggregation type
  string ml_aggr_type = List.get("aggregation: type","Uncoupled");
  CLP.setOption("ml-aggr-type",&ml_aggr_type,"Aggregation scheme");
  string ml_aggr_type_0 = List.get("aggregation: type","Uncoupled");
  CLP.setOption("ml-aggr-type-0",&ml_aggr_type_0,"Aggregation scheme (level 0)");
  string ml_aggr_type_1 = List.get("aggregation: type","Uncoupled");
  CLP.setOption("ml-aggr-type-1",&ml_aggr_type_1,"Aggregation scheme (level 1)");
  string ml_aggr_type_2 = List.get("aggregation: type","Uncoupled");
  CLP.setOption("ml-aggr-type-2",&ml_aggr_type_2,"Aggregation scheme (level 2)");
  string ml_aggr_type_3 = List.get("aggregation: type","Uncoupled");
  CLP.setOption("ml-aggr-type-3",&ml_aggr_type_3,"Aggregation scheme (level 3)");
  string ml_aggr_type_4 = List.get("aggregation: type","Uncoupled");
  CLP.setOption("ml-aggr-type-4",&ml_aggr_type_4,"Aggregation scheme (level 4)");
  // aggregation options
  int ml_aggr_local = List.get("aggregation: local aggregates",16);
  CLP.setOption("ml-aggr-local",&ml_aggr_local,"Local number of aggregates");
  int ml_aggr_global = List.get("aggregation: global aggregates",16);
  CLP.setOption("ml-aggr-global",&ml_aggr_global,"Global number of aggregates");
  // aggregation (nodes per aggregate)
  int ml_aggr_npa = List.get("aggregation: nodes per aggregate",16);
  CLP.setOption("ml-aggr-npa",&ml_aggr_npa,"Number of nodes per aggregate");
  int ml_aggr_npa_0 = List.get("aggregation: nodes per aggregate (level 0)",16);
  CLP.setOption("ml-aggr-npa-0",&ml_aggr_npa_0,"Number of nodes per aggregate (level 0)");
  int ml_aggr_npa_1 = List.get("aggregation: nodes per aggregate (level 1)",16);
  CLP.setOption("ml-aggr-npa-1",&ml_aggr_npa_1,"Number of nodes per aggregate (level 1)");
  int ml_aggr_npa_2 = List.get("aggregation: nodes per aggregate (level 2)",16);
  CLP.setOption("ml-aggr-npa-2",&ml_aggr_npa_2,"Number of nodes per aggregate (level 2)");
  int ml_aggr_npa_3 = List.get("aggregation: nodes per aggregate (level 3)",16);
  CLP.setOption("ml-aggr-npa-3",&ml_aggr_npa_3,"Number of nodes per aggregate (level 3)");
  int ml_aggr_npa_4 = List.get("aggregation: nodes per aggregate (level 4)",16);
  CLP.setOption("ml-aggr-npa-4",&ml_aggr_npa_4,"Number of nodes per aggregate (level 4)");
  // coarse type
  string ml_coarse_type = List.get("coarse: type","Amesos-KLU");
  CLP.setOption("ml-coarse-type",&ml_coarse_type,"Coarse type");
  // prolongator smoother
  string prol_damping = List.get("R and P damping: damping", "classic");
  CLP.setOption("ml-prol-damping",&prol_damping);

  CLP.throwExceptions(false);
  // allow users to specify other options for other packages
  CLP.recogniseAllOptions(false);
  CLP.parse(argc,argv);

  // now insert the read values into `List'
  if (ml_defaults != "not-set")
    SetDefaults(ml_defaults,List);
  List.set("max levels", ml_levels);
  List.set("output", ml_output);
  List.set("smoother: type", ml_smoother_type);
  List.set("smoother: type (level 0)", ml_smoother_type_0);
  List.set("smoother: type (level 1)", ml_smoother_type_1);
  List.set("smoother: type (level 2)", ml_smoother_type_2);
  List.set("smoother: type (level 3)", ml_smoother_type_3);
  List.set("smoother: type (level 4)", ml_smoother_type_4);
  List.set("smoother: sweeps", ml_smoother_sweeps);
  List.set("smoother: damping factor", ml_smoother_damping);
  List.set("smoother: damping factor (level 0)", ml_smoother_damping_0);
  List.set("smoother: damping factor (level 1)", ml_smoother_damping_1);
  List.set("smoother: damping factor (level 2)", ml_smoother_damping_2);
  List.set("smoother: damping factor (level 3)", ml_smoother_damping_3);
  List.set("smoother: damping factor (level 4)", ml_smoother_damping_4);
  List.set("smoother: pre or post", ml_smoother_pre_or_post);
  List.set("aggregation: type", ml_aggr_type);
  List.set("aggregation: type (level 0)", ml_aggr_type_0);
  List.set("aggregation: type (level 1)", ml_aggr_type_1);
  List.set("aggregation: type (level 2)", ml_aggr_type_2);
  List.set("aggregation: type (level 3)", ml_aggr_type_3);
  List.set("aggregation: type (level 4)", ml_aggr_type_4);
  List.set("aggregation: local aggregates", ml_aggr_local);
  List.set("aggregation: nodes per aggregate", ml_aggr_npa);
  List.set("aggregation: nodes per aggregate (level 0)", ml_aggr_npa_0); 
  List.set("aggregation: nodes per aggregate (level 1)", ml_aggr_npa_1);
  List.set("aggregation: nodes per aggregate (level 2)", ml_aggr_npa_2);
  List.set("aggregation: nodes per aggregate (level 3)", ml_aggr_npa_3);
  List.set("aggregation: nodes per aggregate (level 4)", ml_aggr_npa_4);
  List.set("aggregation: damping factor", ml_damping);
  List.set("coarse: type", ml_coarse_type);
  if (prol_damping == "classic")
    List.set("R and P smoothing: type", "classic");
  else {
    List.set("R and P smoothing: type", "advanced");
    List.set("R and P smoothing: damping", prol_damping);
    List.set("aggregation: compute field of values", true);
  }

  return(0);
  
}

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/

