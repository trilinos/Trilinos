/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_common.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

#include "Epetra_CrsMatrix.h"
#include "ml_ValidateParameters.h"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#ifdef HAVE_ML_IFPACK
#include "Ifpack_config.h"
#endif

# ifdef HAVE_ML_TekoSmoothers
namespace Teko {
class InverseLibrary;
}
#endif

using namespace ML_Epetra;

bool ML_Epetra::ValidateMLPParameters(const Teuchos::ParameterList &inList, int depth){
  using Teuchos::ParameterList;
  using Teuchos::Exceptions::InvalidParameterName;
  using Teuchos::Exceptions::InvalidParameterType;
  using Teuchos::Exceptions::InvalidParameterValue;
  using std::cout;
  using std::endl;
  using std::string;

  ParameterList List,*validList;
  bool rv=true;

  /* Build a copy of the list to be validated. */
  for (ParameterList::ConstIterator param = inList.begin(); param != inList.end(); ++param) {
    const std::string pname=inList.name(param);
    if (pname.find("user-defined function",0) == std::string::npos) {
      List.setEntry(pname,inList.entry(param));
    }
  }

  List.setName(inList.name());

  /* Get Defaults + Validate */
  try {
    validList = GetValidMLPParameters();
  }
  catch(...) {
    std::cout << "Error in GetValidMLPParameters.  Please report this bug to the ML "
      "developers." << std::endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_FAILURE);
  }
  try {
    List.validateParameters (*validList, depth, Teuchos::VALIDATE_USED_ENABLED,
			     Teuchos::VALIDATE_DEFAULTS_DISABLED);
  }
#ifdef HAVE_IFPACK_DYNAMIC_FACTORY
  catch(InvalidParameterName &excpt)  {/*rv=false; std::cout<<excpt.what()<<std::endl;*/}
#else
  catch(InvalidParameterName &excpt)  {rv=false; std::cout<<excpt.what()<<std::endl;}
#endif
  catch(InvalidParameterType &excpt)  {rv=false; std::cout<<excpt.what()<<std::endl;}
  catch(InvalidParameterValue &excpt) {rv=false; std::cout<<excpt.what()<<std::endl;}
  catch(...) {rv=false;}
  delete validList;
  return rv;
}

/******************************************************************************/

void ML_Epetra::SetValidSmooParams(Teuchos::ParameterList *PL, Teuchos::Array<std::string> &smoothers)
{
  using Teuchos::AnyNumberParameterEntryValidator;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::setIntParameter;
  using Teuchos::setDoubleParameter;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple; // Different than std::tuple.
  using std::string;
  using std::vector;

  ParameterList dummy;
  // prevent Teuchos from converting parameter types
  AnyNumberParameterEntryValidator::AcceptedTypes intParam(false),
                                                  dblParam(false),
                                                  strParam(false);
  intParam.allowInt(true);
  dblParam.allowDouble(true);
  strParam.allowString(true);

  /* Smoothing Options (Section 6.4.4) */
  setStringToIntegralParameter<int>("smoother: type", std::string("Chebyshev"),
				    "Smoothing algorithm",smoothers,PL);
  setIntParameter("smoother: sweeps", 2, "Number of smoothing sweeps", PL, intParam);
  setDoubleParameter("smoother: line detection threshold",-1.0,"Smoother line detection threshold",PL,dblParam);
  setIntParameter("smoother: line direction nodes", -1, "Number of mesh points in z direction", PL, intParam);
  setDoubleParameter("smoother: damping factor",1.0,"Smoother damping factor",PL,dblParam);
  setDoubleParameter("smoother: Chebyshev eig boost", 1.1, "factor to scale eig_max when using Cheby smoothing", PL,dblParam);
  setStringToIntegralParameter<int>("smoother: pre or post","both","Smooth before/after coarse correction, or both",tuple<std::string>("pre","post","both"),PL);
  setStringToIntegralParameter<int>("smoother: line orientation","use coordinates","indicates grid points are numbered either horizontally or veritcally for extruded meshes", tuple<std::string>("horizontal","vertical","use coordinates"),PL);
  setStringToIntegralParameter<int>("smoother: line GS Type","symmetric","use symmetric GS, forward GS, or pre-forward and post-backward GS",tuple<std::string>("symmetric","standard","efficient symmetric"),PL);
  setStringToIntegralParameter<int>("smoother: line group dofs","separate","group or separate dofs per node", tuple<std::string>("group","separate"),PL);


#ifdef HAVE_ML_AZTECOO
  RCP<std::vector<int> > options = rcp(new std::vector<int>(AZ_OPTIONS_SIZE));
  RCP<std::vector<double> > params = rcp(new std::vector<double>(AZ_PARAMS_SIZE));
  PL->set("smoother: Aztec options",options);
  PL->set("smoother: Aztec params",params);
#endif
  PL->set("smoother: Aztec as solver",false);
  setDoubleParameter("smoother: Chebyshev alpha",20.0,"Damping radius for Chebyshev",PL,dblParam);
  setDoubleParameter("smoother: MLS alpha",20.0,"Damping radius for Chebyshev",PL,dblParam);
  PL->set("smoother: user-defined function",(int(*)(ML_Smoother*,int,double*,int,double*))0);
  PL->set("smoother: user-defined name",std::string("User-defined"));
  PL->set("smoother: Hiptmair efficient symmetric",true);
  setStringToIntegralParameter<int>("subsmoother: type","Chebyshev","Subsmoother algorithm (Maxwell)",tuple<std::string>("Chebyshev","symmetric Gauss-Seidel","MLS"),PL);
  setDoubleParameter("subsmoother: Chebyshev alpha",20.0,"Damping radius for Chebshev",PL,dblParam);
  setDoubleParameter("subsmoother: MLS alpha",20.0,"Damping radius for Chebshev",PL,dblParam);
  setDoubleParameter("subsmoother: damping factor",1.333,"Damping factor for symmetric Gauss-Seidel",PL,dblParam);
  setIntParameter("subsmoother: edge sweeps",4,"Number of edge smoothing sweeps",PL,intParam);
  setIntParameter("subsmoother: node sweeps",4,"Number of node smoothing sweeps",PL,intParam);
# ifdef HAVE_ML_PETSC
  void *petscKSP;
  PL->set("smoother: petsc ksp",petscKSP);
# endif
# ifdef HAVE_ML_TekoSmoothers
  {
    RCP<ParameterList> nullList;
    PL->set("smoother: teko filename",std::string("teko_smoother.xml"));
    PL->set<RCP<const Teko::InverseLibrary> >("smoother: teko inverse library",Teuchos::null);
    PL->set("smoother: teko parameter list",nullList);
    PL->set("smoother: teko inverse",std::string("Amesos"));
    PL->set("smoother: teko is blocked",0);
  }
# endif
  /* Unlisted Options */
  setIntParameter("smoother: self overlap",0,"experimental option",PL,intParam);
  /* Unlisted Options that should probably be listed */
  PL->set("smoother: self list",dummy);
  PL->sublist("smoother: self list").disableRecursiveValidation();
  setDoubleParameter("coarse: add to diag", 0.0,"Unlisted option",PL,dblParam);
  // From ml_Multilevel_Smoothers.cpp:
  setIntParameter("smoother: ParaSails matrix",0,"Unlisted option",PL,intParam);
  setIntParameter("smoother: ParaSails levels",0,"Unlisted option",PL,intParam);
  setDoubleParameter("smoother: ParaSails threshold",0.01,"Unlisted option",PL,dblParam);
  setDoubleParameter("smoother: ParaSails filter",0.05,"Unlisted option",PL,dblParam);
  setDoubleParameter("smoother: ParaSails load balancing",0,"Unlisted option",PL,dblParam);
  setIntParameter("smoother: ParaSails factorized",0,"Unlisted option",PL,intParam);
  PL->set("smoother: ifpack list",dummy);
  PL->sublist("smoother: ifpack list").disableRecursiveValidation();
  PL->set("smoother: ifpack type",std::string(""));
  setIntParameter("smoother: ifpack overlap",0,"Unlisted option",PL,intParam);
  setDoubleParameter("smoother: ifpack level-of-fill",0.0,"Unlisted option",PL,dblParam);
  setDoubleParameter("smoother: ifpack relative threshold",1.0,"Unlisted option",PL,dblParam);
  setDoubleParameter("smoother: ifpack absolute threshold",0.0,"Unlisted option",PL,dblParam);
  /* Unlisted options that should probably go away */
  setIntParameter("smoother: polynomial order",2,"Unlisted option",PL,intParam);
  setIntParameter("smoother: MLS polynomial order",2,"Unlisted option",PL,intParam);
  setIntParameter("coarse: polynomial order",2,"Unlisted option",PL,intParam);
  setIntParameter("coarse: MLS polynomial order",2,"Unlisted option",PL,intParam);

  /*
     Coarse grid options.

     NOTE:  "inList" has already been run through ML_CreateSublist().  In
     particular, all coarse level options are now in a sublist called "coarse:
     list".  Furthermore, in the coares list any coarse level option, i.e.,
     starts with "coarse:", has been converted to a smoother option of the
     same name.  So for example, "coarse: type" is now "smoother: type".
     Admittedly, this does create some peculiar options, such as
     "smoother: max size".

     The upshot is that all valid coarse level options from this point in the
     code execution start with "smoother:".
  */

  /* Coarsest Grid Options (Section 6.4.5) */
  setIntParameter("smoother: max size",75,"Max coarse grid size",PL,intParam);
  //setStringToIntegralParameter<int>("smoother: type","Chebyshev","Coarse solver algorithm",smoothers,PL);
  //setStringToIntegralParameter<int>("smoother: pre or post","post","When to smooth on coarse grid",tuple<std::string>("pre","post","both"),PL);
  //setIntParameter("smoother: sweeps",2,"Number of coarse sweeps",PL,intParam);
  //PL->set("smoother: user-defined function",(int(*)(ML_Smoother*,int,double*,int,double*))0);
  //PL->set("smoother: user-defined name",std::string("User-defined"));
  //setDoubleParameter("smoother: damping factor",1.0,"Coarse smoother damping factor",PL,dblParam);
  setStringToIntegralParameter<int>("smoother: subsmoother type","Chebyshev","Coarse grid subsmoother (Maxwell)",tuple<std::string>("Chebyshev","symmetric Gauss-Seidel","MLS"),PL);
  setIntParameter("smoother: edge sweeps",4,"Number of coarse edge smoothing sweeps",PL,intParam);
  setIntParameter("smoother: node sweeps",4,"Number of coarse node smoothing sweeps",PL,intParam);
  //setDoubleParameter("smoother: Chebyshev alpha",2.0,"Damping radius for Chebshev",PL,dblParam);
  //setDoubleParameter("smoother: MLS alpha",2.0,"Damping radius for Chebshev",PL,dblParam);
  setIntParameter("smoother: max processes",-1,"Maximum number of procs for coarse solve (Superludist/MUMPS)",PL,intParam);


  /* EXPERIMENTAL options*/
  PL->set("smoother: Gauss-Seidel efficient symmetric",false);
  setIntParameter("smoother: Block Chebyshev number of blocks",-1,"Number of blocks to use with Block Chebyshev",PL,intParam);
  PL->set("smoother: Block Chebyshev block starts",(int*)0);
  PL->set("smoother: Block Chebyshev block list",(int*)0);
  PL->set("smoother: chebyshev solve normal equations",false);
  PL->set("smoother: use l1 Gauss-Seidel",false);
  PL->set("use crs matrix storage",false);


  /* Coarse IFPACK Solvers - experimental */
  //PL->set("smoother: ifpack list",dummy);
  //PL->set("smoother: ifpack type",std::string(""));
  //setIntParameter("smoother: ifpack overlap",0,"Unlisted option",PL,intParam);
  //setDoubleParameter("smoother: ifpack level-of-fill",0.0,"Unlisted option",PL,dblParam);
  //setDoubleParameter("smoother: ifpack relative threshold",1.0,"Unlisted option",PL,dblParam);
  //setDoubleParameter("smoother: ifpack absolute threshold",0.0,"Unlisted option",PL,dblParam);

} //SetValidSmooParams()

/******************************************************************************/

void ML_Epetra::SetValidAggrParams(Teuchos::ParameterList *PL)
{
  using Teuchos::AnyNumberParameterEntryValidator;
  using Teuchos::setIntParameter;
  using Teuchos::setDoubleParameter;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;
  using std::string;

  // prevent Teuchos from converting parameter types
  AnyNumberParameterEntryValidator::AcceptedTypes intParam(false),
                                                  dblParam(false),
                                                  strParam(false);
  intParam.allowInt(true);
  dblParam.allowDouble(true);
  strParam.allowString(true);

  /* Aggregation and Prolongator Options (Section 6.4.3) */
  setStringToIntegralParameter<int>("aggregation: type", "Uncoupled", "Aggregation algorithm", tuple<std::string>("Uncoupled","Coupled","MIS","Uncoupled-MIS","METIS","ParMETIS","Zoltan","user"),PL);
  setDoubleParameter("aggregation: threshold",0.0,"Dropping for aggregation",PL,dblParam);
  setDoubleParameter("aggregation: damping factor",1.3333,"Damping factor for smoothed aggregation",PL,dblParam);
  setIntParameter("aggregation: smoothing sweeps",1,"Number of sweeps for prolongator smoother",PL,intParam);
  setIntParameter("aggregation: global aggregates",0,"Number of global aggregates (METIS/ParMETIS)",PL,intParam);
  setIntParameter("aggregation: local aggregates",0,"Number of local aggregates (METIS/ParMETIS)",PL,intParam);
  setIntParameter("aggregation: nodes per aggregate",0,"Number of nodes per aggregate (METIS/ParMETIS)",PL,intParam);
  setIntParameter("aggregation: next-level aggregates per process",128,"Number of next-level rows / process (ParMETIS)",PL,intParam);
  PL->set("aggregation: use tentative restriction",false);
  PL->set("aggregation: symmetrize",false);

  /* Highly experimental */
  PL->set("aggregation: respect materials",false);
  PL->set("aggregation: material type",(int*)0);

}//SetValidAggrParams()

/******************************************************************************/

Teuchos::ParameterList * ML_Epetra::GetValidMLPParameters(){
  using Teuchos::AnyNumberParameterEntryValidator;
  using Teuchos::Array;
  using Teuchos::ParameterList;
  using Teuchos::setIntParameter;
  using Teuchos::setDoubleParameter;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;
  using std::string;

  ParameterList dummy;
  ParameterList * PL = new ParameterList;

  // prevent Teuchos from converting parameter types
  AnyNumberParameterEntryValidator::AcceptedTypes intParam(false),
                                                  dblParam(false),
                                                  strParam(false);
  intParam.allowInt(true);
  dblParam.allowDouble(true);
  strParam.allowString(true);

  /* Allocate List for Smoothing Options */
# if defined(HAVE_ML_PETSC) && defined(HAVE_ML_SUPERLU4_0)
  const int num_smoothers=32;
# elif defined(HAVE_ML_PETSC) || defined(HAVE_ML_SUPERLU4_0)
  const int num_smoothers=31;
#elif defined(HAVE_ML_TekoSmoothers)
  const int num_smoothers=31; // won't work with SUPERLU or PETSC!
# else
  const int num_smoothers=30;
# endif
  const char* smoother_strings[num_smoothers]={"Aztec","IFPACK","Jacobi",
   "ML symmetric Gauss-Seidel","symmetric Gauss-Seidel","ML Gauss-Seidel",
   "Gauss-Seidel","block Gauss-Seidel","symmetric block Gauss-Seidel",
   "Chebyshev","MLS","Hiptmair","Amesos-KLU","Amesos-Superlu",
   "Amesos-UMFPACK","Amesos-Superludist","Amesos-MUMPS","user-defined",
   "SuperLU","IFPACK-Chebyshev","self","do-nothing","IC","ICT","ILU","ILUT",
   "Block Chebyshev","IFPACK-Block Chebyshev","line Jacobi","line Gauss-Seidel"
#  ifdef HAVE_ML_PETSC
   ,"petsc"
#  endif
#  ifdef HAVE_ML_SUPERLU4_0
   ,"SILU"
//#else
//#error "No SuperLU for you!"
#  endif
#  ifdef HAVE_ML_TekoSmoothers
   ,"teko"
#  endif
   };
  Array<std::string> smoothers(num_smoothers);
  for(int i = 0; i<num_smoothers; i++) {
    smoothers[i] = smoother_strings[i];
  }

  /* General Options (Section 6.4.1) */
  setIntParameter("ML output",0,"Output Level",PL,intParam);
  setIntParameter("print unused",-2,"Print unused parameters",PL,intParam);
  setIntParameter("ML print initial list",-2,"Print initial list supplied to constructor",PL,intParam);
  setIntParameter("ML print final list",-2,"Print final list used by constructor",PL,intParam);
  setIntParameter("PDE equations",1,"# of PDE equations per node",PL,intParam);
  setStringToIntegralParameter<int>("eigen-analysis: type","cg","Scheme to compute spectral radius",
                               tuple<std::string>("cg","Anorm","power-method"),PL);
  setIntParameter("eigen-analysis: iterations",10,"# iterations of eigen-anaysis",PL,intParam);
  PL->set("ML label","dummy std::string");
  setIntParameter("print hierarchy",-2,"Print hierarchy.  0 or greater prints individual levels.",PL,intParam);

  /* Multigrid Cycle Options (Section 6.4.2) */
  setIntParameter("cycle applications",1,"# MG cycles",PL,intParam);
  setIntParameter("max levels",10,"Max # of levels",PL,intParam);
  setStringToIntegralParameter<int>("increasing or decreasing", "increasing", "Level numbering",tuple<std::string>("increasing","decreasing"),PL);
  setStringToIntegralParameter<int>("prec type", "MGV","Multigrid cycle type",tuple<std::string>("MGV","MGW","full-MGV","one-level-postsmoothing","two-level-additive","two-level-hybrid","two-level-hybrid2","projected MGV"),PL);
  PL->set("projected modes",(double**)0);
  setIntParameter("number of projected modes",0,"# of modes to be projected out before and after the V-cycle",PL,intParam);

  /* Aggregation and Prolongator Options (Section 6.4.3) */
  SetValidAggrParams(PL);
  for (int i = 0; i < 10; ++i) {
    char param[32];
    sprintf(param,"aggregation: list (level %d)",i);
    SetValidAggrParams(&(PL->sublist(param)));
  }

  PL->set("energy minimization: enable",false);
  setIntParameter("energy minimization: type",2,"Norm to use for energy minimization",PL,intParam);
  setDoubleParameter("energy minimization: droptol",0.0,"Drop tolerance for energy minimization",PL,dblParam);
  PL->set("energy minimization: cheap",false);
  setIntParameter("semicoarsen: number of levels",-1,"number of levels to apply semi-coarsening",PL,intParam);
  setIntParameter("semicoarsen: coarsen rate",-1,"Coarsening rate for structured semi-coarsening",PL,intParam);
  setIntParameter("semicoarsen: line direction nodes", -1,"Number of mesh points in z direction", PL, intParam);
  setStringToIntegralParameter<int>("semicoarsen: line orientation","use coordinates","indicates grid points are numbered either horizontally or veritcally for extruded meshes", tuple<std::string>("horizontal","vertical","use coordinates"),PL);


  /* Smoothing Options (Section 6.4.4) */
  SetValidSmooParams(PL,smoothers);
  for (int i = 0; i < 10; ++i) {
    char param[32];
    sprintf(param,"smoother: list (level %d)",i);
    SetValidSmooParams(&(PL->sublist(param)),smoothers);
  }
  SetValidSmooParams(&(PL->sublist("coarse: list")),smoothers);

  /* Load-balancing Options (Section 6.4.6) */
  setIntParameter("repartition: enable",0,"Enable repartitioning",PL,intParam);
  setStringToIntegralParameter<int>("repartition: partitioner","Zoltan","Repartitioning method",tuple<std::string>("Zoltan","ParMETIS"),PL);
  setDoubleParameter("repartition: max min ratio",1.3,"Specifies desired maximum imbalance ratio",PL,dblParam);
  setIntParameter("repartition: min per proc",512,"Specifies minimum # rows / processor",PL,intParam);
  setIntParameter("repartition: put on single proc",5000,"Specifies max global problem to be put on one processor",PL,intParam);
  setDoubleParameter("repartition: node max min ratio",1.3,"Specifies desired maximum imbalance for nodal heirarchy (Maxwell)",PL,dblParam);
  setIntParameter("repartition: node min per proc",170,"Specifies minimum number of nodes per proc (Maxwell)",PL,intParam);
  setIntParameter("repartition: Zoltan dimensions",0,"Dimension of problem",PL,intParam);
  setIntParameter("repartition: start level",1,"Suppress repartitioning until this level",PL,intParam);

  /* Analysis Options (Section 6.4.7) */
  PL->set("analyze memory",false);
  PL->set("viz: enable",false);
  setStringToIntegralParameter<int>("viz: output format","vtk","Visualization format",tuple<std::string>("vtk","xyz","openx"),PL);
  PL->set("viz: print starting solution",false);
  setIntParameter("viz: equation to plot",-1,"Equation number to print",PL,intParam);

  /* Miscellaneous Options (Section 6.4.8) */
  PL->set("x-coordinates",(double*)0);
  PL->set("y-coordinates",(double*)0);
  PL->set("z-coordinates",(double*)0);
  PL->set("node: x-coordinates",(double*)0);
  PL->set("node: y-coordinates",(double*)0);
  PL->set("node: z-coordinates",(double*)0);
  PL->set("read XML",true);
  PL->set("XML input file","ml_ParameterList.xml",std::string(""));

  /* Smoothed Aggregation and the Null Space (Section 6.4.9) */
  setStringToIntegralParameter<int>("null space: type","default vectors","Type of null space to use",tuple<std::string>("pre-computed","enriched","default vectors","elasticity from coordinates"),PL);
  PL->set("null space: vectors",(double*)0);
  setIntParameter("null space: dimension",0,"Number of user-supplied null space vectors",PL,intParam);
  setIntParameter("null space: vectors to compute",1,"Number of vectors to compute",PL,intParam);
  PL->set("null space: add default vectors",true);

  /* Aggregation Strategies (Section 6.4.10) */
  PL->set("aggregation: aux: enable",false);
  setDoubleParameter("aggregation: aux: threshold",0.0,"Dropping threshold for auxillary matrix",PL,dblParam);

  /* Unlisted Options */
  PL->set("ML debug mode",false);
  setStringToIntegralParameter<int>("default values","SA","Internal Option",tuple<std::string>("SA","DD","DD-ML","maxwell","NSSA","RefMaxwell","DD-ML-LU"),PL);
  PL->set("ML validate parameter list",true);
  setIntParameter("ML validate depth",0,"Internal option to control validation depth",PL,intParam);
  PL->set("ResetList",true);
  setStringToIntegralParameter<int>("SetDefaults","not-set","Internal Option",tuple<std::string>("not-set","SA","DD","DD-ML","maxwell","NSSA","RefMaxwell"),PL);
  setIntParameter("ML node id",-1,"Experimental option to identify the processor node (vis-a-vis core) id",PL,intParam);

  /* Unlisted Options that should probably be listed */
  setIntParameter("aggregation: aux: max levels",10,"Unlisted option",PL,intParam);
  PL->set("low memory usage",false);
  setDoubleParameter("aggregation: edge prolongator drop threshold",0.0,"Unlisted option",PL,dblParam);
  PL->set("zero starting solution",true);
  PL->set("aggregation: block scaling",false);
  setIntParameter("profile: operator iterations",0,"Unlisted option",PL,intParam);
  setDoubleParameter("subsmoother: edge alpha",20.0,"alpha for edge Chebyshev polynomial in Hiptmair",PL,dblParam);
  setDoubleParameter("subsmoother: node alpha",20.0,"alpha for node Chebyshev polynomial in Hiptmair",PL,dblParam);
  PL->set("reuse: enable",false);

  /* Unlisted options that should probably go away */
  setIntParameter("output",0,"Output Level",PL,intParam);

  /* Hightly experimental */
  PL->set("repartition: output timings",false);
  setIntParameter("repartition: estimated iterations",0,"Estimated number of iterations",PL,intParam);
  setStringToIntegralParameter<int>("repartition: Zoltan type","RCB","Type of repartitioner to use",tuple<std::string>("RCB","hypergraph","fast hypergraph"),PL);

  /* EXPERIMENTAL - Half-GS Smoothing */
  PL->set("smoother: Gauss-Seidel efficient symmetric",false);

  /* Coarse IFPACK Solvers - experimental */
  PL->set("coarse: ifpack list",dummy);
  PL->set("coarse: ifpack type",std::string(""));
  setIntParameter("coarse: ifpack overlap",0,"Unlisted option",PL,intParam);
  setDoubleParameter("coarse: ifpack level-of-fill",0.0,"Unlisted option",PL,dblParam);
  setDoubleParameter("coarse: ifpack relative threshold",1.0,"Unlisted option",PL,dblParam);
  setDoubleParameter("coarse: ifpack absolute threshold",0.0,"Unlisted option",PL,dblParam);

  /* EXPERIMENTAL - RefMaxwell block parallelization */
  PL->set("partitioner: options",dummy);
  PL->sublist("partitioner: options").disableRecursiveValidation();

  /* EXPERIMENTAL - node aware code */
  setIntParameter("ML node id", -1, "Unlisted option", PL, intParam);

  return PL;
}


/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/

bool ML_Epetra::ValidateRefMaxwellParameters(const Teuchos::ParameterList &inList){
  using Teuchos::ParameterList;
  using Teuchos::Exceptions::InvalidParameterName;
  using Teuchos::Exceptions::InvalidParameterType;
  using Teuchos::Exceptions::InvalidParameterValue;
  using std::cout;
  using std::endl;
  using std::string;

  ParameterList List,*validList;
  bool rv=true;

  /* Build a list with level-specific stuff stripped */
  //TODO this should be fixed
  for (ParameterList::ConstIterator param = inList.begin(); param != inList.end(); ++param) {
    const std::string pname=inList.name(param);
    if (pname.find("(level",0) == std::string::npos) {
      List.setEntry(pname,inList.entry(param));
    }
  }

  List.setName(inList.name());

  /* Get Defaults + Validate */
  try{
    validList = GetValidRefMaxwellParameters();
  }
  catch(...) {
    std::cout << "Error in GetValidMLPParameters.  Please report this bug to the ML "
      "developers." << std::endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_FAILURE);
  }
  try {
    List.validateParameters(*validList, 0, Teuchos::VALIDATE_USED_DISABLED,
			    Teuchos::VALIDATE_DEFAULTS_DISABLED);
  }
  catch(InvalidParameterName &excpt)  {rv=false; std::cout<<excpt.what();}
  catch(InvalidParameterType &excpt)  {rv=false; std::cout<<excpt.what();}
  catch(InvalidParameterValue &excpt) {rv=false; std::cout<<excpt.what();}
  catch(...) {rv=false;}
  delete validList;
  return rv;
}


Teuchos::ParameterList * ML_Epetra::GetValidRefMaxwellParameters(){
  using Teuchos::ParameterList;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;
  using std::string;

  ParameterList dummy,List11;
  ParameterList * PL = GetValidMLPParameters();

  /* RefMaxwell Options */
  setStringToIntegralParameter<int>("refmaxwell: 11solver","edge matrix free","(1,1) Block Solver",tuple<std::string>("edge matrix free"),PL);
  setStringToIntegralParameter<int>("refmaxwell: 22solver","multilevel","(2,2) Block Solver",tuple<std::string>("multilevel"),PL);
  setStringToIntegralParameter<int>("refmaxwell: mode","additive","Mode for RefMaxwell",tuple<std::string>("additive","212","121"),PL);
  PL->set("edge matrix free: coarse",dummy);
  List11.set("aggregation: aux: user matrix",(Epetra_CrsMatrix*)0);
  PL->set("refmaxwell: 11list",List11);
  PL->set("refmaxwell: 22list",dummy);


  // HAQ - clobber the smoother list to avoid validation.  MUST FIX LATER
  PL->set("smoother: type","IFPACK");


  /* RefMaxwell Options - Unsupported */
  PL->set("refmaxwell: aggregate with sigma",false);
  PL->set("refmaxwell: lump m1",false);
  PL->set("refmaxwell: disable addon",false);
  PL->set("refmaxwell: normalize prolongator",false);
  PL->set("refmaxwell: parallelize blocks",false);
  PL->set("refmaxwell: local nodal list",dummy);
  PL->set("refmaxwell: enable local nodal solver",false);
  PL->set("refmaxwell: global to local nodal transfer matrix",(Epetra_CrsMatrix*)0);
  PL->set("refmaxwell: drop nodal correction",false);//HAQ
  return PL;
}


#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
