
#include "ml_common.h"
#include "ml_common.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

#include "ml_ValidateParameters.h"
#include "Teuchos_StandardParameterEntryValidators.hpp"

using namespace Teuchos;
using namespace ML_Epetra;


bool ML_Epetra::ValidateMLPParameters(const Teuchos::ParameterList &inList){
  Teuchos::ParameterList List,*validList;
  bool rv=true;
    
  /* Build a list with level-specific stuff stripped */
  for(ParameterList::ConstIterator param=inList.begin(); param!=inList.end(); param++){
    const std::string pname=inList.name(param);
    if(pname.find("(level",0) == string::npos &&
       pname.find("user-defined function",0) == string::npos)
      List.setEntry(pname,inList.entry(param));
  }

  /* Get Defaults + Validate */
  try{
  validList=GetValidMLPParameters();
  }
  catch(...) {
    std::cout<<"Error in GetValidMLPParameters: The developers messed something up.  Sorry."<<std::endl;
#   ifdef HAVE_MPI
    MPI_Finalize();
#   endif
    exit(EXIT_FAILURE);
  }
  try{
    List.validateParameters(*validList,0,VALIDATE_USED_DISABLED,VALIDATE_DEFAULTS_DISABLED);
  }
  catch(Exceptions::InvalidParameterName &excpt)  {rv=false; std::cout<<excpt.what()<<std::endl;}
  catch(Exceptions::InvalidParameterType &excpt)  {rv=false; std::cout<<excpt.what()<<std::endl;}
  catch(Exceptions::InvalidParameterValue &excpt) {rv=false; std::cout<<excpt.what()<<std::endl;}
  catch(...) {rv=false;}
  delete validList;
  return rv;
}


Teuchos::ParameterList * ML_Epetra::GetValidMLPParameters(){
  Teuchos::ParameterList dummy;
  ParameterList * PL = new ParameterList;

  // prevent Teuchos from converting parameter types
  AnyNumberParameterEntryValidator::AcceptedTypes intParam(false),
                                                  dblParam(false),
                                                  strParam(false);
  intParam.allowInt(true); 
  dblParam.allowDouble(true);
  strParam.allowString(true); 

  /* Allocate List for Smoothing Options */
  const int num_smoothers=17;
  const char* smoother_strings[num_smoothers]={"Aztec","IFPACK","Jacobi","ML symmetric Gauss-Seidel","symmetric Gauss-Seidel","ML Gauss-Seidel","Gauss-Seidel","Chebyshev","MLS","Hiptmair","Amesos-KLU","Amesos-Superlu","Amesos-UMFPACK","Amesos-Superludist","Amesos-MUMPS","user-defined","SuperLU"};
  Array<std::string> smoothers(num_smoothers);
  for(int i=0;i<num_smoothers;i++) smoothers[i] = smoother_strings[i];

  /* General Options (Section 6.4.1) */
  setIntParameter("ML output",0,"Output Level",PL,intParam);
  setIntParameter("print unused",-2,"Print unused parameters",PL,intParam);
  setIntParameter("ML print initial list",-2,"Print initial list supplied to constructor",PL,intParam);
  setIntParameter("ML print final list",-2,"Print final list used by constructor",PL,intParam);
  setIntParameter("PDE equations",1,"# of PDE equations per node",PL,intParam);
  setStringToIntegralParameter<int>("eigen-analysis: type","cg","Scheme to compute spectral radius",
                               tuple<std::string>("cg","Anorm","power-method"),PL);
  setIntParameter("eigen-analysis: iterations",10,"# iterations of eigen-anaysis",PL,intParam);

  /* Multigrid Cycle Options (Section 6.4.2) */
  setIntParameter("cycle applications",1,"# MG cycles",PL,intParam);
  setIntParameter("max levels",10,"Max # of levels",PL,intParam);
  setStringToIntegralParameter<int>("increasing or decreasing", "increasing", "Level numbering",tuple<std::string>("increasing","decreasing"),PL);
  setStringToIntegralParameter<int>("prec type", "MGV","Multigrid cycle type",tuple<std::string>("MGV","MGW","full-MGV","one-level-postsmoothing","two-level-additive","two-level-hybrid","two-level-hybrid2","projected MGV"),PL);
  PL->set("projected mode",(double**)0);
  setIntParameter("number of projected modes",0,"# of modes to be projected out before and after the V-cycle",PL,intParam);

  /* Aggregation and Prolongator Options (Section 6.4.3) */
  setStringToIntegralParameter<int>("aggregation: type", "Uncoupled", "Aggregation algorithm",tuple<std::string>("Uncoupled","Coupled","MIS","Uncoupled-MIS","METIS","ParMETIS","Zoltan","user"),PL);
  setDoubleParameter("aggregation: threshold",0.0,"Dropping for aggregation",PL,dblParam);
  setDoubleParameter("aggregation: damping factor",1.3333,"Damping factor for smoothed aggregation",PL,dblParam);
  setIntParameter("aggregation: smoothing sweeps",1,"Number of sweeps for prolongator smoother",PL,intParam);
  setIntParameter("aggregation: global aggregates",0,"Number of global aggregates (METIS/ParMETIS)",PL,intParam);
  setIntParameter("aggregation: local aggregates",0,"Number of local aggregates (METIS/ParMETIS)",PL,intParam);
  setIntParameter("aggregation: nodes per aggregate",0,"Number of nodes per aggregate (METIS/ParMETIS)",PL,intParam);
  setIntParameter("aggregation: next-level aggregates per process",128,"Number of next-level rows / process (ParMETIS)",PL,intParam);
  PL->set("aggregation: use tentative restriction",false);
  PL->set("aggregation: symmetrize",false);
  PL->set("energy minimization: enable",false);
  setIntParameter("energy minimization: type",2,"Norm to use for energy minimization",PL,intParam);  
  setDoubleParameter("energy minimization: droptol",0.0,"Drop tolerance for energy minimization",PL,dblParam);
  PL->set("energy minimization: cheap",false);

  /* Smoothing Options (Section 6.4.4) */
  setStringToIntegralParameter<int>("smoother: type",std::string("Chebyshev"),"Smoothing algorithm",smoothers,PL);
  setIntParameter("smoother: sweeps",2,"Number of smoothing sweeps",PL,intParam);
  setDoubleParameter("smoother: damping factor",1.0,"Smoother damping factor",PL,dblParam);
  setStringToIntegralParameter<int>("smoother: pre or post","both","Smooth before/after coarse correction, or both",tuple<std::string>("pre","post","both"),PL);

  RCP<std::vector<int> > options = rcp(new std::vector<int>(AZ_OPTIONS_SIZE));
  RCP<std::vector<double> > params = rcp(new std::vector<double>(AZ_PARAMS_SIZE));
  PL->set("smoother: Aztec options",options);
  PL->set("smoother: Aztec params",params);
  PL->set("smoother: Aztec as solver",false);
  setDoubleParameter("smoother: Chebyshev alpha",20.0,"Damping radius for Chebyshev",PL,dblParam);
  setDoubleParameter("smoother: MLS alpha",20.0,"Damping radius for Chebyshev",PL,dblParam);
  PL->set("smoother: user-defined function",(int(*)(ML_Smoother*,int,double*,int,double*))0);
  PL->set("smoother: user-defined name",std::string("User-defined"));
  PL->set("smoother: Hiptmair efficient symmetric",true); 
  setStringToIntegralParameter<int>("subsmoother: type","Chebyshev","Subsmoother algorithm (Maxwell)",tuple<std::string>("Chebyshev","symmetric Gauss-Seidel","MLS"),PL);
  setDoubleParameter("subsmoother: Chebyshev alpha",20.0,"Damping radius for Chebshev",PL,dblParam);
  setDoubleParameter("subsmoother: MLS alpha",20.0,"Damping radius for Chebshev",PL,dblParam);
  setDoubleParameter("subsmoother: SGS damping factor",1.333,"Damping factor for symmetric Gauss-Seidel",PL,dblParam); 
  setIntParameter("subsmoother: edge sweeps",4,"Number of edge smoothing sweeps",PL,intParam);
  setIntParameter("subsmoother: node sweeps",4,"Number of node smoothing sweeps",PL,intParam);   
  
  /* Coarsest Grid Options (Section 6.4.5) */
  setIntParameter("coarse: max size",75,"Max coarse grid size",PL,intParam);
  setStringToIntegralParameter<int>("coarse: type","Chebyshev","Coarse solver algorithm",smoothers,PL);
  setStringToIntegralParameter<int>("coarse: pre or post","post","When to smooth on coarse grid",tuple<std::string>("pre","post","both"),PL);
  setIntParameter("coarse: sweeps",2,"Number of coarse sweeps",PL,intParam);  
  PL->set("coarse: user-defined function",(int(*)(ML_Smoother*,int,double*,int,double*))0);
  PL->set("coarse: user-defined name",std::string("User-defined"));
  setDoubleParameter("coarse: damping factor",1.0,"Coarse smoother damping factor",PL,dblParam); 
  setStringToIntegralParameter<int>("coarse: subsmoother type","Chebyshev","Coarse grid subsmoother (Maxwell)",tuple<std::string>("Chebyshev","symmetric Gauss-Seidel","MLS"),PL);
  setIntParameter("coarse: edge sweeps",4,"Number of coarse edge smoothing sweeps",PL,intParam);
  setIntParameter("coarse: node sweeps",4,"Number of coarse node smoothing sweeps",PL,intParam);
  setDoubleParameter("coarse: Chebyshev alpha",2.0,"Damping radius for Chebshev",PL,dblParam);
  setDoubleParameter("coarse: MLS alpha",2.0,"Damping radius for Chebshev",PL,dblParam);
  setIntParameter("coarse: max processes",-1,"Maximum number of procs for coarse solve (Superludist/MUMPS)",PL,intParam);  

  /* Load-balancing Options (Section 6.4.6) */
  setIntParameter("repartition: enable",0,"Enable repartitioning",PL,intParam);
  setStringToIntegralParameter<int>("repartition: partitioner","Zoltan","Repartitioning method",tuple<std::string>("Zoltan","ParMETIS"),PL);
  setDoubleParameter("repartition: max min ratio",1.3,"Specifies desired maximum imbalance ratio",PL,dblParam);
  setIntParameter("repartition: min per proc",512,"Specifies minimum # rows / processor",PL,intParam);
  setDoubleParameter("repartition: node max min ratio",1.3,"Specifies desired maximum imbalance for nodal heirarchy (Maxwell)",PL,dblParam);
  setIntParameter("repartition: node min per proc",170,"Specifies minimum number of nodes per proc (Maxwell)",PL,intParam);
  setIntParameter("repartition: Zoltan dimensions",0,"Dimension of problem",PL,intParam);

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
  setStringToIntegralParameter<int>("null space: type","default vectors","Type of null space to use",tuple<std::string>("pre-computed","enriched","default vectors"),PL);
  PL->set("null space: vectors",(double*)0); 
  setIntParameter("null space: dimension",0,"Number of user-supplied null space vectors",PL,intParam);
  setIntParameter("null space: vectors to compute",1,"Number of vectors to compute",PL,intParam);
  PL->set("null space: add default vectors",true);

  /* Aggregation Strategies (Section 6.4.10) */
  PL->set("aggregation: aux: enable",false);
  setDoubleParameter("aggregation: aux: threshold",0.0,"Dropping threshold for auxillary matrix",PL,dblParam);  

  /* Unlisted Options */ 
  PL->set("ML debug mode",false);
  setStringToIntegralParameter<int>("default values","SA","Internal Option",tuple<std::string>("SA","DD","DD-ML","maxwell","NSSA","RefMaxwell"),PL);
  PL->set("ML validate parameter list",true);
  PL->set("ResetList",true); 
  setStringToIntegralParameter<int>("SetDefaults","not-set","Internal Option",tuple<std::string>("not-set","SA","DD","DD-ML","maxwell","NSSA","RefMaxwell"),PL);

  /* Unlisted Options that should probably be listed */
  setIntParameter("aggregation: aux: max levels",10,"Unlisted option",PL,intParam);
  PL->set("low memory usage",false);
  setDoubleParameter("aggregation: edge prolongator drop threshold",0.0,"Unlisted option",PL,dblParam);
  PL->set("zero starting solution",true);
  PL->set("print hierarchy",false);  
  PL->set("smoother: self list",dummy);
  PL->set("aggregation: block scaling",false);
  setIntParameter("profile: operator iterations",0,"Unlisted option",PL,intParam);
  
  // From ml_Multilevel_Smoothers.cpp:
  setIntParameter("smoother: ParaSails matrix",0,"Unlisted option",PL,intParam);
  setIntParameter("smoother: ParaSails levels",0,"Unlisted option",PL,intParam);
  setDoubleParameter("smoother: ParaSails threshold",0.01,"Unlisted option",PL,dblParam);
  setDoubleParameter("smoother: ParaSails filter",0.05,"Unlisted option",PL,dblParam);
  setDoubleParameter("smoother: ParaSails load balancing",0,"Unlisted option",PL,dblParam);
  setIntParameter("smoother: ParaSails factorized",0,"Unlisted option",PL,intParam);
  PL->set("smoother: ifpack list",dummy); 
  PL->set("smoother: ifpack type",std::string(""));
  setIntParameter("smoother: ifpack overlap",0,"Unlisted option",PL,intParam);
  setDoubleParameter("smoother: ifpack level-of-fill",0.0,"Unlisted option",PL,dblParam);
  setDoubleParameter("smoother: ifpack relative threshold",1.0,"Unlisted option",PL,dblParam);
  setDoubleParameter("smoother: ifpack absolute threshold",0.0,"Unlisted option",PL,dblParam);
  
  /* Unlisted options that should probably go away */
  setIntParameter("output",0,"Output Level",PL,intParam);
  setIntParameter("smoother: polynomial order",2,"Unlisted option",PL,intParam);
  setIntParameter("smoother: MLS polynomial order",2,"Unlisted option",PL,intParam);  
  setIntParameter("coarse: polynomial order",2,"Unlisted option",PL,intParam);
  setIntParameter("coarse: MLS polynomial order",2,"Unlisted option",PL,intParam);

  /* Hightly experimental */
  PL->set("aggregation: respect materials",false);
  PL->set("aggregation: material type",(int*)0); 



  return PL;
}


/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/

bool ML_Epetra::ValidateRefMaxwellParameters(const Teuchos::ParameterList &inList){
  Teuchos::ParameterList List,*validList;
  bool rv=true;
  
  /* Build a list with level-specific stuff stripped */
  for(ParameterList::ConstIterator param=inList.begin(); param!=inList.end(); param++){
    const std::string pname=inList.name(param);
    if(pname.find("(level",0) == std::string::npos)
      List.setEntry(pname,inList.entry(param));
  }

  /* Get Defaults + Validate */
  try{
  validList=GetValidRefMaxwellParameters();
  }
  catch(...) {
    std::cout<<"Error in GetValidMLPParameters: The developers messed something up.  Sorry."<<std::endl;
#   ifdef HAVE_MPI
    MPI_Finalize();
#   endif
    exit(EXIT_FAILURE);
  }
  try{
    List.validateParameters(*validList,0,VALIDATE_USED_DISABLED,VALIDATE_DEFAULTS_DISABLED);
  }
  catch(Exceptions::InvalidParameterName &excpt)  {rv=false; std::cout<<excpt.what()<<std::endl;}
  catch(Exceptions::InvalidParameterType &excpt)  {rv=false; std::cout<<excpt.what()<<std::endl;}
  catch(Exceptions::InvalidParameterValue &excpt) {rv=false; std::cout<<excpt.what()<<std::endl;}
  catch(...) {rv=false;}
  delete validList;
  return rv;
}


Teuchos::ParameterList * ML_Epetra::GetValidRefMaxwellParameters(){
  Teuchos::ParameterList dummy;
  ParameterList * PL = GetValidMLPParameters();

  /* RefMaxwell Options - This should get added to the manual */
  setStringToIntegralParameter<int>("refmaxwell: 11solver","edge matrix free","(1,1) Block Solver",tuple<std::string>("edge matrix free"),PL);
  setStringToIntegralParameter<int>("refmaxwell: 22solver","multilevel","(2,2) Block Solver",tuple<std::string>("multilevel"),PL);
  setStringToIntegralParameter<int>("refmaxwell: mode","additive","Mode for RefMaxwell",tuple<std::string>("additive","212","121"),PL);
  PL->set("refmaxwell: 11list",dummy);
  PL->set("refmaxwell: 22list",dummy);

  /* RefMaxwell Options - Unsupported */
  PL->set("refmaxwell: aggregate with sigma",false);
  PL->set("refmaxwell: lump m1",false);
  PL->set("refmaxwell: disable addon",false); 
  PL->set("refmaxwell: normalize prolongator",false);
  return PL;
}


#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
