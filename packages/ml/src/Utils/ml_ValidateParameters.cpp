/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#include "ml_common.h"
#include "ml_common.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

#include "Epetra_CrsMatrix.h"
#include "ml_ValidateParameters.h"
#include "Teuchos_StandardParameterEntryValidators.hpp"

using namespace Teuchos;
using namespace ML_Epetra;
using namespace std;


bool ML_Epetra::ValidateMLPParameters(const Teuchos::ParameterList &inList,int depth){
  Teuchos::ParameterList List,*validList;
  bool rv=true;

  /* Build a copy of the list to be validated. */
  for(ParameterList::ConstIterator param=inList.begin(); param!=inList.end(); param++){
    const string pname=inList.name(param);
    if( pname.find("user-defined function",0) == string::npos )
    {
      List.setEntry(pname,inList.entry(param));
    }
  }

  List.setName(inList.name());

  /* Get Defaults + Validate */
  try{
  validList=GetValidMLPParameters();
  }
  catch(...) {
    cout<<"Error in GetValidMLPParameters: The developers messed something up.  Sorry."<<endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_FAILURE);
  }
  try{
    List.validateParameters(*validList,depth,VALIDATE_USED_ENABLED,VALIDATE_DEFAULTS_DISABLED);
  }
  catch(Exceptions::InvalidParameterName &excpt)  {rv=false; cout<<excpt.what()<<endl;}
  catch(Exceptions::InvalidParameterType &excpt)  {rv=false; cout<<excpt.what()<<endl;}
  catch(Exceptions::InvalidParameterValue &excpt) {rv=false; cout<<excpt.what()<<endl;}
  catch(...) {rv=false;}
  delete validList;
  return rv;
}

/******************************************************************************/

void ML_Epetra::SetValidSmooParams(ParameterList *PL, Array<string> &smoothers)
{
  Teuchos::ParameterList dummy;
  // prevent Teuchos from converting parameter types
  AnyNumberParameterEntryValidator::AcceptedTypes intParam(false),
                                                  dblParam(false),
                                                  strParam(false);
  intParam.allowInt(true); 
  dblParam.allowDouble(true);
  strParam.allowString(true); 

  /* Smoothing Options (Section 6.4.4) */
  setStringToIntegralParameter<int>("smoother: type",string("Chebyshev"),"Smoothing algorithm",smoothers,PL);
  setIntParameter("smoother: sweeps",2,"Number of smoothing sweeps",PL,intParam);
  setDoubleParameter("smoother: damping factor",1.0,"Smoother damping factor",PL,dblParam);
  setStringToIntegralParameter<int>("smoother: pre or post","both","Smooth before/after coarse correction, or both",tuple<string>("pre","post","both"),PL);
#ifdef HAVE_ML_AZTECOO
  RCP<vector<int> > options = rcp(new vector<int>(AZ_OPTIONS_SIZE));
  RCP<vector<double> > params = rcp(new vector<double>(AZ_PARAMS_SIZE));
  PL->set("smoother: Aztec options",options);
  PL->set("smoother: Aztec params",params);
#endif
  PL->set("smoother: Aztec as solver",false);
  setDoubleParameter("smoother: Chebyshev alpha",20.0,"Damping radius for Chebyshev",PL,dblParam);
  setDoubleParameter("smoother: MLS alpha",20.0,"Damping radius for Chebyshev",PL,dblParam);
  PL->set("smoother: user-defined function",(int(*)(ML_Smoother*,int,double*,int,double*))0);
  PL->set("smoother: user-defined name",string("User-defined"));
  PL->set("smoother: Hiptmair efficient symmetric",true); 
  setStringToIntegralParameter<int>("subsmoother: type","Chebyshev","Subsmoother algorithm (Maxwell)",tuple<string>("Chebyshev","symmetric Gauss-Seidel","MLS"),PL);
  setDoubleParameter("subsmoother: Chebyshev alpha",20.0,"Damping radius for Chebshev",PL,dblParam);
  setDoubleParameter("subsmoother: MLS alpha",20.0,"Damping radius for Chebshev",PL,dblParam);
  setDoubleParameter("subsmoother: damping factor",1.333,"Damping factor for symmetric Gauss-Seidel",PL,dblParam); 
  setIntParameter("subsmoother: edge sweeps",4,"Number of edge smoothing sweeps",PL,intParam);
  setIntParameter("subsmoother: node sweeps",4,"Number of node smoothing sweeps",PL,intParam);   
# ifdef HAVE_PETSC
  void *petscKSP;
  PL->set("smoother: petsc ksp",petscKSP);
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
  PL->set("smoother: ifpack type",string(""));
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
  //setStringToIntegralParameter<int>("smoother: pre or post","post","When to smooth on coarse grid",tuple<string>("pre","post","both"),PL);
  //setIntParameter("smoother: sweeps",2,"Number of coarse sweeps",PL,intParam);  
  //PL->set("smoother: user-defined function",(int(*)(ML_Smoother*,int,double*,int,double*))0);
  //PL->set("smoother: user-defined name",string("User-defined"));
  //setDoubleParameter("smoother: damping factor",1.0,"Coarse smoother damping factor",PL,dblParam); 
  setStringToIntegralParameter<int>("smoother: subsmoother type","Chebyshev","Coarse grid subsmoother (Maxwell)",tuple<string>("Chebyshev","symmetric Gauss-Seidel","MLS"),PL);
  setIntParameter("smoother: edge sweeps",4,"Number of coarse edge smoothing sweeps",PL,intParam);
  setIntParameter("smoother: node sweeps",4,"Number of coarse node smoothing sweeps",PL,intParam);
  //setDoubleParameter("smoother: Chebyshev alpha",2.0,"Damping radius for Chebshev",PL,dblParam);
  //setDoubleParameter("smoother: MLS alpha",2.0,"Damping radius for Chebshev",PL,dblParam);
  setIntParameter("smoother: max processes",-1,"Maximum number of procs for coarse solve (Superludist/MUMPS)",PL,intParam);  


  /* EXPERIMENTAL - Half-GS Smoothing */
  PL->set("smoother: Gauss-Seidel efficient symmetric",false);
  setIntParameter("smoother: Block Chebyshev number of blocks",-1,"Number of blocks to use with Block Chebyshev",PL,intParam);    
  PL->set("smoother: Block Chebyshev block starts",(int*)0);
  PL->set("smoother: Block Chebyshev block list",(int*)0);
  
  /* Coarse IFPACK Solvers - experimental */
  //PL->set("smoother: ifpack list",dummy);
  //PL->set("smoother: ifpack type",std::string(""));
  //setIntParameter("smoother: ifpack overlap",0,"Unlisted option",PL,intParam);
  //setDoubleParameter("smoother: ifpack level-of-fill",0.0,"Unlisted option",PL,dblParam);
  //setDoubleParameter("smoother: ifpack relative threshold",1.0,"Unlisted option",PL,dblParam);
  //setDoubleParameter("smoother: ifpack absolute threshold",0.0,"Unlisted option",PL,dblParam);

} //SetValidSmooParams()

/******************************************************************************/

void ML_Epetra::SetValidAggrParams(ParameterList *PL)
{
  // prevent Teuchos from converting parameter types
  AnyNumberParameterEntryValidator::AcceptedTypes intParam(false),
                                                  dblParam(false),
                                                  strParam(false);
  intParam.allowInt(true); 
  dblParam.allowDouble(true);
  strParam.allowString(true); 

  /* Aggregation and Prolongator Options (Section 6.4.3) */
  setStringToIntegralParameter<int>("aggregation: type", "Uncoupled", "Aggregation algorithm", tuple<string>("Uncoupled","Coupled","MIS","Uncoupled-MIS","METIS","ParMETIS","Zoltan","user"),PL);
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
# ifdef HAVE_PETSC
  const int num_smoothers=29;
# else
  const int num_smoothers=28;
# endif
  const char* smoother_strings[num_smoothers]={"Aztec","IFPACK","Jacobi",
   "ML symmetric Gauss-Seidel","symmetric Gauss-Seidel","ML Gauss-Seidel",
   "Gauss-Seidel","block Gauss-Seidel","symmetric block Gauss-Seidel",
   "Chebyshev","MLS","Hiptmair","Amesos-KLU","Amesos-Superlu",
   "Amesos-UMFPACK","Amesos-Superludist","Amesos-MUMPS","user-defined",
   "SuperLU","IFPACK-Chebyshev","self","do-nothing","IC","ICT","ILU","ILUT",
   "Block Chebyshev","IFPACK-Block Chebyshev"
#  ifdef HAVE_PETSC
   ,"petsc"
#  endif
   };
  Array<string> smoothers(num_smoothers);
  for(int i=0;i<num_smoothers;i++) smoothers[i] = smoother_strings[i];

  /* General Options (Section 6.4.1) */
  setIntParameter("ML output",0,"Output Level",PL,intParam);
  setIntParameter("print unused",-2,"Print unused parameters",PL,intParam);
  setIntParameter("ML print initial list",-2,"Print initial list supplied to constructor",PL,intParam);
  setIntParameter("ML print final list",-2,"Print final list used by constructor",PL,intParam);
  setIntParameter("PDE equations",1,"# of PDE equations per node",PL,intParam);
  setStringToIntegralParameter<int>("eigen-analysis: type","cg","Scheme to compute spectral radius",
                               tuple<string>("cg","Anorm","power-method"),PL);
  setIntParameter("eigen-analysis: iterations",10,"# iterations of eigen-anaysis",PL,intParam);
  PL->set("ML label","dummy string");
  setIntParameter("print hierarchy",-2,"Print hierarchy.  0 or greater prints individual levels.",PL,intParam);

  /* Multigrid Cycle Options (Section 6.4.2) */
  setIntParameter("cycle applications",1,"# MG cycles",PL,intParam);
  setIntParameter("max levels",10,"Max # of levels",PL,intParam);
  setStringToIntegralParameter<int>("increasing or decreasing", "increasing", "Level numbering",tuple<string>("increasing","decreasing"),PL);
  setStringToIntegralParameter<int>("prec type", "MGV","Multigrid cycle type",tuple<string>("MGV","MGW","full-MGV","one-level-postsmoothing","two-level-additive","two-level-hybrid","two-level-hybrid2","projected MGV"),PL);
  PL->set("projected mode",(double**)0);
  setIntParameter("number of projected modes",0,"# of modes to be projected out before and after the V-cycle",PL,intParam);

  /* Aggregation and Prolongator Options (Section 6.4.3) */
  SetValidAggrParams(PL);
  for (int i=0; i<10;i++) {
    char param[32];
    sprintf(param,"aggregation: list (level %d)",i); 
    SetValidAggrParams(&(PL->sublist(param)));
  }
  
  PL->set("energy minimization: enable",false);
  setIntParameter("energy minimization: type",2,"Norm to use for energy minimization",PL,intParam);  
  setDoubleParameter("energy minimization: droptol",0.0,"Drop tolerance for energy minimization",PL,dblParam);
  PL->set("energy minimization: cheap",false);

  /* Smoothing Options (Section 6.4.4) */
  SetValidSmooParams(PL,smoothers);
  for (int i=0; i<10;i++) {
    char param[32];
    sprintf(param,"smoother: list (level %d)",i); 
    SetValidSmooParams(&(PL->sublist(param)),smoothers);
  }
  SetValidSmooParams(&(PL->sublist("coarse: list")),smoothers);

  /* Load-balancing Options (Section 6.4.6) */
  setIntParameter("repartition: enable",0,"Enable repartitioning",PL,intParam);
  setStringToIntegralParameter<int>("repartition: partitioner","Zoltan","Repartitioning method",tuple<string>("Zoltan","ParMETIS"),PL);
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
  setStringToIntegralParameter<int>("viz: output format","vtk","Visualization format",tuple<string>("vtk","xyz","openx"),PL);
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
  PL->set("XML input file","ml_ParameterList.xml",string(""));

  /* Smoothed Aggregation and the Null Space (Section 6.4.9) */
  setStringToIntegralParameter<int>("null space: type","default vectors","Type of null space to use",tuple<string>("pre-computed","enriched","default vectors","from coordinates"),PL);
  PL->set("null space: vectors",(double*)0); 
  setIntParameter("null space: dimension",0,"Number of user-supplied null space vectors",PL,intParam);
  setIntParameter("null space: vectors to compute",1,"Number of vectors to compute",PL,intParam);
  PL->set("null space: add default vectors",true);

  /* Aggregation Strategies (Section 6.4.10) */
  PL->set("aggregation: aux: enable",false);
  setDoubleParameter("aggregation: aux: threshold",0.0,"Dropping threshold for auxillary matrix",PL,dblParam);  

  /* Unlisted Options */ 
  PL->set("ML debug mode",false);
  setStringToIntegralParameter<int>("default values","SA","Internal Option",tuple<string>("SA","DD","DD-ML","maxwell","NSSA","RefMaxwell"),PL);
  PL->set("ML validate parameter list",true);
  setIntParameter("ML validate depth",0,"Internal option to control validation depth",PL,intParam);
  PL->set("ResetList",true); 
  setStringToIntegralParameter<int>("SetDefaults","not-set","Internal Option",tuple<string>("not-set","SA","DD","DD-ML","maxwell","NSSA","RefMaxwell"),PL);
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
  setStringToIntegralParameter<int>("repartition: Zoltan type","RCB","Type of repartitioner to use",tuple<string>("RCB","hypergraph","fast hypergraph"),PL);

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
  setIntParameter("ML node id",-1,"Unlisted option",PL,intParam);

  return PL;
}


/***************************************************************************************************/
/***************************************************************************************************/
/***************************************************************************************************/

bool ML_Epetra::ValidateRefMaxwellParameters(const Teuchos::ParameterList &inList){
  Teuchos::ParameterList List,*validList;
  bool rv=true;
  
  /* Build a list with level-specific stuff stripped */
  //TODO this should be fixed
  for(ParameterList::ConstIterator param=inList.begin(); param!=inList.end(); param++){
    const string pname=inList.name(param);
    if(pname.find("(level",0) == string::npos)
      List.setEntry(pname,inList.entry(param));
  }

  List.setName(inList.name());

  /* Get Defaults + Validate */
  try{
  validList=GetValidRefMaxwellParameters();
  }
  catch(...) {
    cout<<"Error in GetValidMLPParameters: The developers messed something up.  Sorry."<<endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_FAILURE);
  }
  try{
    List.validateParameters(*validList,0,VALIDATE_USED_DISABLED,VALIDATE_DEFAULTS_DISABLED);
  }
  catch(Exceptions::InvalidParameterName &excpt)  {rv=false; cout<<excpt.what();}
  catch(Exceptions::InvalidParameterType &excpt)  {rv=false; cout<<excpt.what();}
  catch(Exceptions::InvalidParameterValue &excpt) {rv=false; cout<<excpt.what();}
  catch(...) {rv=false;}
  delete validList;
  return rv;
}


Teuchos::ParameterList * ML_Epetra::GetValidRefMaxwellParameters(){
  Teuchos::ParameterList dummy;
  ParameterList * PL = GetValidMLPParameters();

  /* RefMaxwell Options */
  setStringToIntegralParameter<int>("refmaxwell: 11solver","edge matrix free","(1,1) Block Solver",tuple<string>("edge matrix free"),PL);
  setStringToIntegralParameter<int>("refmaxwell: 22solver","multilevel","(2,2) Block Solver",tuple<string>("multilevel"),PL);
  setStringToIntegralParameter<int>("refmaxwell: mode","additive","Mode for RefMaxwell",tuple<string>("additive","212","121"),PL);
  PL->set("edge matrix free: coarse",dummy);
  PL->set("refmaxwell: 11list",dummy);
  PL->set("refmaxwell: 22list",dummy);

  /* RefMaxwell Options - Unsupported */
  PL->set("refmaxwell: aggregate with sigma",false);
  PL->set("refmaxwell: lump m1",false);
  PL->set("refmaxwell: disable addon",false); 
  PL->set("refmaxwell: normalize prolongator",false);
  PL->set("refmaxwell: parallelize blocks",false);
  PL->set("refmaxwell: local nodal list",dummy);
  PL->set("refmaxwell: enable local nodal solver",false);
  PL->set("refmaxwell: global to local nodal transfer matrix",(Epetra_CrsMatrix*)0);  
  
  return PL;
}


#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/
