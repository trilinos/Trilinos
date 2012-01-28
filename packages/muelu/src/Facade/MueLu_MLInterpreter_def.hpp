#ifndef MUELU_MLINTERPRETER_DEF_HPP
#define MUELU_MLINTERPRETER_DEF_HPP

#include <Xpetra_Operator.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_MLInterpreter_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_FactoryManager.hpp"

#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_HierarchyHelpers.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"

// #include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
MLInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MLInterpreter() { }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
MLInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~MLInterpreter() { }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > MLInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(const Teuchos::ParameterList & params, const RCP<Operator> & A, const RCP<MultiVector> nsp) {
#include "MueLu_UseShortNames.hpp" // needed because some classes are not forward declared in _decl.hpp

  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  // Read in common parameters
  int    maxLevels             = 10;   // multigrid parameters
  int    verbosityLevel        = 10;   // verbosity level
  int    maxCoarseSize         = 50;
  int    nDofsPerNode          = 1;    // coalesce and drop parameters
  double agg_threshold         = 0.0;  // aggregation parameters
  double agg_damping           = 4/3;
  int    agg_smoothingsweeps   = 1;
  int    minPerAgg             = 3;    // optimal for 2d
  int    maxNbrAlreadySelected = 0;
  std::string agg_type         = "Uncoupled";
  bool   bEnergyMinimization = false;  // PGAMG

  if(params.isParameter("max levels"))                      maxLevels           = params.get<int>("max levels");
  if(params.isParameter("ML output"))                       verbosityLevel      = params.get<int>("ML output");
  if(params.isParameter("coarse: max size"))                maxCoarseSize       = params.get<int>("coarse: max size");
  if(params.isParameter("PDE equations"))                   nDofsPerNode        = params.get<int>("PDE equations");
  if(params.isParameter("aggregation: threshold"))          agg_threshold       = params.get<double>("aggregation: threshold");
  if(params.isParameter("aggregation: damping factor"))     agg_damping         = params.get<double>("aggregation: damping factor");
  if(params.isParameter("aggregation: smoothing sweeps"))   agg_smoothingsweeps = params.get<int>   ("aggregation: smoothing sweeps");
  if(params.isParameter("aggregation: type"))               agg_type            = params.get<std::string> ("aggregation: type");
  //if(params.isParameter("aggregation: nodes per aggregate")) minPerAgg        = params.get<int>("aggregation: nodes per aggregate");
  if(params.isParameter("energy minimization: enable"))     bEnergyMinimization = params.get<bool>("energy minimization: enable");

  // Translate verbosity parameter
  Teuchos::EVerbosityLevel eVerbLevel = Teuchos::VERB_NONE;
  if(verbosityLevel == 0)  eVerbLevel = Teuchos::VERB_NONE;
  if(verbosityLevel >  0 ) eVerbLevel = Teuchos::VERB_LOW;
  if(verbosityLevel >  4 ) eVerbLevel = Teuchos::VERB_MEDIUM;
  if(verbosityLevel >  7 ) eVerbLevel = Teuchos::VERB_HIGH;
  if(verbosityLevel >  9 ) eVerbLevel = Teuchos::VERB_EXTREME;

  TEUCHOS_TEST_FOR_EXCEPTION(agg_type != "Uncoupled", Exceptions::RuntimeError, "MueLu::MLInterpreter::Setup(): parameter \"aggregation: type\": only 'Uncoupled' aggregation is supported.");

  // Create MueLu factories
  RCP<NullspaceFactory>     nspFact = rcp(new NullspaceFactory());
  RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory(/*Teuchos::null,nspFact*/));
  //dropFact->SetVerbLevel(toMueLuVerbLevel(eVerbLevel));
  dropFact->SetFixedBlockSize(nDofsPerNode);

  RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory(dropFact));
  if(verbosityLevel > 3) { // TODO fix me: Setup is a static function: we cannot use GetOStream without an object...
  *out << "========================= Aggregate option summary  =========================" << std::endl;
  *out << "min Nodes per aggregate :               " << minPerAgg << std::endl;
  *out << "min # of root nbrs already aggregated : " << maxNbrAlreadySelected << std::endl;
  *out << "aggregate ordering :                    NATURAL" << std::endl;
  *out << "=============================================================================" << std::endl;
  }
  UCAggFact->SetMinNodesPerAggregate(minPerAgg); //TODO should increase if run anything other than 1D
  UCAggFact->SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
  UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
  UCAggFact->SetPhase3AggCreation(0.5);

  RCP<PFactory> PFact;
  RCP<RFactory> RFact;

  if (agg_damping == 0.0 && bEnergyMinimization == false) {
    // tentative prolongation operator (PA-AMG)
    PFact = rcp(new TentativePFactory(UCAggFact/*,nspFact*/));
    RFact = rcp( new TransPFactory(PFact) );
  } else if(agg_damping != 0.0 && bEnergyMinimization == false) {
    // smoothed aggregation (SA-AMG)
    RCP<PFactory>  PtentFact = rcp(new TentativePFactory(UCAggFact/*,nspFact*/));
    PFact  = rcp( new SaPFactory(PtentFact) );
    RCP<SaPFactory> SaPFact = Teuchos::rcp_dynamic_cast<SaPFactory>(PFact);
    SaPFact->SetDampingFactor(agg_damping);
    RFact  = rcp( new TransPFactory(PFact) );
  } else if(bEnergyMinimization == true) {
    // Petrov Galerkin PG-AMG smoothed aggregation (energy minimization in ML)
    RCP<PFactory> PtentFact = rcp(new TentativePFactory(UCAggFact/*,nspFact*/));
    PFact  = rcp( new PgPFactory(PtentFact) );
    RFact  = rcp( new GenericRFactory(PFact) );
  }

  RCP<RAPFactory> AcFact = rcp( new RAPFactory(PFact, RFact) );
  //AcFact->setVerbLevel(eVerbLevel); //(toMueLuVerbLevel(eVerbLevel)); TODO: check me

  RCP<SmootherFactory> coarsestSmooFact;
  coarsestSmooFact = GetCoarsestSolverFactory(params);

  ///////////////////////////////////////////////////

  // Fill hierarchy
  RCP<Hierarchy> hierarchy = Teuchos::rcp(new Hierarchy(A));
  hierarchy->SetDefaultVerbLevel(toMueLuVerbLevel(eVerbLevel));
  hierarchy->SetMaxCoarseSize(Teuchos::as<Xpetra::global_size_t>(maxCoarseSize));

  ///////////////////////////////////////////////////////////

  // Set fine level nullspace
  // use given fine level null space or extract pre-computed nullspace from ML parameter list
  if (nsp != Teuchos::null) {
    RCP<MueLu::Level> Finest = hierarchy->GetLevel();  // get finest level
    Finest->Set("Nullspace",nsp);                      // set user given null space
  } else if(params.isParameter("null space: type")) {
    std::string type = "";
    type = params.get<std::string>("null space: type");
    TEUCHOS_TEST_FOR_EXCEPTION(type != "pre-computed", Exceptions::RuntimeError, "MueLu::Interpreter: no valid nullspace (no pre-computed null space). error.");
    int dimns = -1;
    if(params.isParameter("null space: dimension")) dimns = params.get<int>("null space: dimension");
    TEUCHOS_TEST_FOR_EXCEPTION(dimns == -1, Exceptions::RuntimeError, "MueLu::Interpreter: no valid nullspace (nullspace dim = -1). error.");
    
    const RCP<const Map> rowMap = A->getRowMap();
    RCP<MultiVector> nspVector = MultiVectorFactory::Build(rowMap,dimns,true);
    double* nsdata = NULL;
    if(params.isParameter("null space: vectors")) nsdata = params.get<double*>("null space: vectors");
    TEUCHOS_TEST_FOR_EXCEPTION(nsdata == NULL, Exceptions::RuntimeError, "MueLu::Interpreter: no valid nullspace (nsdata = NULL). error.");
    
    for ( size_t i=0; i < Teuchos::as<size_t>(dimns); i++) {
      Teuchos::ArrayRCP<Scalar> nspVectori = nspVector->getDataNonConst(i);
      const size_t myLength = nspVector->getLocalLength();
      for(size_t j=0; j<myLength; j++) {
        nspVectori[j] = nsdata[i*myLength+j];
      }
      }
    
      RCP<MueLu::Level> Finest = hierarchy->GetLevel();  // get finest level
      Finest->Set("Nullspace",nspVector);                // set user given null space
      //Finest->setDefaultVerbLevel(eVerbLevel);
  }

  ////////////////////////////////////

  // Prepare factory managers

  bool bIsLastLevel = false;
  std::vector<RCP<FactoryManager> > vecManager(maxLevels);
  for(int i=0; i < maxLevels; i++) {
    RCP<SmootherFactory> SmooFactFine   = GetSmootherFactory(params, i);

    vecManager[i] = rcp(new FactoryManager());
    if(SmooFactFine != Teuchos::null)
    	vecManager[i]->SetFactory("Smoother" ,  SmooFactFine);    // Hierarchy.Setup uses TOPSmootherFactory, that only needs "Smoother"
    vecManager[i]->SetFactory("CoarseSolver", coarsestSmooFact);
    vecManager[i]->SetFactory("A", AcFact);                       // same RAP factory
    vecManager[i]->SetFactory("P", PFact);                        // same prolongator and restrictor factories
    vecManager[i]->SetFactory("R", RFact);                        // same prolongator and restrictor factories
    vecManager[i]->SetFactory("Nullspace", nspFact);              // use same nullspace factory throughout all multigrid levels
  }

  // use new Hierarchy::Setup routine
  if(maxLevels == 1) {
	  bIsLastLevel = hierarchy->Setup(0, Teuchos::null, vecManager[0].ptr(), Teuchos::null);
  }
  else
  {
	  bIsLastLevel = hierarchy->Setup(0, Teuchos::null, vecManager[0].ptr(), vecManager[1].ptr()); // true, false because first level
	  for(int i=1; i < maxLevels-1; i++) {
		if(bIsLastLevel == true) break;
		bIsLastLevel = hierarchy->Setup(i, vecManager[i-1].ptr(), vecManager[i].ptr(), vecManager[i+1].ptr());
	  }
	  if(bIsLastLevel == false) {
		bIsLastLevel = hierarchy->Setup(maxLevels-1, vecManager[maxLevels-2].ptr(), vecManager[maxLevels-1].ptr(), Teuchos::null);
	  }
  }
  //*out << *hierarchy << std::endl;


  //for(int i=0; i<hierarchy->GetNumLevels(); i++) {
    /*RCP<Level> l = hierarchy->GetLevel(i);
    l->print(*out, Teuchos::VERB_EXTREME);*/
    /*RCP<SmootherBase> preSmoo = l->Get< RCP<SmootherBase> >("PreSmoother");
    RCP<TrilinosSmoother> tPreSmoo = Teuchos::rcp_dynamic_cast<TrilinosSmoother>(preSmoo);
    if(tPreSmoo!=Teuchos::null) tPreSmoo->print(*out, MueLu::toMueLuVerbLevel(Teuchos::VERB_EXTREME));*/
  //}

  return hierarchy;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<MueLu::SmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > MLInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetCoarsestSolverFactory(const Teuchos::ParameterList & params) {
#include "MueLu_UseShortNames.hpp" // needed because some classes are not forward declared in _decl.hpp

#if   defined(HAVE_MUELU_AMESOS2)
  std::string type = "Amesos-Superlu"; // propose SuperLU as coarsest solver if AMESOS2 is enabled // TODO: check if SuperLU is available in Amesos2
#else
  std::string type = "Amesos-KLU";     // propose KLU as coarsest solver // TODO: only available with Epetra? 
  // TODO: what if neither KLU nor SuperLU are available??
#endif
  
  if(params.isParameter("coarse: type")) type = params.get<std::string>("coarse: type");

  RCP<SmootherPrototype> smooProto;
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  RCP<SmootherFactory> SmooFact;

  if(type == "Jacobi") {
    if(params.isParameter("coarse: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", params.get<int>("coarse: sweeps"));
    if(params.get<double>("coarse: damping factor"))
      ifpackList.set("relaxation: damping factor", params.get<double>("coarse: damping factor"));
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Jacobi");
    smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList) );
  } else if(type == "Gauss-Seidel") {
    if(params.isParameter("coarse: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", params.get<int>("coarse: sweeps"));
    if(params.get<double>("coarse: damping factor"))
      ifpackList.set("relaxation: damping factor", params.get<double>("coarse: damping factor"));
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Gauss-Seidel");
    smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList) );
  } else if (type == "symmetric Gauss-Seidel") {
    if(params.isParameter("coarse: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", params.get<int>("coarse: sweeps"));
    if(params.get<double>("coarse: damping factor"))
      ifpackList.set("relaxation: damping factor", params.get<double>("coarse: damping factor"));
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList) );
  } else if (type == "Chebyshev") {
    ifpackType = "CHEBYSHEV";
    if(params.isParameter("coarse: sweeps"))
      ifpackList.set("chebyshev: degree", params.get<int>("coarse: sweeps"));
    if(params.isParameter("coarse: Chebyshev alpha"))
      ifpackList.set("chebyshev: alpha", params.get<double>("coarse: Chebyshev alpha"));
    smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList) );
  } else if(type == "IFPACK") {
    // TODO change to TrilnosSmoother as soon as Ifpack2 supports all preconditioners from Ifpack
    ifpackType = params.get<std::string>("coarse: ifpack type");
    if(ifpackType == "ILU") {
      ifpackList.set<int>("fact: level-of-fill", (int)params.get<double>("coarse: ifpack level-of-fill"));
      ifpackList.set("partitioner: overlap", params.get<int>("coarse: ifpack overlap"));
      smooProto = rcp( new IfpackSmoother(ifpackType, ifpackList, params.get<int>("coarse: ifpack overlap")) );
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Interpreter: unknown ML smoother type " + type + " (IFPACK) not supported by MueLu. Only ILU is supported.");
  } else if(type == "Amesos-Superlu") {
    Teuchos::ParameterList coarsestSmooList;
    smooProto = Teuchos::rcp( new DirectSolver("Superlu", coarsestSmooList) );
  } else if(type == "Amesos-Superludist") {
    Teuchos::ParameterList coarsestSmooList;
    smooProto = Teuchos::rcp( new DirectSolver("Superludist", coarsestSmooList) );
  } else if(type == "Amesos-KLU") {
    Teuchos::ParameterList coarsestSmooList;
    smooProto = Teuchos::rcp( new DirectSolver("Klu", coarsestSmooList) );
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Interpreter: unknown coarsest solver type. not supported by MueLu.");
  }

  // create smoother factory
  TEUCHOS_TEST_FOR_EXCEPTION(smooProto == Teuchos::null, Exceptions::RuntimeError, "MueLu::Interpreter: no smoother prototype. fatal error.");
  SmooFact = rcp( new SmootherFactory(smooProto) );

  // check if pre- and postsmoothing is set
  std::string preorpost = "both";
  if(params.isParameter("coarse: pre or post")) preorpost = params.get<std::string>("coarse: pre or post");

  if (preorpost == "pre") {
    SmooFact->SetSmootherPrototypes(smooProto, Teuchos::null);
  } else if(preorpost == "post") {
    SmooFact->SetSmootherPrototypes(Teuchos::null, smooProto);
  }

  return SmooFact;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
RCP<MueLu::SmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > MLInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetSmootherFactory(const Teuchos::ParameterList & params, int level) {
#include "MueLu_UseShortNames.hpp" // needed because some classes are not forward declared in _decl.hpp
  char levelchar[11];
  sprintf(levelchar,"(level %d)",level);
  std::string levelstr(levelchar);

  if(params.isSublist("smoother: list " + levelstr)==false)
	  return Teuchos::null;
  //TEUCHOS_TEST_FOR_EXCEPTION(params.isSublist("smoother: list " + levelstr)==false, Exceptions::RuntimeError, "MueLu::Interpreter: no ML smoother parameter list for level. error.");

  std::string type = params.sublist("smoother: list " + levelstr).get<std::string>("smoother: type");
  TEUCHOS_TEST_FOR_EXCEPTION(type.empty(), Exceptions::RuntimeError, "MueLu::Interpreter: no ML smoother type for level. error.");

  const Teuchos::ParameterList smolevelsublist = params.sublist("smoother: list " + levelstr);
  RCP<SmootherPrototype> smooProto;
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  RCP<SmootherFactory> SmooFact;

  if(type == "Jacobi") {
    if(smolevelsublist.isParameter("smoother: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", smolevelsublist.get<int>("smoother: sweeps"));
    if(smolevelsublist.get<double>("smoother: damping factor"))
      ifpackList.set("relaxation: damping factor", smolevelsublist.get<double>("smoother: damping factor"));
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Jacobi");
    smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList) );
  } else if(type == "Gauss-Seidel") {
    if(smolevelsublist.isParameter("smoother: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", smolevelsublist.get<int>("smoother: sweeps"));
    if(smolevelsublist.get<double>("smoother: damping factor"))
      ifpackList.set("relaxation: damping factor", smolevelsublist.get<double>("smoother: damping factor"));
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Gauss-Seidel");
    smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList) );
  } else if (type == "symmetric Gauss-Seidel") {
    if(smolevelsublist.isParameter("smoother: sweeps"))
      ifpackList.set<int>("relaxation: sweeps", smolevelsublist.get<int>("smoother: sweeps"));
    if(smolevelsublist.get<double>("smoother: damping factor"))
      ifpackList.set("relaxation: damping factor", smolevelsublist.get<double>("smoother: damping factor"));
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
    smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList) );
  } else if (type == "Chebyshev") {
    ifpackType = "CHEBYSHEV";
    if(smolevelsublist.isParameter("smoother: sweeps"))
      ifpackList.set("chebyshev: degree", smolevelsublist.get<int>("smoother: sweeps"));
    smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList) );
    // TODO what about the other parameters
  } else if(type == "IFPACK") {
        // TODO change to TrilnosSmoother as soon as Ifpack2 supports all preconditioners from Ifpack
        ifpackType = params.sublist("smoother: list " + levelstr).get<std::string>("smoother: ifpack type");
        if(ifpackType == "ILU") {
            ifpackList.set<int>("fact: level-of-fill", (int)smolevelsublist.get<double>("smoother: ifpack level-of-fill"));
            ifpackList.set("partitioner: overlap", smolevelsublist.get<int>("smoother: ifpack overlap"));
            smooProto = rcp( new IfpackSmoother(ifpackType, ifpackList, smolevelsublist.get<int>("smoother: ifpack overlap")) );
        }
        else
          TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Interpreter: unknown ML smoother type " + type + " (IFPACK) not supported by MueLu. Only ILU is supported.");
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Interpreter: unknown ML smoother type " + type + " not supported by MueLu.");
  }

  // create smoother factory
  //smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList) );
  SmooFact = rcp( new SmootherFactory(smooProto) );

  // check if pre- and postsmoothing is set
  std::string preorpost = "both";
  if(smolevelsublist.isParameter("smoother: pre or post")) preorpost = smolevelsublist.get<std::string>("smoother: pre or post");

  if (preorpost == "pre") {
    SmooFact->SetSmootherPrototypes(smooProto, Teuchos::null);
  } else if(preorpost == "post") {
    SmooFact->SetSmootherPrototypes(Teuchos::null, smooProto);
  }

  // DEBUG only
  /*RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  RCP<TrilinosSmoother> tSmooProto = Teuchos::rcp_dynamic_cast<TrilinosSmoother>(smooProto);
  std::cout << "SmooProto " << level << std::endl;
  tSmooProto->print(*out, toMueLuVerbLevel(Teuchos::VERB_EXTREME));*/
  //

  return SmooFact;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void MLInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::FillMLParameterList(Teuchos::ParameterList & params) {

  params.set("PDE equations",2);
  params.set("aggregation: damping factor", 1.33);
  params.set("aggregation: nodes per aggregate", 27);
  params.set("aggregation: threshold", 0);
  params.set("aggregation: type", "Uncoupled");
  params.set("coarse: max size", 50);
  params.set("coarse: pre or post", "post");
  params.set("coarse: sweeps", 1);
  params.set("coarse: type", "Amesos-Superludist");
  params.set("max levels", 7);
  params.set("null space: add default vectors", 0);
  params.set("null space: dimension", 3);
  params.set("null space: type", "pre-computed");
  params.set("prec type","MGV");

  Teuchos::ParameterList & l0 = params.sublist("smoother: list (level 0)");
  Teuchos::ParameterList & l1 = params.sublist("smoother: list (level 1)");
  Teuchos::ParameterList & l2 = params.sublist("smoother: list (level 2)");
  Teuchos::ParameterList & l3 = params.sublist("smoother: list (level 3)");
  Teuchos::ParameterList & l4 = params.sublist("smoother: list (level 4)");
  Teuchos::ParameterList & l5 = params.sublist("smoother: list (level 5)");
  Teuchos::ParameterList & l6 = params.sublist("smoother: list (level 6)");

  l0.set("smoother: damping factor", 0.9);
  l0.set("smoother: sweeps", 1);
  l0.set("smoother: type", "symmetric Gauss-Seidel");
  l1.set("smoother: damping factor", 0.9);
  l1.set("smoother: sweeps", 1);
  l1.set("smoother: type", "symmetric Gauss-Seidel");
  l2.set("smoother: damping factor", 0.9);
  l2.set("smoother: sweeps", 1);
  l2.set("smoother: type", "symmetric Gauss-Seidel");
  l3.set("smoother: damping factor", 0.9);
  l3.set("smoother: sweeps", 1);
  l3.set("smoother: type", "symmetric Gauss-Seidel");
  l4.set("smoother: damping factor", 0.89);
  l4.set("smoother: sweeps", 1);
  l4.set("smoother: type", "Jacobi");
  l5.set("smoother: damping factor", 0.89);
  l5.set("smoother: sweeps", 12);
  l5.set("smoother: type", "Jacobi");
  l6.set("smoother: damping factor", 0.89);
  l6.set("smoother: sweeps", 14);
  l6.set("smoother: type", "Jacobi");

}

} // namespace MueLu

#define MUELU_MLINTERPRETER_SHORT
#endif /* MUELU_INTERPRETER_DEF_HPP */
