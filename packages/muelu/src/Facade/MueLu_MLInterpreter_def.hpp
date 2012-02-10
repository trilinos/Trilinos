#ifndef MUELU_MLINTERPRETER_DEF_HPP
#define MUELU_MLINTERPRETER_DEF_HPP

#include <Teuchos_XMLParameterListHelpers.hpp>

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
  MLInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MLInterpreter(Teuchos::ParameterList & paramList) : nullspace_(NULL) {
    SetParameterList(paramList);
  }
  
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  MLInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MLInterpreter(const std::string & xmlFileName) : nullspace_(NULL) {
    Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile(xmlFileName);
    SetParameterList(*paramList);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void MLInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetParameterList(const Teuchos::ParameterList & paramList) {

    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)); // TODO: use internal out (GetOStream())

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

    if(paramList.isParameter("max levels"))                      maxLevels           = paramList.get<int>("max levels");
    if(paramList.isParameter("ML output"))                       verbosityLevel      = paramList.get<int>("ML output");
    if(paramList.isParameter("coarse: max size"))                maxCoarseSize       = paramList.get<int>("coarse: max size");
    if(paramList.isParameter("PDE equations"))                   nDofsPerNode        = paramList.get<int>("PDE equations");
    if(paramList.isParameter("aggregation: threshold"))          agg_threshold       = paramList.get<double>("aggregation: threshold");
    if(paramList.isParameter("aggregation: damping factor"))     agg_damping         = paramList.get<double>("aggregation: damping factor");
    if(paramList.isParameter("aggregation: smoothing sweeps"))   agg_smoothingsweeps = paramList.get<int>   ("aggregation: smoothing sweeps");
    if(paramList.isParameter("aggregation: type"))               agg_type            = paramList.get<std::string> ("aggregation: type");
    //if(paramList.isParameter("aggregation: nodes per aggregate")) minPerAgg        = paramList.get<int>("aggregation: nodes per aggregate");
    if(paramList.isParameter("energy minimization: enable"))     bEnergyMinimization = paramList.get<bool>("energy minimization: enable");

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
      PFact = rcp( new TentativePFactory(UCAggFact/*,nspFact*/) );
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
    coarsestSmooFact = GetCoarsestSolverFactory(paramList);

    //
    //
    //

    // Set fine level nullspace
    // extract pre-computed nullspace from ML parameter list
    if(paramList.isParameter("null space: type")) {
      std::string type = "";
      type = paramList.get<std::string>("null space: type");
      TEUCHOS_TEST_FOR_EXCEPTION(type != "pre-computed", Exceptions::RuntimeError, "MueLu::Interpreter: no valid nullspace (no pre-computed null space). error.");

      int dimns = -1;
      if(paramList.isParameter("null space: dimension")) dimns = paramList.get<int>("null space: dimension");
      TEUCHOS_TEST_FOR_EXCEPTION(dimns == -1, Exceptions::RuntimeError, "MueLu::Interpreter: no valid nullspace (nullspace dim = -1). error.");

      double* nsdata = NULL;
      if(paramList.isParameter("null space: vectors")) nsdata = paramList.get<double*>("null space: vectors");
      TEUCHOS_TEST_FOR_EXCEPTION(nsdata == NULL, Exceptions::RuntimeError, "MueLu::Interpreter: no valid nullspace (nsdata = NULL). error.");

      nullspaceDim_ = dimns;
      nullspace_    = nsdata;
    }

    //
    //
    //

    // Hierarchy options
    this->verbLevel_       = toMueLuVerbLevel(eVerbLevel);
    this->numDesiredLevel_ = maxLevels;
    this->maxCoarseSize_   = maxCoarseSize;

    // Prepare factory managers
    for(int levelID=0; levelID < maxLevels; levelID++) {
      RCP<SmootherFactory> SmooFactFine   = GetSmootherFactory(paramList, levelID);
      RCP<FactoryManager> manager = rcp(new FactoryManager());
      if(SmooFactFine != Teuchos::null)
    	manager->SetFactory("Smoother" ,  SmooFactFine);    // Hierarchy.Setup uses TOPSmootherFactory, that only needs "Smoother"
      manager->SetFactory("CoarseSolver", coarsestSmooFact);
      manager->SetFactory("A", AcFact);                     // same RAP factory
      manager->SetFactory("P", PFact);                      // same prolongator and restrictor factories
      manager->SetFactory("R", RFact);                      // same prolongator and restrictor factories
      manager->SetFactory("Nullspace", nspFact);            // use same nullspace factory throughout all multigrid levels

      this->AddFactoryManager(levelID, 1, manager);
    }

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void MLInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetupHierarchy(Hierarchy & H) const {

    if (nullspace_ != NULL) {
      
      RCP<Level> fineLevel = H.GetLevel(0);
      const RCP<const Map> rowMap = fineLevel->Get< RCP<Operator> >("A")->getRowMap();
      RCP<MultiVector> nullspace = MultiVectorFactory::Build(rowMap, nullspaceDim_, true);
      
      for ( size_t i=0; i < Teuchos::as<size_t>(nullspaceDim_); i++) {
        Teuchos::ArrayRCP<Scalar> nullspacei = nullspace->getDataNonConst(i);
        const size_t              myLength   = nullspace->getLocalLength();

        for(size_t j = 0; j < myLength; j++) {
          nullspacei[j] = nullspace_[i*myLength + j];
        }
      }

      fineLevel->Set("Nullspace", nullspace_);
    }
    
    HierarchyManager::SetupHierarchy(H);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > MLInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetCoarsestSolverFactory(const Teuchos::ParameterList & paramList) {

#if (defined(HAVE_MUELU_EPETRA) && defined( HAVE_MUELU_AMESOS)) || (defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2)) // FIXME: test is wrong (ex: compiled with Epetra&&Tpetra&&Amesos2 but without Amesos => error running Epetra problem)
    std::string type = ""; // use default defined by AmesosSmoother or Amesos2Smoother
#else
    std::string type = "symmetric Gauss-Seidel"; // use a sym Gauss-Seidel (with no damping) as fallback "coarse solver" (TODO: needs Ifpack(2))
#endif

    if(paramList.isParameter("coarse: type")) type = paramList.get<std::string>("coarse: type");

    RCP<SmootherPrototype> smooProto;
    std::string ifpackType;
    Teuchos::ParameterList ifpackList;
    RCP<SmootherFactory> SmooFact;

    if(type == "Jacobi") {
      if(paramList.isParameter("coarse: sweeps"))
        ifpackList.set<int>("relaxation: sweeps", paramList.get<int>("coarse: sweeps"));
      else ifpackList.set<int>("relaxation: sweeps", 1);
      if(paramList.isParameter("coarse: damping factor"))
        ifpackList.set("relaxation: damping factor", paramList.get<double>("coarse: damping factor"));
      else ifpackList.set("relaxation: damping factor", 1.0);
      ifpackType = "RELAXATION";
      ifpackList.set("relaxation: type", "Jacobi");
      smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList) );
    } else if(type == "Gauss-Seidel") {
      if(paramList.isParameter("coarse: sweeps"))
        ifpackList.set<int>("relaxation: sweeps", paramList.get<int>("coarse: sweeps"));
      else ifpackList.set<int>("relaxation: sweeps", 1);
      if(paramList.isParameter("coarse: damping factor"))
        ifpackList.set("relaxation: damping factor", paramList.get<double>("coarse: damping factor"));
      else ifpackList.set("relaxation: damping factor", 1.0);
      ifpackType = "RELAXATION";
      ifpackList.set("relaxation: type", "Gauss-Seidel");
      smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList) );
    } else if (type == "symmetric Gauss-Seidel") {
      if(paramList.isParameter("coarse: sweeps"))
        ifpackList.set<int>("relaxation: sweeps", paramList.get<int>("coarse: sweeps"));
      else ifpackList.set<int>("relaxation: sweeps", 1);
      if(paramList.isParameter("coarse: damping factor"))
        ifpackList.set("relaxation: damping factor", paramList.get<double>("coarse: damping factor"));
      else ifpackList.set("relaxation: damping factor", 1.0);
      ifpackType = "RELAXATION";
      ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
      smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList) );
    } else if (type == "Chebyshev") {
      ifpackType = "CHEBYSHEV";
      if(paramList.isParameter("coarse: sweeps"))
        ifpackList.set("chebyshev: degree", paramList.get<int>("coarse: sweeps"));
      if(paramList.isParameter("coarse: Chebyshev alpha"))
        ifpackList.set("chebyshev: alpha", paramList.get<double>("coarse: Chebyshev alpha"));
      smooProto = rcp( new TrilinosSmoother(ifpackType, ifpackList) );
    } else if(type == "IFPACK") {
#ifdef HAVE_MUELU_IFPACK
      // TODO change to TrilinosSmoother as soon as Ifpack2 supports all preconditioners from Ifpack
      ifpackType = paramList.get<std::string>("coarse: ifpack type");
      if(ifpackType == "ILU") {
        ifpackList.set<int>("fact: level-of-fill", (int)paramList.get<double>("coarse: ifpack level-of-fill"));
        ifpackList.set("partitioner: overlap", paramList.get<int>("coarse: ifpack overlap"));
        smooProto = rcp( new IfpackSmoother(ifpackType, ifpackList, paramList.get<int>("coarse: ifpack overlap")) );
      }
      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Interpreter: unknown ML smoother type " + type + " (IFPACK) not supported by MueLu. Only ILU is supported.");
#else // HAVE_MUELU_IFPACK
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Interpreter: MueLu compiled without Ifpack support");
#endif // HAVE_MUELU_IFPACK
    } else if(type == "Amesos-Superlu") {
      smooProto = Teuchos::rcp( new DirectSolver("Superlu") );
    } else if(type == "Amesos-Superludist") {
      smooProto = Teuchos::rcp( new DirectSolver("Superludist") );
    } else if(type == "Amesos-KLU") {
      smooProto = Teuchos::rcp( new DirectSolver("Klu") );
    } else if(type == "") {
      smooProto = Teuchos::rcp( new DirectSolver("") );
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Interpreter: unknown coarsest solver type. '" << type << "' not supported by MueLu.");
    }

    // create smoother factory
    TEUCHOS_TEST_FOR_EXCEPTION(smooProto == Teuchos::null, Exceptions::RuntimeError, "MueLu::Interpreter: no smoother prototype. fatal error.");
    SmooFact = rcp( new SmootherFactory(smooProto) );

    // check if pre- and postsmoothing is set
    std::string preorpost = "both";
    if(paramList.isParameter("coarse: pre or post")) preorpost = paramList.get<std::string>("coarse: pre or post");

    if (preorpost == "pre") {
      SmooFact->SetSmootherPrototypes(smooProto, Teuchos::null);
    } else if(preorpost == "post") {
      SmooFact->SetSmootherPrototypes(Teuchos::null, smooProto);
    }

    return SmooFact;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > MLInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetSmootherFactory(const Teuchos::ParameterList & paramList, int level) {

    char levelchar[11];
    sprintf(levelchar,"(level %d)",level);
    std::string levelstr(levelchar);

    if(paramList.isSublist("smoother: list " + levelstr)==false)
      return Teuchos::null;
    //TEUCHOS_TEST_FOR_EXCEPTION(paramList.isSublist("smoother: list " + levelstr)==false, Exceptions::RuntimeError, "MueLu::Interpreter: no ML smoother parameter list for level. error.");

    std::string type = paramList.sublist("smoother: list " + levelstr).get<std::string>("smoother: type");
    TEUCHOS_TEST_FOR_EXCEPTION(type.empty(), Exceptions::RuntimeError, "MueLu::Interpreter: no ML smoother type for level. error.");

    const Teuchos::ParameterList smolevelsublist = paramList.sublist("smoother: list " + levelstr);
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
#ifdef HAVE_MUELU_IFPACK
      // TODO change to TrilinosSmoother as soon as Ifpack2 supports all preconditioners from Ifpack
      ifpackType = paramList.sublist("smoother: list " + levelstr).get<std::string>("smoother: ifpack type");
      if(ifpackType == "ILU") {
        ifpackList.set<int>("fact: level-of-fill", (int)smolevelsublist.get<double>("smoother: ifpack level-of-fill"));
        ifpackList.set("partitioner: overlap", smolevelsublist.get<int>("smoother: ifpack overlap"));
        smooProto = rcp( new IfpackSmoother(ifpackType, ifpackList, smolevelsublist.get<int>("smoother: ifpack overlap")) );
      }
      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Interpreter: unknown ML smoother type " + type + " (IFPACK) not supported by MueLu. Only ILU is supported.");
#else // HAVE_MUELU_IFPACK
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Interpreter: MueLu compiled without Ifpack support");
#endif // HAVE_MUELU_IFPACK
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
  void MLInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::FillMLParameterList(Teuchos::ParameterList & paramList) {

    paramList.set("PDE equations",2);
    paramList.set("aggregation: damping factor", 1.33);
    paramList.set("aggregation: nodes per aggregate", 27);
    paramList.set("aggregation: threshold", 0);
    paramList.set("aggregation: type", "Uncoupled");
    paramList.set("coarse: max size", 50);
    paramList.set("coarse: pre or post", "post");
    paramList.set("coarse: sweeps", 1);
    paramList.set("coarse: type", "Amesos-Superludist");
    paramList.set("max levels", 7);
//     paramList.set("null space: add default vectors", 0);
//     paramList.set("null space: dimension", 3);
//     paramList.set("null space: type", "pre-computed");
    paramList.set("prec type","MGV");

    Teuchos::ParameterList & l0 = paramList.sublist("smoother: list (level 0)");
    Teuchos::ParameterList & l1 = paramList.sublist("smoother: list (level 1)");
    Teuchos::ParameterList & l2 = paramList.sublist("smoother: list (level 2)");
    Teuchos::ParameterList & l3 = paramList.sublist("smoother: list (level 3)");
    Teuchos::ParameterList & l4 = paramList.sublist("smoother: list (level 4)");
    Teuchos::ParameterList & l5 = paramList.sublist("smoother: list (level 5)");
    Teuchos::ParameterList & l6 = paramList.sublist("smoother: list (level 6)");

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
