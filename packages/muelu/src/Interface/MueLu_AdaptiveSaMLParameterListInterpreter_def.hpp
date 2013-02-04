/*
 * MueLu_AdaptiveSaMLParamterListInterpreter_def.hpp
 *
 *  Created on: Jan 28, 2013
 *      Author: tobias
 */

#ifndef MUELU_ADAPTIVESAMLPARAMETERLISTINTERPRETER_DEF_HPP_
#define MUELU_ADAPTIVESAMLPARAMETERLISTINTERPRETER_DEF_HPP_

#include <Teuchos_XMLParameterListHelpers.hpp>

#include "MueLu_ConfigDefs.hpp"
#if defined(HAVE_MUELU_ML) && defined(HAVE_MUELU_EPETRA)
#include <ml_ValidateParameters.h>
#endif

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_AdaptiveSaMLParameterListInterpreter_decl.hpp"

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
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_ParameterListUtils.hpp"
#include "MueLu_MLParameterListInterpreter.hpp"

#include "MueLu_Utilities.hpp"

#include "MueLu_DisableMultipleCallCheck.hpp"

// Note: do not add options that are only recognized by MueLu.

// TODO: this parameter list interpreter should force MueLu to use default ML parameters
// - Ex: smoother sweep=2 by default for ML

// Read a parameter value from a parameter list and store it into a variable named 'varName'
#define MUELU_READ_PARAM(paramList, paramStr, varType, defaultValue, varName) \
  varType varName = defaultValue; if (paramList.isParameter(paramStr)) varName = paramList.get<varType>(paramStr);

// Read a parameter value from a paraeter list and copy it into a new parameter list (with another parameter name)
#define MUELU_COPY_PARAM(paramList, paramStr, varType, defaultValue, outParamList, outParamStr) \
  if (paramList.isParameter(paramStr))                                  \
    outParamList.set<varType>(outParamStr, paramList.get<varType>(paramStr)); \
  else outParamList.set<varType>(outParamStr, defaultValue);            \

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  AdaptiveSaMLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AdaptiveSaMLParameterListInterpreter(Teuchos::ParameterList & paramList,Teuchos::RCP<MultiVector> nspVector, std::vector<RCP<FactoryBase> > factoryList) : TransferFacts_(factoryList), blksize_(1) {
    SetParameterList(paramList);
    //nspVector_ = nspVector;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  AdaptiveSaMLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AdaptiveSaMLParameterListInterpreter(const std::string & xmlFileName, std::vector<RCP<FactoryBase> > factoryList) : nullspace_(NULL), TransferFacts_(factoryList), blksize_(1) {
    Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile(xmlFileName);
    SetParameterList(*paramList);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AdaptiveSaMLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetParameterList(const Teuchos::ParameterList & paramList_in) {
    Teuchos::ParameterList paramList = paramList_in;

    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)); // TODO: use internal out (GetOStream())

    //
    // Read top-level of the parameter list
    //

    // hard-coded default values == ML defaults according to the manual
    MUELU_READ_PARAM(paramList, "ML output",                                int,                   0,       verbosityLevel);
    MUELU_READ_PARAM(paramList, "max levels",                               int,                  10,       maxLevels);
    MUELU_READ_PARAM(paramList, "PDE equations",                            int,                   1,       nDofsPerNode);

    MUELU_READ_PARAM(paramList, "coarse: max size",                         int,                 128,       maxCoarseSize);

    MUELU_READ_PARAM(paramList, "aggregation: type",                std::string,         "Uncoupled",       agg_type);
    MUELU_READ_PARAM(paramList, "aggregation: threshold",                double,                 0.0,       agg_threshold);
    MUELU_READ_PARAM(paramList, "aggregation: damping factor",           double, (double)4/(double)3,       agg_damping);
    MUELU_READ_PARAM(paramList, "aggregation: smoothing sweeps",            int,                   1,       agg_smoothingsweeps);
    MUELU_READ_PARAM(paramList, "aggregation: nodes per aggregate",         int,                   1,       minPerAgg);

    MUELU_READ_PARAM(paramList, "null space: type",                 std::string,   "default vectors",       nullspaceType);
    MUELU_READ_PARAM(paramList, "null space: dimension",                    int,                  -1,       nullspaceDim); // TODO: ML default not in documentation
    MUELU_READ_PARAM(paramList, "null space: vectors",                   double*,               NULL,       nullspaceVec); // TODO: ML default not in documentation

    MUELU_READ_PARAM(paramList, "energy minimization: enable",             bool,               false,       bEnergyMinimization);


    //
    // Move smoothers/aggregation/coarse parameters to sublists
    //

    // ML allows to have level-specific smoothers/aggregation/coarse parameters at the top level of the list or/and defined in sublists:
    // See also: ML Guide section 6.4.1, MueLu::CreateSublists, ML_CreateSublists
    ParameterList paramListWithSubList;
    MueLu::CreateSublists(paramList, paramListWithSubList);
    paramList = paramListWithSubList; // swap

    // std::cout << std::endl << "Parameter list after CreateSublists" << std::endl;
    // std::cout << paramListWithSubList << std::endl;

    //
    // Validate parameter list
    //

    {
      bool validate = paramList.get("ML validate parameter list", true); /* true = default in ML */
      if (validate) {

#if defined(HAVE_MUELU_ML) && defined(HAVE_MUELU_EPETRA)
        // Validate parameter list using ML validator
        int depth = paramList.get("ML validate depth", 5); /* 5 = default in ML */
        TEUCHOS_TEST_FOR_EXCEPTION(! ML_Epetra::ValidateMLPParameters(paramList, depth), Exceptions::RuntimeError,
                                   "ERROR: ML's Teuchos::ParameterList contains incorrect parameter!");
#else
        // If no validator available: issue a warning and set parameter value to false in the output list
        *out << "Warning: MueLu_ENABLE_ML=OFF. The parameter list cannot be validated." << std::endl;
        paramList.set("ML validate parameter list", false);

#endif // HAVE_MUELU_ML
      } // if(validate)
    } // scope


    //
    //
    //

    int maxNbrAlreadySelected = 0;

    // Matrix option
    this->blksize_ = nDofsPerNode;

    // Translate verbosity parameter
    Teuchos::EVerbosityLevel eVerbLevel = Teuchos::VERB_NONE;
    if (verbosityLevel == 0) eVerbLevel = Teuchos::VERB_NONE;
    if (verbosityLevel >  0) eVerbLevel = Teuchos::VERB_LOW;
    if (verbosityLevel >  4) eVerbLevel = Teuchos::VERB_MEDIUM;
    if (verbosityLevel >  7) eVerbLevel = Teuchos::VERB_HIGH;
    if (verbosityLevel >  9) eVerbLevel = Teuchos::VERB_EXTREME;

    TEUCHOS_TEST_FOR_EXCEPTION(agg_type != "Uncoupled" && agg_type != "Coupled", Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter::Setup(): parameter \"aggregation: type\": only 'Uncoupled' or 'Coupled' aggregation is supported.");

    // Create MueLu factories
    // RCP<NullspaceFactory>     nspFact = rcp(new NullspaceFactory());
    RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
    //dropFact->SetVerbLevel(toMueLuVerbLevel(eVerbLevel));

    RCP<FactoryBase> CoupledAggFact = Teuchos::null;
    if(agg_type == "Uncoupled") {
      // Uncoupled aggregation
      RCP<UncoupledAggregationFactory> CoupledAggFact2 = rcp(new UncoupledAggregationFactory());
      CoupledAggFact2->SetMinNodesPerAggregate(minPerAgg); //TODO should increase if run anything other than 1D
      CoupledAggFact2->SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
      CoupledAggFact2->SetOrdering(MueLu::AggOptions::NATURAL);
      CoupledAggFact = CoupledAggFact2;
    } else {
      // Coupled Aggregation (default)
      RCP<CoupledAggregationFactory> CoupledAggFact2 = rcp(new CoupledAggregationFactory());
      CoupledAggFact2->SetMinNodesPerAggregate(minPerAgg); //TODO should increase if run anything other than 1D
      CoupledAggFact2->SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
      CoupledAggFact2->SetOrdering(MueLu::AggOptions::NATURAL);
      CoupledAggFact2->SetPhase3AggCreation(0.5);
      CoupledAggFact = CoupledAggFact2;
    }
    if (verbosityLevel > 3) { // TODO fix me: Setup is a static function: we cannot use GetOStream without an object...
      *out << "========================= Aggregate option summary  =========================" << std::endl;
      *out << "min Nodes per aggregate :               " << minPerAgg << std::endl;
      *out << "min # of root nbrs already aggregated : " << maxNbrAlreadySelected << std::endl;
      *out << "aggregate ordering :                    NATURAL" << std::endl;
      *out << "=============================================================================" << std::endl;
    }

    RCP<Factory> PFact;
    RCP<Factory> RFact;
    RCP<Factory> PtentFact = rcp( new TentativePFactory() );
    if (agg_damping == 0.0 && bEnergyMinimization == false) {
      // tentative prolongation operator (PA-AMG)
      PFact = PtentFact;
      RFact = rcp( new TransPFactory() );
    } else if (agg_damping != 0.0 && bEnergyMinimization == false) {
      // smoothed aggregation (SA-AMG)
      RCP<SaPFactory> SaPFact =  rcp( new SaPFactory() );
      SaPFact->SetDampingFactor(agg_damping);
      PFact  = SaPFact;
      RFact  = rcp( new TransPFactory() );
    } else if (bEnergyMinimization == true) {
      // Petrov Galerkin PG-AMG smoothed aggregation (energy minimization in ML)
      PFact  = rcp( new PgPFactory() );
      RFact  = rcp( new GenericRFactory() );
    }

    RCP<RAPFactory> AcFact = rcp( new RAPFactory() );
    for (size_t i = 0; i<TransferFacts_.size(); i++) {
      AcFact->AddTransferFactory(TransferFacts_[i]); // THIS WILL BE REPLACED with a call to the MLParamterListInterpreter
    }

    //
    // Nullspace factory
    //

    // Set fine level nullspace
    // extract pre-computed nullspace from ML parameter list
    // store it in nullspace_ and nullspaceDim_
    if (nullspaceType != "default vectors") {
      TEUCHOS_TEST_FOR_EXCEPTION(nullspaceType != "pre-computed", Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: no valid nullspace (no pre-computed null space). error.");
      TEUCHOS_TEST_FOR_EXCEPTION(nullspaceDim  == -1,             Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: no valid nullspace (nullspace dim == -1). error.");
      TEUCHOS_TEST_FOR_EXCEPTION(nullspaceVec  == NULL,           Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: no valid nullspace (nullspace == NULL). You have to provide a valid fine-level nullspace in \'null space: vectors\'");

      nullspaceDim_ = nullspaceDim;
      nullspace_    = nullspaceVec;
    }

    Teuchos::RCP<NullspaceFactory> nspFact = Teuchos::rcp(new NullspaceFactory());
    nspFact->SetFactory("Nullspace", PtentFact);

    //
    // Hierarchy + FactoryManager
    //

    // Hierarchy options
    this->SetVerbLevel(toMueLuVerbLevel(eVerbLevel));
    this->numDesiredLevel_ = maxLevels;
    this->maxCoarseSize_   = maxCoarseSize;

    //
    // Coarse Smoother
    //
    ParameterList& coarseList = paramList.sublist("coarse: list ");
    //    coarseList.get("smoother: type", "Amesos-KLU"); // set default
    //RCP<SmootherFactory> coarseFact = this->GetSmootherFactory(coarseList);
    RCP<SmootherFactory> coarseFact = MLParameterListInterpreter::GetSmootherFactory(coarseList);

    // Smoothers Top Level Parameters

    RCP<ParameterList> topLevelSmootherParam = ExtractSetOfParameters(paramList, "smoother");
    // std::cout << std::endl << "Top level smoother parameters:" << std::endl;
    // std::cout << *topLevelSmootherParam << std::endl;

    //

    // Prepare factory managers
    // TODO: smootherFact can be reuse accross level if same parameters/no specific parameterList

    for (int levelID=0; levelID < maxLevels; levelID++) {

      //
      // Level FactoryManager
      //

      RCP<FactoryManager> manager = rcp(new FactoryManager());
      RCP<FactoryManager> initmanager = rcp(new FactoryManager());

      //
      // Smoothers
      //

      {
        // Merge level-specific parameters with global parameters. level-specific parameters takes precedence.
        // TODO: unit-test this part alone

        ParameterList levelSmootherParam = GetMLSubList(paramList, "smoother", levelID); // copy
        MergeParameterList(*topLevelSmootherParam, levelSmootherParam, false); /* false = do no overwrite levelSmootherParam parameters by topLevelSmootherParam parameters */
        // std::cout << std::endl << "Merged List for level  " << levelID << std::endl;
        // std::cout << levelSmootherParam << std::endl;

        //RCP<SmootherFactory> smootherFact = this->GetSmootherFactory(levelSmootherParam); // TODO: missing AFact input arg.
        RCP<SmootherFactory> smootherFact = MLParameterListInterpreter::GetSmootherFactory(levelSmootherParam); // TODO: missing AFact input arg.
        manager->SetFactory("Smoother", smootherFact);
        smootherFact->DisableMultipleCallCheck();

        std::string ifpackType = "RELAXATION";
        Teuchos::ParameterList smootherParamList;
        smootherParamList.set("relaxation: type", "Gauss-Seidel");
        smootherParamList.set("smoother: sweeps", 1);
        smootherParamList.set("smoother: damping factor", 1.0); // not 1.0 since then dirichlet bcs are zeroed out...
        RCP<SmootherPrototype> smooProto = rcp( new TrilinosSmoother(ifpackType, smootherParamList, 0) );

        RCP<SmootherFactory> initSmootherFact = rcp( new SmootherFactory() );
        initSmootherFact->SetSmootherPrototypes(smooProto, smooProto);

        initmanager->SetFactory("Smoother", initSmootherFact);
        initmanager->SetFactory("CoarseSolver", initSmootherFact);
        initSmootherFact->DisableMultipleCallCheck();

      }

      //
      // Misc
      //

      Teuchos::rcp_dynamic_cast<PFactory>(PFact)->DisableMultipleCallCheck();
      Teuchos::rcp_dynamic_cast<PFactory>(PtentFact)->DisableMultipleCallCheck();
      Teuchos::rcp_dynamic_cast<TwoLevelFactoryBase>(RFact)->DisableMultipleCallCheck();
      Teuchos::rcp_dynamic_cast<SingleLevelFactoryBase>(coarseFact)->DisableMultipleCallCheck();
      Teuchos::rcp_dynamic_cast<SingleLevelFactoryBase>(dropFact)->DisableMultipleCallCheck();
      Teuchos::rcp_dynamic_cast<SingleLevelFactoryBase>(CoupledAggFact)->DisableMultipleCallCheck();
      Teuchos::rcp_dynamic_cast<TwoLevelFactoryBase>(AcFact)->DisableMultipleCallCheck();
      Teuchos::rcp_dynamic_cast<SingleLevelFactoryBase>(nspFact)->DisableMultipleCallCheck();

      manager->SetFactory("CoarseSolver", coarseFact); // TODO: should not be done in the loop
      manager->SetFactory("Graph", dropFact);
      manager->SetFactory("Aggregates", CoupledAggFact);
      manager->SetFactory("DofsPerNode", dropFact);
      manager->SetFactory("A", AcFact);
      manager->SetFactory("P", PFact);
      manager->SetFactory("Ptent", PtentFact);
      manager->SetFactory("R", RFact);
      manager->SetFactory("Nullspace", nspFact);

      //initmanager->SetFactory("CoarseSolver", coarseFact);
      initmanager->SetFactory("Graph", dropFact);
      initmanager->SetFactory("Aggregates", CoupledAggFact);
      initmanager->SetFactory("DofsPerNode", dropFact);
      initmanager->SetFactory("A", AcFact);
      initmanager->SetFactory("P", PtentFact); // use nonsmoothed transfers
      initmanager->SetFactory("Ptent", PtentFact);
      initmanager->SetFactory("R", RFact);
      initmanager->SetFactory("Nullspace", nspFact);

      this->AddFactoryManager(levelID, 1, manager);
      this->AddInitFactoryManager(levelID, 1, initmanager);
    } // for (level loop)
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AdaptiveSaMLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetupInitHierarchy(Hierarchy & H) const {
    TEUCHOS_TEST_FOR_EXCEPTION(!H.GetLevel(0)->IsAvailable("A"), Exceptions::RuntimeError, "No fine level operator");

    RCP<Level> l = H.GetLevel(0);
    RCP<Matrix> Op = l->Get<RCP<Matrix> >("A");
    SetupMatrix(*Op); // use overloaded SetupMatrix routine
    this->SetupExtra(H);

    // Setup Hierarchy
    H.SetMaxCoarseSize(this->maxCoarseSize_); // TODO

    int  levelID     = 0;
    int  lastLevelID = this->numDesiredLevel_ - 1;
    bool isLastLevel = false;

    while(!isLastLevel) {
      bool r = H.Setup(levelID,
          InitLvlMngr(levelID-1, lastLevelID),
          InitLvlMngr(levelID,   lastLevelID),
          InitLvlMngr(levelID+1, lastLevelID));

      isLastLevel = r || (levelID == lastLevelID);
      levelID++;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AdaptiveSaMLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetupHierarchy(Hierarchy & H) const {

    // set fine level null space
    // usually this null space is provided from outside (by the user) using
    // the ML parameter lists.
    if (this->nullspace_ != NULL) {
      RCP<Level> fineLevel = H.GetLevel(0);
      const RCP<const Map> rowMap = fineLevel->Get< RCP<Matrix> >("A")->getRowMap();
      RCP<MultiVector> nullspace = MultiVectorFactory::Build(rowMap, nullspaceDim_, true);

      for ( size_t i=0; i < Teuchos::as<size_t>(nullspaceDim_); i++) {
        Teuchos::ArrayRCP<Scalar> nullspacei = nullspace->getDataNonConst(i);
        const size_t              myLength   = nullspace->getLocalLength();

        for (size_t j = 0; j < myLength; j++) {
          nullspacei[j] = nullspace_[i*myLength + j];
        }
      }

      fineLevel->Set("Nullspace", nullspace);
    }

    // keep aggregates
    H.Keep("Aggregates", HierarchyManager::GetFactoryManager(0)->GetFactory("Aggregates").get());

    ///////////////////////////////

    // build hierarchy for initialization
    SetupInitHierarchy(H);

    {
      // do some iterations with the built hierarchy to improve the null space
      Teuchos::RCP<MueLu::Level> Finest = H.GetLevel(0);  // get finest level,MueLu::NoFactory::get()
      Teuchos::RCP<MultiVector> nspVector2 = Finest->Get<Teuchos::RCP<MultiVector> >("Nullspace");

      MueLu::Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Write("orig_nsp.vec", *nspVector2);

      Teuchos::RCP<MultiVector> homogRhsVec = MultiVectorFactory::Build(nspVector2->getMap(),nspVector2->getNumVectors(),true);
      //homogRhsVec->putScalar(0.0);

      // do 1 multigrid cycle for improving the null space by "solving"
      //     A B_f = 0
      // where A is the system matrix and B_f the fine level null space vectors
      H.Iterate(*homogRhsVec, 1, *nspVector2, false);

      // store improved fine level null space
      Finest->Set("Nullspace",nspVector2);

      MueLu::Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Write("new_nsp.vec", *nspVector2);

    }

    {
      // do some clean up.
      // remove all old default factories. Build new ones for the second build.
      // this is a little bit tricky to understand
      for(size_t k=0; k < HierarchyManager::getNumFactoryManagers(); k++) {
        HierarchyManager::GetFactoryManager(k)->Clean();
        //Teuchos::rcp_dynamic_cast<const SingleLevelFactoryBase>(HierarchyManager::GetFactoryManager(k)->GetFactory("Smoother"))->DisableMultipleCallCheck(); // after changing to MLParamterListInterpreter functions
      }
      // not sure about this. i only need it if Smoother is defined explicitely (not using default smoother)
      // need this: otherwise RAPFactory::Build is complaining on level 0
      //            and TentativePFactory::Build is complaining on level 1
      Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(HierarchyManager::GetFactoryManager(0)->GetFactory("A"))->DisableMultipleCallCheck();
      Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(HierarchyManager::GetFactoryManager(1)->GetFactory("P"))->DisableMultipleCallCheck();

      HierarchyManager::SetupHierarchy(H);
    }

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AdaptiveSaMLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddTransferFactory(const RCP<FactoryBase>& factory) {
    // check if it's a TwoLevelFactoryBase based transfer factory
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast, "Transfer factory is not derived from TwoLevelFactoryBase. Since transfer factories will be handled by the RAPFactory they have to be derived from TwoLevelFactoryBase!");
    TransferFacts_.push_back(factory);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t AdaptiveSaMLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::NumTransferFactories() const {
    return TransferFacts_.size();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AdaptiveSaMLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetupMatrix(Matrix & Op) const {
    Op.SetFixedBlockSize(blksize_);
  }

} // namespace MueLu


#endif /* MUELU_ADAPTIVESAMLPARAMETERLISTINTERPRETER_DEF_HPP_ */
