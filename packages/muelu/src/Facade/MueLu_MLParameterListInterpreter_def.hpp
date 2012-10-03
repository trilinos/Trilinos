// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, 
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_MLPARAMETERLISTINTERPRETER_DEF_HPP
#define MUELU_MLPARAMETERLISTINTERPRETER_DEF_HPP

#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_MLParameterListInterpreter_decl.hpp"

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
  MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MLParameterListInterpreter(Teuchos::ParameterList & paramList, std::vector<RCP<FactoryBase> > factoryList) : nullspace_(NULL), TransferFacts_(factoryList), blksize_(1) {
    SetParameterList(paramList);
  }
  
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MLParameterListInterpreter(const std::string & xmlFileName, std::vector<RCP<FactoryBase> > factoryList) : nullspace_(NULL), TransferFacts_(factoryList), blksize_(1) {
    Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::getParametersFromXmlFile(xmlFileName);
    SetParameterList(*paramList);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetParameterList(const Teuchos::ParameterList & paramList) {

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
    //
    //

    int maxNbrAlreadySelected = 0;

    // Matrix option
    blksize_ = nDofsPerNode;

    // Translate verbosity parameter
    Teuchos::EVerbosityLevel eVerbLevel = Teuchos::VERB_NONE;
    if (verbosityLevel == 0) eVerbLevel = Teuchos::VERB_NONE;
    if (verbosityLevel >  0) eVerbLevel = Teuchos::VERB_LOW;
    if (verbosityLevel >  4) eVerbLevel = Teuchos::VERB_MEDIUM;
    if (verbosityLevel >  7) eVerbLevel = Teuchos::VERB_HIGH;
    if (verbosityLevel >  9) eVerbLevel = Teuchos::VERB_EXTREME;

    TEUCHOS_TEST_FOR_EXCEPTION(agg_type != "Uncoupled", Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter::Setup(): parameter \"aggregation: type\": only 'Uncoupled' aggregation is supported.");

    // Create MueLu factories
    // RCP<NullspaceFactory>     nspFact = rcp(new NullspaceFactory());
    RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
    //dropFact->SetVerbLevel(toMueLuVerbLevel(eVerbLevel));

    RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
    if (verbosityLevel > 3) { // TODO fix me: Setup is a static function: we cannot use GetOStream without an object...
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

    RCP<PFactory> PtentFact = rcp(new TentativePFactory(UCAggFact));

    if (agg_damping == 0.0 && bEnergyMinimization == false) {
      // tentative prolongation operator (PA-AMG)
      PFact = PtentFact; //rcp( new TentativePFactory() );
      RFact = rcp( new TransPFactory() );
    } else if (agg_damping != 0.0 && bEnergyMinimization == false) {
      // smoothed aggregation (SA-AMG)
      RCP<SaPFactory> SaPFact =  rcp( new SaPFactory(PtentFact) );
      SaPFact->SetDampingFactor(agg_damping);
      PFact  = SaPFact;
      RFact  = rcp( new TransPFactory() );
    } else if (bEnergyMinimization == true) {
      // Petrov Galerkin PG-AMG smoothed aggregation (energy minimization in ML)
      PFact  = rcp( new PgPFactory(PtentFact) );
      RFact  = rcp( new GenericRFactory() );
    }

    RCP<RAPFactory> AcFact = rcp( new RAPFactory() );
    for (size_t i = 0; i<TransferFacts_.size(); i++) {
      AcFact->AddTransferFactory(TransferFacts_[i]);
    }
    
    std::string t("coarse");
    RCP<SmootherFactory> coarsestSmooFact = GetSmootherFactory(paramList, std::string("coarse")); // TODO: missing AFact input arg.

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

    Teuchos::RCP<NullspaceFactory> nspFact = Teuchos::rcp(new NullspaceFactory("Nullspace", PtentFact));

    //
    // Hierarchy + FactoryManager
    //

    // Hierarchy options
    this->verbLevel_       = toMueLuVerbLevel(eVerbLevel);
    this->numDesiredLevel_ = maxLevels;
    this->maxCoarseSize_   = maxCoarseSize;

    // Prepare factory managers
    for (int levelID=0; levelID < maxLevels; levelID++) {
      RCP<SmootherFactory> SmooFactFine   = GetSmootherFactory(paramList, std::string("smoother"), levelID); // TODO: missing AFact input arg.
      RCP<FactoryManager> manager = rcp(new FactoryManager());
      if (SmooFactFine != Teuchos::null)
    	manager->SetFactory("Smoother" ,  SmooFactFine);    // Hierarchy.Setup uses TOPSmootherFactory, that only needs "Smoother"
      manager->SetFactory("CoarseSolver", coarsestSmooFact);
      manager->SetFactory("Graph", dropFact);
      manager->SetFactory("Aggregates", UCAggFact);
      manager->SetFactory("DofsPerNode", dropFact);
      manager->SetFactory("A", AcFact);                     // same RAP factory
      manager->SetFactory("P", PFact);                      // same prolongator and restrictor factories
      manager->SetFactory("Ptent", PtentFact);              // same prolongator and restrictor factories
      manager->SetFactory("R", RFact);                      // same prolongator and restrictor factories
      manager->SetFactory("Nullspace", nspFact);            // use same nullspace factory throughout all multigrid levels

      this->AddFactoryManager(levelID, 1, manager);
    }

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetupHierarchy(Hierarchy & H) const {

    // if nullspace_ has already been extracted from ML parameter list
    if (nullspace_ != NULL) {
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
    HierarchyManager::SetupHierarchy(H);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetSmootherFactory(const Teuchos::ParameterList & topLevelParamList, const std::string & smootherTypeStr, const int levelID, const RCP<FactoryBase> & AFact) {

    //
    // Select the parameter list to be parsed 
    //

    Teuchos::ParameterList paramList; // parameter list that will be parsed
    std::string levelStr;             // string uses to identify the sub parameter list that will be parsed

    if (smootherTypeStr == "smoother") {
      TEUCHOS_TEST_FOR_EXCEPT(levelID == -1);

      char levelChar[11];
      sprintf(levelChar, "(level %d)", levelID);
      levelStr = levelChar;
    
      if (paramList.isSublist(smootherTypeStr + ": list " + levelStr) == false)
        return Teuchos::null;
      //TEUCHOS_TEST_FOR_EXCEPTION(paramList.isSublist("smoother: list " + levelStr)==false, Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: no ML smoother parameter list for level. error.");
      
      paramList = paramList.sublist("smoother: list " + levelStr);
  
    } else if (smootherTypeStr == "coarse") {

      levelStr = "";
      paramList = topLevelParamList;

    } else {

      TEUCHOS_TEST_FOR_EXCEPT(true);

    }

    //
    // Get 'type'
    //

    //TODO: fix defaults!!

    // Default coarse grid smoother
    std::string type;
    if (smootherTypeStr == "coarse") {
#if (defined(HAVE_MUELU_EPETRA) && defined( HAVE_MUELU_AMESOS)) || (defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2)) // FIXME: test is wrong (ex: compiled with Epetra&&Tpetra&&Amesos2 but without Amesos => error running Epetra problem)
      type = ""; // use default defined by AmesosSmoother or Amesos2Smoother
#else
      type = "symmetric Gauss-Seidel"; // use a sym Gauss-Seidel (with no damping) as fallback "coarse solver" (TODO: needs Ifpack(2))
#endif
    } else {
      // TODO: default smoother?
      type = "";
    }

    if (paramList.isParameter(smootherTypeStr + ": type")) type = paramList.get<std::string>(smootherTypeStr + ": type");
    TEUCHOS_TEST_FOR_EXCEPTION(type.empty(), Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: no ML smoother type for level. error.");

    //
    // Create the smoother prototype
    //

    RCP<SmootherPrototype> smooProto;
    std::string ifpackType;
    Teuchos::ParameterList smootherParamList;

    if (type == "Jacobi" || type == "Gauss-Seidel" || type == "symmetric Gauss-Seidel") {
      if (type == "symmetric Gauss-Seidel") type = "Symmetric Gauss-Seidel"; // FIXME

      ifpackType = "RELAXATION";
      smootherParamList.set("relaxation: type", type);

      MUELU_COPY_PARAM(paramList, smootherTypeStr + ": sweeps",            int,   2, smootherParamList, "relaxation: sweeps");
      MUELU_COPY_PARAM(paramList, smootherTypeStr + ": damping factor", double, 1.0, smootherParamList, "relaxation: damping factor");

      smooProto = rcp( new TrilinosSmoother(ifpackType, smootherParamList, 0, AFact) );

    } else if (type == "Chebyshev") {

      ifpackType = "CHEBYSHEV";

      MUELU_COPY_PARAM(paramList, smootherTypeStr + ": sweeps",          int, 2,  smootherParamList, "chebyshev: degree");
      MUELU_COPY_PARAM(paramList, smootherTypeStr + ": Chebyshev alpha", int, 30, smootherParamList, "chebyshev: alpha");

      smooProto = rcp( new TrilinosSmoother(ifpackType, smootherParamList, 0, AFact) );

    } else if (type == "IFPACK") { // TODO: this option is not described in the ML Guide v5.0

#ifdef HAVE_MUELU_IFPACK
      ifpackType = paramList.get<std::string>(smootherTypeStr + ": ifpack type");

      if (ifpackType == "ILU") {
        MUELU_COPY_PARAM(paramList, smootherTypeStr + ": ifpack level-of-fill", int, 2,  smootherParamList, "fact: level-of-fill");
        MUELU_COPY_PARAM(paramList, smootherTypeStr + ": ifpack overlap",       int, 2,  smootherParamList, "partitioner: overlap");

        // TODO change to TrilinosSmoother as soon as Ifpack2 supports all preconditioners from Ifpack
        smooProto = MueLu::GetIfpackSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(ifpackType, smootherParamList, paramList.get<int>(smootherTypeStr + ": ifpack overlap"), AFact);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: unknown ML smoother type " + type + " (IFPACK) not supported by MueLu. Only ILU is supported.");
      }
#else // HAVE_MUELU_IFPACK
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: MueLu compiled without Ifpack support");
#endif // HAVE_MUELU_IFPACK

    } else if (type.length() > strlen("Amesos") && type.substr(0, strlen("Amesos")) == "Amesos") {  /* catch Amesos-* */
      std::string solverType = type.substr(strlen("Amesos")+1);  /* ("Amesos-KLU" -> "KLU") */

      // Validator: following upper/lower case is what is allowed by ML
      bool valid = false;
      const int  validatorSize = 5;
      std::string validator[validatorSize] = {"Superlu", "Superludist", "KLU", "UMFPACK"}; /* TODO: should "" be allowed? */
      for (int i=0; i < validatorSize; i++) { if (validator[i] == solverType) valid = true; }
      TEUCHOS_TEST_FOR_EXCEPTION(!valid, Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: unknown smoother type. '" << type << "' not supported.");

      // FIXME: MueLu should accept any Upper/Lower case. Not the case for the moment
      std::transform(solverType.begin()+1, solverType.end(), solverType.begin()+1, ::tolower);

      smooProto = Teuchos::rcp( new DirectSolver(solverType, Teuchos::ParameterList(), AFact) );

    } else {

      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: unknown smoother type. '" << type << "' not supported by MueLu.");

    }
    TEUCHOS_TEST_FOR_EXCEPTION(smooProto == Teuchos::null, Exceptions::RuntimeError, "MueLu::MLParameterListInterpreter: no smoother prototype. fatal error.");

    //
    // Create the smoother factory
    //

    RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory() );

    // Set parameters of the smoother factory
    MUELU_READ_PARAM(paramList, smootherTypeStr + ": pre or post", std::string, "both", preOrPost);
    if (preOrPost == "both") {
      SmooFact->SetSmootherPrototypes(smooProto, smooProto);
    } else if (preOrPost == "pre") {
      SmooFact->SetSmootherPrototypes(smooProto, Teuchos::null);
    } else if (preOrPost == "post") {
      SmooFact->SetSmootherPrototypes(Teuchos::null, smooProto);
    }

    return SmooFact;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddTransferFactory(const RCP<FactoryBase>& factory) {
    // check if it's a TwoLevelFactoryBase based transfer factory
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast, "Transfer factory is not derived from TwoLevelFactoryBase. Since transfer factories will be handled by the RAPFactory they have to be derived from TwoLevelFactoryBase!");
    TransferFacts_.push_back(factory);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::NumTransferFactories() const {
    return TransferFacts_.size();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetupMatrix(Matrix & Op) const {
      Op.SetFixedBlockSize(blksize_);
  }

} // namespace MueLu

#define MUELU_MLPARAMETERLISTINTERPRETER_SHORT
#endif /* MUELU_MLPARAMETERLISTINTERPRETER_DEF_HPP */

//TODO: see if it can be factorized with ML interpreter (ex: generation of Ifpack param list)
