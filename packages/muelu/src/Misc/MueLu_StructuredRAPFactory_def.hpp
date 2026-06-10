// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_STRUCTUREDRAPFACTORY_DEF_HPP
#define MUELU_STRUCTUREDRAPFACTORY_DEF_HPP

#include <sstream>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <stdexcept>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_TripleMatrixMultiply.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_CrsGraph.hpp>

#include "MueLu_StructuredRAPFactory_decl.hpp"

#include "MueLu_Utilities.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_Behavior.hpp"
#include "Teuchos_TestForException.hpp"
#include "MueLu_Behavior.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::StructuredRAPFactory()
  : hasDeclaredInput_(false) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~StructuredRAPFactory() = default;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("transpose: use implicit");     // should always be used right now, guarantees symmetry
  SET_VALID_ENTRY("rap: triple product");         // in the long term this has to be the only option for multiplication
  SET_VALID_ENTRY("rap: structure type");         // set by user to define matrix structure (e.g. Laplace2D)
  SET_VALID_ENTRY("rap: fix zero diagonals");
  SET_VALID_ENTRY("rap: fix zero diagonals threshold");
  SET_VALID_ENTRY("rap: fix zero diagonals replacement");
  SET_VALID_ENTRY("rap: relative diagonal floor");
#undef SET_VALID_ENTRY
  validParamList->set<RCP<const FactoryBase> >("A", null, "Generating factory of the matrix A used during the prolongator smoothing process");
  validParamList->set<RCP<const FactoryBase> >("P", null, "Prolongator factory");
  validParamList->set<RCP<const FactoryBase> >("R", null, "Restrictor factory");
  validParamList->set< RCP<const FactoryBase> >("numDimensions",      null, "Number of spacial dimensions in the problem.");
  validParamList->set< RCP<const FactoryBase> >("lCoarseNodesPerDim", null, "Number of nodes per spatial dimension on the coarse grid.");

  validParamList->set<bool>("CheckMainDiagonal", false, "Check main diagonal for zeros");
  validParamList->set<bool>("RepairMainDiagonal", false, "Repair zeros on main diagonal");

  // Make sure we don't recursively validate options for the matrixmatrix kernels
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  const Teuchos::ParameterList& pL = GetParameterList();
  if (!pL.get<bool>("transpose: use implicit"))
    Input(coarseLevel, "R");

  Input(fineLevel, "A");
  Input(coarseLevel, "P");

  // get structure information
  Input(fineLevel, "numDimensions");
  Input(fineLevel, "lCoarseNodesPerDim");

  // call DeclareInput of all user-given transfer factories
  for (std::vector<RCP<const FactoryBase> >::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it)
    (*it)->CallDeclareInput(coarseLevel);

  hasDeclaredInput_ = true;
}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetLaplace1D(RCP<Matrix>& Ac, RCP<Matrix> P,
                                                                              	     Teuchos::Array<LocalOrdinal> lCoarseNodesPerDim) const
  {         
    // Define some containers for the compressed row storage
    int ncoarse = lCoarseNodesPerDim[0];
    int maxNnzOnRow = 3; // nnzOnRow is here ~3 for Laplace1D
    const RCP<ParameterList> paramList = Teuchos::null;
  
    // get graph container  
    auto rowMap = P->getDomainMap();
    auto colMap = P->getColMap();
    RCP<CrsGraph> myGraph = CrsGraphFactory::Build(rowMap, colMap, maxNnzOnRow, paramList);
    
    int end = 4+(ncoarse-2)*3;
    const ArrayRCP<LO> colind(end);
    const ArrayRCP<size_t> rowptr(ncoarse+1);
    rowptr[0] = 0;

    // set the crs pattern into the graph with local Indices
    rowptr[1] = rowptr[0]+2;    
    colind[0] = 0;
    colind[1] = 1;
    for(int rowIdx=1; rowIdx<ncoarse-1; rowIdx++) {
      int k = rowptr[rowIdx];
      rowptr[rowIdx+1] = rowptr[rowIdx]+3;
      colind[k]   = rowIdx-1;
      colind[k+1] = rowIdx;
      colind[k+2] = rowIdx+1;
    }
    rowptr[ncoarse] = rowptr[ncoarse-1]+2;
    colind[end-2] = ncoarse-2;
    colind[end-1] = ncoarse-1;
    
    GetOStream(Statistics2) << "StructuredRAP: Graph indices created!\n";
    myGraph->setAllIndices(rowptr, colind);

    GetOStream(Statistics2) << "StructuredRAP: Graph is created and filled!\n";
    myGraph->fillComplete();
 
    // build Ac with static graph pattern ...
    Ac = MatrixFactory::Build(myGraph, paramList); 
  }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  const bool doTranspose       = true;
  const bool doFillComplete    = true;
  const bool doOptimizeStorage = true;
  RCP<Matrix> Ac;

  {
         
      Teuchos::Array<LocalOrdinal> lCoarseNodesPerDim(3);
      lCoarseNodesPerDim = Get<Teuchos::Array<LocalOrdinal>>(fineLevel, "lCoarseNodesPerDim");
  
    FactoryMonitor m(*this, "Computing Ac", coarseLevel);
      std::ostringstream levelstr;
      levelstr << coarseLevel.GetLevelID();
      std::string labelstr = FormattingHelper::getColonLabel(coarseLevel.getObjectLabel());

    TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_ == false, Exceptions::RuntimeError,
                               "MueLu::RAPFactory::Build(): CallDeclareInput has not been called before Build!");

      const Teuchos::ParameterList& pL = GetParameterList();
      RCP<Matrix> A = Get< RCP<Matrix> >(fineLevel,   "A");
      RCP<Matrix> P = Get< RCP<Matrix> >(coarseLevel, "P");
      RCP<Matrix> AP;

      // bool isTpetra = A->getRowMap()->lib() == Xpetra::UseTpetra;
#ifdef KOKKOS_ENABLE_CUDA
      bool isCuda = typeid(Node).name() == typeid(Kokkos::Compat::KokkosCudaWrapperNode).name();
#else
      bool isCuda = false;
#endif

      if (pL.get<bool>("rap: triple product") == false || isCuda) {
        
#ifdef KOKKOS_ENABLE_CUDA
          GetOStream(Warnings1) << "Triple product R x A x P has not been implemented for Cuda.\n";
#endif
        GetOStream(Warnings1) << "For the structured RAP Factory switch to triple product R x A x P.\n";
        // For now, force hard exit
        exit(1);

      } 
      else {
            
        RCP<ParameterList> RAPparams = rcp(new ParameterList);
        if(pL.isSublist("matrixmatrix: kernel params"))
          RAPparams->sublist("matrixmatrix: kernel params") = pL.sublist("matrixmatrix: kernel params");

    if (coarseLevel.IsAvailable("RAP reuse data", this)) {
      GetOStream(static_cast<MsgType>(Runtime0 | Test)) << "Reusing previous RAP data" << std::endl;
      
          RAPparams = coarseLevel.Get< RCP<ParameterList> >("RAP reuse data", this);

          if (RAPparams->isParameter("graph"))
            Ac = RAPparams->get< RCP<Matrix> >("graph");

          // Some eigenvalue may have been cached with the matrix in the previous run.
          // As the matrix values will be updated, we need to reset the eigenvalue.
          Ac->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<SC>::one());
        
    } else {
      // if reuse data not available, try to get sparse fill graph via the knowledge of the matrix structure
          // Here we should technically also ask for the corsening rate / interpolation order
          std::string structureType = pL.get<std::string>("rap: structure type");
          if(structureType=="Laplace1D")
          {
            GetOStream(Statistics2) << "StructuredRAP: Using Laplace1D pattern determination routine.\n";
            GetLaplace1D(Ac, P, lCoarseNodesPerDim);
          }
        }
        
        // We *always* need global constants for the RAP, but not for the temps
        RAPparams->set("compute global constants: temporaries",RAPparams->get("compute global constants: temporaries",false));
        RAPparams->set("compute global constants",true);

        if (pL.get<bool>("transpose: use implicit") == true) {

          SubFactoryMonitor m2(*this, "MxMxM: R x A x P (implicit)", coarseLevel);

          Xpetra::TripleMatrixMultiply<SC,LO,GO,NO>::
            MultiplyRAP(*P, doTranspose, *A, !doTranspose, *P, !doTranspose, *Ac, doFillComplete,
                        doOptimizeStorage, labelstr+std::string("MueLu::R*A*P-implicit-")+levelstr.str(),
                        RAPparams);
            
        } else {
          
          RCP<Matrix> R = Get< RCP<Matrix> >(coarseLevel, "R");

          SubFactoryMonitor m2(*this, "MxMxM: R x A x P (explicit)", coarseLevel);

          Xpetra::TripleMatrixMultiply<SC,LO,GO,NO>::
            MultiplyRAP(*R, !doTranspose, *A, !doTranspose, *P, !doTranspose, *Ac, doFillComplete,
                        doOptimizeStorage, labelstr+std::string("MueLu::R*A*P-explicit-")+levelstr.str(),
                        RAPparams);
        }
      
        Teuchos::ArrayView<const double> relativeFloor = pL.get<Teuchos::Array<double> >("rap: relative diagonal floor")();
        if(relativeFloor.size() > 0) {
          Xpetra::MatrixUtils<SC,LO,GO,NO>::RelativeDiagonalBoost(Ac, relativeFloor,GetOStream(Statistics2));
        }

        bool repairZeroDiagonals = pL.get<bool>("RepairMainDiagonal") || pL.get<bool>("rap: fix zero diagonals");
        bool checkAc             = pL.get<bool>("CheckMainDiagonal")|| pL.get<bool>("rap: fix zero diagonals"); ;
        if (checkAc || repairZeroDiagonals) {
          using magnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
          magnitudeType threshold;
          if (pL.isType<magnitudeType>("rap: fix zero diagonals threshold"))
            threshold = pL.get<magnitudeType>("rap: fix zero diagonals threshold");
          else
            threshold = Teuchos::as<magnitudeType>(pL.get<double>("rap: fix zero diagonals threshold"));
          Scalar replacement = Teuchos::as<Scalar>(pL.get<double>("rap: fix zero diagonals replacement"));
          Xpetra::MatrixUtils<SC,LO,GO,NO>::CheckRepairMainDiagonal(Ac, repairZeroDiagonals, GetOStream(Warnings1), threshold, replacement);
      }


        if (IsPrint(Statistics2)) {
          RCP<ParameterList> params = rcp(new ParameterList());;
          params->set("printLoadBalancingInfo", true);
          params->set("printCommInfo",          true);
          GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*Ac, "Ac", params);
        }

        if(!Ac.is_null()) {std::ostringstream oss; oss << "A_" << coarseLevel.GetLevelID(); Ac->setObjectLabel(oss.str());}
        Set(coarseLevel, "A",         Ac);

        RAPparams->set("graph", Ac);
        Set(coarseLevel, "RAP reuse data", RAPparams);
      }
    }

  if (Behavior::debug())
    MatrixUtils::checkLocalRowMapMatchesColMap(*Ac);

  if (transferFacts_.begin() != transferFacts_.end()) {
    SubFactoryMonitor m(*this, "Projections", coarseLevel);

    // call Build of all user-given transfer factories
    for (std::vector<RCP<const FactoryBase> >::const_iterator it = transferFacts_.begin(); it != transferFacts_.end(); ++it) {
      RCP<const FactoryBase> fac = *it;
      GetOStream(Runtime0) << "RAPFactory: call transfer factory: " << fac->description() << std::endl;
      fac->CallBuild(coarseLevel);
      // Coordinates transfer is marginally different from all other operations
      // because it is *optional*, and not required. For instance, we may need
      // coordinates only on level 4 if we start repartitioning from that level,
      // but we don't need them on level 1,2,3. As our current Hierarchy setup
      // assumes propagation of dependencies only through three levels, this
      // means that we need to rely on other methods to propagate optional data.
      //
      // The method currently used is through RAP transfer factories, which are
      // simply factories which are called at the end of RAP with a single goal:
      // transfer some fine data to coarser level. Because these factories are
      // kind of outside of the mainline factories, they behave different. In
      // particular, we call their Build method explicitly, rather than through
      // Get calls. This difference is significant, as the Get call is smart
      // enough to know when to release all factory dependencies, and Build is
      // dumb. This led to the following CoordinatesTransferFactory sequence:
      // 1. Request level 0
      // 2. Request level 1
      // 3. Request level 0
      // 4. Release level 0
      // 5. Release level 1
      //
      // The problem is missing "6. Release level 0". Because it was missing,
      // we had outstanding request on "Coordinates", "Aggregates" and
      // "CoarseMap" on level 0.
      //
      // This was fixed by explicitly calling Release on transfer factories in
      // RAPFactory. I am still unsure how exactly it works, but now we have
      // clear data requests for all levels.
      coarseLevel.Release(*fac);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void StructuredRAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::AddTransferFactory(const RCP<const FactoryBase>& factory) {
  // check if it's a TwoLevelFactoryBase based transfer factory
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const TwoLevelFactoryBase>(factory) == Teuchos::null, Exceptions::BadCast,
                             "MueLu::RAPFactory::AddTransferFactory: Transfer factory is not derived from TwoLevelFactoryBase. "
                             "This is very strange. (Note: you can remove this exception if there's a good reason for)");
  TEUCHOS_TEST_FOR_EXCEPTION(hasDeclaredInput_, Exceptions::RuntimeError, "MueLu::RAPFactory::AddTransferFactory: Factory is being added after we have already declared input");
  transferFacts_.push_back(factory);
}

}  // namespace MueLu

#define MUELU_STRUCTUREDRAPFACTORY_SHORT
#endif // MUELU_STRUCTUREDRAPFACTORY_DEF_HPP
