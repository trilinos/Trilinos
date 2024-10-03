// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_COALESCEDROPFACTORY_DEF_HPP
#define MUELU_COALESCEDROPFACTORY_DEF_HPP

#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_StridedMap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>

#include <Xpetra_IO.hpp>

#include "MueLu_CoalesceDropFactory_decl.hpp"

#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_LWGraph.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_LWGraph.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"
#include "MueLu_Utilities.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Tpetra_CrsGraphTransposer.hpp"
#endif

#include <algorithm>
#include <cstdlib>
#include <string>

// If defined, read environment variables.
// Should be removed once we are confident that this works.
//#define DJS_READ_ENV_VARIABLES

namespace MueLu {

namespace Details {
template <class real_type, class LO>
struct DropTol {
  DropTol()               = default;
  DropTol(DropTol const&) = default;
  DropTol(DropTol&&)      = default;

  DropTol& operator=(DropTol const&) = default;
  DropTol& operator=(DropTol&&)      = default;

  DropTol(real_type val_, real_type diag_, LO col_, bool drop_)
    : val{val_}
    , diag{diag_}
    , col{col_}
    , drop{drop_} {}

  real_type val{Teuchos::ScalarTraits<real_type>::zero()};
  real_type diag{Teuchos::ScalarTraits<real_type>::zero()};
  LO col{Teuchos::OrdinalTraits<LO>::invalid()};
  bool drop{true};

  // CMS: Auxillary information for debugging info
  //      real_type aux_val {Teuchos::ScalarTraits<real_type>::nan()};
};
}  // namespace Details

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("aggregation: drop tol");
  SET_VALID_ENTRY("aggregation: use ml scaling of drop tol");
  SET_VALID_ENTRY("aggregation: Dirichlet threshold");
  SET_VALID_ENTRY("aggregation: greedy Dirichlet");
  SET_VALID_ENTRY("aggregation: row sum drop tol");
  SET_VALID_ENTRY("aggregation: drop scheme");
  SET_VALID_ENTRY("aggregation: block diagonal: interleaved blocksize");
  SET_VALID_ENTRY("aggregation: distance laplacian directional weights");
  SET_VALID_ENTRY("aggregation: dropping may create Dirichlet");

  {
    // "signed classical" is the Ruge-Stuben style (relative to max off-diagonal), "sign classical sa" is the signed version of the sa criterion (relative to the diagonal values)
    validParamList->getEntry("aggregation: drop scheme").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("signed classical sa", "classical", "distance laplacian", "signed classical", "block diagonal", "block diagonal classical", "block diagonal distance laplacian", "block diagonal signed classical", "block diagonal colored signed classical"))));
  }
  SET_VALID_ENTRY("aggregation: distance laplacian algo");
  SET_VALID_ENTRY("aggregation: classical algo");
  SET_VALID_ENTRY("aggregation: coloring: localize color graph");
#undef SET_VALID_ENTRY
  validParamList->set<bool>("lightweight wrap", true, "Experimental option for lightweight graph access");

  validParamList->set<RCP<const FactoryBase>>("A", Teuchos::null, "Generating factory of the matrix A");
  validParamList->set<RCP<const FactoryBase>>("UnAmalgamationInfo", Teuchos::null, "Generating factory for UnAmalgamationInfo");
  validParamList->set<RCP<const FactoryBase>>("Coordinates", Teuchos::null, "Generating factory for Coordinates");
  validParamList->set<RCP<const FactoryBase>>("BlockNumber", Teuchos::null, "Generating factory for BlockNUmber");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CoalesceDropFactory()
  : predrop_(Teuchos::null) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
  Input(currentLevel, "UnAmalgamationInfo");

  const ParameterList& pL = GetParameterList();
  if (pL.get<bool>("lightweight wrap") == true) {
    std::string algo = pL.get<std::string>("aggregation: drop scheme");
    if (algo == "distance laplacian" || algo == "block diagonal distance laplacian") {
      Input(currentLevel, "Coordinates");
    }
    if (algo == "signed classical sa")
      ;
    else if (algo.find("block diagonal") != std::string::npos || algo.find("signed classical") != std::string::npos) {
      Input(currentLevel, "BlockNumber");
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Build", currentLevel);

  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;
  typedef Xpetra::MultiVectorFactory<real_type, LO, GO, NO> RealValuedMultiVectorFactory;

  if (predrop_ != Teuchos::null)
    GetOStream(Parameters0) << predrop_->description();

  RCP<Matrix> realA              = Get<RCP<Matrix>>(currentLevel, "A");
  RCP<AmalgamationInfo> amalInfo = Get<RCP<AmalgamationInfo>>(currentLevel, "UnAmalgamationInfo");
  const ParameterList& pL        = GetParameterList();
  bool doExperimentalWrap        = pL.get<bool>("lightweight wrap");

  GetOStream(Parameters0) << "lightweight wrap = " << doExperimentalWrap << std::endl;
  std::string algo                         = pL.get<std::string>("aggregation: drop scheme");
  const bool aggregationMayCreateDirichlet = pL.get<bool>("aggregation: dropping may create Dirichlet");

  RCP<RealValuedMultiVector> Coords;
  RCP<Matrix> A;

  bool use_block_algorithm   = false;
  LO interleaved_blocksize   = as<LO>(pL.get<int>("aggregation: block diagonal: interleaved blocksize"));
  bool useSignedClassicalRS  = false;
  bool useSignedClassicalSA  = false;
  bool generateColoringGraph = false;

  // NOTE:  If we're doing blockDiagonal, we'll not want to do rowSum twice (we'll do it
  // in the block diagonalization). So we'll clobber the rowSumTol with -1.0 in this case
  typename STS::magnitudeType rowSumTol = as<typename STS::magnitudeType>(pL.get<double>("aggregation: row sum drop tol"));

  RCP<LocalOrdinalVector> ghostedBlockNumber;

  if (algo == "distance laplacian") {
    // Grab the coordinates for distance laplacian
    Coords = Get<RCP<RealValuedMultiVector>>(currentLevel, "Coordinates");
    A      = realA;
  } else if (algo == "signed classical sa") {
    useSignedClassicalSA = true;
    algo                 = "classical";
    A                    = realA;
  } else if (algo == "signed classical" || algo == "block diagonal colored signed classical" || algo == "block diagonal signed classical") {
    useSignedClassicalRS = true;
    //      if(realA->GetFixedBlockSize() > 1) {
    RCP<LocalOrdinalVector> BlockNumber = Get<RCP<LocalOrdinalVector>>(currentLevel, "BlockNumber");
    // Ghost the column block numbers if we need to
    RCP<const Import> importer = realA->getCrsGraph()->getImporter();
    if (!importer.is_null()) {
      SubFactoryMonitor m1(*this, "Block Number import", currentLevel);
      ghostedBlockNumber = Xpetra::VectorFactory<LO, LO, GO, NO>::Build(importer->getTargetMap());
      ghostedBlockNumber->doImport(*BlockNumber, *importer, Xpetra::INSERT);
    } else {
      ghostedBlockNumber = BlockNumber;
    }
    //      }
    if (algo == "block diagonal colored signed classical")
      generateColoringGraph = true;
    algo = "classical";
    A    = realA;

  } else if (algo == "block diagonal") {
    // Handle the "block diagonal" filtering and then leave
    BlockDiagonalize(currentLevel, realA, false);
    return;
  } else if (algo == "block diagonal classical" || algo == "block diagonal distance laplacian") {
    // Handle the "block diagonal" filtering, and then continue onward
    use_block_algorithm        = true;
    RCP<Matrix> filteredMatrix = BlockDiagonalize(currentLevel, realA, true);
    if (algo == "block diagonal distance laplacian") {
      // We now need to expand the coordinates by the interleaved blocksize
      RCP<RealValuedMultiVector> OldCoords = Get<RCP<RealValuedMultiVector>>(currentLevel, "Coordinates");
      if (OldCoords->getLocalLength() != realA->getLocalNumRows()) {
        LO dim = (LO)OldCoords->getNumVectors();
        Coords = RealValuedMultiVectorFactory::Build(realA->getRowMap(), dim);
        for (LO k = 0; k < dim; k++) {
          ArrayRCP<const real_type> old_vec = OldCoords->getData(k);
          ArrayRCP<real_type> new_vec       = Coords->getDataNonConst(k);
          for (LO i = 0; i < (LO)OldCoords->getLocalLength(); i++) {
            LO new_base = i * dim;
            for (LO j = 0; j < interleaved_blocksize; j++)
              new_vec[new_base + j] = old_vec[i];
          }
        }
      } else {
        Coords = OldCoords;
      }
      algo = "distance laplacian";
    } else if (algo == "block diagonal classical") {
      algo = "classical";
    }
    // All cases
    A         = filteredMatrix;
    rowSumTol = -1.0;
  } else {
    A = realA;
  }

  // Distance Laplacian weights
  Array<double> dlap_weights = pL.get<Array<double>>("aggregation: distance laplacian directional weights");
  enum { NO_WEIGHTS = 0,
         SINGLE_WEIGHTS,
         BLOCK_WEIGHTS };
  int use_dlap_weights = NO_WEIGHTS;
  if (algo == "distance laplacian") {
    LO dim = (LO)Coords->getNumVectors();
    // If anything isn't 1.0 we need to turn on the weighting
    bool non_unity = false;
    for (LO i = 0; !non_unity && i < (LO)dlap_weights.size(); i++) {
      if (dlap_weights[i] != 1.0) {
        non_unity = true;
      }
    }
    if (non_unity) {
      LO blocksize = use_block_algorithm ? as<LO>(pL.get<int>("aggregation: block diagonal: interleaved blocksize")) : 1;
      if ((LO)dlap_weights.size() == dim)
        use_dlap_weights = SINGLE_WEIGHTS;
      else if ((LO)dlap_weights.size() == blocksize * dim)
        use_dlap_weights = BLOCK_WEIGHTS;
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(1, Exceptions::RuntimeError,
                                   "length of 'aggregation: distance laplacian directional weights' must equal the coordinate dimension OR the coordinate dimension times the blocksize");
      }
      if (GetVerbLevel() & Statistics1)
        GetOStream(Statistics1) << "Using distance laplacian weights: " << dlap_weights << std::endl;
    }
  }

  // decide wether to use the fast-track code path for standard maps or the somewhat slower
  // code path for non-standard maps
  /*bool bNonStandardMaps = false;
  if (A->IsView("stridedMaps") == true) {
    Teuchos::RCP<const Map> myMap = A->getRowMap("stridedMaps");
    Teuchos::RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(myMap);
    TEUCHOS_TEST_FOR_EXCEPTION(strMap == null, Exceptions::RuntimeError, "Map is not of type StridedMap");
    if (strMap->getStridedBlockId() != -1 || strMap->getOffset() > 0)
      bNonStandardMaps = true;
  }*/

  if (doExperimentalWrap) {
    TEUCHOS_TEST_FOR_EXCEPTION(predrop_ != null && algo != "classical", Exceptions::RuntimeError, "Dropping function must not be provided for \"" << algo << "\" algorithm");
    TEUCHOS_TEST_FOR_EXCEPTION(algo != "classical" && algo != "distance laplacian" && algo != "signed classical", Exceptions::RuntimeError, "\"algorithm\" must be one of (classical|distance laplacian|signed classical)");

    SC threshold;
    // If we're doing the ML-style halving of the drop tol at each level, we do that here.
    if (pL.get<bool>("aggregation: use ml scaling of drop tol"))
      threshold = pL.get<double>("aggregation: drop tol") / pow(2.0, currentLevel.GetLevelID());
    else
      threshold = as<SC>(pL.get<double>("aggregation: drop tol"));

    std::string distanceLaplacianAlgoStr = pL.get<std::string>("aggregation: distance laplacian algo");
    std::string classicalAlgoStr         = pL.get<std::string>("aggregation: classical algo");
    real_type realThreshold              = STS::magnitude(threshold);  // CMS: Rename this to "magnitude threshold" sometime

    ////////////////////////////////////////////////////
    // Remove this bit once we are confident that cut-based dropping works.
#ifdef HAVE_MUELU_DEBUG
    int distanceLaplacianCutVerbose = 0;
#endif
#ifdef DJS_READ_ENV_VARIABLES
    if (getenv("MUELU_DROP_TOLERANCE_MODE")) {
      distanceLaplacianAlgoStr = std::string(getenv("MUELU_DROP_TOLERANCE_MODE"));
    }

    if (getenv("MUELU_DROP_TOLERANCE_THRESHOLD")) {
      auto tmp      = atoi(getenv("MUELU_DROP_TOLERANCE_THRESHOLD"));
      realThreshold = 1e-4 * tmp;
    }

#ifdef HAVE_MUELU_DEBUG
    if (getenv("MUELU_DROP_TOLERANCE_VERBOSE")) {
      distanceLaplacianCutVerbose = atoi(getenv("MUELU_DROP_TOLERANCE_VERBOSE"));
    }
#endif
#endif
    ////////////////////////////////////////////////////

    enum decisionAlgoType { defaultAlgo,
                            unscaled_cut,
                            scaled_cut,
                            scaled_cut_symmetric };

    decisionAlgoType distanceLaplacianAlgo = defaultAlgo;
    decisionAlgoType classicalAlgo         = defaultAlgo;
    if (algo == "distance laplacian") {
      if (distanceLaplacianAlgoStr == "default")
        distanceLaplacianAlgo = defaultAlgo;
      else if (distanceLaplacianAlgoStr == "unscaled cut")
        distanceLaplacianAlgo = unscaled_cut;
      else if (distanceLaplacianAlgoStr == "scaled cut")
        distanceLaplacianAlgo = scaled_cut;
      else if (distanceLaplacianAlgoStr == "scaled cut symmetric")
        distanceLaplacianAlgo = scaled_cut_symmetric;
      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "\"aggregation: distance laplacian algo\" must be one of (default|unscaled cut|scaled cut), not \"" << distanceLaplacianAlgoStr << "\"");
      GetOStream(Runtime0) << "algorithm = \"" << algo << "\" distance laplacian algorithm = \"" << distanceLaplacianAlgoStr << "\": threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << std::endl;
    } else if (algo == "classical") {
      if (classicalAlgoStr == "default")
        classicalAlgo = defaultAlgo;
      else if (classicalAlgoStr == "unscaled cut")
        classicalAlgo = unscaled_cut;
      else if (classicalAlgoStr == "scaled cut")
        classicalAlgo = scaled_cut;
      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "\"aggregation: classical algo\" must be one of (default|unscaled cut|scaled cut), not \"" << classicalAlgoStr << "\"");
      GetOStream(Runtime0) << "algorithm = \"" << algo << "\" classical algorithm = \"" << classicalAlgoStr << "\": threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << std::endl;

    } else
      GetOStream(Runtime0) << "algorithm = \"" << algo << "\": threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << std::endl;
    Set<bool>(currentLevel, "Filtering", (threshold != STS::zero()));

    const typename STS::magnitudeType dirichletThreshold = STS::magnitude(as<SC>(pL.get<double>("aggregation: Dirichlet threshold")));

    // NOTE: We don't support signed classical RS or SA with cut drop at present
    TEUCHOS_TEST_FOR_EXCEPTION(useSignedClassicalRS && classicalAlgo != defaultAlgo, Exceptions::RuntimeError, "\"aggregation: classical algo\" != default is not supported for scalled classical aggregation");
    TEUCHOS_TEST_FOR_EXCEPTION(useSignedClassicalSA && classicalAlgo != defaultAlgo, Exceptions::RuntimeError, "\"aggregation: classical algo\" != default is not supported for scalled classical sa aggregation");

    GO numDropped = 0, numTotal = 0;
    std::string graphType = "unamalgamated";  // for description purposes only

    /* NOTE: storageblocksize (from GetStorageBlockSize()) is the size of a block in the chosen storage scheme.
     BlockSize is the number of storage blocks that must kept together during the amalgamation process.

     Both of these quantities may be different than numPDEs (from GetFixedBlockSize()), but the following must always hold:

     numPDEs = BlockSize * storageblocksize.

     If numPDEs==1
       Matrix is point storage (classical CRS storage).  storageblocksize=1 and BlockSize=1
       No other values makes sense.

     If numPDEs>1
       If matrix uses point storage, then storageblocksize=1  and BlockSize=numPDEs.
       If matrix uses block storage, with block size of n, then storageblocksize=n, and BlockSize=numPDEs/n.
       Thus far, only storageblocksize=numPDEs and BlockSize=1 has been tested.
    */
    TEUCHOS_TEST_FOR_EXCEPTION(A->GetFixedBlockSize() % A->GetStorageBlockSize() != 0, Exceptions::RuntimeError, "A->GetFixedBlockSize() needs to be a multiple of A->GetStorageBlockSize()");
    const LO BlockSize = A->GetFixedBlockSize() / A->GetStorageBlockSize();

    /************************** RS or SA-style Classical Dropping (and variants) **************************/
    if (algo == "classical") {
      if (predrop_ == null) {
        // ap: this is a hack: had to declare predrop_ as mutable
        predrop_ = rcp(new PreDropFunctionConstVal(threshold));
      }

      if (predrop_ != null) {
        RCP<PreDropFunctionConstVal> predropConstVal = rcp_dynamic_cast<PreDropFunctionConstVal>(predrop_);
        TEUCHOS_TEST_FOR_EXCEPTION(predropConstVal == Teuchos::null, Exceptions::BadCast,
                                   "MueLu::CoalesceFactory::Build: cast to PreDropFunctionConstVal failed.");
        // If a user provided a predrop function, it overwrites the XML threshold parameter
        SC newt = predropConstVal->GetThreshold();
        if (newt != threshold) {
          GetOStream(Warnings0) << "switching threshold parameter from " << threshold << " (list) to " << newt << " (user function" << std::endl;
          threshold = newt;
        }
      }
      // At this points we either have
      //     (predrop_ != null)
      // Therefore, it is sufficient to check only threshold
      if (BlockSize == 1 && threshold == STS::zero() && !useSignedClassicalRS && !useSignedClassicalSA && A->hasCrsGraph()) {
        // Case 1:  scalar problem, no dropping => just use matrix graph
        RCP<LWGraph> graph = rcp(new LWGraph(A->getCrsGraph(), "graph of A"));
        // Detect and record rows that correspond to Dirichlet boundary conditions
        auto boundaryNodes = MueLu::Utilities<SC, LO, GO, NO>::DetectDirichletRows_kokkos_host(*A, dirichletThreshold);
        if (rowSumTol > 0.)
          Utilities::ApplyRowSumCriterionHost(*A, rowSumTol, boundaryNodes);

        graph->SetBoundaryNodeMap(boundaryNodes);
        numTotal = A->getLocalNumEntries();

        if (GetVerbLevel() & Statistics1) {
          GO numLocalBoundaryNodes  = 0;
          GO numGlobalBoundaryNodes = 0;
          for (size_t i = 0; i < boundaryNodes.size(); ++i)
            if (boundaryNodes[i])
              numLocalBoundaryNodes++;
          RCP<const Teuchos::Comm<int>> comm = A->getRowMap()->getComm();
          MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
          GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
        }

        Set(currentLevel, "DofsPerNode", 1);
        Set(currentLevel, "Graph", graph);

      } else if ((BlockSize == 1 && threshold != STS::zero()) ||
                 (BlockSize == 1 && threshold == STS::zero() && !A->hasCrsGraph()) ||
                 (BlockSize == 1 && useSignedClassicalRS) ||
                 (BlockSize == 1 && useSignedClassicalSA)) {
        // Case 2:  scalar problem with dropping => record the column indices of undropped entries, but still use original
        //                                          graph's map information, e.g., whether index is local
        // OR a matrix without a CrsGraph

        // allocate space for the local graph
        typename LWGraph::row_type::non_const_type rows("rows", A->getLocalNumRows() + 1);
        typename LWGraph::entries_type::non_const_type columns("columns", A->getLocalNumEntries());

        using MT = typename STS::magnitudeType;
        RCP<Vector> ghostedDiag;
        ArrayRCP<const SC> ghostedDiagVals;
        ArrayRCP<const MT> negMaxOffDiagonal;
        // RS style needs the max negative off-diagonal, SA style needs the diagonal
        if (useSignedClassicalRS) {
          if (ghostedBlockNumber.is_null()) {
            auto negMaxOffDiagonalVec = MueLu::Utilities<SC, LO, GO, NO>::GetMatrixMaxMinusOffDiagonal(*A);
            negMaxOffDiagonal         = negMaxOffDiagonalVec->getData(0);
            if (GetVerbLevel() & Statistics1)
              GetOStream(Statistics1) << "Calculated max point off-diagonal" << std::endl;
          } else {
            auto negMaxOffDiagonalVec = MueLu::Utilities<SC, LO, GO, NO>::GetMatrixMaxMinusOffDiagonal(*A, *ghostedBlockNumber);
            negMaxOffDiagonal         = negMaxOffDiagonalVec->getData(0);
            if (GetVerbLevel() & Statistics1)
              GetOStream(Statistics1) << "Calculating max block off-diagonal" << std::endl;
          }
        } else {
          ghostedDiag     = MueLu::Utilities<SC, LO, GO, NO>::GetMatrixOverlappedDiagonal(*A);
          ghostedDiagVals = ghostedDiag->getData(0);
        }
        auto boundaryNodes = MueLu::Utilities<SC, LO, GO, NO>::DetectDirichletRows_kokkos_host(*A, dirichletThreshold);
        if (rowSumTol > 0.) {
          if (ghostedBlockNumber.is_null()) {
            if (GetVerbLevel() & Statistics1)
              GetOStream(Statistics1) << "Applying point row sum criterion." << std::endl;
            Utilities::ApplyRowSumCriterionHost(*A, rowSumTol, boundaryNodes);
          } else {
            if (GetVerbLevel() & Statistics1)
              GetOStream(Statistics1) << "Applying block row sum criterion." << std::endl;
            Utilities::ApplyRowSumCriterionHost(*A, *ghostedBlockNumber, rowSumTol, boundaryNodes);
          }
        }

        ArrayRCP<const LO> g_block_id;
        if (!ghostedBlockNumber.is_null())
          g_block_id = ghostedBlockNumber->getData(0);

        LO realnnz = 0;
        rows(0)    = 0;
        for (LO row = 0; row < Teuchos::as<LO>(A->getRowMap()->getLocalNumElements()); ++row) {
          size_t nnz          = A->getNumEntriesInLocalRow(row);
          bool rowIsDirichlet = boundaryNodes[row];
          ArrayView<const LO> indices;
          ArrayView<const SC> vals;
          A->getLocalRowView(row, indices, vals);

          if (classicalAlgo == defaultAlgo) {
            // FIXME the current predrop function uses the following
            // FIXME    if(std::abs(vals[k]) > std::abs(threshold_) || grow == gcid )
            // FIXME but the threshold doesn't take into account the rows' diagonal entries
            // FIXME For now, hardwiring the dropping in here

            LO rownnz = 0;
            if (useSignedClassicalRS) {
              // Signed classical RS style
              for (LO colID = 0; colID < Teuchos::as<LO>(nnz); colID++) {
                LO col         = indices[colID];
                MT max_neg_aik = realThreshold * STS::real(negMaxOffDiagonal[row]);
                MT neg_aij     = -STS::real(vals[colID]);
                /*                  if(row==1326) printf("A(%d,%d) = %6.4e, block = (%d,%d) neg_aij = %6.4e max_neg_aik = %6.4e\n",row,col,vals[colID],
                                     g_block_id.is_null() ? -1 :  g_block_id[row],
                                     g_block_id.is_null() ? -1 :  g_block_id[col],
                                     neg_aij, max_neg_aik);*/
                if ((!rowIsDirichlet && (g_block_id.is_null() || g_block_id[row] == g_block_id[col]) && neg_aij > max_neg_aik) || row == col) {
                  columns[realnnz++] = col;
                  rownnz++;
                } else
                  numDropped++;
              }
              rows(row + 1) = realnnz;
            } else if (useSignedClassicalSA) {
              // Signed classical SA style
              for (LO colID = 0; colID < Teuchos::as<LO>(nnz); colID++) {
                LO col = indices[colID];

                bool is_nonpositive = STS::real(vals[colID]) <= 0;
                MT aiiajj           = STS::magnitude(threshold * threshold * ghostedDiagVals[col] * ghostedDiagVals[row]);                        // eps^2*|a_ii|*|a_jj|
                MT aij              = is_nonpositive ? STS::magnitude(vals[colID] * vals[colID]) : (-STS::magnitude(vals[colID] * vals[colID]));  // + |a_ij|^2, if a_ij < 0, - |a_ij|^2 if a_ij >=0
                /*
                if(row==1326) printf("A(%d,%d) = %6.4e, raw_aij = %6.4e aij = %6.4e aiiajj = %6.4e\n",row,col,vals[colID],
                                     vals[colID],aij, aiiajj);
                */

                if ((!rowIsDirichlet && aij > aiiajj) || row == col) {
                  columns(realnnz++) = col;
                  rownnz++;
                } else
                  numDropped++;
              }
              rows[row + 1] = realnnz;
            } else {
              // Standard abs classical
              for (LO colID = 0; colID < Teuchos::as<LO>(nnz); colID++) {
                LO col    = indices[colID];
                MT aiiajj = STS::magnitude(threshold * threshold * ghostedDiagVals[col] * ghostedDiagVals[row]);  // eps^2*|a_ii|*|a_jj|
                MT aij    = STS::magnitude(vals[colID] * vals[colID]);                                            // |a_ij|^2

                if ((!rowIsDirichlet && aij > aiiajj) || row == col) {
                  columns(realnnz++) = col;
                  rownnz++;
                } else
                  numDropped++;
              }
              rows(row + 1) = realnnz;
            }
          } else {
            /* Cut Algorithm */
            // CMS
            using DropTol = Details::DropTol<real_type, LO>;
            std::vector<DropTol> drop_vec;
            drop_vec.reserve(nnz);
            const real_type zero = Teuchos::ScalarTraits<real_type>::zero();
            const real_type one  = Teuchos::ScalarTraits<real_type>::one();
            LO rownnz            = 0;
            // NOTE: This probably needs to be fixed for rowsum

            // find magnitudes
            for (LO colID = 0; colID < (LO)nnz; colID++) {
              LO col = indices[colID];
              if (row == col) {
                drop_vec.emplace_back(zero, one, colID, false);
                continue;
              }

              // Don't aggregate boundaries
              if (boundaryNodes[colID]) continue;
              typename STS::magnitudeType aiiajj = STS::magnitude(threshold * threshold * ghostedDiagVals[col] * ghostedDiagVals[row]);  // eps^2*|a_ii|*|a_jj|
              typename STS::magnitudeType aij    = STS::magnitude(vals[colID] * vals[colID]);                                            // |a_ij|^2
              drop_vec.emplace_back(aij, aiiajj, colID, false);
            }

            const size_t n = drop_vec.size();

            if (classicalAlgo == unscaled_cut) {
              std::sort(drop_vec.begin(), drop_vec.end(), [](DropTol const& a, DropTol const& b) {
                return a.val > b.val;
              });

              bool drop = false;
              for (size_t i = 1; i < n; ++i) {
                if (!drop) {
                  auto const& x = drop_vec[i - 1];
                  auto const& y = drop_vec[i];
                  auto a        = x.val;
                  auto b        = y.val;
                  if (a > realThreshold * b) {
                    drop = true;
#ifdef HAVE_MUELU_DEBUG
                    if (distanceLaplacianCutVerbose) {
                      std::cout << "DJS: KEEP, N, ROW:  " << i + 1 << ", " << n << ", " << row << std::endl;
                    }
#endif
                  }
                }
                drop_vec[i].drop = drop;
              }
            } else if (classicalAlgo == scaled_cut) {
              std::sort(drop_vec.begin(), drop_vec.end(), [](DropTol const& a, DropTol const& b) {
                return a.val / a.diag > b.val / b.diag;
              });
              bool drop = false;
              //                  printf("[%d] Scaled Cut: ",(int)row);
              //                  printf("%3d(%4s) ",indices[drop_vec[0].col],"keep");
              for (size_t i = 1; i < n; ++i) {
                if (!drop) {
                  auto const& x = drop_vec[i - 1];
                  auto const& y = drop_vec[i];
                  auto a        = x.val / x.diag;
                  auto b        = y.val / y.diag;
                  if (a > realThreshold * b) {
                    drop = true;

#ifdef HAVE_MUELU_DEBUG
                    if (distanceLaplacianCutVerbose) {
                      std::cout << "DJS: KEEP, N, ROW:  " << i + 1 << ", " << n << ", " << row << std::endl;
                    }
#endif
                  }
                  //                      printf("%3d(%4s) ",indices[drop_vec[i].col],drop?"drop":"keep");
                }
                drop_vec[i].drop = drop;
              }
              //                  printf("\n");
            }
            std::sort(drop_vec.begin(), drop_vec.end(), [](DropTol const& a, DropTol const& b) {
              return a.col < b.col;
            });

            for (LO idxID = 0; idxID < (LO)drop_vec.size(); idxID++) {
              LO col = indices[drop_vec[idxID].col];
              // don't drop diagonal
              if (row == col) {
                columns[realnnz++] = col;
                rownnz++;
                continue;
              }

              if (!drop_vec[idxID].drop) {
                columns[realnnz++] = col;
                rownnz++;
              } else {
                numDropped++;
              }
            }
            // CMS
            rows[row + 1] = realnnz;
          }
        }  // end for row

        numTotal = A->getLocalNumEntries();

        if (aggregationMayCreateDirichlet) {
          // If the only element remaining after filtering is diagonal, mark node as boundary
          for (LO row = 0; row < Teuchos::as<LO>(A->getRowMap()->getLocalNumElements()); ++row) {
            if (rows[row + 1] - rows[row] <= 1)
              boundaryNodes[row] = true;
          }
        }

        RCP<LWGraph> graph = rcp(new LWGraph(rows, Kokkos::subview(columns, Kokkos::make_pair(0, realnnz)), A->getRowMap(), A->getColMap(), "thresholded graph of A"));
        graph->SetBoundaryNodeMap(boundaryNodes);
        if (GetVerbLevel() & Statistics1) {
          GO numLocalBoundaryNodes  = 0;
          GO numGlobalBoundaryNodes = 0;
          for (size_t i = 0; i < boundaryNodes.size(); ++i)
            if (boundaryNodes(i))
              numLocalBoundaryNodes++;
          RCP<const Teuchos::Comm<int>> comm = A->getRowMap()->getComm();
          MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
          GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
        }
        Set(currentLevel, "Graph", graph);
        Set(currentLevel, "DofsPerNode", 1);

        // If we're doing signed classical, we might want to block-diagonalize *after* the dropping
        if (generateColoringGraph) {
          RCP<LWGraph> colorGraph;
          RCP<const Import> importer = A->getCrsGraph()->getImporter();
          BlockDiagonalizeGraph(graph, ghostedBlockNumber, colorGraph, importer);
          Set(currentLevel, "Coloring Graph", colorGraph);
          // #define CMS_DUMP
#ifdef CMS_DUMP
          {
            Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("m_regular_graph." + std::to_string(currentLevel.GetLevelID()), *rcp_dynamic_cast<LWGraph>(graph)->GetCrsGraph());
            Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write("m_color_graph." + std::to_string(currentLevel.GetLevelID()), *rcp_dynamic_cast<LWGraph>(colorGraph)->GetCrsGraph());
            // int rank = graph->GetDomainMap()->getComm()->getRank();
            // {
            //   std::ofstream ofs(std::string("m_color_graph_") + std::to_string(currentLevel.GetLevelID())+std::string("_") + std::to_string(rank) + std::string(".dat"),std::ofstream::out);
            //   RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(ofs));
            //   colorGraph->print(*fancy,Debug);
            // }
            // {
            //   std::ofstream ofs(std::string("m_regular_graph_") + std::to_string(currentLevel.GetLevelID())+std::string("_") + std::to_string(rank) + std::string(".dat"),std::ofstream::out);
            //   RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(ofs));
            //   graph->print(*fancy,Debug);
            // }
          }
#endif
        }  // end generateColoringGraph
      } else if (BlockSize > 1 && threshold == STS::zero()) {
        // Case 3:  Multiple DOF/node problem without dropping
        const RCP<const Map> rowMap = A->getRowMap();
        const RCP<const Map> colMap = A->getColMap();

        graphType = "amalgamated";

        // build node row map (uniqueMap) and node column map (nonUniqueMap)
        // the arrays rowTranslation and colTranslation contain the local node id
        // given a local dof id. The data is calculated by the AmalgamationFactory and
        // stored in the variable container "UnAmalgamationInfo"
        RCP<const Map> uniqueMap    = amalInfo->getNodeRowMap();
        RCP<const Map> nonUniqueMap = amalInfo->getNodeColMap();
        Array<LO> rowTranslation    = *(amalInfo->getRowTranslation());
        Array<LO> colTranslation    = *(amalInfo->getColTranslation());

        // get number of local nodes
        LO numRows = Teuchos::as<LocalOrdinal>(uniqueMap->getLocalNumElements());

        // Allocate space for the local graph
        typename LWGraph::row_type::non_const_type rows("rows", numRows + 1);
        typename LWGraph::entries_type::non_const_type columns("columns", A->getLocalNumEntries());

        typename LWGraph::boundary_nodes_type amalgBoundaryNodes("amalgBoundaryNodes", numRows);
        Kokkos::deep_copy(amalgBoundaryNodes, false);

        // Detect and record rows that correspond to Dirichlet boundary conditions
        // TODO If we use ArrayRCP<LO>, then we can record boundary nodes as usual.  Size
        // TODO the array one bigger than the number of local rows, and the last entry can
        // TODO hold the actual number of boundary nodes.  Clever, huh?
        ArrayRCP<bool> pointBoundaryNodes;
        pointBoundaryNodes = Teuchos::arcp_const_cast<bool>(MueLu::Utilities<SC, LO, GO, NO>::DetectDirichletRows(*A, dirichletThreshold));
        if (rowSumTol > 0.)
          Utilities::ApplyRowSumCriterion(*A, rowSumTol, pointBoundaryNodes);

        // extract striding information
        LO blkSize     = A->GetFixedBlockSize();  //< the full block size (number of dofs per node in strided map)
        LO blkId       = -1;                      //< the block id within the strided map (or -1 if it is a full block map)
        LO blkPartSize = A->GetFixedBlockSize();  //< stores the size of the block within the strided map
        if (A->IsView("stridedMaps") == true) {
          Teuchos::RCP<const Map> myMap         = A->getRowMap("stridedMaps");
          Teuchos::RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(myMap);
          TEUCHOS_TEST_FOR_EXCEPTION(strMap == null, Exceptions::RuntimeError, "Map is not of type StridedMap");
          blkSize = Teuchos::as<const LO>(strMap->getFixedBlockSize());
          blkId   = strMap->getStridedBlockId();
          if (blkId > -1)
            blkPartSize = Teuchos::as<LO>(strMap->getStridingData()[blkId]);
        }

        // loop over all local nodes
        LO realnnz = 0;
        rows(0)    = 0;
        Array<LO> indicesExtra;
        for (LO row = 0; row < numRows; row++) {
          ArrayView<const LO> indices;
          indicesExtra.resize(0);

          // The amalgamated row is marked as Dirichlet iff all point rows are Dirichlet
          // Note, that pointBoundaryNodes lives on the dofmap (and not the node map).
          // Therefore, looping over all dofs is fine here. We use blkPartSize as we work
          // with local ids.
          // TODO: Here we have different options of how to define a node to be a boundary (or Dirichlet)
          // node.
          bool isBoundary = false;
          if (pL.get<bool>("aggregation: greedy Dirichlet") == true) {
            for (LO j = 0; j < blkPartSize; j++) {
              if (pointBoundaryNodes[row * blkPartSize + j]) {
                isBoundary = true;
                break;
              }
            }
          } else {
            isBoundary = true;
            for (LO j = 0; j < blkPartSize; j++) {
              if (!pointBoundaryNodes[row * blkPartSize + j]) {
                isBoundary = false;
                break;
              }
            }
          }

          // Merge rows of A
          // The array indicesExtra contains local column node ids for the current local node "row"
          if (!isBoundary)
            MergeRows(*A, row, indicesExtra, colTranslation);
          else
            indicesExtra.push_back(row);
          indices = indicesExtra;
          numTotal += indices.size();

          // add the local column node ids to the full columns array which
          // contains the local column node ids for all local node rows
          LO nnz = indices.size(), rownnz = 0;
          for (LO colID = 0; colID < nnz; colID++) {
            LO col             = indices[colID];
            columns(realnnz++) = col;
            rownnz++;
          }

          if (rownnz == 1) {
            // If the only element remaining after filtering is diagonal, mark node as boundary
            // FIXME: this should really be replaced by the following
            //    if (indices.size() == 1 && indices[0] == row)
            //        boundaryNodes[row] = true;
            // We do not do it this way now because there is no framework for distinguishing isolated
            // and boundary nodes in the aggregation algorithms
            amalgBoundaryNodes[row] = true;
          }
          rows(row + 1) = realnnz;
        }  // for (LO row = 0; row < numRows; row++)

        RCP<LWGraph> graph = rcp(new LWGraph(rows, Kokkos::subview(columns, Kokkos::make_pair(0, realnnz)), uniqueMap, nonUniqueMap, "amalgamated graph of A"));
        graph->SetBoundaryNodeMap(amalgBoundaryNodes);

        if (GetVerbLevel() & Statistics1) {
          GO numLocalBoundaryNodes  = 0;
          GO numGlobalBoundaryNodes = 0;

          for (size_t i = 0; i < amalgBoundaryNodes.size(); ++i)
            if (amalgBoundaryNodes(i))
              numLocalBoundaryNodes++;

          RCP<const Teuchos::Comm<int>> comm = A->getRowMap()->getComm();
          MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
          GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes
                                  << " agglomerated Dirichlet nodes" << std::endl;
        }

        Set(currentLevel, "Graph", graph);
        Set(currentLevel, "DofsPerNode", blkSize);  // full block size

      } else if (BlockSize > 1 && threshold != STS::zero()) {
        // Case 4:  Multiple DOF/node problem with dropping
        const RCP<const Map> rowMap = A->getRowMap();
        const RCP<const Map> colMap = A->getColMap();
        graphType                   = "amalgamated";

        // build node row map (uniqueMap) and node column map (nonUniqueMap)
        // the arrays rowTranslation and colTranslation contain the local node id
        // given a local dof id. The data is calculated by the AmalgamationFactory and
        // stored in the variable container "UnAmalgamationInfo"
        RCP<const Map> uniqueMap    = amalInfo->getNodeRowMap();
        RCP<const Map> nonUniqueMap = amalInfo->getNodeColMap();
        Array<LO> rowTranslation    = *(amalInfo->getRowTranslation());
        Array<LO> colTranslation    = *(amalInfo->getColTranslation());

        // get number of local nodes
        LO numRows = Teuchos::as<LocalOrdinal>(uniqueMap->getLocalNumElements());

        // Allocate space for the local graph
        typename LWGraph::row_type::non_const_type rows("rows", numRows + 1);
        typename LWGraph::entries_type::non_const_type columns("columns", A->getLocalNumEntries());

        typename LWGraph::boundary_nodes_type amalgBoundaryNodes("amalgBoundaryNodes", numRows);
        Kokkos::deep_copy(amalgBoundaryNodes, false);

        // Detect and record rows that correspond to Dirichlet boundary conditions
        // TODO If we use ArrayRCP<LO>, then we can record boundary nodes as usual.  Size
        // TODO the array one bigger than the number of local rows, and the last entry can
        // TODO hold the actual number of boundary nodes.  Clever, huh?
        auto pointBoundaryNodes = MueLu::Utilities<SC, LO, GO, NO>::DetectDirichletRows_kokkos_host(*A, dirichletThreshold);
        if (rowSumTol > 0.)
          Utilities::ApplyRowSumCriterionHost(*A, rowSumTol, pointBoundaryNodes);

        // extract striding information
        LO blkSize     = A->GetFixedBlockSize();  //< the full block size (number of dofs per node in strided map)
        LO blkId       = -1;                      //< the block id within the strided map (or -1 if it is a full block map)
        LO blkPartSize = A->GetFixedBlockSize();  //< stores the size of the block within the strided map
        if (A->IsView("stridedMaps") == true) {
          Teuchos::RCP<const Map> myMap         = A->getRowMap("stridedMaps");
          Teuchos::RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(myMap);
          TEUCHOS_TEST_FOR_EXCEPTION(strMap == null, Exceptions::RuntimeError, "Map is not of type StridedMap");
          blkSize = Teuchos::as<const LO>(strMap->getFixedBlockSize());
          blkId   = strMap->getStridedBlockId();
          if (blkId > -1)
            blkPartSize = Teuchos::as<LO>(strMap->getStridingData()[blkId]);
        }

        // extract diagonal data for dropping strategy
        RCP<Vector> ghostedDiag                  = MueLu::Utilities<SC, LO, GO, NO>::GetMatrixOverlappedDiagonal(*A);
        const ArrayRCP<const SC> ghostedDiagVals = ghostedDiag->getData(0);

        // loop over all local nodes
        LO realnnz = 0;
        rows[0]    = 0;
        Array<LO> indicesExtra;
        for (LO row = 0; row < numRows; row++) {
          ArrayView<const LO> indices;
          indicesExtra.resize(0);

          // The amalgamated row is marked as Dirichlet iff all point rows are Dirichlet
          // Note, that pointBoundaryNodes lives on the dofmap (and not the node map).
          // Therefore, looping over all dofs is fine here. We use blkPartSize as we work
          // with local ids.
          // TODO: Here we have different options of how to define a node to be a boundary (or Dirichlet)
          // node.
          bool isBoundary = false;
          if (pL.get<bool>("aggregation: greedy Dirichlet") == true) {
            for (LO j = 0; j < blkPartSize; j++) {
              if (pointBoundaryNodes[row * blkPartSize + j]) {
                isBoundary = true;
                break;
              }
            }
          } else {
            isBoundary = true;
            for (LO j = 0; j < blkPartSize; j++) {
              if (!pointBoundaryNodes[row * blkPartSize + j]) {
                isBoundary = false;
                break;
              }
            }
          }

          // Merge rows of A
          // The array indicesExtra contains local column node ids for the current local node "row"
          if (!isBoundary)
            MergeRowsWithDropping(*A, row, ghostedDiagVals, threshold, indicesExtra, colTranslation);
          else
            indicesExtra.push_back(row);
          indices = indicesExtra;
          numTotal += indices.size();

          // add the local column node ids to the full columns array which
          // contains the local column node ids for all local node rows
          LO nnz = indices.size(), rownnz = 0;
          for (LO colID = 0; colID < nnz; colID++) {
            LO col             = indices[colID];
            columns[realnnz++] = col;
            rownnz++;
          }

          if (rownnz == 1) {
            // If the only element remaining after filtering is diagonal, mark node as boundary
            // FIXME: this should really be replaced by the following
            //    if (indices.size() == 1 && indices[0] == row)
            //        boundaryNodes[row] = true;
            // We do not do it this way now because there is no framework for distinguishing isolated
            // and boundary nodes in the aggregation algorithms
            amalgBoundaryNodes[row] = true;
          }
          rows[row + 1] = realnnz;
        }  // for (LO row = 0; row < numRows; row++)
        // columns.resize(realnnz);

        RCP<LWGraph> graph = rcp(new LWGraph(rows, Kokkos::subview(columns, Kokkos::make_pair(0, realnnz)), uniqueMap, nonUniqueMap, "amalgamated graph of A"));
        graph->SetBoundaryNodeMap(amalgBoundaryNodes);

        if (GetVerbLevel() & Statistics1) {
          GO numLocalBoundaryNodes  = 0;
          GO numGlobalBoundaryNodes = 0;

          for (size_t i = 0; i < amalgBoundaryNodes.size(); ++i)
            if (amalgBoundaryNodes(i))
              numLocalBoundaryNodes++;

          RCP<const Teuchos::Comm<int>> comm = A->getRowMap()->getComm();
          MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
          GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes
                                  << " agglomerated Dirichlet nodes" << std::endl;
        }

        Set(currentLevel, "Graph", graph);
        Set(currentLevel, "DofsPerNode", blkSize);  // full block size
      }

    } else if (algo == "distance laplacian") {
      LO blkSize   = A->GetFixedBlockSize();
      GO indexBase = A->getRowMap()->getIndexBase();
      // [*0*] : FIXME
      // ap: somehow, if I move this line to [*1*], Belos throws an error
      // I'm not sure what's going on. Do we always have to Get data, if we did
      // DeclareInput for it?
      //        RCP<RealValuedMultiVector> Coords = Get< RCP<RealValuedMultiVector > >(currentLevel, "Coordinates");

      // Detect and record rows that correspond to Dirichlet boundary conditions
      // TODO If we use ArrayRCP<LO>, then we can record boundary nodes as usual.  Size
      // TODO the array one bigger than the number of local rows, and the last entry can
      // TODO hold the actual number of boundary nodes.  Clever, huh?
      auto pointBoundaryNodes = MueLu::Utilities<SC, LO, GO, NO>::DetectDirichletRows_kokkos_host(*A, dirichletThreshold);
      if (rowSumTol > 0.)
        Utilities::ApplyRowSumCriterionHost(*A, rowSumTol, pointBoundaryNodes);

      if ((blkSize == 1) && (threshold == STS::zero())) {
        // Trivial case: scalar problem, no dropping. Can return original graph
        RCP<LWGraph> graph = rcp(new LWGraph(A->getCrsGraph(), "graph of A"));
        graph->SetBoundaryNodeMap(pointBoundaryNodes);
        graphType = "unamalgamated";
        numTotal  = A->getLocalNumEntries();

        if (GetVerbLevel() & Statistics1) {
          GO numLocalBoundaryNodes  = 0;
          GO numGlobalBoundaryNodes = 0;
          for (size_t i = 0; i < pointBoundaryNodes.size(); ++i)
            if (pointBoundaryNodes(i))
              numLocalBoundaryNodes++;
          RCP<const Teuchos::Comm<int>> comm = A->getRowMap()->getComm();
          MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
          GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
        }

        Set(currentLevel, "DofsPerNode", blkSize);
        Set(currentLevel, "Graph", graph);

      } else {
        // ap: We make quite a few assumptions here; general case may be a lot different,
        // but much much harder to implement. We assume that:
        //  1) all maps are standard maps, not strided maps
        //  2) global indices of dofs in A are related to dofs in coordinates in a simple arithmetic
        //     way: rows i*blkSize, i*blkSize+1, ..., i*blkSize + (blkSize-1) correspond to node i
        //
        // NOTE: Potentially, some of the code below could be simplified with UnAmalgamationInfo,
        // but as I totally don't understand that code, here is my solution

        // [*1*]: see [*0*]

        // Check that the number of local coordinates is consistent with the #rows in A
        TEUCHOS_TEST_FOR_EXCEPTION(A->getRowMap()->getLocalNumElements() / blkSize != Coords->getLocalLength(), Exceptions::Incompatible,
                                   "Coordinate vector length (" << Coords->getLocalLength() << ") is incompatible with number of rows in A (" << A->getRowMap()->getLocalNumElements() << ") by modulo block size (" << blkSize << ").");

        const RCP<const Map> colMap = A->getColMap();
        RCP<const Map> uniqueMap, nonUniqueMap;
        Array<LO> colTranslation;
        if (blkSize == 1) {
          uniqueMap    = A->getRowMap();
          nonUniqueMap = A->getColMap();
          graphType    = "unamalgamated";

        } else {
          uniqueMap = Coords->getMap();
          TEUCHOS_TEST_FOR_EXCEPTION(uniqueMap->getIndexBase() != indexBase, Exceptions::Incompatible,
                                     "Different index bases for matrix and coordinates");

          AmalgamationFactory::AmalgamateMap(*(A->getColMap()), *A, nonUniqueMap, colTranslation);

          graphType = "amalgamated";
        }
        LO numRows = Teuchos::as<LocalOrdinal>(uniqueMap->getLocalNumElements());

        RCP<RealValuedMultiVector> ghostedCoords;
        RCP<Vector> ghostedLaplDiag;
        Teuchos::ArrayRCP<SC> ghostedLaplDiagData;
        if (threshold != STS::zero()) {
          // Get ghost coordinates
          RCP<const Import> importer;
          {
            SubFactoryMonitor m1(*this, "Import construction", currentLevel);
            if (blkSize == 1 && realA->getCrsGraph()->getImporter() != Teuchos::null) {
              GetOStream(Warnings1) << "Using existing importer from matrix graph" << std::endl;
              importer = realA->getCrsGraph()->getImporter();
            } else {
              GetOStream(Warnings0) << "Constructing new importer instance" << std::endl;
              importer = ImportFactory::Build(uniqueMap, nonUniqueMap);
            }
          }  // subtimer
          ghostedCoords = Xpetra::MultiVectorFactory<real_type, LO, GO, NO>::Build(nonUniqueMap, Coords->getNumVectors());
          {
            SubFactoryMonitor m1(*this, "Coordinate import", currentLevel);
            ghostedCoords->doImport(*Coords, *importer, Xpetra::INSERT);
          }  // subtimer

          // Construct Distance Laplacian diagonal
          RCP<Vector> localLaplDiag = VectorFactory::Build(uniqueMap);
          Array<LO> indicesExtra;
          Teuchos::Array<Teuchos::ArrayRCP<const real_type>> coordData;
          if (threshold != STS::zero()) {
            const size_t numVectors = ghostedCoords->getNumVectors();
            coordData.reserve(numVectors);
            for (size_t j = 0; j < numVectors; j++) {
              Teuchos::ArrayRCP<const real_type> tmpData = ghostedCoords->getData(j);
              coordData.push_back(tmpData);
            }
          }
          {
            SubFactoryMonitor m1(*this, "Laplacian local diagonal", currentLevel);
            ArrayRCP<SC> localLaplDiagData = localLaplDiag->getDataNonConst(0);
            for (LO row = 0; row < numRows; row++) {
              ArrayView<const LO> indices;

              if (blkSize == 1) {
                ArrayView<const SC> vals;
                A->getLocalRowView(row, indices, vals);

              } else {
                // Merge rows of A
                indicesExtra.resize(0);
                MergeRows(*A, row, indicesExtra, colTranslation);
                indices = indicesExtra;
              }

              LO nnz               = indices.size();
              bool haveAddedToDiag = false;
              for (LO colID = 0; colID < nnz; colID++) {
                const LO col = indices[colID];

                if (row != col) {
                  if (use_dlap_weights == SINGLE_WEIGHTS) {
                    /*printf("[%d,%d] Unweighted Distance = %6.4e Weighted Distance = %6.4e\n",row,col,
                           MueLu::Utilities<real_type,LO,GO,NO>::Distance2(coordData, row, col),
                           MueLu::Utilities<real_type,LO,GO,NO>::Distance2(dlap_weights(),coordData, row, col));*/
                    localLaplDiagData[row] += STS::one() / MueLu::Utilities<real_type, LO, GO, NO>::Distance2(dlap_weights(), coordData, row, col);
                  } else if (use_dlap_weights == BLOCK_WEIGHTS) {
                    int block_id    = row % interleaved_blocksize;
                    int block_start = block_id * interleaved_blocksize;
                    localLaplDiagData[row] += STS::one() / MueLu::Utilities<real_type, LO, GO, NO>::Distance2(dlap_weights(block_start, interleaved_blocksize), coordData, row, col);
                  } else {
                    //                    printf("[%d,%d] Unweighted Distance = %6.4e\n",row,col,MueLu::Utilities<real_type,LO,GO,NO>::Distance2(coordData, row, col));
                    localLaplDiagData[row] += STS::one() / MueLu::Utilities<real_type, LO, GO, NO>::Distance2(coordData, row, col);
                  }
                  haveAddedToDiag = true;
                }
              }
              // Deal with the situation where boundary conditions have only been enforced on rows, but not on columns.
              // We enforce dropping of these entries by assigning a very large number to the diagonal entries corresponding to BCs.
              if (!haveAddedToDiag)
                localLaplDiagData[row] = STS::rmax();
            }
          }  // subtimer
          {
            SubFactoryMonitor m1(*this, "Laplacian distributed diagonal", currentLevel);
            ghostedLaplDiag = VectorFactory::Build(nonUniqueMap);
            ghostedLaplDiag->doImport(*localLaplDiag, *importer, Xpetra::INSERT);
            ghostedLaplDiagData = ghostedLaplDiag->getDataNonConst(0);
          }  // subtimer

        } else {
          GetOStream(Runtime0) << "Skipping distance laplacian construction due to 0 threshold" << std::endl;
        }

        // NOTE: ghostedLaplDiagData might be zero if we don't actually calculate the laplacian

        // allocate space for the local graph
        typename LWGraph::row_type::non_const_type rows("rows", numRows + 1);
        typename LWGraph::entries_type::non_const_type columns("columns", A->getLocalNumEntries());

#ifdef HAVE_MUELU_DEBUG
        // DEBUGGING
        for (LO i = 0; i < (LO)columns.size(); i++) columns[i] = -666;
#endif

        // Extra array for if we're allowing symmetrization with cutting
        ArrayRCP<LO> rows_stop;
        bool use_stop_array = threshold != STS::zero() && distanceLaplacianAlgo == scaled_cut_symmetric;
        if (use_stop_array)
          // rows_stop = typename LWGraph::row_type::non_const_type("rows_stop", numRows);
          rows_stop.resize(numRows);

        typename LWGraph::boundary_nodes_type amalgBoundaryNodes("amalgBoundaryNodes", numRows);
        Kokkos::deep_copy(amalgBoundaryNodes, false);

        LO realnnz = 0;
        rows(0)    = 0;

        Array<LO> indicesExtra;
        {
          SubFactoryMonitor m1(*this, "Laplacian dropping", currentLevel);
          Teuchos::Array<Teuchos::ArrayRCP<const real_type>> coordData;
          if (threshold != STS::zero()) {
            const size_t numVectors = ghostedCoords->getNumVectors();
            coordData.reserve(numVectors);
            for (size_t j = 0; j < numVectors; j++) {
              Teuchos::ArrayRCP<const real_type> tmpData = ghostedCoords->getData(j);
              coordData.push_back(tmpData);
            }
          }

          ArrayView<const SC> vals;  // CMS hackery
          for (LO row = 0; row < numRows; row++) {
            ArrayView<const LO> indices;
            indicesExtra.resize(0);
            bool isBoundary = false;

            if (blkSize == 1) {
              //	      ArrayView<const SC>     vals;//CMS uncomment
              A->getLocalRowView(row, indices, vals);
              isBoundary = pointBoundaryNodes[row];
            } else {
              // The amalgamated row is marked as Dirichlet iff all point rows are Dirichlet
              for (LO j = 0; j < blkSize; j++) {
                if (!pointBoundaryNodes[row * blkSize + j]) {
                  isBoundary = false;
                  break;
                }
              }

              // Merge rows of A
              if (!isBoundary)
                MergeRows(*A, row, indicesExtra, colTranslation);
              else
                indicesExtra.push_back(row);
              indices = indicesExtra;
            }
            numTotal += indices.size();

            LO nnz = indices.size(), rownnz = 0;

            if (use_stop_array) {
              rows(row + 1) = rows(row) + nnz;
              realnnz       = rows(row);
            }

            if (threshold != STS::zero()) {
              // default
              if (distanceLaplacianAlgo == defaultAlgo) {
                /* Standard Distance Laplacian */
                for (LO colID = 0; colID < nnz; colID++) {
                  LO col = indices[colID];

                  if (row == col) {
                    columns(realnnz++) = col;
                    rownnz++;
                    continue;
                  }

                  // We do not want the distance Laplacian aggregating boundary nodes
                  if (isBoundary) continue;

                  SC laplVal;
                  if (use_dlap_weights == SINGLE_WEIGHTS) {
                    laplVal = STS::one() / MueLu::Utilities<real_type, LO, GO, NO>::Distance2(dlap_weights(), coordData, row, col);
                  } else if (use_dlap_weights == BLOCK_WEIGHTS) {
                    int block_id    = row % interleaved_blocksize;
                    int block_start = block_id * interleaved_blocksize;
                    laplVal         = STS::one() / MueLu::Utilities<real_type, LO, GO, NO>::Distance2(dlap_weights(block_start, interleaved_blocksize), coordData, row, col);
                  } else {
                    laplVal = STS::one() / MueLu::Utilities<real_type, LO, GO, NO>::Distance2(coordData, row, col);
                  }
                  real_type aiiajj = STS::magnitude(realThreshold * realThreshold * ghostedLaplDiagData[row] * ghostedLaplDiagData[col]);
                  real_type aij    = STS::magnitude(laplVal * laplVal);

                  if (aij > aiiajj) {
                    columns(realnnz++) = col;
                    rownnz++;
                  } else {
                    numDropped++;
                  }
                }
              } else {
                /* Cut Algorithm */
                using DropTol = Details::DropTol<real_type, LO>;
                std::vector<DropTol> drop_vec;
                drop_vec.reserve(nnz);
                const real_type zero = Teuchos::ScalarTraits<real_type>::zero();
                const real_type one  = Teuchos::ScalarTraits<real_type>::one();

                // find magnitudes
                for (LO colID = 0; colID < nnz; colID++) {
                  LO col = indices[colID];

                  if (row == col) {
                    drop_vec.emplace_back(zero, one, colID, false);
                    continue;
                  }
                  // We do not want the distance Laplacian aggregating boundary nodes
                  if (isBoundary) continue;

                  SC laplVal;
                  if (use_dlap_weights == SINGLE_WEIGHTS) {
                    laplVal = STS::one() / MueLu::Utilities<real_type, LO, GO, NO>::Distance2(dlap_weights(), coordData, row, col);
                  } else if (use_dlap_weights == BLOCK_WEIGHTS) {
                    int block_id    = row % interleaved_blocksize;
                    int block_start = block_id * interleaved_blocksize;
                    laplVal         = STS::one() / MueLu::Utilities<real_type, LO, GO, NO>::Distance2(dlap_weights(block_start, interleaved_blocksize), coordData, row, col);
                  } else {
                    laplVal = STS::one() / MueLu::Utilities<real_type, LO, GO, NO>::Distance2(coordData, row, col);
                  }

                  real_type aiiajj = STS::magnitude(ghostedLaplDiagData[row] * ghostedLaplDiagData[col]);
                  real_type aij    = STS::magnitude(laplVal * laplVal);

                  drop_vec.emplace_back(aij, aiiajj, colID, false);
                }

                const size_t n = drop_vec.size();

                if (distanceLaplacianAlgo == unscaled_cut) {
                  std::sort(drop_vec.begin(), drop_vec.end(), [](DropTol const& a, DropTol const& b) {
                    return a.val > b.val;
                  });

                  bool drop = false;
                  for (size_t i = 1; i < n; ++i) {
                    if (!drop) {
                      auto const& x = drop_vec[i - 1];
                      auto const& y = drop_vec[i];
                      auto a        = x.val;
                      auto b        = y.val;
                      if (a > realThreshold * b) {
                        drop = true;
#ifdef HAVE_MUELU_DEBUG
                        if (distanceLaplacianCutVerbose) {
                          std::cout << "DJS: KEEP, N, ROW:  " << i + 1 << ", " << n << ", " << row << std::endl;
                        }
#endif
                      }
                    }
                    drop_vec[i].drop = drop;
                  }
                } else if (distanceLaplacianAlgo == scaled_cut || distanceLaplacianAlgo == scaled_cut_symmetric) {
                  std::sort(drop_vec.begin(), drop_vec.end(), [](DropTol const& a, DropTol const& b) {
                    return a.val / a.diag > b.val / b.diag;
                  });

                  bool drop = false;
                  for (size_t i = 1; i < n; ++i) {
                    if (!drop) {
                      auto const& x = drop_vec[i - 1];
                      auto const& y = drop_vec[i];
                      auto a        = x.val / x.diag;
                      auto b        = y.val / y.diag;
                      if (a > realThreshold * b) {
                        drop = true;
#ifdef HAVE_MUELU_DEBUG
                        if (distanceLaplacianCutVerbose) {
                          std::cout << "DJS: KEEP, N, ROW:  " << i + 1 << ", " << n << ", " << row << std::endl;
                        }
#endif
                      }
                    }
                    drop_vec[i].drop = drop;
                  }
                }

                std::sort(drop_vec.begin(), drop_vec.end(), [](DropTol const& a, DropTol const& b) {
                  return a.col < b.col;
                });

                for (LO idxID = 0; idxID < (LO)drop_vec.size(); idxID++) {
                  LO col = indices[drop_vec[idxID].col];

                  // don't drop diagonal
                  if (row == col) {
                    columns(realnnz++) = col;
                    rownnz++;
                    //		    printf("(%d,%d) KEEP %13s matrix = %6.4e\n",row,row,"DIAGONAL",drop_vec[idxID].aux_val);
                    continue;
                  }

                  if (!drop_vec[idxID].drop) {
                    columns(realnnz++) = col;
                    //		    printf("(%d,%d) KEEP dlap = %6.4e matrix = %6.4e\n",row,col,drop_vec[idxID].val/drop_vec[idxID].diag,drop_vec[idxID].aux_val);
                    rownnz++;
                  } else {
                    //		    printf("(%d,%d) DROP dlap = %6.4e matrix = %6.4e\n",row,col,drop_vec[idxID].val/drop_vec[idxID].diag,drop_vec[idxID].aux_val);
                    numDropped++;
                  }
                }
              }
            } else {
              // Skip laplace calculation and threshold comparison for zero threshold
              for (LO colID = 0; colID < nnz; colID++) {
                LO col             = indices[colID];
                columns(realnnz++) = col;
                rownnz++;
              }
            }

            if (rownnz == 1) {
              // If the only element remaining after filtering is diagonal, mark node as boundary
              // FIXME: this should really be replaced by the following
              //    if (indices.size() == 1 && indices[0] == row)
              //        boundaryNodes[row] = true;
              // We do not do it this way now because there is no framework for distinguishing isolated
              // and boundary nodes in the aggregation algorithms
              amalgBoundaryNodes[row] = true;
            }

            if (use_stop_array)
              rows_stop[row] = rownnz + rows[row];
            else
              rows[row + 1] = realnnz;
          }  // for (LO row = 0; row < numRows; row++)

        }  // subtimer

        if (use_stop_array) {
          // Do symmetrization of the cut matrix
          // NOTE: We assume nested row/column maps here
          for (LO row = 0; row < numRows; row++) {
            for (LO colidx = rows[row]; colidx < rows_stop[row]; colidx++) {
              LO col = columns[colidx];
              if (col >= numRows) continue;

              bool found = false;
              for (LO t_col = rows(col); !found && t_col < rows_stop[col]; t_col++) {
                if (columns[t_col] == row)
                  found = true;
              }
              // We didn't find the transpose buddy, so let's symmetrize, unless we'd be symmetrizing
              // into a Dirichlet unknown.  In that case don't.
              if (!found && !pointBoundaryNodes[col] && Teuchos::as<typename LWGraph::row_type::value_type>(rows_stop[col]) < rows[col + 1]) {
                LO new_idx = rows_stop[col];
                //		  printf("(%d,%d) SYMADD entry\n",col,row);
                columns[new_idx] = row;
                rows_stop[col]++;
                numDropped--;
              }
            }
          }

          // Condense everything down
          LO current_start = 0;
          for (LO row = 0; row < numRows; row++) {
            LO old_start = current_start;
            for (LO col = rows(row); col < rows_stop[row]; col++) {
              if (current_start != col) {
                columns(current_start) = columns(col);
              }
              current_start++;
            }
            rows[row] = old_start;
          }
          rows(numRows) = realnnz = current_start;
        }

        RCP<LWGraph> graph;
        {
          SubFactoryMonitor m1(*this, "Build amalgamated graph", currentLevel);
          graph = rcp(new LWGraph(rows, Kokkos::subview(columns, Kokkos::make_pair(0, realnnz)), uniqueMap, nonUniqueMap, "amalgamated graph of A"));
          graph->SetBoundaryNodeMap(amalgBoundaryNodes);
        }  // subtimer

        if (GetVerbLevel() & Statistics1) {
          GO numLocalBoundaryNodes  = 0;
          GO numGlobalBoundaryNodes = 0;

          for (size_t i = 0; i < amalgBoundaryNodes.size(); ++i)
            if (amalgBoundaryNodes(i))
              numLocalBoundaryNodes++;

          RCP<const Teuchos::Comm<int>> comm = A->getRowMap()->getComm();
          MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
          GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " agglomerated Dirichlet nodes"
                                  << " using threshold " << dirichletThreshold << std::endl;
        }

        Set(currentLevel, "Graph", graph);
        Set(currentLevel, "DofsPerNode", blkSize);
      }
    }

    if ((GetVerbLevel() & Statistics1) && !(A->GetFixedBlockSize() > 1 && threshold != STS::zero())) {
      RCP<const Teuchos::Comm<int>> comm = A->getRowMap()->getComm();
      GO numGlobalTotal, numGlobalDropped;
      MueLu_sumAll(comm, numTotal, numGlobalTotal);
      MueLu_sumAll(comm, numDropped, numGlobalDropped);
      GetOStream(Statistics1) << "Number of dropped entries in " << graphType << " matrix graph: " << numGlobalDropped << "/" << numGlobalTotal;
      if (numGlobalTotal != 0)
        GetOStream(Statistics1) << " (" << 100 * Teuchos::as<double>(numGlobalDropped) / Teuchos::as<double>(numGlobalTotal) << "%)";
      GetOStream(Statistics1) << std::endl;
    }

  } else {
    // what Tobias has implemented

    SC threshold = as<SC>(pL.get<double>("aggregation: drop tol"));
    // GetOStream(Runtime0) << "algorithm = \"" << algo << "\": threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << std::endl;
    GetOStream(Runtime0) << "algorithm = \""
                         << "failsafe"
                         << "\": threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << std::endl;
    Set<bool>(currentLevel, "Filtering", (threshold != STS::zero()));

    RCP<const Map> rowMap = A->getRowMap();
    RCP<const Map> colMap = A->getColMap();

    LO blockdim  = 1;                       // block dim for fixed size blocks
    GO indexBase = rowMap->getIndexBase();  // index base of maps
    GO offset    = 0;

    // 1) check for blocking/striding information
    if (A->IsView("stridedMaps") &&
        Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap("stridedMaps")) != Teuchos::null) {
      Xpetra::viewLabel_t oldView  = A->SwitchToView("stridedMaps");  // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
      RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap());
      TEUCHOS_TEST_FOR_EXCEPTION(strMap == Teuchos::null, Exceptions::BadCast, "MueLu::CoalesceFactory::Build: cast to strided row map failed.");
      blockdim = strMap->getFixedBlockSize();
      offset   = strMap->getOffset();
      oldView  = A->SwitchToView(oldView);
      GetOStream(Statistics1) << "CoalesceDropFactory::Build():"
                              << " found blockdim=" << blockdim << " from strided maps. offset=" << offset << std::endl;
    } else
      GetOStream(Statistics1) << "CoalesceDropFactory::Build(): no striding information available. Use blockdim=1 with offset=0" << std::endl;

    // 2) get row map for amalgamated matrix (graph of A)
    //    with same distribution over all procs as row map of A
    RCP<const Map> nodeMap = amalInfo->getNodeRowMap();
    GetOStream(Statistics1) << "CoalesceDropFactory: nodeMap " << nodeMap->getLocalNumElements() << "/" << nodeMap->getGlobalNumElements() << " elements" << std::endl;

    // 3) create graph of amalgamated matrix
    RCP<CrsGraph> crsGraph = CrsGraphFactory::Build(nodeMap, A->getLocalMaxNumRowEntries() * blockdim);

    LO numRows  = A->getRowMap()->getLocalNumElements();
    LO numNodes = nodeMap->getLocalNumElements();
    typename LWGraph::boundary_nodes_type amalgBoundaryNodes("amalgBoundaryNodes", numNodes);
    Kokkos::deep_copy(amalgBoundaryNodes, false);
    const ArrayRCP<int> numberDirichletRowsPerNode(numNodes, 0);  // helper array counting the number of Dirichlet nodes associated with node
    bool bIsDiagonalEntry = false;                                // boolean flag stating that grid==gcid

    // 4) do amalgamation. generate graph of amalgamated matrix
    //    Note, this code is much more inefficient than the leightwight implementation
    //    Most of the work has already been done in the AmalgamationFactory
    for (LO row = 0; row < numRows; row++) {
      // get global DOF id
      GO grid = rowMap->getGlobalElement(row);

      // reinitialize boolean helper variable
      bIsDiagonalEntry = false;

      // translate grid to nodeid
      GO nodeId = AmalgamationFactory::DOFGid2NodeId(grid, blockdim, offset, indexBase);

      size_t nnz = A->getNumEntriesInLocalRow(row);
      Teuchos::ArrayView<const LO> indices;
      Teuchos::ArrayView<const SC> vals;
      A->getLocalRowView(row, indices, vals);

      RCP<std::vector<GO>> cnodeIds = Teuchos::rcp(new std::vector<GO>);  // global column block ids
      LO realnnz                    = 0;
      for (LO col = 0; col < Teuchos::as<LO>(nnz); col++) {
        GO gcid = colMap->getGlobalElement(indices[col]);  // global column id

        if (vals[col] != STS::zero()) {
          GO cnodeId = AmalgamationFactory::DOFGid2NodeId(gcid, blockdim, offset, indexBase);
          cnodeIds->push_back(cnodeId);
          realnnz++;  // increment number of nnz in matrix row
          if (grid == gcid) bIsDiagonalEntry = true;
        }
      }

      if (realnnz == 1 && bIsDiagonalEntry == true) {
        LO lNodeId = nodeMap->getLocalElement(nodeId);
        numberDirichletRowsPerNode[lNodeId] += 1;             // increment Dirichlet row counter associated with lNodeId
        if (numberDirichletRowsPerNode[lNodeId] == blockdim)  // mark full Dirichlet nodes
          amalgBoundaryNodes[lNodeId] = true;
      }

      Teuchos::ArrayRCP<GO> arr_cnodeIds = Teuchos::arcp(cnodeIds);

      if (arr_cnodeIds.size() > 0)
        crsGraph->insertGlobalIndices(nodeId, arr_cnodeIds());
    }
    // fill matrix graph
    crsGraph->fillComplete(nodeMap, nodeMap);

    // 5) create MueLu Graph object
    RCP<LWGraph> graph = rcp(new LWGraph(crsGraph, "amalgamated graph of A"));

    // Detect and record rows that correspond to Dirichlet boundary conditions
    graph->SetBoundaryNodeMap(amalgBoundaryNodes);

    if (GetVerbLevel() & Statistics1) {
      GO numLocalBoundaryNodes  = 0;
      GO numGlobalBoundaryNodes = 0;
      for (size_t i = 0; i < amalgBoundaryNodes.size(); ++i)
        if (amalgBoundaryNodes(i))
          numLocalBoundaryNodes++;
      RCP<const Teuchos::Comm<int>> comm = A->getRowMap()->getComm();
      MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
      GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
    }

    // 6) store results in Level
    // graph->SetBoundaryNodeMap(gBoundaryNodeMap);
    Set(currentLevel, "DofsPerNode", blockdim);
    Set(currentLevel, "Graph", graph);

  }  // if (doExperimentalWrap) ... else ...

}  // Build

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MergeRows(const Matrix& A, const LO row, Array<LO>& cols, const Array<LO>& translation) const {
  typedef typename ArrayView<const LO>::size_type size_type;

  // extract striding information
  LO blkSize = A.GetFixedBlockSize();  //< stores the size of the block within the strided map
  if (A.IsView("stridedMaps") == true) {
    Teuchos::RCP<const Map> myMap         = A.getRowMap("stridedMaps");
    Teuchos::RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(myMap);
    TEUCHOS_TEST_FOR_EXCEPTION(strMap == null, Exceptions::RuntimeError, "Map is not of type StridedMap");
    if (strMap->getStridedBlockId() > -1)
      blkSize = Teuchos::as<LO>(strMap->getStridingData()[strMap->getStridedBlockId()]);
  }

  // count nonzero entries in all dof rows associated with node row
  size_t nnz = 0, pos = 0;
  for (LO j = 0; j < blkSize; j++)
    nnz += A.getNumEntriesInLocalRow(row * blkSize + j);

  if (nnz == 0) {
    cols.resize(0);
    return;
  }

  cols.resize(nnz);

  // loop over all local dof rows associated with local node "row"
  ArrayView<const LO> inds;
  ArrayView<const SC> vals;
  for (LO j = 0; j < blkSize; j++) {
    A.getLocalRowView(row * blkSize + j, inds, vals);
    size_type numIndices = inds.size();

    if (numIndices == 0)  // skip empty dof rows
      continue;

    // cols: stores all local node ids for current local node id "row"
    cols[pos++] = translation[inds[0]];
    for (size_type k = 1; k < numIndices; k++) {
      LO nodeID = translation[inds[k]];
      // Here we try to speed up the process by reducing the size of an array
      // to sort. This works if the column nonzeros belonging to the same
      // node are stored consequently.
      if (nodeID != cols[pos - 1])
        cols[pos++] = nodeID;
    }
  }
  cols.resize(pos);
  nnz = pos;

  // Sort and remove duplicates
  std::sort(cols.begin(), cols.end());
  pos = 0;
  for (size_t j = 1; j < nnz; j++)
    if (cols[j] != cols[pos])
      cols[++pos] = cols[j];
  cols.resize(pos + 1);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MergeRowsWithDropping(const Matrix& A, const LO row, const ArrayRCP<const SC>& ghostedDiagVals, SC threshold, Array<LO>& cols, const Array<LO>& translation) const {
  typedef typename ArrayView<const LO>::size_type size_type;
  typedef Teuchos::ScalarTraits<SC> STS;

  // extract striding information
  LO blkSize = A.GetFixedBlockSize();  //< stores the size of the block within the strided map
  if (A.IsView("stridedMaps") == true) {
    Teuchos::RCP<const Map> myMap         = A.getRowMap("stridedMaps");
    Teuchos::RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(myMap);
    TEUCHOS_TEST_FOR_EXCEPTION(strMap == null, Exceptions::RuntimeError, "Map is not of type StridedMap");
    if (strMap->getStridedBlockId() > -1)
      blkSize = Teuchos::as<LO>(strMap->getStridingData()[strMap->getStridedBlockId()]);
  }

  // count nonzero entries in all dof rows associated with node row
  size_t nnz = 0, pos = 0;
  for (LO j = 0; j < blkSize; j++)
    nnz += A.getNumEntriesInLocalRow(row * blkSize + j);

  if (nnz == 0) {
    cols.resize(0);
    return;
  }

  cols.resize(nnz);

  // loop over all local dof rows associated with local node "row"
  ArrayView<const LO> inds;
  ArrayView<const SC> vals;
  for (LO j = 0; j < blkSize; j++) {
    A.getLocalRowView(row * blkSize + j, inds, vals);
    size_type numIndices = inds.size();

    if (numIndices == 0)  // skip empty dof rows
      continue;

    // cols: stores all local node ids for current local node id "row"
    LO prevNodeID = -1;
    for (size_type k = 0; k < numIndices; k++) {
      LO dofID  = inds[k];
      LO nodeID = translation[inds[k]];

      // we avoid a square root by using squared values
      typename STS::magnitudeType aiiajj = STS::magnitude(threshold * threshold * ghostedDiagVals[dofID] * ghostedDiagVals[row * blkSize + j]);  // eps^2 * |a_ii| * |a_jj|
      typename STS::magnitudeType aij    = STS::magnitude(vals[k] * vals[k]);

      // check dropping criterion
      if (aij > aiiajj || (row * blkSize + j == dofID)) {
        // accept entry in graph

        // Here we try to speed up the process by reducing the size of an array
        // to sort. This works if the column nonzeros belonging to the same
        // node are stored consequently.
        if (nodeID != prevNodeID) {
          cols[pos++] = nodeID;
          prevNodeID  = nodeID;
        }
      }
    }
  }
  cols.resize(pos);
  nnz = pos;

  // Sort and remove duplicates
  std::sort(cols.begin(), cols.end());
  pos = 0;
  for (size_t j = 1; j < nnz; j++)
    if (cols[j] != cols[pos])
      cols[++pos] = cols[j];
  cols.resize(pos + 1);

  return;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockDiagonalize(Level& currentLevel, const RCP<Matrix>& A, bool generate_matrix) const {
  typedef Teuchos::ScalarTraits<SC> STS;

  const ParameterList& pL                              = GetParameterList();
  const typename STS::magnitudeType dirichletThreshold = STS::magnitude(as<SC>(pL.get<double>("aggregation: Dirichlet threshold")));
  const typename STS::magnitudeType rowSumTol          = as<typename STS::magnitudeType>(pL.get<double>("aggregation: row sum drop tol"));

  RCP<LocalOrdinalVector> BlockNumber = Get<RCP<LocalOrdinalVector>>(currentLevel, "BlockNumber");
  RCP<LocalOrdinalVector> ghostedBlockNumber;
  GetOStream(Statistics1) << "Using BlockDiagonal Graph before dropping (with provided blocking)" << std::endl;

  // Ghost the column block numbers if we need to
  RCP<const Import> importer = A->getCrsGraph()->getImporter();
  if (!importer.is_null()) {
    SubFactoryMonitor m1(*this, "Block Number import", currentLevel);
    ghostedBlockNumber = Xpetra::VectorFactory<LO, LO, GO, NO>::Build(importer->getTargetMap());
    ghostedBlockNumber->doImport(*BlockNumber, *importer, Xpetra::INSERT);
  } else {
    ghostedBlockNumber = BlockNumber;
  }

  // Accessors for block numbers
  Teuchos::ArrayRCP<const LO> row_block_number = BlockNumber->getData(0);
  Teuchos::ArrayRCP<const LO> col_block_number = ghostedBlockNumber->getData(0);

  // allocate space for the local graph
  typename CrsMatrix::local_matrix_type::row_map_type::HostMirror::non_const_type rows_mat;
  typename LWGraph::row_type::non_const_type rows_graph;
  typename LWGraph::entries_type::non_const_type columns;
  typename CrsMatrix::local_matrix_type::values_type::HostMirror::non_const_type values;
  RCP<CrsMatrixWrap> crs_matrix_wrap;

  if (generate_matrix) {
    crs_matrix_wrap = rcp(new CrsMatrixWrap(A->getRowMap(), A->getColMap(), 0));
    rows_mat        = typename CrsMatrix::local_matrix_type::row_map_type::HostMirror::non_const_type("rows_mat", A->getLocalNumRows() + 1);
  } else {
    rows_graph = typename LWGraph::row_type::non_const_type("rows_graph", A->getLocalNumRows() + 1);
  }
  columns = typename LWGraph::entries_type::non_const_type("columns", A->getLocalNumEntries());
  values  = typename CrsMatrix::local_matrix_type::values_type::HostMirror::non_const_type("values", A->getLocalNumEntries());

  LO realnnz    = 0;
  GO numDropped = 0, numTotal = 0;
  for (LO row = 0; row < Teuchos::as<LO>(A->getRowMap()->getLocalNumElements()); ++row) {
    LO row_block = row_block_number[row];
    size_t nnz   = A->getNumEntriesInLocalRow(row);
    ArrayView<const LO> indices;
    ArrayView<const SC> vals;
    A->getLocalRowView(row, indices, vals);

    LO rownnz = 0;
    for (LO colID = 0; colID < Teuchos::as<LO>(nnz); colID++) {
      LO col       = indices[colID];
      LO col_block = col_block_number[col];

      if (row_block == col_block) {
        if (generate_matrix) values[realnnz] = vals[colID];
        columns[realnnz++] = col;
        rownnz++;
      } else
        numDropped++;
    }
    if (generate_matrix)
      rows_mat[row + 1] = realnnz;
    else
      rows_graph[row + 1] = realnnz;
  }

  auto boundaryNodes = MueLu::Utilities<SC, LO, GO, NO>::DetectDirichletRows_kokkos_host(*A, dirichletThreshold);
  if (rowSumTol > 0.)
    Utilities::ApplyRowSumCriterionHost(*A, rowSumTol, boundaryNodes);

  numTotal = A->getLocalNumEntries();

  if (GetVerbLevel() & Statistics1) {
    GO numLocalBoundaryNodes  = 0;
    GO numGlobalBoundaryNodes = 0;
    for (size_t i = 0; i < boundaryNodes.size(); ++i)
      if (boundaryNodes(i))
        numLocalBoundaryNodes++;
    RCP<const Teuchos::Comm<int>> comm = A->getRowMap()->getComm();
    MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
    GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;

    GO numGlobalTotal, numGlobalDropped;
    MueLu_sumAll(comm, numTotal, numGlobalTotal);
    MueLu_sumAll(comm, numDropped, numGlobalDropped);
    GetOStream(Statistics1) << "Number of dropped entries in block-diagonalized matrix graph: " << numGlobalDropped << "/" << numGlobalTotal;
    if (numGlobalTotal != 0)
      GetOStream(Statistics1) << " (" << 100 * Teuchos::as<double>(numGlobalDropped) / Teuchos::as<double>(numGlobalTotal) << "%)";
    GetOStream(Statistics1) << std::endl;
  }

  Set(currentLevel, "Filtering", true);

  if (generate_matrix) {
    // NOTE: Trying to use A's Import/Export objects will cause the code to segfault back in Build() with errors on the Import
    // if you're using Epetra.  I'm not really sure why. By using the Col==Domain and Row==Range maps, we get null Import/Export objects
    // here, which is legit, because we never use them anyway.
    if constexpr (std::is_same<typename LWGraph::row_type,
                               typename CrsMatrix::local_matrix_type::row_map_type>::value) {
      crs_matrix_wrap->getCrsMatrix()->setAllValues(rows_mat, columns, values);
    } else {
      auto rows_mat2 = typename CrsMatrix::local_matrix_type::row_map_type::non_const_type("rows_mat2", rows_mat.extent(0));
      Kokkos::deep_copy(rows_mat2, rows_mat);
      auto columns2 = typename CrsMatrix::local_graph_type::entries_type::non_const_type("columns2", columns.extent(0));
      Kokkos::deep_copy(columns2, columns);
      auto values2 = typename CrsMatrix::local_matrix_type::values_type::non_const_type("values2", values.extent(0));
      Kokkos::deep_copy(values2, values);
      crs_matrix_wrap->getCrsMatrix()->setAllValues(rows_mat2, columns2, values2);
    }
    crs_matrix_wrap->getCrsMatrix()->expertStaticFillComplete(A->getColMap(), A->getRowMap());
  } else {
    RCP<LWGraph> graph = rcp(new LWGraph(rows_graph, Kokkos::subview(columns, Kokkos::make_pair(0, realnnz)), A->getRowMap(), A->getColMap(), "block-diagonalized graph of A"));
    graph->SetBoundaryNodeMap(boundaryNodes);
    Set(currentLevel, "Graph", graph);
  }

  Set(currentLevel, "DofsPerNode", 1);
  return crs_matrix_wrap;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockDiagonalizeGraph(const RCP<LWGraph>& inputGraph, const RCP<LocalOrdinalVector>& ghostedBlockNumber, RCP<LWGraph>& outputGraph, RCP<const Import>& importer) const {
  TEUCHOS_TEST_FOR_EXCEPTION(ghostedBlockNumber.is_null(), Exceptions::RuntimeError, "BlockDiagonalizeGraph(): ghostedBlockNumber is null.");
  const ParameterList& pL = GetParameterList();

  const bool localizeColoringGraph = pL.get<bool>("aggregation: coloring: localize color graph");

  GetOStream(Statistics1) << "Using BlockDiagonal Graph after Dropping (with provided blocking)";
  if (localizeColoringGraph)
    GetOStream(Statistics1) << ", with localization" << std::endl;
  else
    GetOStream(Statistics1) << ", without localization" << std::endl;

  // Accessors for block numbers
  Teuchos::ArrayRCP<const LO> row_block_number = ghostedBlockNumber->getData(0);
  Teuchos::ArrayRCP<const LO> col_block_number = ghostedBlockNumber->getData(0);

  // allocate space for the local graph
  ArrayRCP<size_t> rows_mat;
  typename LWGraph::row_type::non_const_type rows_graph("rows_graph", inputGraph->GetNodeNumVertices() + 1);
  typename LWGraph::entries_type::non_const_type columns("columns", inputGraph->GetNodeNumEdges());

  LO realnnz    = 0;
  GO numDropped = 0, numTotal = 0;
  const LO numRows = Teuchos::as<LO>(inputGraph->GetDomainMap()->getLocalNumElements());
  if (localizeColoringGraph) {
    for (LO row = 0; row < numRows; ++row) {
      LO row_block = row_block_number[row];
      auto indices = inputGraph->getNeighborVertices(row);

      LO rownnz = 0;
      for (LO colID = 0; colID < Teuchos::as<LO>(indices.length); colID++) {
        LO col       = indices(colID);
        LO col_block = col_block_number[col];

        if ((row_block == col_block) && (col < numRows)) {
          columns(realnnz++) = col;
          rownnz++;
        } else
          numDropped++;
      }
      rows_graph(row + 1) = realnnz;
    }
  } else {
    // ghosting of boundary node map
    auto boundaryNodes       = inputGraph->GetBoundaryNodeMap();
    auto boundaryNodesVector = Xpetra::VectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>::Build(inputGraph->GetDomainMap());
    for (size_t i = 0; i < inputGraph->GetNodeNumVertices(); i++)
      boundaryNodesVector->getDataNonConst(0)[i] = boundaryNodes[i];
    // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("boundary",*boundaryNodesVector);
    auto boundaryColumnVector = Xpetra::VectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>::Build(inputGraph->GetImportMap());
    boundaryColumnVector->doImport(*boundaryNodesVector, *importer, Xpetra::INSERT);
    auto boundaryColumn = boundaryColumnVector->getData(0);

    for (LO row = 0; row < numRows; ++row) {
      LO row_block = row_block_number[row];
      auto indices = inputGraph->getNeighborVertices(row);

      LO rownnz = 0;
      for (LO colID = 0; colID < Teuchos::as<LO>(indices.length); colID++) {
        LO col       = indices(colID);
        LO col_block = col_block_number[col];

        if ((row_block == col_block) && ((row == col) || (boundaryColumn[col] == 0))) {
          columns(realnnz++) = col;
          rownnz++;
        } else
          numDropped++;
      }
      rows_graph(row + 1) = realnnz;
    }
  }

  numTotal = inputGraph->GetNodeNumEdges();

  if (GetVerbLevel() & Statistics1) {
    RCP<const Teuchos::Comm<int>> comm = inputGraph->GetDomainMap()->getComm();
    GO numGlobalTotal, numGlobalDropped;
    MueLu_sumAll(comm, numTotal, numGlobalTotal);
    MueLu_sumAll(comm, numDropped, numGlobalDropped);
    GetOStream(Statistics1) << "Number of dropped entries in block-diagonalized matrix graph: " << numGlobalDropped << "/" << numGlobalTotal;
    if (numGlobalTotal != 0)
      GetOStream(Statistics1) << " (" << 100 * Teuchos::as<double>(numGlobalDropped) / Teuchos::as<double>(numGlobalTotal) << "%)";
    GetOStream(Statistics1) << std::endl;
  }

  if (localizeColoringGraph) {
    outputGraph = rcp(new LWGraph(rows_graph, Kokkos::subview(columns, Kokkos::make_pair(0, realnnz)), inputGraph->GetDomainMap(), inputGraph->GetImportMap(), "block-diagonalized graph of A"));
    outputGraph->SetBoundaryNodeMap(inputGraph->GetBoundaryNodeMap());
  } else {
    TEUCHOS_ASSERT(inputGraph->GetDomainMap()->lib() == Xpetra::UseTpetra);
#ifdef HAVE_XPETRA_TPETRA
    auto outputGraph2 = rcp(new LWGraph(rows_graph, Kokkos::subview(columns, Kokkos::make_pair(0, realnnz)), inputGraph->GetDomainMap(), inputGraph->GetImportMap(), "block-diagonalized graph of A"));

    auto tpGraph     = Xpetra::toTpetra(rcp_const_cast<const CrsGraph>(outputGraph2->GetCrsGraph()));
    auto sym         = rcp(new Tpetra::CrsGraphTransposer<LocalOrdinal, GlobalOrdinal, Node>(tpGraph));
    auto tpGraphSym  = sym->symmetrize();
    auto lclGraphSym = tpGraphSym->getLocalGraphHost();
    auto colIndsSym  = lclGraphSym.entries;

    auto rowsSym = tpGraphSym->getLocalRowPtrsHost();
    typename LWGraph::row_type::non_const_type rows_graphSym("rows_graphSym", rowsSym.size());
    for (size_t row = 0; row < rowsSym.size(); row++)
      rows_graphSym(row) = rowsSym(row);
    outputGraph = rcp(new LWGraph(rows_graphSym, colIndsSym, inputGraph->GetDomainMap(), Xpetra::toXpetra(tpGraphSym->getColMap()), "block-diagonalized graph of A"));
    outputGraph->SetBoundaryNodeMap(inputGraph->GetBoundaryNodeMap());
#endif
  }
}

}  // namespace MueLu

#endif  // MUELU_COALESCEDROPFACTORY_DEF_HPP
