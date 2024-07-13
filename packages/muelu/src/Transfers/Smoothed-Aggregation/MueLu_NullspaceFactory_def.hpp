// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_NULLSPACEFACTORY_DEF_HPP
#define MUELU_NULLSPACEFACTORY_DEF_HPP

#include "MueLu_NullspaceFactory_decl.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_MasterList.hpp"
#include "Xpetra_Access.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("nullspace: calculate rotations");
#undef SET_VALID_ENTRY
  validParamList->set<std::string>("Fine level nullspace", "Nullspace", "Variable name which is used to store null space multi vector on the finest level (default=\"Nullspace\"). For block matrices also \"Nullspace1\" to \"Nullspace9\" are accepted to describe the null space vectors for the (i,i) block (i=1..9).");

  validParamList->set<RCP<const FactoryBase>>("A", Teuchos::null, "Generating factory of the fine level matrix (only needed if default null space is generated)");
  validParamList->set<RCP<const FactoryBase>>("Nullspace", Teuchos::null, "Generating factory of the fine level null space");
  validParamList->set<RCP<const FactoryBase>>("Coordinates", Teuchos::null, "Generating factory of the coordinates");

  // TODO not very elegant.

  // 1/20/2016: we could add a sublist (e.g. "Nullspaces" which is excluded from parameter validation)
  validParamList->set<RCP<const FactoryBase>>("Nullspace1", Teuchos::null, "Generating factory of the fine level null space associated with the (1,1) block in your n x n block matrix.");
  validParamList->set<RCP<const FactoryBase>>("Nullspace2", Teuchos::null, "Generating factory of the fine level null space associated with the (2,2) block in your n x n block matrix.");
  validParamList->set<RCP<const FactoryBase>>("Nullspace3", Teuchos::null, "Generating factory of the fine level null space associated with the (3,3) block in your n x n block matrix.");
  validParamList->set<RCP<const FactoryBase>>("Nullspace4", Teuchos::null, "Generating factory of the fine level null space associated with the (4,4) block in your n x n block matrix.");
  validParamList->set<RCP<const FactoryBase>>("Nullspace5", Teuchos::null, "Generating factory of the fine level null space associated with the (5,5) block in your n x n block matrix.");
  validParamList->set<RCP<const FactoryBase>>("Nullspace6", Teuchos::null, "Generating factory of the fine level null space associated with the (6,6) block in your n x n block matrix.");
  validParamList->set<RCP<const FactoryBase>>("Nullspace7", Teuchos::null, "Generating factory of the fine level null space associated with the (7,7) block in your n x n block matrix.");
  validParamList->set<RCP<const FactoryBase>>("Nullspace8", Teuchos::null, "Generating factory of the fine level null space associated with the (8,8) block in your n x n block matrix.");
  validParamList->set<RCP<const FactoryBase>>("Nullspace9", Teuchos::null, "Generating factory of the fine level null space associated with the (9,9) block in your n x n block matrix.");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  const ParameterList& pL = GetParameterList();
  std::string nspName     = pL.get<std::string>("Fine level nullspace");

  // only request "A" in DeclareInput if
  // 1) there is not nspName (e.g. "Nullspace") is available in Level, AND
  // 2) it is the finest level (i.e. LevelID == 0)
  if (currentLevel.IsAvailable(nspName, NoFactory::get()) == false && currentLevel.GetLevelID() == 0)
    Input(currentLevel, "A");

  if (currentLevel.GetLevelID() == 0 &&
      currentLevel.IsAvailable("Coordinates", NoFactory::get()) &&  // we have coordinates (provided by user app)
      pL.get<bool>("nullspace: calculate rotations")) {             // and we want to calculate rotation modes
    calculateRotations_ = true;
    Input(currentLevel, "Coordinates");
  }

  if (currentLevel.GetLevelID() != 0) {
    // validate nullspaceFact_
    // 1) The factory for "Nullspace" (or nspName) must not be Teuchos::null, since the default factory
    //    for "Nullspace" is a NullspaceFactory
    // 2) The factory for "Nullspace" (or nspName) must be a TentativePFactory or any other TwoLevelFactoryBase derived object
    //    which generates the variable "Nullspace" as output
    TEUCHOS_TEST_FOR_EXCEPTION(GetFactory(nspName).is_null(), Exceptions::RuntimeError,
                               "MueLu::NullspaceFactory::DeclareInput(): You must declare an existing factory which "
                               "produces the variable \"Nullspace\" in the NullspaceFactory (e.g. a TentativePFactory).");
    currentLevel.DeclareInput("Nullspace", GetFactory(nspName).get(), this); /* ! "Nullspace" and nspName mismatch possible here */
  }
}

template <class NullspaceType, class CoordsType, class MeanCoordsType, class LO>
class NullspaceFunctor {
 private:
  NullspaceType nullspace;
  CoordsType coords;
  MeanCoordsType mean;
  LO numPDEs;
  LO nullspaceDim;
  typedef typename NullspaceType::value_type SC;
  typedef Kokkos::ArithTraits<SC> ATS;

 public:
  NullspaceFunctor(NullspaceType nullspace_, CoordsType coords_, MeanCoordsType mean_, LO numPDEs_, LO nullspaceDim_)
    : nullspace(nullspace_)
    , coords(coords_)
    , mean(mean_)
    , numPDEs(numPDEs_)
    , nullspaceDim(nullspaceDim_) {
    static_assert(static_cast<int>(NullspaceType::rank) == 2, "Nullspace needs to be a rank 2 view.");
    static_assert(static_cast<int>(CoordsType::rank) == 2, "Coords needs to be a rank 2 view.");
    static_assert(static_cast<int>(MeanCoordsType::rank) == 1, "Mean needs to be a rank 1 view.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const LO j) const {
    SC one = ATS::one();
    for (LO i = 0; i < numPDEs; i++)
      nullspace(j * numPDEs + i, i) = one;
    if ((nullspaceDim > numPDEs) && (numPDEs > 1)) {
      // xy rotation
      nullspace(j * numPDEs + 0, numPDEs) = -(coords(j, 1) - mean(1));
      nullspace(j * numPDEs + 1, numPDEs) = (coords(j, 0) - mean(0));
    }
    if ((nullspaceDim == numPDEs + 3) && (numPDEs > 2)) {
      // xz rotation
      nullspace(j * numPDEs + 1, numPDEs + 1) = -(coords(j, 2) - mean(2));
      nullspace(j * numPDEs + 2, numPDEs + 1) = (coords(j, 1) - mean(1));

      // yz rotation
      nullspace(j * numPDEs + 0, numPDEs + 2) = -(coords(j, 2) - mean(2));
      nullspace(j * numPDEs + 2, numPDEs + 2) = (coords(j, 0) - mean(0));
    }
  }
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& currentLevel) const {
  FactoryMonitor m(*this, "Nullspace factory", currentLevel);

  RCP<MultiVector> nullspace;

  // TEUCHOS_TEST_FOR_EXCEPTION(currentLevel.GetLevelID() != 0, Exceptions::RuntimeError, "MueLu::NullspaceFactory::Build(): NullspaceFactory can be used for finest level (LevelID == 0) only.");
  const ParameterList& pL = GetParameterList();
  std::string nspName     = pL.get<std::string>("Fine level nullspace");

  // get coordinates and compute mean of coordinates. (or centroid).

  RCP<RealValuedMultiVector> Coords;

  if (currentLevel.GetLevelID() == 0) {
    if (currentLevel.IsAvailable(nspName, NoFactory::get())) {
      // When a fine nullspace have already been defined by user using Set("Nullspace", ...) or
      // Set("Nullspace1", ...), we use it.
      nullspace = currentLevel.Get<RCP<MultiVector>>(nspName, NoFactory::get());
      GetOStream(Runtime1) << "Use user-given nullspace " << nspName << ":"
                           << " nullspace dimension=" << nullspace->getNumVectors()
                           << " nullspace length=" << nullspace->getGlobalLength() << std::endl;

    } else {
      // "Nullspace" (nspName) is not available
      auto A = Get<RCP<Matrix>>(currentLevel, "A");

      // Determine numPDEs
      LO numPDEs = 1;
      if (A->IsView("stridedMaps") == true) {
        Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps");  // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
        TEUCHOS_TEST_FOR_EXCEPTION(rcp_dynamic_cast<const StridedMap>(A->getRowMap()).is_null(), Exceptions::BadCast,
                                   "MueLu::CoalesceFactory::Build: cast to strided row map failed.");
        numPDEs = rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getFixedBlockSize();
        oldView = A->SwitchToView(oldView);
      }

      LO nullspaceDim = numPDEs;

      CoordsType coordsView;
      MeanCoordsType meanView;
      if (calculateRotations_) {
        Coords = Get<RCP<RealValuedMultiVector>>(currentLevel, "Coordinates");
        if (Coords->getNumVectors() > 1) nullspaceDim++;
        if (Coords->getNumVectors() > 2) nullspaceDim += 2;

        meanView = MeanCoordsType("mean coords", Coords->getNumVectors());
        Teuchos::Array<coordinate_type> hostMeans(Coords->getNumVectors());
        Coords->meanValue(hostMeans);
        Kokkos::View<typename RealValuedMultiVector::impl_scalar_type*, Kokkos::HostSpace> hostMeanView(hostMeans.getRawPtr(), hostMeans.size());
        Kokkos::deep_copy(meanView, hostMeanView);
        coordsView = Coords->getDeviceLocalView(Xpetra::Access::ReadOnly);
        GetOStream(Runtime1) << "Generating nullspace with rotations: dimension = " << nullspaceDim << std::endl;
      } else
        GetOStream(Runtime1) << "Generating canonical nullspace: dimension = " << numPDEs << std::endl;

      nullspace = MultiVectorFactory::Build(A->getDomainMap(), nullspaceDim);

      fillNullspaceVector(nullspace, numPDEs, nullspaceDim, coordsView, meanView);
    }

  } else {
    // On coarser levels always use "Nullspace" as variable name, since it is expected by
    // tentative P factory to be "Nullspace"
    nullspace = currentLevel.Get<RCP<MultiVector>>("Nullspace", GetFactory(nspName).get());  // NOTE: "Nullspace" and nspName mismatch possible here
  }

  // provide "Nullspace" variable on current level (used by TentativePFactory)
  Set(currentLevel, "Nullspace", nullspace);

}  // Build

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void NullspaceFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::fillNullspaceVector(const RCP<MultiVector>& nullspace, LocalOrdinal numPDEs, LocalOrdinal nullspaceDim, CoordsType coordsView, MeanCoordsType meanView) const {
  RCP<BlockedMultiVector> bnsp = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(nullspace);
  if (bnsp.is_null() == true) {
    auto nullspaceView = nullspace->getDeviceLocalView(Xpetra::Access::OverwriteAll);

    int numBlocks = nullspace->getLocalLength() / numPDEs;
    if (nullspaceDim > numPDEs)
      TEUCHOS_TEST_FOR_EXCEPTION(numBlocks != coordsView.extent_int(0), Exceptions::RuntimeError, "MueLu::NullspaceFactory::fillNullspaceVector(): number of coordinates does not match  ndofs/numPDEs.");

    NullspaceFunctor<decltype(nullspaceView), decltype(coordsView), decltype(meanView), LO> nullspaceFunctor(nullspaceView, coordsView, meanView, numPDEs, nullspaceDim);
    Kokkos::parallel_for("MueLu:NullspaceF:Build:for", range_type(0, numBlocks), nullspaceFunctor);

    /*
    // Scale columns to match what Galeri does. Not sure that this is necessary as the qr factorizatoin
    // of the tentative prolongator also takes care of scaling issues. I'm leaving the code here
    // just in case.
    if ( (int) nullspaceDim > numPDEs ) {
    Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> norms2(nullspaceDim);
    nullspace->norm2(norms2);
    Teuchos::Array<Scalar> norms2scalar(nullspaceDim);
    for (int i = 0; i < nullspaceDim; i++)
    norms2scalar[i] = norms2[0] / norms2[i];
    nullspace->scale(norms2scalar);
    }
    */

  } else {
    RCP<const BlockedMap> bmap = bnsp->getBlockedMap();
    for (size_t r = 0; r < bmap->getNumMaps(); r++) {
      Teuchos::RCP<MultiVector> part = bnsp->getMultiVector(r);
      fillNullspaceVector(part, numPDEs, nullspaceDim, coordsView, meanView);
    }
  }
}

}  // namespace MueLu

#endif  // MUELU_NULLSPACEFACTORY_DEF_HPP
