// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_EDGEPROLONGATORPATTERNFACTORY_DEF_HPP
#define MUELU_EDGEPROLONGATORPATTERNFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>
// #include <Xpetra_IO.hpp>

#include "MueLu_EdgeProlongatorPatternFactory_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
//#include "MueLu_Utilities.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> EdgeProlongatorPatternFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("FineD0", Teuchos::null, "Generating factory for the fine discrete gradient");
  validParamList->set<RCP<const FactoryBase> >("CoarseD0", Teuchos::null, "Generating factory for the coarse discrete gradient");
  validParamList->set<RCP<const FactoryBase> >("PnodalEmin", Teuchos::null, "Generating factory for the nodal prolongator");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void EdgeProlongatorPatternFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  Input(fineLevel, "D0", "FineD0");
  Input(coarseLevel, "D0", "CoarseD0");
  Input(coarseLevel, "PnodalEmin");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void EdgeProlongatorPatternFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "EdgeProlongatorPattern", coarseLevel);

  RCP<Matrix> D  = Get<RCP<Matrix> >(fineLevel, "D0", "FineD0");
  RCP<Matrix> Dc = Get<RCP<Matrix> >(coarseLevel, "D0", "CoarseD0");
  RCP<Matrix> Pn = Get<RCP<Matrix> >(coarseLevel, "PnodalEmin");

  const auto one = Teuchos::ScalarTraits<Scalar>::one();

  // |FineD| * |Pnodal| * |CoarseD^T|

  auto absD = MatrixFactory::BuildCopy(D);
  absD->setAllToScalar(one);

  auto absPn = MatrixFactory::BuildCopy(Pn);
  absPn->setAllToScalar(one);

  RCP<Matrix> absD_absPn = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*absD, false, *absPn, false, GetOStream(Statistics2), true, true);
  absD_absPn->setAllToScalar(one);

  // If we rebalanced then Dc lives on a smaller communicator than D.
  // Since we need to perform matrix-matrix multiplications with Dc, we construct a version of it that lives on the same communicator.
  auto comm = absD_absPn->getRowMap()->getComm();
  if (Dc.is_null() || Dc->getRowMap()->getComm()->getSize() < comm->getSize()) {
    auto lib = absD_absPn->getRowMap()->lib();
    if (Dc.is_null()) {
      Kokkos::View<GlobalOrdinal*, typename Node::memory_space> dummy("", 0);
      auto big_coarse_nodal_map    = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, -1, dummy, 0, comm);
      auto big_coarse_edge_map     = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, -1, dummy, 0, comm);
      auto big_coarse_nodal_colmap = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, -1, dummy, 0, comm);

      typename Matrix::local_matrix_device_type dummyLocalMatrix;
      Dc = MatrixFactory::Build(dummyLocalMatrix, big_coarse_edge_map, big_coarse_nodal_colmap, big_coarse_nodal_map, big_coarse_edge_map);

    } else {
      auto big_coarse_nodal_map    = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, -1, Dc->getDomainMap()->getMyGlobalIndicesDevice(), 0, comm);
      auto big_coarse_edge_map     = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, -1, Dc->getRangeMap()->getMyGlobalIndicesDevice(), 0, comm);
      auto big_coarse_nodal_colmap = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, -1, Dc->getColMap()->getMyGlobalIndicesDevice(), 0, comm);

      Dc = MatrixFactory::Build(Dc->getLocalMatrixDevice(), big_coarse_edge_map, big_coarse_nodal_colmap, big_coarse_nodal_map, big_coarse_edge_map);
    }
  }
  auto absDc = MatrixFactory::BuildCopy(Dc);
  absDc->setAllToScalar(one);

  RCP<Matrix> absD_absPn_absDcT = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*absD_absPn, false, *absDc, true, GetOStream(Statistics2), true, true);

  RCP<Matrix> filtered;
  {
#if KOKKOS_VERSION >= 40799
    using ATS           = KokkosKernels::ArithTraits<Scalar>;
    using magnitudeType = typename ATS::magnitudeType;
    using magATS        = KokkosKernels::ArithTraits<magnitudeType>;
#else
    using ATS           = Kokkos::ArithTraits<Scalar>;
    using magnitudeType = typename ATS::magnitudeType;
    using magATS        = Kokkos::ArithTraits<magnitudeType>;
#endif
    auto eps = magATS::epsilon();

    RCP<MultiVector> oneVec = MultiVectorFactory::Build(absDc->getDomainMap(), 1);
    oneVec->putScalar(one);
    RCP<MultiVector> singleParent = MultiVectorFactory::Build(absDc->getRowMap(), 1);
    absDc->apply(*oneVec, *singleParent, Teuchos::NO_TRANS);
    // ghost singleParent
    RCP<MultiVector> singleParentGhosted;
    auto importer = absD_absPn_absDcT->getCrsGraph()->getImporter();
    if (importer.is_null()) {
      singleParentGhosted = singleParent;
    } else {
      singleParentGhosted = MultiVectorFactory::Build(importer->getTargetMap(), 1);
      singleParentGhosted->doImport(*singleParent, *importer, Xpetra::INSERT);
    }

    auto lclSingleParent = singleParentGhosted->getLocalViewDevice(Tpetra::Access::ReadOnly);

    // Filter matrix using criterion
    filtered = Xpetra::applyFilter_LID(
        absD_absPn_absDcT,
        KOKKOS_LAMBDA(const LocalOrdinal row,
                      const LocalOrdinal col,
                      const typename Matrix::impl_scalar_type val) {
          return ((ATS::magnitude(val - 2.0) < eps) || ((lclSingleParent(col, 0) == 1.0) && (ATS::magnitude(val - 1.0) < eps)));
        });
  }

  // Xpetra::IO<SC, LO, GO, NO>::Write("pattern.code", *filtered);

  Set(coarseLevel, "Ppattern", filtered->getCrsGraph());
}

}  // namespace MueLu

#endif  // MUELU_EDGEPROLONGATORPATTERNFACTORY_DEF_HPP
