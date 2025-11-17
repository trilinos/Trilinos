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

  auto absDc = MatrixFactory::BuildCopy(Dc);
  absDc->setAllToScalar(one);

  RCP<Matrix> temp1 = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*absD, false, *absPn, false, GetOStream(Statistics2), true, true);
  temp1->setAllToScalar(one);
  RCP<Matrix> temp2 = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*temp1, false, *absDc, true, GetOStream(Statistics2), true, true);

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
    auto importer = absDc->getCrsGraph()->getImporter();
    if (importer.is_null()) {
      singleParentGhosted = singleParent;
    } else {
      singleParentGhosted = MultiVectorFactory::Build(absDc->getColMap(), 1);
      singleParentGhosted->doImport(*singleParent, *importer, Xpetra::INSERT);
    }

    auto lclSingleParent = singleParentGhosted->getLocalViewDevice(Tpetra::Access::ReadOnly);

    // Filter matrix using criterion
    filtered = Xpetra::applyFilter_LID(
        temp2,
        KOKKOS_LAMBDA(const LocalOrdinal row,
                      const LocalOrdinal col,
                      const Matrix::impl_scalar_type val) {
          return ((ATS::magnitude(val - 2.0) < eps) || ((lclSingleParent(col, 0) == 1.0) && (ATS::magnitude(val - 1.0) < eps)));
        });
  }

  // Xpetra::IO<SC, LO, GO, NO>::Write("pattern.code", *filtered);

  Set(coarseLevel, "Ppattern", filtered->getCrsGraph());
}

}  // namespace MueLu

#endif  // MUELU_EDGEPROLONGATORPATTERNFACTORY_DEF_HPP
