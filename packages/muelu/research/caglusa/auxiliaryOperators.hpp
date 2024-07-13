// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef AUXILIARYOPERATORS_HPP
#define AUXILIARYOPERATORS_HPP

#include <Xpetra_IO.hpp>
#include <MueLu_IOhelpers.hpp>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_HierarchicalOperator_decl.hpp>
#include <Xpetra_HierarchicalOperator_def.hpp>

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_MasterList.hpp>
#include <MueLu_HierarchyManager.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_AmalgamationFactory.hpp>
#include <MueLu_CoalesceDropFactory_kokkos.hpp>
#include <MueLu_ThresholdAFilterFactory.hpp>

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
buildDistanceLaplacian(RCP<const Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >& graph,
                       RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> >& coords) {
  const Scalar ONE       = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar ZERO      = Teuchos::ScalarTraits<Scalar>::zero();
  const LocalOrdinal INV = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
  auto distLapl          = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(graph);

  auto rowMap         = graph->getRowMap();
  auto colMap         = graph->getColMap();
  auto ghosted_coords = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(colMap, coords->getNumVectors());
  ghosted_coords->doImport(*coords, *graph->getImporter(), Xpetra::INSERT);

  {
    auto lcl_coords         = coords->getHostLocalView(Xpetra::Access::ReadOnly);
    auto lcl_ghosted_coords = ghosted_coords->getHostLocalView(Xpetra::Access::ReadOnly);
    auto lcl_distLapl       = distLapl->getLocalMatrixHost();

    // TODO: parallel_for
    for (LocalOrdinal rlid = 0; rlid < lcl_distLapl.numRows(); ++rlid) {
      auto row               = lcl_distLapl.row(rlid);
      Scalar diag            = ZERO;
      LocalOrdinal diagIndex = INV;
      for (LocalOrdinal k = 0; k < row.length; ++k) {
        LocalOrdinal clid = row.colidx(k);
        if (rowMap->getGlobalElement(rlid) == colMap->getGlobalElement(clid)) {
          diagIndex = k;
        } else {
          Scalar dist = ZERO;
          for (size_t j = 0; j < lcl_coords.extent(1); j++) {
            auto s = lcl_coords(rlid, j) - lcl_ghosted_coords(clid, j);
            dist += s * s;
          }
          row.value(k) = ONE / std::sqrt(dist);
          diag -= row.value(k);
        }
      }
      TEUCHOS_ASSERT(diagIndex != INV);
      row.value(diagIndex) = diag;
    }
  }
  distLapl->fillComplete();
  return distLapl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
constructAuxiliaryOperator(RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > op,
                           Teuchos::ParameterList& problemParams) {
  using IO        = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using IOhelpers = MueLu::IOhelpers<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  RCP<Xpetra::HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > hop = rcp_dynamic_cast<Xpetra::HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(op);

  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > auxOp;

  const std::string auxOpStr = problemParams.get<std::string>("auxiliary operator");

  if ((auxOpStr == "near") || (auxOpStr == "distanceLaplacian")) {
    if (hop.is_null())
      auxOp = rcp_dynamic_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(op, true);
    else
      auxOp = hop->nearFieldMatrix();
#ifdef MUELU_HIERARCHICAL_DEBUG
    // CoalesceDropFactory_kokkos assumes fitted row and column maps
    Xpetra::MatrixUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::checkLocalRowMapMatchesColMap(*auxOp);
#endif
    // coordinates
    auto coords = Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node>::ReadMultiVector(problemParams.get<std::string>("coordinates"), op->getRangeMap());

    {
      // apply dropping to auxOp
      Level fineLevel;
      fineLevel.SetFactoryManager(Teuchos::null);
      fineLevel.SetLevelID(0);
      fineLevel.Set("A", auxOp);
      fineLevel.Set("Coordinates", coords);
      fineLevel.Set("DofsPerNode", 1);
      fineLevel.setlib(auxOp->getDomainMap()->lib());
      auto amalgFact = rcp(new AmalgamationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>());
      auto dropFact  = rcp(new CoalesceDropFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>());
      dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

      double dropTol = problemParams.get<double>("drop tolerance");
      // double dropTol = 0.1; // 1D
      // double dropTol = 0.03; // 2D
      std::string dropScheme = "classical";
      dropFact->SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(dropTol));
      dropFact->SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(dropScheme));

      fineLevel.Request("A", dropFact.get());
      fineLevel.Get("A", auxOp, dropFact.get());
    }

    {
      // filter out small entries in auxOp
      Level fineLevel;
      fineLevel.SetFactoryManager(Teuchos::null);
      fineLevel.SetLevelID(0);
      fineLevel.Set("A", auxOp);
      auto filterFact = rcp(new ThresholdAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>("A", 1.0e-8, true, -1));
      fineLevel.Request("A", filterFact.get());
      filterFact->Build(fineLevel);
      auxOp = fineLevel.Get<RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >("A", filterFact.get());
    }

    if (auxOpStr == "distanceLaplacian") {
      // build distance Laplacian using graph of auxOp and coordinates
      auto graph = auxOp->getCrsGraph();
      auxOp      = buildDistanceLaplacian<Scalar, LocalOrdinal, GlobalOrdinal, Node>(graph, coords);
    }

  } else {
    const bool readBinary = problemParams.get<bool>("read binary", false);
    const bool readLocal  = problemParams.get<bool>("read local", false);

    // colmap of auxiliary operator
    auto aux_colmap = IO::ReadMap(problemParams.get<std::string>("aux colmap"), op->getRangeMap()->lib(), op->getRangeMap()->getComm(), readBinary);

    auxOp = IOhelpers::Read(auxOpStr, op->getRangeMap(), aux_colmap, op->getRangeMap(), op->getRangeMap(), true, readBinary, readLocal);
  }

  return auxOp;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
constructHierarchyFromAuxiliary(RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > op,
                                RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > auxH,
                                Teuchos::ParameterList& params,
                                Teuchos::FancyOStream& out) {
  using op_type  = Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using hop_type = Xpetra::HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using mat_type = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  params.set("coarse: max size", 1);
  params.set("max levels", auxH->GetNumLevels());

  const bool implicitTranspose = params.get("transpose: use implicit", MueLu::MasterList::getDefault<bool>("transpose: use implicit"));

  auto hop = rcp_dynamic_cast<hop_type>(op);
  if (!hop.is_null())
    op->describe(out, Teuchos::VERB_EXTREME);

  RCP<Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > H = rcp(new Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>());
  RCP<Level> lvl                                               = H->GetLevel(0);
  auto op_mat                                                  = Teuchos::rcp_dynamic_cast<mat_type>(op);
  if (!op_mat.is_null())
    lvl->Set("A", op_mat);
  else {
    lvl->Set("A", Teuchos::rcp_dynamic_cast<op_type>(op, true));
  }
  // lvl->Set("Coordinates", coords);
  for (int lvlNo = 1; lvlNo < auxH->GetNumLevels(); lvlNo++) {
    RCP<Level> fineLvl = H->GetLevel(lvlNo - 1);
    H->AddNewLevel();
    lvl               = H->GetLevel(lvlNo);
    RCP<Level> auxLvl = auxH->GetLevel(lvlNo);
    // auto mgr = auxLvl->GetFactoryManager();
    // auxLvl->print(std::cout, MueLu::Debug);

    RCP<mat_type> P = auxLvl->Get<RCP<mat_type> >("P");
    lvl->Set("P", P);
    params.sublist("level " + std::to_string(lvlNo)).set("P", P);

    if (!implicitTranspose) {
      TEUCHOS_ASSERT(auxLvl->IsAvailable("R"));
      RCP<mat_type> R = auxLvl->Get<RCP<mat_type> >("R");
      lvl->Set("R", R);
      params.sublist("level " + std::to_string(lvlNo)).set("R", R);
    }

    if (fineLvl->IsType<RCP<op_type> >("A")) {
      RCP<op_type> fineAOp = fineLvl->Get<RCP<op_type> >("A");
      auto fineA           = rcp_dynamic_cast<hop_type>(fineAOp);
      auto coarseA         = fineA->restrict(P);

#ifdef MUELU_HIERARCHICAL_DEBUG
      {
        // Test that the Galerkin product worked
        using MagnitudeType     = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
        const Scalar one        = Teuchos::ScalarTraits<Scalar>::one();
        const MagnitudeType tol = 10000 * Teuchos::ScalarTraits<MagnitudeType>::eps();
        auto testLHS            = MultiVectorFactory::Build(coarseA->getDomainMap(), 1);
        auto testRHS_HOp_coarse = MultiVectorFactory::Build(coarseA->getRangeMap(), 1);
        auto testRHS_HOp_fine   = MultiVectorFactory::Build(coarseA->getRangeMap(), 1);
        auto temp1              = MultiVectorFactory::Build(fineA->getDomainMap(), 1);
        auto temp2              = MultiVectorFactory::Build(fineA->getRangeMap(), 1);
        testLHS->putScalar(one);
        coarseA->apply(*testLHS, *testRHS_HOp_coarse);
        P->apply(*testLHS, *temp1);
        fineA->apply(*temp1, *temp2);
        P->apply(*temp2, *testRHS_HOp_fine, Teuchos::TRANS);
        testRHS_HOp_fine->update(one, *testRHS_HOp_coarse, -one);
        auto norm = testRHS_HOp_fine->getVector(0)->norm2();
        out << "|P^T*op_fine*P*1 - op_H_coarse*1| = " << norm << std::endl;
        TEUCHOS_ASSERT(norm < tol);
      }
#endif

      if ((lvlNo + 1 == auxH->GetNumLevels()) || !coarseA->hasFarField() || coarseA->denserThanDenseMatrix()) {
        // coarseA->describe(out, Teuchos::VERB_EXTREME);

        auto matA = coarseA->toMatrix();

#ifdef MUELU_HIERARCHICAL_DEBUG
        {
          // test that the conversion to Crs format worked
          using MagnitudeType     = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
          const Scalar one        = Teuchos::ScalarTraits<Scalar>::one();
          const MagnitudeType tol = 10000 * Teuchos::ScalarTraits<MagnitudeType>::eps();
          auto testLHS            = MultiVectorFactory::Build(coarseA->getDomainMap(), 1);
          auto testRHS_HOp        = MultiVectorFactory::Build(coarseA->getRangeMap(), 1);
          auto testRHS_dense      = MultiVectorFactory::Build(coarseA->getRangeMap(), 1);
          testLHS->putScalar(one);
          coarseA->apply(*testLHS, *testRHS_HOp);
          matA->apply(*testLHS, *testRHS_dense);
          testRHS_dense->update(one, *testRHS_HOp, -one);
          auto norm = testRHS_dense->getVector(0)->norm2();
          out << "|op_dense*1 - op_H*1| = " << norm << std::endl;
          TEUCHOS_ASSERT(norm < tol);
        }
#endif

        using std::endl;
        using std::setw;
        const size_t numRows   = matA->getRowMap()->getGlobalNumElements();
        const size_t nnz       = matA->getGlobalNumEntries();
        const double nnzPerRow = Teuchos::as<double>(nnz) / numRows;
        std::ostringstream oss;
        oss << std::left;
        // oss << setw(9) << "rows"  << setw(12) << "nnz"  << setw(14) << "nnz/row" << setw(12)  << endl;
        oss << setw(9) << numRows << setw(12) << nnz << setw(14) << nnzPerRow << endl;
        out << oss.str();

        lvl->Set("A", matA);
      } else {
        coarseA->describe(out, Teuchos::VERB_EXTREME, /*printHeader=*/false);
        lvl->Set("A", rcp_dynamic_cast<op_type>(coarseA));
      }
    } else {
      // classical RAP
      RCP<mat_type> fineAmat = fineLvl->Get<RCP<mat_type> >("A");
      // auto fineAmat = rcp_dynamic_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(fineAOp, true);
      Level fineLevel, coarseLevel;
      fineLevel.SetFactoryManager(Teuchos::null);
      coarseLevel.SetFactoryManager(Teuchos::null);
      coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));
      fineLevel.SetLevelID(0);
      coarseLevel.SetLevelID(1);
      fineLevel.Set("A", fineAmat);
      coarseLevel.Set("P", P);
      RCP<RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rapFact = rcp(new RAPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>());
      Teuchos::ParameterList rapList                                      = *(rapFact->GetValidParameterList());
      rapList.set("transpose: use implicit", true);
      rapFact->SetParameterList(rapList);
      coarseLevel.Request("A", rapFact.get());
      RCP<mat_type> matA = coarseLevel.Get<RCP<mat_type> >("A", rapFact.get());

      using std::endl;
      using std::setw;
      const size_t numRows   = matA->getRowMap()->getGlobalNumElements();
      const size_t nnz       = matA->getGlobalNumEntries();
      const double nnzPerRow = Teuchos::as<double>(nnz) / numRows;
      std::ostringstream oss;
      oss << std::left;
      oss << setw(9) << "rows" << setw(12) << "nnz" << setw(14) << "nnz/row" << setw(12) << endl;
      oss << setw(9) << numRows << setw(12) << nnz << setw(14) << nnzPerRow << endl;
      out << oss.str();

      lvl->Set("A", matA);
    }
  }

  RCP<HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node> > mueLuFactory = rcp(new ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node>(params, op->getDomainMap()->getComm()));
  H->setlib(op->getDomainMap()->lib());
  H->SetProcRankVerbose(op->getDomainMap()->getComm()->getRank());
  mueLuFactory->SetupHierarchy(*H);
  H->IsPreconditioner(true);

  return H;
}

}  // namespace MueLu

#endif
