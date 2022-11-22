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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Tpetra_RowMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>

#include <Tpetra_BlockedMap_decl.hpp>
#include <Tpetra_BlockedMap_def.hpp>

#include <Tpetra_BlockedMatrix_decl.hpp>
#include <Tpetra_BlockedMatrix_def.hpp>

#include <Tpetra_HierarchicalOperator_decl.hpp>
#include <Tpetra_HierarchicalOperator_def.hpp>

#include <Xpetra_TpetraBlockedMap.hpp>
#include <Xpetra_TpetraBlockedMatrix.hpp>

#include <Xpetra_HierarchicalOperator_decl.hpp>
#include <Xpetra_HierarchicalOperator_def.hpp>


#include <Xpetra_IO.hpp>
#include <Xpetra_MatrixUtils.hpp>

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_MasterList.hpp>
#include "MueLu_RAPFactory.hpp"
#include "MueLu_Exceptions.hpp"
#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_AmalgamationFactory_kokkos.hpp>
#include <MueLu_CoalesceDropFactory_kokkos.hpp>
#include <MueLu_ThresholdAFilterFactory.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosMueLuAdapter.hpp>
#endif

#define MUELU_HIERARCHICAL_DEBUG

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
struct IOhelpers {

  static
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Read(const std::string&   filename,
       const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rowMap,
       RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > colMap,
       const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > domainMap        = Teuchos::null,
       const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rangeMap         = Teuchos::null,
       const bool           callFillComplete = true,
       const bool           binary           = false,
       const bool           readLocal        = false) {
    using IO = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
    Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A;
    if (readLocal)
      A = IO::ReadLocal(filename, rowMap, colMap, domainMap, rangeMap, callFillComplete, binary);
    else
      A = IO::Read(filename, rowMap, colMap, domainMap, rangeMap, callFillComplete, binary);
    return A;
  }

  static
  Teuchos::RCP<Xpetra::HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  Read(Teuchos::ParameterList& hierarchicalParams,
       RCP< const Teuchos::Comm<int> >& comm) {
    using HOp = Xpetra::HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
    using blocked_matrix_type = typename HOp::blocked_matrix_type;
    using blocked_map_type = typename blocked_matrix_type::blocked_map_type;
    using matrix_type = typename HOp::matrix_type;
    using map_type = typename HOp::map_type;
    using lo_vec_type = typename blocked_map_type::lo_vec_type;

    auto  lib = Xpetra::UseTpetra;
    RCP<HOp>                       op;
    RCP<const map_type>            map, near_colmap, clusterCoeffMap, ghosted_clusterCoeffMap, clusterMap, ghosted_clusterMap;
    RCP<matrix_type>               nearField, basisMatrix, kernelApproximations, kernelBlockGraph;

    std::vector<RCP<blocked_matrix_type> > transferMatrices;
    RCP<lo_vec_type>               clusterSizes;
    RCP<blocked_map_type>          blockedClusterMap, ghosted_blockedClusterMap;
    RCP<blocked_matrix_type>       blockKernelApproximations;

    const bool readBinary = hierarchicalParams.get<bool>("read binary", false);
    const bool readLocal = hierarchicalParams.get<bool>("read local", false);

    using IO = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>;

    // row, domain and range map of the operator
    map = IO::ReadMap(hierarchicalParams.get<std::string>("map"), lib, comm, readBinary);
    // colmap of near field
    near_colmap = IO::ReadMap(hierarchicalParams.get<std::string>("near colmap"), lib, comm, readBinary);
    // 1-to-1 map for the cluster coefficients
    clusterCoeffMap = IO::ReadMap(hierarchicalParams.get<std::string>("coefficient map"), lib, comm, readBinary);
    // overlapping map for the cluster coefficients
    ghosted_clusterCoeffMap = IO::ReadMap(hierarchicalParams.get<std::string>("ghosted coefficient map"), lib, comm, readBinary);
    // 1-to-1 map for the clusters
    clusterMap = IO::ReadMap(hierarchicalParams.get<std::string>("cluster map"), lib, comm, readBinary);
    // overlapping map for the clusters
    ghosted_clusterMap = IO::ReadMap(hierarchicalParams.get<std::string>("ghosted cluster map"), lib, comm, readBinary);

    // blocked cluster map
    clusterSizes = Xpetra::IO<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node>::ReadMultiVector(hierarchicalParams.get<std::string>("gid_cluster_to_gid_coeff"), clusterMap)->getVectorNonConst(0);
    blockedClusterMap = rcp(new blocked_map_type(clusterCoeffMap, clusterSizes));

    // near field interactions
    nearField = Read(hierarchicalParams.get<std::string>("near field matrix"), map, near_colmap, map, map, true, readBinary, readLocal);

    // far field basis expansion coefficients
    basisMatrix = IOhelpers::Read(hierarchicalParams.get<std::string>("basis expansion coefficient matrix"), map, clusterCoeffMap, clusterCoeffMap, map, true, readBinary, readLocal);

    // far field interactions
    kernelApproximations = IOhelpers::Read(hierarchicalParams.get<std::string>("far field interaction matrix"), clusterCoeffMap, ghosted_clusterCoeffMap, clusterCoeffMap, clusterCoeffMap, true, readBinary, readLocal);
    // block graph of far field interactions
    kernelBlockGraph = IOhelpers::Read(hierarchicalParams.get<std::string>("far field interaction matrix")+".block", clusterMap, ghosted_clusterMap, clusterMap, clusterMap, true, readBinary, readLocal);

    {
      auto import = kernelBlockGraph->getCrsGraph()->getImporter();
      RCP<lo_vec_type> ghosted_clusterSizes = Xpetra::VectorFactory<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node>::Build(ghosted_clusterMap);
      ghosted_clusterSizes->doImport(*clusterSizes, *import, Xpetra::INSERT);
      ghosted_blockedClusterMap = rcp(new blocked_map_type(ghosted_clusterCoeffMap, ghosted_clusterSizes));
    }

    blockKernelApproximations = rcp(new blocked_matrix_type(kernelApproximations, kernelBlockGraph, blockedClusterMap, ghosted_blockedClusterMap));

    // Transfer matrices
    auto transfersList = hierarchicalParams.sublist("shift coefficient matrices");
    for (int i = 0; i < transfersList.numParams(); i++) {
      std::string filename = transfersList.get<std::string>(std::to_string(i));
      auto transferPoint = IOhelpers::Read(filename, clusterCoeffMap, clusterCoeffMap, clusterCoeffMap, clusterCoeffMap, true, readBinary, readLocal);
      auto transferBlock = IOhelpers::Read(filename+".block", clusterMap, clusterMap, clusterMap, clusterMap, true, readBinary, readLocal);
      auto transfer = rcp(new blocked_matrix_type(transferPoint, transferBlock, blockedClusterMap));
      transferMatrices.push_back(transfer);
    }

    RCP<Teuchos::ParameterList> params;
    if (hierarchicalParams.isSublist("params")) {
      params = Teuchos::rcpFromRef(hierarchicalParams.sublist("params"));
    }
    op = rcp(new HOp(nearField, blockKernelApproximations, basisMatrix, transferMatrices, params));

    return op;
  }

};


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
buildDistanceLaplacian(RCP<const Xpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& graph,
                       RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType,LocalOrdinal,GlobalOrdinal,Node> >& coords)
{
  #include "MueLu_UseShortNames.hpp"

  const Scalar ONE = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  const LocalOrdinal INV = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
  auto distLapl = MatrixFactory::Build(graph);

  auto rowMap = graph->getRowMap();
  auto colMap = graph->getColMap();
  auto ghosted_coords = MultiVectorFactory::Build(colMap, coords->getNumVectors());
  ghosted_coords->doImport(*coords, *graph->getImporter(), Xpetra::INSERT);

  {
    auto lcl_coords = coords->getHostLocalView(Xpetra::Access::ReadOnly);
    auto lcl_ghosted_coords = ghosted_coords->getHostLocalView(Xpetra::Access::ReadOnly);
    auto lcl_distLapl = distLapl->getLocalMatrixHost();

    // TODO: parallel_for
    for (LocalOrdinal rlid = 0; rlid < lcl_distLapl.numRows(); ++rlid) {
      auto row = lcl_distLapl.row(rlid);
      Scalar diag = ZERO;
      LocalOrdinal diagIndex = INV;
      for (LocalOrdinal k = 0; k < row.length; ++k) {
        LocalOrdinal clid = row.colidx(k);
        if (rowMap->getGlobalElement(rlid) == colMap->getGlobalElement(clid)) {
          diagIndex = k;
        } else {
          Scalar dist = ZERO;
          for (size_t j = 0; j < lcl_coords.extent(1); j++) {
            auto s = lcl_coords(rlid,j) - lcl_ghosted_coords(clid,j);
            dist += s*s;
          }
          row.value(k) = ONE/std::sqrt(dist);
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


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
constructHierarchyFromAuxiliary(RCP<Xpetra::HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > op,
                                RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> > auxH,
                                Teuchos::ParameterList& params,
                                Teuchos::FancyOStream& out) {
  #include "MueLu_UseShortNames.hpp"

  params.set("coarse: max size", 1);
  params.set("max levels", auxH->GetNumLevels());

  const bool implicitTranspose = params.get("transpose: use implicit", MueLu::MasterList::getDefault<bool>("transpose: use implicit"));

  op->describe(out, Teuchos::VERB_EXTREME);

  RCP<Hierarchy> H = rcp(new Hierarchy());
  RCP<Level> lvl = H->GetLevel(0);
  lvl->Set("A", rcp_dynamic_cast<Operator>(op));
  // lvl->Set("Coordinates", coords);
  for(int lvlNo = 1; lvlNo < auxH->GetNumLevels(); lvlNo++) {
    RCP<Level> fineLvl = H->GetLevel(lvlNo-1);
    H->AddNewLevel();
    lvl = H->GetLevel(lvlNo);
    RCP<Level> auxLvl = auxH->GetLevel(lvlNo);
    // auto mgr = auxLvl->GetFactoryManager();
    // auxLvl->print(std::cout, MueLu::Debug);

    RCP<Matrix> P = auxLvl->Get<RCP<Matrix> >("P");
    RCP<Operator> fineAOp = fineLvl->Get<RCP<Operator> >("A");
    lvl->Set("P", P);
    params.sublist("level "+std::to_string(lvlNo)).set("P", P);

    if (!implicitTranspose) {
      TEUCHOS_ASSERT(auxLvl->IsAvailable("R"));
      RCP<Matrix> R = auxLvl->Get<RCP<Matrix> >("R");
      lvl->Set("R", R);
      params.sublist("level "+std::to_string(lvlNo)).set("R", R);
    }

    auto fineA = rcp_dynamic_cast<Xpetra::HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(fineAOp);
    if (!fineA.is_null()) {
      auto coarseA = fineA->restrict(P);

#ifdef MUELU_HIERARCHICAL_DEBUG
      {
        // Test that the Galerkin product worked
        using MagnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
        const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
        const MagnitudeType tol = 10000*Teuchos::ScalarTraits<MagnitudeType>::eps();
        auto testLHS            = MultiVectorFactory::Build(coarseA->getDomainMap(), 1);
        auto testRHS_HOp_coarse = MultiVectorFactory::Build(coarseA->getRangeMap(),  1);
        auto testRHS_HOp_fine   = MultiVectorFactory::Build(coarseA->getRangeMap(),  1);
        auto temp1              = MultiVectorFactory::Build(fineA->getDomainMap(),   1);
        auto temp2              = MultiVectorFactory::Build(fineA->getRangeMap(),    1);
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

      if ((lvlNo+1 == auxH->GetNumLevels()) || !coarseA->hasFarField() || coarseA->denserThanDenseMatrix()) {
        // coarseA->describe(out, Teuchos::VERB_EXTREME);

        auto matA = coarseA->toMatrix();

#ifdef MUELU_HIERARCHICAL_DEBUG
        {
          // test that the conversion to Crs format worked
          using MagnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
          const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
          const MagnitudeType tol = 10000*Teuchos::ScalarTraits<MagnitudeType>::eps();
          auto testLHS       = MultiVectorFactory::Build(coarseA->getDomainMap(), 1);
          auto testRHS_HOp   = MultiVectorFactory::Build(coarseA->getRangeMap(),  1);
          auto testRHS_dense = MultiVectorFactory::Build(coarseA->getRangeMap(),  1);
          testLHS->putScalar(one);
          coarseA->apply(*testLHS, *testRHS_HOp);
          matA->apply(*testLHS, *testRHS_dense);
          testRHS_dense->update(one, *testRHS_HOp, -one);
          auto norm = testRHS_dense->getVector(0)->norm2();
          out << "|op_dense*1 - op_H*1| = " << norm << std::endl;
          TEUCHOS_ASSERT(norm < tol);
        }
#endif

        using std::setw;
        using std::endl;
        const size_t numRows = matA->getRowMap()->getGlobalNumElements();
        const size_t nnz = matA->getGlobalNumEntries();
        const double nnzPerRow = Teuchos::as<double>(nnz)/numRows;
        std::ostringstream oss;
        oss << std::left;
        // oss << setw(9) << "rows"  << setw(12) << "nnz"  << setw(14) << "nnz/row" << setw(12)  << endl;
        oss << setw(9) << numRows << setw(12) << nnz << setw(14) << nnzPerRow << endl;
        out << oss.str();

        lvl->Set("A", matA);
      } else {
        coarseA->describe(out, Teuchos::VERB_EXTREME, /*printHeader=*/false);
        lvl->Set("A", rcp_dynamic_cast<Operator>(coarseA));
      }
    } else {
      // classical RAP
      auto fineAmat = rcp_dynamic_cast<Matrix>(fineAOp, true);
      Level fineLevel, coarseLevel;
      fineLevel.SetFactoryManager(Teuchos::null);
      coarseLevel.SetFactoryManager(Teuchos::null);
      coarseLevel.SetPreviousLevel(rcpFromRef(fineLevel));
      fineLevel.SetLevelID(0);
      coarseLevel.SetLevelID(1);
      fineLevel.Set("A", fineAmat);
      coarseLevel.Set("P", P);
      RCP<RAPFactory> rapFact = rcp(new RAPFactory());
      Teuchos::ParameterList rapList = *(rapFact->GetValidParameterList());
      rapList.set("transpose: use implicit", true);
      rapFact->SetParameterList(rapList);
      coarseLevel.Request("A", rapFact.get());
      RCP<Matrix> matA = coarseLevel.Get<RCP<Matrix> >("A", rapFact.get());

      using std::setw;
      using std::endl;
      const size_t numRows = matA->getRowMap()->getGlobalNumElements();
      const size_t nnz = matA->getGlobalNumEntries();
      const double nnzPerRow = Teuchos::as<double>(nnz)/numRows;
      std::ostringstream oss;
      oss << std::left;
      oss << setw(9) << "rows"  << setw(12) << "nnz"  << setw(14) << "nnz/row" << setw(12)  << endl;
      oss << setw(9) << numRows << setw(12) << nnz << setw(14) << nnzPerRow << endl;
      out << oss.str();

      lvl->Set("A", matA);
    }
  }

  RCP<HierarchyManager> mueLuFactory = rcp(new ParameterListInterpreter(params,op->getDomainMap()->getComm()));
  H->setlib(op->getDomainMap()->lib());
  H->SetProcRankVerbose(op->getDomainMap()->getComm()->getRank());
  mueLuFactory->SetupHierarchy(*H);
  H->IsPreconditioner(true);

  return H;
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
  #include "MueLu_UseShortNames.hpp"

  std::string xmlHierarchical = "hierarchical-1d-mm-global.xml"; clp.setOption("xml",      &xmlHierarchical, "XML describing the hierarchical operator");
  std::string xmlBelos        = "belos.xml";                     clp.setOption("xmlBelos", &xmlBelos,        "XML with Belos parameters");
  std::string xmlMueLu        = "muelu.xml";                     clp.setOption("xmlMueLu", &xmlMueLu,        "XML with MueLu parameters");
  std::string xmlAuxHierarchy = "auxiliary.xml";                 clp.setOption("xmlAux",   &xmlAuxHierarchy, "XML with MueLu parameters for the auxiliary hierarchy");
  bool printTimings  = true; clp.setOption("timings", "notimings", &printTimings,  "print timings to screen");
  bool doTests       = true; clp.setOption("tests",   "notests",   &doTests,       "Test operator using known LHS & RHS.");
  bool doUnPrecSolve = true; clp.setOption("unPrec",  "noUnPrec",  &doUnPrecSolve, "Solve unpreconditioned");
  bool doPrecSolve   = true; clp.setOption("prec",    "noPrec",    &doPrecSolve,   "Solve preconditioned with AMG");

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  Teuchos::RCP<Teuchos::StackedTimer> stacked_timer = rcp(new Teuchos::StackedTimer("Hierarchical Driver"));
  Teuchos::RCP<Teuchos::FancyOStream> verbose_out = Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcpFromRef(std::cout)));
  verbose_out->setShowProcRank(true);
  stacked_timer->setVerboseOstream(verbose_out);
  Teuchos::TimeMonitor::setStackedTimer(stacked_timer);

  using HOp = Xpetra::HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using blocked_matrix_type = typename HOp::blocked_matrix_type;
  using blocked_map_type = typename blocked_matrix_type::blocked_map_type;
  using matrix_type = typename HOp::matrix_type;
  using map_type = typename HOp::map_type;
  using mv_type = typename HOp::mv_type;
  using lo_vec_type = typename blocked_map_type::lo_vec_type;
  using coord_mv = Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType,LocalOrdinal,GlobalOrdinal,Node>;
  using MagnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using IO = Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
  using IOhelpers = IOhelpers<Scalar,LocalOrdinal,GlobalOrdinal,Node>;

  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& out = *fancy;
  out.setOutputToRootOnly(0);
  bool success = true;
  const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  const MagnitudeType tol = 10000*Teuchos::ScalarTraits<MagnitudeType>::eps();

  Teuchos::ParameterList hierarchicalParams;
  Teuchos::updateParametersFromXmlFileAndBroadcast(xmlHierarchical, Teuchos::Ptr<Teuchos::ParameterList>(&hierarchicalParams), *comm);

  RCP<HOp> op;
  {
    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Read hierarchical matrix")));
    op = IOhelpers::Read(hierarchicalParams, comm);
  }

  out << "Compression: " << op->getCompression() << " of dense matrix."<< std::endl;

  RCP<const map_type> map = op->getDomainMap();
  RCP<matrix_type>    auxOp;
  RCP<mv_type>        X_ex, RHS, X;
  RCP<coord_mv>       coords;
  {
    // Read in auxiliary stuff

    // coordinates
    coords = Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::coordinateType,LocalOrdinal,GlobalOrdinal,Node>::ReadMultiVector(hierarchicalParams.get<std::string>("coordinates"), map);

    // Auxiliary matrix used for multigrid construction
    const std::string auxOpStr = hierarchicalParams.get<std::string>("auxiliary operator");
    if ((auxOpStr == "near") || (auxOpStr == "distanceLaplacian")) {
      auxOp = op->nearFieldMatrix();
#ifdef HAVE_MUELU_DEBUG
      // CoalesceDropFactory_kokkos assumes fitted row and column maps
      Xpetra::MatrixUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::checkLocalRowMapMatchesColMap(*auxOp);
#endif

      {
        // apply dropping to auxOp
        Level fineLevel;
        fineLevel.SetFactoryManager(Teuchos::null);
        fineLevel.SetLevelID(0);
        fineLevel.Set("A",auxOp);
        fineLevel.Set("Coordinates",coords);
        fineLevel.Set("DofsPerNode",1);
        fineLevel.setlib(auxOp->getDomainMap()->lib());
        auto amalgFact = rcp(new AmalgamationFactory_kokkos());
        auto dropFact = rcp(new CoalesceDropFactory_kokkos());
        dropFact->SetFactory("UnAmalgamationInfo", amalgFact);

        double dropTol = hierarchicalParams.get<double>("drop tolerance");
        // double dropTol = 0.1; // 1D
        // double dropTol = 0.03; // 2D
        std::string dropScheme = "classical";
        dropFact->SetParameter("aggregation: drop tol",Teuchos::ParameterEntry(dropTol));
        dropFact->SetParameter("aggregation: drop scheme",Teuchos::ParameterEntry(dropScheme));

        fineLevel.Request("A",dropFact.get());
        fineLevel.Get("A", auxOp, dropFact.get());
      }

      {
        // filter out small entries in auxOp
        Level fineLevel;
        fineLevel.SetFactoryManager(Teuchos::null);
        fineLevel.SetLevelID(0);
        fineLevel.Set("A",auxOp);
        auto filterFact = rcp(new ThresholdAFilterFactory("A", 1.0e-8, true, -1));
        fineLevel.Request("A",filterFact.get());
        filterFact->Build(fineLevel);
        auxOp = fineLevel.Get< RCP<Matrix> >("A",filterFact.get());
      }

      if (auxOpStr == "distanceLaplacian") {
        // build distance Laplacian using graph of auxOp and coordinates
        auto graph = auxOp->getCrsGraph();
        auxOp = buildDistanceLaplacian<Scalar,LocalOrdinal,GlobalOrdinal,Node>(graph, coords);
      }

    } else {
      const bool readBinary = hierarchicalParams.get<bool>("read binary", false);
      const bool readLocal = hierarchicalParams.get<bool>("read local", false);

      // colmap of auxiliary operator
      RCP<const map_type> aux_colmap = IO::ReadMap(hierarchicalParams.get<std::string>("aux colmap"), lib, comm, readBinary);

      auxOp = IOhelpers::Read(auxOpStr, map, aux_colmap, map, map, true, readBinary, readLocal);
    }

    // known pair of LHS, RHS
    X_ex = IO::ReadMultiVector(hierarchicalParams.get<std::string>("exact solution"), map);
    RHS  = IO::ReadMultiVector(hierarchicalParams.get<std::string>("right-hand side"), map);
    // solution vector
    X    = MultiVectorFactory::Build(map, 1);

  }

  if (doTests) {
    // Some simple apply tests
    Scalar opX_exRHS, MopX_exRHS, MopTX_exRHS;
    {
      op->apply(*X_ex, *X);

      // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("X.mtx", *X);
      X->update(one, *RHS, -one);
      opX_exRHS = X->getVector(0)->norm2();
      out << "|op*X_ex - RHS| = " << opX_exRHS << std::endl;
      // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("diff.mtx", *X);
    }

    {
      op->apply(*X_ex, *X, Teuchos::NO_TRANS, -one);

      // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("X2.mtx", *X);
      X->update(one, *RHS, one);
      MopX_exRHS = X->getVector(0)->norm2();
      out << "|(-op)*X_ex + RHS| = " << MopX_exRHS << std::endl;
      // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("diff2.mtx", *X);
    }

    {
      op->apply(*X_ex, *X, Teuchos::TRANS, -one);

      // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("X2.mtx", *X);
      X->update(one, *RHS, one);
      MopTX_exRHS = X->getVector(0)->norm2();
      out << "|(-op^T)*X_ex + RHS| = " << MopTX_exRHS << std::endl;
      // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("diff2.mtx", *X);
    }

    TEUCHOS_ASSERT(opX_exRHS < tol);
    TEUCHOS_ASSERT(MopX_exRHS < tol);
    TEUCHOS_ASSERT(MopTX_exRHS < tol);
  }

#ifdef HAVE_MUELU_BELOS
  Teuchos::ParameterList belosParams;
  Teuchos::updateParametersFromXmlFileAndBroadcast(xmlBelos, Teuchos::Ptr<Teuchos::ParameterList>(&belosParams), *comm);
  if (doUnPrecSolve) {
    // Solve linear system using unpreconditioned Krylov method
    out << "\n*********************************************************\n";
    out << "Unpreconditioned Krylov method\n";
    out << "*********************************************************\n\n";

    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Unpreconditioned solve")));

    using MV = typename HOp::mv_type;
    using OP = Belos::OperatorT<MV>;

    X->putScalar(zero);
    RCP<OP> belosOp = rcp(new Belos::XpetraOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(op));
    RCP<Belos::LinearProblem<Scalar, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<Scalar, MV, OP>(belosOp, X, RHS));

    std::string belosType = "Pseudoblock CG";
    auto belosSolverList = rcpFromRef(belosParams.sublist(belosType));

    bool set = belosProblem->setProblem();
    if (set == false) {
      throw MueLu::Exceptions::RuntimeError("ERROR:  Belos::LinearProblem failed to set up correctly!");
    }

    // Create an iterative solver manager
    Belos::SolverFactory<Scalar, MV, OP> solverFactory;
    RCP< Belos::SolverManager<Scalar, MV, OP> > solver = solverFactory.create(belosType, belosSolverList);
    solver->setProblem(belosProblem);

    // Perform solve
    Belos::ReturnType ret = solver->solve();
    int numIts = solver->getNumIters();

    // Get the number of iterations for this solve.
    out << "Number of iterations performed for this solve: " << numIts << std::endl;

    // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("X.mtx", *X);
    X->update(one, *X_ex, -one);
    out << "|X-X_ex| = " << X->getVector(0)->norm2() << std::endl << std::endl;

    success &= (ret == Belos::Converged);

  }
#endif // HAVE_MUELU_BELOS

  if (doPrecSolve) {
    // Solve linear system using a AMG preconditioned Krylov method

    RCP<Hierarchy> auxH, H;

    {
      ////////////////////////////////////////////////////////////////
      // Build the auxiliary hierarchy
      out << "\n*********************************************************\n";
      out << "Building the auxiliary hierarchy\n";
      out << "*********************************************************\n\n";

      Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Construct auxiliary hierarchy")));

      Teuchos::ParameterList auxParams;
      Teuchos::updateParametersFromXmlFileAndBroadcast(xmlAuxHierarchy, Teuchos::Ptr<Teuchos::ParameterList>(&auxParams), *comm);
      auxParams.set("hierarchy label", "Auxiliary");
      auxParams.sublist("user data").set("Coordinates", coords);

      auxH = MueLu::CreateXpetraPreconditioner(auxOp, auxParams);
    }

    {
      ////////////////////////////////////////////////////////////////
      // Construct the main hierarchy
      out << "\n*********************************************************\n";
      out << "Building the main hierarchy\n";
      out << "*********************************************************\n\n";

      Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Construct hierarchy")));

      Teuchos::ParameterList params;
      Teuchos::updateParametersFromXmlFileAndBroadcast(xmlMueLu, Teuchos::Ptr<Teuchos::ParameterList>(&params), *comm);
      params.sublist("user data").set("Coordinates", coords);

      H = constructHierarchyFromAuxiliary(op, auxH, params, out);
    }


#ifdef HAVE_MUELU_BELOS
    {
      ////////////////////////////////////////////////////////////////
      // Set up the Krylov solver

      Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Preconditioned solve")));

      using MV = typename HOp::mv_type;
      using OP = Belos::OperatorT<MV>;

      X->putScalar(zero);
      RCP<OP> belosOp = rcp(new Belos::XpetraOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(op));
      RCP<OP> belosPrec = rcp(new Belos::MueLuOp <SC, LO, GO, NO>(H));
      RCP<Belos::LinearProblem<Scalar, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<Scalar, MV, OP>(belosOp, X, RHS));

      std::string belosType = "Pseudoblock CG";
      auto belosSolverList = rcpFromRef(belosParams.sublist(belosType));

      belosProblem->setRightPrec(belosPrec);

      bool set = belosProblem->setProblem();
      if (set == false) {
        throw MueLu::Exceptions::RuntimeError("ERROR:  Belos::LinearProblem failed to set up correctly!");
      }

      // Create an iterative solver manager
      Belos::SolverFactory<Scalar, MV, OP> solverFactory;
      RCP< Belos::SolverManager<Scalar, MV, OP> > solver = solverFactory.create(belosType, belosSolverList);
      solver->setProblem(belosProblem);

      // Perform solve
      Belos::ReturnType ret = solver->solve();
      int numIts = solver->getNumIters();

      // Get the number of iterations for this solve.
      out << "Number of iterations performed for this solve: " << numIts << std::endl;

      // Xpetra::IO<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Write("X.mtx", *X);
      X->update(one, *X_ex, -one);
      out << "|X-X_ex| = " << X->getVector(0)->norm2() << std::endl;

      success &= (ret == Belos::Converged);
    }
#endif // HAVE_MUELU_BELOS

  }

  stacked_timer->stop("Hierarchical Driver");
  Teuchos::StackedTimer::OutputOptions options;
  options.output_fraction = options.output_histogram = options.output_minmax = true;
  if (printTimings)
    stacked_timer->report(out, comm, options);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} //main

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc,argv);
}
