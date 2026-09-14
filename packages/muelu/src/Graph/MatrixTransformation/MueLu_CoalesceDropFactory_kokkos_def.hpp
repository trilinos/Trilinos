// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_COALESCEDROPFACTORY_KOKKOS_DEF_HPP
#define MUELU_COALESCEDROPFACTORY_KOKKOS_DEF_HPP

#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "Xpetra_Matrix.hpp"
#include "Xpetra_MultiVectorFactory.hpp"

#include "MueLu_CoalesceDropFactory_kokkos_decl.hpp"

#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_BoundaryDetection.hpp"
#include "MueLu_DroppingCommon.hpp"
#include "MueLu_MatrixConstruction.hpp"

// The different dropping algorithms are split up over several translation units. This speeds up compilation and also avoids launch latencies on GPU.
#include "MueLu_ScalarDroppingBase.hpp"
#include "MueLu_ScalarDroppingClassical.hpp"
#include "MueLu_ScalarDroppingDistanceLaplacian.hpp"
#include "MueLu_VectorDroppingBase.hpp"
#include "MueLu_VectorDroppingClassical.hpp"
#include "MueLu_VectorDroppingDistanceLaplacian.hpp"

namespace MueLu {

namespace MaterialDistanceDiagnostics {

template <class Scalar>
double realValue(const Scalar& value) {
  using non_const_scalar_type = typename std::remove_const<Scalar>::type;
  return Teuchos::ScalarTraits<non_const_scalar_type>::real(value);
}

struct Tensor2DInfo {
  double a;
  double b;
  double d;
  double lambdaMax;
  double lambdaMin;
  double vx;
  double vy;
};

template <class material_view_type>
Tensor2DInfo getTensor2DInfo(const material_view_type& material, const size_t row) {
  using scalar_type = typename material_view_type::value_type;

  const double a  = realValue<scalar_type>(material(row, 0));
  const double b0 = realValue<scalar_type>(material(row, 1));
  const double b1 = realValue<scalar_type>(material(row, 2));
  const double d  = realValue<scalar_type>(material(row, 3));
  const double b  = 0.5 * (b0 + b1);

  const double delta     = std::sqrt((a - d) * (a - d) + 4.0 * b * b);
  const double lambdaMax = 0.5 * (a + d + delta);
  const double lambdaMin = 0.5 * (a + d - delta);

  double vx = b;
  double vy = lambdaMax - a;
  if (std::abs(vx) + std::abs(vy) <= 100.0 * std::numeric_limits<double>::epsilon()) {
    if (a >= d) {
      vx = 1.0;
      vy = 0.0;
    } else {
      vx = 0.0;
      vy = 1.0;
    }
  }
  const double norm = std::sqrt(vx * vx + vy * vy);

  Tensor2DInfo info;
  info.a         = a;
  info.b         = b;
  info.d         = d;
  info.lambdaMax = lambdaMax;
  info.lambdaMin = lambdaMin;
  info.vx        = vx / norm;
  info.vy        = vy / norm;
  return info;
}

inline double tensorContrast(const Tensor2DInfo& rowTensor, const Tensor2DInfo& colTensor) {
  const double floor    = 100.0 * std::numeric_limits<double>::min();
  const double maxRatio = std::max(rowTensor.lambdaMax, colTensor.lambdaMax) / std::max(std::min(rowTensor.lambdaMax, colTensor.lambdaMax), floor);
  const double minRatio = std::max(rowTensor.lambdaMin, colTensor.lambdaMin) / std::max(std::min(rowTensor.lambdaMin, colTensor.lambdaMin), floor);
  const double da       = rowTensor.a - colTensor.a;
  const double db       = rowTensor.b - colTensor.b;
  const double dd       = rowTensor.d - colTensor.d;
  const double diffNorm = std::sqrt(da * da + 2.0 * db * db + dd * dd);
  const double rowNorm  = std::sqrt(rowTensor.a * rowTensor.a + 2.0 * rowTensor.b * rowTensor.b + rowTensor.d * rowTensor.d);
  const double colNorm  = std::sqrt(colTensor.a * colTensor.a + 2.0 * colTensor.b * colTensor.b + colTensor.d * colTensor.d);
  const double relDiff  = diffNorm / std::max(0.5 * (rowNorm + colNorm), floor);
  return std::max(std::max(maxRatio, minRatio), 1.0 + relDiff);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PrintMaterialDistanceDiagnostics(Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                                      Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node>> coords,
                                      Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> material,
                                      const Kokkos::View<DecisionType*, typename Node::device_type::memory_space>& results,
                                      const double coverageThreshold,
                                      const Factory& factory) {
  using magnitude_type = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using coords_type    = Xpetra::MultiVector<magnitude_type, LocalOrdinal, GlobalOrdinal, Node>;
  using material_type  = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  const size_t spatialDim = coords->getNumVectors();
  if (spatialDim != 2 || material->getNumVectors() != 4) {
    if (factory.IsPrint(Statistics1)) {
      factory.GetOStream(Statistics1) << "Material distance diagnostics skipped: currently implemented for 2D tensor materials only" << std::endl;
    }
    return;
  }

  auto importer = A.getCrsGraph()->getImporter();
  Teuchos::RCP<coords_type> ghostedCoords;
  Teuchos::RCP<material_type> ghostedMaterial;
  if (!importer.is_null()) {
    ghostedCoords = Xpetra::MultiVectorFactory<magnitude_type, LocalOrdinal, GlobalOrdinal, Node>::Build(importer->getTargetMap(), coords->getNumVectors(), false);
    ghostedCoords->doImport(*coords, *importer, Xpetra::INSERT);
    ghostedMaterial = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(importer->getTargetMap(), material->getNumVectors(), false);
    ghostedMaterial->doImport(*material, *importer, Xpetra::INSERT);
  } else {
    ghostedCoords   = coords;
    ghostedMaterial = material;
  }

  auto coordsHost          = coords->getLocalViewHost(Tpetra::Access::ReadOnly);
  auto ghostedCoordsHost   = ghostedCoords->getLocalViewHost(Tpetra::Access::ReadOnly);
  auto materialHost        = material->getLocalViewHost(Tpetra::Access::ReadOnly);
  auto ghostedMaterialHost = ghostedMaterial->getLocalViewHost(Tpetra::Access::ReadOnly);
  auto resultsHost         = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), results);

  const double materialContrastThreshold = factory.GetParameterList().get<double>("aggregation: material distance diagnostics contrast");
  const double noCoverage                = std::numeric_limits<double>::max();

  GlobalOrdinal localNodes                     = 0;
  GlobalOrdinal localPoorCoverageNodes         = 0;
  GlobalOrdinal localInterfaceNodes            = 0;
  GlobalOrdinal localPoorInterfaceNodes        = 0;
  GlobalOrdinal localKeptSameMaterial          = 0;
  GlobalOrdinal localDroppedSameMaterial       = 0;
  GlobalOrdinal localKeptCrossMaterial         = 0;
  GlobalOrdinal localDroppedCrossMaterial      = 0;
  GlobalOrdinal localRowsWithNoKeptOffdiag     = 0;
  GlobalOrdinal localPoorSameMaterialCoverage  = 0;
  GlobalOrdinal localPoorInterfaceSameCoverage = 0;
  double localCoverageSum                      = 0.0;
  double localSameMaterialCoverageSum          = 0.0;
  double localInterfaceCoverageSum             = 0.0;
  double localInterfaceSameMaterialCoverageSum = 0.0;
  double localCoverageMin                      = noCoverage;
  double localSameMaterialCoverageMin          = noCoverage;
  double localInterfaceCoverageMin             = noCoverage;
  double localInterfaceSameMaterialCoverageMin = noCoverage;

  size_t resultOffset = 0;
  for (LocalOrdinal row = 0; row < Teuchos::as<LocalOrdinal>(A.getRowMap()->getLocalNumElements()); ++row) {
    Teuchos::ArrayView<const LocalOrdinal> indices;
    Teuchos::ArrayView<const Scalar> vals;
    A.getLocalRowView(row, indices, vals);

    const auto rowTensor            = getTensor2DInfo(materialHost, row);
    double bestCoverage             = 0.0;
    double bestSameMaterialCoverage = 0.0;
    bool interfaceNode              = false;
    bool hasKeptOffdiag             = false;

    for (LocalOrdinal k = 0; k < indices.size(); ++k) {
      const LocalOrdinal col = indices[k];
      if (row == col)
        continue;

      const double dx         = realValue<magnitude_type>(ghostedCoordsHost(col, 0)) - realValue<magnitude_type>(coordsHost(row, 0));
      const double dy         = realValue<magnitude_type>(ghostedCoordsHost(col, 1)) - realValue<magnitude_type>(coordsHost(row, 1));
      const double h          = std::sqrt(dx * dx + dy * dy);
      const auto colTensor    = getTensor2DInfo(ghostedMaterialHost, col);
      const bool sameMaterial = tensorContrast(rowTensor, colTensor) <= materialContrastThreshold;
      if (h > 0.0) {
        const double coverage = std::abs(dx * rowTensor.vx + dy * rowTensor.vy) / h;
        bestCoverage          = std::max(bestCoverage, coverage);
        if (sameMaterial)
          bestSameMaterialCoverage = std::max(bestSameMaterialCoverage, coverage);
      }

      interfaceNode = interfaceNode || !sameMaterial;

      const bool kept = resultsHost(resultOffset + k) == KEEP;
      hasKeptOffdiag  = hasKeptOffdiag || kept;
      if (sameMaterial) {
        if (kept)
          ++localKeptSameMaterial;
        else
          ++localDroppedSameMaterial;
      } else {
        if (kept)
          ++localKeptCrossMaterial;
        else
          ++localDroppedCrossMaterial;
      }
    }

    ++localNodes;
    localCoverageSum += bestCoverage;
    localSameMaterialCoverageSum += bestSameMaterialCoverage;
    localCoverageMin             = std::min(localCoverageMin, bestCoverage);
    localSameMaterialCoverageMin = std::min(localSameMaterialCoverageMin, bestSameMaterialCoverage);
    if (bestCoverage * bestCoverage < coverageThreshold)
      ++localPoorCoverageNodes;
    if (bestSameMaterialCoverage * bestSameMaterialCoverage < coverageThreshold)
      ++localPoorSameMaterialCoverage;
    if (interfaceNode) {
      ++localInterfaceNodes;
      localInterfaceCoverageSum += bestCoverage;
      localInterfaceSameMaterialCoverageSum += bestSameMaterialCoverage;
      localInterfaceCoverageMin             = std::min(localInterfaceCoverageMin, bestCoverage);
      localInterfaceSameMaterialCoverageMin = std::min(localInterfaceSameMaterialCoverageMin, bestSameMaterialCoverage);
      if (bestCoverage * bestCoverage < coverageThreshold)
        ++localPoorInterfaceNodes;
      if (bestSameMaterialCoverage * bestSameMaterialCoverage < coverageThreshold)
        ++localPoorInterfaceSameCoverage;
    }
    if (!hasKeptOffdiag)
      ++localRowsWithNoKeptOffdiag;

    resultOffset += indices.size();
  }

  auto comm = A.getRowMap()->getComm();
  Teuchos::Array<GlobalOrdinal> localCounts(11);
  localCounts[0]  = localNodes;
  localCounts[1]  = localPoorCoverageNodes;
  localCounts[2]  = localInterfaceNodes;
  localCounts[3]  = localPoorInterfaceNodes;
  localCounts[4]  = localKeptSameMaterial;
  localCounts[5]  = localDroppedSameMaterial;
  localCounts[6]  = localKeptCrossMaterial;
  localCounts[7]  = localDroppedCrossMaterial;
  localCounts[8]  = localRowsWithNoKeptOffdiag;
  localCounts[9]  = localPoorSameMaterialCoverage;
  localCounts[10] = localPoorInterfaceSameCoverage;
  Teuchos::Array<GlobalOrdinal> globalCounts(localCounts.size());
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, Teuchos::as<int>(localCounts.size()), localCounts.getRawPtr(), globalCounts.getRawPtr());

  Teuchos::Array<double> localSums(4);
  localSums[0] = localCoverageSum;
  localSums[1] = localInterfaceCoverageSum;
  localSums[2] = localSameMaterialCoverageSum;
  localSums[3] = localInterfaceSameMaterialCoverageSum;
  Teuchos::Array<double> globalSums(localSums.size());
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, Teuchos::as<int>(localSums.size()), localSums.getRawPtr(), globalSums.getRawPtr());

  double globalCoverageMin                      = noCoverage;
  double globalInterfaceCoverageMin             = noCoverage;
  double globalSameMaterialCoverageMin          = noCoverage;
  double globalInterfaceSameMaterialCoverageMin = noCoverage;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &localCoverageMin, &globalCoverageMin);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &localInterfaceCoverageMin, &globalInterfaceCoverageMin);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &localSameMaterialCoverageMin, &globalSameMaterialCoverageMin);
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &localInterfaceSameMaterialCoverageMin, &globalInterfaceSameMaterialCoverageMin);

  if (factory.IsPrint(Statistics1)) {
    const double globalCoverageMean                      = globalCounts[0] > 0 ? globalSums[0] / Teuchos::as<double>(globalCounts[0]) : 0.0;
    const double globalInterfaceCoverageMean             = globalCounts[2] > 0 ? globalSums[1] / Teuchos::as<double>(globalCounts[2]) : 0.0;
    const double globalSameMaterialCoverageMean          = globalCounts[0] > 0 ? globalSums[2] / Teuchos::as<double>(globalCounts[0]) : 0.0;
    const double globalInterfaceSameMaterialCoverageMean = globalCounts[2] > 0 ? globalSums[3] / Teuchos::as<double>(globalCounts[2]) : 0.0;
    if (globalCoverageMin == noCoverage)
      globalCoverageMin = 0.0;
    if (globalInterfaceCoverageMin == noCoverage)
      globalInterfaceCoverageMin = 0.0;
    if (globalSameMaterialCoverageMin == noCoverage)
      globalSameMaterialCoverageMin = 0.0;
    if (globalInterfaceSameMaterialCoverageMin == noCoverage)
      globalInterfaceSameMaterialCoverageMin = 0.0;

    factory.GetOStream(Statistics1)
        << "Material distance diagnostics:" << std::endl
        << "  directional coverage c_i: min=" << globalCoverageMin
        << ", mean=" << globalCoverageMean
        << ", poor(c_i^2 < " << coverageThreshold << ")=" << globalCounts[1] << "/" << globalCounts[0] << std::endl
        << "  same-material directional coverage: min=" << globalSameMaterialCoverageMin
        << ", mean=" << globalSameMaterialCoverageMean
        << ", poor(c_i^2 < " << coverageThreshold << ")=" << globalCounts[9] << "/" << globalCounts[0] << std::endl
        << "  tensor-contrast interface nodes: " << globalCounts[2]
        << ", poor interface coverage=" << globalCounts[3] << "/" << globalCounts[2]
        << ", poor same-material interface coverage=" << globalCounts[10] << "/" << globalCounts[2]
        << ", interface mean=" << globalInterfaceCoverageMean
        << ", interface same-material mean=" << globalInterfaceSameMaterialCoverageMean
        << ", interface min=" << globalInterfaceCoverageMin
        << ", interface same-material min=" << globalInterfaceSameMaterialCoverageMin << std::endl
        << "  offdiag entries kept same/cross: " << globalCounts[4] << "/" << globalCounts[6]
        << ", dropped same/cross: " << globalCounts[5] << "/" << globalCounts[7] << std::endl
        << "  rows with no kept offdiag entry: " << globalCounts[8] << std::endl;
  }
}

}  // namespace MaterialDistanceDiagnostics

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> CoalesceDropFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("aggregation: use blocking");
  SET_VALID_ENTRY("aggregation: drop tol");
  SET_VALID_ENTRY("aggregation: use ml scaling of drop tol");
  SET_VALID_ENTRY("aggregation: Dirichlet threshold");
  SET_VALID_ENTRY("aggregation: greedy Dirichlet");
  SET_VALID_ENTRY("aggregation: row sum drop tol");
  SET_VALID_ENTRY("aggregation: strength-of-connection: matrix");
  SET_VALID_ENTRY("aggregation: strength-of-connection: measure");
  SET_VALID_ENTRY("aggregation: drop scheme");
  SET_VALID_ENTRY("aggregation: block diagonal: interleaved blocksize");
  SET_VALID_ENTRY("aggregation: distance laplacian metric");
  SET_VALID_ENTRY("aggregation: material distance: interface penalty");
  SET_VALID_ENTRY("aggregation: material distance: interface penalty strength");
  SET_VALID_ENTRY("aggregation: material distance: interface penalty floor");
  SET_VALID_ENTRY("aggregation: material distance: interface penalty shape weight");
  SET_VALID_ENTRY("aggregation: material distance diagnostics");
  SET_VALID_ENTRY("aggregation: material distance diagnostics threshold");
  SET_VALID_ENTRY("aggregation: material distance diagnostics contrast");
  SET_VALID_ENTRY("aggregation: Minv scheme");
  SET_VALID_ENTRY("aggregation: distance laplacian directional weights");
  SET_VALID_ENTRY("aggregation: dropping may create Dirichlet");
#ifdef HAVE_MUELU_COALESCEDROP_ALLOW_OLD_PARAMETERS
  SET_VALID_ENTRY("aggregation: distance laplacian algo");
  SET_VALID_ENTRY("aggregation: classical algo");
#endif
  SET_VALID_ENTRY("aggregation: symmetrize graph after dropping");
  SET_VALID_ENTRY("aggregation: coloring: use color graph");
  SET_VALID_ENTRY("aggregation: coloring: localize color graph");

  SET_VALID_ENTRY("filtered matrix: use lumping");
  SET_VALID_ENTRY("filtered matrix: reuse graph");
  SET_VALID_ENTRY("filtered matrix: reuse eigenvalue");

  SET_VALID_ENTRY("filtered matrix: use root stencil");
  SET_VALID_ENTRY("filtered matrix: lumping choice");
  SET_VALID_ENTRY("filtered matrix: use spread lumping");
  SET_VALID_ENTRY("filtered matrix: spread lumping diag dom growth factor");
  SET_VALID_ENTRY("filtered matrix: spread lumping diag dom cap");
  SET_VALID_ENTRY("filtered matrix: Dirichlet threshold");
  SET_VALID_ENTRY("filtered matrix: count negative diagonals");

#undef SET_VALID_ENTRY
  validParamList->set<bool>("lightweight wrap", true, "Experimental option for lightweight graph access");
#ifndef HAVE_MUELU_COALESCEDROP_ALLOW_OLD_PARAMETERS
  validParamList->getEntry("aggregation: drop scheme").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("point-wise", "cut-drop"))));
#else
  validParamList->getEntry("aggregation: drop scheme").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("point-wise", "cut-drop", "signed classical sa", "classical", "distance laplacian", "signed classical", "block diagonal", "block diagonal classical", "block diagonal distance laplacian", "block diagonal signed classical", "block diagonal colored signed classical", "signed classical distance laplacian", "signed classical sa distance laplacian"))));
  validParamList->getEntry("aggregation: classical algo").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("default", "unscaled cut", "scaled cut", "scaled cut symmetric"))));
  validParamList->getEntry("aggregation: distance laplacian algo").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("default", "unscaled cut", "scaled cut", "scaled cut symmetric"))));
#endif
  validParamList->getEntry("aggregation: strength-of-connection: matrix").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("A", "distance laplacian", "MinvA"))));
  validParamList->getEntry("aggregation: strength-of-connection: measure").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("smoothed aggregation", "signed smoothed aggregation", "signed ruge-stueben", "unscaled"))));
  validParamList->getEntry("aggregation: distance laplacian metric").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("unweighted", "material"))));
  validParamList->getEntry("aggregation: material distance: interface penalty").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("none", "log-frobenius"))));
  validParamList->getEntry("aggregation: Minv scheme").setValidator(rcp(new Teuchos::StringValidator(Teuchos::tuple<std::string>("spai", "fsai"))));

  validParamList->set<RCP<const FactoryBase>>("A", Teuchos::null, "Generating factory of the matrix A");
  validParamList->set<RCP<const FactoryBase>>("UnAmalgamationInfo", Teuchos::null, "Generating factory for UnAmalgamationInfo");
  validParamList->set<RCP<const FactoryBase>>("Coordinates", Teuchos::null, "Generating factory for Coordinates");
  validParamList->set<RCP<const FactoryBase>>("BlockNumber", Teuchos::null, "Generating factory for BlockNumber");
  validParamList->set<RCP<const FactoryBase>>("Material", Teuchos::null, "Generating factory for Material");
  validParamList->set<RCP<const FactoryBase>>("M", Teuchos::null, "Generating factory for M");
  validParamList->set<RCP<const FactoryBase>>("Minv", Teuchos::null, "Generating factory for Minv");
  validParamList->set<RCP<const FactoryBase>>("MinvA", Teuchos::null, "Generating factory for MinvA");

  // Make sure we don't recursively validate options for project auxiliary matrices
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("project auxiliary matrices", norecurse, "matrices that will be projected on coarse levels");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoalesceDropFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
  Input(currentLevel, "A");
  Input(currentLevel, "UnAmalgamationInfo");

  const ParameterList& pL = GetParameterList();

  std::string socUsesMatrix = pL.get<std::string>("aggregation: strength-of-connection: matrix");
  bool needCoords           = (socUsesMatrix == "distance laplacian");
  const bool needM          = (socUsesMatrix == "MinvA");
#ifdef HAVE_MUELU_COALESCEDROP_ALLOW_OLD_PARAMETERS
  std::string droppingMethod = pL.get<std::string>("aggregation: drop scheme");
  needCoords |= (droppingMethod.find("distance laplacian") != std::string::npos);
#endif
  if (needCoords) {
    Input(currentLevel, "Coordinates");
    std::string distLaplMetric = pL.get<std::string>("aggregation: distance laplacian metric");
    if (distLaplMetric == "material")
      Input(currentLevel, "Material");
  }
  if (needM && (currentLevel.GetLevelID() != 0)) {
    if (pL.isSublist("project auxiliary matrices")) {
      auto projectList = pL.sublist("project auxiliary matrices");
      if (projectList.isParameter("MinvA"))
        Input(currentLevel, "MinvA");
      else if (projectList.isParameter("Minv"))
        Input(currentLevel, "Minv");
      else if (projectList.isParameter("M"))
        Input(currentLevel, "M");
    }
  }

  bool useBlocking = pL.get<bool>("aggregation: use blocking");
#ifdef HAVE_MUELU_COALESCEDROP_ALLOW_OLD_PARAMETERS
  useBlocking |= (droppingMethod.find("block diagonal") != std::string::npos);
#endif
  if (useBlocking) {
    Input(currentLevel, "BlockNumber");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoalesceDropFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Build(Level& currentLevel) const {
  auto A = Get<RCP<Matrix>>(currentLevel, "A");
  TEUCHOS_TEST_FOR_EXCEPTION(A->GetFixedBlockSize() % A->GetStorageBlockSize() != 0, Exceptions::RuntimeError, "A->GetFixedBlockSize() needs to be a multiple of A->GetStorageBlockSize()");
  LO blkSize = A->GetFixedBlockSize() / A->GetStorageBlockSize();

  std::tuple<GlobalOrdinal, GlobalOrdinal, boundary_nodes_type> results;
  if (blkSize == 1)
    results = BuildScalar(currentLevel);
  else
    results = BuildVector(currentLevel);

  if (GetVerbLevel() & Statistics1) {
    GlobalOrdinal numDropped = std::get<0>(results);
    GlobalOrdinal numKept    = std::get<1>(results);
    auto boundaryNodes       = std::get<2>(results);

    GO numLocalBoundaryNodes = 0;

    Kokkos::parallel_reduce(
        "MueLu:CoalesceDropF:Build:bnd", range_type(0, boundaryNodes.extent(0)),
        KOKKOS_LAMBDA(const LO i, GO& n) {
          if (boundaryNodes(i))
            n++;
        },
        numLocalBoundaryNodes);

    if (IsPrint(Statistics1)) {
      auto comm = A->getRowMap()->getComm();

      std::vector<GlobalOrdinal> localStats = {numLocalBoundaryNodes, numDropped, numKept};
      std::vector<GlobalOrdinal> globalStats(3);
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 3, localStats.data(), globalStats.data());

      GO numGlobalBoundaryNodes = globalStats[0];
      GO numGlobalDropped       = globalStats[1];
      GO numGlobalKept          = globalStats[2];
      GO numGlobalTotal         = numGlobalKept + numGlobalDropped;

      GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
      if (numGlobalTotal != 0) {
        GetOStream(Statistics1) << "Number of dropped entries: "
                                << numGlobalDropped << "/" << numGlobalTotal
                                << " (" << 100 * Teuchos::as<double>(numGlobalDropped) / Teuchos::as<double>(numGlobalTotal) << "%)" << std::endl;
      }
    }
  }
}

template <class local_matrix_type, class boundary_nodes_view, class... Functors>
void runBoundaryFunctors(local_matrix_type& lclA, boundary_nodes_view& boundaryNodes, Functors&... functors) {
  using local_ordinal_type = typename local_matrix_type::ordinal_type;
  using execution_space    = typename local_matrix_type::execution_space;
  using range_type         = Kokkos::RangePolicy<local_ordinal_type, execution_space>;
  auto range               = range_type(0, boundaryNodes.extent(0));
  auto boundaries          = BoundaryDetection::BoundaryFunctor(lclA, functors...);
  Kokkos::parallel_for("CoalesceDrop::BoundaryDetection", range, boundaries);
}

template <class magnitudeType>
void translateOldAlgoParam(const Teuchos::ParameterList& pL, std::string& droppingMethod, bool& useBlocking, std::string& socUsesMatrix, std::string& socUsesMeasure, bool& symmetrizeDroppedGraph, bool& generateColoringGraph, magnitudeType& threshold, MueLu::MatrixConstruction::lumpingType& lumpingChoice) {
  std::set<std::string> validDroppingMethods = {"piece-wise", "cut-drop"};

  if (!pL.get<bool>("filtered matrix: use lumping")) lumpingChoice = MueLu::MatrixConstruction::no_lumping;

  if (validDroppingMethods.find(droppingMethod) == validDroppingMethods.end()) {
    std::string algo                     = droppingMethod;
    std::string classicalAlgoStr         = pL.get<std::string>("aggregation: classical algo");
    std::string distanceLaplacianAlgoStr = pL.get<std::string>("aggregation: distance laplacian algo");

    // Remove prefix "block diagonal" from algo
    if (algo.find("block diagonal") == 0) {
      useBlocking = true;
      algo        = algo.substr(14);
      if (algo != "") {
        algo = algo.substr(1);
      }
    }

    if ((algo == "classical") || (algo == "signed classical sa") || (algo == "signed classical") || (algo == "colored signed classical")) {
      socUsesMatrix = "A";

      if (algo == "classical") {
        socUsesMeasure = "smoothed aggregation";
      } else if (algo == "signed classical sa") {
        socUsesMeasure = "signed smoothed aggregation";
      } else if (algo == "signed classical") {
        socUsesMeasure = "signed ruge-stueben";
      } else if (algo == "colored signed classical") {
        socUsesMeasure        = "signed ruge-stueben";
        generateColoringGraph = true;
      }

      if (classicalAlgoStr == "default")
        droppingMethod = "point-wise";
      else if (classicalAlgoStr == "unscaled cut") {
        socUsesMeasure = "unscaled";
        droppingMethod = "cut-drop";
      } else if (classicalAlgoStr == "scaled cut") {
        droppingMethod = "cut-drop";
      } else if (classicalAlgoStr == "scaled cut symmetric") {
        droppingMethod         = "cut-drop";
        symmetrizeDroppedGraph = true;
      }
    } else if ((algo == "distance laplacian") || (algo == "signed classical sa distance laplacian") || (algo == "signed classical distance laplacian")) {
      socUsesMatrix = "distance laplacian";

      if (algo == "distance laplacian") {
        socUsesMeasure = "smoothed aggregation";
      } else if (algo == "signed classical sa distance laplacian") {
        socUsesMeasure = "signed smoothed aggregation";
      } else if (algo == "signed classical distance laplacian") {
        socUsesMeasure = "signed ruge-stueben";
      }

      if (distanceLaplacianAlgoStr == "default")
        droppingMethod = "point-wise";
      else if (distanceLaplacianAlgoStr == "unscaled cut") {
        socUsesMeasure = "unscaled";
        droppingMethod = "cut-drop";
      } else if (distanceLaplacianAlgoStr == "scaled cut") {
        droppingMethod = "cut-drop";
      } else if (distanceLaplacianAlgoStr == "scaled cut symmetric") {
        droppingMethod         = "cut-drop";
        symmetrizeDroppedGraph = true;
      }
    } else if (algo == "") {
      // algo was "block diagonal", but we process and remove the "block diagonal" part
      socUsesMatrix = "A";
      threshold     = Teuchos::ScalarTraits<magnitudeType>::zero();
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::tuple<GlobalOrdinal, GlobalOrdinal, typename MueLu::LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Node>::boundary_nodes_type> CoalesceDropFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildScalar(Level& currentLevel) const {
  FactoryMonitor m(*this, "BuildScalar", currentLevel);

  using MatrixType        = Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using GraphType         = Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type = typename MatrixType::local_matrix_device_type;
  using local_graph_type  = typename GraphType::local_graph_device_type;
  using rowptr_type       = typename local_graph_type::row_map_type::non_const_type;
  using entries_type      = typename local_graph_type::entries_type::non_const_type;
  using values_type       = typename local_matrix_type::values_type::non_const_type;
  using device_type       = typename Node::device_type;
  using memory_space      = typename device_type::memory_space;
  using results_view_type = Kokkos::View<DecisionType*, memory_space>;
  using magnitudeType     = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using doubleMultiVector = Xpetra::MultiVector<magnitudeType, LO, GO, NO>;

  typedef Teuchos::ScalarTraits<Scalar> STS;
  const magnitudeType zero = Teuchos::ScalarTraits<magnitudeType>::zero();

  auto A = Get<RCP<Matrix>>(currentLevel, "A");

  const bool needToBuildFilteredA = currentLevel.IsRequested("A", this);

  //////////////////////////////////////////////////////////////////////
  // Process parameterlist
  const ParameterList& pL = GetParameterList();

  // Boundary detection
  const magnitudeType dirichletThreshold       = STS::magnitude(as<SC>(pL.get<double>("aggregation: Dirichlet threshold")));
  const magnitudeType rowSumTol                = as<magnitudeType>(pL.get<double>("aggregation: row sum drop tol"));
  const LocalOrdinal dirichletNonzeroThreshold = 1;

  // Dropping
  bool useBlocking                    = pL.get<bool>("aggregation: use blocking");
  std::string droppingMethod          = pL.get<std::string>("aggregation: drop scheme");
  std::string socUsesMatrix           = pL.get<std::string>("aggregation: strength-of-connection: matrix");
  std::string socUsesMeasure          = pL.get<std::string>("aggregation: strength-of-connection: measure");
  std::string distanceLaplacianMetric = pL.get<std::string>("aggregation: distance laplacian metric");
  std::string MinvScheme              = pL.get<std::string>("aggregation: Minv scheme");
  bool symmetrizeDroppedGraph         = pL.get<bool>("aggregation: symmetrize graph after dropping");
  magnitudeType threshold;
  // If we're doing the ML-style halving of the drop tol at each level, we do that here.
  if (pL.get<bool>("aggregation: use ml scaling of drop tol"))
    threshold = pL.get<double>("aggregation: drop tol") / pow(2.0, currentLevel.GetLevelID());
  else
    threshold = as<magnitudeType>(pL.get<double>("aggregation: drop tol"));
  bool aggregationMayCreateDirichlet = pL.get<bool>("aggregation: dropping may create Dirichlet");

  // Fill
  const bool reuseGraph      = pL.get<bool>("filtered matrix: reuse graph") && needToBuildFilteredA;
  const bool reuseEigenvalue = pL.get<bool>("filtered matrix: reuse eigenvalue");

  const bool useRootStencil                            = pL.get<bool>("filtered matrix: use root stencil");
  const bool useSpreadLumping                          = pL.get<bool>("filtered matrix: use spread lumping");
  const std::string lumpingChoiceString                = pL.get<std::string>("filtered matrix: lumping choice");
  MueLu::MatrixConstruction::lumpingType lumpingChoice = MueLu::MatrixConstruction::no_lumping;
  if (needToBuildFilteredA) {
    if (lumpingChoiceString == "diag lumping")
      lumpingChoice = MueLu::MatrixConstruction::diag_lumping;
    else if (lumpingChoiceString == "distributed lumping")
      lumpingChoice = MueLu::MatrixConstruction::distributed_lumping;
  }

  const magnitudeType filteringDirichletThreshold = as<magnitudeType>(pL.get<double>("filtered matrix: Dirichlet threshold"));

  // coloring graph
  bool generateColoringGraph         = pL.get<bool>("aggregation: coloring: use color graph");
  const bool localizeColoringGraph   = pL.get<bool>("aggregation: coloring: localize color graph");
  const bool symmetrizeColoringGraph = true;

#ifdef HAVE_MUELU_COALESCEDROP_ALLOW_OLD_PARAMETERS
  translateOldAlgoParam(pL, droppingMethod, useBlocking, socUsesMatrix, socUsesMeasure, symmetrizeDroppedGraph, generateColoringGraph, threshold, lumpingChoice);
#endif

  {
    std::stringstream ss;
    ss << "dropping scheme = \"" << droppingMethod << "\", strength-of-connection measure = \"" << socUsesMeasure << "\", strength-of-connection matrix = \"" << socUsesMatrix << "\", ";
    if (socUsesMatrix == "distance laplacian")
      ss << "distance laplacian metric = \"" << distanceLaplacianMetric << "\", ";
    ss << "threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << ", useBlocking = " << useBlocking;
    ss << ", symmetrizeDroppedGraph = " << symmetrizeDroppedGraph << std::endl;

    GetOStream(Runtime0) << ss.str();
  }

  TEUCHOS_ASSERT(!useRootStencil);
  TEUCHOS_ASSERT(!useSpreadLumping);
  TEUCHOS_ASSERT((lumpingChoice != MueLu::MatrixConstruction::distributed_lumping) || !reuseGraph);
  if (droppingMethod == "cut-drop")
    TEUCHOS_TEST_FOR_EXCEPTION(threshold > 1.0, Exceptions::RuntimeError, "For cut-drop algorithms, \"aggregation: drop tol\" = " << threshold << ", needs to be <= 1.0");

  //////////////////////////////////////////////////////////////////////
  // We perform four sweeps over the rows of A:
  // Pass 1: detection of boundary nodes
  // Pass 2: diagonal extraction
  // Pass 3: drop decision for each entry and construction of the rowptr of the filtered matrix
  // Pass 4: fill of the filtered matrix
  //
  // Pass 1 and 3 apply a sequence of criteria to each row of the matrix.

  // TODO: We could merge pass 1 and 2.

  // Compute matrix that is used for dropping
  RCP<Matrix> A_drop;
  if (threshold != zero) {
    if ((socUsesMatrix == "A") || (socUsesMatrix == "MinvA")) {
      // Get matrix used for dropping
      if (socUsesMatrix == "A")
        A_drop = A;
      else if (socUsesMatrix == "MinvA") {
        bool storeMinvOnLevel = false, storeMinvAOnLevel = false;
        if (pL.isSublist("project auxiliary matrices")) {
          auto projectList = pL.sublist("project auxiliary matrices");
          if (projectList.isParameter("MinvA"))
            storeMinvAOnLevel = true;
          else {
            if (projectList.isParameter("Minv"))
              storeMinvOnLevel = true;
            else
              storeMinvAOnLevel = true;  // nothing found in project list, so default to storing MinvA
          }
        } else
          storeMinvAOnLevel = true;  // no project list given, so default to storing MinvA

        if (IsAvailable(currentLevel, "MinvA")) {
          A_drop = Get<RCP<Matrix>>(currentLevel, "MinvA");
        } else {
          RCP<Matrix> Minv;
          bool MueLu_Minv = IsAvailable(currentLevel, "Minv");
          bool User_Minv  = currentLevel.IsAvailable("Minv", NoFactory::get());
          if (MueLu_Minv || User_Minv) {
            if (MueLu_Minv)
              Minv = Get<RCP<Matrix>>(currentLevel, "Minv");
            else
              Minv = currentLevel.Get<RCP<Matrix>>("Minv", NoFactory::get());
          } else {  // get M and create Minv

            RCP<Matrix> M;
            if (currentLevel.GetLevelID() == 0) {
              M = currentLevel.Get<RCP<Matrix>>("M", NoFactory::get());
            } else {
              M = Get<RCP<Matrix>>(currentLevel, "M");
            }
            // Find Dirichlets in A to stick corresponding Dirichlets in M
#ifdef moreExperimentsToSeeIfHelpful
            auto boundaryNodes = MueLu::Utilities<SC, LO, GO, NO>::DetectDirichletRows_kokkos(*A, 8.0 * Teuchos::ScalarTraits<magnitudeType>::eps());
            if (GetVerbLevel() & Statistics1) {
              int nAbc = 0;
              for (size_t iii = 0; iii < A->getLocalNumRows(); iii++)
                if (boundaryNodes[iii]) nAbc++;
              auto MbcBefore = MueLu::Utilities<SC, LO, GO, NO>::DetectDirichletRows_kokkos(*M, 8.0 * Teuchos::ScalarTraits<magnitudeType>::eps());
              int nMbcBefore = 0;
              for (size_t iii = 0; iii < M->getLocalNumRows(); iii++)
                if (MbcBefore[iii]) nMbcBefore++;
              GetOStream(Statistics1) << "MueLu_CoalesceDropFactory_kokkos_def: # A bcs = " << nAbc << ", # M bcs before MueLu enforcement = " << nMbcBefore;
            }
            MueLu::Utilities<SC, LO, GO, NO>::ApplyOAZToMatrixRows(M, boundaryNodes);
            if (GetVerbLevel() & Statistics1) {
              auto MbcAfter = MueLu::Utilities<SC, LO, GO, NO>::DetectDirichletRows_kokkos(*M, 8.0 * Teuchos::ScalarTraits<magnitudeType>::eps());
              int nMbcAfter = 0;
              for (size_t iii = 0; iii < M->getLocalNumRows(); iii++)
                if (MbcAfter[iii]) nMbcAfter++;
              GetOStream(Statistics1) << ", # M bcs after MueLu enforcement = " << nMbcAfter << std::endl;
            }
#endif

            // Create Minv via sparse approximate inverse

            Minv = Utilities::SPAI(M, MinvScheme);

            if (storeMinvOnLevel) currentLevel.Set("Minv", Minv);
          }  // finished if/else (currentLevel.IsAvailable("Minv", *mtf)

          // build MinvA matrix with same sparsity pattern as A.
          A_drop      = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(A);
          auto params = Teuchos::rcp(new Teuchos::ParameterList());
          params->set("MM Throw For Non-Existent Entries", false);
          A_drop = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Minv, false, *A, false, A_drop, GetOStream(Statistics2), true, true, std::string("MinvA"), params);

          // Enforce Dirichlet BCs by zeroing out off-diagonal rows and columns and putting a 1 on diag of MinvA
          auto boundaryNodes = MueLu::Utilities<SC, LO, GO, NO>::DetectDirichletRows_kokkos(*A, 8.0 * Teuchos::ScalarTraits<magnitudeType>::eps());
          auto dirichletCols = MueLu::Utilities<SC, LO, GO, NO>::DetectDirichletCols(*A, boundaryNodes);
          MueLu::Utilities<SC, LO, GO, NO>::ApplyOAZToMatrixRows(A_drop, boundaryNodes);
          MueLu::Utilities<SC, LO, GO, NO>::ZeroDirichletCols(A_drop, dirichletCols, Teuchos::ScalarTraits<SC>::zero(), true);

          if (storeMinvAOnLevel) currentLevel.Set("MinvA", A_drop);
        }  // finished if/else  (currentLevel.IsAvailable("MinvA", NoFactory::get()))
      }    // else if (socUsesMatrix == "MinvA") {
    }
  }

  auto crsA  = toCrsMatrix(A);
  auto lclA  = crsA->getLocalMatrixDevice();
  auto range = range_type(0, lclA.numRows());

  //////////////////////////////////////////////////////////////////////
  // Pass 1: Detect boundary nodes
  //
  // The following criteria are available:
  // - BoundaryDetection::PointDirichletFunctor
  //   Marks rows as Dirichlet based on value threshold and number of off-diagonal entries
  // - BoundaryDetection::RowSumFunctor
  //   Marks rows as Dirichlet bases on row-sum criterion

  // Dirichlet nodes
  auto boundaryNodes = boundary_nodes_type("boundaryNodes", lclA.numRows());  // initialized to false
  {
    SubFactoryMonitor mBoundary(*this, "Boundary detection", currentLevel);

    // macro that applies boundary detection functors
    auto dirichlet_detection = BoundaryDetection::PointDirichletFunctor(lclA, boundaryNodes, dirichletThreshold, dirichletNonzeroThreshold);

    if (rowSumTol <= 0.) {
      runBoundaryFunctors(lclA, boundaryNodes, dirichlet_detection);
    } else {
      auto apply_rowsum = BoundaryDetection::RowSumFunctor(lclA, boundaryNodes, rowSumTol);
      runBoundaryFunctors(lclA, boundaryNodes, dirichlet_detection, apply_rowsum);
    }
  }
  // In what follows, boundaryNodes can still still get modified if aggregationMayCreateDirichlet == true.
  // Otherwise we're now done with it now.

  //////////////////////////////////////////////////////////////////////
  // Pass 2 & 3: Diagonal extraction and determine dropping and construct
  //             rowptr of filtered matrix
  //
  // The following criteria are available:
  // - Misc::PointwiseDropBoundaryFunctor
  //   Drop all rows that have been marked as Dirichlet
  // - Misc::DropOffRankFunctor
  //   Drop all entries that are off-rank
  // - ClassicalDropping::DropFunctor
  //   Classical dropping
  // - DistanceLaplacian::DropFunctor
  //   Distance Laplacian dropping
  // - Misc::KeepDiagonalFunctor
  //   Mark diagonal as KEEP
  // - Misc::MarkSingletonFunctor
  //   Mark singletons after dropping as Dirichlet
  // - Misc::BlockDiagonalizeFunctor
  //   Drop coupling between blocks
  //
  // For the block diagonal variants we first block diagonalized and then apply "blocksize = 1" algorithms.

  // rowptr of filtered A
  auto filtered_rowptr = rowptr_type("filtered_rowptr", lclA.numRows() + 1);
  // Number of nonzeros of filtered A
  LocalOrdinal nnz_filtered = 0;
  // dropping decisions for each entry
  auto results = results_view_type("results", lclA.nnz());  // initialized to UNDECIDED
  {
    SubFactoryMonitor mDropping(*this, "Dropping decisions", currentLevel);

    if (threshold != zero) {
      if ((socUsesMatrix == "A") || (socUsesMatrix == "MinvA")) {
        if (socUsesMeasure == "unscaled") {
          ScalarDroppingClassical<Scalar, LocalOrdinal, GlobalOrdinal, Node, Misc::UnscaledMeasure>::runDroppingFunctors_on_A(*A_drop, results, filtered_rowptr, nnz_filtered, boundaryNodes, droppingMethod, threshold,
                                                                                                                              aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, currentLevel, *this);
        } else if (socUsesMeasure == "smoothed aggregation") {
          ScalarDroppingClassical<Scalar, LocalOrdinal, GlobalOrdinal, Node, Misc::SmoothedAggregationMeasure>::runDroppingFunctors_on_A(*A_drop, results, filtered_rowptr, nnz_filtered, boundaryNodes, droppingMethod, threshold,
                                                                                                                                         aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, currentLevel, *this);
        } else if (socUsesMeasure == "signed ruge-stueben") {
          ScalarDroppingClassical<Scalar, LocalOrdinal, GlobalOrdinal, Node, Misc::SignedRugeStuebenMeasure>::runDroppingFunctors_on_A(*A_drop, results, filtered_rowptr, nnz_filtered, boundaryNodes, droppingMethod, threshold,
                                                                                                                                       aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, currentLevel, *this);
        } else if (socUsesMeasure == "signed smoothed aggregation") {
          ScalarDroppingClassical<Scalar, LocalOrdinal, GlobalOrdinal, Node, Misc::SignedSmoothedAggregationMeasure>::runDroppingFunctors_on_A(*A_drop, results, filtered_rowptr, nnz_filtered, boundaryNodes, droppingMethod, threshold,
                                                                                                                                               aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, currentLevel, *this);
        }
      } else if (socUsesMatrix == "distance laplacian") {
        auto coords = Get<RCP<doubleMultiVector>>(currentLevel, "Coordinates");
        if (socUsesMeasure == "unscaled") {
          ScalarDroppingDistanceLaplacian<Scalar, LocalOrdinal, GlobalOrdinal, Node, Misc::UnscaledMeasure>::runDroppingFunctors_on_dlap(*A, results, filtered_rowptr, nnz_filtered, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, distanceLaplacianMetric, currentLevel, *this);
        } else if (socUsesMeasure == "smoothed aggregation") {
          ScalarDroppingDistanceLaplacian<Scalar, LocalOrdinal, GlobalOrdinal, Node, Misc::SmoothedAggregationMeasure>::runDroppingFunctors_on_dlap(*A, results, filtered_rowptr, nnz_filtered, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, distanceLaplacianMetric, currentLevel, *this);
        } else if (socUsesMeasure == "signed ruge-stueben") {
          ScalarDroppingDistanceLaplacian<Scalar, LocalOrdinal, GlobalOrdinal, Node, Misc::SignedRugeStuebenMeasure>::runDroppingFunctors_on_dlap(*A, results, filtered_rowptr, nnz_filtered, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, distanceLaplacianMetric, currentLevel, *this);
        } else if (socUsesMeasure == "signed smoothed aggregation") {
          ScalarDroppingDistanceLaplacian<Scalar, LocalOrdinal, GlobalOrdinal, Node, Misc::SignedSmoothedAggregationMeasure>::runDroppingFunctors_on_dlap(*A, results, filtered_rowptr, nnz_filtered, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, distanceLaplacianMetric, currentLevel, *this);
        }
      }
    } else {
      Kokkos::deep_copy(results, KEEP);

      if (symmetrizeDroppedGraph) {
        auto drop_boundaries = Misc::PointwiseSymmetricDropBoundaryFunctor(*A, boundaryNodes, results);
        ScalarDroppingBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::template runDroppingFunctors<>(*A, results, filtered_rowptr, nnz_filtered, useBlocking, currentLevel, *this, drop_boundaries);
      } else {
        auto no_op = Misc::NoOpFunctor<LocalOrdinal>();
        ScalarDroppingBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::template runDroppingFunctors<>(*A, results, filtered_rowptr, nnz_filtered, useBlocking, currentLevel, *this, no_op);
      }
    }

    if (symmetrizeDroppedGraph) {
      auto symmetrize = Misc::SymmetrizeFunctor(lclA, results);
      ScalarDroppingBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::template runDroppingFunctors<>(*A, results, filtered_rowptr, nnz_filtered, useBlocking, currentLevel, *this, symmetrize);
    }
  }
  GO numDropped = lclA.nnz() - nnz_filtered;
  GO numGlobalDropped;
  Teuchos::reduceAll(*A->getRowMap()->getComm(), Teuchos::REDUCE_SUM, 1, &numDropped, &numGlobalDropped);

  if (pL.get<bool>("aggregation: material distance diagnostics") && distanceLaplacianMetric == "material" && IsPrint(Statistics1)) {
    auto coords   = Get<RCP<doubleMultiVector>>(currentLevel, "Coordinates");
    auto material = Get<RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>(currentLevel, "Material");
    MaterialDistanceDiagnostics::PrintMaterialDistanceDiagnostics(*A, coords, material, results,
                                                                  pL.get<double>("aggregation: material distance diagnostics threshold"),
                                                                  *this);
  }
  // We now know the number of entries of filtered A and have the final rowptr.

  //////////////////////////////////////////////////////////////////////
  // Pass 4: Create local matrix for filtered A
  //
  // Dropped entries are optionally lumped to the diagonal.

  RCP<LWGraph_kokkos> graph;
  RCP<Matrix> filteredA;
  if (numGlobalDropped > 0) {
    SubFactoryMonitor mFill(*this, "Filtered matrix fill", currentLevel);

    local_graph_type lclGraph;
    local_matrix_type lclFilteredA;
    auto colidx = entries_type("entries", nnz_filtered);
    lclGraph    = local_graph_type(colidx, filtered_rowptr);
    if (reuseGraph) {
      // needToBuildFilteredA is true
      filteredA    = MatrixFactory::BuildCopy(A);
      lclFilteredA = filteredA->getLocalMatrixDevice();
    } else if (needToBuildFilteredA) {
      auto values  = values_type("values", nnz_filtered);
      lclFilteredA = local_matrix_type("filteredA", lclGraph, lclA.numCols());
    }

    if (lumpingChoice != MueLu::MatrixConstruction::no_lumping) {
      if (reuseGraph) {
        auto fillFunctor = MatrixConstruction::PointwiseFillReuseFunctor<local_matrix_type, local_graph_type, true>(lclA, results, lclFilteredA, lclGraph, filteringDirichletThreshold);
        Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_lumped_reuse", range, fillFunctor);
      } else {
        if (lumpingChoice == MueLu::MatrixConstruction::diag_lumping) {
          auto fillFunctor = MatrixConstruction::PointwiseFillNoReuseFunctor<local_matrix_type, MueLu::MatrixConstruction::diag_lumping>(lclA, results, lclGraph, lclFilteredA, filteringDirichletThreshold);
          Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_lumped_noreuse", range, fillFunctor);
        } else if (lumpingChoice == MueLu::MatrixConstruction::distributed_lumping) {
          auto fillFunctor = MatrixConstruction::PointwiseFillNoReuseFunctor<local_matrix_type, MueLu::MatrixConstruction::distributed_lumping>(lclA, results, lclGraph, lclFilteredA, filteringDirichletThreshold);
          Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_lumped_noreuse", range, fillFunctor);
        }
      }
    } else {
      if (reuseGraph) {
        auto fillFunctor = MatrixConstruction::PointwiseFillReuseFunctor<local_matrix_type, local_graph_type, false>(lclA, results, lclFilteredA, lclGraph, filteringDirichletThreshold);
        Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_unlumped_reuse", range, fillFunctor);
      } else {
        if (needToBuildFilteredA) {
          auto fillFunctor = MatrixConstruction::PointwiseFillNoReuseFunctor<local_matrix_type, MueLu::MatrixConstruction::no_lumping>(lclA, results, lclGraph, lclFilteredA, filteringDirichletThreshold);
          Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_unlumped_noreuse", range, fillFunctor);
        } else {
          auto fillFunctor = MatrixConstruction::PointwiseFillNoReuseFunctor<local_matrix_type, MueLu::MatrixConstruction::no_lumping, /*constructFilteredA=*/false>(lclA, results, lclGraph, lclFilteredA, filteringDirichletThreshold);
          Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_unlumped_noreuse", range, fillFunctor);
        }
      }
    }

    if (needToBuildFilteredA) {
      if (!reuseGraph)
        filteredA = MatrixFactory::Build(lclFilteredA, A->getRowMap(), A->getColMap(), A->getDomainMap(), A->getRangeMap());
      filteredA->SetFixedBlockSize(A->GetFixedBlockSize());

      if (reuseEigenvalue) {
        // Reuse max eigenvalue from A
        // It is unclear what eigenvalue is the best for the smoothing, but we already may have
        // the D^{-1}A estimate in A, may as well use it.
        // NOTE: ML does that too
        filteredA->SetMaxEigenvalueEstimate(A->GetMaxEigenvalueEstimate());
      } else {
        filteredA->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<SC>::one());
      }
    }

    graph = rcp(new LWGraph_kokkos(lclGraph, A->getRowMap(), A->getColMap(), "amalgamated graph of A"));
    graph->SetBoundaryNodeMap(boundaryNodes);
  } else {
    if (needToBuildFilteredA)
      filteredA = A;
    graph = rcp(new LWGraph_kokkos(A->getCrsGraph()->getLocalGraphDevice(), A->getRowMap(), A->getColMap(), "amalgamated graph of A"));
    graph->SetBoundaryNodeMap(boundaryNodes);
  }

  // Construct a second graph for coloring
  if (generateColoringGraph) {
    SubFactoryMonitor mColoringGraph(*this, "Construct coloring graph", currentLevel);

    filtered_rowptr = rowptr_type("rowptr_coloring_graph", lclA.numRows() + 1);
    if (localizeColoringGraph) {
      auto drop_offrank = Misc::DropOffRankFunctor(lclA, results);
      ScalarDroppingBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::template runDroppingFunctors<>(*A, results, filtered_rowptr, nnz_filtered, useBlocking, currentLevel, *this, drop_offrank);
    }
    if (symmetrizeColoringGraph) {
      auto symmetrize = Misc::SymmetrizeFunctor(lclA, results);
      ScalarDroppingBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::template runDroppingFunctors<>(*A, results, filtered_rowptr, nnz_filtered, useBlocking, currentLevel, *this, symmetrize);
    }
    auto colidx            = entries_type("entries_coloring_graph", nnz_filtered);
    auto lclGraph          = local_graph_type(colidx, filtered_rowptr);
    auto graphConstruction = MatrixConstruction::GraphConstruction<local_matrix_type, local_graph_type>(lclA, results, lclGraph);
    Kokkos::parallel_for("MueLu::CoalesceDrop::Construct_coloring_graph", range, graphConstruction);

    auto colorGraph = rcp(new LWGraph_kokkos(lclGraph, A->getRowMap(), A->getColMap(), "coloring graph"));
    Set(currentLevel, "Coloring Graph", colorGraph);
  }

  if (needToBuildFilteredA && pL.get<bool>("filtered matrix: count negative diagonals")) {
    // Count the negative diagonals (and display that information)
    GlobalOrdinal neg_count = MueLu::Utilities<SC, LO, GO, NO>::CountNegativeDiagonalEntries(*filteredA);
    GetOStream(Runtime0) << "CoalesceDrop: Negative diagonals: " << neg_count << std::endl;
  }

  LO dofsPerNode = 1;
  Set(currentLevel, "DofsPerNode", dofsPerNode);
  Set(currentLevel, "Graph", graph);
  if (needToBuildFilteredA)
    Set(currentLevel, "A", filteredA);

  return std::make_tuple(numDropped, (GlobalOrdinal)nnz_filtered, boundaryNodes);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::tuple<GlobalOrdinal, GlobalOrdinal, typename MueLu::LWGraph_kokkos<LocalOrdinal, GlobalOrdinal, Node>::boundary_nodes_type> CoalesceDropFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BuildVector(Level& currentLevel) const {
  FactoryMonitor m(*this, "BuildVector", currentLevel);

  using MatrixType        = Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using GraphType         = Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
  using local_matrix_type = typename MatrixType::local_matrix_device_type;
  using local_graph_type  = typename GraphType::local_graph_device_type;
  using rowptr_type       = typename local_graph_type::row_map_type::non_const_type;
  using entries_type      = typename local_graph_type::entries_type::non_const_type;
  using values_type       = typename local_matrix_type::values_type::non_const_type;
  using device_type       = typename Node::device_type;
  using memory_space      = typename device_type::memory_space;
  using results_view_type = Kokkos::View<DecisionType*, memory_space>;
  using magnitudeType     = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  using doubleMultiVector = Xpetra::MultiVector<magnitudeType, LO, GO, NO>;

  typedef Teuchos::ScalarTraits<Scalar> STS;
  const magnitudeType zero = Teuchos::ScalarTraits<magnitudeType>::zero();

  auto A = Get<RCP<Matrix>>(currentLevel, "A");

  const bool needToBuildFilteredA = currentLevel.IsRequested("A", this);

  /* NOTE: storageblocksize (from GetStorageBlockSize()) is the size of a block in the chosen storage scheme.
     blkSize is the number of storage blocks that must kept together during the amalgamation process.

     Both of these quantities may be different than numPDEs (from GetFixedBlockSize()), but the following must always hold:

     numPDEs = blkSize * storageblocksize.

     If numPDEs==1
       Matrix is point storage (classical CRS storage).  storageblocksize=1 and  blkSize=1
       No other values makes sense.

     If numPDEs>1
       If matrix uses point storage, then storageblocksize=1  and blkSize=numPDEs.
       If matrix uses block storage, with block size of n, then storageblocksize=n, and blkSize=numPDEs/n.
       Thus far, only storageblocksize=numPDEs and blkSize=1 has been tested.
    */

  TEUCHOS_TEST_FOR_EXCEPTION(A->GetFixedBlockSize() % A->GetStorageBlockSize() != 0, Exceptions::RuntimeError, "A->GetFixedBlockSize() needs to be a multiple of A->GetStorageBlockSize()");
  LO blkSize = A->GetFixedBlockSize() / A->GetStorageBlockSize();

  auto amalInfo = Get<RCP<AmalgamationInfo>>(currentLevel, "UnAmalgamationInfo");

  const RCP<const Map> rowMap = A->getRowMap();
  const RCP<const Map> colMap = A->getColMap();

  // build a node row map (uniqueMap = non-overlapping) and a node column map
  // (nonUniqueMap = overlapping). The arrays rowTranslation and colTranslation
  // stored in the AmalgamationInfo class container contain the local node id
  // given a local dof id. The data is calculated in the AmalgamationFactory and
  // stored in the variable "UnAmalgamationInfo" (which is of type AmalagamationInfo)
  const RCP<const Map> uniqueMap    = amalInfo->getNodeRowMap();
  const RCP<const Map> nonUniqueMap = amalInfo->getNodeColMap();
  Array<LO> rowTranslationArray     = *(amalInfo->getRowTranslation());  // TAW should be transform that into a View?
  Array<LO> colTranslationArray     = *(amalInfo->getColTranslation());

  Kokkos::View<LO*, Kokkos::MemoryUnmanaged>
      rowTranslationView(rowTranslationArray.getRawPtr(), rowTranslationArray.size());
  Kokkos::View<LO*, Kokkos::MemoryUnmanaged>
      colTranslationView(colTranslationArray.getRawPtr(), colTranslationArray.size());

  // get number of local nodes
  LO numNodes = Teuchos::as<LocalOrdinal>(uniqueMap->getLocalNumElements());
  typedef typename Kokkos::View<LocalOrdinal*, typename Node::device_type> id_translation_type;
  id_translation_type rowTranslation("dofId2nodeId", rowTranslationArray.size());
  id_translation_type colTranslation("ov_dofId2nodeId", colTranslationArray.size());
  Kokkos::deep_copy(rowTranslation, rowTranslationView);
  Kokkos::deep_copy(colTranslation, colTranslationView);

  // extract striding information
  blkSize                  = A->GetFixedBlockSize();  //< the full block size (number of dofs per node in strided map)
  LocalOrdinal blkId       = -1;                      //< the block id within a strided map or -1 if it is a full block map
  LocalOrdinal blkPartSize = A->GetFixedBlockSize();  //< stores block size of part blkId (or the full block size)
  if (A->IsView("stridedMaps") == true) {
    const RCP<const Map> myMap         = A->getRowMap("stridedMaps");
    const RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(myMap);
    TEUCHOS_TEST_FOR_EXCEPTION(strMap.is_null() == true, Exceptions::RuntimeError, "Map is not of type stridedMap");
    blkSize = Teuchos::as<const LocalOrdinal>(strMap->getFixedBlockSize());
    blkId   = strMap->getStridedBlockId();
    if (blkId > -1)
      blkPartSize = Teuchos::as<LocalOrdinal>(strMap->getStridingData()[blkId]);
  }

  TEUCHOS_TEST_FOR_EXCEPTION(A->getRowMap()->getLocalNumElements() % blkPartSize != 0, MueLu::Exceptions::RuntimeError, "MueLu::CoalesceDropFactory: Number of local elements is " << A->getRowMap()->getLocalNumElements() << " but should be a multiple of " << blkPartSize);

  //////////////////////////////////////////////////////////////////////
  // Process parameterlist
  const ParameterList& pL = GetParameterList();

  // Boundary detection
  const magnitudeType dirichletThreshold       = STS::magnitude(as<SC>(pL.get<double>("aggregation: Dirichlet threshold")));
  const magnitudeType rowSumTol                = as<magnitudeType>(pL.get<double>("aggregation: row sum drop tol"));
  const LocalOrdinal dirichletNonzeroThreshold = 1;
  const bool useGreedyDirichlet                = pL.get<bool>("aggregation: greedy Dirichlet");
  TEUCHOS_TEST_FOR_EXCEPTION(rowSumTol > zero, MueLu::Exceptions::RuntimeError, "MueLu::CoalesceDropFactory: RowSum is not implemented for vectorial problems.");

  // Dropping
  bool useBlocking                    = pL.get<bool>("aggregation: use blocking");
  std::string droppingMethod          = pL.get<std::string>("aggregation: drop scheme");
  std::string socUsesMatrix           = pL.get<std::string>("aggregation: strength-of-connection: matrix");
  std::string socUsesMeasure          = pL.get<std::string>("aggregation: strength-of-connection: measure");
  std::string distanceLaplacianMetric = pL.get<std::string>("aggregation: distance laplacian metric");
  bool symmetrizeDroppedGraph         = pL.get<bool>("aggregation: symmetrize graph after dropping");
  magnitudeType threshold;
  // If we're doing the ML-style halving of the drop tol at each level, we do that here.
  if (pL.get<bool>("aggregation: use ml scaling of drop tol"))
    threshold = pL.get<double>("aggregation: drop tol") / pow(2.0, currentLevel.GetLevelID());
  else
    threshold = as<magnitudeType>(pL.get<double>("aggregation: drop tol"));
  bool aggregationMayCreateDirichlet = pL.get<bool>("aggregation: dropping may create Dirichlet");

  // Fill
  const bool reuseGraph      = pL.get<bool>("filtered matrix: reuse graph") && needToBuildFilteredA;
  const bool reuseEigenvalue = pL.get<bool>("filtered matrix: reuse eigenvalue");

  const bool useRootStencil                            = pL.get<bool>("filtered matrix: use root stencil");
  const bool useSpreadLumping                          = pL.get<bool>("filtered matrix: use spread lumping");
  const std::string lumpingChoiceString                = pL.get<std::string>("filtered matrix: lumping choice");
  MueLu::MatrixConstruction::lumpingType lumpingChoice = MueLu::MatrixConstruction::no_lumping;
  if (needToBuildFilteredA) {
    if (lumpingChoiceString == "diag lumping")
      lumpingChoice = MueLu::MatrixConstruction::diag_lumping;
    else if (lumpingChoiceString == "distributed lumping")
      lumpingChoice = MueLu::MatrixConstruction::distributed_lumping;
  }

  const magnitudeType filteringDirichletThreshold = as<magnitudeType>(pL.get<double>("filtered matrix: Dirichlet threshold"));

  // coloring graph
  bool generateColoringGraph         = pL.get<bool>("aggregation: coloring: use color graph");
  const bool localizeColoringGraph   = pL.get<bool>("aggregation: coloring: localize color graph");
  const bool symmetrizeColoringGraph = true;

#ifdef HAVE_MUELU_COALESCEDROP_ALLOW_OLD_PARAMETERS
  translateOldAlgoParam(pL, droppingMethod, useBlocking, socUsesMatrix, socUsesMeasure, symmetrizeDroppedGraph, generateColoringGraph, threshold, lumpingChoice);
#endif
  {
    std::stringstream ss;
    ss << "dropping scheme = \"" << droppingMethod << "\", strength-of-connection measure = \"" << socUsesMeasure << "\", strength-of-connection matrix = \"" << socUsesMatrix << "\", ";
    if (socUsesMatrix == "distance laplacian")
      ss << "distance laplacian metric = \"" << distanceLaplacianMetric << "\", ";
    ss << "threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << ", useBlocking = " << useBlocking;
    ss << ", symmetrizeDroppedGraph = " << symmetrizeDroppedGraph << std::endl;

    GetOStream(Runtime0) << ss.str();
  }

  TEUCHOS_ASSERT(!useRootStencil);
  TEUCHOS_ASSERT(!useSpreadLumping);
  TEUCHOS_ASSERT((lumpingChoice != MueLu::MatrixConstruction::distributed_lumping) || !reuseGraph);
  if (droppingMethod == "cut-drop")
    TEUCHOS_TEST_FOR_EXCEPTION(threshold > 1.0, Exceptions::RuntimeError, "For cut-drop algorithms, \"aggregation: drop tol\" = " << threshold << ", needs to be <= 1.0");

  //////////////////////////////////////////////////////////////////////
  // We perform four sweeps over the rows of A:
  // Pass 1: detection of boundary nodes
  // Pass 2: diagonal extraction
  // Pass 3: drop decision for each entry and construction of the rowptr of the filtered matrix
  // Pass 4: fill of the filtered matrix
  //
  // Pass 1 and 3 apply a sequence of criteria to each row of the matrix.

  // TODO: We could merge pass 1 and 2.

  auto crsA  = toCrsMatrix(A);
  auto lclA  = crsA->getLocalMatrixDevice();
  auto range = range_type(0, numNodes);

  //////////////////////////////////////////////////////////////////////
  // Pass 1: Detect boundary nodes
  //
  // The following criteria are available:
  // - BoundaryDetection::VectorDirichletFunctor
  //   Marks rows as Dirichlet based on value threshold and number of off-diagonal entries

  // Dirichlet nodes
  auto boundaryNodes = boundary_nodes_type("boundaryNodes", numNodes);  // initialized to false
  {
    SubFactoryMonitor mBoundary(*this, "Boundary detection", currentLevel);

    if (useGreedyDirichlet) {
      auto dirichlet_detection = BoundaryDetection::VectorDirichletFunctor<local_matrix_type, true>(lclA, blkPartSize, boundaryNodes, dirichletThreshold, dirichletNonzeroThreshold);
      runBoundaryFunctors(lclA, boundaryNodes, dirichlet_detection);
    } else {
      auto dirichlet_detection = BoundaryDetection::VectorDirichletFunctor<local_matrix_type, false>(lclA, blkPartSize, boundaryNodes, dirichletThreshold, dirichletNonzeroThreshold);
      runBoundaryFunctors(lclA, boundaryNodes, dirichlet_detection);
    }
  }
  // In what follows, boundaryNodes can still still get modified if aggregationMayCreateDirichlet == true.
  // Otherwise we're now done with it now.

  //////////////////////////////////////////////////////////////////////
  // Pass 2 & 3: Diagonal extraction and determine dropping and construct
  //             rowptr of filtered matrix
  //
  // The following criteria are available:
  // - Misc::VectorDropBoundaryFunctor
  //   Drop all rows that have been marked as Dirichlet
  // - Misc::DropOffRankFunctor
  //   Drop all entries that are off-rank
  // - ClassicalDropping::DropFunctor
  //   Classical dropping
  // - DistanceLaplacian::VectorDropFunctor
  //   Distance Laplacian dropping
  // - Misc::KeepDiagonalFunctor
  //   Mark diagonal as KEEP
  // - Misc::MarkSingletonFunctor
  //   Mark singletons after dropping as Dirichlet

  // rowptr of filtered A
  auto filtered_rowptr = rowptr_type("rowptr", lclA.numRows() + 1);
  auto graph_rowptr    = rowptr_type("rowptr", numNodes + 1);
  // Number of nonzeros of filtered A and graph
  Kokkos::pair<LocalOrdinal, LocalOrdinal> nnz = {0, 0};

  // dropping decisions for each entry
  auto results = results_view_type("results", lclA.nnz());  // initialized to UNDECIDED

  RCP<Matrix> mergedA;
  {
    SubFactoryMonitor mDropping(*this, "Dropping decisions", currentLevel);

    {
      // Construct merged A.

      auto merged_rowptr      = rowptr_type("rowptr", numNodes + 1);
      LocalOrdinal nnz_merged = 0;

      auto functor = MatrixConstruction::MergeCountFunctor(lclA, blkPartSize, colTranslation, merged_rowptr);
      Kokkos::parallel_scan("MergeCount", range, functor, nnz_merged);

      local_graph_type lclMergedGraph;
      auto colidx_merged = entries_type("entries", nnz_merged);
      auto values_merged = values_type("values", nnz_merged);

      local_matrix_type lclMergedA = local_matrix_type("mergedA",
                                                       numNodes, nonUniqueMap->getLocalNumElements(),
                                                       nnz_merged,
                                                       values_merged, merged_rowptr, colidx_merged);

      auto fillFunctor = MatrixConstruction::MergeFillFunctor<local_matrix_type>(lclA, blkSize, colTranslation, lclMergedA);
      Kokkos::parallel_for("MueLu::CoalesceDrop::MergeFill", range, fillFunctor);

      mergedA = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(lclMergedA, uniqueMap, nonUniqueMap, uniqueMap, uniqueMap);
    }

    if (threshold != zero) {
      if (socUsesMatrix == "A") {
        if (socUsesMeasure == "unscaled") {
          VectorDroppingClassical<Scalar, LocalOrdinal, GlobalOrdinal, Node, Misc::UnscaledMeasure>::runDroppingFunctors_on_A(*A, *mergedA, blkPartSize, rowTranslation, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, currentLevel, *this);
        } else if (socUsesMeasure == "smoothed aggregation") {
          VectorDroppingClassical<Scalar, LocalOrdinal, GlobalOrdinal, Node, Misc::SmoothedAggregationMeasure>::runDroppingFunctors_on_A(*A, *mergedA, blkPartSize, rowTranslation, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, currentLevel, *this);
        } else if (socUsesMeasure == "signed ruge-stueben") {
          VectorDroppingClassical<Scalar, LocalOrdinal, GlobalOrdinal, Node, Misc::SignedRugeStuebenMeasure>::runDroppingFunctors_on_A(*A, *mergedA, blkPartSize, rowTranslation, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, currentLevel, *this);
        } else if (socUsesMeasure == "signed smoothed aggregation") {
          VectorDroppingClassical<Scalar, LocalOrdinal, GlobalOrdinal, Node, Misc::SignedSmoothedAggregationMeasure>::runDroppingFunctors_on_A(*A, *mergedA, blkPartSize, rowTranslation, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, currentLevel, *this);
        }
      } else if (socUsesMatrix == "distance laplacian") {
        auto coords = Get<RCP<doubleMultiVector>>(currentLevel, "Coordinates");

        Array<double> dlap_weights         = pL.get<Array<double>>("aggregation: distance laplacian directional weights");
        LocalOrdinal interleaved_blocksize = as<LocalOrdinal>(pL.get<int>("aggregation: block diagonal: interleaved blocksize"));
        if (socUsesMeasure == "distance laplacian") {
          LO dim = (LO)coords->getNumVectors();
          // If anything isn't 1.0 we need to turn on the weighting
          bool non_unity = false;
          for (LO i = 0; !non_unity && i < (LO)dlap_weights.size(); i++) {
            if (dlap_weights[i] != 1.0) {
              non_unity = true;
            }
          }
          if (non_unity) {
            if ((LO)dlap_weights.size() == dim) {
              distanceLaplacianMetric = "weighted";
            } else if ((LO)dlap_weights.size() == interleaved_blocksize * dim)
              distanceLaplacianMetric = "block weighted";
            else {
              TEUCHOS_TEST_FOR_EXCEPTION(1, Exceptions::RuntimeError,
                                         "length of 'aggregation: distance laplacian directional weights' must equal the coordinate dimension OR the coordinate dimension times the blocksize");
            }
            if (GetVerbLevel() & Statistics1)
              GetOStream(Statistics1) << "Using distance laplacian weights: " << dlap_weights << std::endl;
          }
        }

        if (socUsesMeasure == "unscaled") {
          VectorDroppingDistanceLaplacian<Scalar, LocalOrdinal, GlobalOrdinal, Node, Misc::UnscaledMeasure>::runDroppingFunctors_on_dlap(*A, *mergedA, blkPartSize, rowTranslation, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, distanceLaplacianMetric, dlap_weights, interleaved_blocksize, currentLevel, *this);
        } else if (socUsesMeasure == "smoothed aggregation") {
          VectorDroppingDistanceLaplacian<Scalar, LocalOrdinal, GlobalOrdinal, Node, Misc::SmoothedAggregationMeasure>::runDroppingFunctors_on_dlap(*A, *mergedA, blkPartSize, rowTranslation, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, distanceLaplacianMetric, dlap_weights, interleaved_blocksize, currentLevel, *this);
        } else if (socUsesMeasure == "signed ruge-stueben") {
          VectorDroppingDistanceLaplacian<Scalar, LocalOrdinal, GlobalOrdinal, Node, Misc::SignedRugeStuebenMeasure>::runDroppingFunctors_on_dlap(*A, *mergedA, blkPartSize, rowTranslation, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, distanceLaplacianMetric, dlap_weights, interleaved_blocksize, currentLevel, *this);
        } else if (socUsesMeasure == "signed smoothed aggregation") {
          VectorDroppingDistanceLaplacian<Scalar, LocalOrdinal, GlobalOrdinal, Node, Misc::SignedSmoothedAggregationMeasure>::runDroppingFunctors_on_dlap(*A, *mergedA, blkPartSize, rowTranslation, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, boundaryNodes, droppingMethod, threshold, aggregationMayCreateDirichlet, symmetrizeDroppedGraph, useBlocking, distanceLaplacianMetric, dlap_weights, interleaved_blocksize, currentLevel, *this);
        }
      }
    } else {
      Kokkos::deep_copy(results, KEEP);

      auto no_op = Misc::NoOpFunctor<LocalOrdinal>();
      VectorDroppingBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::template runDroppingFunctors<>(*A, *mergedA, blkPartSize, rowTranslation, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, useBlocking, currentLevel, *this, no_op);
    }

    if (symmetrizeDroppedGraph) {
      auto symmetrize = Misc::SymmetrizeFunctor(lclA, results);
      VectorDroppingBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::template runDroppingFunctors<>(*A, *mergedA, blkPartSize, rowTranslation, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, useBlocking, currentLevel, *this, symmetrize);
    }
  }
  LocalOrdinal nnz_filtered = nnz.first;
  LocalOrdinal nnz_graph    = nnz.second;
  GO numTotal               = mergedA->getLocalNumEntries();
  GO numDropped             = numTotal - nnz_graph;
  GO numGlobalDropped;
  Teuchos::reduceAll(*A->getRowMap()->getComm(), Teuchos::REDUCE_SUM, 1, &numDropped, &numGlobalDropped);
  // We now know the number of entries of filtered A and have the final rowptr.

  //////////////////////////////////////////////////////////////////////
  // Pass 4: Create local matrix for filtered A
  //
  // Dropped entries are optionally lumped to the diagonal.

  RCP<LWGraph_kokkos> graph;
  RCP<Matrix> filteredA;
  if (numGlobalDropped > 0) {
    SubFactoryMonitor mFill(*this, "Filtered matrix fill", currentLevel);

    local_graph_type lclGraph;
    {
      auto colidx = entries_type("entries", nnz_graph);
      lclGraph    = local_graph_type(colidx, graph_rowptr);
    }

    local_matrix_type lclFilteredA;
    if (reuseGraph) {
      // needToBuildFilteredA is true
      lclFilteredA = local_matrix_type("filteredA", lclA.graph, lclA.numCols());
    } else if (needToBuildFilteredA) {
      auto colidx  = entries_type("entries", nnz_filtered);
      auto values  = values_type("values", nnz_filtered);
      lclFilteredA = local_matrix_type("filteredA",
                                       lclA.numRows(), lclA.numCols(),
                                       nnz_filtered,
                                       values, filtered_rowptr, colidx);
    }

    if (lumpingChoice != MueLu::MatrixConstruction::no_lumping) {
      if (reuseGraph) {
        auto fillFunctor = MatrixConstruction::VectorFillFunctor<local_matrix_type, true, true>(lclA, blkPartSize, colTranslation, results, lclFilteredA, lclGraph, filteringDirichletThreshold);
        Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_lumped_reuse", range, fillFunctor);
      } else {
        auto fillFunctor = MatrixConstruction::VectorFillFunctor<local_matrix_type, true, false>(lclA, blkPartSize, colTranslation, results, lclFilteredA, lclGraph, filteringDirichletThreshold);
        Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_lumped_noreuse", range, fillFunctor);
      }
    } else {
      if (reuseGraph) {
        auto fillFunctor = MatrixConstruction::VectorFillFunctor<local_matrix_type, false, true>(lclA, blkSize, colTranslation, results, lclFilteredA, lclGraph, filteringDirichletThreshold);
        Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_unlumped_reuse", range, fillFunctor);
      } else {
        if (needToBuildFilteredA) {
          auto fillFunctor = MatrixConstruction::VectorFillFunctor<local_matrix_type, false, false>(lclA, blkSize, colTranslation, results, lclFilteredA, lclGraph, filteringDirichletThreshold);
          Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_unlumped_noreuse", range, fillFunctor);
        } else {
          auto fillFunctor = MatrixConstruction::VectorFillFunctor<local_matrix_type, false, false, /*constructFilteredA=*/false>(lclA, blkSize, colTranslation, results, lclFilteredA, lclGraph, filteringDirichletThreshold);
          Kokkos::parallel_for("MueLu::CoalesceDrop::Fill_unlumped_noreuse", range, fillFunctor);
        }
      }
    }

    if (needToBuildFilteredA) {
      filteredA = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(lclFilteredA, A->getRowMap(), A->getColMap(), A->getDomainMap(), A->getRangeMap());
      filteredA->SetFixedBlockSize(blkSize);

      if (reuseEigenvalue) {
        // Reuse max eigenvalue from A
        // It is unclear what eigenvalue is the best for the smoothing, but we already may have
        // the D^{-1}A estimate in A, may as well use it.
        // NOTE: ML does that too
        filteredA->SetMaxEigenvalueEstimate(A->GetMaxEigenvalueEstimate());
      } else {
        filteredA->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<SC>::one());
      }
    }

    graph = rcp(new LWGraph_kokkos(lclGraph, uniqueMap, nonUniqueMap, "amalgamated graph of A"));
    graph->SetBoundaryNodeMap(boundaryNodes);
  } else {
    if (needToBuildFilteredA)
      filteredA = A;
    graph = rcp(new LWGraph_kokkos(mergedA->getCrsGraph()->getLocalGraphDevice(), uniqueMap, nonUniqueMap, "amalgamated graph of A"));
    graph->SetBoundaryNodeMap(boundaryNodes);
  }

  // Construct a second graph for coloring
  if (generateColoringGraph) {
    SubFactoryMonitor mColoringGraph(*this, "Construct coloring graph", currentLevel);

    filtered_rowptr = rowptr_type("rowptr_coloring_graph", lclA.numRows() + 1);
    graph_rowptr    = rowptr_type("rowptr", numNodes + 1);
    if (localizeColoringGraph) {
      auto drop_offrank = Misc::DropOffRankFunctor(lclA, results);
      VectorDroppingBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::template runDroppingFunctors<>(*A, *mergedA, blkPartSize, rowTranslation, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, useBlocking, currentLevel, *this, drop_offrank);
    }
    if (symmetrizeColoringGraph) {
      auto symmetrize = Misc::SymmetrizeFunctor(lclA, results);
      VectorDroppingBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::template runDroppingFunctors<>(*A, *mergedA, blkPartSize, rowTranslation, colTranslation, results, filtered_rowptr, graph_rowptr, nnz, useBlocking, currentLevel, *this, symmetrize);
    }
    auto colidx            = entries_type("entries_coloring_graph", nnz_filtered);
    auto lclGraph          = local_graph_type(colidx, filtered_rowptr);
    auto graphConstruction = MatrixConstruction::GraphConstruction<local_matrix_type, local_graph_type>(lclA, results, lclGraph);
    Kokkos::parallel_for("MueLu::CoalesceDrop::Construct_coloring_graph", range, graphConstruction);

    auto colorGraph = rcp(new LWGraph_kokkos(lclGraph, A->getRowMap(), A->getColMap(), "coloring graph"));
    Set(currentLevel, "Coloring Graph", colorGraph);
  }

  LO dofsPerNode = blkSize;

  Set(currentLevel, "DofsPerNode", dofsPerNode);
  Set(currentLevel, "Graph", graph);
  if (needToBuildFilteredA)
    Set(currentLevel, "A", filteredA);

  return std::make_tuple(numDropped, (GlobalOrdinal)nnz_graph, boundaryNodes);
}

}  // namespace MueLu
#endif  // MUELU_COALESCEDROPFACTORY_KOKKOS_DEF_HPP
