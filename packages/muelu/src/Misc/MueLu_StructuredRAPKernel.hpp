// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// @HEADER

#ifndef MUELU_STRUCTUREDRAPKERNEL_HPP
#define MUELU_STRUCTUREDRAPKERNEL_HPP

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>

#include <Teuchos_Array.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_VectorFactory.hpp>

namespace MueLu {
namespace Details {

/**
 * Specialized numeric kernel for two- and three-dimensional structured RAP.
 *
 * The implemented path requires two or three dimensions and order-zero interpolation.
 * Its numeric plan is generated from the detected fine stencil and sparse graphs
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class StructuredRAPKernel {
 public:
  using Matrix = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

 private:
  using LO          = LocalOrdinal;
  using device_type = typename Node::device_type;
  using lo_view     = Kokkos::View<LocalOrdinal*, device_type>;
  using size_view   = Kokkos::View<size_t*, device_type>;

  KOKKOS_INLINE_FUNCTION
  static LocalOrdinal coarseAnchor(const LocalOrdinal coarse,
                                   const LocalOrdinal numFine,
                                   const LocalOrdinal numCoarse,
                                   const LocalOrdinal rate) {
    return coarse + 1 == numCoarse ? numFine - 1 : coarse * rate;
  }

  // Inclusive fine-coordinate range belonging to an order-zero coarse node.
  // Keep the midpoint tie on the lower coarse node, including the short end interval.
  KOKKOS_INLINE_FUNCTION
  static Kokkos::Array<LO, 2> aggregateBounds(LO coarse, LO numFine,
                                            LO numCoarse, LO rate) {
    const LO anchor = coarseAnchor(coarse, numFine, numCoarse, rate);
    const LO previous = coarse == 0 ? LO(0) : coarseAnchor(coarse - 1, numFine, numCoarse, rate);
    const LO next = coarse + 1 == numCoarse ? numFine - 1
                                          : coarseAnchor(coarse + 1, numFine, numCoarse, rate);
    return {{coarse == 0 ? LO(0) : previous + (anchor - previous) / 2 + 1,
             coarse + 1 == numCoarse ? numFine - 1 : anchor + (next - anchor) / 2}};
  }

  static size_t boundaryClass(LO coarse, LO numCoarse,
                             size_t leftDepth, size_t rightDepth) {
    const size_t left = static_cast<size_t>(coarse);
    const size_t right = static_cast<size_t>(numCoarse - coarse - 1);
    return left < leftDepth ? left
                           : (right < rightDepth ? leftDepth + rightDepth - right : leftDepth);
  }

  // Keep host mirrors alive until the preparation fence, allowing asynchronous copies.
  template <class View>
  struct PlanArray {
    View device;
    typename View::host_mirror_type host;

    PlanArray(const char* label, size_t size)
        : device(Kokkos::ViewAllocateWithoutInitializing(label), size),
          host(Kokkos::create_mirror_view(device)) {}

    View copyToDevice(const typename device_type::execution_space& executionSpace) const {
      Kokkos::deep_copy(executionSpace, device, host);
      return device;
    }
  };

  struct AxisClassMap {
    lo_view boundaryClasses;
    LO numCoordinates = 0;
    LO leftCount = 0;
    LO rightCount = 0;
    LO interiorClass = 0;

    KOKKOS_INLINE_FUNCTION
    LO operator()(const LO coordinate) const {
      if (coordinate < leftCount)
        return boundaryClasses(coordinate);
      const LO rightBegin = numCoordinates - rightCount;
      if (coordinate >= rightBegin)
        return boundaryClasses(leftCount + coordinate - rightBegin);
      return interiorClass;
    }

    KOKKOS_INLINE_FUNCTION
    bool isInterior(const LO coordinate) const {
      return coordinate >= leftCount &&
             coordinate < numCoordinates - rightCount;
    }
  };

  KOKKOS_INLINE_FUNCTION
  static LocalOrdinal constantInterpolationCoarse(
      const LocalOrdinal fine,
      const LocalOrdinal numFine,
      const LocalOrdinal rate) {
    LocalOrdinal endRate = (numFine - 1) % rate;
    if (endRate == 0)
      endRate = rate;

    LocalOrdinal coarse = fine / rate;
    const LocalOrdinal remainder = fine % rate;
    const LocalOrdinal activeRate =
        fine < numFine - endRate ? rate : endRate;
    if (remainder > activeRate / 2)
      ++coarse;
    return coarse;
  }

  KOKKOS_INLINE_FUNCTION
  static LO coarseRowForFineRow(
      const LO fineRow, const LO dofsPerNode,
      const Kokkos::Array<LO, 3>& fineNodes,
      const Kokkos::Array<LO, 3>& coarseNodes,
      const Kokkos::Array<LO, 3>& rate) {
    const LO fineNode = fineRow / dofsPerNode;
    const LO fineX = fineNode % fineNodes[0];
    const LO fineYZ = fineNode / fineNodes[0];
    const LO coarseX = constantInterpolationCoarse(
        fineX, fineNodes[0], rate[0]);
    const LO coarseY = constantInterpolationCoarse(
        fineYZ % fineNodes[1], fineNodes[1], rate[1]);
    const LO coarseZ = constantInterpolationCoarse(
        fineYZ / fineNodes[1], fineNodes[2], rate[2]);
    return static_cast<LO>(
        (static_cast<size_t>(coarseX) +
         static_cast<size_t>(coarseNodes[0]) *
             (static_cast<size_t>(coarseY) +
              static_cast<size_t>(coarseNodes[1]) * coarseZ)) *
            dofsPerNode + fineRow % dofsPerNode);
  }

 public:
  template <class FineStencil>
  static void Compute(const Matrix& A,
                      const Matrix& P,
                      Matrix& Ac,
                      const FineStencil& fineStencil,
                      const int interpolationOrder,
                      const Teuchos::Array<LocalOrdinal>& lFineNodesPerDim,
                      const Teuchos::Array<LocalOrdinal>& lCoarseNodesPerDim,
                      const Teuchos::Array<int>& structuredCoarseningRate,
                      const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
    using SC              = Scalar;
    using execution_space = typename device_type::execution_space;
    using range_policy    = Kokkos::RangePolicy<execution_space,
                                                Kokkos::IndexType<size_t>>;
    const std::string timerPrefix = "MueLu: StructuredRAPKernel: ";

    const int numDimensions = fineStencil.numDimensions;
    const LO dofsPerNode     = fineStencil.dofsPerNode;

    Kokkos::Array<LO, 3> localFineNodes{{Teuchos::as<LO>(1), Teuchos::as<LO>(1), Teuchos::as<LO>(1)}};
    Kokkos::Array<LO, 3> localCoarseNodes{{Teuchos::as<LO>(1), Teuchos::as<LO>(1), Teuchos::as<LO>(1)}};
    Kokkos::Array<LO, 3> coarseningRate{{Teuchos::as<LO>(1), Teuchos::as<LO>(1), Teuchos::as<LO>(1)}};
    for (int dim = 0; dim < numDimensions; ++dim) {
      localFineNodes[dim] = lFineNodesPerDim[dim];
      localCoarseNodes[dim] = lCoarseNodesPerDim[dim];
      coarseningRate[dim] = static_cast<LO>(
          structuredCoarseningRate.size() == 1
              ? structuredCoarseningRate[0]
              : structuredCoarseningRate[dim]);
    }

    execution_space executionSpace;
    bool structuredDimensionsMatch = true;
    for (int dim = 0; dim < numDimensions && structuredDimensionsMatch; ++dim) {
      const LO expectedCoarseNodes =
          localFineNodes[dim] == 1
              ? LO(1)
              : static_cast<LO>((localFineNodes[dim] - 2) /
                                    coarseningRate[dim] +
                                2);
      structuredDimensionsMatch =
          localCoarseNodes[dim] == expectedCoarseNodes;
    }

    const size_t expectedFineRows =
        static_cast<size_t>(localFineNodes[0]) *
        static_cast<size_t>(localFineNodes[1]) *
        static_cast<size_t>(localFineNodes[2]) *
        static_cast<size_t>(dofsPerNode);
    const size_t expectedCoarseRows =
        static_cast<size_t>(localCoarseNodes[0]) *
        static_cast<size_t>(localCoarseNodes[1]) *
        static_cast<size_t>(localCoarseNodes[2]) *
        static_cast<size_t>(dofsPerNode);
    const bool distributed = A.getRowMap()->getComm()->getSize() > 1;
    const bool useDetectedStencil =
        interpolationOrder == 0 &&
        (numDimensions >= 1 && numDimensions <= 3) &&
        structuredDimensionsMatch &&
        A.GetFixedBlockSize() == dofsPerNode &&
        A.getRowMap()->getLocalNumElements() == expectedFineRows &&
        Ac.getRowMap()->getLocalNumElements() == expectedCoarseRows;

    if (!useDetectedStencil)
      throw std::runtime_error(
          "StructuredRAPKernel supports only 1D, 2D, or 3D order-zero structured RAP.");

    const SC zero  = Teuchos::ScalarTraits<SC>::zero();
    const auto localA = A.getLocalMatrixDevice();
    auto localAc      = Ac.getLocalMatrixDevice();

    using GOVector = Xpetra::Vector<GlobalOrdinal, LO, GlobalOrdinal, Node>;
    Teuchos::RCP<const Xpetra::Map<LO, GlobalOrdinal, Node>> remoteFineMap;
    Teuchos::RCP<GOVector> remoteCoarseGids;
    if (distributed) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          !P.getRowMap()->isSameAs(*A.getRowMap()), std::runtime_error,
          "StructuredRAPKernel requires A and P to have identical row maps.");

      const auto graphImporter = A.getCrsGraph()->getImporter();
      TEUCHOS_TEST_FOR_EXCEPTION(
          graphImporter.is_null(), std::runtime_error,
          "StructuredRAPKernel requires an A graph importer in parallel.");
      const auto remoteLids = graphImporter->getRemoteLIDs();
      lo_view remoteLidsDevice(
          Kokkos::ViewAllocateWithoutInitializing(
              "StructuredRAP: remote fine LIDs"), remoteLids.size());
      auto remoteLidsMirror = Kokkos::create_mirror_view(remoteLidsDevice);
      for (size_t i = 0; i < static_cast<size_t>(remoteLids.size()); ++i)
        remoteLidsMirror(i) = remoteLids[i];
      Kokkos::deep_copy(executionSpace, remoteLidsDevice, remoteLidsMirror);
      Kokkos::View<GlobalOrdinal*, device_type> remoteFineGidsDevice(
          Kokkos::ViewAllocateWithoutInitializing(
              "StructuredRAP: remote fine GIDs"), remoteLids.size());
      const auto fineColumnMap = A.getColMap()->getLocalMap();
      Kokkos::parallel_for(
          "StructuredRAP: map remote fine GIDs",
          range_policy(executionSpace, 0, remoteLids.size()),
          KOKKOS_LAMBDA(const size_t i) {
            remoteFineGidsDevice(i) =
                fineColumnMap.getGlobalElement(remoteLidsDevice(i));
          });
      const auto remoteFineGidsHost = Kokkos::create_mirror_view_and_copy(
          Kokkos::HostSpace(), remoteFineGidsDevice);
      Teuchos::Array<GlobalOrdinal> remoteFineGids(remoteLids.size());
      for (size_t i = 0; i < static_cast<size_t>(remoteLids.size()); ++i)
        remoteFineGids[i] = remoteFineGidsHost(i);
      remoteFineMap = Xpetra::MapFactory<LO, GlobalOrdinal, Node>::Build(
          A.getRowMap()->lib(),
          Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
          remoteFineGids(), A.getRowMap()->getIndexBase(),
          A.getRowMap()->getComm());
      const auto remoteImporter =
          graphImporter->createRemoteOnlyImport(remoteFineMap);

      const auto localP = P.getLocalMatrixDevice();
      auto ownedCoarseGids =
          Xpetra::VectorFactory<GlobalOrdinal, LO, GlobalOrdinal, Node>::Build(
              A.getRowMap(), false);
      auto ownedData = ownedCoarseGids->getLocalViewDevice(
          Tpetra::Access::OverwriteAll);
      const auto coarseColumnMap = P.getColMap()->getLocalMap();
      const auto exportLidsHost = remoteImporter->getExportLIDs();
      PlanArray<lo_view> exportLids(
          "StructuredRAP: halo export LIDs", exportLidsHost.size());
      for (size_t i = 0; i < static_cast<size_t>(exportLidsHost.size()); ++i)
        exportLids.host(i) = exportLidsHost[i];
      const auto exportLidsDevice = exportLids.copyToDevice(executionSpace);
      size_t invalidInterpolationRows = 0;
      Kokkos::parallel_reduce(
          "StructuredRAP: build halo coarse target GIDs",
          range_policy(executionSpace, 0,
                       static_cast<size_t>(exportLidsHost.size())),
          KOKKOS_LAMBDA(const size_t exportIndex, size_t& invalid) {
            const LO row = exportLidsDevice(exportIndex);
            const size_t rowBegin = localP.graph.row_map(row);
            const size_t rowEnd   = localP.graph.row_map(row + 1);
            if (rowEnd - rowBegin != size_t(1)) {
              ++invalid;
              return;
            }
            ownedData(row, 0) = coarseColumnMap.getGlobalElement(
                localP.graph.entries(rowBegin));
          },
          invalidInterpolationRows);
      executionSpace.fence("StructuredRAP: coarse target GIDs ready");
      TEUCHOS_TEST_FOR_EXCEPTION(
          invalidInterpolationRows != 0, std::runtime_error,
          "StructuredRAPKernel requires exactly one interpolation entry per fine row.");

      remoteCoarseGids =
          Xpetra::VectorFactory<GlobalOrdinal, LO, GlobalOrdinal, Node>::Build(
              remoteFineMap, false);
      remoteCoarseGids->doImport(
          *ownedCoarseGids, *remoteImporter, Xpetra::INSERT);
    }

    const size_t numFineStencilOffsets = fineStencil.stencilOffsets.size();
    Kokkos::Array<int, 3> minStencilOffset{{
        fineStencil.stencilOffsets[0].x,
        fineStencil.stencilOffsets[0].y,
        fineStencil.stencilOffsets[0].z}};
    Kokkos::Array<int, 3> maxStencilOffset = minStencilOffset;
    for (size_t offset = 1; offset < numFineStencilOffsets; ++offset) {
      const int values[3] = {fineStencil.stencilOffsets[offset].x,
                             fineStencil.stencilOffsets[offset].y,
                             fineStencil.stencilOffsets[offset].z};
      for (int dim = 0; dim < 3; ++dim) {
        minStencilOffset[dim] = std::min(minStencilOffset[dim], values[dim]);
        maxStencilOffset[dim] = std::max(maxStencilOffset[dim], values[dim]);
      }
    }
    const size_t numCoarseNodes =
        static_cast<size_t>(localCoarseNodes[0]) *
        static_cast<size_t>(localCoarseNodes[1]) *
        static_cast<size_t>(localCoarseNodes[2]);
    const size_t numDofs = static_cast<size_t>(dofsPerNode);

    Kokkos::Array<AxisClassMap, 3> coarseAggregatePlanClass;
    size_view reverseContributionOffsets;
    size_view reverseFineNodeOffsets;
    size_view reverseARowEntryOffsets;
    size_view reverseDirectAValueOffsets;
    Kokkos::Array<size_t, 3> numAggregatePlanClasses;
    size_t maxAcEntriesPerRow = 0;
    {
      Teuchos::TimeMonitor timer(
          *Teuchos::TimeMonitor::getNewTimer(
              timerPrefix + "prepare numeric entry plan"));

      const auto aRowMapHost = Kokkos::create_mirror_view_and_copy(
          Kokkos::HostSpace(), localA.graph.row_map);
      const auto aEntriesHost = Kokkos::create_mirror_view_and_copy(
          Kokkos::HostSpace(), localA.graph.entries);
      const auto acRowMapHost = Kokkos::create_mirror_view_and_copy(
          Kokkos::HostSpace(), localAc.graph.row_map);
      const auto acEntriesHost = Kokkos::create_mirror_view_and_copy(
          Kokkos::HostSpace(), localAc.graph.entries);

      const bool fineStencilHasCenter = std::any_of(
          fineStencil.stencilOffsets.begin(), fineStencil.stencilOffsets.end(),
          [](const auto& offset) {
            return offset.x == 0 && offset.y == 0 && offset.z == 0;
          });
      maxAcEntriesPerRow =
          (numFineStencilOffsets + (fineStencilHasCenter ? size_t(0) : size_t(1))) * numDofs;

      const size_t leftBoundaryDepth[3] = {
          static_cast<size_t>(std::max(0, -minStencilOffset[0])),
          static_cast<size_t>(std::max(0, -minStencilOffset[1])),
          static_cast<size_t>(std::max(0, -minStencilOffset[2]))};
      const size_t rightBoundaryDepth[3] = {
          static_cast<size_t>(std::max(0, maxStencilOffset[0])),
          static_cast<size_t>(std::max(0, maxStencilOffset[1])),
          static_cast<size_t>(std::max(0, maxStencilOffset[2]))};

      struct AxisAggregatePlanClass {
        LO representativeCoarse;
        LO fineBegin;
        LO fineEnd;
      };
      auto buildAxisAggregatePlanClasses = [&](
          const LO numFine, const LO numCoarse, const LO rate,
          const int minOffset, const int maxOffset,
          const size_t leftDepth, const size_t rightDepth,
          std::vector<AxisAggregatePlanClass>& classes) {
        std::vector<std::vector<int>> aggregateSignatures;
        std::vector<int> aggregateSignature;
        auto aggregateClass = [&](const LO coarse) {
          const auto bounds = aggregateBounds(coarse, numFine, numCoarse, rate);
          const LO fineBegin = bounds[0], fineEnd = bounds[1];
          aggregateSignature.clear();
          aggregateSignature.push_back(static_cast<int>(
              boundaryClass(coarse, numCoarse, leftDepth, rightDepth)));
          for (LO fine = fineBegin; fine <= fineEnd; ++fine) {
            const LO fineCoarse = constantInterpolationCoarse(
                fine, numFine, rate);
            aggregateSignature.push_back(static_cast<int>(boundaryClass(
                fineCoarse, numCoarse, leftDepth, rightDepth)));
            for (int offset = minOffset; offset <= maxOffset; ++offset) {
              const LO neighbor = fine + static_cast<LO>(offset);
              const bool valid = neighbor >= LO(0) && neighbor < numFine;
              aggregateSignature.push_back(valid);
              aggregateSignature.push_back(valid ? static_cast<int>(
                  constantInterpolationCoarse(neighbor, numFine, rate) -
                  fineCoarse) : 0);
            }
          }
          const auto found = std::find(
              aggregateSignatures.begin(), aggregateSignatures.end(),
              aggregateSignature);
          if (found != aggregateSignatures.end())
            return static_cast<LO>(found - aggregateSignatures.begin());
          aggregateSignatures.push_back(aggregateSignature);
          classes.push_back({coarse, fineBegin, fineEnd});
          return static_cast<LO>(classes.size() - 1);
        };

        AxisClassMap result;
        result.numCoordinates = numCoarse;
        const LO explicitLimit = static_cast<LO>(
            2 * (leftDepth + rightDepth + size_t(3)));
        std::vector<LO> leftClasses, rightClasses;
        if (numCoarse <= explicitLimit) {
          result.leftCount = numCoarse;
          leftClasses.reserve(static_cast<size_t>(numCoarse));
          for (LO coarse = 0; coarse < numCoarse; ++coarse)
            leftClasses.push_back(aggregateClass(coarse));
        } else {
          const LO representative = numCoarse / LO(2);
          result.interiorClass = aggregateClass(representative);
          for (LO coarse = 0; coarse < representative; ++coarse) {
            const LO id = aggregateClass(coarse);
            if (id == result.interiorClass)
              break;
            leftClasses.push_back(id);
          }
          for (LO coarse = numCoarse - 1; coarse > representative; --coarse) {
            const LO id = aggregateClass(coarse);
            if (id == result.interiorClass)
              break;
            rightClasses.push_back(id);
          }
          std::reverse(rightClasses.begin(), rightClasses.end());
          result.leftCount = static_cast<LO>(leftClasses.size());
          result.rightCount = static_cast<LO>(rightClasses.size());
        }

        result.boundaryClasses = lo_view(
            Kokkos::ViewAllocateWithoutInitializing(
                "StructuredRAP: compact axis aggregate classes"),
            leftClasses.size() + rightClasses.size());
        auto host = Kokkos::create_mirror_view(result.boundaryClasses);
        size_t entry = 0;
        for (const LO id : leftClasses)
          host(entry++) = id;
        for (const LO id : rightClasses)
          host(entry++) = id;
        Kokkos::deep_copy(result.boundaryClasses, host);
        return result;
      };

      std::vector<AxisAggregatePlanClass> aggregatePlanClasses[3];
      for (int dim = 0; dim < 3; ++dim) {
        coarseAggregatePlanClass[dim] = buildAxisAggregatePlanClasses(
            localFineNodes[dim], localCoarseNodes[dim], coarseningRate[dim],
            minStencilOffset[dim], maxStencilOffset[dim],
            leftBoundaryDepth[dim], rightBoundaryDepth[dim],
            aggregatePlanClasses[dim]);
        numAggregatePlanClasses[dim] = aggregatePlanClasses[dim].size();
      }

      const size_t reversePlanRowCount =
          numAggregatePlanClasses[0] * numAggregatePlanClasses[1] *
          numAggregatePlanClasses[2] * numDofs;
      const size_t reverseBucketCount =
          reversePlanRowCount * maxAcEntriesPerRow;
      struct PlanContribution {
        size_t fineNodeOffset, aRowEntryOffset, directAValueOffset;
      };
      std::vector<std::vector<PlanContribution>> reverseBuckets(reverseBucketCount);
      struct PendingContribution {
        size_t reverseRow, fineNodeOffset, aRowEntryOffset, directAValueOffset;
        size_t acRowBegin, acRowEnd;
        LO fineColumn;
      };
      std::vector<PendingContribution> pendingContributions;
      for (size_t zClass = 0; zClass < numAggregatePlanClasses[2]; ++zClass) {
        const auto& zAggregate = aggregatePlanClasses[2][zClass];
        const LO coarseZ = zAggregate.representativeCoarse;
        const size_t fineZBegin = static_cast<size_t>(zAggregate.fineBegin);
        const size_t fineZEnd = static_cast<size_t>(zAggregate.fineEnd);
        for (size_t yClass = 0; yClass < numAggregatePlanClasses[1]; ++yClass) {
          const auto& yAggregate = aggregatePlanClasses[1][yClass];
          const LO coarseY = yAggregate.representativeCoarse;
          const size_t fineYBegin = static_cast<size_t>(yAggregate.fineBegin);
          const size_t fineYEnd = static_cast<size_t>(yAggregate.fineEnd);
          for (size_t xClass = 0; xClass < numAggregatePlanClasses[0]; ++xClass) {
            const auto& xAggregate = aggregatePlanClasses[0][xClass];
            const LO coarseX = xAggregate.representativeCoarse;
            const size_t fineXBegin = static_cast<size_t>(xAggregate.fineBegin);
            const size_t fineXEnd = static_cast<size_t>(xAggregate.fineEnd);
            const size_t coarseNode =
                static_cast<size_t>(coarseX) +
                static_cast<size_t>(localCoarseNodes[0]) *
                    (static_cast<size_t>(coarseY) +
                     static_cast<size_t>(localCoarseNodes[1]) * coarseZ);
            const size_t firstFineNode =
                fineXBegin + static_cast<size_t>(localFineNodes[0]) *
                    (fineYBegin + static_cast<size_t>(localFineNodes[1]) * fineZBegin);
            for (size_t rowDof = 0; rowDof < numDofs; ++rowDof) {
              const size_t reverseRow =
                  ((zClass * numAggregatePlanClasses[1] + yClass) *
                       numAggregatePlanClasses[0] + xClass) * numDofs + rowDof;
              const size_t acRow = coarseNode * numDofs + rowDof;
              const size_t acRowBegin = acRowMapHost(acRow);
              const size_t acRowEnd = acRowMapHost(acRow + 1);
              const size_t firstAValue =
                  aRowMapHost(firstFineNode * numDofs + rowDof);
              for (size_t fineZ = fineZBegin; fineZ <= fineZEnd; ++fineZ) {
                for (size_t fineY = fineYBegin; fineY <= fineYEnd; ++fineY) {
                  for (size_t fineX = fineXBegin; fineX <= fineXEnd; ++fineX) {
                    const size_t sourceFineNode =
                        fineX + static_cast<size_t>(localFineNodes[0]) *
                            (fineY + static_cast<size_t>(localFineNodes[1]) * fineZ);
                    const size_t aRowBegin =
                        aRowMapHost(sourceFineNode * numDofs + rowDof);
                    const size_t aRowEnd =
                        aRowMapHost(sourceFineNode * numDofs + rowDof + 1);
                    for (size_t aEntry = aRowBegin; aEntry < aRowEnd; ++aEntry)
                      pendingContributions.push_back(PendingContribution{
                          reverseRow, sourceFineNode - firstFineNode,
                          aEntry - aRowBegin, aEntry - firstAValue,
                          acRowBegin, acRowEnd,
                          aEntriesHost(aEntry)});
                  }
                }
              }
            }
          }
        }
      }

      std::vector<LO> planFineColumns;
      planFineColumns.reserve(pendingContributions.size());
      for (const auto& contribution : pendingContributions)
        planFineColumns.push_back(contribution.fineColumn);
      std::sort(planFineColumns.begin(), planFineColumns.end());
      planFineColumns.erase(
          std::unique(planFineColumns.begin(), planFineColumns.end()),
          planFineColumns.end());

      PlanArray<lo_view> fineColumns(
          "StructuredRAP: plan fine columns", planFineColumns.size());
      for (size_t column = 0; column < planFineColumns.size(); ++column)
        fineColumns.host(column) = planFineColumns[column];
      const auto fineColumnsDevice = fineColumns.copyToDevice(executionSpace);
      lo_view targetColumnsDevice(
          Kokkos::ViewAllocateWithoutInitializing(
              "StructuredRAP: plan target columns"),
          planFineColumns.size());
      const auto fineColumnMap = A.getColMap()->getLocalMap();
      const auto fineRowMap = A.getRowMap()->getLocalMap();
      const auto coarseRowMap = Ac.getRowMap()->getLocalMap();
      const auto coarseColumnMap = Ac.getColMap()->getLocalMap();
      const LO invalid = Teuchos::OrdinalTraits<LO>::invalid();
      size_t invalidTargets = 0;
      if (distributed) {
        const auto remoteFineLocalMap = remoteFineMap->getLocalMap();
        const auto remoteCoarseData = remoteCoarseGids->getLocalViewDevice(
            Tpetra::Access::ReadOnly);
        Kokkos::parallel_reduce(
            "StructuredRAP: map plan columns", range_policy(
                executionSpace, 0, planFineColumns.size()),
            KOKKOS_LAMBDA(const size_t column, size_t& invalidCount) {
              const GlobalOrdinal fineGid = fineColumnMap.getGlobalElement(
                  fineColumnsDevice(column));
              const LO ownedFine = fineRowMap.getLocalElement(fineGid);
              GlobalOrdinal coarseGid;
              if (ownedFine != invalid) {
                coarseGid = coarseRowMap.getGlobalElement(
                    coarseRowForFineRow(ownedFine, dofsPerNode,
                                        localFineNodes, localCoarseNodes,
                                        coarseningRate));
              } else {
                const LO remoteFine = remoteFineLocalMap.getLocalElement(fineGid);
                if (remoteFine == invalid) {
                  ++invalidCount;
                  return;
                }
                coarseGid = remoteCoarseData(remoteFine, 0);
              }
              const LO target = coarseColumnMap.getLocalElement(coarseGid);
              if (target == invalid)
                ++invalidCount;
              else
                targetColumnsDevice(column) = target;
            }, invalidTargets);
      } else {
        Kokkos::parallel_reduce(
            "StructuredRAP: map plan columns", range_policy(
                executionSpace, 0, planFineColumns.size()),
            KOKKOS_LAMBDA(const size_t column, size_t& invalidCount) {
              const GlobalOrdinal fineGid = fineColumnMap.getGlobalElement(
                  fineColumnsDevice(column));
              const LO ownedFine = fineRowMap.getLocalElement(fineGid);
              if (ownedFine == invalid) {
                ++invalidCount;
                return;
              }
              const GlobalOrdinal coarseGid =
                  coarseRowMap.getGlobalElement(coarseRowForFineRow(
                      ownedFine, dofsPerNode, localFineNodes,
                      localCoarseNodes, coarseningRate));
              const LO target = coarseColumnMap.getLocalElement(coarseGid);
              if (target == invalid)
                ++invalidCount;
              else
                targetColumnsDevice(column) = target;
            }, invalidTargets);
      }
      TEUCHOS_TEST_FOR_EXCEPTION(
          invalidTargets != 0, std::runtime_error,
          "StructuredRAPKernel could not map a fine column to its coarse interpolation target.");
      const auto targetColumnsHost = Kokkos::create_mirror_view_and_copy(
          Kokkos::HostSpace(), targetColumnsDevice);
      for (const auto& contribution : pendingContributions) {
        const size_t column = static_cast<size_t>(std::lower_bound(
            planFineColumns.begin(), planFineColumns.end(),
            contribution.fineColumn) - planFineColumns.begin());
        const LO targetColumn = targetColumnsHost(column);
        for (size_t acEntry = contribution.acRowBegin;
             acEntry < contribution.acRowEnd; ++acEntry) {
          if (acEntriesHost(acEntry) == targetColumn) {
            const size_t bucket = contribution.reverseRow * maxAcEntriesPerRow +
                                  acEntry - contribution.acRowBegin;
            reverseBuckets[bucket].push_back(PlanContribution{
                contribution.fineNodeOffset,
                contribution.aRowEntryOffset,
                contribution.directAValueOffset});
            break;
          }
        }
      }

      size_t numContributions = 0;
      for (const auto& bucket : reverseBuckets)
        numContributions += bucket.size();
      PlanArray<size_view> offsets(
          "StructuredRAP: reverse contribution offsets",
          reverseBucketCount + 1);
      PlanArray<size_view> fineNodeOffsets(
          "StructuredRAP: reverse fine node offsets", numContributions);
      PlanArray<size_view> aRowEntryOffsets(
          "StructuredRAP: reverse A row entry offsets", numContributions);
      PlanArray<size_view> directAValueOffsets(
          "StructuredRAP: reverse direct A value offsets", numContributions);
      size_t contribution = 0;
      for (size_t bucket = 0; bucket < reverseBucketCount; ++bucket) {
        offsets.host(bucket) = contribution;
        for (const auto& entry : reverseBuckets[bucket]) {
          fineNodeOffsets.host(contribution) = entry.fineNodeOffset;
          aRowEntryOffsets.host(contribution) = entry.aRowEntryOffset;
          directAValueOffsets.host(contribution) = entry.directAValueOffset;
          ++contribution;
        }
      }
      offsets.host(reverseBucketCount) = contribution;
      reverseContributionOffsets = offsets.copyToDevice(executionSpace);
      reverseFineNodeOffsets = fineNodeOffsets.copyToDevice(executionSpace);
      reverseARowEntryOffsets = aRowEntryOffsets.copyToDevice(executionSpace);
      reverseDirectAValueOffsets = directAValueOffsets.copyToDevice(executionSpace);

      executionSpace.fence("StructuredRAP: numeric entry plan prepared");
    }

    {
      Teuchos::TimeMonitor timer(
          *Teuchos::TimeMonitor::getNewTimer(
              timerPrefix +
              "numeric (constant stencil)"));
      range_policy policy(executionSpace, 0, numCoarseNodes);
      Kokkos::parallel_for(
          "StructuredRAP: constant stencil P^T A P",
          policy,
          KOKKOS_LAMBDA(const size_t coarseNodeIndex) {
            const SC* KOKKOS_RESTRICT aValues = localA.values.data();
            SC* KOKKOS_RESTRICT acValues = localAc.values.data();
            const size_t* KOKKOS_RESTRICT fineNodeOffsets =
                reverseFineNodeOffsets.data();
            const size_t* KOKKOS_RESTRICT aRowEntryOffsets =
                reverseARowEntryOffsets.data();
            const size_t* KOKKOS_RESTRICT directAValueOffsets =
                reverseDirectAValueOffsets.data();
            const LO coarseNode = static_cast<LO>(coarseNodeIndex);
            const LO coarseX = coarseNode % localCoarseNodes[0];
            const LO coarseYZ = coarseNode / localCoarseNodes[0];
            const LO coarseY = coarseYZ % localCoarseNodes[1];
            const LO coarseZ = coarseYZ / localCoarseNodes[1];
            const size_t firstCoarseRow = coarseNodeIndex * numDofs;

            const auto fineXBounds = aggregateBounds(
                coarseX, localFineNodes[0], localCoarseNodes[0], coarseningRate[0]);
            const auto fineYBounds = aggregateBounds(
                coarseY, localFineNodes[1], localCoarseNodes[1], coarseningRate[1]);
            const auto fineZBounds = aggregateBounds(
                coarseZ, localFineNodes[2], localCoarseNodes[2], coarseningRate[2]);
            const size_t firstFineNode =
                static_cast<size_t>(fineXBounds[0]) +
                static_cast<size_t>(localFineNodes[0]) *
                    (static_cast<size_t>(fineYBounds[0]) +
                     static_cast<size_t>(localFineNodes[1]) * fineZBounds[0]);
            const bool useDirectOffsets =
                coarseAggregatePlanClass[0].isInterior(coarseX) &&
                (numDimensions < 2 ||
                 coarseAggregatePlanClass[1].isInterior(coarseY)) &&
                (numDimensions < 3 ||
                 coarseAggregatePlanClass[2].isInterior(coarseZ));
            const size_t reversePlanNodeBegin =
                ((static_cast<size_t>(coarseAggregatePlanClass[2](coarseZ)) *
                      numAggregatePlanClasses[1] +
                  coarseAggregatePlanClass[1](coarseY)) *
                     numAggregatePlanClasses[0] +
                 coarseAggregatePlanClass[0](coarseX)) * numDofs;

            for (size_t rowDof = 0; rowDof < numDofs; ++rowDof) {
              const size_t acRow = firstCoarseRow + rowDof;
              const size_t acRowBegin = localAc.graph.row_map(acRow);
              const size_t acRowEnd = localAc.graph.row_map(acRow + size_t(1));
              const size_t reversePlanRow = reversePlanNodeBegin + rowDof;
              const size_t firstAValue =
                  localA.graph.row_map(firstFineNode * numDofs + rowDof);
              for (size_t acOffset = 0; acOffset < acRowEnd - acRowBegin; ++acOffset) {
                const size_t bucket =
                    reversePlanRow * maxAcEntriesPerRow + acOffset;
                const size_t contributionBegin = reverseContributionOffsets(bucket);
                const size_t contributionEnd = reverseContributionOffsets(bucket + 1);
                SC sum0 = zero, sum1 = zero, sum2 = zero, sum3 = zero;
                size_t contribution = contributionBegin;
                if (useDirectOffsets) {
                  for (; contribution + 3 < contributionEnd; contribution += 4) {
                    sum0 += aValues[firstAValue + directAValueOffsets[contribution]];
                    sum1 += aValues[firstAValue + directAValueOffsets[contribution + 1]];
                    sum2 += aValues[firstAValue + directAValueOffsets[contribution + 2]];
                    sum3 += aValues[firstAValue + directAValueOffsets[contribution + 3]];
                  }
                  SC sum = (sum0 + sum1) + (sum2 + sum3);
                  for (; contribution < contributionEnd; ++contribution)
                    sum += aValues[firstAValue + directAValueOffsets[contribution]];
                  acValues[acRowBegin + acOffset] = sum;
                } else {
                  for (; contribution + 3 < contributionEnd; contribution += 4) {
                    const size_t aRow0 =
                        (firstFineNode + fineNodeOffsets[contribution]) * numDofs + rowDof;
                    const size_t aRow1 =
                        (firstFineNode + fineNodeOffsets[contribution + 1]) * numDofs + rowDof;
                    const size_t aRow2 =
                        (firstFineNode + fineNodeOffsets[contribution + 2]) * numDofs + rowDof;
                    const size_t aRow3 =
                        (firstFineNode + fineNodeOffsets[contribution + 3]) * numDofs + rowDof;
                    sum0 += aValues[localA.graph.row_map(aRow0) +
                                    aRowEntryOffsets[contribution]];
                    sum1 += aValues[localA.graph.row_map(aRow1) +
                                    aRowEntryOffsets[contribution + 1]];
                    sum2 += aValues[localA.graph.row_map(aRow2) +
                                    aRowEntryOffsets[contribution + 2]];
                    sum3 += aValues[localA.graph.row_map(aRow3) +
                                    aRowEntryOffsets[contribution + 3]];
                  }
                  SC sum = (sum0 + sum1) + (sum2 + sum3);
                  for (; contribution < contributionEnd; ++contribution) {
                    const size_t aRow =
                        (firstFineNode + fineNodeOffsets[contribution]) * numDofs + rowDof;
                    sum += aValues[localA.graph.row_map(aRow) +
                                   aRowEntryOffsets[contribution]];
                  }
                  acValues[acRowBegin + acOffset] = sum;
                }
              }
            }
          });
      executionSpace.fence(
          "StructuredRAP: constant stencil product complete");
    }

    if (!Ac.isFillComplete()) {
      Teuchos::TimeMonitor timer(
          *Teuchos::TimeMonitor::getNewTimer(timerPrefix +
                                         "fillComplete"));
      Ac.fillComplete(P.getDomainMap(), P.getDomainMap(), params);
    }
  }
};

}  // namespace Details
}  // namespace MueLu

#endif  // MUELU_STRUCTUREDRAPKERNEL_HPP
