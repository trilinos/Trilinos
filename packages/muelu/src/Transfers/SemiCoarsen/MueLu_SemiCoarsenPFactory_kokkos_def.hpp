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
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_SEMICOARSENPFACTORY_KOKKOS_DEF_HPP
#define MUELU_SEMICOARSENPFACTORY_KOKKOS_DEF_HPP

#include <stdlib.h>

#include <Kokkos_Core.hpp>

#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_LU_Serial_Impl.hpp>
#include <KokkosBatched_Trsm_Decl.hpp>
#include <KokkosBatched_Trsm_Serial_Impl.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_SemiCoarsenPFactory_kokkos_decl.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal,
          class DeviceType>
RCP<const ParameterList>
SemiCoarsenPFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal,
                           Kokkos::Compat::KokkosDeviceWrapperNode<
                               DeviceType>>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  std::string name = "semicoarsen: coarsen rate";
  validParamList->setEntry(name, MasterList::getEntry(name));
  validParamList->set<RCP<const FactoryBase>>(
      "A", Teuchos::null, "Generating factory of the matrix A");
  validParamList->set<RCP<const FactoryBase>>(
      "Nullspace", Teuchos::null, "Generating factory of the nullspace");
  validParamList->set<RCP<const FactoryBase>>(
      "Coordinates", Teuchos::null, "Generating factory for coordinates");

  validParamList->set<RCP<const FactoryBase>>(
      "LineDetection_VertLineIds", Teuchos::null,
      "Generating factory for LineDetection vertical line ids");
  validParamList->set<RCP<const FactoryBase>>(
      "LineDetection_Layers", Teuchos::null,
      "Generating factory for LineDetection layer ids");
  validParamList->set<RCP<const FactoryBase>>(
      "CoarseNumZLayers", Teuchos::null,
      "Generating factory for number of coarse z-layers");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal,
          class DeviceType>
void SemiCoarsenPFactory_kokkos<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::
    DeclareInput(Level &fineLevel, Level & /* coarseLevel */) const {
  Input(fineLevel, "A");
  Input(fineLevel, "Nullspace");

  Input(fineLevel, "LineDetection_VertLineIds");
  Input(fineLevel, "LineDetection_Layers");
  Input(fineLevel, "CoarseNumZLayers");

  // check whether fine level coordinate information is available.
  // If yes, request the fine level coordinates and generate coarse coordinates
  // during the Build call
  if (fineLevel.GetLevelID() == 0) {
    if (fineLevel.IsAvailable("Coordinates", NoFactory::get())) {
      fineLevel.DeclareInput("Coordinates", NoFactory::get(), this);
      bTransferCoordinates_ = true;
    }
  } else if (bTransferCoordinates_ == true) {
    // on coarser levels we check the default factory providing "Coordinates"
    // or the factory declared to provide "Coordinates"
    // first, check which factory is providing coordinate information
    RCP<const FactoryBase> myCoordsFact = GetFactory("Coordinates");
    if (myCoordsFact == Teuchos::null) {
      myCoordsFact = fineLevel.GetFactoryManager()->GetFactory("Coordinates");
    }
    if (fineLevel.IsAvailable("Coordinates", myCoordsFact.get())) {
      fineLevel.DeclareInput("Coordinates", myCoordsFact.get(), this);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal,
          class DeviceType>
void SemiCoarsenPFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal,
                                Kokkos::Compat::KokkosDeviceWrapperNode<
                                    DeviceType>>::Build(Level &fineLevel,
                                                        Level &coarseLevel)
    const {
  return BuildP(fineLevel, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal,
          class DeviceType>
void SemiCoarsenPFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal,
                                Kokkos::Compat::KokkosDeviceWrapperNode<
                                    DeviceType>>::BuildP(Level &fineLevel,
                                                         Level &coarseLevel)
    const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  // obtain general variables
  RCP<Matrix> A = Get<RCP<Matrix>>(fineLevel, "A");
  RCP<MultiVector> fineNullspace =
      Get<RCP<MultiVector>>(fineLevel, "Nullspace");

  // get user-provided coarsening rate parameter (constant over all levels)
  const ParameterList &pL = GetParameterList();
  LO CoarsenRate = as<LO>(pL.get<int>("semicoarsen: coarsen rate"));
  TEUCHOS_TEST_FOR_EXCEPTION(
      CoarsenRate < 2, Exceptions::RuntimeError,
      "semicoarsen: coarsen rate must be greater than 1");

  // collect general input data
  LO BlkSize = A->GetFixedBlockSize();
  RCP<const Map> rowMap = A->getRowMap();
  LO Ndofs = rowMap->getLocalNumElements();
  LO Nnodes = Ndofs / BlkSize;

  // collect line detection information generated by the LineDetectionFactory
  // instance
  LO FineNumZLayers = Get<LO>(fineLevel, "CoarseNumZLayers");
  Teuchos::ArrayRCP<LO> TVertLineId =
      Get<Teuchos::ArrayRCP<LO>>(fineLevel, "LineDetection_VertLineIds");
  Teuchos::ArrayRCP<LO> TLayerId =
      Get<Teuchos::ArrayRCP<LO>>(fineLevel, "LineDetection_Layers");

  // compute number of coarse layers
  TEUCHOS_TEST_FOR_EXCEPTION(FineNumZLayers < 2, Exceptions::RuntimeError,
                             "Cannot coarsen further");
  using coordT = typename Teuchos::ScalarTraits<Scalar>::coordinateType;
  LO CoarseNumZLayers =
      (LO)floor(((coordT)(FineNumZLayers + 1)) / ((coordT)CoarsenRate) - 1.0);
  if (CoarseNumZLayers < 1)
    CoarseNumZLayers = 1;

  // generate transfer operator with semicoarsening
  RCP<Matrix> P;
  RCP<MultiVector> coarseNullspace;
  BuildSemiCoarsenP(coarseLevel, Ndofs, Nnodes, BlkSize, FineNumZLayers,
                    CoarseNumZLayers, TLayerId, TVertLineId, A, fineNullspace, P,
                    coarseNullspace);

  // Store number of coarse z-layers on the coarse level container
  // This information is used by the LineDetectionAlgorithm
  // TODO get rid of the NoFactory
  coarseLevel.Set("NumZLayers", Teuchos::as<LO>(CoarseNumZLayers),
                  MueLu::NoFactory::get());

  // store semicoarsening transfer on coarse level
  Set(coarseLevel, "P", P);
  Set(coarseLevel, "Nullspace", coarseNullspace);

  // transfer coordinates
  if (bTransferCoordinates_) {
    SubFactoryMonitor m2(*this, "TransferCoordinates", coarseLevel);
    typedef Xpetra::MultiVector<
        typename Teuchos::ScalarTraits<Scalar>::coordinateType, LO, GO, NO>
        xdMV;
    RCP<xdMV> fineCoords = Teuchos::null;
    if (fineLevel.GetLevelID() == 0 &&
        fineLevel.IsAvailable("Coordinates", NoFactory::get())) {
      fineCoords = fineLevel.Get<RCP<xdMV>>("Coordinates", NoFactory::get());
    } else {
      RCP<const FactoryBase> myCoordsFact = GetFactory("Coordinates");
      if (myCoordsFact == Teuchos::null) {
        myCoordsFact = fineLevel.GetFactoryManager()->GetFactory("Coordinates");
      }
      if (fineLevel.IsAvailable("Coordinates", myCoordsFact.get())) {
        fineCoords =
            fineLevel.Get<RCP<xdMV>>("Coordinates", myCoordsFact.get());
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION(fineCoords == Teuchos::null,
                               Exceptions::RuntimeError,
                               "No Coordinates found provided by the user.");

    TEUCHOS_TEST_FOR_EXCEPTION(fineCoords->getNumVectors() != 3,
                               Exceptions::RuntimeError,
                               "Three coordinates arrays must be supplied if "
                               "line detection orientation not given.");
    ArrayRCP<typename Teuchos::ScalarTraits<Scalar>::coordinateType> x =
        fineCoords->getDataNonConst(0);
    ArrayRCP<typename Teuchos::ScalarTraits<Scalar>::coordinateType> y =
        fineCoords->getDataNonConst(1);
    ArrayRCP<typename Teuchos::ScalarTraits<Scalar>::coordinateType> z =
        fineCoords->getDataNonConst(2);

    // determine the maximum and minimum z coordinate value on the current
    // processor.
    typename Teuchos::ScalarTraits<Scalar>::coordinateType zval_max =
        -Teuchos::ScalarTraits<
            typename Teuchos::ScalarTraits<Scalar>::coordinateType>::one() /
        Teuchos::ScalarTraits<
            typename Teuchos::ScalarTraits<Scalar>::coordinateType>::sfmin();
    typename Teuchos::ScalarTraits<Scalar>::coordinateType zval_min =
        Teuchos::ScalarTraits<
            typename Teuchos::ScalarTraits<Scalar>::coordinateType>::one() /
        Teuchos::ScalarTraits<
            typename Teuchos::ScalarTraits<Scalar>::coordinateType>::sfmin();
    for (auto it = z.begin(); it != z.end(); ++it) {
      if (*it > zval_max)
        zval_max = *it;
      if (*it < zval_min)
        zval_min = *it;
    }

    LO myCoarseZLayers = Teuchos::as<LO>(CoarseNumZLayers);

    ArrayRCP<typename Teuchos::ScalarTraits<Scalar>::coordinateType>
        myZLayerCoords = Teuchos::arcp<
            typename Teuchos::ScalarTraits<Scalar>::coordinateType>(
            myCoarseZLayers);
    if (myCoarseZLayers == 1) {
      myZLayerCoords[0] = zval_min;
    } else {
      typename Teuchos::ScalarTraits<Scalar>::coordinateType dz =
          (zval_max - zval_min) / (myCoarseZLayers - 1);
      for (LO k = 0; k < myCoarseZLayers; ++k) {
        myZLayerCoords[k] = k * dz;
      }
    }

    // Note, that the coarse level node coordinates have to be in vertical
    // ordering according to the numbering of the vertical lines

    // number of vertical lines on current node:
    LO numVertLines = Nnodes / FineNumZLayers;
    LO numLocalCoarseNodes = numVertLines * myCoarseZLayers;

    RCP<const Map> coarseCoordMap = MapFactory::Build(
        fineCoords->getMap()->lib(),
        Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
        Teuchos::as<size_t>(numLocalCoarseNodes),
        fineCoords->getMap()->getIndexBase(), fineCoords->getMap()->getComm());
    RCP<xdMV> coarseCoords = Xpetra::MultiVectorFactory<
        typename Teuchos::ScalarTraits<Scalar>::coordinateType, LO, GO,
        NO>::Build(coarseCoordMap, fineCoords->getNumVectors());
    coarseCoords->putScalar(-1.0);
    ArrayRCP<typename Teuchos::ScalarTraits<Scalar>::coordinateType> cx =
        coarseCoords->getDataNonConst(0);
    ArrayRCP<typename Teuchos::ScalarTraits<Scalar>::coordinateType> cy =
        coarseCoords->getDataNonConst(1);
    ArrayRCP<typename Teuchos::ScalarTraits<Scalar>::coordinateType> cz =
        coarseCoords->getDataNonConst(2);

    // loop over all vert line indices (stop as soon as possible)
    LO cntCoarseNodes = 0;
    for (LO vt = 0; vt < TVertLineId.size(); ++vt) {
      // vertical line id in *vt
      LO curVertLineId = TVertLineId[vt];

      if (cx[curVertLineId * myCoarseZLayers] == -1.0 &&
          cy[curVertLineId * myCoarseZLayers] == -1.0) {
        // loop over all local myCoarseZLayers
        for (LO n = 0; n < myCoarseZLayers; ++n) {
          cx[curVertLineId * myCoarseZLayers + n] = x[vt];
          cy[curVertLineId * myCoarseZLayers + n] = y[vt];
          cz[curVertLineId * myCoarseZLayers + n] = myZLayerCoords[n];
        }
        cntCoarseNodes += myCoarseZLayers;
      }

      TEUCHOS_TEST_FOR_EXCEPTION(cntCoarseNodes > numLocalCoarseNodes,
                                 Exceptions::RuntimeError,
                                 "number of coarse nodes is inconsistent.");
      if (cntCoarseNodes == numLocalCoarseNodes)
        break;
    }

    // set coarse level coordinates
    Set(coarseLevel, "Coordinates", coarseCoords);
  } /* end bool bTransferCoordinates */
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal,
          class DeviceType>
void SemiCoarsenPFactory_kokkos<
    Scalar, LocalOrdinal, GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::
    BuildSemiCoarsenP(Level &coarseLevel, const LO NFRows, const LO NFNodes,
                      const LO DofsPerNode, const LO NFLayers,
                      const LO NCLayers, const ArrayRCP<LO> LayerId,
                      const ArrayRCP<LO> VertLineId, const RCP<Matrix> &Amat,
                      const RCP<MultiVector> fineNullspace, RCP<Matrix> &P,
                      RCP<MultiVector> &coarseNullspace) const {
  SubFactoryMonitor m2(*this, "BuildSemiCoarsenP", coarseLevel);
  using impl_SC = typename Kokkos::ArithTraits<SC>::val_type;
  using impl_ATS = Kokkos::ArithTraits<impl_SC>;
  using LOView1D = Kokkos::View<LO *, DeviceType>;
  using LOView2D = Kokkos::View<LO **, DeviceType>;

  // Construct a map from fine level column to layer ids (including ghost nodes)
  // Note: this is needed to sum all couplings within a layer
  const auto FCol2LayerVector =
      Xpetra::VectorFactory<LO, LO, GO, NO>::Build(Amat->getColMap());
  const auto localTemp =
      Xpetra::VectorFactory<LO, LO, GO, NO>::Build(Amat->getDomainMap());
  RCP<const Import> importer = Amat->getCrsGraph()->getImporter();
  if (importer == Teuchos::null)
    importer = ImportFactory::Build(Amat->getDomainMap(), Amat->getColMap());
  {
    // Fill local temp with layer ids and fill ghost nodes
    const auto localTempHost = localTemp->getHostLocalView(Xpetra::Access::ReadWrite);
    for (int row = 0; row < NFRows; row++)
      localTempHost(row, 0) = LayerId[row / DofsPerNode];
    const auto localTempView = localTemp->getDeviceLocalView(Xpetra::Access::ReadWrite);
    Kokkos::deep_copy(localTempView, localTempHost);
    FCol2LayerVector->doImport(*localTemp, *(importer), Xpetra::INSERT);
  }
  const auto FCol2LayerView = FCol2LayerVector->getDeviceLocalView(Xpetra::Access::ReadOnly);
  const auto FCol2Layer = Kokkos::subview(FCol2LayerView, Kokkos::ALL(), 0);

  // Construct a map from fine level column to local dof per node id (including
  // ghost nodes) Note: this is needed to sum all couplings within a layer
  const auto FCol2DofVector =
      Xpetra::VectorFactory<LO, LO, GO, NO>::Build(Amat->getColMap());
  {
    // Fill local temp with local dof per node ids and fill ghost nodes
    const auto localTempHost = localTemp->getHostLocalView(Xpetra::Access::ReadWrite);
    for (int row = 0; row < NFRows; row++)
      localTempHost(row, 0) = row % DofsPerNode;
    const auto localTempView = localTemp->getDeviceLocalView(Xpetra::Access::ReadWrite);
    Kokkos::deep_copy(localTempView, localTempHost);
    FCol2DofVector->doImport(*localTemp, *(importer), Xpetra::INSERT);
  }
  const auto FCol2DofView = FCol2DofVector->getDeviceLocalView(Xpetra::Access::ReadOnly);
  const auto FCol2Dof = Kokkos::subview(FCol2DofView, Kokkos::ALL(), 0);

  // Compute NVertLines
  // TODO: Read this from line detection factory
  int NVertLines = -1;
  if (NFNodes != 0)
    NVertLines = VertLineId[0];
  for (int node = 1; node < NFNodes; ++node)
    if (VertLineId[node] > NVertLines)
      NVertLines = VertLineId[node];
  NVertLines++;

  // Construct a map from Line, Layer ids to fine level node
  LOView2D LineLayer2Node("LineLayer2Node", NVertLines, NFLayers);
  typename LOView2D::HostMirror LineLayer2NodeHost =
      Kokkos::create_mirror_view(LineLayer2Node);
  for (int node = 0; node < NFNodes; ++node)
    LineLayer2NodeHost(VertLineId[node], LayerId[node]) = node;
  Kokkos::deep_copy(LineLayer2Node, LineLayer2NodeHost);

  // Construct a map from coarse layer id to fine layer id
  LOView1D CLayer2FLayer("CLayer2FLayer", NCLayers);
  typename LOView1D::HostMirror CLayer2FLayerHost =
      Kokkos::create_mirror_view(CLayer2FLayer);
  using coordT = typename Teuchos::ScalarTraits<Scalar>::coordinateType;
  const LO FirstStride =
      (LO)ceil(((coordT)(NFLayers + 1)) / ((coordT)(NCLayers + 1)));
  const coordT RestStride =
      ((coordT)(NFLayers - FirstStride + 1)) / ((coordT)NCLayers);
  const LO NCpts =
      (LO)floor((((coordT)(NFLayers - FirstStride + 1)) / RestStride) + .00001);
  TEUCHOS_TEST_FOR_EXCEPTION(NCLayers != NCpts, Exceptions::RuntimeError,
                             "sizes do not match.");
  coordT stride = (coordT)FirstStride;
  for (int clayer = 0; clayer < NCLayers; ++clayer) {
    CLayer2FLayerHost(clayer) = (LO)floor(stride) - 1;
    stride += RestStride;
  }
  Kokkos::deep_copy(CLayer2FLayer, CLayer2FLayerHost);

  // Compute start layer and stencil sizes for layer interpolation at each
  // coarse layer
  int MaxStencilSize = 1;
  LOView1D CLayer2StartLayer("CLayer2StartLayer", NCLayers);
  LOView1D CLayer2StencilSize("CLayer2StencilSize", NCLayers);
  typename LOView1D::HostMirror CLayer2StartLayerHost =
      Kokkos::create_mirror_view(CLayer2StartLayer);
  typename LOView1D::HostMirror CLayer2StencilSizeHost =
      Kokkos::create_mirror_view(CLayer2StencilSize);
  for (int clayer = 0; clayer < NCLayers; ++clayer) {
    const int startLayer = (clayer > 0) ? CLayer2FLayerHost(clayer - 1) + 1 : 0;
    const int stencilSize = (clayer < NCLayers - 1)
                                ? CLayer2FLayerHost(clayer + 1) - startLayer
                                : NFLayers - startLayer;

    if (MaxStencilSize < stencilSize)
      MaxStencilSize = stencilSize;
    CLayer2StartLayerHost(clayer) = startLayer;
    CLayer2StencilSizeHost(clayer) = stencilSize;
  }
  Kokkos::deep_copy(CLayer2StartLayer, CLayer2StartLayerHost);
  Kokkos::deep_copy(CLayer2StencilSize, CLayer2StencilSizeHost);

  // Allocate storage for the coarse layer interpolation matrices on all
  // vertical lines Note: Contributions to each matrix are collapsed to vertical
  // lines. Thus, each vertical line gives rise to a block tridiagonal matrix.
  // Here we store the full matrix to be compatible with kokkos kernels batch LU
  // and tringular solve.
  int Nmax = MaxStencilSize * DofsPerNode;
  Kokkos::View<impl_SC ***, DeviceType> BandMat(
      "BandMat", NVertLines, Nmax, Nmax);
  Kokkos::View<impl_SC ***, DeviceType> BandSol(
      "BandSol", NVertLines, Nmax, DofsPerNode);

  // Precompute number of nonzeros in prolongation matrix and allocate P views
  // Note: Each coarse dof (NVertLines*NCLayers*DofsPerNode) contributes an
  // interpolation stencil (StencilSize*DofsPerNode)
  int NnzP = 0;
  for (int clayer = 0; clayer < NCLayers; ++clayer)
    NnzP += CLayer2StencilSizeHost(clayer);
  NnzP *= NVertLines * DofsPerNode * DofsPerNode;
  Kokkos::View<impl_SC *, DeviceType> Pvals("Pvals", NnzP);
  Kokkos::View<LO *, DeviceType> Pcols("Pcols", NnzP);

  // Precompute Pptr
  // Note: Each coarse layer stencil dof contributes DofsPerNode to the
  // corresponding row in P
  Kokkos::View<size_t *, DeviceType> Pptr("Pptr", NFRows + 1);
  typename Kokkos::View<size_t *, DeviceType>::HostMirror PptrHost =
      Kokkos::create_mirror_view(Pptr);
  Kokkos::deep_copy(PptrHost, 0);
  for (int line = 0; line < NVertLines; ++line) {
    for (int clayer = 0; clayer < NCLayers; ++clayer) {
      const int stencilSize = CLayer2StencilSizeHost(clayer);
      const int startLayer = CLayer2StartLayerHost(clayer);
      for (int snode = 0; snode < stencilSize; ++snode) {
        for (int dofi = 0; dofi < DofsPerNode; ++dofi) {
          const int layer = startLayer + snode;
          const int AmatBlkRow = LineLayer2NodeHost(line, layer);
          const int AmatRow = AmatBlkRow * DofsPerNode + dofi;
          PptrHost(AmatRow + 1) += DofsPerNode;
        }
      }
    }
  }
  for (int i = 2; i < NFRows + 1; ++i)
    PptrHost(i) += PptrHost(i - 1);
  TEUCHOS_TEST_FOR_EXCEPTION(NnzP != (int)PptrHost(NFRows),
                             Exceptions::RuntimeError,
                             "Number of nonzeros in P does not match");
  Kokkos::deep_copy(Pptr, PptrHost);

  // Precompute Pptr offsets
  // Note: These are used to determine the nonzero index in Pvals and Pcols
  Kokkos::View<LO *, Kokkos::DefaultHostExecutionSpace> layerBuckets(
      "layerBuckets", NFLayers);
  Kokkos::deep_copy(layerBuckets, 0);
  LOView2D CLayerSNode2PptrOffset("CLayerSNode2PptrOffset", NCLayers,
                                  MaxStencilSize);
  typename LOView2D::HostMirror CLayerSNode2PptrOffsetHost =
      Kokkos::create_mirror_view(CLayerSNode2PptrOffset);
  for (int clayer = 0; clayer < NCLayers; ++clayer) {
    const int stencilSize = CLayer2StencilSizeHost(clayer);
    const int startLayer = CLayer2StartLayerHost(clayer);
    for (int snode = 0; snode < stencilSize; ++snode) {
      const int layer = startLayer + snode;
      CLayerSNode2PptrOffsetHost(clayer, snode) = layerBuckets(layer);
      layerBuckets(layer)++;
    }
  }
  Kokkos::deep_copy(CLayerSNode2PptrOffset, CLayerSNode2PptrOffsetHost);

  { // Fill P - fill and solve each block tridiagonal system and fill P views
    SubFactoryMonitor m3(*this, "Fill P", coarseLevel);

    const auto localAmat = Amat->getLocalMatrixDevice();
    const auto zero = impl_ATS::zero();
    const auto one = impl_ATS::one();

    using range_policy = Kokkos::RangePolicy<execution_space>;
    Kokkos::parallel_for(
        "MueLu::SemiCoarsenPFactory_kokkos::BuildSemiCoarsenP Fill P",
        range_policy(0, NVertLines), KOKKOS_LAMBDA(const int line) {
          for (int clayer = 0; clayer < NCLayers; ++clayer) {

            // Initialize BandSol
            auto bandSol =
                Kokkos::subview(BandSol, line, Kokkos::ALL(), Kokkos::ALL());
            for (int row = 0; row < Nmax; ++row)
              for (int dof = 0; dof < DofsPerNode; ++dof)
                bandSol(row, dof) = zero;

            // Initialize BandMat (set unused row diagonal to 1.0)
            const int stencilSize = CLayer2StencilSize(clayer);
            const int N = stencilSize * DofsPerNode;
            auto bandMat =
                Kokkos::subview(BandMat, line, Kokkos::ALL(), Kokkos::ALL());
            for (int row = 0; row < Nmax; ++row)
              for (int col = 0; col < Nmax; ++col)
                bandMat(row, col) =
                    (row == col && row >= N) ? one : zero;

            // Loop over layers in stencil and fill banded matrix and rhs
            const int flayer = CLayer2FLayer(clayer);
            const int startLayer = CLayer2StartLayer(clayer);
            for (int snode = 0; snode < stencilSize; ++snode) {

              const int layer = startLayer + snode;
              if (layer == flayer) { // If layer in stencil is a coarse layer
                for (int dof = 0; dof < DofsPerNode; ++dof) {
                  const int row = snode * DofsPerNode + dof;
                  bandMat(row, row) = one;
                  bandSol(row, dof) = one;
                }
              } else { // Not a coarse layer
                const int AmatBlkRow = LineLayer2Node(line, layer);
                for (int dofi = 0; dofi < DofsPerNode; ++dofi) {

                  // Get Amat row info
                  const int AmatRow = AmatBlkRow * DofsPerNode + dofi;
                  const auto localAmatRow = localAmat.rowConst(AmatRow);
                  const int AmatRowLeng = localAmatRow.length;

                  const int row = snode * DofsPerNode + dofi;
                  for (int dofj = 0; dofj < DofsPerNode; ++dofj) {
                    const int col = snode * DofsPerNode + dofj;

                    // Sum values along row which correspond to stencil
                    // layer/dof and fill bandMat
                    auto val = zero;
                    for (int i = 0; i < AmatRowLeng; ++i) {
                      const int colidx = localAmatRow.colidx(i);
                      if (FCol2Layer(colidx) == layer &&
                          FCol2Dof(colidx) == dofj)
                        val += localAmatRow.value(i);
                    }
                    bandMat(row, col) = val;

                    if (snode > 0) {
                      // Sum values along row which correspond to stencil
                      // layer/dof below and fill bandMat
                      val = zero;
                      for (int i = 0; i < AmatRowLeng; ++i) {
                        const int colidx = localAmatRow.colidx(i);
                        if (FCol2Layer(colidx) == layer - 1 &&
                            FCol2Dof(colidx) == dofj)
                          val += localAmatRow.value(i);
                      }
                      bandMat(row, col - DofsPerNode) = val;
                    }

                    if (snode < stencilSize - 1) {
                      // Sum values along row which correspond to stencil
                      // layer/dof above and fill bandMat
                      val = zero;
                      for (int i = 0; i < AmatRowLeng; ++i) {
                        const int colidx = localAmatRow.colidx(i);
                        if (FCol2Layer(colidx) == layer + 1 &&
                            FCol2Dof(colidx) == dofj)
                          val += localAmatRow.value(i);
                      }
                      bandMat(row, col + DofsPerNode) = val;
                    }
                  }
                }
              }
            }

            // Batch LU and triangular solves
            namespace KB = KokkosBatched;
            using lu_type = typename KB::SerialLU<KB::Algo::LU::Unblocked>;
            lu_type::invoke(bandMat);
            using trsv_l_type =
                typename KB::SerialTrsm<KB::Side::Left, KB::Uplo::Lower,
                                        KB::Trans::NoTranspose, KB::Diag::Unit,
                                        KB::Algo::Trsm::Unblocked>;
            trsv_l_type::invoke(one, bandMat, bandSol);
            using trsv_u_type = typename KB::SerialTrsm<
                KB::Side::Left, KB::Uplo::Upper, KB::Trans::NoTranspose,
                KB::Diag::NonUnit, KB::Algo::Trsm::Unblocked>;
            trsv_u_type::invoke(one, bandMat, bandSol);

            // Fill prolongation views with solution
            for (int snode = 0; snode < stencilSize; ++snode) {
              for (int dofj = 0; dofj < DofsPerNode; ++dofj) {
                for (int dofi = 0; dofi < DofsPerNode; ++dofi) {
                  const int layer = startLayer + snode;
                  const int AmatBlkRow = LineLayer2Node(line, layer);
                  const int AmatRow = AmatBlkRow * DofsPerNode + dofi;

                  const int pptrOffset = CLayerSNode2PptrOffset(clayer, snode);
                  const int pptr =
                      Pptr(AmatRow) + pptrOffset * DofsPerNode + dofj;

                  const int col =
                      line * NCLayers + clayer; // coarse node (block row) index
                  Pcols(pptr) = col * DofsPerNode + dofj;
                  Pvals(pptr) = bandSol(snode * DofsPerNode + dofi, dofj);
                }
              }
            }
          }
        });
  } // Fill P

  // Build P
  RCP<const Map> rowMap = Amat->getRowMap();
  Xpetra::global_size_t GNdofs = rowMap->getGlobalNumElements();
  Xpetra::global_size_t itemp = GNdofs / NFLayers;
  std::vector<size_t> stridingInfo_{(size_t)DofsPerNode};
  RCP<const Map> coarseMap = StridedMapFactory::Build(
      rowMap->lib(), NCLayers * itemp, NCLayers * NVertLines * DofsPerNode, 0,
      stridingInfo_, rowMap->getComm(), -1, 0);
  P = rcp(new CrsMatrixWrap(rowMap, coarseMap, 0));
  RCP<CrsMatrix> PCrs = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();
  PCrs->setAllValues(Pptr, Pcols, Pvals);
  PCrs->expertStaticFillComplete(coarseMap, Amat->getDomainMap());

  // set StridingInformation of P
  if (Amat->IsView("stridedMaps") == true)
    P->CreateView("stridedMaps", Amat->getRowMap("stridedMaps"), coarseMap);
  else
    P->CreateView("stridedMaps", P->getRangeMap(), coarseMap);

  // Construct coarse nullspace and inject fine nullspace
  coarseNullspace =
      MultiVectorFactory::Build(coarseMap, fineNullspace->getNumVectors());
  const int numVectors = fineNullspace->getNumVectors();
  const auto fineNullspaceView = fineNullspace->getDeviceLocalView(Xpetra::Access::ReadOnly);
  const auto coarseNullspaceView = coarseNullspace->getDeviceLocalView(Xpetra::Access::ReadWrite);
  using range_policy = Kokkos::RangePolicy<execution_space>;
  Kokkos::parallel_for(
      "MueLu::SemiCoarsenPFactory_kokkos::BuildSemiCoarsenP Inject Nullspace",
      range_policy(0, NVertLines), KOKKOS_LAMBDA(const int line) {
        for (int clayer = 0; clayer < NCLayers; ++clayer) {
          const int layer = CLayer2FLayer(clayer);
          const int AmatBlkRow =
              LineLayer2Node(line, layer); // fine node (block row) index
          const int col =
              line * NCLayers + clayer; // coarse node (block row) index
          for (int k = 0; k < numVectors; ++k) {
            for (int dofi = 0; dofi < DofsPerNode; ++dofi) {
              const int fRow = AmatBlkRow * DofsPerNode + dofi;
              const int cRow = col * DofsPerNode + dofi;
              coarseNullspaceView(cRow, k) = fineNullspaceView(fRow, k);
            }
          }
        }
      });
}

} // namespace MueLu

#define MUELU_SEMICOARSENPFACTORY_KOKKOS_SHORT
#endif // MUELU_SEMICOARSENPFACTORY_KOKKOS_DEF_HPP
