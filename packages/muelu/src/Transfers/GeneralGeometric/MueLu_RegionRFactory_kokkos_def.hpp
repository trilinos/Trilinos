// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REGIONRFACTORY_KOKKOS_DEF_HPP
#define MUELU_REGIONRFACTORY_KOKKOS_DEF_HPP

#include "Kokkos_UnorderedMap.hpp"

#include "MueLu_RegionRFactory_kokkos_decl.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsGraphFactory.hpp>

#include "MueLu_Types.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> RegionRFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null,
                                               "Generating factory of the matrix A");
  validParamList->set<RCP<const FactoryBase> >("numDimensions", Teuchos::null,
                                               "Number of spatial dimensions in the problem.");
  validParamList->set<RCP<const FactoryBase> >("lNodesPerDim", Teuchos::null,
                                               "Number of local nodes per spatial dimension on the fine grid.");
  validParamList->set<RCP<const FactoryBase> >("Nullspace", Teuchos::null,
                                               "Fine level nullspace used to construct the coarse level nullspace.");
  validParamList->set<RCP<const FactoryBase> >("Coordinates", Teuchos::null,
                                               "Fine level coordinates used to construct piece-wise linear prolongator and coarse level coordinates.");
  validParamList->set<bool>("keep coarse coords", false, "Flag to keep coordinates for special coarse grid solve");

  return validParamList;
}  // GetValidParameterList()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RegionRFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    DeclareInput(Level& fineLevel, Level& /* coarseLevel */) const {
  Input(fineLevel, "A");
  Input(fineLevel, "numDimensions");
  Input(fineLevel, "lNodesPerDim");
  Input(fineLevel, "Nullspace");
  Input(fineLevel, "Coordinates");

}  // DeclareInput()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RegionRFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Build(Level& fineLevel, Level& coarseLevel) const {
  // Set debug outputs based on environment variable
  RCP<Teuchos::FancyOStream> out;
  if (const char* dbg = std::getenv("MUELU_REGIONRFACTORY_DEBUG")) {
    out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
  } else {
    out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
  }

  *out << "Starting RegionRFactory_kokkos::Build." << std::endl;

  // First get the inputs from the fineLevel
  const int numDimensions = Get<int>(fineLevel, "numDimensions");
  Array<LO> lFineNodesPerDim(3, Teuchos::OrdinalTraits<LO>::one());
  {
    Array<LO> lNodesPerDim = Get<Array<LO> >(fineLevel, "lNodesPerDim");
    for (int dim = 0; dim < numDimensions; ++dim) {
      lFineNodesPerDim[dim] = lNodesPerDim[dim];
    }
  }
  *out << "numDimensions " << numDimensions << " and lFineNodesPerDim: " << lFineNodesPerDim
       << std::endl;

  // Let us check that the inputs verify our assumptions
  if (numDimensions < 1 || numDimensions > 3) {
    throw std::runtime_error("numDimensions must be 1, 2 or 3!");
  }
  for (int dim = 0; dim < numDimensions; ++dim) {
    if (lFineNodesPerDim[dim] % 3 != 1) {
      throw std::runtime_error("The number of fine node in each direction need to be 3n+1");
    }
  }
  Array<LO> lCoarseNodesPerDim(3, Teuchos::OrdinalTraits<LO>::one());

  const RCP<Matrix> A = Get<RCP<Matrix> >(fineLevel, "A");

  RCP<realvaluedmultivector_type> fineCoordinates, coarseCoordinates;
  fineCoordinates = Get<RCP<realvaluedmultivector_type> >(fineLevel, "Coordinates");
  if (static_cast<int>(fineCoordinates->getNumVectors()) != numDimensions) {
    throw std::runtime_error("The number of vectors in the coordinates is not equal to numDimensions!");
  }

  // Let us create R and pass it down to the
  // appropriate specialization and see what we
  // get back!
  RCP<Matrix> R;

  if (numDimensions == 1) {
    throw std::runtime_error("RegionRFactory_kokkos no implemented for 1D case yet.");
  } else if (numDimensions == 2) {
    throw std::runtime_error("RegionRFactory_kokkos no implemented for 2D case yet.");
  } else if (numDimensions == 3) {
    Build3D(numDimensions, lFineNodesPerDim, A, fineCoordinates,
            R, coarseCoordinates, lCoarseNodesPerDim);
  }

  const Teuchos::ParameterList& pL = GetParameterList();

  // Reuse pattern if available (multiple solve)
  RCP<ParameterList> Tparams;
  if (pL.isSublist("matrixmatrix: kernel params"))
    Tparams = rcp(new ParameterList(pL.sublist("matrixmatrix: kernel params")));
  else
    Tparams = rcp(new ParameterList);

  // R->describe(*out, Teuchos::VERB_EXTREME);
  *out << "Compute P=R^t" << std::endl;
  // By default, we don't need global constants for transpose
  Tparams->set("compute global constants: temporaries", Tparams->get("compute global constants: temporaries", false));
  Tparams->set("compute global constants", Tparams->get("compute global constants", false));
  std::string label = "MueLu::RegionR-transR" + Teuchos::toString(coarseLevel.GetLevelID());
  RCP<Matrix> P     = Utilities::Transpose(*R, true, label, Tparams);

  *out << "Compute coarse nullspace" << std::endl;
  RCP<MultiVector> fineNullspace   = Get<RCP<MultiVector> >(fineLevel, "Nullspace");
  RCP<MultiVector> coarseNullspace = MultiVectorFactory::Build(R->getRowMap(),
                                                               fineNullspace->getNumVectors());
  R->apply(*fineNullspace, *coarseNullspace, Teuchos::NO_TRANS, Teuchos::ScalarTraits<SC>::one(),
           Teuchos::ScalarTraits<SC>::zero());

  *out << "Set data on coarse level" << std::endl;
  Set(coarseLevel, "numDimensions", numDimensions);
  Set(coarseLevel, "lNodesPerDim", lCoarseNodesPerDim);
  Set(coarseLevel, "Nullspace", coarseNullspace);
  Set(coarseLevel, "Coordinates", coarseCoordinates);
  if (pL.get<bool>("keep coarse coords")) {
    coarseLevel.Set<RCP<realvaluedmultivector_type> >("Coordinates2", coarseCoordinates, NoFactory::get());
  }

  R->SetFixedBlockSize(A->GetFixedBlockSize());
  P->SetFixedBlockSize(A->GetFixedBlockSize());

  Set(coarseLevel, "R", R);
  Set(coarseLevel, "P", P);

}  // Build()

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RegionRFactory_kokkos<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Build3D(const int numDimensions,
            Teuchos::Array<LocalOrdinal>& lFineNodesPerDim,
            const RCP<Matrix>& A,
            const RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> >& fineCoordinates,
            RCP<Matrix>& R,
            RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> >& coarseCoordinates,
            Teuchos::Array<LocalOrdinal>& lCoarseNodesPerDim) const {
  using local_matrix_type = typename CrsMatrix::local_matrix_type;
  using local_graph_type  = typename local_matrix_type::staticcrsgraph_type;
  using row_map_type      = typename local_matrix_type::row_map_type::non_const_type;
  using entries_type      = typename local_matrix_type::index_type::non_const_type;
  using values_type       = typename local_matrix_type::values_type::non_const_type;
  using impl_scalar_type  = typename Kokkos::ArithTraits<Scalar>::val_type;

  // Set debug outputs based on environment variable
  RCP<Teuchos::FancyOStream> out;
  if (const char* dbg = std::getenv("MUELU_REGIONRFACTORY_DEBUG")) {
    out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
  } else {
    out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
  }

  // Now compute number of coarse grid points
  for (int dim = 0; dim < numDimensions; ++dim) {
    lCoarseNodesPerDim[dim] = lFineNodesPerDim[dim] / 3 + 1;
  }
  *out << "lCoarseNodesPerDim " << lCoarseNodesPerDim << std::endl;

  // Grab the block size here and multiply all existing offsets by it
  const LO blkSize = A->GetFixedBlockSize();
  *out << "blkSize " << blkSize << std::endl;

  // Based on lCoarseNodesPerDim and lFineNodesPerDim
  // we can compute numRows, numCols and NNZ for R
  const LO numRows = blkSize * lCoarseNodesPerDim[0] * lCoarseNodesPerDim[1] * lCoarseNodesPerDim[2];
  const LO numCols = blkSize * lFineNodesPerDim[0] * lFineNodesPerDim[1] * lFineNodesPerDim[2];

  // Create the coarse coordinates multivector
  // so we can fill it on the fly while computing
  // the restriction operator
  RCP<Map> rowMap = MapFactory::Build(A->getRowMap()->lib(),
                                      Teuchos::OrdinalTraits<GO>::invalid(),
                                      numRows,
                                      A->getRowMap()->getIndexBase(),
                                      A->getRowMap()->getComm());

  RCP<Map> coordRowMap = MapFactory::Build(A->getRowMap()->lib(),
                                           Teuchos::OrdinalTraits<GO>::invalid(),
                                           lCoarseNodesPerDim[0] * lCoarseNodesPerDim[1] * lCoarseNodesPerDim[2],
                                           A->getRowMap()->getIndexBase(),
                                           A->getRowMap()->getComm());

  coarseCoordinates = Xpetra::MultiVectorFactory<real_type, LO, GO, NO>::Build(coordRowMap,
                                                                               numDimensions);

  // Get device views of coordinates
  auto fineCoordsView   = fineCoordinates->getDeviceLocalView(Xpetra::Access::ReadOnly);
  auto coarseCoordsView = coarseCoordinates->getDeviceLocalView(Xpetra::Access::OverwriteAll);

  Array<ArrayRCP<const real_type> > fineCoordData(numDimensions);
  Array<ArrayRCP<real_type> > coarseCoordData(numDimensions);
  for (int dim = 0; dim < numDimensions; ++dim) {
    fineCoordData[dim]   = fineCoordinates->getData(dim);
    coarseCoordData[dim] = coarseCoordinates->getDataNonConst(dim);
  }

  // Let us set some parameter that will be useful
  // while constructing R

  // Length of interpolation stencils based on geometry
  const LO cornerStencilLength   = 27;
  const LO edgeStencilLength     = 45;
  const LO faceStencilLength     = 75;
  const LO interiorStencilLength = 125;

  // Number of corner, edge, face and interior nodes
  const LO numCorners   = 8;
  const LO numEdges     = 4 * (lCoarseNodesPerDim[0] - 2) + 4 * (lCoarseNodesPerDim[1] - 2) + 4 * (lCoarseNodesPerDim[2] - 2);
  const LO numFaces     = 2 * (lCoarseNodesPerDim[0] - 2) * (lCoarseNodesPerDim[1] - 2) + 2 * (lCoarseNodesPerDim[0] - 2) * (lCoarseNodesPerDim[2] - 2) + 2 * (lCoarseNodesPerDim[1] - 2) * (lCoarseNodesPerDim[2] - 2);
  const LO numInteriors = (lCoarseNodesPerDim[0] - 2) * (lCoarseNodesPerDim[1] - 2) * (lCoarseNodesPerDim[2] - 2);

  const LO nnz = (numCorners * cornerStencilLength + numEdges * edgeStencilLength + numFaces * faceStencilLength + numInteriors * interiorStencilLength) * blkSize;

  // Having the number of rows and columns we can genrate
  // the appropriate maps for our operator.

  *out << "R statistics:" << std::endl
       << "  -numRows= " << numRows << std::endl
       << "  -numCols= " << numCols << std::endl
       << "  -nnz=     " << nnz << std::endl;

  row_map_type row_map(Kokkos::ViewAllocateWithoutInitializing("row_map"), numRows + 1);
  typename row_map_type::HostMirror row_map_h = Kokkos::create_mirror_view(row_map);

  entries_type entries(Kokkos::ViewAllocateWithoutInitializing("entries"), nnz);
  typename entries_type::HostMirror entries_h = Kokkos::create_mirror_view(entries);

  values_type values(Kokkos::ViewAllocateWithoutInitializing("values"), nnz);
  typename values_type::HostMirror values_h = Kokkos::create_mirror_view(values);

  // Compute the basic interpolation
  // coefficients for 1D rate of 3
  // coarsening.
  Array<SC> coeffs({1.0 / 3.0, 2.0 / 3.0, 1.0, 2.0 / 3.0, 1.0 / 3.0});
  row_map_h(0) = 0;

  // Define some offsets that
  // will be needed often later on
  const LO edgeLineOffset     = 2 * cornerStencilLength + (lCoarseNodesPerDim[0] - 2) * edgeStencilLength;
  const LO faceLineOffset     = 2 * edgeStencilLength + (lCoarseNodesPerDim[0] - 2) * faceStencilLength;
  const LO interiorLineOffset = 2 * faceStencilLength + (lCoarseNodesPerDim[0] - 2) * interiorStencilLength;

  const LO facePlaneOffset     = 2 * edgeLineOffset + (lCoarseNodesPerDim[1] - 2) * faceLineOffset;
  const LO interiorPlaneOffset = 2 * faceLineOffset + (lCoarseNodesPerDim[1] - 2) * interiorLineOffset;

  // Let us take care of the corners
  // first since we always have
  // corners to deal with!
  {
    // Corner 1
    LO coordRowIdx = 0, rowIdx = 0, coordColumnOffset = 0, columnOffset = 0, entryOffset = 0;
    for (LO l = 0; l < blkSize; ++l) {
      for (LO k = 0; k < 3; ++k) {
        for (LO j = 0; j < 3; ++j) {
          for (LO i = 0; i < 3; ++i) {
            entries_h(entryOffset + k * 9 + j * 3 + i + cornerStencilLength * l) = columnOffset + (k * lFineNodesPerDim[1] * lFineNodesPerDim[0] + j * lFineNodesPerDim[0] + i) * blkSize + l;
            values_h(entryOffset + k * 9 + j * 3 + i + cornerStencilLength * l)  = coeffs[k + 2] * coeffs[j + 2] * coeffs[i + 2];
          }
        }
      }
    }
    for (LO l = 0; l < blkSize; ++l) {
      row_map_h(rowIdx + 1 + l) = entryOffset + cornerStencilLength * (l + 1);
    }
    for (int dim = 0; dim < numDimensions; ++dim) {
      coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
    }

    // Corner 5
    coordRowIdx += (lCoarseNodesPerDim[2] - 1) * lCoarseNodesPerDim[1] * lCoarseNodesPerDim[0];
    rowIdx = coordRowIdx * blkSize;
    coordColumnOffset += (lFineNodesPerDim[2] - 1) * lFineNodesPerDim[1] * lFineNodesPerDim[0];
    columnOffset = coordColumnOffset * blkSize;
    entryOffset += (facePlaneOffset + (lCoarseNodesPerDim[2] - 2) * interiorPlaneOffset) * blkSize;
    for (LO l = 0; l < blkSize; ++l) {
      for (LO k = 0; k < 3; ++k) {
        for (LO j = 0; j < 3; ++j) {
          for (LO i = 0; i < 3; ++i) {
            entries_h(entryOffset + k * 9 + j * 3 + i + cornerStencilLength * l) = columnOffset + ((k - 2) * lFineNodesPerDim[1] * lFineNodesPerDim[0] + j * lFineNodesPerDim[0] + i) * blkSize + l;
            values_h(entryOffset + k * 9 + j * 3 + i + cornerStencilLength * l)  = coeffs[k] * coeffs[j + 2] * coeffs[i + 2];
          }
        }
      }
    }
    for (LO l = 0; l < blkSize; ++l) {
      row_map_h(rowIdx + 1 + l) = entryOffset + cornerStencilLength * (l + 1);
    }
    for (int dim = 0; dim < numDimensions; ++dim) {
      coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
    }

    // Corner 2
    coordRowIdx       = (lCoarseNodesPerDim[0] - 1);
    rowIdx            = coordRowIdx * blkSize;
    coordColumnOffset = (lFineNodesPerDim[0] - 1);
    columnOffset      = coordColumnOffset * blkSize;
    entryOffset       = (cornerStencilLength + (lCoarseNodesPerDim[0] - 2) * edgeStencilLength) * blkSize;
    for (LO l = 0; l < blkSize; ++l) {
      for (LO k = 0; k < 3; ++k) {
        for (LO j = 0; j < 3; ++j) {
          for (LO i = 0; i < 3; ++i) {
            entries_h(entryOffset + k * 9 + j * 3 + i + cornerStencilLength * l) = columnOffset + (k * lFineNodesPerDim[1] * lFineNodesPerDim[0] + j * lFineNodesPerDim[0] + (i - 2)) * blkSize + l;
            values_h(entryOffset + k * 9 + j * 3 + i + cornerStencilLength * l)  = coeffs[k + 2] * coeffs[j + 2] * coeffs[i];
          }
        }
      }
    }
    for (LO l = 0; l < blkSize; ++l) {
      row_map_h(rowIdx + 1 + l) = entryOffset + cornerStencilLength * (l + 1);
    }
    for (int dim = 0; dim < numDimensions; ++dim) {
      coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
    }

    // Corner 6
    coordRowIdx += (lCoarseNodesPerDim[2] - 1) * lCoarseNodesPerDim[1] * lCoarseNodesPerDim[0];
    rowIdx = coordRowIdx * blkSize;
    coordColumnOffset += (lFineNodesPerDim[2] - 1) * lFineNodesPerDim[1] * lFineNodesPerDim[0];
    columnOffset = coordColumnOffset * blkSize;
    entryOffset += (facePlaneOffset + (lCoarseNodesPerDim[2] - 2) * interiorPlaneOffset) * blkSize;
    for (LO l = 0; l < blkSize; ++l) {
      for (LO k = 0; k < 3; ++k) {
        for (LO j = 0; j < 3; ++j) {
          for (LO i = 0; i < 3; ++i) {
            entries_h(entryOffset + k * 9 + j * 3 + i + cornerStencilLength * l) = columnOffset + ((k - 2) * lFineNodesPerDim[1] * lFineNodesPerDim[0] + j * lFineNodesPerDim[0] + i - 2) * blkSize + l;
            values_h(entryOffset + k * 9 + j * 3 + i + cornerStencilLength * l)  = coeffs[k] * coeffs[j + 2] * coeffs[i];
          }
        }
      }
    }
    for (LO l = 0; l < blkSize; ++l) {
      row_map_h(rowIdx + 1 + l) = entryOffset + cornerStencilLength * (l + 1);
    }
    for (int dim = 0; dim < numDimensions; ++dim) {
      coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
    }

    // Corner 3
    coordRowIdx       = (lCoarseNodesPerDim[1] - 1) * lCoarseNodesPerDim[0];
    rowIdx            = coordRowIdx * blkSize;
    coordColumnOffset = (lFineNodesPerDim[1] - 1) * lFineNodesPerDim[0];
    columnOffset      = coordColumnOffset * blkSize;
    entryOffset       = (edgeLineOffset + (lCoarseNodesPerDim[1] - 2) * faceLineOffset) * blkSize;
    for (LO l = 0; l < blkSize; ++l) {
      for (LO k = 0; k < 3; ++k) {
        for (LO j = 0; j < 3; ++j) {
          for (LO i = 0; i < 3; ++i) {
            entries_h(entryOffset + k * 9 + j * 3 + i + cornerStencilLength * l) = columnOffset + (k * lFineNodesPerDim[1] * lFineNodesPerDim[0] + (j - 2) * lFineNodesPerDim[0] + i) * blkSize + l;
            values_h(entryOffset + k * 9 + j * 3 + i + cornerStencilLength * l)  = coeffs[k + 2] * coeffs[j] * coeffs[i + 2];
          }
        }
      }
    }
    for (LO l = 0; l < blkSize; ++l) {
      row_map_h(rowIdx + 1 + l) = entryOffset + cornerStencilLength * (l + 1);
    }
    for (int dim = 0; dim < numDimensions; ++dim) {
      coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
    }

    // Corner 7
    coordRowIdx += (lCoarseNodesPerDim[2] - 1) * lCoarseNodesPerDim[1] * lCoarseNodesPerDim[0];
    rowIdx = coordRowIdx * blkSize;
    coordColumnOffset += (lFineNodesPerDim[2] - 1) * lFineNodesPerDim[1] * lFineNodesPerDim[0];
    columnOffset = coordColumnOffset * blkSize;
    entryOffset += (facePlaneOffset + (lCoarseNodesPerDim[2] - 2) * interiorPlaneOffset) * blkSize;
    for (LO l = 0; l < blkSize; ++l) {
      for (LO k = 0; k < 3; ++k) {
        for (LO j = 0; j < 3; ++j) {
          for (LO i = 0; i < 3; ++i) {
            entries_h(entryOffset + k * 9 + j * 3 + i + cornerStencilLength * l) = columnOffset + ((k - 2) * lFineNodesPerDim[1] * lFineNodesPerDim[0] + (j - 2) * lFineNodesPerDim[0] + i) * blkSize + l;
            values_h(entryOffset + k * 9 + j * 3 + i + cornerStencilLength * l)  = coeffs[k] * coeffs[j] * coeffs[i + 2];
          }
        }
      }
    }
    for (LO l = 0; l < blkSize; ++l) {
      row_map_h(rowIdx + 1 + l) = entryOffset + cornerStencilLength * (l + 1);
    }
    for (int dim = 0; dim < numDimensions; ++dim) {
      coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
    }

    // Corner 4
    coordRowIdx       = (lCoarseNodesPerDim[1] * lCoarseNodesPerDim[0] - 1);
    rowIdx            = coordRowIdx * blkSize;
    coordColumnOffset = (lFineNodesPerDim[1] * lFineNodesPerDim[0] - 1);
    columnOffset      = coordColumnOffset * blkSize;
    entryOffset       = (edgeLineOffset + (lCoarseNodesPerDim[1] - 2) * faceLineOffset +
                   cornerStencilLength + (lCoarseNodesPerDim[0] - 2) * edgeStencilLength) *
                  blkSize;
    for (LO l = 0; l < blkSize; ++l) {
      for (LO k = 0; k < 3; ++k) {
        for (LO j = 0; j < 3; ++j) {
          for (LO i = 0; i < 3; ++i) {
            entries_h(entryOffset + k * 9 + j * 3 + i + cornerStencilLength * l) = columnOffset + (k * lFineNodesPerDim[1] * lFineNodesPerDim[0] + (j - 2) * lFineNodesPerDim[0] + (i - 2)) * blkSize + l;
            values_h(entryOffset + k * 9 + j * 3 + i + cornerStencilLength * l)  = coeffs[k + 2] * coeffs[j] * coeffs[i];
          }
        }
      }
    }
    for (LO l = 0; l < blkSize; ++l) {
      row_map_h(rowIdx + 1 + l) = entryOffset + cornerStencilLength * (l + 1);
    }
    for (int dim = 0; dim < numDimensions; ++dim) {
      coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
    }

    // Corner 8
    coordRowIdx += (lCoarseNodesPerDim[2] - 1) * lCoarseNodesPerDim[1] * lCoarseNodesPerDim[0];
    rowIdx = coordRowIdx * blkSize;
    coordColumnOffset += (lFineNodesPerDim[2] - 1) * lFineNodesPerDim[1] * lFineNodesPerDim[0];
    columnOffset = coordColumnOffset * blkSize;
    entryOffset += (facePlaneOffset + (lCoarseNodesPerDim[2] - 2) * interiorPlaneOffset) * blkSize;
    for (LO l = 0; l < blkSize; ++l) {
      for (LO k = 0; k < 3; ++k) {
        for (LO j = 0; j < 3; ++j) {
          for (LO i = 0; i < 3; ++i) {
            entries_h(entryOffset + k * 9 + j * 3 + i + cornerStencilLength * l) = columnOffset + ((k - 2) * lFineNodesPerDim[1] * lFineNodesPerDim[0] + (j - 2) * lFineNodesPerDim[0] + (i - 2)) * blkSize + l;
            values_h(entryOffset + k * 9 + j * 3 + i + cornerStencilLength * l)  = coeffs[k] * coeffs[j] * coeffs[i];
          }
        }
      }
    }
    for (LO l = 0; l < blkSize; ++l) {
      row_map_h(rowIdx + 1 + l) = entryOffset + cornerStencilLength * (l + 1);
    }
    for (int dim = 0; dim < numDimensions; ++dim) {
      coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
    }
  }  // Corners are done!

  // Edges along 0 direction
  if (lCoarseNodesPerDim[0] - 2 > 0) {
    LO coordRowIdx, rowIdx, coordColumnOffset, columnOffset, entryOffset;
    for (LO edgeIdx = 0; edgeIdx < lCoarseNodesPerDim[0] - 2; ++edgeIdx) {
      // Edge 0
      coordRowIdx       = (edgeIdx + 1);
      rowIdx            = coordRowIdx * blkSize;
      coordColumnOffset = (edgeIdx + 1) * 3;
      columnOffset      = coordColumnOffset * blkSize;
      entryOffset       = (cornerStencilLength + edgeIdx * edgeStencilLength) * blkSize;
      for (LO l = 0; l < blkSize; ++l) {
        for (LO k = 0; k < 3; ++k) {
          for (LO j = 0; j < 3; ++j) {
            for (LO i = 0; i < 5; ++i) {
              entries_h(entryOffset + k * 15 + j * 5 + i + edgeStencilLength * l) = columnOffset + (k * lFineNodesPerDim[1] * lFineNodesPerDim[0] + j * lFineNodesPerDim[0] + i - 2) * blkSize + l;
              values_h(entryOffset + k * 15 + j * 5 + i + edgeStencilLength * l)  = coeffs[k + 2] * coeffs[j + 2] * coeffs[i];
            }
          }
        }
      }
      for (LO l = 0; l < blkSize; ++l) {
        row_map_h(rowIdx + 1 + l) = entryOffset + edgeStencilLength * (l + 1);
      }
      for (int dim = 0; dim < numDimensions; ++dim) {
        coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
      }

      // Edge 1
      coordRowIdx       = ((lCoarseNodesPerDim[1] - 1) * lCoarseNodesPerDim[0] + edgeIdx + 1);
      rowIdx            = coordRowIdx * blkSize;
      coordColumnOffset = ((lFineNodesPerDim[1] - 1) * lFineNodesPerDim[0] + (edgeIdx + 1) * 3);
      columnOffset      = coordColumnOffset * blkSize;
      entryOffset       = (edgeLineOffset + (lCoarseNodesPerDim[1] - 2) * faceLineOffset + cornerStencilLength + edgeIdx * edgeStencilLength) * blkSize;
      for (LO l = 0; l < blkSize; ++l) {
        for (LO k = 0; k < 3; ++k) {
          for (LO j = 0; j < 3; ++j) {
            for (LO i = 0; i < 5; ++i) {
              entries_h(entryOffset + k * 15 + j * 5 + i + edgeStencilLength * l) = columnOffset + (k * lFineNodesPerDim[1] * lFineNodesPerDim[0] + (j - 2) * lFineNodesPerDim[0] + i - 2) * blkSize + l;
              values_h(entryOffset + k * 15 + j * 5 + i + edgeStencilLength * l)  = coeffs[k + 2] * coeffs[j] * coeffs[i];
            }
          }
        }
      }
      for (LO l = 0; l < blkSize; ++l) {
        row_map_h(rowIdx + 1 + l) = entryOffset + edgeStencilLength * (l + 1);
      }
      for (int dim = 0; dim < numDimensions; ++dim) {
        coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
      }

      // Edge 2
      coordRowIdx       = ((lCoarseNodesPerDim[2] - 1) * lCoarseNodesPerDim[1] * lCoarseNodesPerDim[0] + edgeIdx + 1);
      rowIdx            = coordRowIdx * blkSize;
      coordColumnOffset = ((lFineNodesPerDim[2] - 1) * lFineNodesPerDim[1] * lFineNodesPerDim[0] + (edgeIdx + 1) * 3);
      columnOffset      = coordColumnOffset * blkSize;
      entryOffset       = (facePlaneOffset + (lCoarseNodesPerDim[2] - 2) * interiorPlaneOffset + cornerStencilLength + edgeIdx * edgeStencilLength) * blkSize;
      for (LO l = 0; l < blkSize; ++l) {
        for (LO k = 0; k < 3; ++k) {
          for (LO j = 0; j < 3; ++j) {
            for (LO i = 0; i < 5; ++i) {
              entries_h(entryOffset + k * 15 + j * 5 + i + edgeStencilLength * l) = columnOffset + ((k - 2) * lFineNodesPerDim[1] * lFineNodesPerDim[0] + j * lFineNodesPerDim[0] + i - 2) * blkSize + l;
              values_h(entryOffset + k * 15 + j * 5 + i + edgeStencilLength * l)  = coeffs[k] * coeffs[j + 2] * coeffs[i];
            }
          }
        }
      }
      for (LO l = 0; l < blkSize; ++l) {
        row_map_h(rowIdx + 1 + l) = entryOffset + edgeStencilLength * (l + 1);
      }
      for (int dim = 0; dim < numDimensions; ++dim) {
        coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
      }

      // Edge 3
      coordRowIdx       = ((lCoarseNodesPerDim[2] - 1) * lCoarseNodesPerDim[1] * lCoarseNodesPerDim[0] + (lCoarseNodesPerDim[1] - 1) * lCoarseNodesPerDim[0] + edgeIdx + 1);
      rowIdx            = coordRowIdx * blkSize;
      coordColumnOffset = ((lFineNodesPerDim[2] - 1) * lFineNodesPerDim[1] * lFineNodesPerDim[0] + (lFineNodesPerDim[1] - 1) * lFineNodesPerDim[0] + (edgeIdx + 1) * 3);
      columnOffset      = coordColumnOffset * blkSize;
      entryOffset       = (facePlaneOffset + (lCoarseNodesPerDim[2] - 2) * interiorPlaneOffset + edgeLineOffset + (lCoarseNodesPerDim[1] - 2) * faceLineOffset + cornerStencilLength + edgeIdx * edgeStencilLength) * blkSize;
      for (LO l = 0; l < blkSize; ++l) {
        for (LO k = 0; k < 3; ++k) {
          for (LO j = 0; j < 3; ++j) {
            for (LO i = 0; i < 5; ++i) {
              entries_h(entryOffset + k * 15 + j * 5 + i + edgeStencilLength * l) = columnOffset + ((k - 2) * lFineNodesPerDim[1] * lFineNodesPerDim[0] + (j - 2) * lFineNodesPerDim[0] + i - 2) * blkSize + l;
              values_h(entryOffset + k * 15 + j * 5 + i + edgeStencilLength * l)  = coeffs[k] * coeffs[j] * coeffs[i];
            }
          }
        }
      }
      for (LO l = 0; l < blkSize; ++l) {
        row_map_h(rowIdx + 1 + l) = entryOffset + edgeStencilLength * (l + 1);
      }
      for (int dim = 0; dim < numDimensions; ++dim) {
        coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
      }
    }
  }

  // Edges along 1 direction
  if (lCoarseNodesPerDim[1] - 2 > 0) {
    LO coordRowIdx, rowIdx, coordColumnOffset, columnOffset, entryOffset;
    for (LO edgeIdx = 0; edgeIdx < lCoarseNodesPerDim[1] - 2; ++edgeIdx) {
      // Edge 0
      coordRowIdx       = (edgeIdx + 1) * lCoarseNodesPerDim[0];
      rowIdx            = coordRowIdx * blkSize;
      coordColumnOffset = (edgeIdx + 1) * 3 * lFineNodesPerDim[0];
      columnOffset      = coordColumnOffset * blkSize;
      entryOffset       = (edgeLineOffset + edgeIdx * faceLineOffset) * blkSize;
      for (LO l = 0; l < blkSize; ++l) {
        for (LO k = 0; k < 3; ++k) {
          for (LO j = 0; j < 5; ++j) {
            for (LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k * 15 + j * 3 + i + edgeStencilLength * l) = columnOffset + (k * lFineNodesPerDim[1] * lFineNodesPerDim[0] + (j - 2) * lFineNodesPerDim[0] + i) * blkSize + l;
              values_h(entryOffset + k * 15 + j * 3 + i + edgeStencilLength * l)  = coeffs[k + 2] * coeffs[j] * coeffs[i + 2];
            }
          }
        }
      }
      for (LO l = 0; l < blkSize; ++l) {
        row_map_h(rowIdx + 1 + l) = entryOffset + edgeStencilLength * (l + 1);
      }
      for (int dim = 0; dim < numDimensions; ++dim) {
        coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
      }

      // Edge 1
      coordRowIdx       = ((edgeIdx + 1) * lCoarseNodesPerDim[0] + lCoarseNodesPerDim[0] - 1);
      rowIdx            = coordRowIdx * blkSize;
      coordColumnOffset = ((edgeIdx + 1) * 3 * lFineNodesPerDim[0] + lFineNodesPerDim[0] - 1);
      columnOffset      = coordColumnOffset * blkSize;
      entryOffset       = (edgeLineOffset + edgeIdx * faceLineOffset + edgeStencilLength + (lCoarseNodesPerDim[0] - 2) * faceStencilLength) * blkSize;
      for (LO l = 0; l < blkSize; ++l) {
        for (LO k = 0; k < 3; ++k) {
          for (LO j = 0; j < 5; ++j) {
            for (LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k * 15 + j * 3 + i + edgeStencilLength * l) = columnOffset + (k * lFineNodesPerDim[1] * lFineNodesPerDim[0] + (j - 2) * lFineNodesPerDim[0] + i - 2) * blkSize + l;
              values_h(entryOffset + k * 15 + j * 3 + i + edgeStencilLength * l)  = coeffs[k + 2] * coeffs[j] * coeffs[i];
            }
          }
        }
      }
      for (LO l = 0; l < blkSize; ++l) {
        row_map_h(rowIdx + 1 + l) = entryOffset + edgeStencilLength * (l + 1);
      }
      for (int dim = 0; dim < numDimensions; ++dim) {
        coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
      }

      // Edge 2
      coordRowIdx       = ((lCoarseNodesPerDim[2] - 1) * lCoarseNodesPerDim[1] * lCoarseNodesPerDim[0] + (edgeIdx + 1) * lCoarseNodesPerDim[0]);
      rowIdx            = coordRowIdx * blkSize;
      coordColumnOffset = ((lFineNodesPerDim[2] - 1) * lFineNodesPerDim[1] * lFineNodesPerDim[0] + (edgeIdx + 1) * 3 * lFineNodesPerDim[0]);
      columnOffset      = coordColumnOffset * blkSize;
      entryOffset       = (facePlaneOffset + (lCoarseNodesPerDim[2] - 2) * interiorPlaneOffset + edgeLineOffset + edgeIdx * faceLineOffset) * blkSize;
      for (LO l = 0; l < blkSize; ++l) {
        for (LO k = 0; k < 3; ++k) {
          for (LO j = 0; j < 5; ++j) {
            for (LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k * 15 + j * 3 + i + edgeStencilLength * l) = columnOffset + ((k - 2) * lFineNodesPerDim[1] * lFineNodesPerDim[0] + (j - 2) * lFineNodesPerDim[0] + i) * blkSize + l;
              values_h(entryOffset + k * 15 + j * 3 + i + edgeStencilLength * l)  = coeffs[k] * coeffs[j] * coeffs[i + 2];
            }
          }
        }
      }
      for (LO l = 0; l < blkSize; ++l) {
        row_map_h(rowIdx + 1 + l) = entryOffset + edgeStencilLength * (l + 1);
      }
      for (int dim = 0; dim < numDimensions; ++dim) {
        coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
      }

      // Edge 3
      coordRowIdx       = ((lCoarseNodesPerDim[2] - 1) * lCoarseNodesPerDim[1] * lCoarseNodesPerDim[0] + (edgeIdx + 1) * lCoarseNodesPerDim[0] + lCoarseNodesPerDim[0] - 1);
      rowIdx            = coordRowIdx * blkSize;
      coordColumnOffset = ((lFineNodesPerDim[2] - 1) * lFineNodesPerDim[1] * lFineNodesPerDim[0] + (edgeIdx + 1) * 3 * lFineNodesPerDim[0] + lFineNodesPerDim[0] - 1);
      columnOffset      = coordColumnOffset * blkSize;
      entryOffset       = (facePlaneOffset + (lCoarseNodesPerDim[2] - 2) * interiorPlaneOffset + edgeLineOffset + edgeIdx * faceLineOffset + edgeStencilLength + (lCoarseNodesPerDim[0] - 2) * faceStencilLength) * blkSize;
      for (LO l = 0; l < blkSize; ++l) {
        for (LO k = 0; k < 3; ++k) {
          for (LO j = 0; j < 5; ++j) {
            for (LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k * 15 + j * 3 + i + edgeStencilLength * l) = columnOffset + ((k - 2) * lFineNodesPerDim[1] * lFineNodesPerDim[0] + (j - 2) * lFineNodesPerDim[0] + i - 2) * blkSize + l;
              values_h(entryOffset + k * 15 + j * 3 + i + edgeStencilLength * l)  = coeffs[k] * coeffs[j] * coeffs[i];
            }
          }
        }
      }
      for (LO l = 0; l < blkSize; ++l) {
        row_map_h(rowIdx + 1 + l) = entryOffset + edgeStencilLength * (l + 1);
      }
      for (int dim = 0; dim < numDimensions; ++dim) {
        coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
      }
    }
  }

  // Edges along 2 direction
  if (lCoarseNodesPerDim[2] - 2 > 0) {
    LO coordRowIdx, rowIdx, coordColumnOffset, columnOffset, entryOffset;
    for (LO edgeIdx = 0; edgeIdx < lCoarseNodesPerDim[2] - 2; ++edgeIdx) {
      // Edge 0
      coordRowIdx       = (edgeIdx + 1) * lCoarseNodesPerDim[1] * lCoarseNodesPerDim[0];
      rowIdx            = coordRowIdx * blkSize;
      coordColumnOffset = (edgeIdx + 1) * 3 * lFineNodesPerDim[1] * lFineNodesPerDim[0];
      columnOffset      = coordColumnOffset * blkSize;
      entryOffset       = (facePlaneOffset + edgeIdx * interiorPlaneOffset) * blkSize;
      for (LO l = 0; l < blkSize; ++l) {
        for (LO k = 0; k < 5; ++k) {
          for (LO j = 0; j < 3; ++j) {
            for (LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k * 9 + j * 3 + i + edgeStencilLength * l) = columnOffset + ((k - 2) * lFineNodesPerDim[1] * lFineNodesPerDim[0] + j * lFineNodesPerDim[0] + i) * blkSize + l;
              values_h(entryOffset + k * 9 + j * 3 + i + edgeStencilLength * l)  = coeffs[k] * coeffs[j + 2] * coeffs[i + 2];
            }
          }
        }
      }
      for (LO l = 0; l < blkSize; ++l) {
        row_map_h(rowIdx + 1 + l) = entryOffset + edgeStencilLength * (l + 1);
      }
      for (int dim = 0; dim < numDimensions; ++dim) {
        coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
      }

      // Edge 1
      coordRowIdx       = ((edgeIdx + 1) * lCoarseNodesPerDim[1] * lCoarseNodesPerDim[0] + lCoarseNodesPerDim[0] - 1);
      rowIdx            = coordRowIdx * blkSize;
      coordColumnOffset = ((edgeIdx + 1) * 3 * lFineNodesPerDim[1] * lFineNodesPerDim[0] + lFineNodesPerDim[0] - 1);
      columnOffset      = coordColumnOffset * blkSize;
      entryOffset       = (facePlaneOffset + faceLineOffset - edgeStencilLength + edgeIdx * interiorPlaneOffset) * blkSize;
      for (LO l = 0; l < blkSize; ++l) {
        for (LO k = 0; k < 5; ++k) {
          for (LO j = 0; j < 3; ++j) {
            for (LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k * 9 + j * 3 + i + edgeStencilLength * l) = columnOffset + ((k - 2) * lFineNodesPerDim[1] * lFineNodesPerDim[0] + j * lFineNodesPerDim[0] + i - 2) * blkSize + l;
              values_h(entryOffset + k * 9 + j * 3 + i + edgeStencilLength * l)  = coeffs[k] * coeffs[j + 2] * coeffs[i];
            }
          }
        }
      }
      for (LO l = 0; l < blkSize; ++l) {
        row_map_h(rowIdx + 1 + l) = entryOffset + edgeStencilLength * (l + 1);
      }
      for (int dim = 0; dim < numDimensions; ++dim) {
        coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
      }

      // Edge 2
      coordRowIdx       = ((edgeIdx + 1) * lCoarseNodesPerDim[1] * lCoarseNodesPerDim[0] + (lCoarseNodesPerDim[1] - 1) * lCoarseNodesPerDim[0]);
      rowIdx            = coordRowIdx * blkSize;
      coordColumnOffset = ((edgeIdx + 1) * 3 * lFineNodesPerDim[1] * lFineNodesPerDim[0] + (lFineNodesPerDim[1] - 1) * lFineNodesPerDim[0]);
      columnOffset      = coordColumnOffset * blkSize;
      entryOffset       = (facePlaneOffset + edgeIdx * interiorPlaneOffset + faceLineOffset + (lCoarseNodesPerDim[1] - 2) * interiorLineOffset) * blkSize;
      for (LO l = 0; l < blkSize; ++l) {
        for (LO k = 0; k < 5; ++k) {
          for (LO j = 0; j < 3; ++j) {
            for (LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k * 9 + j * 3 + i + edgeStencilLength * l) = columnOffset + ((k - 2) * lFineNodesPerDim[1] * lFineNodesPerDim[0] + (j - 2) * lFineNodesPerDim[0] + i) * blkSize + l;
              values_h(entryOffset + k * 9 + j * 3 + i + edgeStencilLength * l)  = coeffs[k] * coeffs[j] * coeffs[i + 2];
            }
          }
        }
      }
      for (LO l = 0; l < blkSize; ++l) {
        row_map_h(rowIdx + 1 + l) = entryOffset + edgeStencilLength * (l + 1);
      }
      for (int dim = 0; dim < numDimensions; ++dim) {
        coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
      }

      // Edge 3
      coordRowIdx       = ((edgeIdx + 2) * lCoarseNodesPerDim[1] * lCoarseNodesPerDim[0] - 1);
      rowIdx            = coordRowIdx * blkSize;
      coordColumnOffset = ((edgeIdx + 1) * 3 * lFineNodesPerDim[1] * lFineNodesPerDim[0] + lFineNodesPerDim[1] * lFineNodesPerDim[0] - 1);
      columnOffset      = coordColumnOffset * blkSize;
      entryOffset       = (facePlaneOffset + (edgeIdx + 1) * interiorPlaneOffset - edgeStencilLength) * blkSize;
      for (LO l = 0; l < blkSize; ++l) {
        for (LO k = 0; k < 5; ++k) {
          for (LO j = 0; j < 3; ++j) {
            for (LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k * 9 + j * 3 + i + edgeStencilLength * l) = columnOffset + ((k - 2) * lFineNodesPerDim[1] * lFineNodesPerDim[0] + (j - 2) * lFineNodesPerDim[0] + i - 2) * blkSize + l;
              values_h(entryOffset + k * 9 + j * 3 + i + edgeStencilLength * l)  = coeffs[k] * coeffs[j] * coeffs[i];
            }
          }
        }
      }
      for (LO l = 0; l < blkSize; ++l) {
        row_map_h(rowIdx + 1 + l) = entryOffset + edgeStencilLength * (l + 1);
      }
      for (int dim = 0; dim < numDimensions; ++dim) {
        coarseCoordData[dim][coordRowIdx] = fineCoordData[dim][coordColumnOffset];
      }
    }
  }

  // TODO: KOKKOS parallel_for used from here. Not sure if it should be used for edges.
  Kokkos::deep_copy(row_map, row_map_h);
  Kokkos::deep_copy(entries, entries_h);
  Kokkos::deep_copy(values, values_h);

  // Create views on device for nodes per dim
  LOTupleView lFineNodesPerDim_d("lFineNodesPerDim");
  LOTupleView lCoarseNodesPerDim_d("lCoarseNodesPerDim");

  typename Kokkos::View<LO[3], device_type>::HostMirror lCoarseNodesPerDim_h = Kokkos::create_mirror_view(lCoarseNodesPerDim_d);
  typename Kokkos::View<LO[3], device_type>::HostMirror lFineNodesPerDim_h   = Kokkos::create_mirror_view(lFineNodesPerDim_d);

  for (int dim = 0; dim < numDimensions; ++dim) {
    lCoarseNodesPerDim_h(dim) = lCoarseNodesPerDim[dim];
    lFineNodesPerDim_h(dim)   = lFineNodesPerDim[dim];
  }

  Kokkos::deep_copy(lCoarseNodesPerDim_d, lCoarseNodesPerDim_h);
  Kokkos::deep_copy(lFineNodesPerDim_d, lFineNodesPerDim_h);

  // Faces in 0-1 plane
  if ((lCoarseNodesPerDim[0] - 2 > 0) && (lCoarseNodesPerDim[1] - 2 > 0)) {
    Kokkos::parallel_for(
        "Faces in 0-1 plane region R",
        Kokkos::RangePolicy<typename device_type::execution_space>(0, (lCoarseNodesPerDim[1] - 2) * (lCoarseNodesPerDim[0] - 2)),
        KOKKOS_LAMBDA(const LO faceIdx) {
          LO coordRowIdx, rowIdx, coordColumnOffset, columnOffset, entryOffset;
          LO gridIdx[3]                = {0, 0, 0};
          impl_scalar_type coeffs_d[5] = {1.0 / 3.0, 2.0 / 3.0, 1.0, 2.0 / 3.0, 1.0 / 3.0};
          // Last step in the loop
          // update the grid indices
          // for next grid point
          for (LO i = 0; i < faceIdx; i++) {
            ++gridIdx[0];
            if (gridIdx[0] == lCoarseNodesPerDim_d(0) - 2) {
              gridIdx[0] = 0;
              ++gridIdx[1];
            }
          }

          // Face 0
          coordRowIdx       = ((gridIdx[1] + 1) * lCoarseNodesPerDim_d(0) + gridIdx[0] + 1);
          rowIdx            = coordRowIdx * blkSize;
          coordColumnOffset = 3 * ((gridIdx[1] + 1) * lFineNodesPerDim_d(0) + gridIdx[0] + 1);
          columnOffset      = coordColumnOffset * blkSize;
          entryOffset       = (edgeLineOffset + edgeStencilLength + gridIdx[1] * faceLineOffset + gridIdx[0] * faceStencilLength) * blkSize;
          for (LO l = 0; l < blkSize; ++l) {
            for (LO k = 0; k < 3; ++k) {
              for (LO j = 0; j < 5; ++j) {
                for (LO i = 0; i < 5; ++i) {
                  entries(entryOffset + k * 25 + j * 5 + i + faceStencilLength * l) = columnOffset + (k * lFineNodesPerDim_d(1) * lFineNodesPerDim_d(0) + (j - 2) * lFineNodesPerDim_d(0) + i - 2) * blkSize + l;
                  values(entryOffset + k * 25 + j * 5 + i + faceStencilLength * l)  = coeffs_d[k + 2] * coeffs_d[j] * coeffs_d[i];
                }
              }
            }
          }
          for (LO l = 0; l < blkSize; ++l) {
            row_map(rowIdx + 1 + l) = entryOffset + faceStencilLength * (l + 1);
          }
          for (int dim = 0; dim < numDimensions; ++dim) {
            coarseCoordsView(coordRowIdx, dim) = fineCoordsView(coordColumnOffset, dim);
          }

          // Face 1
          coordRowIdx += (lCoarseNodesPerDim_d(2) - 1) * lCoarseNodesPerDim_d(1) * lCoarseNodesPerDim_d(0);
          rowIdx = coordRowIdx * blkSize;
          coordColumnOffset += (lFineNodesPerDim_d(2) - 1) * lFineNodesPerDim_d(1) * lFineNodesPerDim_d(0);
          columnOffset = coordColumnOffset * blkSize;
          entryOffset += (facePlaneOffset + (lCoarseNodesPerDim_d(2) - 2) * interiorPlaneOffset) * blkSize;
          for (LO l = 0; l < blkSize; ++l) {
            for (LO k = 0; k < 3; ++k) {
              for (LO j = 0; j < 5; ++j) {
                for (LO i = 0; i < 5; ++i) {
                  entries(entryOffset + k * 25 + j * 5 + i + faceStencilLength * l) = columnOffset + ((k - 2) * lFineNodesPerDim_d(1) * lFineNodesPerDim_d(0) + (j - 2) * lFineNodesPerDim_d(0) + i - 2) * blkSize + l;
                  values(entryOffset + k * 25 + j * 5 + i + faceStencilLength * l)  = coeffs_d[k] * coeffs_d[j] * coeffs_d[i];
                }
              }
            }
          }
          for (LO l = 0; l < blkSize; ++l) {
            row_map(rowIdx + 1 + l) = entryOffset + faceStencilLength * (l + 1);
          }
          for (int dim = 0; dim < numDimensions; ++dim) {
            coarseCoordsView(coordRowIdx, dim) = fineCoordsView(coordColumnOffset, dim);
          }
        });  // parallel_for faces in 0-1 plane
  }

  // Faces in 0-2 plane
  if ((lCoarseNodesPerDim[0] - 2 > 0) && (lCoarseNodesPerDim[2] - 2 > 0)) {
    Kokkos::parallel_for(
        "Faces in 0-2 plane region R",
        Kokkos::RangePolicy<typename device_type::execution_space>(0, (lCoarseNodesPerDim[2] - 2) * (lCoarseNodesPerDim[0] - 2)),
        KOKKOS_LAMBDA(const LO faceIdx) {
          LO coordRowIdx, rowIdx, coordColumnOffset, columnOffset, entryOffset;
          LO gridIdx[3]                = {0, 0, 0};
          impl_scalar_type coeffs_d[5] = {1.0 / 3.0, 2.0 / 3.0, 1.0, 2.0 / 3.0, 1.0 / 3.0};
          // Last step in the loop
          // update the grid indices
          // for next grid point
          for (LO i = 0; i < faceIdx; i++) {
            ++gridIdx[0];
            if (gridIdx[0] == lCoarseNodesPerDim_d(0) - 2) {
              gridIdx[0] = 0;
              ++gridIdx[2];
            }
          }

          // Face 0
          coordRowIdx       = ((gridIdx[2] + 1) * lCoarseNodesPerDim_d(1) * lCoarseNodesPerDim_d(0) + (gridIdx[0] + 1));
          rowIdx            = coordRowIdx * blkSize;
          coordColumnOffset = ((gridIdx[2] + 1) * lFineNodesPerDim_d(1) * lFineNodesPerDim_d(0) + gridIdx[0] + 1) * 3;
          columnOffset      = coordColumnOffset * blkSize;
          entryOffset       = (facePlaneOffset + gridIdx[2] * interiorPlaneOffset + edgeStencilLength + gridIdx[0] * faceStencilLength) * blkSize;
          for (LO l = 0; l < blkSize; ++l) {
            for (LO k = 0; k < 5; ++k) {
              for (LO j = 0; j < 3; ++j) {
                for (LO i = 0; i < 5; ++i) {
                  entries(entryOffset + k * 15 + j * 5 + i + faceStencilLength * l) = columnOffset + ((k - 2) * lFineNodesPerDim_d(1) * lFineNodesPerDim_d(0) + j * lFineNodesPerDim_d(0) + i - 2) * blkSize + l;
                  values(entryOffset + k * 15 + j * 5 + i + faceStencilLength * l)  = coeffs_d[k] * coeffs_d[j + 2] * coeffs_d[i];
                }
              }
            }
          }
          for (LO l = 0; l < blkSize; ++l) {
            row_map(rowIdx + 1 + l) = entryOffset + faceStencilLength * (l + 1);
          }
          for (int dim = 0; dim < numDimensions; ++dim) {
            coarseCoordsView(coordRowIdx, dim) = fineCoordsView(coordColumnOffset, dim);
          }

          // Face 1
          coordRowIdx += (lCoarseNodesPerDim_d(1) - 1) * lCoarseNodesPerDim_d(0);
          rowIdx = coordRowIdx * blkSize;
          coordColumnOffset += (lFineNodesPerDim_d(1) - 1) * lFineNodesPerDim_d(0);
          columnOffset = coordColumnOffset * blkSize;
          entryOffset += (faceLineOffset + (lCoarseNodesPerDim_d(1) - 2) * interiorLineOffset) * blkSize;
          for (LO l = 0; l < blkSize; ++l) {
            for (LO k = 0; k < 5; ++k) {
              for (LO j = 0; j < 3; ++j) {
                for (LO i = 0; i < 5; ++i) {
                  entries(entryOffset + k * 15 + j * 5 + i + faceStencilLength * l) = columnOffset + ((k - 2) * lFineNodesPerDim_d(1) * lFineNodesPerDim_d(0) + (j - 2) * lFineNodesPerDim_d(0) + i - 2) * blkSize + l;
                  values(entryOffset + k * 15 + j * 5 + i + faceStencilLength * l)  = coeffs_d[k] * coeffs_d[j] * coeffs_d[i];
                }
              }
            }
          }
          for (LO l = 0; l < blkSize; ++l) {
            row_map(rowIdx + 1 + l) = entryOffset + faceStencilLength * (l + 1);
          }
          for (int dim = 0; dim < numDimensions; ++dim) {
            coarseCoordsView(coordRowIdx, dim) = fineCoordsView(coordColumnOffset, dim);
          }
        });  // parallel_for faces in 0-2 plane
  }

  // Faces in 1-2 plane
  if ((lCoarseNodesPerDim[1] - 2 > 0) && (lCoarseNodesPerDim[2] - 2 > 0)) {
    Kokkos::parallel_for(
        "Faces in 1-2 plane region R",
        Kokkos::RangePolicy<typename device_type::execution_space>(0, (lCoarseNodesPerDim[2] - 2) * (lCoarseNodesPerDim[1] - 2)),
        KOKKOS_LAMBDA(const LO faceIdx) {
          LO coordRowIdx, rowIdx, coordColumnOffset, columnOffset, entryOffset;
          LO gridIdx[3]                = {0, 0, 0};
          impl_scalar_type coeffs_d[5] = {1.0 / 3.0, 2.0 / 3.0, 1.0, 2.0 / 3.0, 1.0 / 3.0};
          // Last step in the loop
          // update the grid indices
          // for next grid point
          for (LO i = 0; i < faceIdx; i++) {
            ++gridIdx[1];
            if (gridIdx[1] == lCoarseNodesPerDim_d(1) - 2) {
              gridIdx[1] = 0;
              ++gridIdx[2];
            }
          }

          // Face 0
          coordRowIdx       = ((gridIdx[2] + 1) * lCoarseNodesPerDim_d(1) * lCoarseNodesPerDim_d(0) + (gridIdx[1] + 1) * lCoarseNodesPerDim_d(0));
          rowIdx            = coordRowIdx * blkSize;
          coordColumnOffset = ((gridIdx[2] + 1) * lFineNodesPerDim_d(1) * lFineNodesPerDim_d(0) + (gridIdx[1] + 1) * lFineNodesPerDim_d(0)) * 3;
          columnOffset      = coordColumnOffset * blkSize;
          entryOffset       = (facePlaneOffset + gridIdx[2] * interiorPlaneOffset + faceLineOffset + gridIdx[1] * interiorLineOffset) * blkSize;
          for (LO l = 0; l < blkSize; ++l) {
            for (LO k = 0; k < 5; ++k) {
              for (LO j = 0; j < 5; ++j) {
                for (LO i = 0; i < 3; ++i) {
                  entries(entryOffset + k * 15 + j * 3 + i + faceStencilLength * l) = columnOffset + ((k - 2) * lFineNodesPerDim_d(1) * lFineNodesPerDim_d(0) + (j - 2) * lFineNodesPerDim_d(0) + i) * blkSize + l;
                  values(entryOffset + k * 15 + j * 3 + i + faceStencilLength * l)  = coeffs_d[k] * coeffs_d[j] * coeffs_d[i + 2];
                }
              }
            }
          }
          for (LO l = 0; l < blkSize; ++l) {
            row_map(rowIdx + 1 + l) = entryOffset + faceStencilLength * (l + 1);
          }
          for (int dim = 0; dim < numDimensions; ++dim) {
            coarseCoordsView(coordRowIdx, dim) = fineCoordsView(coordColumnOffset, dim);
          }

          // Face 1
          coordRowIdx += (lCoarseNodesPerDim_d(0) - 1);
          rowIdx = coordRowIdx * blkSize;
          coordColumnOffset += (lFineNodesPerDim_d(0) - 1);
          columnOffset = coordColumnOffset * blkSize;
          entryOffset += (faceStencilLength + (lCoarseNodesPerDim_d(0) - 2) * interiorStencilLength) * blkSize;
          for (LO l = 0; l < blkSize; ++l) {
            for (LO k = 0; k < 5; ++k) {
              for (LO j = 0; j < 5; ++j) {
                for (LO i = 0; i < 3; ++i) {
                  entries(entryOffset + k * 15 + j * 3 + i + faceStencilLength * l) = columnOffset + ((k - 2) * lFineNodesPerDim_d(1) * lFineNodesPerDim_d(0) + (j - 2) * lFineNodesPerDim_d(0) + i - 2) * blkSize + l;
                  values(entryOffset + k * 15 + j * 3 + i + faceStencilLength * l)  = coeffs_d[k] * coeffs_d[j] * coeffs_d[i];
                }
              }
            }
          }
          for (LO l = 0; l < blkSize; ++l) {
            row_map(rowIdx + 1 + l) = entryOffset + faceStencilLength * (l + 1);
          }
          for (int dim = 0; dim < numDimensions; ++dim) {
            coarseCoordsView(coordRowIdx, dim) = fineCoordsView(coordColumnOffset, dim);
          }
        });  // parallel_for faces in 1-2 plane
  }

  if (numInteriors > 0) {
    // Allocate and compute arrays
    // containing column offsets
    // and values associated with
    // interior points
    LO countRowEntries = 0;
    Kokkos::View<LO[125]> coordColumnOffsets_d("coordColOffset");
    auto coordColumnOffsets_h = Kokkos::create_mirror_view(coordColumnOffsets_d);

    for (LO k = -2; k < 3; ++k) {
      for (LO j = -2; j < 3; ++j) {
        for (LO i = -2; i < 3; ++i) {
          coordColumnOffsets_h(countRowEntries) = k * lFineNodesPerDim[1] * lFineNodesPerDim[0] + j * lFineNodesPerDim[0] + i;
          ++countRowEntries;
        }
      }
    }
    Kokkos::deep_copy(coordColumnOffsets_d, coordColumnOffsets_h);

    LO countValues = 0;
    Kokkos::View<impl_scalar_type*> interiorValues_d("interiorValues", 125);
    auto interiorValues_h = Kokkos::create_mirror_view(interiorValues_d);
    for (LO k = 0; k < 5; ++k) {
      for (LO j = 0; j < 5; ++j) {
        for (LO i = 0; i < 5; ++i) {
          interiorValues_h(countValues) = coeffs[k] * coeffs[j] * coeffs[i];
          ++countValues;
        }
      }
    }
    Kokkos::deep_copy(interiorValues_d, interiorValues_h);

    Kokkos::parallel_for(
        "interior idx region R", Kokkos::RangePolicy<typename device_type::execution_space>(0, numInteriors),
        KOKKOS_LAMBDA(const LO interiorIdx) {
          LO coordRowIdx, rowIdx, coordColumnOffset, columnOffset, entryOffset;
          LO gridIdx[3];
          gridIdx[0] = 0;
          gridIdx[1] = 0;
          gridIdx[2] = 0;
          // First step in the loop
          // update the grid indices
          // for the grid point
          for (LO i = 0; i < interiorIdx; i++) {
            ++gridIdx[0];
            if (gridIdx[0] == lCoarseNodesPerDim_d(0) - 2) {
              gridIdx[0] = 0;
              ++gridIdx[1];
              if (gridIdx[1] == lCoarseNodesPerDim_d(1) - 2) {
                gridIdx[1] = 0;
                ++gridIdx[2];
              }
            }
          }

          coordRowIdx       = ((gridIdx[2] + 1) * lCoarseNodesPerDim_d(0) * lCoarseNodesPerDim_d(1) + (gridIdx[1] + 1) * lCoarseNodesPerDim_d(0) + gridIdx[0] + 1);
          rowIdx            = coordRowIdx * blkSize;
          coordColumnOffset = ((gridIdx[2] + 1) * 3 * lFineNodesPerDim_d(1) * lFineNodesPerDim_d(0) + (gridIdx[1] + 1) * 3 * lFineNodesPerDim_d(0) + (gridIdx[0] + 1) * 3);
          columnOffset      = coordColumnOffset * blkSize;
          entryOffset       = (facePlaneOffset + faceLineOffset + faceStencilLength + gridIdx[2] * interiorPlaneOffset + gridIdx[1] * interiorLineOffset + gridIdx[0] * interiorStencilLength) * blkSize;
          for (LO l = 0; l < blkSize; ++l) {
            row_map(rowIdx + 1 + l) = entryOffset + interiorStencilLength * (l + 1);
          }
          // Fill the column indices
          // and values in the approproate
          // views.
          for (LO l = 0; l < blkSize; ++l) {
            for (LO entryIdx = 0; entryIdx < interiorStencilLength; ++entryIdx) {
              entries(entryOffset + entryIdx + interiorStencilLength * l) = columnOffset + coordColumnOffsets_d(entryIdx) * blkSize + l;
              values(entryOffset + entryIdx + interiorStencilLength * l)  = interiorValues_d(entryIdx);
            }
          }
          for (int dim = 0; dim < numDimensions; ++dim) {
            coarseCoordsView(coordRowIdx, dim) = fineCoordsView(coordColumnOffset, dim);
          }
        });  // Kokkos::parallel_for interior idx
    //
  }

  local_graph_type localGraph(entries, row_map);
  local_matrix_type localR("R", numCols, values, localGraph);

  R = MatrixFactory::Build(localR,          // the local data
                           rowMap,          // rowMap
                           A->getRowMap(),  // colMap
                           A->getRowMap(),  // domainMap == colMap
                           rowMap,          // rangeMap  == rowMap
                           Teuchos::null);  // params for optimized construction

}  // Build3D()

}  // namespace MueLu

#define MUELU_REGIONRFACTORY_KOKKOS_SHORT
#endif  // MUELU_REGIONRFACTORY_DEF_HPP
