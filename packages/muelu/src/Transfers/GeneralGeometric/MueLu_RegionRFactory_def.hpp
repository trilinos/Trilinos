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
#ifndef MUELU_REGIONRFACTORY_DEF_HPP
#define MUELU_REGIONRFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsGraphFactory.hpp>

// #include "MueLu_PFactory.hpp"
// #include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_MasterList.hpp"

#include "MueLu_RegionRFactory_decl.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> RegionRFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set<RCP<const FactoryBase> >("A",                            Teuchos::null,
                                                 "Generating factory of the matrix A");
    validParamList->set<RCP<const FactoryBase> >("numDimensions",                Teuchos::null,
                                                 "Number of spatial dimensions in the problem.");
    validParamList->set<RCP<const FactoryBase> >("lNodesPerDim",                 Teuchos::null,
                                                 "Number of local nodes per spatial dimension on the fine grid.");
    validParamList->set<RCP<const FactoryBase> >("Nullspace",                    Teuchos::null,
                                                 "Fine level nullspace used to construct the coarse level nullspace.");
    validParamList->set<RCP<const FactoryBase> >("Coordinates",                  Teuchos::null,
                                                 "Fine level coordinates used to construct piece-wise linear prolongator and coarse level coordinates.");
    validParamList->set<bool>                   ("keep coarse coords",           false, "Flag to keep coordinates for special coarse grid solve");

    return validParamList;
  } // GetValidParameterList()

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void RegionRFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  DeclareInput(Level& fineLevel, Level& /* coarseLevel */) const {

    Input(fineLevel, "A");
    Input(fineLevel, "numDimensions");
    Input(fineLevel, "lNodesPerDim");
    Input(fineLevel, "Nullspace");
    Input(fineLevel, "Coordinates");

  } // DeclareInput()

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void RegionRFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  Build(Level& fineLevel, Level& coarseLevel) const {

    // Set debug outputs based on environment variable
    RCP<Teuchos::FancyOStream> out;
    if(const char* dbg = std::getenv("MUELU_REGIONRFACTORY_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setShowAllFrontMatter(false).setShowProcRank(true);
    } else {
      out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
    }

    *out << "Starting RegionRFactory::Build." << std::endl;

    // First get the inputs from the fineLevel
    const int numDimensions = Get<int>(fineLevel, "numDimensions");
    Array<LO> lFineNodesPerDim(3, Teuchos::OrdinalTraits<LO>::one());
    {
      Array<LO> lNodesPerDim = Get<Array<LO> >(fineLevel, "lNodesPerDim");
      for(int dim = 0; dim < numDimensions; ++dim) {
        lFineNodesPerDim[dim] = lNodesPerDim[dim];
      }
    }
    *out << "numDimensions " << numDimensions << " and lFineNodesPerDim: " << lFineNodesPerDim
         << std::endl;

    // Let us check that the inputs verify our assumptions
    if(numDimensions < 1 || numDimensions > 3) {
      throw std::runtime_error("numDimensions must be 1, 2 or 3!");
    }
    for(int dim = 0; dim < numDimensions; ++dim) {
      if(lFineNodesPerDim[dim] % 3 != 1) {
        throw std::runtime_error("The number of fine node in each direction need to be 3n+1");
      }
    }
    Array<LO> lCoarseNodesPerDim(3, Teuchos::OrdinalTraits<LO>::one());

    const RCP<Matrix> A = Get<RCP<Matrix> >(fineLevel, "A");

    RCP<realvaluedmultivector_type> fineCoordinates, coarseCoordinates;
    fineCoordinates = Get< RCP<realvaluedmultivector_type> >(fineLevel, "Coordinates");
    if(static_cast<int>(fineCoordinates->getNumVectors()) != numDimensions) {
      throw std::runtime_error("The number of vectors in the coordinates is not equal to numDimensions!");
    }

    // Let us create R and pass it down to the
    // appropriate specialization and see what we
    // get back!
    RCP<Matrix> R;

    if(numDimensions == 1) {
      throw std::runtime_error("RegionRFactory no implemented for 1D case yet.");
    } else if(numDimensions == 2) {
      throw std::runtime_error("RegionRFactory no implemented for 2D case yet.");
    } else if(numDimensions == 3) {
      Build3D(numDimensions, lFineNodesPerDim, A, fineCoordinates,
              R, coarseCoordinates, lCoarseNodesPerDim);
    }

    const Teuchos::ParameterList& pL = GetParameterList();

    // Reuse pattern if available (multiple solve)
    RCP<ParameterList> Tparams;
    if(pL.isSublist("matrixmatrix: kernel params"))
        Tparams=rcp(new ParameterList(pL.sublist("matrixmatrix: kernel params")));
      else
        Tparams= rcp(new ParameterList);

    // R->describe(*out, Teuchos::VERB_EXTREME);
    *out << "Compute P=R^t" << std::endl;
    // By default, we don't need global constants for transpose
    Tparams->set("compute global constants: temporaries",Tparams->get("compute global constants: temporaries", false));
    Tparams->set("compute global constants", Tparams->get("compute global constants",false));
    std::string label = "MueLu::RegionR-transR" + Teuchos::toString(coarseLevel.GetLevelID());
    RCP<Matrix> P = Utilities::Transpose(*R, true, label, Tparams);

    *out << "Compute coarse nullspace" << std::endl;
    RCP<MultiVector> fineNullspace   = Get<RCP<MultiVector> >(fineLevel, "Nullspace");
    RCP<MultiVector> coarseNullspace = MultiVectorFactory::Build(R->getRowMap(),
                                                                 fineNullspace->getNumVectors());
    R->apply(*fineNullspace, *coarseNullspace, Teuchos::NO_TRANS, Teuchos::ScalarTraits<SC>::one(),
             Teuchos::ScalarTraits<SC>::zero());

    *out << "Set data on coarse level" << std::endl;
    Set(coarseLevel, "numDimensions", numDimensions);
    Set(coarseLevel, "lNodesPerDim",  lCoarseNodesPerDim);
    Set(coarseLevel, "Nullspace",     coarseNullspace);
    Set(coarseLevel, "Coordinates",   coarseCoordinates);
    if(pL.get<bool>("keep coarse coords")) {
      coarseLevel.Set<RCP<realvaluedmultivector_type> >("Coordinates2", coarseCoordinates, NoFactory::get());
    }
    Set(coarseLevel, "R", R);
    Set(coarseLevel, "P", P);

  } // Build()

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void RegionRFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
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

    // Set debug outputs based on environment variable
    RCP<Teuchos::FancyOStream> out;
    if(const char* dbg = std::getenv("MUELU_REGIONRFACTORY_DEBUG")) {
      out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setShowAllFrontMatter(false).setShowProcRank(true);
    } else {
      out = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
    }


    // Now compute number of coarse grid points
    for(int dim = 0; dim < 3; ++dim) {
      lCoarseNodesPerDim[dim] = lFineNodesPerDim[dim] / 3 + 1;
    }
    *out << "lCoarseNodesPerDim " << lCoarseNodesPerDim << std::endl;

    // Based on lCoarseNodesPerDim and lFineNodesPerDim
    // we can compute numRows, numCols and NNZ for R
    const LO numRows = lCoarseNodesPerDim[0]*lCoarseNodesPerDim[1]*lCoarseNodesPerDim[2];
    const LO numCols = lFineNodesPerDim[0]*lFineNodesPerDim[1]*lFineNodesPerDim[2];

    // Create the coarse coordinates multivector
    // so we can fill it on the fly while computing
    // the restriction operator
    RCP<Map> rowMap = MapFactory::Build(A->getRowMap()->lib(),
                                        Teuchos::OrdinalTraits<GO>::invalid(),
                                        numRows,
                                        A->getRowMap()->getIndexBase(),
                                        A->getRowMap()->getComm());

    coarseCoordinates = Xpetra::MultiVectorFactory<real_type, LO, GO, NO>::Build(rowMap,
                                                                                 numDimensions);
    Array<ArrayRCP<const real_type> > fineCoordData(numDimensions);
    Array<ArrayRCP<real_type> > coarseCoordData(numDimensions);
    for(int dim = 0; dim < numDimensions; ++dim) {
      fineCoordData[dim] = fineCoordinates->getData(dim);
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
    const LO numEdges     = 4*(lCoarseNodesPerDim[0] - 2)
      + 4*(lCoarseNodesPerDim[1] - 2)
      + 4*(lCoarseNodesPerDim[2] - 2);
    const LO numFaces     = 2*(lCoarseNodesPerDim[0] - 2)*(lCoarseNodesPerDim[1] - 2)
      + 2*(lCoarseNodesPerDim[0] - 2)*(lCoarseNodesPerDim[2] - 2)
      + 2*(lCoarseNodesPerDim[1] - 2)*(lCoarseNodesPerDim[2] - 2);
    const LO numInteriors = (lCoarseNodesPerDim[0] - 2)*(lCoarseNodesPerDim[1] - 2)
      *(lCoarseNodesPerDim[2] - 2);

    const LO nnz = numCorners*cornerStencilLength + numEdges*edgeStencilLength
      + numFaces*faceStencilLength + numInteriors*interiorStencilLength;

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
    Array<SC> coeffs({1.0/3.0, 2.0/3.0, 1.0, 2.0/3.0, 1.0/3.0});
    row_map_h(0) = 0;

    // Define some offsets that
    // will be needed often later on
    const LO edgeLineOffset = 2*cornerStencilLength + (lCoarseNodesPerDim[0] - 2)*edgeStencilLength;
    const LO faceLineOffset = 2*edgeStencilLength + (lCoarseNodesPerDim[0] - 2)*faceStencilLength;
    const LO interiorLineOffset = 2*faceStencilLength
      + (lCoarseNodesPerDim[0] - 2)*interiorStencilLength;

    const LO facePlaneOffset = 2*edgeLineOffset + (lCoarseNodesPerDim[1] - 2)*faceLineOffset;
    const LO interiorPlaneOffset = 2*faceLineOffset + (lCoarseNodesPerDim[1] - 2)*interiorLineOffset;

    // Let us take care of the corners
    // first since we always have
    // corners to deal with!
    {
      // Corner 1
      LO rowIdx = 0, columnOffset = 0, entryOffset = 0;
      row_map_h(rowIdx + 1) = entryOffset + cornerStencilLength;
      for(LO k = 0; k < 3; ++k) {
        for(LO j = 0; j < 3; ++j) {
          for(LO i = 0; i < 3; ++i) {
            entries_h(entryOffset + k*9 + j*3 + i) = columnOffset
              + k*lFineNodesPerDim[1]*lFineNodesPerDim[0] + j*lFineNodesPerDim[0] + i;
            values_h(entryOffset + k*9 + j*3 + i)  = coeffs[k + 2]*coeffs[j + 2]*coeffs[i + 2];
          }
        }
      }
      for(int dim = 0; dim <numDimensions; ++dim) {
        coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
      }

      // Corner 5
      rowIdx += (lCoarseNodesPerDim[2] - 1)*lCoarseNodesPerDim[1]*lCoarseNodesPerDim[0];
      columnOffset += (lFineNodesPerDim[2] - 1)*lFineNodesPerDim[1]*lFineNodesPerDim[0];
      entryOffset += facePlaneOffset + (lCoarseNodesPerDim[2] - 2)*interiorPlaneOffset;
      row_map_h(rowIdx + 1) = entryOffset + cornerStencilLength;
      for(LO k = 0; k < 3; ++k) {
        for(LO j = 0; j < 3; ++j) {
          for(LO i = 0; i < 3; ++i) {
            entries_h(entryOffset + k*9 + j*3 + i) = columnOffset
              + (k - 2)*lFineNodesPerDim[1]*lFineNodesPerDim[0] + j*lFineNodesPerDim[0] + i;
            values_h(entryOffset + k*9 + j*3 + i)  = coeffs[k]*coeffs[j + 2]*coeffs[i + 2];
          }
        }
      }
      for(int dim = 0; dim <numDimensions; ++dim) {
        coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
      }

      // Corner 2
      rowIdx = lCoarseNodesPerDim[0] - 1;
      columnOffset = lFineNodesPerDim[0] - 1;
      entryOffset = cornerStencilLength + (lCoarseNodesPerDim[0] - 2)*edgeStencilLength;

      row_map_h(rowIdx + 1) = entryOffset + cornerStencilLength;
      for(LO k = 0; k < 3; ++k) {
        for(LO j = 0; j < 3; ++j) {
          for(LO i = 0; i < 3; ++i) {
            entries_h(entryOffset + k*9 + j*3 + i) = columnOffset
              + k*lFineNodesPerDim[1]*lFineNodesPerDim[0]
              + j*lFineNodesPerDim[0] + (i - 2);
            values_h(entryOffset + k*9 + j*3 + i)  = coeffs[k + 2]*coeffs[j + 2]*coeffs[i];
          }
        }
      }
      for(int dim = 0; dim <numDimensions; ++dim) {
        coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
      }

      // Corner 6
      rowIdx += (lCoarseNodesPerDim[2] - 1)*lCoarseNodesPerDim[1]*lCoarseNodesPerDim[0];
      columnOffset += (lFineNodesPerDim[2] - 1)*lFineNodesPerDim[1]*lFineNodesPerDim[0];
      entryOffset += facePlaneOffset + (lCoarseNodesPerDim[2] - 2)*interiorPlaneOffset;
      row_map_h(rowIdx + 1) = entryOffset + cornerStencilLength;
      for(LO k = 0; k < 3; ++k) {
        for(LO j = 0; j < 3; ++j) {
          for(LO i = 0; i < 3; ++i) {
            entries_h(entryOffset + k*9 + j*3 + i) = columnOffset
              + (k - 2)*lFineNodesPerDim[1]*lFineNodesPerDim[0] + j*lFineNodesPerDim[0] + i - 2;
            values_h(entryOffset + k*9 + j*3 + i)  = coeffs[k]*coeffs[j + 2]*coeffs[i];
          }
        }
      }
      for(int dim = 0; dim <numDimensions; ++dim) {
        coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
      }

      // Corner 3
      rowIdx = (lCoarseNodesPerDim[1] - 1)*lCoarseNodesPerDim[0];
      columnOffset = (lFineNodesPerDim[1] - 1)*lFineNodesPerDim[0];
      entryOffset = edgeLineOffset + (lCoarseNodesPerDim[1] - 2)*faceLineOffset;

      row_map_h(rowIdx + 1) = entryOffset + cornerStencilLength;
      for(LO k = 0; k < 3; ++k) {
        for(LO j = 0; j < 3; ++j) {
          for(LO i = 0; i < 3; ++i) {
            entries_h(entryOffset + k*9 + j*3 + i) = columnOffset
              + k*lFineNodesPerDim[1]*lFineNodesPerDim[0]
              + (j - 2)*lFineNodesPerDim[0] + i;
            values_h(entryOffset + k*9 + j*3 + i)  = coeffs[k + 2]*coeffs[j]*coeffs[i + 2];
          }
        }
      }
      for(int dim = 0; dim <numDimensions; ++dim) {
        coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
      }

      // Corner 7
      rowIdx += (lCoarseNodesPerDim[2] - 1)*lCoarseNodesPerDim[1]*lCoarseNodesPerDim[0];
      columnOffset += (lFineNodesPerDim[2] - 1)*lFineNodesPerDim[1]*lFineNodesPerDim[0];
      entryOffset += facePlaneOffset + (lCoarseNodesPerDim[2] - 2)*interiorPlaneOffset;
      row_map_h(rowIdx + 1) = entryOffset + cornerStencilLength;
      for(LO k = 0; k < 3; ++k) {
        for(LO j = 0; j < 3; ++j) {
          for(LO i = 0; i < 3; ++i) {
            entries_h(entryOffset + k*9 + j*3 + i) = columnOffset
              + (k - 2)*lFineNodesPerDim[1]*lFineNodesPerDim[0] + (j - 2)*lFineNodesPerDim[0] + i;
            values_h(entryOffset + k*9 + j*3 + i)  = coeffs[k]*coeffs[j]*coeffs[i + 2];
          }
        }
      }
      for(int dim = 0; dim <numDimensions; ++dim) {
        coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
      }

      // Corner 4
      rowIdx = lCoarseNodesPerDim[1]*lCoarseNodesPerDim[0] - 1;
      columnOffset = lFineNodesPerDim[1]*lFineNodesPerDim[0] - 1;
      entryOffset = edgeLineOffset + (lCoarseNodesPerDim[1] - 2)*faceLineOffset +
        cornerStencilLength + (lCoarseNodesPerDim[0] - 2)*edgeStencilLength;

      row_map_h(rowIdx + 1) = entryOffset + cornerStencilLength;
      for(LO k = 0; k < 3; ++k) {
        for(LO j = 0; j < 3; ++j) {
          for(LO i = 0; i < 3; ++i) {
            entries_h(entryOffset + k*9 + j*3 + i) = columnOffset
              + k*lFineNodesPerDim[1]*lFineNodesPerDim[0]
              + (j - 2)*lFineNodesPerDim[0] + (i - 2);
            values_h(entryOffset + k*9 + j*3 + i)  = coeffs[k + 2]*coeffs[j]*coeffs[i];
          }
        }
      }
      for(int dim = 0; dim <numDimensions; ++dim) {
        coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
      }

      // Corner 8
      rowIdx += (lCoarseNodesPerDim[2] - 1)*lCoarseNodesPerDim[1]*lCoarseNodesPerDim[0];
      columnOffset += (lFineNodesPerDim[2] - 1)*lFineNodesPerDim[1]*lFineNodesPerDim[0];
      entryOffset += facePlaneOffset + (lCoarseNodesPerDim[2] - 2)*interiorPlaneOffset;
      row_map_h(rowIdx + 1) = entryOffset + cornerStencilLength;
      for(LO k = 0; k < 3; ++k) {
        for(LO j = 0; j < 3; ++j) {
          for(LO i = 0; i < 3; ++i) {
            entries_h(entryOffset + k*9 + j*3 + i) = columnOffset
              + (k - 2)*lFineNodesPerDim[1]*lFineNodesPerDim[0] + (j - 2)*lFineNodesPerDim[0] + (i - 2);
            values_h(entryOffset + k*9 + j*3 + i)  = coeffs[k]*coeffs[j]*coeffs[i];
          }
        }
      }
      for(int dim = 0; dim <numDimensions; ++dim) {
        coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
      }
    } // Corners are done!

    // Edges along 0 direction
    if(lCoarseNodesPerDim[0] - 2 > 0) {

      LO rowIdx, columnOffset, entryOffset;
      for(LO edgeIdx = 0; edgeIdx < lCoarseNodesPerDim[0] - 2; ++edgeIdx) {

        // Edge 0
        rowIdx = edgeIdx + 1;
        columnOffset = (edgeIdx + 1)*3;
        entryOffset  = cornerStencilLength + edgeIdx*edgeStencilLength;
        row_map_h(rowIdx + 1) = entryOffset + edgeStencilLength;
        for(LO k = 0; k < 3; ++k) {
          for(LO j = 0; j < 3; ++j) {
            for(LO i = 0; i < 5; ++i) {
              entries_h(entryOffset + k*15 + j*5 + i) = columnOffset
                + k*lFineNodesPerDim[1]*lFineNodesPerDim[0] + j*lFineNodesPerDim[0] + i - 2;
              values_h(entryOffset + k*15 + j*5 + i)  = coeffs[k + 2]*coeffs[j + 2]*coeffs[i];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }

        // Edge 1
        rowIdx = (lCoarseNodesPerDim[1] - 1)*lCoarseNodesPerDim[0] + edgeIdx + 1;
        columnOffset = (lFineNodesPerDim[1] - 1)*lFineNodesPerDim[0] + (edgeIdx + 1)*3;
        entryOffset  = edgeLineOffset + (lCoarseNodesPerDim[1] - 2)*faceLineOffset
          + cornerStencilLength + edgeIdx*edgeStencilLength;
        row_map_h(rowIdx + 1) = entryOffset + edgeStencilLength;
        for(LO k = 0; k < 3; ++k) {
          for(LO j = 0; j < 3; ++j) {
            for(LO i = 0; i < 5; ++i) {
              entries_h(entryOffset + k*15 + j*5 + i) = columnOffset
                + k*lFineNodesPerDim[1]*lFineNodesPerDim[0] + (j - 2)*lFineNodesPerDim[0] + i - 2;
              values_h(entryOffset + k*15 + j*5 + i)  = coeffs[k + 2]*coeffs[j]*coeffs[i];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }

        // Edge 2
        rowIdx = (lCoarseNodesPerDim[2] - 1)*lCoarseNodesPerDim[1]*lCoarseNodesPerDim[0]
          + edgeIdx + 1;
        columnOffset = (lFineNodesPerDim[2] - 1)*lFineNodesPerDim[1]*lFineNodesPerDim[0]
          + (edgeIdx + 1)*3;
        entryOffset  = facePlaneOffset + (lCoarseNodesPerDim[2] - 2)*interiorPlaneOffset
          + cornerStencilLength + edgeIdx*edgeStencilLength;
        row_map_h(rowIdx + 1) = entryOffset + edgeStencilLength;
        for(LO k = 0; k < 3; ++k) {
          for(LO j = 0; j < 3; ++j) {
            for(LO i = 0; i < 5; ++i) {
              entries_h(entryOffset + k*15 + j*5 + i) = columnOffset
                + (k - 2)*lFineNodesPerDim[1]*lFineNodesPerDim[0] + j*lFineNodesPerDim[0] + i - 2;
              values_h(entryOffset + k*15 + j*5 + i)  = coeffs[k]*coeffs[j + 2]*coeffs[i];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }

        // Edge 3
        rowIdx = (lCoarseNodesPerDim[2] - 1)*lCoarseNodesPerDim[1]*lCoarseNodesPerDim[0]
          + (lCoarseNodesPerDim[1] - 1)*lCoarseNodesPerDim[0] + edgeIdx + 1;
        columnOffset = (lFineNodesPerDim[2] - 1)*lFineNodesPerDim[1]*lFineNodesPerDim[0]
          + (lFineNodesPerDim[1] - 1)*lFineNodesPerDim[0] + (edgeIdx + 1)*3;
        entryOffset  =  facePlaneOffset + (lCoarseNodesPerDim[2] - 2)*interiorPlaneOffset
          + edgeLineOffset + (lCoarseNodesPerDim[1] - 2)*faceLineOffset
          + cornerStencilLength + edgeIdx*edgeStencilLength;
        row_map_h(rowIdx + 1) = entryOffset + edgeStencilLength;
        for(LO k = 0; k < 3; ++k) {
          for(LO j = 0; j < 3; ++j) {
            for(LO i = 0; i < 5; ++i) {
              entries_h(entryOffset + k*15 + j*5 + i) = columnOffset
                + (k - 2)*lFineNodesPerDim[1]*lFineNodesPerDim[0]
                + (j - 2)*lFineNodesPerDim[0] + i - 2;
              values_h(entryOffset + k*15 + j*5 + i)  = coeffs[k]*coeffs[j]*coeffs[i];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }
      }
    }

    // Edges along 1 direction
    if(lCoarseNodesPerDim[1] - 2 > 0) {

      LO rowIdx, columnOffset, entryOffset;
      for(LO edgeIdx = 0; edgeIdx < lCoarseNodesPerDim[1] - 2; ++edgeIdx) {

        // Edge 0
        rowIdx = (edgeIdx + 1)*lCoarseNodesPerDim[0];
        columnOffset = (edgeIdx + 1)*3*lFineNodesPerDim[0];
        entryOffset  = edgeLineOffset + edgeIdx*faceLineOffset;
        row_map_h(rowIdx + 1) = entryOffset + edgeStencilLength;
        for(LO k = 0; k < 3; ++k) {
          for(LO j = 0; j < 5; ++j) {
            for(LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k*15 + j*3 + i) = columnOffset
                + k*lFineNodesPerDim[1]*lFineNodesPerDim[0] + (j - 2)*lFineNodesPerDim[0] + i;
              values_h(entryOffset + k*15 + j*3 + i)  = coeffs[k + 2]*coeffs[j]*coeffs[i + 2];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }

        // Edge 1
        rowIdx = (edgeIdx + 1)*lCoarseNodesPerDim[0] + lCoarseNodesPerDim[0] - 1;
        columnOffset = (edgeIdx + 1)*3*lFineNodesPerDim[0] + lFineNodesPerDim[0] - 1;
        entryOffset  = edgeLineOffset + edgeIdx*faceLineOffset
          + edgeStencilLength + (lCoarseNodesPerDim[0] - 2)*faceStencilLength;
        row_map_h(rowIdx + 1) = entryOffset + edgeStencilLength;
        for(LO k = 0; k < 3; ++k) {
          for(LO j = 0; j < 5; ++j) {
            for(LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k*15 + j*3 + i) = columnOffset
                + k*lFineNodesPerDim[1]*lFineNodesPerDim[0] + (j - 2)*lFineNodesPerDim[0] + i - 2;
              values_h(entryOffset + k*15 + j*3 + i)  = coeffs[k + 2]*coeffs[j]*coeffs[i];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }

        // Edge 2
        rowIdx = (lCoarseNodesPerDim[2] - 1)*lCoarseNodesPerDim[1]*lCoarseNodesPerDim[0]
          + (edgeIdx + 1)*lCoarseNodesPerDim[0];
        columnOffset = (lFineNodesPerDim[2] - 1)*lFineNodesPerDim[1]*lFineNodesPerDim[0]
          + (edgeIdx + 1)*3*lFineNodesPerDim[0];
        entryOffset  = facePlaneOffset + (lCoarseNodesPerDim[2] - 2)*interiorPlaneOffset
          + edgeLineOffset + edgeIdx*faceLineOffset;
        row_map_h(rowIdx + 1) = entryOffset + edgeStencilLength;
        for(LO k = 0; k < 3; ++k) {
          for(LO j = 0; j < 5; ++j) {
            for(LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k*15 + j*3 + i) = columnOffset
                + (k - 2)*lFineNodesPerDim[1]*lFineNodesPerDim[0] + (j - 2)*lFineNodesPerDim[0] + i;
              values_h(entryOffset + k*15 + j*3 + i)  = coeffs[k]*coeffs[j]*coeffs[i + 2];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }

        // Edge 3
        rowIdx = (lCoarseNodesPerDim[2] - 1)*lCoarseNodesPerDim[1]*lCoarseNodesPerDim[0]
          + (edgeIdx + 1)*lCoarseNodesPerDim[0] + lCoarseNodesPerDim[0] - 1;
        columnOffset = (lFineNodesPerDim[2] - 1)*lFineNodesPerDim[1]*lFineNodesPerDim[0]
          + (edgeIdx + 1)*3*lFineNodesPerDim[0] + lFineNodesPerDim[0] - 1;
        entryOffset  = facePlaneOffset + (lCoarseNodesPerDim[2] - 2)*interiorPlaneOffset
          + edgeLineOffset + edgeIdx*faceLineOffset
          + edgeStencilLength + (lCoarseNodesPerDim[0] - 2)*faceStencilLength;
        row_map_h(rowIdx + 1) = entryOffset + edgeStencilLength;
        for(LO k = 0; k < 3; ++k) {
          for(LO j = 0; j < 5; ++j) {
            for(LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k*15 + j*3 + i) = columnOffset
                + (k - 2)*lFineNodesPerDim[1]*lFineNodesPerDim[0]
                + (j - 2)*lFineNodesPerDim[0] + i - 2;
              values_h(entryOffset + k*15 + j*3 + i)  = coeffs[k]*coeffs[j]*coeffs[i];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }
      }
    }

    // Edges along 2 direction
    if(lCoarseNodesPerDim[2] - 2 > 0) {

      LO rowIdx, columnOffset, entryOffset;
      for(LO edgeIdx = 0; edgeIdx < lCoarseNodesPerDim[2] - 2; ++edgeIdx) {

        // Edge 0
        rowIdx = (edgeIdx + 1)*lCoarseNodesPerDim[1]*lCoarseNodesPerDim[0];
        columnOffset = (edgeIdx + 1)*3*lFineNodesPerDim[1]*lFineNodesPerDim[0];
        entryOffset  = facePlaneOffset + edgeIdx*interiorPlaneOffset;
        row_map_h(rowIdx + 1) = entryOffset + edgeStencilLength;
        for(LO k = 0; k < 5; ++k) {
          for(LO j = 0; j < 3; ++j) {
            for(LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k*9 + j*3 + i) = columnOffset
                + (k - 2)*lFineNodesPerDim[1]*lFineNodesPerDim[0] + j*lFineNodesPerDim[0] + i;
              values_h(entryOffset + k*9 + j*3 + i)  = coeffs[k]*coeffs[j + 2]*coeffs[i + 2];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }

        // Edge 1
        rowIdx = (edgeIdx + 1)*lCoarseNodesPerDim[1]*lCoarseNodesPerDim[0]
          + lCoarseNodesPerDim[0] - 1;
        columnOffset = (edgeIdx + 1)*3*lFineNodesPerDim[1]*lFineNodesPerDim[0]
          + lFineNodesPerDim[0] - 1;
        entryOffset  = facePlaneOffset + faceLineOffset - edgeStencilLength
          + edgeIdx*interiorPlaneOffset;
        row_map_h(rowIdx + 1) = entryOffset + edgeStencilLength;
        for(LO k = 0; k < 5; ++k) {
          for(LO j = 0; j < 3; ++j) {
            for(LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k*9 + j*3 + i) = columnOffset
                + (k - 2)*lFineNodesPerDim[1]*lFineNodesPerDim[0] + j*lFineNodesPerDim[0] + i - 2;
              values_h(entryOffset + k*9 + j*3 + i)  = coeffs[k]*coeffs[j + 2]*coeffs[i];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }

        // Edge 2
        rowIdx = (edgeIdx + 1)*lCoarseNodesPerDim[1]*lCoarseNodesPerDim[0]
          + (lCoarseNodesPerDim[1] - 1)*lCoarseNodesPerDim[0];
        columnOffset = (edgeIdx + 1)*3*lFineNodesPerDim[1]*lFineNodesPerDim[0]
          + (lFineNodesPerDim[1] - 1)*lFineNodesPerDim[0];
        entryOffset  = facePlaneOffset + edgeIdx*interiorPlaneOffset + faceLineOffset
          + (lCoarseNodesPerDim[1] - 2)*interiorLineOffset;
        row_map_h(rowIdx + 1) = entryOffset + edgeStencilLength;
        for(LO k = 0; k < 5; ++k) {
          for(LO j = 0; j < 3; ++j) {
            for(LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k*9 + j*3 + i) = columnOffset
                + (k - 2)*lFineNodesPerDim[1]*lFineNodesPerDim[0] + (j - 2)*lFineNodesPerDim[0] + i;
              values_h(entryOffset + k*9 + j*3 + i)  = coeffs[k]*coeffs[j]*coeffs[i + 2];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }

        // Edge 3
        rowIdx = (edgeIdx + 2)*lCoarseNodesPerDim[1]*lCoarseNodesPerDim[0] - 1;
        columnOffset = (edgeIdx + 1)*3*lFineNodesPerDim[1]*lFineNodesPerDim[0]
          + lFineNodesPerDim[1]*lFineNodesPerDim[0] - 1;
        entryOffset  = facePlaneOffset + (edgeIdx + 1)*interiorPlaneOffset - edgeStencilLength;
        row_map_h(rowIdx + 1) = entryOffset + edgeStencilLength;
        for(LO k = 0; k < 5; ++k) {
          for(LO j = 0; j < 3; ++j) {
            for(LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k*9 + j*3 + i) = columnOffset
                + (k - 2)*lFineNodesPerDim[1]*lFineNodesPerDim[0]
                + (j - 2)*lFineNodesPerDim[0] + i - 2;
              values_h(entryOffset + k*9 + j*3 + i)  = coeffs[k]*coeffs[j]*coeffs[i];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }
      }
    }

    // Faces in 0-1 plane
    if((lCoarseNodesPerDim[0] - 2 > 0) && (lCoarseNodesPerDim[1] - 2 > 0)) {

      Array<LO> gridIdx(3);
      LO rowIdx, columnOffset, entryOffset;
      for(LO faceIdx=0; faceIdx < (lCoarseNodesPerDim[1]-2)*(lCoarseNodesPerDim[0]-2); ++faceIdx) {

        // Face 0
        rowIdx = (gridIdx[1] + 1)*lCoarseNodesPerDim[0] + gridIdx[0] + 1;
        columnOffset = 3*((gridIdx[1] + 1)*lFineNodesPerDim[0] + gridIdx[0] + 1);
        entryOffset  = edgeLineOffset + edgeStencilLength
          + gridIdx[1]*faceLineOffset + gridIdx[0]*faceStencilLength;
        row_map_h(rowIdx + 1) = entryOffset + faceStencilLength;
        for(LO k = 0; k < 3; ++k) {
          for(LO j = 0; j < 5; ++j) {
            for(LO i = 0; i < 5; ++i) {
              entries_h(entryOffset + k*25 + j*5 + i) = columnOffset
                + k*lFineNodesPerDim[1]*lFineNodesPerDim[0] + (j - 2)*lFineNodesPerDim[0] + i - 2;
              values_h(entryOffset + k*25 + j*5 + i)  = coeffs[k + 2]*coeffs[j]*coeffs[i];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }

        // Face 1
        rowIdx += (lCoarseNodesPerDim[2] - 1)*lCoarseNodesPerDim[1]*lCoarseNodesPerDim[0];
        columnOffset += (lFineNodesPerDim[2] - 1)*lFineNodesPerDim[1]*lFineNodesPerDim[0];
        entryOffset  += facePlaneOffset + (lCoarseNodesPerDim[2] - 2)*interiorPlaneOffset;
        row_map_h(rowIdx + 1) = entryOffset + faceStencilLength;
        for(LO k = 0; k < 3; ++k) {
          for(LO j = 0; j < 5; ++j) {
            for(LO i = 0; i < 5; ++i) {
              entries_h(entryOffset + k*25 + j*5 + i) = columnOffset
                + (k - 2)*lFineNodesPerDim[1]*lFineNodesPerDim[0]
                + (j - 2)*lFineNodesPerDim[0] + i - 2;
              values_h(entryOffset + k*25 + j*5 + i)  = coeffs[k]*coeffs[j]*coeffs[i];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }

        // Last step in the loop
        // update the grid indices
        // for next grid point
        ++gridIdx[0];
        if(gridIdx[0] == lCoarseNodesPerDim[0] - 2) {
          gridIdx[0] = 0;
          ++gridIdx[1];
        }
      }
    }

    // Faces in 0-2 plane
    if((lCoarseNodesPerDim[0] - 2 > 0) && (lCoarseNodesPerDim[2] - 2 > 0)) {

      Array<LO> gridIdx(3);
      LO rowIdx, columnOffset, entryOffset;
      for(LO faceIdx=0; faceIdx < (lCoarseNodesPerDim[2]-2)*(lCoarseNodesPerDim[0]-2); ++faceIdx) {

        // Face 0
        rowIdx = (gridIdx[2] + 1)*lCoarseNodesPerDim[1]*lCoarseNodesPerDim[0] + (gridIdx[0] + 1);
        columnOffset = ((gridIdx[2] + 1)*lFineNodesPerDim[1]*lFineNodesPerDim[0]
          + gridIdx[0] + 1)*3;
        entryOffset  = facePlaneOffset + gridIdx[2]*interiorPlaneOffset + edgeStencilLength
          + gridIdx[0]*faceStencilLength;
        row_map_h(rowIdx + 1) = entryOffset + faceStencilLength;
        for(LO k = 0; k < 5; ++k) {
          for(LO j = 0; j < 3; ++j) {
            for(LO i = 0; i < 5; ++i) {
              entries_h(entryOffset + k*15 + j*5 + i) = columnOffset
                + (k - 2)*lFineNodesPerDim[1]*lFineNodesPerDim[0] + j*lFineNodesPerDim[0] + i - 2;
              values_h(entryOffset + k*15 + j*5 + i)  = coeffs[k]*coeffs[j + 2]*coeffs[i];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }

        // Face 1
        rowIdx += (lCoarseNodesPerDim[1] - 1)*lCoarseNodesPerDim[0];
        columnOffset += (lFineNodesPerDim[1] - 1)*lFineNodesPerDim[0];
        entryOffset  += faceLineOffset + (lCoarseNodesPerDim[1] - 2)*interiorLineOffset;
        row_map_h(rowIdx + 1) = entryOffset + faceStencilLength;
        for(LO k = 0; k < 5; ++k) {
          for(LO j = 0; j < 3; ++j) {
            for(LO i = 0; i < 5; ++i) {
              entries_h(entryOffset + k*15 + j*5 + i) = columnOffset
                + (k - 2)*lFineNodesPerDim[1]*lFineNodesPerDim[0]
                + (j - 2)*lFineNodesPerDim[0] + i - 2;
              values_h(entryOffset + k*15 + j*5 + i)  = coeffs[k]*coeffs[j]*coeffs[i];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }

        // Last step in the loop
        // update the grid indices
        // for next grid point
        ++gridIdx[0];
        if(gridIdx[0] == lCoarseNodesPerDim[0] - 2) {
          gridIdx[0] = 0;
          ++gridIdx[2];
        }
      }
    }

    // Faces in 1-2 plane
    if((lCoarseNodesPerDim[1] - 2 > 0) && (lCoarseNodesPerDim[2] - 2 > 0)) {

      Array<LO> gridIdx(3);
      LO rowIdx, columnOffset, entryOffset;
      for(LO faceIdx=0; faceIdx < (lCoarseNodesPerDim[2]-2)*(lCoarseNodesPerDim[1]-2); ++faceIdx) {

        // Face 0
        rowIdx = (gridIdx[2] + 1)*lCoarseNodesPerDim[1]*lCoarseNodesPerDim[0]
          + (gridIdx[1] + 1)*lCoarseNodesPerDim[0];
        columnOffset = ((gridIdx[2] + 1)*lFineNodesPerDim[1]*lFineNodesPerDim[0]
                        + (gridIdx[1] + 1)*lFineNodesPerDim[0])*3;
        entryOffset  = facePlaneOffset + gridIdx[2]*interiorPlaneOffset + faceLineOffset
          + gridIdx[1]*interiorLineOffset;
        row_map_h(rowIdx + 1) = entryOffset + faceStencilLength;
        for(LO k = 0; k < 5; ++k) {
          for(LO j = 0; j < 5; ++j) {
            for(LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k*15 + j*3 + i) = columnOffset
                + (k - 2)*lFineNodesPerDim[1]*lFineNodesPerDim[0] + (j - 2)*lFineNodesPerDim[0] + i;
              values_h(entryOffset + k*15 + j*3 + i)  = coeffs[k]*coeffs[j]*coeffs[i + 2];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }

        // Face 1
        rowIdx += lCoarseNodesPerDim[0] - 1;
        columnOffset += lFineNodesPerDim[0] - 1;
        entryOffset  += faceStencilLength + (lCoarseNodesPerDim[0] - 2)*interiorStencilLength;
        row_map_h(rowIdx + 1) = entryOffset + faceStencilLength;
        for(LO k = 0; k < 5; ++k) {
          for(LO j = 0; j < 5; ++j) {
            for(LO i = 0; i < 3; ++i) {
              entries_h(entryOffset + k*15 + j*3 + i) = columnOffset
                + (k - 2)*lFineNodesPerDim[1]*lFineNodesPerDim[0]
                + (j - 2)*lFineNodesPerDim[0] + i - 2;
              values_h(entryOffset + k*15 + j*3 + i)  = coeffs[k]*coeffs[j]*coeffs[i];
            }
          }
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }

        // Last step in the loop
        // update the grid indices
        // for next grid point
        ++gridIdx[1];
        if(gridIdx[1] == lCoarseNodesPerDim[1] - 2) {
          gridIdx[1] = 0;
          ++gridIdx[2];
        }
      }
    }

    if(numInteriors > 0) {
      // Allocate and compute arrays
      // containing column offsets
      // and values associated with
      // interior points
      LO countRowEntries = 0;
      Array<LO> columnOffsets(125);
      for(LO k = -2; k < 3; ++k) {
        for(LO j = -2; j < 3; ++j) {
          for(LO i = -2; i < 3; ++i) {
            columnOffsets[countRowEntries] = k*lFineNodesPerDim[1]*lFineNodesPerDim[0]
              + j*lFineNodesPerDim[0] + i;
            ++countRowEntries;
          }
        }
      }

      LO countValues = 0;
      Array<SC> interiorValues(125);
      for(LO k = 0; k < 5; ++k) {
        for(LO j = 0; j < 5; ++j) {
          for(LO i = 0; i < 5; ++i) {
            interiorValues[countValues] = coeffs[k]*coeffs[j]*coeffs[i];
            ++countValues;
          }
        }
      }

      LO rowIdx, columnOffset, entryOffset;
      Array<LO> gridIdx(3);
      for(LO interiorIdx = 0; interiorIdx < numInteriors; ++interiorIdx) {
        rowIdx = (gridIdx[2] + 1)*lCoarseNodesPerDim[0]*lCoarseNodesPerDim[1]
          + (gridIdx[1] + 1)*lCoarseNodesPerDim[0]
          + gridIdx[0] + 1;
        columnOffset = (gridIdx[2] + 1)*3*lFineNodesPerDim[1]*lFineNodesPerDim[0]
          + (gridIdx[1] + 1)*3*lFineNodesPerDim[0] + (gridIdx[0] + 1)*3;

        entryOffset = facePlaneOffset + faceLineOffset + faceStencilLength
          + gridIdx[2]*interiorPlaneOffset + gridIdx[1]*interiorLineOffset
          + gridIdx[0]*interiorStencilLength;
        row_map_h(rowIdx + 1) = entryOffset + interiorStencilLength;

        // Fill the column indices
        // and values in the approproate
        // views.
        for(LO entryIdx = 0; entryIdx < interiorStencilLength; ++entryIdx) {
          entries_h(entryOffset + entryIdx) = columnOffset + columnOffsets[entryIdx];
          values_h(entryOffset + entryIdx) = interiorValues[entryIdx];
        }
        for(int dim = 0; dim <numDimensions; ++dim) {
          coarseCoordData[dim][rowIdx] = fineCoordData[dim][columnOffset];
        }

        // Last step in the loop
        // update the grid indices
        // for next grid point
        ++gridIdx[0];
        if(gridIdx[0] == lCoarseNodesPerDim[0] - 2) {
          gridIdx[0] = 0;
          ++gridIdx[1];
          if(gridIdx[1] == lCoarseNodesPerDim[1] - 2) {
            gridIdx[1] = 0;
            ++gridIdx[2];
          }
        }
      }
    }

    Kokkos::deep_copy(row_map, row_map_h);
    Kokkos::deep_copy(entries, entries_h);
    Kokkos::deep_copy(values,  values_h);

    local_graph_type localGraph(entries, row_map);
    local_matrix_type localR("R", numCols, values, localGraph);

    R = MatrixFactory::Build(localR,            // the local data
                             rowMap,            // rowMap
                             A->getRowMap(),    // colMap
                             A->getRowMap(),    // domainMap == colMap
                             rowMap,            // rangeMap  == rowMap
                             Teuchos::null);    // params for optimized construction

  } // Build3D()

} //namespace MueLu

#define MUELU_REGIONRFACTORY_SHORT
#endif // MUELU_REGIONRFACTORY_DEF_HPP
