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
#ifndef MUELU_COALESCEDROPFACTORY_KOKKOS_DEF_HPP
#define MUELU_COALESCEDROPFACTORY_KOKKOS_DEF_HPP

#ifdef HAVE_MUELU_KOKKOS_REFACTOR
#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

#include "Xpetra_Matrix.hpp"

#include "MueLu_CoalesceDropFactory_kokkos_decl.hpp"

#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_LWGraph_kokkos.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities_kokkos.hpp"

namespace MueLu {


  namespace { // anonymous

    template<class LO, class RowType>
    class ScanFunctor {
    public:
      ScanFunctor(RowType rows_) : rows(rows_) { }

      KOKKOS_INLINE_FUNCTION
      void operator()(const LO i, LO& upd, const bool& final) const {
        upd += rows(i);
        if (final)
          rows(i) = upd;
      }

    private:
      RowType rows;
    };

    template<class LO, class GhostedViewType>
    class ClassicalDropFunctor {
    private:
      typedef typename GhostedViewType::value_type      SC;
      typedef          Kokkos::ArithTraits<SC>          ATS;
      typedef typename ATS::magnitudeType               magnitudeType;

      GhostedViewType   diag;   // corresponds to overlapped diagonal multivector (2D View)
      magnitudeType     eps;

    public:
      ClassicalDropFunctor(GhostedViewType ghostedDiag, magnitudeType threshold) :
          diag(ghostedDiag),
          eps(threshold)
      { }

      // Return true if we drop, false if not
      KOKKOS_FORCEINLINE_FUNCTION
      bool operator()(LO row, LO col, SC val) const {
        // We avoid square root by using squared values
        auto aiiajj = ATS::magnitude(diag(row, 0)) * ATS::magnitude(diag(col, 0));   // |a_ii|*|a_jj|
        auto aij2   = ATS::magnitude(val)          * ATS::magnitude(val);            // |a_ij|^2

        return (aij2 <= eps*eps * aiiajj);
      }
    };

    template<class LO, class CoordsType>
    class DistanceFunctor {
    private:
      typedef typename CoordsType::value_type           SC;
      typedef          Kokkos::ArithTraits<SC>          ATS;
      typedef typename ATS::magnitudeType               magnitudeType;

    public:
      typedef SC value_type;

    public:
      DistanceFunctor(CoordsType coords_) : coords(coords_) { }

      KOKKOS_INLINE_FUNCTION
      magnitudeType distance2(LO row, LO col) const {
        SC d = ATS::zero(), s;
        for (size_t j = 0; j < coords.extent(1); j++) {
          s = coords(row,j) - coords(col,j);
          d += s*s;
        }
        return ATS::magnitude(d);
      }
    private:
      CoordsType coords;
    };

    template<class LO, class GhostedViewType, class DistanceFunctor>
    class DistanceLaplacianDropFunctor {
    private:
      typedef typename GhostedViewType::value_type      SC;
      typedef          Kokkos::ArithTraits<SC>          ATS;
      typedef typename ATS::magnitudeType               magnitudeType;

    public:
      DistanceLaplacianDropFunctor(GhostedViewType ghostedLaplDiag, DistanceFunctor distFunctor_, magnitudeType threshold) :
          diag(ghostedLaplDiag),
          distFunctor(distFunctor_),
          eps(threshold)
      { }

      // Return true if we drop, false if not
      KOKKOS_INLINE_FUNCTION
      bool operator()(LO row, LO col, SC val) const {
        // We avoid square root by using squared values

        // We ignore incoming value of val as we operate on an auxiliary
        // distance Laplacian matrix
        typedef typename DistanceFunctor::value_type      dSC;
        typedef          Kokkos::ArithTraits<dSC>         dATS;
        auto fval = dATS::one() / distFunctor.distance2(row, col);

        auto aiiajj = ATS::magnitude(diag(row, 0)) * ATS::magnitude(diag(col, 0));   // |a_ii|*|a_jj|
        auto aij2   = ATS::magnitude(fval)         * ATS::magnitude(fval);           // |a_ij|^2

        return (aij2 <= eps*eps * aiiajj);
      }

    private:
      GhostedViewType   diag;   // corresponds to overlapped diagonal multivector (2D View)
      DistanceFunctor   distFunctor;
      magnitudeType     eps;
    };

    template<class SC, class LO, class MatrixType, class BndViewType, class DropFunctorType>
    class ScalarFunctor {
    private:
      typedef typename MatrixType::StaticCrsGraphType   graph_type;
      typedef typename graph_type::row_map_type         rows_type;
      typedef typename graph_type::entries_type         cols_type;
      typedef typename MatrixType::values_type          vals_type;
      typedef          Kokkos::ArithTraits<SC>          ATS;
      typedef typename ATS::magnitudeType               magnitudeType;

    public:
      ScalarFunctor(MatrixType A_, BndViewType bndNodes_, DropFunctorType dropFunctor_,
                    typename rows_type::non_const_type rows_,
                    typename cols_type::non_const_type colsAux_,
                    typename vals_type::non_const_type valsAux_,
                    bool reuseGraph_, bool lumping_, SC threshold_) :
          A(A_),
          bndNodes(bndNodes_),
          dropFunctor(dropFunctor_),
          rows(rows_),
          colsAux(colsAux_),
          valsAux(valsAux_),
          reuseGraph(reuseGraph_),
          lumping(lumping_)
      {
        rowsA = A.graph.row_map;
        zero = ATS::zero();
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const LO row, LO& nnz) const {
        auto rowView = A.rowConst(row);
        auto length  = rowView.length;
        auto offset  = rowsA(row);

        SC diag = zero;
        LO rownnz = 0;
        LO diagID = -1;
        for (decltype(length) colID = 0; colID < length; colID++) {
          LO col = rowView.colidx(colID);
          SC val = rowView.value (colID);

          if (!dropFunctor(row, col, rowView.value(colID)) || row == col) {
            colsAux(offset+rownnz) = col;

            LO valID = (reuseGraph ? colID : rownnz);
            valsAux(offset+valID) = val;
            if (row == col)
              diagID = valID;

            rownnz++;

          } else {
            // Rewrite with zeros (needed for reuseGraph)
            valsAux(offset+colID) = zero;
            diag += val;
          }
        }
        // How to assert on the device?
        // assert(diagIndex != -1);
        rows(row+1) = rownnz;
        // if (lumping && diagID != -1) {
        if (lumping) {
          // Add diag to the diagonal

          // NOTE_KOKKOS: valsAux was allocated with
          // ViewAllocateWithoutInitializing. This is not a problem here
          // because we explicitly set this value above.
          valsAux(offset+diagID) += diag;
        }

        // If the only element remaining after filtering is diagonal, mark node as boundary
        // FIXME: this should really be replaced by the following
        //    if (indices.size() == 1 && indices[0] == row)
        //        boundaryNodes[row] = true;
        // We do not do it this way now because there is no framework for distinguishing isolated
        // and boundary nodes in the aggregation algorithms
        bndNodes(row) = (rownnz == 1);

        nnz += rownnz;
      }

    private:
      MatrixType                            A;
      BndViewType                           bndNodes;
      DropFunctorType                       dropFunctor;

      rows_type                             rowsA;

      typename rows_type::non_const_type    rows;
      typename cols_type::non_const_type    colsAux;
      typename vals_type::non_const_type    valsAux;

      bool                                  reuseGraph;
      bool                                  lumping;
      SC                                    zero;
    };

    // collect number nonzeros of blkSize rows in nnz_(row+1)
    template<class MatrixType, class NnzType, class blkSizeType>
    class Stage1aVectorFunctor {
    private:
      typedef typename MatrixType::ordinal_type LO;

    public:
      Stage1aVectorFunctor(MatrixType kokkosMatrix_, NnzType nnz_, blkSizeType blkSize_) :
        kokkosMatrix(kokkosMatrix_),
        nnz(nnz_),
        blkSize(blkSize_) { }

      KOKKOS_INLINE_FUNCTION
      void operator()(const LO row, LO& totalnnz) const {

        // the following code is more or less what MergeRows is doing
        // count nonzero entries in all dof rows associated with node row
        LO nodeRowMaxNonZeros = 0;
        for (LO j = 0; j < blkSize; j++) {
          auto rowView = kokkosMatrix.row(row * blkSize + j);
          nodeRowMaxNonZeros += rowView.length;
        }
        nnz(row + 1) = nodeRowMaxNonZeros;
        totalnnz += nodeRowMaxNonZeros;
      }


    private:
      MatrixType kokkosMatrix; //< local matrix part
      NnzType nnz;             //< View containing number of nonzeros for current row
      blkSizeType blkSize;     //< block size (or partial block size in strided maps)
    };


    // build the dof-based column map containing the local dof ids belonging to blkSize rows in matrix
    // the DofIds may not be sorted.
    template<class MatrixType, class NnzType, class blkSizeType, class ColDofType>
    class Stage1bVectorFunctor {
    private:
      typedef typename MatrixType::ordinal_type LO;
      //typedef typename MatrixType::value_type   SC;
      //typedef typename MatrixType::device_type  NO;

    private:
      MatrixType kokkosMatrix; //< local matrix part
      NnzType nnz;             //< View containing dof offsets for dof columns
      blkSizeType blkSize;     //< block size (or partial block size in strided maps)
      ColDofType coldofs;      //< view containing the local dof ids associated with columns for the blkSize rows (not sorted)

    public:
      Stage1bVectorFunctor(MatrixType kokkosMatrix_,
                           NnzType nnz_,
                           blkSizeType blkSize_,
                           ColDofType coldofs_) :
        kokkosMatrix(kokkosMatrix_),
        nnz(nnz_),
        blkSize(blkSize_),
        coldofs(coldofs_) {
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const LO rowNode) const {

        LO pos = nnz(rowNode);
        for (LO j = 0; j < blkSize; j++) {
          auto rowView = kokkosMatrix.row(rowNode * blkSize + j);
          auto numIndices = rowView.length;

          for (decltype(numIndices) k = 0; k < numIndices; k++) {
            auto dofID = rowView.colidx(k);
            coldofs(pos) = dofID;
            pos ++;
          }
        }
      }

    };

    // sort column ids
    // translate them into (unique) node ids
    // count the node column ids per node row
    template<class MatrixType, class ColDofNnzType, class ColDofType, class Dof2NodeTranslationType, class BdryNodeType>
    class Stage1cVectorFunctor {
    private:
      typedef typename MatrixType::ordinal_type LO;

    private:
      ColDofNnzType coldofnnz; //< view containing start and stop indices for subviews
      ColDofType coldofs;      //< view containing the local dof ids associated with columns for the blkSize rows (not sorted)
      Dof2NodeTranslationType dof2node; //< view containing the local node id associated with the local dof id
      ColDofNnzType colnodennz; //< view containing number of column nodes for each node row
      BdryNodeType  bdrynode;  //< view containing with numNodes booleans. True if node is (full) dirichlet boundardy node.

    public:
      Stage1cVectorFunctor(ColDofNnzType coldofnnz_, ColDofType coldofs_, Dof2NodeTranslationType dof2node_, ColDofNnzType colnodennz_, BdryNodeType bdrynode_) :
        coldofnnz(coldofnnz_),
        coldofs(coldofs_),
        dof2node(dof2node_),
        colnodennz(colnodennz_),
        bdrynode(bdrynode_) {
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const LO rowNode, LO& nnz) const {
        LO begin = coldofnnz(rowNode);
        LO end   = coldofnnz(rowNode+1);
        LO n     = end - begin;
        for (LO i = 0; i < (n-1); i++) {
          for (LO j = 0; j < (n-i-1); j++) {
            if (coldofs(j+begin) > coldofs(j+begin+1)) {
              LO temp = coldofs(j+begin);
              coldofs(j+begin) = coldofs(j+begin+1);
              coldofs(j+begin+1) = temp;
            }
          }
        }
        size_t cnt = 0;
        LO lastNodeID = -1;
        for (LO i = 0; i < n; i++) {
          LO dofID  = coldofs(begin + i);
          LO nodeID = dof2node(dofID);
          if(nodeID != lastNodeID) {
            lastNodeID = nodeID;
            coldofs(begin+cnt) = nodeID;
            cnt++;
          }
        }
        if(cnt == 1)
          bdrynode(rowNode) = true;
        else
          bdrynode(rowNode) = false;
        colnodennz(rowNode+1) = cnt;
        nnz += cnt;
      }

    };

    // fill column node id view
    template<class MatrixType, class ColDofNnzType, class ColDofType, class ColNodeNnzType, class ColNodeType>
    class Stage1dVectorFunctor {
    private:
      typedef typename MatrixType::ordinal_type LO;
      typedef typename MatrixType::value_type SC;

    private:
      ColDofType     coldofs;    //< view containing mixed node and dof indices (only input)
      ColDofNnzType  coldofnnz;  //< view containing the start and stop indices for subviews (dofs)
      ColNodeType colnodes;      //< view containing the local node ids associated with columns
      ColNodeNnzType colnodennz; //< view containing start and stop indices for subviews

    public:
      Stage1dVectorFunctor(ColDofType coldofs_, ColDofNnzType coldofnnz_, ColNodeType colnodes_, ColNodeNnzType colnodennz_) :
        coldofs(coldofs_),
        coldofnnz(coldofnnz_),
        colnodes(colnodes_),
        colnodennz(colnodennz_) {
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const LO rowNode) const {
        auto dofbegin  = coldofnnz(rowNode);
        auto nodebegin = colnodennz(rowNode);
        auto nodeend   = colnodennz(rowNode+1);
        auto n = nodeend - nodebegin;

        for (decltype(nodebegin) i = 0; i < n; i++) {
          colnodes(nodebegin + i) = coldofs(dofbegin + i);
        }
      }
    };


  } // namespace

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RCP<const ParameterList> CoalesceDropFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("aggregation: drop tol");
    SET_VALID_ENTRY("aggregation: Dirichlet threshold");
    SET_VALID_ENTRY("aggregation: drop scheme");
    SET_VALID_ENTRY("filtered matrix: use lumping");
    SET_VALID_ENTRY("filtered matrix: reuse graph");
    SET_VALID_ENTRY("filtered matrix: reuse eigenvalue");
    {
      typedef Teuchos::StringToIntegralParameterEntryValidator<int> validatorType;
      validParamList->getEntry("aggregation: drop scheme").setValidator(
        rcp(new validatorType(Teuchos::tuple<std::string>("classical", "distance laplacian"), "aggregation: drop scheme")));
    }
#undef  SET_VALID_ENTRY
    validParamList->set< bool >                  ("lightweight wrap",            true, "Experimental option for lightweight graph access");

    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Generating factory for UnAmalgamationInfo");
    validParamList->set< RCP<const FactoryBase> >("Coordinates",        Teuchos::null, "Generating factory for Coordinates");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CoalesceDropFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "A");
    Input(currentLevel, "UnAmalgamationInfo");

    const ParameterList& pL = GetParameterList();
    if (pL.get<bool>("lightweight wrap") == true) {
      if (pL.get<std::string>("aggregation: drop scheme") == "distance laplacian")
        Input(currentLevel, "Coordinates");
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CoalesceDropFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::
  Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    typedef Teuchos::ScalarTraits<SC> STS;
    const SC zero    = STS::zero();

    auto A         = Get< RCP<Matrix> >(currentLevel, "A");
    LO   blkSize   = A->GetFixedBlockSize();

    auto amalInfo = Get< RCP<AmalgamationInfo> >(currentLevel, "UnAmalgamationInfo");

    const ParameterList& pL = GetParameterList();

    bool doLightWeightWrap = pL.get<bool>("lightweight wrap");
    GetOStream(Warnings0) << "lightweight wrap is deprecated" << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(!doLightWeightWrap, Exceptions::RuntimeError,
                               "MueLu KokkosRefactor only supports \"lightweight wrap\"=\"true\"");

    std::string algo = pL.get<std::string>("aggregation: drop scheme");

    double threshold = pL.get<double>("aggregation: drop tol");
    GetOStream(Runtime0) << "algorithm = \"" << algo << "\": threshold = " << threshold
                         << ", blocksize = " << A->GetFixedBlockSize() << std::endl;

    const typename STS::magnitudeType dirichletThreshold =
        STS::magnitude(as<SC>(pL.get<double>("aggregation: Dirichlet threshold")));

    GO numDropped = 0, numTotal = 0;

    RCP<LWGraph_kokkos> graph;
    LO                  dofsPerNode = -1;

    typedef typename LWGraph_kokkos::boundary_nodes_type boundary_nodes_type;
    boundary_nodes_type boundaryNodes;

    RCP<Matrix> filteredA;
    if (blkSize == 1 && threshold == zero) {
      // Scalar problem without dropping

      // Detect and record rows that correspond to Dirichlet boundary conditions
      boundaryNodes = Utilities_kokkos::DetectDirichletRows(*A, dirichletThreshold);

      // Trivial LWGraph construction
      graph = rcp(new LWGraph_kokkos(A->getLocalMatrix().graph, A->getRowMap(), A->getColMap(), "graph of A"));
      graph->SetBoundaryNodeMap(boundaryNodes);

      numTotal = A->getNodeNumEntries();
      dofsPerNode = 1;

      filteredA = A;

    } else if (blkSize == 1 && threshold != zero) {
      // Scalar problem with dropping

      typedef typename Matrix::local_matrix_type          local_matrix_type;
      typedef typename LWGraph_kokkos::local_graph_type   kokkos_graph_type;
      typedef typename kokkos_graph_type::row_map_type::non_const_type    rows_type;
      typedef typename kokkos_graph_type::entries_type::non_const_type    cols_type;
      typedef typename local_matrix_type::values_type::non_const_type     vals_type;

      LO   numRows      = A->getNodeNumRows();
      auto kokkosMatrix = A->getLocalMatrix();
      auto nnzA  = kokkosMatrix.nnz();
      auto rowsA = kokkosMatrix.graph.row_map;

      typedef Kokkos::ArithTraits<SC>     ATS;
      typedef typename ATS::magnitudeType magnitudeType;

      bool reuseGraph = pL.get<bool>("filtered matrix: reuse graph");
      bool lumping    = pL.get<bool>("filtered matrix: use lumping");
      if (lumping)
        GetOStream(Runtime0) << "Lumping dropped entries" << std::endl;

      // FIXME_KOKKOS: replace by ViewAllocateWithoutInitializing + setting a single value
      rows_type rows   ("FA_rows",     numRows+1);
      cols_type colsAux(Kokkos::ViewAllocateWithoutInitializing("FA_aux_cols"), nnzA);
      vals_type valsAux;
      if (reuseGraph) {
        SubFactoryMonitor m2(*this, "CopyMatrix", currentLevel);

        // Share graph with the original matrix
        filteredA = MatrixFactory::Build(A->getCrsGraph());

        // Do a no-op fill-complete
        RCP<ParameterList> fillCompleteParams(new ParameterList);
        fillCompleteParams->set("No Nonlocal Changes", true);
        filteredA->fillComplete(fillCompleteParams);

        // No need to reuseFill, just modify in place
        valsAux = filteredA->getLocalMatrix().values;

      } else {
        // Need an extra array to compress
        valsAux = vals_type(Kokkos::ViewAllocateWithoutInitializing("FA_aux_vals"), nnzA);
      }

      typename boundary_nodes_type::non_const_type bndNodes(Kokkos::ViewAllocateWithoutInitializing("boundaryNodes"), numRows);

      LO nnzFA = 0;
      {
        if (algo == "classical") {
          // Construct overlapped matrix diagonal
          RCP<Vector> ghostedDiag;
          {
            SubFactoryMonitor m2(*this, "Ghosted diag construction", currentLevel);
            ghostedDiag     = Utilities_kokkos::GetMatrixOverlappedDiagonal(*A);
          }

          // Filter out entries
          {
            SubFactoryMonitor m2(*this, "MainLoop", currentLevel);

            auto ghostedDiagView = ghostedDiag->template getLocalView<DeviceType>();

            ClassicalDropFunctor<LO, decltype(ghostedDiagView)> dropFunctor(ghostedDiagView, threshold);
            ScalarFunctor<SC, LO, local_matrix_type, decltype(bndNodes), decltype(dropFunctor)>
                scalarFunctor(kokkosMatrix, bndNodes, dropFunctor, rows, colsAux, valsAux, reuseGraph, lumping, threshold);

            Kokkos::parallel_reduce("MueLu:CoalesceDropF:Build:scalar_filter:main_loop", range_type(0,numRows),
                                    scalarFunctor, nnzFA);
          }

        } else if (algo == "distance laplacian") {
          typedef Xpetra::MultiVector<double,LO,GO,NO> doubleMultiVector;
          auto coords = Get<RCP<doubleMultiVector> >(currentLevel, "Coordinates");

          auto uniqueMap    = A->getRowMap();
          auto nonUniqueMap = A->getColMap();

          // Construct ghosted coordinates
          RCP<const Import> importer;
          {
            SubFactoryMonitor m2(*this, "Coords Import construction", currentLevel);
            importer = ImportFactory::Build(uniqueMap, nonUniqueMap);
          }
          RCP<doubleMultiVector> ghostedCoords;
          {
            SubFactoryMonitor m2(*this, "Ghosted coords construction", currentLevel);
            ghostedCoords = Xpetra::MultiVectorFactory<double,LO,GO,NO>::Build(nonUniqueMap, coords->getNumVectors());
            ghostedCoords->doImport(*coords, *importer, Xpetra::INSERT);
          }

          auto ghostedCoordsView = ghostedCoords->template getLocalView<DeviceType>();
          DistanceFunctor<LO, decltype(ghostedCoordsView)> distFunctor(ghostedCoordsView);

          // Construct Laplacian diagonal
          RCP<Vector> localLaplDiag;
          {
            SubFactoryMonitor m2(*this, "Local Laplacian diag construction", currentLevel);

            localLaplDiag = VectorFactory::Build(uniqueMap);

            auto localLaplDiagView = localLaplDiag->template getLocalView<DeviceType>();
            auto kokkosGraph = kokkosMatrix.graph;

            Kokkos::parallel_for("MueLu:CoalesceDropF:Build:scalar_filter:laplacian_diag", range_type(0,numRows),
              KOKKOS_LAMBDA(const LO row) {
                auto rowView = kokkosGraph.rowConst(row);
                auto length  = rowView.length;

                SC d = ATS::zero();
                for (decltype(length) colID = 0; colID < length; colID++) {
                  auto col = rowView(colID);
                  if (row != col)
                    d += ATS::one()/distFunctor.distance2(row, col);
                }
                localLaplDiagView(row,0) = d;
              });
          }

          // Construct ghosted Laplacian diagonal
          RCP<Vector> ghostedLaplDiag;
          {
            SubFactoryMonitor m2(*this, "Ghosted Laplacian diag construction", currentLevel);
            ghostedLaplDiag = VectorFactory::Build(nonUniqueMap);
            ghostedLaplDiag->doImport(*localLaplDiag, *importer, Xpetra::INSERT);
          }

          // Filter out entries
          {
            SubFactoryMonitor m2(*this, "MainLoop", currentLevel);

            auto ghostedLaplDiagView = ghostedLaplDiag->template getLocalView<DeviceType>();

            DistanceLaplacianDropFunctor<LO, decltype(ghostedLaplDiagView), decltype(distFunctor)>
                dropFunctor(ghostedLaplDiagView, distFunctor, threshold);
            ScalarFunctor<SC, LO, local_matrix_type, decltype(bndNodes), decltype(dropFunctor)>
                scalarFunctor(kokkosMatrix, bndNodes, dropFunctor, rows, colsAux, valsAux, reuseGraph, lumping, threshold);

            Kokkos::parallel_reduce("MueLu:CoalesceDropF:Build:scalar_filter:main_loop", range_type(0,numRows),
                                    scalarFunctor, nnzFA);
          }
        }

      }
      numDropped = nnzA - nnzFA;

      boundaryNodes = bndNodes;

      {
        SubFactoryMonitor m2(*this, "CompressRows", currentLevel);

        // parallel_scan (exclusive)
        Kokkos::parallel_scan("MueLu:CoalesceDropF:Build:scalar_filter:compress_rows", range_type(0,numRows+1),
          KOKKOS_LAMBDA(const LO i, LO& update, const bool& final_pass) {
            update += rows(i);
            if (final_pass)
              rows(i) = update;
          });
      }

      // Compress cols (and optionally vals)
      // We use a trick here: we moved all remaining elements to the beginning
      // of the original row in the main loop, so we don't need to check for
      // INVALID here, and just stop when achieving the new number of elements
      // per row.
      cols_type cols(Kokkos::ViewAllocateWithoutInitializing("FA_cols"), nnzFA);
      vals_type vals;
      if (reuseGraph) {
        // Only compress cols
        SubFactoryMonitor m2(*this, "CompressColsAndVals", currentLevel);

        Kokkos::parallel_for("MueLu:TentativePF:Build:compress_cols", range_type(0,numRows),
          KOKKOS_LAMBDA(const LO i) {
            // Is there Kokkos memcpy?
            LO rowStart  = rows(i);
            LO rowAStart = rowsA(i);
            size_t rownnz = rows(i+1) - rows(i);
            for (size_t j = 0; j < rownnz; j++)
              cols(rowStart+j) = colsAux(rowAStart+j);
          });
      } else {
        // Compress cols and vals
        SubFactoryMonitor m2(*this, "CompressColsAndVals", currentLevel);

        vals = vals_type(Kokkos::ViewAllocateWithoutInitializing("FA_vals"), nnzFA);

        Kokkos::parallel_for("MueLu:TentativePF:Build:compress_cols", range_type(0,numRows),
          KOKKOS_LAMBDA(const LO i) {
            LO rowStart  = rows(i);
            LO rowAStart = rowsA(i);
            size_t rownnz = rows(i+1) - rows(i);
            for (size_t j = 0; j < rownnz; j++) {
              cols(rowStart+j) = colsAux(rowAStart+j);
              vals(rowStart+j) = colsAux(rowAStart+j);
            }
          });
      }

      kokkos_graph_type kokkosGraph(cols, rows);

      {
        SubFactoryMonitor m2(*this, "LWGraph construction", currentLevel);

        graph = rcp(new LWGraph_kokkos(kokkosGraph, A->getRowMap(), A->getColMap(), "filtered graph of A"));
        graph->SetBoundaryNodeMap(boundaryNodes);
      }

      numTotal = A->getNodeNumEntries();

      dofsPerNode = 1;

      if (!reuseGraph) {
        SubFactoryMonitor m2(*this, "LocalMatrix+FillComplete", currentLevel);

        local_matrix_type localFA = local_matrix_type("A", numRows, kokkosMatrix.numCols(), nnzFA, vals, rows, cols);
#if 1
        // FIXME: this should be gone once Xpetra propagate "local matrix + 4 maps" constructor
        // FIXME: we don't need to rebuild importers/exporters. We should reuse A's
        auto filteredACrs = CrsMatrixFactory::Build(A->getRowMap(), A->getColMap(), localFA);
        filteredACrs->resumeFill();  // we need that for rectangular matrices
        filteredACrs->expertStaticFillComplete(A->getDomainMap(), A->getRangeMap());
#endif
        filteredA = rcp(new CrsMatrixWrap(filteredACrs));
      }

      filteredA->SetFixedBlockSize(A->GetFixedBlockSize());

      if (pL.get<bool>("filtered matrix: reuse eigenvalue")) {
        // Reuse max eigenvalue from A
        // It is unclear what eigenvalue is the best for the smoothing, but we already may have
        // the D^{-1}A estimate in A, may as well use it.
        // NOTE: ML does that too
        filteredA->SetMaxEigenvalueEstimate(A->GetMaxEigenvalueEstimate());
      } else {
        filteredA->SetMaxEigenvalueEstimate(-Teuchos::ScalarTraits<SC>::one());
      }

    } else if (blkSize > 1 && threshold == zero) {
      // Case 3:  block problem without filtering
      //
      // FIXME_KOKKOS: this code is completely unoptimized. It really should do
      // a very simple thing: merge rows and produce nodal graph. But the code
      // seems very complicated. Can we do better?

      TEUCHOS_TEST_FOR_EXCEPTION(A->getRowMap()->getNodeNumElements() % blkSize != 0, MueLu::Exceptions::RuntimeError, "MueLu::CoalesceDropFactory: Number of local elements is " << A->getRowMap()->getNodeNumElements() << " but should be a multiply of " << blkSize);

      const RCP<const Map> rowMap = A->getRowMap();
      const RCP<const Map> colMap = A->getColMap();

      // build a node row map (uniqueMap = non-overlapping) and a node column map
      // (nonUniqueMap = overlapping). The arrays rowTranslation and colTranslation
      // stored in the AmalgamationInfo class container contain the local node id
      // given a local dof id. The data is calculated in the AmalgamationFactory and
      // stored in the variable "UnAmalgamationInfo" (which is of type AmalagamationInfo)
      const RCP<const Map> uniqueMap = amalInfo->getNodeRowMap();
      const RCP<const Map> nonUniqueMap = amalInfo->getNodeColMap();
      Array<LO> rowTranslationArray = *(amalInfo->getRowTranslation()); // TAW should be transform that into a View?
      Array<LO> colTranslationArray = *(amalInfo->getColTranslation());

      // get number of local nodes
      LO numNodes = Teuchos::as<LocalOrdinal>(uniqueMap->getNodeNumElements());
      typedef typename Kokkos::View<LocalOrdinal*, DeviceType> id_translation_type;
      id_translation_type rowTranslation("dofId2nodeId",rowTranslationArray.size());
      id_translation_type colTranslation("ov_dofId2nodeId",colTranslationArray.size());

      // TODO change this to lambdas
      for (decltype(rowTranslationArray.size()) i = 0; i < rowTranslationArray.size(); ++i)
        rowTranslation(i) = rowTranslationArray[i];
      for (decltype(colTranslationArray.size()) i = 0; i < colTranslationArray.size(); ++i)
        colTranslation(i) = colTranslationArray[i];
      // extract striding information
      blkSize = A->GetFixedBlockSize();  //< the full block size (number of dofs per node in strided map)
      LocalOrdinal blkId   = -1;         //< the block id within a strided map or -1 if it is a full block map
      LocalOrdinal blkPartSize = A->GetFixedBlockSize(); //< stores block size of part blkId (or the full block size)
      if(A->IsView("stridedMaps") == true) {
        const RCP<const Map> myMap = A->getRowMap("stridedMaps");
        const RCP<const StridedMap> strMap = Teuchos::rcp_dynamic_cast<const StridedMap>(myMap);
        TEUCHOS_TEST_FOR_EXCEPTION(strMap.is_null() == true, Exceptions::RuntimeError, "Map is not of type stridedMap");
        blkSize = Teuchos::as<const LocalOrdinal>(strMap->getFixedBlockSize());
        blkId   = strMap->getStridedBlockId();
        if (blkId > -1)
          blkPartSize = Teuchos::as<LocalOrdinal>(strMap->getStridingData()[blkId]);
      }
      auto kokkosMatrix = A->getLocalMatrix(); // access underlying kokkos data

      //
      typedef typename Matrix::local_matrix_type          local_matrix_type;
      typedef typename LWGraph_kokkos::local_graph_type   kokkos_graph_type;
      typedef typename kokkos_graph_type::row_map_type    row_map_type;
      //typedef typename row_map_type::HostMirror           row_map_type_h;
      typedef typename kokkos_graph_type::entries_type    entries_type;

      // Stage 1c: get number of dof-nonzeros per blkSize node rows
      typename row_map_type::non_const_type dofNnz("nnz_map", numNodes + 1);
      LO numDofCols = 0;
      Stage1aVectorFunctor<decltype(kokkosMatrix), decltype(dofNnz), decltype(blkPartSize)> stage1aFunctor(kokkosMatrix, dofNnz, blkPartSize);
      Kokkos::parallel_reduce("MueLu:CoalesceDropF:Build:scalar_filter:stage1a", range_type(0,numNodes), stage1aFunctor, numDofCols);
      // parallel_scan (exclusive)
      ScanFunctor<LO,decltype(dofNnz)> scanFunctor(dofNnz);
      Kokkos::parallel_scan("MueLu:CoalesceDropF:Build:scalar_filter:stage1_scan", range_type(0,numNodes+1), scanFunctor);

      typename entries_type::non_const_type dofcols("dofcols", numDofCols/*dofNnz(numNodes)*/); // why does dofNnz(numNodes) work? should be a parallel reduce, i guess
      Stage1bVectorFunctor <decltype(kokkosMatrix), decltype(dofNnz), decltype(blkPartSize), decltype(dofcols)> stage1bFunctor(kokkosMatrix, dofNnz, blkPartSize, dofcols);
      Kokkos::parallel_for("MueLu:CoalesceDropF:Build:scalar_filter:stage1b", range_type(0,numNodes), stage1bFunctor);

      // we have dofcols and dofids from Stage1dVectorFunctor
      LO numNodeCols = 0;
      typename row_map_type::non_const_type rows("nnz_nodemap", numNodes + 1);
      typename boundary_nodes_type::non_const_type bndNodes("boundaryNodes", numNodes);
      Stage1cVectorFunctor <decltype(kokkosMatrix), decltype(dofNnz), decltype(dofcols), decltype(colTranslation), decltype(bndNodes)> stage1cFunctor(dofNnz, dofcols, colTranslation,rows,bndNodes);
      Kokkos::parallel_reduce("MueLu:CoalesceDropF:Build:scalar_filter:stage1c", range_type(0,numNodes), stage1cFunctor,numNodeCols);

      // parallel_scan (exclusive)
      ScanFunctor<LO,decltype(rows)> scanNodeFunctor(rows);
      Kokkos::parallel_scan("MueLu:CoalesceDropF:Build:scalar_filter:stage1_scan", range_type(0,numNodes+1), scanNodeFunctor);

      // create column node view
      typename entries_type::non_const_type cols("nodecols", numNodeCols);


      Stage1dVectorFunctor <decltype(kokkosMatrix), decltype(dofNnz), decltype(dofcols), decltype(rows), decltype(cols)> stage1dFunctor(dofcols, dofNnz, cols, rows);
      Kokkos::parallel_for("MueLu:CoalesceDropF:Build:scalar_filter:stage1c", range_type(0,numNodes), stage1dFunctor);
      kokkos_graph_type kokkosGraph(cols, rows);

      // create LW graph
      graph = rcp(new LWGraph_kokkos(kokkosGraph, uniqueMap, nonUniqueMap, "amalgamated graph of A"));

      boundaryNodes = bndNodes;
      graph->SetBoundaryNodeMap(boundaryNodes);
      numTotal = A->getNodeNumEntries();

      dofsPerNode = blkSize;

      filteredA = A;

    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu: CoalesceDropFactory_kokkos: Block filtering is not implemented");
    }

    if (GetVerbLevel() & Statistics1) {
      GO numLocalBoundaryNodes  = 0;
      GO numGlobalBoundaryNodes = 0;

      Kokkos::parallel_reduce("MueLu:CoalesceDropF:Build:bnd", range_type(0, boundaryNodes.extent(0)),
        KOKKOS_LAMBDA(const LO i, GO& n) {
          if (boundaryNodes(i))
            n++;
        }, numLocalBoundaryNodes);

      auto comm = A->getRowMap()->getComm();
      MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
      GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
    }

    if ((GetVerbLevel() & Statistics1) && threshold != zero) {
      auto comm = A->getRowMap()->getComm();

      GO numGlobalTotal, numGlobalDropped;
      MueLu_sumAll(comm, numTotal,   numGlobalTotal);
      MueLu_sumAll(comm, numDropped, numGlobalDropped);

      if (numGlobalTotal != 0) {
        GetOStream(Statistics1) << "Number of dropped entries: "
            << numGlobalDropped << "/" << numGlobalTotal
            << " (" << 100*Teuchos::as<double>(numGlobalDropped)/Teuchos::as<double>(numGlobalTotal) << "%)" << std::endl;
      }
    }

    Set(currentLevel, "DofsPerNode",  dofsPerNode);
    Set(currentLevel, "Graph",        graph);
    Set(currentLevel, "A",            filteredA);
  }
}
#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_COALESCEDROPFACTORY_KOKKOS_DEF_HPP
