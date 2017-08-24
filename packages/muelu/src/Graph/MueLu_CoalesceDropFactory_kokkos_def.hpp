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

    template<class MatrixType, class GhostedDiagType, class RowType>
    class Stage1ScalarFunctor {
    private:
      typedef typename MatrixType::ordinal_type LO;
      typedef typename MatrixType::value_type   SC;
      typedef Kokkos::ArithTraits<SC>           ATS;
      typedef typename ATS::magnitudeType       magnitudeType;

    public:
      Stage1ScalarFunctor(MatrixType kokkosMatrix_, double threshold_, GhostedDiagType ghostedDiag_, RowType rows_) :
        kokkosMatrix(kokkosMatrix_),
        threshold(threshold_),
        ghostedDiag(ghostedDiag_),
        rows(rows_)
      { }

      KOKKOS_INLINE_FUNCTION
      void operator()(const LO row, LO& nnz) const {
        auto rowView = kokkosMatrix.row (row);
        auto length  = rowView.length;

        LO rownnz = 0;
        for (decltype(length) colID = 0; colID < length; colID++) {
          LO col = rowView.colidx(colID);

          // Avoid square root by using squared values
          magnitudeType aiiajj = threshold*threshold * ATS::magnitude(ghostedDiag(row, 0))*ATS::magnitude(ghostedDiag(col, 0));   // eps^2*|a_ii|*|a_jj|
          magnitudeType aij2   = ATS::magnitude(rowView.value(colID)) * ATS::magnitude(rowView.value(colID));                     // |a_ij|^2

          if (aij2 > aiiajj || row == col)
            rownnz++;
        }
        rows(row+1) = rownnz;
        nnz += rownnz;
      }

    private:
      MatrixType        kokkosMatrix;
      double            threshold;
      GhostedDiagType   ghostedDiag;
      RowType           rows;
    };

    template<class MatrixType, class GhostedDiagType, class RowType, class ColType, class BndNodesType>
    class Stage2ScalarFunctor {
    private:
      typedef typename MatrixType::ordinal_type LO;
      typedef typename MatrixType::size_type    GO;
      typedef typename MatrixType::value_type   SC;
      typedef Kokkos::ArithTraits<SC>           ATS;
      typedef typename ATS::magnitudeType       magnitudeType;

    public:
      Stage2ScalarFunctor(MatrixType kokkosMatrix_, GhostedDiagType ghostedDiag_, RowType rows_, ColType cols_, BndNodesType bndNodes_, double threshold_) :
        kokkosMatrix(kokkosMatrix_),
        ghostedDiag(ghostedDiag_),
        rows(rows_),
        cols(cols_),
        bndNodes(bndNodes_),
        threshold(threshold_)
      { }

      KOKKOS_INLINE_FUNCTION
      void operator()(const LO row, GO& dropped) const {
        auto rowView = kokkosMatrix.row (row);
        auto length = rowView.length;

        LO rownnz = 0;
        for (decltype(length) colID = 0; colID < length; colID++) {
          LO col = rowView.colidx(colID);

          // Avoid square root by using squared values
          magnitudeType aiiajj = threshold*threshold * ATS::magnitude(ghostedDiag(row, 0))*ATS::magnitude(ghostedDiag(col, 0));   // eps^2*|a_ii|*|a_jj|
          magnitudeType aij2   = ATS::magnitude(rowView.value(colID)) * ATS::magnitude(rowView.value(colID));                     // |a_ij|^2

          if (aij2 > aiiajj || row == col) {
            cols(rows(row) + rownnz) = col;
            rownnz++;
          } else {
            dropped++;
          }
          if (rownnz == 1) {
            // If the only element remaining after filtering is diagonal, mark node as boundary
            // FIXME: this should really be replaced by the following
            //    if (indices.size() == 1 && indices[0] == row)
            //        boundaryNodes[row] = true;
            // We do not do it this way now because there is no framework for distinguishing isolated
            // and boundary nodes in the aggregation algorithms
            bndNodes(row) = true;
          }
        }
      }

    private:
      MatrixType        kokkosMatrix;
      GhostedDiagType   ghostedDiag;
      RowType           rows;
      ColType           cols;
      BndNodesType      bndNodes;
      double            threshold;
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
  void CoalesceDropFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);

    typedef Teuchos::ScalarTraits<SC> STS;
    SC zero = STS::zero();

    auto A         = Get< RCP<Matrix> >(currentLevel, "A");
    LO   blkSize   = A->GetFixedBlockSize();
    GO   indexBase = A->getRowMap()->getIndexBase();

    auto amalInfo = Get< RCP<AmalgamationInfo> >(currentLevel, "UnAmalgamationInfo");

    const ParameterList& pL = GetParameterList();

    // bool doLightWeightWrap = pL.get<bool>("lightweight wrap");
    GetOStream(Warnings0) << "lightweight wrap is deprecated" << std::endl;

    std::string algo = pL.get<std::string>("aggregation: drop scheme");

    double threshold = pL.get<double>("aggregation: drop tol");
    GetOStream(Runtime0) << "algorithm = \"" << algo << "\": threshold = " << threshold << ", blocksize = " << A->GetFixedBlockSize() << std::endl;

    Set(currentLevel, "Filtering", (threshold != zero));

    const typename STS::magnitudeType dirichletThreshold = STS::magnitude(as<SC>(pL.get<double>("aggregation: Dirichlet threshold")));

    GO numDropped = 0, numTotal = 0;

    RCP<LWGraph_kokkos> graph;
    LO                  dofsPerNode = -1;

    typedef typename LWGraph_kokkos::boundary_nodes_type boundary_nodes_type;
    boundary_nodes_type boundaryNodes;

    if (blkSize == 1 && threshold == zero) {
      // Case 1:  scalar problem without dropping

      graph = rcp(new LWGraph_kokkos(A->getLocalMatrix().graph, A->getRowMap(), A->getColMap(), "graph of A"));

      // Detect and record rows that correspond to Dirichlet boundary conditions
      boundaryNodes = Utilities_kokkos::DetectDirichletRows(*A, dirichletThreshold);
      graph->SetBoundaryNodeMap(boundaryNodes);
      numTotal = A->getNodeNumEntries();
      dofsPerNode = 1;
    } else if (blkSize > 1 && threshold == zero) {
      // Case 3:  block problem without filtering

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
    } else if (algo == "classical") {

      if (blkSize == 1 && threshold != zero) {
        // Case 2:  scalar problem with filtering

        RCP<Vector> ghostedDiagVec = Utilities_kokkos::GetMatrixOverlappedDiagonal(*A);
        auto        ghostedDiag    = ghostedDiagVec->template getLocalView<typename NO::device_type>();

        typedef typename Matrix::local_matrix_type          local_matrix_type;
        typedef typename LWGraph_kokkos::local_graph_type   kokkos_graph_type;
        typedef typename kokkos_graph_type::row_map_type    row_map_type;
        typedef typename kokkos_graph_type::entries_type    entries_type;

        LO   numRows      = A->getNodeNumRows();
        auto kokkosMatrix = A->getLocalMatrix();

        typedef Kokkos::ArithTraits<SC>     ATS;
        typedef typename ATS::magnitudeType magnitudeType;

        // Stage 1: calculate the number of remaining entries per row
        typename row_map_type::non_const_type rows("row_map", numRows+1);       // rows(0) = 0 automatically

        LO realnnz = 0;

        Stage1ScalarFunctor<decltype(kokkosMatrix), decltype(ghostedDiag), decltype(rows)> stage1Functor(kokkosMatrix, threshold, ghostedDiag, rows);
        Kokkos::parallel_reduce("MueLu:CoalesceDropF:Build:scalar_filter:stage1_reduce", range_type(0,numRows), stage1Functor, realnnz);

        // parallel_scan (exclusive)
        ScanFunctor<LO,decltype(rows)> scanFunctor(rows);
        Kokkos::parallel_scan("MueLu:CoalesceDropF:Build:scalar_filter:stage1_scan", range_type(0,numRows+1), scanFunctor);


        // Stage 2: fill in the column indices
        typename boundary_nodes_type::non_const_type bndNodes("boundaryNodes", numRows);
        typename entries_type::non_const_type        cols    ("entries",       realnnz);

        typename decltype(kokkosMatrix)::size_type    kokkos_numDropped;

        Stage2ScalarFunctor<decltype(kokkosMatrix), decltype(ghostedDiag), decltype(rows), decltype(cols), decltype(bndNodes)>
          stage2Functor(kokkosMatrix, ghostedDiag, rows, cols, bndNodes, threshold);
        Kokkos::parallel_reduce("MueLu:CoalesceDropF:Build:scalar_filter:stage2_reduce", range_type(0,numRows), stage2Functor, kokkos_numDropped);
        numDropped = kokkos_numDropped;

        boundaryNodes = bndNodes;

        kokkos_graph_type kokkosGraph(cols, rows);

        graph = rcp(new LWGraph_kokkos(kokkosGraph, A->getRowMap(), A->getColMap(), "filtered graph of A"));
        graph->SetBoundaryNodeMap(boundaryNodes);

        numTotal = A->getNodeNumEntries();

        dofsPerNode = 1;

      } else if (blkSize > 1 && threshold != zero) {
        // Case 4:  block problem with filtering

        throw Exceptions::RuntimeError("Block systems with filtering are not implemented");

        // Detect and record rows that correspond to Dirichlet boundary conditions
        // boundary_nodes_type pointBoundaryNodes = Utilities_kokkos::DetectDirichletRows(*A, dirichletThreshold);
      }

    } else if (algo == "distance laplacian") {
      typedef Xpetra::MultiVector<double,LO,GO,NO> doubleMultiVector;

      auto coords = Get<RCP<doubleMultiVector> >(currentLevel, "Coordinates");

      if (blkSize == 1 && threshold != zero) {
        // Case 2:  scalar problem with filtering

        throw Exceptions::RuntimeError("not implemented");

      } else if (blkSize > 1 && threshold != zero) {
        // Case 4:  block problem with filtering

        throw Exceptions::RuntimeError("Block systems with filtering are not implemented");
      }
    }
    if (GetVerbLevel() & Statistics1) {
      GO numLocalBoundaryNodes  = 0;
      GO numGlobalBoundaryNodes = 0;
#if 0
      // Convert to functors later
      Kokkos::parallel_reduce("MueLu:CoalesceDropF:Build:bnd", range_type(0,boundaryNodes.dimension_0()), KOKKOS_LAMBDA(const LO i, GO& n) {
        if (boundaryNodes(i))
          n++;
      }, numLocalBoundaryNodes);
#endif

      RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
      MueLu_sumAll(comm, numLocalBoundaryNodes, numGlobalBoundaryNodes);
      GetOStream(Statistics1) << "Detected " << numGlobalBoundaryNodes << " Dirichlet nodes" << std::endl;
    }
    if ((GetVerbLevel() & Statistics1) && threshold != zero) {
      RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
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
  }
}
#endif // HAVE_MUELU_KOKKOS_REFACTOR
#endif // MUELU_COALESCEDROPFACTORY_KOKKOS_DEF_HPP
