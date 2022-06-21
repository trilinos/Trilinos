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

#include <Xpetra_IO.hpp>

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include "MueLu_Exceptions.hpp"
#include <MueLu_CreateXpetraPreconditioner.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosMueLuAdapter.hpp>
#endif

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

namespace Tpetra {


  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  removeSmallEntries(Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A,
                     typename Teuchos::ScalarTraits<Scalar>::magnitudeType tol) {

    using crs_matrix = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using row_ptr_type = typename crs_matrix::local_graph_device_type::row_map_type::non_const_type;
    using col_idx_type = typename crs_matrix::local_graph_device_type::entries_type::non_const_type;
    using vals_type = typename crs_matrix::local_matrix_device_type::values_type;

    auto lclA = A->getLocalMatrixHost();

    auto rowptr = row_ptr_type("rowptr", lclA.numRows()+1);

    for (LocalOrdinal rlid = 0; rlid < lclA.numRows(); ++rlid) {
      auto row = lclA.row(rlid);
      for (size_t k = 0; k<row.length; ++k) {
        if (Teuchos::ScalarTraits<Scalar>::magnitude(row.value(k)) > tol) {
          rowptr(rlid+1) += 1;
        }
      }
    }
    for (LocalOrdinal rlid = 0; rlid < lclA.numRows(); ++rlid) {
      rowptr(rlid+1) += rowptr(rlid);
    }

    auto nnz = rowptr(lclA.numRows());
    auto idx = col_idx_type("idx", nnz);
    auto vals = vals_type("vals", nnz);

    for (LocalOrdinal rlid = 0; rlid < lclA.numRows(); ++rlid) {
      auto row = lclA.row(rlid);
      auto I = rowptr(rlid);
      for (size_t k = 0; k<row.length; ++k) {
        if (Teuchos::ScalarTraits<Scalar>::magnitude(row.value(k)) > tol) {
          idx(I) = row.colidx(k);
          vals(I) = row.value(k);
          I += 1;
        }
      }
    }

    auto newA =  Teuchos::rcp(new crs_matrix(A->getRowMap(), A->getColMap(), rowptr, idx, vals));
    newA->fillComplete(A->getDomainMap(),
                       A->getRangeMap());
    return newA;
  }


  /*
    A standard map that is broken up into blocks of consecutive LIDs.
   */

  template <class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class BlockedMap {

  public:
    using map_type = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
    using lo_vec_type = Tpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node>;

    BlockedMap(const RCP<const map_type>& pointMap,
               const RCP<lo_vec_type>& blockSizes)
      :
      pointMap_(pointMap),
      blockMap_(blockSizes->getMap()),
      blockSizes_(blockSizes)
    {
      auto lclBlockSizes = blockSizes_->getLocalViewHost(Tpetra::Access::ReadOnly);
      LocalOrdinal minClusterSize = Teuchos::OrdinalTraits<LocalOrdinal>::max();
      LocalOrdinal maxClusterSize = 0;
      offsets_ = Kokkos::View<size_t*>("offsets", blockMap_->getLocalNumElements()+1);
      offsets_(0) = 0;
      for (size_t blockNum = 0; blockNum < blockMap_->getLocalNumElements(); ++blockNum) {
        offsets_(blockNum+1) = offsets_(blockNum) + lclBlockSizes(blockNum, 0);
        minClusterSize = std::min(minClusterSize, lclBlockSizes(blockNum, 0));
        maxClusterSize = std::max(maxClusterSize, lclBlockSizes(blockNum, 0));
      }
      TEUCHOS_ASSERT_EQUALITY(offsets_(blockMap_->getLocalNumElements()), pointMap->getLocalNumElements());
      minClusterSize_ = minClusterSize;
      maxClusterSize_ = maxClusterSize;
    }

    // private:
    RCP<const map_type> pointMap_;
    RCP<const map_type> blockMap_;
    RCP<lo_vec_type>    blockSizes_;
    Kokkos::View<size_t*> offsets_;
    LocalOrdinal minClusterSize_;
    LocalOrdinal maxClusterSize_;

  };

  template <class Scalar = Tpetra::Operator<>::scalar_type,
            class LocalOrdinal = typename Tpetra::Operator<Scalar>::local_ordinal_type,
            class GlobalOrdinal = typename Tpetra::Operator<Scalar, LocalOrdinal>::global_ordinal_type,
            class Node = typename Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
  class BlockedMatrix {

  public:
    using matrix_type = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using blocked_map_type = BlockedMap<LocalOrdinal,GlobalOrdinal,Node>;
    using lo_vec_type = typename blocked_map_type::lo_vec_type;

    BlockedMatrix(const RCP<matrix_type>& pointA,
                  const RCP<matrix_type>& blockA,
                  const RCP<blocked_map_type>& blockMap,
                  const RCP<blocked_map_type>& ghosted_blockMap=Teuchos::null)
      :
      pointA_(pointA),
      blockA_(blockA),
      blockMap_(blockMap),
      ghosted_blockMap_(ghosted_blockMap)
    {
      TEUCHOS_ASSERT(blockA_->getDomainMap()->isSameAs(*blockA_->getRangeMap()));
      TEUCHOS_ASSERT(blockA_->getDomainMap()->isSameAs(*blockA_->getRowMap()));
      TEUCHOS_ASSERT(blockA_->getDomainMap()->isSameAs(*blockMap_->blockMap_));

      TEUCHOS_ASSERT(pointA_->getDomainMap()->isSameAs(*pointA_->getRangeMap()));
      TEUCHOS_ASSERT(pointA_->getDomainMap()->isSameAs(*pointA_->getRowMap()));
      TEUCHOS_ASSERT(pointA_->getDomainMap()->isSameAs(*blockMap_->pointMap_));

      // {
      // auto lcl_block = blockA_->getLocalMatrixHost();
      // auto lcl_point = pointA_->getLocalMatrixHost();
      // auto lcl_offsets = blockMap_->offsets_;
      // auto lcl_ghosted_offsets = ghosted_blockMap_->offsets_;

      // for (LocalOrdinal brlid = 0; brlid < ; ++brlid) {
      //   auto brow = lcl_block.row(brlid);
      //   for (LocalOrdinal n = 0; n < brow.length; ++n) {
      //     auto bclid = brow.colidx(n);

      //     const LocalOrdinal row_start = lcl_offsets(brlid);
      //     const LocalOrdinal row_end = lcl_offsets(brlid+1);
      //     const LocalOrdinal col_start = lcl_ghosted_offsets(bclid);
      //     const LocalOrdinal col_end = lcl_ghosted_offsets(bclid+1);

      //     for (LocalOrdinal rlid = row_start; rlid < row_end; ++rlid) {
      //       auto row = lcl_point.row(rlid);

      //     }
      //   }
      // }

      // }

    }

    void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
               Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
               Teuchos::ETransp mode = Teuchos::NO_TRANS,
               Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
               Scalar beta  = Teuchos::ScalarTraits<Scalar>::zero()) const {
      pointA_->apply(X, Y, mode, alpha, beta);
    }

    void localApply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                    Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                    Teuchos::ETransp mode = Teuchos::NO_TRANS,
                    Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                    Scalar beta  = Teuchos::ScalarTraits<Scalar>::zero()) const {
      pointA_->localApply(X, Y, mode, alpha, beta);
    }

  // private:
    RCP<matrix_type>      pointA_;
    RCP<matrix_type>      blockA_;
    RCP<blocked_map_type> blockMap_;
    RCP<blocked_map_type> ghosted_blockMap_;

  };

  /*
    The unknowns of the kernel approximations are collected in the clusterMap.
    For H-matrices, this is just a concatenation.
    For H2-matrices, the map also contains the intermediate clusters that might be needed in upward/downward pass.


    H = nearField
        + basisMatrix *
          ((I+transferMatrices[K-1]^T) * ... * (I+transferMatrices[0]^T)) *
          kernelApproximations *
          ((I+transferMatrices[0]) * ... * (I+transferMatrices[K-1])) *
          basisMatrix^T

   nearField and basisMatrix are CRS matrices.
   kernelApproximations and transferMatrices[.] are blocked CRS matrices

   Maps:
   map (standard): domain and range of H;
                   domain, range, row of nearField
   clusterMap (blocked map): domain, range, row, column of transferMatrices;
                             domain, range, row of kernelApproximations
   ghosted_clusterMap (blocked map): column of kernelApproximations


   For H-matrices:
   K = 0, i.e. no transfer matrices

   For H2-matrices:
   upward and downward pass in the cluster hierarchy are encoded in transfer matrices
   */

  template <class Scalar = Tpetra::Operator<>::scalar_type,
            class LocalOrdinal = typename Tpetra::Operator<Scalar>::local_ordinal_type,
            class GlobalOrdinal = typename Tpetra::Operator<Scalar, LocalOrdinal>::global_ordinal_type,
            class Node = typename Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
  class HierarchicalOperator : public Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

  public:
    using matrix_type = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using mv_type = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
    using map_type = Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>;

    //! The RowMatrix representing the base class of CrsMatrix
    using row_matrix_type = RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    using impl_scalar_type = typename row_matrix_type::impl_scalar_type;
    using mag_type = typename Kokkos::ArithTraits<impl_scalar_type>::mag_type;

    using local_inds_device_view_type =
          typename row_matrix_type::local_inds_device_view_type;
    using local_inds_host_view_type =
          typename row_matrix_type::local_inds_host_view_type;
    using nonconst_local_inds_host_view_type =
          typename row_matrix_type::nonconst_local_inds_host_view_type;

    using global_inds_device_view_type =
          typename row_matrix_type::global_inds_device_view_type;
    using global_inds_host_view_type =
          typename row_matrix_type::global_inds_host_view_type;
    using nonconst_global_inds_host_view_type =
          typename row_matrix_type::nonconst_global_inds_host_view_type;

    using values_device_view_type =
          typename row_matrix_type::values_device_view_type;
    using values_host_view_type =
          typename row_matrix_type::values_host_view_type;
    using nonconst_values_host_view_type =
          typename row_matrix_type::nonconst_values_host_view_type;

    using blocked_matrix_type = BlockedMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using blocked_map_type = BlockedMap<LocalOrdinal, GlobalOrdinal, Node>;

    //! @name Constructor/Destructor
    //@{

    //! Constructor
    HierarchicalOperator(const RCP<matrix_type>& nearField,
                         const RCP<blocked_matrix_type>& kernelApproximations,
                         const RCP<matrix_type>& basisMatrix,
                         std::vector<RCP<blocked_matrix_type> >& transferMatrices)
      :
      nearField_(nearField),
      kernelApproximations_(kernelApproximations),
      basisMatrix_(basisMatrix),
      transferMatrices_(transferMatrices)
    {
      auto map = nearField_->getDomainMap();
      clusterCoeffMap_ = basisMatrix_->getDomainMap();

      const bool doDebugChecks = true;

      if (doDebugChecks) {
        // near field matrix lives on map and is nonlocal
        TEUCHOS_ASSERT(map->isSameAs(*nearField_->getRangeMap()));
        TEUCHOS_ASSERT(map->isSameAs(*nearField_->getRowMap()));

        // basis matrix is entirely local and maps from clusterCoeffMap_ to map.
        TEUCHOS_ASSERT(map->isSameAs(*basisMatrix->getRangeMap()));
        TEUCHOS_ASSERT(map->isSameAs(*basisMatrix->getRowMap()));
        // TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*basisMatrix->getDomainMap()));
        TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*basisMatrix->getColMap()));

        // kernel approximations live on clusterCoeffMap and are nonlocal
        TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*kernelApproximations_->pointA_->getDomainMap()));
        TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*kernelApproximations_->pointA_->getRangeMap()));
        TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*kernelApproximations_->pointA_->getRowMap()));

        for (size_t i = 0; i<transferMatrices_.size(); i++) {
          // transfer matrices are entirely local, block diagonal on clusterCoeffMap
          TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*transferMatrices_[i]->pointA_->getDomainMap()));
          TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*transferMatrices_[i]->pointA_->getColMap()));
          TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*transferMatrices_[i]->pointA_->getRowMap()));
          TEUCHOS_ASSERT(clusterCoeffMap_->isSameAs(*transferMatrices_[i]->pointA_->getRangeMap()));
        }
      }

      // Set the two importers to Isend
      RCP<Teuchos::ParameterList> distParams = rcp(new Teuchos::ParameterList());
      distParams->set("Send type", "Isend");
      {
        RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > nearFieldImporter = nearField_->getGraph()->getImporter();
        nearFieldImporter->getDistributor().setParameterList(distParams);
        auto revDistor = nearFieldImporter->getDistributor().getReverse(false);
        if (!revDistor.is_null())
          revDistor->setParameterList(distParams);
      }

      {
        RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > kernelApproximationsImporter = kernelApproximations_->pointA_->getGraph()->getImporter();
        kernelApproximationsImporter->getDistributor().setParameterList(distParams);
        auto revDistor = kernelApproximationsImporter->getDistributor().getReverse(false);
        if (!revDistor.is_null())
          revDistor->setParameterList(distParams);
      }

      // Allocate memory for apply with vectors
      allocateMemory(1);
    }

    //! Returns the Tpetra::Map object associated with the domain of this operator.
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const {
      return nearField_->getDomainMap();
    }

    //! Returns the Tpetra::Map object associated with the range of this operator.
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const {
      return nearField_->getRangeMap();
    }

    //! Returns in Y the result of a Tpetra::Operator applied to a Tpetra::MultiVector X.
    /*!
      \param[in]  X - Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
      \param[out] Y -Tpetra::MultiVector of dimension NumVectors containing result.
    */
    void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
               Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
               Teuchos::ETransp mode = Teuchos::NO_TRANS,
               Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
               Scalar beta  = Teuchos::ScalarTraits<Scalar>::zero()) const {
      const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
      const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
      bool flip = true;

      allocateMemory(X.getNumVectors());

      // near field - part 1
      RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > nearFieldImporter = nearField_->getGraph()->getImporter();
      {
        if (mode == Teuchos::NO_TRANS) {
          X_colmap_->beginImport(X, *nearFieldImporter, INSERT);
        } else if (mode == Teuchos::TRANS) {
          nearField_->localApply(X, *X_colmap_, mode, alpha, zero);
          Y.scale (beta);
          Y.beginExport(*X_colmap_, *nearFieldImporter, ADD_ASSIGN);
        }
      }

      // upward pass
      {
        Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("upward pass")));

        basisMatrix_->localApply(X, *coefficients_, Teuchos::TRANS);

        for (int i = Teuchos::as<int>(transferMatrices_.size())-1; i>=0; i--)
          if (flip) {
            coefficients2_->assign(*coefficients_);
            transferMatrices_[i]->localApply(*coefficients_, *coefficients2_, Teuchos::NO_TRANS, one, one);
            flip = false;
          } else {
            coefficients_->assign(*coefficients2_);
            transferMatrices_[i]->localApply(*coefficients2_, *coefficients_, Teuchos::NO_TRANS, one, one);
            flip = true;
          }
      }

      // far field interactions - part 1
      {
        Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("far field 1")));

        RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > kernelApproximationsImporter = kernelApproximations_->pointA_->getGraph()->getImporter();
        if (flip) {
          if (mode == Teuchos::NO_TRANS) {
            coefficients_colmap_->beginImport(*coefficients_, *kernelApproximationsImporter, INSERT);
          } else if (mode == Teuchos::TRANS) {
            kernelApproximations_->localApply(*coefficients_, *coefficients_colmap_, mode, alpha);
            coefficients2_->putScalar(zero);
            coefficients2_->beginExport(*coefficients_colmap_, *kernelApproximationsImporter, ADD_ASSIGN);
          }
        } else {
          if (mode == Teuchos::NO_TRANS) {
            coefficients_colmap_->beginImport(*coefficients2_, *kernelApproximationsImporter, INSERT);
          } else if (mode == Teuchos::TRANS) {
            kernelApproximations_->localApply(*coefficients2_, *coefficients_colmap_, mode, alpha);
            coefficients_->putScalar(zero);
            coefficients_->beginExport(*coefficients_colmap_, *kernelApproximationsImporter, ADD_ASSIGN);
          }
        }
      }

      // near field - part 2
      {
        Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("near field 2")));

        if (mode == Teuchos::NO_TRANS) {
          X_colmap_->endImport(X, *nearFieldImporter, INSERT);
          nearField_->localApply(*X_colmap_, Y, mode, alpha, beta);
        } else if (mode == Teuchos::TRANS) {
          Y.endExport(*X_colmap_, *nearFieldImporter, ADD_ASSIGN);
        }
      }

      // far field interactions - part 2
      {
        Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("far field 2")));

        RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > kernelApproximationsImporter = kernelApproximations_->pointA_->getGraph()->getImporter();
        if (flip) {
          if (mode == Teuchos::NO_TRANS) {
            coefficients_colmap_->endImport(*coefficients_, *kernelApproximationsImporter, INSERT);
            kernelApproximations_->localApply(*coefficients_colmap_, *coefficients2_, mode, alpha);
          } else if (mode == Teuchos::TRANS) {
            coefficients2_->endExport(*coefficients_colmap_, *kernelApproximationsImporter, ADD_ASSIGN);
          }
        } else {
          if (mode == Teuchos::NO_TRANS) {
            coefficients_colmap_->endImport(*coefficients2_, *kernelApproximationsImporter, INSERT);
            kernelApproximations_->localApply(*coefficients_colmap_, *coefficients_, mode, alpha);
          } else if (mode == Teuchos::TRANS) {
            coefficients_->endExport(*coefficients_colmap_, *kernelApproximationsImporter, ADD_ASSIGN);
          }
        }
      }

      // downward pass
      {
        Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("downward pass")));

        for (size_t i = 0; i<transferMatrices_.size(); i++)
          if (flip) {
            coefficients_->assign(*coefficients2_);
            transferMatrices_[i]->localApply(*coefficients2_, *coefficients_, Teuchos::TRANS, one, one);
            flip = false;
          } else {
            coefficients2_->assign(*coefficients_);
            transferMatrices_[i]->localApply(*coefficients_, *coefficients2_, Teuchos::TRANS, one, one);
            flip = true;
          }
        if (flip)
          basisMatrix_->localApply(*coefficients2_, Y, Teuchos::NO_TRANS, one, one);
        else
          basisMatrix_->localApply(*coefficients_, Y, Teuchos::NO_TRANS, one, one);
      }
    }


    RCP<HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > restrict(const RCP<matrix_type>& P) {
      // H_c = P^T * H * P
      using lo_vec_type = typename blocked_map_type::lo_vec_type;
      using vec_type = typename Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
      using MagnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;

      // P^T * basisMatrix
      RCP<matrix_type> newBasisMatrix = rcp(new matrix_type(P->getDomainMap(), clusterCoeffMap_, 0));
      MatrixMatrix::Multiply(*P, true, *basisMatrix_, false, *newBasisMatrix);

      //
      auto clusterMap = kernelApproximations_->blockA_->getRowMap();
      auto clusterSizes = kernelApproximations_->blockMap_->blockSizes_;
      auto ghosted_clusterMap = kernelApproximations_->blockA_->getColMap();
      auto ghosted_clusterSizes = kernelApproximations_->ghosted_blockMap_->blockSizes_;
      auto lcl_clusterSizes = clusterSizes->getLocalViewHost(Tpetra::Access::ReadOnly);
      auto lcl_ghosted_clusterSizes = ghosted_clusterSizes->getLocalViewHost(Tpetra::Access::ReadOnly);

      // get number of unknowns associated with cluster via new basisMatrix
      RCP<vec_type> numUnknownsPerCluster = rcp(new vec_type(clusterMap, false));
      {
        auto lcl_numUnknownsPerCluster = numUnknownsPerCluster->getLocalViewHost(Tpetra::Access::OverwriteAll);
        // Compute the transpose of the newBasisMatrix.
        RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > newBasisMatrixT;
        Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer(newBasisMatrix);
        RCP<Teuchos::ParameterList> transposeParams = rcp(new Teuchos::ParameterList);
        newBasisMatrixT = transposer.createTranspose(transposeParams);

        auto rowptr = newBasisMatrixT->getLocalRowPtrsHost();
        LocalOrdinal clusterStart = 0;
        LocalOrdinal clusterEnd = 0;
        for (LocalOrdinal cluster = 0; cluster < lcl_clusterSizes.extent(0); ++cluster) {
          clusterStart = clusterEnd;
          clusterEnd += lcl_clusterSizes(cluster, 0);
          LocalOrdinal maxEntries = 0;
          for (LocalOrdinal row = clusterStart; row < clusterEnd; ++row) {
            LocalOrdinal numEntriesPerRow = rowptr(row+1)-rowptr(row);
            maxEntries = std::max(maxEntries, numEntriesPerRow);
          }
          lcl_numUnknownsPerCluster(cluster, 0) = maxEntries;
        }
        TEUCHOS_ASSERT_EQUALITY(clusterEnd+1, rowptr.extent(0));
      }

      const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
      for (int i = Teuchos::as<int>(transferMatrices_.size())-1; i>=0; i--)
        transferMatrices_[i]->blockA_->localApply(*numUnknownsPerCluster, *numUnknownsPerCluster, Teuchos::NO_TRANS, one, one);

      // get ghosted numUnknownsPerCluster
      RCP<vec_type> ghosted_numUnknownsPerCluster = rcp(new vec_type(ghosted_clusterMap, false));
      auto import = kernelApproximations_->blockA_->getCrsGraph()->getImporter();
      ghosted_numUnknownsPerCluster->doImport(*numUnknownsPerCluster, *import, Tpetra::INSERT);

      // coarse cluster pair graph
      RCP<matrix_type> newKernelBlockGraph = rcp(new matrix_type(kernelApproximations_->blockA_->getCrsGraph()));
      newKernelBlockGraph->resumeFill();
      // point entries of cluster pairs that should be moved to the near field
      RCP<matrix_type> diffKernelApprox = rcp(new matrix_type(kernelApproximations_->pointA_->getCrsGraph()));

      {
        auto lcl_BlockGraph = kernelApproximations_->blockA_->getLocalMatrixHost();
        auto lcl_newBlockGraph = newKernelBlockGraph->getLocalMatrixHost();
        auto lcl_KernelApprox = kernelApproximations_->pointA_->getLocalMatrixHost();
        auto lcl_diffKernelApprox = diffKernelApprox->getLocalMatrixHost();
        auto lcl_numUnknownsPerCluster = numUnknownsPerCluster->getLocalViewHost(Tpetra::Access::ReadOnly);
        auto lcl_ghosted_numUnknownsPerCluster = ghosted_numUnknownsPerCluster->getLocalViewHost(Tpetra::Access::ReadOnly);
        auto lcl_offsets = kernelApproximations_->blockMap_->offsets_;
        auto lcl_ghosted_offsets = kernelApproximations_->ghosted_blockMap_->offsets_;

        int dropped = 0;
        int kept = 0;
        int ignored = 0;
        // loop over cluster pairs
        for (LocalOrdinal brlid = 0; brlid < lcl_BlockGraph.numRows(); ++brlid) {
          size_t brsize = lcl_clusterSizes(brlid, 0);
          auto brow = lcl_BlockGraph.row(brlid);
          auto new_brow = lcl_newBlockGraph.row(brlid);
          for (LocalOrdinal k = 0; k < brow.length; ++k) {
            if (brow.value(k) > 0.5) {
              LocalOrdinal bclid = brow.colidx(k);
              size_t bcsize = lcl_ghosted_clusterSizes(bclid, 0);
              size_t bsize = brsize*bcsize;

              if (brsize * bcsize >= lcl_numUnknownsPerCluster(brlid, 0) * lcl_ghosted_numUnknownsPerCluster(bclid, 0)) {
                ++dropped;
                new_brow.value(k) = 0;

                const LocalOrdinal row_start = lcl_offsets(brlid);
                const LocalOrdinal row_end = lcl_offsets(brlid+1);
                const LocalOrdinal col_start = lcl_ghosted_offsets(bclid);
                const LocalOrdinal col_end = lcl_ghosted_offsets(bclid+1);
                TEUCHOS_ASSERT_EQUALITY(row_end-row_start, brsize);
                TEUCHOS_ASSERT_EQUALITY(col_end-col_start, bcsize);
                // loop over rows of kernelApproximations in pointwise indexing
                for (LocalOrdinal rlid = row_start; rlid <  row_end; ++rlid) {
                  auto diff_row = lcl_diffKernelApprox.row(rlid);
                  auto row = lcl_KernelApprox.row(rlid);
                  int removed = 0;
                  for (LocalOrdinal n = 0; n < row.length; ++n) {
                    if ((col_start <= row.colidx(n)) && (col_end > row.colidx(n))) {
                      diff_row.value(n) = row.value(n);
                      ++removed;
                    }
                  }
                  if (removed != bcsize) {
                    std::ostringstream oss;
                    oss << "brlid " << brlid << " row " << rlid << std::endl;
                    oss << "col_start "<< col_start << " col_end " << col_end << std::endl;
                    for (LocalOrdinal n = 0; n < row.length; ++n) {
                      oss << row.colidx(n) << " " << row.value(n) << std::endl;
                    }
                    std::cout << oss.str();
                  }
                  TEUCHOS_ASSERT_EQUALITY(removed, bcsize);
                }

              } else {
                ++kept;
                new_brow.value(k) = brow.value(k);
              }
            } else {
              ++ignored;
              new_brow.value(k) = brow.value(k);
            }
          }
        }

        // std::cout << "dropped " << dropped << " kept " << kept << " ignored " << ignored << std::endl;
      }
      newKernelBlockGraph->fillComplete(kernelApproximations_->blockA_->getDomainMap(),
                                        kernelApproximations_->blockA_->getRangeMap());
      newKernelBlockGraph = removeSmallEntries(newKernelBlockGraph, Teuchos::ScalarTraits<MagnitudeType>::eps());
      diffKernelApprox->fillComplete(clusterCoeffMap_,
                                     clusterCoeffMap_);

      std::vector<RCP<blocked_matrix_type> > newTransferMatrices;
      {
        const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
        RCP<mv_type> mv_temp = rcp(new mv_type(newKernelBlockGraph->getDomainMap(), 1));
        mv_temp->putScalar(one);
        RCP<mv_type> clusterUseCount = rcp(new mv_type(newKernelBlockGraph->getDomainMap(), 1));
        clusterUseCount->putScalar(zero);
        newKernelBlockGraph->apply(*mv_temp, *clusterUseCount, Teuchos::NO_TRANS);
        newKernelBlockGraph->apply(*mv_temp, *clusterUseCount, Teuchos::TRANS, one, one);

        std::vector<int> keepTransfers;
        auto comm = getComm();
        // std::ostringstream oss;
        for (int i = Teuchos::as<int>(transferMatrices_.size())-1; i>=0; i--) {
          transferMatrices_[i]->blockA_->localApply(*clusterUseCount, *mv_temp, Teuchos::TRANS);
          auto lcl_counts = mv_temp->getLocalViewHost(Tpetra::Access::ReadOnly);
          Scalar lcl_count = zero, gbl_count = zero;
          for (LocalOrdinal n = 0; n < lcl_counts.extent(0); ++n) {
            lcl_count += lcl_counts(n, 0);
          }
          Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &lcl_count, &gbl_count);
          // oss << "Transfer " << i << " count " << gbl_count << std::endl;

          if (gbl_count < 0.5) {
            RCP<matrix_type> temp = rcp(new matrix_type(P->getDomainMap(), clusterCoeffMap_, 0));
            MatrixMatrix::Multiply(*newBasisMatrix, false, *transferMatrices_[i]->pointA_, true, *temp);
            newBasisMatrix = temp;
          } else
            keepTransfers.push_back(i);
        }
        // std::cout << oss.str();

        for (auto it = keepTransfers.begin(); it != keepTransfers.end(); ++it) {
          newTransferMatrices.insert(newTransferMatrices.begin(), transferMatrices_[*it]);
        }
      }

      // coarse point matrix of cluster pairs
      RCP<matrix_type> newKernelApprox;
      {
        RCP<matrix_type> temp = MatrixMatrix::add(one, false, *kernelApproximations_->pointA_, -one, false, *diffKernelApprox);
        newKernelApprox = removeSmallEntries(temp, Teuchos::ScalarTraits<MagnitudeType>::eps());
      }

      // static int lvlNo = 1;
      // Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::writeSparseFile("kernel"+std::to_string(lvlNo), *newKernelApprox);
      // Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::writeSparseFile("diffKernel"+std::to_string(lvlNo), *diffKernelApprox);
      // Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::writeSparseFile("kernelGraph"+std::to_string(lvlNo), *newKernelBlockGraph);
      // Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::writeDenseFile("numUnknownsPerCluster"+std::to_string(lvlNo), numUnknownsPerCluster);
      // ++lvlNo;

      RCP<blocked_matrix_type> newBlockedKernelApproximation = rcp(new blocked_matrix_type(newKernelApprox, newKernelBlockGraph, kernelApproximations_->blockMap_, kernelApproximations_->ghosted_blockMap_));

      // construct identity on clusterCoeffMap_
      RCP<matrix_type> identity = rcp(new matrix_type(clusterCoeffMap_, 1));
      Teuchos::ArrayView<const GlobalOrdinal> gblRows = clusterCoeffMap_->getLocalElementList ();
      for (auto it = gblRows.begin (); it != gblRows.end (); ++it) {
        Teuchos::Array<GlobalOrdinal> col (1, *it);
        Teuchos::Array<Scalar> val (1, one);
        identity->insertGlobalValues (*it, col (), val ());
      }
      identity->fillComplete ();

      // transfer = basisMatrix_ * (identity + transferMatrices_[0]) * ... * (identity + transferMatrices_[n-1])
      RCP<matrix_type> transfer = rcp(new matrix_type(*newBasisMatrix));
      for (size_t i = 0; i<transferMatrices_.size(); i++) {
        RCP<matrix_type> temp = MatrixMatrix::add(one, false, *identity, one, false, *transferMatrices_[i]->pointA_);
        RCP<matrix_type> temp2 = rcp(new matrix_type(newBasisMatrix->getRowMap(), 0));
        MatrixMatrix::Multiply(*transfer, false, *temp, true, *temp2);
        transfer = temp2;
      }

      // diffFarField = transfer * diffKernelApprox * transfer^T
      RCP<matrix_type> diffFarField;
      {
        RCP<matrix_type> temp = rcp(new matrix_type(newBasisMatrix->getRowMap(), 0));
        MatrixMatrix::Multiply(*transfer, false, *diffKernelApprox, false, *temp);
        diffFarField = rcp(new matrix_type(newBasisMatrix->getRowMap(), 0));
        MatrixMatrix::Multiply(*temp, false, *transfer, true, *diffFarField);
      }

      // P^T * nearField * P
      RCP<matrix_type> newNearField;
      {
        RCP<matrix_type> temp = rcp(new matrix_type(nearField_->getRowMap(), 0));
        MatrixMatrix::Multiply(*nearField_, false, *P, false, *temp);
        RCP<matrix_type> temp2 = rcp(new matrix_type(P->getDomainMap(), 0));
        MatrixMatrix::Multiply(*P, true, *temp, false, *temp2);
        newNearField = MatrixMatrix::add(one, false, *temp2, one, false, *diffFarField);
      }

      return rcp(new HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>(newNearField, newBlockedKernelApproximation, newBasisMatrix, newTransferMatrices));
    }

    RCP<matrix_type> toMatrix() {
      const Scalar one = Teuchos::ScalarTraits<Scalar>::one();

      // construct identity on clusterCoeffMap_
      RCP<matrix_type> identity = rcp(new matrix_type(clusterCoeffMap_, 1));
      Teuchos::ArrayView<const GlobalOrdinal> gblRows = clusterCoeffMap_->getLocalElementList ();
      for (auto it = gblRows.begin (); it != gblRows.end (); ++it) {
        Teuchos::Array<GlobalOrdinal> col (1, *it);
        Teuchos::Array<Scalar> val (1, one);
        identity->insertGlobalValues (*it, col (), val ());
      }
      identity->fillComplete ();

      // transfer = basisMatrix_ * (identity + transferMatrices_[0]) * ... * (identity + transferMatrices_[n-1])
      RCP<matrix_type> transfer = rcp(new matrix_type(*basisMatrix_));
      for (size_t i = 0; i<transferMatrices_.size(); i++) {
        RCP<matrix_type> temp = MatrixMatrix::add(one, false, *identity, one, false, *transferMatrices_[i]->pointA_);
        RCP<matrix_type> temp2 = rcp(new matrix_type(basisMatrix_->getRowMap(), 0));
        MatrixMatrix::Multiply(*transfer, false, *temp, true, *temp2);
        transfer = temp2;
      }

      // farField = transfer * kernelApproximations_ * transfer^T
      RCP<matrix_type> temp = rcp(new matrix_type(basisMatrix_->getRowMap(), 0));
      MatrixMatrix::Multiply(*transfer, false, *kernelApproximations_->pointA_, false, *temp);
      RCP<matrix_type> farField = rcp(new matrix_type(basisMatrix_->getRowMap(), 0));
      MatrixMatrix::Multiply(*temp, false, *transfer, true, *farField);

      // nearField_ + farField
      return MatrixMatrix::add(one, false, *nearField_, one, false, *farField);

    }

    double getCompression() {
      size_t nnz = (nearField_->getGlobalNumEntries() +
                    kernelApproximations_->pointA_->getGlobalNumEntries() +
                    basisMatrix_->getGlobalNumEntries());
      for (size_t i = 0; i < transferMatrices_.size(); i++)
        nnz += transferMatrices_[i]->pointA_->getGlobalNumEntries();
      return Teuchos::as<double>(nnz) / (getDomainMap()->getGlobalNumElements()*getDomainMap()->getGlobalNumElements());
    }

    RCP<matrix_type> nearFieldMatrix() {
      return nearField_;
    }

    // Fake RowMatrix interface
    Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getRowMap() const {
      return nearField_->getRowMap();
    }

    Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getColMap() const {
      return nearField_->getColMap();
    }

    Teuchos::RCP<const Teuchos::Comm<int> > getComm() const {
      return nearField_->getDomainMap()->getComm();
    }

    Teuchos::RCP<const RowGraph<LocalOrdinal,GlobalOrdinal,Node> > getGraph() const {
      return nearField_->getCrsGraph();
    }

    global_size_t getGlobalNumRows() const {
      return nearField_->getGlobalNumRows();
    }

    global_size_t getGlobalNumCols() const {
      return nearField_->getGlobalNumCols();
    }

    size_t getLocalNumRows() const {
      return nearField_->getLocalNumRows();
    }

    size_t getLocalNumCols() const {
      return nearField_->getLocalNumCols();
    }

    GlobalOrdinal getIndexBase() const {
      return nearField_->getIndexBase();
    }

    global_size_t getGlobalNumEntries() const {
      return nearField_->getGlobalNumEntries();
    }

    size_t getLocalNumEntries() const {
      return nearField_->getLocalNumEntries();
    }

    size_t getNumEntriesInGlobalRow (GlobalOrdinal globalRow) const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

    size_t getNumEntriesInLocalRow (LocalOrdinal localRow) const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

    size_t getGlobalMaxNumRowEntries () const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

    size_t getLocalMaxNumRowEntries () const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

    bool hasColMap () const {
      return false;
    }

    bool isLocallyIndexed() const {
      return true;
    }

    bool isGloballyIndexed() const {
      return true;
    }

    bool isFillComplete() const {
      return true;
    }

    bool supportsRowViews() const {
      return false;
    }

    void
    getGlobalRowCopy (GlobalOrdinal GlobalRow,
                      nonconst_global_inds_host_view_type &Indices,
                      nonconst_values_host_view_type &Values,
                      size_t& NumEntries) const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

    void
    getLocalRowCopy (LocalOrdinal LocalRow,
                     nonconst_local_inds_host_view_type &Indices,
                     nonconst_values_host_view_type &Values,
                     size_t& NumEntries) const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

    void
    getGlobalRowView (GlobalOrdinal GlobalRow,
                      global_inds_host_view_type &indices,
                      values_host_view_type &values) const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

    void
    getLocalRowView (LocalOrdinal LocalRow,
                     local_inds_host_view_type & indices,
                     values_host_view_type & values) const {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

    void getLocalDiagCopy (Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &diag) const {
      nearField_->getLocalDiagCopy(diag);
    }

    void leftScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

    void rightScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
      throw MueLu::Exceptions::RuntimeError("Not implemented.");
    }

    mag_type getFrobeniusNorm() const {
      return 0.;
    }

    void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
      using std::setw;
      using std::endl;
      const size_t numRows = nearField_->getRowMap()->getGlobalNumElements();
      const size_t nnzNearField = nearField_->getGlobalNumEntries();
      const double nnzNearPerRow = Teuchos::as<double>(nnzNearField)/numRows;
      const size_t nnzKernelApprox = kernelApproximations_->pointA_->getGlobalNumEntries();
      const size_t numClusterPairs = kernelApproximations_->blockA_->getGlobalNumEntries();
      const size_t nnzBasis = basisMatrix_->getGlobalNumEntries();
      size_t nnzTransfer = 0;
      for (size_t i = 0; i<transferMatrices_.size(); i++)
        nnzTransfer += transferMatrices_[i]->pointA_->getGlobalNumEntries();
      const size_t nnzTotal = nnzNearField+nnzKernelApprox+nnzBasis+nnzTransfer;
      const double nnzTotalPerRow = Teuchos::as<double>(nnzTotal)/numRows;
      std::ostringstream oss;
      oss << std::left;
      oss << setw(9) << "rows"  << setw(12) << "nnz(near)"  << setw(14) << "nnz(near)/row" << setw(12) << "nnz(basis)" << setw(15) << "#cluster pairs" << setw(12)<< "nnz(kernel)"    << setw(14) << "nnz(transfer)" << setw(12) << "nnz(total)" << setw(14) << "nnz(total)/row" << endl;
      oss << setw(9) << numRows << setw(12) << nnzNearField << setw(14) << nnzNearPerRow   << setw(12) << nnzBasis     << setw(15) << numClusterPairs  << setw(12) << nnzKernelApprox << setw(14) << nnzTransfer     << setw(12) << nnzTotal     << setw(14) << nnzTotalPerRow   << endl;
      out << oss.str();
    }

  private:

    void allocateMemory(size_t numVectors) const {
      if (coefficients_.is_null() || coefficients_->getNumVectors() != numVectors) {
        coefficients_  = rcp(new mv_type(clusterCoeffMap_, numVectors));
        coefficients2_ = rcp(new mv_type(clusterCoeffMap_, numVectors));
        X_colmap_ = rcp(new mv_type(nearField_->getColMap(), numVectors));
        coefficients_colmap_  = rcp(new mv_type(kernelApproximations_->pointA_->getColMap(), numVectors));
      }
    }

    RCP<matrix_type> nearField_;
    RCP<blocked_matrix_type> kernelApproximations_;
    RCP<matrix_type> basisMatrix_;
    std::vector<RCP<blocked_matrix_type> > transferMatrices_;
    RCP<const map_type> clusterCoeffMap_;
    mutable RCP<mv_type> coefficients_, coefficients2_;
    mutable RCP<mv_type> X_colmap_, coefficients_colmap_;
  };

}

namespace Xpetra {

  template <class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class MyBlockedMap {

  public:
    using map_type = Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
    using lo_vec_type = Xpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node>;
    using tpetra_blocked_map_type = Tpetra::BlockedMap<LocalOrdinal,GlobalOrdinal,Node>;

    MyBlockedMap(const RCP<const map_type>& pointMap,
               const RCP<lo_vec_type>& blockSizes)
    {
      using TpLOVec = TpetraVector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node>;
      tpBlockedMap_ = rcp(new tpetra_blocked_map_type(toTpetra(pointMap),
                                                      rcp_dynamic_cast<TpLOVec>(blockSizes)->getTpetra_Vector()));
    }

    RCP<tpetra_blocked_map_type> getTpetra_BlockedMap() const {
      return tpBlockedMap_;
    }

  private:
    RCP<tpetra_blocked_map_type> tpBlockedMap_;

  };


  template <class Scalar = Tpetra::Operator<>::scalar_type,
            class LocalOrdinal = typename Tpetra::Operator<Scalar>::local_ordinal_type,
            class GlobalOrdinal = typename Tpetra::Operator<Scalar, LocalOrdinal>::global_ordinal_type,
            class Node = typename Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
  class BlockedMatrix {

  public:
    using matrix_type = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using blocked_map_type = MyBlockedMap<LocalOrdinal,GlobalOrdinal,Node>;
    using tpetra_blocked_matrix_type = Tpetra::BlockedMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;

    BlockedMatrix(const RCP<matrix_type>& pointA,
                  const RCP<matrix_type>& blockA,
                  const RCP<blocked_map_type>& blockMap,
                  const RCP<blocked_map_type>& ghosted_blockMap=Teuchos::null)
    {
      using TpCrs = TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
      using CrsWrap = CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
      RCP<typename blocked_map_type::tpetra_blocked_map_type> tp_ghosted_blockMap;
      if (!ghosted_blockMap.is_null())
        tp_ghosted_blockMap = ghosted_blockMap->getTpetra_BlockedMap();
      blockMatrix_ = rcp(new tpetra_blocked_matrix_type(rcp_dynamic_cast<TpCrs>(rcp_dynamic_cast<CrsWrap>(pointA)->getCrsMatrix(), true)->getTpetra_CrsMatrixNonConst(),
                                                        rcp_dynamic_cast<TpCrs>(rcp_dynamic_cast<CrsWrap>(blockA)->getCrsMatrix(), true)->getTpetra_CrsMatrixNonConst(),
                                                        blockMap->getTpetra_BlockedMap(),
                                                        tp_ghosted_blockMap));
    }

    void apply(const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
               Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
               Teuchos::ETransp mode = Teuchos::NO_TRANS,
               Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
               Scalar beta  = Teuchos::ScalarTraits<Scalar>::zero()) const {
      blockMatrix_->apply(Xpetra::toTpetra(X), Xpetra::toTpetra(Y), mode, alpha, beta);
    }

    void localApply(const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                    Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                    Teuchos::ETransp mode = Teuchos::NO_TRANS,
                    Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                    Scalar beta  = Teuchos::ScalarTraits<Scalar>::zero()) const {
      blockMatrix_->localApply(Xpetra::toTpetra(X), Xpetra::toTpetra(Y), mode, alpha, beta);
    }

    RCP<tpetra_blocked_matrix_type> getTpetra_BlockedMatrix() const {
      return blockMatrix_;
    }

  private:
    RCP<tpetra_blocked_matrix_type> blockMatrix_;

  };

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class HierarchicalOperator : public TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

  public:
    using tHOp = Tpetra::HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
    using map_type = Map<LocalOrdinal,GlobalOrdinal,Node>;
    using mv_type = MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
    using matrix_type = Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
    using blocked_matrix_type = BlockedMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;

    //! @name Constructor/Destructor
    //@{

    //! Constructor
    HierarchicalOperator(const RCP<tHOp>& op) : op_(op) { }

    HierarchicalOperator(const RCP<matrix_type>& nearField,
                         const RCP<blocked_matrix_type>& kernelApproximations,
                         const RCP<matrix_type>& basisMatrix,
                         std::vector<RCP<blocked_matrix_type> >& transferMatrices) {
      using TpCrs = TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
      using TpGOVec = TpetraMultiVector<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node>;
      using CrsWrap = CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>;

      std::vector<RCP<typename tHOp::blocked_matrix_type> > tTransferMatrices;
      for (size_t i = 0; i<transferMatrices.size(); i++) {
        auto transferT = transferMatrices[i]->getTpetra_BlockedMatrix();
        tTransferMatrices.push_back(transferT);
      }

      op_ = rcp(new tHOp(rcp_dynamic_cast<TpCrs>(rcp_dynamic_cast<CrsWrap>(nearField)->getCrsMatrix(), true)->getTpetra_CrsMatrixNonConst(),
                         kernelApproximations->getTpetra_BlockedMatrix(),
                         rcp_dynamic_cast<TpCrs>(rcp_dynamic_cast<CrsWrap>(basisMatrix)->getCrsMatrix(), true)->getTpetra_CrsMatrixNonConst(),
                         tTransferMatrices));
    }

    //! Returns the Tpetra::Map object associated with the domain of this operator.
    Teuchos::RCP<const map_type> getDomainMap() const {
      return toXpetra(op_->getDomainMap());
    }

    //! Returns the Tpetra::Map object associated with the range of this operator.
    Teuchos::RCP<const map_type> getRangeMap() const {
      return toXpetra(op_->getRangeMap());
    }

    //! \brief Computes the operator-multivector application.
    /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{mode}} \cdot X + \beta \cdot Y\f$. However, the details of operation
        vary according to the values of \c alpha and \c beta. Specifically
        - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
        - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
     */
    void apply (const mv_type& X, mv_type& Y,
                Teuchos::ETransp mode = Teuchos::NO_TRANS,
                Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const {
      op_->apply(Xpetra::toTpetra(X), Xpetra::toTpetra(Y), mode, alpha, beta);
    }

    //! Compute a residual R = B - (*this) * X
    void residual(const mv_type & X,
                  const mv_type & B,
                  mv_type& R) const {
      Tpetra::Details::residual(*op_, toTpetra(X), toTpetra(B), toTpetra(R));
    }

    RCP<HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > restrict(const RCP<matrix_type>& P) {
      using TpCrs = TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
      using CrsWrap = CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
      return rcp(new HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>(op_->restrict(rcp_dynamic_cast<TpCrs>(rcp_dynamic_cast<CrsWrap>(P)->getCrsMatrix(), true)->getTpetra_CrsMatrixNonConst())));
    }

    RCP<matrix_type> toMatrix() {
      using TpCrs = TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
      using CrsWrap = CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
      auto tpMat = rcp(new TpCrs(op_->toMatrix()));
      return rcp(new CrsWrap(rcp_dynamic_cast<CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(tpMat)));
    }

    double getCompression() {
      return op_->getCompression();
    }

    RCP<matrix_type> nearFieldMatrix() {
      auto tpMat = rcp(new TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(op_->nearFieldMatrix()));
      return rcp(new CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rcp_dynamic_cast<CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(tpMat)));
    }

    //! Gets the operator out
    RCP<Tpetra::Operator< Scalar, LocalOrdinal, GlobalOrdinal, Node> > getOperator() { return op_; }

    RCP<const Tpetra::Operator< Scalar, LocalOrdinal, GlobalOrdinal, Node> > getOperatorConst() const { return op_; }

    void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
      op_->describe(out, verbLevel);
    }

  private:
    RCP<tHOp> op_;
  };
}


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
    map = IO::ReadMap(hierarchicalParams.get<std::string>("map"), lib, comm);
    // colmap of near field
    near_colmap = IO::ReadMap(hierarchicalParams.get<std::string>("near colmap"), lib, comm);
    // 1-to-1 map for the cluster coefficients
    clusterCoeffMap = IO::ReadMap(hierarchicalParams.get<std::string>("coefficient map"), lib, comm);
    // overlapping map for the cluster coefficients
    ghosted_clusterCoeffMap = IO::ReadMap(hierarchicalParams.get<std::string>("ghosted coefficient map"), lib, comm);
    // 1-to-1 map for the clusters
    clusterMap = IO::ReadMap(hierarchicalParams.get<std::string>("cluster map"), lib, comm);
    // overlapping map for the clusters
    ghosted_clusterMap = IO::ReadMap(hierarchicalParams.get<std::string>("ghosted cluster map"), lib, comm);

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

    op = rcp(new HOp(nearField, blockKernelApproximations, basisMatrix, transferMatrices));

    return op;
  }

};


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
  #include "MueLu_UseShortNames.hpp"

  std::string xmlHierarchical  = "hierarchical-1d-mm-global.xml";
  std::string xmlMueLu        = "muelu.xml";
  std::string xmlAuxHierarchy = "aux.xml";
  clp.setOption("xml",    &xmlHierarchical);
  clp.setOption("xmlMueLu", &xmlMueLu);
  clp.setOption("xmlAux", &xmlAuxHierarchy);
  bool        printTimings      = true;               clp.setOption("timings", "notimings",  &printTimings,      "print timings to screen");

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

  Teuchos::ParameterList         hierarchicalParams;
  Teuchos::updateParametersFromXmlFileAndBroadcast(xmlHierarchical, Teuchos::Ptr<Teuchos::ParameterList>(&hierarchicalParams), *comm);

  RCP<HOp> op;
  {
    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("Read hierarchical matrix")));
    op = IOhelpers::Read(hierarchicalParams, comm);
  }

  out << "Compression: " << op->getCompression() << " of dense matrix."<< std::endl;

  RCP<const map_type> map = op->getDomainMap();
  RCP<matrix_type>    auxOp;
  RCP<mv_type>       X_ex, RHS, X;
  RCP<coord_mv>       coords;
  {
    // Read in auxiliary stuff

    const bool readBinary = hierarchicalParams.get<bool>("read binary", false);
    const bool readLocal = hierarchicalParams.get<bool>("read local", false);

    // colmap of auxiliary operator
    RCP<const map_type> aux_colmap = IO::ReadMap(hierarchicalParams.get<std::string>("aux colmap"), lib, comm);

    // Auxiliary matrix used for multigrid construction
    const std::string auxOpStr = hierarchicalParams.get<std::string>("auxiliary operator");
    if (auxOpStr == "near")
      auxOp = op->nearFieldMatrix();
    else
      auxOp  = IOhelpers::Read(auxOpStr, map, aux_colmap, map, map, true, readBinary, readLocal);

    X_ex = IO::ReadMultiVector(hierarchicalParams.get<std::string>("exact solution"), map);
    RHS  = IO::ReadMultiVector(hierarchicalParams.get<std::string>("right-hand side"), map);
    X    = MultiVectorFactory::Build(map, 1);

    coords = Xpetra::IO<typename Teuchos::ScalarTraits<Scalar>::coordinateType,LocalOrdinal,GlobalOrdinal,Node>::ReadMultiVector(hierarchicalParams.get<std::string>("coordinates"), map);
  }

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

#ifdef HAVE_MUELU_BELOS
  {
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
    RCP<Teuchos::ParameterList> belosList = Teuchos::parameterList();
    belosList->set("Maximum Iterations",    1000); // Maximum number of iterations allowed
    belosList->set("Convergence Tolerance", 1e-5);    // Relative convergence tolerance requested
    belosList->set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    belosList->set("Output Frequency",      1);
    belosList->set("Output Style",          Belos::Brief);

    bool set = belosProblem->setProblem();
    if (set == false) {
      throw MueLu::Exceptions::RuntimeError("ERROR:  Belos::LinearProblem failed to set up correctly!");
    }

    // Create an iterative solver manager
    Belos::SolverFactory<Scalar, MV, OP> solverFactory;
    RCP< Belos::SolverManager<Scalar, MV, OP> > solver = solverFactory.create(belosType, belosList);
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

  {
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
      auxParams.sublist("user data").set("Coordinates", coords);
      TEUCHOS_ASSERT_EQUALITY(auxParams.get("multigrid algorithm", "unsmoothed"), "unsmoothed");

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
      params.set("coarse: max size", 1);
      params.set("max levels", auxH->GetNumLevels());
      const std::string multigridAlgo = params.get("multigrid algorithm", "unsmoothed");

      op->describe(out, Teuchos::VERB_EXTREME);

      H = rcp(new Hierarchy());
      RCP<Level> lvl = H->GetLevel(0);
      lvl->Set("A", rcp_dynamic_cast<Operator>(op));
      lvl->Set("Coordinates", coords);
      for(int lvlNo = 1; lvlNo < auxH->GetNumLevels(); lvlNo++) {
        H->AddNewLevel();
        RCP<Level> auxLvl = auxH->GetLevel(lvlNo);
        // auto mgr = auxLvl->GetFactoryManager();
        // auxLvl->print(std::cout, MueLu::Debug);
        RCP<Level> fineLvl = H->GetLevel(lvlNo-1);
        lvl = H->GetLevel(lvlNo);
        auto P = auxLvl->Get<RCP<Matrix> >("P");
        auto fineA = rcp_dynamic_cast<HOp>(fineLvl->Get<RCP<Operator> >("A"));

        if (multigridAlgo == "sa") {
          // MueLu::Level lvl0, lvl1;
          // lvl0.SetFactoryManager(Teuchos::null);
          // lvl1.SetFactoryManager(Teuchos::null);
          // lvl0.setlib(P->getDomainMap()->lib());
          // lvl1.setlib(P->getDomainMap()->lib());
          // lvl0.SetLevelID(0);
          // lvl1.SetLevelID(1);
          // lvl1.SetPreviousLevel(rcpFromRef(lvl0));
          // lvl0.Set("A", fineA->nearFieldMatrix());

          // lvl1.Set("Ptent", P);
          // auto sapFact = rcp(new MueLu::SaPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>());
          // lvl1.Request("A", sapFact.get());

          // P = lvl1.Get< RCP<Matrix> >("P", sapFact.get());

        }

        lvl->Set("P", P);
        params.sublist("level "+std::to_string(lvlNo)).set("P", P);


        auto coarseA = fineA->restrict(P);
        coarseA->describe(out, Teuchos::VERB_EXTREME);
        if (lvlNo+1 == auxH->GetNumLevels())
          lvl->Set("A", coarseA->toMatrix());
        else
          lvl->Set("A", rcp_dynamic_cast<Operator>(coarseA));
      }

      RCP<HierarchyManager> mueLuFactory = rcp(new ParameterListInterpreter(params,op->getDomainMap()->getComm()));
      H->setlib(op->getDomainMap()->lib());
      H->SetProcRankVerbose(op->getDomainMap()->getComm()->getRank());
      mueLuFactory->SetupHierarchy(*H);
      H->IsPreconditioner(true);
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
      RCP<Teuchos::ParameterList> belosList = Teuchos::parameterList();
      belosList->set("Maximum Iterations",    1000); // Maximum number of iterations allowed
      belosList->set("Convergence Tolerance", 1e-5);    // Relative convergence tolerance requested
      belosList->set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
      belosList->set("Output Frequency",      1);
      belosList->set("Output Style",          Belos::Brief);

      belosProblem->setRightPrec(belosPrec);

      bool set = belosProblem->setProblem();
      if (set == false) {
        throw MueLu::Exceptions::RuntimeError("ERROR:  Belos::LinearProblem failed to set up correctly!");
      }

      // Create an iterative solver manager
      Belos::SolverFactory<Scalar, MV, OP> solverFactory;
      RCP< Belos::SolverManager<Scalar, MV, OP> > solver = solverFactory.create(belosType, belosList);
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

    stacked_timer->stop("Hierarchical Driver");
    Teuchos::StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    if (printTimings)
      stacked_timer->report(out, comm, options);

#endif // HAVE_MUELU_BELOS
  }

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} //main

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc,argv);
}
