// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_GAUSS_SEIDEL_HPP
#define IFPACK2_GAUSS_SEIDEL_HPP

namespace Ifpack2
{
namespace Details
{
  template<typename Scalar, typename LO, typename GO, typename NT>
  struct GaussSeidel
  {
    using crs_matrix_type = Tpetra::CrsMatrix<Scalar, LO, GO, NT>;
    using bcrs_matrix_type = Tpetra::BlockCrsMatrix<Scalar, LO, GO, NT>;
    using row_matrix_type = Tpetra::RowMatrix<Scalar, LO, GO, NT>;
    using local_matrix_device_type = typename crs_matrix_type::local_matrix_device_type;
    using vector_type = Tpetra::Vector<Scalar, LO, GO, NT>;
    using multivector_type = Tpetra::MultiVector<Scalar, LO, GO, NT>;
    using block_multivector_type = Tpetra::BlockMultiVector<Scalar, LO, GO, NT>;
    using mem_space_t = typename local_matrix_device_type::memory_space;
    using rowmap_t = typename local_matrix_device_type::row_map_type::HostMirror;
    using entries_t = typename local_matrix_device_type::index_type::HostMirror;
    using values_t = typename local_matrix_device_type::values_type::HostMirror;
    using const_rowmap_t = typename rowmap_t::const_type;
    using const_entries_t = typename entries_t::const_type;
    using const_values_t = typename values_t::const_type;
    using Offset = typename rowmap_t::non_const_value_type;
    using IST = typename crs_matrix_type::impl_scalar_type;
    using KAT = Kokkos::ArithTraits<IST>;
    //Type of view representing inverse diagonal blocks, and its HostMirror.
    using InverseBlocks = Kokkos::View<IST***, typename bcrs_matrix_type::device_type>;
    using InverseBlocksHost = typename InverseBlocks::HostMirror;

    typedef typename crs_matrix_type::nonconst_global_inds_host_view_type nonconst_global_inds_host_view_type;
    typedef typename crs_matrix_type::nonconst_local_inds_host_view_type nonconst_local_inds_host_view_type;
    typedef typename crs_matrix_type::nonconst_values_host_view_type nonconst_values_host_view_type;

    //Setup for CrsMatrix
    GaussSeidel(Teuchos::RCP<const crs_matrix_type> A_, Teuchos::RCP<vector_type>& inverseDiagVec_, Teuchos::ArrayRCP<LO>& applyRows_, Scalar omega_)
    {
      A = A_;
      numRows = A_->getLocalNumRows();
      haveRowMatrix = false;
      inverseDiagVec = inverseDiagVec_;
      applyRows = applyRows_;
      blockSize = 1;
      omega = omega_;
    }

    GaussSeidel(Teuchos::RCP<const row_matrix_type> A_, Teuchos::RCP<vector_type>& inverseDiagVec_, Teuchos::ArrayRCP<LO>& applyRows_, Scalar omega_)
    {
      A = A_;
      numRows = A_->getLocalNumRows();
      haveRowMatrix = true;
      inverseDiagVec = inverseDiagVec_;
      applyRows = applyRows_;
      blockSize = 1;
      omega = omega_;
      //Here, need to make a deep CRS copy to avoid slow access via getLocalRowCopy
      rowMatrixRowmap = rowmap_t(Kokkos::ViewAllocateWithoutInitializing("Arowmap"), numRows + 1);
      rowMatrixEntries = entries_t(Kokkos::ViewAllocateWithoutInitializing("Aentries"), A_->getLocalNumEntries());
      rowMatrixValues = values_t(Kokkos::ViewAllocateWithoutInitializing("Avalues"), A_->getLocalNumEntries());
      size_t maxDegree = A_->getLocalMaxNumRowEntries();
      nonconst_values_host_view_type rowValues("rowValues", maxDegree);
      nonconst_local_inds_host_view_type rowEntries("rowEntries", maxDegree);
      size_t accum = 0;
      for(LO i = 0; i <= numRows; i++)
      {
        rowMatrixRowmap(i) = accum;
        if(i == numRows)
          break;
        size_t degree;
        A_->getLocalRowCopy(i, rowEntries, rowValues, degree);
        accum += degree;
        size_t rowBegin = rowMatrixRowmap(i);
        for(size_t j = 0; j < degree; j++)
        {
          rowMatrixEntries(rowBegin + j) = rowEntries[j];
          rowMatrixValues(rowBegin + j) = rowValues[j];
        }
      }
    }

    GaussSeidel(Teuchos::RCP<const bcrs_matrix_type> A_, const InverseBlocks& inverseBlockDiag_, Teuchos::ArrayRCP<LO>& applyRows_, Scalar omega_)
    {
      A = A_;
      numRows = A_->getLocalNumRows();
      haveRowMatrix = false;
      //note: next 2 lines are no-ops if inverseBlockDiag_ is already host-accessible
      inverseBlockDiag = Kokkos::create_mirror_view(inverseBlockDiag_);
      Kokkos::deep_copy(inverseBlockDiag, inverseBlockDiag_);
      applyRows = applyRows_;
      omega = omega_;
      blockSize = A_->getBlockSize();
    }

    template<bool useApplyRows, bool multipleRHS, bool omegaNotOne>
    void applyImpl(const const_values_t& Avalues, const const_rowmap_t& Arowmap, const const_entries_t& Aentries, multivector_type& x, const multivector_type& b, Tpetra::ESweepDirection direction)
    {
      //note: direction is either Forward or Backward (Symmetric is handled in apply())
      LO numApplyRows = useApplyRows ? (LO) applyRows.size() : numRows;
      //note: inverseDiagMV always has only one column
      auto inverseDiag = Kokkos::subview(inverseDiagVec->getLocalViewHost(Tpetra::Access::ReadOnly), Kokkos::ALL(), 0);
      bool forward = direction == Tpetra::Forward;
      if(multipleRHS)
      {
        LO numVecs = x.getNumVectors();
        Teuchos::Array<IST> accum(numVecs);
        auto xlcl = x.get2dViewNonConst();
        auto blcl = b.get2dView();
        for(LO i = 0; i < numApplyRows; i++)
        {
          LO row;
          if(forward)
            row = useApplyRows ? applyRows[i] : i;
          else
            row = useApplyRows ? applyRows[numApplyRows - 1 - i] : numApplyRows - 1 - i;
          for(LO k = 0; k < numVecs; k++)
          {
            accum[k] = KAT::zero();
          }
          Offset rowBegin = Arowmap(row);
          Offset rowEnd = Arowmap(row + 1);
          for(Offset j = rowBegin; j < rowEnd; j++)
          {
            LO col = Aentries(j);
            IST val = Avalues(j);
            for(LO k = 0; k < numVecs; k++)
            {
              accum[k] += val * IST(xlcl[k][col]);
            }
          }
          //Update x
          IST dinv = inverseDiag(row);
          for(LO k = 0; k < numVecs; k++)
          {
            if(omegaNotOne)
              xlcl[k][row] += Scalar(omega * dinv * (IST(blcl[k][row]) - accum[k]));
            else
              xlcl[k][row] += Scalar(dinv * (IST(blcl[k][row]) - accum[k]));
          }
        }
      }
      else
      {
        auto xlcl = Kokkos::subview(x.getLocalViewHost(Tpetra::Access::ReadWrite), Kokkos::ALL(), 0);
        auto blcl = Kokkos::subview(b.getLocalViewHost(Tpetra::Access::ReadOnly), Kokkos::ALL(), 0);
        auto dlcl = Kokkos::subview(inverseDiagVec->getLocalViewHost(Tpetra::Access::ReadOnly), Kokkos::ALL(), 0);
        for(LO i = 0; i < numApplyRows; i++)
        {
          LO row;
          if(forward)
            row = useApplyRows ? applyRows[i] : i;
          else
            row = useApplyRows ? applyRows[numApplyRows - 1 - i] : numApplyRows - 1 - i;
          IST accum = KAT::zero();
          Offset rowBegin = Arowmap(row);
          Offset rowEnd = Arowmap(row + 1);
          for(Offset j = rowBegin; j < rowEnd; j++)
          {
            accum += Avalues(j) * xlcl(Aentries(j));
          }
          //Update x
          IST dinv = dlcl(row);
          if(omegaNotOne)
            xlcl(row) += omega * dinv * (blcl(row) - accum);
          else
            xlcl(row) += dinv * (blcl(row) - accum);
        }
      }
    }

    void applyBlock(block_multivector_type& x, const block_multivector_type& b, Tpetra::ESweepDirection direction)
    {
      if(direction == Tpetra::Symmetric)
      {
        applyBlock(x, b, Tpetra::Forward);
        applyBlock(x, b, Tpetra::Backward);
        return;
      }
      auto Abcrs = Teuchos::rcp_dynamic_cast<const bcrs_matrix_type>(A);
      if(Abcrs.is_null())
        throw std::runtime_error("Ifpack2::Details::GaussSeidel::applyBlock: A must be a BlockCrsMatrix");
      auto Avalues = Abcrs->getValuesHost();
      auto AlclGraph = Abcrs->getCrsGraph().getLocalGraphHost();
      auto Arowmap = AlclGraph.row_map;
      auto Aentries = AlclGraph.entries;
      //Number of scalars in Avalues per block entry.
      Offset bs2 = blockSize * blockSize;
      LO numVecs = x.getNumVectors();
      Kokkos::View<IST**, Kokkos::LayoutLeft, Kokkos::HostSpace> accum(
          Kokkos::ViewAllocateWithoutInitializing("BlockGaussSeidel Accumulator"), blockSize, numVecs);
      Kokkos::View<IST**, Kokkos::LayoutLeft, Kokkos::HostSpace> dinv_accum(
          Kokkos::ViewAllocateWithoutInitializing("BlockGaussSeidel A_ii^-1*Accumulator"), blockSize, numVecs);
      bool useApplyRows = !applyRows.is_null();
      bool forward = direction == Tpetra::Forward;
      LO numApplyRows = useApplyRows ? applyRows.size() : numRows;
      for(LO i = 0; i < numApplyRows; i++)
      {
        LO row;
        if(forward)
          row = useApplyRows ? applyRows[i] : i;
        else
          row = useApplyRows ? applyRows[numApplyRows - 1 - i] : numApplyRows - 1 - i;
        for(LO v = 0; v < numVecs; v++)
        {
          auto bRow = b.getLocalBlockHost (row, v, Tpetra::Access::ReadOnly);
          for(LO k = 0; k < blockSize; k++)
          {
            accum(k, v) = KAT::zero();
          }
        }
        Offset rowBegin = Arowmap(row);
        Offset rowEnd = Arowmap(row + 1);
        for(Offset j = rowBegin; j < rowEnd; j++)
        {
          LO col = Aentries(j);
          const IST* blk = &Avalues(j * bs2);
          for(LO v = 0; v < numVecs; v++)
          {
            auto xCol = x.getLocalBlockHost (col, v, Tpetra::Access::ReadOnly);
            for(LO br = 0; br < blockSize; br++)
            {
              for(LO bc = 0; bc < blockSize; bc++)
              {
                IST Aval = blk[br * blockSize + bc];
                accum(br, v) += Aval * xCol(bc);
              }
            }
          }
        }
        //Update x: term is omega * Aii^-1 * accum, where omega is scalar, Aii^-1 is bs*bs, and accum is bs*nv
        auto invBlock = Kokkos::subview(inverseBlockDiag, row, Kokkos::ALL(), Kokkos::ALL());
        Kokkos::deep_copy(dinv_accum, KAT::zero());
        for(LO v = 0; v < numVecs; v++)
        {
          auto bRow = b.getLocalBlockHost (row, v, Tpetra::Access::ReadOnly);
          for(LO br = 0; br < blockSize; br++)
          {
            accum(br, v) = bRow(br) - accum(br, v);
          }
        }
        for(LO v = 0; v < numVecs; v++)
        {
          for(LO br = 0; br < blockSize; br++)
          {
            for(LO bc = 0; bc < blockSize; bc++)
            {
              dinv_accum(br, v) += invBlock(br, bc) * accum(bc, v);
            }
          }
        }
        //Update x
        for(LO v = 0; v < numVecs; v++)
        {
          auto xRow = x.getLocalBlockHost (row, v, Tpetra::Access::ReadWrite);
          for(LO k = 0; k < blockSize; k++)
          {
            xRow(k) += omega * dinv_accum(k, v);
          }
        }
      }
    }

    //Version of apply for CrsMatrix/RowMatrix (for BlockCrsMatrix, call applyBlock)
    void apply(multivector_type& x, const multivector_type& b, Tpetra::ESweepDirection direction)
    {
      if(direction == Tpetra::Symmetric)
      {
        apply(x, b, Tpetra::Forward);
        apply(x, b, Tpetra::Backward);
        return;
      }
      const_values_t Avalues;
      const_rowmap_t Arowmap;
      const_entries_t Aentries;
      if(haveRowMatrix)
      {
        Avalues = rowMatrixValues;
        Arowmap = rowMatrixRowmap;
        Aentries = rowMatrixEntries;
      }
      else
      {
        auto Acrs = Teuchos::rcp_dynamic_cast<const crs_matrix_type>(A);
        if(Acrs.is_null())
          throw std::runtime_error("Ifpack2::Details::GaussSeidel::apply: either haveRowMatrix, or A is CrsMatrix");
        auto Alcl = Acrs->getLocalMatrixHost();
        Avalues = Alcl.values;
        Arowmap = Alcl.graph.row_map;
        Aentries = Alcl.graph.entries;
      }
      if(applyRows.is_null())
      {
        if(x.getNumVectors() > 1)
          this->template applyImpl<false, true, true>(Avalues, Arowmap, Aentries, x, b, direction);
        else
        {
          //Optimize for the all-rows, single vector, omega = 1 case
          if(omega == KAT::one())
            this->template applyImpl<false, false, false>(Avalues, Arowmap, Aentries, x, b, direction);
          else
            this->template applyImpl<false, false, true>(Avalues, Arowmap, Aentries, x, b, direction);
        }
      }
      else
      {
        this->template applyImpl<true, true, true>(Avalues, Arowmap, Aentries, x, b, direction);
      }
    }

    Teuchos::RCP<const row_matrix_type> A;
    //For efficiency, if input is a RowMatrix, keep a persistent copy of the CRS formatted local matrix.
    bool haveRowMatrix;
    values_t rowMatrixValues;
    rowmap_t rowMatrixRowmap;
    entries_t rowMatrixEntries;
    LO numRows;
    IST omega;
    //If set up with BlockCrsMatrix, the block size. Otherwise 1.
    LO blockSize;
    //If using a non-block matrix, the inverse diagonal.
    Teuchos::RCP<vector_type> inverseDiagVec;
    //If using a block matrix, the inverses of all diagonal blocks.
    InverseBlocksHost inverseBlockDiag;
    //If null, apply over all rows in natural order. Otherwise, apply for each row listed in order.
    Teuchos::ArrayRCP<LO> applyRows;
  };
}
}

#endif
