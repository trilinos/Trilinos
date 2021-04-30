/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

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
    using local_matrix_type = typename crs_matrix_type::local_matrix_type;
    using vector_type = Tpetra::Vector<Scalar, LO, GO, NT>;
    using multivector_type = Tpetra::MultiVector<Scalar, LO, GO, NT>;
    using block_multivector_type = Tpetra::BlockMultiVector<Scalar, LO, GO, NT>;
    using mem_space_t = typename local_matrix_type::memory_space;
    using rowmap_t = typename local_matrix_type::row_map_type::HostMirror;
    using entries_t = typename local_matrix_type::index_type::HostMirror;
    using values_t = typename local_matrix_type::values_type::HostMirror;
    using Offset = typename rowmap_t::non_const_value_type;
    using IST = typename crs_matrix_type::impl_scalar_type;
    using KAT = Kokkos::ArithTraits<IST>;
    //Type of view representing inverse diagonal blocks, and its HostMirror.
    using InverseBlocks = Kokkos::View<IST***, typename bcrs_matrix_type::device_type>;
    using InverseBlocksHost = typename InverseBlocks::HostMirror;

    //Setup for CrsMatrix
    GaussSeidel(const crs_matrix_type& A, Teuchos::RCP<vector_type>& inverseDiagVec_, Teuchos::ArrayRCP<LO>& applyRows_, Scalar omega_)
    {
      numRows = A.getNodeNumRows();
      inverseDiagVec = inverseDiagVec_;
      applyRows = applyRows_;
      blockSize = 1;
      omega = omega_;
      auto Alocal = A.getLocalMatrix();
      Arowmap = Kokkos::create_mirror_view(Alocal.graph.row_map);
      Aentries = Kokkos::create_mirror_view(Alocal.graph.entries);
      Avalues = Kokkos::create_mirror_view(Alocal.values);
      Kokkos::deep_copy(Arowmap, Alocal.graph.row_map);
      Kokkos::deep_copy(Aentries, Alocal.graph.entries);
      Kokkos::deep_copy(Avalues, Alocal.values);
    }

    GaussSeidel(const row_matrix_type& A, Teuchos::RCP<vector_type>& inverseDiagVec_, Teuchos::ArrayRCP<LO>& applyRows_, Scalar omega_)
    {
      numRows = A.getNodeNumRows();
      inverseDiagVec = inverseDiagVec_;
      applyRows = applyRows_;
      blockSize = 1;
      omega = omega_;
      //Here, need to make a deep CRS copy to avoid slow access via getLocalRowCopy
      Arowmap = rowmap_t(Kokkos::ViewAllocateWithoutInitializing("Arowmap"), numRows + 1);
      Aentries = entries_t(Kokkos::ViewAllocateWithoutInitializing("Aentries"), A.getNodeNumEntries());
      Avalues = values_t(Kokkos::ViewAllocateWithoutInitializing("Avalues"), A.getNodeNumEntries());
      size_t maxDegree = A.getNodeMaxNumRowEntries();
      Teuchos::Array<Scalar> rowValues(maxDegree);
      Teuchos::Array<LO> rowEntries(maxDegree);
      size_t accum = 0;
      for(LO i = 0; i <= numRows; i++)
      {
        Arowmap(i) = accum;
        if(i == numRows)
          break;
        size_t degree;
        A.getLocalRowCopy(i, rowEntries(), rowValues(), degree);
        accum += degree;
        size_t rowBegin = Arowmap(i);
        for(size_t j = 0; j < degree; j++)
        {
          Aentries(rowBegin + j) = rowEntries[j];
          Avalues(rowBegin + j) = rowValues[j];
        }
      }
    }

    GaussSeidel(const bcrs_matrix_type& A, const InverseBlocks& inverseBlockDiag_, Teuchos::ArrayRCP<LO>& applyRows_, Scalar omega_)
    {
      numRows = A.getNodeNumRows();
      //note: next 2 lines are no-ops if inverseBlockDiag_ is already host-accessible
      inverseBlockDiag = Kokkos::create_mirror_view(inverseBlockDiag_);
      Kokkos::deep_copy(inverseBlockDiag, inverseBlockDiag_);
      applyRows = applyRows_;
      omega = omega_;
      auto AlocalGraph = A.getCrsGraph().getLocalGraph();
      //A.sync_host();  //note: this only syncs values, not graph
      Avalues = A.getValuesHost();
      Arowmap = Kokkos::create_mirror_view(AlocalGraph.row_map);
      Aentries = Kokkos::create_mirror_view(AlocalGraph.entries);
      Kokkos::deep_copy(Arowmap, AlocalGraph.row_map);
      Kokkos::deep_copy(Aentries, AlocalGraph.entries);
      blockSize = A.getBlockSize();
    }

    template<bool useApplyRows, bool multipleRHS, bool omegaNotOne>
    void applyImpl(multivector_type& x, const multivector_type& b, Tpetra::ESweepDirection direction)
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
          auto bRow = b.getLocalBlock (row, v, Tpetra::Access::ReadOnly);
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
          IST* blk = &Avalues(j * bs2);
          for(LO v = 0; v < numVecs; v++)
          {
            auto xCol = x.getLocalBlock (col, v, Tpetra::Access::ReadOnly);
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
          auto bRow = b.getLocalBlock (row, v, Tpetra::Access::ReadOnly);
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
          auto xRow = x.getLocalBlock (row, v, Tpetra::Access::ReadWrite);
          for(LO k = 0; k < blockSize; k++)
          {
            xRow(k) += omega * dinv_accum(k, v);
          }
        }
      }
    }

    void apply(multivector_type& x, const multivector_type& b, Tpetra::ESweepDirection direction)
    {
      if(direction == Tpetra::Symmetric)
      {
        apply(x, b, Tpetra::Forward);
        apply(x, b, Tpetra::Backward);
        return;
      }
      else if(applyRows.is_null())
      {
        if(x.getNumVectors() > 1)
          this->template applyImpl<false, true, true>(x, b, direction);
        else
        {
          //Optimize for the all-rows, single vector, omega = 1 case
          if(omega == KAT::one())
            this->template applyImpl<false, false, false>(x, b, direction);
          else
            this->template applyImpl<false, false, true>(x, b, direction);
        }
      }
      else
      {
        this->template applyImpl<true, true, true>(x, b, direction);
      }
    }

    values_t Avalues; //length = Aentries.extent(0) * blockSize * blockSize
    rowmap_t Arowmap;
    entries_t Aentries;
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
