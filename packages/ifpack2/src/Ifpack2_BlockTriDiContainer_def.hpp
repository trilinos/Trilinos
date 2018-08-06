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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_BLOCKTRIDICONTAINER_DEF_HPP
#define IFPACK2_BLOCKTRIDICONTAINER_DEF_HPP

#include <Teuchos_Details_MpiTypeTraits.hpp>

#include <Tpetra_Distributor.hpp>
#include <Tpetra_BlockMultiVector.hpp>

#include <Kokkos_ArithTraits.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Vector.hpp>
#include <KokkosBatched_AddRadial_Decl.hpp>
#include <KokkosBatched_AddRadial_Impl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Gemm_Serial_Impl.hpp>
#include <KokkosBatched_Gemv_Decl.hpp>
#include <KokkosBatched_Gemv_Serial_Impl.hpp>
#include <KokkosBatched_Trsm_Decl.hpp>
#include <KokkosBatched_Trsm_Serial_Impl.hpp>
#include <KokkosBatched_Trsv_Decl.hpp>
#include <KokkosBatched_Trsv_Serial_Impl.hpp>
#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_LU_Serial_Impl.hpp>

#include "Ifpack2_BlockTriDiContainer_decl.hpp"
#include "Ifpack2_BlockTriDiContainer_impl.hpp"

#include <memory>


namespace Ifpack2 {
  
  ///
  /// BlockTriDiContainer, ImplSimdTag
  ///

  template <typename MatrixType>
  void
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::initInternal (const Teuchos::RCP<const row_matrix_type>& matrix,
                  const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                  const Teuchos::RCP<const import_type>& importer,
                  const bool overlapCommAndComp,
                  const bool useSeqMethod) 
  {
    // create pointer of impl
    impl_ = Teuchos::rcp(new BlockTriDiContainerDetails::ImplObject<MatrixType>());

    using impl_type = BlockTriDiContainerDetails::ImplType<MatrixType>;
    using block_crs_matrix_type = typename impl_type::tpetra_block_crs_matrix_type;

    impl_->A = Teuchos::rcp_dynamic_cast<const block_crs_matrix_type>(matrix);
    TEUCHOS_TEST_FOR_EXCEPT_MSG
      (impl_->A.is_null(), "BlockTriDiContainer currently supports Tpetra::BlockCrsMatrix only.");

    impl_->tpetra_importer = Teuchos::null;
    impl_->async_importer  = Teuchos::null;
    
    if (importer.is_null()) // there is no given importer, then create one
      if (useSeqMethod) 
        impl_->tpetra_importer = BlockTriDiContainerDetails::createBlockCrsTpetraImporter<MatrixType>(impl_->A);
      else 
        impl_->async_importer = BlockTriDiContainerDetails::createBlockCrsAsyncImporter<MatrixType>(impl_->A);
    else 
      impl_->tpetra_importer = importer; // if there is a given importer, use it

    // as a result, there are 
    // 1) tpetra_importer is     null , async_importer is     null (no need for importer)
    // 2) tpetra_importer is NOT null , async_importer is     null (sequential method is used)
    // 3) tpetra_importer is     null , async_importer is NOT null (async method is used)

    // temporary disabling 
    impl_->overlap_communication_and_computation = overlapCommAndComp;
    
    impl_->Z = typename impl_type::tpetra_multivector_type();

    impl_->part_interface  = BlockTriDiContainerDetails::createPartInterface<MatrixType>(impl_->A, partitions);
    impl_->block_tridiags  = BlockTriDiContainerDetails::createBlockTridiags<MatrixType>(impl_->part_interface);
    impl_->norm_manager    = BlockTriDiContainerDetails::NormManager<MatrixType>(impl_->A->getComm());
  }

  template <typename MatrixType>
  void
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::clearInternal ()
  {
    using impl_type = BlockTriDiContainerDetails::ImplType<MatrixType>;
    using part_interface_type = BlockTriDiContainerDetails::PartInterface<MatrixType>;
    using block_tridiags_type = BlockTriDiContainerDetails::BlockTridiags<MatrixType>;
    using amd_type = BlockTriDiContainerDetails::AmD<MatrixType>;
    using norm_manager_type = BlockTriDiContainerDetails::NormManager<MatrixType>;
    
    impl_->A = Teuchos::null;
    impl_->tpetra_importer = Teuchos::null;
    impl_->async_importer  = Teuchos::null;

    impl_->Z = typename impl_type::tpetra_multivector_type();

    impl_->part_interface  = part_interface_type();
    impl_->block_tridiags  = block_tridiags_type();
    impl_->a_minus_d       = amd_type();
    impl_->work            = typename impl_type::vector_type_1d_view();
    impl_->norm_manager    = norm_manager_type();

    impl_ = Teuchos::null;
  }

  template <typename MatrixType>
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::BlockTriDiContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                       const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                       const Teuchos::RCP<const import_type>& importer,
                       int OverlapLevel, scalar_type DampingFactor)
    : Container<MatrixType>(matrix, partitions, importer, OverlapLevel, DampingFactor)
  {
    const bool useSeqMethod = false;
    const bool overlapCommAndComp = false;
    initInternal(matrix, partitions, importer, overlapCommAndComp, useSeqMethod);
  }

  template <typename MatrixType>
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::BlockTriDiContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                       const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                       const bool overlapCommAndComp, const bool useSeqMethod)
    : Container<MatrixType>(matrix, partitions, Teuchos::null, 0,
                            Kokkos::ArithTraits<magnitude_type>::one())
  {
    initInternal(matrix, partitions, Teuchos::null, overlapCommAndComp, useSeqMethod);
  }

  template <typename MatrixType>
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::~BlockTriDiContainer ()
  {
  }

  template <typename MatrixType>
  bool 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::isInitialized () const
  {
    return IsInitialized_;
  }

  template <typename MatrixType>
  bool 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::isComputed () const
  {
    return IsComputed_;
  }

  template <typename MatrixType>
  void 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::setParameters (const Teuchos::ParameterList& /* List */)
  {
    // the solver doesn't currently take any parameters
  }

  template <typename MatrixType>
  void 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::initialize ()
  {
    IsInitialized_ = true;
    // We assume that if you called this method, you intend to recompute
    // everything.
    IsComputed_ = false;
    TEUCHOS_ASSERT(!impl_->A.is_null()); // when initInternal is called, A_ must be set
    {
      BlockTriDiContainerDetails::performSymbolicPhase<MatrixType>
        (impl_->A, 
         impl_->part_interface, impl_->block_tridiags, 
         impl_->a_minus_d, 
         impl_->overlap_communication_and_computation);    
    }
  }

  template <typename MatrixType>
  void 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::compute ()
  {
    IsComputed_ = false;
    if (!this->isInitialized())
      this->initialize();
    {
      BlockTriDiContainerDetails::performNumericPhase<MatrixType>
        (impl_->A, 
         impl_->part_interface, impl_->block_tridiags, 
         Kokkos::ArithTraits<magnitude_type>::zero());
    }
    IsComputed_ = true;
  }

  template <typename MatrixType>
  void 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::clearBlocks ()
  {
    clearInternal();
    IsInitialized_ = false;
    IsComputed_ = false;
  }

  template <typename MatrixType>
  void 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::applyInverseJacobi (const mv_type& X, mv_type& Y, bool zeroStartingSolution,
                        int numSweeps) const
  {
    const magnitude_type tol = Kokkos::ArithTraits<magnitude_type>::zero();
    const int check_tol_every = 1;

    BlockTriDiContainerDetails::applyInverseJacobi<MatrixType>
      (impl_->A,
       impl_->tpetra_importer, 
       impl_->async_importer, 
       impl_->overlap_communication_and_computation,
       X, Y, impl_->Z,
       impl_->part_interface, impl_->block_tridiags, impl_->a_minus_d,
       impl_->work,
       impl_->norm_manager,
       this->DampingFactor_,
       zeroStartingSolution,
       numSweeps,
       tol,
       check_tol_every);
  }

  template <typename MatrixType>
  typename BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>::ComputeParameters
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::createDefaultComputeParameters () const
  {
    return ComputeParameters();
  }

  template <typename MatrixType>
  void 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::compute (const ComputeParameters& in)
  {
    IsComputed_ = false;
    if (!this->isInitialized())
      this->initialize();
    {
      BlockTriDiContainerDetails::performNumericPhase<MatrixType>
        (impl_->A, 
         impl_->part_interface, impl_->block_tridiags, 
         in.addRadiallyToDiagonal);
    }
    IsComputed_ = true;
  }

  template <typename MatrixType>
  typename BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>::ApplyParameters
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::createDefaultApplyParameters () const
  {
    ApplyParameters in;
    in.dampingFactor = this->DampingFactor_;
    return in;
  }

  template <typename MatrixType>
  int 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::applyInverseJacobi (const mv_type& X, mv_type& Y, 
                        const ApplyParameters& in) const
  {
    int r_val = 0;
    {
      r_val = BlockTriDiContainerDetails::applyInverseJacobi<MatrixType>
        (impl_->A,
         impl_->tpetra_importer, 
         impl_->async_importer,
         impl_->overlap_communication_and_computation,
         X, Y, impl_->Z,
         impl_->part_interface, impl_->block_tridiags, impl_->a_minus_d,
         impl_->work,
         impl_->norm_manager,
         in.dampingFactor,
         in.zeroStartingSolution,
         in.maxNumSweeps,
         in.tolerance,
         in.checkToleranceEvery);
    }
    return r_val;
  }

  template <typename MatrixType>
  const Teuchos::ArrayRCP<const typename BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>::magnitude_type>
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::getNorms0 () const {
    const auto p = impl_->norm_manager.getNorms0();
    return Teuchos::arcp(p, 0, p ? impl_->norm_manager.getNumVectors() : 0, false);
  }

  template <typename MatrixType>
  const Teuchos::ArrayRCP<const typename BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>::magnitude_type>
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::getNormsFinal () const {
    const auto p = impl_->norm_manager.getNormsFinal();
    return Teuchos::arcp(p, 0, p ? impl_->norm_manager.getNumVectors() : 0, false);
  }

  template <typename MatrixType>
  void 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::apply (HostView& /* X */, HostView& /* Y */, int /* blockIndex */, int /* stride */, Teuchos::ETransp /* mode */,
           scalar_type /* alpha */, scalar_type /* beta */) const
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "BlockTriDiContainer::apply is not implemented. You may have reached this message "
                                << "because you want to use this container's performance-portable Jacobi iteration. In "
                                << "that case, set \"relaxation: type\" to \"MT Split Jacobi\" rather than \"Jacobi\".");
  }

  template <typename MatrixType>
  void 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::weightedApply (HostView& /* X */, HostView& /* Y */, HostView& /* D */, int /* blockIndex */, int /* stride */,
                   Teuchos::ETransp /* mode */, scalar_type /* alpha */, scalar_type /* beta */) const
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "BlockTriDiContainer::weightedApply is not implemented.");
  }

  template <typename MatrixType>
  std::ostream& 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::print (std::ostream& os) const
  {
    Teuchos::FancyOStream fos(Teuchos::rcp(&os,false));
    fos.setOutputToRootOnly(0);
    describe(fos);
    return os;
  }

  template <typename MatrixType>
  std::string 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::description () const
  {
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    if (isInitialized()) {
      if (isComputed()) {
        oss << "{status = initialized, computed";
      }
      else {
        oss << "{status = initialized, not computed";
      }
    }
    else {
      oss << "{status = not initialized, not computed";
    }

    oss << "}";
    return oss.str();
  }

  template <typename MatrixType>
  void
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>::
  describe (Teuchos::FancyOStream& os,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using std::endl;
    if(verbLevel==Teuchos::VERB_NONE) return;
    os << "================================================================================" << endl
       << "Ifpack2::BlockTriDiContainer" << endl
       << "Number of blocks        = " << this->numBlocks_ << endl
       << "isInitialized()         = " << IsInitialized_ << endl
       << "isComputed()            = " << IsComputed_ << endl
       << "================================================================================" << endl
       << endl;
  }

  template <typename MatrixType>
  std::string 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::getName() { return "Ifpack2::BlockTriDiContainer::ImplSimdTag"; }

#define IFPACK2_BLOCKTRIDICONTAINER_INSTANT(S,LO,GO,N)                  \
  template class Ifpack2::BlockTriDiContainer< Tpetra::RowMatrix<S, LO, GO, N> >;
}
#endif
