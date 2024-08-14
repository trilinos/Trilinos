// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_BLOCKTRIDICONTAINER_DEF_HPP
#define IFPACK2_BLOCKTRIDICONTAINER_DEF_HPP

#include <Teuchos_Details_MpiTypeTraits.hpp>

#include <Tpetra_Distributor.hpp>
#include <Tpetra_BlockMultiVector.hpp>
#include <Tpetra_BlockCrsMatrix_Helpers.hpp>

#include <Kokkos_ArithTraits.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Vector.hpp>
#include <KokkosBatched_AddRadial_Decl.hpp>
#include <KokkosBatched_AddRadial_Impl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Gemm_Serial_Impl.hpp>
#include <KokkosBatched_Gemv_Decl.hpp>
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
                  const Teuchos::RCP<const import_type>& importer,
                  const bool overlapCommAndComp,
                  const bool useSeqMethod,
                  const int block_size,
                  const bool explicitConversion) 
  {
    IFPACK2_BLOCKHELPER_TIMER_WITH_FENCE("BlockTriDiContainer::initInternal", initInternal, typename BlockHelperDetails::ImplType<MatrixType>::execution_space);

    // create pointer of impl
    {
      IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::createImpl", createImpl);
      impl_ = Teuchos::rcp(new BlockTriDiContainerDetails::ImplObject<MatrixType>());
      IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
    }

    using impl_type = BlockHelperDetails::ImplType<MatrixType>;
    // using block_crs_matrix_type = typename impl_type::tpetra_block_crs_matrix_type;

    {
      IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::setA", setA);
      if (explicitConversion) {
        impl_->A = Teuchos::rcp_dynamic_cast<const block_crs_matrix_type>(matrix);
        if (impl_->A.is_null()) {
          TEUCHOS_TEST_FOR_EXCEPT_MSG
            (block_size == -1, "A pointwise matrix and block_size = -1 were given as inputs.");
          {
            IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::setA::convertToBlockCrsMatrix", convertToBlockCrsMatrix);
            impl_->A = Tpetra::convertToBlockCrsMatrix(*Teuchos::rcp_dynamic_cast<const crs_matrix_type>(matrix), block_size, true);
            IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
          }
        }
      }
      else {
        impl_->A = matrix;
      }
      IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
    }

    impl_->tpetra_importer = Teuchos::null;
    impl_->async_importer  = Teuchos::null;
    
    if (useSeqMethod)
    {
      IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::createBlockCrsTpetraImporter useSeqMethod", useSeqMethod);
      if (importer.is_null()) // there is no given importer, then create one
        impl_->tpetra_importer = BlockTriDiContainerDetails::createBlockCrsTpetraImporter<MatrixType>(impl_->A);
      else
        impl_->tpetra_importer = importer; // if there is a given importer, use it
      IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
    }
    else
    {
      IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::createBlockCrsTpetraImporter", createBlockCrsTpetraImporter);
      //Leave tpetra_importer null even if user provided an importer.
      //It is not used in the performant codepath (!useSeqMethod)
      impl_->async_importer = BlockTriDiContainerDetails::createBlockCrsAsyncImporter<MatrixType>(impl_->A);
      IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
    }

    // as a result, there are 
    // 1) tpetra_importer is     null , async_importer is     null (no need for importer)
    // 2) tpetra_importer is NOT null , async_importer is     null (sequential method is used)
    // 3) tpetra_importer is     null , async_importer is NOT null (async method is used)

    // temporary disabling 
    impl_->overlap_communication_and_computation = overlapCommAndComp;

    {
      IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::createZ", createZ);
      impl_->Z = typename impl_type::tpetra_multivector_type();
      IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
    }
    {
      IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::createW", createW);
      impl_->W = typename impl_type::impl_scalar_type_1d_view();
      IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
    }

    IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
  }

  template <typename MatrixType>
  void
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::clearInternal ()
  {
    IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::clearInternal", clearInternal);
    using impl_type = BlockHelperDetails::ImplType<MatrixType>;
    using part_interface_type = BlockHelperDetails::PartInterface<MatrixType>;
    using block_tridiags_type = BlockTriDiContainerDetails::BlockTridiags<MatrixType>;
    using amd_type = BlockHelperDetails::AmD<MatrixType>;
    using norm_manager_type = BlockHelperDetails::NormManager<MatrixType>;
    
    impl_->A = Teuchos::null;
    impl_->tpetra_importer = Teuchos::null;
    impl_->async_importer  = Teuchos::null;

    impl_->Z = typename impl_type::tpetra_multivector_type();
    impl_->W = typename impl_type::impl_scalar_type_1d_view();

    impl_->part_interface  = part_interface_type();
    impl_->block_tridiags  = block_tridiags_type();
    impl_->a_minus_d       = amd_type();
    impl_->work            = typename impl_type::vector_type_1d_view();
    impl_->norm_manager    = norm_manager_type();

    impl_ = Teuchos::null;
    IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
  }

  template <typename MatrixType>
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::BlockTriDiContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                       const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                       const Teuchos::RCP<const import_type>& importer,
                       bool pointIndexed)
    : Container<MatrixType>(matrix, partitions, pointIndexed), partitions_(partitions)
  {
    IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::BlockTriDiContainer", BlockTriDiContainer);
    const bool useSeqMethod = false;
    const bool overlapCommAndComp = false;
    initInternal(matrix, importer, overlapCommAndComp, useSeqMethod);
    n_subparts_per_part_ = -1;
    block_size_ = -1;
    IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
  }

  template <typename MatrixType>
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::BlockTriDiContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                       const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                       const int n_subparts_per_part,
                       const bool overlapCommAndComp, 
                       const bool useSeqMethod,
                       const int block_size,
                       const bool explicitConversion)
    : Container<MatrixType>(matrix, partitions, false), partitions_(partitions)
  {
    IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::BlockTriDiContainer", BlockTriDiContainer);
    initInternal(matrix, Teuchos::null, overlapCommAndComp, useSeqMethod, block_size, explicitConversion);
    n_subparts_per_part_ = n_subparts_per_part;
    block_size_ = block_size;
    IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
  }

  template <typename MatrixType>
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::~BlockTriDiContainer ()
  {
  }

  template <typename MatrixType>
  void 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::setParameters (const Teuchos::ParameterList& List)
  {
    if (List.isType<int>("partitioner: subparts per part"))
      n_subparts_per_part_ = List.get<int>("partitioner: subparts per part");
    if (List.isType<int>("partitioner: block size"))
      block_size_ = List.get<int>("partitioner: block size");
  }

  template <typename MatrixType>
  void 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::initialize ()
  {
    IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::initialize", initialize);
    this->IsInitialized_ = true;
    {
      auto bA = Teuchos::rcp_dynamic_cast<const block_crs_matrix_type>(impl_->A);
      if (bA.is_null()) {
        TEUCHOS_TEST_FOR_EXCEPT_MSG
          (block_size_ == -1, "A pointwise matrix and block_size = -1 were given as inputs.");
        {
          IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::initialize::getBlockCrsGraph", getBlockCrsGraph);
          auto A = Teuchos::rcp_dynamic_cast<const crs_matrix_type>(impl_->A);
          impl_->blockGraph = Tpetra::getBlockCrsGraph(*A, block_size_, true);
          IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
        }
      }
      else {
        impl_->blockGraph = Teuchos::rcpFromRef(bA->getCrsGraph());
      }
    }

    {
      IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::createPartInterfaceBlockTridiagsNormManager", createPartInterfaceBlockTridiagsNormManager);
      impl_->part_interface  = BlockTriDiContainerDetails::createPartInterface<MatrixType>(impl_->A, impl_->blockGraph, partitions_, n_subparts_per_part_);
      impl_->block_tridiags  = BlockTriDiContainerDetails::createBlockTridiags<MatrixType>(impl_->part_interface);
      impl_->norm_manager    = BlockHelperDetails::NormManager<MatrixType>(impl_->A->getComm());
      IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
    }

    // We assume that if you called this method, you intend to recompute
    // everything.
    this->IsComputed_ = false;
    TEUCHOS_ASSERT(!impl_->A.is_null()); // when initInternal is called, A_ must be set
    {
      BlockTriDiContainerDetails::performSymbolicPhase<MatrixType>
        (impl_->A, 
         impl_->blockGraph, 
         impl_->part_interface, 
         impl_->block_tridiags, 
         impl_->a_minus_d, 
         impl_->overlap_communication_and_computation);    
    }
    IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
  }

  template <typename MatrixType>
  void 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::compute ()
  {
    IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::compute", compute);
    this->IsComputed_ = false;
    if (!this->isInitialized())
      this->initialize();
    {
      BlockTriDiContainerDetails::performNumericPhase<MatrixType>
        (impl_->A,
         impl_->blockGraph, 
         impl_->part_interface, impl_->block_tridiags, 
         Kokkos::ArithTraits<magnitude_type>::zero());
    }
    this->IsComputed_ = true;
    IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
  }

  template <typename MatrixType>
  void 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::clearBlocks ()
  {
    IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::clearBlocks", clearBlocks);
    clearInternal();
    this->IsInitialized_ = false;
    this->IsComputed_ = false;
    Container<MatrixType>::clearBlocks();
    IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
  }

  template <typename MatrixType>
  void 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::applyInverseJacobi (const mv_type& X, mv_type& Y, scalar_type dampingFactor,
                        bool zeroStartingSolution, int numSweeps) const
  {
    IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::applyInverseJacobi", applyInverseJacobi);
    const magnitude_type tol = Kokkos::ArithTraits<magnitude_type>::zero();
    const int check_tol_every = 1;

    BlockTriDiContainerDetails::applyInverseJacobi<MatrixType>
      (impl_->A,
       impl_->blockGraph,
       impl_->tpetra_importer, 
       impl_->async_importer, 
       impl_->overlap_communication_and_computation,
       X, Y, impl_->Z, impl_->W,
       impl_->part_interface, impl_->block_tridiags, impl_->a_minus_d,
       impl_->work,
       impl_->norm_manager,
       dampingFactor,
       zeroStartingSolution,
       numSweeps,
       tol,
       check_tol_every);
    IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
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
    IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::compute", compute);
    this->IsComputed_ = false;
    if (!this->isInitialized())
      this->initialize();
    {
      BlockTriDiContainerDetails::performNumericPhase<MatrixType>
        (impl_->A,
         impl_->blockGraph, 
         impl_->part_interface, impl_->block_tridiags, 
         in.addRadiallyToDiagonal);
    }
    this->IsComputed_ = true;
    IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
  }

  template <typename MatrixType>
  typename BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>::ApplyParameters
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::createDefaultApplyParameters () const
  {
    ApplyParameters in;
    in.dampingFactor = Teuchos::ScalarTraits<scalar_type>::one();
    return in;
  }

  template <typename MatrixType>
  int 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::applyInverseJacobi (const mv_type& X, mv_type& Y, 
                        const ApplyParameters& in) const
  {
    IFPACK2_BLOCKHELPER_TIMER("BlockTriDiContainer::applyInverseJacobi", applyInverseJacobi);
    int r_val = 0;
    {
      r_val = BlockTriDiContainerDetails::applyInverseJacobi<MatrixType>
        (impl_->A,
         impl_->blockGraph,
         impl_->tpetra_importer, 
         impl_->async_importer,
         impl_->overlap_communication_and_computation,
         X, Y, impl_->Z, impl_->W,
         impl_->part_interface, impl_->block_tridiags, impl_->a_minus_d,
         impl_->work,
         impl_->norm_manager,
         in.dampingFactor,
         in.zeroStartingSolution,
         in.maxNumSweeps,
         in.tolerance,
         in.checkToleranceEvery);
    }
    IFPACK2_BLOCKHELPER_TIMER_FENCE(typename BlockHelperDetails::ImplType<MatrixType>::execution_space)
    return r_val;
  }

  template <typename MatrixType>
  const typename BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>::magnitude_type
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::getNorms0 () const {
    return impl_->norm_manager.getNorms0();
  }

  template <typename MatrixType>
  const typename BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>::magnitude_type
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::getNormsFinal () const {
    return impl_->norm_manager.getNormsFinal();
  }

  template <typename MatrixType>
  void 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::apply (ConstHostView /* X */, HostView /* Y */, int /* blockIndex */, Teuchos::ETransp /* mode */,
           scalar_type /* alpha */, scalar_type /* beta */) const
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "BlockTriDiContainer::apply is not implemented. You may have reached this message "
                                << "because you want to use this container's performance-portable Jacobi iteration. In "
                                << "that case, set \"relaxation: type\" to \"MT Split Jacobi\" rather than \"Jacobi\".");
  }

  template <typename MatrixType>
  void 
  BlockTriDiContainer<MatrixType, BlockTriDiContainerDetails::ImplSimdTag>
  ::weightedApply (ConstHostView /* X */, HostView /* Y */, ConstHostView /* D */, int /* blockIndex */,
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
    if (this->isInitialized()) {
      if (this->isComputed()) {
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
       << "isInitialized()         = " << this->IsInitialized_ << endl
       << "isComputed()            = " << this->IsComputed_ << endl
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
