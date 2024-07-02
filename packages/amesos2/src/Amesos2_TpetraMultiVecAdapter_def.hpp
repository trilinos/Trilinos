// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
  \file   Amesos2_TpetraMultiVecAdapter_def.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Wed May 26 19:48:32 CDT 2010

  \brief  Amesos2::MultiVecAdapter specialization for the
          Tpetra::MultiVector class.
*/

#ifndef AMESOS2_TPETRA_MULTIVEC_ADAPTER_DEF_HPP
#define AMESOS2_TPETRA_MULTIVEC_ADAPTER_DEF_HPP

#include <type_traits>
#include "Amesos2_TpetraMultiVecAdapter_decl.hpp"
#include "Amesos2_Kokkos_View_Copy_Assign.hpp"


namespace Amesos2 {

  using Tpetra::MultiVector;

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  MultiVecAdapter<
    MultiVector<Scalar,
                LocalOrdinal,
                GlobalOrdinal,
                Node> >::MultiVecAdapter( const Teuchos::RCP<multivec_t>& m )
  : mv_(m)
  {}

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  Teuchos::RCP<
    MultiVector<Scalar,
                LocalOrdinal,
                GlobalOrdinal,
                Node> >
  MultiVecAdapter<
    MultiVector<Scalar,
                LocalOrdinal,
                GlobalOrdinal,
                Node> >::clone() const
  {
      using MV = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
      Teuchos::RCP<MV> Y (new MV (mv_->getMap(), mv_->getNumVectors(), false));
      Y->setCopyOrView (Teuchos::View);
      return Y;
  }

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  typename MultiVecAdapter<
    MultiVector<Scalar,
                LocalOrdinal,
                GlobalOrdinal,
                Node> >::multivec_t::impl_scalar_type *
  MultiVecAdapter<
    MultiVector<Scalar,
                LocalOrdinal,
                GlobalOrdinal,
                Node> >::getMVPointer_impl() const
  {
  TEUCHOS_TEST_FOR_EXCEPTION( this->getGlobalNumVectors() != 1,
		      std::invalid_argument,
		      "Amesos2_TpetraMultiVectorAdapter: getMVPointer_impl should only be called for case with a single vector and single MPI process" );

    auto contig_local_view_2d = mv_->getLocalViewHost(Tpetra::Access::ReadWrite);
    auto contig_local_view_1d = Kokkos::subview (contig_local_view_2d, Kokkos::ALL (), 0);
    return contig_local_view_1d.data();
  }

  // TODO Proper type handling: 
  // Consider a MultiVectorTraits class
  // typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> multivector_type
  // NOTE: In this class, above already has a typedef multivec_t
  // typedef typename multivector_type::impl_scalar_type return_scalar_type; // this is the POD type the dual_view_type is templated on
  // Traits class needed to do this generically for the general MultiVectorAdapter interface


  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  void
  MultiVecAdapter<
    MultiVector<Scalar,
                LocalOrdinal,
                GlobalOrdinal,
                Node> >::get1dCopy(const Teuchos::ArrayView<scalar_t>& av,
                                   size_t lda,
                                   Teuchos::Ptr<
                                     const Tpetra::Map<LocalOrdinal,
                                                       GlobalOrdinal,
                                                       Node> > distribution_map,
                                                       EDistribution distribution) const
  {
    using Teuchos::as;
    using Teuchos::RCP;
    typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
    const size_t num_vecs = getGlobalNumVectors ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      distribution_map.getRawPtr () == NULL, std::invalid_argument,
      "Amesos2::MultiVecAdapter::get1dCopy: distribution_map argument is null.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      mv_.is_null (), std::logic_error,
      "Amesos2::MultiVecAdapter::get1dCopy: mv_ is null.");
    // Check mv_ before getMap(), because the latter calls mv_->getMap().
    TEUCHOS_TEST_FOR_EXCEPTION(
      this->getMap ().is_null (), std::logic_error,
      "Amesos2::MultiVecAdapter::get1dCopy: this->getMap() returns null.");

#ifdef HAVE_AMESOS2_DEBUG
    const size_t requested_vector_length = distribution_map->getLocalNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      lda < requested_vector_length, std::invalid_argument,
      "Amesos2::MultiVecAdapter::get1dCopy: On process " <<
      distribution_map->getComm ()->getRank () << " of the distribution Map's "
      "communicator, the given stride lda = " << lda << " is not large enough "
      "for the local vector length " << requested_vector_length << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      as<size_t> (av.size ()) < as<size_t> ((num_vecs - 1) * lda + requested_vector_length),
      std::invalid_argument, "Amesos2::MultiVector::get1dCopy: MultiVector "
      "storage not large enough given leading dimension and number of vectors." );
#endif // HAVE_AMESOS2_DEBUG

    // Special case when number vectors == 1 and single MPI process
    if ( num_vecs == 1 && this->getComm()->getRank() == 0 && this->getComm()->getSize() == 1 ) {
      mv_->get1dCopy (av, lda);
    }
    else {

      // (Re)compute the Export object if necessary.  If not, then we
      // don't need to clone distribution_map; we can instead just get
      // the previously cloned target Map from the Export object.
      RCP<const map_type> distMap;
      if (exporter_.is_null () ||
          ! exporter_->getSourceMap ()->isSameAs (* (this->getMap ())) ||
          ! exporter_->getTargetMap ()->isSameAs (* distribution_map)) {
        // Since we're caching the Export object, and since the Export
        // needs to keep the distribution Map, we have to make a copy of
        // the latter in order to ensure that it will stick around past
        // the scope of this function call.  (Ptr is not reference
        // counted.)  
        distMap = rcp(new map_type(*distribution_map));
        // (Re)create the Export object.
        exporter_ = rcp (new export_type (this->getMap (), distMap));
      }
      else {
        distMap = exporter_->getTargetMap ();
      }

      multivec_t redist_mv (distMap, num_vecs);

      // Redistribute the input (multi)vector.
      redist_mv.doExport (*mv_, *exporter_, Tpetra::REPLACE);

      if ( distribution != CONTIGUOUS_AND_ROOTED ) {
        // Do this if GIDs contiguous - existing functionality
        // Copy the imported (multi)vector's data into the ArrayView.
        redist_mv.get1dCopy (av, lda);
      }
      else {
        // Do this if GIDs not contiguous...
        // sync is needed for example if mv was updated on device, but will be passed through Amesos2 to solver running on host
        //redist_mv.template sync < host_execution_space > ();
        auto contig_local_view_2d = redist_mv.getLocalViewHost(Tpetra::Access::ReadOnly);
        if ( redist_mv.isConstantStride() ) {
          for ( size_t j = 0; j < num_vecs; ++j) {
            auto av_j = av(lda*j, lda);
            for ( size_t i = 0; i < lda; ++i ) {
              av_j[i] = contig_local_view_2d(i,j); //lda may not be correct if redist_mv is not constant stride...
            }
          }
        }
        else {
          // ... lda should come from Teuchos::Array* allocation,
          // not the MultiVector, since the MultiVector does NOT
          // have constant stride in this case.
          // TODO lda comes from X->getGlobalLength() in solve_impl - should this be changed???
          const size_t lclNumRows = redist_mv.getLocalLength();
          for (size_t j = 0; j < redist_mv.getNumVectors(); ++j) {
            auto av_j = av(lda*j, lclNumRows);
            auto X_lcl_j_2d = redist_mv.getLocalViewHost(Tpetra::Access::ReadOnly);
            auto X_lcl_j_1d = Kokkos::subview (X_lcl_j_2d, Kokkos::ALL (), j);

            using val_type = typename std::remove_const<typename decltype( X_lcl_j_1d )::value_type>::type;
            Kokkos::View<val_type*, Kokkos::HostSpace> umavj ( const_cast< val_type* > ( reinterpret_cast<const val_type*> ( av_j.getRawPtr () ) ), av_j.size () );
            Kokkos::deep_copy (umavj, X_lcl_j_1d);
          }
        }
      }
    }
  }

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  template <typename KV>
  bool
  MultiVecAdapter<
    MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>
  >::get1dCopy_kokkos_view(
    bool bInitialize,
    KV& kokkos_view,
    [[maybe_unused]] size_t lda,
    Teuchos::Ptr<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > distribution_map,
    EDistribution distribution
  ) const {
    using Teuchos::as;
    using Teuchos::RCP;
    typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
    const size_t num_vecs = getGlobalNumVectors ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      distribution_map.getRawPtr () == NULL, std::invalid_argument,
      "Amesos2::MultiVecAdapter::get1dCopy_kokkos_view: distribution_map argument is null.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      mv_.is_null (), std::logic_error,
      "Amesos2::MultiVecAdapter::get1dCopy_kokkos_view: mv_ is null.");
    // Check mv_ before getMap(), because the latter calls mv_->getMap().
    TEUCHOS_TEST_FOR_EXCEPTION(
      this->getMap ().is_null (), std::logic_error,
      "Amesos2::MultiVecAdapter::get1dCopy_kokkos_view: this->getMap() returns null.");

#ifdef HAVE_AMESOS2_DEBUG
    const size_t requested_vector_length = distribution_map->getLocalNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      lda < requested_vector_length, std::invalid_argument,
      "Amesos2::MultiVecAdapter::get1dCopy_kokkos_view: On process " <<
      distribution_map->getComm ()->getRank () << " of the distribution Map's "
      "communicator, the given stride lda = " << lda << " is not large enough "
      "for the local vector length " << requested_vector_length << ".");

    // Note do not check size since deep_copy_or_assign_view below will allocate
    // if necessary - but may just assign ptrs.
#endif // HAVE_AMESOS2_DEBUG

    // Special case when number vectors == 1 and single MPI process
    if ( num_vecs == 1 && this->getComm()->getRank() == 0 && this->getComm()->getSize() == 1 ) {
      if(mv_->isConstantStride()) {
        bool bAssigned;
        //deep_copy_or_assign_view(bInitialize, kokkos_view, mv_->getLocalViewDevice(Tpetra::Access::ReadOnly), bAssigned);
        deep_copy_only(bInitialize, kokkos_view, mv_->getLocalViewDevice(Tpetra::Access::ReadOnly), bAssigned);
        return bAssigned; // if bAssigned is true we are accessing the mv data directly without a copy
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Resolve handling for non-constant stride.");
      }
    }
    else {

      // (Re)compute the Export object if necessary.  If not, then we
      // don't need to clone distribution_map; we can instead just get
      // the previously cloned target Map from the Export object.
      RCP<const map_type> distMap;
      if (exporter_.is_null () ||
          ! exporter_->getSourceMap ()->isSameAs (* (this->getMap ())) ||
          ! exporter_->getTargetMap ()->isSameAs (* distribution_map)) {
        // Since we're caching the Export object, and since the Export
        // needs to keep the distribution Map, we have to make a copy of
        // the latter in order to ensure that it will stick around past
        // the scope of this function call.  (Ptr is not reference
        // counted.)
        distMap = rcp(new map_type(*distribution_map));
        // (Re)create the Export object.
        exporter_ = rcp (new export_type (this->getMap (), distMap));
      }
      else {
        distMap = exporter_->getTargetMap ();
      }

      multivec_t redist_mv (distMap, num_vecs);

      // Redistribute the input (multi)vector.
      redist_mv.doExport (*mv_, *exporter_, Tpetra::REPLACE);

      if ( distribution != CONTIGUOUS_AND_ROOTED ) {
        // Do this if GIDs contiguous - existing functionality
        // Copy the imported (multi)vector's data into the Kokkos View.
        bool bAssigned;
        deep_copy_or_assign_view(bInitialize, kokkos_view, redist_mv.getLocalViewDevice(Tpetra::Access::ReadOnly), bAssigned);
        return false; // do not return bAssigned because redist_mv was already copied so always return false
      }
      else {
        if(redist_mv.isConstantStride()) {
          bool bAssigned; // deep_copy_or_assign_view sets true if assigned (no deep copy)
          deep_copy_or_assign_view(bInitialize, kokkos_view, redist_mv.getLocalViewDevice(Tpetra::Access::ReadOnly), bAssigned);
          return false; // do not return bAssigned because redist_mv was already copied so always return false
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Kokkos adapter non-constant stride not imlemented.");
        }
      }
    }
  }

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  Teuchos::ArrayRCP<Scalar>
  MultiVecAdapter<
    MultiVector<Scalar,
                LocalOrdinal,
                GlobalOrdinal,
                Node> >::get1dViewNonConst (bool local)
  {
    // FIXME (mfh 22 Jan 2014) When I first found this routine, all of
    // its code was commented out, and it didn't return anything.  The
    // latter is ESPECIALLY dangerous, given that the return value is
    // an ArrayRCP.  Thus, I added the exception throw below.
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Amesos2::MultiVecAdapter::get1dViewNonConst: "
      "Not implemented.");

    // if( local ){
    //   this->localize();
    //   /* Use the global element list returned by
    //    * mv_->getMap()->getLocalElementList() to get a subCopy of mv_ which we
    //    * assign to l_l_mv_, then return get1dViewNonConst() of l_l_mv_
    //    */
    //   if(l_l_mv_.is_null() ){
    //  Teuchos::ArrayView<const GlobalOrdinal> nodeElements_go
    //    = mv_->getMap()->getLocalElementList();
    //  Teuchos::Array<size_t> nodeElements_st(nodeElements_go.size());

    //  // Convert the node element to a list of size_t type objects
    //  typename Teuchos::ArrayView<const GlobalOrdinal>::iterator it_go;
    //  Teuchos::Array<size_t>::iterator it_st = nodeElements_st.begin();
    //  for( it_go = nodeElements_go.begin(); it_go != nodeElements_go.end(); ++it_go ){
    //    *(it_st++) = Teuchos::as<size_t>(*it_go);
    //  }

    //  // To be consistent with the globalize() function, get a view of the local mv
    //  l_l_mv_ = l_mv_->subViewNonConst(nodeElements_st);

    //  return(l_l_mv_->get1dViewNonConst());
    //   } else {
    //  // We need to re-import values to the local, since changes may have been
    //  // made to the global structure that are not reflected in the local
    //  // view.
    //  Teuchos::ArrayView<const GlobalOrdinal> nodeElements_go
    //    = mv_->getMap()->getLocalElementList();
    //  Teuchos::Array<size_t> nodeElements_st(nodeElements_go.size());

    //  // Convert the node element to a list of size_t type objects
    //  typename Teuchos::ArrayView<const GlobalOrdinal>::iterator it_go;
    //  Teuchos::Array<size_t>::iterator it_st = nodeElements_st.begin();
    //  for( it_go = nodeElements_go.begin(); it_go != nodeElements_go.end(); ++it_go ){
    //    *(it_st++) = Teuchos::as<size_t>(*it_go);
    //  }

    //  return l_l_mv_->get1dViewNonConst();
    //   }
    // } else {
    //   if( mv_->isDistributed() ){
    //  this->localize();

    //  return l_mv_->get1dViewNonConst();
    //   } else {                      // not distributed, no need to import
    //  return mv_->get1dViewNonConst();
    //   }
    // }
  }


  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node>
  void
  MultiVecAdapter<
    MultiVector<Scalar,
                LocalOrdinal,
                GlobalOrdinal,
                Node> >::put1dData(const Teuchos::ArrayView<const scalar_t>& new_data,
                                   size_t lda,
                                   Teuchos::Ptr<
                                     const Tpetra::Map<LocalOrdinal,
                                                       GlobalOrdinal,
                                                       Node> > source_map,
                                                       EDistribution distribution )
  {
    using Teuchos::RCP;
    typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

    TEUCHOS_TEST_FOR_EXCEPTION(
      source_map.getRawPtr () == NULL, std::invalid_argument,
      "Amesos2::MultiVecAdapter::put1dData: source_map argument is null.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      mv_.is_null (), std::logic_error,
      "Amesos2::MultiVecAdapter::put1dData: the internal MultiVector mv_ is null.");
    // getMap() calls mv_->getMap(), so test first whether mv_ is null.
    TEUCHOS_TEST_FOR_EXCEPTION(
      this->getMap ().is_null (), std::logic_error,
      "Amesos2::MultiVecAdapter::put1dData: this->getMap() returns null.");

    const size_t num_vecs = getGlobalNumVectors ();

    // Special case when number vectors == 1 and single MPI process
    if ( num_vecs == 1 && this->getComm()->getRank() == 0 && this->getComm()->getSize() == 1 ) {
      // num_vecs = 1; stride does not matter
      auto mv_view_to_modify_2d = mv_->getLocalViewHost(Tpetra::Access::OverwriteAll);
      for ( size_t i = 0; i < lda; ++i ) {
        mv_view_to_modify_2d(i,0) = new_data[i]; // Only one vector
      }
    }
    else {

      // (Re)compute the Import object if necessary.  If not, then we
      // don't need to clone source_map; we can instead just get the
      // previously cloned source Map from the Import object.
      RCP<const map_type> srcMap;
      if (importer_.is_null () ||
          ! importer_->getSourceMap ()->isSameAs (* source_map) ||
          ! importer_->getTargetMap ()->isSameAs (* (this->getMap ()))) {

        // Since we're caching the Import object, and since the Import
        // needs to keep the source Map, we have to make a copy of the
        // latter in order to ensure that it will stick around past the
        // scope of this function call.  (Ptr is not reference counted.)
        srcMap = rcp(new map_type(*source_map));
        importer_ = rcp (new import_type (srcMap, this->getMap ()));
      }
      else {
        srcMap = importer_->getSourceMap ();
      }

      if ( distribution != CONTIGUOUS_AND_ROOTED ) {
        // Do this if GIDs contiguous - existing functionality
        // Redistribute the output (multi)vector.
        const multivec_t source_mv (srcMap, new_data, lda, num_vecs);
        mv_->doImport (source_mv, *importer_, Tpetra::REPLACE);
      }
      else {
        multivec_t redist_mv (srcMap, num_vecs); // unused for ROOTED case
        if ( redist_mv.isConstantStride() ) {
          auto contig_local_view_2d = redist_mv.getLocalViewHost(Tpetra::Access::OverwriteAll);
          for ( size_t j = 0; j < num_vecs; ++j) {
            auto av_j = new_data(lda*j, lda);
            for ( size_t i = 0; i < lda; ++i ) {
              contig_local_view_2d(i,j) = av_j[i];
            }
          }
        }
        else {
          // ... lda should come from Teuchos::Array* allocation,
          // not the MultiVector, since the MultiVector does NOT
          // have constant stride in this case.
          // TODO lda comes from X->getGlobalLength() in solve_impl - should this be changed???
          const size_t lclNumRows = redist_mv.getLocalLength();
          for (size_t j = 0; j < redist_mv.getNumVectors(); ++j) {
            auto av_j = new_data(lda*j, lclNumRows);
            auto X_lcl_j_2d = redist_mv.getLocalViewHost(Tpetra::Access::ReadOnly);
            auto X_lcl_j_1d = Kokkos::subview (X_lcl_j_2d, Kokkos::ALL (), j);

            using val_type = typename std::remove_const<typename decltype( X_lcl_j_1d )::value_type>::type;
            Kokkos::View<val_type*, Kokkos::HostSpace> umavj ( const_cast< val_type* > ( reinterpret_cast<const val_type*> ( av_j.getRawPtr () ) ), av_j.size () );
            Kokkos::deep_copy (umavj, X_lcl_j_1d);
          }
        }

        //typedef typename multivec_t::node_type::memory_space memory_space;
        //redist_mv.template sync <memory_space> ();

        mv_->doImport (redist_mv, *importer_, Tpetra::REPLACE);
      }
    }

  }

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node>
  template <typename KV>
  void
  MultiVecAdapter<
    MultiVector<Scalar,
                LocalOrdinal,
                GlobalOrdinal,
                Node> >::put1dData_kokkos_view(KV& kokkos_new_data,
                                   size_t lda,
                                   Teuchos::Ptr<
                                     const Tpetra::Map<LocalOrdinal,
                                                       GlobalOrdinal,
                                                       Node> > source_map,
                                                       EDistribution distribution )
  {
    using Teuchos::RCP;
    typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

    TEUCHOS_TEST_FOR_EXCEPTION(
      source_map.getRawPtr () == NULL, std::invalid_argument,
      "Amesos2::MultiVecAdapter::put1dData_kokkos_view: source_map argument is null.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      mv_.is_null (), std::logic_error,
      "Amesos2::MultiVecAdapter::put1dData_kokkos_view: the internal MultiVector mv_ is null.");
    // getMap() calls mv_->getMap(), so test first whether mv_ is null.
    TEUCHOS_TEST_FOR_EXCEPTION(
      this->getMap ().is_null (), std::logic_error,
      "Amesos2::MultiVecAdapter::put1dData_kokkos_view: this->getMap() returns null.");

    const size_t num_vecs = getGlobalNumVectors ();

    // Special case when number vectors == 1 and single MPI process
    if ( num_vecs == 1 && this->getComm()->getRank() == 0 && this->getComm()->getSize() == 1 ) {

      // num_vecs = 1; stride does not matter

      // If this is the optimized path then kokkos_new_data will be the dst
      auto mv_view_to_modify_2d = mv_->getLocalViewDevice(Tpetra::Access::OverwriteAll);
      //deep_copy_or_assign_view(mv_view_to_modify_2d, kokkos_new_data);
      deep_copy_only(mv_view_to_modify_2d, kokkos_new_data);
    }
    else {

      // (Re)compute the Import object if necessary.  If not, then we
      // don't need to clone source_map; we can instead just get the
      // previously cloned source Map from the Import object.
      RCP<const map_type> srcMap;
      if (importer_.is_null () ||
          ! importer_->getSourceMap ()->isSameAs (* source_map) ||
          ! importer_->getTargetMap ()->isSameAs (* (this->getMap ()))) {

        // Since we're caching the Import object, and since the Import
        // needs to keep the source Map, we have to make a copy of the
        // latter in order to ensure that it will stick around past the
        // scope of this function call.  (Ptr is not reference counted.)
        srcMap = rcp(new map_type(*source_map));
        importer_ = rcp (new import_type (srcMap, this->getMap ()));
      }
      else {
        srcMap = importer_->getSourceMap ();
      }

      if ( distribution != CONTIGUOUS_AND_ROOTED ) {
        // Use View scalar type, not MV Scalar because we want Kokkos::complex, not
        // std::complex to avoid a Kokkos::complex<double> to std::complex<float>
        // conversion which would require a double copy and fail here. Then we'll be
        // setup to safely reinterpret_cast complex to std if necessary.
        typedef typename multivec_t::dual_view_type::t_host::value_type tpetra_mv_view_type;
        Kokkos::View<tpetra_mv_view_type**,typename KV::array_layout,
          Kokkos::HostSpace> convert_kokkos_new_data;
        deep_copy_or_assign_view(convert_kokkos_new_data, kokkos_new_data);
#ifdef HAVE_TEUCHOS_COMPLEX
        // convert_kokkos_new_data may be Kokkos::complex and Scalar could be std::complex
        auto pData = reinterpret_cast<Scalar*>(convert_kokkos_new_data.data());
#else
        auto pData = convert_kokkos_new_data.data();
#endif

        const multivec_t source_mv (srcMap, Teuchos::ArrayView<const scalar_t>(
          pData, kokkos_new_data.size()), lda, num_vecs);
        mv_->doImport (source_mv, *importer_, Tpetra::REPLACE);
      }
      else {
        multivec_t redist_mv (srcMap, num_vecs); // unused for ROOTED case
        // Cuda solvers won't currently use this path since they are just serial
        // right now, so this mirror should be harmless (and not strictly necessary).
        // Adding it for future possibilities though we may then refactor this
        // for better efficiency if the kokkos_new_data view is on device.
        auto host_kokkos_new_data = Kokkos::create_mirror_view(kokkos_new_data);
        Kokkos::deep_copy(host_kokkos_new_data, kokkos_new_data);
        if ( redist_mv.isConstantStride() ) {
          auto contig_local_view_2d = redist_mv.getLocalViewHost(Tpetra::Access::OverwriteAll);
          for ( size_t j = 0; j < num_vecs; ++j) {
            auto av_j = Kokkos::subview(host_kokkos_new_data, Kokkos::ALL, j);
            for ( size_t i = 0; i < lda; ++i ) {
              contig_local_view_2d(i,j) = av_j(i);
            }
          }
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Kokkos adapter "
            "CONTIGUOUS_AND_ROOTED not implemented for put1dData_kokkos_view "
            "with non constant stride.");
        }

        //typedef typename multivec_t::node_type::memory_space memory_space;
        //redist_mv.template sync <memory_space> ();

        mv_->doImport (redist_mv, *importer_, Tpetra::REPLACE);
      }
    }

  }

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  std::string
  MultiVecAdapter<
    MultiVector<Scalar,
                LocalOrdinal,
                GlobalOrdinal,
                Node> >::description() const
  {
    std::ostringstream oss;
    oss << "Amesos2 adapter wrapping: ";
    oss << mv_->description();
    return oss.str();
  }


  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  void
  MultiVecAdapter<
    MultiVector<Scalar,
                LocalOrdinal,
                GlobalOrdinal,
                Node> >::describe (Teuchos::FancyOStream& os,
                                   const Teuchos::EVerbosityLevel verbLevel) const
  {
    mv_->describe (os, verbLevel);
  }


  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  const char* MultiVecAdapter<
    MultiVector<Scalar,
                LocalOrdinal,
                GlobalOrdinal,
                Node> >::name = "Amesos2 adapter for Tpetra::MultiVector";

} // end namespace Amesos2

#endif // AMESOS2_TPETRA_MULTIVEC_ADAPTER_DEF_HPP
