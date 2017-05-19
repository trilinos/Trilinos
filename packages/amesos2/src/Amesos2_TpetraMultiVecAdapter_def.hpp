// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
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

#include "Amesos2_TpetraMultiVecAdapter_decl.hpp"


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

    typedef typename multivec_t::dual_view_type dual_view_type;
    typedef typename dual_view_type::host_mirror_space host_execution_space;
    mv_->template sync<host_execution_space> ();
    auto contig_local_view_2d = mv_->template getLocalView<host_execution_space>();
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
    const size_t requested_vector_length = distribution_map->getNodeNumElements ();
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
        // counted.)  Map's clone() method suffices, even though it only
        // makes a shallow copy of some of Map's data, because Map is
        // immutable and those data are reference-counted (e.g.,
        // ArrayRCP or RCP).
        distMap = distribution_map->template clone<Node> (distribution_map->getNode ());

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
        typedef typename multivec_t::dual_view_type dual_view_type;
        typedef typename dual_view_type::host_mirror_space host_execution_space;
        redist_mv.template sync < host_execution_space > ();

        auto contig_local_view_2d = redist_mv.template getLocalView<host_execution_space>();
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
            auto X_j = redist_mv.getVector(j);
            auto X_lcl_j_2d = redist_mv.template getLocalView<host_execution_space> ();
            auto X_lcl_j_1d = Kokkos::subview (X_lcl_j_2d, Kokkos::ALL (), j);
            for ( size_t i = 0; i < lclNumRows; ++i ) {
              av_j[i] = X_lcl_j_1d(i);
            }
          }
        }

        auto global_contiguous_size = distMap->getGlobalNumElements(); //maybe use getGlobalLength() from the mv
        auto local_contiguous_size = (distMap->getComm()->getRank() == 0) ? global_contiguous_size : 0;
        RCP<const map_type> contigMap = rcp( new map_type(global_contiguous_size, local_contiguous_size, 0, distMap->getComm() ));

        typedef Tpetra::Export<local_ordinal_t, global_ordinal_t, node_t> contiguous_export_t;
        RCP<contiguous_export_t> contig_exporter = rcp( new contiguous_export_t(redist_mv.getMap(), contigMap) );
        multivec_t contig_mv( contigMap, num_vecs);
        contig_mv.doExport(redist_mv, *contig_exporter, Tpetra::INSERT);
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
    //    * mv_->getMap()->getNodeElementList() to get a subCopy of mv_ which we
    //    * assign to l_l_mv_, then return get1dViewNonConst() of l_l_mv_
    //    */
    //   if(l_l_mv_.is_null() ){
    //  Teuchos::ArrayView<const GlobalOrdinal> nodeElements_go
    //    = mv_->getMap()->getNodeElementList();
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
    //    = mv_->getMap()->getNodeElementList();
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
      typedef typename multivec_t::dual_view_type::host_mirror_space host_execution_space;
      // num_vecs = 1; stride does not matter
      auto mv_view_to_modify_2d = mv_->template getLocalView<host_execution_space>();
      for ( size_t i = 0; i < lda; ++i ) {
        mv_view_to_modify_2d(i,0) = new_data[i];
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
        // Map's clone() method suffices, even though it only makes a
        // shallow copy of some of Map's data, because Map is immutable
        // and those data are reference-counted (e.g., ArrayRCP or RCP).
        srcMap = source_map->template clone<Node> (source_map->getNode ());
        importer_ = rcp (new import_type (srcMap, this->getMap ()));
      }
      else {
        srcMap = importer_->getSourceMap ();
      }

      multivec_t redist_mv (srcMap, num_vecs);

      if ( distribution != CONTIGUOUS_AND_ROOTED ) {
        // Do this if GIDs contiguous - existing functionality
        // Redistribute the output (multi)vector.
        const multivec_t source_mv (srcMap, new_data, lda, num_vecs);
        mv_->doImport (source_mv, *importer_, Tpetra::REPLACE);
      }
      else {
        typedef typename multivec_t::dual_view_type dual_view_type;
        typedef typename dual_view_type::host_mirror_space host_execution_space;
        redist_mv.template modify< host_execution_space > ();

        if ( redist_mv.isConstantStride() ) {
          auto contig_local_view_2d = redist_mv.template getLocalView<host_execution_space>();
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
            auto X_j = redist_mv.getVector(j);
            auto X_lcl_j_2d = redist_mv.template getLocalView<host_execution_space> ();
            auto X_lcl_j_1d = Kokkos::subview (X_lcl_j_2d, Kokkos::ALL (), j);
            for ( size_t i = 0; i < lclNumRows; ++i ) {
              X_lcl_j_1d(i) = av_j[i];
            }
          }
        }

        typedef typename multivec_t::node_type::memory_space memory_space;
        redist_mv.template sync <memory_space> ();

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
