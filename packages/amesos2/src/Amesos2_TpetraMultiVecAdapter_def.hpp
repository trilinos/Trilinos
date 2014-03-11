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
                                                       Node> > distribution_map ) const
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

    // Copy the imported (multi)vector's data into the ArrayView.
    redist_mv.get1dCopy (av, lda);
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
                                                       Node> > source_map)
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

    const multivec_t source_mv (srcMap, new_data, lda, num_vecs);

    // Redistribute the output (multi)vector.
    mv_->doImport (source_mv, *importer_, Tpetra::REPLACE);
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
