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
  {
    mv_map_ = this->getMap();
  }

  
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
    using Teuchos::rcpFromPtr;
    using Teuchos::as;
    
    size_t num_vecs = getGlobalNumVectors();

#ifdef HAVE_AMESOS2_DEBUG
    size_t requested_vector_length = distribution_map->getNodeNumElements();
    TEST_FOR_EXCEPTION( lda < requested_vector_length,
			std::invalid_argument,
			"Given stride is not large enough for local vector length" );
    TEST_FOR_EXCEPTION( as<size_t>(av.size()) < as<size_t>((num_vecs-1) * lda + requested_vector_length),
			std::invalid_argument,
			"MultiVector storage not large enough given leading dimension "
			"and number of vectors" );
#endif

    multivec_t redist_mv(rcpFromPtr(distribution_map), num_vecs);

    const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> importer(mv_map_, rcpFromPtr(distribution_map));
    redist_mv.doImport(*mv_, importer, Tpetra::REPLACE);

    // do copy
    redist_mv.get1dCopy(av, lda);
  }

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  Teuchos::ArrayRCP<Scalar>
  MultiVecAdapter<
    MultiVector<Scalar,
		LocalOrdinal,
		GlobalOrdinal,
		Node> >::get1dViewNonConst(bool local)
  {
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
    using Teuchos::rcpFromPtr;

    const size_t num_vecs  = getGlobalNumVectors();
    const multivec_t source_mv(rcpFromPtr(source_map), new_data, lda, num_vecs);
    const Tpetra::Import<local_ordinal_t, global_ordinal_t, node_t> importer(rcpFromPtr(source_map), mv_map_);

    mv_->doImport(source_mv, importer, Tpetra::REPLACE);
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


  /// Print a description of this adapter to the Fancy Output Stream.
  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  void
  MultiVecAdapter<
    MultiVector<Scalar,
		LocalOrdinal,
		GlobalOrdinal,
		Node> >::describe(Teuchos::FancyOStream& os,
				  const Teuchos::EVerbosityLevel verbLevel=
				  Teuchos::Describable::verbLevel_default) const
  {

  }


  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  const char* MultiVecAdapter<
    MultiVector<Scalar,
		LocalOrdinal,
		GlobalOrdinal,
		Node> >::name
  = "Amesos2 adapter for Tpetra::MultiVector";


} // end namespace Amesos2

#endif // AMESOS2_TPETRA_MULTIVEC_ADAPTER_DEF_HPP
