// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2010 Sandia Corporation
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

namespace Amesos {

  using Tpetra::MultiVector;

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  void
  MultiVecAdapter<
    MultiVector<Scalar,
		LocalOrdinal,
		GlobalOrdinal,
		Node> >::get1dCopy(const Teuchos::ArrayView<Scalar>& av,
				   size_t lda,
				   bool global_copy) const
  {
    // size_t local_length = getLocalLength();
    global_size_type global_length = getGlobalLength();
    size_t num_vecs = getGlobalNumVectors();

    Teuchos::RCP<const Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type > > o_map, l_map;
    o_map = getMap();

    if( global_copy ){
      // l_map = Tpetra::createLocalMapWithNode
      //        <local_ordinal_type,
      //        global_ordinal_type,
      //        node_type>( Teuchos::as<size_t>(global_length),
      //                    getComm(),
      //                    o_map->getNode() );

      size_t num_my_elements;
      if( getComm()->getRank() == 0 ){
	num_my_elements = Teuchos::as<size_t>(global_length);
      } else {
	num_my_elements = Teuchos::OrdinalTraits<size_t>::zero();
      }

      l_map = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal>( global_length,
									num_my_elements,
									o_map->getIndexBase(),
									o_map->getComm(),
									o_map->getNode()));
    } else {
      l_map = o_map;
    }

    Teuchos::RCP<multivec_type> l_mv;
    l_mv = Teuchos::rcp(new multivec_type(l_map, num_vecs));

    Teuchos::RCP<Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> > importer;
    importer = Teuchos::rcp(new Tpetra::Import<LocalOrdinal,GlobalOrdinal>(o_map, l_map));

    l_mv->doImport(*mv_, *importer, Tpetra::REPLACE);

    // Finally, do copy
    l_mv->get1dCopy(av, lda);
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


  // Implementation is almost identical to get1dCopy()
  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  void
  MultiVecAdapter<
    MultiVector<Scalar,
		LocalOrdinal,
		GlobalOrdinal,
		Node> >::get2dCopy(Teuchos::ArrayView<const Teuchos::ArrayView<scalar_type> > A) const
  {
    // TEST_FOR_EXCEPTION(Teuchos::as<size_t>(A.size()) != this->getGlobalNumVectors(),
    //                 std::length_error,
    //                 "get2dCopy() : Length of A must be equal to the number of vectors");

    // if( mv_->isDistributed() ){
    //   this->localize();

    //   l_mv_->get2dCopy(A);
    // } else {
    //   mv_->get2dCopy(A);
    // }
  }


  // Implementation is almost identical to get1dViewNonConst()
  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >
  MultiVecAdapter<
    MultiVector<Scalar,
		LocalOrdinal,
		GlobalOrdinal,
		Node> >::get2dViewNonConst( bool local )
  {
    // if( local ){
    //   this->localize();
    //   /* Use the global element list returned by
    //    * mv_->getMap()->getNodeElementList() to get a subCopy of mv_ which we
    //    * assign to l_mv_, then return get2dViewNonConst() of l_mv_
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

    //  return(l_l_mv_->get2dViewNonConst());
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

    //  return l_l_mv_->get2dViewNonConst();
    //   }
    // } else {
    //   if( mv_->isDistributed() ){
    //  this->localize();

    //  return l_mv_->get2dViewNonConst();
    //   } else {                      // not distributed, no need to import
    //  return mv_->get2dViewNonConst();
    //   }
    // }
  }


  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  void
  MultiVecAdapter<
    MultiVector<Scalar,
		LocalOrdinal,
		GlobalOrdinal,
		Node> >::globalize(int root)
  {
    if( mv_->isDistributed() ){
      mv_->doImport(*l_mv_, *exporter_, Tpetra::REPLACE);
    }
  }


  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node>
  template <typename Value_t>
  void
  MultiVecAdapter<
    MultiVector<Scalar,
		LocalOrdinal,
		GlobalOrdinal,
		Node> >::globalize(const Teuchos::ArrayView<Value_t>& newVals,
				   int root)
  {
#ifdef HAVE_AMESOS2_DEBUG
    std::cout << "globalizing : " << newVals.toString() << std::endl;
#endif
    size_t num_vecs  = getGlobalNumVectors();
    global_size_type global_num_elems = getGlobalLength();
    size_t num_elems;
    const Teuchos::RCP<const Teuchos::Comm<int> > comm = getComm();

    if( comm->getRank() == root ){ // root distributes its values
      num_elems = Teuchos::as<size_t>(global_num_elems);
    } else {
      num_elems = Teuchos::OrdinalTraits<global_size_type>::zero();
    }
    Teuchos::RCP<const Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type > > o_map, l_map;
    o_map = getMap();
    l_map = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal>(global_num_elems,
								     num_elems,
								     o_map->getIndexBase(),
								     getComm(),
								     o_map->getNode()));

    // Do conversion to Scalar values
    Teuchos::Array<Scalar> converted_vals(newVals.size());
    typedef typename Teuchos::Array<Scalar>::iterator scalar_it_t;
    typedef typename Teuchos::ArrayView<Value_t>::iterator value_it_t;
    value_it_t val_it = newVals.begin();
    scalar_it_t it = converted_vals.begin();
    scalar_it_t end = converted_vals.end();
    for( ; it != end; ++it )
      {
	*it = Teuchos::as<Scalar>(*val_it++);
      }

    Teuchos::RCP<multivec_type> l_mv;
    l_mv = Teuchos::rcp(new multivec_type(l_map,
					  converted_vals(),
					  Teuchos::as<size_t>(global_num_elems),
					  num_vecs));

    Teuchos::RCP<Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> > importer;
    importer = Teuchos::rcp(new Tpetra::Import<LocalOrdinal,GlobalOrdinal>(l_map, o_map));

    mv_->doImport(*l_mv, *importer, Tpetra::REPLACE);
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
		Node> >::describe(
				  Teuchos::FancyOStream& os,
				  const Teuchos::EVerbosityLevel verbLevel=
				  Teuchos::Describable::verbLevel_default) const
  {

  }

  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  void
  MultiVecAdapter<
    MultiVector<Scalar,
		LocalOrdinal,
		GlobalOrdinal,
		Node> >::localize(bool root) const
  {
    if( l_mv_.is_null() ){
      // create local multivector, maps, and do import
      o_map_ = mv_->getMap();
      const Teuchos::RCP<const Teuchos::Comm<int> > comm = o_map_->getComm();

      size_t num_vecs  = this->getGlobalNumVectors();
      global_size_type global_num_elems = this->getGlobalLength();
      size_t num_elems;
      if( root ){ // root gets all the vectors
	num_elems = Teuchos::as<size_t>(global_num_elems);
      } else {
	num_elems = Teuchos::OrdinalTraits<global_size_type>::zero();
      }
      // l_map_ = Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal>(
      //   numVectors,
      //   comm,
      //   o_map_->getNode());
      l_map_ = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal>(global_num_elems,
									num_elems,
									o_map_->getIndexBase(),
									comm,
									o_map_->getNode()));

      l_mv_ = Teuchos::rcp(new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>(l_map_,num_vecs));

      importer_ = Teuchos::rcp(new Tpetra::Import<LocalOrdinal,GlobalOrdinal>(o_map_, l_map_));
      exporter_ = Teuchos::rcp(new Tpetra::Export<LocalOrdinal,GlobalOrdinal>(o_map_, l_map_));

      l_mv_->doImport(*mv_, *importer_, Tpetra::REPLACE);
    } else {
      // Just update local values
      l_mv_->doImport(*mv_, *importer_, Tpetra::REPLACE);
    }
  }


  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
  const char* MultiVecAdapter<
    MultiVector<Scalar,
		LocalOrdinal,
		GlobalOrdinal,
		Node> >::name
  = "Amesos2 adapter for Tpetra::MultiVector";


} // end namespace Amesos

#endif // AMESOS2_TPETRA_MULTIVEC_ADAPTER_DEF_HPP
