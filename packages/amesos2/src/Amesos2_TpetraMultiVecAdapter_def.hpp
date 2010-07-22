/**
  \file   Amesos2_TpetraMultiVecAdapter_def.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Wed May 26 19:48:32 CDT 2010

  \brief  Amesos2::MultiVecAdapter specialization for the
          Tpetra::MultiVector class.
*/

#ifndef AMESOS2_TPETRA_MULTIVEC_ADAPTER_DEF_HPP
#define AMESOS2_TPETRA_MULTIVEC_ADAPTER_DEF_HPP

using Tpetra::MultiVector;

namespace Amesos {


template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
void
MultiVecAdapter<
  MultiVector<Scalar,
              LocalOrdinal,
              GlobalOrdinal,
              Node> >::get1dCopy(
                const Teuchos::ArrayView<Scalar>& A,
                size_t lda) const
{
  if( mv_->isDistributed() ){
    this->localize();

    l_mv_->get1dCopy(A,lda);
  } else {
    mv_->get1dCopy(A,lda);
  }
}


template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
Teuchos::ArrayRCP<Scalar>
MultiVecAdapter<
  MultiVector<Scalar,
              LocalOrdinal,
              GlobalOrdinal,
              Node> >::get1dViewNonConst(bool local)
{
  if( local ){
    localize();
    /* Use the global element list returned by
     * mv_->getMap()->getNodeElementList() to get a subCopy of mv_ which we
     * assign to l_l_mv_, then return get1dViewNonConst() of l_l_mv_
     */
    if(l_l_mv_.is_null() ){
      Teuchos::ArrayView<const GlobalOrdinal> nodeElements_go
        = mv_->getMap()->getNodeElementList();
      Teuchos::Array<size_t> nodeElements_st(nodeElements_go.size());

      // Convert the node element to a list of size_t type objects
      typename Teuchos::ArrayView<const GlobalOrdinal>::iterator it_go;
      Teuchos::Array<size_t>::iterator it_st = nodeElements_st.begin();
      for( it_go = nodeElements_go.begin(); it_go != nodeElements_go.end(); ++it_go ){
        *(it_st++) = Teuchos::as<size_t>(*it_go);
      }

      // To be consistent with the globalize() function, get a view of the local mv
      l_l_mv_ = l_mv_->subViewNonConst(nodeElements_st);

      return(l_l_mv_->get1dViewNonConst());
    } else {
      // We need to re-import values to the local, since changes may have been
      // made to the global structure that are not reflected in the local
      // view.
      Teuchos::ArrayView<const GlobalOrdinal> nodeElements_go
        = mv_->getMap()->getNodeElementList();
      Teuchos::Array<size_t> nodeElements_st(nodeElements_go.size());

      // Convert the node element to a list of size_t type objects
      typename Teuchos::ArrayView<const GlobalOrdinal>::iterator it_go;
      Teuchos::Array<size_t>::iterator it_st = nodeElements_st.begin();
      for( it_go = nodeElements_go.begin(); it_go != nodeElements_go.end(); ++it_go ){
        *(it_st++) = Teuchos::as<size_t>(*it_go);
      }

      return l_l_mv_->get1dViewNonConst();
    }
  } else {
    if( mv_->isDistributed() ){
      this->localize();

      return l_mv_->get1dViewNonConst();
    } else {                      // not distributed, no need to import
      return mv_->get1dViewNonConst();
    }
  }
}


// Implementation is almost identical to get1dCopy()
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
void
MultiVecAdapter<
  MultiVector<Scalar,
              LocalOrdinal,
              GlobalOrdinal,
              Node> >::get2dCopy(
                Teuchos::ArrayView<const Teuchos::ArrayView<scalar_type> > A) const
{
  TEST_FOR_EXCEPTION(
    Teuchos::as<size_t>(A.size()) != this->getGlobalNumVectors(),
    std::length_error,
    "get2dCopy() : Length of A must be equal to the number of vectors");
  
  if( mv_->isDistributed() ){
    this->localize();

    l_mv_->get2dCopy(A);
  } else {
    mv_->get2dCopy(A);
  }
}


// Implementation is almost identical to get1dViewNonConst()
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >
MultiVecAdapter<
  MultiVector<Scalar,
              LocalOrdinal,
              GlobalOrdinal,
              Node> >::get2dViewNonConst( bool local ) const
{
  if( local ){
    localize();
    /* Use the global element list returned by
     * mv_->getMap()->getNodeElementList() to get a subCopy of mv_ which we
     * assign to l_mv_, then return get2dViewNonConst() of l_mv_
     */
    if(l_l_mv_.is_null() ){
      Teuchos::ArrayView<const GlobalOrdinal> nodeElements_go
        = mv_->getMap()->getNodeElementList();
      Teuchos::Array<size_t> nodeElements_st(nodeElements_go.size());

      // Convert the node element to a list of size_t type objects
      typename Teuchos::ArrayView<const GlobalOrdinal>::iterator it_go;
      Teuchos::Array<size_t>::iterator it_st = nodeElements_st.begin();
      for( it_go = nodeElements_go.begin(); it_go != nodeElements_go.end(); ++it_go ){
        *(it_st++) = Teuchos::as<size_t>(*it_go);
      }

      // To be consistent with the globalize() function, get a view of the local mv
      l_l_mv_ = l_mv_->subViewNonConst(nodeElements_st);

      return(l_l_mv_->get2dViewNonConst());
    } else {
      // We need to re-import values to the local, since changes may have been
      // made to the global structure that are not reflected in the local
      // view.
      Teuchos::ArrayView<const GlobalOrdinal> nodeElements_go
        = mv_->getMap()->getNodeElementList();
      Teuchos::Array<size_t> nodeElements_st(nodeElements_go.size());

      // Convert the node element to a list of size_t type objects
      typename Teuchos::ArrayView<const GlobalOrdinal>::iterator it_go;
      Teuchos::Array<size_t>::iterator it_st = nodeElements_st.begin();
      for( it_go = nodeElements_go.begin(); it_go != nodeElements_go.end(); ++it_go ){
        *(it_st++) = Teuchos::as<size_t>(*it_go);
      }

      return l_l_mv_->get2dViewNonConst();
    }
  } else {
    if( mv_->isDistributed() ){
      this->localize();

      return l_mv_->get2dViewNonConst();
    } else {                      // not distributed, no need to import
      return mv_->get2dViewNonConst();
    }
  }
}


template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node >
void
MultiVecAdapter<
  MultiVector<Scalar,
              LocalOrdinal,
              GlobalOrdinal,
              Node> >::globalize()
{
  if( mv_->isDistributed() ){
    l_mv_->doExport(*mv_, *importer_, Tpetra::REPLACE);
  }
}


template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, class Node>
template <typename Value_t>
void
MultiVecAdapter<
  MultiVector<Scalar,
              LocalOrdinal,
              GlobalOrdinal,
              Node> >::globalize(const Teuchos::ArrayView<Value_t>& newVals)
{
  if( mv_->isDistributed() ){
    if( l_mv_.is_null() ){
      localize();
    }
    Teuchos::ArrayRCP<Scalar> l_ptr = l_mv_->get1dViewNonConst();

    typename Teuchos::ArrayRCP<Scalar>::iterator it = l_ptr.begin();
    typename Teuchos::ArrayRCP<Scalar>::iterator end = l_ptr.end();

    typename Teuchos::ArrayView<Value_t>::iterator val_it = newVals.begin();

    // Update local view, doing conversion
    for( ; it != end; ++it ){
      *it = Teuchos::as<Scalar>(*val_it++);
    }

    l_mv_->doExport(*mv_, *importer_, Tpetra::REPLACE);
  } else {
    Teuchos::ArrayRCP<Scalar> ptr = mv_->get1dViewNonConst();

    typename Teuchos::ArrayRCP<Scalar>::iterator it = ptr.begin();
    typename Teuchos::ArrayRCP<Scalar>::iterator end = ptr.end();

    typename Teuchos::ArrayView<Value_t>::iterator val_it = newVals.begin();

    // Update view, doing conversion
    for( ; it != end; ++it ){
      *it = Teuchos::as<Scalar>(*val_it++);
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
              Node> >::localize()
{
  if( l_mv_.is_null() ){
    // create local multivector, maps, and do import
    o_map_ = mv_->getMap();
    Teuchos::RCP<const Teuchos::Comm<int> > comm = o_map_->getComm();

    size_t numVectors = getGlobalNumVectors();
    //l_map_ = Tpetra::createLocalMap<LocalOrdinal,GlobalOrdinal>(numVectors, comm);
    l_map_ = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal>(
        (Tpetra::global_size_t)numVectors,
        Teuchos::OrdinalTraits<GlobalOrdinal>::zero(),
        comm,
        Tpetra::LocallyReplicated,
        o_map_->getNode()));

    l_mv_ = Teuchos::rcp(
      new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>(l_map_,numVectors));

    importer_ = Teuchos::rcp(
      new Tpetra::Import<LocalOrdinal,GlobalOrdinal>(o_map_, l_map_));

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
