#ifndef AMESOS2_TPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP
#define AMESOS2_TPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP

namespace Amesos {

  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    typename LocalMatOps >
  ConcreteMatrixAdapter<
    Tpetra::CrsMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node,
		      LocalMatOps>
    >::ConcreteMatrixAdapter(RCP<matrix_t> m) 
      : AbstractConcreteMatrixAdapter<Tpetra::RowMatrix<Scalar,
							LocalOrdinal,
							GlobalOrdinal,
							Node>,
				      Tpetra::CrsMatrix<Scalar,
							LocalOrdinal,
							GlobalOrdinal,
							Node,
							LocalMatOps> >(m) // with implicit cast
    {}

  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    typename LocalMatOps >
  Teuchos::RCP<const MatrixAdapter<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > >
  ConcreteMatrixAdapter<
    Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>
    >::get_impl(EDistribution d) const
    {
      Teuchos::RCP<const Tpetra::Map<local_ordinal_t, global_ordinal_t, node_t > > o_map, l_map;
      o_map = this->getRowMap();
      global_size_t numrows = this->getGlobalNumRows();
      
      if( d == Util::Rooted ){
	if( Teuchos::size(*(this->comm_)) == 1 ){ // then mat_ (and this adapter) is already "rooted"
	  return( rcpFromRef(*this) );
	}

	global_size_t num_my_elements;
	if( this->getComm()->getRank() == 0 ){
	  num_my_elements = numrows;
	} else {
	  num_my_elements = Teuchos::OrdinalTraits<global_size_t>::zero();
	}
	l_map = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>( numrows,
									       num_my_elements,
									       o_map->getIndexBase(),
									       this->getComm(),
									       o_map->getNode()));
      } else if( d == Util::Globally_Replicated ){
	l_map = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>( numrows,
									       o_map->getIndexBase(),
									       this->getComm(),
									       Tpetra::LocallyReplicated,
									       o_map->getNode()));
      }
      
      Teuchos::RCP<matrix_t> l_mat;
      l_mat = Teuchos::rcp(new matrix_t(l_map, this->getMaxRowNNZ()));
      
      Teuchos::RCP<Tpetra::Import<local_ordinal_t, global_ordinal_t, node_t> > importer;
      importer = Teuchos::rcp(new Tpetra::Import<local_ordinal_t, global_ordinal_t, node_t>(o_map, l_map));
      
      l_mat->doImport(*(this->mat_), *importer, Tpetra::REPLACE);
      
      return( rcp(new ConcreteMatrixAdapter<matrix_t>(l_mat)) );
    }

} // end namespace Amesos

#endif	// AMESOS2_TPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP
