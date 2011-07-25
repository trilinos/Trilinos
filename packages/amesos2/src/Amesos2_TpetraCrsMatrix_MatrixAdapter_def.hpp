#ifndef AMESOS2_TPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP
#define AMESOS2_TPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP

#include "Amesos2_TpetraCrsMatrix_MatrixAdapter_decl.hpp"
#include "Amesos2_TpetraRowMatrix_AbstractMatrixAdapter_def.hpp"
#include "Amesos2_MatrixAdapter_def.hpp"

namespace Amesos2 {

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
    >::ConcreteMatrixAdapter(Teuchos::RCP<matrix_t> m) 
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
    >::get_impl(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map) const
    {
      typedef Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> map_t;
      Teuchos::RCP<matrix_t> t_mat;
      t_mat = Teuchos::rcp(new matrix_t(Teuchos::rcpFromPtr(map), this->getMaxRowNNZ()));
      
      Teuchos::RCP<Tpetra::Import<local_ordinal_t, global_ordinal_t, node_t> > importer;
      importer = Teuchos::rcp(new Tpetra::Import<local_ordinal_t, global_ordinal_t, node_t>(this->getRowMap(), Teuchos::rcpFromPtr(map)));
      
      t_mat->doImport(*(this->mat_), *importer, Tpetra::REPLACE);
      
      return( rcp(new ConcreteMatrixAdapter<matrix_t>(t_mat)) );
    }

} // end namespace Amesos2

#endif	// AMESOS2_TPETRACRSMATRIX_MATRIXADAPTER_DEF_HPP
