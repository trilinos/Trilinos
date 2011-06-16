#ifndef AMESOS2_TPETRAROWMATRIX_MATRIXADAPTER_DECL_HPP
#define AMESOS2_TPETRAROWMATRIX_MATRIXADAPTER_DECL_HPP

#include "Amesos2_config.h"

#include <Teuchos_ArrayView.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "Amesos2_AbstractConcreteMatrixAdapter.hpp"
#include "Amesos2_MatrixAdapter_decl.hpp"
#include "Amesos2_Util.cpp"

namespace Amesos {

  template <typename Scalar,
	    typename LocalOrdinal,
	    typename GlobalOrdinal,
	    typename Node,
	    class DerivedMat >
  class AbstractConcreteMatrixAdapter< Tpetra::RowMatrix<Scalar,
							 LocalOrdinal,
							 GlobalOrdinal,
							 Node >,
				       DerivedMat >
    : public MatrixAdapter< DerivedMat > {
  public:
    // Give our base class access to our private implementation functions
    friend class MatrixAdapter< DerivedMat >;

    typedef Tpetra::RowMatrix<Scalar,
			      LocalOrdinal,
			      GlobalOrdinal,
			      Node >              matrix_t;
						 
    typedef Scalar                                scalar_t;
    typedef LocalOrdinal                   local_ordinal_t;
    typedef GlobalOrdinal                 global_ordinal_t;
    typedef Node                                    node_t;

  private:
    typedef MatrixAdapter< DerivedMat >            super_t;
    
  public:
    typedef typename super_t::global_size_t  global_size_t;

    typedef AbstractConcreteMatrixAdapter<matrix_t,
					  DerivedMat> type;

    typedef Util::no_special_impl             get_crs_spec;
    typedef Util::no_special_impl             get_ccs_spec;
    typedef Util::row_access                  major_access;

    AbstractConcreteMatrixAdapter(RCP<matrix_t> m);

  public:			// these functions should technically be private

    // implementation functions
    void getGlobalRowCopy_impl(global_ordinal_t row,
			       const Teuchos::ArrayView<global_ordinal_t>& indices,
			       const Teuchos::ArrayView<scalar_t>& vals,
			       size_t& nnz) const;
    
    void getGlobalColCopy_impl(global_ordinal_t col,
			       const Teuchos::ArrayView<global_ordinal_t>& indices,
			       const Teuchos::ArrayView<scalar_t>& vals,
			       size_t& nnz) const;

    global_size_t getGlobalNNZ_impl() const;

    size_t getMaxRowNNZ_impl() const;

    size_t getMaxColNNZ_impl() const;

    size_t getGlobalRowNNZ_impl(global_ordinal_t row) const;

    size_t getLocalRowNNZ_impl(local_ordinal_t row) const;
    
    size_t getGlobalColNNZ_impl(global_ordinal_t col) const;

    size_t getLocalColNNZ_impl(local_ordinal_t col) const;

    // Brunt of the work is put on the implementation for converting
    // their maps to a Tpetra::Map
    const RCP<const Tpetra::Map<local_ordinal_t,
				global_ordinal_t,
				node_t> >
    getRowMap_impl() const;

    const RCP<const Tpetra::Map<local_ordinal_t,
				global_ordinal_t,
				node_t> >
    getColMap_impl() const;

    const RCP<const Teuchos::Comm<int> > getComm_impl() const;

    bool isLocallyIndexed_impl() const;

    bool isGloballyIndexed_impl() const;

    // Because instantiation of the subclasses could be wildly
    // different (cf subclasses of Epetra_RowMatrix), this method
    // hands off implementation to the adapter for the subclass
    RCP<const super_t> get_impl(EDistribution d) const;

  };

} // end namespace Amesos

#endif	// AMESOS2_TPETRAROWMATRIX_MATRIXADAPTER_DECL_HPP
