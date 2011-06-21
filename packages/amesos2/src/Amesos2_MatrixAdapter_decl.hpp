#ifndef AMESOS2_MATRIXADAPTER_DECL_HPP
#define AMESOS2_MATRIXADAPTER_DECL_HPP

#include "Amesos2_config.h"

#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_VerbosityLevel.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <Tpetra_ConfigDefs.hpp>	// for global_size_t

#include "Amesos2_ConcreteMatrixAdapter.hpp"
#include "Amesos2_Util.hpp"
#include "Amesos2_MatrixTraits.hpp"

namespace Amesos {

  using Amesos::Util::has_special_impl;
  using Amesos::Util::no_special_impl;
  using Amesos::Util::row_access;
  using Amesos::Util::col_access;
  
  using Teuchos::RCP;

  typedef enum {
    Distributed,                /**< no processor has a view of the entire matrix, only local pieces */
    Distributed_No_Overlap,     /**< no row or column may be present on more than one processor */
    Globally_Replicated,        /**< each processor has a view of the entire matrix */
    Rooted                     /**< only \c rank=0 has a full view, all others have nothing. */
    // SameDistribution            /**< Use whatever distribution the matrix adapter currently has */
  } EDistribution;

  typedef enum {
    Sorted_Indices,             /**< row/col indices need to appear in sorted order */
    Arbitrary                   /**< index order can be arbitrary */
  } EStorage_Ordering;

  /**
   * \brief A Matrix adapter interface for Amesos2.
   *
   * All Amesos2 solver interfaces are expected to use this matrix
   * adapter interface to make their lives easier.  The methods have
   * been chosen to cater to a wide variety of third-party direct
   * sparse solvers' needs.
   */
  template < class Matrix >
  class MatrixAdapter {
  public:

    typedef typename MatrixTraits<Matrix>::scalar_t                 scalar_t;
    typedef typename MatrixTraits<Matrix>::local_ordinal_t   local_ordinal_t;
    typedef typename MatrixTraits<Matrix>::global_ordinal_t global_ordinal_t;
    typedef typename MatrixTraits<Matrix>::node_t                     node_t;
    // typedef typename MatrixTraits<Matrix>::local_mat_ops_t   local_mat_ops_t;
    typedef Tpetra::global_size_t                              global_size_t;

    typedef MatrixAdapter<Matrix> type;
    typedef ConcreteMatrixAdapter<Matrix> adapter_t;

    MatrixAdapter(RCP<Matrix> m);


    /**
     * \brief Gets a compressed-row storage summary of \c this
     *
     * Extracts a compressed-row storage format of the matrix and stores the
     * information in the user-supplied containers.
     *
     * \param [out] nzval will hold the values of the nonzero entries of \c this
     * \param [out] colind will hold the column indices of \c this for each row.
     * \param [out] rowptr is of size <tt>nrow + 1</tt> and <tt>rowptr[j]</tt>
     *              stores the location in \c nzval and \c colind which starts
     *              row \c j of \c this.  <tt>rowptr[nrow] = nnz</tt>, where \c
     *              nrow is the number of rows in this matrix.
     * \param [out] nnz is the number of nonzero entries in this matrix.
     * \param [in]  distribution
     * \param [in]  ordering
     *
     * \exception std::length_error Thrown if \c nzval or \c colind is not
     * large enough to hold the global number of nonzero values.
     *
     * \exception std::length_error Thrown if \c rowptr is not at least
     * <tt>nrow + 1</tt> in size, the required size.
     *
     * \exception std::runtime_error Thrown if there is an error while extracting
     * row values from the underlying matrix.
     */
    void getCrs(const Teuchos::ArrayView<scalar_t> nzval,
		const Teuchos::ArrayView<global_ordinal_t> colind,
		const Teuchos::ArrayView<global_size_t> rowptr,
		global_size_t& nnz,
		EDistribution distribution,
		EStorage_Ordering ordering) const;

    /**
     * \brief Gets a compressed-column storage summary of \c this
     *
     * Extracts a compressed-column storage format of the matrix and stores the
     * information in the user-supplied containers.
     *
     * \param [out] nzval will hold the values of the nonzero entries of \c this
     * \param [out] rowind will hold the row indices of \c this for each column.
     * \param [out] colptr is of size <tt>ncol + 1</tt> and <tt>colptr[j]</tt>
     *              stores the location in \c nzval and \c rowind which starts
     *              column \c j of \c this.  <tt>colptr[ncol] = nnz</tt>, where \c
     *              ncol is the number of columns in this matrix.
     * \param [out] nnz is the number of nonzero entries in this matrix.
     * \param [in]  distribution
     * \param [in]  ordering
     *
     * \exception std::length_error Thrown if \c nzval or \c rowind is not
     * large enough to hold the global number of nonzero values.
     *
     * \exception std::length_error Thrown if \c colptr is not at least
     * <tt>ncol + 1</tt> in size, the required size.
     *
     * \exception std::runtime_error Thrown if there is an error while extracting
     * row values from the underlying matrix.
     */
    void getCcs(const Teuchos::ArrayView<scalar_t> nzval,
		const Teuchos::ArrayView<global_ordinal_t> rowind,
		const Teuchos::ArrayView<global_size_t> colptr,
		global_size_t& nnz,
		EDistribution distribution,
		EStorage_Ordering ordering) const;


    /// Returns the Teuchos::Comm object associated with this matrix.
    const RCP<const Teuchos::Comm<int> > getComm() const
    {
      return comm_;
    }

    /// Get the number of rows in this matrix
    global_size_t getGlobalNumRows() const;

    /// Get the number of columns in this matrix
    global_size_t getGlobalNumCols() const;

    /// Get the global number of non-zeros in this sparse matrix
    global_size_t getGlobalNNZ() const;

    /// Returns a short description of this Solver
    std::string description() const;

    /// Describes of this matrix adapter with some level of verbosity.
    void describe(Teuchos::FancyOStream &out,
		  const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;


  private:

    void help_getCrs(const Teuchos::ArrayView<scalar_t> nzval,
		     const Teuchos::ArrayView<global_ordinal_t> colind,
		     const Teuchos::ArrayView<global_size_t> rowptr,
		     global_size_t& nnz,
		     EDistribution distribution,
		     EStorage_Ordering ordering,
		     has_special_impl hsi) const;

    void help_getCrs(const Teuchos::ArrayView<scalar_t> nzval,
		     const Teuchos::ArrayView<global_ordinal_t> colind,
		     const Teuchos::ArrayView<global_size_t> rowptr,
		     global_size_t& nnz,
		     EDistribution distribution,
		     EStorage_Ordering ordering,
		     no_special_impl nsi) const;

    void do_getCrs(const Teuchos::ArrayView<scalar_t> nzval,
		   const Teuchos::ArrayView<global_ordinal_t> colind,
		   const Teuchos::ArrayView<global_size_t> rowptr,
		   global_size_t& nnz,
		   EDistribution distribution,
		   EStorage_Ordering ordering,
		   row_access ra) const;

    void do_getCrs(const Teuchos::ArrayView<scalar_t> nzval,
		   const Teuchos::ArrayView<global_ordinal_t> colind,
		   const Teuchos::ArrayView<global_size_t> rowptr,
		   global_size_t& nnz,
		   EDistribution distribution,
		   EStorage_Ordering ordering,
		   col_access ca) const;

    void help_getCcs(const Teuchos::ArrayView<scalar_t> nzval,
		     const Teuchos::ArrayView<global_ordinal_t> rowind,
		     const Teuchos::ArrayView<global_size_t> colptr,
		     global_size_t& nnz,
		     EDistribution distribution,
		     EStorage_Ordering ordering,
		     has_special_impl hsi) const;

    void help_getCcs(const Teuchos::ArrayView<scalar_t> nzval,
		     const Teuchos::ArrayView<global_ordinal_t> rowind,
		     const Teuchos::ArrayView<global_size_t> colptr,
		     global_size_t& nnz,
		     EDistribution distribution,
		     EStorage_Ordering ordering,
		     no_special_impl nsi) const;

    void do_getCcs(const Teuchos::ArrayView<scalar_t> nzval,
		   const Teuchos::ArrayView<global_ordinal_t> rowind,
		   const Teuchos::ArrayView<global_size_t> colptr,
		   global_size_t& nnz,
		   EDistribution distribution,
		   EStorage_Ordering ordering,
		   row_access ra) const;

    void do_getCcs(const Teuchos::ArrayView<scalar_t> nzval,
		   const Teuchos::ArrayView<global_ordinal_t> rowind,
		   const Teuchos::ArrayView<global_size_t> colptr,
		   global_size_t& nnz,
		   EDistribution distribution,
		   EStorage_Ordering ordering,
		   col_access ca) const;

  protected:
    // These methods will link to concrete implementations, and may
    // also be used by them
    void getGlobalRowCopy(global_ordinal_t row,
			  const Teuchos::ArrayView<global_ordinal_t>& indices,
			  const Teuchos::ArrayView<scalar_t>& vals,
			  size_t& nnz) const;

    void getGlobalColCopy(global_ordinal_t col,
			  const Teuchos::ArrayView<global_ordinal_t>& indices,
			  const Teuchos::ArrayView<scalar_t>& vals,
			  size_t& nnz) const;

    size_t getMaxRowNNZ() const;

    size_t getMaxColNNZ() const;

    size_t getGlobalRowNNZ(global_ordinal_t row) const;

    size_t getLocalRowNNZ(local_ordinal_t row) const;

    size_t getGlobalColNNZ(global_ordinal_t col) const;

    size_t getLocalColNNZ(local_ordinal_t col) const;

    RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t> >
    getRowMap() const {
      return row_map_;
    }

    RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t> >
    getColMap() const {
      return col_map_;
    }

    bool isLocallyIndexed() const;

    bool isGloballyIndexed() const;

    RCP<const type> get(EDistribution d) const;

  protected:
    const RCP<const Matrix> mat_;

    // only need to be mutable for the initial assignment, is there
    // another way to do this?
    mutable RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > row_map_;

    mutable RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > col_map_;

    mutable RCP<const Teuchos::Comm<int> > comm_;
    
  };				// end class MatrixAdapter


  // Factory creation method
  template<typename Matrix>
  Teuchos::RCP<MatrixAdapter<Matrix> >
  createMatrixAdapter(Teuchos::RCP<Matrix> m){
    return( rcp(new ConcreteMatrixAdapter<Matrix>(m)) );
  }

} // end namespace Amesos

#endif	// AMESOS2_MATRIXADAPTER_DECL_HPP
