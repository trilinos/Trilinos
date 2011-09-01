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


#ifndef AMESOS2_MATRIXADAPTER_DECL_HPP
#define AMESOS2_MATRIXADAPTER_DECL_HPP

#include "Amesos2_config.h"

#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_VerbosityLevel.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <Tpetra_ConfigDefs.hpp>	// for global_size_t

// #include "Amesos2_ConcreteMatrixAdapter_decl.hpp"
#include "Amesos2_Util.hpp"
#include "Amesos2_MatrixTraits.hpp"

namespace Amesos2 {

  template <class M> class ConcreteMatrixAdapter;
  
  /**
   * \brief A Matrix adapter interface for Amesos2.
   *
   * All Amesos2 solver interfaces are expected to use this matrix
   * adapter interface to make their lives easier.  The methods have
   * been chosen to cater to a wide variety of third-party direct
   * sparse solvers' needs.
   *
   * \ingroup amesos2_matrix_adapters
   */
  template < class Matrix >
  class MatrixAdapter {

  public:

    typedef typename MatrixTraits<Matrix>::scalar_t                 scalar_t;
    typedef typename MatrixTraits<Matrix>::local_ordinal_t   local_ordinal_t;
    typedef typename MatrixTraits<Matrix>::global_ordinal_t global_ordinal_t;
    typedef typename MatrixTraits<Matrix>::node_t                     node_t;
    typedef Tpetra::global_size_t                              global_size_t;

    typedef Matrix                                                  matrix_t;
    typedef MatrixAdapter<Matrix>                                       type;
    typedef ConcreteMatrixAdapter<Matrix>                          adapter_t;

    // template<typename S, typename GO, typename GS, typename Op>
    // friend class Util::get_cxs_helper<MatrixAdapter<Matrix>,S,GO,GS,Op>;
    // template<class M, typename S, typename GO, typename GS, typename Op>
    // friend class Util::get_cxs_helper;

    MatrixAdapter(Teuchos::RCP<Matrix> m);


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
     * \param [out] nnz is the local number of nonzero entries in the
     *              representation of this matrix.  For example, if this
     *              processor has 12 of the 23 nzvals in its rows, then \c nnz
     *              will be 12 on exit.
     * \param [in]  rowmap A Tpetra::Map describing the desired distribution of
     *              the rows of the CRS representation on the calling processors.
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
		const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap,
		EStorage_Ordering ordering=ARBITRARY) const;

    /**
     * Convenience overload for the getCrs function that uses an enum
     * to describe some of the most basic distributions that could be
     * desired.
     */
    void getCrs(const Teuchos::ArrayView<scalar_t> nzval,
		const Teuchos::ArrayView<global_ordinal_t> colind,
		const Teuchos::ArrayView<global_size_t> rowptr,
		global_size_t& nnz,
		EDistribution distribution,
		EStorage_Ordering ordering=ARBITRARY) const;

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
     * \param [in]  colmap A Tpetra::Map describing the desired distribution of
     *              the columns of the CCS representation on the calling processors.
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
		const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > colmap,
		EStorage_Ordering ordering=ARBITRARY) const;

    /**
     * Convenience overload for the getCcs function that uses an enum
     * to describe some of the most basic distributions that could be
     * desired.
     */
    void getCcs(const Teuchos::ArrayView<scalar_t> nzval,
		const Teuchos::ArrayView<global_ordinal_t> rowind,
		const Teuchos::ArrayView<global_size_t> colptr,
		global_size_t& nnz,
		EDistribution distribution,
		EStorage_Ordering ordering=ARBITRARY) const;


    /// Returns the Teuchos::Comm object associated with this matrix.
    const Teuchos::RCP<const Teuchos::Comm<int> > getComm() const
    {
      return comm_;
    }

    /// Get the number of rows in this matrix
    global_size_t getGlobalNumRows() const;

    /// Get the number of columns in this matrix
    global_size_t getGlobalNumCols() const;

    /// Get the global number of non-zeros in this sparse matrix
    global_size_t getGlobalNNZ() const;

    /// Get the number of rows local to the calling process
    size_t getLocalNumRows() const;

    /// Get the number of columns local to the calling process
    size_t getLocalNumCols() const;

    /// Get the local number of non-zeros on this processor
    size_t getLocalNNZ() const;

    Teuchos::RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t> >
    getRowMap() const {
      return row_map_;
    }

    Teuchos::RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t> >
    getColMap() const {
      return col_map_;
    }

    Teuchos::RCP<const type> get(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map) const;

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
		     const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap,
		     EStorage_Ordering ordering,
		     has_special_impl hsi) const;

    void help_getCrs(const Teuchos::ArrayView<scalar_t> nzval,
		     const Teuchos::ArrayView<global_ordinal_t> colind,
		     const Teuchos::ArrayView<global_size_t> rowptr,
		     global_size_t& nnz,
		     const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap,
		     EStorage_Ordering ordering,
		     no_special_impl nsi) const;

    void do_getCrs(const Teuchos::ArrayView<scalar_t> nzval,
		   const Teuchos::ArrayView<global_ordinal_t> colind,
		   const Teuchos::ArrayView<global_size_t> rowptr,
		   global_size_t& nnz,
		   const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap,
		   EStorage_Ordering ordering,
		   row_access ra) const;

    void do_getCrs(const Teuchos::ArrayView<scalar_t> nzval,
		   const Teuchos::ArrayView<global_ordinal_t> colind,
		   const Teuchos::ArrayView<global_size_t> rowptr,
		   global_size_t& nnz,
		   const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap,
		   EStorage_Ordering ordering,
		   col_access ca) const;

    void help_getCcs(const Teuchos::ArrayView<scalar_t> nzval,
		     const Teuchos::ArrayView<global_ordinal_t> rowind,
		     const Teuchos::ArrayView<global_size_t> colptr,
		     global_size_t& nnz,
		     const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > colmap,
		     EStorage_Ordering ordering,
		     has_special_impl hsi) const;

    void help_getCcs(const Teuchos::ArrayView<scalar_t> nzval,
		     const Teuchos::ArrayView<global_ordinal_t> rowind,
		     const Teuchos::ArrayView<global_size_t> colptr,
		     global_size_t& nnz,
		     const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > colmap,
		     EStorage_Ordering ordering,
		     no_special_impl nsi) const;

    void do_getCcs(const Teuchos::ArrayView<scalar_t> nzval,
		   const Teuchos::ArrayView<global_ordinal_t> rowind,
		   const Teuchos::ArrayView<global_size_t> colptr,
		   global_size_t& nnz,
		   const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > colmap,
		   EStorage_Ordering ordering,
		   row_access ra) const;

    void do_getCcs(const Teuchos::ArrayView<scalar_t> nzval,
		   const Teuchos::ArrayView<global_ordinal_t> rowind,
		   const Teuchos::ArrayView<global_size_t> colptr,
		   global_size_t& nnz,
		   const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > colmap,
		   EStorage_Ordering ordering,
		   col_access ca) const;

  protected:
    // These methods will link to concrete implementations, and may
    // also be used by them

    /**
     * \param [out] row the global matrix row
     * \param [out] indices global column indices
     * \param [out] vals the non-zero values in row \c row
     * \param [out] nnz the number of nonzeros extracted from row \c row
     */
    void getGlobalRowCopy(global_ordinal_t row,
			  const Teuchos::ArrayView<global_ordinal_t>& indices,
			  const Teuchos::ArrayView<scalar_t>& vals,
			  size_t& nnz) const;

    /**
     * \param [out] col the global matrix col
     * \param [out] indices global column indices
     * \param [out] vals the non-zero values in row \c row
     * \param [out] nnz the number of nonzeros extracted from row \c row
     */
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

    bool isLocallyIndexed() const;

    bool isGloballyIndexed() const;

  protected:
    const Teuchos::RCP<const Matrix> mat_;

    // only need to be mutable for the initial assignment, is there
    // another way to do this?
    mutable Teuchos::RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > row_map_;

    mutable Teuchos::RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > col_map_;

    mutable Teuchos::RCP<const Teuchos::Comm<int> > comm_;
    
  };				// end class MatrixAdapter


  // Factory creation method
  template <class Matrix>
  Teuchos::RCP<MatrixAdapter<Matrix> >
  createMatrixAdapter(Teuchos::RCP<Matrix> m);
    
  template <class Matrix>
  Teuchos::RCP<const MatrixAdapter<Matrix> >
  createConstMatrixAdapter(Teuchos::RCP<const Matrix> m);

} // end namespace Amesos2

#endif	// AMESOS2_MATRIXADAPTER_DECL_HPP
