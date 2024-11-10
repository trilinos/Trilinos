// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef AMESOS2_MATRIXADAPTER_DECL_HPP
#define AMESOS2_MATRIXADAPTER_DECL_HPP

#include "Amesos2_config.h"

#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_VerbosityLevel.hpp>
#include <Teuchos_FancyOStream.hpp>

#include <Tpetra_ConfigDefs.hpp>        // for global_size_t

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

    typedef typename MatrixTraits<Matrix>::global_host_idx_type  global_host_idx_t;
    typedef typename MatrixTraits<Matrix>::global_host_val_type  global_host_val_t;

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
     * \param [in]  distribution (optional: default = ROOTED 
     *              - only CONTIGUOUS_AND_ROOTED has an effect behavior)
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
    template<typename KV_S, typename KV_GO, typename KV_GS>
    void getCrs_kokkos_view(KV_S & nzval,
                            KV_GO & colind,
                            KV_GS & rowptr,
                            global_size_t& nnz,
                            const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap,
                            EStorage_Ordering ordering=ARBITRARY,
                            EDistribution distribution=ROOTED) const; // This was placed as last argument to preserve API


    /**
     * Convenience overload for the getCrs function that uses an enum
     * to describe some of the most basic distributions that could be
     * desired.
     */
    template<typename KV_S, typename KV_GO, typename KV_GS>
    void getCrs_kokkos_view(KV_S & nzval,
                            KV_GO & colind,
                            KV_GS & rowptr,
                            global_size_t& nnz,
                            EDistribution distribution=ROOTED,
                            EStorage_Ordering ordering=ARBITRARY) const; // This was placed as last argument to preserve API

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
     * \param [in]  distribution (optional: default = ROOTED 
     *              - only CONTIGUOUS_AND_ROOTED has an effect behavior)
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
    template<typename KV_S, typename KV_GO, typename KV_GS>
    void getCcs_kokkos_view(KV_S & nzval,
                            KV_GO & rowind,
                            KV_GS & colptr,
                            global_size_t& nnz,
                            const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > colmap,
                            EStorage_Ordering ordering=ARBITRARY,
                            EDistribution distribution=ROOTED) const; // This was placed as last argument to preserve API

    /**
     * Convenience overload for the getCcs function that uses an enum
     * to describe some of the most basic distributions that could be
     * desired.
     */
    template<typename KV_S, typename KV_GO, typename KV_GS>
    void getCcs_kokkos_view(KV_S & nzval,
                            KV_GO & rowind,
                            KV_GS & colptr,
                            global_size_t& nnz,
                            EDistribution distribution=ROOTED,
                            EStorage_Ordering ordering=ARBITRARY) const; // This was placed as last argument to preserve API


    /// Returns the Teuchos::Comm object associated with this matrix.
    const Teuchos::RCP<const Teuchos::Comm<int> > getComm() const
    {
      return comm_;
    }

    /// Get the number of rows in this matrix
    global_size_t getGlobalNumRows() const;

    /// Get the number of columns in this matrix
    global_size_t getGlobalNumCols() const;

    /// Get the indexbase for the row map
    global_size_t getRowIndexBase() const;

    /// Get the indexbase for the column map
    global_size_t getColumnIndexBase() const;

    /// Get the global number of non-zeros in this sparse matrix
    global_size_t getGlobalNNZ() const;

    /// Get the number of rows local to the calling process
    size_t getLocalNumRows() const;

    /// Get the number of columns local to the calling process
    size_t getLocalNumCols() const;

    /// Get the local number of non-zeros on this processor
    size_t getLocalNNZ() const;

    Teuchos::RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> >
    getMap() const {
      return static_cast<const adapter_t*>(this)->getMap_impl();
    }

    Teuchos::RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> >
    getRowMap() const {
      return row_map_;
    }

    Teuchos::RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> >
    getColMap() const {
      return col_map_;
    }

    Teuchos::RCP<const type> get(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map, EDistribution distribution = ROOTED) const;

    /// Returns a short description of this Solver
    std::string description() const;

    /// Describes of this matrix adapter with some level of verbosity.
    void describe(Teuchos::FancyOStream &out,
                  const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

    /// Return kokkos view of CRS row pointer of matrixA_
    template<typename KV>
    void returnRowPtr_kokkos_view(KV & view) const;

    /// Return kokkos view of CRS column indices of matrixA_
    template<typename KV>
    void returnColInd_kokkos_view(KV & view) const;

    /// Return kokkos view of CRS values of matrixA_
    template<typename KV>
    void returnValues_kokkos_view(KV & view) const;


  private:
    template<typename KV_S, typename KV_GO, typename KV_GS>
    void help_getCrs_kokkos_view(KV_S & nzval,
         KV_GO & colind,
         KV_GS & rowptr,
         global_size_t& nnz,
         const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap,
         EDistribution distribution,
         EStorage_Ordering ordering,
         no_special_impl nsi) const;

    template<typename KV_S, typename KV_GO, typename KV_GS>
    void do_getCrs_kokkos_view(KV_S & nzval,
        KV_GO & colind,
        KV_GS & rowptr,
        global_size_t& nnz,
        const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap,
        EDistribution distribution,
        EStorage_Ordering ordering,
        row_access ra) const;

    template<typename KV_S, typename KV_GO, typename KV_GS>
    void help_getCcs_kokkos_view(KV_S & nzval,
         KV_GO & colind,
         KV_GS & rowptr,
         global_size_t& nnz,
         const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap,
         EDistribution distribution,
         EStorage_Ordering ordering,
         no_special_impl nsi) const;

    template<typename KV_S, typename KV_GO, typename KV_GS>
    void do_getCcs_kokkos_view(KV_S & nzval,
        KV_GO & rowind,
        KV_GS & colptr,
        global_size_t& nnz,
        const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap,
        EDistribution distribution,
        EStorage_Ordering ordering,
        row_access ra) const;

  protected:
    // These methods will link to concrete implementations, and may
    // also be used by them

    /**
     * \param [out] row the global matrix row
     * \param [out] indices global column indices
     * \param [out] vals the non-zero values in row \c row
     * \param [out] nnz the number of nonzeros extracted from row \c row
     */
    template<typename KV_GO, typename KV_S>
    void getGlobalRowCopy_kokkos_view(global_ordinal_t row,
                                      KV_GO & indices,
                                      KV_S & vals,
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

    mutable Teuchos::RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > row_map_;

    mutable Teuchos::RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > col_map_;

    mutable Teuchos::RCP<const Teuchos::Comm<int> > comm_;
  };                                // end class MatrixAdapter


  // Factory creation method
  template <class Matrix>
  Teuchos::RCP<MatrixAdapter<Matrix> >
  createMatrixAdapter(Teuchos::RCP<Matrix> m);
    
  template <class Matrix>
  Teuchos::RCP<const MatrixAdapter<Matrix> >
  createConstMatrixAdapter(Teuchos::RCP<const Matrix> m);

} // end namespace Amesos2

#endif        // AMESOS2_MATRIXADAPTER_DECL_HPP
