// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/**
 * \file   Amesos2_TpetraCrsMatrix_MatrixAdapter_decl.hpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Mon Jun 13 11:39:24 2011
 *
 * \brief Specialization of the ConcreteMatrixAdapter for
 * Tpetra::CrsMatrix.  Inherits all its functionality from the
 * Tpetra::RowMatrix specialization of \c AbstractConcreteMatrixAdapter.
 */

#ifndef AMESOS2_TPETRACRSMATRIX_MATRIXADAPTER_DECL_HPP
#define AMESOS2_TPETRACRSMATRIX_MATRIXADAPTER_DECL_HPP

#include "Amesos2_config.h"

#include "Amesos2_TpetraRowMatrix_AbstractMatrixAdapter_decl.hpp"
#include "Amesos2_MatrixAdapter_decl.hpp"

namespace Amesos2 {

  /**
   * \brief MatrixAdapter definitions for Tpetra::CrsMatrix objects.
   *
   * Defines only the get_impl() method, which returns an instance of
   * a Amesos2::MatrixAdapter whose underlying matrix has the given
   * distribution based on the Tpetra::Map.
   *
   * All other significant functionality is inherited from this
   * class's superclass.
   *
   * \ingroup amesos2_matrix_adapters
   */
  template <typename Scalar,
            typename LocalOrdinal,
            typename GlobalOrdinal,
            typename Node>
  class ConcreteMatrixAdapter<Tpetra::CrsMatrix<Scalar,
                                                LocalOrdinal,
                                                GlobalOrdinal,
                                                Node> >
    : public AbstractConcreteMatrixAdapter<Tpetra::RowMatrix<Scalar,
                                                             LocalOrdinal,
                                                             GlobalOrdinal,
                                                             Node>,
                                           Tpetra::CrsMatrix<Scalar,
                                                             LocalOrdinal,
                                                             GlobalOrdinal,
                                                             Node> >
  {
    // Give our matrix adapter class access to our private
    // implementation functions
    friend class MatrixAdapter<Tpetra::RowMatrix<Scalar,
                                                 LocalOrdinal,
                                                 GlobalOrdinal,
                                                 Node> >;
  public:
    typedef Tpetra::CrsMatrix<Scalar,
                              LocalOrdinal,
                              GlobalOrdinal,
                              Node>                  matrix_t;
  private:
    typedef AbstractConcreteMatrixAdapter<
      Tpetra::RowMatrix<Scalar,
                        LocalOrdinal,
                        GlobalOrdinal,
                        Node>, matrix_t>                super_t;
  public:
    // 'import' superclass types
    typedef typename super_t::scalar_t                        scalar_t;
    typedef typename super_t::local_ordinal_t          local_ordinal_t;
    typedef typename super_t::global_ordinal_t        global_ordinal_t;
    typedef typename super_t::node_t                            node_t;
    typedef typename super_t::global_size_t              global_size_t;

    typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>       map_t;
    typedef ConcreteMatrixAdapter<matrix_t>                       type;

    typedef Kokkos::DefaultHostExecutionSpace                 HostExecSpaceType;

    ConcreteMatrixAdapter(RCP<matrix_t> m);

    RCP<const MatrixAdapter<matrix_t> > get_impl(const Teuchos::Ptr<const map_t> map, EDistribution distribution = ROOTED) const;
    RCP<const MatrixAdapter<matrix_t> > reindex_impl(Teuchos::RCP<const map_t> &contigRowMap, Teuchos::RCP<const map_t> &contigColMap, const EPhase current_phase) const;

    template<typename KV_S, typename KV_GO, typename KV_GS, typename host_ordinal_type_array, typename host_scalar_type_array>
    LocalOrdinal gather_impl(KV_S& nzvals, KV_GO& indices, KV_GS& pointers,
                             host_ordinal_type_array &perm_g2l,
                             host_ordinal_type_array &recvCountRows, host_ordinal_type_array &recvDisplRows,
                             host_ordinal_type_array &recvCounts, host_ordinal_type_array &recvDispls,
                             host_ordinal_type_array &transpose_map, host_scalar_type_array &nzvals_t,
                             bool column_major, EPhase current_phase) const;

    //! Print a description of this adapter to the given output stream
    void
    describe (Teuchos::FancyOStream& os,
              const Teuchos::EVerbosityLevel verbLevel =
              Teuchos::Describable::verbLevel_default) const;
  };

} // end namespace Amesos2

#endif  // AMESOS2_TPETRACRSMATRIX_MATRIXADAPTER_DECL_HPP
