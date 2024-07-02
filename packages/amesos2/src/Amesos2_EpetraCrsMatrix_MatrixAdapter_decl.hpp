// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/**
 * \file   Amesos2_EpetraCrsMatrix_MatrixAdapter_decl.hpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Tue Jun 14 17:17:00 MDT 2011
 * 
 * \brief Specialization of the ConcreteMatrixAdapter for
 * Epetra_CrsMatrix.  Inherits all its functionality from the
 * Epetra_RowMatrix specialization of \c AbstractConcreteMatrixAdapter.
 */

#ifndef AMESOS2_EPETRACRSMATRIX_MATRIXADAPTER_DECL_HPP
#define AMESOS2_EPETRACRSMATRIX_MATRIXADAPTER_DECL_HPP

#include "Amesos2_config.h"

#include <Epetra_CrsMatrix.h>

#include "Amesos2_EpetraRowMatrix_AbstractMatrixAdapter_decl.hpp"
#include "Amesos2_MatrixAdapter_decl.hpp"

namespace Amesos2 {

  // template <class M, class D> class AbstractConcreteMatrixAdapter;

  /**
   * \brief MatrixAdapter definitions for Epetra_CrsMatrix objects.
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
  template <>
  class ConcreteMatrixAdapter< Epetra_CrsMatrix >
    : public AbstractConcreteMatrixAdapter< Epetra_RowMatrix, Epetra_CrsMatrix >
  {
    // Give our matrix adapter class access to our private
    // implementation functions
    friend class MatrixAdapter< Epetra_RowMatrix >;
  public:
    typedef Epetra_CrsMatrix                               matrix_t;
  private:
    typedef AbstractConcreteMatrixAdapter<Epetra_RowMatrix,
					  Epetra_CrsMatrix> super_t;
  public:
    // 'import' superclass types
    typedef super_t::scalar_t                              scalar_t;
    typedef super_t::local_ordinal_t                local_ordinal_t;
    typedef super_t::global_ordinal_t              global_ordinal_t;
    typedef super_t::node_t                                  node_t;
    typedef super_t::global_size_t                    global_size_t;
    
    typedef ConcreteMatrixAdapter<matrix_t>                    type;
    
    ConcreteMatrixAdapter(RCP<matrix_t> m);
    
    RCP<const MatrixAdapter<matrix_t> > get_impl(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map, EDistribution distribution = ROOTED) const;
    
  };

} // end namespace Amesos2

#endif	// AMESOS2_EPETRACRSMATRIX_MATRIXADAPTER_DECL_HPP
