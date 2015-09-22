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
    
    RCP<const MatrixAdapter<matrix_t> > get_impl(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map) const;
    
  };

} // end namespace Amesos2

#endif	// AMESOS2_EPETRACRSMATRIX_MATRIXADAPTER_DECL_HPP
