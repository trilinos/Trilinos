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
 * \file   Amesos2_EpetraRowMatrix_AbstractMatrixAdapter_decl.hpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Tue Jun 14 17:17:23 2011
 * 
 * \brief  Provides the Epetra_RowMatrix abstraction for the concrete
 *         Epetra row matric adapters.
 */


#ifndef AMESOS2_EPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DECL_HPP
#define AMESOS2_EPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DECL_HPP

#include "Amesos2_config.h"

#include <Teuchos_ArrayView.hpp>

#include <Epetra_RowMatrix.h>

#include "Amesos2_AbstractConcreteMatrixAdapter.hpp"
#include "Amesos2_MatrixAdapter_decl.hpp"
#include "Amesos2_MatrixTraits.hpp"
#include "Amesos2_Util.hpp"

namespace Amesos2 {

  using Teuchos::RCP;
  
  /**
   * \brief Amesos2::MatrixAdapter definitions for objects deriving
   * from Epetra_RowMatrix.
   *
   * This class provides definitions for classes that derive
   * from/implement the Epetra_RowMatrix interface.  Most methods
   * required for compliance with the Amesos2::MatrixAdapter interface
   * are defined here.  The only method that derived class must define
   * is the get() method, which relies on each derived object knowing
   * how to construct an instance of itself (something which the
   * abstract base class cannot know).
   *
   * \ingroup amesos2_matrix_adapters
   */
  template < class DerivedMat >
  class AbstractConcreteMatrixAdapter< Epetra_RowMatrix, DerivedMat >
    : public MatrixAdapter< DerivedMat > {
  public:
    // Give our base class access to our private implementation functions
    friend class MatrixAdapter< DerivedMat >;
						 
    typedef MatrixTraits<Epetra_RowMatrix>::scalar_t                 scalar_t;
    typedef MatrixTraits<Epetra_RowMatrix>::local_ordinal_t   local_ordinal_t;
    typedef MatrixTraits<Epetra_RowMatrix>::global_ordinal_t global_ordinal_t;
    typedef MatrixTraits<Epetra_RowMatrix>::node_t                     node_t;

    typedef DerivedMat                                               matrix_t;

  private:
    typedef MatrixAdapter< DerivedMat >                               super_t;
    
  public:
    typedef typename super_t::global_size_t                     global_size_t;

    typedef AbstractConcreteMatrixAdapter<matrix_t, DerivedMat>          type;

    // subclasses should override these typedef's in case of specialization
    typedef no_special_impl                                      get_crs_spec;
    typedef no_special_impl                                      get_ccs_spec;
    typedef MatrixTraits<Epetra_RowMatrix>::major_access         major_access;

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

    size_t getLocalNNZ_impl() const;

    global_size_t getGlobalNumRows_impl() const;

    global_size_t getGlobalNumCols_impl() const;

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
    RCP<const super_t> get_impl(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map) const;

  };

} // end namespace Amesos2

#endif	// AMESOS2_EPETRAROWMATRIX_ABSTRACTMATRIXADAPTER_DECL_HPP
