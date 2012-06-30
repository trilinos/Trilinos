//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef KOKKOS_CRSMATRIXBASE_HPP
#define KOKKOS_CRSMATRIXBASE_HPP

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsGraphBase.hpp"

#include <Teuchos_ScalarTraits.hpp>

namespace Kokkos {

  /*! @class CrsMatrixBase
      @brief An abstract base class providing a template for Kokkos-level sparse matrix objects.
      @ingroup kokkos_crs_ops

      The matrix data is packed into a single contiguous array. The array is of length num-total-nonzeros.  Two auxiliary arrays
      give indices to the beginning of each row, according to the standard compressed sparse row (CSR) format.

      \tparam Scalar  Defines the type of the matrix valeus;
                      same as the Scalar template parameter of the encapsulating Tpetra::CrsMatrix object.
      \tparam Ordinal Defines the type of the column indices;
                      same as the LocalOrdinal template parameter of the encapsulating Tpetra::CrsMatrix object.
      \tparam Node    Kokkos Node type; same as the Node template parameter of the encapsulating Tpetra::CrsMatrix object.
  */
  template <class Scalar,
            class Ordinal,
            class Node>
  class CrsMatrixBase {
  public:

    typedef Scalar                        scalar_type;
    typedef Ordinal                       ordinal_type;
    typedef Node                          node_type;

    //! @name Constructors/Destructor
    //@{

    //! Constructor
    CrsMatrixBase(const RCP<const CrsGraphBase<Ordinal,Node> > &graph,
                  const RCP<ParameterList> &params);

    //! CrsMatrixBase Destructor
    virtual ~CrsMatrixBase();

    //@}

    //! @name Data entry and accessor methods.
    //@{

    //! Submit the matrix values.
    /**
        Must be congruous with the associated graph.

        This is used by Tpetra.
     */
    virtual void setValues(const ArrayRCP<const Scalar> &vals) = 0;

    //@}

  private:
    RCP<const CrsGraphBase<Ordinal,Node> > graph_;
    RCP<ParameterList> params_;
  };

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  CrsMatrixBase<Scalar,Ordinal,Node>::CrsMatrixBase(
      const RCP<const CrsGraphBase<Ordinal,Node> > &graph,
      const RCP<ParameterList> &params)
  : graph_(graph)
  , params_(params)
  {}

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  CrsMatrixBase<Scalar,Ordinal,Node>::~CrsMatrixBase() {
  }

} // namespace Kokkos

#endif /* KOKKOS_CRSMATRIXBASE_HPP */
