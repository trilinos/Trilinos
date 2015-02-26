// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_TPETRA_UTILITIES_HPP
#define STOKHOS_TPETRA_UTILITIES_HPP

#include "Stokhos_Tpetra_UQ_PCE.hpp"
#include "Stokhos_Tpetra_MP_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

namespace Stokhos {

  //! Get mean values matrix for mean-based preconditioning
  /*! Default implementation for all scalar types where "mean" is the same
   * as the scalar type.
   */
  template <class ViewType>
  class GetMeanValsFunc {
  public:
    typedef ViewType MeanViewType;
    typedef typename ViewType::execution_space execution_space;
    typedef typename ViewType::size_type size_type;

    GetMeanValsFunc(const ViewType& vals) {
      mean_vals = ViewType("mean-values", vals.dimension_0());
      Kokkos::deep_copy( mean_vals, vals );
    }

    MeanViewType getMeanValues() const { return mean_vals; }

  private:
    MeanViewType mean_vals;
  };

  //! Get mean values matrix for mean-based preconditioning
  /*! Specialization for Sacado::UQ::PCE
   */
  template <class Storage, class Layout, class Memory, class Device>
  class GetMeanValsFunc< Kokkos::View< Sacado::UQ::PCE<Storage>*,
                                       Layout, Memory, Device > > {
  public:
    typedef Sacado::UQ::PCE<Storage> Scalar;
    typedef Kokkos::View< Scalar*, Layout, Memory, Device > ViewType;
    typedef ViewType MeanViewType;
    typedef typename ViewType::execution_space execution_space;
    typedef typename ViewType::size_type size_type;

    GetMeanValsFunc(const ViewType& vals_) : vals(vals_) {
      const size_type nnz = vals.dimension_0();
      mean_vals = ViewType("mean-values", vals.cijk(), nnz, 1);
      Kokkos::parallel_for( nnz, *this );
    }

    KOKKOS_INLINE_FUNCTION
    void operator() (const size_type i) const {
      mean_vals(i) = vals(i).fastAccessCoeff(0);
    }

    MeanViewType getMeanValues() const { return mean_vals; }

  private:
    MeanViewType mean_vals;
    ViewType vals;
  };

  //! Get mean values matrix for mean-based preconditioning
  /*! Specialization for Sacado::MP::Vector
   */
  template <class Storage, class Layout, class Memory, class Device>
  class GetMeanValsFunc< Kokkos::View< Sacado::MP::Vector<Storage>*,
                                       Layout, Memory, Device > > {
  public:
    typedef Sacado::MP::Vector<Storage> Scalar;
    typedef Kokkos::View< Scalar*, Layout, Memory, Device > ViewType;
    typedef ViewType MeanViewType;
    typedef typename ViewType::execution_space execution_space;
    typedef typename ViewType::size_type size_type;

    GetMeanValsFunc(const ViewType& vals_) :
      vals(vals_), vec_size(vals.sacado_size())
    {
      const size_type nnz = vals.dimension_0();
      mean_vals = ViewType("mean-values", nnz, 1);
      Kokkos::parallel_for( nnz, *this );
    }

    KOKKOS_INLINE_FUNCTION
    void operator() (const size_type i) const
    {
      typename Scalar::value_type s = 0.0;
      for (size_type j=0; j<vec_size; ++j)
        s += vals(i).fastAccessCoeff(j);
      mean_vals(i) = s;
    }

    MeanViewType getMeanValues() const { return mean_vals; }

  private:
    MeanViewType mean_vals;
    ViewType vals;
    const size_type vec_size;
  };

  template <typename Scalar, typename LO, typename GO, typename N>
  Teuchos::RCP< Tpetra::CrsMatrix<Scalar,LO,GO,N> >
  build_mean_matrix(const Tpetra::CrsMatrix<Scalar,LO,GO,N>& A)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef Tpetra::CrsMatrix<Scalar,LO,GO,N> MatrixType;
    typedef Tpetra::Map<LO,GO,N> Map;

    typedef Kokkos::CrsMatrix<Scalar, LO, typename N::execution_space, void, size_t> KokkosMatrixType;

    typedef typename KokkosMatrixType::StaticCrsGraphType KokkosGraphType;
    typedef typename KokkosMatrixType::values_type KokkosMatrixValuesType;

    RCP< const Map > rmap = A.getRowMap();
    RCP< const Map > cmap = A.getColMap();

    KokkosMatrixType kokkos_matrix = A.getLocalMatrix();
    KokkosGraphType kokkos_graph = kokkos_matrix.graph;
    KokkosMatrixValuesType matrix_values = kokkos_matrix.values;
    const size_t ncols = kokkos_matrix.numCols();
    typedef GetMeanValsFunc <KokkosMatrixValuesType > MeanFunc;
    typedef typename MeanFunc::MeanViewType KokkosMeanMatrixValuesType;
    MeanFunc meanfunc(matrix_values);
    KokkosMeanMatrixValuesType mean_matrix_values = meanfunc.getMeanValues();

    // From here on we are assuming that
    // KokkosMeanMatrixValuesType == KokkosMatrixValuestype

    KokkosMatrixType mean_kokkos_matrix(
      "mean-matrix", ncols, mean_matrix_values, kokkos_graph);
    RCP < MatrixType > mean_matrix =
      rcp( new MatrixType(rmap, cmap, mean_kokkos_matrix) );
    return mean_matrix;
  }

}

#endif // STOKHOS_TPETRA_UTILITIES_HPP
