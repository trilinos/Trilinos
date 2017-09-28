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

#ifndef STOKHOS_MEAN_BASED_PRECONDITIONER_HPP
#define STOKHOS_MEAN_BASED_PRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "SGPreconditioner.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Tpetra_MultiVector.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

namespace Kokkos {
namespace Example {

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

    GetMeanValsFunc(const ViewType& vals)
    {
      mean_vals = ViewType("mean-values", vals.dimension_0());
      Kokkos::deep_copy( mean_vals, vals );
    }

    MeanViewType getMeanValues() const { return mean_vals; }

  private:
    MeanViewType mean_vals;
  };

  /*!
   * \brief A stochastic preconditioner based on applying the inverse of the
   * mean.
   */
  template<class Scalar, class LO, class GO, class N>
  class MeanBasedPreconditioner : public SGPreconditioner<Scalar, LO, GO, N> {
  public:

    //! Constructor
    MeanBasedPreconditioner() {}

    //! Destructor
    virtual ~MeanBasedPreconditioner() {}

    //! Setup preconditioner
    virtual
    Teuchos::RCP<Tpetra::Operator<Scalar,LO,GO,N> >
    setupPreconditioner(
      const Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO,N> >& A,
      const Teuchos::RCP<Teuchos::ParameterList>& precParams,
      const Teuchos::RCP<Tpetra::MultiVector<double,LO,GO,N> >& coords)
    {
      using Teuchos::ArrayView;
      using Teuchos::Array;
      typedef Tpetra::CrsMatrix<Scalar,LO,GO,N> MatrixType;
      typedef Tpetra::Map<LO,GO,N> Map;
      typedef Tpetra::Operator<Scalar,LO,GO,N> OperatorType;
      typedef MueLu::TpetraOperator<Scalar,LO,GO,N> PreconditionerType;

      typedef typename MatrixType::local_matrix_type KokkosMatrixType;

      typedef typename KokkosMatrixType::StaticCrsGraphType KokkosGraphType;
      typedef typename KokkosMatrixType::values_type KokkosMatrixValuesType;

      Teuchos::RCP< const Map > rmap = A->getRowMap();
      Teuchos::RCP< const Map > cmap = A->getColMap();

      KokkosMatrixType  kokkos_matrix = A->getLocalMatrix();
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
      Teuchos::RCP < MatrixType > M =
          Teuchos::rcp( new MatrixType(rmap, cmap, mean_kokkos_matrix) );

      Teuchos::RCP< OperatorType > M_op = M;
      Teuchos::RCP< PreconditionerType > mueluPreconditioner;
      mueluPreconditioner =
        MueLu::CreateTpetraPreconditioner(M_op,*precParams,coords);
      return mueluPreconditioner;
    };

  private:

    //! Private to prohibit copying
    MeanBasedPreconditioner(const MeanBasedPreconditioner&);

    //! Private to prohibit copying
    MeanBasedPreconditioner& operator=(const MeanBasedPreconditioner&);


  }; // class MeanBasedPreconditioner

}
}

#endif // STOKHOS_MEAN_BASED_PRECONDITIONER_HPP
