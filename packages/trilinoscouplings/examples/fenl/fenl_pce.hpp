/*
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
*/

#include "Stokhos_Tpetra_UQ_PCE.hpp"
#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"

#if defined( HAVE_STOKHOS_BELOS )
#include "Belos_TpetraAdapter_UQ_PCE.hpp"
#endif

#if defined( HAVE_STOKHOS_MUELU )
#include "Stokhos_MueLu_UQ_PCE.hpp"
#endif

#include <Kokkos_Core.hpp>
#include <HexElement.hpp>
#include <fenl.hpp>
#include <fenl_functors_pce.hpp>
#include <fenl_impl.hpp>
#include <MeanBasedPreconditioner.hpp>

namespace Kokkos {
namespace Example {

#if defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

  //! Get mean values matrix for mean-based preconditioning
  /*! Specialization for Sacado::UQ::PCE
   */
  template <class Storage, class ... P>
  class GetMeanValsFunc< Kokkos::View< Sacado::UQ::PCE<Storage>*,
                                       P... > > {
  public:
    typedef Sacado::UQ::PCE<Storage> Scalar;
    typedef Kokkos::View< Scalar*, P... > ViewType;
    typedef ViewType MeanViewType;
    typedef typename ViewType::execution_space execution_space;
    typedef typename ViewType::size_type size_type;

    GetMeanValsFunc(const ViewType& vals_) : vals(vals_)
    {
      const size_type nnz = vals.dimension_0();
      typename Scalar::cijk_type mean_cijk =
        Stokhos::create_mean_based_product_tensor<execution_space, typename Storage::ordinal_type, typename Storage::value_type>();
      mean_vals = Kokkos::make_view<ViewType>("mean-values", mean_cijk, nnz, 1);
      Kokkos::parallel_for( nnz, *this );
    }

    KOKKOS_INLINE_FUNCTION
    void operator() (const size_type i) const
    {
      mean_vals(i) = vals(i).fastAccessCoeff(0);
    }

     MeanViewType getMeanValues() const { return mean_vals; }

  private:
    MeanViewType mean_vals;
    ViewType vals;
  };

#else

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

    GetMeanValsFunc(const ViewType& vals_) : vals(vals_)
    {
      const size_type nnz = vals.dimension_0();
      mean_vals = ViewType("mean-values", vals.cijk(), nnz, 1);
      Kokkos::parallel_for( nnz, *this );
    }

    KOKKOS_INLINE_FUNCTION
    void operator() (const size_type i) const
    {
      mean_vals(i) = vals(i).fastAccessCoeff(0);
    }

     MeanViewType getMeanValues() const { return mean_vals; }

  private:
    MeanViewType mean_vals;
    ViewType vals;
  };

#endif

  /*!
   * \brief A stochastic preconditioner based on applying the inverse of the
   * mean.  Specialized for UQ::PCE
   */
  template<class Storage, class LO, class GO, class N>
  class MeanBasedPreconditioner<Sacado::UQ::PCE<Storage>,LO,GO,N> :
    public SGPreconditioner<Sacado::UQ::PCE<Storage>, LO, GO, N> {
  public:

    typedef Sacado::UQ::PCE<Storage> Scalar;

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
      typedef typename Scalar::cijk_type Cijk;

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

      Cijk cijk = Kokkos::getGlobalCijkTensor<Cijk>();
      Cijk mean_cijk =
        Stokhos::create_mean_based_product_tensor<typename Storage::execution_space,typename Storage::ordinal_type,typename Storage::value_type>();
      Kokkos::setGlobalCijkTensor(mean_cijk);

      KokkosMatrixType mean_kokkos_matrix(
        "mean-matrix", ncols, mean_matrix_values, kokkos_graph);
      Teuchos::RCP < MatrixType > M =
          Teuchos::rcp( new MatrixType(rmap, cmap, mean_kokkos_matrix) );

      Teuchos::RCP< OperatorType > M_op = M;
      Teuchos::RCP< PreconditionerType > mueluPreconditioner;
      mueluPreconditioner =
        MueLu::CreateTpetraPreconditioner(M_op,*precParams,coords);

      Kokkos::setGlobalCijkTensor(cijk);

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
