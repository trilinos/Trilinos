// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
      mean_vals = ViewType("mean-values", vals.extent(0));
      Kokkos::deep_copy( mean_vals, vals );
    }

    MeanViewType getMeanValues() const { return mean_vals; }

  private:
    MeanViewType mean_vals;
  };

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

    GetMeanValsFunc(const ViewType& vals_) : vals(vals_) {
      const size_type nnz = vals.extent(0);
      typename Scalar::cijk_type mean_cijk =
        Stokhos::create_mean_based_product_tensor<execution_space, typename Storage::ordinal_type, typename Storage::value_type>();
      mean_vals = Kokkos::make_view<ViewType>("mean-values", mean_cijk, nnz, 1);
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
  template <class Storage, class ... P>
  class GetMeanValsFunc< Kokkos::View< Sacado::MP::Vector<Storage>*,
                                       P... > > {
  public:
    typedef Sacado::MP::Vector<Storage> Scalar;
    typedef Kokkos::View< Scalar*, P... > ViewType;
    typedef ViewType MeanViewType;
    typedef typename ViewType::execution_space execution_space;
    typedef typename ViewType::size_type size_type;

    GetMeanValsFunc(const ViewType& vals_) :
      vals(vals_), vec_size(Kokkos::dimension_scalar(vals))
    {
      const size_type nnz = vals.extent(0);
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

  //! Get mean values matrix for mean-based preconditioning
  /*! Default implementation for all scalar types where "mean" is the same
   * as the scalar type.
   */
  template <class ViewType>
  class GetScalarMeanValsFunc {
  public:
    typedef ViewType MeanViewType;
    typedef typename ViewType::execution_space execution_space;
    typedef typename ViewType::size_type size_type;

    GetScalarMeanValsFunc(const ViewType& vals) {
      mean_vals = ViewType("mean-values", vals.extent(0));
      Kokkos::deep_copy( mean_vals, vals );
    }

    MeanViewType getMeanValues() const { return mean_vals; }

  private:
    MeanViewType mean_vals;
  };

  //! Get mean values matrix for mean-based preconditioning
  /*! Specialization for Sacado::UQ::PCE
   */
  template <class Storage, class ... P>
  class GetScalarMeanValsFunc< Kokkos::View< Sacado::UQ::PCE<Storage>*,
                                             P... > > {
  public:
    typedef Sacado::UQ::PCE<Storage> Scalar;
    typedef typename Scalar::value_type MeanScalar;
    typedef Kokkos::View< Scalar*, P... > ViewType;
    typedef Kokkos::View< MeanScalar*, P... > MeanViewType;
    typedef typename ViewType::execution_space execution_space;
    typedef typename ViewType::size_type size_type;

    GetScalarMeanValsFunc(const ViewType& vals_) : vals(vals_) {
      const size_type nnz = vals.extent(0);
      mean_vals = MeanViewType("mean-values", nnz);
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
  template <class Storage, class ... P>
  class GetScalarMeanValsFunc< Kokkos::View< Sacado::MP::Vector<Storage>*,
                                             P... > > {
  public:
    typedef Sacado::MP::Vector<Storage> Scalar;
    typedef typename Scalar::value_type MeanScalar;
    typedef Kokkos::View< Scalar*, P... > ViewType;
    typedef Kokkos::View< MeanScalar*, P... > MeanViewType;
    typedef typename ViewType::execution_space execution_space;
    typedef typename ViewType::size_type size_type;

    GetScalarMeanValsFunc(const ViewType& vals_) :
      vals(vals_), vec_size(Kokkos::dimension_scalar(vals))
    {
      const size_type nnz = vals.extent(0);
      mean_vals = ViewType("mean-values", nnz);
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
    typedef typename MatrixType::local_matrix_device_type KokkosMatrixType;
    typedef typename KokkosMatrixType::values_type KokkosMatrixValuesType;

    KokkosMatrixType kokkos_matrix = A.getLocalMatrixDevice();
    KokkosMatrixValuesType matrix_values = kokkos_matrix.values;
    const size_t ncols = kokkos_matrix.numCols();
    typedef GetMeanValsFunc <KokkosMatrixValuesType > MeanFunc;
    typedef typename MeanFunc::MeanViewType KokkosMeanMatrixValuesType;
    MeanFunc meanfunc(matrix_values);
    KokkosMeanMatrixValuesType mean_matrix_values = meanfunc.getMeanValues();

    // From here on we are assuming that
    // KokkosMeanMatrixValuesType == KokkosMatrixValuestype

    RCP < MatrixType > mean_matrix =
      rcp( new MatrixType(A.getCrsGraph(), mean_matrix_values) );
    mean_matrix->fillComplete();
    return mean_matrix;
  }

  template <typename Scalar, typename LO, typename GO, typename N>
  Teuchos::RCP< Tpetra::CrsMatrix<typename Scalar::value_type,LO,GO,N> >
  build_mean_scalar_matrix(const Tpetra::CrsMatrix<Scalar,LO,GO,N>& A)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef typename Scalar::value_type BaseScalar;
    typedef Tpetra::CrsMatrix<Scalar,LO,GO,N> MatrixType;
    typedef Tpetra::CrsMatrix<BaseScalar,LO,GO,N> ScalarMatrixType;
    typedef typename MatrixType::local_matrix_device_type KokkosMatrixType;
    typedef typename ScalarMatrixType::local_matrix_device_type ScalarKokkosMatrixType;
    typedef typename KokkosMatrixType::values_type KokkosMatrixValuesType;

    KokkosMatrixType kokkos_matrix = A.getLocalMatrixDevice();
    KokkosMatrixValuesType matrix_values = kokkos_matrix.values;
    typedef GetScalarMeanValsFunc <KokkosMatrixValuesType > MeanFunc;
    typedef typename MeanFunc::MeanViewType KokkosMeanMatrixValuesType;
    MeanFunc meanfunc(matrix_values);
    KokkosMeanMatrixValuesType mean_matrix_values = meanfunc.getMeanValues();

    // From here on we are assuming that
    // KokkosMeanMatrixValuesType == ScalarKokkosMatrixValuesType

    RCP < ScalarMatrixType > mean_matrix =
      rcp( new ScalarMatrixType(A.getCrsGraph(), mean_matrix_values) );
    mean_matrix->fillComplete();
    return mean_matrix;
  }

  namespace Impl {

  // Functor for copying a PCE view to a scalar view
  // (Assumes view is rank 2, LayoutLeft)
  template <typename ExecSpace>
  struct CopyPCE2Scalar {
    typedef ExecSpace exec_space;
    template <typename DstView, typename SrcView>
    CopyPCE2Scalar(const DstView& dst, const SrcView& src) {
      impl(dst,src);
    }

    template <typename DstView, typename SrcView>
    void impl(const DstView& dst, const SrcView& src) {
      typedef typename SrcView::non_const_value_type Scalar;
      const size_t m = src.extent(0);
      const size_t n = src.extent(1);
      const size_t p = Kokkos::dimension_scalar(src);
      Kokkos::RangePolicy<exec_space> policy(0,m);
      Kokkos::parallel_for( policy, KOKKOS_LAMBDA(const size_t i)
      {
        for (size_t j=0; j<n; ++j) {
          const Scalar& s = src(i,j);
          for (size_t k=0; k<p; ++k)
            dst(i,j*p+k) = s.fastAccessCoeff(k);
        }
      });
    }
  };

  // Functor for copying a scalar view to a PCE view
  // (Assumes view is rank 2, LayoutLeft)
  template <typename ExecSpace>
  struct CopyScalar2PCE {
    typedef ExecSpace exec_space;

    template <typename DstView, typename SrcView>
    CopyScalar2PCE(const DstView& dst, const SrcView& src) {
      impl(dst,src);
    }

    template <typename DstView, typename SrcView>
    void impl(const DstView& dst, const SrcView& src) {
      typedef typename DstView::non_const_value_type Scalar;
      const size_t m = dst.extent(0);
      const size_t n = dst.extent(1);
      const size_t p = Kokkos::dimension_scalar(dst);

      Kokkos::RangePolicy<exec_space> policy(0,m);
      Kokkos::parallel_for( policy, KOKKOS_LAMBDA(const size_t i)
      {
        for (size_t j=0; j<n; ++j) {
          Scalar& d = dst(i,j);
          for (size_t k=0; k<p; ++k)
            d.fastAccessCoeff(k) = src(i,j*p+k);
        }
      });
    }
  };

#ifdef KOKKOS_ENABLE_CUDA
  // Specialization for CopyPCE2Scalar specifically for Cuda that ensures
  // coalesced reads and writes
  template <>
  struct CopyPCE2Scalar<Kokkos::Cuda> {
    typedef Kokkos::Cuda exec_space;

    template <typename DstView, typename SrcView>
    CopyPCE2Scalar(const DstView& dst, const SrcView& src) {
      impl(dst,src);
    }

    template <typename DstView, typename SrcView>
    void impl(const DstView& dst, const SrcView& src) {
      typedef typename DstView::non_const_value_type Scalar;
      typedef Kokkos::TeamPolicy<exec_space> Policy;
      typedef typename Policy::member_type Member;

      const size_t m = src.extent(0);
      const size_t n = src.extent(1);
      const size_t p = Kokkos::dimension_scalar(src);

      const size_t ChunkSize = 16;
      const size_t M = (m+ChunkSize-1)/ChunkSize;

      Policy policy(M,ChunkSize,ChunkSize);
      Kokkos::parallel_for( policy, [=] __device__(const Member& team)
      {
        __shared__ Scalar tmp[ChunkSize][ChunkSize];
        const size_t i_block = blockIdx.x*ChunkSize;

        for (size_t j=0; j<n; ++j) {
          for (size_t k_block=0; k_block<p; k_block+=ChunkSize) {

            // Make sure previous iteration has completed before overwriting tmp
            __syncthreads();

            // Read ChunkSize x ChunkSize block (coalesced on k)
            size_t i = i_block + threadIdx.y;
            size_t k = k_block + threadIdx.x;
            if (i < m && k < p)
              tmp[threadIdx.y][threadIdx.x] = src(i,j).fastAccessCoeff(k);

            // Wait for all threads to finish
            __syncthreads();

            // Write ChunkSize x ChunkSize block (coalesced on i for LayoutLeft)
            i = i_block + threadIdx.x;
            k = k_block + threadIdx.y;
            if (i < m && k < p)
                dst(i,j*p+k) = tmp[threadIdx.x][threadIdx.y];

          }
        }
      });
    }
  };

  // Specialization for Scalar2PCE specifically for Cuda that ensures
  // coalesced reads and writes
  template <>
  struct CopyScalar2PCE<Kokkos::Cuda> {
    typedef Kokkos::Cuda exec_space;

    template <typename DstView, typename SrcView>
    CopyScalar2PCE(const DstView& dst, const SrcView& src) {
      impl(dst,src);
    }

    template <typename DstView, typename SrcView>
    void impl(const DstView& dst, const SrcView& src) {
      typedef typename SrcView::non_const_value_type Scalar;
      typedef Kokkos::TeamPolicy<exec_space> Policy;
      typedef typename Policy::member_type Member;

      const size_t m = dst.extent(0);
      const size_t n = dst.extent(1);
      const size_t p = Kokkos::dimension_scalar(dst);

      const size_t ChunkSize = 16;
      const size_t M = (m+ChunkSize-1)/ChunkSize;

      Policy policy(M,ChunkSize,ChunkSize);
      Kokkos::parallel_for( policy, [=] __device__(const Member& team)
      {
        __shared__ Scalar tmp[ChunkSize][ChunkSize];
        const size_t i_block = blockIdx.x*ChunkSize;

        for (size_t j=0; j<n; ++j) {
          for (size_t k_block=0; k_block<p; k_block+=ChunkSize) {

            // Make sure previous iteration has completed before overwriting tmp
            __syncthreads();

            // Read ChunkSize x ChunkSize block (coalesced on i for LayoutLeft)
            size_t i = i_block + threadIdx.x;
            size_t k = k_block + threadIdx.y;
            if (i < m && k < p)
              tmp[threadIdx.x][threadIdx.y] = src(i,j*p+k);

            // Wait for all threads to finish
            __syncthreads();

            // Write ChunkSize x ChunkSize block (coalesced on k)
            i = i_block + threadIdx.y;
            k = k_block + threadIdx.x;
            if (i < m && k < p)
                dst(i,j).fastAccessCoeff(k) = tmp[threadIdx.y][threadIdx.x];

          }
        }
      });
    }
  };
#endif

  } // Impl

  template <typename DstView, typename SrcView>
  void copy_pce_to_scalar(const DstView& dst, const SrcView& src)
  {
    Impl::CopyPCE2Scalar<typename DstView::execution_space>(dst,src);
  }

  template <typename DstView, typename SrcView>
  void copy_scalar_to_pce(const DstView& dst, const SrcView& src)
  {
    Impl::CopyScalar2PCE<typename DstView::execution_space>(dst,src);
  }

  // Tpetra operator wrapper allowing a mean0-based operator (with double
  // scalar type) to be applied to a UQ::PCE multi-vector
  template <typename Scalar,
            typename LocalOrdinal,
            typename GlobalOrdinal,
            typename Node>
  class MeanBasedTpetraOperator :
    virtual public Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  public:
    typedef Scalar scalar_type;
    typedef LocalOrdinal local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node node_type;
    typedef typename scalar_type::value_type base_scalar_type;
    typedef Tpetra::Operator<base_scalar_type,LocalOrdinal,GlobalOrdinal,Node> scalar_op_type;

    MeanBasedTpetraOperator(const Teuchos::RCP<const scalar_op_type>& mb_op_) :
      mb_op(mb_op_) {}

    virtual ~MeanBasedTpetraOperator() {}

    virtual Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
    getDomainMap() const {
      return mb_op->getDomainMap();
    }

    virtual Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
    getRangeMap() const {
      return mb_op->getRangeMap();
    }

    virtual void
    apply (const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
           Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
           Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const
    {
      typedef typename scalar_mv_type::device_type device_type;

      auto xv = X.getLocalViewDevice(Tpetra::Access::ReadOnly);
      auto yv = Y.getLocalViewDevice(Tpetra::Access::ReadWrite);
      const size_t pce_size = Kokkos::dimension_scalar(xv);
      if (X_s == Teuchos::null ||
          X_s->getNumVectors() != X.getNumVectors()*pce_size)
        X_s = Teuchos::rcp(new scalar_mv_type(X.getMap(),
                                              X.getNumVectors()*pce_size));
      if (Y_s == Teuchos::null ||
          Y_s->getNumVectors() != Y.getNumVectors()*pce_size)
        Y_s = Teuchos::rcp(new scalar_mv_type(Y.getMap(),
                                              Y.getNumVectors()*pce_size));
      base_scalar_type alpha_s = alpha.fastAccessCoeff(0);
      base_scalar_type beta_s = beta.fastAccessCoeff(0);

      {
        auto xv_s = X_s->getLocalViewDevice(Tpetra::Access::ReadWrite);
        auto yv_s = Y_s->getLocalViewDevice(Tpetra::Access::ReadWrite);

        copy_pce_to_scalar(xv_s, xv);
        if (beta_s != 0.0)
          copy_pce_to_scalar(yv_s, yv);
      }

      mb_op->apply(*X_s, *Y_s, mode, alpha_s, beta_s);

      {
        auto yv_s = Y_s->getLocalViewDevice(Tpetra::Access::ReadOnly);
        copy_scalar_to_pce(yv, yv_s);
      }
    }

    virtual bool hasTransposeApply() const {
      return mb_op->hasTransposeApply();
    }

  private:

    typedef Tpetra::MultiVector<base_scalar_type,LocalOrdinal,GlobalOrdinal,Node> scalar_mv_type;
    mutable Teuchos::RCP<scalar_mv_type> X_s, Y_s;
    Teuchos::RCP<const scalar_op_type> mb_op;

  };

}

#endif // STOKHOS_TPETRA_UTILITIES_HPP
