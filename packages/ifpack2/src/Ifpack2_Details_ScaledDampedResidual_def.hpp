// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DETAILS_SCALEDDAMPEDRESIDUAL_DEF_HPP
#define IFPACK2_DETAILS_SCALEDDAMPEDRESIDUAL_DEF_HPP

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Export_decl.hpp"
#include "Tpetra_Import_decl.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Teuchos_Assert.hpp"
#include <type_traits>
#include "KokkosSparse_spmv_impl.hpp"

namespace Ifpack2 {
namespace Details {
namespace Impl {

/// \brief Functor for computing W := alpha * D * (B - A*X) + beta * W.
///
/// This is an implementation detail of scaled_damped_residual_vector,
/// which in turn is an implementation detail of ScaledDampedResidual.
template<class WVector,
         class DVector,
         class BVector,
         class AMatrix,
         class XVector,
         class Scalar,
         bool use_beta>
struct ScaledDampedResidualVectorFunctor {
  static_assert (static_cast<int> (WVector::rank) == 1,
                 "WVector must be a rank 1 View.");
  static_assert (static_cast<int> (DVector::rank) == 1,
                 "DVector must be a rank 1 View.");
  static_assert (static_cast<int> (BVector::rank) == 1,
                 "BVector must be a rank 1 View.");
  static_assert (static_cast<int> (XVector::rank) == 1,
                 "XVector must be a rank 1 View.");

  using execution_space = typename AMatrix::execution_space;
  using LO = typename AMatrix::non_const_ordinal_type;
  using value_type = typename AMatrix::non_const_value_type;
  using team_policy = typename Kokkos::TeamPolicy<execution_space>;
  using team_member = typename team_policy::member_type;
  using ATV = Kokkos::ArithTraits<value_type>;

  const Scalar alpha;
  WVector m_w;
  DVector m_d;
  BVector m_b;
  AMatrix m_A;
  XVector m_x;
  const Scalar beta;

  const LO rows_per_team;

  ScaledDampedResidualVectorFunctor (const Scalar& alpha_,
                                     const WVector& m_w_,
                                     const DVector& m_d_,
                                     const BVector& m_b_,
                                     const AMatrix& m_A_,
                                     const XVector& m_x_,
                                     const Scalar& beta_,
                                     const int rows_per_team_) :
    alpha (alpha_),
    m_w (m_w_),
    m_d (m_d_),
    m_b (m_b_),
    m_A (m_A_),
    m_x (m_x_),
    beta (beta_),
    rows_per_team (rows_per_team_)
  {
    const size_t numRows = m_A.numRows ();
    const size_t numCols = m_A.numCols ();

    TEUCHOS_ASSERT( m_w.extent (0) == m_d.extent (0) );
    TEUCHOS_ASSERT( m_w.extent (0) == m_b.extent (0) );
    TEUCHOS_ASSERT( numRows == size_t (m_w.extent (0)) );
    TEUCHOS_ASSERT( numCols <= size_t (m_x.extent (0)) );
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const team_member& dev) const
  {
    using residual_value_type = typename BVector::non_const_value_type;
    using KAT = Kokkos::ArithTraits<residual_value_type>;

    Kokkos::parallel_for
      (Kokkos::TeamThreadRange (dev, 0, rows_per_team),
       [&] (const LO& loop) {
         const LO lclRow =
           static_cast<LO> (dev.league_rank ()) * rows_per_team + loop;
         if (lclRow >= m_A.numRows ()) {
           return;
         }
         const auto A_row = m_A.rowConst(lclRow);
         const LO row_length = static_cast<LO> (A_row.length);
         residual_value_type A_x = KAT::zero ();

         Kokkos::parallel_reduce
           (Kokkos::ThreadVectorRange (dev, row_length),
            [&] (const LO iEntry, residual_value_type& lsum) {
              const auto A_val = A_row.value(iEntry);
              lsum += A_val * m_x(A_row.colidx(iEntry));
            }, A_x);

         Kokkos::single
           (Kokkos::PerThread(dev),
            [&] () {
              const auto alpha_D_res =
                alpha * m_d(lclRow) * (m_b(lclRow) - A_x);
              if (use_beta) {
                m_w(lclRow) = beta * m_w(lclRow) + alpha_D_res;
              }
              else {
                m_w(lclRow) = alpha_D_res;
              }
            });
       });
  }

};


// W := alpha * D * (B - A*X) + beta * W.
template<class WVector,
         class DVector,
         class BVector,
         class AMatrix,
         class XVector,
         class Scalar>
static void
scaled_damped_residual_vector
(const Scalar& alpha,
 const WVector& w,
 const DVector& d,
 const BVector& b,
 const AMatrix& A,
 const XVector& x,
 const Scalar& beta)
{
  using execution_space = typename AMatrix::execution_space;

  if (A.numRows () == 0) {
    return;
  }

  int team_size = -1;
  int vector_length = -1;
  int64_t rows_per_thread = -1;

  const int64_t rows_per_team = KokkosSparse::Impl::spmv_launch_parameters<execution_space>(A.numRows(), A.nnz(), rows_per_thread, team_size, vector_length);
  int64_t worksets = (b.extent (0) + rows_per_team - 1) / rows_per_team;

  using Kokkos::Dynamic;
  using Kokkos::Static;
  using Kokkos::Schedule;
  using Kokkos::TeamPolicy;
  using policy_type_dynamic = TeamPolicy<execution_space, Schedule<Dynamic> >;
  using policy_type_static = TeamPolicy<execution_space, Schedule<Static> >;
  const char kernel_label[] = "scaled_damped_residual_vector";
  policy_type_dynamic policyDynamic (1, 1);
  policy_type_static  policyStatic (1, 1);
  if (team_size < 0) {
    policyDynamic = policy_type_dynamic (worksets, Kokkos::AUTO, vector_length);
    policyStatic  = policy_type_static  (worksets, Kokkos::AUTO, vector_length);
  }
  else {
    policyDynamic = policy_type_dynamic (worksets, team_size, vector_length);
    policyStatic  = policy_type_static  (worksets, team_size, vector_length);
  }

  // Canonicalize template arguments to avoid redundant instantiations.
  using w_vec_type = typename WVector::non_const_type;
  using d_vec_type = typename DVector::const_type;
  using b_vec_type = typename BVector::const_type;
  using matrix_type = AMatrix;
  using x_vec_type = typename XVector::const_type;
  using scalar_type = typename Kokkos::ArithTraits<Scalar>::val_type;

  if (beta == Kokkos::ArithTraits<Scalar>::zero ()) {
    constexpr bool use_beta = false;
    using functor_type =
      ScaledDampedResidualVectorFunctor<w_vec_type, d_vec_type,
                                        b_vec_type, matrix_type,
                                        x_vec_type, scalar_type,
                                        use_beta>;
    functor_type func (alpha, w, d, b, A, x, beta, rows_per_team);
    if(A.nnz()>10000000)
      Kokkos::parallel_for (kernel_label, policyDynamic, func);
    else
      Kokkos::parallel_for (kernel_label, policyStatic, func);
  }
  else {
    constexpr bool use_beta = true;
    using functor_type =
      ScaledDampedResidualVectorFunctor<w_vec_type, d_vec_type,
                                        b_vec_type, matrix_type,
                                        x_vec_type, scalar_type,
                                        use_beta>;
    functor_type func (alpha, w, d, b, A, x, beta, rows_per_team);
    if(A.nnz()>10000000)
      Kokkos::parallel_for (kernel_label, policyDynamic, func);
    else
      Kokkos::parallel_for (kernel_label, policyStatic, func);
  }
}

} // namespace Impl

template<class TpetraOperatorType>
ScaledDampedResidual<TpetraOperatorType>::
ScaledDampedResidual (const Teuchos::RCP<const operator_type>& A)
{
  setMatrix (A);
}

template<class TpetraOperatorType>
void
ScaledDampedResidual<TpetraOperatorType>::
setMatrix (const Teuchos::RCP<const operator_type>& A)
{
  if (A_op_.get () != A.get ()) {
    A_op_ = A;

    // We'll (re)allocate these on demand.
    X_colMap_ = std::unique_ptr<vector_type> (nullptr);
    V1_ = std::unique_ptr<multivector_type> (nullptr);

    using Teuchos::rcp_dynamic_cast;
    Teuchos::RCP<const crs_matrix_type> A_crs =
      rcp_dynamic_cast<const crs_matrix_type> (A);
    if (A_crs.is_null ()) {
      A_crs_ = Teuchos::null;
      imp_ = Teuchos::null;
      exp_ = Teuchos::null;
    }
    else {
      TEUCHOS_ASSERT( A_crs->isFillComplete () );
      A_crs_ = A_crs;
      auto G = A_crs->getCrsGraph ();
      imp_ = G->getImporter ();
      exp_ = G->getExporter ();
    }
  }
}

template<class TpetraOperatorType>
void
ScaledDampedResidual<TpetraOperatorType>::
compute (multivector_type& W,
         const SC& alpha,
         vector_type& D_inv,
         multivector_type& B,
         multivector_type& X,
         const SC& beta)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  if (canFuse (B)) {
    // "nonconst" here has no effect other than on the return type.
    W_vec_ = W.getVectorNonConst (0);
    B_vec_ = B.getVectorNonConst (0);
    X_vec_ = X.getVectorNonConst (0);
    TEUCHOS_ASSERT( ! A_crs_.is_null () );
    fusedCase (*W_vec_, alpha, D_inv, *B_vec_, *A_crs_, *X_vec_, beta);
  }
  else {
    TEUCHOS_ASSERT( ! A_op_.is_null () );
    unfusedCase (W, alpha, D_inv, B, *A_op_, X, beta);
  }
}

template<class TpetraOperatorType>
typename ScaledDampedResidual<TpetraOperatorType>::vector_type&
ScaledDampedResidual<TpetraOperatorType>::
importVector (vector_type& X_domMap)
{
  if (imp_.is_null ()) {
    return X_domMap;
  }
  else {
    if (X_colMap_.get () == nullptr) {
      using V = vector_type;
      X_colMap_ = std::unique_ptr<V> (new V (imp_->getTargetMap ()));
    }
    X_colMap_->doImport (X_domMap, *imp_, Tpetra::REPLACE);
    return *X_colMap_;
  }
}

template<class TpetraOperatorType>
bool
ScaledDampedResidual<TpetraOperatorType>::
canFuse (const multivector_type& B) const
{
  return B.getNumVectors () == size_t (1) &&
    ! A_crs_.is_null () &&
    exp_.is_null ();
}

template<class TpetraOperatorType>
void
ScaledDampedResidual<TpetraOperatorType>::
unfusedCase (multivector_type& W,
             const SC& alpha,
             vector_type& D_inv,
             multivector_type& B,
             const operator_type& A,
             multivector_type& X,
             const SC& beta)
{
  if (V1_.get () == nullptr) {
    using MV = multivector_type;
    const size_t numVecs = B.getNumVectors ();
    V1_ = std::unique_ptr<MV> (new MV (B.getMap (), numVecs));
  }
  const SC one = Teuchos::ScalarTraits<SC>::one ();

  // V1 = B - A*X
  Tpetra::deep_copy (*V1_, B);
  A.apply (X, *V1_, Teuchos::NO_TRANS, -one, one);

  // W := alpha * D_inv * V1 + beta * W
  W.elementWiseMultiply (alpha, D_inv, *V1_, beta);
}

template<class TpetraOperatorType>
void
ScaledDampedResidual<TpetraOperatorType>::
fusedCase (vector_type& W,
           const SC& alpha,
           vector_type& D_inv,
           vector_type& B,
           const crs_matrix_type& A,
           vector_type& X,
           const SC& beta)
{
  vector_type& X_colMap = importVector (X);

  using Impl::scaled_damped_residual_vector;
  using STS = Teuchos::ScalarTraits<SC>;

  auto A_lcl = A.getLocalMatrixDevice ();
  auto Dinv_lcl = Kokkos::subview(D_inv.getLocalViewDevice(Tpetra::Access::ReadOnly), Kokkos::ALL(), 0);
  auto B_lcl = Kokkos::subview(B.getLocalViewDevice(Tpetra::Access::ReadOnly), Kokkos::ALL(), 0);
  auto X_lcl = Kokkos::subview(X_colMap.getLocalViewDevice(Tpetra::Access::ReadOnly), Kokkos::ALL(), 0);
  if (beta == STS::zero ()) {
    auto W_lcl = Kokkos::subview(W.getLocalViewDevice(Tpetra::Access::OverwriteAll), Kokkos::ALL(), 0);
    scaled_damped_residual_vector (alpha, W_lcl, Dinv_lcl,
        B_lcl, A_lcl, X_lcl, beta);
  }
  else { // need to read _and_ write W if beta != 0
    auto W_lcl = Kokkos::subview(W.getLocalViewDevice(Tpetra::Access::ReadWrite), Kokkos::ALL(), 0);
    scaled_damped_residual_vector (alpha, W_lcl, Dinv_lcl,
        B_lcl, A_lcl, X_lcl, beta);
  }
}

} // namespace Details
} // namespace Ifpack2

#define IFPACK2_DETAILS_SCALEDDAMPEDRESIDUAL_INSTANT(SC,LO,GO,NT) \
  template class Ifpack2::Details::ScaledDampedResidual<Tpetra::Operator<SC, LO, GO, NT> >;

#endif // IFPACK2_DETAILS_SCALEDDAMPEDRESIDUAL_DEF_HPP
