#ifndef IFPACK2_DETAILS_SCALEDDAMPEDRESIDUAL_HPP
#define IFPACK2_DETAILS_SCALEDDAMPEDRESIDUAL_HPP

#include "Ifpack2_config.h"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_withLocalAccess_MultiVector.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Import.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"
#include <memory>
#include <type_traits>

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
  static_assert (static_cast<int> (WVector::Rank) == 1,
                 "WVector must be a rank 1 View.");
  static_assert (static_cast<int> (DVector::Rank) == 1,
                 "DVector must be a rank 1 View.");
  static_assert (static_cast<int> (BVector::Rank) == 1,
                 "BVector must be a rank 1 View.");
  static_assert (static_cast<int> (XVector::Rank) == 1,
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

template<class ExecutionSpace>
int64_t
scaled_damped_residual_vector_launch_parameters (int64_t numRows,
                                                 int64_t nnz,
                                                 int64_t rows_per_thread,
                                                 int& team_size,
                                                 int& vector_length)
{
  using execution_space = typename ExecutionSpace::execution_space;

  int64_t rows_per_team;
  int64_t nnz_per_row = nnz/numRows;

  if (nnz_per_row < 1) {
    nnz_per_row = 1;
  }

  if (vector_length < 1) {
    vector_length = 1;
    while (vector_length<32 && vector_length*6 < nnz_per_row) {
      vector_length *= 2;
    }
  }

  // Determine rows per thread
  if (rows_per_thread < 1) {
#ifdef KOKKOS_ENABLE_CUDA
    if (std::is_same<Kokkos::Cuda, execution_space>::value) {
      rows_per_thread = 1;
    }
    else
#endif
    {
      if (nnz_per_row < 20 && nnz > 5000000) {
        rows_per_thread = 256;
      }
      else {
        rows_per_thread = 64;
      }
    }
  }

#ifdef KOKKOS_ENABLE_CUDA
  if (team_size < 1) {
    if (std::is_same<Kokkos::Cuda,execution_space>::value) {
      team_size = 256/vector_length;
    }
    else {
      team_size = 1;
    }
  }
#endif

  rows_per_team = rows_per_thread * team_size;

  if (rows_per_team < 0) {
    int64_t nnz_per_team = 4096;
    int64_t conc = execution_space::concurrency ();
    while ((conc * nnz_per_team * 4 > nnz) &&
           (nnz_per_team > 256)) {
      nnz_per_team /= 2;
    }
    rows_per_team = (nnz_per_team + nnz_per_row - 1) / nnz_per_row;
  }

  return rows_per_team;
}

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

  const int64_t rows_per_team =
    scaled_damped_residual_vector_launch_parameters<execution_space>
      (A.numRows (), A.nnz (), rows_per_thread, team_size, vector_length);
  int64_t worksets = (b.extent (0) + rows_per_team - 1) / rows_per_team;

  using Kokkos::Dynamic;
  using Kokkos::Schedule;
  using Kokkos::TeamPolicy;
  using policy_type = TeamPolicy<execution_space, Schedule<Dynamic>>;
  const char kernel_label[] = "scaled_damped_residual_vector";
  policy_type policy (1, 1);
  if (team_size < 0) {
    policy = policy_type (worksets, Kokkos::AUTO, vector_length);
  }
  else {
    policy = policy_type (worksets, team_size, vector_length);
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
    Kokkos::parallel_for (kernel_label, policy, func);
  }
  else {
    constexpr bool use_beta = true;
    using functor_type =
      ScaledDampedResidualVectorFunctor<w_vec_type, d_vec_type,
                                        b_vec_type, matrix_type,
                                        x_vec_type, scalar_type,
                                        use_beta>;
    functor_type func (alpha, w, d, b, A, x, beta, rows_per_team);
    Kokkos::parallel_for (kernel_label, policy, func);
  }
}

} // namespace Impl

/// \brief Compute scaled damped residual for Chebyshev.
///
/// This is an implementation detail of Ifpack2::Chebyshev.  Given a
/// linear system A*X=B and an "inverse diagonal" matrix (stored as a
/// vector) D_inv, it computes a "scaled damped" residual vector W :=
/// alpha*D_inv*(B-A*X) + beta*W.
///
/// \tparam TpetraOperatorType Specialization of Tpetra::Operator.
///
/// \note To Ifpack2 developers: We can't fuse this with X := X + W,
///   because data dependencies in the input X are not elementwise
///   (unless A is diagonal).
template<class TpetraOperatorType>
class ScaledDampedResidual {
private:
  using SC = typename TpetraOperatorType::scalar_type;
  using LO = typename TpetraOperatorType::local_ordinal_type;
  using GO = typename TpetraOperatorType::global_ordinal_type;
  using NT = typename TpetraOperatorType::node_type;

  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using multivector_type = Tpetra::MultiVector<SC, LO, GO, NT>;
  using operator_type = Tpetra::Operator<SC, LO, GO, NT>;
  using vector_type = Tpetra::Vector<SC, LO, GO, NT>;

public:
  ScaledDampedResidual (const Teuchos::RCP<const operator_type>& A)
  {
    setMatrix (A);
  }

  void
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

  void
  compute (multivector_type& W,
           const SC& alpha,
           vector_type& D_inv,
           multivector_type& B,
           multivector_type& X,
           const SC& beta)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    const size_t numVecs = B.getNumVectors ();

    if (canFuse (B)) {
      // "nonconst" here has no effect other than on the return type.
      // We need to name these return types because Intel 17 isn't so
      // good at type deduction.
      RCP<vector_type> W_vec = W.getVectorNonConst (0);
      RCP<vector_type> B_vec = B.getVectorNonConst (0);
      RCP<vector_type> X_vec = X.getVectorNonConst (0);
      fusedCase (*W_vec, alpha, D_inv, *B_vec, *A_crs_, *X_vec, beta);
    }
    else {
      unfusedCase (W, alpha, D_inv, B, *A_op_, X, beta);
    }
  }

private:
  using import_type = Tpetra::Import<LO, GO, NT>;
  using export_type = Tpetra::Export<LO, GO, NT>;

  Teuchos::RCP<const operator_type> A_op_;
  Teuchos::RCP<const crs_matrix_type> A_crs_;
  Teuchos::RCP<const import_type> imp_;
  Teuchos::RCP<const export_type> exp_;
  std::unique_ptr<vector_type> X_colMap_;
  std::unique_ptr<multivector_type> V1_;

  // Do the Import, if needed, and return the column Map version of X.
  vector_type&
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

  bool canFuse (const multivector_type& B) const
  {
    return B.getNumVectors () == size_t (1) &&
      ! A_crs_.is_null () &&
      exp_.is_null ();
  }

  void
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

  void
  fusedCase (vector_type& W,
             const SC& alpha,
             vector_type& D_inv,
             vector_type& B,
             const crs_matrix_type& A,
             vector_type& X,
             const SC& beta)
  {
    vector_type& X_colMap = importVector (X);

    // Only need these aliases because we lack C++14 generic lambdas.
    using Tpetra::with_local_access_function_argument_type;
    using ro_lcl_vec_type =
      with_local_access_function_argument_type<
        decltype (readOnly (B))>;
    using wo_lcl_vec_type =
      with_local_access_function_argument_type<
        decltype (writeOnly (B))>;
    using rw_lcl_vec_type =
      with_local_access_function_argument_type<
        decltype (readWrite (B))>;

    using Tpetra::withLocalAccess;
    using Tpetra::readOnly;
    using Tpetra::readWrite;
    using Tpetra::writeOnly;
    using Impl::scaled_damped_residual_vector;
    using STS = Teuchos::ScalarTraits<SC>;

    auto A_lcl = A.getLocalMatrix ();
    if (beta == STS::zero ()) {
      withLocalAccess
        ([&] (const wo_lcl_vec_type& W_lcl,
              const ro_lcl_vec_type& D_lcl,
              const ro_lcl_vec_type& B_lcl,
              const ro_lcl_vec_type& X_lcl) {
           scaled_damped_residual_vector (alpha, W_lcl, D_lcl,
                                          B_lcl, A_lcl, X_lcl, beta);
         },
         writeOnly (W),
         readOnly (D_inv),
         readOnly (B),
         readOnly (X_colMap));
    }
    else { // need to read _and_ write W if beta != 0
      withLocalAccess
        ([&] (const rw_lcl_vec_type& W_lcl,
              const ro_lcl_vec_type& D_lcl,
              const ro_lcl_vec_type& B_lcl,
              const ro_lcl_vec_type& X_lcl) {
           scaled_damped_residual_vector (alpha, W_lcl, D_lcl,
                                          B_lcl, A_lcl, X_lcl, beta);
         },
         readWrite (W),
         readOnly (D_inv),
         readOnly (B),
         readOnly (X_colMap));
    }
  }
};

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_SCALEDDAMPEDRESIDUAL_HPP
