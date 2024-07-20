// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _TEUCHOS_BLAS_MP_VECTOR_HPP_
#define _TEUCHOS_BLAS_MP_VECTOR_HPP_

#include "Teuchos_BLAS.hpp"
#include "Sacado_MP_Vector.hpp"

namespace Teuchos
{

namespace details
{

template <typename Storage>
class GivensRotator<Sacado::MP::Vector<Storage>, false>
{
public:
  typedef Sacado::MP::Vector<Storage> ScalarType;
  typedef ScalarType c_type;

  void
  ROTG(ScalarType *da,
       ScalarType *db,
       ScalarType *c,
       ScalarType *s) const
  {
    typedef ScalarTraits<ScalarType> STS;

    ScalarType r, roe, scale, z, da_scaled, db_scaled;
    auto m_da = (STS::magnitude(*da) > STS::magnitude(*db));
    mask_assign(m_da, roe) = {*da, *db};

    scale = STS::magnitude(*da) + STS::magnitude(*db);

    auto m_scale = scale != STS::zero();

    da_scaled = *da;
    db_scaled = *db;

    *c = *da;
    *s = *db;

    ScalarType tmp = STS::one();
    mask_assign(m_scale, tmp) /= scale;

    mask_assign(m_scale, da_scaled) *= tmp;
    mask_assign(m_scale, db_scaled) *= tmp;

    r = scale * STS::squareroot(da_scaled * da_scaled + db_scaled * db_scaled);
    auto m_roe = roe < 0;
    mask_assign(m_roe, r) = -r;

    tmp = STS::one();
    mask_assign(m_scale, tmp) /= r;

    mask_assign(m_scale, *c) *= tmp;
    mask_assign(m_scale, *s) *= tmp;

    mask_assign(!m_scale, *c) = STS::one();
    mask_assign(!m_scale, *s) = STS::zero();

    mask_assign(*c != STS::zero(), z) /= {STS::one(), *c, STS::zero()};
    mask_assign(!m_scale, z) = STS::zero();
    mask_assign(m_da, z) = *s;

    *da = r;
    *db = z;
  }
};
} // namespace details
} // namespace Teuchos

//namespace Sacado {
//  namespace MP {
namespace Teuchos
{
//! Vector specializations for Teuchos::BLAS wrappers
template <typename OrdinalType, typename Storage>
class BLAS<OrdinalType, Sacado::MP::Vector<Storage>> : public Teuchos::DefaultBLASImpl<OrdinalType, Sacado::MP::Vector<Storage>>
{

  typedef typename Teuchos::ScalarTraits<Sacado::MP::Vector<Storage>>::magnitudeType MagnitudeType;
  typedef typename Sacado::ValueType<Sacado::MP::Vector<Storage>>::type ValueType;
  typedef typename Sacado::ScalarType<Sacado::MP::Vector<Storage>>::type scalar_type;
  typedef typename Sacado::MP::Vector<Storage> ScalarType;
  typedef Teuchos::DefaultBLASImpl<OrdinalType, Sacado::MP::Vector<Storage>> BLASType;

public:
  template <typename alpha_type, typename A_type>
  void TRSM(Teuchos::ESide side, Teuchos::EUplo uplo,
            Teuchos::ETransp transa, Teuchos::EDiag diag,
            const OrdinalType m, const OrdinalType n,
            const alpha_type &alpha,
            const A_type *A, const OrdinalType lda,
            ScalarType *B, const OrdinalType ldb) const
  {
    OrdinalType izero = OrdinalTraits<OrdinalType>::zero();
    OrdinalType ione = OrdinalTraits<OrdinalType>::one();
    alpha_type alpha_zero = ScalarTraits<alpha_type>::zero();
    A_type A_zero = ScalarTraits<A_type>::zero();
    ScalarType B_zero = ScalarTraits<ScalarType>::zero();
    alpha_type alpha_one = ScalarTraits<alpha_type>::one();
    ScalarType B_one = ScalarTraits<ScalarType>::one();
    ScalarType temp;
    OrdinalType NRowA = m;
    bool BadArgument = false;
    bool noUnit = (EDiagChar[diag] == 'N');
    bool noConj = (ETranspChar[transa] == 'T');

    if (!(ESideChar[side] == 'L'))
    {
      NRowA = n;
    }

    // Quick return.
    if (n == izero || m == izero)
    {
      return;
    }
    if (m < izero)
    {
      std::cout << "BLAS::TRSM Error: M == " << m << std::endl;
      BadArgument = true;
    }
    if (n < izero)
    {
      std::cout << "BLAS::TRSM Error: N == " << n << std::endl;
      BadArgument = true;
    }
    if (lda < NRowA)
    {
      std::cout << "BLAS::TRSM Error: LDA < " << NRowA << std::endl;
      BadArgument = true;
    }
    if (ldb < m)
    {
      std::cout << "BLAS::TRSM Error: LDB < MAX(1,M)" << std::endl;
      BadArgument = true;
    }

    if (!BadArgument)
    {
      int i, j, k;
      // Set the solution to the zero vector.
      auto alpha_is_zero = (alpha == alpha_zero);
      for (j = izero; j < n; j++)
      {
        for (i = izero; i < m; i++)
        {
          mask_assign(alpha_is_zero, B[j * ldb + i]) = B_zero;
        }
      }

      auto alpha_is_not_one = (alpha != alpha_one);

      { // Start the operations.
        if (ESideChar[side] == 'L')
        {
          //
          // Perform computations for OP(A)*X = alpha*B
          //
          if (ETranspChar[transa] == 'N')
          {
            //
            //  Compute B = alpha*inv( A )*B
            //
            if (EUploChar[uplo] == 'U')
            {
              // A is upper triangular.
              for (j = izero; j < n; j++)
              {
                // Perform alpha*B if alpha is not 1.
                for (i = izero; i < m; i++)
                {
                  mask_assign(alpha_is_not_one, B[j * ldb + i]) *= alpha;
                }

                // Perform a backsolve for column j of B.
                for (k = (m - ione); k > -ione; k--)
                {
                  // If this entry is zero, we don't have to do anything.
                  auto B_is_not_zero = (B[j * ldb + k] != B_zero);

                  if (noUnit)
                  {
                    mask_assign(B_is_not_zero, B[j * ldb + k]) /= A[k * lda + k];
                  }
                  for (i = izero; i < k; i++)
                  {
                    mask_assign(B_is_not_zero, B[j * ldb + i]) -= B[j * ldb + k] * A[k * lda + i];
                  }
                }
              }
            }
            else
            { // A is lower triangular.
              for (j = izero; j < n; j++)
              {
                // Perform alpha*B if alpha is not 1.
                for (i = izero; i < m; i++)
                {
                  mask_assign(alpha_is_not_one, B[j * ldb + i]) *= alpha;
                }
                // Perform a forward solve for column j of B.
                for (k = izero; k < m; k++)
                {
                  // If this entry is zero, we don't have to do anything.
                  auto B_is_not_zero = (B[j * ldb + k] != B_zero);
                  if (noUnit)
                  {
                    mask_assign(B_is_not_zero, B[j * ldb + k]) /= A[k * lda + k];
                  }
                  for (i = k + ione; i < m; i++)
                  {
                    mask_assign(B_is_not_zero, B[j * ldb + i]) -= B[j * ldb + k] * A[k * lda + i];
                  }
                }
              }
            } // end if (uplo == 'U')
          }   // if (transa =='N')
          else
          {
            //
            //  Compute B = alpha*inv( A' )*B
            //  or      B = alpha*inv( conj(A') )*B
            //
            if (EUploChar[uplo] == 'U')
            {
              // A is upper triangular.
              for (j = izero; j < n; j++)
              {
                for (i = izero; i < m; i++)
                {
                  temp = alpha * B[j * ldb + i];
                  if (noConj)
                  {
                    for (k = izero; k < i; k++)
                    {
                      temp -= A[i * lda + k] * B[j * ldb + k];
                    }
                    if (noUnit)
                    {
                      temp /= A[i * lda + i];
                    }
                  }
                  else
                  {
                    for (k = izero; k < i; k++)
                    {
                      temp -= ScalarTraits<A_type>::conjugate(A[i * lda + k]) * B[j * ldb + k];
                    }
                    if (noUnit)
                    {
                      temp /= ScalarTraits<A_type>::conjugate(A[i * lda + i]);
                    }
                  }
                  B[j * ldb + i] = temp;
                }
              }
            }
            else
            { // A is lower triangular.
              for (j = izero; j < n; j++)
              {
                for (i = (m - ione); i > -ione; i--)
                {
                  temp = alpha * B[j * ldb + i];
                  if (noConj)
                  {
                    for (k = i + ione; k < m; k++)
                    {
                      temp -= A[i * lda + k] * B[j * ldb + k];
                    }
                    if (noUnit)
                    {
                      temp /= A[i * lda + i];
                    }
                  }
                  else
                  {
                    for (k = i + ione; k < m; k++)
                    {
                      temp -= ScalarTraits<A_type>::conjugate(A[i * lda + k]) * B[j * ldb + k];
                    }
                    if (noUnit)
                    {
                      temp /= ScalarTraits<A_type>::conjugate(A[i * lda + i]);
                    }
                  }
                  B[j * ldb + i] = temp;
                }
              }
            }
          }
        } // if (side == 'L')
        else
        {
          // side == 'R'
          //
          // Perform computations for X*OP(A) = alpha*B
          //
          if (ETranspChar[transa] == 'N')
          {
            //
            //  Compute B = alpha*B*inv( A )
            //
            if (EUploChar[uplo] == 'U')
            {
              // A is upper triangular.
              // Perform a backsolve for column j of B.
              for (j = izero; j < n; j++)
              {
                // Perform alpha*B if alpha is not 1.
                for (i = izero; i < m; i++)
                {
                  mask_assign(alpha_is_not_one, B[j * ldb + i]) *= alpha;
                }
                for (k = izero; k < j; k++)
                {
                  // If this entry is zero, we don't have to do anything.
                  auto A_is_not_zero = (A[j * lda + k] != A_zero);
                  for (i = izero; i < m; i++)
                  {
                    mask_assign(A_is_not_zero, B[j * ldb + i]) -= A[j * lda + k] * B[k * ldb + i];
                  }
                }
                if (noUnit)
                {
                  temp = B_one / A[j * lda + j];
                  for (i = izero; i < m; i++)
                  {
                    B[j * ldb + i] *= temp;
                  }
                }
              }
            }
            else
            { // A is lower triangular.
              for (j = (n - ione); j > -ione; j--)
              {
                // Perform alpha*B if alpha is not 1.
                for (i = izero; i < m; i++)
                {
                  mask_assign(alpha_is_not_one, B[j * ldb + i]) *= alpha;
                }

                // Perform a forward solve for column j of B.
                for (k = j + ione; k < n; k++)
                {
                  // If this entry is zero, we don't have to do anything.
                  auto A_is_not_zero = (A[j * lda + k] != A_zero);
                  for (i = izero; i < m; i++)
                  {
                    mask_assign(A_is_not_zero, B[j * ldb + i]) -= A[j * lda + k] * B[k * ldb + i];
                  }
                }
                if (noUnit)
                {
                  temp = B_one / A[j * lda + j];
                  for (i = izero; i < m; i++)
                  {
                    B[j * ldb + i] *= temp;
                  }
                }
              }
            } // end if (uplo == 'U')
          }   // if (transa =='N')
          else
          {
            //
            //  Compute B = alpha*B*inv( A' )
            //  or      B = alpha*B*inv( conj(A') )
            //
            if (EUploChar[uplo] == 'U')
            {
              // A is upper triangular.
              for (k = (n - ione); k > -ione; k--)
              {
                if (noUnit)
                {
                  if (noConj)
                    temp = B_one / A[k * lda + k];
                  else
                    temp = B_one / ScalarTraits<A_type>::conjugate(A[k * lda + k]);
                  for (i = izero; i < m; i++)
                  {
                    B[k * ldb + i] *= temp;
                  }
                }
                for (j = izero; j < k; j++)
                {
                  auto A_is_not_zero = (A[k * lda + j] != A_zero);
                  if (noConj)
                    mask_assign(A_is_not_zero, temp) = A[k * lda + j];
                  else
                    mask_assign(A_is_not_zero, temp) = ScalarTraits<A_type>::conjugate(A[k * lda + j]);
                  for (i = izero; i < m; i++)
                  {
                    mask_assign(A_is_not_zero, B[j * ldb + i]) -= temp * B[k * ldb + i];
                  }
                }
                for (i = izero; i < m; i++)
                {
                  mask_assign(alpha_is_not_one, B[k * ldb + i]) *= alpha;
                }
              }
            }
            else
            { // A is lower triangular.
              for (k = izero; k < n; k++)
              {
                if (noUnit)
                {
                  if (noConj)
                    temp = B_one / A[k * lda + k];
                  else
                    temp = B_one / ScalarTraits<A_type>::conjugate(A[k * lda + k]);
                  for (i = izero; i < m; i++)
                  {
                    B[k * ldb + i] *= temp;
                  }
                }
                for (j = k + ione; j < n; j++)
                {
                  auto A_is_not_zero = (A[k * lda + j] != A_zero);
                  if (noConj)
                    mask_assign(A_is_not_zero, temp) = A[k * lda + j];
                  else
                    mask_assign(A_is_not_zero, temp) = ScalarTraits<A_type>::conjugate(A[k * lda + j]);
                  for (i = izero; i < m; i++)
                  {
                    mask_assign(A_is_not_zero, B[j * ldb + i]) -= temp * B[k * ldb + i];
                  }
                }
                for (i = izero; i < m; i++)
                {
                  mask_assign(alpha_is_not_one, B[k * ldb + i]) *= alpha;
                }
              }
            }
          }
        }
      }
    }
  }
}; // class BLAS
//  }
//}

} // namespace Teuchos

#endif // _TEUCHOS_BLAS__MP_VECTOR_HPP_
