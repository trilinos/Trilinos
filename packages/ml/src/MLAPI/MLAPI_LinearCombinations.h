#ifndef ML_LINEARCOMBINATION_H
#define ML_LINEARCOMBINATION_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_LinearCombination.h

\brief Concrete implementations of BaseLinearCombination.

\author Marzio Sala, D-INFK/ETHZ.

\date Last updated on Mar-06.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "MLAPI_BaseLinearCombination.h"

namespace MLAPI {

class BaseOperator;
class MultiVector;

// ============================================================================
class LinearCombinationAdd : public BaseLinearCombination
{
public:
  LinearCombinationAdd(const BaseLinearCombination& left,
                       const BaseLinearCombination& right) :
    left_(left),
    right_(right)
  {}

  inline const Space GetVectorSpace() const
  {
    return(left_.GetVectorSpace());
  }

  inline void Update(MultiVector& v) const
  {
    left_.Update(v);
    right_.Update(v);
  }

  inline void Set(MultiVector& v) const
  {
    left_.Set(v);
    right_.Update(v);
  }

private:
  const BaseLinearCombination& left_;
  const BaseLinearCombination& right_;
};

// ============================================================================
class LinearCombinationMixed : public BaseLinearCombination
{
public:

  LinearCombinationMixed(const BaseLinearCombination& left,
                         const MultiVector& right, double alpha) :
    left_(left),
    right_(right),
    alpha_(alpha)
  {}

  const Space GetVectorSpace() const;

  void Update(MultiVector& v) const;

  void Set(MultiVector& v) const;

private:
  const BaseLinearCombination& left_;
  const MultiVector            right_;
  double                       alpha_;
};

// ============================================================================
class LinearCombinationScaled : public BaseLinearCombination
{
public:
  LinearCombinationScaled(const BaseLinearCombination& left, double scalar) :
    left_(left),
    scalar_(scalar)
  {}

  const Space GetVectorSpace() const;

  void Set(MultiVector& v) const;

  void Update(MultiVector& v) const;

private:
  const BaseLinearCombination& left_;
  double                       scalar_;
};

// ============================================================================
// scaled vector, ScaledMultiVector = alpha * MultiVector
// ============================================================================
class MultiVectorScaled : public BaseLinearCombination
{
public:
  MultiVectorScaled(const MultiVector& vector, const double alpha) :
    vector_(vector),
    alpha_(alpha)
  {}

  const Space GetVectorSpace() const;

  const MultiVector& GetMultiVector() const
  {
    return(vector_);
  }

  double GetScalar() const
  {
    return(alpha_);
  }

  void Update(MultiVector& v) const;

  void Set(MultiVector& v) const;

private:
  const MultiVector vector_;
  double            alpha_;
};

// ============================================================================
class MultiVectorCombination : public BaseLinearCombination
{
public:
  MultiVectorCombination(const double alpha,
                         const MultiVector & x,
                         const double beta,
                         const MultiVector & y) :
    x_(x),
    y_(y),
    alpha_(alpha),
    beta_(beta)
  {}

  const Space GetVectorSpace() const;

  const MultiVector GetLeftMultiVector() const
  {
    return(x_);
  }

  inline double GetLeftScalar() const
  {
    return(alpha_);
  }

  inline const MultiVector GetRightMultiVector() const
  {
    return(y_);
  }

  inline double GetRightScalar() const
  {
    return(beta_);
  }

  void Update(MultiVector& v) const;

  void Set(MultiVector& v) const;

private:
  const MultiVector x_;
  const MultiVector y_;
  double alpha_, beta_;
};

// ============================================================================
// v = A * x
// ============================================================================
class BaseOperatorTimesMultiVector : public BaseLinearCombination
{
public:
  BaseOperatorTimesMultiVector(const BaseOperator& A,
                               const MultiVector& x) :
    A_(A),
    x_(x)
  {}

  const Space GetVectorSpace() const;

  inline const BaseOperator& GetBaseOperator() const
  {
    return(A_);
  }

  inline const MultiVector& GetMultiVector() const
  {
    return(x_);
  }

  void Update(MultiVector& v) const;

  void Set(MultiVector& v) const;

private:
  const BaseOperator& A_;
  const MultiVector   x_;
};

// ============================================================================
// v += alpha * b + beta * A * x
// ============================================================================
class Residual : public BaseLinearCombination
{
public:

  Residual(double alpha, const MultiVector& b, double beta,
           const BaseOperator& A, const MultiVector& x) :
    A_(A),
    b_(b),
    x_(x),
    alpha_(alpha),
    beta_(beta)
  {}

  const Space GetVectorSpace() const;

  void Update(MultiVector& v) const;

  void Set(MultiVector& v) const;

private:

  const BaseOperator& A_;
  const MultiVector  b_;
  const MultiVector  x_;
  double             alpha_;
  double             beta_;
};

} // namespace MLAPI

#endif
