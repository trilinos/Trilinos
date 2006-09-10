#include "ml_common.h"
#ifdef HAVE_ML_MLAPI
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

#include "MLAPI_MultiVector.h"
#include "MLAPI_MultiVector_Utils.h"
#include "MLAPI_BaseOperator.h"
#include "MLAPI_BaseLinearCombination.h"
#include "MLAPI_LinearCombinations.h"

namespace MLAPI {

const Space LinearCombinationMixed::GetVectorSpace() const
{
  return(right_.GetVectorSpace());
}

void LinearCombinationMixed::Update(MultiVector& v) const
{
  StackPush();

  // FIXME THIS IS A DUPLICATE
  MultiVector tmp(v.GetVectorSpace());
  left_.Update(tmp);
  v.Update(1.0, tmp, 1.0);
  if (v.IsAlias(right_))
    v.Scale(1.0 + alpha_);
  else
    v.Update(1.0, right_, alpha_);

  StackPop();
}

void LinearCombinationMixed::Set(MultiVector& v) const
{
  StackPush();

  if (v.GetVectorSpace() != left_.GetVectorSpace() ||
      v.IsAlias(right_))
    v.Reshape(left_.GetVectorSpace());

  left_.Set(v);
  v.Update(1.0, right_, alpha_);

  StackPop();
}

const Space MultiVectorScaled::GetVectorSpace() const
{
  return(vector_.GetVectorSpace());
}

const Space LinearCombinationScaled::GetVectorSpace() const
{
  return(left_.GetVectorSpace());
}

void LinearCombinationScaled::Update(MultiVector& v) const
{
  StackPush();

  MultiVector tmp(v.GetVectorSpace());
  left_.Update(tmp);
  v.Update(scalar_, tmp, 1.0);

  StackPop();
}

void LinearCombinationScaled::Set(MultiVector& v) const
{
  StackPush();

  left_.Set(v);
  v.Scale(scalar_);

  StackPop();
}

void MultiVectorScaled::Update(MultiVector& v) const
{
  StackPush();

  // FIXME THIS IS A DUPLICATE
  if (v.GetVectorSpace() != vector_.GetVectorSpace() ||
      v.IsAlias(vector_))
    v.Reshape(vector_.GetVectorSpace());

  v.Update(alpha_, vector_, 1.0);

  StackPop();
}

void MultiVectorScaled::Set(MultiVector& v) const
{
  StackPush();

  if (v.GetVectorSpace() != vector_.GetVectorSpace() ||
      v.IsAlias(vector_))
    v.Reshape(vector_.GetVectorSpace());

  v.Update(alpha_, vector_);

  StackPop();
}

const Space MultiVectorCombination::GetVectorSpace() const
{
  return(x_.GetVectorSpace());
}

void MultiVectorCombination::Update(MultiVector& v) const
{
  StackPush();

  if (v.GetVectorSpace() != x_.GetVectorSpace() ||
      v.IsAlias(x_) || v.IsAlias(y_))
    v = Duplicate(v);

  v.Update(alpha_, x_, 1.0);
  v.Update(beta_, y_, 1.0);

  StackPop();
}

void MultiVectorCombination::Set(MultiVector& v) const
{
  StackPush();

  if (v.GetVectorSpace() != x_.GetVectorSpace() ||
      v.IsAlias(x_) || v.IsAlias(y_))
    v.Reshape(x_.GetVectorSpace());

  v.Update(alpha_, x_, beta_, y_);

  StackPop();
}

const Space BaseOperatorTimesMultiVector::GetVectorSpace() const
{
  return(x_.GetVectorSpace());
}

void BaseOperatorTimesMultiVector::Update(MultiVector& v) const
{
  StackPush();

  // FIXME THIS IS DUPLICATE
  if (v.GetVectorSpace() != x_.GetVectorSpace() ||
      v.IsAlias(x_))
    v.Reshape(x_.GetVectorSpace());
  
  MultiVector tmp(v.GetVectorSpace());
  A_.Apply(x_, tmp);
  v.Update(1.0, tmp, 1.0);

  StackPop();
}

void BaseOperatorTimesMultiVector::Set(MultiVector& v) const
{
  StackPush();

  if (v.GetVectorSpace() != A_.GetOperatorRangeSpace() ||
      v.IsAlias(x_))
  {
    v.Reshape(A_.GetOperatorRangeSpace(), x_.GetNumVectors());
  }

  v = 0.0;

  A_.Apply(x_, v);

  StackPop();
}

const Space Residual::GetVectorSpace() const
{
  return(b_.GetVectorSpace());
}

void Residual::Update(MultiVector& v) const
{
  StackPush();
  
  // FIXME: THIS IS A DUPLICATE
  if (v.IsAlias(b_) || v.IsAlias(x_))
  {
    v.Reshape(b_.GetVectorSpace(), b_.GetNumVectors());
  }

  MultiVector tmp(v.GetVectorSpace());
  A_.Apply(x_, tmp);
  tmp.Update(alpha_, b_, beta_);
  v.Update(1.0, tmp, 1.0);

  StackPop();
}

void Residual::Set(MultiVector& v) const
{
  StackPush();

  if (v.GetVectorSpace() != A_.GetOperatorRangeSpace() ||
      v.IsAlias(b_) || v.IsAlias(x_))
  {
    v.Reshape(b_.GetVectorSpace(), b_.GetNumVectors());
  }

  A_.Apply(x_, v);
  v.Update(alpha_, b_, beta_);

  StackPop();
}

} // namespace MLAPI
#endif // defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#endif // ifdef HAVE_ML_MLAPI
