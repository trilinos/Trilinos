/*
// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#ifndef DUMMY_TEST_MODEL_EVALUATOR_DEF_HPP
#define DUMMY_TEST_MODEL_EVALUATOR_DEF_HPP


#include "Thyra_DummyTestModelEvaluator_decl.hpp"
#include "Thyra_SimpleDenseLinearOp.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"


namespace Thyra {


// Nonmember constuctors


template<class Scalar>
Teuchos::RCP<DummyTestModelEvaluator<Scalar> >
dummyTestModelEvaluator()
{
  return Teuchos::rcp(new DummyTestModelEvaluator<Scalar>);
}


// Initializers/Accessors


template<class Scalar>
DummyTestModelEvaluator<Scalar>::DummyTestModelEvaluator()
  : x_space_(Thyra::defaultSpmdVectorSpace<Scalar>(2)),
    f_space_(x_space_),
    W_factory_(Thyra::defaultSerialDenseLinearOpWithSolveFactory<Scalar>()),
    p_(x_space_->dim(), Teuchos::ScalarTraits<Scalar>::zero())
{

  using Teuchos::RCP;
  using Thyra::VectorBase;
  using Thyra::createMember;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x);
  prototypeInArgs_ = inArgs;
  
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);
  outArgs.setSupports(MEB::OUT_ARG_W_prec);
  prototypeOutArgs_ = outArgs;

  nominalValues_ = inArgs;
  x0_ = createMember(x_space_);
  V_S(x0_.ptr(), ST::zero());
  nominalValues_.set_x(x0_);

  set_p(Teuchos::tuple<Scalar>(2.0, 0.0)());
  set_x0(Teuchos::tuple<Scalar>(1.0, 1.0)());

}


template<class Scalar>
void DummyTestModelEvaluator<Scalar>::set_p(const Teuchos::ArrayView<const Scalar> &p)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(p_.size(), p.size());
#endif
  p_().assign(p);
}


template<class Scalar>
void DummyTestModelEvaluator<Scalar>::set_x0(const Teuchos::ArrayView<const Scalar> &x0_in)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(x_space_->dim(), x0_in.size());
#endif
  Thyra::DetachedVectorView<Scalar> x0(x0_);
  x0.sv().values()().assign(x0_in);
}


// Public functions overridden from ModelEvaulator


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
DummyTestModelEvaluator<Scalar>::get_x_space() const
{
  return x_space_;
}


template<class Scalar>
Teuchos::RCP<const VectorSpaceBase<Scalar> >
DummyTestModelEvaluator<Scalar>::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> >
DummyTestModelEvaluator<Scalar>::get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
DummyTestModelEvaluator<Scalar>::get_f_space() const
{
  return f_space_;
}


template<class Scalar>
Teuchos::RCP<const VectorSpaceBase<Scalar> >
DummyTestModelEvaluator<Scalar>::get_g_space(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
DummyTestModelEvaluator<Scalar>::getNominalValues() const
{
  return nominalValues_;
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DummyTestModelEvaluator<Scalar>::getLowerBounds() const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
  return ModelEvaluatorBase::InArgs<Scalar>();
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DummyTestModelEvaluator<Scalar>::getUpperBounds() const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
  return ModelEvaluatorBase::InArgs<Scalar>();
}


template<class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
DummyTestModelEvaluator<Scalar>::create_W_op() const
{
  return createNonconstSimpleDenseLinearOp<Scalar>(
    createMembers<Scalar>(f_space_, x_space_->dim())
    );
}


template<class Scalar>
Teuchos::RCP<Thyra::PreconditionerBase<Scalar> >
DummyTestModelEvaluator<Scalar>::create_W_prec() const
{
  return nonconstUnspecifiedPrec<Scalar>(
    createNonconstSimpleDenseLinearOp<Scalar>(
      createMembers<Scalar>(f_space_, x_space_->dim())
      )
    );
}


template<class Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
DummyTestModelEvaluator<Scalar>::get_W_factory() const
{
  return W_factory_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
DummyTestModelEvaluator<Scalar>::createInArgs() const
{
  return prototypeInArgs_;
}


template<class Scalar>
void DummyTestModelEvaluator<Scalar>::reportFinalPoint(
  const ModelEvaluatorBase::InArgs<Scalar> &finalPoint,
  const bool wasSolved
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
DummyTestModelEvaluator<Scalar>::createOutArgsImpl() const
{
  return prototypeOutArgs_;
}


template<class Scalar>
void DummyTestModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement to just copy inArgs and outArgs!
  
}


} // namespace Thyra


//
// Explicit instantiation macro
//
// Must be expanded from within the global namespace!
//

#define DUMMY_TEST_MODEL_EVALUATOR_INSTANT(SCALAR) \
  \
  template class DummyTestModelEvaluator<SCALAR >; \
  \
  template Teuchos::RCP<DummyTestModelEvaluator<SCALAR > > \
  dummyTestModelEvaluator(); \


#endif // DUMMY_TEST_MODEL_EVALUATOR_DEF_HPP
