/*
// @HEADER
// ***********************************************************************
// 
//    OptiPack: Collection of simple Thyra-based Optimization ANAs
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#include "Diagonal_ThyraROME.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_Assert.hpp"


namespace TriKota {


//
// DiagonalROME
//


// Constructors, Initialization, Misc.
using Teuchos::RCP;
using Teuchos::rcp;

template<class Scalar>
DiagonalROME<Scalar>::DiagonalROME(
  const int localDim,
  const RCP<const Teuchos::Comm<Thyra::Ordinal> > &comm
  )
  :Np_(1), Ng_(1), comm_(comm), localDim_(localDim),
   nonlinearTermFactor_(0.0), g_offset_(0.0)
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Thyra::createMember;

  TEUCHOS_ASSERT( localDim > 0 );

  // Get the comm
  if (is_null(comm_)) {
    comm_ = Teuchos::DefaultComm<Thyra::Ordinal>::getComm();
  }

  // Locally replicated space for g
  g_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(comm_, 1, 1);

  // Distributed space for p
  p_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(comm_, localDim, -1);

  // Default solution
  const RCP<Thyra::VectorBase<Scalar> > ps = createMember<Scalar>(p_space_);
  V_S(ps.ptr(), ST::zero());
  ps_ = ps;

  // Default diagonal
  const RCP<Thyra::VectorBase<Scalar> > diag = createMember<Scalar>(p_space_);
  V_S(diag.ptr(), ST::one());
  diag_ = diag;
  diag_bar_ = diag;

  // Default inner product
  const RCP<Thyra::VectorBase<Scalar> > s_bar = createMember<Scalar>(p_space_);
  V_S(s_bar.ptr(), ST::one());
  s_bar_ = s_bar;

  // Default response offset
  g_offset_ = ST::zero();

}


template<class Scalar>
void DiagonalROME<Scalar>::setSolutionVector(
  const RCP<const Thyra::VectorBase<Scalar> > &ps)
{
  ps_ = ps.assert_not_null();
}


template<class Scalar>
const RCP<const Thyra::VectorBase<Scalar> >
DiagonalROME<Scalar>::getSolutionVector() const
{
  return ps_;
}


template<class Scalar>
void DiagonalROME<Scalar>::setDiagonalVector(
  const RCP<const Thyra::VectorBase<Scalar> > &diag)
{
  diag_ = diag;
  diag_bar_ = diag;
}


template<class Scalar>
void DiagonalROME<Scalar>::setDiagonalBarVector(
  const RCP<const Thyra::VectorBase<Scalar> > &diag_bar)
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::createMember;
  using Thyra::ele_wise_divide;
  using Thyra::V_S;

  diag_bar_ = diag_bar.assert_not_null();

  // Reset the scalar product for p_space!

  RCP<Thyra::VectorBase<Scalar> > s_bar = createMember<Scalar>(p_space_);

  // s_bar[i] = diag[i] / diag_bar[i] 
  V_S( s_bar.ptr(), ST::zero() );
  ele_wise_divide( ST::one(), *diag_, *diag_bar_, s_bar.ptr() );
  s_bar_ = s_bar;
  
  const RCP<Thyra::ScalarProdVectorSpaceBase<Scalar> > sp_p_space =
    rcp_dynamic_cast<Thyra::ScalarProdVectorSpaceBase<Scalar> >(p_space_, true);
  //sp_p_space->setScalarProd(diagonalScalarProd<Scalar>(s_bar_));

}


template<class Scalar>
void DiagonalROME<Scalar>::setNonlinearTermFactor(
  const Scalar &nonlinearTermFactor)
{
  nonlinearTermFactor_ = nonlinearTermFactor;
}


template<class Scalar>
void DiagonalROME<Scalar>::setScalarOffset(
  const Scalar &g_offset)
{
  g_offset_ = g_offset;
}


// Public functions overridden from ModelEvaulator


template<class Scalar>
int DiagonalROME<Scalar>::Np() const
{
  return Np_;
}


template<class Scalar>
int DiagonalROME<Scalar>::Ng() const
{
  return Ng_;
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
DiagonalROME<Scalar>::get_p_space(int l) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, Np_ );
#endif
  return p_space_;
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
DiagonalROME<Scalar>::get_g_space(int j) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( j, 0, Ng_ );
#endif
  return g_space_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
DiagonalROME<Scalar>::getNominalValues() const
{
  Thyra::ModelEvaluatorBase::InArgs<Scalar> initialGuess =
    this->createInArgs();
  RCP<Thyra::VectorBase<Scalar> > p_init =
    Thyra::createMember<Scalar>(p_space_);
  Thyra::V_S( p_init.ptr(), 1.5 );
  initialGuess.set_p(0, p_init);
  return initialGuess;
}

template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
DiagonalROME<Scalar>::createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(Np_);
  return inArgs;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
DiagonalROME<Scalar>::createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef MEB::DerivativeSupport DS;
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(Np_,Ng_);
  outArgs.setSupports(MEB::OUT_ARG_DgDp, 0 ,0, MEB::DERIV_TRANS_MV_BY_ROW);
  return outArgs;
}


template<class Scalar>
void DiagonalROME<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs
  ) const
{

  using Teuchos::as;
  using Teuchos::outArg;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Thyra::get_mv;
  using Thyra::ConstDetachedSpmdVectorView;
  using Thyra::DetachedSpmdVectorView;
  typedef Thyra::Ordinal Ordinal;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef MEB::DerivativeMultiVector<Scalar> DMV;

  const ConstDetachedSpmdVectorView<Scalar> p(inArgs.get_p(0));
  const ConstDetachedSpmdVectorView<Scalar> ps(ps_);
  const ConstDetachedSpmdVectorView<Scalar> diag(diag_);
  const ConstDetachedSpmdVectorView<Scalar> s_bar(s_bar_);

  // g(p)
  if (!is_null(outArgs.get_g(0))) {
    Scalar g_val = ST::zero();
    for (Ordinal i = 0; i < p.subDim(); ++i) {
      const Scalar p_ps = p[i] - ps[i];
      g_val += diag[i] * p_ps*p_ps;
      if (nonlinearTermFactor_ != ST::zero()) {
        g_val += nonlinearTermFactor_ * p_ps * p_ps * p_ps;
      }
    }
    Scalar global_g_val;
    Teuchos::reduceAll<Ordinal, Scalar>(*comm_, Teuchos::REDUCE_SUM, 
      g_val, outArg(global_g_val) );
    DetachedSpmdVectorView<Scalar>(outArgs.get_g(0))[0] =
      as<Scalar>(0.5) * global_g_val + g_offset_;
  }

  // DgDp[i]
  if (!outArgs.get_DgDp(0,0).isEmpty()) {
    const RCP<Thyra::MultiVectorBase<Scalar> > DgDp_trans_mv =
      get_mv<Scalar>(outArgs.get_DgDp(0,0), "DgDp^T", MEB::DERIV_TRANS_MV_BY_ROW);
    const DetachedSpmdVectorView<Scalar> DgDp_grad(DgDp_trans_mv->col(0));
    for (Thyra::Ordinal i = 0; i < p.subDim(); ++i) {
      const Scalar p_ps = p[i] - ps[i];
      Scalar DgDp_grad_i = diag[i] * p_ps;
      if (nonlinearTermFactor_ != ST::zero()) {
        DgDp_grad_i += as<Scalar>(1.5) * nonlinearTermFactor_ * p_ps * p_ps;
      }
      DgDp_grad[i] = DgDp_grad_i / s_bar[i];

    }
  }
  
}

// External constructor
template<class Scalar>
const Teuchos::RCP<TriKota::DiagonalROME<Scalar> >
createModel(
  const int globalDim,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType &g_offset
  )
{
  using Teuchos::RCP;

  const RCP<const Teuchos::Comm<Thyra::Ordinal> > comm =
    Teuchos::DefaultComm<Thyra::Ordinal>::getComm();

  const int numProcs = comm->getSize();
  TEST_FOR_EXCEPT_MSG( numProcs > globalDim,
    "Error, the number of processors can not be greater than the global"
    " dimension of the vectors!." );
  const int localDim = globalDim / numProcs;
  const int localDimRemainder = globalDim % numProcs;
  TEST_FOR_EXCEPT_MSG( localDimRemainder != 0,
    "Error, the number of processors must divide into the global number"
    " of elements exactly for now!." );

  const RCP<TriKota::DiagonalROME<Scalar> > model =
    Teuchos::rcp(new TriKota::DiagonalROME<Scalar>(localDim));
  const RCP<const Thyra::VectorSpaceBase<Scalar> > p_space = model->get_p_space(0);
  const RCP<Thyra::VectorBase<Scalar> > ps = createMember(p_space);
  const Scalar ps_val = 2.0;
  Thyra::V_S(ps.ptr(), ps_val);
  model->setSolutionVector(ps);
  model->setScalarOffset(g_offset);

  return model;
}


}

