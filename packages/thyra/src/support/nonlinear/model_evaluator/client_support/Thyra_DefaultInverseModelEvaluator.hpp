// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAUL_INVERSE_MODEL_EVALUATOR_HPP
#define THYRA_DEFAUL_INVERSE_MODEL_EVALUATOR_HPP

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_SpmdMultiVectorFileIO.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_Time.hpp"

namespace Thyra {

/** \brief This class wraps any ModelEvaluator object and adds a simple, but
 * fairly general, inverse-type of response function.
 *
 * The following response function is added to the end of the supported
 * response functions:

 \verbatim

  g_(getUnderlyingModel()->Ng())(x,p,...)
    = stateMultiplier*(1/nx)*sum((xs(i)*(x(i)-xt(i)))^2,i=0...nx-1)
    + paramMultiplier*(1/np)*sum((ps(i)*(p(i)-pt(i)))^2,i=0...np-1)

 \endverbatim 

 * where <tt>xt</tt> is the target vector for the states, <tt>xs</tt> is a
 * scaling vector (where some components may be zero), <tt.stateMultiplier</tt> is a
 * scalar for the state matching term, <tt>p</tt> is one of the parameter
 * subvectors supported by the underlying model, <tt>pt</tt> is a nomial
 * parameter vector for which violations are penalized, <tt>ps</tt> is the
 * scaling vector for the parameter change violations, and <tt>paramMultiplier</tt>
 * is a scalar for the parameter change violations.
 *
 * Note that <tt>this->Ng() == getUnderlyingModel()->Ng() + 1</tt>.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class DefaultInverseModelEvaluator
  : virtual public ModelEvaluatorDelegatorBase<Scalar>
  , virtual public Teuchos::ParameterListAcceptor
{
public:

  /** \brief State target vector <tt>xt</tt>. */
  STANDARD_COMPOSITION_MEMBERS( VectorBase<Scalar>, stateTarget );

  /** \brief Parameter base vector <tt>pt</tt>. */
  STANDARD_COMPOSITION_MEMBERS( VectorBase<Scalar>, parameterBase );

  /** \name Constructors/initializers/accessors/utilities. */
  //@{

  /** \brief . */
  DefaultInverseModelEvaluator();

  /** \brief . */
  DefaultInverseModelEvaluator(
    const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >       &thyraModel
    );

  /** \brief . */
  void initialize(
    const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >       &thyraModel
    );

  /** \brief . */
  void uninitialize(
    Teuchos::RefCountPtr<ModelEvaluator<Scalar> >             *thyraModel
    );

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList);
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> getParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RefCountPtr<const Teuchos::ParameterList> getValidParameters() const;

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  int Ng() const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_g_space(int j) const;
  /** \brief . */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const;
  /** \brief . */
  void evalModel(
    const ModelEvaluatorBase::InArgs<Scalar>    &inArgs
    ,const ModelEvaluatorBase::OutArgs<Scalar>  &outArgs
    ) const;

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  Teuchos::RefCountPtr<Teuchos::ParameterList>  paramList_;

  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > inv_g_space_;

  int            p_idx_;
  int            Ng_;

  double         stateMultiplier_; 
  std::string    stateTargetFileNameBase_;
  std::string    stateScaleFileNameBase_;
  double         parameterMultiplier_; 
  std::string    parameterBaseFileNameBase_;
  std::string    parameterScaleFileNameBase_;


  static const std::string ParameterSubvectorIndex_name_;
  static const int ParameterSubvectorIndex_default_;

  static const std::string StateMultiplier_name_;
  static const double StateMultiplier_default_;

  static const std::string StateTargetFileNameBase_name_;
  static const std::string StateTargetFileNameBase_default_;

  static const std::string StateScaleFileNameBase_name_;
  static const std::string StateScaleFileNameBase_default_;

  static const std::string ParameterMultiplier_name_;
  static const double ParameterMultiplier_default_;

  static const std::string ParameterBaseFileNameBase_name_;
  static const std::string ParameterBaseFileNameBase_default_;

  static const std::string ParameterScaleFileNameBase_name_;
  static const std::string ParameterScaleFileNameBase_default_;

};

// /////////////////////////////////
// Implementations

// Static data members

template<class Scalar>
const std::string
DefaultInverseModelEvaluator<Scalar>::ParameterSubvectorIndex_name_
= "Parameter Subvector Index";

template<class Scalar>
const int
DefaultInverseModelEvaluator<Scalar>::ParameterSubvectorIndex_default_
= 0;

template<class Scalar>
const std::string
DefaultInverseModelEvaluator<Scalar>::StateMultiplier_name_
= "State Multiplier";

template<class Scalar>
const double
DefaultInverseModelEvaluator<Scalar>::StateMultiplier_default_
= 1.0;

template<class Scalar>
const std::string
DefaultInverseModelEvaluator<Scalar>::StateTargetFileNameBase_name_
= "State Target File Name Base";

template<class Scalar>
const std::string
DefaultInverseModelEvaluator<Scalar>::StateTargetFileNameBase_default_
= "";

template<class Scalar>
const std::string
DefaultInverseModelEvaluator<Scalar>::StateScaleFileNameBase_name_
= "State Scale File Name Base";

template<class Scalar>
const std::string
DefaultInverseModelEvaluator<Scalar>::StateScaleFileNameBase_default_
= "";

template<class Scalar>
const std::string
DefaultInverseModelEvaluator<Scalar>::ParameterMultiplier_name_
= "Parameter Multiplier";

template<class Scalar>
const double
DefaultInverseModelEvaluator<Scalar>::ParameterMultiplier_default_
= 1e-6;

template<class Scalar>
const std::string
DefaultInverseModelEvaluator<Scalar>::ParameterBaseFileNameBase_name_
= "Parameter Base File Name Base";

template<class Scalar>
const std::string
DefaultInverseModelEvaluator<Scalar>::ParameterBaseFileNameBase_default_
= "";

template<class Scalar>
const std::string
DefaultInverseModelEvaluator<Scalar>::ParameterScaleFileNameBase_name_
= "Parameter Scale File Name Base";

template<class Scalar>
const std::string
DefaultInverseModelEvaluator<Scalar>::ParameterScaleFileNameBase_default_
= "";

// Constructors/initializers/accessors/utilities

template<class Scalar>
DefaultInverseModelEvaluator<Scalar>::DefaultInverseModelEvaluator()
  :p_idx_(0)
{}

template<class Scalar>
DefaultInverseModelEvaluator<Scalar>::DefaultInverseModelEvaluator(
  const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >   &thyraModel
  )
  :p_idx_(0)
{
  initialize(thyraModel);
}

template<class Scalar>
void DefaultInverseModelEvaluator<Scalar>::initialize(
  const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >   &thyraModel
  )
{
  this->ModelEvaluatorDelegatorBase<Scalar>::initialize(thyraModel);
  Ng_ = thyraModel->Ng()+1;
  inv_g_space_= thyraModel->get_x_space()->smallVecSpcFcty()->createVecSpc(1);
}

template<class Scalar>
void DefaultInverseModelEvaluator<Scalar>::uninitialize(
  Teuchos::RefCountPtr<ModelEvaluator<Scalar> >  *thyraModel
  )
{
  if(thyraModel) *thyraModel = this->getUnderlyingModel();
  this->ModelEvaluatorDelegatorBase<Scalar>::uninitialize();
}

// Overridden from Teuchos::ParameterListAcceptor

template<class Scalar>
void DefaultInverseModelEvaluator<Scalar>::setParameterList(
  Teuchos::RefCountPtr<Teuchos::ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(paramList.get()==NULL);
  paramList->validateParameters(*getValidParameters());
  paramList_ = paramList;
  p_idx_ = paramList_->get(
    ParameterSubvectorIndex_name_,ParameterSubvectorIndex_default_);
  stateMultiplier_ = paramList_->get(
    StateMultiplier_name_,StateMultiplier_default_);
  stateTargetFileNameBase_ = paramList_->get(
    StateTargetFileNameBase_name_,StateTargetFileNameBase_default_);
  stateScaleFileNameBase_ = paramList_->get(
    StateScaleFileNameBase_name_,StateScaleFileNameBase_default_);
  parameterMultiplier_ = paramList_->get(
    ParameterMultiplier_name_,ParameterMultiplier_default_);
  parameterBaseFileNameBase_ = paramList_->get(
    ParameterBaseFileNameBase_name_,ParameterBaseFileNameBase_default_);
  parameterScaleFileNameBase_ = paramList_->get(
    ParameterScaleFileNameBase_name_,ParameterScaleFileNameBase_default_);
  //
  Thyra::SpmdMultiVectorFileIO<Scalar> fileIO;
  if(stateTargetFileNameBase_.length()) {
    stateTarget_ = fileIO.readVectorFromFile(
      stateTargetFileNameBase_,this->get_x_space() );
  }
  if(stateScaleFileNameBase_.length()) {
    TEST_FOR_EXCEPT(true);
  }
  if(parameterBaseFileNameBase_.length()) {
    parameterBase_ = fileIO.readVectorFromFile(
      parameterBaseFileNameBase_,this->get_p_space(p_idx_) );
  }
  if(parameterScaleFileNameBase_.length()) {
    TEST_FOR_EXCEPT(true);
  }
}

template<class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList>
DefaultInverseModelEvaluator<Scalar>::getParameterList()
{
  return paramList_;
}

template<class Scalar>
Teuchos::RefCountPtr<Teuchos::ParameterList>
DefaultInverseModelEvaluator<Scalar>::unsetParameterList()
{
  Teuchos::RefCountPtr<Teuchos::ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}

template<class Scalar>
Teuchos::RefCountPtr<const Teuchos::ParameterList>
DefaultInverseModelEvaluator<Scalar>::getParameterList() const
{
  return paramList_;
}

template<class Scalar>
Teuchos::RefCountPtr<const Teuchos::ParameterList>
DefaultInverseModelEvaluator<Scalar>::getValidParameters() const
{
  static Teuchos::RefCountPtr<Teuchos::ParameterList> pl;
  if(pl.get()==NULL) {
    pl = Teuchos::rcp(new Teuchos::ParameterList());
    pl->set(
      StateMultiplier_name_,StateMultiplier_default_
      ,"stateMultiplier"
      );
    pl->set(
      StateTargetFileNameBase_name_,StateTargetFileNameBase_default_
      ,"Base-name of file(s) that contain state matching vector xt." 
      );
    pl->set(
      StateScaleFileNameBase_name_,StateScaleFileNameBase_default_
      ,"Base-name of file(s) that contain state scaling vector xs." 
      );
    pl->set(
      ParameterMultiplier_name_,ParameterMultiplier_default_
      ,"parameterMultiplier"
      );
    pl->set(
      ParameterBaseFileNameBase_name_,ParameterBaseFileNameBase_default_
      ,"Base-name of file(s) that contain parameter base vector pt." 
      );
    pl->set(
      ParameterScaleFileNameBase_name_,ParameterScaleFileNameBase_default_
      ,"Base-name of file(s) that contain parameter scaling vector ps." 
      );
  }
  return pl;
}

// Overridden from ModelEvaulator.

template<class Scalar>
int DefaultInverseModelEvaluator<Scalar>::Ng() const
{
  return Ng_;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
DefaultInverseModelEvaluator<Scalar>::get_g_space(int j) const
{
  if(j==Ng_-1)
    return inv_g_space_;
  return this->getUnderlyingModel()->get_g_space(j);
}

template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
DefaultInverseModelEvaluator<Scalar>::createOutArgs() const
{
  typedef ModelEvaluatorBase MEB;
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  const MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel->createOutArgs();
  const int Np = wrappedOutArgs.Np(), Ng = wrappedOutArgs.Ng();
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(Np,Ng+1);
  outArgs.setSupports(wrappedOutArgs);
  outArgs.setSupports(MEB::OUT_ARG_DgDx,Ng,MEB::DERIV_TRANS_MV_BY_ROW);
  outArgs.setSupports(MEB::OUT_ARG_DgDp,Ng,p_idx_,MEB::DERIV_TRANS_MV_BY_ROW);
  return outArgs;
}

template<class Scalar>
void DefaultInverseModelEvaluator<Scalar>::evalModel(
  const ModelEvaluatorBase::InArgs<Scalar>     &inArgs
  ,const ModelEvaluatorBase::OutArgs<Scalar>   &outArgs
  ) const
{
  typedef ModelEvaluatorBase MEB;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  typedef typename ST::magnitudeType ScalarMag;

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_BEGIN(
    "Thyra::DefaultInverseModelEvaluator",inArgs,outArgs
    );

  //
  // Compute all of the desired functions in the base model
  //

  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nComputing the base point ...\n";

  MEB::InArgs<Scalar>  wrappedInArgs = inArgs;
  MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel->createOutArgs();
  wrappedOutArgs.setArgs(outArgs,true);
  thyraModel->evalModel(wrappedInArgs,wrappedOutArgs);
  
  bool failed = wrappedOutArgs.isFailed();
  
  if(!failed) {
    //
    // Compute the inverse response function and its derivatives if asked to
    // do so.
    //
    const VectorBase<Scalar>
      &x = *inArgs.get_x(),
      &p = *inArgs.get_p(p_idx_);
    //
    VectorBase<Scalar>
      *g_out = outArgs.get_g(outArgs.Ng()-1).get();
    MultiVectorBase<Scalar>
      *DgDx_trans_out = get_mv(
        outArgs.get_DgDx(outArgs.Ng()-1),"DgDx",MEB::DERIV_TRANS_MV_BY_ROW
        ).get();
    MultiVectorBase<Scalar>
      *DgDp_trans_out = get_mv(
        outArgs.get_DgDp(outArgs.Ng()-1,p_idx_),"DgDp",MEB::DERIV_TRANS_MV_BY_ROW
        ).get();
    //
    const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
      x_space = this->get_x_space(),
      p_space = this->get_p_space(p_idx_);
    const Index
      nx = x_space->dim(),
      np = p_space->dim();
    Teuchos::RefCountPtr<VectorBase<Scalar> > x_minus_xt, p_minus_pt;
    if( g_out || DgDx_trans_out ) {
      x_minus_xt = createMember(x_space);
      V_VmV( &*x_minus_xt, x, *stateTarget_ );
    }
    if( g_out || DgDp_trans_out ) {
      p_minus_pt = createMember(p_space);
      V_VmV( &*p_minus_pt, p, *parameterBase_ );
    }
    if(g_out) {
      if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
        *out << "\nComputing inverse response function g ...\n";
      set_ele(
        0
        ,stateMultiplier_*(1.0/nx)*dot(*x_minus_xt,*x_minus_xt)
        +parameterMultiplier_*(1.0/np)*dot(*p_minus_pt,*p_minus_pt)
        ,g_out
        );
    }
    if(DgDx_trans_out) {
      // trans(DgDx) = stateMultiplier*(2/nx)*(x-xt)
      if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
        *out << "\nComputing inverse response function derivative DgDx_trans ...\n";
      V_StV( &*DgDx_trans_out->col(0), Scalar(stateMultiplier_*(2.0/nx)), *x_minus_xt );
    }
    if(DgDp_trans_out) {
      // trans(DgDp) = parameterMultiplier*(2/np)*(p-pt)
      if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
        *out << "\nComputing inverse response function derivative DgDp_trans ...\n";
      V_StV( &*DgDp_trans_out->col(0), Scalar(parameterMultiplier_*(2.0/np)), *p_minus_pt );
    }
  }

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();
  
}

// Public functions overridden from Teuchos::Describable

template<class Scalar>
std::string DefaultInverseModelEvaluator<Scalar>::description() const
{
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  std::ostringstream oss;
  oss << "Thyra::DefaultInverseModelEvaluator{";
  oss << "thyraModel=";
  if(thyraModel.get())
    oss << "\'"<<thyraModel->description()<<"\'";
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}

} // namespace Thyra

#endif // THYRA_DEFAUL_INVERSE_MODEL_EVALUATOR_HPP
