// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_MODEL_EVALUATOR_BASE_DEF_HPP
#define THYRA_MODEL_EVALUATOR_BASE_DEF_HPP


#include "Thyra_ModelEvaluatorBase_decl.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"


namespace Thyra {


namespace ModelEvaluatorHelperPack {


template<class Scalar>
inline
RCP<const Thyra::VectorBase<Scalar> >
condCloneVec(
  const RCP<const Thyra::VectorBase<Scalar> > &vec,
  bool cloneObject
  )
{
  if(cloneObject)
    return vec->clone_v();
  return vec;
}

template<class Scalar>
inline
RCP<const Thyra::MultiVectorBase<Scalar> >
condCloneMultiVec(
  const RCP<const Thyra::MultiVectorBase<Scalar> > &vec,
  bool cloneObject
  )
{
  if(cloneObject)
    return vec->clone_mv();
  return vec;
}

inline
RCP<const Stokhos::ProductEpetraVector >
condCloneVec_mp(
  const RCP<const Stokhos::ProductEpetraVector > &vec,
  bool cloneObject
  )
{
  if(cloneObject)
  {
    printf("Warning:  clone_v not implemented for ProductEpetraVector: %s %d\n",__FILE__,__LINE__);
    //return vec->clone_v();  //JF clone_v not implemented for ProductEpetraVector
  }
  return vec;
}

} // namespace ModelEvaluatorHelperPack


//
// ModelEvaluatorBase::InArgs
//


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>::InArgs()
  :modelEvalDescription_("WARNING!  THIS INARGS OBJECT IS UNINITIALIZED!")
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Teuchos::ScalarTraits<typename ST::magnitudeType> SMT;
  std::fill_n(&supports_[0],NUM_E_IN_ARGS_MEMBERS,false);
  t_     = SMT::zero();
  alpha_ = ST::zero();
  beta_  = ST::zero();
  step_size_ = ST::zero();
  stage_number_ = ST::one();
}


template<class Scalar>
int ModelEvaluatorBase::InArgs<Scalar>::Np() const
{ return p_.size(); }

template<class Scalar>
int ModelEvaluatorBase::InArgs<Scalar>::Ng() const
{ return g_multiplier_.size(); }

template<class Scalar>
bool ModelEvaluatorBase::InArgs<Scalar>::supports(EInArgsMembers arg) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_IN_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"model = \'"<<modelEvalDescription_
    <<"\': Error, arg="<<toString(arg)<<" is invalid!"
    );
  return supports_[arg];
}

template<class Scalar>
bool ModelEvaluatorBase::InArgs<Scalar>::supports(EInArgs_p_mp /* arg */, int l) const
{
  assert_l(l);
  return supports_p_mp_[l];
}

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_x_dot_dot(
  const RCP<const VectorBase<Scalar> > &x_dot_dot
  )
{ assert_supports(IN_ARG_x_dot_dot); x_dot_dot_ = x_dot_dot; }

template<class Scalar>
RCP<const VectorBase<Scalar> >
ModelEvaluatorBase::InArgs<Scalar>::get_x_dot_dot() const
{ assert_supports(IN_ARG_x_dot_dot); return x_dot_dot_; }

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_x_dot(
  const RCP<const VectorBase<Scalar> > &x_dot
  )
{ assert_supports(IN_ARG_x_dot); x_dot_ = x_dot; }


template<class Scalar>
RCP<const VectorBase<Scalar> >
ModelEvaluatorBase::InArgs<Scalar>::get_x_dot() const
{ assert_supports(IN_ARG_x_dot); return x_dot_; }


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_x(
  const RCP<const VectorBase<Scalar> > &x
  )
{ assert_supports(IN_ARG_x); x_ = x; }


template<class Scalar>
RCP<const VectorBase<Scalar> >
ModelEvaluatorBase::InArgs<Scalar>::get_x() const
{ assert_supports(IN_ARG_x); return x_; }


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_x_direction(
  const RCP<const MultiVectorBase<Scalar> > &x_direction
  )
{ assert_supports(IN_ARG_x); x_direction_ = x_direction; }


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_p_direction(
  int l, const RCP<const MultiVectorBase<Scalar> > &p_direction_l
  )
{ assert_l(l); p_direction_[l] = p_direction_l; }


template<class Scalar>
RCP<const MultiVectorBase<Scalar> >
ModelEvaluatorBase::InArgs<Scalar>::get_x_direction() const
{ assert_supports(IN_ARG_x); return x_direction_; }


template<class Scalar>
RCP<const MultiVectorBase<Scalar> >
ModelEvaluatorBase::InArgs<Scalar>::get_p_direction(int l) const
{ assert_l(l); return p_direction_[l]; }


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_x_dot_mp(
  const RCP<const Stokhos::ProductEpetraVector > &x_dot_mp
  )
{ assert_supports(IN_ARG_x_dot_mp); x_dot_mp_ = x_dot_mp; }


template<class Scalar>
RCP<const Stokhos::ProductEpetraVector >
ModelEvaluatorBase::InArgs<Scalar>::get_x_dot_mp() const
{ assert_supports(IN_ARG_x_dot_mp); return x_dot_mp_; }


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_x_mp(
  const RCP<const Stokhos::ProductEpetraVector > &x_mp
  )
{ assert_supports(IN_ARG_x_mp); x_mp_ = x_mp; }


template<class Scalar>
RCP<const Stokhos::ProductEpetraVector >
ModelEvaluatorBase::InArgs<Scalar>::get_x_mp() const
{ assert_supports(IN_ARG_x_mp); return x_mp_; }

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_f_multiplier(
  const RCP<const VectorBase<Scalar> > &f_multiplier
  )
{ assert_supports(IN_ARG_x); f_multiplier_ = f_multiplier; }

template<class Scalar>
RCP<const VectorBase<Scalar> >
ModelEvaluatorBase::InArgs<Scalar>::get_f_multiplier() const
{ assert_supports(IN_ARG_x); return f_multiplier_; }

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_g_multiplier(
  int j, const RCP<const VectorBase<Scalar> > &g_multiplier
  )
{
  assert_j(j);
  assert_supports(IN_ARG_x);
  g_multiplier_[j] = g_multiplier;
}

template<class Scalar>
RCP<const VectorBase<Scalar> >
ModelEvaluatorBase::InArgs<Scalar>::get_g_multiplier(int j) const
{
  assert_j(j);
  assert_supports(IN_ARG_x);
  return g_multiplier_[j];
}

#ifdef HAVE_THYRA_ME_POLYNOMIAL

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_x_dot_poly(
  const RCP<const Teuchos::Polynomial< VectorBase<Scalar> > > &x_dot_poly
  )
{ assert_supports(IN_ARG_x_dot_poly); x_dot_poly_ = x_dot_poly; }


template<class Scalar>
RCP<const Teuchos::Polynomial< VectorBase<Scalar> > >
ModelEvaluatorBase::InArgs<Scalar>::get_x_dot_poly() const
{ assert_supports(IN_ARG_x_dot_poly); return x_dot_poly_; }


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_x_poly(
  const RCP<const Teuchos::Polynomial< VectorBase<Scalar> > > &x_poly
  )
{ assert_supports(IN_ARG_x_poly); x_poly_ = x_poly; }


template<class Scalar>
RCP<const Teuchos::Polynomial< VectorBase<Scalar> > >
ModelEvaluatorBase::InArgs<Scalar>::get_x_poly() const
{ assert_supports(IN_ARG_x_poly); return x_poly_; }


#endif // HAVE_THYRA_ME_POLYNOMIAL

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_p(
  int l, const RCP<const VectorBase<Scalar> > &p_l
  )
{ assert_l(l); p_[l] = p_l; }


template<class Scalar>
RCP<const VectorBase<Scalar> >
ModelEvaluatorBase::InArgs<Scalar>::get_p(int l) const
{ assert_l(l); return p_[l]; }


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_p_mp(
  int l, const RCP<const Stokhos::ProductEpetraVector > &p_mp_l
  )
{ assert_supports(IN_ARG_p_mp, l); p_mp_[l] = p_mp_l; }

template<class Scalar>
RCP<const Stokhos::ProductEpetraVector >
ModelEvaluatorBase::InArgs<Scalar>::get_p_mp(int l) const
{ assert_supports(IN_ARG_p_mp, l); return p_mp_[l]; }


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_t( ScalarMag t )
{ assert_supports(IN_ARG_t); t_ = t; }


template<class Scalar>
typename ModelEvaluatorBase::InArgs<Scalar>::ScalarMag
ModelEvaluatorBase::InArgs<Scalar>::get_t() const
{ assert_supports(IN_ARG_t); return t_; }


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_alpha( Scalar alpha )
{ assert_supports(IN_ARG_alpha); alpha_ = alpha; }


template<class Scalar>
Scalar ModelEvaluatorBase::InArgs<Scalar>::get_alpha() const
{ assert_supports(IN_ARG_alpha); return alpha_; }


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_beta( Scalar beta )
{ assert_supports(IN_ARG_beta); beta_ = beta; }


template<class Scalar>
Scalar ModelEvaluatorBase::InArgs<Scalar>::get_beta() const
{ assert_supports(IN_ARG_beta); return beta_; }

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_W_x_dot_dot_coeff( Scalar W_x_dot_dot_coeff )
{ assert_supports(IN_ARG_W_x_dot_dot_coeff); W_x_dot_dot_coeff_ = W_x_dot_dot_coeff; }


template<class Scalar>
Scalar ModelEvaluatorBase::InArgs<Scalar>::get_W_x_dot_dot_coeff() const
{ assert_supports(IN_ARG_W_x_dot_dot_coeff); return W_x_dot_dot_coeff_; }

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_step_size( Scalar step_size)
{ assert_supports(IN_ARG_step_size); step_size_ = step_size; }

template<class Scalar>
Scalar ModelEvaluatorBase::InArgs<Scalar>::get_step_size() const
{ assert_supports(IN_ARG_step_size); return step_size_; }

template<class Scalar>
Scalar ModelEvaluatorBase::InArgs<Scalar>::get_stage_number() const
{ assert_supports(IN_ARG_stage_number); return stage_number_; }


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::set_stage_number( Scalar stage_number)
{ assert_supports(IN_ARG_stage_number); stage_number_ = stage_number; }


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::setArgs(
  const InArgs<Scalar>& inArgs, bool ignoreUnsupported, bool cloneObjects
  )
{
  using ModelEvaluatorHelperPack::condCloneVec;
  using ModelEvaluatorHelperPack::condCloneMultiVec;
  using ModelEvaluatorHelperPack::condCloneVec_mp;
  if( inArgs.supports(IN_ARG_x_dot_dot) && nonnull(inArgs.get_x_dot_dot()) ) {
    if(supports(IN_ARG_x_dot_dot) || !ignoreUnsupported)
      set_x_dot_dot(condCloneVec(inArgs.get_x_dot_dot(),cloneObjects));
  }
  if( inArgs.supports(IN_ARG_x_dot) && nonnull(inArgs.get_x_dot()) ) {
    if(supports(IN_ARG_x_dot) || !ignoreUnsupported)
      set_x_dot(condCloneVec(inArgs.get_x_dot(),cloneObjects));
  }
  if( inArgs.supports(IN_ARG_x_dot_mp) && nonnull(inArgs.get_x_dot_mp()) ) {
    if(supports(IN_ARG_x_dot_mp) || !ignoreUnsupported)
      set_x_dot_mp(condCloneVec_mp(inArgs.get_x_dot_mp(),cloneObjects));
  }
  if( inArgs.supports(IN_ARG_x) && nonnull(inArgs.get_x()) ) {
    if(supports(IN_ARG_x) || !ignoreUnsupported)
      set_x(condCloneVec(inArgs.get_x(),cloneObjects));
  }
  if( inArgs.supports(IN_ARG_x_mp) && nonnull(inArgs.get_x_mp()) ) {
    if(supports(IN_ARG_x_mp) || !ignoreUnsupported)
      set_x_mp(condCloneVec_mp(inArgs.get_x_mp(),cloneObjects));
  }
  if( inArgs.supports(IN_ARG_x) && nonnull(inArgs.get_x_direction()) ) {
    if(supports(IN_ARG_x) || !ignoreUnsupported)
      set_x_direction(condCloneMultiVec(inArgs.get_x_direction(),cloneObjects));
  }
  if( inArgs.supports(IN_ARG_x) && nonnull(inArgs.get_f_multiplier()) ) {
    if(supports(IN_ARG_x) || !ignoreUnsupported)
      set_f_multiplier(condCloneVec(inArgs.get_f_multiplier(),cloneObjects));
  }
  const int min_Ng = TEUCHOS_MIN(this->Ng(),inArgs.Ng());
  for (int j = 0; j < min_Ng; ++j) {
    if (nonnull(inArgs.get_g_multiplier(j)))
      set_g_multiplier(j,condCloneVec(inArgs.get_g_multiplier(j),cloneObjects));
  }
#ifdef HAVE_THYRA_ME_POLYNOMIAL
  if( inArgs.supports(IN_ARG_x_dot_poly) && nonnull(inArgs.get_x_dot_poly()) ) {
    if(supports(IN_ARG_x_dot_poly) || !ignoreUnsupported) {
      TEUCHOS_TEST_FOR_EXCEPT(
        cloneObjects && "Have not implemented cloning for x_dot_poly yet!" );
      set_x_dot_poly(inArgs.get_x_dot_poly());
    }
  }
  if( inArgs.supports(IN_ARG_x_poly) && nonnull(inArgs.get_x_poly()) ) {
    if(supports(IN_ARG_x_poly) || !ignoreUnsupported) {
      TEUCHOS_TEST_FOR_EXCEPT(
        cloneObjects && "Have not implemented cloning for x_poly yet!" );
      set_x_poly(inArgs.get_x_poly());
    }
  }
#endif // HAVE_THYRA_ME_POLYNOMIAL
  const int min_Np = TEUCHOS_MIN(this->Np(),inArgs.Np());
  for (int l = 0; l < min_Np; ++l) {
    if (nonnull(inArgs.get_p(l)))
      set_p(l,condCloneVec(inArgs.get_p(l),cloneObjects));
  }
  for (int l = 0; l < min_Np; ++l) {
    if (nonnull(inArgs.get_p_direction(l)))
      set_p_direction(l,condCloneMultiVec(inArgs.get_p_direction(l),cloneObjects));
  }
  for (int l = 0; l < min_Np; ++l) {
    if (inArgs.supports(IN_ARG_p_mp,l)) {
      if (nonnull(inArgs.get_p_mp(l)))
        set_p_mp(l,condCloneVec_mp(inArgs.get_p_mp(l),cloneObjects));
    }
  }
  if (inArgs.supports(IN_ARG_t)) {
    if(supports(IN_ARG_t) || !ignoreUnsupported)
      set_t(inArgs.get_t());
  }
  if (inArgs.supports(IN_ARG_alpha)) {
    if(supports(IN_ARG_alpha) || !ignoreUnsupported)
      set_alpha(inArgs.get_alpha());
  }
  if (inArgs.supports(IN_ARG_beta)) {
    if(supports(IN_ARG_beta) || !ignoreUnsupported)
      set_beta(inArgs.get_beta());
  }
  if (inArgs.supports(IN_ARG_W_x_dot_dot_coeff)) {
    if(supports(IN_ARG_W_x_dot_dot_coeff) || !ignoreUnsupported)
      set_W_x_dot_dot_coeff(inArgs.get_W_x_dot_dot_coeff());
  }
  if (inArgs.supports(IN_ARG_step_size)) {
    if(supports(IN_ARG_step_size) || !ignoreUnsupported)
      set_step_size(inArgs.get_step_size());
  }
  if (inArgs.supports(IN_ARG_stage_number)) {
    if(supports(IN_ARG_stage_number) || !ignoreUnsupported)
      set_stage_number(inArgs.get_stage_number());
  }
  // Extended inArgs
  if (extended_inargs_.size() > 0)
    TEUCHOS_TEST_FOR_EXCEPTION(cloneObjects,
                               std::runtime_error,
                               "Extended InArgs does not support cloning!");
  this->extended_inargs_ = inArgs.extended_inargs_;
}


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::assertSameSupport(
  const InArgs<Scalar> &inArgs
  ) const
{
  for ( int inArg_i = 0; inArg_i < NUM_E_IN_ARGS_MEMBERS; ++inArg_i ) {
    const EInArgsMembers inArg_arg = static_cast<EInArgsMembers>(inArg_i);
    const std::string inArg_name = toString(inArg_arg);
    TEUCHOS_TEST_FOR_EXCEPTION(
      supports(inArg_arg) != inArgs.supports(inArg_arg), std::logic_error,
      "Error, the input argument "<<inArg_name<<" with support "<<inArgs.supports(inArg_arg)<<"\n"
      "in the InArgs object for the model:\n\n"
      "  "<<inArgs.modelEvalDescription()<<"\n\n"
      "is not the same the argument "<<inArg_name<<" with support "<<supports(inArg_arg)<<"\n"
      "in the InArgs object for the model:\n\n"
      "  "<<modelEvalDescription()<<"\n\n"
      "and these two InArgs objects are not compatible!"
      );
  }
  TEUCHOS_ASSERT_EQUALITY( this->Np(), inArgs.Np() );
}


template<class Scalar>
std::string ModelEvaluatorBase::InArgs<Scalar>::modelEvalDescription() const
{
  return modelEvalDescription_;
}


template<class Scalar>
std::string ModelEvaluatorBase::InArgs<Scalar>::description() const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  std::ostringstream oss;
  oss
    << "Thyra::ModelEvaluatorBase::InArgs<"<<ST::name()<<">"
    << "{"
    << "model="<<modelEvalDescription_
    << ",Np="<<Np()
    << "}";
  return oss.str();
}


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::describe(
  Teuchos::FancyOStream &out_arg, const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using std::endl;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  using Teuchos::OSTab;
  using Teuchos::describe;
  using Teuchos::includesVerbLevel;
  typedef RCP<const VectorBase<Scalar> > CV_ptr;

  if(verbLevel == Teuchos::VERB_NONE)
    return;

  RCP<Teuchos::FancyOStream>
    out = Teuchos::rcp(&out_arg,false);
  const bool dump_x = includesVerbLevel(verbLevel,Teuchos::VERB_HIGH);
  const Teuchos::EVerbosityLevel x_verbLevel =
    dump_x?Teuchos::VERB_EXTREME:verbLevel;
  const bool print_x_nrm = includesVerbLevel(verbLevel,Teuchos::VERB_LOW);
  const bool dump_p = includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM);
  const Teuchos::EVerbosityLevel p_verbLevel =
    dump_p?Teuchos::VERB_EXTREME:verbLevel;
  const bool print_p_nrm = includesVerbLevel(verbLevel,Teuchos::VERB_LOW);
  OSTab tab(out);

  *out <<"Thyra::ModelEvaluatorBase::InArgs<"<<ST::name()<<">:\n";
  tab.incrTab();

  *out <<"model = " << modelEvalDescription_ << "\n";
  *out <<"Np = " << Np() << "\n";

  CV_ptr x_dot_dot;
  if ( this->supports(IN_ARG_x_dot_dot) && !is_null(x_dot_dot=get_x_dot_dot()) ) {
    *out << "x_dot_dot = " << Teuchos::describe(*x_dot_dot,x_verbLevel);
    if (print_x_nrm)
      *out << "||x_dot_dot|| = " << norm(*x_dot_dot) << endl;
  }

  CV_ptr x_dot;
  if ( this->supports(IN_ARG_x_dot) && !is_null(x_dot=get_x_dot()) ) {
    *out << "x_dot = " << Teuchos::describe(*x_dot,x_verbLevel);
    if (print_x_nrm)
      *out << "||x_dot|| = " << norm(*x_dot) << endl;
  }

  CV_ptr x;
  if ( this->supports(IN_ARG_x) && !is_null(x=get_x()) ) {
    *out << "x = " << Teuchos::describe(*x,x_verbLevel);
    if (print_x_nrm)
      *out << "||x|| = " << norm(*x) << endl;
  }

  if (print_x_nrm) {
    for( int l = 0; l < Np(); ++l ) {
      CV_ptr p_l;
      if ( !is_null(p_l = this->get_p(l)) ) {
        *out << "p("<<l<<") = " << Teuchos::describe(*p_l,p_verbLevel);
        if (print_p_nrm)
          *out << "||p("<<l<<")|| = " << norm(*p_l) << endl;
      }
    }
  }

  if (includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM)) {
    if (this->supports(IN_ARG_t)) {
      *out << "t = " << t_ << endl;
    }
    if (this->supports(IN_ARG_alpha)) {
      *out << "alpha = " << alpha_ << endl;
    }
    if (this->supports(IN_ARG_beta)) {
      *out << "beta = " << beta_ << endl;
    }
    if (this->supports(IN_ARG_W_x_dot_dot_coeff)) {
      *out << "W_x_dot_dot_coeff = " << W_x_dot_dot_coeff_ << endl;
    }
    if (this->supports(IN_ARG_step_size)) {
      *out << "step_size = " << step_size_ << endl;
    }
    if (this->supports(IN_ARG_stage_number)) {
      *out << "stage_number = " << stage_number_ << endl;
    }
  }

}


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::_setModelEvalDescription(
  const std::string &modelEvalDescription_in
  )
{
  modelEvalDescription_ = modelEvalDescription_in;
}


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::_set_Np_Ng(int Np_in, int Ng_in)
{
  p_.resize(Np_in);
  p_direction_.resize(Np_in);
  p_mp_.resize(Np_in);
  supports_p_mp_.resize(Np_in);
  g_multiplier_.resize(Ng_in);
}


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::_setSupports(
  EInArgsMembers arg, bool supports_in
  )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_IN_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"model = \'"<<modelEvalDescription_
    <<"\': Error, arg="<<toString(arg)<<" is invalid!");
  supports_[arg] = supports_in;
}

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::_setSupports(
  EInArgs_p_mp /* arg */, int l, bool supports_in
  )
{
  assert_l(l);
  supports_p_mp_[l] = supports_in;
}


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::_setSupports(
  const InArgs<Scalar>& inArgs, const int Np_in
  )
{
  std::copy(
    &inArgs.supports_[0],
    &inArgs.supports_[0] + NUM_E_IN_ARGS_MEMBERS, &supports_[0] );
  this->_set_Np_Ng( Np_in >= 0 ? Np_in : inArgs.Np(), inArgs.Ng() );
}


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::_setUnsupportsAndRelated(
  EInArgsMembers arg
  )
{
  switch(arg) {
    case IN_ARG_x: {
      this->_setSupports(IN_ARG_x_dot_dot,false);
      this->_setSupports(IN_ARG_x_dot,false);
      this->_setSupports(IN_ARG_x_dot_poly,false);
      this->_setSupports(IN_ARG_alpha,false);
      this->_setSupports(IN_ARG_beta,false);
      this->_setSupports(IN_ARG_W_x_dot_dot_coeff,false);
      this->_setSupports(IN_ARG_step_size,false);
      this->_setSupports(IN_ARG_stage_number,false);
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(
        true ,std::logic_error,
        "Error, can not handle args other than IN_ARG_x yet!"
        );
  }
  this->_setSupports(arg,false);
}


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::assert_supports(
  EInArgsMembers arg
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !supports_[arg], std::logic_error
    ,"Thyra::ModelEvaluatorBase::InArgs<"
    << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_supports(arg): "
    "model = \'"<<modelEvalDescription_<<"\': Error, "
    "The argument arg = " << toString(arg) << " is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::assert_supports(
  EInArgs_p_mp /* arg */, int l
  ) const
{
  assert_l(l);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !supports_p_mp_[l], std::logic_error
    ,"Thyra::ModelEvaluatorBase::InArgs<"
    << Teuchos::ScalarTraits<Scalar>::name() <<">::assert_supports(IN_ARG_p_mp,1): "
    "model = \'"<<modelEvalDescription_<<"\': Error, "
    "The argument p_mp(l) with index l = " << l << " is not supported!"
    );
}


template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::assert_l(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !( 0 <= l && l < Np() ), std::logic_error
    ,"Thyra::ModelEvaluatorBase::InArgs<Scalar>::assert_l(l):\n\n"
    " model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The parameter l = " << l << " is not in the range [0,"<<Np()<<")!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::InArgs<Scalar>::assert_j(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !( 0 <= j && j < Ng() ), std::logic_error
    ,"Thyra::ModelEvaluatorBase::InArgs<Scalar>::assert_j(j):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The auxiliary function g("<<j<<")"
    " is not in the range [0,"<<Ng()<<")!"
    );
}

//
// ModelEvaluatorBase::DerivativeMultiVector
//


template<class Scalar>
std::string ModelEvaluatorBase::DerivativeMultiVector<Scalar>::description() const
{
  using std::endl;
  std::ostringstream oss;
  oss << "DerivativeMultiVector{";
  if (is_null(getMultiVector())) {
    oss << "NULL";
  }
  else {
    oss
      << "multiVec=" << getMultiVector()->description()
      << ",orientation=" << toString(getOrientation());
  }
  oss << "}";
  return oss.str();
}


template<class Scalar>
void ModelEvaluatorBase::DerivativeMultiVector<Scalar>::describe(
  Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using std::endl;
  using Teuchos::describe;
  Teuchos::OSTab tab1(out);
  out << "DerivativeMultiVector\n";
  Teuchos::OSTab tab2(out);
  out
    << "multiVec = "
    << describe(*getMultiVector(),verbLevel)
    << "orientation = "
    << toString(getOrientation()) << endl;
}


// 2007/06/12: rabartl: The above description() and describe(...) functions
// have to be defined here and not in the class DerivativeMultiVector since it
// relies on the non-member function
// toString(ModelEvaluatorBase::EDerivativeMultiVectorOrientation) which is
// defined after the class definition for ModelEvaluatorBase.  This was caught
// by the intel compiler.  I am not sure why this worked with gcc.


//
// ModelEvaluatorBase::Derivative
//


template<class Scalar>
std::string
ModelEvaluatorBase::Derivative<Scalar>::description() const
{
  using std::endl;
  std::ostringstream oss;
  oss << "Derivative{";
  if (isEmpty()) {
    oss << "NULL";
  }
  else if (!is_null(getLinearOp())) {
    oss << "linearOp=" << getLinearOp()->description();
  }
  else {
    oss << "derivMultiVec=" << getDerivativeMultiVector().description();
  }
  oss << "}";
  return oss.str();
}


template<class Scalar>
void ModelEvaluatorBase::Derivative<Scalar>::describe(
  Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  using std::endl;
  using Teuchos::describe;
  Teuchos::OSTab tab1(out);
  out << "Derivative:";
  if (isEmpty()) {
    out << " NULL\n";
  }
  else if (!is_null(getLinearOp())) {
    out
      << endl
      << "linearOp = " << describe(*getLinearOp(),verbLevel);
  }
  else {
    out
      << endl
      << "derivMultiVec = ";
    getDerivativeMultiVector().describe(out,verbLevel);
  }
}


//
// ModelEvaluatorBase::OutArgs
//


template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>::OutArgs()
  :modelEvalDescription_("WARNING!  THIS OUTARGS OBJECT IS UNINITIALIZED!"),
   isFailed_(false)
{
  std::fill_n(&supports_[0],NUM_E_OUT_ARGS_MEMBERS,false);
#ifdef Thyra_BUILD_HESSIAN_SUPPORT
  this->_setHessianSupports(false);
#endif
}


template<class Scalar>
int ModelEvaluatorBase::OutArgs<Scalar>::Np() const
{ return DfDp_.size(); }


template<class Scalar>
int ModelEvaluatorBase::OutArgs<Scalar>::Ng() const
{ return g_.size(); }


template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgsMembers arg
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"model = \'"<<modelEvalDescription_
    <<"\': Error, arg="<<toString(arg)<<" is invalid!"
    );
  return supports_[arg];
}


template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgsDfDp /* arg */, int l
  ) const
{
  assert_l(l);
  return supports_DfDp_[l];
}


template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgsDgDx_dot /* arg */, int j
  ) const
{
  assert_j(j);
  return supports_DgDx_dot_[j];
}


template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgsDgDx /* arg */, int j
  ) const
{
  assert_j(j);
  return supports_DgDx_[j];
}


template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgsDgDp /* arg */, int j, int l
  ) const
{
  assert_j(j);
  assert_l(l);
  return supports_DgDp_[ j*Np() + l ];
}


#ifdef Thyra_BUILD_HESSIAN_SUPPORT

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_hess_vec_prod_f_xx /* arg */
  ) const
{
  return supports_hess_vec_prod_f_xx_;
}

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_hess_vec_prod_f_xp /* arg */, int l
  ) const
{
  assert_l(l);
  return supports_hess_vec_prod_f_xp_[l];
}

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_hess_vec_prod_f_px /* arg */, int l
  ) const
{
  assert_l(l);
  return supports_hess_vec_prod_f_px_[l];
}

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_hess_vec_prod_f_pp /* arg */, int l1, int l2
  ) const
{
  assert_l(l1);
  assert_l(l2);
  return supports_hess_vec_prod_f_pp_[ l1*Np() + l2 ];
}

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_hess_vec_prod_g_xx /* arg */, int j
  ) const
{
  assert_j(j);
  return supports_hess_vec_prod_g_xx_[j];
}

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_hess_vec_prod_g_xp /* arg */, int j, int l
  ) const
{
  assert_j(j);
  assert_l(l);
  return supports_hess_vec_prod_g_xp_[ j*Np() + l ];
}

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_hess_vec_prod_g_px /* arg */, int j, int l
  ) const
{
  assert_j(j);
  assert_l(l);
  return supports_hess_vec_prod_g_px_[ j*Np() + l ];
}

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_hess_vec_prod_g_pp /* arg */, int j, int l1, int l2
  ) const
{
  assert_j(j);
  assert_l(l1);
  assert_l(l2);
  return supports_hess_vec_prod_g_pp_[ j*Np()*Np() + l1*Np() + l2 ];
}

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_hess_f_xx /* arg */
  ) const
{
  return supports_hess_f_xx_;
}

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_hess_f_xp /* arg */, int l
  ) const
{
  assert_l(l);
  return supports_hess_f_xp_[l];
}

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_hess_f_pp /* arg */, int l1, int l2
  ) const
{
  assert_l(l1);
  assert_l(l2);
  return supports_hess_f_pp_[ l1*Np() + l2 ];
}

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_H_xx /* arg */
  ) const
{
  return supports_H_xx_;
}

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_H_xp /* arg */, int l
  ) const
{
  assert_l(l);
  return supports_H_xp_[l];
}

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_H_pp /* arg */, int l1, int l2
  ) const
{
  assert_l(l1);
  assert_l(l2);
  return supports_H_pp_[ l1*Np() + l2 ];
}

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_hess_g_xx /* arg */, int j
  ) const
{
  assert_j(j);
  return supports_hess_g_xx_[j];
}

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_hess_g_xp /* arg */, int j, int l
  ) const
{
  assert_j(j);
  assert_l(l);
  return supports_hess_g_xp_[ j*Np() + l ];
}

template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_hess_g_pp /* arg */, int j, int l1, int l2
  ) const
{
  assert_j(j);
  assert_l(l1);
  assert_l(l2);
  return supports_hess_g_pp_[ j*Np()*Np() + l1*Np() + l2 ];
}

#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT


template<class Scalar>
bool
ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgs_g_mp /* arg */, int j
  ) const
{
  assert_j(j);
  return supports_g_mp_[j];
}

template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgsDfDp_mp /* arg */, int l
  ) const
{
  assert_l(l);
  return supports_DfDp_mp_[l];
}


template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgsDgDx_dot_mp /* arg */, int j
  ) const
{
  assert_j(j);
  return supports_DgDx_dot_mp_[j];
}


template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgsDgDx_mp /* arg */, int j
  ) const
{
  assert_j(j);
  return supports_DgDx_mp_[j];
}


template<class Scalar>
const ModelEvaluatorBase::DerivativeSupport&
ModelEvaluatorBase::OutArgs<Scalar>::supports(
  EOutArgsDgDp_mp /* arg */, int j, int l
  ) const
{
  assert_j(j);
  assert_l(l);
  return supports_DgDp_mp_[ j*Np() + l ];
}



template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_f(
  const Evaluation<VectorBase<Scalar> > &f
  )
{
  assert_supports(OUT_ARG_f);
  f_ = f;
}


template<class Scalar>
ModelEvaluatorBase::Evaluation<VectorBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_f() const
{
  assert_supports(OUT_ARG_f);
  return f_;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_g(
  int j, const Evaluation<VectorBase<Scalar> > &g_j
  )
{
  assert_j(j);
  g_[j] = g_j;
}


template<class Scalar>
ModelEvaluatorBase::Evaluation<VectorBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_g(int j) const
{
  assert_j(j);
  return g_[j];
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_f_mp(
  const RCP<Stokhos::ProductEpetraVector > &f_mp
  )
{
  assert_supports(OUT_ARG_f_mp);
  f_mp_ = f_mp;
}


template<class Scalar>
RCP<Stokhos::ProductEpetraVector >
ModelEvaluatorBase::OutArgs<Scalar>::get_f_mp() const
{
  assert_supports(OUT_ARG_f_mp);
  return f_mp_;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_g_mp(
  int j, const RCP<Stokhos::ProductEpetraVector> &g_mp_j
  )
{
  assert_supports(OUT_ARG_g_mp,j);
  g_mp_[j] = g_mp_j;
}


template<class Scalar>
RCP<Stokhos::ProductEpetraVector>
ModelEvaluatorBase::OutArgs<Scalar>::get_g_mp(int j) const
{
  assert_supports(OUT_ARG_g_mp,j);
  return g_mp_[j];
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_W(
  const RCP<LinearOpWithSolveBase<Scalar> > &W
  )
{
  assert_supports(OUT_ARG_W);
  W_ = W;
}


template<class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_W() const
{
  assert_supports(OUT_ARG_W);
  return W_;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_W_mp(
  const RCP<Stokhos::ProductEpetraOperator> &W_mp
  )
{
  assert_supports(OUT_ARG_W_mp);
  W_mp_ = W_mp;
}


template<class Scalar>
RCP<Stokhos::ProductEpetraOperator>
ModelEvaluatorBase::OutArgs<Scalar>::get_W_mp() const
{
  assert_supports(OUT_ARG_W_mp);
  return W_mp_;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_W_op(
  const RCP<LinearOpBase<Scalar> > &W_op
  )
{
  assert_supports(OUT_ARG_W_op);
  W_op_ = W_op;
}


template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_W_op() const
{
  assert_supports(OUT_ARG_W_op);
  return W_op_;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_W_prec(
  const RCP<PreconditionerBase<Scalar> > &W_prec
  )
{
  assert_supports(OUT_ARG_W_prec);
  W_prec_ = W_prec;
}


template<class Scalar>
RCP<PreconditionerBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_W_prec() const
{
  assert_supports(OUT_ARG_W_prec);
  return W_prec_;
}


template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_W_properties() const
{
  assert_supports(OUT_ARG_f);
  return W_properties_;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_DfDp(
  int l, const Derivative<Scalar> &DfDp_l
  )
{
  assert_supports(OUT_ARG_DfDp,l,DfDp_l);
  DfDp_[l] = DfDp_l;
}


template<class Scalar>
ModelEvaluatorBase::Derivative<Scalar>
ModelEvaluatorBase::OutArgs<Scalar>::get_DfDp(int l) const
{
  assert_supports(OUT_ARG_DfDp,l);
  return DfDp_[l];
}


template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_DfDp_properties(int l) const
{
  assert_supports(OUT_ARG_DfDp,l);
  return DfDp_properties_[l];
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_DfDp_mp(
  int l, const MPDerivative &DfDp_mp_l
  )
{
  assert_supports(OUT_ARG_DfDp_mp,l,DfDp_mp_l);
  DfDp_mp_[l] = DfDp_mp_l;
}


template<class Scalar>
ModelEvaluatorBase::MPDerivative
ModelEvaluatorBase::OutArgs<Scalar>::get_DfDp_mp(int l) const
{
  assert_supports(OUT_ARG_DfDp_mp,l);
  return DfDp_mp_[l];
}


template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_DfDp_mp_properties(int l) const
{
  assert_supports(OUT_ARG_DfDp_mp,l);
  return DfDp_mp_properties_[l];
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_DgDx_dot(
  int j, const Derivative<Scalar> &DgDx_dot_j
  )
{
  assert_supports(OUT_ARG_DgDx_dot,j,DgDx_dot_j);
  DgDx_dot_[j] = DgDx_dot_j;
}


template<class Scalar>
ModelEvaluatorBase::Derivative<Scalar>
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDx_dot(int j) const
{
  assert_supports(OUT_ARG_DgDx_dot,j);
  return DgDx_dot_[j];
}


template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDx_dot_properties(int j) const
{
  assert_supports(OUT_ARG_DgDx_dot,j);
  return DgDx_dot_properties_[j];
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_DgDx_dot_mp(
  int j, const MPDerivative &DgDx_dot_mp_j
  )
{
  assert_supports(OUT_ARG_DgDx_dot_mp,j,DgDx_dot_mp_j);
  DgDx_dot_mp_[j] = DgDx_dot_mp_j;
}


template<class Scalar>
ModelEvaluatorBase::MPDerivative
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDx_dot_mp(int j) const
{
  assert_supports(OUT_ARG_DgDx_dot_mp,j);
  return DgDx_dot_mp_[j];
}


template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDx_dot_mp_properties(int j) const
{
  assert_supports(OUT_ARG_DgDx_dot_mp,j);
  return DgDx_dot_mp_properties_[j];
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_DgDx(
  int j, const Derivative<Scalar> &DgDx_j
  )
{
  assert_supports(OUT_ARG_DgDx,j,DgDx_j);
  DgDx_[j] = DgDx_j;
}


template<class Scalar>
ModelEvaluatorBase::Derivative<Scalar>
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDx(int j) const
{
  assert_supports(OUT_ARG_DgDx,j);
  return DgDx_[j];
}


template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDx_properties(int j) const
{
  assert_supports(OUT_ARG_DgDx,j);
  return DgDx_properties_[j];
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_DgDx_mp(
  int j, const MPDerivative &DgDx_mp_j
  )
{
  assert_supports(OUT_ARG_DgDx_mp,j,DgDx_mp_j);
  DgDx_mp_[j] = DgDx_mp_j;
}


template<class Scalar>
ModelEvaluatorBase::MPDerivative
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDx_mp(int j) const
{
  assert_supports(OUT_ARG_DgDx_mp,j);
  return DgDx_mp_[j];
}


template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDx_mp_properties(int j) const
{
  assert_supports(OUT_ARG_DgDx_mp,j);
  return DgDx_mp_properties_[j];
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_DgDp(
  int j, int l, const Derivative<Scalar> &DgDp_j_l
  )
{
  assert_supports(OUT_ARG_DgDp,j,l,DgDp_j_l);
  DgDp_[ j*Np() + l ] = DgDp_j_l;
}


template<class Scalar>
ModelEvaluatorBase::Derivative<Scalar>
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDp(int j, int l) const
{
  assert_supports(OUT_ARG_DgDp,j,l);
  return DgDp_[ j*Np() + l ];
}


template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDp_properties(int j, int l) const
{
  assert_supports(OUT_ARG_DgDp,j,l);
  return DgDp_properties_[ j*Np() + l ];
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_DgDp_mp(
  int j, int l, const MPDerivative &DgDp_mp_j_l
  )
{
  assert_supports(OUT_ARG_DgDp_mp,j,l,DgDp_mp_j_l);
  DgDp_mp_[ j*Np() + l ] = DgDp_mp_j_l;
}


template<class Scalar>
ModelEvaluatorBase::MPDerivative
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDp_mp(int j, int l) const
{
  assert_supports(OUT_ARG_DgDp_mp,j,l);
  return DgDp_mp_[ j*Np() + l ];
}


template<class Scalar>
ModelEvaluatorBase::DerivativeProperties
ModelEvaluatorBase::OutArgs<Scalar>::get_DgDp_mp_properties(int j, int l) const
{
  assert_supports(OUT_ARG_DgDp_mp,j,l);
  return DgDp_mp_properties_[ j*Np() + l ];
}


#ifdef Thyra_BUILD_HESSIAN_SUPPORT

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_hess_vec_prod_f_xx(
  const RCP<MultiVectorBase<Scalar> > &hess_vec_prod_f_xx
  )
{
  assert_supports(OUT_ARG_hess_vec_prod_f_xx);
  hess_vec_prod_f_xx_ = hess_vec_prod_f_xx;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_hess_vec_prod_f_xp(
  int l, const RCP<MultiVectorBase<Scalar> > &hess_vec_prod_f_xp_l
  )
{
  assert_supports(OUT_ARG_hess_vec_prod_f_xp,l);
  hess_vec_prod_f_xp_[l] = hess_vec_prod_f_xp_l;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_hess_vec_prod_f_px(
  int l, const RCP<MultiVectorBase<Scalar> > &hess_vec_prod_f_px_l
  )
{
  assert_supports(OUT_ARG_hess_vec_prod_f_px,l);
  hess_vec_prod_f_px_[l] = hess_vec_prod_f_px_l;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_hess_vec_prod_f_pp(
  int l1, int l2, const RCP<MultiVectorBase<Scalar> > &hess_vec_prod_f_pp_l1_l2
  )
{
  assert_supports(OUT_ARG_hess_vec_prod_f_pp,l1,l2);
  hess_vec_prod_f_pp_[ l1*Np() + l2 ] = hess_vec_prod_f_pp_l1_l2;
}

template<class Scalar>
RCP<MultiVectorBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_hess_vec_prod_f_xx() const
{
  assert_supports(OUT_ARG_hess_vec_prod_f_xx);
  return hess_vec_prod_f_xx_;
}

template<class Scalar>
RCP<MultiVectorBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_hess_vec_prod_f_xp(int l) const
{
  assert_supports(OUT_ARG_hess_vec_prod_f_xp,l);
  return hess_vec_prod_f_xp_[l];
}

template<class Scalar>
RCP<MultiVectorBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_hess_vec_prod_f_px(int l) const
{
  assert_supports(OUT_ARG_hess_vec_prod_f_px,l);
  return hess_vec_prod_f_px_[l];
}

template<class Scalar>
RCP<MultiVectorBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_hess_vec_prod_f_pp(int l1, int l2) const
{
  assert_supports(OUT_ARG_hess_vec_prod_f_pp,l1,l2);
  return hess_vec_prod_f_pp_[ l1*Np() + l2 ];
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_hess_vec_prod_g_xx(
  int j, const RCP<MultiVectorBase<Scalar> > &hess_vec_prod_g_xx_j
  )
{
  assert_supports(OUT_ARG_hess_vec_prod_g_xx,j);
  hess_vec_prod_g_xx_[j] = hess_vec_prod_g_xx_j;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_hess_vec_prod_g_xp(
  int j, int l, const RCP<MultiVectorBase<Scalar> > &hess_vec_prod_g_xp_j_l
  )
{
  assert_supports(OUT_ARG_hess_vec_prod_g_xp,j,l);
  hess_vec_prod_g_xp_[ j*Np() + l ] = hess_vec_prod_g_xp_j_l;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_hess_vec_prod_g_px(
  int j, int l, const RCP<MultiVectorBase<Scalar> > &hess_vec_prod_g_px_j_l
  )
{
  assert_supports(OUT_ARG_hess_vec_prod_g_px,j,l);
  hess_vec_prod_g_px_[ j*Np() + l ] = hess_vec_prod_g_px_j_l;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_hess_vec_prod_g_pp(
  int j, int l1, int l2, const RCP<MultiVectorBase<Scalar> > &hess_vec_prod_g_pp_j_l1_l2
  )
{
  assert_supports(OUT_ARG_hess_vec_prod_g_pp,j,l1,l2);
  hess_vec_prod_g_pp_[ j*Np()*Np() + l1*Np() + l2 ] = hess_vec_prod_g_pp_j_l1_l2;
}

template<class Scalar>
RCP<MultiVectorBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_hess_vec_prod_g_xx(int j) const
{
  assert_supports(OUT_ARG_hess_vec_prod_g_xx,j);
  return hess_vec_prod_g_xx_[j];
}

template<class Scalar>
RCP<MultiVectorBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_hess_vec_prod_g_xp(int j, int l) const
{
  assert_supports(OUT_ARG_hess_vec_prod_g_xp,j,l);
  return hess_vec_prod_g_xp_[ j*Np() + l ];
}

template<class Scalar>
RCP<MultiVectorBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_hess_vec_prod_g_px(int j, int l) const
{
  assert_supports(OUT_ARG_hess_vec_prod_g_px,j,l);
  return hess_vec_prod_g_px_[ j*Np() + l ];
}

template<class Scalar>
RCP<MultiVectorBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_hess_vec_prod_g_pp(int j, int l1, int l2) const
{
  assert_supports(OUT_ARG_hess_vec_prod_g_pp,j,l1,l2);
  return hess_vec_prod_g_pp_[ j*Np()*Np() + l1*Np() + l2 ];
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_hess_f_xx(
  const RCP<LinearOpBase<Scalar> > &hess_f_xx
  )
{
  assert_supports(OUT_ARG_hess_f_xx);
  hess_f_xx_ = hess_f_xx;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_hess_f_xp(
  int l, const RCP<LinearOpBase<Scalar> > &hess_f_xp_l
  )
{
  assert_supports(OUT_ARG_hess_f_xp,l);
  hess_f_xp_[l] = hess_f_xp_l;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_hess_f_pp(
  int l1, int l2, const RCP<LinearOpBase<Scalar> > &hess_f_pp_l1_l2
  )
{
  assert_supports(OUT_ARG_hess_f_pp,l1,l2);
  hess_f_pp_[ l1*Np() + l2 ] = hess_f_pp_l1_l2;
}

template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_hess_f_xx() const
{
  assert_supports(OUT_ARG_hess_f_xx);
  return hess_f_xx_;
}

template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_hess_f_xp(int l) const
{
  assert_supports(OUT_ARG_hess_f_xp,l);
  return hess_f_xp_[l];
}

template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_hess_f_pp(int l1, int l2) const
{
  assert_supports(OUT_ARG_hess_f_pp,l1,l2);
  return hess_f_pp_[ l1*Np() + l2 ];
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_hess_g_xx(
  int j, const RCP<LinearOpBase<Scalar> > &hess_g_xx_j
  )
{
  assert_supports(OUT_ARG_hess_g_xx,j);
  hess_g_xx_[j] = hess_g_xx_j;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_hess_g_xp(
  int j, int l, const RCP<LinearOpBase<Scalar> > &hess_g_xp_j_l
  )
{
  assert_supports(OUT_ARG_hess_g_xp,j,l);
  hess_g_xp_[ j*Np() + l ] = hess_g_xp_j_l;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_hess_g_pp(
  int j, int l1, int l2, const RCP<LinearOpBase<Scalar> > &hess_g_pp_j_l1_l2
  )
{
  assert_supports(OUT_ARG_hess_g_pp,j,l1,l2);
  hess_g_pp_[ j*Np()*Np() + l1*Np() + l2 ] = hess_g_pp_j_l1_l2;
}

template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_hess_g_xx(int j) const
{
  assert_supports(OUT_ARG_hess_g_xx,j);
  return hess_g_xx_[j];
}

template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_hess_g_xp(int j, int l) const
{
  assert_supports(OUT_ARG_hess_g_xp,j,l);
  return hess_g_xp_[ j*Np() + l ];
}

template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_hess_g_pp(int j, int l1, int l2) const
{
  assert_supports(OUT_ARG_hess_g_pp,j,l1,l2);
  return hess_g_pp_[ j*Np()*Np() + l1*Np() + l2 ];
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_H_xx(
  const RCP<LinearOpBase<Scalar> > &H_xx
  )
{
  assert_supports(OUT_ARG_H_xx);
  H_xx_ = H_xx;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_H_xp(
  int l, const RCP<LinearOpBase<Scalar> > &H_xp_l
  )
{
  assert_supports(OUT_ARG_H_xp,l);
  H_xp_[l] = H_xp_l;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_H_pp(
  int l1, int l2, const RCP<LinearOpBase<Scalar> > &H_pp_l1_l2
  )
{
  assert_supports(OUT_ARG_H_pp,l1,l2);
  H_pp_[ l1*Np() + l2 ] = H_pp_l1_l2;
}

template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_H_xx() const
{
  assert_supports(OUT_ARG_H_xx);
  return H_xx_;
}

template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_H_xp(int l) const
{
  assert_supports(OUT_ARG_H_xp,l);
  return H_xp_[l];
}

template<class Scalar>
RCP<LinearOpBase<Scalar> >
ModelEvaluatorBase::OutArgs<Scalar>::get_H_pp(int l1, int l2) const
{
  assert_supports(OUT_ARG_H_pp,l1,l2);
  return H_pp_[ l1*Np() + l2 ];
}

#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT


#ifdef HAVE_THYRA_ME_POLYNOMIAL


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::set_f_poly(
  const RCP<Teuchos::Polynomial< VectorBase<Scalar> > > &f_poly
  )
{
  f_poly_ = f_poly;
}


template<class Scalar>
RCP<Teuchos::Polynomial< VectorBase<Scalar> > >
ModelEvaluatorBase::OutArgs<Scalar>::get_f_poly() const
{
  return f_poly_;
}


#endif // HAVE_THYRA_ME_POLYNOMIAL


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::setArgs(
  const OutArgs<Scalar>& inputOutArgs, bool ignoreUnsupported
  )
{
  typedef ModelEvaluatorBase MEB;
  const int min_Np = TEUCHOS_MIN(this->Np(),inputOutArgs.Np());
  const int min_Ng = TEUCHOS_MIN(this->Ng(),inputOutArgs.Ng());
  // f
  if ( inputOutArgs.supports(OUT_ARG_f) && nonnull(inputOutArgs.get_f()) ) {
    if ( supports(OUT_ARG_f) || !ignoreUnsupported )
      set_f(inputOutArgs.get_f());
  }
  if ( inputOutArgs.supports(OUT_ARG_f_mp) && nonnull(inputOutArgs.get_f_mp()) ) {
    if ( supports(OUT_ARG_f_mp) || !ignoreUnsupported )
      set_f_mp(inputOutArgs.get_f_mp());
  }
#ifdef HAVE_THYRA_ME_POLYNOMIAL
  // f_poly
  if ( inputOutArgs.supports(OUT_ARG_f_poly) && nonnull(inputOutArgs.get_f_poly()) ) {
    if ( supports(OUT_ARG_f_poly) || !ignoreUnsupported )
      set_f_poly(inputOutArgs.get_f_poly());
  }
#endif // HAVE_THYRA_ME_POLYNOMIAL
  // g(j)
  for ( int j = 0; j < min_Ng; ++j ) {
    if ( nonnull(inputOutArgs.get_g(j)) )
      set_g(j,inputOutArgs.get_g(j));
  }
  for ( int j = 0; j < min_Ng; ++j ) {
    if ( inputOutArgs.supports(OUT_ARG_g_mp,j) && nonnull(inputOutArgs.get_g_mp(j)) ) {
    if ( supports(OUT_ARG_g_mp,j) || !ignoreUnsupported )
      set_g_mp(j,inputOutArgs.get_g_mp(j));
    }
  }
  // W
  if( inputOutArgs.supports(OUT_ARG_W) && nonnull(inputOutArgs.get_W()) ) {
    if ( supports(OUT_ARG_W) || !ignoreUnsupported )
      set_W(inputOutArgs.get_W());
  }
  if( inputOutArgs.supports(OUT_ARG_W_mp) && nonnull(inputOutArgs.get_W_mp()) ) {
    if ( supports(OUT_ARG_W_mp) || !ignoreUnsupported )
      set_W_mp(inputOutArgs.get_W_mp());
  }
  // W_op
  if( inputOutArgs.supports(OUT_ARG_W_op) && nonnull(inputOutArgs.get_W_op()) ) {
    if ( supports(OUT_ARG_W_op) || !ignoreUnsupported )
      set_W_op(inputOutArgs.get_W_op());
  }
  // W_prec
  if( inputOutArgs.supports(OUT_ARG_W_prec) && nonnull(inputOutArgs.get_W_prec()) ) {
    if ( supports(OUT_ARG_W_prec) || !ignoreUnsupported )
      set_W_prec(inputOutArgs.get_W_prec());
  }
  // DfDp(l)
  for ( int l = 0; l < min_Np; ++l ) {
    MEB::Derivative<Scalar> DfDp_l;
    if ( !inputOutArgs.supports(OUT_ARG_DfDp,l).none()
      && !(DfDp_l=inputOutArgs.get_DfDp(l)).isEmpty() )
    {
      if ( DfDp_l.isSupportedBy(supports(OUT_ARG_DfDp,l)) || !ignoreUnsupported )
        set_DfDp(l,DfDp_l);
    }
  }
  for ( int l = 0; l < min_Np; ++l ) {
    MEB::MPDerivative DfDp_mp_l;
    if ( !inputOutArgs.supports(OUT_ARG_DfDp_mp,l).none()
      && !(DfDp_mp_l=inputOutArgs.get_DfDp_mp(l)).isEmpty() )
    {
      if ( DfDp_mp_l.isSupportedBy(supports(OUT_ARG_DfDp_mp,l)) || !ignoreUnsupported )
        set_DfDp_mp(l,DfDp_mp_l);
    }
  }
  // DgDx_dot(j) and DgDx(j)
  for ( int j = 0; j < min_Ng; ++j ) {
    // DgDx_dot(j)
    MEB::Derivative<Scalar> DgDx_dot_j;
    if ( !inputOutArgs.supports(OUT_ARG_DgDx_dot,j).none()
      && !(DgDx_dot_j=inputOutArgs.get_DgDx_dot(j)).isEmpty() )
    {
      if( DgDx_dot_j.isSupportedBy(supports(OUT_ARG_DgDx_dot,j)) || !ignoreUnsupported )
        set_DgDx_dot(j,DgDx_dot_j);
    }
    // DgDx(j)
    MEB::Derivative<Scalar> DgDx_j;
    if ( !inputOutArgs.supports(OUT_ARG_DgDx,j).none()
      && !(DgDx_j=inputOutArgs.get_DgDx(j)).isEmpty() ) {
      if ( DgDx_j.isSupportedBy(supports(OUT_ARG_DgDx,j)) || !ignoreUnsupported )
        set_DgDx(j,DgDx_j);
    }
  }
  for ( int j = 0; j < min_Ng; ++j ) {
    // DgDx_dot(j)
    MEB::MPDerivative DgDx_dot_mp_j;
    if ( !inputOutArgs.supports(OUT_ARG_DgDx_dot_mp,j).none()
      && !(DgDx_dot_mp_j=inputOutArgs.get_DgDx_dot_mp(j)).isEmpty() )
    {
      if( DgDx_dot_mp_j.isSupportedBy(supports(OUT_ARG_DgDx_dot_mp,j)) || !ignoreUnsupported )
        set_DgDx_dot_mp(j,DgDx_dot_mp_j);
    }
    // DgDx(j)
    MEB::MPDerivative DgDx_mp_j;
    if ( !inputOutArgs.supports(OUT_ARG_DgDx_mp,j).none()
      && !(DgDx_mp_j=inputOutArgs.get_DgDx_mp(j)).isEmpty() ) {
      if ( DgDx_mp_j.isSupportedBy(supports(OUT_ARG_DgDx_mp,j)) || !ignoreUnsupported )
        set_DgDx_mp(j,DgDx_mp_j);
    }
  }
  // DgDp(j,l)
  for ( int l = 0; l < min_Np; ++l ) {
    for ( int j = 0; j < min_Ng; ++j ) {
      MEB::Derivative<Scalar> DgDp_j_l;
      if ( !inputOutArgs.supports(OUT_ARG_DgDp,j,l).none()
        && !(DgDp_j_l=inputOutArgs.get_DgDp(j,l)).isEmpty() )
      {
        if ( DgDp_j_l.isSupportedBy(supports(OUT_ARG_DgDp,j,l)) || !ignoreUnsupported )
          set_DgDp(j,l,DgDp_j_l);
      }
    }
  }
  for ( int l = 0; l < min_Np; ++l ) {
    for ( int j = 0; j < min_Ng; ++j ) {
      MEB::MPDerivative DgDp_mp_j_l;
      if ( !inputOutArgs.supports(OUT_ARG_DgDp_mp,j,l).none()
        && !(DgDp_mp_j_l=inputOutArgs.get_DgDp_mp(j,l)).isEmpty() )
      {
        if ( DgDp_mp_j_l.isSupportedBy(supports(OUT_ARG_DgDp_mp,j,l)) || !ignoreUnsupported )
          set_DgDp_mp(j,l,DgDp_mp_j_l);
      }
    }
  }

#ifdef Thyra_BUILD_HESSIAN_SUPPORT

  // hess_vec_prod_f
  if( inputOutArgs.supports(OUT_ARG_hess_vec_prod_f_xx) && nonnull(inputOutArgs.get_hess_vec_prod_f_xx()) ) {
    if ( supports(OUT_ARG_hess_vec_prod_f_xx) || !ignoreUnsupported )
      set_hess_vec_prod_f_xx(inputOutArgs.get_hess_vec_prod_f_xx());
  }
  for ( int l1 = 0; l1 < min_Np; ++l1 ) {
    if( inputOutArgs.supports(OUT_ARG_hess_vec_prod_f_xp,l1) && nonnull(inputOutArgs.get_hess_vec_prod_f_xp(l1)) ) {
      if ( supports(OUT_ARG_hess_vec_prod_f_xp,l1) || !ignoreUnsupported )
        set_hess_vec_prod_f_xp(l1,inputOutArgs.get_hess_vec_prod_f_xp(l1));
    }
    if( inputOutArgs.supports(OUT_ARG_hess_vec_prod_f_px,l1) && nonnull(inputOutArgs.get_hess_vec_prod_f_px(l1)) ) {
      if ( supports(OUT_ARG_hess_vec_prod_f_px,l1) || !ignoreUnsupported )
        set_hess_vec_prod_f_px(l1,inputOutArgs.get_hess_vec_prod_f_px(l1));
    }
    for ( int l2 = 0; l2 < min_Np; ++l2 ) {
      if( inputOutArgs.supports(OUT_ARG_hess_vec_prod_f_pp,l1,l2) && nonnull(inputOutArgs.get_hess_vec_prod_f_pp(l1,l2)) ) {
        if ( supports(OUT_ARG_hess_vec_prod_f_pp,l1,l2) || !ignoreUnsupported )
          set_hess_vec_prod_f_pp(l1,l2,inputOutArgs.get_hess_vec_prod_f_pp(l1,l2));
      }
    }
  }

  // hess_vec_prod_g
  for ( int j = 0; j < min_Ng; ++j ) {
    if( inputOutArgs.supports(OUT_ARG_hess_vec_prod_g_xx,j) && nonnull(inputOutArgs.get_hess_vec_prod_g_xx(j)) ) {
      if ( supports(OUT_ARG_hess_vec_prod_g_xx,j) || !ignoreUnsupported )
        set_hess_vec_prod_g_xx(j,inputOutArgs.get_hess_vec_prod_g_xx(j));
    }
    for ( int l1 = 0; l1 < min_Np; ++l1 ) {
      if( inputOutArgs.supports(OUT_ARG_hess_vec_prod_g_xp,j,l1) && nonnull(inputOutArgs.get_hess_vec_prod_g_xp(j,l1)) ) {
        if ( supports(OUT_ARG_hess_vec_prod_g_xp,j,l1) || !ignoreUnsupported )
          set_hess_vec_prod_g_xp(j,l1,inputOutArgs.get_hess_vec_prod_g_xp(j,l1));
      }
      if( inputOutArgs.supports(OUT_ARG_hess_vec_prod_g_px,j,l1) && nonnull(inputOutArgs.get_hess_vec_prod_g_px(j,l1)) ) {
        if ( supports(OUT_ARG_hess_vec_prod_g_px,j,l1) || !ignoreUnsupported )
          set_hess_vec_prod_g_px(j,l1,inputOutArgs.get_hess_vec_prod_g_px(j,l1));
      }
      for ( int l2 = 0; l2 < min_Np; ++l2 ) {
        if( inputOutArgs.supports(OUT_ARG_hess_vec_prod_g_pp,j,l1,l2) && nonnull(inputOutArgs.get_hess_vec_prod_g_pp(j,l1,l2)) ) {
          if ( supports(OUT_ARG_hess_vec_prod_g_pp,j,l1,l2) || !ignoreUnsupported )
            set_hess_vec_prod_g_pp(j,l1,l2,inputOutArgs.get_hess_vec_prod_g_pp(j,l1,l2));
        }
      }
    }
  }

  // hess_f
  if( inputOutArgs.supports(OUT_ARG_hess_f_xx) && nonnull(inputOutArgs.get_hess_f_xx()) ) {
    if ( supports(OUT_ARG_hess_f_xx) || !ignoreUnsupported )
      set_hess_f_xx(inputOutArgs.get_hess_f_xx());
  }
  for ( int l1 = 0; l1 < min_Np; ++l1 ) {
    if( inputOutArgs.supports(OUT_ARG_hess_f_xp,l1) && nonnull(inputOutArgs.get_hess_f_xp(l1)) ) {
      if ( supports(OUT_ARG_hess_f_xp,l1) || !ignoreUnsupported )
        set_hess_f_xp(l1,inputOutArgs.get_hess_f_xp(l1));
    }
    for ( int l2 = 0; l2 < min_Np; ++l2 ) {
      if( inputOutArgs.supports(OUT_ARG_hess_f_pp,l1,l2) && nonnull(inputOutArgs.get_hess_f_pp(l1,l2)) ) {
        if ( supports(OUT_ARG_hess_f_pp,l1,l2) || !ignoreUnsupported )
          set_hess_f_pp(l1,l2,inputOutArgs.get_hess_f_pp(l1,l2));
      }
    }
  }

  // hess_g
  for ( int j = 0; j < min_Ng; ++j ) {
    if( inputOutArgs.supports(OUT_ARG_hess_g_xx,j) && nonnull(inputOutArgs.get_hess_g_xx(j)) ) {
      if ( supports(OUT_ARG_hess_g_xx,j) || !ignoreUnsupported )
        set_hess_g_xx(j,inputOutArgs.get_hess_g_xx(j));
    }
    for ( int l1 = 0; l1 < min_Np; ++l1 ) {
      if( inputOutArgs.supports(OUT_ARG_hess_g_xp,j,l1) && nonnull(inputOutArgs.get_hess_g_xp(j,l1)) ) {
        if ( supports(OUT_ARG_hess_g_xp,j,l1) || !ignoreUnsupported )
          set_hess_g_xp(j,l1,inputOutArgs.get_hess_g_xp(j,l1));
      }
      for ( int l2 = 0; l2 < min_Np; ++l2 ) {
        if( inputOutArgs.supports(OUT_ARG_hess_g_pp,j,l1,l2) && nonnull(inputOutArgs.get_hess_g_pp(j,l1,l2)) ) {
          if ( supports(OUT_ARG_hess_g_pp,j,l1,l2) || !ignoreUnsupported )
            set_hess_g_pp(j,l1,l2,inputOutArgs.get_hess_g_pp(j,l1,l2));
        }
      }
    }
  }

  // H
  if( inputOutArgs.supports(OUT_ARG_H_xx) && nonnull(inputOutArgs.get_H_xx()) ) {
    if ( supports(OUT_ARG_H_xx) || !ignoreUnsupported )
      set_H_xx(inputOutArgs.get_H_xx());
  }
  for ( int l1 = 0; l1 < min_Np; ++l1 ) {
    if( inputOutArgs.supports(OUT_ARG_H_xp,l1) && nonnull(inputOutArgs.get_H_xp(l1)) ) {
      if ( supports(OUT_ARG_H_xp,l1) || !ignoreUnsupported )
        set_H_xp(l1,inputOutArgs.get_H_xp(l1));
    }
    for ( int l2 = 0; l2 < min_Np; ++l2 ) {
      if( inputOutArgs.supports(OUT_ARG_H_pp,l1,l2) && nonnull(inputOutArgs.get_H_pp(l1,l2)) ) {
        if ( supports(OUT_ARG_H_pp,l1,l2) || !ignoreUnsupported )
          set_H_pp(l1,l2,inputOutArgs.get_H_pp(l1,l2));
      }
    }
  }

#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

  // Extended outArgs
  this->extended_outargs_ = inputOutArgs.extended_outargs_;

  // ToDo: Add more args as needed!
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::setFailed() const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  isFailed_ = true;
  if( this->supports(OUT_ARG_f) && nonnull(this->get_f()) ) {
    assign(this->get_f().ptr(),ST::nan());
  }
  for( int j = 0; j < this->Ng(); ++j ) {
    if (nonnull(this->get_g(j)))
      assign(this->get_g(j).ptr(),ST::nan());
  }
  // ToDo: Set other objects to NaN as well!
}


template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::isFailed() const
{
  return isFailed_;
}


template<class Scalar>
bool ModelEvaluatorBase::OutArgs<Scalar>::isEmpty() const
{
  if (!is_null(f_))
    return false;
  if (!is_null(W_))
    return false;
  if (!is_null(W_op_))
    return false;
  for ( int l = 0; l < Np(); ++l ) {
    if (!DfDp_[l].isEmpty())
      return false;
  }
#ifdef HAVE_THYRA_ME_POLYNOMIAL
  if (!is_null(f_poly_))
    return false;
#endif // HAVE_THYRA_ME_POLYNOMIAL
  for ( int j = 0; j < Ng(); ++j ) {
    if (!is_null(g_[j]))
      return false;
    if (!DgDx_dot_[j].isEmpty())
      return false;
    if (!DgDx_[j].isEmpty())
      return false;
    for ( int l = 0; l < Np(); ++l ) {
      if (!DgDp_[j*Np()+l].isEmpty())
        return false;
    }
  }
  return true;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assertSameSupport(
  const OutArgs<Scalar> &outArgs
  ) const
{

  for ( int outArg_i = 0; outArg_i < NUM_E_OUT_ARGS_MEMBERS; ++outArg_i ) {
    const EOutArgsMembers outArg_arg = static_cast<EOutArgsMembers>(outArg_i);
    const std::string outArg_name = toString(outArg_arg);
    TEUCHOS_TEST_FOR_EXCEPTION(
      supports(outArg_arg) != outArgs.supports(outArg_arg), std::logic_error,
      "Error, the output argument "<<outArg_name<<" with support "<<outArgs.supports(outArg_arg)<<"\n"
      "in the OutArgs object for the model:\n\n"
      "  "<<outArgs.modelEvalDescription()<<"\n\n"
      "is not the same the argument "<<outArg_name<<" with support "<<supports(outArg_arg)<<"\n"
      "in the OutArgs object for the model:\n\n"
      "  "<<modelEvalDescription()<<"\n\n"
      "and these two OutArgs objects are not compatible!"
      );
  }

  const int l_Np = this->Np();
  const int l_Ng = this->Ng();
  TEUCHOS_ASSERT_EQUALITY( l_Np, outArgs.Np() );
  TEUCHOS_ASSERT_EQUALITY( l_Ng, outArgs.Ng() );

  if (supports(OUT_ARG_f)) {
    for ( int l = 0; l < l_Np; ++l ) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        !supports(OUT_ARG_DfDp,l).isSameSupport(outArgs.supports(OUT_ARG_DfDp,l)),
        std::logic_error,
        "Error, the support for DfDp("<<l<<") is not the same for the models\n\n"
        "  "<<outArgs.modelEvalDescription()<<"\n\n"
        "and:\n\n"
        "  "<<modelEvalDescription()<<"\n\n"
        "and these two OutArgs objects are not compatible!"
        );
    }
  }

  for ( int j = 0; j < l_Ng; ++j ) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      !supports(OUT_ARG_DgDx_dot,j).isSameSupport(outArgs.supports(OUT_ARG_DgDx_dot,j)),
      std::logic_error,
      "Error, the support for DgDx_dot("<<j<<") is not the same for the models\n\n"
      "  "<<outArgs.modelEvalDescription()<<"\n\n"
      "and:\n\n"
      "  "<<modelEvalDescription()<<"\n\n"
      "and these two OutArgs objects are not compatible!"
      );
    TEUCHOS_TEST_FOR_EXCEPTION(
      !supports(OUT_ARG_DgDx,j).isSameSupport(outArgs.supports(OUT_ARG_DgDx,j)),
      std::logic_error,
      "Error, the support for DgDx("<<j<<") is not the same for the models\n\n"
      "  "<<outArgs.modelEvalDescription()<<"\n\n"
      "and:\n\n"
      "  "<<modelEvalDescription()<<"\n\n"
      "and these two OutArgs objects are not compatible!"
      );
    for ( int l = 0; l < l_Np; ++l ) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        !supports(OUT_ARG_DgDp,j,l).isSameSupport(outArgs.supports(OUT_ARG_DgDp,j,l)),
        std::logic_error,
        "Error, the support for DgDp("<<j<<","<<l<<") is not the same for the models\n\n"
        "  "<<outArgs.modelEvalDescription()<<"\n\n"
        "and:\n\n"
        "  "<<modelEvalDescription()<<"\n\n"
        "and these two OutArgs objects are not compatible!"
        );
    }
  }
}


template<class Scalar>
std::string ModelEvaluatorBase::OutArgs<Scalar>::modelEvalDescription() const
{
  return modelEvalDescription_;
}


template<class Scalar>
std::string ModelEvaluatorBase::OutArgs<Scalar>::description() const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  std::ostringstream oss;
  oss
    << "Thyra::ModelEvaluatorBase::OutArgs<"<<ST::name()<<">"
    << "{"
    << "model="<<modelEvalDescription_
    << ",Np="<<Np()
    << ",Ng="<<Ng()
    << "}";
  return oss.str();
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::describe(
  Teuchos::FancyOStream &out_arg, const Teuchos::EVerbosityLevel verbLevel
  ) const
{

  using Teuchos::OSTab;
  using Teuchos::describe;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef RCP<const VectorBase<Scalar> > CV_ptr;
  typedef RCP<const LinearOpBase<Scalar> > CLO_ptr;
  typedef RCP<const LinearOpWithSolveBase<Scalar> > CLOWS_ptr;
  typedef ModelEvaluatorBase MEB;
  typedef MEB::Derivative<Scalar> Deriv;

  if( verbLevel == Teuchos::VERB_NONE && verbLevel == Teuchos::VERB_DEFAULT )
    return;

  RCP<Teuchos::FancyOStream>
    out = Teuchos::rcp(&out_arg,false);
  OSTab tab(out);

  *out <<"Thyra::ModelEvaluatorBase::OutArgs<"<<ST::name()<<">:\n";
  tab.incrTab();

  *out <<"model = " << modelEvalDescription_ << "\n";
  *out <<"Np = " << Np() << "\n";
  *out <<"Ng = " << Ng() << "\n";

  CV_ptr f;
  if (this->supports(OUT_ARG_f) && !is_null(f=get_f()) ) {
    *out << "f = " << Teuchos::describe(*f,verbLevel);
  }

  for( int j = 0; j < Ng(); ++j ) {
    CV_ptr g_j;
    if (!is_null(g_j=this->get_g(j)))
      *out << "g("<<j<<") = " << Teuchos::describe(*g_j,verbLevel);
  }

  CLOWS_ptr W;
  if ( this->supports(OUT_ARG_W) && !is_null(W=get_W()) ) {
    *out << "W = " << Teuchos::describe(*W,verbLevel);
  }

  CLO_ptr W_op;
  if ( this->supports(OUT_ARG_W_op) && !is_null(W_op=get_W_op()) ) {
    *out << "W_op = " << Teuchos::describe(*W_op,verbLevel);
  }

  for( int l = 0; l < Np(); ++l ) {
    Deriv DfDp_l;
    if (
      !this->supports(OUT_ARG_DfDp,l).none()
      && !(DfDp_l=get_DfDp(l)).isEmpty()
      )
    {
      *out << "DfDp("<<l<<") = ";
      DfDp_l.describe(*out,verbLevel);
    }
  }

  for( int j = 0; j < Ng(); ++j ) {

    Deriv DgDx_dot_j;
    if (
      !this->supports(OUT_ARG_DgDx_dot,j).none()
      && !(DgDx_dot_j=get_DgDx_dot(j)).isEmpty()
      )
    {
      *out << "DgDx_dot("<<j<<") = ";
      DgDx_dot_j.describe(*out,verbLevel);
    }

    Deriv DgDx_j;
    if (
      !this->supports(OUT_ARG_DgDx,j).none()
      && !(DgDx_j=get_DgDx(j)).isEmpty()
      )
    {
      *out << "DgDx("<<j<<") = ";
      DgDx_j.describe(*out,verbLevel);
    }

    for( int l = 0; l < Np(); ++l ) {

      Deriv DgDp_j_l;
      if (
        !this->supports(OUT_ARG_DgDp,j,l).none()
        && !(DgDp_j_l=get_DgDp(j,l)).isEmpty()
        )
      {
        *out << "DgDp("<<j<<","<<l<<") = ";
        DgDp_j_l.describe(*out,verbLevel);
      }
    }

  }

  // ToDo: Add output for more objects?

}


// protected


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setModelEvalDescription(
  const std::string &modelEvalDescription_in
  )
{
  modelEvalDescription_ = modelEvalDescription_in;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_Np_Ng(int Np_in, int Ng_in)
{
  if(Np_in) {
    supports_DfDp_.resize(Np_in);
    DfDp_.resize(Np_in); std::fill_n(DfDp_.begin(),Np_in,Derivative<Scalar>());
    DfDp_properties_.resize(Np_in); std::fill_n(DfDp_properties_.begin(),Np_in,DerivativeProperties());

    supports_DfDp_mp_.resize(Np_in);
    DfDp_mp_.resize(Np_in); std::fill_n(DfDp_mp_.begin(),Np_in,MPDerivative());
    DfDp_mp_properties_.resize(Np_in); std::fill_n(DfDp_mp_properties_.begin(),Np_in,DerivativeProperties());
  }
  if(Ng_in) {
    g_.resize(Ng_in); std::fill_n(g_.begin(),Ng_in,Teuchos::null);
    supports_DgDx_dot_.resize(Ng_in);
    DgDx_dot_.resize(Ng_in); std::fill_n(DgDx_dot_.begin(),Ng_in,Derivative<Scalar>());
    DgDx_dot_properties_.resize(Ng_in); std::fill_n(DgDx_dot_properties_.begin(),Ng_in,DerivativeProperties());
    supports_DgDx_.resize(Ng_in);
    DgDx_.resize(Ng_in); std::fill_n(DgDx_.begin(),Ng_in,Derivative<Scalar>());
    DgDx_properties_.resize(Ng_in); std::fill_n(DgDx_properties_.begin(),Ng_in,DerivativeProperties());

#ifdef Thyra_BUILD_HESSIAN_SUPPORT
    supports_hess_vec_prod_g_xx_.resize(Ng_in);
    hess_vec_prod_g_xx_.resize(Ng_in); std::fill_n(hess_vec_prod_g_xx_.begin(),Ng_in,RCP<MultiVectorBase<Scalar> >());

    supports_hess_g_xx_.resize(Ng_in);
    hess_g_xx_.resize(Ng_in); std::fill_n(hess_g_xx_.begin(),Ng_in,RCP<LinearOpBase<Scalar> >());
#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

    g_mp_.resize(Ng_in); std::fill_n(g_mp_.begin(),Ng_in,Teuchos::null);
    supports_g_mp_.resize(Ng_in);
    supports_DgDx_dot_mp_.resize(Ng_in);
    DgDx_dot_mp_.resize(Ng_in); std::fill_n(DgDx_dot_mp_.begin(),Ng_in,MPDerivative());
    DgDx_dot_mp_properties_.resize(Ng_in); std::fill_n(DgDx_dot_mp_properties_.begin(),Ng_in,DerivativeProperties());
    supports_DgDx_mp_.resize(Ng_in);
    DgDx_mp_.resize(Ng_in); std::fill_n(DgDx_mp_.begin(),Ng_in,MPDerivative());
    DgDx_mp_properties_.resize(Ng_in); std::fill_n(DgDx_mp_properties_.begin(),Ng_in,DerivativeProperties());
  }
#ifdef Thyra_BUILD_HESSIAN_SUPPORT
  if(Np_in) {
    const int Np = Np_in;
    const int NpNp = Np_in*Np_in;

    supports_hess_vec_prod_f_xp_.resize(Np);
    hess_vec_prod_f_xp_.resize(Np); std::fill_n(hess_vec_prod_f_xp_.begin(),Np,RCP<MultiVectorBase<Scalar> >());

    supports_hess_vec_prod_f_px_.resize(Np);
    hess_vec_prod_f_px_.resize(Np); std::fill_n(hess_vec_prod_f_px_.begin(),Np,RCP<MultiVectorBase<Scalar> >());

    supports_hess_vec_prod_f_pp_.resize(NpNp);
    hess_vec_prod_f_pp_.resize(NpNp); std::fill_n(hess_vec_prod_f_pp_.begin(),NpNp,RCP<MultiVectorBase<Scalar> >());

    supports_hess_f_xp_.resize(Np);
    hess_f_xp_.resize(Np); std::fill_n(hess_f_xp_.begin(),Np,RCP<LinearOpBase<Scalar> >());

    supports_hess_f_pp_.resize(NpNp);
    hess_f_pp_.resize(NpNp); std::fill_n(hess_f_pp_.begin(),NpNp,RCP<LinearOpBase<Scalar> >());

    supports_H_xp_.resize(Np);
    H_xp_.resize(Np); std::fill_n(H_xp_.begin(),Np,RCP<LinearOpBase<Scalar> >());

    supports_H_pp_.resize(NpNp);
    H_pp_.resize(NpNp); std::fill_n(H_pp_.begin(),NpNp,RCP<LinearOpBase<Scalar> >());
  }
#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT
  if(Np_in && Ng_in) {
    const int NpNg = Np_in*Ng_in;
    const int NpNpNg = Np_in*Np_in*Ng_in;
    supports_DgDp_.resize(NpNg);
    DgDp_.resize(NpNg); std::fill_n(DgDp_.begin(),NpNg,Derivative<Scalar>());
    DgDp_properties_.resize(NpNg); std::fill_n(DgDp_properties_.begin(),NpNg,DerivativeProperties());

#ifdef Thyra_BUILD_HESSIAN_SUPPORT
    supports_hess_vec_prod_g_xp_.resize(NpNg);
    hess_vec_prod_g_xp_.resize(NpNg); std::fill_n(hess_vec_prod_g_xp_.begin(),NpNg,RCP<MultiVectorBase<Scalar> >());

    supports_hess_vec_prod_g_px_.resize(NpNg);
    hess_vec_prod_g_px_.resize(NpNg); std::fill_n(hess_vec_prod_g_px_.begin(),NpNg,RCP<MultiVectorBase<Scalar> >());

    supports_hess_vec_prod_g_pp_.resize(NpNpNg);
    hess_vec_prod_g_pp_.resize(NpNpNg); std::fill_n(hess_vec_prod_g_pp_.begin(),NpNpNg,RCP<MultiVectorBase<Scalar> >());

    supports_hess_g_xp_.resize(NpNg);
    hess_g_xp_.resize(NpNg); std::fill_n(hess_g_xp_.begin(),NpNg,RCP<LinearOpBase<Scalar> >());

    supports_hess_g_pp_.resize(NpNpNg);
    hess_g_pp_.resize(NpNpNg); std::fill_n(hess_g_pp_.begin(),NpNpNg,RCP<LinearOpBase<Scalar> >());
#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

    supports_DgDp_mp_.resize(NpNg);
    DgDp_mp_.resize(NpNg); std::fill_n(DgDp_mp_.begin(),NpNg,MPDerivative());
    DgDp_mp_properties_.resize(NpNg); std::fill_n(DgDp_mp_properties_.begin(),NpNg,DerivativeProperties());
  }
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgsMembers arg, bool supports_in )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    int(arg)>=NUM_E_OUT_ARGS_MEMBERS || int(arg) < 0,std::logic_error
    ,"model = \'"<<modelEvalDescription_
    <<"\': Error, arg="<<toString(arg)<<" is invalid!"
    );
  supports_[arg] = supports_in;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgsDfDp /* arg */, int l, const DerivativeSupport& supports_in
  )
{
  assert_supports(OUT_ARG_f);
  assert_l(l);
  supports_DfDp_[l] = supports_in;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgsDgDx_dot /* arg */, int j, const DerivativeSupport& supports_in
  )
{
  assert_j(j);
  supports_DgDx_dot_[j] = supports_in;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgsDgDx /* arg */, int j, const DerivativeSupport& supports_in
  )
{
  assert_j(j);
  supports_DgDx_[j] = supports_in;
}


#ifdef Thyra_BUILD_HESSIAN_SUPPORT

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_hess_vec_prod_f_xx /* arg */, bool supports_in
  )
{
  supports_hess_vec_prod_f_xx_ = supports_in;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_hess_vec_prod_f_xp /* arg */, int l, bool supports_in
  )
{
  assert_l(l);
  supports_hess_vec_prod_f_xp_[l] = supports_in;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_hess_vec_prod_f_px /* arg */, int l, bool supports_in
  )
{
  assert_l(l);
  supports_hess_vec_prod_f_px_[l] = supports_in;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_hess_vec_prod_f_pp /* arg */, int l1, int l2, bool supports_in
  )
{
  assert_l(l1);
  assert_l(l2);
  supports_hess_vec_prod_f_pp_[ l1*Np()+ l2 ] = supports_in;
}

#endif // ifdef Thyra_BUILD_HESSIAN_SUPPORT


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgsDgDp /* arg */, int j, int l, const DerivativeSupport& supports_in
  )
{
  assert_j(j);
  assert_l(l);
  supports_DgDp_[ j*Np()+ l ] = supports_in;
}


#ifdef Thyra_BUILD_HESSIAN_SUPPORT

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_hess_vec_prod_g_xx /* arg */, int j, bool supports_in
  )
{
  assert_j(j);
  supports_hess_vec_prod_g_xx_[j] = supports_in;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_hess_vec_prod_g_xp /* arg */, int j, int l, bool supports_in
  )
{
  assert_j(j);
  assert_l(l);
  supports_hess_vec_prod_g_xp_[ j*Np()+ l ] = supports_in;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_hess_vec_prod_g_px /* arg */, int j, int l, bool supports_in
  )
{
  assert_j(j);
  assert_l(l);
  supports_hess_vec_prod_g_px_[ j*Np()+ l ] = supports_in;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_hess_vec_prod_g_pp /* arg */, int j, int l1, int l2, bool supports_in
  )
{
  assert_j(j);
  assert_l(l1);
  assert_l(l2);
  supports_hess_vec_prod_g_pp_[ j*Np()*Np()+ l1*Np() + l2 ] = supports_in;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_hess_f_xx /* arg */, bool supports_in
  )
{
  supports_hess_f_xx_ = supports_in;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_hess_f_xp /* arg */, int l, bool supports_in
  )
{
  assert_l(l);
  supports_hess_f_xp_[l] = supports_in;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_hess_f_pp /* arg */, int l1, int l2, bool supports_in
  )
{
  assert_l(l1);
  assert_l(l2);
  supports_hess_f_pp_[ l1*Np()+ l2 ] = supports_in;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_hess_g_xx /* arg */, int j, bool supports_in
  )
{
  assert_j(j);
  supports_hess_g_xx_[j] = supports_in;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_hess_g_xp /* arg */, int j, int l, bool supports_in
  )
{
  assert_j(j);
  assert_l(l);
  supports_hess_g_xp_[ j*Np()+ l ] = supports_in;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_hess_g_pp /* arg */, int j, int l1, int l2, bool supports_in
  )
{
  assert_j(j);
  assert_l(l1);
  assert_l(l2);
  supports_hess_g_pp_[ j*Np()*Np()+ l1*Np() + l2 ] = supports_in;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_H_xx /* arg */, bool supports_in
  )
{
  supports_H_xx_ = supports_in;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_H_xp /* arg */, int l, bool supports_in
  )
{
  assert_l(l);
  supports_H_xp_[l] = supports_in;
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_H_pp /* arg */, int l1, int l2, bool supports_in
  )
{
  assert_l(l1);
  assert_l(l2);
  supports_H_pp_[ l1*Np()+ l2 ] = supports_in;
}

#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgs_g_mp /* arg */, int j, bool supports_in
  )
{
  //assert_supports(OUT_ARG_g_mp,j);
  assert_j(j);
  supports_g_mp_[j] = supports_in;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgsDfDp_mp /* arg */, int l, const DerivativeSupport& supports_in
  )
{
  assert_supports(OUT_ARG_f_mp); //JF this is not in epetraext ME
  assert_l(l);
  supports_DfDp_mp_[l] = supports_in;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgsDgDx_dot_mp /* arg */, int j, const DerivativeSupport& supports_in
  )
{
  assert_j(j);
  supports_DgDx_dot_mp_[j] = supports_in;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgsDgDx_mp /* arg */, int j, const DerivativeSupport& supports_in
  )
{
  assert_j(j);
  supports_DgDx_mp_[j] = supports_in;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  EOutArgsDgDp_mp /* arg */, int j, int l, const DerivativeSupport& supports_in
  )
{
  assert_j(j);
  assert_l(l);
  supports_DgDp_mp_[ j*Np()+ l ] = supports_in;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_W_properties(
  const DerivativeProperties &properties
  )
{
  W_properties_ = properties;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_DfDp_properties(
  int l, const DerivativeProperties &properties
  )
{
  assert_supports(OUT_ARG_DfDp,l);
  DfDp_properties_[l] = properties;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_DgDx_dot_properties(
  int j, const DerivativeProperties &properties
  )
{
  assert_supports(OUT_ARG_DgDx_dot,j);
  DgDx_dot_properties_[j] = properties;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_DgDx_properties(
  int j, const DerivativeProperties &properties
  )
{
  assert_supports(OUT_ARG_DgDx,j);
  DgDx_properties_[j] = properties;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_DgDp_properties(
  int j, int l, const DerivativeProperties &properties
  )
{
  assert_supports(OUT_ARG_DgDp,j,l);
  DgDp_properties_[ j*Np()+ l ] = properties;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_DfDp_mp_properties(
  int l, const DerivativeProperties &properties
  )
{
  assert_supports(OUT_ARG_DfDp_mp,l);
  DfDp_mp_properties_[l] = properties;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_DgDx_dot_mp_properties(
  int j, const DerivativeProperties &properties
  )
{
  assert_supports(OUT_ARG_DgDx_dot_mp,j);
  DgDx_dot_mp_properties_[j] = properties;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_DgDx_mp_properties(
  int j, const DerivativeProperties &properties
  )
{
  assert_supports(OUT_ARG_DgDx_mp,j);
  DgDx_mp_properties_[j] = properties;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_set_DgDp_mp_properties(
  int j, int l, const DerivativeProperties &properties
  )
{
  assert_supports(OUT_ARG_DgDp_mp,j,l);
  DgDp_mp_properties_[ j*Np()+ l ] = properties;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setSupports(
  const OutArgs<Scalar>& inputOutArgs
  )
{
  typedef ModelEvaluatorBase MEB;
  const int l_Np = TEUCHOS_MIN(this->Np(),inputOutArgs.Np());
  const int l_Ng = TEUCHOS_MIN(this->Ng(),inputOutArgs.Ng());
  std::copy(
    &inputOutArgs.supports_[0],
    &inputOutArgs.supports_[0] + NUM_E_OUT_ARGS_MEMBERS, &supports_[0] );
  for( int l = 0; l < l_Np; ++l ) {
    DerivativeSupport ds = inputOutArgs.supports(MEB::OUT_ARG_DfDp,l);
    if (!ds.none()) {
      this->_setSupports(MEB::OUT_ARG_DfDp,l,ds);
      this->_set_DfDp_properties(l,inputOutArgs.get_DfDp_properties(l));
    }
  }
  for( int l = 0; l < l_Np; ++l ) {
    DerivativeSupport ds = inputOutArgs.supports(MEB::OUT_ARG_DfDp_mp,l);
    if (!ds.none()) {
      this->_setSupports(MEB::OUT_ARG_DfDp_mp,l,ds);
      this->_set_DfDp_mp_properties(l,inputOutArgs.get_DfDp_mp_properties(l));
    }
  }
  for( int j = 0; j < l_Ng; ++j ) {
    DerivativeSupport ds = inputOutArgs.supports(MEB::OUT_ARG_DgDx_dot,j);
    this->_setSupports(MEB::OUT_ARG_DgDx_dot,j,ds);
    if(!ds.none()) this->_set_DgDx_dot_properties(j,inputOutArgs.get_DgDx_dot_properties(j));
  }
  for( int j = 0; j < l_Ng; ++j ) {
    DerivativeSupport ds = inputOutArgs.supports(MEB::OUT_ARG_DgDx_dot_mp,j);
    this->_setSupports(MEB::OUT_ARG_DgDx_dot_mp,j,ds);
    if(!ds.none()) this->_set_DgDx_dot_mp_properties(j,inputOutArgs.get_DgDx_dot_mp_properties(j));
  }
  for( int j = 0; j < l_Ng; ++j ) {
    DerivativeSupport ds = inputOutArgs.supports(MEB::OUT_ARG_DgDx,j);
    this->_setSupports(MEB::OUT_ARG_DgDx,j,ds);
    if(!ds.none()) this->_set_DgDx_properties(j,inputOutArgs.get_DgDx_properties(j));
  }
  for( int j = 0; j < l_Ng; ++j ) {
    DerivativeSupport ds = inputOutArgs.supports(MEB::OUT_ARG_DgDx_mp,j);
    this->_setSupports(MEB::OUT_ARG_DgDx_mp,j,ds);
    if(!ds.none()) this->_set_DgDx_mp_properties(j,inputOutArgs.get_DgDx_mp_properties(j));
  }
  for( int j = 0; j < l_Ng; ++j ) for( int l = 0; l < l_Np; ++l ) {
    DerivativeSupport ds = inputOutArgs.supports(MEB::OUT_ARG_DgDp,j,l);
    this->_setSupports(MEB::OUT_ARG_DgDp,j,l,ds);
    if(!ds.none()) this->_set_DgDp_properties(j,l,inputOutArgs.get_DgDp_properties(j,l));
  }
  for( int j = 0; j < l_Ng; ++j ) for( int l = 0; l < l_Np; ++l ) {
    DerivativeSupport ds = inputOutArgs.supports(MEB::OUT_ARG_DgDp_mp,j,l);
    this->_setSupports(MEB::OUT_ARG_DgDp_mp,j,l,ds);
    if(!ds.none()) this->_set_DgDp_mp_properties(j,l,inputOutArgs.get_DgDp_mp_properties(j,l));
  }
  if(this->supports(OUT_ARG_W) || this->supports(OUT_ARG_W_op))
    this->_set_W_properties(inputOutArgs.get_W_properties());
  if(this->supports(OUT_ARG_W_mp))
    this->_set_W_properties(inputOutArgs.get_W_properties());  //JF should this be W_mp_properties?

#ifdef Thyra_BUILD_HESSIAN_SUPPORT
  this->_setSupports(MEB::OUT_ARG_hess_vec_prod_f_xx,inputOutArgs.supports(OUT_ARG_hess_vec_prod_f_xx));
  for( int l1 = 0; l1 < l_Np; ++l1 ) {
    this->_setSupports(MEB::OUT_ARG_hess_vec_prod_f_xp,l1,inputOutArgs.supports(MEB::OUT_ARG_hess_vec_prod_f_xp,l1));
    this->_setSupports(MEB::OUT_ARG_hess_vec_prod_f_px,l1,inputOutArgs.supports(MEB::OUT_ARG_hess_vec_prod_f_px,l1));
    for( int l2 = 0; l2 < l_Np; ++l2 ) {
      this->_setSupports(MEB::OUT_ARG_hess_vec_prod_f_pp,l1,l2,inputOutArgs.supports(MEB::OUT_ARG_hess_vec_prod_f_pp,l1,l2));
    }
  }
  for( int j = 0; j < l_Ng; ++j ) {
    this->_setSupports(MEB::OUT_ARG_hess_vec_prod_g_xx,j,inputOutArgs.supports(MEB::OUT_ARG_hess_vec_prod_g_xx,j));
    for( int l1 = 0; l1 < l_Np; ++l1 ) {
      this->_setSupports(MEB::OUT_ARG_hess_vec_prod_g_xp,j,l1,inputOutArgs.supports(MEB::OUT_ARG_hess_vec_prod_g_xp,j,l1));
      this->_setSupports(MEB::OUT_ARG_hess_vec_prod_g_px,j,l1,inputOutArgs.supports(MEB::OUT_ARG_hess_vec_prod_g_px,j,l1));
      for( int l2 = 0; l2 < l_Np; ++l2 ) {
        this->_setSupports(MEB::OUT_ARG_hess_vec_prod_g_pp,j,l1,l2,inputOutArgs.supports(MEB::OUT_ARG_hess_vec_prod_g_pp,j,l1,l2));
      }
    }
  }
  this->_setSupports(MEB::OUT_ARG_hess_f_xx,inputOutArgs.supports(OUT_ARG_hess_f_xx));
  for( int l1 = 0; l1 < l_Np; ++l1 ) {
    this->_setSupports(MEB::OUT_ARG_hess_f_xp,l1,inputOutArgs.supports(MEB::OUT_ARG_hess_f_xp,l1));
    for( int l2 = 0; l2 < l_Np; ++l2 ) {
      this->_setSupports(MEB::OUT_ARG_hess_f_pp,l1,l2,inputOutArgs.supports(MEB::OUT_ARG_hess_f_pp,l1,l2));
    }
  }
  for( int j = 0; j < l_Ng; ++j ) {
    this->_setSupports(MEB::OUT_ARG_hess_g_xx,j,inputOutArgs.supports(MEB::OUT_ARG_hess_g_xx,j));
    for( int l1 = 0; l1 < l_Np; ++l1 ) {
      this->_setSupports(MEB::OUT_ARG_hess_g_xp,j,l1,inputOutArgs.supports(MEB::OUT_ARG_hess_g_xp,j,l1));
      for( int l2 = 0; l2 < l_Np; ++l2 ) {
        this->_setSupports(MEB::OUT_ARG_hess_g_pp,j,l1,l2,inputOutArgs.supports(MEB::OUT_ARG_hess_g_pp,j,l1,l2));
      }
    }
  }
  this->_setSupports(MEB::OUT_ARG_H_xx,inputOutArgs.supports(OUT_ARG_H_xx));
  for( int l1 = 0; l1 < l_Np; ++l1 ) {
    this->_setSupports(MEB::OUT_ARG_H_xp,l1,inputOutArgs.supports(MEB::OUT_ARG_H_xp,l1));
    for( int l2 = 0; l2 < l_Np; ++l2 ) {
      this->_setSupports(MEB::OUT_ARG_H_pp,l1,l2,inputOutArgs.supports(MEB::OUT_ARG_H_pp,l1,l2));
    }
  }
#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

  extended_outargs_ = inputOutArgs.extended_outargs_;
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setUnsupportsAndRelated(
  EInArgsMembers arg
  )
{
  switch(arg) {
    case IN_ARG_x: {
      const int l_Ng = this->Ng();
      for( int j = 0; j < l_Ng; ++j ) {
        this->_setSupports(OUT_ARG_DgDx_dot,j,DerivativeSupport());
        this->_setSupports(OUT_ARG_DgDx,j,DerivativeSupport());
      }
      break;
    }
    case IN_ARG_x_mp: {
      const int l_Ng = this->Ng();
      for( int j = 0; j < l_Ng; ++j ) {
        this->_setSupports(OUT_ARG_DgDx_dot_mp,j,DerivativeSupport());
        this->_setSupports(OUT_ARG_DgDx_mp,j,DerivativeSupport());
      }
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(
        true ,std::logic_error,
        "Error, can not handle args other than IN_ARG_x yet!"
        );
  }
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setUnsupportsAndRelated(
  EOutArgsMembers arg
  )
{
  switch(arg) {
    case OUT_ARG_f: {
      this->_setSupports(OUT_ARG_W,false);
      this->_setSupports(OUT_ARG_W_op,false);
      this->_setSupports(OUT_ARG_f_poly,false);
      const int l_Np = this->Np();
      for( int l = 0; l < l_Np; ++l )
        this->_setSupports(OUT_ARG_DfDp,l,DerivativeSupport());
      break;
    }
    case OUT_ARG_f_mp: {
      this->_setSupports(OUT_ARG_W_mp,false);
      this->_setSupports(OUT_ARG_W_op,false);
      this->_setSupports(OUT_ARG_f_poly,false);
      const int l_Np = this->Np();
      for( int l = 0; l < l_Np; ++l )
        this->_setSupports(OUT_ARG_DfDp_mp,l,DerivativeSupport());
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(
        true ,std::logic_error,
        "Error, can not handle args other than OUT_ARG_f yet!"
        );
  }
  this->_setSupports(arg,false);
}

#ifdef Thyra_BUILD_HESSIAN_SUPPORT
template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::_setHessianSupports(
  const bool supports
  )
{
  typedef ModelEvaluatorBase MEB;
  this->_setSupports(MEB::OUT_ARG_hess_vec_prod_f_xx,supports);
  this->_setSupports(MEB::OUT_ARG_hess_f_xx,supports);
  this->_setSupports(MEB::OUT_ARG_H_xx,supports);
  for (int l1=0; l1<this->Np(); ++l1)
  {
    this->_setSupports(MEB::OUT_ARG_hess_vec_prod_f_xp,l1,supports);
    this->_setSupports(MEB::OUT_ARG_hess_vec_prod_f_px,l1,supports);
    this->_setSupports(MEB::OUT_ARG_hess_f_xp,l1,supports);
    this->_setSupports(MEB::OUT_ARG_H_xp,l1,supports);
    for (int l2=0; l2<this->Np(); ++l2)
    {
      this->_setSupports(MEB::OUT_ARG_hess_vec_prod_f_pp,l1,l2,supports);
      this->_setSupports(MEB::OUT_ARG_hess_f_pp,l1,l2,supports);
      this->_setSupports(MEB::OUT_ARG_H_pp,l1,l2,supports);
    }
  }

  for (int j=0; j<this->Ng(); ++j)
  {
    this->_setSupports(MEB::OUT_ARG_hess_g_xx,j,supports);
    this->_setSupports(MEB::OUT_ARG_hess_vec_prod_g_xx,j,supports);
    for (int l1=0; l1<this->Np(); ++l1)
    {
      this->_setSupports(MEB::OUT_ARG_hess_vec_prod_g_xp,j,l1,supports);
      this->_setSupports(MEB::OUT_ARG_hess_vec_prod_g_px,j,l1,supports);
      this->_setSupports(MEB::OUT_ARG_hess_g_xp,j,l1,supports);
      for (int l2=0; l2<this->Np(); ++l2)
      {
        this->_setSupports(MEB::OUT_ARG_hess_vec_prod_g_pp,j,l1,l2,supports);
        this->_setSupports(MEB::OUT_ARG_hess_g_pp,j,l1,l2,supports);
      }
    }
  }
}
#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

// private


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(EOutArgsMembers arg) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !this->supports(arg), std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(arg):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument arg = " << toString(arg) << " is not supported!"
    );
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgsDfDp arg, int l, const Derivative<Scalar> &deriv
  ) const
{
  const DerivativeSupport derivSupport = this->supports(arg,l);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !deriv.isSupportedBy(derivSupport), std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_DfDp,l):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument DfDp("<<l<<") = " << deriv.description() << "\n"
    "is not supported!\n\n"
    "The supported types include " << derivSupport.description() << "!"
    );
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgsDgDx_dot arg, int j, const Derivative<Scalar> &deriv
  ) const
{
  const DerivativeSupport derivSupport = this->supports(arg,j);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !deriv.isSupportedBy(derivSupport), std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_DgDx_dot,j):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument DgDx_dot("<<j<<") = " << deriv.description() << "\n"
    "is not supported!\n\n"
    "The supported types include " << derivSupport.description() << "!"
    );
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgsDgDx arg, int j, const Derivative<Scalar> &deriv
  ) const
{
  const DerivativeSupport derivSupport = this->supports(arg,j);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !deriv.isSupportedBy(derivSupport), std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_DgDx,j):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument DgDx("<<j<<") = " << deriv.description() << "\n"
    "is not supported!\n\n"
    "The supported types include " << derivSupport.description() << "!"
    );
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgsDgDp arg, int j, int l, const Derivative<Scalar> &deriv
  ) const
{
  const DerivativeSupport derivSupport = this->supports(arg,j,l);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !deriv.isSupportedBy(derivSupport), std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_DgDp,j,l):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument DgDp("<<j<<","<<l<<") = " << deriv.description() << "\n"
    "is not supported!\n\n"
    "The supported types include " << derivSupport.description() << "!"
    );
}


#ifdef Thyra_BUILD_HESSIAN_SUPPORT

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_hess_vec_prod_f_xx arg
  ) const
{
  const bool support = this->supports(arg);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !support, std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_hess_vec_prod_f_xx):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument hess_vec_prod_f_xx() is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_hess_vec_prod_f_xp arg, int l
  ) const
{
  const bool support = this->supports(arg,l);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !support, std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_hess_vec_prod_f_xp,l):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument hess_vec_prod_f_xp("<<l<<") is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_hess_vec_prod_f_px arg, int l
  ) const
{
  const bool support = this->supports(arg,l);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !support, std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_hess_vec_prod_f_px,l):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument hess_vec_prod_f_px("<<l<<") is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_hess_vec_prod_f_pp arg, int l1, int l2
  ) const
{
  const bool support = this->supports(arg,l1,l2);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !support, std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_hess_vec_prod_f_pp,l1,l2):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument hess_vec_prod_f_pp("<<l1<<","<<l2<<") is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_hess_vec_prod_g_xx arg, int j
  ) const
{
  const bool support = this->supports(arg,j);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !support, std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_hess_vec_prod_g_xx,j):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument hess_vec_prod_g_xx("<<j<<") is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_hess_vec_prod_g_xp arg, int j, int l
  ) const
{
  const bool support = this->supports(arg,j,l);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !support, std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_hess_vec_prod_g_xp,j,l):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument hess_vec_prod_g_xp("<<j<<","<<l<<") is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_hess_vec_prod_g_px arg, int j, int l
  ) const
{
  const bool support = this->supports(arg,j,l);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !support, std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_hess_vec_prod_g_px,j,l):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument hess_vec_prod_g_px("<<j<<","<<l<<") is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_hess_vec_prod_g_pp arg, int j, int l1, int l2
  ) const
{
  const bool support = this->supports(arg,j,l1,l2);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !support, std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_hess_vec_prod_g_pp,j,l1,l2):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument hess_vec_prod_g_pp("<<j<<","<<l1<<","<<l2<<") is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_hess_f_xx arg
  ) const
{
  const bool support = this->supports(arg);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !support, std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_hess_f_xx):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument hess_f_xx() is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_hess_f_xp arg, int l
  ) const
{
  const bool support = this->supports(arg,l);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !support, std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_hess_f_xp,l):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument hess_f_xp() is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_hess_f_pp arg, int l1, int l2
  ) const
{
  const bool support = this->supports(arg,l1,l2);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !support, std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_hess_f_pp,l1,l2):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument hess_f_pp() is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_hess_g_xx arg, int j
  ) const
{
  const bool support = this->supports(arg,j);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !support, std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_hess_g_xx,j):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument hess_g_xx() is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_hess_g_xp arg, int j, int l
  ) const
{
  const bool support = this->supports(arg,j,l);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !support, std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_hess_g_xp,j,l):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument hess_g_xp() is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_hess_g_pp arg, int j, int l1, int l2
  ) const
{
  const bool support = this->supports(arg,j,l1,l2);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !support, std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_hess_g_pp,j,l1,l2):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument hess_g_pp() is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_H_xx arg
  ) const
{
  const bool support = this->supports(arg);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !support, std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_H_xx):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument H_xx() is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_H_xp arg, int l
  ) const
{
  const bool support = this->supports(arg,l);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !support, std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_H_xp,l):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument H_xp() is not supported!"
    );
}

template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_H_pp arg, int l1, int l2
  ) const
{
  const bool support = this->supports(arg,l1,l2);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !support, std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_H_pp,l1,l2):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument H_pp() is not supported!"
    );
}

#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgs_g_mp /* arg */, int j
  ) const
{
  assert_j(j);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !supports_g_mp_[j], std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_g_mp,j):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument g_mp("<<j<<") \n"
    "is not supported!\n\n"
    );
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgsDfDp_mp arg, int l, const MPDerivative &deriv
  ) const
{
  const DerivativeSupport derivSupport = this->supports(arg,l);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !deriv.isSupportedBy(derivSupport), std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_DfDp_mp,l):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument DfDp_mp("<<l<<") = " << deriv.description() << "\n"
    "is not supported!\n\n"
    "The supported types include " << derivSupport.description() << "!"
    );
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgsDgDx_dot_mp arg, int j, const MPDerivative &deriv
  ) const
{
  const DerivativeSupport derivSupport = this->supports(arg,j);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !deriv.isSupportedBy(derivSupport), std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_DgDx_dot_mp,j):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument DgDx_dot_mp("<<j<<") = " << deriv.description() << "\n"
    "is not supported!\n\n"
    "The supported types include " << derivSupport.description() << "!"
    );
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgsDgDx_mp arg, int j, const MPDerivative &deriv
  ) const
{
  const DerivativeSupport derivSupport = this->supports(arg,j);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !deriv.isSupportedBy(derivSupport), std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_DgDx_mp,j):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument DgDx_mp("<<j<<") = " << deriv.description() << "\n"
    "is not supported!\n\n"
    "The supported types include " << derivSupport.description() << "!"
    );
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(
  EOutArgsDgDp_mp arg, int j, int l, const MPDerivative &deriv
  ) const
{
  const DerivativeSupport derivSupport = this->supports(arg,j,l);
  TEUCHOS_TEST_FOR_EXCEPTION(
    !deriv.isSupportedBy(derivSupport), std::logic_error,
    "Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_supports(OUT_ARG_DgDp_mp,j,l):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The argument DgDp_mp("<<j<<","<<l<<") = " << deriv.description() << "\n"
    "is not supported!\n\n"
    "The supported types include " << derivSupport.description() << "!"
    );
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_l(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !( 0 <= l && l < Np() ), std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_l(l):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error,  The parameter subvector p("<<l<<")"
    " is not in the range [0,"<<Np()<<")!"
    );
}


template<class Scalar>
void ModelEvaluatorBase::OutArgs<Scalar>::assert_j(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !( 0 <= j && j < Ng() ), std::logic_error
    ,"Thyra::ModelEvaluatorBase::OutArgs<Scalar>::assert_j(j):\n\n"
    "model = \'"<<modelEvalDescription_<<"\':\n\n"
    "Error, The auxiliary function g("<<j<<")"
    " is not in the range [0,"<<Ng()<<")!"
    );
}


//
// ModelEvaluatorBase::InArgsSetup
//


template<class Scalar>
ModelEvaluatorBase::InArgsSetup<Scalar>::InArgsSetup()
{}


template<class Scalar>
ModelEvaluatorBase::InArgsSetup<Scalar>::InArgsSetup( const InArgs<Scalar>& inArgs )
  :InArgs<Scalar>(inArgs)
{}


template<class Scalar>
void ModelEvaluatorBase::InArgsSetup<Scalar>::setModelEvalDescription(
  const std::string &modelEvalDescription_in )
{
  this->_setModelEvalDescription(modelEvalDescription_in);
}


template<class Scalar>
void ModelEvaluatorBase::InArgsSetup<Scalar>::set_Np(int Np_in)
{ this->_set_Np_Ng(Np_in, 0); }

template<class Scalar>
void ModelEvaluatorBase::InArgsSetup<Scalar>::set_Np_Ng(int Np_in, int Ng_in)
{ this->_set_Np_Ng(Np_in, Ng_in); }

template<class Scalar>
void ModelEvaluatorBase::InArgsSetup<Scalar>::setSupports( EInArgsMembers arg, bool supports_in )
{ this->_setSupports(arg,supports_in); }


template<class Scalar>
void ModelEvaluatorBase::InArgsSetup<Scalar>::setSupports( EInArgs_p_mp arg, int l, bool supports_in )
{ this->_setSupports(arg,l,supports_in); }


template<class Scalar>
void ModelEvaluatorBase::InArgsSetup<Scalar>::setSupports(
  const InArgs<Scalar>& inArgs, const int Np_in
  )
{
  this->_setSupports(inArgs, Np_in);
}


template<class Scalar>
void ModelEvaluatorBase::InArgsSetup<Scalar>::setUnsupportsAndRelated(
  EInArgsMembers arg
  )
{
  this->_setUnsupportsAndRelated(arg);
}


//
// ModelEvaluatorBase::OutArgsSetup
//


template<class Scalar>
ModelEvaluatorBase::OutArgsSetup<Scalar>::OutArgsSetup()
{}


template<class Scalar>
ModelEvaluatorBase::OutArgsSetup<Scalar>::OutArgsSetup(
  const OutArgs<Scalar>& inputOutArgs
  )
  :OutArgs<Scalar>(inputOutArgs)
{}


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setModelEvalDescription(
  const std::string &modelEvalDescription_in
  )
{ this->_setModelEvalDescription(modelEvalDescription_in); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_Np_Ng(int Np_in, int Ng_in)
{ this->_set_Np_Ng(Np_in, Ng_in); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgsMembers arg, bool supports_in
  )
{ this->_setSupports(arg,supports_in); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgsDfDp arg, int l, const DerivativeSupport& supports_in
  )
{ this->_setSupports(arg,l,supports_in); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgsDgDx_dot arg, int j, const DerivativeSupport& supports_in
  )
{ this->_setSupports(arg,j,supports_in); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgsDgDx arg, int j, const DerivativeSupport& supports_in
  )
{ this->_setSupports(arg,j,supports_in); }


#ifdef Thyra_BUILD_HESSIAN_SUPPORT

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_hess_vec_prod_f_xx arg, bool supports_in
  )
{ this->_setSupports(arg,supports_in); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_hess_vec_prod_f_xp arg, int l, bool supports_in
  )
{ this->_setSupports(arg,l,supports_in); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_hess_vec_prod_f_px arg, int l, bool supports_in
  )
{ this->_setSupports(arg,l,supports_in); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_hess_vec_prod_f_pp arg, int l1, int l2, bool supports_in
  )
{ this->_setSupports(arg,l1,l2,supports_in); }

#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgsDgDp arg, int j, int l, const DerivativeSupport& supports_in
  )
{ this->_setSupports(arg,j,l,supports_in); }


#ifdef Thyra_BUILD_HESSIAN_SUPPORT

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_hess_vec_prod_g_xx arg, int j, bool supports_in
  )
{ this->_setSupports(arg,j,supports_in); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_hess_vec_prod_g_xp arg, int j, int l, bool supports_in
  )
{ this->_setSupports(arg,j,l,supports_in); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_hess_vec_prod_g_px arg, int j, int l, bool supports_in
  )
{ this->_setSupports(arg,j,l,supports_in); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_hess_vec_prod_g_pp arg, int j, int l1, int l2, bool supports_in
  )
{ this->_setSupports(arg,j,l1,l2,supports_in); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_hess_f_xx arg, bool supports_in
  )
{ this->_setSupports(arg,supports_in); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_hess_f_xp arg, int l, bool supports_in
  )
{ this->_setSupports(arg,l,supports_in); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_hess_f_pp arg, int l1, int l2, bool supports_in
  )
{ this->_setSupports(arg,l1,l2,supports_in); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_hess_g_xx arg, int j, bool supports_in
  )
{ this->_setSupports(arg,j,supports_in); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_hess_g_xp arg, int j, int l, bool supports_in
  )
{ this->_setSupports(arg,j,l,supports_in); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_hess_g_pp arg, int j, int l1, int l2, bool supports_in
  )
{ this->_setSupports(arg,j,l1,l2,supports_in); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_H_xx arg, bool supports_in
  )
{ this->_setSupports(arg,supports_in); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_H_xp arg, int l, bool supports_in
  )
{ this->_setSupports(arg,l,supports_in); }

template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_H_pp arg, int l1, int l2, bool supports_in
  )
{ this->_setSupports(arg,l1,l2,supports_in); }

#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgs_g_mp arg, int j, bool supports_in
  )
{ this->_setSupports(arg,j,supports_in); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgsDfDp_mp arg, int l, const DerivativeSupport& supports_in
  )
{ this->_setSupports(arg,l,supports_in); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgsDgDx_dot_mp arg, int j, const DerivativeSupport& supports_in
  )
{ this->_setSupports(arg,j,supports_in); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgsDgDx_mp arg, int j, const DerivativeSupport& supports_in
  )
{ this->_setSupports(arg,j,supports_in); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  EOutArgsDgDp_mp arg, int j, int l, const DerivativeSupport& supports_in
  )
{ this->_setSupports(arg,j,l,supports_in); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_W_properties(
  const DerivativeProperties &properties
  )
{ this->_set_W_properties(properties); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_DfDp_properties(
  int l, const DerivativeProperties &properties
  )
{ this->_set_DfDp_properties(l,properties); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_DgDx_dot_properties(
  int j, const DerivativeProperties &properties
  )
{ this->_set_DgDx_dot_properties(j,properties); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_DgDx_properties(
  int j, const DerivativeProperties &properties
  )
{ this->_set_DgDx_properties(j,properties); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_DgDp_properties(
  int j, int l, const DerivativeProperties &properties
  )
{ this->_set_DgDp_properties(j,l,properties); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_DfDp_mp_properties(
  int l, const DerivativeProperties &properties
  )
{ this->_set_DfDp_mp_properties(l,properties); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_DgDx_dot_mp_properties(
  int j, const DerivativeProperties &properties
  )
{ this->_set_DgDx_dot_mp_properties(j,properties); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_DgDx_mp_properties(
  int j, const DerivativeProperties &properties
  )
{ this->_set_DgDx_mp_properties(j,properties); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::set_DgDp_mp_properties(
  int j, int l, const DerivativeProperties &properties
  )
{ this->_set_DgDp_mp_properties(j,l,properties); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setSupports(
  const OutArgs<Scalar>& inputOutArgs
  )
{ this->_setSupports(inputOutArgs); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setUnsupportsAndRelated(
  EInArgsMembers arg
  )
{ this->_setUnsupportsAndRelated(arg); }


template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setUnsupportsAndRelated(
  EOutArgsMembers arg
  )
{ this->_setUnsupportsAndRelated(arg); }


#ifdef Thyra_BUILD_HESSIAN_SUPPORT
template<class Scalar>
void ModelEvaluatorBase::OutArgsSetup<Scalar>::setHessianSupports(
  const bool supports
  )
{ this->_setHessianSupports(supports); }
#endif  // ifdef Thyra_BUILD_HESSIAN_SUPPORT

} // namespace Thyra



//
// Explicit instantiation macro
//
// Must be expanded from within the Thyra namespace!
//


#define THYRA_MODEL_EVALUATOR_BASE_INSTANT(SCALAR) \
  \
  template class ModelEvaluatorBase::InArgs<SCALAR >; \
  \
  template std::string \
  ModelEvaluatorBase::Derivative<SCALAR >::description() const; \
  \
  template \
  void ModelEvaluatorBase::Derivative<SCALAR >::describe( \
    Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel \
    ) const; \
  \
  template class ModelEvaluatorBase::OutArgs<SCALAR >; \
  \
  template class ModelEvaluatorBase::InArgsSetup<SCALAR >; \
  \
  template class ModelEvaluatorBase::OutArgsSetup<SCALAR >;


#endif // THYRA_MODEL_EVALUATOR_BASE_DEF_HPP
