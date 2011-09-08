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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DIRECTIONAL_FINITE_DIFF_CALCULATOR_DEF_HPP
#define THYRA_DIRECTIONAL_FINITE_DIFF_CALCULATOR_DEF_HPP


#include "Thyra_DirectionalFiniteDiffCalculator_decl.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"


namespace Thyra {


namespace DirectionalFiniteDiffCalculatorTypes {


//
// Undocumented utility class for setting support for new derivative objects
// on an OutArgs object!  Warning, users should not attempt to play these
// tricks on their own!
//
// Note that because of the design of the OutArgs and OutArgsSetup classes,
// you can only change the list of supported arguments in a subclass of
// ModelEvaluatorBase since OutArgsSetup is a protected type.  The fact that
// the only way to do this is convoluted is a feature!
//
template<class Scalar>
class OutArgsCreator : public StateFuncModelEvaluatorBase<Scalar>
{
public:
  // Public functions overridden from ModelEvaulator.
  RCP<const VectorSpaceBase<Scalar> > get_x_space() const
    { TEST_FOR_EXCEPT(true); return Teuchos::null; }
  RCP<const VectorSpaceBase<Scalar> > get_f_space() const
    { TEST_FOR_EXCEPT(true); return Teuchos::null; }
  ModelEvaluatorBase::InArgs<Scalar> createInArgs() const
    { TEST_FOR_EXCEPT(true); return ModelEvaluatorBase::InArgs<Scalar>(); }
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const
    { TEST_FOR_EXCEPT(true); return ModelEvaluatorBase::OutArgs<Scalar>(); }
  void evalModel(
    const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const
    { TEST_FOR_EXCEPT(true); }
  // Static function that does the magic!
  static ModelEvaluatorBase::OutArgs<Scalar> createOutArgs(
    const ModelEvaluator<Scalar> &model,
    const SelectedDerivatives &fdDerivatives
    )
    {
      
      typedef ModelEvaluatorBase MEB;

      const MEB::OutArgs<Scalar> wrappedOutArgs = model.createOutArgs();
      const int Np = wrappedOutArgs.Np(), Ng = wrappedOutArgs.Ng();
      MEB::OutArgsSetup<Scalar> outArgs;

      outArgs.setModelEvalDescription(
        "DirectionalFiniteDiffCalculator: " + model.description()
        );

      // Start off by supporting everything that the underlying model supports
      // computing!

      outArgs.set_Np_Ng(Np,Ng);
      outArgs.setSupports(wrappedOutArgs);

      // Add support for finite difference DfDp(l) if asked

      const SelectedDerivatives::supports_DfDp_t
        &supports_DfDp = fdDerivatives.supports_DfDp_;
      for(
        SelectedDerivatives::supports_DfDp_t::const_iterator
          itr = supports_DfDp.begin();
        itr != supports_DfDp.end();
        ++itr
        )
      {
        const int l = *itr;
        assert_p_space(model,l);
        outArgs.setSupports(MEB::OUT_ARG_DfDp,l,MEB::DERIV_MV_BY_COL);
      }

      // Add support for finite difference DgDp(j,l) if asked
      
      const SelectedDerivatives::supports_DgDp_t
        &supports_DgDp = fdDerivatives.supports_DgDp_;
      for(
        SelectedDerivatives::supports_DgDp_t::const_iterator
          itr = supports_DgDp.begin();
        itr != supports_DgDp.end();
        ++itr
        )
      {
        const int j = itr->first;
        const int l = itr->second;
        assert_p_space(model,l);
        outArgs.setSupports(MEB::OUT_ARG_DgDp,j,l,MEB::DERIV_MV_BY_COL);
      }

      return outArgs;

    }

private:

  static void assert_p_space( const ModelEvaluator<Scalar> &model, const int l )
    {
#ifdef TEUCHOS_DEBUG
      const bool p_space_l_is_in_core = model.get_p_space(l)->hasInCoreView();
      TEST_FOR_EXCEPTION(
        !p_space_l_is_in_core, std::logic_error,
        "Error, for the model " << model.description()
        << ", the space p_space("<<l<<") must be in-core so that they can"
        " act as the domain space for the multi-vector derivative!"
        );
#endif
    }

};


} // namespace DirectionalFiniteDiffCalculatorTypes


// Private static data members


template<class Scalar>
const std::string DirectionalFiniteDiffCalculator<Scalar>::FDMethod_name = "FD Method";
template<class Scalar>
const RCP<
  Teuchos::StringToIntegralParameterEntryValidator<
    Thyra::DirectionalFiniteDiffCalculatorTypes::EFDMethodType
  >
>
DirectionalFiniteDiffCalculator<Scalar>::fdMethodValidator
= Teuchos::rcp(
  new Teuchos::StringToIntegralParameterEntryValidator<Thyra::DirectionalFiniteDiffCalculatorTypes::EFDMethodType>(
    Teuchos::tuple<std::string>(
      "order-one"
      ,"order-two"
      ,"order-two-central"
      ,"order-two-auto"
      ,"order-four"
      ,"order-four-central"
      ,"order-four-auto"
      )
    ,Teuchos::tuple<std::string>(
      "Use O(eps) one sided finite differences (cramped bounds)"
      ,"Use O(eps^2) one sided finite differences (cramped bounds)"
      ,"Use O(eps^2) two sided central finite differences"
      ,"Use \"order-two-central\" when not cramped by bounds, otherwise use \"order-two\""
      ,"Use O(eps^4) one sided finite differences (cramped bounds)"
      ,"Use O(eps^4) two sided central finite differences"
      ,"Use \"order-four-central\" when not cramped by bounds, otherwise use \"order-four\""
      )
    ,Teuchos::tuple<Thyra::DirectionalFiniteDiffCalculatorTypes::EFDMethodType>(
      Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_ONE
      ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_TWO
      ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_TWO_CENTRAL
      ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_TWO_AUTO
      ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_FOUR
      ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_FOUR_CENTRAL
      ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_ORDER_FOUR_AUTO
      )
    ,""
    )
  );
template<class Scalar>
const std::string DirectionalFiniteDiffCalculator<Scalar>::FDMethod_default = "order-one";

template<class Scalar>
const std::string DirectionalFiniteDiffCalculator<Scalar>::FDStepLength_name = "FD Step Length";
template<class Scalar>
const double DirectionalFiniteDiffCalculator<Scalar>::FDStepLength_default = -1.0;

template<class Scalar>
const std::string DirectionalFiniteDiffCalculator<Scalar>::FDStepSelectType_name = "FD Step Select Type";
template<class Scalar>
const RCP<
  Teuchos::StringToIntegralParameterEntryValidator<
    Thyra::DirectionalFiniteDiffCalculatorTypes::EFDStepSelectType
    >
  >
DirectionalFiniteDiffCalculator<Scalar>::fdStepSelectTypeValidator
= Teuchos::rcp(
  new Teuchos::StringToIntegralParameterEntryValidator<Thyra::DirectionalFiniteDiffCalculatorTypes::EFDStepSelectType>(
    Teuchos::tuple<std::string>(
      "Absolute"
      ,"Relative"
      )
    ,Teuchos::tuple<std::string>(
      "Use absolute step size \""+FDStepLength_name+"\""
      ,"Use relative step size \""+FDStepLength_name+"\"*||xo||inf"
      )
    ,Teuchos::tuple<Thyra::DirectionalFiniteDiffCalculatorTypes::EFDStepSelectType>(
      Thyra::DirectionalFiniteDiffCalculatorTypes::FD_STEP_ABSOLUTE
      ,Thyra::DirectionalFiniteDiffCalculatorTypes::FD_STEP_RELATIVE
      )
    ,""
    )
  );
template<class Scalar>
const std::string DirectionalFiniteDiffCalculator<Scalar>::FDStepSelectType_default = "Absolute";


// Constructors/initializer


template<class Scalar>
DirectionalFiniteDiffCalculator<Scalar>::DirectionalFiniteDiffCalculator(
  EFDMethodType               fd_method_type_in
  ,EFDStepSelectType          fd_step_select_type_in
  ,ScalarMag                  fd_step_size_in
  ,ScalarMag                  fd_step_size_min_in
  )
  :fd_method_type_(fd_method_type_in)
  ,fd_step_select_type_(fd_step_select_type_in)
  ,fd_step_size_(fd_step_size_in)
  ,fd_step_size_min_(fd_step_size_min_in)
{}


// Overriden from ParameterListAcceptor


template<class Scalar>
void DirectionalFiniteDiffCalculator<Scalar>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(paramList.get()==0);
  paramList->validateParameters(*getValidParameters());
  paramList_ = paramList;
  fd_method_type_ = fdMethodValidator->getIntegralValue(
    *paramList_,FDMethod_name,FDMethod_default);
  fd_step_select_type_ = fdStepSelectTypeValidator->getIntegralValue(
    *paramList_,FDStepSelectType_name,FDStepSelectType_default);
  fd_step_size_ = paramList_->get(
    FDStepLength_name,FDStepLength_default);
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
}


template<class Scalar>
RCP<ParameterList>
DirectionalFiniteDiffCalculator<Scalar>::getNonconstParameterList()
{
  return paramList_;
}


template<class Scalar>
RCP<ParameterList>
DirectionalFiniteDiffCalculator<Scalar>::unsetParameterList()
{
  RCP<ParameterList> _paramList = paramList_;
  paramList_ = Teuchos::null;
  return _paramList;
}


template<class Scalar>
RCP<const ParameterList>
DirectionalFiniteDiffCalculator<Scalar>::getParameterList() const
{
  return paramList_;
}


template<class Scalar>
RCP<const ParameterList>
DirectionalFiniteDiffCalculator<Scalar>::getValidParameters() const
{
  using Teuchos::rcp_implicit_cast;
  typedef Teuchos::ParameterEntryValidator PEV;
  static RCP<ParameterList> pl;
  if(pl.get()==NULL) {
    pl = Teuchos::parameterList();
    pl->set(
      FDMethod_name, FDMethod_default,
      "The method used to compute the finite differences.",
      rcp_implicit_cast<const PEV>(fdMethodValidator)
      );
    pl->set(
      FDStepSelectType_name,FDStepSelectType_default,
      "Method used to select the finite difference step length.",
      rcp_implicit_cast<const PEV>(fdStepSelectTypeValidator)
      );
    pl->set(
      FDStepLength_name,FDStepLength_default
      ,"The length of the finite difference step to take.\n"
      "A value of < 0.0 means that the step length will be determined automatically."
      );
    Teuchos::setupVerboseObjectSublist(&*pl);
  }
  return pl;
}


template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
DirectionalFiniteDiffCalculator<Scalar>::createOutArgs(
  const ModelEvaluator<Scalar> &model,
  const SelectedDerivatives   &fdDerivatives
  )
{
  return DirectionalFiniteDiffCalculatorTypes::OutArgsCreator<Scalar>::createOutArgs(
    model, fdDerivatives );
}


template<class Scalar>
void DirectionalFiniteDiffCalculator<Scalar>::calcVariations(
  const ModelEvaluator<Scalar> &model,
  const ModelEvaluatorBase::InArgs<Scalar> &bp, // basePoint
  const ModelEvaluatorBase::InArgs<Scalar> &dir, // directions
  const ModelEvaluatorBase::OutArgs<Scalar> &bfunc, // baseFunctionValues
  const ModelEvaluatorBase::OutArgs<Scalar> &var // variations
  ) const
{

  using std::string;

#ifdef THYRA_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR(
    string("Thyra::DirectionalFiniteDiffCalculator<")+ST::name()+">::calcVariations(...)"
    );
#endif
  
  using std::setw;
  using std::endl;
  using std::right;
  using Teuchos::as;
  typedef ModelEvaluatorBase MEB;
  namespace DFDCT = DirectionalFiniteDiffCalculatorTypes;
  typedef VectorBase<Scalar> V;
  typedef RCP<VectorBase<Scalar> > VectorPtr;
  
  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  const bool trace = (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_MEDIUM));
  Teuchos::OSTab tab(out);

  if(out.get() && trace)
    *out << "\nEntering DirectionalFiniteDiffCalculator<Scalar>::calcVariations(...)\n";

  if(out.get() && trace)
    *out
      << "\nbasePoint=\n" << describe(bp,verbLevel)
      << "\ndirections=\n" << describe(dir,verbLevel)
      << "\nbaseFunctionValues=\n" << describe(bfunc,verbLevel)
      << "\nvariations=\n" << describe(var,Teuchos::VERB_LOW);

#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPTION(
    var.isEmpty(), std::logic_error,
    "Error, all of the variations can not be null!"
    );
#endif
  
  //
  // To illustrate the theory behind this implementation consider
  // the generic multi-variable function h(z) <: R^n -> R.  Now let's
  // consider we have the base point zo and the vector v to
  // perturb h(z) along.  First form the function h(zo+a*v) and then
  // let's compute d(h)/d(a) at a = 0:
  // 
  // (1) d(h(zo+a*v))/d(a) at a = 0
  //         = sum( d(h)/d(x(i)) * d(x(i))/d(a), i = 1...n)
  //         = sum( d(h)/d(x(i)) * v(i), i = 1...n)
  //         = d(h)/d(a) * v
  //
  // Now we can approximate (1) using, for example, central differences as:
  // 
  // (2) d(h(zo+a*v))/d(a) at a = 0
  //          \approx ( h(zo+h*v) - h(zo+h*v) ) / (2*h)
  //
  // If we equate (1) and (2) we have the approximation:
  // 
  // (3) d(h)/d(a) * v \approx ( h(zo+h*v) - g(zo+h*v) ) / (2*h)
  // 
  // 

  // /////////////////////////////////////////
  // Validate the input

  // ToDo: Validate input!

  switch(this->fd_method_type()) {
    case DFDCT::FD_ORDER_ONE:
      if(out.get()&&trace) *out<<"\nUsing one-sided, first-order finite differences ...\n";
      break;
    case DFDCT::FD_ORDER_TWO:
      if(out.get()&&trace) *out<<"\nUsing one-sided, second-order finite differences ...\n";
      break;
    case DFDCT::FD_ORDER_TWO_CENTRAL:
      if(out.get()&&trace) *out<<"\nUsing second-order central finite differences ...\n";
      break;
    case DFDCT::FD_ORDER_TWO_AUTO:
      if(out.get()&&trace) *out<<"\nUsing auto selection of some second-order finite difference method ...\n";
      break;
    case DFDCT::FD_ORDER_FOUR:
      if(out.get()&&trace) *out<<"\nUsing one-sided, fourth-order finite differences ...\n";
      break;
    case DFDCT::FD_ORDER_FOUR_CENTRAL:
      if(out.get()&&trace) *out<<"\nUsing fourth-order central finite differences ...\n";
      break;
    case DFDCT::FD_ORDER_FOUR_AUTO:
      if(out.get()&&trace) *out<<"\nUsing auto selection of some fourth-order finite difference method ...\n";
      break;
    default:
      TEST_FOR_EXCEPT(true); // Should not get here!
  }

  // ////////////////////////
  // Find the step size

  //
  // Get defaults for the optimal step sizes
  //

  const ScalarMag
    sqrt_epsilon = SMT::squareroot(SMT::eps()),
    u_optimal_1  = sqrt_epsilon,
    u_optimal_2  = SMT::squareroot(sqrt_epsilon),
    u_optimal_4  = SMT::squareroot(u_optimal_2);

  ScalarMag
    bp_norm = SMT::zero();
  // ToDo above: compute a reasonable norm somehow based on the base-point vector(s)!
  
  ScalarMag
    uh_opt = 0.0;
  switch(this->fd_method_type()) {
    case DFDCT::FD_ORDER_ONE:
      uh_opt = u_optimal_1 * ( fd_step_select_type() == DFDCT::FD_STEP_ABSOLUTE ? 1.0 : bp_norm + 1.0 );
      break;
    case DFDCT::FD_ORDER_TWO:
    case DFDCT::FD_ORDER_TWO_CENTRAL:
    case DFDCT::FD_ORDER_TWO_AUTO:
      uh_opt = u_optimal_2 * ( fd_step_select_type() == DFDCT::FD_STEP_ABSOLUTE ? 1.0 : bp_norm + 1.0 );
      break;
    case DFDCT::FD_ORDER_FOUR:
    case DFDCT::FD_ORDER_FOUR_CENTRAL:
    case DFDCT::FD_ORDER_FOUR_AUTO:
      uh_opt = u_optimal_4 * ( fd_step_select_type() == DFDCT::FD_STEP_ABSOLUTE ? 1.0 : bp_norm + 1.0 );
      break;
    default:
      TEST_FOR_EXCEPT(true); // Should not get here!
  }

  if(out.get()&&trace) *out<<"\nDefault optimal step length uh_opt = " << uh_opt << " ...\n";

  //
  // Set the step sizes used.
  //

  ScalarMag
    uh      = this->fd_step_size();

  if( uh < 0 )
    uh = uh_opt;
  else if( fd_step_select_type() == DFDCT::FD_STEP_RELATIVE )
    uh *= (bp_norm + 1.0);

  if(out.get()&&trace) *out<<"\nStep size to be used uh="<<uh<<"\n";

  //
  // Find step lengths that stay in bounds!
  //

  // ToDo: Add logic for variable bounds when needed!

  //
  // Set the actual method being used
  //
  // ToDo: Consider cramped bounds and method order.
  //
  
  DFDCT::EFDMethodType l_fd_method_type = this->fd_method_type();
  switch(l_fd_method_type) {
    case DFDCT::FD_ORDER_TWO_AUTO:
      l_fd_method_type = DFDCT::FD_ORDER_TWO_CENTRAL;
      break;
    case DFDCT::FD_ORDER_FOUR_AUTO:
      l_fd_method_type = DFDCT::FD_ORDER_FOUR_CENTRAL;
      break;
    default:
      break; // Okay
  }

  //if(out.get()&&trace) *out<<"\nStep size to fit in bounds: uh="<<uh"\n";

  int p_saved = -1;
  if(out.get())
    p_saved = out->precision();

  // ///////////////////////////////////////////////
  // Compute the directional derivatives

  const int Np = var.Np(), Ng = var.Ng();
  
  // Setup storage for perturbed variables
  VectorPtr per_x, per_x_dot;
  std::vector<VectorPtr> per_p(Np);
  MEB::InArgs<Scalar> pp = model.createInArgs();
  pp.setArgs(bp); // Set all args to start with
  if( bp.supports(MEB::IN_ARG_x) ) {
    if( dir.get_x().get() )
      pp.set_x(per_x=createMember(model.get_x_space()));
    else
      pp.set_x(bp.get_x());
  }
  if( bp.supports(MEB::IN_ARG_x_dot) ) {
    if( dir.get_x_dot().get() )
      pp.set_x_dot(per_x_dot=createMember(model.get_x_space()));
    else
      pp.set_x_dot(bp.get_x_dot());
  }
  for( int l = 0; l < Np; ++l ) {
    if( dir.get_p(l).get() )
      pp.set_p(l,per_p[l]=createMember(model.get_p_space(l)));
    else
      pp.set_p(l,bp.get_p(l));
  }
  if(out.get() && trace)
    *out
      << "\nperturbedPoint after initial setup (with some unintialized vectors) =\n"
      << describe(pp,verbLevel);
  
  // Setup storage for perturbed functions
  bool all_funcs_at_base_computed = true;
  MEB::OutArgs<Scalar> pfunc = model.createOutArgs();
  {
    VectorPtr f;
    if( var.supports(MEB::OUT_ARG_f) && (f=var.get_f()).get() ) {
      pfunc.set_f(createMember(model.get_f_space()));
      assign(f.ptr(),ST::zero());
      if(!bfunc.get_f().get()) all_funcs_at_base_computed = false;
    }
    for( int j = 0; j < Ng; ++j ) {
      VectorPtr g_j;
      if( (g_j=var.get_g(j)).get() ) {
        pfunc.set_g(j,createMember(model.get_g_space(j)));
        assign(g_j.ptr(),ST::zero());
        if(!bfunc.get_g(j).get()) all_funcs_at_base_computed = false;
      }
    }
  }
  if(out.get() && trace)
    *out
      << "\nperturbedFunctions after initial setup (with some unintialized vectors) =\n"
      << describe(pfunc,verbLevel);
  
  const int dbl_p = 15;
  if(out.get())
    *out << std::setprecision(dbl_p);
    
  //
  // Compute the weighted sum of the terms
  //
    
  int num_evals = 0;
  ScalarMag dwgt = SMT::zero();
  switch(l_fd_method_type) {
    case DFDCT::FD_ORDER_ONE: // may only need one eval if f(xo) etc is passed in
      num_evals = 2;
      dwgt      = ScalarMag(1.0);
      break;
    case DFDCT::FD_ORDER_TWO: // may only need two evals if c(xo) etc is passed in
      num_evals = 3;
      dwgt      = ScalarMag(2.0);
      break;
    case DFDCT::FD_ORDER_TWO_CENTRAL:
      num_evals = 2;
      dwgt      = ScalarMag(2.0);
      break;
    case DFDCT::FD_ORDER_FOUR:
      num_evals = 5;
      dwgt      = ScalarMag(12.0);
      break;
    case DFDCT::FD_ORDER_FOUR_CENTRAL:
      num_evals = 4;
      dwgt      = ScalarMag(12.0);
      break;
    default:
      TEST_FOR_EXCEPT(true); // Should not get here!
  }
  for( int eval_i = 1; eval_i <= num_evals; ++eval_i ) {
    // Set the step constant and the weighting constant
    ScalarMag
      uh_i   = 0.0,
      wgt_i  = 0.0;
    switch(l_fd_method_type) {
      case DFDCT::FD_ORDER_ONE: {
        switch(eval_i) {
          case 1:
            uh_i  = +0.0;
            wgt_i = -1.0;
            break;
          case 2:
            uh_i  = +1.0;
            wgt_i = +1.0;
            break;
        }
        break;
      }
      case DFDCT::FD_ORDER_TWO: {
        switch(eval_i) {
          case 1:
            uh_i  = +0.0;
            wgt_i = -3.0;
            break;
          case 2:
            uh_i  = +1.0;
            wgt_i = +4.0;
            break;
          case 3:
            uh_i  = +2.0;
            wgt_i = -1.0;
            break;
        }
        break;
      }
      case DFDCT::FD_ORDER_TWO_CENTRAL: {
        switch(eval_i) {
          case 1:
            uh_i  = -1.0;
            wgt_i = -1.0;
            break;
          case 2:
            uh_i  = +1.0;
            wgt_i = +1.0;
            break;
        }
        break;
      }
      case DFDCT::FD_ORDER_FOUR: {
        switch(eval_i) {
          case 1:
            uh_i  = +0.0;
            wgt_i = -25.0;
            break;
          case 2:
            uh_i  = +1.0;
            wgt_i = +48.0;
            break;
          case 3:
            uh_i  = +2.0;
            wgt_i = -36.0;
            break;
          case 4:
            uh_i  = +3.0;
            wgt_i = +16.0;
            break;
          case 5:
            uh_i  = +4.0;
            wgt_i = -3.0;
            break;
        }
        break;
      }
      case DFDCT::FD_ORDER_FOUR_CENTRAL: {
        switch(eval_i) {
          case 1:
            uh_i  = -2.0;
            wgt_i = +1.0;
            break;
          case 2:
            uh_i  = -1.0;
            wgt_i = -8.0;
            break;
          case 3:
            uh_i  = +1.0;
            wgt_i = +8.0;
            break;
          case 4:
            uh_i  = +2.0;
            wgt_i = -1.0;
            break;
        }
        break;
      }
      case DFDCT::FD_ORDER_TWO_AUTO:
      case DFDCT::FD_ORDER_FOUR_AUTO:
        break; // Okay
      default:
        TEST_FOR_EXCEPT(true);
    }

    if(out.get() && trace)
      *out << "\neval_i="<<eval_i<<", uh_i="<<uh_i<<", wgt_i="<<wgt_i<<"\n";
    Teuchos::OSTab tab2(out);
    
    // Compute the weighted term and add it to the sum
    if(uh_i == ST::zero()) {
      MEB::OutArgs<Scalar> bfuncall;
      if(!all_funcs_at_base_computed) {
        // Compute missing functions at the base point
        bfuncall = model.createOutArgs();
        if( pfunc.supports(MEB::OUT_ARG_f) && pfunc.get_f().get() && !bfunc.get_f().get() ) {
          bfuncall.set_f(createMember(model.get_f_space()));
        }
        for( int j = 0; j < Ng; ++j ) {
          if( pfunc.get_g(j).get() && !bfunc.get_g(j).get() ) {
            bfuncall.set_g(j,createMember(model.get_g_space(j)));
          }
        }
        model.evalModel(bp,bfuncall);
        bfuncall.setArgs(bfunc);
      }
      else {
        bfuncall = bfunc;
      }
      // Use functions at the base point
      if(out.get() && trace)
        *out << "\nSetting variations = wgt_i * basePoint ...\n";
      VectorPtr f;
      if( pfunc.supports(MEB::OUT_ARG_f) && (f=var.get_f()).get() ) {
        V_StV<Scalar>(f.ptr(), wgt_i, *bfuncall.get_f());
      }
      for( int j = 0; j < Ng; ++j ) {
        VectorPtr g_j;
        if( (g_j=var.get_g(j)).get() ) {
          V_StV<Scalar>(g_j.ptr(), wgt_i, *bfuncall.get_g(j));
        }
      }
    }
    else {
      if(out.get() && trace)
        *out << "\nSetting perturbedPoint = basePoint + uh_i*uh*direction ...\n";
      // z = zo + uh_i*uh*v
      {
        if ( dir.supports(MEB::IN_ARG_x) && dir.get_x().get() )
          V_StVpV(per_x.ptr(),as<Scalar>(uh_i*uh),*dir.get_x(),*bp.get_x());
        if ( dir.supports(MEB::IN_ARG_x_dot) && dir.get_x_dot().get() )
          V_StVpV(per_x_dot.ptr(), as<Scalar>(uh_i*uh), *dir.get_x_dot(), *bp.get_x_dot());
        for ( int l = 0; l < Np; ++l ) {
          if( dir.get_p(l).get() )
            V_StVpV(per_p[l].ptr(), as<Scalar>(uh_i*uh), *dir.get_p(l), *bp.get_p(l));
        }
      }
      if(out.get() && trace)
        *out << "\nperturbedPoint=\n" << describe(pp,verbLevel);
      // Compute perturbed function values h(zo+uh_i*uh)
      if(out.get() && trace)
        *out << "\nCompute perturbedFunctions at perturbedPoint...\n";
      model.evalModel(pp,pfunc);
      if(out.get() && trace)
        *out << "\nperturbedFunctions=\n" << describe(pfunc,verbLevel);
      // Sum perturbed function values into total variation
      {
        // var_h += wgt_i*perturbed_h
        if(out.get() && trace)
          *out << "\nComputing variations += wgt_i*perturbedfunctions ...\n";
        VectorPtr f;
        if( pfunc.supports(MEB::OUT_ARG_f) && (f=pfunc.get_f()).get() )
          Vp_StV<Scalar>(var.get_f().ptr(), wgt_i, *f);
        for( int j = 0; j < Ng; ++j ) {
          VectorPtr g_j;
          if( (g_j=pfunc.get_g(j)).get() )
            Vp_StV<Scalar>(var.get_g(j).ptr(), wgt_i, *g_j);
        }
      }
    }
    if(out.get() && trace)
      *out << "\nvariations=\n" << describe(var,verbLevel);
  }

  //
  // Multiply by the scaling factor!
  //
  
  {
    // var_h *= 1.0/(dwgt*uh)
    const Scalar alpha = ST::one()/(dwgt*uh);
    if(out.get() && trace)
      *out
        << "\nComputing variations *= (1.0)/(dwgt*uh),"
        << " where (1.0)/(dwgt*uh) = (1.0)/("<<dwgt<<"*"<<uh<<") = "<<alpha<<" ...\n";
    VectorPtr f;
    if( var.supports(MEB::OUT_ARG_f) && (f=var.get_f()).get() )
      Vt_S(f.ptr(),alpha);
    for( int j = 0; j < Ng; ++j ) {
      VectorPtr g_j;
      if( (g_j=var.get_g(j)).get() )
        Vt_S(g_j.ptr(),alpha);
    }
    if(out.get() && trace)
      *out << "\nFinal variations=\n" << describe(var,verbLevel);
  }
  
  if(out.get())
    *out << std::setprecision(p_saved);

  if(out.get() && trace)
    *out << "\nLeaving DirectionalFiniteDiffCalculator<Scalar>::calcVariations(...)\n";
  
}


template<class Scalar>
void DirectionalFiniteDiffCalculator<Scalar>::calcDerivatives(
  const ModelEvaluator<Scalar> &model,
  const ModelEvaluatorBase::InArgs<Scalar> &bp, // basePoint
  const ModelEvaluatorBase::OutArgs<Scalar> &bfunc, // baseFunctionValues
  const ModelEvaluatorBase::OutArgs<Scalar> &deriv // derivatives
  ) const
{

  using std::string;
  //typedef Teuchos::ScalarTraits<Scalar> ST;

#ifdef THYRA_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR(
    string("Thyra::DirectionalFiniteDiffCalculator<")+ST::name()+">::calcDerivatives(...)"
    );
#endif

  typedef ModelEvaluatorBase MEB;
  typedef RCP<VectorBase<Scalar> > VectorPtr;
  typedef RCP<MultiVectorBase<Scalar> > MultiVectorPtr;

  RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  const bool trace = (static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_MEDIUM));
  Teuchos::OSTab tab(out);

  if(out.get() && trace)
    *out << "\nEntering DirectionalFiniteDiffCalculator<Scalar>::calcDerivatives(...)\n";

  if(out.get() && trace)
    *out
      << "\nbasePoint=\n" << describe(bp,verbLevel)
      << "\nbaseFunctionValues=\n" << describe(bfunc,verbLevel)
      << "\nderivatives=\n" << describe(deriv,Teuchos::VERB_LOW);
  
  //
  // We will only compute finite differences w.r.t. p(l) for now
  //
  const int Np = bp.Np(), Ng = bfunc.Ng();
  MEB::InArgs<Scalar> dir = model.createInArgs();
  MEB::OutArgs<Scalar> var = model.createOutArgs();
  MultiVectorPtr DfDp_l;
  std::vector<MEB::DerivativeMultiVector<Scalar> > DgDp_l(Ng);
  std::vector<VectorPtr> var_g(Ng);
  for( int l = 0; l < Np; ++l ) {
    if(out.get() && trace)
      *out << "\nComputing derivatives for parameter subvector p("<<l<<") ...\n";
    Teuchos::OSTab tab2(out);
    //
    // Load up OutArgs var object of derivative variations to compute
    //
    bool hasDerivObject = false;
    // DfDp(l)
    if (
      !deriv.supports(MEB::OUT_ARG_DfDp,l).none()
      && !deriv.get_DfDp(l).isEmpty()
      )
    {
      hasDerivObject = true;
      std::ostringstream name; name << "DfDp("<<l<<")";
      DfDp_l = get_mv(deriv.get_DfDp(l),name.str(),MEB::DERIV_MV_BY_COL);
    }
    else {
      DfDp_l = Teuchos::null;
    }
    // DgDp(j=1...Ng,l)
    for ( int j = 0; j < Ng; ++j ) {
      if (
        !deriv.supports(MEB::OUT_ARG_DgDp,j,l).none()
        &&
        !deriv.get_DgDp(j,l).isEmpty()
        )
      {
        hasDerivObject = true;
        std::ostringstream name; name << "DgDp("<<j<<","<<l<<")";
        DgDp_l[j] = get_dmv(deriv.get_DgDp(j,l),name.str());
        if( DgDp_l[j].getMultiVector().get() && !var_g[j].get() )
        {
          // Need a temporary vector for the variation for g(j)
          var_g[j] = createMember(model.get_g_space(j));
        }
      }
      else{
        DgDp_l[j] = MEB::DerivativeMultiVector<Scalar>();
        var_g[j] = Teuchos::null;
      }
    }
    //
    // Compute derivative variations by finite differences
    //
    if (hasDerivObject) {
      VectorPtr e_i = createMember(model.get_p_space(l));
      dir.set_p(l,e_i);
      assign(e_i.ptr(),ST::zero());
      const int np_l = e_i->space()->dim();
      for( int i = 0 ; i < np_l; ++ i ) {
        if(out.get() && trace)
          *out << "\nComputing derivatives for single variable p("<<l<<")("<<i<<") ...\n";
        Teuchos::OSTab tab3(out);
        if(DfDp_l.get()) var.set_f(DfDp_l->col(i)); // Compute DfDp(l)(i) in place!
        for(int j = 0; j < Ng; ++j) {
          MultiVectorPtr DgDp_j_l;
          if( (DgDp_j_l=DgDp_l[j].getMultiVector()).get() ) {
            var.set_g(j,var_g[j]); // Computes d(g(j))/d(p(l)(i))
          }
        }
        set_ele(i,ST::one(),e_i.ptr());
        this->calcVariations(
          model,bp,dir,bfunc,var
          );
        set_ele(i,ST::zero(),e_i.ptr());
        if (DfDp_l.get()) var.set_f(Teuchos::null);
        for (int j = 0; j < Ng; ++j) {
          MultiVectorPtr DgDp_j_l;
          if ( !is_null(DgDp_j_l=DgDp_l[j].getMultiVector()) ) {
            assign( DgDp_j_l->col(i).ptr(), *var_g[j] );
          }
        }
      }
      dir.set_p(l,Teuchos::null);
    }
  }
    
  if(out.get() && trace)
    *out
      << "\nderivatives=\n" << describe(deriv,verbLevel);
  
  if(out.get() && trace)
    *out << "\nLeaving DirectionalFiniteDiffCalculator<Scalar>::calcDerivatives(...)\n";

}


} // namespace Thyra


#endif	// THYRA_DIRECTIONAL_FINITE_DIFF_CALCULATOR_DEF_HPP
