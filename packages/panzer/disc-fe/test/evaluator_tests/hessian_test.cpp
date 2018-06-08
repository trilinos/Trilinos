// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include "Phalanx_FieldManager.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

namespace panzer {

/** Seed to compute a second directional derivative.This is used
  * to compute (given a function f)
  *
  *    H(x) * v = d/dt(\nabla f(x+t*v))|_{t=0}
  *       
  * where x is the value to differentiate about, v is the direction. This function
  * seeds for one value of the independent variable x[i] in the direction v[i]. What
  * is returned from this function is the independent variable to pass the derivative
  * to.
  *
  * \param[in] num_vars Size of the vector x (and v)
  * \param[in] index Index of variable that this is seeding
  * \param[in] xi Value of the x vector to seed with
  * \param[in] vi Value of the v vector to seed with
  */
inline panzer::Traits::HessianType seed_second_deriv(int num_vars, int index, double /* xi */, double vi)
{
  typedef panzer::Traits::HessianType SecondFadType;

//  SecondFadType x = SecondFadType(1,panzer::Traits::FadType(num_vars,index,xi));
//  x.fastAccessDx(0) = vi;

  Sacado::Fad::SFad<panzer::Traits::RealType,1> xi_fad;
  xi_fad.fastAccessDx(0) = vi;
  SecondFadType x = SecondFadType(num_vars,index,xi_fad);

  return x;
}

//**********************************************************************
template<typename EvalT, typename Traits>
class InputConditionsEvaluator
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    InputConditionsEvaluator(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

public:
  PHX::MDField<ScalarT,panzer::IP> x;
  PHX::MDField<ScalarT,panzer::IP> y;
  PHX::MDField<ScalarT,panzer::IP> dx;
  PHX::MDField<ScalarT,panzer::IP> dy;

}; // end of class InputConditionsEvaluator


//**********************************************************************
template<typename EvalT, typename Traits>
InputConditionsEvaluator<EvalT, Traits>::
InputConditionsEvaluator(
  const Teuchos::ParameterList&  /* p */)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Read from parameters
  const std::string x_name = "X";
  const std::string y_name = "Y";

  RCP<PHX::DataLayout> dl = rcp(new PHX::MDALayout<IP>(5));

  // grab information from quadrature rule
  x = PHX::MDField<ScalarT,IP>(x_name, dl);
  y = PHX::MDField<ScalarT,IP>(y_name, dl);
  dx = PHX::MDField<ScalarT,IP>("d"+x_name, dl);
  dy = PHX::MDField<ScalarT,IP>("d"+y_name, dl);

  this->addEvaluatedField(x);
  this->addEvaluatedField(y);
  this->addEvaluatedField(dx);
  this->addEvaluatedField(dy);
  
  std::string n = "InputConditions evaluator (" + PHX::typeAsString<ScalarT>()+")";
  this->setName(n);
}

//**********************************************************************

template<typename EvalT, typename Traits>
void
InputConditionsEvaluator<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData  /* sd */,
  PHX::FieldManager<Traits>&  fm)
{
  this->utils.setFieldData(x,fm);
  this->utils.setFieldData(y,fm);
  this->utils.setFieldData(dx,fm);
  this->utils.setFieldData(dy,fm);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
InputConditionsEvaluator<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  /* workset */)
{ 
  double x_val = 0.25;
  double y_val = 0.5;
  double dx_val = 2.0;
  double dy_val = 3.0;

  for(int i=0;i<5;i++) {
    dx(i) = ScalarT(dx_val);
    dy(i) = ScalarT(dy_val);
    x(i) = seed_second_deriv(2,0,x_val,dx_val);
    y(i) = seed_second_deriv(2,1,y_val,dy_val);
  }
}

//**********************************************************************
template<typename EvalT, typename Traits>
class HessianTestEvaluator
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    HessianTestEvaluator(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

public:
  PHX::MDField<const ScalarT,panzer::IP> x;
  PHX::MDField<const ScalarT,panzer::IP> y;
  PHX::MDField<ScalarT,panzer::IP> result;

}; // end of class HessianTestEvaluator


//**********************************************************************
template<typename EvalT, typename Traits>
HessianTestEvaluator<EvalT, Traits>::
HessianTestEvaluator(
  const Teuchos::ParameterList&  /* p */)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Read from parameters
  const std::string x_name = "X";
  const std::string y_name = "Y";
  const std::string result_name = "Result";

  RCP<PHX::DataLayout> dl = rcp(new PHX::MDALayout<IP>(5));

  // grab information from quadrature rule
  x = PHX::MDField<const ScalarT,IP>(x_name, dl);
  y = PHX::MDField<const ScalarT,IP>(y_name, dl);
  result = PHX::MDField<ScalarT,IP>(result_name, dl);

  this->addDependentField(x);
  this->addDependentField(y);
  this->addEvaluatedField(result);
  
  std::string n = "Hessian test evaluator (" + PHX::typeAsString<ScalarT>()+")";
  this->setName(n);
}

//**********************************************************************

template<typename EvalT, typename Traits>
void
HessianTestEvaluator<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData  /* sd */,
  PHX::FieldManager<Traits>&  fm)
{
  this->utils.setFieldData(x,fm);
  this->utils.setFieldData(y,fm);
  this->utils.setFieldData(result,fm);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
HessianTestEvaluator<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  /* workset */)
{ 
  // Grad = y * std::cos(x*y)
  //      = x * std::cos(x*y)-0.25*std::sin(y)
  
  // Hess*v = dy * std::cos(x*y) - y * ( dx * sin(x*y) + dy * sin(x*y))
  //        = dx * std::cos(x*y) - x * ( dx * sin(x*y) + dy * sin(x*y)) - 0.25 * dy * std::cos(y)
  //
  for(int ip=0;ip<Teuchos::as<int>(result.dimension_0());++ip)
    result(ip) = std::sin(x(ip)*y(ip))+0.25*std::cos(y(ip));
}

//**********************************************************************
double func(double x, double y)
{
  return std::sin(x*y)+0.25*std::cos(y);
}

std::vector<double> hess_func(double x, double y,double dx, double dy)
{
  std::vector<double> v(2);
  v[0] = dy * std::cos(x*y) - y * ( y*dx * sin(x*y) + x*dy * sin(x*y));
  v[1] = dx * std::cos(x*y) - x * ( y*dx * sin(x*y) + x*dy * sin(x*y)) - 0.25 * dy * std::cos(y);

  return v;
}

//**********************************************************************
TEUCHOS_UNIT_TEST(hessian_test,correctness)
{
  typedef InputConditionsEvaluator<panzer::Traits::Hessian,panzer::Traits> InputCondEval;
  typedef HessianTestEvaluator<panzer::Traits::Hessian,panzer::Traits> HessTestEval;
  typedef panzer::Traits::HessianType ScalarT;
  typedef Sacado::ScalarValue<ScalarT> Value;
 
  using Teuchos::RCP;
  using Teuchos::rcp;


  // the one and only evaluator
  Teuchos::ParameterList empty_pl;
  RCP<InputCondEval> ic_eval = rcp(new InputCondEval(empty_pl));
  RCP<HessTestEval> ht_eval = rcp(new HessTestEval(empty_pl));

  Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm
     = Teuchos::rcp(new PHX::FieldManager<panzer::Traits>); 
  fm->registerEvaluator<panzer::Traits::Hessian>(ic_eval);
  fm->registerEvaluator<panzer::Traits::Hessian>(ht_eval);
  fm->requireField<panzer::Traits::Hessian>(ht_eval->result.fieldTag());

  std::vector<PHX::index_size_type> derivative_dimensions;
  derivative_dimensions.push_back(4);
  fm->setKokkosExtendedDataTypeDimensions<panzer::Traits::Hessian>(derivative_dimensions);

  panzer::Traits::SD setupData;
  fm->postRegistrationSetup(setupData);

  panzer::Workset workset;
  panzer::Traits::PED preEvalData;

  fm->preEvaluate<panzer::Traits::Hessian>(preEvalData);
  fm->evaluateFields<panzer::Traits::Hessian>(workset);
  fm->postEvaluate<panzer::Traits::Hessian>(0);

  for(int i=0;i<5;i++) {
    double x  = Value::eval(ic_eval->x(i));
    double y  = Value::eval(ic_eval->y(i));
    double dx = Value::eval(ic_eval->dx(i));
    double dy = Value::eval(ic_eval->dy(i));
    double f = func(x,y);
    std::vector<double> hess = hess_func(x,y,dx,dy);

    ScalarT r = ht_eval->result(i);

    TEST_EQUALITY(Value::eval(r),f);
    TEST_EQUALITY(r.fastAccessDx(0).fastAccessDx(0),hess[0]);
    TEST_EQUALITY(r.fastAccessDx(1).fastAccessDx(0),hess[1]);
  }
}


}
