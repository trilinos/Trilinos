#ifndef __Step01_LinearFunction_impl_hpp__
#define __Step01_LinearFunction_impl_hpp__

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"

namespace user_app {

//**********************************************************************
template <typename EvalT,typename Traits>
LinearFunction<EvalT,Traits>::LinearFunction(const std::string & name,
                                             double acoeff,double bcoeff,
                                             const panzer::IntegrationRule & ir)
  : acoeff_(acoeff) 
  , bcoeff_(bcoeff) 
  , id_(ir)
{
  using Teuchos::RCP;

  Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_scalar;

  result = PHX::MDField<ScalarT,panzer::Cell,panzer::Point>(name, data_layout);
  this->addEvaluatedField(result);

  this->setName("Linear Function("+name+")");
}

//**********************************************************************
template <typename EvalT,typename Traits>
void LinearFunction<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  const auto & ip_coordinates = workset.getIntegrationValues(id_).ip_coordinates;
  for (panzer::index_t cell = 0; cell < workset.numCells(); ++cell) {
    for (int point = 0; point < result.extent_int(1); ++point) {

      const double& x = ip_coordinates(cell,point,0);
      const double& y = ip_coordinates(cell,point,1);

      result(cell,point) = acoeff_*x + acoeff_*y + bcoeff_;
    }
  }
}

//**********************************************************************
}

#endif
