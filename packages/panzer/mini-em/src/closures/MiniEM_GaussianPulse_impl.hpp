#ifndef MINIEM_GAUSSIAN_PULSE_IMPL_HPP
#define MINIEM_GAUSSIAN_PULSE_IMPL_HPP

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"

namespace mini_em {

//**********************************************************************
template <typename EvalT,typename Traits>
GaussianPulse<EvalT,Traits>::GaussianPulse(const std::string & name,
                                           const panzer::IntegrationRule & ir,
                                           const panzer::FieldLayoutLibrary & fl,
                                           const double & dt)
{
  using Teuchos::RCP;

  Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_vector;
  ir_degree = ir.cubature_degree;
  ir_dim = ir.spatial_dimension;

  current = PHX::MDField<ScalarT,Cell,Point,Dim>(name, data_layout);
  this->addEvaluatedField(current);
  
  alpha = 1.0;
  beta  = 5.0*dt;

  Teuchos::RCP<const panzer::PureBasis> basis = fl.lookupBasis("E_edge");
  const std::string coordName = panzer::GatherBasisCoordinates<EvalT,Traits>::fieldName(basis->name());
  coords = PHX::MDField<const ScalarT,Cell,Point,Dim>(coordName, basis->coordinates); 
  this->addDependentField(coords);

  std::string n = "Gaussian Pulse";
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void GaussianPulse<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{ 
  using panzer::index_t;

  double time = workset.time;

  const ScalarT factor = std::exp(-(time-2.0*beta)*(time-2.0*beta)/beta/beta);
  const ScalarT scale = 1.0/alpha/alpha;
  if (ir_dim == 3) {
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int point = 0; point < current.extent_int(1); ++point) {
        const ScalarT& x = coords(cell,point,0);
        const ScalarT& y = coords(cell,point,1);
        const ScalarT& z = coords(cell,point,2);
        const ScalarT  r2 = (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5);
        current(cell,point,0) = 0.0;
        current(cell,point,1) = 0.0;
        current(cell,point,2) = std::exp(-r2*scale)*factor;
      }
    }
  } else {
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int point = 0; point < current.extent_int(1); ++point) {
        const ScalarT& x = coords(cell,point,0);
        const ScalarT& y = coords(cell,point,1);
        const ScalarT  r2 = (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5);
        current(cell,point,0) = 0.0;
        current(cell,point,1) = std::exp(-r2*scale)*factor;
      }
    }
  }
}

//**********************************************************************
}

#endif
