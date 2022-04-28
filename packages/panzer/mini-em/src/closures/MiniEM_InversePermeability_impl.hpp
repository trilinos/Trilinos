#ifndef MINIEM_INVERSEPERMEABILITY_IMPL_HPP
#define MINIEM_INVERSEPERMEABILITY_IMPL_HPP

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"

namespace mini_em {

//**********************************************************************
template <typename EvalT,typename Traits>
InversePermeability<EvalT,Traits>::InversePermeability(const std::string & name,
                                                       const panzer::IntegrationRule & ir,
                                                       const panzer::FieldLayoutLibrary & fl,
                                                       const double & mu_,
                                                       const std::string& DoF_)
{
  using Teuchos::RCP;

  Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_scalar;
  ir_degree = ir.cubature_degree;

  permeability = PHX::MDField<ScalarT,Cell,Point>(name, data_layout);
  this->addEvaluatedField(permeability);
  
  mu = mu_;

  Teuchos::RCP<const panzer::PureBasis> basis = fl.lookupBasis(DoF_);
  const std::string coordName = panzer::GatherBasisCoordinates<EvalT,Traits>::fieldName(basis->name());
  coords = PHX::MDField<const ScalarT,Cell,Point,Dim>(coordName, basis->coordinates); 
  this->addDependentField(coords);

  std::string n = "InversePermeability";
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void InversePermeability<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{ 
  using panzer::index_t;
  auto inv_mu = 1./mu;
  Kokkos::MDRangePolicy<PHX::exec_space,Kokkos::Rank<2>> policy({0,0},{workset.num_cells,permeability.extent_int(1)});
  Kokkos::parallel_for("panzer:InversePermeability",policy,KOKKOS_LAMBDA (const int cell,const int point) {
      // const ScalarT& x = coords(cell,point,0);
      // const ScalarT& y = coords(cell,point,1);
      // const ScalarT& z = coords(cell,point,2);
      permeability(cell,point) = inv_mu;
    });
  Kokkos::fence();
}

//**********************************************************************
}

#endif
