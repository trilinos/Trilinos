#ifndef MINIEM_TENSORCONDUCTIVITY_IMPL_HPP
#define MINIEM_TENSORCONDUCTIVITY_IMPL_HPP

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"

namespace mini_em {

//**********************************************************************
template <typename EvalT,typename Traits>
TensorConductivity<EvalT,Traits>::TensorConductivity(const std::string & name,
                                         const panzer::IntegrationRule & ir,
                                         const panzer::FieldLayoutLibrary & fl,
                                         const double & sigma_,
                                         const std::string& DoF_)
{
  using Teuchos::RCP;

  Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_tensor;
  ir_degree = ir.cubature_degree;
  ir_dim = ir.spatial_dimension;

  conductivity = PHX::MDField<ScalarT,Cell,Point,Dim,Dim>(name, data_layout);
  this->addEvaluatedField(conductivity);

  sigma = sigma_;

  Teuchos::RCP<const panzer::PureBasis> basis = fl.lookupBasis(DoF_);
  const std::string coordName = panzer::GatherBasisCoordinates<EvalT,Traits>::fieldName(basis->name());
  coords = PHX::MDField<const ScalarT,Cell,Point,Dim>(coordName, basis->coordinates);
  this->addDependentField(coords);

  std::string n = "TensorConductivity";
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void TensorConductivity<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  using panzer::index_t;

  if (ir_dim == 3) {
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int point = 0; point < conductivity.extent_int(1); ++point) {
        // const ScalarT& x = coords(cell,point,0);
        // const ScalarT& y = coords(cell,point,1);
        // const ScalarT& z = coords(cell,point,2);
        conductivity(cell,point,0,0) = sigma;
        conductivity(cell,point,0,1) = 0.;
        conductivity(cell,point,0,2) = 0.;

        conductivity(cell,point,1,0) = 0.;
        conductivity(cell,point,1,1) = sigma;
        conductivity(cell,point,1,2) = 0.;

        conductivity(cell,point,2,0) = 0.;
        conductivity(cell,point,2,1) = 0.;
        conductivity(cell,point,2,2) = sigma;
      }
    }
  } else {
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int point = 0; point < conductivity.extent_int(1); ++point) {
        // const ScalarT& x = coords(cell,point,0);
        // const ScalarT& y = coords(cell,point,1);
        conductivity(cell,point,0,0) = sigma;
        conductivity(cell,point,0,1) = 0.;

        conductivity(cell,point,1,0) = 0.;
        conductivity(cell,point,1,1) = sigma;
      }
    }
  }
}

//**********************************************************************
}

#endif
