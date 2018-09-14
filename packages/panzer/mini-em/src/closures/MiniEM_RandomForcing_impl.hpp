#ifndef MINIEM_RANDOM_FORCING_IMPL_HPP
#define MINIEM_RANDOM_FORCING_IMPL_HPP

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"

namespace mini_em {

//**********************************************************************
template <typename EvalT,typename Traits>
RandomForcing<EvalT,Traits>::RandomForcing(const std::string & name,
                                           const panzer::IntegrationRule & ir,
                                           const panzer::FieldLayoutLibrary & fl,
                                           const unsigned int & seed)
{
  using Teuchos::RCP;

  Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_vector;
  ir_degree = ir.cubature_degree;
  ir_dim = ir.spatial_dimension;

  current = PHX::MDField<ScalarT,Cell,Point,Dim>(name, data_layout);
  this->addEvaluatedField(current);

  Teuchos::RCP<const panzer::PureBasis> basis = fl.lookupBasis("E_edge");
  const std::string coordName = panzer::GatherBasisCoordinates<EvalT,Traits>::fieldName(basis->name());
  coords = PHX::MDField<const ScalarT,Cell,Point,Dim>(coordName, basis->coordinates);
  this->addDependentField(coords);

  std::string n = "Random Forcing";
  this->setName(n);
  std::srand(seed);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void RandomForcing<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  using panzer::index_t;

  // double time = workset.time;


  if (ir_dim == 3) {
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int point = 0; point < current.extent_int(1); ++point) {
        // const ScalarT& x = coords(cell,point,0);
        // const ScalarT& y = coords(cell,point,1);
        // const ScalarT& z = coords(cell,point,2);
        current(cell,point,0) = double(std::rand())/double(RAND_MAX);
        current(cell,point,1) = double(std::rand())/double(RAND_MAX);
        current(cell,point,2) = double(std::rand())/double(RAND_MAX);
      }
    }
  } else {
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int point = 0; point < current.extent_int(1); ++point) {
        // const ScalarT& x = coords(cell,point,0);
        // const ScalarT& y = coords(cell,point,1);
        current(cell,point,0) = double(std::rand())/double(RAND_MAX);
        current(cell,point,1) = double(std::rand())/double(RAND_MAX);
      }
    }
  }
}

//**********************************************************************
}

#endif
