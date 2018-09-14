#ifndef MINIEM_RANDOM_FORCING_DECL_HPP
#define MINIEM_RANDOM_FORCING_DECL_HPP

#include "PanzerAdaptersSTK_config.hpp"

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_FieldLibrary.hpp"

#include <string>

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace mini_em {

  using panzer::Cell;
  using panzer::Point;
  using panzer::Dim;

/** Random current source
  */
template<typename EvalT, typename Traits>
class RandomForcing : public panzer::EvaluatorWithBaseImpl<Traits>,
                      public PHX::EvaluatorDerived<EvalT, Traits>  {

public:
    RandomForcing(const std::string & name,
                  const panzer::IntegrationRule & ir,
                  const panzer::FieldLayoutLibrary & fl,
                  const unsigned int & seed);

    void evaluateFields(typename Traits::EvalData d);


private:
  typedef typename EvalT::ScalarT ScalarT;

  // Simulation source
  PHX::MDField<ScalarT,Cell,Point,Dim> current;
  PHX::MDField<const ScalarT,Cell,Point,Dim> coords;
  int ir_degree, ir_index, ir_dim;
};

}

#include "MiniEM_RandomForcing_impl.hpp"

#endif
