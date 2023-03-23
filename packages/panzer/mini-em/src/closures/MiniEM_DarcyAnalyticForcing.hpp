#ifndef MINIEM_DARCYANALYTICFORCING_HPP
#define MINIEM_DARCYANALYTICFORCING_HPP

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

/** Darcy source with know analytic solution
  */
template<typename EvalT, typename Traits>
class DarcyAnalyticForcing : public panzer::EvaluatorWithBaseImpl<Traits>,
                      public PHX::EvaluatorDerived<EvalT, Traits>  {

public:
    DarcyAnalyticForcing(const std::string & name,
                  const panzer::IntegrationRule & ir,
                  const panzer::FieldLayoutLibrary & fl,
                  const double kappa,
                  const std::string& basisName="u");

    void postRegistrationSetup(typename Traits::SetupData d,
                               PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);


private:
  typedef typename EvalT::ScalarT ScalarT;

  // Simulation source
  PHX::MDField<ScalarT,Cell,Point> source;
  int ir_degree, ir_index, ir_dim;
  double kappa_;

  using device_type = PHX::Device;
};

}

#include "MiniEM_DarcyAnalyticForcing_impl.hpp"

#endif
