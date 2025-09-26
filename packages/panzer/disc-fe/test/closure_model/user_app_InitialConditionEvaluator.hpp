#ifndef USER_APP_INITIAL_CONDITION_EVALUATOR_HPP
#define USER_APP_INITIAL_CONDITION_EVALUATOR_HPP

#include "Phalanx_Evaluator_WithBaseImpl.hpp"

namespace user_app {

  template<typename EvalT, typename Traits>
  class InitialConditionEvaluator
    : public panzer::EvaluatorWithBaseImpl<Traits>,
      public PHX::EvaluatorDerived<EvalT, Traits>
  {
  public:
    using ScalarT = typename EvalT::ScalarT;

    InitialConditionEvaluator(const std::string& field_name,
                              const Teuchos::RCP<PHX::DataLayout>& layout)
    {
      temperature_ = PHX::MDField<ScalarT,panzer::Cell,panzer::Point>(field_name,layout);
      this->addEvaluatedField(temperature_);
      this->setName(std::string("InitialConditionEvaluator: ")+field_name);
    }

    void evaluateFields(typename Traits::EvalData workset)
    {
      const auto& coords = workset.getCellNodes();
      const int num_points = coords.extent(1);
      auto temperature = temperature_;
      auto policy = panzer::HP::inst().teamPolicy<ScalarT,PHX::exec_space>(workset.num_cells);
      Kokkos::parallel_for("Initial Condition Evaluator",policy,KOKKOS_LAMBDA(Kokkos::TeamPolicy<PHX::exec_space>::member_type team)
      {
        const int cell = team.league_rank();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_points), [&] (const int point)
        {
          const double& x = coords(cell,point,0);
          if (x > -0.05)
            temperature(cell,point) = 2.0;
          else
            temperature(cell,point) = 1.0;
        });
      });
    }
    
  private:
    PHX::MDField<ScalarT,panzer::Cell,panzer::Point> temperature_;
  };
}

#endif
