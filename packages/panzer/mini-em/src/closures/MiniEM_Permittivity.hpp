#ifndef MINIEM_PERMITTIVITY_DECL_HPP
#define MINIEM_PERMITTIVITY_DECL_HPP

#include "PanzerAdaptersSTK_config.hpp"

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_FieldLibrary.hpp"

#include <string>

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace mini_em {
    
  using panzer::Cell;
  using panzer::Point;

  template<typename EvalT, typename Traits>
  class Permittivity : public panzer::EvaluatorWithBaseImpl<Traits>,
                       public PHX::EvaluatorDerived<EvalT, Traits>  {

  public:
    Permittivity(const std::string & name,
                 const panzer::IntegrationRule & ir,
                 const panzer::FieldLayoutLibrary & fl,
                 const double & epsilon_,
                 const std::string& DoF_);
                                                                        
    void evaluateFields(typename Traits::EvalData d);


  private:
    typedef typename EvalT::ScalarT ScalarT;

    PHX::MDField<ScalarT,Cell,Point> permittivity;
    PHX::MDField<const ScalarT,Cell,Point,Dim> coords;
    int ir_degree;
    double epsilon;
  };

}

#include "MiniEM_Permittivity_impl.hpp"

#endif
