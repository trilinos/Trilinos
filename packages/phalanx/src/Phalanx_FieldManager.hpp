#ifndef PHX_FIELD_MANAGER_HPP
#define PHX_FIELD_MANAGER_HPP

#include <cstddef>
#include <iostream>
#include <vector>
#include <algorithm>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_Field.hpp"
#include "Phalanx_EvaluationContainer_TemplateManager.hpp"

namespace PHX {

  template<typename Traits>
  class FieldManager {
    
  public:

    typedef typename PHX::EvaluationContainer_TemplateManager<Traits>::iterator iterator;

    FieldManager();

    ~FieldManager();
    
    void requireFieldForAllTypes(const PHX::FieldTag& v);
    
    template<typename EvalT>
    void requireFieldForScalarType(const PHX::FieldTag& v);

    void registerEvaluatorForAllTypes(const Teuchos::RCP< PHX::Evaluator<Traits> >& p);
    
    template<typename EvalT>
    void registerEvaluatorForScalarType(const Teuchos::RCP< PHX::Evaluator<Traits> >& p);

    void registerEvaluatorForScalarType(typename PHX::FieldManager<Traits>::iterator it,
				 const Teuchos::RCP< PHX::Evaluator<Traits> >& p);
    
    template<typename DataT, typename EvalT> 
    void getFieldData(PHX::Field<DataT>& h);
    
    template<typename DataT, typename EvalT> 
    void getFieldData(const PHX::FieldTag& v, Teuchos::ArrayRCP<DataT>& d);
    
    void postRegistrationSetup(std::size_t max_num_cells);

    template<typename EvalT>
    void evaluateFields(typename Traits::EvalData d);

    template<typename EvalT>
    void preEvaluate(typename Traits::PreEvalData d);

    template<typename EvalT>
    void postEvaluate(typename Traits::PostEvalData d);

    std::size_t getMaxNumCells() const;

    //! Return iterator to first EvaluationContainer
    typename FieldManager::iterator begin();

    //! Return iterator to last EvaluationContainer
    typename FieldManager::iterator end();

    void print(std::ostream& os) const;

  private:

    typedef PHX::EvaluationContainer_TemplateManager<Traits> SCTM;

    std::size_t m_num_evaluation_types;

    PHX::EvaluationContainer_TemplateManager<Traits> m_eval_containers;

    std::size_t m_max_num_cells;

  };

  template<typename Traits>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::FieldManager<Traits>& vm);

} 

#include "Phalanx_FieldManager_Def.hpp"

#endif 
