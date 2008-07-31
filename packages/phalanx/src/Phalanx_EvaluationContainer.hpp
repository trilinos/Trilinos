#ifndef PHX_SCALAR_CONTAINER_HPP
#define PHX_SCALAR_CONTAINER_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_EvaluationContainer_Base.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_DataContainer_TemplateManager.hpp"

namespace PHX {

  /*! \brief Container that holds all data associated with a scalar type.


  */
  template <typename ScalarT, typename Traits>
  class EvaluationContainer : public PHX::EvaluationContainerBase<Traits> {
    
  public:
    
    EvaluationContainer();
    
    ~EvaluationContainer();
    
    //! Requests that the container must compute this field.
    void requireField(const PHX::FieldTag& f);

    void 
    registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& p);

    template <typename DataT> 
    Teuchos::ArrayRCP<DataT> getFieldData(const PHX::FieldTag& f);

    void postRegistrationSetup(std::size_t max_num_cells,
			       PHX::FieldManager<Traits>& fm);

    void evaluateFields(typename Traits::EvalData d);

    void preEvaluate(typename Traits::PreEvalData d);

    void postEvaluate(typename Traits::PostEvalData d);

    void print(std::ostream& os) const;

  protected:

    typedef PHX::DataContainer_TemplateManager<ScalarT, Traits> DCTM;

    PHX::DataContainer_TemplateManager<ScalarT, Traits> 
    data_container_template_manager_;
    
    typename Traits::Allocator allocator_;

  };
  
} 

#include "Phalanx_EvaluationContainer_Def.hpp"

#endif 
