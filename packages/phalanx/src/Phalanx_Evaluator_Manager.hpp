// @HEADER
// @HEADER

#ifndef PHX_FIELD_EVALUATOR_MANAGER_HPP
#define PHX_FIELD_EVALUATOR_MANAGER_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_FieldTag_Comparison.hpp"
#include "Phalanx_Evaluator.hpp"
#include "Phalanx_DebugStrings.hpp"

namespace PHX {
  
  template<typename Traits> class FieldManager;

  /*! @brief Class to sort which Evaluators should be called and the order in which to call them such that all dependencies are met.
   */
  template<typename Traits>
  class EvaluatorManager {

  public:

    EvaluatorManager(const std::string& evaluator_type_name = "???");
    
    ~EvaluatorManager();
    
    //! Require a variable to be evaluated.
    void requireField(const PHX::FieldTag& v);
    
    //! Registers a variable provider with the manager.
    void 
    registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& p);
    
    /*! Sets up all field dependencies.  This should only be called
      once all variables and DOFs have been added and all providers
      have been registered.  Sorts variable and creates dependency
      lists and evaluation order
    */
    void sortAndOrderEvaluators();
    
    /*! Calls post registration setup on all variable providers.
    */
    void postRegistrationSetup(PHX::FieldManager<Traits>& vm);
    
    //! Compute the required variables for the fill on the specific element.
    void evaluateFields(typename Traits::EvalData d);
    
    /*! \brief This routine is called before each residual/Jacobian fill.
      
        This routine is called ONCE on the provider before the fill
        loop over elements is started.  This allows us to reset global
        objects between each fill.  An example is to reset a provider
        that monitors the maximum grid peclet number in a cell.  This
        call would zero out the maximum for a new fill.
    */
    void preEvaluate(typename Traits::PreEvalData d);
    
    /*! \brief This routine is called after each residual/Jacobian fill.
      
        This routine is called ONCE on the provider after the fill
        loop over elements is completed.  This allows us to evaluate
        any post fill data.  An example is to print out some
        statistics such as the maximum grid peclet number in a cell.
    */
    void postEvaluate(typename Traits::PostEvalData d);
    
    void setEvaluationTypeName(const std::string& evaluation_type_name);
    
    const std::vector< Teuchos::RCP<PHX::FieldTag> >& getFieldTags();

    bool sortingCalled() const;

    //! Printing
    void print(std::ostream& os) const;
    
  protected:
    
    /*! @brief Create and arrange the dependency list in the correct
     *  order it should be evaluated.
     */
    void createProviderEvaluationOrder();
    
  protected:
    
    //! Fields required by the user.
    std::vector< Teuchos::RCP<PHX::FieldTag> > fields_;
    
    //@{
    /*!
      @name Evaluator Objects
      @brief Stores information about variable provider objects.
    */    
    std::vector< Teuchos::RCP<PHX::Evaluator<Traits> > > 
    varProviders;
    
    std::vector< std::vector< Teuchos::RCP<PHX::FieldTag> > > 
    providerVariables;

    std::vector< std::vector< Teuchos::RCP<PHX::FieldTag> > > 
    providerRequirements;

    std::vector<std::string> providerNames;
    //@}

    
    //@{
    /*! @name Evaluation Order Objects
      
        Stores information about the order that providers need to be
        called to evaluate fields correctly.
    */
    std::vector<int> providerEvalOrderIndex;
    //@}
    
    std::string evaluation_type_name_;

    //! Flag to tell the setup has been called.
    bool sorting_called_;
    
  };
  
  template<typename Traits>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::EvaluatorManager<Traits>& m);

}

#include "Phalanx_Evaluator_Manager_Def.hpp"

#endif
