#ifndef PANZER_ASSEMBLY_ENGINE_HPP
#define PANZER_ASSEMBLY_ENGINE_HPP

#include "Panzer_Base.hpp"
#include "Teuchos_RCP.hpp"
#include "Panzer_BC.hpp"

namespace panzer {
  template <typename LO, typename GO> class FieldManagerBuilder;
  struct AssemblyEngineInArgs;
}

namespace panzer {

  //! Class for the matrix and residual fill.
  template <typename EvalT, typename LO, typename GO>
    class AssemblyEngine : public panzer::Base {

  public:    
    
    AssemblyEngine(const Teuchos::RCP<panzer::FieldManagerBuilder<LO,GO> >& 
		   fmb);
    
    void evaluate(const panzer::AssemblyEngineInArgs& input_arguments);

    Teuchos::RCP<panzer::FieldManagerBuilder<LO,GO> > getManagerBuilder()
      { return m_field_manager_builder; }
    
  protected:
      
    void evaluateBCs(const panzer::BCType bc_type, 
		     const panzer::AssemblyEngineInArgs& input_arguments);

  protected:
    
      Teuchos::RCP<panzer::FieldManagerBuilder<LO,GO> > 
      m_field_manager_builder;
    
  };
  
}

#include "Panzer_AssemblyEngineT.hpp"

#endif
