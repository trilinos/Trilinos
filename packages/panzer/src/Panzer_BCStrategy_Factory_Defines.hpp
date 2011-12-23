#include <iostream>
#include "Teuchos_RCP.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_BCStrategy_TemplateManager.hpp"

#undef PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER
#define PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(key, fClass, fType)				        \
									\
  struct fType ## _TemplateBuilder {					\
    const panzer::BC& m_bc;						\
    const Teuchos::RCP<panzer::GlobalData> m_global_data;               \
    fType ## _TemplateBuilder(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) : m_bc(bc), m_global_data(global_data) {} \
									\
    template<typename EvalT>						\
    Teuchos::RCP<panzer::BCStrategyBase> build() const {           	\
      fClass <EvalT>* ptr = new fClass <EvalT>(m_bc,m_global_data);   	\
      return Teuchos::rcp(ptr);						\
    }									\
    									\
  };

#undef PANZER_BUILD_BCSTRATEGY_OBJECTS
#define PANZER_BUILD_BCSTRATEGY_OBJECTS(key, fClass, fType)	        \
  if (bc.strategy() == key) {						\
    fType ## _TemplateBuilder builder(bc, global_data);			\
      bcs_tm->buildObjects(builder);				        \
      found = true;                                                     \
    }
