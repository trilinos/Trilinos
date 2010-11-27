#include <iostream>
#include "Teuchos_RCP.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_BCStrategy_TemplateManager.hpp"

#undef PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER
#define PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(key, fClass, fType)				        \
									\
  struct fType ## _TemplateBuilder {					\
    const panzer::BC& m_bc;						\
    fType ## _TemplateBuilder(const panzer::BC& bc) : m_bc(bc) {}	\
									\
    template<typename EvalT>						\
    Teuchos::RCP<panzer::BCStrategyBase> build() const {           	\
      fClass <EvalT>* ptr = new fClass <EvalT>(m_bc);             	\
      return Teuchos::rcp(ptr);						\
    }									\
    									\
  };

#undef PANZER_BUILD_BCSTRATEGY_OBJECTS
#define PANZER_BUILD_BCSTRATEGY_OBJECTS(key, fClass, fType)	        \
  if (bc.strategy() == key) {						\
      fType ## _TemplateBuilder builder(bc);				\
      bcs_tm->buildObjects(builder);				        \
      found = true;                                                     \
    }
