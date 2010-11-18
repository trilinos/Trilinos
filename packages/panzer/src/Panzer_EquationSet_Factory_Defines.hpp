#include <iostream>
#include "Teuchos_TestForException.hpp"
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_TemplateManager.hpp"

#undef PANZER_DECLARE_EQSET_TEMPLATE_BUILDER
#define PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(key, fClass, fType)	\
  									\
  struct fType ## _TemplateBuilder {					\
    const panzer::InputEquationSet& m_input_eq_set;			\
    const panzer::CellData& m_cell_data;                                \
    fType ## _TemplateBuilder(const panzer::InputEquationSet& ies, const panzer::CellData& cd) : m_input_eq_set(ies), m_cell_data(cd) {} \
									\
    template<typename EvalT>						\
    Teuchos::RCP<panzer::EquationSetBase> build() const {           	\
      fClass <EvalT>* ptr = new fClass <EvalT>(m_input_eq_set, m_cell_data); 	\
      return Teuchos::rcp(ptr);						\
    }									\
    									\
  };

#undef PANZER_BUILD_EQSET_OBJECTS
#define PANZER_BUILD_EQSET_OBJECTS(key, fClass, fType)                  \
    if (ies.name == key) {                                              \
      fType ## _TemplateBuilder builder(ies, cell_data);		\
      eq_set->buildObjects(builder);				        \
      found = true;                                                     \
    }
