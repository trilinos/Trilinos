#ifndef PIKE_STATUS_TEST_ABSTRACT_FACTORY_HPP
#define PIKE_STATUS_TEST_ABSTRACT_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include <string>

namespace Teuchos { class ParameterList; }

namespace pike {

  class StatusTest;

  class StatusTestAbstractFactory {
    
  public:
    
    virtual ~StatusTestAbstractFactory() {}

    //! Returns true if the status test type can be built by this factory.
    virtual bool supportsType(const std::string& type) const = 0;

    //! Builds the status test object from the parameter list.
    virtual
    Teuchos::RCP<pike::StatusTest> 
    buildStatusTests(const Teuchos::RCP<Teuchos::ParameterList>& p) const = 0;

  };

}

#endif
