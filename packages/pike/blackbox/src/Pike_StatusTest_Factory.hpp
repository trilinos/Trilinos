#ifndef PIKE_STATUS_TEST_FACTORY
#define PIKE_STATUS_TEST_FACTORY

#include "Teuchos_RCP.hpp"

namespace Teuchos { class ParameterList; }

namespace pike {

  class StatusTest;

  Teuchos::RCP<pike::StatusTest> buildStatusTests(const Teuchos::ParameterList& p);

}

#endif
