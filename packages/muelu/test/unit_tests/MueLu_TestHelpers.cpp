#include <Teuchos_UnitTestRepository.hpp>

#include "MueLu_TestHelpers.hpp"

// static members initialization of the class MueLu::TestHelpers::Parameters
Xpetra::Parameters MueLu::TestHelpers::Parameters::xpetraParameters = Xpetra::Parameters(Teuchos::UnitTestRepository::getCLP());
