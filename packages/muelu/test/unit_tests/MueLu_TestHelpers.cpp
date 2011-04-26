#include <Teuchos_UnitTestRepository.hpp>

#include "MueLu_TestHelpers.hpp"

// static members initialization of the class MueLu::TestHelpers::Parameters
Cthulhu::Parameters MueLu::TestHelpers::Parameters::cthulhuParameters = Cthulhu::Parameters(Teuchos::UnitTestRepository::getCLP());
