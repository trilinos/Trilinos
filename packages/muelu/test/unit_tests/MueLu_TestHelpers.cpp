#include <Teuchos_UnitTestRepository.hpp>

#include "MueLu_TestHelpers.hpp"

namespace MueLuTests {

  // static members initialization of the class TestHelpers::Parameters
  Xpetra::Parameters TestHelpers::Parameters::xpetraParameters = Xpetra::Parameters(Teuchos::UnitTestRepository::getCLP());

}
