#ifndef MUELU_TEST_HELPERS_SMOOTHERS_H
#define MUELU_TEST_HELPERS_SMOOTHERS_H

#include <Teuchos_FancyOStream.hpp>
#include <Xpetra_Operator.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_UseDefaultTypes.hpp"

// Helper functions to test if derived classes conforms to the SmootherBase and SmootherPrototype interfaces

namespace MueLuTests {

  namespace Smoother {
#include "MueLu_UseShortNames.hpp"

    //! Test if Apply() throw an exception when called before Setup()
    void testApplyNoSetup(SmootherPrototype & smoother, Teuchos::FancyOStream & out, bool & success);

    //! Apply smoother with Poisson1D(125), RHS=0 and X initialized to 1
    //  This method calls Setup() internally
    ST::magnitudeType testApply_A125_X1_RHS0(SmootherPrototype & smoother, Teuchos::FancyOStream & out, bool & success);
    
    //!
    //! Apply smoother with Poisson1D(125) RHS choosed randomly and X initialized to 0
    ST::magnitudeType testApply_A125_X0_RandomRHS(SmootherPrototype & smoother, Teuchos::FancyOStream & out, bool & success);

    //! Test if a smoother reduces effectively the error on a simple problem
    // TODO

    //! Test direct solver
    void testDirectSolver(SmootherPrototype & smoother, Teuchos::FancyOStream & out, bool & success);

  } // namespace Smoother
  
} // MueLuTests

#endif // MUELU_TEST_HELPERS_SMOOTHERS_H
