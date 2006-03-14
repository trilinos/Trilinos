
#ifndef TEST_SINGLE_BELOS_THYRA_SOLVER_HPP
#define TEST_SINGLE_BELOS_THYRA_SOLVER_HPP

#include "BelosConfigDefs.hpp"
#include "Teuchos_FancyOStream.hpp"

namespace Teuchos { class ParameterList; }

namespace Thyra {

/** \brief Testing function for a single belos solver with a single matrix.
 *
 */
bool test_single_belos_thyra_solver(
  const std::string                       matrixFile
  ,const bool                             testTranspose
  ,const int                              numRandomVectors
  ,const double                           maxFwdError
  ,const int                              maxIterations
  ,const double                           maxResid
  ,const double                           maxSolutionError
  ,const bool                             showAllTests
  ,const bool                             dumpAll
  ,Teuchos::ParameterList                 *fwdSolveParamList
  ,Teuchos::ParameterList                 *adjSolveParamList
  ,Teuchos::FancyOStream                  *out
  );

} // namespace Thyra

#endif // TEST_SINGLE_BELOS_THYRA_SOLVER_HPP
