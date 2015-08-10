
/// \file   Ifpack2_Details_LinearSolverFactory_def.hpp
/// \author Mark Hoemmen
/// \brief  Definition of Ifpack2::Details::LinearSolverFactory.

#ifndef IFPACK2_DETAILS_LINEARSOLVERFACTORY_DEF_HPP
#define IFPACK2_DETAILS_LINEARSOLVERFACTORY_DEF_HPP

#include <Ifpack2_Details_LinearSolver.hpp>
#include <Trilinos_Details_LinearSolverFactory.hpp>
#include <type_traits> // std::is_same

namespace Ifpack2 {
namespace Details {

template<class SC, class LO, class GO, class NT>
Teuchos::RCP<typename LinearSolverFactory<SC, LO, GO, NT>::solver_type>
LinearSolverFactory<SC, LO, GO, NT>::
getLinearSolver (const std::string& solverName)
{
  return Teuchos::rcp (new Ifpack2::Details::LinearSolver<SC, LO, GO, NT> (solverName));
}

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_LINEARSOLVERFACTORY_DEF_HPP
