
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

template<class MV, class OP>
Teuchos::RCP<Trilinos::Details::LinearSolver<MV, OP> >
LinearSolverFactory<MV, OP>::getLinearSolver (const std::string& solverName)
{
  typedef typename MV::scalar_type SC;
  typedef typename MV::local_ordinal_type LO;
  typedef typename MV::global_ordinal_type GO;
  typedef typename MV::node_type NT;

  static_assert(std::is_same<MV, Tpetra::MultiVector<SC, LO, GO, NT> >::value,
                "Ifpack2::Details::LinearSolverFactory::getLinearSolver: "
                "The MV template parameter must be a Tpetra::MultiVector "
                "specialization.  The most likely reason for seeing this error "
                "is that the two template parameters MV and OP got mixed up.");
  static_assert(std::is_same<OP, Tpetra::Operator<SC, LO, GO, NT> >::value,
                "Ifpack2::Details::LinearSolverFactory::getLinearSolver: "
                "The MV template parameter must be a Tpetra::MultiVector "
                "specialization.  The most likely reason for seeing this error "
                "is that the two template parameters MV and OP got mixed up.");

  return Teuchos::rcp (new Ifpack2::Details::LinearSolver<SC, LO, GO, NT> (solverName));
}

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_LINEARSOLVERFACTORY_DEF_HPP
