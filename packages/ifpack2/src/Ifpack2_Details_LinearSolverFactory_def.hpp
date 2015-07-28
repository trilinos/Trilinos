
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

template<class SC, class LO, class GO, class NT>
RegisterLinearSolverFactory<SC, LO, GO, NT>::
RegisterLinearSolverFactory () {
  typedef Tpetra::MultiVector<SC, LO, GO, NT> MV;
  typedef Tpetra::Operator<SC, LO, GO, NT> OP;
  typedef typename MV::mag_type mag_type;
  typedef Trilinos::Details::LinearSolverFactory<MV, OP, mag_type> factory_base_type;
  typedef Ifpack2::Details::LinearSolverFactory<SC, LO, GO, NT> factory_impl_type;

#ifdef HAVE_TEUCHOSCORE_CXX11
  typedef std::shared_ptr<factory_base_type> base_ptr_type;
  typedef std::shared_ptr<factory_impl_type> impl_ptr_type;
#else
  typedef Teuchos::RCP<factory_base_type> base_ptr_type;
  typedef Teuchos::RCP<factory_impl_type> impl_ptr_type;
#endif // HAVE_TEUCHOSCORE_CXX11

  impl_ptr_type factory (new factory_impl_type ());
  base_ptr_type factoryBase = factory; // implicit cast to base class

  TEUCHOS_TEST_FOR_EXCEPTION
    (factoryBase.get () == NULL, std::logic_error, "Factory is null!  This "
     "should never happen!  Please report this bug to the Ifpack2 developers.");

#ifdef HAVE_TEUCHOS_DEBUG
  {
    using std::cerr;
    using std::endl;
    using Teuchos::TypeNameTraits;
    cerr << "Registering Ifpack2 LinearSolverFactory for"
         << " SC = " << TypeNameTraits<SC>::name ()
         << ", LO = " << TypeNameTraits<LO>::name ()
         << ", GO = " << TypeNameTraits<GO>::name ()
         << ", NT = " << TypeNameTraits<NT>::name ()
         << ", and mag_type = " << TypeNameTraits<mag_type>::name ()
         << endl;
  }
#endif // HAVE_TEUCHOS_DEBUG
  Trilinos::Details::registerLinearSolverFactory<MV, OP, mag_type> ("Ifpack2", factoryBase);
}

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_LINEARSOLVERFACTORY_DEF_HPP
