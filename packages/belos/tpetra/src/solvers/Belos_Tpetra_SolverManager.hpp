// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_TPETRA_SOLVER_MANAGER_HPP
#define BELOS_TPETRA_SOLVER_MANAGER_HPP

#include "Belos_Tpetra_Krylov.hpp"
#include "Belos_Tpetra_SolverManagerBase.hpp"

namespace BelosTpetra {
namespace Impl {

/// \brief Concrete Belos::SolverManager subclass that users get when
///   they ask for one of the methods that a subclass of Krylov
///   implements.
template<class SC, class MV, class OP,
	 template<class, class, class> class KrylovSubclassType>
class SolverManager :
    public SolverManagerBase<SC, MV, OP>
{
private:
  using solver_impl_type = KrylovSubclassType<SC, MV, OP>;  
  using solver_base_type = Krylov<SC, MV, OP>;
  
  static Teuchos::RCP<solver_base_type>
  makeSolverImplementation ()
  {
    return Teuchos::RCP<solver_base_type> (static_cast<solver_base_type*> (new solver_impl_type));
  }
  
public:
  SolverManager (const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) :
    SolverManagerBase<SC, MV, OP>::SolverManagerBase (makeSolverImplementation (),
						      params)
  {}

  virtual ~SolverManager () = default;

  virtual Teuchos::RCP<Belos::SolverManager<SC, MV, OP> >
  clone () const override
  {
    using this_type = SolverManager<SC, MV, OP, KrylovSubclassType>;
    Teuchos::RCP<this_type> solver (new this_type);
    return Teuchos::rcp_implicit_cast<Belos::SolverManager<SC, MV, OP>> (solver);
  }
};

} // namespace Impl
} // namespace BelosTpetra

#endif // BELOS_TPETRA_SOLVER_MANAGER_HPP

