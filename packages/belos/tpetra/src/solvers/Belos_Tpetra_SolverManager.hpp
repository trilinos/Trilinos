//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
//@HEADER

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

