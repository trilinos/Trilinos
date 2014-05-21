/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_EXAMPLE_BELOS_SOLVE
#define KOKKOS_EXAMPLE_BELOS_SOLVE

// Tpetra
#include "Kokkos_ArithTraits.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "CGSolve.hpp" // for result_struct

#include "TrilinosCouplings_config.h"
#if defined(HAVE_TRILINOSCOUPLINGS_BELOS) && defined(HAVE_TRILINOSCOUPLINGS_MUELU)

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "BelosTpetraAdapter.hpp"

// MueLu
#include "MueLu_CreateTpetraPreconditioner.hpp"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {

template <class S, class LO, class GO, class N, class LMO>
result_struct
belos_solve(
  const Teuchos::RCP<Tpetra::CrsMatrix<S,LO,GO,N,LMO> >& A,
  const Teuchos::RCP<Tpetra::Vector<S,LO,GO,N> >& b,
  const Teuchos::RCP<Tpetra::Vector<S,LO,GO,N> >& x,
  const int use_muelu,
  const unsigned max_iter = 200,
  const typename Kokkos::Details::ArithTraits<S>::mag_type tolerance =
    Kokkos::Details::ArithTraits<S>::epsilon())
{
  typedef Tpetra::Operator<S,LO,GO,N> OperatorType;
  typedef Tpetra::MultiVector<S,LO,GO,N> VectorType;
  typedef typename VectorType::dot_type BelosScalarType;
  typedef Belos::LinearProblem<BelosScalarType, VectorType, OperatorType> ProblemType;
  typedef Belos::PseudoBlockCGSolMgr<BelosScalarType, VectorType, OperatorType> SolverType;
  typedef MueLu::TpetraOperator<S,LO,GO,N> PreconditionerType;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  //--------------------------------
  // Create preconditioner
  RCP<PreconditionerType> mueluPreconditioner;
  if (use_muelu) {
    std::string xmlFileName="muelu.xml";
    mueluPreconditioner = MueLu::CreateTpetraPreconditioner(A,xmlFileName);
  }

  //--------------------------------
  // Set up linear solver
  RCP<ParameterList> belosParams = Teuchos::parameterList();
  belosParams->set("Maximum Iterations", Teuchos::as<int>(max_iter));
  belosParams->set("Convergence Tolerance", tolerance);
  belosParams->set("Output Frequency", 1);
  belosParams->set("Output Style", Belos::Brief);
  RCP<SolverType> solver = rcp(new SolverType);
  RCP<ProblemType> problem = rcp(new ProblemType(A, x, b));
  if (use_muelu) problem->setRightPrec(mueluPreconditioner);
  const bool isSet = problem->setProblem();
  TEUCHOS_TEST_FOR_EXCEPTION(!isSet, std::runtime_error,
                             "Belos failed to set problem correctly.");
  solver->setParameters(belosParams);
  solver->setProblem(problem);

  //--------------------------------
  // Solve for nonlinear update
  Belos::ReturnType result = solver->solve();
  bool converged = (result == Belos::Converged);
  TEUCHOS_TEST_FOR_EXCEPTION(!converged, std::runtime_error,
                             "Belos solver did not converge!");

  result_struct cgsolve;
  cgsolve.iteration = solver->getNumIters();
  //cgsolve.norm_res = solver->achievedTol();

  Teuchos::Array< RCP<Teuchos::Time> > timers = solver->getTimers();
  for (int i=0; i<timers.size(); ++i) {
    if (timers[i]->name() == "Belos: PseudoBlockCGSolMgr total solve time")
      cgsolve.iter_time = timers[i]->totalElapsedTime();
    else if (timers[i]->name() == "Belos: Operation Op*x")
      cgsolve.matvectime = timers[i]->totalElapsedTime();
  }

  return cgsolve;
}


} // namespace Example
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#else

namespace Kokkos {
namespace Example {

template <class S, class LO, class GO, class N, class LMO>
result_struct
belos_solve(
  const Teuchos::RCP<Tpetra::CrsMatrix<S,LO,GO,N,LMO> >& A,
  const Teuchos::RCP<Tpetra::Vector<S,LO,GO,N> >& b,
  const Teuchos::RCP<Tpetra::Vector<S,LO,GO,N> >& x,
  const unsigned max_iter = 200,
  const typename Kokkos::Details::ArithTraits<S>::mag_type tolerance =
    Kokkos::Details::ArithTraits<S>::epsilon())
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
                             "Belos solver requested but not compiled!");
}

} // namespace Example
} // namespace Kokkos

#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXAMPLE_BELOS_SOLVE */
