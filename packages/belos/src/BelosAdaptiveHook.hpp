// BelosAdaptiveHook.hpp — header-only post-solve callback registry.
//
// After a flexible GMRES solve converges, BelosBlockGmresSolMgr invokes
// invoke() here.  The Teko package registers an actual implementation via
// registerHook() from a static initialiser in Teko_KrylovSurrogateInit.cpp.
// Everything is specialised to double so that Belos stays free of Teko
// headers and there is no circular dependency.
//
// registerHook()/invoke() are declared here but DEFINED (non-inline) in
// BelosAdaptiveHook.cpp, which is compiled once into libbelos.so.  This is
// deliberate: BlockGmresSolMgr<double, Thyra::MultiVectorBase<double>,
// Thyra::LinearOpBase<double>> is instantiated in more than one shared
// library (libteko.so and the pybind11 extension module).  A header-only
// inline function with a function-local static (the usual Meyers-singleton
// pattern) would give each of those libraries its own private copy of the
// registry whenever they are built with different visibility settings (e.g.
// -fvisibility=hidden), so registerHook() in one library would silently not
// be seen by invoke() in another. Routing both through a single
// out-of-line definition in libbelos.so avoids that entirely.
#pragma once

#include <functional>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
// BelosGmresIteration.hpp's GmresIterationState stores
// Teuchos::RCP<const Teuchos::SerialDenseMatrix<...>> members but does not
// include Teuchos_SerialDenseMatrix.hpp itself (it relies on whatever
// translation unit includes it having already pulled that in transitively).
// Include it explicitly here since this header is the entry point.
#include "Teuchos_SerialDenseMatrix.hpp"
#include "BelosGmresIteration.hpp"
#include "BelosLinearProblem.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"

namespace Belos {
namespace AdaptiveHook {

using SC      = double;
using MV      = Thyra::MultiVectorBase<SC>;
using OP      = Thyra::LinearOpBase<SC>;
using State   = GmresIterationState<SC, MV>;
using Problem = LinearProblem<SC, MV, OP>;

// Diagnostics about the solve that just converged, captured by
// BelosBlockGmresSolMgr::solve() at the point it invokes the hook. These let
// the hook (Teko::KrylovSurrogate::adaptiveLoop) record convergence
// diagnostics for the first solve alongside those of any follow-up solve it
// runs itself.
//   wall_time_sec — time spent inside solve()'s main iteration loop (timer
//                   started at the top of solve()), excluding any
//                   preconditioner construction the caller did beforehand.
//   achieved_tol  — the same relative-residual quantity solve() will store in
//                   achievedTol_, computed a few lines earlier so it is
//                   available at hook time.
struct SolveMetrics {
    double wall_time_sec;
    double achieved_tol;
};

// Function type called after convergence of a flexible GMRES solve.
//   state       — final Krylov state (V, Z, curDim, H, R)
//   problem     — the linear problem (operator, prec, LHS, RHS)
//   blocked_op  — problem operator cast to BlockedLinearOpBase (may be null
//                 if the operator is not blocked; hook must handle this)
//   params      — the resolved parameter list of the solver that just
//                 converged (Maximum Iterations, Num Blocks, Convergence
//                 Tolerance, Flexible Gmres, Verbosity, etc.), so the hook
//                 can reuse the same settings for any follow-up solve.
//   metrics     — wall time and achieved tolerance for this solve.
using HookFn = std::function<void(
    const State&                                              state,
    Teuchos::RCP<Problem>                                     problem,
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<SC>>        blocked_op,
    Teuchos::RCP<const Teuchos::ParameterList>                params,
    const SolveMetrics&                                       metrics)>;

// Called once from Teko_KrylovSurrogateInit.cpp static initialiser.
void registerHook(HookFn fn);

// Called from BelosBlockGmresSolMgr::solve() after a flexible solve converges.
void invoke(
    const State&                                              state,
    Teuchos::RCP<Problem>                                     problem,
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<SC>>        blocked_op,
    Teuchos::RCP<const Teuchos::ParameterList>                params,
    const SolveMetrics&                                       metrics);

} // namespace AdaptiveHook
} // namespace Belos
