// BelosAdaptiveHook.cpp — out-of-line storage for the adaptive-hook registry
// declared in BelosAdaptiveHook.hpp.  See that header for why this lives in
// exactly one translation unit (compiled into libbelos.so) instead of being
// a header-only inline function.

#include "BelosAdaptiveHook.hpp"

namespace Belos {
namespace AdaptiveHook {

namespace {
HookFn& globalHook() {
    static HookFn hook;
    return hook;
}
} // namespace

void registerHook(HookFn fn) { globalHook() = std::move(fn); }

void invoke(
    const State&                                              state,
    Teuchos::RCP<Problem>                                     problem,
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<SC>>        blocked_op)
{
    auto& h = globalHook();
    if (h) h(state, problem, blocked_op);
}

} // namespace AdaptiveHook
} // namespace Belos
