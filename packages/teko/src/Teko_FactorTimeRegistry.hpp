// Teko_FactorTimeRegistry.hpp — process-wide accumulator of wall time spent
// factoring inverse operators.
//
// The free functions Teko::buildInverse / Teko::rebuildInverse (defined in
// Teko_InverseFactory.cpp) add their elapsed wall time here. The Krylov
// surrogate adaptive loop (Teko_KrylovSurrogate.hpp) drains the accumulator at
// hook entry to recover the setup (factorization) cost of the first solve —
// work the application performed before the solve, which is otherwise
// invisible to the hook. This keeps the application-facing interface
// unchanged: any application that builds its block preconditioner through
// Teko's buildInverse gets timed automatically; one that factors some other
// way simply reports zero setup time.
//
// The accessors are deliberately defined out-of-line in libteko.so (NOT
// inline here) so that every shared library agrees on a single accumulator —
// see the visibility discussion at the top of BelosAdaptiveHook.hpp.
#pragma once

namespace Teko {
namespace FactorTimeRegistry {

// Accumulate seconds of factorization wall time (thread-safe).
void add(double seconds);

// Current accumulated total in seconds.
double read();

// Zero the accumulator.
void reset();

}  // namespace FactorTimeRegistry
}  // namespace Teko
