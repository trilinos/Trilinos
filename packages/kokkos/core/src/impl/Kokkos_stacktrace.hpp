#ifndef KOKKOS_STACKTRACE_HPP
#define KOKKOS_STACKTRACE_HPP

#include <ostream>
#include <string>

namespace Kokkos {
namespace Impl {

/// \brief Return the demangled version of the input symbol, or the
///   original input if demangling is not possible.
std::string demangle (const std::string& name);

/// \brief Save the current stacktrace.
///
/// You may only save one stacktrace at a time.  If you call this
/// twice, the second call will overwrite the result of the first
/// call.
void save_stacktrace ();

/// \brief Print the currently saved stacktrace, if any, to the given
///   output stream.
void print_saved_stacktrace (std::ostream& out);

void print_demangled_saved_stacktrace (std::ostream& out);
  
} // namespace Impl
} // namespace Kokkos

#endif // KOKKOS_STACKTRACE_HPP
