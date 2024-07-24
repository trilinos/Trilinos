// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_initializeKokkos.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Kokkos_Core.hpp"
#include "Tpetra_Details_checkLaunchBlocking.hpp"
#include "Tpetra_Details_KokkosTeuchosTimerInjection.hpp"
#include <cstdlib> // std::atexit
#include <string>
#include <vector>

namespace Tpetra {
namespace Details {

void finalizeKokkosIfNeeded() {
  if(!Kokkos::is_finalized()) {
    Kokkos::finalize();
  }
}
  
void
initializeKokkos ()
{
  if (! Kokkos::is_initialized ()) {
    std::vector<std::string> args = Teuchos::GlobalMPISession::getArgv ();
    int narg = static_cast<int> (args.size ()); // must be nonconst

    std::vector<char*> args_c;
    std::vector<std::unique_ptr<char[]>> args_;
    for (auto const& x : args) {
      args_.emplace_back(new char[x.size() + 1]);
      char* ptr = args_.back().get();
      strcpy(ptr, x.c_str());
      args_c.push_back(ptr);
    }
    args_c.push_back(nullptr);

    Kokkos::initialize (narg, narg == 0 ? nullptr : args_c.data ());
    checkOldCudaLaunchBlocking();

    std::atexit (finalizeKokkosIfNeeded);

  }
  // Add Kokkos calls to the TimeMonitor if the environment says so
  Tpetra::Details::AddKokkosDeepCopyToTimeMonitor();
  Tpetra::Details::AddKokkosFenceToTimeMonitor();
  Tpetra::Details::AddKokkosFunctionsToTimeMonitor();
}

} // namespace Details
} // namespace Tpetra

