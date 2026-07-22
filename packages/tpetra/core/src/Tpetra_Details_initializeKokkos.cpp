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
#include "Teuchos_oblackholestream.hpp"
#include "Kokkos_Core.hpp"
#include "Tpetra_Details_checkLaunchBlocking.hpp"
#include "Tpetra_Details_KokkosTeuchosTimerInjection.hpp"
#include <cstdlib>  // std::atexit
#include <string>
#include <vector>
#include "KokkosKernels_EagerInitialize.hpp"

namespace Tpetra {
namespace Details {

class HideOutputExceptOnProcess0 {
 public:
  HideOutputExceptOnProcess0(std::ostream& stream,
                             const int myRank)
    : stream_(stream)
    , originalBuffer_(stream.rdbuf()) {
    if (myRank != 0) {
      stream.rdbuf(blackHole_.rdbuf());
    }
  }

  ~HideOutputExceptOnProcess0() {
    stream_.rdbuf(originalBuffer_);
  }

 private:
  std::ostream& stream_;
  decltype(std::cout.rdbuf()) originalBuffer_;
  Teuchos::oblackholestream blackHole_;
};

void finalizeKokkosIfNeeded() {
  if (!Kokkos::is_finalized()) {
    Kokkos::finalize();
  }
}

void initializeKokkos(int* argc, char*** argv, int myRank) {
  if (!Kokkos::is_initialized()) {
    HideOutputExceptOnProcess0 hideCerr(std::cerr, myRank);
    HideOutputExceptOnProcess0 hideCout(std::cout, myRank);

    if (!argc) {
      // If there are no args, we try to find them via Teuchos
      std::vector<std::string> args = Teuchos::GlobalMPISession::getArgv();
      int narg                      = static_cast<int>(args.size());  // must be nonconst

      std::vector<char*> args_c;
      std::vector<std::unique_ptr<char[]>> args_;
      for (auto const& x : args) {
        args_.emplace_back(new char[x.size() + 1]);
        char* ptr = args_.back().get();
        strcpy(ptr, x.c_str());
        args_c.push_back(ptr);
      }
      args_c.push_back(nullptr);

      Kokkos::initialize(narg, narg == 0 ? nullptr : args_c.data());
    } else {
      // If there are args, we use them
      Kokkos::initialize(*argc, *argv);
    }

    checkOldCudaLaunchBlocking();

    std::atexit(finalizeKokkosIfNeeded);
  }
  // Add Kokkos calls to the TimeMonitor if the environment says so
  Tpetra::Details::AddKokkosDeepCopyToTimeMonitor();
  Tpetra::Details::AddKokkosFenceToTimeMonitor();
  Tpetra::Details::AddKokkosFunctionsToTimeMonitor();

  // Now that the Kokkos backend(s) are initialized,
  // initialize all KokkosKernels TPLs.
  KokkosKernels::eager_initialize();
}

}  // namespace Details
}  // namespace Tpetra
