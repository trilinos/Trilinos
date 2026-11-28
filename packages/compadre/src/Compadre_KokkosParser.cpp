// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#include "Compadre_KokkosParser.hpp"
#include <iostream>
#include <mutex>
#include <stdexcept>

using namespace Compadre;

namespace {

// the Kokkos runtime started by a KokkosParser, shared by every KokkosParser and by every holder of
// KokkosParser::runtime(); its ScopeGuard finalizes Kokkos when the last of them releases it
std::weak_ptr<Kokkos::ScopeGuard> started_runtime;
std::mutex started_runtime_mutex;

void throwIfFinalized() {
    if (Kokkos::is_finalized()) {
        throw std::runtime_error("KokkosParser: Kokkos has already been finalized, and it cannot be initialized "
                "again in the same process.");
    }
}

// for the constructors given arguments or settings: initializes Kokkos with start(), and throws if Kokkos is
// already initialized, since it would then be initialized twice with arguments
template <typename Start>
std::shared_ptr<void> startRuntime(Start&& start) {
    std::lock_guard<std::mutex> lock(started_runtime_mutex);
    throwIfFinalized();
    if (Kokkos::is_initialized()) {
        throw std::runtime_error("KokkosParser: Kokkos is already initialized, so it cannot be initialized again "
                "with these arguments. Code that runs after Kokkos is initialized, or without access to the "
                "command line, should use KokkosParser() to join the running Kokkos.");
    }
    std::shared_ptr<Kokkos::ScopeGuard> runtime = start();
    started_runtime = runtime;
    return runtime;
}

// for KokkosParser(): joins the running Kokkos (sharing its runtime if a KokkosParser started it); if Kokkos is
// not initialized, warns that command-line arguments are ignored and initializes it from the environment only
std::shared_ptr<void> joinOrStartFromEnvironment() {
    std::lock_guard<std::mutex> lock(started_runtime_mutex);
    if (auto runtime = started_runtime.lock()) return runtime;
    if (Kokkos::is_initialized()) return nullptr; // initialized outside Compadre, which also finalizes it
    throwIfFinalized();
    std::cerr << "Compadre::KokkosParser: WARNING: initializing Kokkos without the command line, so Kokkos "
                 "command-line arguments are ignored; only the KOKKOS_* environment variables are used.\n";
    std::shared_ptr<Kokkos::ScopeGuard> runtime = std::make_shared<Kokkos::ScopeGuard>(KokkosInitArguments());
    started_runtime = runtime;
    return runtime;
}

} // anonymous namespace

std::shared_ptr<void> KokkosParser::runtime() {
    std::lock_guard<std::mutex> lock(started_runtime_mutex);
    return started_runtime.lock();
}

KokkosParser::KokkosParser(KokkosInitArguments args) {
    _runtime = startRuntime([&]() { return std::make_shared<Kokkos::ScopeGuard>(args); });
}

// copies the arguments, since Kokkos removes the arguments it recognizes from the argv it is given
KokkosParser::KokkosParser(int narg, char* args[])
        : KokkosParser(std::vector<std::string>(args, args + narg)) {}

KokkosParser::KokkosParser(std::vector<std::string> stdvec_args) {
    // Kokkos takes argv as char* entries followed by a nullptr, and removes the entries it recognizes
    std::vector<char*> char_args;
    for (auto& arg : stdvec_args) {
        char_args.push_back(arg.data());
    }
    char_args.push_back(nullptr);
    int narg = (int)stdvec_args.size();

    _runtime = startRuntime([&]() { return std::make_shared<Kokkos::ScopeGuard>(narg, char_args.data()); });
}

KokkosParser::KokkosParser() {
    _runtime = joinOrStartFromEnvironment();
}

std::string KokkosParser::status() {
    std::stringstream stream;
    Kokkos::print_configuration(stream, true);
    std::string status = stream.str();
    return status;
}
