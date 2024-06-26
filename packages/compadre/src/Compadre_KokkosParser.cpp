// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#include "Compadre_KokkosParser.hpp"

using namespace Compadre;

// for InitArguments, pass them directly in to Kokkos
KokkosParser::KokkosParser(KokkosInitArguments args, bool print_status) {
    this->ksg = !Kokkos::is_initialized()
#ifdef COMPADRE_KOKKOS_GREATEREQUAL_3_7
                && !Kokkos::is_finalized()
#endif
                ?
                new Kokkos::ScopeGuard(args) : nullptr;
    if (print_status) this->status();
}

// for command line arguments, pass them directly in to Kokkos
KokkosParser::KokkosParser(int narg, char* args[], bool print_status) {
    this->ksg = !Kokkos::is_initialized()
#ifdef COMPADRE_KOKKOS_GREATEREQUAL_3_7
                && !Kokkos::is_finalized()
#endif
                ?
                new Kokkos::ScopeGuard(narg, args) : nullptr;
    if (print_status) this->status();
}

KokkosParser::KokkosParser(std::vector<std::string> stdvec_args, bool print_status) {
    std::vector<char*> char_args;
    for (const auto& arg : stdvec_args) {
        char_args.push_back((char*)arg.data());
    }
    char_args.push_back(nullptr);
    int narg = (int)stdvec_args.size();

    this->ksg = !Kokkos::is_initialized()
#ifdef COMPADRE_KOKKOS_GREATEREQUAL_3_7
                && !Kokkos::is_finalized()
#endif
                ?
                new Kokkos::ScopeGuard(narg, char_args.data()) : nullptr;
    if (print_status) this->status();
}

KokkosParser::KokkosParser(bool print_status) : KokkosParser(KokkosInitArguments(), print_status) {}

std::string KokkosParser::status() {
    std::stringstream stream;
    Kokkos::print_configuration(stream, true);
    std::string status = stream.str();
    return status;
}
