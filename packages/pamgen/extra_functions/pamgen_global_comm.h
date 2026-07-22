// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef PAMGEN_GLOBAL_COMM_HPP
#define PAMGEN_GLOBAL_COMM_HPP

#ifdef HAVE_MPI

#include <mpi.h>
#include <mutex>

namespace PAMGEN_NEVADA {

static std::mutex mpi_mutex;
static MPI_Comm Global_MPI_Comm = MPI_COMM_WORLD; // CHECK: ALLOW MPI_COMM_WORLD

inline void initialize_global_comm(MPI_Comm comm) {
    std::lock_guard<std::mutex> guard(mpi_mutex);
    Global_MPI_Comm = comm;
}

inline MPI_Comm get_global_comm() {
    std::lock_guard<std::mutex> guard(mpi_mutex);
    return Global_MPI_Comm;
}

}

#endif // HAVE_MPI
#endif // PAMGEN_GLOBAL_COMM_HPP
