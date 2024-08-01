// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_GLOBAL_COMM_HPP
#define BELOS_GLOBAL_COMM_HPP

#ifdef HAVE_MPI

#include <mpi.h>
#include <mutex>

namespace Belos {

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
#endif // BELOS_GLOBAL_COMM_HPP
