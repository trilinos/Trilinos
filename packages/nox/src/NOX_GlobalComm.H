// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_GLOBAL_COMM_H
#define NOX_GLOBAL_COMM_H

#ifdef HAVE_MPI

#include <mpi.h>
#include <mutex>

namespace NOX {

static std::mutex mpi_mutex;
static MPI_Comm Global_MPI_Comm = MPI_COMM_WORLD;

inline void initialize_global_comm(MPI_Comm comm) {
    std::lock_guard<std::mutex> guard(mpi_mutex);
    Global_MPI_Comm = comm;
}

inline MPI_Comm get_global_comm() {
    std::lock_guard<std::mutex> guard(mpi_mutex);
    return Global_MPI_Comm;
}

}

#endif /* HAVE_MPI */
#endif /* NOX_GLOBAL_COMM_H */