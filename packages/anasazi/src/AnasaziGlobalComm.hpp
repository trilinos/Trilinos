// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ANASAZI_GLOBAL_COMM_HPP
#define ANASAZI_GLOBAL_COMM_HPP

#include "AnasaziConfigDefs.hpp"

#ifdef HAVE_MPI

#include <mpi.h>
#include <mutex>

namespace Anasazi {

static std::mutex mpi_mutex;
static MPI_Comm Global_Anasazi_Comm = MPI_COMM_WORLD;

inline void initialize_global_comm(MPI_Comm comm) {
  std::lock_guard<std::mutex> guard(mpi_mutex);
  Global_Anasazi_Comm = comm;
}

inline MPI_Comm get_global_comm() {
  std::lock_guard<std::mutex> guard(mpi_mutex);
  return Global_Anasazi_Comm;
}
} // namespace Anasazi

#endif // HAVE_MPI

#endif // ANASAZI_GLOBAL_COMM_HPP