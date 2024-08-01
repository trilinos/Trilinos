/*
//@HEADER
// *****************************************************************************
//                        Adelus
//
// Copyright 2020 NTESS and the Adelus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER
*/

#ifndef __ADELUS_GLOBALCOMM_H__
#define __ADELUS_GLOBALCOMM_H__

#include <mpi.h>
#include <mutex>

namespace Adelus {

static std::mutex mpi_mutex;
static MPI_Comm Global_Adelus_Comm = MPI_COMM_WORLD;

inline void initialize_global_comm(MPI_Comm comm) {
  std::lock_guard<std::mutex> guard(mpi_mutex);
  Global_Adelus_Comm = comm;
}

inline MPI_Comm get_global_comm() {
  std::lock_guard<std::mutex> guard(mpi_mutex);
  return Global_Adelus_Comm;
}

} // namespace Adelus

#endif // __ADELUS_GLOBALCOMM_H__
