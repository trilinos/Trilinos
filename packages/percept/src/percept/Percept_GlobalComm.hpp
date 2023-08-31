// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef percept_global_comm_hpp
#define percept_global_comm_hpp

#include <mpi.h>
#include <mutex>

namespace percept {

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

#endif /* percept_global_comm_hpp */
