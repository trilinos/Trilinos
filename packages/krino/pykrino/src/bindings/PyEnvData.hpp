// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_PYKRINO_SRC_PYENV_DATA_HPP_
#define KRINO_PYKRINO_SRC_PYENV_DATA_HPP_

#include <pybind11/pybind11.h>
#include <stk_util/parallel/Parallel.hpp>

namespace pykrino {

// Process-wide setup needed before building any STK mesh: Kokkos and the stk::EnvData communicator
// (a duplicate of comm).  MPI must already be initialized (pykrino imports mpi4py, which does it).
// Later calls are no-ops until finalize.  Mesh constructors call this, so users need not.
void ensure_initialized(MPI_Comm comm = stk::parallel_machine_world());
// ensure_initialized, then the communicator to build meshes on.
MPI_Comm get_initialized_comm();

// Frees the communicator and finalizes Kokkos if pykrino initialized it.  Kokkos can be finalized only once per
// process and STK meshes cannot be destroyed after it, so Kokkos is left alone while any mesh is alive (meshes still
// alive at interpreter exit are destroyed later, during teardown).  Registered with atexit; safe to call twice.
void finalize();

// Count of live STK meshes created by pykrino, maintained by PyMesh.
void increment_num_live_meshes();
void decrement_num_live_meshes();

// Compatibility with scripts that construct pykrino.EnvData(comm) before building meshes: initializes with comm.
// No longer owns anything, so it may be deleted at any time.
class PyEnvData {
public:
    PyEnvData();
    explicit PyEnvData(pybind11::object comm);

    int rank() const;
    int size() const;
};

void bind_env_data(pybind11::module_& m);

} // namespace pykrino

#endif /* KRINO_PYKRINO_SRC_PYENV_DATA_HPP_ */
