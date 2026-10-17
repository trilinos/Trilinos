// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <PyEnvData.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <pybind11/pybind11.h>
#include <Kokkos_Core.hpp>

#include <stdexcept>
#include <string>

namespace pykrino {

namespace py = pybind11;

static void check_mpi(int ierr)
{
    if (ierr != MPI_SUCCESS)
    {
        char errstr[MPI_MAX_ERROR_STRING];
        int len = 0;

        MPI_Error_string(ierr, errstr, &len);

        throw std::runtime_error(
            std::string("MPI error: ") + std::string(errstr, len));
    }
}

namespace {

struct ProcessState
{
    bool isInitialized = false;
    bool doIOwnKokkos = false;
    MPI_Comm comm = MPI_COMM_NULL;
    int numLiveMeshes = 0;
};

ProcessState & process_state()
{
    static ProcessState state;
    return state;
}

MPI_Comm get_comm_from_mpi4py(pybind11::object comm)
{
    py::object fcomm_obj = comm.attr("py2f")();
    MPI_Fint fcomm = fcomm_obj.cast<MPI_Fint>();
    MPI_Comm inputComm = MPI_Comm_f2c(fcomm);
    if (inputComm == MPI_COMM_NULL)
    {
        throw std::runtime_error("Received MPI_COMM_NULL");
    }
    return inputComm;
}

}

void ensure_initialized(MPI_Comm comm)
{
    ProcessState & state = process_state();
    if (state.isInitialized)
        return;

    int mpiInitialized = 0;
    check_mpi(MPI_Initialized(&mpiInitialized));
    if (!mpiInitialized)
    {
        throw std::runtime_error("MPI has not been initialized (import mpi4py.MPI before using pykrino).");
    }

    if (Kokkos::is_finalized())
    {
        throw std::runtime_error("Kokkos has already been finalized in this process, so no more meshes can be built.");
    }
    if (!Kokkos::is_initialized())
    {
        Kokkos::initialize();
        state.doIOwnKokkos = true;
    }

    check_mpi(MPI_Comm_dup(comm, &state.comm));
    stk::EnvData::instance().initialize(state.comm);
    state.isInitialized = true;
}

MPI_Comm get_initialized_comm()
{
    ensure_initialized();
    return stk::EnvData::parallel_comm();
}

void finalize()
{
    ProcessState & state = process_state();
    if (!state.isInitialized)
        return;

    int mpiFinalized = 0;
    MPI_Finalized(&mpiFinalized);
    if (!mpiFinalized && state.comm != MPI_COMM_NULL)
    {
        MPI_Comm_free(&state.comm);
        state.comm = MPI_COMM_NULL;
    }

    if (state.doIOwnKokkos && state.numLiveMeshes == 0 && Kokkos::is_initialized())
    {
        Kokkos::finalize();
        state.doIOwnKokkos = false;
    }
    state.isInitialized = false;
}

void increment_num_live_meshes()
{
    ++process_state().numLiveMeshes;
}

void decrement_num_live_meshes()
{
    --process_state().numLiveMeshes;
}

PyEnvData::PyEnvData()
{
    ensure_initialized(MPI_COMM_SELF);
}

PyEnvData::PyEnvData(pybind11::object comm)
{
    ensure_initialized(get_comm_from_mpi4py(comm));
}

int PyEnvData::rank() const
{
    return stk::EnvData::parallel_rank();
}

int PyEnvData::size() const
{
    return stk::EnvData::parallel_size();
}

void bind_env_data(py::module_& m) {
    py::class_<PyEnvData>(m, "EnvData")
        .def(py::init<>(), "Initialize pykrino on MPI_COMM_SELF (kept for compatibility; see pykrino.initialize)")
        .def(py::init<py::object>(),
             py::arg("comm"),
             "Initialize pykrino on an mpi4py communicator (kept for compatibility; see pykrino.initialize)")
        .def("rank", &PyEnvData::rank)
        .def("size", &PyEnvData::size);

    m.def("initialize",
          [](py::object comm)
          {
              if (comm.is_none())
                  ensure_initialized();
              else
                  ensure_initialized(get_comm_from_mpi4py(comm));
          },
          py::arg("comm") = py::none(),
          "Set up Kokkos and the STK communicator (a duplicate of comm, an mpi4py communicator, default COMM_WORLD).  "
          "Optional: done automatically by the first mesh.  Needed only to use another communicator, and must "
          "then be called before the first mesh.");
    m.def("finalize", &finalize,
          "Release pykrino's communicator and finalize Kokkos, if no mesh is alive.  Called automatically at exit.");
    m.def("parallel_rank", []() { ensure_initialized(); return stk::EnvData::parallel_rank(); },
          "Rank of this processor in pykrino's communicator.");
    m.def("parallel_size", []() { ensure_initialized(); return stk::EnvData::parallel_size(); },
          "Number of processors in pykrino's communicator.");
}

} // namespace pykrino
