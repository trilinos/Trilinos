// bindings.cpp - Main entry point for pykrino Python module
// This file is referenced by the CMakeLists.txt and serves as the entry point
// for the pybind11 module. The actual bindings are defined in PySurface.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <PySurface.hpp>
#include <PyEnvData.hpp>
#include <PyMesh.hpp>
#include <PySmoothMesh.hpp>
#include <PyLogRedirecter.hpp>
#include <PyRefinement.hpp>

namespace py = pybind11;

PYBIND11_MODULE(_pykrino, m)
{
    pykrino::bind_surfaces(m);
    pykrino::bind_env_data(m);
    pykrino::bind_mesh(m);
    pykrino::bind_smooth_mesh(m);
    pykrino::bind_log_redirecter(m);
    pykrino::bind_refinement(m);
}
