// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <PyRefinement.hpp>
#include <PyMesh.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_MeshSizeMetric.hpp>
#include <Akri_MetricEdgeMarking.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_util/diag/Timer.hpp>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <string>
#include <stdexcept>

namespace pykrino {

namespace py = pybind11;

PyRefinement::PyRefinement(PyMesh &mesh, bool force64Bit, bool assert32Bit)
: myMesh(mesh)
{
    stk::mesh::Part * activePart = nullptr;
    myRefinement = std::make_shared<krino::Refinement>(mesh.get().meta_data(), activePart, force64Bit, assert32Bit, sierra::Diag::sierraTimer());
    
    if (!myRefinement) {
        throw std::runtime_error("Failed to create krino::Refinement object");
    }
}

bool PyRefinement::do_uniform_refinement(int numUniformRefinementLevels)
{
    if (!myRefinement) {
        throw std::runtime_error("Refinement object is not initialized");
    }
    
    return myRefinement->do_uniform_refinement(numUniformRefinementLevels);
}

void PyRefinement::delete_parent_elements()
{
    if (!myRefinement) {
        throw std::runtime_error("Refinement object is not initialized");
    }

    myRefinement->delete_parent_elements();
}

static std::vector<double> to_point(const stk::math::Vector3d & x, const unsigned dim)
{
    return std::vector<double>(x.data(), x.data()+dim);
}

size_t PyRefinement::refine_edges_too_long(const krino::MeshSizeMetric & sizeMetric, double maxRelativeEdgeLength, unsigned maxNumPasses)
{
    if (!myRefinement) {
        throw std::runtime_error("Refinement object is not initialized");
    }

    const stk::mesh::BulkData & mesh = myMesh.get().bulk_data();
    const krino::FieldRef coordsField(mesh.mesh_meta_data().coordinate_field());
    return krino::refine_edges_too_long(mesh, *myRefinement, coordsField, sizeMetric, maxRelativeEdgeLength, maxNumPasses);
}

size_t PyRefinement::refine_edges_too_long_in_metric(const MetricFunction & metricFunction, double maxRelativeEdgeLength, unsigned maxNumPasses)
{
    const stk::mesh::BulkData & mesh = myMesh.get().bulk_data();
    const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
    const krino::MetricTensorFunctionMeshSizeMetric sizeMetric(mesh, krino::FieldRef(mesh.mesh_meta_data().coordinate_field()),
        [&metricFunction, dim](const stk::math::Vector3d & x, double * metric)
        {
            const std::vector<double> values = metricFunction(to_point(x, dim));
            if (values.size() != dim*dim)
                throw std::runtime_error("The metric function must return " + std::to_string(dim*dim) + " values, not " + std::to_string(values.size()) + ".");
            std::copy(values.begin(), values.end(), metric);
        });
    return refine_edges_too_long(sizeMetric, maxRelativeEdgeLength, maxNumPasses);
}

size_t PyRefinement::refine_edges_too_long_in_size(const SizeFunction & sizeFunction, double maxRelativeEdgeLength, unsigned maxNumPasses)
{
    const stk::mesh::BulkData & mesh = myMesh.get().bulk_data();
    const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
    const krino::IsotropicSizeFunctionMeshSizeMetric sizeMetric(mesh, krino::FieldRef(mesh.mesh_meta_data().coordinate_field()),
        [&sizeFunction, dim](const stk::math::Vector3d & x) { return sizeFunction(to_point(x, dim)); });
    return refine_edges_too_long(sizeMetric, maxRelativeEdgeLength, maxNumPasses);
}

void bind_refinement(py::module_ &m)
{
    py::class_<PyRefinement, std::shared_ptr<PyRefinement>>(m, "Refinement")
        .def(
            py::init<
                PyMesh &,
                bool,
                bool
            >(),
            py::arg("mesh"),
            py::arg("force_64bit") = false,
            py::arg("assert_32bit") = false,
            "Create a Refinement object from a mesh"
        )
        .def("do_uniform_refinement",
             &PyRefinement::do_uniform_refinement,
             py::arg("num_uniform_refinement_levels"),
             "Perform uniform refinement")
        .def("delete_parent_elements",
             &PyRefinement::delete_parent_elements,
             "Delete parent elements after refinement")
        .def("refine_edges_too_long_in_metric",
             &PyRefinement::refine_edges_too_long_in_metric,
             py::arg("metric_function"),
             py::arg("max_relative_edge_length") = std::sqrt(2.),
             py::arg("max_num_passes") = 10,
             "Refine every edge that is too long in a metric, repeating until none is or max_num_passes passes have "
             "been done.  metric_function(x) gives the metric M at the point x (spatial_dimension coordinates) as "
             "spatial_dimension*spatial_dimension values in row-major order; M must be symmetric positive definite.  "
             "An edge is too long if its length sqrt(e^T M e), averaged over its two nodes, exceeds "
             "max_relative_edge_length.  Returns the number of edges still too long (0 if converged).  Parent "
             "elements are kept, so call delete_parent_elements() before smoothing.  The function must not raise "
             "on only some processors.")
        .def("refine_edges_too_long_in_size",
             &PyRefinement::refine_edges_too_long_in_size,
             py::arg("size_function"),
             py::arg("max_relative_edge_length") = std::sqrt(2.),
             py::arg("max_num_passes") = 10,
             "As refine_edges_too_long_in_metric, for an isotropic size: size_function(x) gives the target edge "
             "length h at the point x, and an edge is too long if its length over h exceeds max_relative_edge_length.");
}

} // namespace pykrino
