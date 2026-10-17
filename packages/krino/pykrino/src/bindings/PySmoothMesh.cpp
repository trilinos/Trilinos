// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <PySmoothMesh.hpp>
#include <PyMesh.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <Akri_Optimize.hpp>
#include <Akri_ROLOptimize.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_KinematicMeshUtils.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_ModifiedNeoHookeanModel.hpp>
#include <Akri_Objective_KinematicMesh.hpp>
#include <Akri_Objective_KinematicModel.hpp>
#include <Akri_Objective_MeanRatioNorm.hpp>
#include <Akri_SethHillModel.hpp>
#include <Akri_SizeShapeModel.hpp>

#include <string>
#include <stdexcept>
#include <memory>

namespace pykrino {

namespace py = pybind11;

// PyGradientFilter implementation
PyGradientFilter::PyGradientFilter() {
    // Default empty filter (does nothing)
    myFilter = [](const stk::mesh::Entity /*node*/,
                  const stk::math::Vector3d & /*nodeCoords*/,
                  stk::math::Vector3d & /*dir*/) {
        // No-op by default
    };
}

PyGradientFilter PyGradientFilter::create_cube_boundary_filter(double min_coord, double max_coord) {
    PyGradientFilter filter;

    filter.myFilter = [min_coord, max_coord](const stk::mesh::Entity /*node*/,
                                            const stk::math::Vector3d &nodeCoords,
                                            stk::math::Vector3d &dir) {
        // Zero out movement directions for nodes on cube faces
        for (unsigned i = 0; i < 3; ++i) {
            if (nodeCoords[i] <= min_coord + 1e-6 || nodeCoords[i] >= max_coord - 1e-6) {
                dir[i] = 0.0;
            }
        }
    };

    return filter;
}

krino::GradientFilter PyGradientFilter::get() const {
    return myFilter;
}

// PyMeshObjectiveBase implementation extension
void PyMeshObjectiveBase::set_gradient_filter(const PyGradientFilter & filter) {
    if (myMeshObj) {
        myMeshObj->set_gradient_filter(filter.get());
    }
}

PyOptimizer::PyOptimizer(const Type type,
    const double stepTol,
    const double gradTol,
    const unsigned maxIter,
    const bool doOutput)
{
  myOptFn = build_optimizer_function(type, stepTol, gradTol, maxIter, doOutput);
}

krino::OptimizerFunction PyOptimizer::build_optimizer_function(const Type type,
    const double stepTol,
    const double gradTol,
    const unsigned maxIter,
    const bool doOutput)
{
    switch (type) {
        case Type::LBFGS:
          return [=](const krino::ObjectiveInterface &objFn, krino::DistributedVector &x)
          {
              const double initialObjective = objFn.compute_value(x);
              krino::lbfgs(objFn, x, stepTol, gradTol, initialObjective, maxIter, krino::OptimizationDefaults::maxLevelsBFGS, doOutput);
          };
        case Type::SteepestDescent:
          return [=](const krino::ObjectiveInterface &objFn, krino::DistributedVector &x)
          {
              const double initialObjective = objFn.compute_value(x);
              krino::steepest_descent(objFn, x, stepTol, gradTol, initialObjective, maxIter, doOutput);
          };
        case Type::ROL:
          return [=](const krino::ObjectiveInterface &objFn, krino::DistributedVector &x)
          {
            const double initialObjective = objFn.compute_value(x);
              krino::rol_optimize(objFn, x, stepTol, gradTol, initialObjective, maxIter, doOutput);
          };
        default:
            throw std::runtime_error("Unknown optimizer type");
    }
}

PyElementObjective::PyElementObjective(std::shared_ptr<krino::ElementObjective> obj)
  : myElemObj(std::move(obj))
{
    if (!myElemObj) {
        throw std::runtime_error("Cannot construct PyElementObjective from null krino::ElementObjective");
    }
}

const krino::ElementObjective & PyElementObjective::get() const
{
    return *myElemObj;
}

PyKinematicElementObjective::PyKinematicElementObjective(std::shared_ptr<krino::KinematicElementObjective> elemObj)
  : myElemObj(std::move(elemObj))
{
    if (!myElemObj) {
        throw std::runtime_error("Cannot construct PyKinematicElementObjective from null krino::KinematicElementObjective");
    }
}

const krino::KinematicElementObjective & PyKinematicElementObjective::get() const
{
    return *myElemObj;
}

PyMeshObjectiveBase::PyMeshObjectiveBase(PyMesh mesh, std::shared_ptr<krino::MeshObjectiveBase> meshObj)
  : myMesh(std::move(mesh)),
    myMeshObj(std::move(meshObj))
{
    if (!myMeshObj) {
        throw std::runtime_error("Cannot construct PyMeshObjectiveBase from null krino::MeshObjectiveBase");
    }
}

static krino::CoordinatesFieldRef get_coordinates_field(const PyMesh & mesh)
{
    const stk::mesh::MetaData & meta = mesh.get().meta_data();
    return krino::CoordinatesFieldRef(meta.coordinate_field(), meta.spatial_dimension());
}

const krino::MeshObjectiveBase & PyMeshObjectiveBase::get() const
{
    return *myMeshObj;
}

std::vector<stk::mesh::EntityId> PyMeshObjectiveBase::get_owned_element_ids() const
{
    const stk::mesh::BulkData & bulk = myMesh.get().bulk_data();
    std::vector<stk::mesh::EntityId> elemIds;
    elemIds.reserve(myMeshObj->get_owned_elements().size());
    for (const stk::mesh::Entity elem : myMeshObj->get_owned_elements())
        elemIds.push_back(bulk.identifier(elem));
    return elemIds;
}

std::pair<std::vector<stk::mesh::EntityId>, std::vector<double>> PyMeshObjectiveBase::get_node_ids_and_coordinates() const
{
    const stk::mesh::BulkData & bulk = myMesh.get().bulk_data();
    const krino::CoordinatesFieldRef coordsField = get_coordinates_field(myMesh);
    const unsigned dim = coordsField.dim();
    const auto & nodes = myMeshObj->get_nodes();

    std::vector<stk::mesh::EntityId> nodeIds;
    std::vector<double> coords;
    nodeIds.reserve(nodes.size());
    coords.reserve(nodes.size()*dim);
    for (size_t n = 0; n < nodes.size(); ++n)
    {
        nodeIds.push_back(bulk.identifier(nodes[n]));
        const stk::math::Vector3d nodeCoords = krino::get_vector_field(bulk, coordsField, nodes[n], dim);
        for (unsigned d = 0; d < dim; ++d)
            coords.push_back(nodeCoords[d]);
    }
    return {nodeIds, coords};
}

PyMeanRatioElementObjective::PyMeanRatioElementObjective()
: PyElementObjective(std::make_shared<krino::MeanRatioNormObjective>())
{}

// SethHillParams implementation
PySethHillParams::PySethHillParams(double mu, double kappa, int m, int n)
  : myParams{mu, kappa, m, n}
{
}

krino::SethHillParams PySethHillParams::get() const {
  return myParams;
}

// SethHillElementObjective implementation
PySethHillElementObjective::PySethHillElementObjective()
: PyKinematicElementObjective(std::make_shared<krino::KinematicModelObjective<krino::SimplifiedSethHillModel>>())
{
}

PySethHillElementObjective::PySethHillElementObjective(const PySethHillParams & params)
: PyKinematicElementObjective(std::make_shared<krino::KinematicModelObjective<krino::SethHillModel>>(krino::SethHillModel(params.get())))
{
}

// ModifiedNeoHookeanElementObjective implementation
PyModifiedNeoHookeanElementObjective::PyModifiedNeoHookeanElementObjective()
: PyKinematicElementObjective(std::make_shared<krino::KinematicModelObjective<krino::ModifiedNeoHookeanModel>>())
{
}

PyMeshObjective::PyMeshObjective(
    PyMesh mesh,
    PyElementObjective elemObj)
: PyMeshObjectiveBase(mesh, std::make_shared<krino::MeshObjective>(
    mesh.get().bulk_data(),
    get_coordinates_field(mesh),
    mesh.get().meta_data().universal_part(),
    elemObj.get())),
  myElemObj(elemObj)
{}

PyKinematicMeshObjective::PyKinematicMeshObjective(
    PyMesh mesh,
    std::shared_ptr<krino::KinematicMeshObjective> meshObj,
    PyKinematicElementObjective elemObj)
: PyMeshObjectiveBase(std::move(mesh), std::move(meshObj)),
  myElemObj(std::move(elemObj))
{}

PyKinematicConstantIsotropicMeshObjective::PyKinematicConstantIsotropicMeshObjective(
    PyMesh mesh,
    PyKinematicElementObjective elemObj,
    double refSize)
: PyKinematicMeshObjective(mesh,
    std::make_shared<krino::KinematicConstantIsotropicMeshObjective>(
        mesh.get().bulk_data(), get_coordinates_field(mesh), mesh.get().meta_data().universal_part(), elemObj.get(), refSize),
    elemObj)
{}

PyKinematicIsotropicMeshObjective::PyKinematicIsotropicMeshObjective(
    PyMesh mesh,
    PyKinematicElementObjective elemObj)
: PyKinematicIsotropicMeshObjective(mesh,
    std::make_shared<krino::KinematicIsotropicMeshObjective>(
        mesh.get().bulk_data(), get_coordinates_field(mesh), mesh.get().meta_data().universal_part(), elemObj.get()),
    elemObj)
{}

PyKinematicIsotropicMeshObjective::PyKinematicIsotropicMeshObjective(
    PyMesh mesh,
    PyKinematicElementObjective elemObj,
    krino::SmoothReference reference)
: PyKinematicIsotropicMeshObjective(mesh, elemObj)
{
    myIsotropicMeshObj->set_element_ref_sizes(krino::compute_element_ref_sizes_from_reference(
        mesh.get().bulk_data(), get_coordinates_field(mesh), myIsotropicMeshObj->get_owned_elements(), reference));
}

PyKinematicIsotropicMeshObjective::PyKinematicIsotropicMeshObjective(
    PyMesh mesh,
    std::shared_ptr<krino::KinematicIsotropicMeshObjective> meshObj,
    PyKinematicElementObjective elemObj)
: PyKinematicMeshObjective(std::move(mesh), meshObj, std::move(elemObj)),
  myIsotropicMeshObj(std::move(meshObj))
{}

void PyKinematicIsotropicMeshObjective::set_element_ref_sizes(const std::vector<double> & elemRefSizes)
{
    myIsotropicMeshObj->set_element_ref_sizes(elemRefSizes);
}

PyKinematicAnisotropicMeshObjective::PyKinematicAnisotropicMeshObjective(
    PyMesh mesh,
    PyKinematicElementObjective elemObj)
: PyKinematicAnisotropicMeshObjective(mesh,
    std::make_shared<krino::KinematicAnisotropicMeshObjective>(
        mesh.get().bulk_data(), get_coordinates_field(mesh), mesh.get().meta_data().universal_part(), elemObj.get()),
    elemObj)
{}

PyKinematicAnisotropicMeshObjective::PyKinematicAnisotropicMeshObjective(
    PyMesh mesh,
    PyKinematicElementObjective elemObj,
    const std::string & nodalMetricFieldName)
: PyKinematicAnisotropicMeshObjective(mesh, elemObj)
{
    set_element_metric_map_data_from_nodal_metric_field(nodalMetricFieldName);
}

PyKinematicAnisotropicMeshObjective::PyKinematicAnisotropicMeshObjective(
    PyMesh mesh,
    std::shared_ptr<krino::KinematicAnisotropicMeshObjective> meshObj,
    PyKinematicElementObjective elemObj)
: PyKinematicMeshObjective(std::move(mesh), meshObj, std::move(elemObj)),
  myAnisotropicMeshObj(std::move(meshObj))
{}

void PyKinematicAnisotropicMeshObjective::set_element_metric_map_data(const std::vector<double> & elemMetricMapData)
{
    myAnisotropicMeshObj->set_element_metric_map_data(elemMetricMapData);
}

void PyKinematicAnisotropicMeshObjective::set_element_metric_map_data_from_nodal_metric_field(const std::string & nodalMetricFieldName)
{
    const stk::mesh::MetaData & meta = myMesh.get().meta_data();
    const krino::AuxMetaData & auxMeta = krino::AuxMetaData::get(meta);
    if (!auxMeta.has_field(stk::topology::NODE_RANK, nodalMetricFieldName))
        throw std::runtime_error("There is no nodal field named " + nodalMetricFieldName + " for the metric.");
    const krino::FieldRef metricField = auxMeta.get_field(stk::topology::NODE_RANK, nodalMetricFieldName);
    myAnisotropicMeshObj->set_element_metric_map_data(krino::compute_element_metric_mappings_from_nodal_metric_field(
        myMesh.get().bulk_data(), metricField, myAnisotropicMeshObj->get_owned_elements()));
}

void PyKinematicAnisotropicMeshObjective::set_element_metric_map_data_from_nodal_metrics(const std::vector<double> & nodalMetrics)
{
    myAnisotropicMeshObj->set_element_metric_map_data_from_nodal_metrics(nodalMetrics);
}

PySmoothMesh::PySmoothMesh() {}

void PySmoothMesh::smooth_mesh(
    PyMesh mesh,
    PyOptimizer optimizer,
    PyMeshObjectiveBase meshObj,
    bool print_diagnostics)
{
    auto optimizer_func = optimizer.get_optimizer_function();
    
    const stk::mesh::MetaData & meta = mesh.get().meta_data();
    const stk::mesh::BulkData & bulk = mesh.get().bulk_data();

    krino::smooth_mesh_by_optimization(
        bulk,
        krino::CoordinatesFieldRef(meta.coordinate_field(), meta.spatial_dimension()),
        optimizer_func,
        meshObj.get(),
        print_diagnostics
    );
}

void bind_smooth_mesh(py::module_ &m)
{
    py::class_<PyOptimizer> optimizer(m, "Optimizer");

    py::enum_<pykrino::PyOptimizer::Type>(optimizer, "Type")
        .value("LBFGS", pykrino::PyOptimizer::Type::LBFGS)
        .value("SteepestDescent", pykrino::PyOptimizer::Type::SteepestDescent)
        .value("ROL", pykrino::PyOptimizer::Type::ROL)
        .export_values();

    optimizer.def(py::init<PyOptimizer::Type, double, double, unsigned, bool>(),
        py::arg("type"),
        py::arg("stepTol") = krino::OptimizationDefaults::stepTol,
        py::arg("gradTol") = krino::OptimizationDefaults::gradTol,
        py::arg("maxIter") = krino::OptimizationDefaults::maxIter,
        py::arg("doOutput") = krino::OptimizationDefaults::doOutput);

    py::class_<PyElementObjective, std::shared_ptr<PyElementObjective>>(
        m, "Element Objective");

    py::class_<PyMeanRatioElementObjective,
               PyElementObjective,
               std::shared_ptr<PyMeanRatioElementObjective>>(
        m, "MeanRatioElementObjective")
        .def(py::init<>(), "Mean ratio norm objective");

    py::class_<PySethHillParams>(m, "SethHillParams")
        .def(py::init<double, double, int, int>(),
             py::arg("mu") = 0.0,
             py::arg("kappa") = 0.0,
             py::arg("m") = 0,
             py::arg("n") = 0,
             "Create SethHill parameters");

    py::class_<PyKinematicElementObjective,
               std::shared_ptr<PyKinematicElementObjective>>(
        m, "Element Kinematic Objective");

    py::class_<PySethHillElementObjective,
    PyKinematicElementObjective,
               std::shared_ptr<PySethHillElementObjective>>(
        m, "SethHillElementObjective")
        .def(py::init<>(), "SethHill objective with default parameters")
        .def(py::init<const PySethHillParams &>(),
             py::arg("params"),
             "SethHill objective with custom parameters");

    py::class_<PyModifiedNeoHookeanElementObjective,
               PyKinematicElementObjective,
               std::shared_ptr<PyModifiedNeoHookeanElementObjective>>(
        m, "ModifiedNeoHookeanElementObjective")
        .def(py::init<>(), "Modified Neo-Hookean objective for mesh smoothing");

    py::class_<PyMeshObjectiveBase, std::shared_ptr<PyMeshObjectiveBase>>(
        m, "Mesh Objective Base Class")
        .def("set_gradient_filter",
             &PyMeshObjectiveBase::set_gradient_filter,
             py::arg("filter"),
             "Set a gradient filter on this mesh objective")
        .def("get_owned_element_ids",
             &PyMeshObjectiveBase::get_owned_element_ids,
             "Global ids of the locally owned elements, in the order used for per-element data")
        .def("get_node_ids_and_coordinates",
             &PyMeshObjectiveBase::get_node_ids_and_coordinates,
             "Global ids of the nodes used by the locally owned elements (owned, then shared but not owned) "
             "and their coordinates (spatial dimension values per node), in the order used for per-node data.");

    py::class_<PyGradientFilter>(m, "GradientFilter")
        .def(py::init<>(), "Create an empty gradient filter")
        .def_static("create_cube_boundary_filter",
            &PyGradientFilter::create_cube_boundary_filter,
            py::arg("min_coord") = -1.0,
            py::arg("max_coord") = 1.0,
            "Create a cube boundary gradient filter");

    py::class_<PyMeshObjective,
               PyMeshObjectiveBase,
               std::shared_ptr<PyMeshObjective>>(
        m, "MeshObjective")
        .def(py::init<PyMesh, PyElementObjective>(),
             py::arg("mesh"),
             py::arg("element_objective"),
             "Create a mesh objective from a mesh and element objective.");

    py::enum_<krino::SmoothReference>(m, "SmoothReference")
        .value("EqualVolume", krino::SmoothReference::EqualVolume)
        .value("AvgEdgeLength", krino::SmoothReference::AvgEdgeLength)
        .value("Max", krino::SmoothReference::Max);

    py::class_<PyKinematicMeshObjective,
               PyMeshObjectiveBase,
               std::shared_ptr<PyKinematicMeshObjective>>(
        m, "KinematicMeshObjective");

    py::class_<PyKinematicConstantIsotropicMeshObjective,
               PyKinematicMeshObjective,
               std::shared_ptr<PyKinematicConstantIsotropicMeshObjective>>(
        m, "KinematicConstantIsotropicMeshObjective")
        .def(py::init<PyMesh, PyKinematicElementObjective, double>(),
             py::arg("mesh"),
             py::arg("kinematic_element_objective"),
             py::arg("ref_size"),
             "Kinematic mesh objective with the same ideal element size for every element.");

    py::class_<PyKinematicIsotropicMeshObjective,
               PyKinematicMeshObjective,
               std::shared_ptr<PyKinematicIsotropicMeshObjective>>(
        m, "KinematicIsotropicMeshObjective")
        .def(py::init<PyMesh, PyKinematicElementObjective>(),
             py::arg("mesh"),
             py::arg("kinematic_element_objective"),
             "Kinematic mesh objective with per-element ideal sizes, to be set with set_element_ref_sizes.")
        .def(py::init<PyMesh, PyKinematicElementObjective, krino::SmoothReference>(),
             py::arg("mesh"),
             py::arg("kinematic_element_objective"),
             py::arg("reference"),
             "Kinematic mesh objective with per-element ideal sizes computed from the current mesh.")
        .def("set_element_ref_sizes",
             &PyKinematicIsotropicMeshObjective::set_element_ref_sizes,
             py::arg("elem_ref_sizes"),
             "Set the ideal sizes of the locally owned elements, in the order of get_owned_element_ids.");

    py::class_<PyKinematicAnisotropicMeshObjective,
               PyKinematicMeshObjective,
               std::shared_ptr<PyKinematicAnisotropicMeshObjective>>(
        m, "KinematicAnisotropicMeshObjective")
        .def(py::init<PyMesh, PyKinematicElementObjective>(),
             py::arg("mesh"),
             py::arg("kinematic_element_objective"),
             "Kinematic mesh objective with a per-element anisotropic metric, to be set with one of the "
             "set_element_metric_map_data methods.")
        .def(py::init<PyMesh, PyKinematicElementObjective, const std::string &>(),
             py::arg("mesh"),
             py::arg("kinematic_element_objective"),
             py::arg("nodal_metric_field_name"),
             "Kinematic mesh objective with the metric of each element from a nodal metric field of the mesh "
             "(see set_element_metric_map_data_from_nodal_metric_field).")
        .def("set_element_metric_map_data",
             &PyKinematicAnisotropicMeshObjective::set_element_metric_map_data,
             py::arg("elem_metric_map_data"),
             "Set the metric mapping F_M (with metric M = F_M^T F_M) of the locally owned elements, in the order of "
             "get_owned_element_ids, spatial_dimension*spatial_dimension values per element in row-major order.")
        .def("set_element_metric_map_data_from_nodal_metric_field",
             &PyKinematicAnisotropicMeshObjective::set_element_metric_map_data_from_nodal_metric_field,
             py::arg("nodal_metric_field_name"),
             "Set the metric of each element from a nodal field holding the metric M at each node, as a "
             "spatial_dimension*spatial_dimension matrix in row-major order.  The metric of each element is the "
             "log-Euclidean mean of the metrics at its nodes.  Each metric must be symmetric positive definite.")
        .def("set_element_metric_map_data_from_nodal_metrics",
             &PyKinematicAnisotropicMeshObjective::set_element_metric_map_data_from_nodal_metrics,
             py::arg("nodal_metrics"),
             "Set the metric of each element from the metric M at each node of get_node_ids_and_coordinates, in "
             "that order, spatial_dimension*spatial_dimension values per node in row-major order.  The metric "
             "of each element is the log-Euclidean mean of the metrics at its nodes.  Each metric must be "
             "symmetric positive definite.");

    py::class_<PySmoothMesh>(m, "SmoothMesh")
        .def(py::init<>(), "Create a mesh smoothing object")
        .def("smooth_mesh",
             &PySmoothMesh::smooth_mesh,
             py::arg("mesh"),
             py::arg("optimizer"),
             py::arg("mesh_objective"),
             py::arg("print_diagnostics") = false,
             "Smooth a mesh using specified optimizer and objective");
}

} // namespace pykrino
