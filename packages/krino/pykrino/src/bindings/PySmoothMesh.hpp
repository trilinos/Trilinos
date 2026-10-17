// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#ifndef KRINO_PYKRINO_SRC_PYSMOOTHMESH_HPP_
#define KRINO_PYKRINO_SRC_PYSMOOTHMESH_HPP_

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <string>
#include <utility>
#include <vector>
#include <PyMesh.hpp>
#include <Akri_SmoothMesh.hpp>
#include <Akri_Objective_KinematicMesh.hpp>
#include <Akri_Objective_MeanRatioNorm.hpp>
#include <Akri_ModifiedNeoHookeanModel.hpp>
#include <Akri_Objective_KinematicModel.hpp>
#include <Akri_SethHillModel.hpp>
#include <Akri_SizeShapeModel.hpp>
#include <Akri_Optimize.hpp>
#include <Akri_SethHillConfig.hpp>

namespace pykrino {

// Wrapper for optimizer functions
class PyOptimizer {
public:
    enum class Type { LBFGS, SteepestDescent, ROL };

    PyOptimizer(const Type type,
        const double stepTol = krino::OptimizationDefaults::stepTol,
        const double gradTol = krino::OptimizationDefaults::gradTol,
        const unsigned maxIter = krino::OptimizationDefaults::maxIter,
        const bool doOutput = krino::OptimizationDefaults::doOutput);

    krino::OptimizerFunction get_optimizer_function() const { return myOptFn; }

private:
    static krino::OptimizerFunction build_optimizer_function(const Type type,
        const double stepTol = krino::OptimizationDefaults::stepTol,
        const double gradTol = krino::OptimizationDefaults::gradTol,
        const unsigned maxIter = krino::OptimizationDefaults::maxIter,
        const bool doOutput = krino::OptimizationDefaults::doOutput);
    krino::OptimizerFunction myOptFn;
};

class PyElementObjective {
public:
    explicit PyElementObjective(std::shared_ptr<krino::ElementObjective> obj);
    const krino::ElementObjective & get() const;

protected:
    std::shared_ptr<krino::ElementObjective> myElemObj;
};

class PyKinematicElementObjective {
public:
    explicit PyKinematicElementObjective(std::shared_ptr<krino::KinematicElementObjective> obj);
    const krino::KinematicElementObjective & get() const;

protected:
    std::shared_ptr<krino::KinematicElementObjective> myElemObj;
};

class PyMeanRatioElementObjective : public PyElementObjective {
public:
    PyMeanRatioElementObjective();
};

// SethHill parameter structure
class PySethHillParams {
public:
    PySethHillParams(double mu = 0.0, double kappa = 0.0, int m = 0, int n = 0);

    krino::SethHillParams get() const;

private:
    krino::SethHillParams myParams;
};

class PySethHillElementObjective : public PyKinematicElementObjective {
public:
    PySethHillElementObjective();
    PySethHillElementObjective(const PySethHillParams & params);
};

class PyModifiedNeoHookeanElementObjective : public PyKinematicElementObjective {
public:
    PyModifiedNeoHookeanElementObjective();
};

class PyGradientFilter {
public:
    PyGradientFilter();

    // Create a cube boundary filter (like build_node_search_direction_filter_for_2x2x2_cube)
    static PyGradientFilter create_cube_boundary_filter(double min_coord = -1.0, double max_coord = 1.0);

    krino::GradientFilter get() const;

private:
    krino::GradientFilter myFilter;
};

class PyMeshObjectiveBase {
public:
    PyMeshObjectiveBase(PyMesh mesh, std::shared_ptr<krino::MeshObjectiveBase> meshObj);
    const krino::MeshObjectiveBase & get() const;

    void set_gradient_filter(const PyGradientFilter & filter);
    std::vector<stk::mesh::EntityId> get_owned_element_ids() const;
    // Global ids of the nodes of the objective (locally owned, then shared but not owned), and their
    // coordinates, spatial_dimension per node, in the order used for per-node data.
    std::pair<std::vector<stk::mesh::EntityId>, std::vector<double>> get_node_ids_and_coordinates() const;

protected:
    PyMesh myMesh;
    std::shared_ptr<krino::MeshObjectiveBase> myMeshObj;
};

class PyMeshObjective : public PyMeshObjectiveBase {
public:
    PyMeshObjective(
        PyMesh mesh,
        PyElementObjective elemObj);

private:
    PyElementObjective myElemObj;
};

class PyKinematicMeshObjective : public PyMeshObjectiveBase {
protected:
    PyKinematicMeshObjective(
        PyMesh mesh,
        std::shared_ptr<krino::KinematicMeshObjective> meshObj,
        PyKinematicElementObjective elemObj);

private:
    PyKinematicElementObjective myElemObj;
};

class PyKinematicConstantIsotropicMeshObjective : public PyKinematicMeshObjective {
public:
    PyKinematicConstantIsotropicMeshObjective(
        PyMesh mesh,
        PyKinematicElementObjective elemObj,
        double refSize);
};

class PyKinematicIsotropicMeshObjective : public PyKinematicMeshObjective {
public:
    PyKinematicIsotropicMeshObjective(
        PyMesh mesh,
        PyKinematicElementObjective elemObj);
    PyKinematicIsotropicMeshObjective(
        PyMesh mesh,
        PyKinematicElementObjective elemObj,
        krino::SmoothReference reference);

    void set_element_ref_sizes(const std::vector<double> & elemRefSizes);

private:
    PyKinematicIsotropicMeshObjective(
        PyMesh mesh,
        std::shared_ptr<krino::KinematicIsotropicMeshObjective> meshObj,
        PyKinematicElementObjective elemObj);

    std::shared_ptr<krino::KinematicIsotropicMeshObjective> myIsotropicMeshObj;
};

class PyKinematicAnisotropicMeshObjective : public PyKinematicMeshObjective {
public:
    PyKinematicAnisotropicMeshObjective(
        PyMesh mesh,
        PyKinematicElementObjective elemObj);
    PyKinematicAnisotropicMeshObjective(
        PyMesh mesh,
        PyKinematicElementObjective elemObj,
        const std::string & nodalMetricFieldName);

    void set_element_metric_map_data(const std::vector<double> & elemMetricMapData);
    void set_element_metric_map_data_from_nodal_metric_field(const std::string & nodalMetricFieldName);
    void set_element_metric_map_data_from_nodal_metrics(const std::vector<double> & nodalMetrics);

private:
    PyKinematicAnisotropicMeshObjective(
        PyMesh mesh,
        std::shared_ptr<krino::KinematicAnisotropicMeshObjective> meshObj,
        PyKinematicElementObjective elemObj);

    std::shared_ptr<krino::KinematicAnisotropicMeshObjective> myAnisotropicMeshObj;
};

// Main smoothing wrapper class
class PySmoothMesh {
public:
    PySmoothMesh();

    void smooth_mesh(
        PyMesh mesh,
        PyOptimizer optimizer,
        PyMeshObjectiveBase meshObj,
        bool print_diagnostics = false
    );

private:
    // Helper function to create gradient filter
    static krino::GradientFilter create_gradient_filter(pybind11::function filter_func);

    // Helper function to create solution projector
    static krino::SolutionProjector create_solution_projector(pybind11::function projector_func);
};

void bind_smooth_mesh(pybind11::module_ &m);

} // namespace pykrino

#endif /* KRINO_PYKRINO_SRC_PYSMOOTHMESH_HPP_ */
