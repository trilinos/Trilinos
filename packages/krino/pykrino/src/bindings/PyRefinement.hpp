// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#ifndef KRINO_PYKRINO_SRC_PYREFINEMENT_HPP_
#define KRINO_PYKRINO_SRC_PYREFINEMENT_HPP_

#include <pybind11/pybind11.h>
#include <Akri_Refinement.hpp>
#include <PyMesh.hpp>
#include <functional>
#include <memory>
#include <vector>

namespace krino { class MeshSizeMetric; }

namespace pykrino {

class PyRefinement
{
public:
    // Metric M (spatial_dimension*spatial_dimension values, row-major) or isotropic size h at a point.
    using MetricFunction = std::function<std::vector<double>(const std::vector<double> &)>;
    using SizeFunction = std::function<double(const std::vector<double> &)>;

    PyRefinement(PyMesh &mesh, bool force64Bit = false, bool assert32Bit = false);
    virtual ~PyRefinement() = default;
    
    const krino::Refinement & get() const { return *myRefinement; }
    
    bool do_uniform_refinement(int numUniformRefinementLevels);
    void delete_parent_elements();

    // Refines until no edge is too long in the metric or size, or maxNumPasses; returns the number still too long.
    size_t refine_edges_too_long_in_metric(const MetricFunction & metricFunction, double maxRelativeEdgeLength, unsigned maxNumPasses);
    size_t refine_edges_too_long_in_size(const SizeFunction & sizeFunction, double maxRelativeEdgeLength, unsigned maxNumPasses);
    
protected:
    size_t refine_edges_too_long(const krino::MeshSizeMetric & sizeMetric, double maxRelativeEdgeLength, unsigned maxNumPasses);

    PyMesh myMesh;
    std::shared_ptr<krino::Refinement> myRefinement;
};

void bind_refinement(pybind11::module_ &m);

} // namespace pykrino

#endif /* KRINO_PYKRINO_SRC_PYREFINEMENT_HPP_ */
