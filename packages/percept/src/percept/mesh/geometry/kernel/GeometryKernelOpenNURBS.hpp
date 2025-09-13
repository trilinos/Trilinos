// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef GEOMETRYKERNEL_OPENNURBS_HPP
#define GEOMETRYKERNEL_OPENNURBS_HPP

#if HAVE_OPENNURBS

#include <opennurbs.h>
#include <percept/mesh/geometry/kernel/GeometryKernel.hpp>

namespace percept {

class GeometryKernelOpenNURBS : public GeometryKernel
{
public:
    GeometryKernelOpenNURBS();
    virtual ~GeometryKernelOpenNURBS();

    virtual bool read_file(const std::string& file_name,
                           std::vector<GeometryHandle>& geometry_entities) override;
    virtual bool debug_dump_file(const std::string& file_name) override;

    virtual void snap_to(KernelPoint& point, GeometryHandle geom,
                         double *converged_tolerance = NULL,
                         double *uvw_computed = NULL,
                         double *uvw_hint = NULL, void *extra_hint = NULL) override;

    virtual void normal_at(KernelPoint& point, GeometryHandle geom, std::vector<double>& normal, void *extra_hint = NULL) override;

  bool debug_is_curve(int geom) const;

  bool debug_is_surface(int geom) const;

private:
    ONX_Model onModel;
};
}

#endif

#endif // GEOMETRYKERNEL_OPENNURBS_HPP
