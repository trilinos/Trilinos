// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef GEOMETRYFACTORY_HPP
#define GEOMETRYFACTORY_HPP

#include <string>
#include <percept/mesh/geometry/kernel/GeometryKernel.hpp>
#include <percept/mesh/geometry/kernel/MeshGeometry.hpp>
#include <percept/PerceptMesh.hpp>

namespace percept {

class GeometryFactory
{
public:
    GeometryFactory(GeometryKernel* kernel, MeshGeometry* geometry);
    virtual ~GeometryFactory();

    bool read_file(const std::string& filename, stk::mesh::MetaData* meta_data);
    bool read_file(const std::string& filename, PerceptMesh *mesh)
    {
      return read_file(filename, mesh->get_fem_meta_data());
    }

protected:
    GeometryKernel* geomKernel;
    MeshGeometry* geomDatabase;
};
}
#endif // GEOMETRYFACTORY_HPP
