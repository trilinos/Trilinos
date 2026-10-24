// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#ifndef KRINO_PYKRINO_SRC_PYMESH_HPP_
#define KRINO_PYKRINO_SRC_PYMESH_HPP_

#include <pybind11/pybind11.h>
#include <Akri_MeshInterface.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <memory>

namespace pykrino {

class PyMesh
{
public:
    explicit PyMesh(std::shared_ptr<krino::MeshInterface> mesh);
    virtual ~PyMesh() = default;
    void populate_mesh();
    const krino::MeshInterface & get() const { return *myMesh; }
    krino::MeshInterface & get() { return *myMesh; }
    void write_mesh(const std::string & fileName, bool compose = true);
protected:
    std::shared_ptr<krino::MeshInterface> myMesh;
    // Shared by the copies of this PyMesh; counts the mesh as live until the last copy is gone (see finalize).
    // Objects that hold STK data of the mesh (Refinement, mesh objectives) keep a PyMesh copy, so they are counted.
    std::shared_ptr<void> myLiveMeshCount;
};

class PyMeshFromFile : public PyMesh
{
public:
    PyMeshFromFile(
        const std::string & fileName,
        const std::string & decompMethod = "",
        const bool useAllSidesForShells = true,
        const bool readFields = true);

    PyMeshFromFile(const PyMeshFromFile &) = delete;
    PyMeshFromFile &operator=(const PyMeshFromFile &) = delete;
};

class PyBoundingBoxMesh : public PyMesh
{
public:
    PyBoundingBoxMesh(
        const std::string & elementType,
        const std::vector<double> & minCoords,
        const std::vector<double> & maxCoords,
        double meshSize,
        int padCells = 0
    );

    PyBoundingBoxMesh(const PyBoundingBoxMesh &) = delete;
    PyBoundingBoxMesh &operator=(const PyBoundingBoxMesh &) = delete;

    void set_mesh_structure_type(const std::string & structureType);

private:
    explicit PyBoundingBoxMesh(std::shared_ptr<krino::BoundingBoxMesh> bboxMesh);
    stk::topology get_element_topology(const std::string & elementType) const;
    static krino::BoundingBoxMeshStructureType get_structure_type(const std::string & structureType);

    std::shared_ptr<krino::BoundingBoxMesh> myBboxMesh;
};

void bind_mesh(pybind11::module_ &m);

} // namespace pykrino

#endif /* KRINO_PYKRINO_SRC_PYMESH_HPP_ */
