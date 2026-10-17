// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <PyMesh.hpp>
#include <PyEnvData.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <Akri_MeshFromFile.hpp>
#include <Akri_OutputUtils.hpp>

#include <string>
#include <stdexcept>

namespace pykrino {

namespace py = pybind11;

PyMesh::PyMesh(std::shared_ptr<krino::MeshInterface> mesh)
  : myMesh(std::move(mesh))
{
    if (!myMesh) {
        throw std::runtime_error("Cannot construct PyMesh from null krino::MeshInterface");
    }
    increment_num_live_meshes();
    myLiveMeshCount = std::shared_ptr<void>(nullptr, [](void *) { decrement_num_live_meshes(); });
}

void PyMesh::populate_mesh()
{
    if (!myMesh) {
        throw std::runtime_error("Mesh object is not initialized");
    }

    myMesh->populate_mesh();
}

void PyMesh::write_mesh(const std::string & filename, bool compose)
{
  if (compose) {
    krino::output_composed_mesh_with_fields(
        myMesh->bulk_data(), myMesh->meta_data().universal_part(),
        filename, 1, 0.);
  } else {
    krino::output_mesh_with_fields(
        myMesh->bulk_data(), myMesh->meta_data().universal_part(),
        filename, 1, 0.);
  }
}

PyMeshFromFile::PyMeshFromFile(
    const std::string & fileName,
    const std::string & decompMethod,
    const bool useAllSidesForShells,
    const bool readFields)
: PyMesh(std::make_shared<krino::MeshFromFile>(
        fileName,
        get_initialized_comm(),
        decompMethod,
        useAllSidesForShells,
        readFields))
{
}

PyBoundingBoxMesh::PyBoundingBoxMesh(
    const std::string & elementType,
    const std::vector<double> & minCoords,
    const std::vector<double> & maxCoords,
    double meshSize,
    int padCells)
: PyBoundingBoxMesh(std::make_shared<krino::BoundingBoxMesh>(
      get_element_topology(elementType),
      get_initialized_comm()))
{
    if (minCoords.size() != 3 || maxCoords.size() != 3) {
        throw std::runtime_error("minCoords and maxCoords must each have 3 elements");
    }

    krino::BoundingBox_T<double, 3> bbox(
        stk::math::Vector3d(minCoords[0], minCoords[1], minCoords[2]),
        stk::math::Vector3d(maxCoords[0], maxCoords[1], maxCoords[2])
    );

    auto bboxMesh = std::dynamic_pointer_cast<krino::BoundingBoxMesh>(myMesh);
    if (!bboxMesh) {
        throw std::runtime_error("Failed to cast to BoundingBoxMesh");
    }

    myBboxMesh->set_domain(bbox, meshSize, padCells);
}

PyBoundingBoxMesh::PyBoundingBoxMesh(
    std::shared_ptr<krino::BoundingBoxMesh> bboxMesh)
: PyMesh(bboxMesh),
  myBboxMesh(bboxMesh)
{}

stk::topology PyBoundingBoxMesh::get_element_topology(const std::string & elementType) const
{
    static const std::map<std::string, stk::topology> topologyMap = {
        {"TRIANGLE_3_2D", stk::topology::TRIANGLE_3_2D},
        {"QUADRILATERAL_4_2D", stk::topology::QUADRILATERAL_4_2D},
        {"TETRAHEDRON_4", stk::topology::TETRAHEDRON_4},
        {"HEXAHEDRON_8", stk::topology::HEXAHEDRON_8},
        {"TRI_3", stk::topology::TRIANGLE_3_2D},
        {"QUAD_4", stk::topology::QUADRILATERAL_4_2D},
        {"TET_4", stk::topology::TETRAHEDRON_4},
        {"HEX_8", stk::topology::HEXAHEDRON_8}
    };

    auto it = topologyMap.find(elementType);
    if (it == topologyMap.end()) {
        throw std::runtime_error("Unsupported element type: " + elementType);
    }
    return it->second;
}

krino::BoundingBoxMeshStructureType
PyBoundingBoxMesh::get_structure_type(const std::string & structureType)
{
    static const std::map<std::string, krino::BoundingBoxMeshStructureType> structureMap = {
        {"CUBIC", krino::CUBIC_BOUNDING_BOX_MESH},
        {"BCC", krino::BCC_BOUNDING_BOX_MESH},
        {"FLAT_WALLED_BCC", krino::FLAT_WALLED_BCC_BOUNDING_BOX_MESH},
        {"TRIANGULAR_LATTICE", krino::TRIANGULAR_LATTICE_BOUNDING_BOX_MESH},
        {"FLAT_WALLED_TRIANGULAR_LATTICE", krino::FLAT_WALLED_TRIANGULAR_LATTICE_BOUNDING_BOX_MESH}
    };

    auto it = structureMap.find(structureType);
    if (it == structureMap.end()) {
        throw std::runtime_error("Unsupported mesh structure type: " + structureType);
    }

    return it->second;
}

void PyBoundingBoxMesh::set_mesh_structure_type(const std::string & structureType)
{
    myBboxMesh->set_mesh_structure_type(get_structure_type(structureType));
}

void bind_mesh(py::module_ &m)
{
    py::class_<PyMesh, std::shared_ptr<PyMesh>>(m, "Mesh")
        .def("write_mesh",
             &PyMesh::write_mesh,
             py::arg("filename"),
             py::arg("compose") = true,
             "write mesh to file"
            )
        .def(
            "populate_mesh",
            &PyMeshFromFile::populate_mesh,
            "Populate the mesh data structures"
        );
        
    py::class_<PyMeshFromFile, PyMesh, std::shared_ptr<PyMeshFromFile>>(m, "MeshFromFile")
        .def(
            py::init<
                const std::string &,
                const std::string &,
                bool,
                bool
            >(),
            py::arg("filename"),
            py::arg("decomp_method") = "",
            py::arg("use_all_sides_for_shells") = true,
            py::arg("read_fields") = true,
            "Create a MeshFromFile object from an input mesh file"
        );

    py::class_<PyBoundingBoxMesh, PyMesh, std::shared_ptr<PyBoundingBoxMesh>>(m, "BoundingBoxMesh")
        .def(
            py::init<
                const std::string &,
                const std::vector<double> &,
                const std::vector<double> &,
                double,
                int
            >(),
            py::arg("element_type"),
            py::arg("min_coords"),
            py::arg("max_coords"),
            py::arg("mesh_size"),
            py::arg("pad_cells") = 0,
            "Create a BoundingBoxMesh object"
        )
        .def(
            "set_mesh_structure_type",
            &PyBoundingBoxMesh::set_mesh_structure_type,
            py::arg("structure_type"),
            "Set the mesh structure type (CUBIC, BCC, FLAT_WALLED_BCC, TRIANGULAR_LATTICE, FLAT_WALLED_TRIANGULAR_LATTICE)"
        );
}

} // namespace pykrino
