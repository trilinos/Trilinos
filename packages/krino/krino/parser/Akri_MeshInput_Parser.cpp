// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_MeshInput_Parser.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_MeshInputOptions.hpp>
#include <Akri_Phase_Parser.hpp>

#include <stk_util/environment/RuntimeDoomed.hpp>
#include "Akri_Parser.hpp"

namespace krino {

void
MeshInput_Parser::parse(const Parser::Node & base_node)
{
  const Parser::Node fem_nodes = base_node.get_sequence_if_present("finite_element_models");
  if (!fem_nodes) return;

  for ( auto && fem_node : fem_nodes )
  {
    std::string model_name;
    fem_node.get_if_present("name", model_name);
    if (model_name.empty())
    {
      stk::RuntimeDoomedAdHoc() << "Missing finite element model name.\n";
    }
    std::shared_ptr<MeshInputOptions> options = MeshInputOptions::get_or_create(model_name);

    const bool has_generated_mesh = parse_generated_mesh(fem_node, *options);

    std::string mesh_name;
    fem_node.get_if_present("mesh", mesh_name);
    if (!mesh_name.empty() && !has_generated_mesh)
    {
      options->set_filename(mesh_name);
    }
    else if (!has_generated_mesh)
    {
      stk::RuntimeDoomedAdHoc() << "Must specify input mesh name or generated mesh options.\n";
    }

    std::string decomposition_method;
    if (fem_node.get_if_present("decomposition_method", decomposition_method))
    {
      // no checking here for validity
      std::transform(decomposition_method.begin(), decomposition_method.end(), decomposition_method.begin(), ::toupper);
      options->set_decomposition_method(decomposition_method);
    }

    bool useAllSidesForShells = true;
    if (fem_node.get_if_present("use_all_sides_for_shells", useAllSidesForShells))
    {
      options->set_use_all_sides_for_shells(useAllSidesForShells);
    }

    Phase_Parser::parse(fem_node, model_name);
  }
}

bool
MeshInput_Parser::parse_generated_mesh(const Parser::Node & fem_node, MeshInputOptions & options)
{
  const Parser::Node generated_mesh_node = fem_node.get_map_if_present("generated_mesh");
  if (!generated_mesh_node) return false;

  double mesh_size = 0.0;
  if (!generated_mesh_node.get_if_present("mesh_size", mesh_size) || mesh_size < 0.0)
  {
    stk::RuntimeDoomedAdHoc() << "Missing or invalid mesh_size for generated_mesh.\n";
  }
  options.set_generated_mesh_size(mesh_size);

  std::vector<double> domain;
  if (generated_mesh_node.get_if_present("domain", domain))
  {
    if (options.get_generated_mesh_domain_type() != MeshInputOptions::NO_GENERATED_MESH)
    {
      stk::RuntimeDoomedAdHoc() << "Must only specify domain or interface_bounding_box_with_dimension options.\n";
    }
    options.set_generated_mesh_domain_type(MeshInputOptions::GENERATED_MESH_FOR_SPECIFIED_DOMAIN);
    if (domain.size() != 4 && domain.size() != 6)
    {
      stk::RuntimeDoomedAdHoc() << "Domain must be a vector of length 4 for 2D or 6 for 3D (xmin,ymin,zmin, xmax,ymax,zmax).\n";
    }
    options.set_generated_mesh_domain(domain);
  }

  int interface_spatial_dim = 0;
  if (generated_mesh_node.get_if_present("interface_bounding_box_with_dimension", interface_spatial_dim))
  {
    if (options.get_generated_mesh_domain_type() != MeshInputOptions::NO_GENERATED_MESH)
    {
      stk::RuntimeDoomedAdHoc() << "Must only specify domain or interface_bounding_box options.\n";
    }
    if (interface_spatial_dim != 2 && interface_spatial_dim != 3)
    {
      stk::RuntimeDoomedAdHoc() << "interface_bounding_box_with_dimension only support 2 or 3 dimensions.\n";
    }
    const MeshInputOptions::GeneratedMeshDomainType mesh_type =
        interface_spatial_dim == 2 ?
        MeshInputOptions::GENERATED_2D_MESH_FOR_INTERFACE_BOUNDING_BOX :
        MeshInputOptions::GENERATED_3D_MESH_FOR_INTERFACE_BOUNDING_BOX;
    options.set_generated_mesh_domain_type(mesh_type);
  }

  std::string generated_mesh_element_type_string;
  if (generated_mesh_node.get_if_present("element_type", generated_mesh_element_type_string))
  {
    std::transform(generated_mesh_element_type_string.begin(), generated_mesh_element_type_string.end(), generated_mesh_element_type_string.begin(), ::toupper);
    static std::map<std::string, stk::topology> valid_entries =
      { {"TRIANGLE", stk::topology::TRIANGLE_3_2D},
        {"TRI", stk::topology::TRIANGLE_3_2D},
        {"QUADRILATERAL", stk::topology::QUADRILATERAL_4_2D},
        {"QUAD", stk::topology::QUADRILATERAL_4_2D},
        {"HEXAHEDRON", stk::topology::HEXAHEDRON_8},
        {"HEX", stk::topology::HEXAHEDRON_8},
        {"TETRAHEDRON", stk::topology::TETRAHEDRON_4},
        {"TET", stk::topology::TETRAHEDRON_4} };
    auto it = valid_entries.find(generated_mesh_element_type_string);
    if (it == valid_entries.end())
    {
      stk::RuntimeDoomedAdHoc() << "Invalid cdfem_simplex_generation_method type: " << generated_mesh_node.info();
    }
    else
    {
      options.set_generated_mesh_element_type( it->second );
      if (options.get_generated_mesh_spatial_dimension() != (int)options.get_generated_mesh_element_type().dimension())
      {
        stk::RuntimeDoomedAdHoc() << "Mismatch in spatial dimension for generated mesh element type and domain specification. ";
      }
    }
  }
  else
  {
    stk::topology default_topology =
        (options.get_generated_mesh_spatial_dimension() == 2) ?
        stk::topology::TRIANGLE_3_2D :
        stk::topology::TETRAHEDRON_4;
    options.set_generated_mesh_element_type( default_topology );
  }

  std::string generated_mesh_type_string;
  if (generated_mesh_node.get_if_present("mesh_type", generated_mesh_type_string))
  {
    std::transform(generated_mesh_type_string.begin(), generated_mesh_type_string.end(), generated_mesh_type_string.begin(), ::toupper);
    static std::map<std::string, BoundingBoxMeshStructureType> valid_entries =
      { {"CUBIC", BoundingBoxMeshStructureType::CUBIC_BOUNDING_BOX_MESH},
        {"BCC", BoundingBoxMeshStructureType::BCC_BOUNDING_BOX_MESH},
        {"FLAT_WALLED_BCC", BoundingBoxMeshStructureType::FLAT_WALLED_BCC_BOUNDING_BOX_MESH},
        {"TRIANGULAR_LATTICE", BoundingBoxMeshStructureType::TRIANGULAR_LATTICE_BOUNDING_BOX_MESH},
        {"FLAT_WALLED_TRIANGULAR_LATTICE", BoundingBoxMeshStructureType::FLAT_WALLED_TRIANGULAR_LATTICE_BOUNDING_BOX_MESH}};
    auto it = valid_entries.find(generated_mesh_type_string);
    if (it == valid_entries.end())
    {
      stk::RuntimeDoomedAdHoc() << "Invalid mesh_type: " << generated_mesh_node.info();
    }
    else
    {
      options.set_generated_mesh_structure_type(it->second);
      if ((BoundingBoxMeshStructureType::BCC_BOUNDING_BOX_MESH == it->second || BoundingBoxMeshStructureType::FLAT_WALLED_BCC_BOUNDING_BOX_MESH == it->second) && options.get_generated_mesh_spatial_dimension() != 3)
      {
        stk::RuntimeDoomedAdHoc() << "BCC meshes only supported in 3D";
      }
      if ((BoundingBoxMeshStructureType::TRIANGULAR_LATTICE_BOUNDING_BOX_MESH == it->second || BoundingBoxMeshStructureType::FLAT_WALLED_TRIANGULAR_LATTICE_BOUNDING_BOX_MESH == it->second) && options.get_generated_mesh_spatial_dimension() != 2)
      {
        stk::RuntimeDoomedAdHoc() << "BCC meshes only supported in 2D";
      }
    }
  }

  if (options.get_generated_mesh_domain_type() == MeshInputOptions::NO_GENERATED_MESH)
  {
    stk::RuntimeDoomedAdHoc() << "Must specify domain or interface_bounding_box_with_dimension for generated_mesh.\n";
  }
  return true;
}


} // namespace krino
