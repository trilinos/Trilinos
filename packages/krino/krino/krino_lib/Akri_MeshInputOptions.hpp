// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_MeshInputOptions_h
#define Akri_MeshInputOptions_h

#include <Akri_BoundingBoxMesh.hpp>
#include <string>
#include <Ioss_PropertyManager.h>
#include <memory>
#include <set>
#include <stk_topology/topology.hpp>

namespace krino {

  class MeshInputOptions {

  public:
    static std::shared_ptr<MeshInputOptions> get_or_create(const std::string & model_name);
    static MeshInputOptions * get(const std::string & model_name);

    enum GeneratedMeshDomainType
    {
      NO_GENERATED_MESH=0,
      GENERATED_MESH_FOR_SPECIFIED_DOMAIN,
      GENERATED_2D_MESH_FOR_INTERFACE_BOUNDING_BOX,
      GENERATED_3D_MESH_FOR_INTERFACE_BOUNDING_BOX
    };

    MeshInputOptions(const std::string & name) :
      my_name(name),
      my_filename("%B.g"),
      my_filetype("exodusII"),
      my_generated_mesh_domain_type(NO_GENERATED_MESH),
      myGeneratedMeshStructureType(BoundingBoxMeshStructureType::CUBIC_BOUNDING_BOX_MESH),
      my_generated_mesh_size(-1.0) {}

    const std::string & get_name() const { return my_name; }
    bool is_valid() const
    {
      return !my_name.empty() && ((! my_filetype.empty() && ! my_filename.empty()) || use_generated_mesh()) ;
    }

    void set_generated_mesh_structure_type(BoundingBoxMeshStructureType type) { myGeneratedMeshStructureType = type; }
    BoundingBoxMeshStructureType get_generated_mesh_structure_type() const { return myGeneratedMeshStructureType; }

    bool use_generated_mesh() const { return my_generated_mesh_domain_type != NO_GENERATED_MESH; }

    void set_generated_mesh_domain_type(GeneratedMeshDomainType type) { my_generated_mesh_domain_type = type; }
    GeneratedMeshDomainType get_generated_mesh_domain_type() const { return my_generated_mesh_domain_type; }

    int get_generated_mesh_spatial_dimension() const;

    void set_generated_mesh_size(double size) { my_generated_mesh_size = size; }
    double get_generated_mesh_size() const { return my_generated_mesh_size; }

    void set_generated_mesh_element_type(stk::topology top) { my_generated_mesh_element_type = top; }
    stk::topology get_generated_mesh_element_type() const { return my_generated_mesh_element_type; }

    void set_generated_mesh_domain(const std::vector<double> & domain) { my_generated_mesh_domain = domain; }
    const std::vector<double> & get_generated_mesh_domain() const { return my_generated_mesh_domain; }

    void set_filename(std::string filename) { my_filename = filename; }
    const std::string & get_filename() const { return my_filename; }

    void set_filetype(std::string type) { my_filetype = type; }
    const std::string & get_filetype() const { return my_filetype; }

    void set_use_all_sides_for_shells(bool useAllSidesForShells) { myUseAllSidesForShells = useAllSidesForShells; }
    bool get_use_all_sides_for_shells() const { return myUseAllSidesForShells; }

    void set_decomposition_method(std::string decomposition_method) { my_decomposition_method = decomposition_method; }
    const std::string & get_decomposition_method() const { return my_decomposition_method; }

    void set_coordinate_system(std::string coordinate_system) { my_coordinate_system = coordinate_system; }
    const std::string & get_coordinate_system() const { return my_coordinate_system; }

    void add_property(const Ioss::Property & property) { my_properties.add(property); }
    Ioss::PropertyManager & get_properties() { return my_properties; }

  private:
    struct compare_ptr {
      bool operator()( std::shared_ptr<MeshInputOptions> const lhs ,
                       std::shared_ptr<MeshInputOptions> const rhs ) const
      {
        return lhs.get() == nullptr ? rhs.get() != nullptr :
          ( rhs.get() == nullptr ? false : lhs->get_name() < rhs->get_name() );
      }
    };

    typedef std::set<std::shared_ptr<MeshInputOptions>, compare_ptr> Registry;
    static Registry the_registry;

  private:
    std::string my_name;
    std::string my_filename;
    std::string my_filetype;
    bool myUseAllSidesForShells{true};
    GeneratedMeshDomainType my_generated_mesh_domain_type;
    BoundingBoxMeshStructureType myGeneratedMeshStructureType;
    std::vector<double> my_generated_mesh_domain;
    stk::topology my_generated_mesh_element_type = stk::topology::INVALID_TOPOLOGY;
    double my_generated_mesh_size;
    std::string my_decomposition_method;
    std::string my_coordinate_system;
    Ioss::PropertyManager my_properties;
  };

} // namespace krino

#endif /* Akri_MeshInputOptions_h */
