/*
 * GeometryKernelPGEOM.hpp
 *
 *  Created on: Jun 16, 2016
 *      Author: madbrew
 */

#if HAVE_CUBIT

#ifndef PERCEPT_SRC_PERCEPT_MESH_GEOMETRY_KERNEL_GEOMETRYKERNELPGEOM_HPP_
#define PERCEPT_SRC_PERCEPT_MESH_GEOMETRY_KERNEL_GEOMETRYKERNELPGEOM_HPP_

#include "GeometryKernel.hpp"

class PGeom;

namespace percept {

// Ensures that a vector of nodes is both able to be seen by a given process and
// that each node is either owned by or shared to the current process.
// If a single node cannot be validated, then the function stops checking and returns false
// This function is used to optimize mesh-to-geometry association in parallel
bool validate_node_ownership (const stk::mesh::BulkData *mesh, const std::vector<int> &nodeIDsToCheck);
bool validate_element_ownership (const stk::mesh::BulkData *mesh, const int &elementIDToCheck);

bool get_node_from_id(const stk::mesh::BulkData *mesh, int node_id,
                      stk::mesh::Entity &node);
  
bool get_beam_from_ids(const stk::mesh::BulkData *mesh,
                       std::vector<int> edge_node_ids, stk::mesh::Entity &edge);

bool get_shell_from_ids(const stk::mesh::BulkData *mesh,
                        std::vector<int> face_node_ids, stk::mesh::Entity &quad);

// UNUSED FUNCTION WILL PROBABLY DISCARD
bool get_shell_from_ids_LIBERAL(const stk::mesh::BulkData *mesh,
                                std::vector<int> face_node_ids, stk::mesh::Entity &shell);

bool get_hex_from_id(const stk::mesh::BulkData *mesh, int hex_id,
                     stk::mesh::Entity &hex);

bool get_hexes_from_node_ids(const stk::mesh::BulkData *mesh, std::vector<stk::mesh::Entity> &elemsReturned, const std::vector<stk::mesh::Entity> &nodeToSearch);

class GeometryKernelPGEOM: public GeometryKernel 
{

public:
  GeometryKernelPGEOM();
  GeometryKernelPGEOM(PerceptMesh &mesh, std::string geomToImp);
  virtual ~GeometryKernelPGEOM();
  
  bool create_assoc(PerceptMesh &mesh, std::string geomToImp);
  
  virtual bool read_file(const std::string& file_name,
                         std::vector<GeometryHandle>& geometry_entities); /* creates a vector of geometry handles populated by reading info of from
			 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 geometry file. These will be read by get_attribute to create names
			 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 for geom. pieces in the GeometryFactory. */

  virtual void snap_to(KernelPoint& point, GeometryHandle geom,	/* will snap a given point to a given piece of geometry*/
                       double *converged_tolerance = NULL, double *uvw_computed = NULL,
                       double *uvw_hint = NULL, void *extra_hint = NULL);
  
  virtual void normal_at(KernelPoint& point, GeometryHandle geom,
                         std::vector<double>& normal, void *extra_hint = NULL); /* tells where a given point has a normal with a given piece of geometry (!?)
			 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 what if it isn't normal to that piece of geometry?*/

private:  
  PGeom* m_pg;
};

}
#endif /* PERCEPT_SRC_PERCEPT_MESH_GEOMETRY_KERNEL_GEOMETRYKERNELPGEOM_HPP_ */

#endif
