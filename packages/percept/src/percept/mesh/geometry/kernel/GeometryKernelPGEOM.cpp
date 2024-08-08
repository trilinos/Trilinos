
/*
 * GeometryKernelPGEOM.cpp
 *
 *  Created on: Jun 16, 2016
 *      Author: madbrew
 */

#if HAVE_CUBIT // entire file

#include <iostream>
#include <typeinfo>
#include <stdexcept>

#include "GeometryKernelPGEOM.hpp"
#include <string>

#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>

#include <PGeom.hpp>
#include <PGeomAssoc.hpp>

#ifdef HAVE_ACIS
#include <PGeomACIS.hpp>
#endif

#include <CubitVector.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

namespace percept {

  bool validate_node_ownership (const stk::mesh::BulkData *mesh, const std::vector<int> &nodeIDsToCheck)
  {
    stk::mesh::Selector localOrShared(mesh->mesh_meta_data().locally_owned_part());
    localOrShared |=mesh->mesh_meta_data().globally_shared_part();

    for(unsigned i=0;i<nodeIDsToCheck.size(); i++)
      {
        stk::mesh::Entity nodeToCheck = mesh->get_entity(stk::topology::NODE_RANK, (stk::mesh::EntityId) nodeIDsToCheck[i]);
        if (nodeToCheck==stk::mesh::Entity::InvalidEntity){
          return false; //First validate that you can find the entity at all on the current process
        }

        stk::mesh::Bucket &node_bucket = mesh->bucket(nodeToCheck);
        if (!localOrShared(node_bucket))
          return false; //then check to see if you can find it in the local or shared parts
      }

    return true; //if all the nodes pass validation, return true
  }

  bool validate_element_ownership (const stk::mesh::BulkData *mesh, const int &elementIDToCheck)
  {
    stk::mesh::Selector localOrShared(mesh->mesh_meta_data().locally_owned_part());
    localOrShared |=mesh->mesh_meta_data().globally_shared_part();

    stk::mesh::Entity elemToCheck = mesh->get_entity(stk::topology::ELEMENT_RANK, (stk::mesh::EntityId) elementIDToCheck);
    if (elemToCheck==stk::mesh::Entity::InvalidEntity){
      return false;
    }

    stk::mesh::Bucket &elem_bucket = mesh->bucket(elemToCheck);
    if (!localOrShared(elem_bucket))
      return false;

    return true;
  }

  bool get_node_from_id(const stk::mesh::BulkData *mesh, int node_id,
                        stk::mesh::Entity &node) {
    node = mesh->get_entity(stk::topology::NODE_RANK,
                            (stk::mesh::EntityId) node_id);
    return (node == stk::mesh::Entity::InvalidEntity) ? false : true;
  }

  bool get_beam_from_ids(const stk::mesh::BulkData *mesh,
                         const std::vector<int>& edge_node_ids, stk::mesh::Entity &beam)
  {
    std::vector<stk::mesh::Entity> edgeNodes(edge_node_ids.size());

    for (unsigned i = 0; i < edge_node_ids.size(); i++) {
      get_node_from_id(mesh, edge_node_ids[i], edgeNodes[i]);
    }

    unsigned numElements = mesh->num_elements(edgeNodes[0]);
    const stk::mesh::Entity * elements = mesh->begin_elements(edgeNodes[0]);

    for(unsigned e = 0; e < numElements; ++e) { //iterate over the elements in this mesh only
      stk::mesh::Entity beamToInspect = elements[e];
      stk::mesh::Bucket & node_bucket = mesh->bucket(beamToInspect);
      if (node_bucket.topology()!=stk::topology::BEAM_2) //if the current element isn't a beam, move to the next element
        continue;					 //	QUESTION: is there a more efficient means of checking an element's topology than returning an entire bucket?

      const MyPairIterRelation beam_nodes(*mesh, beamToInspect, stk::topology::NODE_RANK);
      bool found = true;
      for (unsigned jj = 0; jj < beam_nodes.size(); jj++) { //simply check if the nodes match. This check doesn't take permutations into account, only whether you can find all the nodes in
        if (beam_nodes[jj].entity() != edgeNodes[0]			// a given entity
            && beam_nodes[jj].entity() != edgeNodes[1]) {
          found = false; //if you can't match the nodes, move to the next element
          break;
        }
      }

      if (found) {
        beam = beamToInspect;
        return true;
      }
    }

    beam = stk::mesh::Entity::InvalidEntity;
    return false;
  } 

  bool get_shell_from_ids(const stk::mesh::BulkData *mesh,
                          const std::vector<int>& face_node_ids, stk::mesh::Entity &shell) 
  {
    std::vector<stk::mesh::Entity> faceNodes(face_node_ids.size());

    for (unsigned i = 0; i < face_node_ids.size(); i++) {
      get_node_from_id(mesh, face_node_ids[i], faceNodes[i]);
    }

    unsigned numElements1 = mesh->num_elements(faceNodes[0]);
    const stk::mesh::Entity * elements_1 = mesh->begin_elements(faceNodes[0]);

    for(unsigned e = 0; e < numElements1; ++e) { //loop over elements
      stk::mesh::Entity element = elements_1[e];

      stk::mesh::Bucket & node_bucket = mesh->bucket(element);
      if (node_bucket.topology()!=stk::topology::SHELL_QUAD_4 && node_bucket.topology()!=stk::topology::SHELL_TRI_3)
        continue; //if the current element isn't a shell, go to the next element
      //	QUESTION: is there a more efficient means of checking an element's topology than returning an entire bucket?

      const MyPairIterRelation face_nodes(*mesh, element, stk::topology::NODE_RANK);

      bool found = true;
      for (unsigned jj = 0; jj < face_nodes.size(); jj++) {
        //simply check if the nodes match. This check doesn't take permutations into account, only whether you can find all the nodes in
        // a given entity
        for (unsigned kk = 0; kk < face_nodes.size(); kk++) {
          if (face_nodes[jj].entity() != faceNodes[kk]) {
            found = false;
          }
          if (!found) break;
        }
      }

      if (found) {
        shell = element;
        return true;
      }
    }

    shell = stk::mesh::Entity::InvalidEntity;
    return false;
  }

  bool get_hex_from_id(const stk::mesh::BulkData *mesh, int hex_id,
                       stk::mesh::Entity &hex) 
  {
    hex = mesh->get_entity(stk::topology::ELEM_RANK,
                           (stk::mesh::EntityId) hex_id);
    return (hex == stk::mesh::Entity::InvalidEntity) ? false : true;
  } 

  bool get_hexes_from_node_ids(const stk::mesh::BulkData *mesh, std::vector<stk::mesh::Entity> &elemsReturned, const std::vector<stk::mesh::Entity> &nodesToSearch)
  {
    stk::mesh::Selector localAndShared(mesh->mesh_meta_data().locally_owned_part());
    localAndShared |= mesh->mesh_meta_data().globally_shared_part();

    for(unsigned i=0;i<nodesToSearch.size();i++) {
      stk::mesh::Entity node = nodesToSearch[i];
      
      stk::mesh::Bucket &node_bucket = mesh->bucket(node);
      if (!localAndShared(node_bucket))
        continue;
      
      unsigned numElements = mesh->num_elements(node);
      const stk::mesh::Entity * elements = mesh->begin_elements(node);
      
      for (unsigned e=0;e<numElements;++e) {
        stk::mesh::Entity element = elements[e];
        stk::mesh::Bucket &node_bucket2 = mesh->bucket(element);
        
        if (node_bucket2.topology()==stk::topology::SHELL_QUAD_4 || node_bucket2.topology() == stk::topology::SHELL_TRI_3)
          continue;
        
        if (elemsReturned.size()==0)
          elemsReturned.push_back(element);
        else {
          bool isIn = false;
          for (unsigned ee=0;ee<elemsReturned.size();ee++) {
            if( mesh->identifier(element) == mesh->identifier(elemsReturned[ee]) )
              isIn=true;
          }
          if (!isIn)
            elemsReturned.push_back(element); //only if you cannot find the element, add it to the elemntsreturned
        }
      }
    }
    
    if (elemsReturned.size()>0)
      return true;
    return false;
  }

  GeometryKernelPGEOM::GeometryKernelPGEOM() {
#ifdef HAVE_ACIS
    m_pg = new PGeomACIS;
#else
    m_pg = new PGeom;
#endif
  }

  GeometryKernelPGEOM::~GeometryKernelPGEOM(){delete m_pg;}

  bool GeometryKernelPGEOM::read_file(const std::string& file_name,
                                      std::vector<GeometryHandle>& geometry_entities) 
  {
    bool file_imported = false;
    if (file_name.find(".sat") != std::string::npos) {
#if HAVE_ACIS
      m_pg->initialize(ACIS_GEOMETRY_ENGINE);
      file_imported = m_pg->import_acis_file(file_name.c_str());
#endif
    }
    else if (file_name.find(".3dm") != std::string::npos) {
      m_pg->initialize(OPENNURBS_GEOMETRY_ENGINE);
      file_imported = m_pg->import_open_nurbs_file(file_name.c_str());
    }


    if (file_imported) {
      const std::string surfa = "geom_surface_tri_";
      const std::string surfb = "geom_surface_quad_";
      const std::string curv = "geom_curve_";

      std::string attribute;
      int id;

      std::vector<int> allSurfs;
      m_pg->get_surfaces(allSurfs);
      for (unsigned i = 0; i < allSurfs.size(); i++) {
        id=allSurfs[i];
        attribute = "#" + surfa + std::to_string(id) + " " + surfb + std::to_string(id);
        geometry_entities.push_back(GeometryHandle(id, SURFACE, attribute));
      }
      std::vector<int> allCurves;
      m_pg->get_curves(allCurves);
      for (unsigned i = 0; i < allCurves.size(); i++) {
        id=allCurves[i];
        attribute = curv + std::to_string(id);
        geometry_entities.push_back(GeometryHandle(id, CURVE, attribute));
      }
      
      if (allSurfs.size() == 0 && allCurves.size() == 0) {
        throw std::runtime_error("Your geometry file had no curves or surfaces. Unusable by GeometryKernelPGEOM");
        return false;
      }
      
      return true;
    }
    
    throw std::runtime_error("NO such file exists");
    
    return false;
  }

  void GeometryKernelPGEOM::snap_to(KernelPoint& point, GeometryHandle geom, /* will snap a given point to a given piece of geometry*/
                                    double *converged_tolerance, double *uvw_computed,
                                    double *uvw_hint, void *extra_hint)
  { //TODO: ADD AN OPTION FOR 2-D
    if (is_curve(geom)) {
      CubitVector old_locations(point[0],point[1],point[2]);
      CubitVector new_locations = m_pg->get_curv_closest_point(geom.m_id, old_locations);
      point[0] = new_locations.x();
      point[1] = new_locations.y();
      point[2] = new_locations.z();
    }
    else if (is_surface(geom)) {
      CubitVector old_locations(point[0], point[1], point[2]);
      CubitVector new_locations = m_pg->get_surf_closest_point(geom.m_id, old_locations);
      
      point[0] = new_locations.x();
      point[1] = new_locations.y();
      point[2] = new_locations.z();
    }
    
  }

  void GeometryKernelPGEOM::normal_at(KernelPoint& point, GeometryHandle geom,
                                      std::vector<double>& normal, void *extra_hint)
  {
    // TODO implement this!
  }

}

#endif // HAVE_CUBIT

