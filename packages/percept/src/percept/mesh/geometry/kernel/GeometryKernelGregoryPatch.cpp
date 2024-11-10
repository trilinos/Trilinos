// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <iostream>
#include <typeinfo>
#include <stdexcept>

#include <percept/mesh/geometry/stk_geom/3D/EvaluateGregoryPatch.hpp>
#include <percept/mesh/geometry/stk_geom/3D/FitGregoryPatches.hpp>
#include <percept/math/Math.hpp>

#include "GeometryKernelGregoryPatch.hpp"
#include <string>

#define DEBUG_GEOM_ON 0

namespace percept {

static bool useLinearOnly = true;

GeometryKernelGregoryPatch::~GeometryKernelGregoryPatch()
{
  delete m_geometryMesh;
}

bool GeometryKernelGregoryPatch::debug_dump_file(const std::string& file_name)
{
  return true;
}


bool GeometryKernelGregoryPatch::read_file
(
 const std::string& file_name,
 std::vector<GeometryHandle>& geometry_entities
 )
{
  m_geometryMesh = new percept::PerceptMesh();
  std::string options = "";
  bool auto_decomp = get_property("auto_decomp") == "true";
  bool exo_large = get_property("exo_large") == "true";
  if (auto_decomp && exo_large)
    options += "large,auto-decomp:yes";
  else if (exo_large)
    options += "large";
  else if (auto_decomp)
    options += "auto-decomp:yes";
  if (auto_decomp || exo_large) m_geometryMesh->set_ioss_read_options(options);
  m_geometryMesh->open(file_name);
  percept::FitGregoryPatches fg(*m_geometryMesh);
  bool doRegister = false;
  fg.register_or_set_fields(doRegister);
  m_geometryMesh->commit();
  if (!m_geometryMesh->m_gregory_control_points_field && !m_geometryMesh->m_gregory_control_points_field_shell)
    {
      throw std::runtime_error("GeometryKernelGregoryPatch:: no Gregory data found");
      return false;
    }
  geometry_entities.resize(0);

  const stk::mesh::PartVector& parts = m_geometryMesh->get_fem_meta_data()->get_parts();
  unsigned nparts = parts.size();
  for (unsigned ipart=0; ipart < nparts; ipart++)
    {
      stk::mesh::Part& part = *parts[ipart];
      bool percept_auto_part = part.attribute<percept::AutoPart>() != 0;
      if (stk::mesh::is_auto_declared_part(part) || (part.name().find("oldElem") != std::string::npos) || percept_auto_part)
        continue;

      if (percept::FitGregoryPatches::is_surface_topology(part.topology()))
        {
          percept::GregoryControlPointsType *field = 0;
          if (part.primary_entity_rank() == m_geometryMesh->side_rank())
            field = m_geometryMesh->m_gregory_control_points_field;
          else
            field = m_geometryMesh->m_gregory_control_points_field_shell;
          if(percept::PerceptMesh::field_is_defined_on_part(field, part))
            {
              geometry_entities.push_back(GeometryHandle(static_cast<int>(part.mesh_meta_data_ordinal()), SURFACE, part.name()));
            }
        }
      else if (part.name() == "edgeseams")
        {
          geometry_entities.push_back(GeometryHandle(static_cast<int>(part.mesh_meta_data_ordinal()), CURVE, "edgeseams"));
        }
    }

  return true;
}

void GeometryKernelGregoryPatch::snap_to
(
 KernelPoint& point,
 GeometryHandle geom,
 double *converged_tolerance,
 double *uvw_computed,
 double *uvw_hint,
 void *hint
 )
{
  unsigned ordinal = static_cast<unsigned>(geom.m_id);
  stk::mesh::Part * part_ptr = &m_geometryMesh->get_fem_meta_data()->get_part(ordinal);
  VERIFY_OP_ON(part_ptr, !=, 0, "bad part");

  VERIFY_OP_ON(hint, !=, 0, "bad hint");
  stk::mesh::Entity node_hint = *static_cast<stk::mesh::Entity *>(hint);

  VERIFY_OP_ON(m_nodeMesh.m_unprojected_coordinates, !=, 0, "must use m_unprojected_coordinates field");
  double *unprojected_coordinates = stk::mesh::field_data(*m_nodeMesh.m_unprojected_coordinates, node_hint);

  bool debug = m_debug || get_debug();
  if (debug)
    {
      std::cout << "P[" << m_nodeMesh.get_rank() << "] GKGP::out1 part= " << part_ptr->name()
                << " node= " << m_nodeMesh.id(node_hint) << " pv: "
                << m_nodeMesh.print_entity_parts_string(node_hint)
                << "\n==> input_point= " << percept::Math::print_3d(point)
                << "\n==>          uc= " << percept::Math::print_3d(unprojected_coordinates) << ", " << unprojected_coordinates[3]
                << std::endl;
    }

  stk::mesh::Selector inMyParts = stk::mesh::selectUnion(m_geometryMeshActiveParts);

  stk::mesh::Entity closest_face = stk::mesh::Entity();
  std::set<stk::mesh::Entity> neighbors, neighbors_local;
  if (!m_geometryMesh->is_valid(closest_face))
    {
      GeometryKernelGregoryPatch::EntitySet& faceSet = m_meshTransfer->meshB()->m_searchMap[node_hint];
      neighbors = faceSet;
    }

  double closest_point[3] = {0,0,0}, closest_uv[2] = {0,0};
  bool allowError = false;
  double point_copy[3];
  bool useUC = true;
  std::string pp = get_property("GKGP:use_unprojected_coordinates");
  useLinearOnly = true;
  if (get_property("GKGP:use_unprojected_coordinates") == "false")
    {
      useLinearOnly = false;
      useUC = false;
    }

  for (int ii=0; ii < m_geometryMesh->get_spatial_dim(); ++ii)
    {
      if (useUC)
        {
          point_copy[ii] = unprojected_coordinates[ii];
        }
      else
        {
          point_copy[ii] = point[ii];
        }
    }

  bool err = findClosestPointInternal(point_copy, neighbors, node_hint, closest_point, closest_uv, allowError, debug);

  // brute force
  if (err)
    {
      err = findClosestPointInternal(point_copy, neighbors, node_hint, closest_point, closest_uv, allowError, true);

      if (1)
        {
          std::cout << "P[" << m_nodeMesh.get_rank() << "] GKGP::out6a WARNING - brute force search... errSearch input point = " << percept::Math::print_3d(point)
                    << " closest_point = " << percept::Math::print_3d(closest_point)
                    << " closest_uv = " << percept::Math::print_2d(closest_uv)
                    << " node_hint= " << m_nodeMesh.id(node_hint)
                    << " neighbors.size= " << neighbors.size()
                    << " err= " << err
                    << std::endl;
        }
      VERIFY_OP_ON(err, == , false, "no err allowed");

      std::cout << "P[" << m_nodeMesh.get_rank() << "] GKGP::out6a WARNING 1 - brute force search... errSearch input point = " << percept::Math::print_3d(point)
                << " node= " << m_nodeMesh.id(node_hint)
                << std::endl;

      std::vector<stk::mesh::Entity> vecFaces, vecShells;

      stk::mesh::get_selected_entities(inMyParts, m_geometryMesh->get_bulk_data()->buckets(m_geometryMesh->side_rank()), vecFaces);
      stk::mesh::get_selected_entities(inMyParts, m_geometryMesh->get_bulk_data()->buckets(m_geometryMesh->element_rank()), vecShells);
      vecFaces.insert(vecFaces.end(), vecShells.begin(), vecShells.end());

      std::set<stk::mesh::Entity> allFaces;
      for (size_t ii = 0; ii < vecFaces.size(); ++ii)
        {
          stk::mesh::Entity face = vecFaces[ii];

          if (percept::FitGregoryPatches::is_surface_topology(m_geometryMesh->bucket(face).topology()))
            {
              stk::mesh::Entity root = m_geometryMesh->rootOfTree(face);
              allFaces.insert(root);
            }
        }

      //  allowError = true;
      err = findClosestPointInternal(point_copy, allFaces, node_hint, closest_point, closest_uv, allowError, debug);

      if (debug)
        {
          std::cout << "P[" << m_nodeMesh.get_rank() << "] GKGP::out6a WARNING - brute force search... errSearch input point = " << percept::Math::print_3d(point)
                    << " closest_point = " << percept::Math::print_3d(closest_point)
                    << " closest_uv = " << percept::Math::print_2d(closest_uv)
                    << " node_hint= " << m_nodeMesh.id(node_hint)
                    << " neighbors.size= " << neighbors.size()
                    << " err= " << err
                    << std::endl;
        }
    }

  if (err)
    {
      static int entered = 0;
      std::cout << "P[" << m_nodeMesh.get_rank() << "] GKGP::out6 errSearch input point = " << percept::Math::print_3d(point)
                << " closest_point = " << percept::Math::print_3d(closest_point)
                << " closest_uv = " << percept::Math::print_2d(closest_uv)
                << " node_hint= " << m_nodeMesh.id(node_hint)
                << " neighbors.size= " << neighbors.size()
                << std::endl;
      m_debug = true;
      if (entered == 0)
        {
          ++entered;
          snap_to(point, geom, converged_tolerance, closest_uv, uvw_hint, hint);
        }
      throw std::runtime_error("error in GeometryKernelGregoryPatch, couldn't project point");
    }
  VERIFY_OP_ON(m_geometryMesh->is_valid(m_found_face), ==, true, "bad found_face");

  if(uvw_computed)
    {
      uvw_computed[0] = closest_uv[0];
      uvw_computed[1] = closest_uv[1];
    }

  percept::Math::copy_3d(point, closest_point);
}

bool GeometryKernelGregoryPatch::findClosestPointInternal(const double *point,
                                                          std::set<stk::mesh::Entity>& neighbors,
                                                          stk::mesh::Entity node_hint,
                                                          double *closest_point,
                                                          double *uvw_computed,
                                                          bool allowError,
                                                          bool debug)
{
  stk::mesh::Entity found_face = stk::mesh::Entity();
  double dmin = std::numeric_limits<double>::max();
  double closest_point_local[3] = {0,0,0}, found_uv_local[2]= {0,0};
  double closest_uv[2] = {0,0};

  bool ldebug = debug;

  bool lmdebug = get_debug();

  size_t nSurfaceElements = 0;
  for (std::set<stk::mesh::Entity>::iterator iter = neighbors.begin(); iter != neighbors.end(); ++iter)
    {
      stk::mesh::Entity neigh = *iter;
      if (!percept::FitGregoryPatches::is_surface_topology(m_geometryMesh->bucket(neigh).topology()))
        continue;
      ++nSurfaceElements;

      stk::mesh::Entity root = m_geometryMesh->rootOfTree(neigh);

      if (debug || ldebug)
        {
          std::cout << "\nP[" << m_nodeMesh.get_rank() << "] GKGP::out7 visiting root= " << m_geometryMesh->id(root) << " node= " << m_nodeMesh.id(node_hint)
                    << "\n==> input_point= " << percept::Math::print_3d(point)
                    << std::endl;
        }

      percept::EvaluateGregoryPatch evgp(*m_geometryMesh, ldebug);

      bool linearOnly = useLinearOnly;
      bool err = evgp.findClosestPoint(point, root, closest_point_local, found_uv_local, linearOnly);
      if (ldebug)
        {
          if (1) std::cout << "\nP[" << m_nodeMesh.get_rank() << "] GKGP::out8 visiting face= " << m_geometryMesh->id(root) << " found_uv_local= " << percept::Math::print_2d(found_uv_local)
                           << " closest_point_local= " << percept::Math::print_3d(closest_point_local)
                           << "==> input_point= " << percept::Math::print_3d(point)
                           << " linearOnly= " << linearOnly << " err= " << err
                           << std::endl;
        }
      double d = percept::Math::distance_squared_3d(point, closest_point_local);

      if (lmdebug) std::cout << "\nP[" << m_nodeMesh.get_rank() << "] GKGP::out8a err= " << err << " face= " << m_geometryMesh->id(root) << " closest_point = " << percept::Math::print_3d(closest_point_local)
                             << "==> input_point= " << percept::Math::print_3d(point)
                            << " found_face= " << found_face << " d= " << d << " dmin= " << dmin
                            << std::endl;

      if (debug)
        {
          double cp_loc[3];
          percept::EvaluateGregoryPatch evgp1(*m_geometryMesh, false);

          evgp1.evaluateGregoryPatch(found_uv_local, root, cp_loc);
          double d0 = percept::Math::distance_squared_3d(cp_loc, closest_point_local);

          std::cout << "\n\nP[" << m_geometryMesh->get_rank()
                    << "] GKGP::out9 point = " << point[0] << ", " << point[1] << ", " << point[2]
                    << " cpoint = " << closest_point_local[0] << ", " << closest_point_local[1] << ", " << closest_point_local[2]
                    << "\nroot= " << m_geometryMesh->id(root) //<< " parts= " << m_geometryMesh->print_entity_parts_string(root)
                    << " rootPr= " << m_geometryMesh->print_entity_compact(root)
                    << " neigh= " << m_geometryMesh->id(neigh) //<< " parts= " << m_geometryMesh->print_entity_parts_string(neigh)
                    << " d= " << d << " dmin= " << dmin << " d0= " << d0
                    << " iv= " << m_geometryMesh->is_valid(found_face)
                    << " err= " << err
                    << " errMsg= " << evgp.error_message()
                    << "\n\n"
                    << std::endl;
        }
      if (err && !allowError)
        {
          if (ldebug)
            std::cout << "P[" << m_geometryMesh->get_rank() << "] GKGP::out10 err in findClosestPoint, face= " << m_geometryMesh->id(root) << " point = " << point[0] << ", "
                      << point[1] << ", " << point[2] << std::endl;
          continue;
        }
      if (!m_geometryMesh->is_valid(found_face) || d < dmin)
        {
          dmin = d;
          found_face = root;
          percept::Math::copy_3d(closest_point, closest_point_local);
          percept::Math::copy_2d(closest_uv, found_uv_local);
          if (lmdebug) std::cout << "\nP[" << m_nodeMesh.get_rank() << "] GKGP::out11 no err, face= " << m_geometryMesh->id(root) << " closest_point = " << percept::Math::print_3d(closest_point) << std::endl;
        }
    }

  if (!nSurfaceElements)
    throw std::runtime_error("couldn't find any surface elements");

  m_found_face = found_face;
  if (!m_geometryMesh->is_valid(found_face))
    {
      return true;
    }

  if (useLinearOnly)
    {
      percept::EvaluateGregoryPatch evgp1(*m_geometryMesh, false);
      evgp1.evaluateGregoryPatch(closest_uv, found_face, closest_point);
    }

  if (1)
    {
      m_nodeToFoundFaceMap[node_hint] = found_face;
    }

  if (lmdebug) std::cout << "\nP[" << m_nodeMesh.get_rank() << "] GKGP:: final found_face= " << m_geometryMesh->id(found_face)
                         << " face= " << m_geometryMesh->print_entity_compact(found_face)
                         << " node_hint= " << m_nodeMesh.id(node_hint)
                         << " closest_uv= " << percept::Math::print_2d(closest_uv)
                         << " closest_point= " << percept::Math::print_3d(closest_point)
                         << " dmin= " << dmin
                         << " d_cp_p= " << percept::Math::distance_squared_3d(point, closest_point)
                         << "\n node= " << m_nodeMesh.print_entity_compact(node_hint)
                         << "\n nodeUC= " << m_nodeMesh.print_entity_compact(node_hint, m_nodeMesh.m_unprojected_coordinates)
                         << std::endl;

  if(uvw_computed)
    {
      uvw_computed[0] = closest_uv[0];
      uvw_computed[1] = closest_uv[1];
    }
  return false;
}

static  void
faceNormal(percept::PerceptMesh& m_eMesh, stk::mesh::Entity face, double * normal)
{
  const percept::MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), face, stk::topology::NODE_RANK );

  double *n0 = static_cast<double*>(stk::mesh::field_data( *m_eMesh.get_coordinates_field() , face_nodes[0].entity() ));
  double *n1 = static_cast<double*>(stk::mesh::field_data( *m_eMesh.get_coordinates_field() , face_nodes[1].entity() ));
  double *n2 = static_cast<double*>(stk::mesh::field_data( *m_eMesh.get_coordinates_field() , face_nodes[2].entity() ));
  // works for triangles, quads, and nonplanar quads
  double *n3 = face_nodes.size() == 3 ? n0 : static_cast<double*>(stk::mesh::field_data( *m_eMesh.get_coordinates_field() , face_nodes[3].entity() ));
  double a[3], b[3];
  for (unsigned jj=0; jj < 3; ++jj)
    {
      a[jj] = n2[jj] - n0[jj];
      b[jj] = n3[jj] - n1[jj];
    }

  percept::Math::cross_3d(a, b, normal);
  percept::Math::normalize_3d(normal);
}

void GeometryKernelGregoryPatch::normal_at(KernelPoint& point, GeometryHandle geom, std::vector<double>& normal, void *hint)
{
  VERIFY_OP_ON(hint, !=, 0, "bad hint");

  double linearNormal[3] = {0,0,0};
  if (1)
    {
      stk::mesh::Entity node_hint = *static_cast<stk::mesh::Entity *>(hint);
      stk::mesh::Entity face = m_nodeToFoundFaceMap[node_hint];
      VERIFY_OP_ON(m_geometryMesh->is_valid(face), ==, true, "bad face");
      faceNormal(*m_geometryMesh, face, &linearNormal[0]);
    }
  double copy_point[3];
  percept::Math::copy_3d(copy_point, point);
  double uv[2] = {0,0};
  double tol = 1.e-6;

  double *cp = &copy_point[0];
  double *cuv = &uv[0];
  snap_to(cp, geom, &tol, cuv, 0, hint);
  bool debug = false;
  percept::EvaluateGregoryPatch evgp(*m_geometryMesh, debug);
  evgp.normalGregoryPatch(uv, m_found_face, &normal[0]);
  double dn = percept::Math::dot_3d(&normal[0], &linearNormal[0]);
  if (dn <= 0.0)
    {
      percept::Math::copy_3d(&normal[0], &linearNormal[0]);
    }
}

void GeometryKernelGregoryPatch::
pre_process(percept::PerceptMesh* eMeshWithNodes , const stk::mesh::PartVector& nodeParts)
{
  VERIFY_OP_ON(eMeshWithNodes, ==, &m_nodeMesh, "bad mesh");
  if (m_nodeMesh.get_rank() == 0)
    std::cout << "GeometryKernelGregoryPatch::pre_process start..." << std::endl;
  m_nodeMeshActiveParts = nodeParts;
  m_geometryMeshActiveParts.resize(0);
  unsigned nparts = nodeParts.size();
  for (unsigned ipart=0; ipart < nparts; ipart++)
    {
      stk::mesh::Part * part = m_geometryMesh->get_fem_meta_data()->get_part(nodeParts[ipart]->name());
      m_geometryMeshActiveParts.push_back(part);
      if (m_debug)
        std::cout << "P[" << m_nodeMesh.get_rank() << "] GKGP::pre_process nodeParts[" << ipart << "] = "
                  << nodeParts[ipart]->name() << " part= " << part << std::endl;
    }

  std::vector<stk::mesh::FieldBase *> fromFields;
  if (m_geometryMesh->m_gregory_control_points_field_shell)
    fromFields.push_back(m_geometryMesh->m_gregory_control_points_field_shell);
  if (m_geometryMesh->m_gregory_control_points_field)
    fromFields.push_back(m_geometryMesh->m_gregory_control_points_field);
  std::vector<stk::mesh::FieldBase *> toFields;
  toFields.push_back(m_nodeMesh.get_coordinates_field());

  if (m_nodeMesh.get_rank() == 0)
    std::cout << "GeometryKernelGregoryPatch::pre_process buildGPSTKMeshTransfer..." << std::endl;

  m_meshTransfer =
    percept::buildGPSTKMeshTransfer(*m_geometryMesh, m_geometryMeshActiveParts, fromFields,
                                    m_nodeMesh, nodeParts, toFields, "name");

  if (m_nodeMesh.get_rank() == 0)
    std::cout << "GeometryKernelGregoryPatch::pre_process buildGPSTKMeshTransfer...done" << std::endl;

  initializeSTKMeshTransfer(&*m_meshTransfer);

  if (m_nodeMesh.get_rank() == 0)
    std::cout << "GeometryKernelGregoryPatch::pre_process initializeSTKMeshTransfer...done" << std::endl;

  m_meshTransfer->apply();

  if (m_nodeMesh.get_rank() == 0)
    std::cout << "GeometryKernelGregoryPatch::pre_process done" << std::endl;

}

}
