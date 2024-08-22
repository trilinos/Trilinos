// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "GregoryPatch.hpp"
#include "GregoryPatchData.hpp"
#include "FitGregoryPatches.hpp"
#include <percept/Util.hpp>
#include <stk_util/parallel/CommSparse.hpp>


#define APRINTLN(a) do { if (debug_print) std::cout << #a << " = " << a << std::endl; } while(0)
#define APRINTLN2(a,b) do { if (debug_print) std::cout << #a << " = " << a << " " << #b << " = " << b << std::endl; } while(0)
#define ID(n) (m_eMesh.id(n))
#define KEY(n) (m_eMesh.entity_key(n))
#define OWNED(n) (m_eMesh.owned(n))
#define SHARED(n) (m_eMesh.shared(n))
#define OWNER(n) (m_eMesh.owner_rank(n))
#define TOP(n) (m_eMesh.topology(n))

namespace percept {

  const double FitGregoryPatches::m_globalAngleCriterionDefault = 135.0;

  static bool checkNAN(MDArray& c)
  {
    double* dd = c.data();
    for (unsigned ii=0; ii < c.size(); ++ii)
      {
        if (dd[ii] != dd[ii])
          return true;
      }
    return false;
  }
  // static
  bool FitGregoryPatches::is_tri_topology(stk::topology topo)
  {
    switch (topo.value())
      {
      case stk::topology::SHELL_TRI_3:
      case stk::topology::TRI_3:
        return true;
      default:
        return false;
      }
    return false;
  }
  bool FitGregoryPatches::is_quad_topology(stk::topology topo)
  {
    switch (topo.value())
      {
      case stk::topology::SHELL_QUAD_4:
      case stk::topology::QUAD_4:
        return true;
      default:
        return false;
      }
    return false;
  }
  bool FitGregoryPatches::is_surface_topology(stk::topology topo)
  {
    return is_tri_topology(topo) || is_quad_topology(topo);
  }

  // template<class Field>
  // static Field *get_or_declare_field
  bool FitGregoryPatches::in_surface_sets(const std::string& partToTest)
  {
    if (m_surfaceSets.size() == 0)
      return true;

    for (SurfaceSets::iterator it = m_surfaceSets.begin(); it != m_surfaceSets.end(); ++it)
      {
        //const std::string& surfaceSetName = it->first;
        std::vector<std::string>& sset = it->second;

        for (unsigned iSurface = 0; iSurface < sset.size(); ++iSurface)
          {
            std::string part_name = sset[iSurface];
            if (part_name == partToTest)
              return true;
          }
      }
    return false;
  }

  void FitGregoryPatches::register_or_set_fields(bool doRegister)
  {
    if (!m_eMesh.m_gregory_control_points_field_set)
      {
        m_eMesh.m_gregory_control_points_field_set = true;
        // we register more than needed for triangles
        // we lay out the control points as {x0,x1,...xn, y0, ... yn, z0,z1, ... zn}
        // Cp(i,j): i=0,n, j=0,3: Cp(i + n*j)
        const int numControlPoints = 3*MaxControlPoints();

        m_eMesh.m_gregory_control_points_field = &(m_eMesh.get_fem_meta_data()->declare_field<GregoryControlPointsType::value_type>(m_eMesh.side_rank(), "gregory_control_points"));
        m_eMesh.m_gregory_control_points_field_shell = &(m_eMesh.get_fem_meta_data()->declare_field<GregoryControlPointsType::value_type>(m_eMesh.element_rank(), "gregory_control_points_shell"));
        m_eMesh.m_node_normals = &(m_eMesh.get_fem_meta_data()->declare_field<NormalsFieldType::value_type>(stk::topology::NODE_RANK, "node_normals"));

        if (doRegister)
          {
            stk::io::set_field_role(*m_eMesh.m_gregory_control_points_field, Ioss::Field::TRANSIENT);
            stk::io::set_field_role(*m_eMesh.m_gregory_control_points_field_shell, Ioss::Field::TRANSIENT);
            stk::io::set_field_role(*m_eMesh.m_node_normals, Ioss::Field::TRANSIENT);

            const stk::mesh::PartVector& pv = m_eMesh.get_fem_meta_data()->get_parts();
            for (unsigned ii=0; ii < pv.size(); ++ii)
              {
                stk::mesh::Part& part = *pv[ii];
                bool percept_auto_part = part.attribute<AutoPart>() != 0;

                if (stk::mesh::is_auto_declared_part(part) || percept_auto_part)
                  continue;

                // 3 slots for "left" of an edge, 3 for "right"
                if (is_surface_topology(part.topology()))
                  {
                    stk::mesh::put_field_on_mesh(*m_eMesh.m_node_normals, part, 3, nullptr);
                    stk::io::set_field_output_type(*m_eMesh.m_node_normals, stk::io::FieldOutputType::VECTOR_3D);
                  }

                if (!in_surface_sets(part.name()))
                  continue;

                if (pv[ii]->primary_entity_rank() == m_eMesh.side_rank())
                  {
                    stk::mesh::put_field_on_mesh(*m_eMesh.m_gregory_control_points_field, part, numControlPoints, nullptr);
                  }
                else if (pv[ii]->primary_entity_rank() == m_eMesh.element_rank()
                         && (pv[ii]->topology() == stk::topology::SHELL_QUAD_4 || pv[ii]->topology() == stk::topology::SHELL_TRI_3))
                  {
                    stk::mesh::put_field_on_mesh(*m_eMesh.m_gregory_control_points_field_shell, part, numControlPoints, nullptr);
                  }
              }
          }
        if (1)
          {
            stk::mesh::Part * part = m_eMesh.get_fem_meta_data()->get_part("edgeseams");
            if (!part)
              {
                part = &m_eMesh.get_fem_meta_data()->declare_part_with_topology("edgeseams", stk::topology::BEAM_2);
                stk::io::put_io_part_attribute(*part);
                if (m_debug) std::cout << "tmp srk edgeseams part.topology= " << part->topology() << std::endl;
              }
            m_edgeSeamsPart = part;
          }
      }
  }

  // computes control points
  /**
   * Algorithm:
   * 1. define sets of surfaces to be treated together
   *
   * for each set:
   * 2. find normals at all face nodes and average to get node normals
   * 3. visit faces and do cubic fit to edges
   *    3a. compare face normal at end nodes - if differs by more than specified
   *           angle, use face normal, else node normal
   *    3b. inject back into Cp
   * 4. visit faces and edge neighbors
   *    4a. form ribbons for left and right p,q,r
   *    4b. fit ribbons, inject into Cp array
   *
   * 5. eval faces with point cloud
   */

  class MyString {
  public:
    std::string m_string;
    MyString(const std::string& str) : m_string(str) {}
  };

  void FitGregoryPatches::
  getCurrentParts(std::vector<stk::mesh::PartVector>& currentParts)
  {
    if (m_surfaceSets.size() == 0)
      {
        currentParts.resize(1);
        const stk::mesh::PartVector& pv = m_eMesh.get_fem_meta_data()->get_parts();
        for (unsigned ii=0; ii < pv.size(); ++ii)
          {
            stk::mesh::Part& part = *pv[ii];
            bool percept_auto_part = part.attribute<AutoPart>() != 0;

            if (stk::mesh::is_auto_declared_part(part) || percept_auto_part)
              continue;

            if (
                (part.primary_entity_rank() == m_eMesh.side_rank()
                 && (part.topology() == stk::topology::QUAD_4 || part.topology() == stk::topology::TRI_3))
                ||
                (part.primary_entity_rank() == m_eMesh.element_rank()
                 && (part.topology() == stk::topology::SHELL_QUAD_4 || part.topology() == stk::topology::SHELL_TRI_3))
                )
              {
                const MyString *existingAttribute = part.attribute<MyString>();
                if (existingAttribute)
                  {
                    m_eMesh.get_fem_meta_data()->remove_attribute(part, existingAttribute);
                    delete existingAttribute;
                  }
                MyString * myStringAttribute = new MyString("all");
                m_eMesh.get_fem_meta_data()->declare_attribute_with_delete(part, myStringAttribute);
                currentParts[0].push_back(pv[ii]);
              }
          }
        if (currentParts[0].size() == 0)
          {
            throw std::runtime_error("No surface parts found");
          }
      }
    else
      {
        currentParts.resize(m_surfaceSets.size());
        int iSurfaceSet = 0;
        for (SurfaceSets::iterator it = m_surfaceSets.begin(); it != m_surfaceSets.end(); ++it)
          {
            const std::string& surfaceSetName = it->first;
            std::vector<std::string>& sset = it->second;

            currentParts[iSurfaceSet].resize(0);
            for (unsigned iSurface = 0; iSurface < sset.size(); ++iSurface)
              {
                std::string part_name = sset[iSurface];
                stk::mesh::Part *found_part = m_eMesh.get_fem_meta_data()->get_part(part_name);
                if (!found_part)
                  {
                    std::cout << "Parts = " << m_eMesh.print_part_vector_string(m_eMesh.get_fem_meta_data()->get_parts(), "\n") << std::endl;
                  }
                VERIFY_OP_ON(found_part, !=, 0, "bad part in FitGregoryPatches: "+part_name);

                const MyString *existingAttribute = found_part->attribute<MyString>();
                if (existingAttribute)
                  {
                    m_eMesh.get_fem_meta_data()->remove_attribute(*found_part, existingAttribute);
                    delete existingAttribute;
                  }
                MyString * myStringAttribute = new MyString(surfaceSetName);
                m_eMesh.get_fem_meta_data()->declare_attribute_with_delete(*found_part, myStringAttribute);

                currentParts[iSurfaceSet].push_back(found_part);
              }
            ++iSurfaceSet;
          }
      }
  }

  void FitGregoryPatches::computeControlPoints(bool doGetNormals, bool createEdgeSeamsPart)
  {
    std::vector<stk::mesh::PartVector> currentParts;
    getCurrentParts(currentParts);

    bool doSeams = true;
    if (m_eMesh.getProperty("FitGregoryPatches::noEdgeFitting") == "true")
      doSeams = false;

    for (unsigned iSurfaceSet = 0; iSurfaceSet < currentParts.size(); ++iSurfaceSet)
      {
        if (doSeams)
          {
            findSeams(currentParts[iSurfaceSet], createEdgeSeamsPart);
            processSeams(currentParts[iSurfaceSet], createEdgeSeamsPart);
          }

        if (doGetNormals)
          {
            getNormals(&currentParts[iSurfaceSet]);
          }

        fitCubics(currentParts[iSurfaceSet]);

        fitRibbons(currentParts[iSurfaceSet]);
      }

    if (m_edgeSeamsPart) addEdgesToMesh(m_edgeSetAll);
  }

  int FitGregoryPatches::
  orient(stk::mesh::Entity face, stk::mesh::Entity node, stk::mesh::Selector& sel, bool debug)
  {
    const MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), face, stk::topology::NODE_RANK );
    bool found = false;
    for (unsigned ii=0; ii < face_nodes.size(); ++ii)
      {
        if (face_nodes[ii].entity() == node)
          {
            found = true;
            break;
          }
      }
    if (!found)
      {
        if (debug) std::cout << "FGP::orient: case A " << std::endl;
        return -1;
      }
    EntitySet visited;
    int val = orientRecurse(face, node, visited, sel, debug);

    if (debug)
      {
        int ii=0;
        std::ostringstream str;
        str << "val = " << val << " visited, face= " << m_eMesh.print_entity_compact(face) << "\n " << m_eMesh.print_entity_compact(node);
        for (EntitySet::iterator it = visited.begin(); it != visited.end(); ++it, ++ii)
          {
            str << "\n" << m_eMesh.print_entity_compact(*it);
          }
        str << "\nneighs= ";

        EntitySet neighbors;
        getEdgeNeighborsSharingNode(face, sel, node, neighbors);
        for (EntitySet::iterator it = neighbors.begin(); it != neighbors.end(); ++it, ++ii)
          {
            str << "\n" << m_eMesh.print_entity_compact(*it);
          }
        str << "Face Parts = " << m_eMesh.print_part_vector_string(m_eMesh.bucket(face).supersets()) << std::endl;
        str << "Node Parts = " << m_eMesh.print_part_vector_string(m_eMesh.bucket(node).supersets()) << std::endl;
        std::cout << str.str() << std::endl;
        if (1)
          {
            EntitySet faceSet;
            faceSet.insert(face);

            m_eMesh.dump_vtk("face.vtk", true, &faceSet);
            m_eMesh.dump_vtk("visited.vtk", true, &visited);
            m_eMesh.dump_vtk("neigh.vtk", true, &neighbors);
          }
      }
    return val;
  }

  int FitGregoryPatches::
  orientRecurse(stk::mesh::Entity face, stk::mesh::Entity node, EntitySet& visited, stk::mesh::Selector& sel, bool debug)
  {
    int val = orientLowLevel(face, node, debug);
    if (val != -1)
      return val;

    EntitySet neighbors;
    getEdgeNeighborsSharingNode(face, sel, node, neighbors);

    for (EntitySet::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
      {
        stk::mesh::Entity neigh = *it;
        if (visited.find(neigh) != visited.end())
          continue;
        visited.insert(neigh);
        int valR = orientRecurse(neigh, node, visited, sel, debug);
        if (valR != -1)
          return valR;
      }
    return -1;
  }

  void FitGregoryPatches::
  getEdgeNeighborsSharingNode(stk::mesh::Entity face, stk::mesh::Selector& sel, stk::mesh::Entity node, EntitySet& neighbors_in)
  {
    EdgeVector& edges = m_nodeToEdgeMap[node];

    neighbors_in.clear();
    EntitySet neighbors, shell_neighbors;
    m_eMesh.get_node_neighbors(face, neighbors, sel, m_eMesh.side_rank());
    m_eMesh.get_node_neighbors(face, shell_neighbors, sel, m_eMesh.element_rank());
    neighbors.insert(shell_neighbors.begin(), shell_neighbors.end());
    for (EntitySet::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
      {
        stk::mesh::Entity neigh = *it;
        if (neigh == face)
          continue;

        if (!is_surface_topology(m_eMesh.topology(neigh)))
          continue;

        int edge_0=0, edge_1=0;
        bool isEdgeN = m_eMesh.is_edge_neighbor(face, neigh, &edge_0, &edge_1);
        if (isEdgeN)
          {
            const MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), face, stk::topology::NODE_RANK );
            stk::mesh::Entity n0 = face_nodes[edge_0].entity();
            stk::mesh::Entity n1 = face_nodes[(edge_0 + 1) % face_nodes.size()].entity();

            // avoid crossing into another quadrant
            bool crossing = false;
            for (unsigned ie = 0; ie < edges.size(); ++ie)
              {
                if ((edges[ie].first == n0 && edges[ie].second == n1)
                    || (edges[ie].first == n1 && edges[ie].second == n0))
                  {
                    crossing = true;
                    break;
                  }
              }

            if (!crossing && (n0 == node || n1 == node))
              {
                neighbors_in.insert(neigh);
              }
          }
      }
  }

  int FitGregoryPatches::
  orientLowLevel(stk::mesh::Entity face, stk::mesh::Entity ndIn, bool debug)
  {
    EdgeVector edges = m_nodeToEdgeMap[ndIn];
    if (edges.size() < 2)
      {
        if (debug) std::cout << "FGP::orientLowLevel: case 0" << std::endl;
        return -1;
      }

    if (debug) std::cout << "FGP::orientLowLevel: face= " << ID(face) << " node= " << ID(ndIn) << " edges.size= " << edges.size() << std::endl;

    const MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), face, stk::topology::NODE_RANK );
    bool found = false;
    for (unsigned ii=0; ii < face_nodes.size(); ++ii)
      {
        if (face_nodes[ii].entity() == ndIn)
          {
            found = true;
            break;
          }
      }
    if (!found)
      {
        if (debug) std::cout << "FGP::orientLowLevel: case 1" << std::endl;
        return -1;
      }

    // find the base edge of the "quadrant" the face/node pair is in
    for (unsigned iedge=0; iedge < edges.size(); ++iedge)
      {
        stk::mesh::Entity nodes[2];
        nodes[0] = ndIn;
        nodes[1] = edges[iedge].first == ndIn ? edges[iedge].second : edges[iedge].first;

        if (debug)
          {
            std::cout << "FGP::orientLowLevel: nodes= " << ID(nodes[0]) << " " << ID(nodes[1]) << std::endl;
          }
        VERIFY_OP_ON(nodes[0], !=, nodes[1], "bad nodes");

        for (unsigned ii=0; ii < face_nodes.size(); ++ii)
          {
            stk::mesh::Entity FN0 = face_nodes[ii].entity();
            stk::mesh::Entity FN1 = face_nodes[(ii + 1) % face_nodes.size()].entity();

            if (debug) std::cout << "FGP::orientLowLevel: FN0= " << ID(FN0) << " FN1= " << ID(FN1) << std::endl;

            /**
             *
             * Example with 3 edges at a corner - should return 0 since
             *   edge e0 is touched "positively" by FN0,FN1
             *
             *                        o n3
             *                       /
             *                      e2
             *           n1        /
             *            o---e0--o---e1--- n2
             *                   n0 (corner)
             *
             *  Note: nodes[0] == n0, nodes[1] = n1
             *
             * nodes[1] ==
             *   n1 == FN1 o-----o ndIn (from arg list) == FN0 == n0 == nodes[0]
             *             |     |
             *             |face |
             *             o-----o
             *
             *
             * Example 2: should return 0 also (but only through orientRecurse)
             *     - the base edge of that quadrant is e0 - we return -1
             *         to keep recursion going until it finds the face
             *         positively touching e0
             *
             *      nodes[0] == ndIn ==
             *          n0 == FN1 o-----o FN0 == n2 == nodes[1]
             *                    |     |
             *                    |face |
             *                    o-----o
             *
             * Example 3: face is above e0, and oriented negatively, return 2
             *    through orientRecurse (orientLowLevel returns -1)
             *
             *             o-----o
             *             |     |
             * nodes[1] == |face |
             *   n1 == FN0 o-----o ndIn (from arg list) == FN1 == n0 == nodes[0]
             */

            if (nodes[0] == FN0 && nodes[1] == FN1)
              {
                if (debug) std::cout << "FGP::orientLowLevel: pos case 2, return= " << iedge << std::endl;
                return iedge;
              }

          }
      }
    if (debug) std::cout << "FGP::orientLowLevel: case 3" << std::endl;
    return -1;
  }

  void FitGregoryPatches::
  faceNormal(stk::mesh::Entity face, double normal[3])
  {
    const MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), face, stk::topology::NODE_RANK );

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

    Math::cross_3d(a, b, normal);
    Math::normalize_3d(normal);
  }

  void FitGregoryPatches::
  getNormals(stk::mesh::PartVector* parts)
  {
    stk::mesh::Selector sel = m_eMesh.get_fem_meta_data()->universal_part();
    if (parts) sel = sel & stk::mesh::selectUnion(*parts);

    m_nodeNormalsMap.clear();
    m_nodeNormalsSetMap.clear();

    std::vector<stk::mesh::Entity> vecNodes, vecFaces, vecShells;
    stk::mesh::get_selected_entities(sel , m_eMesh.get_bulk_data()->buckets(m_eMesh.node_rank()), vecNodes);
    for (unsigned ii=0; ii < vecNodes.size(); ++ii)
      {
        stk::mesh::Entity node = vecNodes[ii];
        double *normals_data = stk::mesh::field_data( *m_eMesh.m_node_normals , node );
        for (unsigned jj=0; jj < 3; ++jj)
          {
            normals_data[jj] = 0.0;
          }
      }

    // sum to nodes
    stk::mesh::get_selected_entities(sel , m_eMesh.get_bulk_data()->buckets(m_eMesh.side_rank()), vecFaces);
    stk::mesh::get_selected_entities(sel , m_eMesh.get_bulk_data()->buckets(m_eMesh.element_rank()), vecShells);
    vecFaces.insert(vecFaces.end(), vecShells.begin(), vecShells.end());
    if (m_debug) std::cout << "FitGregoryPatches::getNormals vecNodes= " << vecNodes.size() << " vecFaces= " << vecFaces.size()  << std::endl;
    for (unsigned ii=0; ii < vecFaces.size(); ++ii)
      {
        stk::mesh::Entity face = vecFaces[ii];
        const MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), face, stk::topology::NODE_RANK );
        double normal[3];
        faceNormal(face, normal);
        for (unsigned kn = 0; kn < face_nodes.size(); ++kn)
          {
            stk::mesh::Entity node = face_nodes[kn].entity();
            if (m_eMesh.aura(node))
              continue;
            double *normals_data = stk::mesh::field_data( *m_eMesh.m_node_normals , node );

            int orientation = orient(face, node, sel);

            PointVector& normals = m_nodeNormalsMap[node];
            std::vector<int>& normalsSet = m_nodeNormalsSetMap[node];

            if (orientation >= 0 and (orientation + 1 > int(normals.size())))
              {
                normals.resize(orientation + 1);
                normalsSet.resize(orientation + 1);
              }
            for (unsigned jj=0; jj < 3; ++jj)
              {
                normals_data[jj] += normal[jj];  // globally averaged normals

                if (orientation >= 0)
                  {
                    normals[orientation][jj] += normal[jj];
                    normalsSet[orientation] = 1;
                  }
              }
          }
      }

    // normalize
    for (unsigned ii=0; ii < vecNodes.size(); ++ii)
      {
        stk::mesh::Entity node = vecNodes[ii];
        if (m_eMesh.aura(node))
          continue;
        double *normals_data = stk::mesh::field_data( *m_eMesh.m_node_normals , node );
        if (Math::norm_3d(normals_data) == 0.0)
          {
            normals_data[0] = 1.0; // arbitrary
            normals_data[1] = 0.0; // arbitrary
            normals_data[2] = 0.0; // arbitrary
          }
        Math::normalize_3d(normals_data);
        if (m_debug) std::cout << "P[" << m_eMesh.get_rank() << "] node= " << m_eMesh.id(node)
                               << " norm= " << Math::norm_3d(normals_data)
                               << std::endl;
        PointVector& normals = m_nodeNormalsMap[node];
        std::vector<int>& normalsSet = m_nodeNormalsSetMap[node];
        for (unsigned jj=0; jj < normals.size(); ++jj)
          {
            if (normalsSet[jj] && Math::norm_3d(normals[jj].data()) > 0.0)
              {
                Math::normalize_3d(normals[jj].data());
              }
            else
              {
                normals[jj][0] = 1.0; // arbitrary
                normals[jj][1] = 0.0;
                normals[jj][2] = 0.0;
              }
          }
      }

    std::vector<const stk::mesh::FieldBase *> fields;
    fields.push_back(m_eMesh.m_node_normals);
    stk::mesh::copy_owned_to_shared  (  *m_eMesh.get_bulk_data(), fields);
  }

  struct Comp
  {
    PerceptMesh& m_eMesh;
    Comp(PerceptMesh& eMesh) : m_eMesh(eMesh) {}
    bool operator()(stk::mesh::Entity& n0, stk::mesh::Entity& n1) const
    {
      return m_eMesh.id(n0)  < m_eMesh.id(n1);
    }
  };

  FitGregoryPatches::Edge FitGregoryPatches::
  create_edge(stk::mesh::Entity face, unsigned edge)
  {
    Comp comp(m_eMesh);
    const MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), face, stk::topology::NODE_RANK );
    stk::mesh::Entity n0 = face_nodes[edge].entity();
    stk::mesh::Entity n1 = face_nodes[(edge+1) % face_nodes.size()].entity();
    return create_edge<stk::mesh::Entity>(n0, n1, comp);
  }

  void FitGregoryPatches::
  findSeams(stk::mesh::PartVector& parts, bool createEdgeSeamsPart)
  {
    stk::mesh::Selector sel = stk::mesh::selectUnion(parts);
    VERIFY_OP_ON(parts.size(), >, 0, "parts size = 0");
    const MyString *surfaceSetNameMS = parts[0]->attribute<MyString>();
    VERIFY_OP_ON(surfaceSetNameMS, !=, 0, "bad surface set name attribute");
    std::string surfaceSetName = surfaceSetNameMS->m_string;
    double angleCriterion = m_globalAngleCriterion;
    if (m_angleMap.find(surfaceSetName) != m_angleMap.end())
      {
        angleCriterion = m_angleMap[surfaceSetName];
      }

    std::vector<stk::mesh::Entity> vecFaces, vecShells;
    stk::mesh::get_selected_entities(sel , m_eMesh.get_bulk_data()->buckets(m_eMesh.side_rank()), vecFaces);
    stk::mesh::get_selected_entities(sel , m_eMesh.get_bulk_data()->buckets(m_eMesh.element_rank()), vecShells);
    vecFaces.insert(vecFaces.end(), vecShells.begin(), vecShells.end());

    m_edgeSeamsMap.clear();
    m_edgeSet.clear();

    // first, find seams due to edges of surfaces
    for (unsigned ii=0; ii < vecFaces.size(); ++ii)
      {
        stk::mesh::Entity face = vecFaces[ii];
        const MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), face, stk::topology::NODE_RANK );
        std::vector<int>& emap = m_edgeSeamsMap[face];
        if (!emap.size())
          {
            emap.resize(face_nodes.size());
            emap.assign(face_nodes.size(), 0);
          }
        for (unsigned edge=0; edge < face_nodes.size(); ++edge)
          {
            emap[edge] = 1;
          }
      }

    for (unsigned ii=0; ii < vecFaces.size(); ++ii)
      {
        stk::mesh::Entity face = vecFaces[ii];
        EntitySet neighbors, shell_neighbors;
        m_eMesh.get_node_neighbors(face, neighbors, sel, m_eMesh.side_rank());
        m_eMesh.get_node_neighbors(face, shell_neighbors, sel, m_eMesh.element_rank());
        neighbors.insert(shell_neighbors.begin(), shell_neighbors.end());
        for (EntitySet::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
          {
            stk::mesh::Entity neigh = *it;
            if (neigh == face)
              continue;

            if (!is_surface_topology(m_eMesh.topology(neigh)))
              continue;

            int edge_0=0, edge_1=0;
            bool isEdgeN = m_eMesh.is_edge_neighbor(face, neigh, &edge_0, &edge_1);
            if (isEdgeN)
              {
                m_edgeSeamsMap[face][edge_0] = 0;
                m_edgeSeamsMap[neigh][edge_1] = 0;
              }
          }

        if (!m_eMesh.owned(face))
          {
            const MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), face, stk::topology::NODE_RANK );
            for (unsigned iedge=0; iedge < face_nodes.size(); ++iedge)
              {
                bool found = false;
                for (EntitySet::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
                  {
                    stk::mesh::Entity neigh = *it;
                    if (neigh == face)
                      continue;

                    if (!is_surface_topology(m_eMesh.topology(neigh)))
                      continue;
                    int edge_0=0, edge_1=0;
                    bool isEdgeN = m_eMesh.is_edge_neighbor(face, neigh, &edge_0, &edge_1);
                    if (isEdgeN && edge_0 == int(iedge))
                      {
                        found = true;
                        break;
                      }
                  }
                if (!found)
                  {
                    m_edgeSeamsMap[face][iedge] = 0;
                  }
              }
          }
      }

    bool debug1 = false;
    typedef std::vector<FitGregoryPatches::Edge> VecEdge;
    VecEdge vecEdge;

    if (1)
      {
        for (unsigned ii=0; ii < vecFaces.size(); ++ii)
          {
            stk::mesh::Entity face = vecFaces[ii];
            std::vector<int>& emap = m_edgeSeamsMap[face];
            const MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), face, stk::topology::NODE_RANK );

            for (unsigned iedge=0; iedge < emap.size(); ++iedge)
              {
                if (emap[iedge])
                  {
                    Edge edge = create_edge(face, iedge);
                    if (debug1) vecEdge.push_back(edge);
                    if (m_edgeSet.find(edge) == m_edgeSet.end())
                      {
                        m_edgeSet.insert(edge);
                      }
                  }
              }
          }
        if (debug1)
          {
            Util::makeUnique(vecEdge);
            std::ostringstream ostr;
            for (unsigned ii=0; ii < vecEdge.size(); ++ii)
              {
                Edge& edge = vecEdge[ii];
                ostr << "P[" << m_eMesh.get_rank() << "] t_edge= "
                     << m_eMesh.id(edge.first) << ", " << m_eMesh.id(edge.second) << "\n";
              }
            std::cout << ostr.str() << std::endl;
            vecEdge.resize(0);
          }
      }

    // now, find seams based on geometry
    for (unsigned ii=0; ii < vecFaces.size(); ++ii)
      {
        stk::mesh::Entity face = vecFaces[ii];

        const MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), face, stk::topology::NODE_RANK );

        typedef std::set<stk::mesh::Entity> EntitySet;
        EntitySet neighbors, shell_neighbors;
        m_eMesh.get_node_neighbors(face, neighbors, sel, m_eMesh.side_rank());
        m_eMesh.get_node_neighbors(face, shell_neighbors, sel, m_eMesh.element_rank());
        neighbors.insert(shell_neighbors.begin(), shell_neighbors.end());
        for (EntitySet::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
          {
            stk::mesh::Entity neigh = *it;
            if (neigh == face)
              continue;
            if (!is_surface_topology(m_eMesh.topology(neigh)))
              continue;

            int edge_0=0, edge_1=0;
            bool isEdgeN = m_eMesh.is_edge_neighbor(face, neigh, &edge_0, &edge_1);
            if (isEdgeN)
              {
                bool seam = isSeam(face, neigh, angleCriterion);
                if (seam)
                  {
                    m_edgeSeamsMap[face][edge_0] = 1;
                    m_edgeSeamsMap[neigh][edge_1] = 1;
                  }
              }
          }
      }

    if (1)
      {
        for (unsigned ii=0; ii < vecFaces.size(); ++ii)
          {
            stk::mesh::Entity face = vecFaces[ii];
            std::vector<int>& emap = m_edgeSeamsMap[face];
            const MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), face, stk::topology::NODE_RANK );
            for (unsigned iedge=0; iedge < emap.size(); ++iedge)
              {
                if (emap[iedge])
                  {
                    Edge edge = create_edge(face, iedge);
                    if (debug1) vecEdge.push_back(edge);

                    if (m_edgeSet.find(edge) == m_edgeSet.end())
                      {
                        m_edgeSet.insert(edge);
                      }
                  }
              }
          }
        if (debug1)
          {
            Util::makeUnique(vecEdge);
            std::ostringstream ostr;
            for (unsigned ii=0; ii < vecEdge.size(); ++ii)
              {
                Edge& edge = vecEdge[ii];
                ostr << "P[" << m_eMesh.get_rank() << "] g_edge= "
                     << m_eMesh.id(edge.first) << ", " << m_eMesh.id(edge.second) << "\n";
              }
            std::cout << ostr.str() << std::endl;
          }
      }

    // parallel share edges
    findGhostEdges(m_edgeSet, m_nodeToEdgeMap);

    m_edgeSetAll.insert(m_edgeSet.begin(), m_edgeSet.end());

    if (createEdgeSeamsPart && m_edgeSeamsPart)
      {
        m_eMesh.get_bulk_data()->modification_begin();
        stk::mesh::EntityId nedges = m_edgeSet.size();
        std::vector<stk::mesh::Entity> edges;
        m_eMesh.createEntities(m_eMesh.element_rank(), nedges, edges);
        nedges=0;
        for (EdgeSet::iterator it = m_edgeSet.begin(); it != m_edgeSet.end(); ++it)
          {
            m_eMesh.get_bulk_data()->declare_relation(edges[nedges], it->first, 0);
            m_eMesh.get_bulk_data()->declare_relation(edges[nedges], it->second, 1);
            stk::mesh::PartVector add(1, m_edgeSeamsPart), remove;
            m_eMesh.get_bulk_data()->change_entity_parts(edges[nedges], add, remove);
            ++nedges;
          }
        stk::mesh::fixup_ghosted_to_shared_nodes(*m_eMesh.get_bulk_data());
        m_eMesh.get_bulk_data()->modification_end();
      }
  }

  struct CompEdges {
    PerceptMesh& m_eMesh;
    stk::mesh::Entity m_centerNode;
    CompEdges(PerceptMesh& eMesh, const stk::mesh::Entity& centerNode) : m_eMesh(eMesh), m_centerNode(centerNode) {}
    bool operator()(const FitGregoryPatches::Edge& e0, const FitGregoryPatches::Edge& e1)
    {
      VERIFY_OP_ON((e0.first == m_centerNode), ||, (e0.second == m_centerNode), "bad node");
      VERIFY_OP_ON((e1.first == m_centerNode), ||, (e1.second == m_centerNode), "bad node");
      stk::mesh::Entity n0 = e0.first == m_centerNode ? e0.second : e0.first;
      stk::mesh::Entity n1 = e1.first == m_centerNode ? e1.second : e1.first;
      return ID(n0) < ID(n1);
    }
  };

  void FitGregoryPatches::
  sortEdges(NodeToEdgeMap& nodeToEdgeMap)
  {
    for (NodeToEdgeMap::iterator it = nodeToEdgeMap.begin(); it != nodeToEdgeMap.end(); ++it)
      {
        EdgeVector& edges = it->second;
        std::sort(edges.begin(), edges.end(), CompEdges(m_eMesh, it->first));
      }
  }

  void FitGregoryPatches::
  processSeams(stk::mesh::PartVector& parts, bool createEdgeSeamsParts)
  {
    m_contiguousEdgeSets.clear();
    m_nodeToEdgeMap.clear();

    findContiguousEdgeSets(m_contiguousEdgeSets, m_nodeToEdgeMap, m_edgeSet);

    sortEdges(m_nodeToEdgeMap);

    findEdgeSetCorners(m_corners, m_globalAngleCriterion, m_contiguousEdgeSets, m_nodeToEdgeMap);

    if (1)
      {
        size_t nc = m_corners.size();
        stk::all_reduce( m_eMesh.get_bulk_data()->parallel() , stk::ReduceSum<1>( & nc ) );
        if (m_eMesh.get_rank() == 0)
          std::cout << "P[" << m_eMesh.get_rank() << "] FitGregoryPatches::processSeams:  findEdgeSetCorners global size= " << nc << std::endl;
      }

    if (m_debug)
      {
        for (EntitySet::iterator it = m_corners.begin(); it != m_corners.end(); ++it)
          {
            std::cout << "P[" << m_eMesh.get_rank() << "] corner= " << m_eMesh.id(*it) << std::endl;
          }
      }

    findTangentVectors(m_corners, m_contiguousEdgeSets, m_nodeToEdgeMap);

    m_contiguousEdgeSetsAll.insert(m_contiguousEdgeSetsAll.end(), m_contiguousEdgeSets.begin(), m_contiguousEdgeSets.end());

    if (createEdgeSeamsParts)
      {
        addEdgesToQAMesh(m_contiguousEdgeSets);
      }
  }

  void FitGregoryPatches::
  findGhostEdges(EdgeSet& mainEdgeSet, NodeToEdgeMap& nodeToEdgeMap)
  {
    stk::CommSparse commAll (m_eMesh.parallel());
    std::vector<int>  procs;
    unsigned proc_size = m_eMesh.get_parallel_size();
    EdgeSet newEdgeSet;

    for (unsigned istage=0; istage < 2; ++istage)
      {
        // pack
        for (EdgeSet::iterator it = mainEdgeSet.begin(); it != mainEdgeSet.end(); ++it)
          {
            const Edge& edge = *it;
            stk::mesh::Entity nodes[2] { edge.first,  edge.second};
            for (unsigned ii=0; ii < 2; ++ii)
              {
                if (m_eMesh.shared(nodes[ii]))
                  {
                    m_eMesh.get_bulk_data()->comm_procs(nodes[ii], procs);
                    for (unsigned jj=0; jj < procs.size(); ++jj)
                      {
                        commAll.send_buffer( procs[jj] ).pack< stk::mesh::EntityId > (ID(nodes[ii]));
                        commAll.send_buffer( procs[jj] ).pack< stk::mesh::EntityId > (ID(nodes[(ii == 0 ? 1 : 0)]));
                      }
                  }
              }
          }

        if (istage == 0)
          {
            commAll.allocate_buffers();
          }
        else
          {

            commAll.communicate();

            // unpack
            for(unsigned from_proc = 0; from_proc < proc_size; ++from_proc )
              {
                stk::CommBuffer & recv_buffer = commAll.recv_buffer( from_proc );

                while ( recv_buffer.remaining() )
                  {
                    stk::mesh::EntityId nodeIds[2] = {0,0};
                    stk::mesh::Entity nodes[2];
                    recv_buffer.unpack< stk::mesh::EntityId >( nodeIds[0] );
                    recv_buffer.unpack< stk::mesh::EntityId >( nodeIds[1] );
                    nodes[0] = m_eMesh.get_bulk_data()->get_entity(m_eMesh.node_rank(), nodeIds[0]);
                    nodes[1] = m_eMesh.get_bulk_data()->get_entity(m_eMesh.node_rank(), nodeIds[1]);
                    if (m_eMesh.is_valid(nodes[0]) && m_eMesh.is_valid(nodes[1]))
                      {
                        Comp comp(m_eMesh);
                        Edge edge = create_edge<stk::mesh::Entity>(nodes[0], nodes[1], comp);
                        newEdgeSet.insert(edge);
                      }
                  }
              }
          }
      }
    mainEdgeSet.insert(newEdgeSet.begin(), newEdgeSet.end());
    nodeToEdgeMap.clear();

    for (EdgeSet::iterator it = mainEdgeSet.begin(); it != mainEdgeSet.end(); ++it)
      {
        nodeToEdgeMap[it->first].push_back( *it);
        nodeToEdgeMap[it->second].push_back( *it);
      }
  }


  void FitGregoryPatches::
  findContiguousEdgeSets(std::vector<EdgeSet>& contiguousEdgeSets, NodeToEdgeMap& nodeToEdgeMap,
                         const EdgeSet& mainEdgeSet)
  {
      {
#if HAVE_BOOST_GRAPH
        BGraphExternal::findContiguousEdgeSets(m_eMesh, contiguousEdgeSets, nodeToEdgeMap, mainEdgeSet);
#endif
        return;
      }

    nodeToEdgeMap.clear();

    EntitySet allNodeSet;
    for (EdgeSet::iterator it = mainEdgeSet.begin(); it != mainEdgeSet.end(); ++it)
      {
        nodeToEdgeMap[it->first].push_back( *it);
        nodeToEdgeMap[it->second].push_back( *it);
        allNodeSet.insert(it->first);
        allNodeSet.insert(it->second);
      }

    EntitySet visitedNodeSet;
    stk::mesh::Entity current_node = mainEdgeSet.begin()->first;
    bool continueEdgeSets = true;
    while(continueEdgeSets)
      {
        EdgeSet newEdgeSet;
        bool continueEdges = true;
        EntitySet newPossibleNodes;
        while(continueEdges)
          {
            visitedNodeSet.insert(current_node);
            allNodeSet.erase(current_node);

            EdgeVector& edgeVec = nodeToEdgeMap[current_node];
            for (unsigned ii = 0; ii < edgeVec.size(); ++ii)
              {
                Edge& edge = edgeVec[ii];
                if (newEdgeSet.find(edge) == newEdgeSet.end())
                  {
                    newEdgeSet.insert(edge);
                    if (visitedNodeSet.find(edge.first) == visitedNodeSet.end())
                      newPossibleNodes.insert(edge.first);
                    if (visitedNodeSet.find(edge.second) == visitedNodeSet.end())
                      newPossibleNodes.insert(edge.second);
                  }
              }
            if (newPossibleNodes.size() == 0)
              {
                contiguousEdgeSets.push_back(newEdgeSet);
                continueEdges = false;
              }
            else
              {
                current_node = *newPossibleNodes.begin();
                newPossibleNodes.erase(current_node);
              }
          }
        if (allNodeSet.size() == 0)
          {
            continueEdgeSets = false;
          }
        else
          {
            current_node = *allNodeSet.begin();
            allNodeSet.erase(current_node);
          }
      }
  }

  double FitGregoryPatches::edgeAngle(stk::mesh::Entity node, const Edge& e0, const Edge& e1)
  {
    stk::mesh::Entity nodes[2] = {e0.first == node ? e0.second : e0.first,
                                  e1.first == node ? e1.second : e1.first};
    double *nd0 = static_cast<double*>(stk::mesh::field_data( *m_eMesh.get_coordinates_field() , nodes[0] ));
    double *nd1 = static_cast<double*>(stk::mesh::field_data( *m_eMesh.get_coordinates_field() , nodes[1] ));
    double *nd = static_cast<double*>(stk::mesh::field_data( *m_eMesh.get_coordinates_field() , node ));

    double c0[3], c1[3], c[3];
    Math::copy_3d(c0, nd0);
    Math::copy_3d(c1, nd1);
    Math::copy_3d(c, nd);
    Math::subtract_3d(c1, c);
    Math::subtract_3d(c, c0);
    Math::normalize_3d(c);
    Math::normalize_3d(c1);
    double dot = Math::dot_3d(c, c1);
    double angle = M_PI - std::acos(dot);
    return angle*180.0/M_PI;
  }

  void FitGregoryPatches::
  findEdgeSetCorners(EntitySet& corners, double angleCriterion, const std::vector<EdgeSet>& contiguousEdgeSets, NodeToEdgeMap& nodeToEdgeMap)
  {
    EntitySet nodeSet;
    corners.clear();
    for (unsigned ii = 0; ii < contiguousEdgeSets.size(); ++ii)
      {
        const EdgeSet& edgeSet = contiguousEdgeSets[ii];

        for (EdgeSet::iterator it = edgeSet.begin(); it != edgeSet.end(); ++it)
          {
            const Edge& edge = *it;
            nodeSet.insert(edge.first);
            nodeSet.insert(edge.second);
          }
      }
    for (EntitySet::iterator it = nodeSet.begin(); it != nodeSet.end(); ++it)
      {
        const stk::mesh::Entity node = *it;
        EdgeVector& edgeVec = nodeToEdgeMap[node];
        if (edgeVec.size() == 2)
          {
            double ea = edgeAngle(node, edgeVec[0], edgeVec[1]);
            if (m_debug) std::cout << "corner for node: " << ID(node) << " ea= " << ea << " angleCriterion= " << angleCriterion << std::endl;
            if (ea < angleCriterion)
              {
                if (m_debug) std::cout << "corner for node: " << ID(node) << " ea= " << ea << " angleCriterion= " << angleCriterion << std::endl;
                corners.insert(node);
              }
          }
        else if (!m_eMesh.aura(node))
          {
            corners.insert(node);
            if (m_debug) std::cout << "corner for node: " << ID(node) << " edgeVecsize= " << edgeVec.size() << std::endl;
          }
      }
  }

  FitGregoryPatches::Point FitGregoryPatches::
  getTangent(stk::mesh::Entity n0, stk::mesh::Entity n1)
  {
    Point tangent;
    double *nd0 = static_cast<double*>(stk::mesh::field_data( *m_eMesh.get_coordinates_field() , n0 ));
    double *nd1 = static_cast<double*>(stk::mesh::field_data( *m_eMesh.get_coordinates_field() , n1 ));
    for (unsigned j=0; j < 3; ++j)
      {
        tangent[j] = nd1[j] - nd0[j];
      }
    Math::normalize_3d(tangent.data());
    return tangent;
  }

  std::string FitGregoryPatches::printEdge(const Edge& edge)
  {
    return toString(ID(edge.first))+" "+toString(ID(edge.second));
  }

  void FitGregoryPatches::
  findTangentVectors(const Edge& edge, const EntitySet& corners, const NodeToEdgeMap& nodeToEdgeMapIn)
  {
    bool debug = false;

    NodeToEdgeMap& nodeToEdgeMap = const_cast<NodeToEdgeMap&>(nodeToEdgeMapIn);
    /**
     *    n0E ---- n0 ---- n1 ---- n1E
     */
    stk::mesh::Entity n0 = edge.first;
    stk::mesh::Entity n1 = edge.second;

    if (m_eMesh.aura(n0) || m_eMesh.aura(n1))
      return;

    if (debug) std::cout << "nodeToEdgeMap.size= " << nodeToEdgeMap[n0].size() << " n1= " <<  nodeToEdgeMap[n1].size()
                         << " edge00= " << printEdge( nodeToEdgeMap[n0][0])
                         << " edge01= " << printEdge( nodeToEdgeMap[n0][1])
                         << " edge10= " << printEdge( nodeToEdgeMap[n1][0])
                         << " edge11= " << printEdge( nodeToEdgeMap[n1][1])
                         << std::endl;

    VERIFY_OP_ON(m_eMesh.id(n0), <, m_eMesh.id(n1), "bad id on edge");
    stk::mesh::Entity n0E = n0;
    stk::mesh::Entity n1E = n1;
    if (corners.find(n0) == corners.end())
      {
        VERIFY_OP_ON(nodeToEdgeMap[n0].size(), ==, 2, "bad nodeToEdgeMap");
        Edge e0E = nodeToEdgeMap[n0][0];
        if (e0E.first == n1 || e0E.second == n1)
          e0E = nodeToEdgeMap[n0][1];
        n0E = e0E.first;
        if (e0E.first == n0)
          n0E = e0E.second;
      }
    if (corners.find(n1) == corners.end())
      {
        VERIFY_OP_ON(nodeToEdgeMap[n1].size(), ==, 2, "bad nodeToEdgeMap");
        Edge e1E = nodeToEdgeMap[n1][0];
        if (e1E.first == n0 || e1E.second == n0)
          e1E = nodeToEdgeMap[n1][1];
        n1E = e1E.first;
        if (e1E.first == n1)
          n1E = e1E.second;
      }
    Point t0 = getTangent(n0, n1);
    if (n0E != n0)
      {
        Point t0E0 = getTangent(n0E, n0);
        for (unsigned j=0; j < 3; ++j)
          {
            t0[j] = 0.5*(t0[j] + t0E0[j]);
          }
      }
    Point t1 = getTangent(n0, n1);
    if (n1E != n1)
      {
        Point t1E1 = getTangent(n1, n1E);
        for (unsigned j=0; j < 3; ++j)
          {
            t1[j] = 0.5*(t1[j] + t1E1[j]);
          }
      }
    m_tangentVectors[edge].resize(2);
    m_tangentVectors[edge][0] = t0;
    m_tangentVectors[edge][1] = t1;
    if (debug)
      std::cout << "n0E .. n1E= "
                << ID(n0E) << " -- "
                << ID(n0) << " -- "
                << ID(n1) << " -- "
                << ID(n1E) << std::endl;

  }

  void FitGregoryPatches::
  findTangentVectors(const EntitySet& corners, const std::vector<EdgeSet>& contiguousEdgeSets, const NodeToEdgeMap& nodeToEdgeMap)
  {
    m_tangentVectors.clear();

    for (unsigned ii = 0; ii < contiguousEdgeSets.size(); ++ii)
      {
        const EdgeSet& edgeSet = contiguousEdgeSets[ii];

        for (EdgeSet::iterator it = edgeSet.begin(); it != edgeSet.end(); ++it)
          {
            findTangentVectors(*it, corners, nodeToEdgeMap);
          }
      }
  }

  void FitGregoryPatches::
  addEdgesToQAMesh(std::vector<EdgeSet>& contiguousEdgeSets)
  {
    typedef std::pair<stk::mesh::EntityId, stk::mesh::EntityId> EdgeId;
    typedef std::set<EdgeId> EdgeIdSet;
    std::vector<EdgeIdSet > contiguousEdgeIdSets(contiguousEdgeSets.size());

    for (unsigned ii=0; ii < contiguousEdgeSets.size(); ++ii)
      {
        EdgeSet& edgeSet = contiguousEdgeSets[ii];

        for (EdgeSet::iterator it = edgeSet.begin(); it != edgeSet.end(); ++it)
          {
            const Edge& edge = *it;
            EdgeId edgeId = create_edge<stk::mesh::EntityId>(m_eMesh.id(edge.first), m_eMesh.id(edge.second));
            contiguousEdgeIdSets[ii].insert(edgeId);
          }
      }

    m_eMesh.reopen();

    stk::mesh::PartVector edgeSeamsParts;
    for (unsigned ii=0; ii < contiguousEdgeIdSets.size(); ++ii)
      {
        stk::mesh::Part& part = m_eMesh.get_fem_meta_data()->declare_part_with_topology("edgeSeams_"+toString(ii), stk::topology::BEAM_2);
        stk::io::put_io_part_attribute(part);
        edgeSeamsParts.push_back(&part);
      }
    m_eMesh.commit();

    m_eMesh.get_bulk_data()->modification_begin();
    for (unsigned ii = 0; ii < contiguousEdgeIdSets.size(); ++ii)
      {
        EdgeIdSet& edgeSet = contiguousEdgeIdSets[ii];
        std::vector<stk::mesh::Entity> edges;
        m_eMesh.createEntities(m_eMesh.element_rank(), edgeSet.size(), edges);
        size_t nedges=0;
        for (EdgeIdSet::iterator it = edgeSet.begin(); it != edgeSet.end(); ++it)
          {
            const EdgeId& edge = *it;
            stk::mesh::Entity n0 = m_eMesh.get_bulk_data()->get_entity(m_eMesh.node_rank(), edge.first);
            stk::mesh::Entity n1 = m_eMesh.get_bulk_data()->get_entity(m_eMesh.node_rank(), edge.second);

            VERIFY_OP_ON(m_eMesh.is_valid(n0), ==, true, "bad n0");
            VERIFY_OP_ON(m_eMesh.is_valid(n1), ==, true, "bad n1");

            m_eMesh.get_bulk_data()->declare_relation(edges[nedges], n0, 0);
            m_eMesh.get_bulk_data()->declare_relation(edges[nedges], n1, 1);
            stk::mesh::PartVector add(1, edgeSeamsParts[ii]), remove;
            m_eMesh.get_bulk_data()->change_entity_parts(edges[nedges], add, remove);
            ++nedges;
          }
      }
    stk::mesh::fixup_ghosted_to_shared_nodes(*m_eMesh.get_bulk_data());
    m_eMesh.get_bulk_data()->modification_end();
  }

  void FitGregoryPatches::
  addEdgesToMesh(EdgeSet& edgeSet)
  {
    m_eMesh.get_bulk_data()->modification_begin();
    std::vector<stk::mesh::Entity> edges;
    m_eMesh.createEntities(m_eMesh.element_rank(), edgeSet.size(), edges);
    size_t nedges=0;
    for (EdgeSet::iterator it = edgeSet.begin(); it != edgeSet.end(); ++it)
      {
        const Edge& edge = *it;
        stk::mesh::Entity n0 = edge.first;
        stk::mesh::Entity n1 = edge.second;

        VERIFY_OP_ON(m_eMesh.is_valid(n0), ==, true, "bad n0");
        VERIFY_OP_ON(m_eMesh.is_valid(n1), ==, true, "bad n1");

        m_eMesh.get_bulk_data()->declare_relation(edges[nedges], n0, 0);
        m_eMesh.get_bulk_data()->declare_relation(edges[nedges], n1, 1);
        stk::mesh::PartVector add(1, m_edgeSeamsPart), remove;
        m_eMesh.get_bulk_data()->change_entity_parts(edges[nedges], add, remove);
        ++nedges;
      }
    stk::mesh::fixup_ghosted_to_shared_nodes(*m_eMesh.get_bulk_data());
    m_eMesh.get_bulk_data()->modification_end();
  }

  /** given two points on an edge, @param {pi,pj} and their associated
   *    @param tangents  {ti,tj}, fit a cubic with tangent to the tangents
   *     passing through the points, returning the control points in @param c
   */

  void  FitGregoryPatches::
  fitCubicWithTangents(MDArray& c, const MDArray& pi, const MDArray& pj, const Point& ti, const Point& tj)
  {
    double len = Math::distance_3d(pi.data(), pj.data());

    VERIFY_OP_ON(Math::norm_3d(&ti[0]), > , 1.e-8, "bad tangent");
    VERIFY_OP_ON(Math::norm_3d(&tj[0]), > , 1.e-8, "bad tangent j");
    VERIFY_OP_ON(len, > , 1.e-8, "bad tangent j");

    c(0,0) = pi(0);
    c(0,1) = pi(1);
    c(0,2) = pi(2);

    c(1,0) = pi(0) + ti[0]*len/3.0;
    c(1,1) = pi(1) + ti[1]*len/3.0;
    c(1,2) = pi(2) + ti[2]*len/3.0;

    c(2,0) = pj(0) - tj[0]*len/3.0;
    c(2,1) = pj(1) - tj[1]*len/3.0;
    c(2,2) = pj(2) - tj[2]*len/3.0;

    c(3,0) = pj(0);
    c(3,1) = pj(1);
    c(3,2) = pj(2);

    if (checkNAN(c))
      {
        VERIFY_MSG("bad fitCubicWithTangents");
      }

  }

  void FitGregoryPatches::
  fitCubics(stk::mesh::PartVector& parts)
  {
    stk::mesh::Selector sel = stk::mesh::selectUnion(parts) & m_eMesh.get_fem_meta_data()->locally_owned_part();
    stk::mesh::Selector selAll = stk::mesh::selectUnion(parts) & m_eMesh.get_fem_meta_data()->universal_part();

    std::vector<stk::mesh::Entity> vecFaces, vecShells;
    stk::mesh::get_selected_entities(sel , m_eMesh.get_bulk_data()->buckets(m_eMesh.side_rank()), vecFaces);
    stk::mesh::get_selected_entities(sel , m_eMesh.get_bulk_data()->buckets(m_eMesh.element_rank()), vecShells);
    vecFaces.insert(vecFaces.end(), vecShells.begin(), vecShells.end());
    for (unsigned ii=0; ii < vecFaces.size(); ++ii)
      {
        stk::mesh::Entity face = vecFaces[ii];
        double *Cp = (m_eMesh.entity_rank(face) == m_eMesh.side_rank()
                      ? stk::mesh::field_data( *m_eMesh.m_gregory_control_points_field, face)
                      : stk::mesh::field_data( *m_eMesh.m_gregory_control_points_field_shell, face));

        const MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), face, stk::topology::NODE_RANK );
        bool isTri = (face_nodes.size() == 3);

        for (unsigned kn = 0; kn < face_nodes.size(); ++kn)
          {
            stk::mesh::Entity node = face_nodes[kn].entity();
            stk::mesh::Entity nodep = face_nodes[(kn+1)%face_nodes.size()].entity();
            bool debug = false; // (m_eMesh.id(node) == 31 && m_eMesh.id(nodep) == 5) || (m_eMesh.id(node) == 5 && m_eMesh.id(nodep) == 31);
            if (debug)
              {
                std::cout << "kn= " << kn << " face= " << m_eMesh.print_entity_compact(face) << std::endl;
              }
            int dim = 3;
            double *normals_data_0 = stk::mesh::field_data( *m_eMesh.m_node_normals , node );
            double *normals_data_p0 = stk::mesh::field_data( *m_eMesh.m_node_normals , nodep );
            MDArray n_0 (normals_data_0, 1, dim);
            MDArray np_0 (normals_data_p0, 1, dim);

            // maybe modified below
            MDArray n("n",n_0.layout());
            Kokkos::deep_copy(n, n_0); 
            MDArray np("np",np_0.layout());
            Kokkos::deep_copy(np, np_0); 
            MDArray c (static_cast<double*>(stk::mesh::field_data( *m_eMesh.get_coordinates_field() , node )), 1, dim);
            MDArray cp (static_cast<double*>(stk::mesh::field_data( *m_eMesh.get_coordinates_field() , nodep )), 1, dim);
            MDArray cf("cf",4,3);

            // if (not a seam, but one of my nodes on the edge is on a seam, use face normal/special average using only the faces across my edge)
            // prep here for extrapolation...

            if (m_edgeSeamsMap[face].size() && m_edgeSeamsMap[face][kn])
              {
                Comp comp(m_eMesh);
                Edge edge = create_edge<stk::mesh::Entity>(node, nodep, comp);
                std::vector<Point> tangents = m_tangentVectors[edge];
                VERIFY_OP_ON(tangents.size(), ==, 2, "bad tangent");
                if (m_eMesh.id(nodep) < m_eMesh.id(node))
                  {
                    VERIFY_OP_ON(edge.first, ==, nodep, "bad nodep");
                    VERIFY_OP_ON(edge.second, ==, node, "bad node");
                    MDArray cfr("cfr",4,3);
                    fitCubicWithTangents(cfr, cp, c, tangents[0], tangents[1]);
                    for (unsigned icp = 0; icp < 4; ++icp)
                      {
                        cf(icp, 0) = cfr(3-icp, 0);
                        cf(icp, 1) = cfr(3-icp, 1);
                        cf(icp, 2) = cfr(3-icp, 2);
                      }
                  }
                else
                  {
                    fitCubicWithTangents(cf, c, cp, tangents[0], tangents[1]);
                  }
                if (checkNAN(cf) || debug)
                  std::cout << "cf=\n " << printContainer(cf) << "\ntangents= "
                            << Math::print_3d(&tangents[0][0]) << " "
                            << Math::print_3d(&tangents[1][0]) << std::endl;
              }
            else
              {
                if (1)
                  {
                    int orient0 = orient(face, node, selAll, debug);
                    int orient1 = orient(face, nodep, selAll, debug);

                    if (debug) {
                      std::cout << "P[" << m_eMesh.get_rank() << "] orient0= " << orient0 << " orient1= " << orient1
                                << " m_nodeToEdgeMap[node]= " << m_nodeToEdgeMap[node].size() << " nodep= " << m_nodeToEdgeMap[nodep].size()
                                << std::endl;
                    }

                    // set the node normal to be the face normal if its not oriented, is a corner, etc
                    PointVector normals0 = m_nodeNormalsMap[node];
                    PointVector normals1 = m_nodeNormalsMap[nodep];
                    VERIFY_OP_ON(orient0, <, int(normals0.size()), "bad normals");
                    VERIFY_OP_ON(orient1, <, int(normals1.size()), "bad normals");

                    if (orient0 >= 0)
                      {
                        Math::copy_3d(n.data(), normals0[orient0].data());
                        VERIFY_OP_ON(Math::norm_3d(n.data()), >, 1.e-8, "bad norm");
                      }
                    if (orient1 >= 0)
                      {
                        Math::copy_3d(np.data(), normals1[orient1].data());
                        VERIFY_OP_ON(Math::norm_3d(np.data()), >, 1.e-8, "bad normp");
                      }

                  }
                if (m_debug) std::cout << "P[" << m_eMesh.get_rank() << "] node= " << m_eMesh.id(node)
                                       << " norm= " << Math::norm_3d(n.data())
                                       << " normp= " << Math::norm_3d(np.data()) << std::endl;
                VERIFY_OP_ON(Math::norm_3d(n.data()), > , 1.e-8, "bad norm");
                VERIFY_OP_ON(Math::norm_3d(np.data()), > , 1.e-8, "bad normp");

                GregoryPatch::fitCubic(cf, c, cp, n, np);
                if (checkNAN(cf))
                  {
                    std::cout << "cf= \n" << printContainer(cf) << " n=\n" << printContainer(n) << " np=\n" << printContainer(np) << std::endl;
                  }

                if (debug)
                  {
                    double c0[3], c1[3], c2[3], c3[3];
                    for (unsigned jj=0; jj< 3; jj++)
                      {
                        c0[jj] = cf(0,jj);
                        c1[jj] = cf(1,jj);
                        c2[jj] = cf(2,jj);
                        c3[jj] = cf(3,jj);
                      }
                    if (Math::distance_3d(c0,c1)< 1.e-8)
                      {
                        std::cout << "P[" << m_eMesh.get_rank() << "] face= " << m_eMesh.id(face) << " kn= " << kn
                                  << " node= " << m_eMesh.id(node) << " nodep= " << m_eMesh.id(nodep)
                                  << "\nn=\n" << printContainer(n) << "\nnp=\n" << printContainer(np) << "\nc=\n" << printContainer(c) << "\ncp=\n" << printContainer(cp) << std::endl;
                      }
                    VERIFY_OP_ON(Math::distance_3d(c0,c1), >, 1.e-8, "bad cf");
                    VERIFY_OP_ON(Math::distance_3d(c2,c3), >, 1.e-8, "bad cf2");
                  }

              }
            if (isTri)
              {
                MDArray qcf("qcf",5,3);
                GregoryPatch::degree_elevate(cf, qcf);
                for (unsigned ip=0; ip < 5; ++ip)
                  {
                    for (unsigned jc=0; jc < 3; ++jc)
                      {
                        double cpi = qcf(ip, jc);
                        Cp[percept::gregory_patch::tri_edges[kn][ip] + jc*MaxControlPoints()] = cpi;
                      }
                  }
              }
            else
              {
                for (unsigned ip=0; ip < 4; ++ip)
                  {
                    for (unsigned jc=0; jc < 3; ++jc)
                      {
                        double cpi = cf(ip, jc);
                        Cp[percept::gregory_patch::quad_edges[kn][ip] + jc*MaxControlPoints()] = cpi;
                      }
                  }
              }
          }
      }

    std::vector<const stk::mesh::FieldBase *> fields;
    fields.push_back(m_eMesh.m_gregory_control_points_field_shell);
    fields.push_back(m_eMesh.m_gregory_control_points_field);
    stk::mesh::communicate_field_data( m_eMesh.get_bulk_data()->aura_ghosting() ,    fields);
  }

  void FitGregoryPatches::
  extractRibbon(stk::mesh::Entity face, int edge, bool reverse, MDArray& p, MDArray& q)
  {
    const MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), face, stk::topology::NODE_RANK );
    bool isTri = (face_nodes.size() == 3);
    double *Cp = (m_eMesh.entity_rank(face) == m_eMesh.side_rank()
                  ? stk::mesh::field_data( *m_eMesh.m_gregory_control_points_field, face)
                  : stk::mesh::field_data( *m_eMesh.m_gregory_control_points_field_shell, face));

    for (unsigned ip=0; ip < 4; ++ip)
      {
        if (isTri)
          {
            int *ribbon = percept::gregory_patch::tri_ribbons[edge][0][ip];
            int *ribbon_l = percept::gregory_patch::tri_ribbons[edge][1][ip];
            for (unsigned jc=0; jc < 3; ++jc)
              {
                if (ribbon_l[0] == ribbon_l[1])
                  {
                    p(ip,jc) = Cp[ribbon_l[0] + jc*MaxControlPoints()];
                  }
                else
                  {
                    // Degree-lowering
                    p(ip,jc) = 1./3.*(4*Cp[ribbon_l[1] + jc*MaxControlPoints()] - Cp[ribbon_l[0] + jc*MaxControlPoints()]);
                  }
                if (ribbon[0] == ribbon[1])
                  {
                    q(ip,jc) = Cp[ribbon[0] + jc*MaxControlPoints()];
                  }
                else
                  {
                    // Degree-lowering
                    q(ip,jc) = 1./3.*(4*Cp[ribbon[1] + jc*MaxControlPoints()] - Cp[ribbon[0] + jc*MaxControlPoints()]);
                  }
              }
          }
        else
          {
            int ribbon = percept::gregory_patch::quad_ribbons[edge][0][ip];
            int ribbon_l = percept::gregory_patch::quad_ribbons[edge][1][ip];
            for (unsigned jc=0; jc < 3; ++jc)
              {
                p(ip,jc) = Cp[ribbon_l + jc*MaxControlPoints()];
                q(ip,jc) = Cp[ribbon + jc*MaxControlPoints()];
              }
          }
      }
    if (reverse)
      {
        MDArray pr("pr",p.layout()), qr("qr",q.layout());
        Kokkos::deep_copy(pr,p);
        Kokkos::deep_copy(qr,q);
        for (unsigned ip=0; ip < 4; ++ip)
          {
            for (unsigned jc=0; jc < 3; ++jc)
              {
                p(3 - ip, jc) = pr(ip, jc);
                q(3 - ip, jc) = qr(ip, jc);
              }
          }
      }
  }

  void FitGregoryPatches::
  putRibbon(stk::mesh::Entity face, int edge, bool reverse, MDArray& p_in)
  {
    const MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), face, stk::topology::NODE_RANK );
    bool isTri = (face_nodes.size() == 3);
    //double *Cp = stk::mesh::field_data( *m_eMesh.m_gregory_control_points_field, face);
    double *Cp = (m_eMesh.entity_rank(face) == m_eMesh.side_rank()
                  ? stk::mesh::field_data( *m_eMesh.m_gregory_control_points_field, face)
                  : stk::mesh::field_data( *m_eMesh.m_gregory_control_points_field_shell, face));
    MDArray p("ploc",p_in.layout());
    Kokkos::deep_copy(p,p_in);
    if (reverse)
      {
        for (unsigned ip=0; ip < 4; ++ip)
          {
            for (unsigned jc=0; jc < 3; ++jc)
              {
                p(3 - ip, jc) = p_in(ip, jc);
              }
          }
      }

    for (unsigned ip=1; ip < 3; ++ip)
      {
        if (isTri)
          {
            int *ribbon_l = percept::gregory_patch::tri_ribbons[edge][1][ip];
            for (unsigned jc=0; jc < 3; ++jc)
              {
                VERIFY_OP_ON(ribbon_l[0], ==, ribbon_l[1], "bad ribbon data");
                Cp[ribbon_l[0] + jc*MaxControlPoints()] = p(ip, jc);
              }
          }
        else
          {
            int ribbon_l = percept::gregory_patch::quad_ribbons[edge][1][ip];
            for (unsigned jc=0; jc < 3; ++jc)
              {
                Cp[ribbon_l + jc*MaxControlPoints()] = p(ip,jc);
              }
          }
      }
  }

  void FitGregoryPatches::
  fitRibbons(stk::mesh::PartVector& parts)
  {
    bool reverseAll = m_reverseAll;
    stk::mesh::Selector sel = stk::mesh::selectUnion(parts) & m_eMesh.get_fem_meta_data()->locally_owned_part();
    stk::mesh::Selector selAll = stk::mesh::selectUnion(parts);

    std::vector<stk::mesh::Entity> vecFaces, vecShells;
    stk::mesh::get_selected_entities(sel , m_eMesh.get_bulk_data()->buckets(m_eMesh.side_rank()), vecFaces);
    stk::mesh::get_selected_entities(sel , m_eMesh.get_bulk_data()->buckets(m_eMesh.element_rank()), vecShells);
    vecFaces.insert(vecFaces.end(), vecShells.begin(), vecShells.end());

    for (unsigned ii=0; ii < vecFaces.size(); ++ii)
      {
        stk::mesh::Entity face = vecFaces[ii];
        bool debug = false; // m_eMesh.id(face) == 654;

        if (debug) std::cout << "P[" << m_eMesh.get_rank() << " FGP:: data for face= " << m_eMesh.identifier(face) << "\n" << printForMathematica(face) << std::endl;

        const MyPairIterRelation face_nodes(*m_eMesh.get_bulk_data(), face, stk::topology::NODE_RANK );
        bool isTri = (face_nodes.size() == 3);
        typedef std::set<stk::mesh::Entity> EntitySet;
        EntitySet neighbors, shell_neighbors;
        std::vector<int> edge_visited(face_nodes.size(), 0);
        std::vector<int> edge_is_seam(face_nodes.size(), 0);

        for (unsigned edge=0; edge < face_nodes.size(); ++edge)
          {
            std::vector<int>& emap = m_edgeSeamsMap[face];
            int seam = 0;
            if (emap.size())
              {
                VERIFY_OP_ON(edge, <, emap.size(), "bad m_edgeSeamsMap");
                seam = m_edgeSeamsMap[face][edge];
                if (seam)
                  {
                    edge_visited[edge] = 1;
                    edge_is_seam[edge] = 1;
                  }
              }
          }

        m_eMesh.get_node_neighbors(face, neighbors, selAll, m_eMesh.side_rank());
        m_eMesh.get_node_neighbors(face, shell_neighbors, selAll, m_eMesh.element_rank());
        neighbors.insert(shell_neighbors.begin(), shell_neighbors.end());
        if (debug) std::cout << "P[" << m_eMesh.get_rank() << "] neighbors.size= " << neighbors.size() << std::endl;
        for (EntitySet::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
          {
            stk::mesh::Entity neigh = *it;
            if (!is_surface_topology(m_eMesh.topology(neigh)))
              continue;

            const MyPairIterRelation neigh_nodes(*m_eMesh.get_bulk_data(), neigh, stk::topology::NODE_RANK );
            bool neighIsTri = (neigh_nodes.size() == 3);
            int edge_0=0, edge_1=0;
            bool isEdgeN = m_eMesh.is_edge_neighbor(face, neigh, &edge_0, &edge_1);
            if (isEdgeN)
              {
                if (int(m_edgeSeamsMap[face].size()) > edge_0 && int(m_edgeSeamsMap[neigh].size()) > edge_1)
                  {
                    VERIFY_OP_ON(m_edgeSeamsMap[face][edge_0], ==, m_edgeSeamsMap[neigh][edge_1], "bad m_edgeSeamsMap");
                    VERIFY_OP_ON(m_edgeSeamsMap[face][edge_0], ==, edge_is_seam[edge_0], "bad m_edgeSeamsMap");
                  }
              }

            if (isEdgeN && !edge_visited[edge_0] && !edge_is_seam[edge_0])
              {
                if (debug) {
                  std::cout << "P[" << m_eMesh.get_rank() << " FGP:: doing neigh ribbons face= " << m_eMesh.identifier(face) << " neigh= " << m_eMesh.identifier(neigh) << std::endl;
                }
                edge_visited[edge_0] = 1;
                MDArray p("p",4,3), q("q",4,3), r("r",4,3), qh("qh",5,3), qcheck("qcheck",4,3);
                extractRibbon(face, edge_0, reverseAll, p, q);
                extractRibbon(neigh, edge_1, !reverseAll, r, qcheck);
                if (checkNAN(p) || checkNAN(q) || checkNAN(r))
                  {
                    std::cout << "p=\n" << printContainer(p) << " q=\n" << printContainer(q) << " r=\n" << printContainer(r) << std::endl;
                    VERIFY_MSG("p bad");
                  }
                GregoryPatch::fitRibbon(p, q, r, qh, isTri, neighIsTri);
                putRibbon(face, edge_0, reverseAll, p);
                putRibbon(neigh, edge_1, !reverseAll, r);
              }
          }

        // check for no neighbors across edge
        for (unsigned edge=0; edge < face_nodes.size(); ++edge)
          {
            if (!edge_visited[edge] || edge_is_seam[edge])
              {
                if (debug) std::cout << "P[" << m_eMesh.get_rank() << " FGP:: doing non-neigh ribbons face= " << m_eMesh.identifier(face) << " edge= " << edge << std::endl;
                MDArray p("p",4,3), q("q",4,3), qh("qh",5,3);
                extractRibbon(face, edge, reverseAll, p, q);
                MDArray pex("pex",p.layout());
                Kokkos::deep_copy(pex,p);
                GregoryPatch::fitRibbonNoNeighbor(p, q, qh, isTri);
                putRibbon(face, edge, reverseAll, p);
              }
          }

        if (debug) std::cout << "P[" << m_eMesh.get_rank() << " FGP:: post ribbon data for face= " << m_eMesh.identifier(face) << "\n" << printForMathematica(face) << std::endl;

      }
    std::vector<const stk::mesh::FieldBase *> fields;
    fields.push_back(m_eMesh.m_gregory_control_points_field_shell);
    fields.push_back(m_eMesh.m_gregory_control_points_field);
    stk::mesh::communicate_field_data( m_eMesh.get_bulk_data()->aura_ghosting() ,    fields);
  }

  //static
  std::string FitGregoryPatches::
  printForMathematica(PerceptMesh& m_eMesh, stk::mesh::Entity face, bool convert, bool printHeader, unsigned lineLength)
  {
    unsigned nn = m_eMesh.get_bulk_data()->num_nodes(face);
    double *Cp = (m_eMesh.entity_rank(face) == m_eMesh.side_rank()
                  ? stk::mesh::field_data( *m_eMesh.m_gregory_control_points_field, face)
                  : stk::mesh::field_data( *m_eMesh.m_gregory_control_points_field_shell, face));

    int npts = (nn == 3 ? NumTriControlPoints() : NumQuadControlPoints());
    std::ostringstream str;

    if (!convert) str << std::setprecision(10);

    if (printHeader)
      str << "face[" << m_eMesh.identifier(face) << "]= ";
    str << "{";
    std::ostringstream istr;
    for (int ii=0; ii < npts; ++ii)
      {
        istr << "{";
        for (int jc=0; jc < 3; ++jc)
          {
            if (convert)
              istr << Util::convert_to_mm(Cp[ii + jc*MaxControlPoints()]);
            else
              istr << Cp[ii + jc*MaxControlPoints()];
            istr << (jc == 2 ? "}" : ",");
          }
        istr << (ii == npts-1 ? (printHeader ? "};" : "},") : ",");
        if (istr.str().length() > lineLength || ii == npts-1)
          {
            istr << "\n";
            str << istr.str();
            istr.clear();
            istr.str("");
          }
      }
    return str.str();
  }

  void FitGregoryPatches::parse(const std::string& file_name)
  {
    if (file_name == "sample.yaml")
      {
        std::ofstream file(file_name.c_str());
        if (!file.good())
          {
            throw std::runtime_error("couldn't open file: "+file_name);
          }
        file
          << "debug: true\n"
          << "globalAngleCriterion: 135 # specified in degrees - this value is used unless overridden below\n"
          << "# the following can be commented out in which case all surfaces will be in one set\n"
          << "surface_sets:\n"
          << "  - set_1: [surface_1, surface_2]\n"
          << "  - set_2: [surface_3, surface_4, surface_10]\n"
          << "# the following can be commented out in which case all surfaces will use globalAngleCriterion\n"
          << "angle_map:\n"
          << "  set_1: 10.0\n"
          << "  surface_3: 9.0\n"
          << "#Note: if QA.activate is set then the surfaces will be\n"
          << "#  converted to shells, and all seams will be detected and put in a\n"
          << "#  special edge part and saved to the somePrefix.e specified - this is\n"
          << "#  for QA'ing of the input data before doing the actual fit.\n"
          << "QA:\n"
          << "  activate: false\n"
          << "  file: somePrefix\n"
          << "  num_divisions: 5 # Note: if > 0, creates a div_someFile_surfName.e showing resulting geometry\n"
          << "#To fire up a visualizer to see the fit (serial runs only), do one of:\n"
          << "  #visualizer_command_prefix: \"ensight \" # note the space at the end\n"
          << "  #visualizer_command_prefix: \"paraview --data=\"\n"
          << std::endl;
        return;
      }
    m_surfaceSets.clear();
    m_angleMap.clear();
    std::ifstream file(file_name.c_str());
    if (!file.good())
      {
        throw std::runtime_error("couldn't open file: "+file_name);
      }

#if HAVE_YAML
    try {
      m_node = YAML::Load(file);
      parse(m_node);
      if (m_debug)
        emit(m_node);
    }
    catch(YAML::ParserException& e) {
      std::cout << e.what() << " input= " << file_name << "\n";
    }
#endif
  }

#if HAVE_YAML
  void FitGregoryPatches::emit(const YAML::Node& node)
    {
      if (m_eMesh.get_rank() != 0) return;

      std::string emit_file = "emit.yaml";
      if (emit_file.length())
        {
          std::ofstream fout(emit_file.c_str());
          if (fout.good())
            {
              YamlUtils::emit(fout, node);
            }
        }
    }

  void FitGregoryPatches::parse(const YAML::Node& node)
    {
#define SIP2(a, Default, node, This) do { set_if_present(node, #a, This.m_ ## a, Default); \
        if (m_debug && m_eMesh.get_rank()==0) std::cout <<  EXPAND_AND_QUOTE(TOKENPASTE(m_, a) ) << " = " << This. TOKENPASTE(m_,a) << std::endl; } while (0)
#define SIP(a, Default) SIP2(a, Default, node, (*this))

      SIP(debug, bool(false)); // m_debug
      // done twice, once to turn on debugging, once to see the messages
      SIP(debug, bool(false));
      SIP(globalAngleCriterion, double(m_globalAngleCriterionDefault));

      const YAML::Node y_QA = node["QA"];
      if (y_QA)
        {
          SIP2(activate, bool(false), (y_QA), (this->m_QA) );
          SIP2(file, std::string(""), (y_QA), (this->m_QA) );
          SIP2(num_divisions, int(5), (y_QA), (this->m_QA) );
        }

      const YAML::Node y_surface_sets = node["surface_sets"];

      if (y_surface_sets)
        {
          VERIFY_OP_ON(y_surface_sets.Type(), ==, YAML::NodeType::Sequence, "bad surface_sets data in yaml file");
          //set_if_present(*y_wedge, "activate", wedge_boundary_layer_special_refinement, false);
          for (unsigned iSurfaceSet = 0; iSurfaceSet < y_surface_sets.size(); ++iSurfaceSet)
            {
              const YAML::Node & y_surface_set = y_surface_sets[iSurfaceSet];
              VERIFY_OP_ON(y_surface_set.Type(), ==, YAML::NodeType::Map, "bad surface_set data");
              for (YAML::const_iterator i = y_surface_set.begin(); i != y_surface_set.end(); ++i)
                {
                  const YAML::Node key   = i->first;
                  const YAML::Node value = i->second;
                  std::string v_key;
                  v_key = key.as<std::string>();
                  VERIFY_OP_ON(value.Type(), ==, YAML::NodeType::Sequence, "bad surface_set value data in [surfaceSetName: [s1,s2...]]");
                  for (unsigned jj=0; jj < value.size(); ++jj)
                    {
                      std::string ss;
                      ss = value[jj].as<std::string>();
                      m_surfaceSets[v_key].push_back(ss);
                    }
                }
            }
        }

      const YAML::Node y_angle_map = node["angle_map"];

      if (y_angle_map)
        {
          VERIFY_OP_ON(y_angle_map.Type(), ==, YAML::NodeType::Map, "bad angle_map data in yaml file");
          for (YAML::const_iterator i = y_angle_map.begin(); i != y_angle_map.end(); ++i)
            {
              const YAML::Node key   = i->first;
              const YAML::Node value = i->second;
              std::string v_key = key.as<std::string>();
              double v_value = value.as<double>();
              m_angleMap[v_key] = v_value;
            }
        }

#undef SIP
#undef SIP2
    }
#endif

  static double angle_deg(double *n0, double *n1)
  {
    double dot = n0[0]*n1[0] + n0[1]*n1[1] + n0[2]*n1[2];
    double angle = M_PI - std::acos(dot);
    return angle*180.0/M_PI;
  }

  bool FitGregoryPatches::
  isSeam(stk::mesh::Entity face_0, stk::mesh::Entity face_1, double angleCriterion)
  {
    double normal_0[3], normal_1[3];
    faceNormal(face_0, normal_0);
    faceNormal(face_1, normal_1);
    double angle =  angle_deg(normal_0, normal_1);
    if (m_debug) std::cout << "face_0= " << m_eMesh.id(face_0) << " face_1= " << m_eMesh.id(face_1) << " angle = " << angle << " angleCriterion= " << angleCriterion << std::endl;
    if (angleCriterion == 0.0)
      return false;
    return angle < angleCriterion;
  }

  void FitGregoryPatches::
  createQA(const std::string& QA_file)
  {
    std::vector<stk::mesh::PartVector> currentParts;
    getCurrentParts(currentParts);

    m_QA.m_filenames.resize(0);

    for (unsigned iSurfaceSet = 0; iSurfaceSet < currentParts.size(); ++iSurfaceSet)
      {
        const MyString *surfaceSetNameMS = currentParts[iSurfaceSet][0]->attribute<MyString>();
        VERIFY_OP_ON(surfaceSetNameMS, !=, 0, "null MyString");
        std::string surfaceSetName = surfaceSetNameMS->m_string;
        std::string file = QA_file+"_"+(surfaceSetName)+".e";
        std::string file_div = "div_"+file;
        m_QA.m_filenames.push_back(file);

        if (m_QA.m_num_divisions > 0)
          PerceptMesh::create_refined_mesh(m_eMesh, file_div, m_QA.m_num_divisions, &currentParts[iSurfaceSet], false);

        PerceptMesh eMesh;
        eMesh.openEmpty(false);
        m_eMesh.convertSurfacesToShells1_meta(eMesh, currentParts[iSurfaceSet]);

        SurfaceSets localSurfaceSet;
        for (unsigned ii=0; ii < currentParts[iSurfaceSet].size(); ++ii)
          {
            localSurfaceSet[surfaceSetName].push_back(currentParts[iSurfaceSet][ii]->name());
          }
        AngleMap localAngleMap;
        localAngleMap[surfaceSetName] = m_angleMap[surfaceSetName];
        FitGregoryPatches tmpFitter(eMesh, localSurfaceSet, localAngleMap, m_globalAngleCriterion);
        tmpFitter.m_QA = m_QA;
        tmpFitter.register_or_set_fields();
        eMesh.commit();

        m_eMesh.convertSurfacesToShells1_bulk(eMesh, currentParts[iSurfaceSet]);

        std::vector<stk::mesh::FieldBase *> nodal_fields_to_copy(1, m_eMesh.m_node_normals);
        std::vector<stk::mesh::FieldBase *> nodal_fields_to_copy_to(1, eMesh.m_node_normals);
        PerceptMesh::copy_nodal_fields(m_eMesh, eMesh, nodal_fields_to_copy, nodal_fields_to_copy_to);

        std::vector<stk::mesh::PartVector> currentPartsNew;
        tmpFitter.getCurrentParts(currentPartsNew);
        VERIFY_OP_ON(currentPartsNew.size(), ==, 1, "bad size");
        tmpFitter.findSeams(currentPartsNew[0], true);
        tmpFitter.processSeams(currentPartsNew[0], true);

        eMesh.save_as(file);
      }
  }

}
