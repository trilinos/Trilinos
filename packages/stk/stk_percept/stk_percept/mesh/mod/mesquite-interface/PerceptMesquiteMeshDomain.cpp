//-------------------------------------------------------------------------
// Filename      : PerceptMesquiteMeshDomain.cpp
//
// Purpose       : mesh domain interface for using Mesquite 
//
// Description   : subclass of Mesquite::MeshDomain
//
// Creator       : Steve Kennon, derived from Steve Owen's work
//
// Creation Date : Dec 2011
//
// Owner         : Steve Kennon
//-------------------------------------------------------------------------

#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)

#include "PerceptMesquiteMeshDomain.hpp"
//#include <mesquite/MsqVertex.hpp>
#include <MsqVertex.hpp>

#include <map>
#include <algorithm>

namespace stk {
  namespace percept {

#define PRINT_ERROR(a) do { std::cout << "PerceptMesquite::Mesh::ERROR: " << a << std::endl; } while (0)

    //! Modifies "coordinate" so that it lies on the
    //! domain to which "entity_handle" is constrained.
    //! The handle determines the domain.  The coordinate
    //! is the proposed new position on that domain.
    void PerceptMesquiteMeshDomain::
    snap_to(Mesquite::Mesh::VertexHandle entity_handle,
            Mesquite::Vector3D &coordinate) const {

      if (!m_meshGeometry) return;
      stk::mesh::Entity* node_ptr = reinterpret_cast<stk::mesh::Entity *>(entity_handle);
      stk::mesh::FieldBase* field = m_eMesh->get_coordinates_field();
      double *f_data = PerceptMesh::field_data(field, *node_ptr);

      double f_data_save[3] = {f_data[0], f_data[1], 0};
      if (m_eMesh->get_spatial_dim() > 2) f_data_save[2] = f_data[2];

//       if (node_ptr->identifier() == 584)
//         {
//           std::cout << "tmp snap_to: node= " << node_ptr->identifier() << std::endl;
//         }

      f_data[0] = coordinate[0];
      f_data[1] = coordinate[1];
      if (m_eMesh->get_spatial_dim() > 2) 
        f_data[2] = coordinate[2];

      static std::vector<stk::mesh::Entity *> nodes(1);
      nodes[0] = node_ptr;
      m_meshGeometry->snap_points_to_geometry(m_eMesh, nodes);
      coordinate[0] = f_data[0];
      coordinate[1] = f_data[1];
      if (m_eMesh->get_spatial_dim() > 2)
        coordinate[2] = f_data[2];

      //if (node_ptr->identifier() == 584)
      if (0)
        {
          std::cout << "tmp snap_to: node= " << node_ptr->identifier() << " orig= " 
                    << f_data_save[0] << " "
                    << f_data_save[1] << " "
                    << f_data_save[2] << " "
                    << " new= " 
                    << f_data[0] << " "
                    << f_data[1] << " "
                    << f_data[2] << " diff1= " 
                    << (std::fabs(f_data[0] - f_data_save[0])+
                        std::fabs(f_data[1] - f_data_save[1])+
                        std::fabs(f_data[2] - f_data_save[2]))
                    << std::endl;
        }

      f_data[0] = f_data_save[0];
      f_data[1] = f_data_save[1];
      if (m_eMesh->get_spatial_dim() > 2)
        f_data[2] = f_data_save[2];

    }
    
    //! Returns the normal of the domain to which
    //! "entity_handle" is constrained.  For non-planar surfaces,
    //! the normal is calculated at the point on the domain that
    //! is closest to the passed in value of "coordinate".  If the
    //! domain does not have a normal, or the normal cannot
    //! be determined, "coordinate" is set to (0,0,0).  Otherwise,
    //! "coordinate" is set to the domain's normal at the
    //! appropriate point.
    //! In summary, the handle determines the domain.  The coordinate
    //! determines the point of interest on that domain.
    //!
    //! User should see also PatchData::get_domain_normal_at_vertex and
    //! PatchData::get_domain_normal_at_element .
    void PerceptMesquiteMeshDomain::
    vertex_normal_at(Mesquite::Mesh::VertexHandle entity_handle,
                     Mesquite::Vector3D &coordinate) const {
      // FIXME srk
      coordinate[0] = 0.0;
      coordinate[1] = 0.0;
      coordinate[2] = 1.0;
      if (m_eMesh->get_spatial_dim() == 3)
        {
          stk::mesh::Entity* node_ptr = reinterpret_cast<stk::mesh::Entity *>(entity_handle);
          //stk::mesh::FieldBase* field = m_eMesh->get_coordinates_field();
          //double *f_data = PerceptMesh::field_data(field, *node_ptr);
          if (!m_meshGeometry)
            {
              return;
            }

          std::vector<double> normal(3,0.0);
          m_meshGeometry->normal_at(m_eMesh, node_ptr, normal);
          coordinate[0] = normal[0];
          coordinate[1] = normal[1];
          coordinate[2] = normal[2];
        }
    }

    void PerceptMesquiteMeshDomain::
    element_normal_at(Mesquite::Mesh::ElementHandle entity_handle,
                      Mesquite::Vector3D &coordinate) const {
      // FIXME srk
      coordinate[0] = 0.0;
      coordinate[1] = 0.0;
      if (m_eMesh->get_spatial_dim() > 2)
        {
          coordinate[2] = 0.0;
        }
      else
        {
          coordinate[2] = 1.0;
          return;
        }
      if (!m_meshGeometry) return;

      stk::mesh::Entity* element_ptr = reinterpret_cast<stk::mesh::Entity*>(entity_handle);

      stk::mesh::PairIterRelation nodes = element_ptr->relations(m_eMesh->node_rank());
      double nodes_size = nodes.size();

      for (unsigned inode=0; inode < nodes.size(); inode++)
        {
          stk::mesh::Entity& node = *nodes[inode].entity();

          std::vector<double> normal(3,0.0);
          m_meshGeometry->normal_at(m_eMesh, &node, normal);
          coordinate[0] += normal[0]/nodes_size;
          coordinate[1] += normal[1]/nodes_size;
          coordinate[2] += normal[2]/nodes_size;
        }
    }
                          
    /**\brief evaluate surface normals
     *
     * Returns normals for a domain.
     *
     *\param handles       The domain evaluated is the one in which
     *                     this mesh entity is constrained.
     *\param coordinates   As input, a list of positions at which to
     *                     evaluate the domain.  As output, the resulting
     *                     domain normals.
     *\param count         The length of the coordinates array.
     */
    void PerceptMesquiteMeshDomain::
    vertex_normal_at( const Mesquite::Mesh::VertexHandle* handles,
                      Mesquite::Vector3D coordinates[],
                      unsigned count,
                      Mesquite::MsqError& err ) const {
      for (unsigned i=0; i < count; i++)
        {
          vertex_normal_at(handles[i], coordinates[i]);
        }
    }
                            
    /**\brief evaluate closest point and normal
     *
     * Given a position in space, return the closest 
     * position in the domain and the domain normal
     * at that point.
     *
     *\param entity_handle Evaluate the subset of the domain contianing
     *                     this entity
     *\param position      Input position for which to evaluate
     *\param closest       Closest position in the domain.
     *\param normal        Domain normal at the location of 'closest'
     */
    void PerceptMesquiteMeshDomain::
    closest_point( Mesquite::Mesh::VertexHandle handle,
                   const Mesquite::Vector3D& position,
                   Mesquite::Vector3D& closest,
                   Mesquite::Vector3D& normal,
                   Mesquite::MsqError& err ) const {

      // this could be more efficient if the MeshGeometry interface supported an
      // auxiliary coordinate (e.g. the "position" arg)
      stk::mesh::Entity* node_ptr = reinterpret_cast<stk::mesh::Entity *>(handle);
      stk::mesh::FieldBase* field = m_eMesh->get_coordinates_field();
      double *f_data = PerceptMesh::field_data(field, *node_ptr);
      if (!m_meshGeometry)
        {
          for (int isd=0; isd < m_eMesh->get_spatial_dim(); isd++)
            {
              closest[isd] = f_data[isd];
              normal[isd]=0;  // FIXME
            }
          normal[2]=1.0;
          return;
        }

      // save coordinates, set to "position", project, copy to closest, set back to saved
      double save_coords[3] = {f_data[0], f_data[1], 0};
      if (m_eMesh->get_spatial_dim() > 2) save_coords[2] = f_data[2];
      f_data[0] = position[0];
      f_data[1] = position[1];
      if (m_eMesh->get_spatial_dim() > 2) 
        f_data[2] = position[2];

      static std::vector<stk::mesh::Entity *> nodes(1);
      nodes[0] = node_ptr;
      m_meshGeometry->snap_points_to_geometry(m_eMesh, nodes);

      closest[0] = f_data[0];
      closest[1] = f_data[1];
      if (m_eMesh->get_spatial_dim() > 2) 
        closest[2] = f_data[2];
      
      // FIXME srk
      normal[0] = 0;
      normal[1] = 0;
      normal[2] = 1;
      if (m_eMesh->get_spatial_dim() > 2) 
        {
          std::vector<double> norm(3,0.0);
          m_meshGeometry->normal_at(m_eMesh, node_ptr, norm);
          normal[0] = norm[0];
          normal[1] = norm[1];
          normal[2] = norm[2];
        }

      if (0)
        {
          std::cout << "tmp closest_point: orig= " 
                    << save_coords[0] << " "
                    << save_coords[1] << " "
                    << save_coords[2] << " "
                    << " new= " 
                    << f_data[0] << " "
                    << f_data[1] << " "
                    << f_data[2] 
                    << std::endl;
        }

      f_data[0] = save_coords[0];
      f_data[1] = save_coords[1];
      if (m_eMesh->get_spatial_dim() > 2) 
        f_data[2] = save_coords[2];
    }
                                
    int PerceptMesquiteMeshDomain::
    classify_node(stk::mesh::Entity& node, size_t& curveOrSurfaceEvaluator) const
    {
      int dof =0;
      if (m_meshGeometry)
        dof = m_meshGeometry->classify_node(node, curveOrSurfaceEvaluator);
      else
        dof = m_eMesh->get_spatial_dim();
      return dof;
    }

    /**\brief Get degrees of freedom in vertex movement.
     *
     * Given a vertex, return how the domain constrains the
     * location of that vertex as the number of degrees of
     * freedom in the motion of the vertex.  If the domain
     * is a geometric domain, the degrees of freedom for a
     * vertex is the dimension of the geometric entity the
     * vertex is constrained to lie on (e.g. point = 0, curve = 1,
     * surface = 2, volume = 3.)
     */
    void PerceptMesquiteMeshDomain::
    domain_DoF( const Mesquite::Mesh::EntityHandle* handle_array,
                unsigned short* dof_array,
                size_t num_handles,
                Mesquite::MsqError& err ) const {

      size_t curveOrSurfaceEvaluator;

      for (size_t i = 0; i < num_handles; i++)
        {
          stk::mesh::Entity* node_ptr = reinterpret_cast<stk::mesh::Entity *>(handle_array[i]);
          int dof = classify_node(*node_ptr, curveOrSurfaceEvaluator);

          if (dof < 0)
            {
              PRINT_ERROR("dof < 0");
              MSQ_SETERR(err)("PerceptMesquiteMeshDomain::domain_DoF classify returned -1", Mesquite::MsqError::INVALID_STATE);
              return;
            }
          dof_array[i] = (unsigned short)dof;
        }
    }
                             

    //static bool myPerceptMesquiteMesh_cpp = true;


  } // namespace percept
} // namespace stk

#endif // STK_BUILT_IN_SIERRA
