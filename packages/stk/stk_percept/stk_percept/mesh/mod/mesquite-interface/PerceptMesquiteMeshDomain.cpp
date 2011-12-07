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

#ifdef STK_BUILT_IN_SIERRA

#include "PerceptMesquiteMeshDomain.hpp"
//#include <mesquite/MsqVertex.hpp>
#include <MsqVertex.hpp>

#include <map>
#include <algorithm>

namespace stk {
  namespace percept {

#if 0

#define PRINT_ERROR(a) do { std::cout << "PerceptMesquite::Mesh::ERROR: " << a << std::endl; } while (0)

    //! Modifies "coordinate" so that it lies on the
    //! domain to which "entity_handle" is constrained.
    //! The handle determines the domain.  The coordinate
    //! is the proposed new position on that domain.
    void PerceptMesquiteMeshDomain::
    snap_to(Mesquite::Mesh::VertexHandle entity_handle,
            Mesquite::Vector3D &coordinate) const {

      stk::mesh::Entity* node_ptr = reinterpret_cast<stk::mesh::Entity *>(entity_handle);
      stk::mesh::FieldBase* field = m_eMesh->getCoordinatesField();
      double *f_data = PerceptMesh::field_data(field, *node_ptr);

      static std::vector<stk::mesh::Entity *> nodes(1);
      nodes[0] = node_ptr;
      m_meshGeometry->snap_points_to_geometry(m_eMesh, nodes);
      coordinate[0] = f_data[0];
      coordinate[1] = f_data[1];
      if (m_eMesh->getSpatialDim() > 2)
        coordinate[2] = f_data[2];
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
      if (m_eMesh->getSpatialDim() > 2)
        coordinate[2] = 0.0;
    }
    void PerceptMesquiteMeshDomain::
    element_normal_at(Mesquite::Mesh::ElementHandle entity_handle,
                      Mesquite::Vector3D &coordinate) const {
      // FIXME srk
      coordinate[0] = 0.0;
      coordinate[1] = 0.0;
      if (m_eMesh->getSpatialDim() > 2)
        coordinate[2] = 0.0;
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
      stk::mesh::FieldBase* field = m_eMesh->getCoordinatesField();
      double *f_data = PerceptMesh::field_data(field, *node_ptr);
      // save coordinates, set to "position", project, copy to closest, set back to saved
      double save_coords[3] = {f_data[0], f_data[1], 0};
      if (m_eMesh->getSpatialDim() > 2) save_coords[2] = f_data[2];
      f_data[0] = position[0];
      f_data[1] = position[1];
      if (m_eMesh->getSpatialDim() > 2) 
        f_data[2] = position[2];

      static std::vector<stk::mesh::Entity *> nodes(1);
      nodes[0] = node_ptr;
      m_meshGeometry->snap_points_to_geometry(m_eMesh, nodes);
      closest[0] = f_data[0];
      closest[1] = f_data[1];
      if (m_eMesh->getSpatialDim() > 2) 
        closest[2] = f_data[2];
      
      // FIXME srk
      normal[0] = 0;
      normal[1] = 0;
      if (m_eMesh->getSpatialDim() > 2) 
        normal[2] = 0;

      f_data[0] = save_coords[0];
      f_data[1] = save_coords[1];
      if (m_eMesh->getSpatialDim() > 2) 
        f_data[2] = save_coords[2];
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

      static std::vector<size_t> curveEvaluators;
      static std::vector<size_t> surfEvaluators;

      for (size_t i = 0; i < num_handles; i++)
        {
          stk::mesh::Entity* node_ptr = reinterpret_cast<stk::mesh::Entity *>(handle_array[i]);
          curveEvaluators.resize(0);
          surfEvaluators.resize(0);
          int dof = m_meshGeometry->classify_node(node_ptr, curveEvaluators, surfEvaluators);
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

#endif

  } // namespace percept
} // namespace stk

#endif // STK_BUILT_IN_SIERRA
