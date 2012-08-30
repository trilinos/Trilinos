//-------------------------------------------------------------------------
// Filename      : PerceptMesquiteMeshDomain.hpp
//
// Purpose       : mesh domain (geometry) interface for using Mesquite 
//
// Description   : subclass of Mesquite::Mesh::Domain
//
// Creator       : Steve Kennon, derived from Steve Owen's work
//
// Creation Date : Dec 2011
//
// Owner         : Steve Kennon
//-------------------------------------------------------------------------
#ifndef PERCEPT_MESQUITE_MESH_DOMAIN_HPP
#define PERCEPT_MESQUITE_MESH_DOMAIN_HPP

#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)


#include <MeshInterface.hpp>
#include <MsqError.hpp>
#include <MsqGeomPrim.hpp>

#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/mesh/geometry/kernel/MeshGeometry.hpp>
#include <map>

class MeshGeometry;

namespace stk {
  namespace percept {

    class PerceptMesquiteMeshDomain : public Mesquite::MeshDomain
    {
      PerceptMesh *m_eMesh;
      MeshGeometry *m_meshGeometry;

      //       stk::mesh::Selector *m_boundarySelector;
      //       std::map<stk::mesh::Entity *, std::pair<stk::mesh::EntityId, unsigned char> > m_mesquiteNodeDataMap;

    public:

      PerceptMesquiteMeshDomain(PerceptMesh *eMesh, MeshGeometry *meshGeometry=0) : m_eMesh(eMesh), m_meshGeometry(meshGeometry) 
      {
        if (m_meshGeometry) 
          m_meshGeometry->m_cache_classify_bucket_is_active = true;
      }

      //       void init(PerceptMesh *eMesh);
      //       int setup();
      //       bool is_on_my_patch_boundary(stk::mesh::Entity *node_ptr);
      //       void clean_out();

    public:
      virtual ~PerceptMesquiteMeshDomain() {}
  
      
      //! Modifies "coordinate" so that it lies on the
      //! domain to which "entity_handle" is constrained.
      //! The handle determines the domain.  The coordinate
      //! is the proposed new position on that domain.
      virtual void snap_to(Mesquite::Mesh::VertexHandle entity_handle,
                           Mesquite::Vector3D &coordinate) const ;
    
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
      virtual void vertex_normal_at(Mesquite::Mesh::VertexHandle entity_handle,
                                    Mesquite::Vector3D &coordinate) const ;
      virtual void element_normal_at(Mesquite::Mesh::ElementHandle entity_handle,
                                     Mesquite::Vector3D &coordinate) const ;
                          
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
      virtual void vertex_normal_at( const Mesquite::Mesh::VertexHandle* handles,
                                     Mesquite::Vector3D coordinates[],
                                     unsigned count,
                                     Mesquite::MsqError& err ) const ;
                            
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
      virtual void closest_point( Mesquite::Mesh::VertexHandle handle,
                                  const Mesquite::Vector3D& position,
                                  Mesquite::Vector3D& closest,
                                  Mesquite::Vector3D& normal,
                                  Mesquite::MsqError& err ) const ;
                                
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
      virtual void domain_DoF( const Mesquite::Mesh::EntityHandle* handle_array,
                               unsigned short* dof_array,
                               size_t num_handles,
                               Mesquite::MsqError& err ) const ;
                             
      int classify_node(stk::mesh::Entity& node, size_t& curveOrSurfaceEvaluator) const;

    };


    //static bool myPerceptMesquiteMeshDomain_hpp = true; 

  } // namespace percept
} // namespace stk

#endif // STK_BUILT_IN_SIERRA
#endif //has file been included
