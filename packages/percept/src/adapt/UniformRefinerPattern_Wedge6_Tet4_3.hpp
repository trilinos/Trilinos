// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Wedge6_Tet4_3_hpp
#define adapt_UniformRefinerPattern_Wedge6_Tet4_3_hpp

#include "UniformRefinerPattern.hpp"
#include "UniformRefinerPattern_Quad4_Tri3_2.hpp"

  namespace percept {

    template <>
    class UniformRefinerPattern<shards::Wedge<6>, shards::Tetrahedron<4>, 3 > : public URP<shards::Wedge<6>, shards::Tetrahedron<4> >
    {
    private:
      UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 2 > *m_face_breaker;

      // permutations from rotating the element so that the first (0) node 
      // has the min ID. The first index is the index of the min node ID, 
      // the second index is the new index for each node after rotation.
      const unsigned min_vertex_renumber[6][6]={{0,1,2,3,4,5},
                                                {1,2,0,4,5,3},
                                                {2,0,1,5,3,4},
                                                {3,5,4,0,2,1},
                                                {4,3,5,1,0,2},
                                                {5,4,3,2,1,0}};
      
    public:

      virtual bool edgeMarkIsEnough() override { return false; }

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Wedge<6>, shards::Tetrahedron<4> >(eMesh)
      {
        EXCEPTWATCH;

        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, false);

        m_face_breaker = new UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 2 > (eMesh, block_names);
      }

      ~UniformRefinerPattern()
      {
        if (m_face_breaker) delete m_face_breaker;
      }

      virtual void doBreak() override {}

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
      {
        needed_entities.resize(1);
        needed_entities[0].first = stk::topology::ELEMENT_RANK;
        setToOne(needed_entities);
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
      {
        EXCEPTWATCH;
        bp.resize(2);
        bp[0] = this;
        bp[1] = m_face_breaker;
      }

      virtual unsigned getNumNewElemPerElem() override { return 3u; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& /*nodeRegistry*/,
                        stk::mesh::Entity element,  NewSubEntityNodesType& /*new_sub_entity_nodes*/, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        EXCEPTWATCH;
        typedef std::array<stk::mesh::EntityId, 4> tet_tuple_type;
        vector<tet_tuple_type> new_elements(3);
        
        stk::mesh::EntityId element_globalIds[6] = {0,0,0,0,0,0};
        const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK );
        
        stk::mesh::EntityId minVal = m_eMesh.identifier(elem_nodes[0].entity());
        int indxMinVal = 0;

        for (int inode=0; inode < 6; inode++) {

          stk::mesh::Entity node = elem_nodes[inode].entity();
          element_globalIds[inode] = m_eMesh.identifier(node);

          if (element_globalIds[inode] < minVal) {
            minVal = element_globalIds[inode];
            indxMinVal = inode;
          }
        }
        
        unsigned opposite_quad_face_nodes[4] = {min_vertex_renumber[indxMinVal][1],
                                                min_vertex_renumber[indxMinVal][2],
                                                min_vertex_renumber[indxMinVal][5],
                                                min_vertex_renumber[indxMinVal][4]};
        
        stk::mesh::EntityId minValFace = element_globalIds[opposite_quad_face_nodes[0]];
        int indxMinValFace = 0;

        for (int inode=1; inode < 4; inode++) {

          const stk::mesh::EntityId next_id = element_globalIds[opposite_quad_face_nodes[inode]];
          
          if (next_id < minValFace) {
            minValFace = next_id;
            indxMinValFace = inode;
          }
        }

        // tet connected to min ID vertex but not opposite quad faces
        new_elements[0] = {element_globalIds[min_vertex_renumber[indxMinVal][0]],
                           element_globalIds[min_vertex_renumber[indxMinVal][4]],
                           element_globalIds[min_vertex_renumber[indxMinVal][5]],
                           element_globalIds[min_vertex_renumber[indxMinVal][3]]};
        
        // tets formed from splitting opposite quad faces using its min ID vertex
        new_elements[1] = {element_globalIds[indxMinVal],
                           element_globalIds[opposite_quad_face_nodes[(0+indxMinValFace)%4]],
                           element_globalIds[opposite_quad_face_nodes[(1+indxMinValFace)%4]],
                           element_globalIds[opposite_quad_face_nodes[(2+indxMinValFace)%4]]};

        new_elements[2] = {element_globalIds[indxMinVal],
                           element_globalIds[opposite_quad_face_nodes[(0+indxMinValFace)%4]],
                           element_globalIds[opposite_quad_face_nodes[(2+indxMinValFace)%4]],
                           element_globalIds[opposite_quad_face_nodes[(3+indxMinValFace)%4]]};
                                         
        for (unsigned ielem=0; ielem < new_elements.size(); ielem++)
          {
            stk::mesh::Entity newElement = *element_pool;

            if (proc_rank_field)
              {
                double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
                fdata[0] = double(eMesh.owner_rank(newElement));
              }

            change_entity_parts(eMesh, element, newElement);

            unsigned nchild = new_elements.size();

            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(new_elements[ielem][0]), 0);
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(new_elements[ielem][1]), 1);
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(new_elements[ielem][2]), 2);
            eMesh.get_bulk_data()->declare_relation(newElement, eMesh.createOrGetNode(new_elements[ielem][3]), 3);

            set_parent_child_relations(eMesh, element, newElement, *ft_element_pool, ielem, &nchild);

            std::vector<stk::mesh::Entity> elements(1,element);
            eMesh.prolongateElementFields( elements, newElement);

            ft_element_pool++;
            element_pool++;
          }
      }

    };

  }

#endif
