// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Pyramid5_Tet4_2_hpp
#define adapt_UniformRefinerPattern_Pyramid5_Tet4_2_hpp

#include "UniformRefinerPattern.hpp"
#include "UniformRefinerPattern_Quad4_Tri3_2.hpp"

  namespace percept {

  /**
   *
   *   Convention: input is the element's nodes and the marks on the edges.  Output is an array
   *     of "elements" defined as local id's of nodes forming those elements, where {0,..,4} represent
   *     the original vertices and {5..12} are the edges, {13} are faces and {14} is the centroid used
   *     to join to faces to discretize the pyramid
   *     (these are numberings used by Shards for the extra/quadratic nodes, and so we map these
   *      numbers back to the local edge/face/centroid using tables; for example, quadratic node 9
   *      corresponds to edge {0,4} see table "edges" and edge_map_rev below)
   *
   *
   *   typedef
   *     MakeTypeList< IndexList< 0 , 1 ,   5 > ,
   *                   IndexList< 1 , 2 ,   6 > ,
   *                   IndexList< 2 , 3 ,   7 > ,
   *                   IndexList< 3 , 0 ,   8 > ,
   *                   IndexList< 0 , 4 ,   9 > ,
   *                   IndexList< 1 , 4 ,  10 > ,
   *                   IndexList< 2 , 4 ,  11 > ,
   *                   IndexList< 3 , 4 ,  12 > >::type
   *     PyramidEdgeNodeMap ;
   *
   *   typedef
   *     MakeTypeList< IndexList< 0, 1, 4,     5, 10,  9 > ,
   *                   IndexList< 1, 2, 4,     6, 11, 10 > ,
   *                   IndexList< 2, 3, 4,     7, 12, 11 > ,
   *                   IndexList< 3, 0, 4,     8,  9, 12 > ,
   *                   IndexList< 0, 3, 2, 1,  8,  7,  6,  5,  13 > >::type
   *     PyramidFaceNodeMap ;
   *
   *
   *   Linear 5-Node Pyramid node locations.
   *
   *  4  - top vertex of pyramid
   *  13 - centroid of quad face
   *
   *
   *                       7
   *            3 o--------o---------o 2
   *             / \ 12          /  /
   *            /   o       o 11   /
   *           /     \    /       /
   *          /      o 4         /
   *         /     /  \         /
   *      8 o  9 /  o \       o 6
   *       /   o   13  o 10   /
   *      /  /          \    /
   *     / /            \   /
   *    //               \ /
   *   o--------o---------o
   *  0         5          1
   *
   *
   *
   */

    template <>
    class UniformRefinerPattern<shards::Pyramid<5>, shards::Tetrahedron<4>, 2 > : public URP<shards::Pyramid<5>, shards::Tetrahedron<4> >
    {
    private:
      UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, 2 > *m_face_breaker;

    public:

      virtual bool edgeMarkIsEnough() override { return false; }

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Pyramid<5>, shards::Tetrahedron<4> >(eMesh)
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

      /// NOTE: we create additional un-used elements if the Pyramid5 can be broken into 6 tets
      virtual unsigned getNumNewElemPerElem() override { return 2u; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& /*nodeRegistry*/,
                        stk::mesh::Entity element,  NewSubEntityNodesType& /*new_sub_entity_nodes*/, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0) override
      {
        EXCEPTWATCH;
        typedef std::array<stk::mesh::EntityId, 4> tet_tuple_type;
        vector<tet_tuple_type> new_elements(2);
        
        static unsigned element_globalIds[5] = {0,0,0,0,0};
        const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK );
        
        for (int inode=0; inode < 5; inode++)
          {
            stk::mesh::Entity node = elem_nodes[inode].entity();
            element_globalIds[inode] = m_eMesh.identifier(node);
          }
        
        unsigned globalIqf  = element_globalIds[0];
        unsigned minVal     = globalIqf;
        unsigned indxMinVal = 0;
        for (unsigned iFaceNodeOrd=1; iFaceNodeOrd < 4; iFaceNodeOrd++)
          {
            globalIqf = element_globalIds[iFaceNodeOrd];
            if (globalIqf < minVal)
              {
                minVal = globalIqf;
                indxMinVal = iFaceNodeOrd;
              }
          }
        unsigned istart = indxMinVal;

        new_elements[0] = {element_globalIds[(0+istart)%4],
                           element_globalIds[(1+istart)%4],
                           element_globalIds[(2+istart)%4],
                           element_globalIds[4]};
                                         
        new_elements[1] = {element_globalIds[(0+istart)%4],
                           element_globalIds[(2+istart)%4],
                           element_globalIds[(3+istart)%4],
                           element_globalIds[4]};
                                         
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
