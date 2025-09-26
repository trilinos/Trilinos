// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinerPattern_Quad4_Het_N_hpp
#define adapt_RefinerPattern_Quad4_Het_N_hpp


#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "RefinerPattern_Line2_Line2_N.hpp"
#include <adapt/TriangulateQuad.hpp>

namespace percept {

#define PRINT_PARTS 0

  extern bool s_do_transition_break;

  // Some explanation: Quad local (non-hanging-node) refinement pattern creates a heterogeneous mesh, so we create two
  //   sub-patterns to deal with the two different resulting topologies (quad and tri).  A third
  //   (parent) pattern is created to refer to the two sub-patterns, akin to URP_Heterogeneous_3D.

  // this pattern is only intended to be used for "transition" element breaking, that is, all hanging
  // nodes are removed by filling the interior of the quad with sub-quads or sub-tris

  // It could also be used for purely local refinement where some tri's get created then propagated
  //================================================================================================================================================================



  //================================================================================================================================================================
  // Just for ensuring triangles are processed by the marker - this pattern does nothing, it doesn't break triangles...
  //================================================================================================================================================================

  struct TriTempPartialNoBreak {};

  // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
  template <>
  class RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1, TriTempPartialNoBreak > : public URP<shards::Triangle<3>, shards::Triangle<3>  >
  {
  public:
    RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Triangle<3>, shards::Triangle<3>  >(eMesh)
    {
      m_primaryEntityRank = m_eMesh.face_rank();
      if (m_eMesh.get_spatial_dim() == 2)
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

      bool sameTopology = true;
      setNeededParts(eMesh, block_names, sameTopology);
      Elem::StdMeshObjTopologies::bootstrap();

    }

    ~RefinerPattern()
    {
    }

    void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh ) override
    {
      bp.resize(1);

      if (eMesh.get_spatial_dim() == 2)
        {
          bp[0] = this;
        }
      else
        bp.resize(0);
    }

    virtual void doBreak() override {}
    void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
    {
      needed_entities.resize(1);
      needed_entities[0].first = m_eMesh.edge_rank();
      setToOne(needed_entities);
    }

    // FIXME - for now, create more than we need (to fix this right we need a change to the Refiner.cpp interface)
    virtual unsigned getNumNewElemPerElem() override { return 4; }

    void
    createNewElements(percept::PerceptMesh& /*eMesh*/, NodeRegistry& /*nodeRegistry*/,
                      stk::mesh::Entity /*element*/,  NewSubEntityNodesType& /*new_sub_entity_nodes*/, vector<stk::mesh::Entity>::iterator& /*element_pool*/,
                      vector<stk::mesh::Entity>::iterator& /*ft_element_pool*/,
                      stk::mesh::FieldBase */*proc_rank_field*/=0) override
    {
    }

  };


  //================================================================================================================================================================
  // Heterogeneous child patterns
  //================================================================================================================================================================

  struct QuadTriPartial {};

  // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
  template <>
  class RefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>, -1, QuadTriPartial > : public URP<shards::Quadrilateral<4>, shards::Triangle<3>  >
  {
    RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > * m_edge_breaker;
    TransitionElementType *m_transition_element_field;
    bool m_transition_element_field_set;
    bool m_use_only_tris;
    bool m_avoid_centroid_node;

  public:


    RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType(), bool use_only_tris=true, bool avoid_centroid_node=false) 
      :  URP<shards::Quadrilateral<4>, shards::Triangle<3>  >(eMesh),
         m_edge_breaker(0),
         m_transition_element_field(0), m_transition_element_field_set(false)
      ,m_use_only_tris(use_only_tris), m_avoid_centroid_node(avoid_centroid_node)
    {
      m_primaryEntityRank = m_eMesh.face_rank();
      if (m_eMesh.get_spatial_dim() == 2)
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

      bool sameTopology = false;
      setNeededParts(eMesh, block_names, sameTopology);
      Elem::StdMeshObjTopologies::bootstrap();

      if (m_eMesh.get_spatial_dim() == 2)
        {
          m_edge_breaker =  new RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > (eMesh, block_names) ;
        }

    }

    ~RefinerPattern()
    {
      if (m_edge_breaker) delete m_edge_breaker;
    }

    void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh ) override
    {
      EXCEPTWATCH;
      bp.resize(2);

      if (eMesh.get_spatial_dim() == 2)
        {
          bp[0] = this;
          bp[1] = m_edge_breaker;
        }
      else
        {
          bp.resize(0);
        }
    }

    virtual void doBreak() override {}
    void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
    {
      needed_entities.resize(2);
      needed_entities[0].first = m_eMesh.edge_rank();
      needed_entities[1].first = (m_eMesh.get_spatial_dim() == 2 ? stk::topology::ELEMENT_RANK : m_eMesh.face_rank());
      setToOne(needed_entities);
    }

    // FIXME - for now, create more than we need (to fix this right we need a change to the Refiner.cpp interface)
    virtual unsigned getNumNewElemPerElem() override { return 8; }

    void
    createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                      stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                      vector<stk::mesh::Entity>::iterator& ft_element_pool,
                      stk::mesh::FieldBase *proc_rank_field=0) override
    {
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
      static vector<quad_to_tri_tuple_type> elems(8);
      static vector<quad_to_quad_tuple_type> elems_quad(8);
      static vector<quad_to_tri_tuple_type_local> elems_local(8);
      static vector<quad_to_quad_tuple_type_local> elems_quad_local(8);
      unsigned num_new_elems=0;

      shards::CellTopology cell_topo(cell_topo_data);
      const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::topology::NODE_RANK);
      //CoordinatesFieldType* coordField = eMesh.get_coordinates_field();

      if (!m_transition_element_field_set)
        {
          m_transition_element_field_set = true;
          m_transition_element_field = eMesh.get_transition_element_field();
          if(m_transition_element_field && m_transition_element_field->name() != "transition_element") m_transition_element_field = nullptr;
        }

      std::vector<stk::mesh::Part*> add_parts;
      std::vector<stk::mesh::Part*> remove_parts;
      add_parts = m_toParts;

      unsigned edge_marks[4] = {0,0,0,0};
      unsigned num_edges_marked=0;
      for (int iedge = 0; iedge < 4; iedge++)
        {
          unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
          if (num_nodes_on_edge)
            {
              edge_marks[iedge] = 1;
              ++num_edges_marked;
            }
        }

      unsigned num_faces_marked = 0;
      stk::mesh::EntityRank rank = (m_eMesh.get_spatial_dim() == 2 ? stk::topology::ELEMENT_RANK : m_eMesh.face_rank());
      if ( new_sub_entity_nodes[rank].size() )
        {
          int iface = 0;
          num_faces_marked = new_sub_entity_nodes[rank][iface].size();
        }

      if (0)
        {
          std::vector<unsigned> em(&edge_marks[0], &edge_marks[0]+4);
          std::cout << "createNewElements:: element= " << m_eMesh.identifier(element)
                    << " num_edges_marked= " << num_edges_marked << " num_faces_marked= " << num_faces_marked
                    << " em= " << em
                    << std::endl;
        }

      if (num_edges_marked == 0)
        return;

      bool avoid_centroid_node = m_eMesh.get_spatial_dim() == 3 || m_avoid_centroid_node;
      bool use_only_tris = m_use_only_tris;
      if (m_eMesh.get_spatial_dim() == 3) use_only_tris = false;
      TriangulateQuad tq(use_only_tris, avoid_centroid_node);
      tq.triangulate_quad_face( edge_marks, elems_local, elems_quad_local);

#define Q_CENTROID_N_EV (useAltNode ? m_eMesh.identifier(new_nodes[0]) : NN(m_primaryEntityRank, 0) )

#define Q_CV_EV(i) ( i < 4 ? VERT_N(i) : (i == 8 ? Q_CENTROID_N_EV : EDGE_N(i-4) ) )

      if (0)
        std::cout << "elem= " << m_eMesh.identifier(element)
          //<< " Q_CENTROID_N_EV = " << Q_CENTROID_N_EV
                  << " new_sub_entity_nodes[m_primaryEntityRank].size()= "
                  << new_sub_entity_nodes[m_primaryEntityRank].size()
                  << " " << (new_sub_entity_nodes[m_primaryEntityRank].size() ?
                       (new_sub_entity_nodes[m_primaryEntityRank][0].size() ?
                        new_sub_entity_nodes[m_primaryEntityRank][0][0] : 123456 ) : 1234567)
                  << std::endl;

      std::vector<stk::mesh::Entity> new_nodes;
      bool useAltNode = false;
      if (!m_avoid_centroid_node && !new_sub_entity_nodes[m_primaryEntityRank][0].size())
        {
          //throw std::runtime_error("no centroid node availabled");
           // std::cout << " bad element: " << std::endl;
           // m_eMesh.print(element);
           // m_eMesh.dump_vtk("bad.vtk");
          if (m_eMesh.get_spatial_dim() == 2)
            {
              if (!m_eMesh.getEntitiesUsingIdServerNewNodes(1, new_nodes))
                {
                  std::cout << "P[" << m_eMesh.get_rank() << "] node_rank entity pool deplenished in RefinerPattern_Quad4_Het_N"
                            << m_eMesh.demangled_stacktrace() << std::endl;
                  throw std::logic_error("node rank entity pool deplenished");
                }
              useAltNode = true;
            }
        }

      if (useAltNode)
        {
          nodeRegistry.forceInMap(new_nodes, NodeRegistry::NR_MARK, element, stk::topology::ELEMENT_RANK, 0u);
          nodeRegistry.prolongateCoords(element, stk::topology::ELEMENT_RANK, 0u);
        }

      num_new_elems = elems_local.size();
      elems.resize(num_new_elems);
      for (unsigned ielem=0; ielem < num_new_elems; ielem++)
        {
          if (0)
            std::cout << "ielem= " << ielem << " num_new_elems= " << num_new_elems << " elems_local= " << elems_local[ielem]
                      << " new_sub_entity_nodes[face].size= " << new_sub_entity_nodes[m_primaryEntityRank].size()
                      << " new_sub_entity_nodes[face][0].size= " << new_sub_entity_nodes[m_primaryEntityRank][0].size()
                      << "\n elem= " << elems[ielem]
                      << std::endl;
          elems[ielem] = { Q_CV_EV(elems_local[ielem][0] ), Q_CV_EV(elems_local[ielem][1] ), Q_CV_EV(elems_local[ielem][2] ) };
        }

      if (0)
        std::cout << "tmp RefinerPattern_QuadTriPartial::num_edges_marked= " << num_edges_marked
                  << " num_new_elems= " << num_new_elems
                  << std::endl;

      bool use_declare_element_side = UniformRefinerPatternBase::USE_DECLARE_ELEMENT_SIDE &&  m_primaryEntityRank == eMesh.side_rank();

      for (unsigned ielem=0; ielem < elems.size(); ielem++)
        {
          stk::mesh::Entity newElement = stk::mesh::Entity();
          if (!use_declare_element_side)
            newElement = *element_pool;


          // 3 nodes of the new tris
          stk::mesh::Entity n0 = eMesh.createOrGetNode(elems[ielem][0]);
          stk::mesh::Entity n1 = eMesh.createOrGetNode(elems[ielem][1]);
          stk::mesh::Entity n2 = eMesh.createOrGetNode(elems[ielem][2]);

          stk::mesh::Entity nodes[3] = {n0, n1, n2};

          create_side_element(eMesh, use_declare_element_side, nodes, 3, newElement);

          if (m_transition_element_field)
            {
              int *transition_element = 0;
              transition_element = stk::mesh::field_data( *m_transition_element_field , newElement );
              if (num_edges_marked == 4)
                {
                  transition_element[0] = 0;
                  // FIXME error???
                  //throw std::runtime_error("logic error in RefinerPattern_Quad4_Het_N");
                }
              else
                {
                  transition_element[0] = 1;
                }
            }

          if (proc_rank_field && eMesh.get_spatial_dim() == 2)
            {
              double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
              if (fdata)
                fdata[0] = double(eMesh.owner_rank(newElement));
            }

          change_entity_parts(eMesh, element, newElement);

          unsigned nchild = eMesh.numChildren(element);
          //set_parent_child_relations(eMesh, element, newElement, ielem);
          set_parent_child_relations(eMesh, element, newElement,  *ft_element_pool, nchild);

          if (0) std::cout << "tmp createNewElements: isParentElement= " << m_eMesh.isParentElement(element, false) << std::endl;
          std::vector<stk::mesh::Entity> elements(1,element);
          eMesh.prolongateElementFields( elements, newElement);

          ++ft_element_pool;
          element_pool++;

        }
    }

  };

  //================================================================================================================================================================
  // Heterogeneous child patterns - quad to quad
  //================================================================================================================================================================

  // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
  struct QuadQuadPartial {};

  template <>
  class RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1, QuadQuadPartial > : public URP<shards::Quadrilateral<4>, shards::Quadrilateral<4>  >
  {
    RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > * m_edge_breaker;
    TransitionElementType *m_transition_element_field;
    bool m_transition_element_field_set;
    bool m_use_only_tris;
    bool m_avoid_centroid_node;

  public:

    RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType(), bool use_only_tris=true, bool avoid_centroid_node=false) 
      :  URP<shards::Quadrilateral<4>, shards::Quadrilateral<4>  >(eMesh),
         m_edge_breaker(0),
         m_transition_element_field(0), m_transition_element_field_set(false)
      ,m_use_only_tris(use_only_tris), m_avoid_centroid_node(avoid_centroid_node)
    {
      m_primaryEntityRank = m_eMesh.face_rank();
      if (m_eMesh.get_spatial_dim() == 2)
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

      bool sameTopology = true;
      setNeededParts(eMesh, block_names, sameTopology);
      Elem::StdMeshObjTopologies::bootstrap();

      if (m_eMesh.get_spatial_dim() == 2)
        {
          m_edge_breaker =  new RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > (eMesh, block_names) ;
        }
    }

    ~RefinerPattern()
    {
      if (m_edge_breaker) delete m_edge_breaker;
    }

    void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh ) override
    {
      EXCEPTWATCH;
      bp.resize(2);
      if (eMesh.get_spatial_dim() == 2)
        {
          bp[0] = this;
          bp[1] = m_edge_breaker;
        }
      else
        {
          bp.resize(0);
        }
    }

    virtual void doBreak() override {}
    void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
    {
      needed_entities.resize(2);
      needed_entities[0].first = m_eMesh.edge_rank();
      needed_entities[1].first = (m_eMesh.get_spatial_dim() == 2 ? stk::topology::ELEMENT_RANK : m_eMesh.face_rank());
      setToOne(needed_entities);
    }

    // FIXME - for now, create more than we need (to fix this right we need a change to the Refiner.cpp interface)
    virtual unsigned getNumNewElemPerElem() override { return 3; }

    void
    createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                      stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                      vector<stk::mesh::Entity>::iterator& ft_element_pool,
                      stk::mesh::FieldBase *proc_rank_field=0) override
    {
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
      static vector<quad_to_tri_tuple_type> elems_tri(8);
      static vector<quad_to_quad_tuple_type> elems(8);
      static vector<quad_to_tri_tuple_type_local> elems_tri_local(8);
      static vector<quad_to_quad_tuple_type_local> elems_local(8);
      unsigned num_new_elems=0;

      shards::CellTopology cell_topo(cell_topo_data);
      const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::topology::NODE_RANK);
      //CoordinatesFieldType* coordField = eMesh.get_coordinates_field();

      if (!m_transition_element_field_set)
        {
          m_transition_element_field_set = true;
          m_transition_element_field = eMesh.get_transition_element_field();
          if(m_transition_element_field->name() != "transition_element") m_transition_element_field = nullptr;
        }

      std::vector<stk::mesh::Part*> add_parts;
      std::vector<stk::mesh::Part*> remove_parts;
      add_parts = m_toParts;

      unsigned edge_marks[4] = {0,0,0,0};
      unsigned num_edges_marked=0;
      for (int iedge = 0; iedge < 4; iedge++)
        {
          unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
          if (num_nodes_on_edge)
            {
              edge_marks[iedge] = 1;
              ++num_edges_marked;
            }
        }

      unsigned num_faces_marked = 0;
      stk::mesh::EntityRank rank = (m_eMesh.get_spatial_dim() == 2 ? stk::topology::ELEMENT_RANK : m_eMesh.face_rank());
      if ( new_sub_entity_nodes[rank].size() )
        {
          int iface = 0;
          num_faces_marked = new_sub_entity_nodes[rank][iface].size();
        }

      if (0)
        std::cout << "createNewElements:: element= " << m_eMesh.identifier(element)
                  << " num_edges_marked= " << num_edges_marked << " num_faces_marked= " << num_faces_marked << std::endl;

      if (num_edges_marked == 0)
        return;

      bool avoid_centroid_node = m_eMesh.get_spatial_dim() == 3 || m_avoid_centroid_node;
      bool use_only_tris = m_use_only_tris;
      if (m_eMesh.get_spatial_dim() == 3) use_only_tris = false;
      TriangulateQuad tq(use_only_tris, avoid_centroid_node);
      tq.triangulate_quad_face(edge_marks, elems_tri_local, elems_local);

      std::vector<stk::mesh::Entity> new_nodes;
      bool useAltNode = false;
      if (!m_avoid_centroid_node && !new_sub_entity_nodes[m_primaryEntityRank][0].size())
        {
          //throw std::runtime_error("no centroid node availabled");
           // std::cout << " bad element: " << std::endl;
           // m_eMesh.print(element);
           // m_eMesh.dump_vtk("bad.vtk");
          if (m_eMesh.get_spatial_dim() == 2)
            {
              if (!m_eMesh.getEntitiesUsingIdServerNewNodes( 1, new_nodes))
                {
                  std::cout << "P[" << m_eMesh.get_rank() << "] node_rank entity pool deplenished in RefinerPattern_Quad4_Het_N" << std::endl;
                  throw std::logic_error("node rank entity pool deplenished");
                }
              useAltNode = true;
            }
        }
      if (useAltNode)
        {
          nodeRegistry.forceInMap(new_nodes, NodeRegistry::NR_MARK, element, stk::topology::ELEMENT_RANK, 0u);
          nodeRegistry.prolongateCoords(element, stk::topology::ELEMENT_RANK, 0u);
        }

      num_new_elems = elems_local.size();
      if (0)
        std::cout << "tmp RefinerPattern_QuadQuadPartial::num_edges_marked= " << num_edges_marked
                  << " num_new_elems= " << num_new_elems
                  << std::endl;

      elems.resize(num_new_elems);
      for (unsigned ielem=0; ielem < num_new_elems; ielem++)
        {
          elems[ielem] = { Q_CV_EV(elems_local[ielem][0] ), Q_CV_EV(elems_local[ielem][1] ), Q_CV_EV(elems_local[ielem][2] ), Q_CV_EV(elems_local[ielem][3] ) };
        }

      bool use_declare_element_side = UniformRefinerPatternBase::USE_DECLARE_ELEMENT_SIDE &&  m_primaryEntityRank == eMesh.side_rank();

      for (unsigned ielem=0; ielem < elems.size(); ielem++)
        {
          stk::mesh::Entity newElement = stk::mesh::Entity();
          if (!use_declare_element_side)
            newElement = *element_pool;

          // 4 nodes of the new quads
          stk::mesh::Entity nodes[4] = {
            eMesh.createOrGetNode(elems[ielem][0]),
            eMesh.createOrGetNode(elems[ielem][1]),
            eMesh.createOrGetNode(elems[ielem][2]),
            eMesh.createOrGetNode(elems[ielem][3]) };

          create_side_element(eMesh, use_declare_element_side, nodes, 4, newElement);

          if (m_transition_element_field)
            {
              int *transition_element = 0;
              transition_element = stk::mesh::field_data( *m_transition_element_field , newElement );
              if (num_edges_marked == 4)
                {
                  transition_element[0] = 0;
                }
              else
                {
                  transition_element[0] = 1;
                }
            }

          if (proc_rank_field && eMesh.get_spatial_dim() == 2)
            {
              double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
              if (fdata)
                fdata[0] = double(eMesh.owner_rank(newElement));
            }

          change_entity_parts(eMesh, element, newElement);


          unsigned nchild = eMesh.numChildren(element);
          //set_parent_child_relations(eMesh, element, newElement, ielem);
          set_parent_child_relations(eMesh, element, newElement, *ft_element_pool, nchild);

          if (0) std::cout << "tmp createNewElements: isParentElement= " << m_eMesh.isParentElement(element, false) << std::endl;
          std::vector<stk::mesh::Entity> elements(1,element);
          eMesh.prolongateElementFields( elements, newElement);

          ft_element_pool++;
          if (!use_declare_element_side)
            element_pool++;

        }
    }

  };


  //================================================================================================================================================================
  // Parent pattern
  //================================================================================================================================================================

  struct QuadHet {};

  template<>
  class RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 8, QuadHet > : public UniformRefinerPatternBase
  {
    std::vector<UniformRefinerPatternBase *> m_bp;
  public:
    RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > * m_edge_breaker;
    std::vector<UniformRefinerPatternBase *> m_bp_exported;

  protected:

    percept::PerceptMesh& m_eMesh;

  public:

    RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType(), bool use_only_tris=true, bool avoid_centroid_node=false) :  m_eMesh(eMesh)
    {
      m_primaryEntityRank = eMesh.face_rank();
      if (m_eMesh.get_spatial_dim() == 2)
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

      Elem::StdMeshObjTopologies::bootstrap();

      m_bp.resize(0);

      m_bp.push_back(new  RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1, QuadQuadPartial > (eMesh, block_names, use_only_tris, avoid_centroid_node) );
      m_bp.push_back(new  RefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>,      -1, QuadTriPartial > (eMesh, block_names, use_only_tris, avoid_centroid_node) );

      m_bp.push_back(new  RefinerPattern<shards::Triangle<3>, shards::Triangle<3>,     -1, TriTempPartialNoBreak > (eMesh, block_names) );

      bool sameTopology = false;
      m_bp[0]->setNeededParts(eMesh, block_names, true); // don't force a new part for quads, even though it's really the same topology
      m_bp[1]->setNeededParts(eMesh, block_names, sameTopology);
      m_bp[2]->setNeededParts(eMesh, block_names, true);

      m_bp_exported.resize(0);
      m_bp_exported.push_back(m_bp[2]);

      if (PRINT_PARTS)
        {
          for (unsigned ibp=0; ibp < m_bp.size(); ibp++)
            {
              std::cout << "tmp Quad_Tri printParts this before=\n" ;
              printParts(m_bp[ibp], false);
            }
        }

      // repeat setNeededParts to catch newly created parts
      bool skipConvertedParts = false;
      m_bp[0]->setNeededParts(eMesh, block_names, true, skipConvertedParts); // don't force a new part for quads, even though it's really the same topology
      m_bp[1]->setNeededParts(eMesh, block_names, sameTopology, skipConvertedParts);
      m_bp[2]->setNeededParts(eMesh, block_names, true, skipConvertedParts);

      if (PRINT_PARTS)
        {
          for (unsigned ibp=0; ibp < m_bp.size(); ibp++)
            {
              std::cout << "tmp Quad_Tri printParts this after=\n" ;
              printParts(m_bp[ibp], false);
            }
        }

      if (1)
        {
          for (unsigned ibp=0; ibp < m_bp.size(); ibp++)
            {
              bool mergeParts = false;
              mergeOrAddParts(m_bp[ibp], this, mergeParts);
            }
        }

      if (PRINT_PARTS)
        {
          std::cout << "tmp Quad_Tri printParts this=\n" ;
          printParts(this, false);
        }

      if (m_eMesh.get_spatial_dim() == 2)
        m_edge_breaker =  new RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > (eMesh, block_names) ;
      else
        m_edge_breaker = 0;

    }

    ~RefinerPattern()
    {
      if (m_edge_breaker) delete m_edge_breaker;
      for (unsigned ibp=0; ibp < m_bp.size(); ibp++)
        {
          if (m_bp[ibp]) delete m_bp[ibp];
        }
    }

    void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh ) override
    {
      EXCEPTWATCH;
      bp.resize(0);

      if (eMesh.get_spatial_dim() == 2)
        {
          bp = m_bp;
          bp.push_back(m_edge_breaker);
        }
      else if (eMesh.get_spatial_dim() == 3)
        {
          bp.resize(0);
          // FIXME
          //             std::cout << "ERROR" ;
          //             exit(1);
        }

    }

    virtual void doBreak() override {
      throw std::runtime_error("shouldn't call RefinerPattern_Quad4_Het_N::doBreak()");
    }
    void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
    {
      needed_entities.resize(2);
      needed_entities[0].first = m_eMesh.edge_rank();
      needed_entities[1].first = (m_eMesh.get_spatial_dim() == 2 ? stk::topology::ELEMENT_RANK : m_eMesh.face_rank());
      setToOne(needed_entities);
    }

    virtual unsigned getNumNewElemPerElem() override { return 8; }

    virtual unsigned getFromTypeKey() override
    {
      return shards::Quadrilateral<4>::key;
    }

    // this is a bit bogus, but need to return something
    virtual unsigned getToTypeKey() override
    {
      return shards::Quadrilateral<4>::key;
    }

    virtual std::string getFromTopoPartName() override {
      shards::CellTopology cell_topo(getFromTopology());
      return cell_topo.getName();
    }
    virtual std::string getToTopoPartName() override {
      shards::CellTopology cell_topo(getToTopology());
      return cell_topo.getName();
    }

    virtual const CellTopologyData * getFromTopology() override { return shards::getCellTopologyData< shards::Quadrilateral<4> >(); }
    virtual const CellTopologyData * getToTopology() override { return shards::getCellTopologyData< shards::Quadrilateral<4> >(); }

    void
    createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                      stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                      stk::mesh::FieldBase *proc_rank_field=0) override
    {
      // quads

      unsigned num_edges_marked=0;
      for (int iedge = 0; iedge < 4; iedge++)
        {
          unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
          if (num_nodes_on_edge)
            {
              ++num_edges_marked;
            }
        }
      unsigned num_faces_marked = 0;
      stk::mesh::EntityRank rank = (m_eMesh.get_spatial_dim() == 2 ? stk::topology::ELEMENT_RANK : m_eMesh.face_rank());
      if ( new_sub_entity_nodes[rank].size() )
        {
          int iface = 0;
          num_faces_marked = new_sub_entity_nodes[rank][iface].size();
        }
      if ( num_edges_marked == 4 && num_faces_marked == 1)
        {
          // vector<stk::mesh::Entity>::iterator start = element_pool;
          // m_bp[0]->createNewElements(eMesh, nodeRegistry, element, new_sub_entity_nodes, element_pool, proc_rank_field);
          // ptrdiff_t n_gen = element_pool - start;
          // if (0)
          //   std::cout << "n_gen = " << n_gen << std::endl;
          throw std::runtime_error("not ready");
        }
      else
        {
          m_bp[0]->createNewElements(eMesh, nodeRegistry, element, new_sub_entity_nodes, element_pool, ft_element_pool, proc_rank_field);
          m_bp[1]->createNewElements(eMesh, nodeRegistry, element, new_sub_entity_nodes, element_pool, ft_element_pool, proc_rank_field);
          m_bp[2]->createNewElements(eMesh, nodeRegistry, element, new_sub_entity_nodes, element_pool, ft_element_pool, proc_rank_field);
        }

    }
  };

}

#undef Q_CENTROID_N_EV
#undef Q_CV_EV
#undef PRINT_PARTS

#endif
