// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinerPattern_Pyr5_Het_N_hpp
#define adapt_RefinerPattern_Pyr5_Het_N_hpp


#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "RefinerPattern_Quad4_Quad4_Transition.hpp"
#include "RefinerPattern_Quad4_Het_N.hpp"

#include <percept/PerceptBoostArray.hpp>
#include <adapt/DiscretizePyr.hpp>

#ifndef NDEBUG
#include <percept/mesh/geometry/volume/VolumeUtil.hpp>
#endif

namespace percept {

#define LOCAL_DEBUG 0
#define LOCAL_DEBUG_PP 0

  extern bool s_do_transition_break;

  // Some explanation: Pyr local (non-hanging-node) refinement pattern creates a heterogeneous mesh, so we create two
  //   sub-patterns to deal with the two different resulting topologies (, tet and pyramid).  A third
  //   (parent) pattern is created to refer to the two sub-patterns, akin to URP_Heterogeneous_3D.

  // this pattern is only intended to be used for "transition" element breaking, that is, all hanging
  // nodes are removed by filling the interior of the pyr with sub-tets or pyramids

  //================================================================================================================================================================


  //================================================================================================================================================================
  // Just for ensuring tets are processed by the marker - this pattern does nothing, it doesn't break tets...
  //================================================================================================================================================================

  /**
  struct TetTempPartialNoBreak {};

  // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
  template <>
  class RefinerPattern<shards::Tetrahedron<4>, shards::Tetrahedron<4>, -1, TetTempPartialNoBreak > : public URP<shards::Tetrahedron<4>, shards::Tetrahedron<4>  >
  */

  //================================================================================================================================================================
  // Just for ensuring pyramids are processed by the marker - this pattern does nothing, it doesn't break pyramids...
  //================================================================================================================================================================

  /**
  struct PyrTempPartialNoBreak {};

  // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
  template <>
  class RefinerPattern<shards::Pyramid<5>, shards::Pyramid<5>, -1, PyrTempPartialNoBreak > : public URP<shards::Pyramid<5>, shards::Pyramid<5>  >
  */


  //================================================================================================================================================================
  // Heterogeneous child patterns
  //================================================================================================================================================================

  // Pyr to Tet
  // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
  struct PyrTetPartial {};
  template <>
  class RefinerPattern<shards::Pyramid<5>, shards::Tetrahedron<4>, -1, PyrTetPartial > : public URP<shards::Pyramid<5>, shards::Tetrahedron<4>  >
  {
    typedef RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1, QuadTransition > FaceBreaker;
    typedef RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1, TriHangingNode >  TriFaceBreakerType;

    FaceBreaker * m_face_breaker;
    TriFaceBreakerType * m_face_breaker_tri;

    TransitionElementType *m_transition_element_field;
    bool m_transition_element_field_set;

  public:

    RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Pyramid<5>, shards::Tetrahedron<4>  >(eMesh),
                                                                                                  m_face_breaker(0),m_face_breaker_tri(0),
                                                                                                  m_transition_element_field(0), m_transition_element_field_set(false)
    {
      m_primaryEntityRank = m_eMesh.element_rank();
      bool sameTopology = false;
      setNeededParts(eMesh, block_names, sameTopology);
      Elem::StdMeshObjTopologies::bootstrap();

      if (LOCAL_DEBUG_PP)
        {
          std::cout << "tmp PyrTetPartial____ printParts after all bp= \n" ;
          printParts(this);
        }

      m_face_breaker =  new FaceBreaker(eMesh, block_names) ;
      m_face_breaker_tri = new TriFaceBreakerType(eMesh, block_names);
    }

    ~RefinerPattern()
    {
      if (m_face_breaker) delete m_face_breaker;
      if (m_face_breaker_tri) delete m_face_breaker_tri;
    }

    void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
    {
      EXCEPTWATCH;
      bp.resize(3);

      bp[0] = this;
      bp[1] = m_face_breaker;
      bp[2] = m_face_breaker_tri;
    }

    virtual void doBreak() override {}
    void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
    {
      needed_entities.resize(3);
      needed_entities[0].first = m_eMesh.edge_rank();
      needed_entities[1].first = m_eMesh.face_rank();
      needed_entities[2].first = m_eMesh.element_rank();
      setToOne(needed_entities);

      int faces[5] = {0,0,0,0,1};
      needed_entities[1].third.resize(5);
      needed_entities[1].third.assign(faces,faces+5);
    }

    // FIXME - for now, create more than we need (to fix this right we need a change to the Refiner.cpp interface)
    virtual unsigned getNumNewElemPerElem() override { return 48; }

    void
    createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                      stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                      vector<stk::mesh::Entity>::iterator& ft_element_pool,
                      stk::mesh::FieldBase *proc_rank_field=0) override
    {
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
      static vector<pyr_to_tet_tuple_type> elems_tet(48);
      static vector<pyr_to_tet_tuple_type_local> elems_tet_local(48);
      static vector<pyr_to_pyr_tuple_type> elems_pyr(24);
      static vector<pyr_to_pyr_tuple_type_local> elems_pyr_local(24);
      unsigned num_new_elems=0;

      shards::CellTopology cell_topo(cell_topo_data);
      const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::topology::NODE_RANK);
      //CoordinatesFieldType* coordField = eMesh.get_coordinates_field();

      if (!m_transition_element_field_set)
        {
          m_transition_element_field_set = true;
          m_transition_element_field = eMesh.get_transition_element_field();
          if(m_transition_element_field && m_transition_element_field->name() != "transition_element_3") m_transition_element_field = nullptr;
        }

      std::vector<stk::mesh::Part*> add_parts;
      std::vector<stk::mesh::Part*> remove_parts;
      add_parts = m_toParts;

      unsigned edge_marks[8] = { 0,0,0,0, 0,0,0,0 };
      unsigned num_edges_marked=0;
      for (int iedge = 0; iedge < 8; iedge++)
        {
          unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
          if (num_nodes_on_edge)
            {
              edge_marks[iedge] = 1;
              ++num_edges_marked;
            }
        }

      unsigned num_faces_marked = 0;
      unsigned face_marks[5] = {0,0,0,0,0};
      stk::mesh::EntityRank rank = m_eMesh.face_rank();

      for (int iface = 0; iface < 5; iface++)
        {
          if ( new_sub_entity_nodes[rank].size() )
            {
              if (new_sub_entity_nodes[rank][iface].size())
                {
                  face_marks[iface] = 1;
                  ++num_faces_marked;
                }
            }
        }

      if (LOCAL_DEBUG)
        {
          std::vector<unsigned> em(&edge_marks[0], &edge_marks[0]+8);
          std::vector<unsigned> fm(&face_marks[0], &face_marks[0]+5);
          std::cout << "PyrTetPartial:createNewElements:: element= " << m_eMesh.identifier(element)
                    << " num_edges_marked= " << num_edges_marked << " num_faces_marked= " << num_faces_marked
                    << "\n em= " << em
                    << "\n fm= " << fm
                    << std::endl;
        }

      if (num_edges_marked == 0 && num_faces_marked == 0)
        return;

      DiscretizePyr::discretize(m_eMesh, element, edge_marks, face_marks,
                                elems_tet_local, elems_pyr_local);

#define Q_CENTROID_N_EV (useAltNode ? m_eMesh.identifier(new_nodes[0]) : NN(m_primaryEntityRank, 0) )

#define Q_CV_EV(i) \
      ( DiscretizePyr::is_mapped_v(i) ? VERT_N(i) :    \
        ( DiscretizePyr::is_mapped_e(i) ? EDGE_N(i - DiscretizePyr::edge_offset) : \
          ( DiscretizePyr::is_mapped_f(i) ? FACE_N(i - DiscretizePyr::face_offset_tri) : Q_CENTROID_N_EV ) ) )

      if (LOCAL_DEBUG == 2)
        std::cout << "PyrTetPartial:createNewElements:: elem= " << m_eMesh.identifier(element)
                  << " m_primaryEntityRank= " << m_primaryEntityRank
          //<< " Q_CENTROID_N_EV = " << Q_CENTROID_N_EV
                  << " new_sub_entity_nodes[m_primaryEntityRank].size()= "
                  << new_sub_entity_nodes[m_primaryEntityRank].size()
                  << " " << (new_sub_entity_nodes[m_primaryEntityRank].size() ?
                       (new_sub_entity_nodes[m_primaryEntityRank][0].size() ?
                        new_sub_entity_nodes[m_primaryEntityRank][0][0] : 123456 ) : 1234567)
                  << std::endl;

      std::vector<stk::mesh::Entity> new_nodes;
      bool useAltNode = false;
      if (new_sub_entity_nodes.size() > static_cast<unsigned>(m_primaryEntityRank) && new_sub_entity_nodes[m_primaryEntityRank].size() && !new_sub_entity_nodes[m_primaryEntityRank][0].size())
        {
          //throw std::runtime_error("no centroid node available");
           // std::cout << " bad element: " << std::endl;
           // m_eMesh.print(element);
           // m_eMesh.dump_vtk("bad.vtk");
          if (1)
            {
              if (!m_eMesh.getEntitiesUsingIdServerNewNodes(1, new_nodes))
                {
                  std::cout << "P[" << m_eMesh.get_rank() << "] node_rank entity pool deplenished in RefinerPattern_Pyr5_Het_N" << std::endl;
                  throw std::logic_error("node rank entity pool deplenished");
                }
              useAltNode = true;
            }
        }

      if (useAltNode)
        {
          nodeRegistry.forceInMap(new_nodes, NodeRegistry::NR_MARK, element, stk::topology::ELEMENT_RANK, 0u);
          nodeRegistry.prolongateCoords(element, stk::topology::ELEMENT_RANK, 0u);
          new_sub_entity_nodes[m_primaryEntityRank][0].resize(1);
          new_sub_entity_nodes[m_primaryEntityRank][0][0] = m_eMesh.identifier(new_nodes[0]);
        }

      num_new_elems = elems_tet_local.size();
      elems_tet.resize(num_new_elems);
      for (unsigned ielem=0; ielem < num_new_elems; ielem++)
        {
          pyr_to_tet_tuple_type ht;
          for (unsigned ii=0; ii < 4; ++ii)
            {
              unsigned ee = elems_tet_local[ielem][ii];
              bool isf = DiscretizePyr::is_mapped_f(ee);
              if (isf)
                {
                  int ifaceOrd = ee - DiscretizePyr::face_offset_tri;

                  if (new_sub_entity_nodes[m_eMesh.face_rank()].size() &&
                      new_sub_entity_nodes[m_eMesh.face_rank()][ifaceOrd].size() &&
                      new_sub_entity_nodes[m_eMesh.face_rank()][ifaceOrd][0])
                    {
                      stk::mesh::EntityId id = new_sub_entity_nodes[m_eMesh.face_rank()][ifaceOrd][0];
                      if (0) std::cout << "no error, id= " << id << std::endl;
                    }
                  else
                    {
                      std::cout << "Pyr_het tet error: ee= " << ee << " isf= " << isf << " ifaceOrd= " << ifaceOrd << std::endl;

                      std::set<stk::mesh::Entity> list;
                      list.insert(element);
                      m_eMesh.dump_vtk("err-face.vtk", false, &list);

                      SubDimCell_SDCEntityType subDimEntity(&m_eMesh);
                      nodeRegistry.getSubDimEntity(subDimEntity, element, m_eMesh.face_rank(), ifaceOrd);
                      std::vector<stk::mesh::Entity> nv;
                      for (unsigned jj=0; jj < subDimEntity.size(); ++jj)
                        {
                          nv.push_back(subDimEntity[jj]);
                        }
                      m_eMesh.dump_vtk(nv, "err-nodes.vtk");
                      throw std::logic_error("Pyr_het tet error");
                    }
                }
              ht[ii] = Q_CV_EV(ee);
            }
          elems_tet[ielem] = ht;

          if (LOCAL_DEBUG)
            std::cout << "PyrTetPartial:: tet ielem= " << ielem << " num_new_elems= " << num_new_elems << " elems_local= " << elems_tet_local[ielem]
                      << " new_sub_entity_nodes[elem_rank].size= " << new_sub_entity_nodes[m_primaryEntityRank].size()
                      << " new_sub_entity_nodes[elem_rank][0].size= " << new_sub_entity_nodes[m_primaryEntityRank][0].size()
                      << " m_primaryEntityRank= " << m_primaryEntityRank
                      << "\n new tet elem= " << elems_tet[ielem]
                      << std::endl;
        }

      if (LOCAL_DEBUG == 2)
        std::cout << "tmp RefinerPattern_Pyr5_Tet4_N::num_edges_marked= " << num_edges_marked
                  << " num_new_elems= " << num_new_elems
                  << std::endl;

      for (unsigned ielem=0; ielem < elems_tet.size(); ielem++)
        {
          stk::mesh::Entity newElement = *element_pool;

          if (m_transition_element_field)
            {
              int *transition_element = 0;
              transition_element = stk::mesh::field_data( *m_transition_element_field , newElement );
              //if (num_edges_marked == 12 && num_faces_marked == 6)
              if (num_edges_marked == 8)
                {
                  transition_element[0] = 0;
                  // FIXME error???
                  //throw std::runtime_error("logic error in RefinerPattern_Pyr5_Het_N");
                }
              else
                {
                  transition_element[0] = 1;
                }
            }

          if (proc_rank_field)
            {
              double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
              if (fdata)
                fdata[0] = double(eMesh.owner_rank(newElement));
            }

          change_entity_parts(eMesh, element, newElement);

          {
            if (!elems_tet[ielem][0])
              {
                std::cout << "P[" << eMesh.get_rank() << "] nid = 0 << " << std::endl;
                //exit(1);
              }
          }

          // 4 nodes of the new tets
          for (unsigned ii=0; ii < 4; ++ii)
            {
              stk::mesh::Entity n0 = eMesh.createOrGetNode(elems_tet[ielem][ii]);
              eMesh.get_bulk_data()->declare_relation(newElement, n0, ii);
            }

          unsigned nchild = eMesh.numChildren(element);
          //set_parent_child_relations(eMesh, element, newElement, ielem);
          set_parent_child_relations(eMesh, element, newElement,  *ft_element_pool, nchild);

          if (0) std::cout << "tmp createNewElements: isParentElement= " << m_eMesh.isParentElement(element, false) << std::endl;
          std::vector<stk::mesh::Entity> elements(1,element);
          eMesh.prolongateElementFields( elements, newElement);

          ++ft_element_pool;
          ++element_pool;

        }
    }

  };

  // Pyramid to Pyramid
  // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
  struct PyrPyrPartial {};
  template <>
  class RefinerPattern<shards::Pyramid<5>, shards::Pyramid<5>, -1, PyrPyrPartial > : public URP< shards::Pyramid<5>, shards::Pyramid<5>  >
  {
    typedef RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1, QuadTransition > QuadFaceBreaker;
    typedef RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1 > TriFaceBreaker;
    QuadFaceBreaker * m_quad_face_breaker;
    TriFaceBreaker * m_tri_face_breaker;

    TransitionElementType *m_transition_element_field;
    bool m_transition_element_field_set;

  public:

    RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Pyramid<5>, shards::Pyramid<5>  >(eMesh),
                                                                                                  m_quad_face_breaker(0), m_tri_face_breaker(0),
                                                                                                  m_transition_element_field(0), m_transition_element_field_set(false)
    {
      m_primaryEntityRank = m_eMesh.element_rank();
      bool sameTopology = true;
      setNeededParts(eMesh, block_names, sameTopology);
      if (LOCAL_DEBUG_PP)
        {
          std::cout << "tmp PyrPyrPartial____ printParts after all bp= \n" ;
          printParts(this);
        }
      Elem::StdMeshObjTopologies::bootstrap();

      m_quad_face_breaker =  new QuadFaceBreaker(eMesh, block_names) ;
      m_tri_face_breaker =  new TriFaceBreaker(eMesh, block_names) ;
    }

    ~RefinerPattern()
    {
      if (m_quad_face_breaker) delete m_quad_face_breaker;
      if (m_tri_face_breaker) delete m_tri_face_breaker;
    }

    void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
    {
      EXCEPTWATCH;
      bp.resize(3);

      bp[0] = this;
      bp[1] = m_quad_face_breaker;
      bp[2] = m_tri_face_breaker;
    }

    virtual void doBreak() override {}
    void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
    {
      needed_entities.resize(3);
      needed_entities[0].first = m_eMesh.edge_rank();
      needed_entities[1].first = m_eMesh.face_rank();
      needed_entities[2].first = m_eMesh.element_rank();
      setToOne(needed_entities);

      int faces[5] = {0,0,0,0,1};
      needed_entities[1].third.resize(5);
      needed_entities[1].third.assign(faces,faces+5);
    }

    // FIXME - for now, create more than we need (to fix this right we need a change to the Refiner.cpp interface)
    virtual unsigned getNumNewElemPerElem() override { return 24; }

    void
    createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                      stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                      vector<stk::mesh::Entity>::iterator& ft_element_pool,
                      stk::mesh::FieldBase *proc_rank_field=0) override
    {
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
      static vector<pyr_to_tet_tuple_type> elems_tet(48);
      static vector<pyr_to_tet_tuple_type_local> elems_tet_local(48);
      static vector<pyr_to_pyr_tuple_type> elems_pyr(48);
      static vector<pyr_to_pyr_tuple_type_local> elems_pyr_local(48);
      unsigned num_new_pyr_elems=0;

      shards::CellTopology cell_topo(cell_topo_data);
      const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::topology::NODE_RANK);
      //CoordinatesFieldType* coordField = eMesh.get_coordinates_field();

      if (!m_transition_element_field_set)
        {
          m_transition_element_field_set = true;
          m_transition_element_field = eMesh.get_transition_element_field();
          if(m_transition_element_field->name() != "transition_element_3") m_transition_element_field = nullptr;
        }

      std::vector<stk::mesh::Part*> add_parts;
      std::vector<stk::mesh::Part*> remove_parts;
      add_parts = m_toParts;

      unsigned edge_marks[8] = { 0,0,0,0, 0,0,0,0,  };
      unsigned num_edges_marked=0;
      for (int iedge = 0; iedge < 8; iedge++)
        {
          unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
          if (num_nodes_on_edge)
            {
              edge_marks[iedge] = 1;
              ++num_edges_marked;
            }
        }

      unsigned face_marks[5] = {0,0,0,0,0};
      unsigned num_faces_marked = 0;
      stk::mesh::EntityRank rank = m_eMesh.face_rank();

      for (int iface = 0; iface < 5; iface++)
        {
          if ( new_sub_entity_nodes[rank].size() )
            {
              if (new_sub_entity_nodes[rank][iface].size())
                {
                  face_marks[iface] = 1;
                  ++num_faces_marked;
                }
            }
        }

      if (LOCAL_DEBUG)
        {
          std::vector<unsigned> em(&edge_marks[0], &edge_marks[0]+8);
          std::vector<unsigned> fm(&face_marks[0], &face_marks[0]+5);
          std::cout << "PyrPyrPartial:createNewElements:: element= " << m_eMesh.identifier(element)
                    << " num_edges_marked= " << num_edges_marked << " num_faces_marked= " << num_faces_marked
                    << "\n em= " << em
                    << "\n fm= " << fm
                    << std::endl;
        }

      if (num_edges_marked == 0 && num_faces_marked == 0)
        return;

      DiscretizePyr::discretize(m_eMesh, element, edge_marks, face_marks,
                                elems_tet_local, elems_pyr_local);


      if (LOCAL_DEBUG == 2)
        std::cout << "PyrPyrPartial:createNewElements:: pyr elem= " << m_eMesh.identifier(element) << " topo= " << m_eMesh.bucket(element).topology()
          //<< " Q_CENTROID_N_EV = " << Q_CENTROID_N_EV
                  << " new_sub_entity_nodes[m_primaryEntityRank].size()= "
                  << new_sub_entity_nodes[m_primaryEntityRank].size()
                  << " " << (new_sub_entity_nodes[m_primaryEntityRank].size() ?
                       (new_sub_entity_nodes[m_primaryEntityRank][0].size() ?
                        new_sub_entity_nodes[m_primaryEntityRank][0][0] : 123456 ) : 1234567)
                  << std::endl;

      std::vector<stk::mesh::Entity> new_nodes;
      bool useAltNode = false;
      if (new_sub_entity_nodes.size() > static_cast<unsigned>(m_primaryEntityRank) && new_sub_entity_nodes[m_primaryEntityRank].size() && !new_sub_entity_nodes[m_primaryEntityRank][0].size())
        {
          //throw std::runtime_error("no centroid node availabled");
           // std::cout << " bad element: " << std::endl;
           // m_eMesh.print(element);
           // m_eMesh.dump_vtk("bad.vtk");
          if (1)
            {
              if (!m_eMesh.getEntitiesUsingIdServerNewNodes( 1, new_nodes))
                {
                  std::cout << "P[" << m_eMesh.get_rank() << "] node_rank entity pool deplenished in RefinerPattern_Pyr5_Het_N" << std::endl;
                  throw std::logic_error("node rank entity pool deplenished");
                }
              useAltNode = true;
            }
        }

      if (useAltNode)
        {
          nodeRegistry.forceInMap(new_nodes, NodeRegistry::NR_MARK, element, stk::topology::ELEMENT_RANK, 0u);
          nodeRegistry.prolongateCoords(element, stk::topology::ELEMENT_RANK, 0u);
          new_sub_entity_nodes[m_primaryEntityRank][0].resize(1);
          new_sub_entity_nodes[m_primaryEntityRank][0][0] = m_eMesh.identifier(new_nodes[0]);
        }

      num_new_pyr_elems = elems_pyr_local.size();
      elems_pyr.resize(num_new_pyr_elems);
      for (unsigned ielem=0; ielem < num_new_pyr_elems; ielem++)
        {
          pyr_to_pyr_tuple_type hp;
          for (unsigned ii=0; ii < 5; ++ii)
            {
              unsigned ee = elems_pyr_local[ielem][ii];

              bool isf = DiscretizePyr::is_mapped_f(ee);
              if (isf)
                {
                  int ifaceOrd = ee - DiscretizePyr::face_offset_tri;
                  if (new_sub_entity_nodes[m_eMesh.face_rank()].size() &&
                      new_sub_entity_nodes[m_eMesh.face_rank()][ifaceOrd].size() &&
                      new_sub_entity_nodes[m_eMesh.face_rank()][ifaceOrd][0])
                    {
                      stk::mesh::EntityId id = new_sub_entity_nodes[m_eMesh.face_rank()][ifaceOrd][0];
                      if (0) std::cout << "no error, id= " << id << std::endl;
                    }
                  else
                    {
                      {
                        // std::vector<unsigned> em(&edge_marks[0], &edge_marks[0]+8);
                        // std::vector<unsigned> fm(&face_marks[0], &face_marks[0]+5);
                        std::cout << "Pyr_het pyr error: ee= " << ee << " isf= " << isf << " ifaceOrd= " << ifaceOrd
                                  << " num_new_pyr_elems= " << num_new_pyr_elems << " num_faces_marked= " << num_faces_marked
                                  << " num_edges_marked= " << num_edges_marked
                                  // << "\n em= " << em
                                  // << "\n fm= " << fm
                                  << std::endl;
                      }
                      std::set<stk::mesh::Entity> list;
                      list.insert(element);
                      m_eMesh.dump_vtk("err-face.vtk", false, &list);

                      SubDimCell_SDCEntityType subDimEntity(&m_eMesh);
                      nodeRegistry.getSubDimEntity(subDimEntity, element, m_eMesh.face_rank(), ifaceOrd);
                      std::vector<stk::mesh::Entity> nv;
                      for (unsigned jj=0; jj < subDimEntity.size(); ++jj)
                        {
                          nv.push_back(subDimEntity[jj]);
                        }
                      m_eMesh.dump_vtk(nv, "err-nodes.vtk");
                      throw std::logic_error("Pyr_het pyr error");
                    }
                }
              hp[ii] = Q_CV_EV(ee);
            }
          elems_pyr[ielem] = hp;

          if (LOCAL_DEBUG)
            std::cout << "PyrPyrPartial:: pyr ielem= " << ielem << " num_new_pyr_elems= " << num_new_pyr_elems << " elems_local= " << elems_pyr_local[ielem]
                      << " new_sub_entity_nodes[elem_rank].size= " << new_sub_entity_nodes[m_primaryEntityRank].size()
                      << " new_sub_entity_nodes[elem_rank][0].size= " << new_sub_entity_nodes[m_primaryEntityRank][0].size()
                      << " m_primaryEntityRank= " << m_primaryEntityRank
                      << "\n new pyr elem= " << elems_pyr[ielem]
                      << std::endl;
        }

      if (LOCAL_DEBUG)
        std::cout << "tmp RefinerPattern_Pyr5_Pyr::num_edges_marked= " << num_edges_marked
                  << " num_new_pyr_elems= " << num_new_pyr_elems << " elems_pyr.size= " << elems_pyr.size()
                  << std::endl;

      bool err = false;
      std::set<stk::mesh::Entity> list;
      list.insert(element);
      for (unsigned ielem=0; ielem < elems_pyr.size(); ielem++)
        {
          stk::mesh::Entity newElement = *element_pool;

          if (m_transition_element_field)
            {
              int *transition_element = 0;
              transition_element = stk::mesh::field_data( *m_transition_element_field , newElement );
              transition_element[0] = 1;
            }

          if (proc_rank_field)
            {
              double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
              if (fdata)
                fdata[0] = double(eMesh.owner_rank(newElement));
            }

          change_entity_parts(eMesh, element, newElement);

          {
            if (!elems_pyr[ielem][0])
              {
                std::cout << "P[" << eMesh.get_rank() << "] nid = 0 << " << std::endl;
                //exit(1);
              }
          }

          // 5 nodes of the new pyramids
          for (unsigned ii=0; ii < 5; ++ii)
            {
              stk::mesh::Entity n0 = eMesh.createOrGetNode(elems_pyr[ielem][ii]);
              eMesh.get_bulk_data()->declare_relation(newElement, n0, ii);
            }

          unsigned nchild = eMesh.numChildren(element);
          //set_parent_child_relations(eMesh, element, newElement, ielem);
          set_parent_child_relations(eMesh, element, newElement,  *ft_element_pool, nchild);

          if (0) std::cout << "tmp PyrPyrPartial createNewElements: isParentElement= " << m_eMesh.isParentElement(element, false) << std::endl;
          if (0) std::cout << "tmp PyrPyrPartial createNewElements: id= " << eMesh.identifier(newElement) << std::endl;
          std::vector<stk::mesh::Entity> elements(1,element);
          eMesh.prolongateElementFields( elements, newElement);

          ++ft_element_pool;
          ++element_pool;

        }
      if (err)
        {
          eMesh.dump_vtk("err-pyr.vtk", false, &list);

          throw std::runtime_error("neg vol PyrPyrPartial");
        }

    }

  };

  //================================================================================================================================================================
  // Parent pattern
  //================================================================================================================================================================

  struct PyrHet {};

  template<>
  class RefinerPattern<shards::Pyramid<5>, shards::Pyramid<5>, 24, PyrHet > : public UniformRefinerPatternBase
  {
    std::vector<UniformRefinerPatternBase *> m_bp;
  public:
    std::vector<UniformRefinerPatternBase *> m_bp_exported;

    typedef RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1, QuadTransition > QuadFaceBreaker;
    typedef RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1 > TriFaceBreaker;
    QuadFaceBreaker * m_quad_face_breaker;
    TriFaceBreaker * m_tri_face_breaker;

  protected:

    percept::PerceptMesh& m_eMesh;

  public:

    RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  m_eMesh(eMesh)
    {
      m_primaryEntityRank = stk::topology::ELEMENT_RANK;

      Elem::StdMeshObjTopologies::bootstrap();

      m_bp.resize(4);

      m_bp[0] = new  RefinerPattern<shards::Pyramid<5>,  shards::Tetrahedron<4>, -1, PyrTetPartial >         (eMesh, block_names) ;
      m_bp[1] = new  RefinerPattern<shards::Pyramid<5>,  shards::Pyramid<5>,     -1, PyrPyrPartial >         (eMesh, block_names) ;
      m_bp[2] = new  RefinerPattern<shards::Tetrahedron<4>, shards::Tetrahedron<4>, -1, TetTempPartialNoBreak > (eMesh, block_names) ;
      m_bp[3] = new  RefinerPattern<shards::Pyramid<5>,     shards::Pyramid<5>,     -1, PyrTempPartialNoBreak > (eMesh, block_names) ;

      m_bp_exported.resize(0);
      m_bp_exported.push_back(m_bp[2]);
      m_bp_exported.push_back(m_bp[3]);

      m_bp[0]->setNeededParts(eMesh, block_names, false);
      m_bp[1]->setNeededParts(eMesh, block_names, true);
      m_bp[2]->setNeededParts(eMesh, block_names, true);
      m_bp[3]->setNeededParts(eMesh, block_names, true);

      // repeat setNeededParts to catch newly created parts
      bool skipConvertedParts = false;
      m_bp[0]->setNeededParts(eMesh, block_names, false, skipConvertedParts);
      m_bp[1]->setNeededParts(eMesh, block_names, true, skipConvertedParts);
      m_bp[2]->setNeededParts(eMesh, block_names, true, skipConvertedParts);
      m_bp[3]->setNeededParts(eMesh, block_names, true, skipConvertedParts);

      if (1)
        {
          for (unsigned ibp=0; ibp < m_bp.size(); ibp++)
            {
              bool mergeParts = false;
              mergeOrAddParts(m_bp[ibp], this, mergeParts);
            }
        }

      if (LOCAL_DEBUG_PP)
        {
          std::cout << "tmp PyrHet____ printParts this= \n" ;
          printParts(this);
        }

      m_quad_face_breaker =  new QuadFaceBreaker(eMesh, block_names) ;
      m_tri_face_breaker =  new TriFaceBreaker(eMesh, block_names) ;
    }

    ~RefinerPattern()
    {
      if (m_quad_face_breaker) delete m_quad_face_breaker;
      if (m_tri_face_breaker) delete m_tri_face_breaker;
      for (unsigned ibp=0; ibp < m_bp.size(); ibp++)
        {
          if (m_bp[ibp]) delete m_bp[ibp];
        }
    }

    void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
    {
      EXCEPTWATCH;
      bp.resize(0);

      bp = m_bp;
      bp.push_back(m_quad_face_breaker);
      bp.push_back(m_tri_face_breaker);
    }

    virtual void doBreak() override {
      throw std::runtime_error("shouldn't call RefinerPattern_Pyr5_Het_N::doBreak()");
    }

    void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
    {
      needed_entities.resize(3);
      needed_entities[0].first = m_eMesh.edge_rank();
      needed_entities[1].first = m_eMesh.face_rank();
      needed_entities[2].first = m_eMesh.element_rank();
      setToOne(needed_entities);

      int faces[5] = {0,0,0,0,1};
      needed_entities[1].third.resize(5);
      needed_entities[1].third.assign(faces,faces+5);

    }

    virtual unsigned getNumNewElemPerElem() override { return 48; }

    virtual unsigned getFromTypeKey() override
    {
      return shards::Pyramid<5>::key;
    }

    // this is a bit bogus, but need to return something
    virtual unsigned getToTypeKey() override
    {
      return shards::Pyramid<5>::key;
    }

    virtual std::string getFromTopoPartName() override {
      shards::CellTopology cell_topo(getFromTopology());
      return cell_topo.getName();
    }
    virtual std::string getToTopoPartName() override {
      shards::CellTopology cell_topo(getToTopology());
      return cell_topo.getName();
    }

    virtual const CellTopologyData * getFromTopology() override { return shards::getCellTopologyData< shards::Pyramid<5> >(); }
    virtual const CellTopologyData * getToTopology() override { return shards::getCellTopologyData< shards::Pyramid<5> >(); }

    void
    createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                      stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                      stk::mesh::FieldBase *proc_rank_field=0) override
    {
      unsigned num_edges_marked=0;
      for (int iedge = 0; iedge < 8; iedge++)
        {
          unsigned num_nodes_on_edge = new_sub_entity_nodes[m_eMesh.edge_rank()][iedge].size();
          if (num_nodes_on_edge)
            {
              ++num_edges_marked;
            }
        }
      unsigned num_faces_marked = 0;
      stk::mesh::EntityRank rank = m_eMesh.face_rank();

      int iface = 4;
      if ( new_sub_entity_nodes[rank].size() )
        {
          if (new_sub_entity_nodes[rank][iface].size())
            ++num_faces_marked;
        }

      if ( num_edges_marked == 8 && num_faces_marked == 1)
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
          for (unsigned ii=0; ii < m_bp.size(); ++ii)
            {
              //m_bp[ii]->m_ep_begin = m_ep_begin;
              //m_bp[ii]->m_ep_end = m_ep_end;
              m_bp[ii]->createNewElements(eMesh, nodeRegistry, element, new_sub_entity_nodes, element_pool, ft_element_pool, proc_rank_field);
            }
        }

    }
  };

}

#undef Q_VERT_N
#undef Q_EDGE_N
#undef Q_CENTROID_N
#undef Q_CENTROID_N_EV
#undef Q_CV_EV

#undef LOCAL_DEBUG
#undef LOCAL_DEBUG_PP

#endif
