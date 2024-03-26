// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_hpp
#define adapt_UniformRefinerPattern_hpp

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <sstream>
#include <cmath>
#include <math.h>
#include <stk_mesh/base/Types.hpp>

#include "Teuchos_RCP.hpp"

#include <percept/stk_mesh.hpp>
#include <percept/PerceptBoostArray.hpp>


#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <percept/Util.hpp>
#include <percept/PerceptMesh.hpp>
#include <adapt/NodeRegistry.hpp>
#include <percept/function/FieldFunction.hpp>

#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

// more efficient fixElementSides implementation using parent/child relations
#define NEW_FIX_ELEMENT_SIDES 1

// set to 0 for doing global (and thus more efficient) computation of node coords and adding to parts
#define STK_ADAPT_URP_LOCAL_NODE_COMPS 0

// set to 1 to turn on some print tracing and cpu/mem tracing
#define FORCE_TRACE_PRINT_ONLY 0
#define TRACE_STAGE_PRINT_ON 0
#define TRACE_STAGE_PRINT (TRACE_STAGE_PRINT_ON && (m_eMesh.get_rank()==0))

#if TRACE_STAGE_PRINT_ON
#  define TRACE_PRINT(a) do { trace_print(a); } while(0)
#  define TRACE_CPU_TIME_AND_MEM_0(a) do { Util::trace_cpu_time_and_mem_0(a); } while(0)
#  define TRACE_CPU_TIME_AND_MEM_1(a) do { Util::trace_cpu_time_and_mem_1(a); } while(0)
#else
#  if FORCE_TRACE_PRINT_ONLY
#    define TRACE_PRINT(a) do { trace_print(a); } while(0)
#  else
#    define TRACE_PRINT(a) do {} while(0)
#  endif
#  define TRACE_CPU_TIME_AND_MEM_0(a) do { } while(0)
#  define TRACE_CPU_TIME_AND_MEM_1(a) do { } while(0)
#endif



  namespace percept {

    extern bool allow_single_refine;

    using std::vector;

    using shards::CellTopology;

    typedef std::vector<std::vector<std::string> > BlockNamesType;
    typedef std::map<std::string, std::string> StringStringMap;

    typedef std::vector<std::vector<std::vector<stk::mesh::EntityId> > > NewSubEntityNodesType;
    typedef Elem::StdMeshObjTopologies::RefTopoX RefTopoX;

    typedef Elem::StdMeshObjTopologies::RefinementTopologyExtraEntry *RefTopoX_arr;

    // useful tools
#define NODE_COORD(node) static_cast<double*>(stk::mesh::field_data( *eMesh.get_coordinates_field() , node ))
#define VERT_COORD(ivert) NODE_COORD(elem_nodes[ivert].entity())
#define EDGE_COORD(iedge,inode) NODE_COORD(elem_nodes[cell_topo_data->edge[iedge].node[inode]].entity())
#define FACE_COORD(iface,inode) NODE_COORD(elem_nodes[cell_topo_data->side[iface].node[inode]].entity())

    /// 2D array new_sub_entity_nodes[entity_rank][ordinal_of_node_on_sub_dim_entity]

#define VERT_N(i) m_eMesh.identifier(elem_nodes[i].entity())

// #define VERT_N(i) elem_nodes[i].entity()->identifier()
// #define EDGE_N(i) new_sub_entity_nodes[m_eMesh.edge_rank()][i][0]
// #define FACE_N(i) new_sub_entity_nodes[m_eMesh.face_rank()][i][0]
// #define NN(i_entity_rank, j_ordinal_on_subDim_entity) new_sub_entity_nodes[i_entity_rank][j_ordinal_on_subDim_entity][0]


#define DEBUG_URP_HPP 0

#if 1

#  define NN(i_entity_rank, j_ordinal_on_subDim_entity) new_sub_entity_nodes[i_entity_rank][j_ordinal_on_subDim_entity][0]

#  define NN_Q(i_entity_rank, j_ordinal_of_entity, k_ordinal_of_node_on_entity) \
       new_sub_entity_nodes[i_entity_rank][j_ordinal_of_entity][k_ordinal_of_node_on_entity]

#  define NN_Q_P(i_entity_rank, j_ordinal_of_entity, k_ordinal_of_node_on_entity, perm) \
       new_sub_entity_nodes[i_entity_rank][j_ordinal_of_entity][perm[k_ordinal_of_node_on_entity]]

#else


#  define NN(i_entity_rank, j_ordinal_on_subDim_entity) \
    ( ( ((unsigned)i_entity_rank < new_sub_entity_nodes.size()) && \
        ((unsigned)j_ordinal_on_subDim_entity < new_sub_entity_nodes[i_entity_rank].size()) \
        && (new_sub_entity_nodes[i_entity_rank][j_ordinal_on_subDim_entity].size() > 0) ) ? \
      (new_sub_entity_nodes[i_entity_rank][j_ordinal_on_subDim_entity][0]) : \
      ( std::cout << new_sub_entity_nodes << " " << __FILE__ << ":" << __LINE__, \
        throw std::runtime_error("error in accessing new_sub_entity_nodes [rank:ord]= "+ \
           toString(i_entity_rank)+":"+toString(j_ordinal_on_subDim_entity)+" F: "+std::string(__FILE__)+":"+toString(__LINE__) ),  0u ) )


#  define NN_Q(i_entity_rank, j_ordinal_of_entity, k_ordinal_of_node_on_entity) \
    new_sub_entity_nodes_check(eMesh, element, new_sub_entity_nodes, i_entity_rank, j_ordinal_of_entity, k_ordinal_of_node_on_entity)

#  define NN_Q_P(i_entity_rank, j_ordinal_of_entity, k_ordinal_of_node_on_entity, perm) \
    new_sub_entity_nodes_check_perm(eMesh, element, new_sub_entity_nodes, i_entity_rank, j_ordinal_of_entity, k_ordinal_of_node_on_entity, perm)

#endif

#define EDGE_N(i) NN(m_eMesh.edge_rank(), i)
#define FACE_N(i) NN(m_eMesh.face_rank(), i)

    //#define EDGE_N_Q(iedge, inode_on_edge) new_sub_entity_nodes[m_eMesh.edge_rank()][iedge][inode_on_edge]
    //#define FACE_N_Q(iface, inode_on_face) new_sub_entity_nodes[m_eMesh.face_rank()][iface][inode_on_face]

#define EDGE_N_Q(iedge, inode_on_edge) NN_Q(m_eMesh.edge_rank(), iedge, inode_on_edge)
#define FACE_N_Q(iface, inode_on_face) NN_Q(m_eMesh.fade_rank(), iface, inode_on_face)


    struct SierraPort {};

    inline int new_sub_entity_nodes_check(PerceptMesh& eMesh, stk::mesh::Entity element, NewSubEntityNodesType& new_sub_entity_nodes, int i_entity_rank, int j_ordinal_of_entity, int k_ordinal_of_node_on_entity)
    {
      try {
        VERIFY_OP_ON((unsigned)i_entity_rank, <, new_sub_entity_nodes.size(), "new_sub_entity_nodes_check 1");
        VERIFY_OP_ON((unsigned)j_ordinal_of_entity, < , new_sub_entity_nodes[i_entity_rank].size(), "new_sub_entity_nodes_check 2");
        VERIFY_OP_ON((unsigned)k_ordinal_of_node_on_entity,  < , new_sub_entity_nodes[i_entity_rank][j_ordinal_of_entity].size(), "new_sub_entity_nodes_check 3");
      }
      catch (const std::exception& exc)
        {
          std::cout << "i_entity_rank= " << i_entity_rank << " j_ordinal_of_entity= " << j_ordinal_of_entity << " k_ordinal_of_node_on_entity= " << k_ordinal_of_node_on_entity << std::endl;
          throw exc;
        }
      return new_sub_entity_nodes[i_entity_rank][j_ordinal_of_entity][k_ordinal_of_node_on_entity];
    }

    inline int new_sub_entity_nodes_check_perm(PerceptMesh& eMesh, stk::mesh::Entity element, NewSubEntityNodesType& new_sub_entity_nodes, int i_entity_rank, int j_ordinal_of_entity, int k_ordinal_of_node_on_entity, const unsigned *perm)
    {
      VERIFY_OP_ON((unsigned)i_entity_rank, <, new_sub_entity_nodes.size(), "new_sub_entity_nodes_check 1");
      VERIFY_OP_ON((unsigned)j_ordinal_of_entity, < , new_sub_entity_nodes[i_entity_rank].size(), "new_sub_entity_nodes_check 2");
      VERIFY_OP_ON((unsigned)k_ordinal_of_node_on_entity,  < , new_sub_entity_nodes[i_entity_rank][j_ordinal_of_entity].size(), "new_sub_entity_nodes_check 3");

      return new_sub_entity_nodes[i_entity_rank][j_ordinal_of_entity][perm[k_ordinal_of_node_on_entity]];
    }



    /// The base class for all refinement patterns
    /// ------------------------------------------------------------------------------------------------------------------------
    //template< typename ToTopology >
    class UniformRefinerPatternBase
    {
    public:
      static const bool USE_DECLARE_ELEMENT_SIDE = false;

    protected:

      enum
        {
          /*
          base_topo_key_hex27      = shards::Hexahedron<27>::key,
          base_topo_key_hex20      = shards::Hexahedron<20>::key,
          base_topo_key_quad8      = shards::Quadrilateral<8>::key,
          base_topo_key_shellquad8 = shards::ShellQuadrilateral<8>::key,
          base_topo_key_shellquad9 = shards::ShellQuadrilateral<9>::key,
          base_topo_key_quad9      = shards::Quadrilateral<9>::key,
          base_topo_key_wedge15    = shards::Wedge<15>::key,
          */

          base_s_beam_2_key       = shards::Beam<2>::key,
          base_s_beam_3_key       = shards::Beam<3>::key,

          base_s_shell_line_2_key = shards::ShellLine<2>::key,
          base_s_shell_line_3_key = shards::ShellLine<3>::key,
          base_s_shell_tri_3_key  = shards::ShellTriangle<3>::key,
          base_s_shell_tri_6_key  = shards::ShellTriangle<6>::key,
          base_s_shell_quad_4_key = shards::ShellQuadrilateral<4>::key,
          base_s_shell_quad_8_key = shards::ShellQuadrilateral<8>::key,
          base_s_shell_quad_9_key = shards::ShellQuadrilateral<9>::key

        };

      stk::mesh::PartVector m_fromParts;
      stk::mesh::PartVector m_toParts;
      static const std::string m_appendConvertString; //="_urpconv_"
      const std::string m_convertSeparatorString; // temporarily use "#" then rename parts to use m_convertSeparatorFinalString
      const std::string m_convertSeparatorFinalString; // "_"
      const std::string m_appendOriginalString; //="_uo_1000"
      static const std::string m_oldElementsPartName;
      stk::mesh::EntityRank m_primaryEntityRank;

      static const unsigned topo_key_hex27      = shards::Hexahedron<27>::key;
      static const unsigned topo_key_hex20      = shards::Hexahedron<20>::key;
      static const unsigned topo_key_quad8      = shards::Quadrilateral<8>::key;
      static const unsigned topo_key_shellquad8 = shards::ShellQuadrilateral<8>::key;
      static const unsigned topo_key_shellquad9 = shards::ShellQuadrilateral<9>::key;
      static const unsigned topo_key_quad9      = shards::Quadrilateral<9>::key;
      static const unsigned topo_key_line3      = shards::Line<3>::key;
      static const unsigned topo_key_shellline3 = shards::ShellLine<3>::key;
      static const unsigned topo_key_wedge15    = shards::Wedge<15>::key;
      static const unsigned topo_key_pyramid13  = shards::Pyramid<13>::key;
      static const unsigned topo_key_pyramid5   = shards::Pyramid<5>::key;
      static const unsigned topo_key_tet4       = shards::Tetrahedron<4>::key;
      static const unsigned topo_key_tet10      = shards::Tetrahedron<10>::key;

      static const unsigned s_shell_line_2_key = shards::ShellLine<2>::key;
      static const unsigned s_shell_line_3_key = shards::ShellLine<3>::key;
      static const unsigned s_shell_tri_3_key  = shards::ShellTriangle<3>::key;
      static const unsigned s_shell_tri_6_key  = shards::ShellTriangle<6>::key;
      static const unsigned s_shell_quad_4_key = shards::ShellQuadrilateral<4>::key;
      static const unsigned s_shell_quad_8_key = shards::ShellQuadrilateral<8>::key;
      static const unsigned s_shell_quad_9_key = shards::ShellQuadrilateral<9>::key;

    public:
      bool m_mark_centroid_always;
      bool m_do_strip_hashes;
      vector<stk::mesh::Entity>::iterator m_ep_begin, m_ep_end;

      UniformRefinerPatternBase() : m_convertSeparatorString("."),
                                    m_convertSeparatorFinalString("_"),
                                    m_appendOriginalString(percept::PerceptMesh::s_omit_part+"_1000"),  // _100000
                                    m_primaryEntityRank(stk::topology::INVALID_RANK),
                                    m_mark_centroid_always(false),
                                    m_do_strip_hashes(true)
      {
        Elem::StdMeshObjTopologies::bootstrap();
      }
      virtual ~UniformRefinerPatternBase() {}

      virtual void doBreak()=0;

      virtual unsigned getFromTypeKey()=0;
      virtual unsigned getToTypeKey()=0;
      virtual const CellTopologyData *  getFromTopology()=0;
      virtual const CellTopologyData *  getToTopology()=0;

      // using mark information, attempt to discern actual number of elements needed
      virtual size_t estimateNumberOfNewElements(percept::PerceptMesh& eMesh, stk::mesh::EntityRank rank, NodeRegistry& nodeRegistry, size_t num_elem_not_ghost);

      // if this is not set, estimateNumberOfNewElements uses any mark of any entity, else only edge marks
      virtual bool edgeMarkIsEnough() { return true; }

      stk::mesh::EntityRank getPrimaryEntityRank() { return m_primaryEntityRank; }
      /// must be provided by derived classes
      /// ------------------------------------------------------------------------------------------------------------------------

      /// supplies the ranks of the sub entities needed during refinement (eg. m_eMesh.face_rank(), m_eMesh.edge_rank(),..)
      /// 10/02/10 and the number of nodes needed for each sub entity
      virtual void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)=0;

      ///
      void addToBreakPatternList(std::set<UniformRefinerPatternBase *>& list, PerceptMesh& eMesh);

      virtual void setNeededParts(PerceptMesh& eMesh, BlockNamesType block_names_ranks,
                                  bool sameTopology=true, bool skipConvertedParts=true);

      void fixSubsets(PerceptMesh& eMesh);
      void addExtraSurfaceParts(PerceptMesh& eMesh);

      void
      genericRefine_createNewElementsBase(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                                          stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                                          vector<stk::mesh::Entity>::iterator& ft_element_pool,

                                          const unsigned fromTopoKey_in, const unsigned toTopoKey_in,
                                          const int ToTopology_node_count,
                                          const int FromTopology_vertex_count,
                                          Elem::CellTopology elem_celltopo_in,
                                          RefTopoX_arr ref_topo_x_in,
                                          const CellTopologyData * const cell_topo_data_toTopo_in,
                                          vector< vector<stk::mesh::EntityId> >& elems,
                                          stk::mesh::FieldBase *proc_rank_field);

      enum { NumNewElements_Enrich = 1 };

#define EXTRA_PRINT_URP_IF 0

      void
      genericEnrich_createNewElementsBase(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                                          stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                                          vector<stk::mesh::Entity>::iterator& ft_element_pool,
                                          const unsigned fromTopoKey_in, const unsigned toTopoKey_in,
                                          const int ToTopology_node_count,
                                          const int FromTopology_vertex_count,
                                          Elem::CellTopology elem_celltopo_in,
                                          const CellTopologyData * const cell_topo_data_toTopo_in,
                                          vector< vector<stk::mesh::EntityId> >& elems,
                                          stk::mesh::FieldBase *proc_rank_field=0);

      /// helper for creating sides
      void create_side_element(PerceptMesh& eMesh, bool use_declare_element_side, stk::mesh::Entity *nodes, unsigned nodes_size, stk::mesh::Entity& newElement);

      /// helpers for interpolating fields, coordinates
      /// ------------------------------------------------------------------------------------------------------------------------

      /// This version uses Intrepid2 for interpolation
      void prolongateFields(percept::PerceptMesh& eMesh, stk::mesh::Entity element, stk::mesh::Entity newElement,  const unsigned *child_nodes,
                             RefTopoX_arr ref_topo_x, stk::mesh::FieldBase *field);

      /// do interpolation for all fields
      /// This version uses Intrepid2
      void prolongateFields(percept::PerceptMesh& eMesh, stk::mesh::Entity element, stk::mesh::Entity newElement, const unsigned *child_nodes,
                             RefTopoX_arr ref_topo_x);

      void prolongateIntrepid2(percept::PerceptMesh& eMesh, stk::mesh::FieldBase* field, shards::CellTopology& cell_topo,
                               MDArray& output_pts, stk::mesh::Entity element, MDArray& input_param_coords, double time_val=0.0);

      stk::mesh::Entity createOrGetNode(NodeRegistry& nodeRegistry, PerceptMesh& eMesh, stk::mesh::EntityId eid);

      void change_entity_parts(percept::PerceptMesh& eMesh, stk::mesh::Entity old_owning_elem, stk::mesh::Entity newElement);

      /*------------------------------------------------------------------------*/
      /*  comments from Shards_CellTopology.hpp with locally added comments
       * \brief  Find the permutation from the expected nodes to the actual nodes,
       *
       *  Find permutation 'p' such that:
       *    actual_node[j] == expected_node[ top.permutation[p].node[j] ]
       *  for all vertices.
       *
       *  So, actual_node[j] is the sub-dim cell; expected_node is the parent cell->subcell[dim][Ord].node[ perm[p][j] ]
       *
       *  Get sub-dim cell's nodes from NodeRegistry (actual_node array); or just sort them into a set
       *  Get parent element's sub-dim cell nodes from element->subcell[dim][ord].node
       *  Get permutation using shards::findPermutation(cell_topo, parent->subcell[dim][ord].node, subdim_cell_sorted_nodes)
       *
       *
       *  Then <b> ParentCell.node(K) == SubCell.node(I) </b> where:
       *  -  SubCellTopology == ParentCellTopology->subcell[dim][Ord].topology
       *  -  K  = ParentCellTopology->subcell[dim][Ord].node[IP]
       *  -  IP = SubCellTopology->permutation[P].node[I]
       *  -  I  = SubCellTopology->permutation_inverse[P].node[IP]

      */

      int getPermutation(PerceptMesh& eMesh, int num_verts, stk::mesh::Entity element, shards::CellTopology& cell_topo, unsigned rank_of_subcell, unsigned ordinal_of_subcell);

      /// supply the number of new elements per element during refinement
      virtual unsigned getNumNewElemPerElem()=0;

      /// given the node database (NodeRegistry), and the newly created nodes, and an iterator for the elements in the element pool,
      ///   create all new sub-elements of the refined element
      virtual void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& element_ft_pool,
                        stk::mesh::FieldBase *proc_rank_field = 0
                        ) = 0;

      /// if numChild is passed in as non-null, use that value, else use getNumNewElemPerElem() as size of child vector
      static void set_parent_child_relations(percept::PerceptMesh& eMesh, stk::mesh::Entity old_owning_elem, stk::mesh::Entity newElement,
                                      stk::mesh::Entity familyTreeNewElement,
                                      unsigned ordinal, unsigned *numChild=0);

      /// optionally overridden (must be overridden if sidesets are to work properly) to provide info on which sub pattern
      /// should be used to refine side sets (and edge sets)
      virtual void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh );

      /// can be overridden in cases where sub-patterns are not to be executed, but only used to
      ///   determine which parts need to be created (in hybrid refine patterns, for example, and
      ///   when a pattern is simply a container for other patterns).
      virtual void setSubPatternsForSetNeededParts( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh );

      /// for i/o to work properly, supply string replacements such as for hex-->tet breaking, you would supply "quad"-->"tri" etc. string maps
      virtual StringStringMap fixSurfaceAndEdgeSetNamesMap();

      /// provided by this class
      /// ------------------------------------------------------------------------------------------------------------------------

      virtual std::string getFromTopoPartName() =0;
      virtual std::string getToTopoPartName() =0;
      virtual std::string getName() { return std::string(); }

      stk::mesh::PartVector& getToParts() { return m_toParts; }
      stk::mesh::PartVector& getFromParts() { return m_fromParts; }
      static const std::string& getAppendConvertString() { return m_appendConvertString; }
      const std::string& getConvertSeparatorString() { return m_convertSeparatorString; }
      const std::string& getConvertSeparatorFinalString() { return m_convertSeparatorFinalString; }
      const std::string& getAppendOriginalString() { return m_appendOriginalString; }
      static const std::string& getOldElementsPartName() { return m_oldElementsPartName; }

      /// utilities
      /// ---------

      /// sets the needed number of nodes on each sub-entity to 1 - this is just a helper - in general, edges and faces have 1 new node
      /// for linear elements, and multiple new nodes in the case of quadratic elements
      void setToOne(std::vector<NeededEntityType>& needed_entities);
      double * midPoint(const double *p1, const double *p2, int spatialDim, double *x);
      double * getCentroid( double* pts[], int len, int spatialDim, double *x);
      static int getTopoDim(shards::CellTopology& cell_topo);

      static Teuchos::RCP<UniformRefinerPatternBase>
      createPattern(std::string refine, std::string enrich, std::string convert, percept::PerceptMesh& eMesh, BlockNamesType& block_names);


      static void mergeOrAddParts(UniformRefinerPatternBase *bp_from, UniformRefinerPatternBase *bp_to, bool merge);
      static void printParts(UniformRefinerPatternBase *bp, bool printAllParts = false);

      static std::string s_convert_options;
      static std::string s_refine_options;
      static std::string s_enrich_options;

      void updateSurfaceBlockMap(percept::PerceptMesh& eMesh, stk::mesh::Part* part, stk::mesh::Part* part_to);
    private:
      void addRefineNewNodesPart(percept::PerceptMesh& eMesh);
      void addActiveParentParts(percept::PerceptMesh& eMesh);
      bool foundIncludeOnlyBlock(percept::PerceptMesh& eMesh, std::vector<std::string>& block_names_include);
      void addOldPart(percept::PerceptMesh& eMesh);
      bool shouldDoThisPart(percept::PerceptMesh& eMesh, BlockNamesType block_names_ranks,
          		bool found_include_only_block, std::vector<std::string>& block_names_include, stk::mesh::Part *  part);
      void setNeededParts_debug1(percept::PerceptMesh& eMesh);
      void setNeededParts_debug2();
    };
    /// Utility intermediate base class providing more support for standard refinement operations
    /// ------------------------------------------------------------------------------------------------------------------------

    template< typename FTopo, typename TTopo  >
    class URP1
    {
    public:
      typedef FTopo FromTopology ;
      typedef TTopo ToTopology ;
    };

    template<typename FromTopology,  typename ToTopology >
    class URP :  public UniformRefinerPatternBase, public URP1<FromTopology, ToTopology>
    {
    public:

      static const unsigned fromTopoKey         = FromTopology::key;
      static const unsigned toTopoKey           = ToTopology::key;

      static const unsigned centroid_node       = (toTopoKey == topo_key_quad9 ? 8 :
                                                   (toTopoKey == topo_key_hex27 ? 20 : 0)
                                                   );

      // return the type of element this pattern can refine
      virtual unsigned getFromTypeKey() { return fromTopoKey; }
      virtual unsigned getToTypeKey() { return toTopoKey; }
      virtual const CellTopologyData * getFromTopology() { return shards::getCellTopologyData< FromTopology >(); }
      virtual const CellTopologyData * getToTopology() { return shards::getCellTopologyData< ToTopology >(); }

      virtual std::string getFromTopoPartName() {
        shards::CellTopology cell_topo(getFromTopology());
        return cell_topo.getName();
      }
      virtual std::string getToTopoPartName() {
        shards::CellTopology cell_topo(getToTopology());
        return cell_topo.getName();
      }

      virtual std::string getName() { return std::string("UniformRefinerPattern_")+getFromTopoPartName()+"_"+getToTopoPartName(); }

    protected:
      percept::PerceptMesh& m_eMesh;
      URP(percept::PerceptMesh& eMesh) : m_eMesh(eMesh) {}

      typedef ToTopology TTopo;
      typedef std::array<stk::mesh::EntityId, ToTopology::node_count > refined_element_type;



      enum { NumNewElements_Enrich = 1 };

      void
      genericEnrich_createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                                      stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes,
                                      vector<stk::mesh::Entity>::iterator& element_pool,
                                      vector<stk::mesh::Entity>::iterator& ft_element_pool,
                                      stk::mesh::FieldBase *proc_rank_field=0)
      {
        Elem::CellTopology elem_celltopo_in = Elem::getCellTopology< FromTopology >();
        const CellTopologyData * const cell_topo_data_toTopo_in = shards::getCellTopologyData< ToTopology >();

        static vector< vector<stk::mesh::EntityId> > elems;
        if (elems.size()==0)
          {
            elems.resize(NumNewElements_Enrich);
            for (unsigned in=0; in < NumNewElements_Enrich; in++)
              {
                elems[in].resize(ToTopology::node_count);
              }
          }
        genericEnrich_createNewElementsBase(eMesh, nodeRegistry, element, new_sub_entity_nodes, element_pool,
                                            ft_element_pool,
                                            fromTopoKey, toTopoKey, ToTopology::node_count, FromTopology::vertex_count,
                                            elem_celltopo_in,
                                            cell_topo_data_toTopo_in,
                                            elems,
                                            proc_rank_field);

      }
      void
      genericRefine_createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                                      stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                                      vector<stk::mesh::Entity>::iterator& ft_element_pool,

                                      stk::mesh::FieldBase *proc_rank_field=0)

      {
        Elem::CellTopology elem_celltopo_in = Elem::getCellTopology< FromTopology >();
        RefTopoX& ref_topo_x_in = Elem::StdMeshObjTopologies::RefinementTopologyExtra< FromTopology > ::refinement_topology;
        const CellTopologyData * const cell_topo_data_toTopo_in = shards::getCellTopologyData< ToTopology >();

        static std::vector<std::vector<stk::mesh::EntityId> > elems;
        if (elems.size()==0)
          {
            elems.resize(getNumNewElemPerElem());
            for (unsigned in=0; in < getNumNewElemPerElem(); in++)
              {
                elems[in].resize(ToTopology::node_count);
              }
          }
        genericRefine_createNewElementsBase(eMesh, nodeRegistry, element, new_sub_entity_nodes, element_pool,
                                            ft_element_pool,
                                            fromTopoKey, toTopoKey, ToTopology::node_count, FromTopology::vertex_count,
                                            elem_celltopo_in,
                                            &ref_topo_x_in[0],
                                            cell_topo_data_toTopo_in,
                                            elems,
                                            proc_rank_field);
      }

      /// utility methods for converting Sierra tables to new format (which groups DOF's on sub-entities)
      static bool on_parent_vertex(unsigned childNodeIdx);
      static bool on_parent_edge(unsigned childNodeIdx, const Elem::RefinementTopology& ref_topo);
      static bool on_parent_face(unsigned childNodeIdx, const Elem::RefinementTopology& ref_topo);
      static bool on_parent_edge_interior(unsigned childNodeIdx, const Elem::RefinementTopology& ref_topo, unsigned& i_edge, unsigned& i_ord, unsigned& n_ord);
      static unsigned renumber_quad_face_interior_nodes(unsigned original_node);

      // not used (yet)
      static unsigned renumber_quad_face_interior_nodes_quad8(unsigned original_node);
      static bool on_parent_face_interior(unsigned childNodeIdx, const Elem::RefinementTopology& ref_topo, unsigned& i_face, unsigned& i_ord, unsigned& n_ord);
      static bool on_parent_volume_interior(unsigned childNodeIdx, const Elem::RefinementTopology& ref_topo, unsigned& i_volume, unsigned& i_ord, unsigned& n_ord);

      /// utility to help convert Sierra tables - this method takes the index of the child node and finds it and adds
      ///   the associated info to the ref_topo_x tables containing the new/additional refinement table information
      static void findRefinedCellTopoInfo(unsigned childNodeIdx,
                                          const Elem::RefinementTopology& ref_topo,
                                          RefTopoX_arr ref_topo_x,  // assumed good for the vertices
                                          unsigned& rank_of_subcell,
                                          unsigned& ordinal_of_subcell,
                                          unsigned& ordinal_of_node_on_subcell,
                                          unsigned& num_node_on_subcell);
      static void findRefinedCellParamCoords(const Elem::RefinementTopology& ref_topo,
                                             RefTopoX_arr ref_topo_x);
      /// continuing in the convert tables theme, this helps to find the new nodes' parametric coordinates
      static void findRefinedCellParamCoordsLinear(const Elem::RefinementTopology& ref_topo,
                                                   RefTopoX_arr ref_topo_x  // assumed good for the vertices
                                                   );

    public:

      /// this is called one time (during code development) to generate and print a table of the extra refinement info

#define DEBUG_PRINT_REF_TOPO_X 1
      static void
      printRefinementTopoX_Table(std::ostream& out = std::cout );

    private:


    public:
      virtual ~URP() {}


    };

    template<typename FromTopology,  typename ToTopology, int NumToCreate, class OptionalTag=void>
    class UniformRefinerPattern : public URP<FromTopology, ToTopology> //, public URP1<FromTopology, ToTopology>
    {
    public:
    };

    // FIXME - temporary for testing only - we need a RefinerPattern with UniformRefinerPattern as a sub-class
    // (i.e. rename UniformRefinerPattern to RefinerPattern, create new UniformRefinerPattern as sub-class of RefinerPattern, specialize)
    template<typename FromTopology,  typename ToTopology, int NumToCreate, class OptionalTag=void>
    class RefinerPattern : public URP<FromTopology, ToTopology> //, public URP1<FromTopology, ToTopology>
    {
    public:
    };


  }

// all the patterns

// homogeneous refine
#include "UniformRefinerPattern_Quad4_Quad4_4.hpp"
#include "UniformRefinerPattern_Line2_Line2_2_sierra.hpp"
#include "UniformRefinerPattern_Beam2_Beam2_2_sierra.hpp"
#include "UniformRefinerPattern_ShellLine2_ShellLine2_2_sierra.hpp"
#include "UniformRefinerPattern_ShellLine3_ShellLine3_2_sierra.hpp"
#include "UniformRefinerPattern_Quad4_Quad4_4_sierra.hpp"
#include "UniformRefinerPattern_Tri3_Tri3_4_sierra.hpp"
#include "UniformRefinerPattern_ShellTri3_ShellTri3_4_sierra.hpp"
#include "UniformRefinerPattern_ShellTri6_ShellTri6_4_sierra.hpp"
#include "UniformRefinerPattern_ShellQuad4_ShellQuad4_4_sierra.hpp"
#include "UniformRefinerPattern_ShellQuad8_ShellQuad8_4_sierra.hpp"


#include "UniformRefinerPattern_Tet4_Tet4_8_sierra.hpp"
#include "UniformRefinerPattern_Hex8_Hex8_8_sierra.hpp"
#include "UniformRefinerPattern_Wedge6_Wedge6_8_sierra.hpp"
#include "UniformRefinerPattern_Pyramid5_Pyramid5_10_sierra.hpp"


#include "UniformRefinerPattern_Line3_Line3_2_sierra.hpp"
#include "UniformRefinerPattern_Beam3_Beam3_2_sierra.hpp"
#include "UniformRefinerPattern_Tri6_Tri6_4_sierra.hpp"
#include "UniformRefinerPattern_Quad8_Quad8_4_sierra.hpp"
#include "UniformRefinerPattern_Quad9_Quad9_4_sierra.hpp"
#include "UniformRefinerPattern_Hex27_Hex27_8_sierra.hpp"
#include "UniformRefinerPattern_Hex20_Hex20_8_sierra.hpp"
#include "UniformRefinerPattern_Tet10_Tet10_8_sierra.hpp"
#include "UniformRefinerPattern_Wedge15_Wedge15_8_sierra.hpp"
#include "UniformRefinerPattern_Wedge18_Wedge18_8_sierra.hpp"
#include "UniformRefinerPattern_Pyramid13_Pyramid13_10_sierra.hpp"

#include "URP_Heterogeneous_3D.hpp"
#include "URP_Heterogeneous_QuadraticRefine_3D.hpp"

// enrich

#include "UniformRefinerPattern_Line2_Line3_1_sierra.hpp"
#include "UniformRefinerPattern_ShellLine2_ShellLine3_1_sierra.hpp"
#include "UniformRefinerPattern_Beam2_Beam3_1_sierra.hpp"

#include "UniformRefinerPattern_Quad4_Quad9_1_sierra.hpp"
#include "UniformRefinerPattern_Quad4_Quad8_1_sierra.hpp"
#include "UniformRefinerPattern_ShellQuad4_ShellQuad8_1_sierra.hpp"
#include "UniformRefinerPattern_ShellQuad4_ShellQuad9_1_sierra.hpp"
#include "UniformRefinerPattern_Tri3_Tri6_1_sierra.hpp"
#include "UniformRefinerPattern_Tet4_Tet10_1_sierra.hpp"
#include "UniformRefinerPattern_Hex8_Hex27_1_sierra.hpp"
#include "UniformRefinerPattern_Hex8_Hex20_1_sierra.hpp"
#include "UniformRefinerPattern_Wedge6_Wedge15_1_sierra.hpp"
#include "UniformRefinerPattern_Wedge6_Wedge18_1_sierra.hpp"
#include "UniformRefinerPattern_Pyramid5_Pyramid13_1_sierra.hpp"

#include "URP_Heterogeneous_Enrich_3D.hpp"

// convert topology
#include "UniformRefinerPattern_Quad4_Tri3_6.hpp"
#include "UniformRefinerPattern_Quad4_Tri3_4.hpp"
#include "UniformRefinerPattern_Quad4_Tri3_2.hpp"
#include "UniformRefinerPattern_Hex8_Tet4_24.hpp"
#include "UniformRefinerPattern_Hex8_Tet4_6_12.hpp"
#include "UniformRefinerPattern_Pyramid5_Tet4_2.hpp"
#include "UniformRefinerPattern_Wedge6_Tet4_3.hpp"

#include "UniformRefinerPattern_Tri3_Quad4_3.hpp"
#include "UniformRefinerPattern_Tet4_Hex8_4.hpp"
#include "UniformRefinerPattern_Wedge6_Hex8_6.hpp"
#include "URP_Tet4_Wedge6_Hex8.hpp"
#include "URP_Hex8_Wedge6_Pyramid5_Tet4.hpp"

// local refinement
#include "RefinerPattern_Tri3_Tri3_2.hpp"
#include "RefinerPattern_Tri3_Tri3_N.hpp"
#include "RefinerPattern_Tet4_Tet4_N.hpp"

#include "RefinerPattern_Tri3_Tri3_HangingNode.hpp"
#include "RefinerPattern_Tet4_Tet4_HangingNode.hpp"
#include "RefinerPattern_Quad4_Quad4_HangingNode.hpp"
#include "RefinerPattern_Quad4_Quad4_Transition.hpp"
#include "RefinerPattern_Hex8_Hex8_HangingNode.hpp"
#include "RefinerPattern_Hex8_Hex8_Transition.hpp"
#include "RefinerPattern_Pyr5_Pyr5_Transition.hpp"
#include "RefinerPattern_Wedge6_Wedge6_Transition.hpp"
#include "RefinerPattern_Quad4_Tri3_Hybrid_Transition.hpp"
#include "RefinerPattern_Hybrid_3D_Transition.hpp"

  namespace percept {

    // refine
    typedef  UniformRefinerPattern<shards::Line<2>,          shards::Line<2>,          2, SierraPort >            Line2_Line2_2;
    typedef  UniformRefinerPattern<shards::Beam<2>,          shards::Beam<2>,          2, SierraPort >            Beam2_Beam2_2;
    typedef  UniformRefinerPattern<shards::ShellLine<2>,     shards::ShellLine<2>,     2, SierraPort >            ShellLine2_ShellLine2_2;
    typedef  UniformRefinerPattern<shards::ShellLine<3>,     shards::ShellLine<3>,     2, SierraPort >            ShellLine3_ShellLine3_2;
    typedef  UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4 >                        Quad4_Quad4_4_Old;
    typedef  UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort >            Quad4_Quad4_4;

    typedef  UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, 4, SierraPort >            Quad4_Quad4_4_Sierra;
    typedef  UniformRefinerPattern<shards::Triangle<3>,      shards::Triangle<3>,      4, SierraPort >            Tri3_Tri3_4;
    typedef  UniformRefinerPattern<shards::ShellTriangle<3>, shards::ShellTriangle<3>, 4, SierraPort >            ShellTri3_ShellTri3_4;
    typedef  UniformRefinerPattern<shards::ShellTriangle<6>, shards::ShellTriangle<6>, 4, SierraPort >            ShellTri6_ShellTri6_4;
    typedef  UniformRefinerPattern<shards::ShellQuadrilateral<4>, shards::ShellQuadrilateral<4>, 4, SierraPort >  ShellQuad4_ShellQuad4_4;
    typedef  UniformRefinerPattern<shards::ShellQuadrilateral<8>, shards::ShellQuadrilateral<8>, 4, SierraPort >  ShellQuad8_ShellQuad8_4;

    typedef  UniformRefinerPattern<shards::Tetrahedron<4>,   shards::Tetrahedron<4>,   8, SierraPort >            Tet4_Tet4_8;
    typedef  UniformRefinerPattern<shards::Hexahedron<8>,    shards::Hexahedron<8>,    8, SierraPort >            Hex8_Hex8_8;
    typedef  UniformRefinerPattern<shards::Wedge<6>,         shards::Wedge<6>,         8, SierraPort >            Wedge6_Wedge6_8;
    typedef  UniformRefinerPattern<shards::Pyramid<5>,       shards::Pyramid<5>,      10, SierraPort >            Pyramid5_Pyramid5_10;

    typedef  UniformRefinerPattern<shards::Line<3>,          shards::Line<3>,          2, SierraPort >            Line3_Line3_2;
    typedef  UniformRefinerPattern<shards::Beam<3>,          shards::Beam<3>,          2, SierraPort >            Beam3_Beam3_2;
    typedef  UniformRefinerPattern<shards::Triangle<6>,      shards::Triangle<6>,      4, SierraPort >            Tri6_Tri6_4;
    typedef  UniformRefinerPattern<shards::Quadrilateral<9>, shards::Quadrilateral<9>, 4, SierraPort >            Quad9_Quad9_4;
    typedef  UniformRefinerPattern<shards::Quadrilateral<8>, shards::Quadrilateral<8>, 4, SierraPort >            Quad8_Quad8_4;
    typedef  UniformRefinerPattern<shards::Hexahedron<27>,   shards::Hexahedron<27>,   8, SierraPort >            Hex27_Hex27_8;
    typedef  UniformRefinerPattern<shards::Hexahedron<20>,   shards::Hexahedron<20>,   8, SierraPort >            Hex20_Hex20_8;
    typedef  UniformRefinerPattern<shards::Tetrahedron<10>,  shards::Tetrahedron<10>,  8, SierraPort >            Tet10_Tet10_8;
    typedef  UniformRefinerPattern<shards::Wedge<15>,        shards::Wedge<15>,        8, SierraPort >            Wedge15_Wedge15_8;
    typedef  UniformRefinerPattern<shards::Wedge<18>,        shards::Wedge<18>,        8, SierraPort >            Wedge18_Wedge18_8;
    typedef  UniformRefinerPattern<shards::Pyramid<13>,      shards::Pyramid<13>,     10, SierraPort >            Pyramid13_Pyramid13_10;

    // enrich
    typedef  UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<9>, 1, SierraPort >            Quad4_Quad9_1;
    typedef  UniformRefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<8>, 1, SierraPort >            Quad4_Quad8_1;
    typedef  UniformRefinerPattern<shards::Beam<2>,          shards::Beam<3>,          1, SierraPort >            Beam2_Beam3_1;

    typedef  UniformRefinerPattern<shards::ShellQuadrilateral<4>, shards::ShellQuadrilateral<8>, 1, SierraPort >  ShellQuad4_ShellQuad8_1;
    typedef  UniformRefinerPattern<shards::ShellQuadrilateral<4>, shards::ShellQuadrilateral<9>, 1, SierraPort >  ShellQuad4_ShellQuad9_1;
    typedef  UniformRefinerPattern<shards::Triangle<3>,      shards::Triangle<6>,      1, SierraPort >            Tri3_Tri6_1;
    typedef  UniformRefinerPattern<shards::Tetrahedron<4>,   shards::Tetrahedron<10>,  1, SierraPort >            Tet4_Tet10_1;
    typedef  UniformRefinerPattern<shards::Hexahedron<8>,    shards::Hexahedron<27>,   1, SierraPort >            Hex8_Hex27_1;
    typedef  UniformRefinerPattern<shards::Hexahedron<8>,    shards::Hexahedron<20>,   1, SierraPort >            Hex8_Hex20_1;
    typedef  UniformRefinerPattern<shards::Wedge<6>,         shards::Wedge<15>,        1, SierraPort >            Wedge6_Wedge15_1;
    typedef  UniformRefinerPattern<shards::Wedge<6>,         shards::Wedge<18>,        1, SierraPort >            Wedge6_Wedge18_1;
    typedef  UniformRefinerPattern<shards::Pyramid<5>,       shards::Pyramid<13>,      1, SierraPort >            Pyramid5_Pyramid13_1;

    // convert
    typedef  UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>,      2 >                        Quad4_Tri3_2;
    typedef  UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>,      4, Specialization >        Quad4_Tri3_4;
    typedef  UniformRefinerPattern<shards::Quadrilateral<4>, shards::Triangle<3>,      6 >                        Quad4_Tri3_6;
    typedef  UniformRefinerPattern<shards::Hexahedron<8>,    shards::Tetrahedron<4>,  24 >                        Hex8_Tet4_24;
    typedef  UniformRefinerPattern<shards::Hexahedron<8>,    shards::Tetrahedron<4>,   6 >                        Hex8_Tet4_6_12;
    typedef  UniformRefinerPattern<shards::Pyramid<5>,       shards::Tetrahedron<4>,   2 >                        Pyramid5_Tet4_2;
    typedef  UniformRefinerPattern<shards::Wedge<6>,         shards::Tetrahedron<4>,   3 >                        Wedge6_Tet4_3;

    typedef  UniformRefinerPattern<shards::Triangle<3>,      shards::Quadrilateral<4>, 3, Specialization >        Tri3_Quad4_3;
    typedef  UniformRefinerPattern<shards::Tetrahedron<4>,   shards::Hexahedron<8>,    4 >                        Tet4_Hex8_4;
    typedef  UniformRefinerPattern<shards::Wedge<6>,         shards::Hexahedron<8>,    6 >                        Wedge6_Hex8_6;
    typedef  URP_Tet4_Wedge6_Hex8 Tet4_Wedge6_Hex8;
    typedef  URP_Hex8_Wedge6_Pyramid5_Tet4 Hex8_Wedge6_Pyramid5_Tet4;

    // local refinement - for testing only right now

    typedef  RefinerPattern<shards::Triangle<3>,      shards::Triangle<3>,      2  >            Local_Tri3_Tri3_2;
    typedef  RefinerPattern<shards::Triangle<3>,      shards::Triangle<3>,     -1  >            Local_Tri3_Tri3_N;

    typedef  RefinerPattern<shards::Tetrahedron<4>,   shards::Tetrahedron<4>,  -1  >            Local_Tet4_Tet4_N;

    typedef  RefinerPattern<shards::Hexahedron<8>,    shards::Hexahedron<8>,  -1  >                Local_Hex8_Hex8_N;
    typedef  RefinerPattern<shards::Hexahedron<8>,    shards::Hexahedron<8>,  -1 , HexTransition > Local_Hex8_Hex8_N_Transition;

    typedef  RefinerPattern<shards::Triangle<3>,      shards::Triangle<3>,-1, TriHangingNode >  Local_Tri3_Tri3_N_HangingNode;
    typedef  RefinerPattern<shards::Tetrahedron<4>,   shards::Tetrahedron<4>,-1, TetHangingNode >  Local_Tet4_Tet4_N_HangingNode;

    typedef  RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>,-1, QuadHangingNode >  Local_Quad4_Quad4_N_HangingNode;
    typedef  RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1, QuadTransition >  Local_Quad4_Quad4_N_Transition;
    typedef  RefinerPattern_Quad4_Tri3_Hybrid_Transition Local_Quad4_Tri3_Hybrid_Transition;

    typedef  RefinerPattern<shards::Pyramid<5>,    shards::Pyramid<5>,  -1 , PyrTransition >          Local_Pyr5_Pyr5_N_Transition;
    typedef  RefinerPattern<shards::Wedge<6>,    shards::Wedge<6>,  -1 , WedgeTransition >          Local_Wedge6_Wedge6_N_Transition;

    typedef  RefinerPattern_Hybrid_3D_Transition                                                    Local_Hybrid_3D;

    //DPM adding enum for use in initializing empty UniformRefiner objects
    enum Pattern
      {
        LINE2_LINE2_2,
        BEAM2_BEAM2_2,
        SHELLLINE2_SHELLLINE2_2,
        SHELLLINE3_SHELLLINE3_2,
        QUAD4_QUAD4_4_OLD,
        QUAD4_QUAD4_4,
        QUAD4_QUAD4_4_SIERRA,
        TRI3_TRI3_4,
        SHELLTRI3_SHELLTRI3_4,
        SHELLTRI6_SHELLTRI6_4,
        SHELLQUAD4_SHELLQUAD4_4,
        SHELLQUAD8_SHELLQUAD8_4,
        TET4_TET4_8,
        HEX8_HEX8_8,
        WEDGE6_WEDGE6_8,
        PYRAMID5_PYRAMID5_10,
        LINE3_LINE3_2,
        BEAM3_BEAM3_2,
        TRI6_TRI6_4,
        QUAD9_QUAD9_4,
        QUAD8_QUAD8_4,
        HEX27_HEX27_8,
        HEX20_HEX20_8,
        TET10_TET10_8,
        WEDGE15_WEDGE15_8,
        WEDGE18_WEDGE18_8,
        PYRAMID13_PYRAMID13_10,
        QUAD4_QUAD9_1,
        QUAD4_QUAD8_1,
        BEAM2_BEAM3_1,
        SHELLQUAD4_SHELLQUAD8_1,
        TRI3_TRI6_1,
        TET4_TET10_1,
        HEX8_HEX27_1,
        HEX8_HEX20_1,
        WEDGE6_WEDGE15_1,
        WEDGE6_WEDGE18_1,
        PYRAMID5_PYRAMID13_1,

        QUAD4_TRI3_2,
        QUAD4_TRI3_4,
        QUAD4_TRI3_6,
        HEX8_TET4_24,
        HEX8_TET4_6_12,

        TRI3_QUAD4_3,
        TET4_HEX8_4,
        WEDGE6_HEX8_6,

        NUM_REF_PATTERNS
      };

  }

#endif
