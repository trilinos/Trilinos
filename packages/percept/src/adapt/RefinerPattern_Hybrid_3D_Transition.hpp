// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinerPattern_Hybrid_3D_Transition_hpp
#define adapt_RefinerPattern_Hybrid_3D_Transition_hpp

#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "RefinerPattern_Hex8_Hex8_Transition.hpp"
#include "RefinerPattern_Pyr5_Pyr5_Transition.hpp"
#include "RefinerPattern_Wedge6_Wedge6_Transition.hpp"
#include "RefinerPattern_Tet4_Tet4_HangingNode.hpp"

#include "RefinerPattern_Quad4_Tri3_Hybrid_Transition.hpp"

#include "RefinerPattern_Quad4_Quad4_Transition.hpp"
#include "RefinerPattern_Tri3_Tri3_HangingNode.hpp"
#include "RefinerPattern_Quad4_Tri3_Hybrid_Transition.hpp"

namespace percept {

  //struct _Hybrid_Transition {};

#define LOCAL_DEBUG_PP 0

  class RefinerPattern_Hybrid_3D_Transition : public UniformRefinerPatternBase
  {
  private:
    std::vector<UniformRefinerPatternBase *> m_bp;

    unsigned m_final_index;

  protected:

    percept::PerceptMesh& m_eMesh;

  public:

    RefinerPattern_Hybrid_3D_Transition(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : m_final_index(0), m_eMesh(eMesh)
    {
      m_primaryEntityRank = stk::topology::ELEMENT_RANK;

      //setNeededParts(eMesh, block_names, true);
      Elem::StdMeshObjTopologies::bootstrap();
      m_bp.resize(11);

      m_bp[0] = new  RefinerPattern<shards::Hexahedron<8>,    shards::Hexahedron<8>,  -1 , HexTransition >    (eMesh, block_names) ;
      m_bp[1] = new  RefinerPattern<shards::Pyramid<5>,       shards::Pyramid<5>,     -1 , PyrTransition >    (eMesh, block_names) ;
      m_bp[2] = new  RefinerPattern<shards::Wedge<6>,         shards::Wedge<6>,       -1 , WedgeTransition >  (eMesh, block_names) ;
      m_bp[3] = new  RefinerPattern<shards::Tetrahedron<4>,   shards::Tetrahedron<4>, -1,  TetHangingNode >   (eMesh, block_names) ;

      //  -----
      //m_bp[3]->m_mark_centroid_always = true;
      //static_cast<RefinerPattern<shards::Tetrahedron<4>,   shards::Tetrahedron<4>, -1,  TetHangingNode > * > ( m_bp[3])->m_transition_breaker->m_mark_centroid_always = true;
      //  -----

      m_bp[4] = new  RefinerPattern<shards::Hexahedron<8>,    shards::Hexahedron<8>,  -1,  HexTempPartialNoBreak >   (eMesh, block_names) ;
      m_bp[5] = new  RefinerPattern<shards::Tetrahedron<4>,   shards::Tetrahedron<4>, -1,  TetTempPartialNoBreak >   (eMesh, block_names) ;

      //  -----
      //m_bp[5]->m_mark_centroid_always = true;
      //  -----

      m_bp[6] = new  RefinerPattern<shards::Pyramid<5>,       shards::Pyramid<5>,     -1,  PyrTempPartialNoBreak >   (eMesh, block_names) ;
      m_bp[7] = new  RefinerPattern<shards::Wedge<6>,         shards::Wedge<6>,       -1,  WedgeTempPartialNoBreak > (eMesh, block_names) ;

      m_final_index = 7;

      //m_bp.push_back(new  RefinerPattern_Quad4_Tri3_Hybrid_Transition (eMesh, block_names));
      m_bp[8] = new  RefinerPattern<shards::Triangle<3>,      shards::Triangle<3>,      -1, TriHangingNode >        (eMesh, block_names) ;
      bool use_only_tris = false;
      bool avoid_centroid_node = true;
      m_bp[9] = new  RefinerPattern<shards::Quadrilateral<4>, shards::Quadrilateral<4>, -1, QuadTransition >        (eMesh, block_names, use_only_tris, avoid_centroid_node) ;
      m_bp[10] = new  RefinerPattern<shards::Triangle<3>,      shards::Triangle<3>,      -1, TriTempPartialNoBreak > (eMesh, block_names) ;

      // bool skipConvertedParts = false;
      // setNeededParts(eMesh, block_names, true, skipConvertedParts);
      // skipConvertedParts = true;
      // setNeededParts(eMesh, block_names, true, skipConvertedParts);

      typedef std::set<UniformRefinerPatternBase *> List;
      List list;
      addToBreakPatternList(list, eMesh);
      bool skipConvertedParts = true;
      if (LOCAL_DEBUG_PP) std::cout << "list size= " << list.size() << std::endl;
      for (List::iterator it = list.begin(); it != list.end(); ++it)
        {
          UniformRefinerPatternBase *bp = *it;
          if (LOCAL_DEBUG_PP) std::cout << "list= " << eMesh.demangle(typeid(*bp).name()) << std::endl;
          bool sameTopology = true;
          if (bp != this)
            {
              sameTopology = (bp->getFromTypeKey() == bp->getToTypeKey());
              bp->setNeededParts(eMesh, block_names, sameTopology, skipConvertedParts);
            }
        }
      skipConvertedParts = false;
      for (List::iterator it = list.begin(); it != list.end(); ++it)
        {
          UniformRefinerPatternBase *bp = *it;
          bool sameTopology = true;
          if (bp != this)
            {
              sameTopology = (bp->getFromTypeKey() == bp->getToTypeKey());
              bp->setNeededParts(eMesh, block_names, sameTopology, skipConvertedParts);
            }
        }

      // for (List::iterator it = list.begin(); it != list.end(); ++it)
      //   {
      //     UniformRefinerPatternBase *bp = *it;
      for (unsigned ibp=0; ibp < m_bp.size(); ibp++)
        {
          bool mergeParts = false;
          mergeOrAddParts(m_bp[ibp], this, mergeParts);
        }

      if (LOCAL_DEBUG_PP)
        {
          std::cout << "tmp Hybrid____ printParts this = \n" ;
          printParts(this);
        }
      m_mark_centroid_always = true;

    }

    ~RefinerPattern_Hybrid_3D_Transition()
    {
      for (unsigned ibp=0; ibp < m_bp.size(); ibp++)
        {
          if (m_bp[ibp] ) delete m_bp[ibp];
        }
    }

    void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
    {
      bp = m_bp;
    }

    virtual void doBreak() override
    {
      throw std::runtime_error("shouldn't call RefinerPattern_Hybrid_3D_Transition::doBreak()");

    }
    virtual unsigned getFromTypeKey() override
    {
      throw std::runtime_error("shouldn't call RefinerPattern_Hybrid_3D_Transition::getFromTypeKey()");

    }
    virtual unsigned getToTypeKey() override
    {
      throw std::runtime_error("shouldn't call RefinerPattern_Hybrid_3D_Transition::getToTypeKey()");

    }

    virtual std::string getFromTopoPartName() override {
      shards::CellTopology cell_topo(getFromTopology());
      return cell_topo.getName();
    }
    virtual std::string getToTopoPartName() override {
      shards::CellTopology cell_topo(getToTopology());
      return cell_topo.getName();
    }

    virtual const CellTopologyData *  getFromTopology() override
    {
      throw std::runtime_error("shouldn't call RefinerPattern_Hybrid_3D_Transition::getFromTopology()");
    }

    virtual const CellTopologyData *  getToTopology() override {
      throw std::runtime_error("shouldn't call RefinerPattern_Hybrid_3D_Transition::getToTopology()");
    }

    void fillNeededEntities(std::vector<NeededEntityType>& /*needed_entities*/) override
    {
      throw std::runtime_error("shouldn't call RefinerPattern_Hybrid_3D_Transition::fillNeededEntities()");
    }

    virtual unsigned getNumNewElemPerElem() override
    {
      throw std::runtime_error("shouldn't call RefinerPattern_Hybrid_3D_Transition::getNumNewElemPerElem()");
      return 8;
    }

    void
    createNewElements(percept::PerceptMesh& /*eMesh*/, NodeRegistry& /*nodeRegistry*/,
                      stk::mesh::Entity /*element*/,  NewSubEntityNodesType& /*new_sub_entity_nodes*/, vector<stk::mesh::Entity>::iterator& /*element_pool*/,
                      vector<stk::mesh::Entity>::iterator& /*ft_element_pool*/,
                      stk::mesh::FieldBase */*proc_rank_field*/=0) override
    {
      throw std::runtime_error("shouldn't call RefinerPattern_Hybrid_3D_Transition::createNewElements()");
    }

  };

}
#undef LOCAL_DEBUG_PP

#endif
