// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinerPattern_Quad4_Tri3_Hybrid_Transition_hpp
#define adapt_RefinerPattern_Quad4_Tri3_Hybrid_Transition_hpp

#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "RefinerPattern_Quad4_Quad4_Transition.hpp"
#include "RefinerPattern_Tri3_Tri3_HangingNode.hpp"

namespace percept {

  //struct Quad4_Tri3_Hybrid_Transition {};
#define PRINT_PARTS 0

  class RefinerPattern_Quad4_Tri3_Hybrid_Transition : public UniformRefinerPatternBase
  {
  private:
    std::vector<UniformRefinerPatternBase *> m_bp;
    std::vector<bool> m_bp_owned;

  protected:

    percept::PerceptMesh& m_eMesh;

  public:

    RefinerPattern_Quad4_Tri3_Hybrid_Transition(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : m_eMesh(eMesh)
    {
      m_primaryEntityRank = eMesh.face_rank();
      if (m_eMesh.get_spatial_dim() == 2)
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

      //setNeededParts(eMesh, block_names, true);
      Elem::StdMeshObjTopologies::bootstrap();
      m_bp.resize(0);
      m_bp_owned.resize(0);

      typedef  RefinerPattern<shards::Quadrilateral<4>,   shards::Quadrilateral<4>,  -1, QuadTransition > QuadTP;
      QuadTP *quadP = new QuadTP(eMesh, block_names);
      m_bp.push_back(quadP);
      m_bp_owned.push_back(true);
      m_bp.push_back(new  RefinerPattern<shards::Triangle<3>,      shards::Triangle<3>,-1, TriHangingNode > (eMesh, block_names));
      m_bp_owned.push_back(true);

      // next pattern is: RefinerPattern<shards::Triangle<3>, shards::Triangle<3>,     -1, TriTempPartialNoBreak > 
      //m_bp.push_back(quadP->m_transition_breaker->m_bp[2]); // FIXME magic number
      for (unsigned ii=0; ii < quadP->m_transition_breaker->m_bp_exported.size(); ++ii)
        {
          m_bp.push_back(quadP->m_transition_breaker->m_bp_exported[ii]);
        }

      m_bp_owned.push_back(false);

      m_bp.push_back(new  RefinerPattern<shards::Line<2>, shards::Line<2>, -1 > (eMesh, block_names));
      m_bp_owned.push_back(true);

      if (PRINT_PARTS)
        {
          std::cout << "tmp printParts this=\n" ;
          printParts(this, false);
          for (unsigned ibp=0; ibp < m_bp.size(); ibp++)
            {
              std::cout << "tmp printParts this after=\n" ;
              printParts(m_bp[ibp], false);
            }
        }

      if (m_eMesh.get_spatial_dim() == 2)
        m_mark_centroid_always = true;

    }

    ~RefinerPattern_Quad4_Tri3_Hybrid_Transition()
    {
      for (unsigned ibp=0; ibp < m_bp.size(); ibp++)
        {
          if (m_bp[ibp] && m_bp_owned[ibp]) delete m_bp[ibp];
        }
    }

    void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
    {
      bp = m_bp;
    }

    virtual void doBreak() override
    {
      throw std::runtime_error("shouldn't call RefinerPattern_Quad4_Tri3_Hybrid_Transition::doBreak()");

    }
    virtual unsigned getFromTypeKey() override
    {
      throw std::runtime_error("shouldn't call RefinerPattern_Quad4_Tri3_Hybrid_Transition::getFromTypeKey()");

    }
    virtual unsigned getToTypeKey() override
    {
      throw std::runtime_error("shouldn't call RefinerPattern_Quad4_Tri3_Hybrid_Transition::getToTypeKey()");

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
      throw std::runtime_error("shouldn't call RefinerPattern_Quad4_Tri3_Hybrid_Transition::getFromTopology()");
    }

    virtual const CellTopologyData *  getToTopology() override {
      throw std::runtime_error("shouldn't call RefinerPattern_Quad4_Tri3_Hybrid_Transition::getToTopology()");
    }

    void fillNeededEntities(std::vector<NeededEntityType>& /*needed_entities*/) override
    {
      throw std::runtime_error("shouldn't call RefinerPattern_Quad4_Tri3_Hybrid_Transition::fillNeededEntities()");
    }

    virtual unsigned getNumNewElemPerElem() override
    {
      throw std::runtime_error("shouldn't call RefinerPattern_Quad4_Tri3_Hybrid_Transition::getNumNewElemPerElem()");
      return 8;
    }

    void
    createNewElements(percept::PerceptMesh& /*eMesh*/, NodeRegistry& /*nodeRegistry*/,
                      stk::mesh::Entity /*element*/,  NewSubEntityNodesType& /*new_sub_entity_nodes*/, vector<stk::mesh::Entity>::iterator& /*element_pool*/,
                      vector<stk::mesh::Entity>::iterator& /*ft_element_pool*/,
                      stk::mesh::FieldBase */*proc_rank_field*/=0) override
    {
      throw std::runtime_error("shouldn't call RefinerPattern_Quad4_Tri3_Hybrid_Transition::createNewElements()");
    }

  };

}
#undef PRINT_PARTS

#endif
