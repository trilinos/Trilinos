// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinerPattern_NoBreak_hpp
#define adapt_RefinerPattern_NoBreak_hpp

#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#include <percept/PerceptBoostArray.hpp>

namespace percept {

  //================================================================================================================================================================
  // Just for ensuring tets are processed by the marker - this pattern does nothing, it doesn't break tets...
  //================================================================================================================================================================

  struct TetTempPartialNoBreak {};

  // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
  template <>
  class RefinerPattern<shards::Tetrahedron<4>, shards::Tetrahedron<4>, -1, TetTempPartialNoBreak > : public URP<shards::Tetrahedron<4>, shards::Tetrahedron<4>  >
  {
  public:

    RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Tetrahedron<4>, shards::Tetrahedron<4>  >(eMesh)
    {
      m_primaryEntityRank = stk::topology::ELEMENT_RANK;

      bool sameTopology = true;
      setNeededParts(eMesh, block_names, sameTopology);
      //getToParts().resize(0);
      Elem::StdMeshObjTopologies::bootstrap();
    }

    ~RefinerPattern()
    {
    }

    void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
    {
      bp.resize(1); // = std::vector<UniformRefinerPatternBase *>(1u, 0);
      bp[0] = this;
    }

    virtual void doBreak() override {}
    void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
    {
      if (!m_mark_centroid_always)
        {
          needed_entities.resize(1);
          needed_entities[0].first = m_eMesh.edge_rank();
          //needed_entities[1].first = m_eMesh.face_rank();
          setToOne(needed_entities);
        }
      else
        {
          needed_entities.resize(2);
          needed_entities[0].first = m_eMesh.edge_rank();
          needed_entities[1].first = m_eMesh.face_rank();
          setToOne(needed_entities);
        }
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
  // Just for ensuring pyramids are processed by the marker - this pattern does nothing, it doesn't break pyramids...
  //================================================================================================================================================================

  struct PyrTempPartialNoBreak {};

  // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
  template <>
  class RefinerPattern<shards::Pyramid<5>, shards::Pyramid<5>, -1, PyrTempPartialNoBreak > : public URP<shards::Pyramid<5>, shards::Pyramid<5>  >
  {
  public:

    RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Pyramid<5>, shards::Pyramid<5>  >(eMesh)
    {
      m_primaryEntityRank = stk::topology::ELEMENT_RANK;

      bool sameTopology = true;
      setNeededParts(eMesh, block_names, sameTopology);
      //getToParts().resize(0);
      Elem::StdMeshObjTopologies::bootstrap();

    }

    ~RefinerPattern()
    {
    }

    void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
    {
      bp.resize(1); // = std::vector<UniformRefinerPatternBase *>(1u, 0);
      bp[0] = this;
    }

    virtual void doBreak() override {}
    void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
    {
      needed_entities.resize(2);
      needed_entities[0].first = m_eMesh.edge_rank();
      needed_entities[1].first = m_eMesh.face_rank();
      setToOne(needed_entities);

      int faces[5] = {0,0,0,0,1};
      needed_entities[1].third.resize(5);
      needed_entities[1].third.assign(faces,faces+5);

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
  // Just for ensuring wedges are processed by the marker - this pattern does nothing, it doesn't break wedges...
  //================================================================================================================================================================

  struct WedgeTempPartialNoBreak {};

  // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
  template <>
  class RefinerPattern<shards::Wedge<6>, shards::Wedge<6>, -1, WedgeTempPartialNoBreak > : public URP<shards::Wedge<6>, shards::Wedge<6>  >
  {
  public:

    RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Wedge<6>, shards::Wedge<6>  >(eMesh)
    {
      m_primaryEntityRank = stk::topology::ELEMENT_RANK;

      bool sameTopology = true;
      setNeededParts(eMesh, block_names, sameTopology);
      //getToParts().resize(0);
      Elem::StdMeshObjTopologies::bootstrap();

    }

    ~RefinerPattern()
    {
    }

    void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
    {
      bp.resize(1); // = std::vector<UniformRefinerPatternBase *>(1u, 0);
      bp[0] = this;
    }

    virtual void doBreak() override {}
    void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
    {
      needed_entities.resize(2);
      needed_entities[0].first = m_eMesh.edge_rank();
      needed_entities[1].first = m_eMesh.face_rank();
      setToOne(needed_entities);

      int faces[5] = {1,1,1,0,0};
      needed_entities[1].third.resize(5);
      needed_entities[1].third.assign(faces,faces+5);

    }

    // FIXME - for now, create more than we need (to fix this right we need a change to the Refiner.cpp interface)
    virtual unsigned getNumNewElemPerElem() override { return 6; }

    void
    createNewElements(percept::PerceptMesh& /*eMesh*/, NodeRegistry& /*nodeRegistry*/,
                      stk::mesh::Entity /*element*/,  NewSubEntityNodesType& /*new_sub_entity_nodes*/, vector<stk::mesh::Entity>::iterator& /*element_pool*/,
                      vector<stk::mesh::Entity>::iterator& /*ft_element_pool*/,
                      stk::mesh::FieldBase */*proc_rank_field*/=0) override
    {
    }

  };

  //================================================================================================================================================================
  // Just for ensuring tets are processed by the marker - this pattern does nothing, it doesn't break tets...
  //================================================================================================================================================================

  struct HexTempPartialNoBreak {};

  // the "-1" here signifies the number of elements created is not fixed, depends on the marking pattern
  template <>
  class RefinerPattern<shards::Hexahedron<8>, shards::Hexahedron<8>, -1, HexTempPartialNoBreak > : public URP<shards::Hexahedron<8>, shards::Hexahedron<8>  >
  {
  public:

    RefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Hexahedron<8>, shards::Hexahedron<8>  >(eMesh)
    {
      m_primaryEntityRank = stk::topology::ELEMENT_RANK;

      bool sameTopology = true;
      setNeededParts(eMesh, block_names, sameTopology);
      //getToParts().resize(0);
      Elem::StdMeshObjTopologies::bootstrap();

    }

    ~RefinerPattern()
    {
    }

    void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& /*eMesh*/ ) override
    {
      bp.resize(1); // = std::vector<UniformRefinerPatternBase *>(1u, 0);
      bp[0] = this;
    }

    virtual void doBreak() override {}
    void fillNeededEntities(std::vector<NeededEntityType>& needed_entities) override
    {
      needed_entities.resize(3);
      needed_entities[0].first = m_eMesh.edge_rank();
      needed_entities[1].first = m_eMesh.face_rank();
      needed_entities[2].first = m_eMesh.element_rank();
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


}

#endif
