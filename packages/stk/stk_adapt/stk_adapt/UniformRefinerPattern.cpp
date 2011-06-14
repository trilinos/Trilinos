
#include <stk_adapt/UniformRefinerPattern.hpp>

namespace stk {
  namespace adapt {

    std::string UniformRefinerPatternBase::s_convert_options = "Hex8_Tet4_24, Hex8_Tet4_6, Quad4_Tri3_2, Quad4_Tri3_6, Quad4_Tri3_4";
    std::string UniformRefinerPatternBase::s_refine_options = "DEFAULT, Quad4_Quad4_4, Tri3_Tri3_4, Tet4_Tet4_8, Hex8_Hex8_8, Wedge6_Wedge6_8, Tri6_Tri6_4, Quad9_Quad9_4, Hex27_Hex27_8, Tet10_Tet10_8, Wedge18_Wedge18_8, ShellTri3_ShellTri3_4, ShellQuad4_ShellQuad4_4";
    std::string UniformRefinerPatternBase::s_enrich_options = "DEFAULT, Quad4_Quad8_1, Quad4_Quad9_1, Tri3_Tri6_1, Tet4_Tet10_1, Hex8_Hex20_1, Hex8_Hex27_1, Wedge6_Wedge15_1, Wedge6_Wedge18_1";
    

#if 0
    Teuchos::RCP<UniformRefinerPatternBase> UniformRefinerPatternBase::
    findDefaultConvert(percept::PerceptMesh& eMesh, BlockNamesType& block_names)
    {
      Teuchos::RCP<UniformRefinerPatternBase> pattern;

      
      //else if (convert == "Quad4_Tri3_6")    pattern  = Teuchos::rcp(new Quad4_Tri3_6(eMesh, block_names));
      pattern  = Teuchos::rcp(new Quad4_Tri3_4(eMesh, block_names));
      else if (convert == "Hex8_Tet4_24")    pattern  = Teuchos::rcp(new Hex8_Tet4_24(eMesh, block_names));
      else if (convert == "Hex8_Tet4_6")     pattern  = Teuchos::rcp(new Hex8_Tet4_6_12(eMesh, block_names));
      
    }
#endif

    Teuchos::RCP<UniformRefinerPatternBase> UniformRefinerPatternBase::
    createPattern(std::string refine, std::string enrich, std::string convert, percept::PerceptMesh& eMesh, BlockNamesType& block_names)
    {
      Teuchos::RCP<UniformRefinerPatternBase> pattern;
      
      // refine
      if      (refine == "DEFAULT")          pattern  = Teuchos::rcp(new URP_Heterogeneous_3D(eMesh, block_names));
      else if (refine == "Quad4_Quad4_4")    pattern  = Teuchos::rcp(new Quad4_Quad4_4_Sierra(eMesh, block_names));
      else if (refine == "Tri3_Tri3_4")      pattern  = Teuchos::rcp(new Tri3_Tri3_4(eMesh, block_names));
      else if (refine == "Tet4_Tet4_8")      pattern  = Teuchos::rcp(new Tet4_Tet4_8(eMesh, block_names));
      else if (refine == "Hex8_Hex8_8")      pattern  = Teuchos::rcp(new Hex8_Hex8_8(eMesh, block_names));
      else if (refine == "Wedge6_Wedge6_8")  pattern  = Teuchos::rcp(new Wedge6_Wedge6_8(eMesh, block_names));

      //    shells
      else if (refine == "ShellTri3_ShellTri3_4")      pattern  = Teuchos::rcp(new ShellTri3_ShellTri3_4(eMesh, block_names));
      else if (refine == "ShellQuad4_ShellQuad4_4")      pattern  = Teuchos::rcp(new ShellQuad4_ShellQuad4_4(eMesh, block_names));

      else if (refine == "Tri6_Tri6_4")      pattern  = Teuchos::rcp(new Tri6_Tri6_4(eMesh, block_names));
      else if (refine == "Quad9_Quad9_4")    pattern  = Teuchos::rcp(new Quad9_Quad9_4(eMesh, block_names));
      else if (refine == "Hex27_Hex27_8")    pattern  = Teuchos::rcp(new Hex27_Hex27_8(eMesh, block_names));
      else if (refine == "Tet10_Tet10_8")    pattern  = Teuchos::rcp(new Tet10_Tet10_8(eMesh, block_names));
      else if (refine == "Wedge15_Wedge15_8") pattern = Teuchos::rcp(new Wedge15_Wedge15_8(eMesh, block_names));
      //else if (refine == "Wedge18_Wedge18_8") pattern = Teuchos::rcp(new Wedge18_Wedge18_8(eMesh, block_names));

      // enrich
      else if (enrich == "DEFAULT")          pattern  = Teuchos::rcp(new URP_Heterogeneous_Enrich_3D(eMesh, block_names));
      else if (enrich == "Quad4_Quad8_1")    pattern  = Teuchos::rcp(new Quad4_Quad8_1(eMesh, block_names));
      else if (enrich == "Quad4_Quad9_1")    pattern  = Teuchos::rcp(new Quad4_Quad9_1(eMesh, block_names));
      else if (enrich == "Tri3_Tri6_1")      pattern  = Teuchos::rcp(new Tri3_Tri6_1(eMesh, block_names));
      else if (enrich == "Tet4_Tet10_1")     pattern  = Teuchos::rcp(new Tet4_Tet10_1(eMesh, block_names));
      else if (enrich == "Hex8_Hex20_1")     pattern  = Teuchos::rcp(new Hex8_Hex20_1(eMesh, block_names));
      else if (enrich == "Hex8_Hex27_1")     pattern  = Teuchos::rcp(new Hex8_Hex27_1(eMesh, block_names));
      else if (enrich == "Wedge6_Wedge15_1") pattern  = Teuchos::rcp(new Wedge6_Wedge15_1(eMesh, block_names));
      else if (enrich == "Wedge6_Wedge18_1") pattern  = Teuchos::rcp(new Wedge6_Wedge18_1(eMesh, block_names));

      // convert
      //else if (convert == "DEFAULT")         pattern  = findDefaultConvert(eMesh, block_names);
      else if (convert == "Quad4_Tri3_6")    pattern  = Teuchos::rcp(new Quad4_Tri3_6(eMesh, block_names));
      else if (convert == "Quad4_Tri3_4")    pattern  = Teuchos::rcp(new Quad4_Tri3_4(eMesh, block_names));
      else if (convert == "Hex8_Tet4_24")    pattern  = Teuchos::rcp(new Hex8_Tet4_24(eMesh, block_names));
      else if (convert == "Hex8_Tet4_6")     pattern  = Teuchos::rcp(new Hex8_Tet4_6_12(eMesh, block_names));
      else
        {
          throw std::invalid_argument( (std::string("UniformRefinerPatternBase::createPattern unknown string: refine= ")+refine+" enrich= "+enrich+
                                        " convert= " + convert).c_str() );
        }

      return pattern;
    }


    /*
      static const SameRankRelationValue * getChildVectorPtr(  SameRankRelation& repo , Entity *parent)
      {
      SameRankRelation::const_iterator i = repo.find( parent );
      if (i != repo.end()) 
      return &i->second;
      else
      return 0;
      }
    */

#if PERCEPT_USE_FAMILY_TREE == 0
    /// if numChild is passed in as non-null, use that value, else use getNumNewElemPerElem() as size of child vector
    void UniformRefinerPatternBase::set_parent_child_relations(percept::PerceptMesh& eMesh, stk::mesh::Entity& parent_elem, stk::mesh::Entity& newElement, 
                                                               unsigned ordinal, unsigned *numChild)
    {
#if NEW_FIX_ELEMENT_SIDES
      VERIFY_OP(ordinal, < , getNumNewElemPerElem(), "logic error in set_parent_child_relations");
      VERIFY_OP(&parent_elem, != , 0, "set_parent_child_relations: parent_elem is null");
      VERIFY_OP(&newElement, != , 0, "set_parent_child_relations: newElement is null");

      if (0 == &parent_elem)
        {
          throw std::logic_error("UniformRefinerPatternBase::set_parent_child_relations parent_elem is null");
        }

      PerceptEntityVector& entity_vector = eMesh.adapt_parent_to_child_relations()[&parent_elem];

      //entity_vector.reserve(getNumNewElemPerElem());
#if 0
      unsigned nchild = getNumNewElemPerElem();
      if (numChild) nchild = *numChild;
      if (entity_vector.size() != nchild)
        {
          entity_vector.resize(nchild);
        }
#else
      if (ordinal + 1 > entity_vector.size())
        {
          entity_vector.resize(ordinal+1);
        }
#endif
      entity_vector[ordinal] = &newElement;

      if (0) std::cout << "tmp here 12 ordinal= " << ordinal << " [ " << getNumNewElemPerElem() << "] newElement_ptr= "<< &newElement<< std::endl;
#endif
    }

#elif PERCEPT_USE_FAMILY_TREE == 1
    /// if numChild is passed in as non-null, use that value, else use getNumNewElemPerElem() as size of child vector
    void UniformRefinerPatternBase::set_parent_child_relations(percept::PerceptMesh& eMesh, stk::mesh::Entity& parent_elem, stk::mesh::Entity& newElement, 
                                                               unsigned ordinal, unsigned *numChild)
    {
#if NEW_FIX_ELEMENT_SIDES

      VERIFY_OP(ordinal, < , getNumNewElemPerElem(), "logic error in set_parent_child_relations");
      VERIFY_OP(&parent_elem, != , 0, "set_parent_child_relations: parent_elem is null");
      VERIFY_OP(&newElement, != , 0, "set_parent_child_relations: newElement is null");

      if (0 == &parent_elem)
        {
          throw std::logic_error("UniformRefinerPatternBase::set_parent_child_relations parent_elem is null");
        }

      const unsigned FAMILY_TREE_RANK = eMesh.element_rank() + 1u;
      stk::mesh::Entity* family_tree = 0;
      mesh::PairIterRelation parent_to_family_tree_relations = parent_elem.relations(FAMILY_TREE_RANK);
      // if this is the first time the parent_elem has been visited, or if the parent_elem is the child of another parent,
      //   (at level 0 only, which is what isChildElement checks), then we need to add a new family tree

      bool ALLOW_MULTI_LEVEL = true;

      if ( (!ALLOW_MULTI_LEVEL && (parent_to_family_tree_relations.size() == 0 ) ) ||
           ( ALLOW_MULTI_LEVEL && (parent_to_family_tree_relations.size() == 0 || (parent_to_family_tree_relations.size() == 1 && eMesh.isChildElement(parent_elem))) )
           )
        {
          stk::mesh::PartVector add(1, &eMesh.getFEM_meta_data()->universal_part());

          if (0 && ALLOW_MULTI_LEVEL)
            {
              std::vector<stk::mesh::Entity *> new_ft;
              eMesh.createEntities( FAMILY_TREE_RANK, 1, new_ft);
              family_tree = new_ft[0];
            }
          else
            {
              // explanation: we want to avoid the above use of BulkData::generate_new_entities due to the parallel comm required, so we
              //   use a bit at the high end of EntityId to differentiate between the two types of family trees, the one that holds
              //   the first level of parent/child (FT_LEVEL_0) and one that holds a child that now is a parent (FT_LEVEL_1)
              // Since we know that the parent_id is unique across processors, we can use it for the family tree (and add the bit for the
              //   second type of family tree) and guarantee uniqueness of family tree id's across processors.
              stk::mesh::EntityId parent_id = parent_elem.identifier();
              stk::mesh::EntityId family_tree_id = parent_id;
              if (parent_to_family_tree_relations.size() == 1) 
                {
                  if (sizeof(stk::mesh::EntityId) < 8u)
                    {
                      throw std::logic_error("UniformRefinerPatternBase::set_parent_child_relations: EntityId not large enough on this platform.");
                    }
                  //stk::mesh::EntityId bit = 1;  // just to avoid having to use a hard-coded constant which may be platform dependent
                  //family_tree_id = family_tree_id | (bit << 60);  // actual
                  // tmp for readability
                  family_tree_id = family_tree_id + 1000000000u;
                  
                  //std::cout << "tmp family_tree_id = " << family_tree_id << " parent_id= " << parent_id << std::endl;
                }
              family_tree = & eMesh.getBulkData()->declare_entity(FAMILY_TREE_RANK, parent_id, add);
            }

          // make the parent be the first relation; children are at the end
          // from->to
          eMesh.getBulkData()->declare_relation(*family_tree, parent_elem, 0u);
          parent_to_family_tree_relations = parent_elem.relations(FAMILY_TREE_RANK);
        }

      if (parent_to_family_tree_relations.size() == 1)
        {
          family_tree = parent_to_family_tree_relations[0].entity();
        }
      else if (ALLOW_MULTI_LEVEL && parent_to_family_tree_relations.size() == 2)
        {
          family_tree = parent_to_family_tree_relations[1].entity();
        }
      else
        {
          throw std::logic_error("UniformRefinerPatternBase::set_parent_child_relations no family_tree");
        }

      //entity_vector.reserve(getNumNewElemPerElem());
      unsigned nchild = getNumNewElemPerElem();
      if (numChild) nchild = *numChild;

      // error check
      if (1)
        {
          mesh::PairIterRelation family_tree_relations = family_tree->relations(parent_elem.entity_rank());
          for (unsigned i = 1; i < family_tree_relations.size(); i++)
            {
              if (family_tree_relations[i].identifier() == (ordinal + 1))
                {
                  throw std::logic_error("UniformRefinerPatternBase::set_parent_child_relations trying to refine a parent element again, or error in ordinal");
                }
            }
        }

      eMesh.getBulkData()->declare_relation(*family_tree, newElement, ordinal + 1);  // the + 1 here is to give space for the parent

      // add all the nodes for ghosting purposes
      if (1)
        {
          
          mesh::PairIterRelation elem_nodes = parent_elem.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
          for (unsigned i = 0; i < elem_nodes.size(); i++)
            {
              bool found = false;
              mesh::PairIterRelation ft_nodes = family_tree->relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
              for (unsigned j = 0; j < ft_nodes.size(); j++)
                {
                  if (ft_nodes[j].entity() == elem_nodes[i].entity())
                    {
                      found = true;
                      break;
                    }
                }
              if (!found)
                {
                  eMesh.getBulkData()->declare_relation(*family_tree, *elem_nodes[i].entity(), ft_nodes.size());
                }
            }

          elem_nodes = newElement.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
          for (unsigned i = 0; i < elem_nodes.size(); i++)
            {
              bool found = false;
              mesh::PairIterRelation ft_nodes = family_tree->relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
              for (unsigned j = 0; j < ft_nodes.size(); j++)
                {
                  if (ft_nodes[j].entity() == elem_nodes[i].entity())
                    {
                      found = true;
                      break;
                    }
                }
              if (!found)
                {
                  eMesh.getBulkData()->declare_relation(*family_tree, *elem_nodes[i].entity(), ft_nodes.size());
                }
            }
        }

      if (0) std::cout << "tmp here 12 ordinal= " << ordinal << " [ " << getNumNewElemPerElem() << "] newElement_ptr= "<< &newElement<< std::endl;
#endif
    }
#endif

  }
}

