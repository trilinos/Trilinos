
#include <stk_adapt/UniformRefinerPattern.hpp>
#include <stk_mesh/base/EntityComm.hpp>

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
#define DEBUG_MULTI_LEVEL 0
#if DEBUG_MULTI_LEVEL
      int parent_to_family_tree_relations_size_0 = parent_to_family_tree_relations.size();
      if (parent_to_family_tree_relations.size() == 1)
        {
          std::cout << " parent_to_family_tree_relations.size()==1 && eMesh.isChildElement(parent_elem) = " <<  eMesh.isChildElement(parent_elem) 
                    << " parent_elem= " << parent_elem
                    << std::endl;
        }
#endif      
      // if this is the first time the parent_elem has been visited, or if the parent_elem is the child of another parent,
      //   (at level 0 only, which is what isChildElement checks), then we need to add a new family tree

      if (parent_to_family_tree_relations.size() == 0 || (parent_to_family_tree_relations.size() == 1 && eMesh.isChildElement(parent_elem) ) )
        {
          stk::mesh::PartVector add(1, &eMesh.getFEM_meta_data()->universal_part());

          // explanation: we want to avoid the above use of BulkData::generate_new_entities due to the parallel comm required, so we
          //   use the parent_id for the familty_tree_id.
          // there are two types of family tree uses, one for 
          //   the first level of parent/child (FAMILY_TREE_LEVEL_0) and one that holds a child that now is a parent (FAMILY_TREE_LEVEL_1)
          //
          // Since we know that the parent_id is unique across processors, we can use it for the family tree
          //   and guarantee uniqueness of family tree id's across processors.
          stk::mesh::EntityId parent_id = parent_elem.identifier();
          stk::mesh::EntityId family_tree_id = parent_id;

          // FIXME
          if (parent_elem.entity_rank() != eMesh.element_rank())
            {
              stk::mesh::EntityId FT_SHIFT_SIDE = 100000000000ull;
              if (family_tree_id > FT_SHIFT_SIDE) 
                throw std::logic_error("FT_SHIFT_SIDE error in set_parent_child_relations");
              family_tree_id += FT_SHIFT_SIDE;
              //std::cout << "tmp family_tree_id = " << family_tree_id << " parent_id= " << parent_id << std::endl;
            }


#if DEBUG_MULTI_LEVEL
          if (parent_to_family_tree_relations.size() == 1) 
            {
              std::cout << "tmp family_tree_id = " << family_tree_id << " parent_id= " << parent_id << std::endl;
            }
#endif
          //unsgined FT_SHIFT = 100000000u;
          unsigned FT_SHIFT = 0u;
          family_tree_id += FT_SHIFT;

          family_tree = & eMesh.getBulkData()->declare_entity(FAMILY_TREE_RANK, family_tree_id, add);

          // make the parent be the first relation; children are at the end
          // from->to
#if DEBUG_MULTI_LEVEL
          std::cout << "tmp super->parent " << family_tree->identifier() << " -> " << parent_elem.identifier() << " " << parent_elem << std::endl;
#endif
          eMesh.getBulkData()->declare_relation(*family_tree, parent_elem, FAMILY_TREE_PARENT);
          //eMesh.getBulkData()->declare_relation( parent_elem, *family_tree, ptft_size-1);
          parent_to_family_tree_relations = parent_elem.relations(FAMILY_TREE_RANK);

#if DEBUG_MULTI_LEVEL
          unsigned ptft_size = parent_to_family_tree_relations.size();
          VERIFY_OP_ON(family_tree, !=, 0,"err1");
          if (1)
            {
              std::cout 
                << "ptft_size= " << ptft_size << " family_tree= " << family_tree
                << " family_tree_id= " << family_tree->identifier() << std::endl;
              for (unsigned k = 0; k < ptft_size; k++)
                {
                  stk::mesh::Entity * ft = parent_to_family_tree_relations[k].entity();
                  std::cout << "ft= " << ft << " ft_id= " << ft->identifier() << std::endl;
                }
            }
          VERIFY_OP_ON(family_tree, ==, parent_to_family_tree_relations[FAMILY_TREE_LEVEL_0].entity(),"err2");
#endif
        }

#if DEBUG_MULTI_LEVEL
      int parent_to_family_tree_relations_size_1 = parent_to_family_tree_relations.size();
#endif      

      if (parent_to_family_tree_relations.size() == 1)
        {
          //VERIFY_OP_ON(family_tree, !=, 0,"err1");
          //VERIFY_OP_ON(family_tree, ==, parent_to_family_tree_relations[FAMILY_TREE_LEVEL_0].entity(),"err2");
          
          unsigned parent_elem_ft_level_0 = eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, parent_elem);
          //family_tree = parent_to_family_tree_relations[FAMILY_TREE_LEVEL_0].entity();
          family_tree = parent_to_family_tree_relations[parent_elem_ft_level_0].entity();
        }
      else if (parent_to_family_tree_relations.size() == 2)
        {
          //VERIFY_OP_ON(family_tree, !=, 0,"err1");
          //VERIFY_OP_ON(family_tree, ==, parent_to_family_tree_relations[FAMILY_TREE_LEVEL_1].entity(),"err2");
          //family_tree = parent_to_family_tree_relations[FAMILY_TREE_LEVEL_1].entity();

          // EXPLANATION:  stk_mesh inserts back-relations in front of existing relations (it uses the std::vector<Relation>::insert method)
          // FIXME - need a unit test to check if this ever breaks in the future (i.e. going to boost::mesh)
          //family_tree = parent_to_family_tree_relations[FAMILY_TREE_LEVEL_0].entity();

          unsigned parent_elem_ft_level_1 = eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_1, parent_elem);
          family_tree = parent_to_family_tree_relations[parent_elem_ft_level_1].entity();

        }
      else
        {
          throw std::logic_error("UniformRefinerPatternBase::set_parent_child_relations no family_tree");
        }

      //entity_vector.reserve(getNumNewElemPerElem());
      unsigned nchild = getNumNewElemPerElem();
      if (numChild) nchild = *numChild;

#if DEBUG_MULTI_LEVEL
      std::cout << "tmp parent_to_family_tree_relations_size_0 = " << parent_to_family_tree_relations_size_0 << " parent_to_family_tree_relations_size_1= " << parent_to_family_tree_relations_size_1 << std::endl;
      std::cout << "tmp super->child " << family_tree->identifier() << " -> " << newElement.identifier() << " [" << (ordinal+1) << "]" << newElement << std::endl;
#endif

      // error check
      if (1)
        {
          mesh::PairIterRelation family_tree_relations = family_tree->relations(parent_elem.entity_rank());
          for (unsigned i = 1; i < family_tree_relations.size(); i++)
            {
              if (family_tree_relations[i].identifier() == (ordinal + 1))
                {
                  std::cout << "UniformRefinerPatternBase::set_parent_child_relations trying to refine a parent element again, or error in ordinal" << std::endl;
                  throw std::logic_error("UniformRefinerPatternBase::set_parent_child_relations trying to refine a parent element again, or error in ordinal");
                }
            }
        }

      eMesh.getBulkData()->declare_relation(*family_tree, newElement, ordinal + 1);  // the + 1 here is to give space for the parent

      // add all the nodes for ghosting purposes
      /** Explanation: child elements can be created in the aura that have nodes in the aura but aren't shared
       *  which doesn't bring over parent/child relations to the other processors.  The code below adds relations
       *  to all nodes of the parent elements thus providing necessary links that the closure code can follow to
       *  gather parent/child relations and send to sharing procs.
       */
      bool workaround_shared_node_issue = true;
      if (workaround_shared_node_issue)
        {
          
          mesh::PairIterRelation parent_elem_nodes = parent_elem.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
          for (unsigned i = 0; i < parent_elem_nodes.size(); i++)
            {
              if (!stk::mesh::in_shared(*parent_elem_nodes[i].entity())) continue;
                
              bool found = false;
              mesh::PairIterRelation ft_nodes = family_tree->relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
              for (unsigned j = 0; j < ft_nodes.size(); j++)
                {
                  if (ft_nodes[j].entity() == parent_elem_nodes[i].entity())
                    {
                      found = true;
                      break;
                    }
                }
              if (!found)
                {
                  eMesh.getBulkData()->declare_relation(*family_tree, *parent_elem_nodes[i].entity(), ft_nodes.size());
                }
            }

          stk::mesh::PairIterRelation child_elem_nodes = newElement.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
          for (unsigned i = 0; i < child_elem_nodes.size(); i++)
            {
              if (!stk::mesh::in_shared(*child_elem_nodes[i].entity())) continue;

              bool found = false;
              mesh::PairIterRelation ft_nodes = family_tree->relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
              for (unsigned j = 0; j < ft_nodes.size(); j++)
                {
                  if (ft_nodes[j].entity() == child_elem_nodes[i].entity())
                    {
                      found = true;
                      break;
                    }
                }
              if (!found)
                {
                  eMesh.getBulkData()->declare_relation(*family_tree, *child_elem_nodes[i].entity(), ft_nodes.size());
                }
            }

          // check for second level and subsequent refinement
          if (parent_to_family_tree_relations.size() == 2)
            {
              unsigned parent_elem_ft_level_0 = eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, parent_elem);
              stk::mesh::Entity* family_tree_level_0 = parent_to_family_tree_relations[parent_elem_ft_level_0].entity();
              
              stk::mesh::PairIterRelation ft_level_0_nodes = family_tree_level_0->relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
              for (unsigned i = 0; i < ft_level_0_nodes.size(); i++)
                {
                  if (!stk::mesh::in_shared(*ft_level_0_nodes[i].entity())) continue;

                  bool found = false;
                  mesh::PairIterRelation ft_nodes = family_tree->relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
                  for (unsigned j = 0; j < ft_nodes.size(); j++)
                    {
                      if (ft_nodes[j].entity() == ft_level_0_nodes[i].entity())
                        {
                          found = true;
                          break;
                        }
                    }
                  if (!found)
                    {
                      eMesh.getBulkData()->declare_relation(*family_tree, *ft_level_0_nodes[i].entity(), ft_nodes.size());
                    }
                }
            }

        }


      if (0) std::cout << "tmp here 12 ordinal= " << ordinal << " [ " << getNumNewElemPerElem() << "] newElement_ptr= "<< &newElement<< std::endl;
#endif
    }
#endif

    void UniformRefinerPatternBase::interpolateElementFields(percept::PerceptMesh& eMesh, stk::mesh::Entity& old_owning_elem, stk::mesh::Entity& newElement)
    {
      // FIXME
//       if (old_owning_elem.entity_rank() != eMesh.element_rank())
//         {
//           return;
//         }
      const stk::mesh::FieldVector & fields = eMesh.getFEM_meta_data()->get_fields();
      unsigned nfields = fields.size();
      for (unsigned ifld = 0; ifld < nfields; ifld++)
        {
          stk::mesh::FieldBase *field = fields[ifld];
          int field_dimension = -1;

          stk::mesh::EntityRank field_rank = stk::mesh::fem::FEMMetaData::NODE_RANK;
          {
            unsigned nfr = field->restrictions().size();
            //if (Util::getFlag(1234)) std::cout << "tmp    number of field restrictions= " << nfr << " for field= " << field->name() <<  std::endl;
            for (unsigned ifr = 0; ifr < nfr; ifr++)
              {
                const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
                //mesh::Part& frpart = metaData.get_part(fr.ordinal());
                
                field_rank = fr.entity_rank();
                field_dimension = fr.dimension() ;
                //if (Util::getFlag(1234)) std::cout << "tmp field_rank= " << field_rank << " field_dimension= " << field_dimension << std::endl;
              }
          }
          if (field_rank == stk::mesh::fem::FEMMetaData::NODE_RANK)
            {
              continue;
            }
          if (field_rank == old_owning_elem.entity_rank())
            {
              unsigned stride_old=0, stride_new=0;
              double *fdata_old = PerceptMesh::field_data_entity(field, old_owning_elem, &stride_old);
              if (!fdata_old) 
                continue;
              if ((int)stride_old != field_dimension)
                {
                  VERIFY_OP_ON((int)stride_old, ==, field_dimension, "interpolateElementFields err1");
                  throw std::runtime_error("interpolateElementFields err1");
                }
              double *fdata_new = PerceptMesh::field_data_entity(field, newElement,  &stride_new);
              if (!fdata_new) 
                continue;
              if ((int)stride_new != field_dimension || stride_new != stride_old)
                {
                  VERIFY_OP_ON((int)stride_new, ==, field_dimension, "interpolateElementFields err2");
                  VERIFY_OP_ON(stride_new, ==, stride_old, "interpolateElementFields err3");
                  throw std::runtime_error("interpolateElementFields err2");
                }
              for (unsigned i = 0; i < stride_old; i++)
                {
                  fdata_new[i] = fdata_old[i];
                }
            }
      
        }
    }


  }
}

