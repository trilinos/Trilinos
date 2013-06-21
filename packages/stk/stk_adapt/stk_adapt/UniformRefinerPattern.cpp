
#include <stk_adapt/UniformRefinerPattern.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>

namespace stk {
  namespace adapt {

    bool allow_single_refine = false;

    STK_Adapt_Auto_Part stk_adapt_auto_part;

    const std::string UniformRefinerPatternBase::m_oldElementsPartName = "urp_oldElements";

    std::string UniformRefinerPatternBase::s_convert_options = "Tet4_Wedge6_Hex8, Wedge6_Hex8_6, Tet4_Hex8_4, Hex8_Tet4_24, Hex8_Tet4_6, Quad4_Tri3_2, Quad4_Tri3_6, Quad4_Tri3_4, Tri3_Quad4_3";
    std::string UniformRefinerPatternBase::s_refine_options = "DEFAULT, Quad4_Quad4_4, Tri3_Tri3_4, Tet4_Tet4_8, Hex8_Hex8_8, Wedge6_Wedge6_8, Pyramid5_Pyramid5_10, "
      " Tri6_Tri6_4, Quad9_Quad9_4, Hex27_Hex27_8, Tet10_Tet10_8, Wedge15_Wedge15_8, Pyramid13_Pyramid13_10, ShellTri3_ShellTri3_4, ShellQuad4_ShellQuad4_4";
    std::string UniformRefinerPatternBase::s_enrich_options = "DEFAULT, Quad4_Quad8_1, Quad4_Quad9_1, Tri3_Tri6_1, Tet4_Tet10_1, Hex8_Hex20_1, Hex8_Hex27_1, "
      " Wedge6_Wedge15_1, Wedge6_Wedge18_1, Pyramid5_Pyramid13_1, Beam2_Beam3_1";


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
      else if (refine == "Pyramid5_Pyramid5_10")  pattern  = Teuchos::rcp(new Pyramid5_Pyramid5_10(eMesh, block_names));

      //    shells
      else if (refine == "ShellTri3_ShellTri3_4")      pattern  = Teuchos::rcp(new ShellTri3_ShellTri3_4(eMesh, block_names));
      else if (refine == "ShellQuad4_ShellQuad4_4")      pattern  = Teuchos::rcp(new ShellQuad4_ShellQuad4_4(eMesh, block_names));

      else if (refine == "Tri6_Tri6_4")      pattern  = Teuchos::rcp(new Tri6_Tri6_4(eMesh, block_names));
      else if (refine == "Quad9_Quad9_4")    pattern  = Teuchos::rcp(new Quad9_Quad9_4(eMesh, block_names));
      else if (refine == "Hex27_Hex27_8")    pattern  = Teuchos::rcp(new Hex27_Hex27_8(eMesh, block_names));
      else if (refine == "Tet10_Tet10_8")    pattern  = Teuchos::rcp(new Tet10_Tet10_8(eMesh, block_names));
      else if (refine == "Wedge15_Wedge15_8") pattern = Teuchos::rcp(new Wedge15_Wedge15_8(eMesh, block_names));
      //else if (refine == "Wedge18_Wedge18_8") pattern = Teuchos::rcp(new Wedge18_Wedge18_8(eMesh, block_names));
      else if (refine == "Pyramid13_Pyramid13_10") pattern = Teuchos::rcp(new Pyramid13_Pyramid13_10(eMesh, block_names));

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
      else if (enrich == "Pyramid5_Pyramid13_1") pattern  = Teuchos::rcp(new Pyramid5_Pyramid13_1(eMesh, block_names));
      else if (enrich == "Beam2_Beam3_1")    pattern  = Teuchos::rcp(new Beam2_Beam3_1(eMesh, block_names));

      // convert
      //else if (convert == "DEFAULT")         pattern  = findDefaultConvert(eMesh, block_names);
      else if (convert == "Quad4_Tri3_2")    pattern  = Teuchos::rcp(new Quad4_Tri3_2(eMesh, block_names));
      else if (convert == "Quad4_Tri3_6")    pattern  = Teuchos::rcp(new Quad4_Tri3_6(eMesh, block_names));
      else if (convert == "Quad4_Tri3_4")    pattern  = Teuchos::rcp(new Quad4_Tri3_4(eMesh, block_names));
      else if (convert == "Tri3_Quad4_3")    pattern  = Teuchos::rcp(new Tri3_Quad4_3(eMesh, block_names));
      else if (convert == "Hex8_Tet4_24")    pattern  = Teuchos::rcp(new Hex8_Tet4_24(eMesh, block_names));
      else if (convert == "Hex8_Tet4_6")     pattern  = Teuchos::rcp(new Hex8_Tet4_6_12(eMesh, block_names));
      else if (convert == "Tet4_Hex8_4")     pattern  = Teuchos::rcp(new Tet4_Hex8_4(eMesh, block_names));
      else if (convert == "Wedge6_Hex8_6")   pattern  = Teuchos::rcp(new Wedge6_Hex8_6(eMesh, block_names));
      else if (convert == "Tet4_Wedge6_Hex8")    pattern  = Teuchos::rcp(new Tet4_Wedge6_Hex8(eMesh, block_names));
      else
        {
          throw std::invalid_argument( (std::string("UniformRefinerPatternBase::createPattern unknown string: refine= ")+refine+" enrich= "+enrich+
                                        " convert= " + convert).c_str() );
        }

      return pattern;
    }


    /*
      static const SameRankRelationValue * getChildVectorPtr(  SameRankRelation& repo , Entity parent)
      {
      SameRankRelation::const_iterator i = repo.find( parent );
      if (i != repo.end())
      return &i->second;
      else
      return 0;
      }
    */


    /// if numChild is passed in as non-null, use that value, else use getNumNewElemPerElem() as size of child vector
    void UniformRefinerPatternBase::set_parent_child_relations(percept::PerceptMesh& eMesh, stk::mesh::Entity parent_elem, stk::mesh::Entity newElement,
                                                               unsigned ordinal, unsigned *numChild)
    {
#if NEW_FIX_ELEMENT_SIDES

      //VERIFY_OP(ordinal, < , getNumNewElemPerElem(), "logic error in set_parent_child_relations");
      VERIFY_OP(parent_elem, != , stk::mesh::Entity(), "set_parent_child_relations: parent_elem is null");
      VERIFY_OP(newElement, != , stk::mesh::Entity(), "set_parent_child_relations: newElement is null");

      if (!eMesh.is_valid(parent_elem))
        {
          throw std::logic_error("UniformRefinerPatternBase::set_parent_child_relations parent_elem is null");
        }

      stk::mesh::FieldBase *refine_level = eMesh.get_field("refine_level");
      if (refine_level)
        {
          double *fdata_new = eMesh.field_data( *static_cast<const ScalarFieldType *>(refine_level) , newElement );
          double *fdata = eMesh.field_data( *static_cast<const ScalarFieldType *>(refine_level) , parent_elem );
          if (fdata && fdata_new)
            fdata_new[0] = fdata[0] + 1.0;
          //std::cout << "fdata= " << fdata << " fdata_new= " << fdata_new[0] << std::endl;
        }

      const unsigned FAMILY_TREE_RANK = stk::mesh::MetaData::ELEMENT_RANK + 1u;
      stk::mesh::Entity family_tree = stk::mesh::Entity();
      percept::MyPairIterRelation parent_to_family_tree_relations (eMesh, parent_elem, FAMILY_TREE_RANK);
      // if this is the first time the parent_elem has been visited, or if the parent_elem is the child of another parent,
      //   (at level 0 only, which is what isChildElement checks), then we need to add a new family tree

      if (parent_to_family_tree_relations.size() == 0 || (parent_to_family_tree_relations.size() == 1 && eMesh.isChildElement(parent_elem) ) )
        {
          stk::mesh::PartVector add(1, &eMesh.get_fem_meta_data()->universal_part());

          // explanation: we want to avoid the above use of BulkData::generate_new_entities due to the parallel comm required, so we
          //   use the parent_id for the familty_tree_id.
          // there are two types of family tree uses, one for
          //   the first level of parent/child (FAMILY_TREE_LEVEL_0) and one that holds a child that now is a parent (FAMILY_TREE_LEVEL_1)
          //
          // Since we know that the parent_id is unique across processors, we can use it for the family tree
          //   and guarantee uniqueness of family tree id's across processors.
          // 06/21/13: but, not with pooling active, as used in unrefinePass2... so, see below for getNextId usage

          stk::mesh::EntityId parent_id = eMesh.identifier(parent_elem);
          stk::mesh::EntityId family_tree_id = parent_id;

          family_tree_id = eMesh.getNextId(FAMILY_TREE_RANK);

          // FIXME
          if (0 && eMesh.entity_rank(parent_elem) != stk::mesh::MetaData::ELEMENT_RANK)
            {
              stk::mesh::EntityId FT_SHIFT_SIDE = 100000000000ull;
              if (family_tree_id > FT_SHIFT_SIDE)
                throw std::logic_error("FT_SHIFT_SIDE error in set_parent_child_relations");
              family_tree_id += FT_SHIFT_SIDE;
              //std::cout << "tmp family_tree_id = " << family_tree_id << " parent_id= " << parent_id << std::endl;
            }

          //stk::mesh::EntityId FT_SHIFT = 100000000000ull;
          stk::mesh::EntityId FT_SHIFT = 0ull;
          family_tree_id += FT_SHIFT;

          family_tree = eMesh.get_bulk_data()->declare_entity(FAMILY_TREE_RANK, family_tree_id, add);

          // make the parent be the first relation; children are at the end
          // from->to
          eMesh.get_bulk_data()->declare_relation(family_tree, parent_elem, FAMILY_TREE_PARENT);
          //eMesh.get_bulk_data()->declare_relation( parent_elem, *family_tree, ptft_size-1);
          percept::MyPairIterRelation new_ptf(eMesh,parent_elem,FAMILY_TREE_RANK);
          parent_to_family_tree_relations = new_ptf;

        }

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
      //
      //unsigned nchild = getNumNewElemPerElem();
      //if (numChild) nchild = *numChild;

      // error check
      if (1)
        {
          percept::MyPairIterRelation family_tree_relations (eMesh, family_tree, eMesh.entity_rank(parent_elem));
          for (unsigned i = 1; i < family_tree_relations.size(); i++)
            {
              if (family_tree_relations[i].relation_ordinal() == (ordinal + 1))
                {
                  std::cout << "UniformRefinerPatternBase::set_parent_child_relations trying to refine a parent element again, or error in ordinal ["
                            << ordinal << "]" << " family_tree_relations.size= " << family_tree_relations.size() 
                            << " parent_elem= " << eMesh.identifier(parent_elem)
                            << std::endl;
                  throw std::logic_error("UniformRefinerPatternBase::set_parent_child_relations trying to refine a parent element again, or error in ordinal");
                }
            }
        }

      eMesh.get_bulk_data()->declare_relation(family_tree, newElement, ordinal + 1);  // the + 1 here is to give space for the parent
      percept::MyPairIterRelation parent_to_family_tree_relations1 (eMesh, parent_elem, FAMILY_TREE_RANK);

      // add all the nodes for ghosting purposes
      /** Explanation: child elements can be created in the aura that have nodes in the aura but aren't shared
       *  which doesn't bring over parent/child relations to the other processors.  The code below adds relations
       *  to all nodes of the parent elements thus providing necessary links that the closure code can follow to
       *  gather parent/child relations and send to sharing procs.
       */
      bool workaround_shared_node_issue = true;
      if (workaround_shared_node_issue)
        {

          std::set<stk::mesh::Entity, stk::mesh::EntityLess> to_add(stk::mesh::EntityLess(*eMesh.get_bulk_data()));

          percept::MyPairIterRelation parent_elem_nodes (eMesh, parent_elem,  stk::mesh::MetaData::NODE_RANK );
          for (unsigned i = 0; i < parent_elem_nodes.size(); i++)
            {
              if (! eMesh.get_bulk_data()->in_shared(eMesh.key(parent_elem_nodes[i].entity()))) continue;

              bool found = false;
              percept::MyPairIterRelation ft_nodes (eMesh, family_tree,  stk::mesh::MetaData::NODE_RANK );
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
                  to_add.insert(parent_elem_nodes[i].entity());
                  VERIFY_OP_ON(eMesh.get_bulk_data()->in_index_range(parent_elem_nodes[i].entity()), ==, true, "parent_elem_nodes bad");
                  //eMesh.get_bulk_data()->declare_relation(family_tree, parent_elem_nodes[i].entity(), ft_nodes.size());
                }
            }

          percept::MyPairIterRelation child_elem_nodes (eMesh, newElement,  stk::mesh::MetaData::NODE_RANK );
          if (child_elem_nodes.size() == 0)
            {
              throw std::runtime_error("child_elem has no nodes");
            }
          for (unsigned i = 0; i < child_elem_nodes.size(); i++)
            {
              if (!eMesh.get_bulk_data()->in_shared(eMesh.key(child_elem_nodes[i].entity()))) continue;

              bool found = false;
              percept::MyPairIterRelation ft_nodes (eMesh, family_tree,  stk::mesh::MetaData::NODE_RANK );
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
                  to_add.insert(child_elem_nodes[i].entity());
                  VERIFY_OP_ON(eMesh.get_bulk_data()->in_index_range(child_elem_nodes[i].entity()), ==, true, "child_elem_nodes bad");
                  //eMesh.get_bulk_data()->declare_relation(family_tree, child_elem_nodes[i].entity(), ft_nodes.size());

                }
            }

          // check for second level and subsequent refinement
          if (parent_to_family_tree_relations1.size() == 2)
            {
              unsigned parent_elem_ft_level_0 = eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, parent_elem);
              stk::mesh::Entity family_tree_level_0 = parent_to_family_tree_relations1[parent_elem_ft_level_0].entity();

              percept::MyPairIterRelation ft_level_0_nodes (eMesh, family_tree_level_0,  stk::mesh::MetaData::NODE_RANK );
              for (unsigned i = 0; i < ft_level_0_nodes.size(); i++)
                {
                  if (!eMesh.get_bulk_data()->in_shared(eMesh.key(ft_level_0_nodes[i].entity()))) continue;

                  bool found = false;
                  percept::MyPairIterRelation ft_nodes (eMesh, family_tree,  stk::mesh::MetaData::NODE_RANK );
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
                      VERIFY_OP_ON(eMesh.get_bulk_data()->in_index_range(ft_level_0_nodes[i].entity()), ==, true, "ft_level_0_nodes bad 0");
                      //eMesh.get_bulk_data()->declare_relation(family_tree, ft_level_0_nodes[i].entity(), ft_nodes.size());
                      to_add.insert(ft_level_0_nodes[i].entity());
                    }
                }
            }

          // add nodes to family_tree
          {
            percept::MyPairIterRelation ft_nodes (eMesh, family_tree,  stk::mesh::MetaData::NODE_RANK );
            unsigned ftns=ft_nodes.size();

            std::vector<stk::mesh::Entity> to_add_vec(to_add.begin(), to_add.end());

            for (unsigned ita=0; ita < to_add_vec.size(); ita++)
              {
                eMesh.get_bulk_data()->declare_relation(family_tree, to_add_vec[ita], ftns+ita);
              }
          }

        }


      if (0) std::cout << "tmp here 12 ordinal= " << ordinal << " [ " << getNumNewElemPerElem() << "] newElement_ptr= "<< &newElement<< std::endl;
      bool foundSide = findSideRelations(eMesh, parent_elem, newElement);
      if (!foundSide && eMesh.entity_rank(parent_elem) < eMesh.element_rank()) {
        //throw std::runtime_error("UniformRefinerPatternBase:: set_parent_child_relations couldn't set child side to elem relations");
      }
#endif
    }

    void UniformRefinerPatternBase::interpolateElementFields(percept::PerceptMesh& eMesh, std::vector<stk::mesh::Entity>& old_owning_elements, stk::mesh::Entity newElement)
    {
      // FIXME
//       if (m_eMesh.entity_rank(old_owning_elem) != stk::mesh::MetaData::ELEMENT_RANK)
//         {
//           return;
//         }
      const stk::mesh::FieldVector & fields = eMesh.get_fem_meta_data()->get_fields();
      double old_owning_elements_size = old_owning_elements.size();
      unsigned nfields = fields.size();
      for (unsigned ifld = 0; ifld < nfields; ifld++)
        {
          stk::mesh::FieldBase *field = fields[ifld];
          if (field->name()=="refine_level") continue;
          int field_dimension = -1;

          stk::mesh::EntityRank field_rank = stk::mesh::MetaData::NODE_RANK;
          {
            unsigned nfr = field->restrictions().size();
            //if (Util::getFlag(1234)) std::cout << "tmp    number of field restrictions= " << nfr << " for field= " << field->name() <<  std::endl;
            for (unsigned ifr = 0; ifr < nfr; ifr++)
              {
                const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
                //mesh::Part& frpart = metaData.get_part(fr.ordinal());

                field_rank = fr . entity_rank();
                field_dimension = fr.dimension() ;
                //if (Util::getFlag(1234)) std::cout << "tmp field_rank= " << field_rank << " field_dimension= " << field_dimension << std::endl;
              }
          }
          if (field_rank == stk::mesh::MetaData::NODE_RANK)
            {
              continue;
            }
          for (unsigned iel=0; iel < old_owning_elements.size(); iel++)
            {
              stk::mesh::Entity old_owning_elem = old_owning_elements[iel];

              if (field_rank == eMesh.entity_rank(old_owning_elem))
                {
                  unsigned stride_old=0, stride_new=0;
                  double *fdata_old = eMesh.field_data(field, old_owning_elem, &stride_old);
                  if (!fdata_old)
                    continue;
                  if ((int)stride_old != field_dimension)
                    {
                      VERIFY_OP_ON((int)stride_old, ==, field_dimension, "interpolateElementFields err1");
                      throw std::runtime_error("interpolateElementFields err1");
                    }
                  double *fdata_new = eMesh.field_data(field, newElement,  &stride_new);
                  if (!fdata_new)
                    continue;
                  if ((int)stride_new != field_dimension || stride_new != stride_old)
                    {
                      VERIFY_OP_ON((int)stride_new, ==, field_dimension, "interpolateElementFields err2");
                      VERIFY_OP_ON(stride_new, ==, stride_old, "interpolateElementFields err3");
                      throw std::runtime_error("interpolateElementFields err2");
                    }
                  if (iel == 0)
                    {
                      for (unsigned i = 0; i < stride_old; i++)
                        {
                          fdata_new[i] = 0.0;
                        }
                    }
                  for (unsigned i = 0; i < stride_old; i++)
                    {
                      fdata_new[i] += fdata_old[i] / old_owning_elements_size;
                    }
                }
            }
        }
    }

    bool UniformRefinerPatternBase::findSideRelations(percept::PerceptMesh& eMesh, stk::mesh::Entity parent, stk::mesh::Entity child)
    {
      VERIFY_OP_ON(eMesh.entity_rank(parent), ==, eMesh.entity_rank(child), "UniformRefinerPatternBase::findSideRelations: bad ranks");
      if (eMesh.entity_rank(parent) == stk::mesh::MetaData::ELEMENT_RANK)
        return true;

      for (unsigned higher_order_rank = eMesh.entity_rank(parent)+1u; higher_order_rank <= stk::mesh::MetaData::ELEMENT_RANK; higher_order_rank++)
        {
          percept::MyPairIterRelation parent_to_elem_rels (eMesh, parent, higher_order_rank);
          VERIFY_OP_ON(parent_to_elem_rels.size(), <=, 1, "UniformRefinerPatternBase::findSideRelations bad number of side to elem relations");
          if (parent_to_elem_rels.size() == 0)
            {
              // nothing to do
              return true;
            }

          for (unsigned i_parent_to_elem=0; i_parent_to_elem < parent_to_elem_rels.size(); i_parent_to_elem++)
            {
              stk::mesh::Entity parents_volume_element = parent_to_elem_rels[i_parent_to_elem].entity();

              std::vector<stk::mesh::Entity> parents_volume_elements_children;
              VERIFY_OP_ON(eMesh.hasFamilyTree(parents_volume_element), == , true, "UniformRefinerPatternBase::findSideRelations parent's volume element has no children.");
              //if (! eMesh.hasFamilyTree(*parents_volume_element) ) return true;
              eMesh.getChildren(parents_volume_element, parents_volume_elements_children);
              for (unsigned i_vol_child=0; i_vol_child < parents_volume_elements_children.size(); i_vol_child++)
                {
                  stk::mesh::Entity parents_volume_elements_child = parents_volume_elements_children[i_vol_child];

                  VERIFY_OP_ON(eMesh.entity_rank(parents_volume_elements_child), ==, higher_order_rank, "UniformRefinerPatternBase::findSideRelations: bad ranks 2");
                  if (connectSides(eMesh, parents_volume_elements_child, child))
                    return true;
                }
            }
        }
      return false;
    }

#define DEBUG_GSPR_1 1

    // if the element (element) has a side that matches  the given side (side_elem), connect them but first delete old connections
    bool UniformRefinerPatternBase::connectSides(percept::PerceptMesh& eMesh, stk::mesh::Entity element, stk::mesh::Entity side_elem)
    {
      EXCEPTWATCH;
      bool debug = false;
      //if (eMesh.identifier(side_elem) == 4348) debug = true;
      //       if (eMesh.identifier(element) == 71896)
      //         {
      //           if (eMesh.identifier(side_elem) == 11174) debug = true;
      //           if (eMesh.identifier(side_elem) == 10190) debug = true;
      //         }
      //       if (eMesh.identifier(side_elem) == 5 && eMesh.identifier(element) == 473) {
      //         debug = true;
      //       }

      shards::CellTopology element_topo(eMesh.get_cell_topology(element));
      unsigned element_nsides = (unsigned)element_topo.getSideCount();

      if (debug) {
        std::cout << "tmp srk connectSides element= "; eMesh.print(element, true, true);
        std::cout << " side= "; eMesh.print(side_elem, true, true);
      }

      // special case for shells
      int topoDim = UniformRefinerPatternBase::getTopoDim(element_topo);

      bool isShell = false;
      if (topoDim < (int)eMesh.entity_rank(element))
        {
          isShell = true;
        }
      int spatialDim = eMesh.get_spatial_dim();
      if (spatialDim == 3 && isShell && eMesh.entity_rank(side_elem) == eMesh.edge_rank())
        {
          element_nsides = (unsigned) element_topo.getEdgeCount();
        }

      if (debug) std::cout << "isShell= " << isShell << std::endl;

      int permIndex = -1;
      int permPolarity = 1;

      unsigned k_element_side = 0;

      // try search
      for (unsigned j_element_side = 0; j_element_side < element_nsides; j_element_side++)
        {
          bool use_coordinate_compare = false;
          eMesh.element_side_permutation(element, side_elem, j_element_side, permIndex, permPolarity, use_coordinate_compare, false);
          if (permIndex >= 0)
            {
              k_element_side = j_element_side;
              break;
            }
        }

      if (permIndex >= 0)
        {
          percept::MyPairIterRelation rels (eMesh, side_elem, eMesh.element_rank());

          // special case for shells
          if (isShell)
            {
              // FIXME for 2D
              if (eMesh.entity_rank(side_elem) == eMesh.face_rank())
                {
                  percept::MyPairIterRelation elem_sides (eMesh, element, eMesh.entity_rank(side_elem));
                  unsigned elem_sides_size= elem_sides.size();
                  if (debug) std::cout << "tmp srk found shell, elem_sides_size= " << elem_sides_size << std::endl;
                  if (elem_sides_size == 1)
                    {
                      stk::mesh::RelationIdentifier rel_id = elem_sides[0].relation_ordinal();
                      if (rel_id > 1)
                        throw std::logic_error("connectSides:: logic 1");
                      k_element_side = (rel_id == 0 ? 1 : 0);
                      if (debug) std::cout << "tmp srk k_element_side= " << k_element_side << " rel_id= " << rel_id << std::endl;
                    }
                }
            }

          int exists=0;
          percept::MyPairIterRelation elem_sides (eMesh, element, eMesh.entity_rank(side_elem));
          unsigned elem_sides_size= elem_sides.size();
          unsigned rel_id = 0;
          for (unsigned iside=0; iside < elem_sides_size; iside++)
            {
              stk::mesh::Entity existing_side = elem_sides[iside].entity();
              if (existing_side == side_elem)
                {
                  ++exists;
                  rel_id = elem_sides[iside].relation_ordinal();
                }

              if (elem_sides[iside].relation_ordinal() == k_element_side ) {
                std::cout << "ERROR: Relation already exists: connectSides element= "; eMesh.print(element, true, true);
                std::cout << " side= " << eMesh.identifier(side_elem); eMesh.print(side_elem, true, true);
                std::cout << " existing_side= " << eMesh.identifier(existing_side); eMesh.print(existing_side, true, true);

                if (DEBUG_GSPR_1)
                  {
                    std::cout << "connectSides:: ERROR: side = " << side_elem << std::endl;
                    eMesh.print(side_elem);
                    std::cout << "\nconnectSides:: ERROR: existing_side = " << existing_side << std::endl;
                    eMesh.print(existing_side);
                    std::cout << "\nelem= " << std::endl;
                    eMesh.print(element);

                    stk::mesh::PartVector side_parts;
                    eMesh.bucket(side_elem).supersets(side_parts);
                    for (unsigned isp = 0; isp < side_parts.size(); isp++)
                      {
                        std::cout << "side parts= " << side_parts[isp]->name() << std::endl;
                      }

                    stk::mesh::PartVector existing_side_parts;
                    eMesh.bucket(existing_side).supersets(existing_side_parts);
                    for (unsigned isp = 0; isp < side_parts.size(); isp++)
                      {
                        std::cout << "existing_side parts= " << existing_side_parts[isp]->name() << std::endl;
                      }

                    stk::mesh::PartVector elem_parts;
                    eMesh.bucket(element).supersets(elem_parts);
                    for (unsigned isp = 0; isp < elem_parts.size(); isp++)
                      {
                        std::cout << "elem parts= " << elem_parts[isp]->name() << std::endl;
                      }

                    VERIFY_OP_ON(elem_sides[iside].relation_ordinal(), !=, k_element_side, "Relation already exists!");

                  }

                VERIFY_OP_ON(elem_sides[iside].relation_ordinal(), !=, k_element_side, "Relation already exists!");
              }

            }
          if (!exists)
            {
              EXCEPTWATCH;
              eMesh.get_bulk_data()->declare_relation(element, side_elem, k_element_side);
            }
          else
            {
              VERIFY_OP_ON(k_element_side, ==, rel_id, "hmmm");
            }
          return true;
        }
      else
        {
          return false;
        }

    }


  }
}

