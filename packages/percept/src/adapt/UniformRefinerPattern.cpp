// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <percept/mesh/geometry/volume/VolumeUtil.hpp>
#include <percept/PerceptUtils.hpp>

#include <adapt/UniformRefinerPattern.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>

#if defined(STK_PERCEPT_USE_INTREPID)
#if !STK_PERCEPT_LITE
#include <percept/element/intrepid/BasisTable.hpp>
#endif
#endif

  namespace percept {

    bool allow_single_refine = false;
    bool s_allow_special_wedge_refine = false;  // for wedge boundary-layer refine

    const std::string UniformRefinerPatternBase::m_oldElementsPartName = "urp_oldElements";
    const std::string UniformRefinerPatternBase::m_appendConvertString = "_urpconv";

    std::string UniformRefinerPatternBase::s_convert_options =
      "Quad4_Tri3_2, Quad4_Tri3_6, Quad4_Tri3_4, Tri3_Quad4_3, Tet4_Wedge6_Hex8, Hex8_Wedge6_Pyramid5_Tet4, Wedge6_Hex8_6, Tet4_Hex8_4, Hex8_Tet4_24, Hex8_Tet4_6, Wedge6_Pyramid5_Tet4_Exclude_Unconnected_Wedges";
    std::string UniformRefinerPatternBase::s_refine_options =
      "DEFAULT, Quad4_Quad4_4, Tri3_Tri3_4, Tet4_Tet4_8, Hex8_Hex8_8, Wedge6_Wedge6_8, Pyramid5_Pyramid5_10, Tri6_Tri6_4, Quad9_Quad9_4, "
      "Hex27_Hex27_8, Tet10_Tet10_8, Wedge15_Wedge15_8, Pyramid13_Pyramid13_10, ShellTri3_ShellTri3_4, ShellQuad4_ShellQuad4_4";
    std::string UniformRefinerPatternBase::s_enrich_options =
      "DEFAULT, Quad4_Quad8_1, Quad4_Quad9_1, Tri3_Tri6_1, Tet4_Tet10_1, Hex8_Hex20_1, Hex8_Hex27_1, "
      "Wedge6_Wedge15_1, Wedge6_Wedge18_1, Pyramid5_Pyramid13_1, Beam2_Beam3_1";


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
      else if (convert == "Tet4_Wedge6_Hex8")pattern  = Teuchos::rcp(new Tet4_Wedge6_Hex8(eMesh, block_names));
      else if (convert == "Hex8_Wedge6_Pyramid5_Tet4")pattern  = Teuchos::rcp(new Hex8_Wedge6_Pyramid5_Tet4(eMesh, block_names));
      else
        {
          throw std::invalid_argument( (std::string("UniformRefinerPatternBase::createPattern unknown string: refine= ")+refine+" enrich= "+enrich+
                                        " convert= " + convert).c_str() );
        }

      return pattern;
    }

    // static
    /// if numChild is passed in as non-null, use that value, else use getNumNewElemPerElem() as size of child vector
    void UniformRefinerPatternBase::set_parent_child_relations(percept::PerceptMesh& eMesh, stk::mesh::Entity parent_elem, stk::mesh::Entity newElement,
                                                               stk::mesh::Entity familyTreeNewElement,
                                                               unsigned ordinal, unsigned *numChild)
    {
      bool debug = false;
      VERIFY_OP(parent_elem, != , stk::mesh::Entity(), "set_parent_child_relations: parent_elem is null");
      VERIFY_OP(newElement, != , stk::mesh::Entity(), "set_parent_child_relations: newElement is null");

      if (eMesh.m_parent_element_field)
        {
          ParentElementType::value_type *fdata_new = NULL;

          if (is_matching_rank(*eMesh.m_parent_element_field, newElement))
            {
              fdata_new = stk::mesh::field_data( *eMesh.m_parent_element_field , newElement );
              if (fdata_new)
                fdata_new[0] = static_cast<ParentElementType::value_type>(eMesh.identifier(parent_elem));
              if (debug && fdata_new)
                {
                  std::cout << "1URP fdata= for entity= " << eMesh.identifier(newElement) << " fdata_new = " << fdata_new[0] << " parent= " << eMesh.identifier(parent_elem) << std::endl;
                }
              if (fdata_new)
                VERIFY_OP_ON(fdata_new[0], ==, static_cast<ParentElementType::value_type>(eMesh.identifier(parent_elem)), "bad parent_field");

            }
          else if (eMesh.m_parent_element_field_side && is_matching_rank(*eMesh.m_parent_element_field_side, newElement))
            {
              VERIFY_OP_ON(eMesh.entity_rank(newElement), ==, eMesh.side_rank(), "hmmm");
              fdata_new = stk::mesh::field_data( *eMesh.m_parent_element_field_side , newElement );
              VERIFY_OP_ON(fdata_new, !=, 0, "bad fdata_new");
              stk::mesh::EntityId predicted_parent_id = 0;
              percept::MyPairIterRelation parent_to_element_relations (eMesh, parent_elem, eMesh.element_rank());
              size_t spt = parent_to_element_relations.size();
              if (spt >= 1)
                {
                  int which_relation = 0; // just pick the first one
                  const stk::mesh::ConnectivityOrdinal parent_ord_conn = parent_to_element_relations[which_relation].relation_ordinal();
                  predicted_parent_id = eMesh.exodus_side_id(eMesh.identifier(parent_to_element_relations[which_relation].entity()), parent_ord_conn);

                  if (fdata_new)
                    fdata_new[0] = static_cast<ParentElementType::value_type>(predicted_parent_id);

                  if (debug && fdata_new)
                    {
                      std::cout << "0URP fdata= for entity= " << eMesh.identifier(newElement) << " fdata_new = " << fdata_new[0] << " parent= " << eMesh.identifier(parent_elem) << std::endl;
                    }
                }
              else
                {
                  //VERIFY_OP_ON(spt, >=, 1, "bad parent_to_element_relations");
                }
            }
          else
            {
              throw std::runtime_error("set_parent_child_relations bad fields");
            }
        }

      if (!eMesh.is_valid(parent_elem))
        {
          throw std::logic_error("UniformRefinerPatternBase::set_parent_child_relations parent_elem is null");
        }

      if (!eMesh.m_refine_level_field_set)
        {
          eMesh.m_refine_level_field_set = true;
          eMesh.m_refine_level_field = eMesh.get_fem_meta_data()->get_field<RefineLevelType::value_type>(stk::topology::ELEMENT_RANK, "refine_level");
        }

      if (eMesh.m_refine_level_field)
        {
          RefineLevelType::value_type *fdata_new = NULL;
          RefineLevelType::value_type *fdata     = NULL;

          if(is_matching_rank(*eMesh.m_refine_level_field, newElement)) {
            fdata_new = stk::mesh::field_data( *eMesh.m_refine_level_field , newElement );
          }
          if(is_matching_rank(*eMesh.m_refine_level_field, parent_elem)) {
            fdata     = stk::mesh::field_data( *eMesh.m_refine_level_field , parent_elem );
          }
          if (fdata && fdata_new)
            fdata_new[0] = fdata[0] + 1;
          //std::cout << "fdata= " << fdata << " fdata_new= " << fdata_new[0] << std::endl;
        }

      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);
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

          if (eMesh.is_valid(familyTreeNewElement))
            {
              family_tree_id = eMesh.identifier(familyTreeNewElement);
            }
          else
            {
              family_tree_id = eMesh.getNextId(FAMILY_TREE_RANK);
            }

          // FIXME
          if (0 && eMesh.entity_rank(parent_elem) != stk::topology::ELEMENT_RANK)
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

          if (eMesh.is_valid(familyTreeNewElement))
            {
              family_tree = familyTreeNewElement;
              //std::cout << "familyTreeNewElement = " << eMesh.identifier(familyTreeNewElement) << std::endl;
            }
          else
            {
              //std::cout << "family_tree_id = " << family_tree_id << std::endl;
              family_tree = eMesh.get_bulk_data()->declare_constraint(family_tree_id, add);
            }

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
          // FIXME - need a unit test to check if this ever breaks in the future
          //family_tree = parent_to_family_tree_relations[FAMILY_TREE_LEVEL_0].entity();

          unsigned parent_elem_ft_level_1 = eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_1, parent_elem);
          family_tree = parent_to_family_tree_relations[parent_elem_ft_level_1].entity();

        }
      else
        {
          throw std::logic_error("UniformRefinerPatternBase::set_parent_child_relations no family_tree");
        }

      if (1)
        {
          percept::MyPairIterRelation family_tree_relations (eMesh, family_tree, eMesh.entity_rank(parent_elem));
          VERIFY_OP_ON(family_tree_relations.size(), ==, ordinal+1, "bad family_tree_relations");
        }
      eMesh.get_bulk_data()->declare_relation(family_tree, newElement, ordinal + 1);  // the + 1 here is to give space for the parent
      percept::MyPairIterRelation parent_to_family_tree_relations1 (eMesh, parent_elem, FAMILY_TREE_RANK);
    }

#define DEBUG_SET_NEEDED_PARTS 0
#define REMOVE_UNDERSCORE_FROM_TOPO_NAME 0

    void UniformRefinerPatternBase::
    addToBreakPatternList(std::set<UniformRefinerPatternBase *>& list, PerceptMesh& eMesh)
    {
      list.insert(this);
      if (DEBUG_SET_NEEDED_PARTS)
        std::cout << "addToBreakPatternList this= " << eMesh.demangle(typeid(*this).name()) << std::endl;
      std::vector<UniformRefinerPatternBase *> bp;
      setSubPatternsForSetNeededParts(bp, eMesh);
      if (DEBUG_SET_NEEDED_PARTS)
        std::cout << "addToBreakPatternList bp.size= " << bp.size() << std::endl;
      for (unsigned ibp=0; ibp < bp.size(); ibp++)
        {
          if (bp[ibp] == this)
            continue;
          if (list.find(bp[ibp]) == list.end())
            {
              list.insert(bp[ibp]);
              if (bp[ibp] == 0)
                {
                  std::cout << "addToBreakPatternList this bad= " << eMesh.demangle(typeid(*this).name()) << std::endl;
                  setSubPatternsForSetNeededParts(bp, eMesh);
                }

#ifndef __clang__
              if (DEBUG_SET_NEEDED_PARTS)
                std::cout << "addToBreakPatternList bp recurse= " << eMesh.demangle(typeid(*bp[ibp]).name()) << std::endl;
#endif
              bp[ibp]->addToBreakPatternList(list, eMesh);
            }
        }
    }

    void UniformRefinerPatternBase::addRefineNewNodesPart(percept::PerceptMesh& eMesh)
    {
    	stk::mesh::Part* new_nodes_part = eMesh.get_non_const_part("refine_new_nodes_part");
    	if (!new_nodes_part)
    	{
    		stk::mesh::Part& part = eMesh.get_fem_meta_data()->declare_part_with_topology("refine_new_nodes_part", stk::topology::NODE);
    		stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
    		meta.declare_attribute_no_delete(part, &percept::auto_part);
    	}
    }

    void UniformRefinerPatternBase::addActiveParentParts(percept::PerceptMesh& eMesh)
    {
    	stk::mesh::EntityRank part_ranks[] = {eMesh.element_rank(), eMesh.side_rank()};
    	for (unsigned irank=0; irank < 2; irank++)
    	{
    		std::string active_part_name = "refine_active_elements_part_"+toString(part_ranks[irank]);
    		std::string inactive_part_name = "refine_inactive_elements_part_"+toString(part_ranks[irank]);
    		stk::mesh::Part* active_elements_part = eMesh.get_non_const_part(active_part_name);
    		if (!active_elements_part)
    		{
    			stk::mesh::Part& part = eMesh.get_fem_meta_data()->declare_part(active_part_name, part_ranks[irank]);
    			stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
    			meta.declare_attribute_no_delete(part, &percept::auto_part);
    		}
    		stk::mesh::Part* inactive_elements_part = eMesh.get_non_const_part(inactive_part_name);
    		if (!inactive_elements_part)
    		{
    			stk::mesh::Part& part = eMesh.get_fem_meta_data()->declare_part(inactive_part_name, part_ranks[irank]);
    			stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
    			meta.declare_attribute_no_delete(part, &percept::auto_part);
    		}
    	}
    }

    bool UniformRefinerPatternBase::foundIncludeOnlyBlock(percept::PerceptMesh& eMesh, std::vector<std::string>& block_names_include)
    {
    	bool found_include_only_block = false;

        stk::mesh::PartVector all_parts = eMesh.get_fem_meta_data()->get_parts();

        for (unsigned ib = 0; ib < block_names_include.size(); ib++)
    	{
    		bool foundPart = false;
    		for (stk::mesh::PartVector::iterator i_part = all_parts.begin(); i_part != all_parts.end(); ++i_part)
    		{
    			stk::mesh::Part * part = *i_part ;

    			std::string bname = block_names_include[ib];
    			if ('+' == bname[0])
    				found_include_only_block = true;
    			bname = bname.substr(1, bname.length()-1);
    			//std::cout << "bname= " << bname << " part= " << part->name() << std::endl;
    			if (eMesh.checkForPartNameWithAliases(*part, bname))
    			{
    				foundPart = true;
    				break;
    			}
    		}
    		if (!foundPart)
    		{
    			std::string msg = "UniformRefinerPattern::setNeededParts unknown block name: " + block_names_include[ib];
    			throw std::runtime_error(msg.c_str());
    		}
    	}

		return found_include_only_block;
    }

    void UniformRefinerPatternBase::addOldPart(percept::PerceptMesh& eMesh)
    {
      stk::mesh::PartVector all_parts = eMesh.get_fem_meta_data()->get_parts();
      std::string oldPartName = UniformRefinerPatternBase::getOldElementsPartName() + toString(m_primaryEntityRank);
      bool foundOldPart = false;
      for (stk::mesh::PartVector::iterator i_part = all_parts.begin(); i_part != all_parts.end(); ++i_part)
        {
          stk::mesh::Part *  part = *i_part ;
          if (oldPartName == part->name())
            {
              foundOldPart = true;
              break;
            }
        }

      if (!foundOldPart)
        {
          if (DEBUG_SET_NEEDED_PARTS) std::cout << "tmp setNeededParts:: declare_part for oldPartName = "
                                                << oldPartName << " rank= " << m_primaryEntityRank << std::endl;
          stk::mesh::Part& part = eMesh.get_fem_meta_data()->declare_part(oldPartName, m_primaryEntityRank);
          stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
          meta.declare_attribute_no_delete(part, &percept::auto_part);
        }
    }

    bool UniformRefinerPatternBase::shouldDoThisPart(percept::PerceptMesh& eMesh, BlockNamesType block_names_ranks,
    		bool found_include_only_block, std::vector<std::string>& block_names_include, stk::mesh::Part *  part)
    {
    	bool doThisPart = (block_names_ranks[eMesh.element_rank()].size() == 0);

    	if (!doThisPart)
      {
        // we found one block with a "+", so this means include only the actual specified list of blocks, except for those excluded with "-"
        if (found_include_only_block)
          {
            doThisPart = false;
            for (unsigned ib = 0; ib < block_names_include.size(); ib++)
              {
                std::string bname = block_names_include[ib];
                if ('+' == bname[0])
                  {
                    bname = bname.substr(1, bname.length()-1);
                    if (eMesh.checkForPartNameWithAliases(*part, bname))
                      {
                        doThisPart = true;
                        break;
                      }
                  }
              }
          }
        else
          // do them all, except for excludes
          {
            doThisPart = true;
          }

        // check for excludes
        if (doThisPart)
          {
            if (block_names_ranks[m_primaryEntityRank].size() == 0)
              doThisPart = false;
            else
              {
                for (unsigned ib = 0; ib < block_names_include.size(); ib++)
                  {
                    std::string bname = block_names_include[ib];
                    if ('-' == bname[0])
                      {
                        bname = bname.substr(1, bname.length()-1);
                        if (eMesh.checkForPartNameWithAliases(*part, bname))
                          {
                            doThisPart = false;
                            break;
                          }
                      }
                  }
              }
          }
      }
    bool isOldElementsPart = ( (part->name()).find(UniformRefinerPatternBase::m_oldElementsPartName) != std::string::npos);
    doThisPart = doThisPart && ( part->primary_entity_rank() == m_primaryEntityRank );
    doThisPart = doThisPart && !isOldElementsPart;

    bool isConvertedPart = ( (part->name()).find(getAppendConvertString()) != std::string::npos);
    if (DEBUG_SET_NEEDED_PARTS)
      {
        std::cout << "isOldElementsPart= " << isOldElementsPart
                  << " m_primaryEntityRank= " << m_primaryEntityRank
                  << " part->primary_entity_rank() = " << part->primary_entity_rank()
                  << " isConvertedPart= " << isConvertedPart
                  << std::endl;
      }

    if (!isOldElementsPart)
      {
        unsigned my_cellTopoKey = getFromTypeKey();
        const CellTopologyData * part_cell_topo_data = eMesh.get_cell_topology(*part);
        if (part_cell_topo_data)
          {
            shards::CellTopology topo(part_cell_topo_data);
            doThisPart = doThisPart && (topo.getKey() == my_cellTopoKey);
          }
        else
          {
            if (part->subsets().size() != 0)
              {
                doThisPart = false;
              }
          }
        if ((DEBUG_SET_NEEDED_PARTS) && doThisPart)
          std::cout << "tmp srk 1 SNP setNeededParts:: "
            //<< eMesh.demangle(typeid(*this).name())
                    << "  part name= " << part->name()
                    << "  doThisPart= " << doThisPart
            // << "  part->primary_entity_rank() = " <<  part->primary_entity_rank()
            // << "  my_cellTopoKey= " << my_cellTopoKey
            // << "  topo.getKey() = " << topo.getKey()
            // << "  topo.getName() = " << topo.getName()
                    << std::endl;
      }

    return doThisPart;
    }

    void UniformRefinerPatternBase::setNeededParts_debug1(percept::PerceptMesh& eMesh)
    {
        if (DEBUG_SET_NEEDED_PARTS)
          {
            stk::mesh::PartVector all_parts = eMesh.get_fem_meta_data()->get_parts();
            std::cout << "\n\n====>\ntmp 0 setNeededParts: for FromTopo= " << getFromTopoPartName() << " ToTopo= " << getToTopoPartName();
            for (stk::mesh::PartVector::iterator i_part = all_parts.begin(); i_part != all_parts.end(); ++i_part)
              {
                stk::mesh::Part *  part = *i_part ;
                //if ( stk::mesh::is_auto_declared_part(*part) )
                //  continue;
                std::cout << " part= " <<  part->name() << std::endl;
              }
            std::cout << std::endl;
          }
    }

    void UniformRefinerPatternBase::setNeededParts_debug2()
    {
    	if (DEBUG_SET_NEEDED_PARTS)
    	{
    		std::cout << "tmp srk SNP m_toParts.size() = " << m_toParts.size() << " m_fromParts.size()= " << m_fromParts.size() << std::endl;
    		for (unsigned i=0; i < m_toParts.size(); i++)
    		{
    			std::cout << "tmp srk SNP m_toParts[" << i << "]= " << m_toParts[i]->name() << std::endl;
    		}
    		for (unsigned i=0; i < m_fromParts.size(); i++)
    		{
    			std::cout << "tmp srk SNP m_fromParts[" << i << "]= " << m_fromParts[i]->name() << std::endl;
    		}
    	}
    }

    void UniformRefinerPatternBase::updateSurfaceBlockMap(percept::PerceptMesh& eMesh, stk::mesh::Part* part,
                                                          stk::mesh::Part* part_to)
    {
      std::vector<const stk::mesh::Part *> surfaces = eMesh.get_fem_meta_data()->get_surfaces_in_surface_to_block_map();

      // case 1: part/part_to are blocks
      if (      part->primary_entity_rank() == stk::topology::ELEMENT_RANK
          && part_to->primary_entity_rank() == stk::topology::ELEMENT_RANK)
      {
        // Add the refined block to all of the mappings that already contain the parent block
        for (auto iter = surfaces.begin(); iter != surfaces.end(); ++iter)
        {
          const stk::mesh::Part * surface = *iter;
          std::vector<const stk::mesh::Part*> blocks = eMesh.get_fem_meta_data()->get_blocks_touching_surface(surface);
          const bool parentBlockInMapping = (find(blocks.begin(), blocks.end(), part) != blocks.end());
          if (parentBlockInMapping)
          {
            const bool blockNotAlreadyInMapping = (find(blocks.begin(), blocks.end(), part_to) == blocks.end());
            if (blockNotAlreadyInMapping)
            {
              std::vector<const stk::mesh::Part*> new_blocks = blocks;
              new_blocks.push_back(part_to);
              eMesh.get_fem_meta_data()->set_surface_to_block_mapping(surface, new_blocks);
            }
          }
        }

        if (part->topology() != part_to->topology()) {
          // If the element topology changes upon refinement, also map the refined block to all
          // surface parts of the appropriate topology that are also subsets of surfaces
          // that are mapped to the "parent" block.  We have to cascade the surface/block
          // mappings down through all topology changes while respecting the original
          // surface/block mapping.

          std::set<stk::topology> faceTopologies;
          for (unsigned i = 0; i < part_to->topology().num_faces(); ++i)
          {
            faceTopologies.insert(part_to->topology().face_topology(i));
          }

          for (const stk::mesh::Part * surface : surfaces)
          {
            std::vector<const stk::mesh::Part*> blocks = eMesh.get_fem_meta_data()->get_blocks_touching_surface(surface);
            const bool parentBlockInMapping = find(blocks.begin(), blocks.end(), part) != blocks.end();
            if (parentBlockInMapping)
            {
              for (const stk::mesh::Part * subsetSurface : surfaces)
              {
                const stk::mesh::PartVector supersetSurfaces = subsetSurface->supersets();
                const bool isSubsetOfAParentBlockSurface = (std::find(supersetSurfaces.begin(), supersetSurfaces.end(), surface) != supersetSurfaces.end());
                const bool subsetSurfaceTopologyMatchesBlock = (std::find(faceTopologies.begin(), faceTopologies.end(), subsetSurface->topology()) != faceTopologies.end());
                if (isSubsetOfAParentBlockSurface && subsetSurfaceTopologyMatchesBlock)
                {
                  std::vector<const stk::mesh::Part*> subsetSurfaceBlocks = eMesh.get_fem_meta_data()->get_blocks_touching_surface(subsetSurface);
                  const bool blockNotAlreadyInMapping = (find(subsetSurfaceBlocks.begin(), subsetSurfaceBlocks.end(), part_to) == subsetSurfaceBlocks.end());
                  if (blockNotAlreadyInMapping)
                  {
                    subsetSurfaceBlocks.push_back(part_to);
                    eMesh.get_fem_meta_data()->set_surface_to_block_mapping(subsetSurface, subsetSurfaceBlocks);
                  }
                }
              }
            }
          }
        }
      }

      // case 2: part/part_to are surfaces
      if (part->primary_entity_rank() == eMesh.get_fem_meta_data()->side_rank() &&
          part_to->primary_entity_rank() == eMesh.get_fem_meta_data()->side_rank()
          && find(surfaces.begin(), surfaces.end(), part) != surfaces.end())
        {
          std::vector<const stk::mesh::Part*> blocks = eMesh.get_fem_meta_data()->get_blocks_touching_surface(part);
          eMesh.get_fem_meta_data()->set_surface_to_block_mapping(part_to, blocks);
        }
    }

    void addDistributionFactorToNewPart(stk::mesh::MetaData & meta, stk::mesh::Part * old_part, stk::mesh::Part * new_part)
    {
      stk::mesh::FieldBase *df_field = const_cast<stk::mesh::FieldBase*>(stk::io::get_distribution_factor_field(*old_part));
      if (df_field) {
          stk::io::set_field_role(*df_field, Ioss::Field::MESH);
          stk::io::set_distribution_factor_field(*new_part, *df_field);
          int side_node_count = new_part->topology().num_nodes();
          stk::mesh::put_field_on_mesh(*df_field, *new_part, side_node_count, nullptr);
      }
    }

    void UniformRefinerPatternBase::setNeededParts(percept::PerceptMesh& eMesh, BlockNamesType block_names_ranks, bool sameTopology, bool skipConvertedParts)
    {
      EXCEPTWATCH;

      if (DEBUG_SET_NEEDED_PARTS)
        std::cout << "\n\n ============= setNeededParts start \n\n " << PerceptMesh::demangle(typeid(*this).name()) << std::endl;

      addRefineNewNodesPart(eMesh);

      addActiveParentParts(eMesh);

      if (block_names_ranks.size() == 0)
      {
        block_names_ranks.resize(percept::EntityRankEnd);
      }

      m_fromParts.resize(0);
      m_toParts.resize(0);

      setNeededParts_debug1(eMesh);

      std::vector<std::string>& block_names_include = block_names_ranks[m_primaryEntityRank];

      stk::mesh::PartVector all_parts = eMesh.get_fem_meta_data()->get_parts();

      bool found_include_only_block = foundIncludeOnlyBlock(eMesh, block_names_include);

      for (stk::mesh::PartVector::iterator i_part = all_parts.begin(); i_part != all_parts.end(); ++i_part)
      {
        stk::mesh::Part *  part = *i_part ;

        if ( stk::mesh::is_auto_declared_part(*part) )
          continue;
        bool is_auto_part = part->attribute<AutoPart>() != 0;
        if (is_auto_part)
          continue;

        bool doThisPart = shouldDoThisPart(eMesh, block_names_ranks,
            found_include_only_block, block_names_include, part);

        if (!doThisPart) continue;

        stk::mesh::EntityRank switch_part_primary_entity_rank  =  part->primary_entity_rank() ;

        if (switch_part_primary_entity_rank == eMesh.edge_rank() ||
            switch_part_primary_entity_rank == stk::topology::ELEMENT_RANK ||
            switch_part_primary_entity_rank == eMesh.face_rank())
        {
          stk::mesh::Part *  block_to=0;
          if (sameTopology)
          {
            block_to = part;
          }
          else
          {
            std::string toTopoPartName = getToTopoPartName();
            if (REMOVE_UNDERSCORE_FROM_TOPO_NAME) Util::replace(toTopoPartName,"_","");
            std::string newPartName = part->name() + getConvertSeparatorString() + toTopoPartName + getConvertSeparatorString() + getAppendConvertString();
            block_to = &eMesh.get_fem_meta_data()->declare_part(newPartName, part->primary_entity_rank());
            if (DEBUG_SET_NEEDED_PARTS) std::cout << "tmp setNeededParts:: declare_part name= " << newPartName
                << " with topo= " << getToTopoPartName() << std::endl;
            stk::mesh::set_topology(*block_to, stk::mesh::get_topology(shards::CellTopology(getToTopology()), eMesh.get_fem_meta_data()->spatial_dimension()));
            eMesh.get_fem_meta_data()->set_part_id(*block_to, part->id());
            if (!stk::io::is_part_io_part(block_to)) {
              stk::io::put_io_part_attribute(*block_to);
            }

            addDistributionFactorToNewPart(*eMesh.get_fem_meta_data(), part, block_to);

            updateSurfaceBlockMap(eMesh, part, block_to);

            if (switch_part_primary_entity_rank == stk::topology::ELEMENT_RANK) {
              // Add the new surface parts needed by topology change upon refinement and make
              // another pass at updating the surface to block mapping for these new parts
              fixSubsets(eMesh);
              addExtraSurfaceParts(eMesh);
              updateSurfaceBlockMap(eMesh, part, block_to);
            }

            stk::mesh::PartVector *pv  = const_cast<stk::mesh::PartVector *>(part->attribute<stk::mesh::PartVector>());
            if (pv == 0)
            {
              pv = new stk::mesh::PartVector;
              eMesh.get_fem_meta_data()->declare_attribute_with_delete(*part, pv);
            }
            pv->push_back(block_to);
            if (DEBUG_SET_NEEDED_PARTS)
            {
              for (unsigned ii=0; ii < pv->size(); ii++)
              {
                std::cout << "tmp srk part.attr = " << part->name() << " block_to= " << (*pv)[ii]->name() << std::endl;
              }
            }
          }

          if (!((part->name()).find(UniformRefinerPatternBase::getOldElementsPartName()) != std::string::npos))
          {
            if (DEBUG_SET_NEEDED_PARTS) std::cout << "tmp setNeededParts:: fromPart = " << part->name() << " toPart = " << block_to->name() << std::endl;
            m_fromParts.push_back(part);
            m_toParts.push_back(block_to);
          }
        }
      }

      if (!sameTopology) fixSubsets(eMesh);
      addExtraSurfaceParts(eMesh);

      addOldPart(eMesh);

      setNeededParts_debug2();
    }

#define DEBUG_fixSubSets 0
    void UniformRefinerPatternBase::fixSubsets(percept::PerceptMesh& eMesh)
    {
      if (DEBUG_fixSubSets)
        printParts(this, false);

      for (unsigned i_fromPart = 0; i_fromPart < m_fromParts.size(); i_fromPart++)
        {
          stk::mesh::Part& fromPart = *m_fromParts[i_fromPart];

          stk::mesh::EntityRank switch_fromPart_primary_entity_rank  = fromPart.primary_entity_rank();

          //case
          if(switch_fromPart_primary_entity_rank == eMesh.edge_rank() ||
             switch_fromPart_primary_entity_rank == eMesh.face_rank() )
            {
              stk::mesh::Part& toPart = *m_toParts[i_fromPart];

              const stk::mesh::PartVector from_subsets = fromPart.subsets();
              const stk::mesh::PartVector from_supersets = fromPart.supersets();
              if (DEBUG_fixSubSets && eMesh.get_rank() == 0)
                std::cout << "fixSubsets:: fromPart= " << fromPart.name() << " " << fromPart.topology()
                          << " toPart= " << toPart.name() << " " << toPart.topology()
                          << " from_subset.size= " << from_subsets.size()
                          << " from_supersets.size= " << from_supersets.size()
                          << std::endl;

              if (from_subsets.size() == 0 && from_supersets.size())
                {
                  for (unsigned i_from_supersets = 0; i_from_supersets < from_supersets.size(); ++i_from_supersets)
                    {
                      stk::mesh::Part& from_superset = *from_supersets[i_from_supersets];
                      if ( stk::mesh::is_auto_declared_part(from_superset) )
                        continue;
                      bool auto_part = from_superset.attribute<AutoPart>() != 0;
                      if (auto_part)
                        continue;

                      if (DEBUG_fixSubSets && eMesh.get_rank() == 0)
                        std::cout << "fixSubsets:: fromPart= " << fromPart.name() << " " << fromPart.topology()
                                  << " toPart= " << toPart.name() << " " << toPart.topology()
                                  << " from_superset= " << from_superset.name() << " " << from_superset.topology()
                                  << std::endl;

                      if (DEBUG_fixSubSets && eMesh.get_rank() == 0)
                        std::cout << "fixSubsets:: declare_part_subset from_superset = " << from_superset.name()  << " toPart.name= " <<  toPart.name() << std::endl;

                      eMesh.get_fem_meta_data()->declare_part_subset(from_superset, toPart );

                    }
                }

              for (unsigned i_from_subset = 0; i_from_subset < from_subsets.size(); i_from_subset++)
                {
                  stk::mesh::Part& from_subset = *from_subsets[i_from_subset];

                  if (DEBUG_fixSubSets && eMesh.get_rank() == 0)
                    std::cout << "fixSubsets:: fromPart= " << fromPart.name() << " " << fromPart.topology()
                              << " toPart= " << toPart.name() << " " << toPart.topology()
                              << " from_subset= " << from_subset.name() << " " << from_subset.topology()
                              << std::endl;

                  const CellTopologyData * from_subset_part_cell_topo_data = eMesh.get_cell_topology(from_subset);
                  const CellTopologyData * to_subset_part_cell_topo_data = eMesh.get_cell_topology(toPart);
                  if (!from_subset_part_cell_topo_data || !to_subset_part_cell_topo_data)
                    continue;
                  if (from_subset_part_cell_topo_data == to_subset_part_cell_topo_data)
                    continue;

                  std::string to_subset_part_cell_topo_data_name = to_subset_part_cell_topo_data->name;
                  if (REMOVE_UNDERSCORE_FROM_TOPO_NAME) Util::replace(to_subset_part_cell_topo_data_name, "_", "");
                  std::string to_subset_name = from_subset.name() + getConvertSeparatorString() + to_subset_part_cell_topo_data_name + getConvertSeparatorString() + getAppendConvertString();

                  stk::mesh::Part* to_subset_p = eMesh.get_fem_meta_data()->get_part(to_subset_name);
                  if (DEBUG_fixSubSets && eMesh.get_rank() == 0)
                    std::cout << "fixSubsets:: declare_part_subset toPart is superset = " << toPart.name()  << " to_subset_name= " <<  to_subset_name << std::endl;
                  if (!to_subset_p)
                    {
                      printParts(this, true);
                    }
                  VERIFY_OP_ON(to_subset_p, !=, 0, std::string("fixSubsets couldn't find part error, part= ")+to_subset_name);
                  stk::mesh::Part& to_subset = *to_subset_p;

                  eMesh.get_fem_meta_data()->declare_part_subset(toPart, to_subset);
                }
            }

        }
    }

    void UniformRefinerPatternBase::addExtraSurfaceParts(percept::PerceptMesh& eMesh)
    {
      if (eMesh.get_spatial_dim() == 2)
        return;

      if (eMesh.getProperty("Refiner_skip_side_part_fixes") == "true")
        return;

      std::vector<const stk::mesh::Part*> surfaces = eMesh.get_fem_meta_data()->get_surfaces_in_surface_to_block_map();
      for (unsigned isu = 0; isu < surfaces.size(); ++isu)
        {
          stk::mesh::Part * from_superset = const_cast<stk::mesh::Part *>(surfaces[isu]);
          VERIFY_OP_ON(from_superset, !=, 0, "error: null surface Part pointer found in surface to block map");
          stk::mesh::PartVector subsets = from_superset->subsets();

          for (auto subset_surface : subsets)
            {
              std::vector<const stk::mesh::Part*> blocks = eMesh.get_fem_meta_data()->get_blocks_touching_surface(subset_surface);
              stk::topology side_topo = subset_surface->topology();
              VERIFY_OP_ON(side_topo, !=, stk::topology::INVALID_TOPOLOGY, "bad side_topo = "+subset_surface->name() + " from_superset= " << from_superset->name());

              std::string ioss_side_topo = "";
              convert_stk_topology_to_ioss_name(side_topo, ioss_side_topo);

              for (auto block : blocks)
                {
                  VERIFY_OP_ON(block, !=, 0, "error: null block Part pointer found in surface to block map for surface=" << from_superset->name());
                  stk::topology elem_topo = block->topology();
                  VERIFY_OP_ON(elem_topo, !=, stk::topology::INVALID_TOPOLOGY, "bad elem_topo");

                  std::string ioss_elem_topo = "";
                  convert_stk_topology_to_ioss_name(elem_topo, ioss_elem_topo);

                  bool is_io_part = stk::io::is_part_io_part(*from_superset);

                  // need to create independent of toPart topology
                  if (is_io_part && from_superset->id()>0)
                    {
                      // Name of form: "name_eltopo_sidetopo_id" or
                      //               "name_block_id_sidetopo_id"
                      // "name" is typically "surface".

                      const std::string id = std::to_string(from_superset->id());

                      const std::string new_part_name = "surface_" + ioss_elem_topo + "_"+ ioss_side_topo +"_" + id;
                      stk::mesh::Part* new_part = eMesh.get_fem_meta_data()->get_part(new_part_name);

                      if (new_part) addDistributionFactorToNewPart(*eMesh.get_fem_meta_data(), from_superset, new_part);

                      // alternate part naming:
                      //       "name_block_id_sidetopo_id"

                      const std::string new_part_name_2 = "surface_" + block->name() + "_" + ioss_side_topo + "_" + id;
                      stk::mesh::Part* new_part_type_2 = eMesh.get_fem_meta_data()->get_part(new_part_name_2);

                      if (new_part_type_2) addDistributionFactorToNewPart(*eMesh.get_fem_meta_data(), from_superset, new_part_type_2);

                      if (!new_part && !new_part_type_2)
                        {
                          new_part = &eMesh.get_fem_meta_data()->declare_part(new_part_name, eMesh.side_rank());
                          stk::mesh::set_topology(*new_part, side_topo);
                          eMesh.get_fem_meta_data()->declare_part_subset(*from_superset, *new_part);
                          eMesh.get_fem_meta_data()->set_part_id(*new_part, from_superset->id());
                          stk::io::put_io_part_attribute(*new_part);

                          addDistributionFactorToNewPart(*eMesh.get_fem_meta_data(), from_superset, new_part);

                          updateSurfaceBlockMap(eMesh, from_superset, new_part);
                        }
                    }
                }
            }
        }
    }

    void UniformRefinerPatternBase::
    genericEnrich_createNewElementsBase(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                                        const unsigned fromTopoKey_in, const unsigned toTopoKey_in,
                                        const int ToTopology_node_count,
                                        const int FromTopology_vertex_count,
                                        Elem::CellTopology elem_celltopo_in,
                                        const CellTopologyData * const cell_topo_data_toTopo_in,
                                        vector< vector<stk::mesh::EntityId> >& elems,
                                        stk::mesh::FieldBase *proc_rank_field)
    {
      percept::PerceptMesh& m_eMesh = eMesh;
      std::vector<NeededEntityType> needed_entities;
      fillNeededEntities(needed_entities);

      const CellTopologyData * const cell_topo_data = eMesh.get_cell_topology(element);

      shards::CellTopology cell_topo(cell_topo_data);
      const percept::MyPairIterRelation elem_nodes (eMesh, element, stk::topology::NODE_RANK);

      std::vector<stk::mesh::Part*> add_parts;
      std::vector<stk::mesh::Part*> remove_parts;

      unsigned n_edges = cell_topo_data->edge_count;

      unsigned n_faces = cell_topo.getFaceCount();
      if (n_faces == 0) n_faces = 1; // 2D face has one "face"

      add_parts = m_toParts;

#if STK_ADAPT_URP_LOCAL_NODE_COMPS
      for (unsigned i_need = 0; i_need < needed_entities.size(); i_need++)
        {
          unsigned nSubDimEntities = 1;
          if (needed_entities[i_need].first == eMesh.edge_rank())
            {
              nSubDimEntities = cell_topo_data->edge_count;
            }
          else if (needed_entities[i_need].first == m_eMesh.face_rank())
            {
              nSubDimEntities = cell_topo_data->side_count;
            }

          for (unsigned iSubDim = 0; iSubDim < nSubDimEntities; iSubDim++)
            {
              nodeRegistry.prolongateCoords(*const_cast<stk::mesh::Entity>(&element), needed_entities[i_need].first, iSubDim);
              nodeRegistry.addToExistingParts(*const_cast<stk::mesh::Entity>(&element), needed_entities[i_need].first, iSubDim);
              nodeRegistry.prolongateFields(*const_cast<stk::mesh::Entity>(&element), needed_entities[i_need].first, iSubDim);
            }
        }
#endif

      //const CellTopologyData * const cell_topo_data_toTopo = shards::getCellTopologyData< ToTopology >();
      const CellTopologyData * const cell_topo_data_toTopo = cell_topo_data_toTopo_in;
      shards::CellTopology cellTopo(cell_topo_data_toTopo);

#define CENTROID_N NN(m_primaryEntityRank,0)

      vector<stk::mesh::EntityId>& EN = elems[0];

      for (int ind = 0; ind < FromTopology_vertex_count; ind++)
        {
          EN[ind] = VERT_N(ind);
        }

      //std::cout << "ToTopology::vertex_count = " << ToTopology::vertex_count << std::endl;
      for (unsigned i_need = 0; i_need < needed_entities.size(); i_need++)
        {
          if (needed_entities[i_need].first == m_eMesh.edge_rank())
            {
              for (unsigned i_edge = 0; i_edge < n_edges; i_edge++)
                {
                  unsigned edge_ord = cell_topo_data_toTopo->edge[i_edge].node[2];
                  stk::mesh::EntityId inode = EDGE_N(i_edge);
                  EN[edge_ord] = inode;
                }
            }
          else if (needed_entities[i_need].first == m_eMesh.face_rank())
            {
              for (unsigned i_face = 0; i_face < n_faces; i_face++)
                {
                  // FIXME assumes face is quadrilateral
                  shards::CellTopology face_topo = cell_topo.getDimension()==2 ? cell_topo : shards::CellTopology(cell_topo.getCellTopologyData( 2, i_face));
                  if (face_topo.getNodeCount() == 4 || toTopoKey_in == topo_key_quad9)
                    {
                      unsigned face_ord = 0;
                      if (toTopoKey_in == topo_key_quad9)
                        {
                          face_ord = 8;
                        }
                      else
                        {
                          face_ord = cell_topo_data_toTopo->side[i_face].node[8];
                        }

                      stk::mesh::EntityId inode = FACE_N(i_face);

                      //std::cout << "tmp P[" << eMesh.get_rank() << "] inode = " << inode << " for i_face = " << i_face << " face_ord= " << face_ord << std::endl;

                      if (!inode)
                        {
                          std::cout << "P[" << eMesh.get_rank() << "] inode = 0 for i_face = " << i_face << " face_ord= " << face_ord << std::endl;
                          //throw std::logic_error("UniformRefinerPatternBase::genericEnrich_createNewElements bad entity id = 0 ");
                        }

                      EN[face_ord] = inode;
                    }
                }
            }
          else if (needed_entities[i_need].first == stk::topology::ELEMENT_RANK)
            {
              const unsigned centroid_node       = (toTopoKey_in == topo_key_quad9 ? 8 :
                                                    (toTopoKey_in == topo_key_hex27 ? 20 : 0)
                                                    );

              EN[ centroid_node ] = CENTROID_N;
            }
        }

#undef CENTROID_N

      bool use_declare_element_side = UniformRefinerPatternBase::USE_DECLARE_ELEMENT_SIDE &&  getPrimaryEntityRank() == eMesh.side_rank();

      std::vector<stk::mesh::Entity> nodes(ToTopology_node_count);

      for (unsigned ielem=0; ielem < elems.size(); ielem++)
        {

          stk::mesh::Entity newElement = stk::mesh::Entity();
          if (!use_declare_element_side)
            {
              newElement = *element_pool;
            }

          for (int inode=0; inode < ToTopology_node_count; inode++)
            {
              stk::mesh::EntityId eid = elems[ielem][inode];
              if (!eid)
                {
                  throw std::logic_error("UniformRefinerPatternBase::genericEnrich_createNewElements bad entity id = 0 ");
                }
              stk::mesh::Entity node = createOrGetNode(nodeRegistry, eMesh, eid);
              nodes[inode] = node;
            }

          create_side_element(eMesh, use_declare_element_side, nodes.data(), static_cast<unsigned>(ToTopology_node_count), newElement);

          // FIXME
          if (m_primaryEntityRank == stk::topology::ELEMENT_RANK &&  proc_rank_field)
            {
              double *fdata = eMesh.field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
              fdata[0] = double(eMesh.owner_rank(newElement));
            }


          change_entity_parts(eMesh, element, newElement);

          if (0 & EXTRA_PRINT_URP_IF)
            {
              std::cout << "tmp newElement: " << std::endl;
              eMesh.print_entity(std::cout, newElement, eMesh.get_coordinates_field() );
            }

          set_parent_child_relations(eMesh, element, newElement, *ft_element_pool, ielem);

          std::vector<stk::mesh::Entity> elements(1,element);
          eMesh.prolongateElementFields( elements, newElement);

          ft_element_pool++;
          if (!use_declare_element_side)
            element_pool++;

        }

    }

    void UniformRefinerPatternBase::
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
                                        stk::mesh::FieldBase *proc_rank_field)
    {
      EXCEPTWATCH;
      percept::PerceptMesh& m_eMesh = eMesh;
      static std::vector<NeededEntityType> needed_entities;
      fillNeededEntities(needed_entities);

      const CellTopologyData * const cell_topo_data = eMesh.get_cell_topology(element);

      VERIFY_OP_ON(cell_topo_data, !=, 0, "bad topology");

      shards::CellTopology cell_topo(cell_topo_data);
      bool isLinearElement = Util::isLinearElement(cell_topo);

      // SPECIAL CASE ALERT  FIXME
      //if (toTopoKey == topo_key_wedge15)
      //  isLinearElement = true;
      //std::cout << "tmp cell_topo= " << cell_topo.getName() << " isLinearElement= " << isLinearElement << std::endl;

      const percept::MyPairIterRelation elem_nodes (eMesh, element, stk::topology::NODE_RANK);

      int topoDim = cell_topo.getDimension();
      //unsigned cell_topo_key = fromTopoKey;
      //getTopoDim(topoDim, cell_topo_key);
      unsigned cellDimension = (unsigned)topoDim;

      unsigned n_edges = cell_topo_data->edge_count;
      if (n_edges == 0) n_edges = 1; // 1D edge has one "edge"
      unsigned n_faces = cell_topo.getFaceCount();
      if (n_faces == 0) n_faces = 1; // 2D face has one "face"

      for (unsigned i_need = 0; i_need < needed_entities.size(); i_need++)
        {
          unsigned nSubDimEntities = 0;
          if (needed_entities[i_need].first == eMesh.edge_rank())
            {
              nSubDimEntities = cell_topo_data->edge_count;
            }
          else if (needed_entities[i_need].first == eMesh.face_rank())
            {
              nSubDimEntities = cell_topo_data->side_count;
            }
          else if (needed_entities[i_need].first == stk::topology::ELEMENT_RANK)
            {
              nSubDimEntities = 1;
            }

          // FIXME - assumes first node on each sub-dim entity is the "linear" one
          for (unsigned iSubDim = 0; iSubDim < nSubDimEntities; iSubDim++)
            {
              //!
#if STK_ADAPT_URP_LOCAL_NODE_COMPS
              nodeRegistry.addToExistingParts(*const_cast<stk::mesh::Entity>(&element), needed_entities[i_need].first, iSubDim);
              if (isLinearElement)
                {
                  nodeRegistry.prolongateFields(*const_cast<stk::mesh::Entity>(&element), needed_entities[i_need].first, iSubDim);
                }
#endif
            }
        }

      const CellTopologyData * const cell_topo_data_toTopo = cell_topo_data_toTopo_in;
      VERIFY_OP_ON(cell_topo_data_toTopo, !=, 0, "bad topo 2");
      shards::CellTopology cellTopo(cell_topo_data_toTopo);

      Elem::CellTopology elem_celltopo = elem_celltopo_in;
      const Elem::RefinementTopology* ref_topo_p = Elem::getRefinementTopology(elem_celltopo); // CHECK
      if (!ref_topo_p)
        throw std::runtime_error("genericRefine_createNewElements:: error, no refinement topology found");
      const Elem::RefinementTopology& ref_topo = *ref_topo_p;

      unsigned num_child = ref_topo.num_child();
      unsigned iChildStart = 0;
      //unsigned iChildEnd = num_child-1;
      // SPECIAL CASE ALERT
      if (fromTopoKey_in == topo_key_pyramid5 || fromTopoKey_in == topo_key_pyramid13)
        {
          num_child = getNumNewElemPerElem();
          if (toTopoKey_in == topo_key_tet4 || toTopoKey_in == topo_key_tet10)
            {
              iChildStart = 6;
              //iChildEnd = 9;
            }
        }
      else
        {
          //VERIFY_OP(num_child, == , getNumNewElemPerElem(), "genericRefine_createNewElements num_child problem");
        }

      // FIXME check if this is a wedge
      //bool homogeneous_child = ref_topo.homogeneous_child();
      //VERIFY_OP(homogeneous_child, ==, true, "genericRefine_createNewElements homogeneous_child");

      //RefTopoX& ref_topo_x = Elem::StdMeshObjTopologies::RefinementTopologyExtra< FromTopology > ::refinement_topology;
      RefTopoX_arr ref_topo_x = ref_topo_x_in;

      for (unsigned iChild = 0; iChild < num_child; iChild++)
        {
          int iChildRefTopo = iChild + iChildStart;

          VERIFY_OP_ON(iChild, <, elems.size(), "bad iChild");

          std::vector<stk::mesh::EntityId>& EN = elems[iChild];

          for (unsigned jNode = 0; jNode < (unsigned)ToTopology_node_count; jNode++)
            {
              unsigned childNodeIdx = ref_topo.child_node(iChildRefTopo)[jNode];

#ifndef NDEBUG
              unsigned childNodeIdxCheck = ref_topo_x[childNodeIdx].ordinal_of_node;
              VERIFY_OP(childNodeIdx, ==, childNodeIdxCheck, "childNodeIdxCheck");
#endif

              stk::mesh::EntityId inode=0;
              unsigned rank_of_subcell            = ref_topo_x[childNodeIdx].rank_of_subcell;
              unsigned ordinal_of_subcell         = ref_topo_x[childNodeIdx].ordinal_of_subcell;
              unsigned ordinal_of_node_on_subcell = ref_topo_x[childNodeIdx].ordinal_of_node_on_subcell;
              unsigned num_nodes_on_subcell       = ref_topo_x[childNodeIdx].num_nodes_on_subcell;

              bool usePerm = true;

              // only need permuation for quadratic elements
              if (num_nodes_on_subcell == 1)
                usePerm = false;

              const unsigned * perm_array = 0;
              if (usePerm)
                {
                  int perm_ord = getPermutation(eMesh,FromTopology_vertex_count, element, cell_topo, rank_of_subcell, ordinal_of_subcell);

                  if (perm_ord < 0)
                    throw std::logic_error("permutation < 0 ");
                  if (1 <= rank_of_subcell && rank_of_subcell <= 2)
                    {
                      perm_array = cell_topo.getCellTopologyData()->subcell[rank_of_subcell][ordinal_of_subcell].topology->permutation[perm_ord].node;
                    }
                }

              // FIXME these ranks are hardcoded because there is a basic assumption in the table generation code that
              //    assumes these ranks are used
              switch (rank_of_subcell)
                {
                case 0:
                  inode = VERT_N(ordinal_of_subcell);
                  break;

                case 1:
                  if (usePerm) // FIXME
                    if (num_nodes_on_subcell > 1)
                      inode = NN_Q_P(eMesh.edge_rank(), ordinal_of_subcell, ordinal_of_node_on_subcell, perm_array);
                    else
                      inode = NN_Q(eMesh.edge_rank(), ordinal_of_subcell, ordinal_of_node_on_subcell);
                  else
                    inode = EDGE_N_Q(ordinal_of_subcell, ordinal_of_node_on_subcell);

                  break;

                case 2:
                  if (cellDimension == 2)
                    {
                      VERIFY_OP(ordinal_of_subcell, == , 0, "createNewElements: ordinal_of_subcell");

                      if (usePerm)
                        {
                          if (num_nodes_on_subcell > 1)
                            {
                              inode = NN_Q_P(m_primaryEntityRank, ordinal_of_subcell, ordinal_of_node_on_subcell, perm_array);
                            }
                          else
                            {
                              inode = NN_Q(m_primaryEntityRank, ordinal_of_subcell, ordinal_of_node_on_subcell);
                            }
                        }
                      else
                        {
                          inode = NN_Q(m_primaryEntityRank, ordinal_of_subcell, ordinal_of_node_on_subcell);
                        }
                    }
                  else
                    {
                      if (usePerm)
                        {
                          if (num_nodes_on_subcell > 1)
                            {
                              inode = NN_Q_P(eMesh.face_rank(), ordinal_of_subcell, ordinal_of_node_on_subcell, perm_array);
                            }
                          else
                            {
                              inode = NN_Q(eMesh.face_rank(), ordinal_of_subcell, ordinal_of_node_on_subcell);
                            }
                        }
                      else
                        {
                          inode = NN_Q(eMesh.face_rank(), ordinal_of_subcell, ordinal_of_node_on_subcell);
                        }
                    }
                  break;
                case 3:
                  inode = NN_Q(m_primaryEntityRank, ordinal_of_subcell, ordinal_of_node_on_subcell);
                  break;
                default:
                  throw std::logic_error("UniformRefinerPattern logic error");
                }
              EN[jNode] = inode;
            }
        }

      bool use_declare_element_side = UniformRefinerPatternBase::USE_DECLARE_ELEMENT_SIDE &&  getPrimaryEntityRank() == eMesh.side_rank();

      std::vector<stk::mesh::Entity> nodes(ToTopology_node_count);

      for (unsigned iChild = 0; iChild < num_child; iChild++)
        {
          int iChildRefTopo = iChild + iChildStart;
          stk::mesh::Entity newElement = stk::mesh::Entity();
          if (!use_declare_element_side)
            {
              newElement = *element_pool;
            }

          for (int inode=0; inode < ToTopology_node_count; inode++)
            {
              VERIFY_OP_ON(inode, <, (int)elems[iChild].size(), "bad inode");
              stk::mesh::EntityId eid = elems[iChild][inode];
              if (!eid)
                {
                  throw std::logic_error("UniformRefinerPatternBase::genericRefine_createNewElements bad entity id = 0 ");
                }

              stk::mesh::Entity node = createOrGetNode(nodeRegistry, eMesh, eid);

              nodes[inode] = node;
            }

          create_side_element(eMesh, use_declare_element_side, nodes.data(), static_cast<unsigned>(ToTopology_node_count), newElement);

          if (m_primaryEntityRank == stk::topology::ELEMENT_RANK &&  proc_rank_field)
            {
              double *fdata = eMesh.field_data( *static_cast<const ScalarFieldType *>(proc_rank_field) , newElement );
              fdata[0] = double(eMesh.owner_rank(newElement));
            }

          change_entity_parts(eMesh, element, newElement);

          if (!isLinearElement)
            {
              prolongateFields(eMesh, element, newElement, ref_topo.child_node(iChildRefTopo),  &ref_topo_x[0], eMesh.get_coordinates_field() );
              prolongateFields(eMesh, element, newElement, ref_topo.child_node(iChildRefTopo),  &ref_topo_x[0]);
            }

          set_parent_child_relations(eMesh, element, newElement, *ft_element_pool, iChildRefTopo);

          std::vector<stk::mesh::Entity> elements(1,element);
          eMesh.prolongateElementFields( elements, newElement);

          ft_element_pool++;
          if (!use_declare_element_side)
            element_pool++;
        }
    }

    /// helpers for interpolating fields, coordinates
    /// ------------------------------------------------------------------------------------------------------------------------

    /// This version uses Intrepid2 for interpolation
    void UniformRefinerPatternBase::
    prolongateFields(percept::PerceptMesh& eMesh, stk::mesh::Entity element, stk::mesh::Entity newElement,  const unsigned *child_nodes,
                      RefTopoX_arr ref_topo_x, stk::mesh::FieldBase *field)
    {
#if STK_PERCEPT_LITE || !defined(STK_PERCEPT_USE_INTREPID)
      VERIFY_MSG("not available in PerceptMeshLite");
#else
      EXCEPTWATCH;

      unsigned *null_u = 0;

      unsigned toTopoKey = getToTypeKey();
      shards::CellTopology cell_topo(eMesh.get_cell_topology(element));

      // FIXME - need topo dimensions here
      int topoDim = getTopoDim(cell_topo);

      int fieldStride = 0;
      stk::mesh::EntityRank fr_type = static_cast<stk::mesh::EntityRank>(field->entity_rank());

      {
        unsigned nfr = field->restrictions().size();
        for (unsigned ifr = 0; ifr < nfr; ifr++)
          {
            const stk::mesh::FieldRestriction& fr = field->restrictions()[ifr];
            fieldStride = fr.num_scalars_per_entity() ;
            stk::mesh::Selector frselect = fr.selector();
          }
        {
          const stk::mesh::FieldBase::Restriction & r =
            stk::mesh::find_restriction(*field, fr_type, stk::mesh::MetaData::get(*field).universal_part());
          fieldStride = r.num_scalars_per_entity();
        }
      }
      // FIXME
      if (!fieldStride || fr_type != stk::topology::NODE_RANK)
        return;

      FieldFunction field_func("tmp", field, eMesh, topoDim, fieldStride);

      MDArray input_pts("input_pts",1, topoDim);
      MDArray input_param_coords("input_param_coords",1, topoDim);
      MDArray output_pts("input_param_coords",1, fieldStride);

      percept::MyPairIterRelation new_elem_nodes (eMesh, newElement, stk::topology::NODE_RANK);
      for (unsigned i_new_node = 0; i_new_node < new_elem_nodes.size(); i_new_node++)
        {
          unsigned childNodeIdx = child_nodes[i_new_node];
          double *param_coord = ref_topo_x[childNodeIdx].parametric_coordinates;

          for (int ip=0; ip < topoDim; ip++)
            {
              input_param_coords(0, ip) = param_coord[ip];
            }
          double time_val=0.0;

          if (toTopoKey == topo_key_wedge15 || toTopoKey == topo_key_quad8 || toTopoKey == topo_key_shellquad8 
              || toTopoKey == topo_key_line3  || toTopoKey == topo_key_shellline3
              || toTopoKey == topo_key_hex20
              || toTopoKey == topo_key_pyramid13 || toTopoKey == topo_key_pyramid5
              || toTopoKey == topo_key_tet10)
            {
              prolongateIntrepid2(eMesh, field, cell_topo, output_pts, element, input_param_coords, time_val);
            }
          else
            {
              field_func(input_pts, output_pts, element, input_param_coords, time_val);
            }

          stk::mesh::Entity new_node = new_elem_nodes[i_new_node].entity();

          {
            double *f_data_new = eMesh.field_data(field, new_node, null_u);
            for (int ifd=0; ifd < fieldStride; ifd++)
              {
                f_data_new[ifd] = output_pts(0, ifd);
              }
          }
        }
#endif
    }

    /// do interpolation for all fields
    /// This version uses Intrepid2
    void UniformRefinerPatternBase::
    prolongateFields(percept::PerceptMesh& eMesh, stk::mesh::Entity element, stk::mesh::Entity newElement, const unsigned *child_nodes,
                      RefTopoX_arr ref_topo_x)
    {
      const stk::mesh::FieldVector & fields = eMesh.get_fem_meta_data()->get_fields();
      unsigned nfields = fields.size();
      //std::cout << "P[" << p_rank << "] info>    Number of fields = " << fields.size() << std::endl;
      for (unsigned ifld = 0; ifld < nfields; ifld++)
        {
          stk::mesh::FieldBase *field = fields[ifld];
          //std::cout << "P[" << eMesh.get_rank() << "] field = " << field->name() << std::endl;

          const stk::mesh::DataTraits & data_traits = field->data_traits();
          if (data_traits.is_floating_point)
            prolongateFields(eMesh, element, newElement, child_nodes, ref_topo_x, field);
        }
    }

    void UniformRefinerPatternBase::
    prolongateIntrepid2(percept::PerceptMesh& eMesh, stk::mesh::FieldBase* field, shards::CellTopology& cell_topo,
                        MDArray& output_pts, stk::mesh::Entity element, MDArray& input_param_coords, double time_val)
    {
#if STK_PERCEPT_LITE
      VERIFY_MSG("not available in PerceptMeshLite");
#else
      int fieldStride = output_pts.extent_int(1);
      unsigned *null_u = 0;

      percept::MyPairIterRelation elem_nodes (eMesh, element, stk::topology::NODE_RANK);
      MDArray basis_val("basis_val",elem_nodes.size(), 1);
      //std::cout << "tmp fieldStride= " << fieldStride << " elem_nodes.size()= " << elem_nodes.size() << std::endl;

      /// special for pyramid
      if (cell_topo.getKey() == topo_key_pyramid13) // || cell_topo.getKey() == topo_key_pyramid5)
        {
          // input is [-1,1]x[-1,1]x[0,1]
          double r = input_param_coords(0,0);
          double s = input_param_coords(0,1);
          double t = input_param_coords(0,2);
          // transform to [0,1]x[0,1]x[0,0.49999]
          r = (1+r)/2;
          s = (1+s)/2;
          t = t*0.49999;
          double bases[] = {-(((-1 + r + t)*(-1 + s + t))/(-1 + 2*t)),
                            ((r - t)*(-1 + s + t))/(-1 + 2*t),
                            -(((r - t)*(s - t))/(-1 + 2*t)),
                            ((s - t)*(-1 + r + t))/(-1 + 2*t),
                            2*t};
          for (unsigned ii=0; ii < 5; ii++)
            {
              basis_val(ii,0) = bases[ii];
            }
        }
      else
        {
          BasisTable::BasisTypeRCP basis = BasisTable::getInstance()->getBasis(cell_topo);
          basis->getValues(basis_val, input_param_coords, Intrepid2::OPERATOR_VALUE);
        }
      if (0)
        std::cout << "\n tmp input_param_coords= "
                  << input_param_coords(0,0) << " "
                  << input_param_coords(0,1) << " "
                  << input_param_coords(0,2) << " " << std::endl;

      for (int i_stride=0; i_stride < fieldStride; i_stride++)
        {
          output_pts(0, i_stride) = 0.0;
        }
      for (unsigned i_node = 0; i_node < elem_nodes.size(); i_node++)
        {
          //std::cout << "tmp basis_val[" << i_node <<"]= " << basis_val(i_node,0) << std::endl;
          stk::mesh::Entity node = elem_nodes[i_node].entity();
          double *f_data = eMesh.field_data(field, node, null_u);
          for (int i_stride=0; i_stride < fieldStride; i_stride++)
            {
              output_pts(0, i_stride) += f_data[i_stride]*basis_val(i_node, 0);
            }
        }
#endif
    }



    stk::mesh::Entity UniformRefinerPatternBase::
    createOrGetNode(NodeRegistry& nodeRegistry, PerceptMesh& eMesh, stk::mesh::EntityId eid)
    {
#if STK_ADAPT_NODEREGISTRY_USE_ENTITY_REPO
      stk::mesh::Entity node_p = nodeRegistry.get_entity_node_Ib(*eMesh.get_bulk_data(), stk::topology::NODE_RANK, eid);
      if (node_p)
        return *node_p;
      else
        return eMesh.createOrGetNode(eid);
#else
      return eMesh.createOrGetNode(eid);
#endif
    }

#define DEBUG_CHANGE_ENTITY_PARTS 0

    void UniformRefinerPatternBase::
    change_entity_parts(percept::PerceptMesh& eMesh, stk::mesh::Entity old_owning_elem, stk::mesh::Entity newElement)
    {
      static std::vector<stk::mesh::Part*> add_parts(1);
      static std::vector<stk::mesh::Part*> remove_parts;

      bool found = false;
      if ( m_fromParts.size() != m_toParts.size())
        {
          std::cout << "this= " << eMesh.demangle(typeid(*this).name())
                    << " m_fromParts.size= " << m_fromParts.size()
                    << " m_toParts.size= " << m_toParts.size() << std::endl;
          std::cout << "tmp change_entity_parts printParts this=\n" ;
          printParts(this, true);
          throw std::runtime_error("bad from/to parts");
        }
      for (unsigned i_part = 0; i_part < m_fromParts.size(); i_part++)
        {
          if (eMesh.bucket(old_owning_elem).member(*m_fromParts[i_part]))
            {
              add_parts[0] = m_toParts[i_part];
              if (DEBUG_CHANGE_ENTITY_PARTS)
                {
                  stk::mesh::PartVector pv = eMesh.bucket(newElement).supersets();
                  std::string str, strOld;
                  for (unsigned ii=0; ii < pv.size(); ++ii)
                    str += pv[ii]->name()+" ";
                  stk::mesh::PartVector pvOld = eMesh.bucket(old_owning_elem).supersets();
                  for (unsigned ii=0; ii < pvOld.size(); ++ii)
                    strOld += pvOld[ii]->name()+" ";
                  if (1)
                    std::cout << "tmp my class= " << eMesh.demangle(typeid(*this).name()) << "\n changing newElement " << eMesh.identifier(newElement)
                              << " rank= " << eMesh.entity_rank(newElement)
                              << " from part= " << m_fromParts[i_part]->name() << " T= " << m_fromParts[i_part]->topology()
                              << " to part=add_parts= " << m_toParts[i_part]->name() << " T= " << m_toParts[i_part]->topology()
                              << " for old elem= " << eMesh.identifier(old_owning_elem)
                              << " rank= " << eMesh.entity_rank(old_owning_elem)
                              << " old supersets= " << strOld
                              << " nsupersets_new size= " << eMesh.bucket(newElement).supersets().size()
                              << " supersets= " << str
                              << std::endl;
                }
              VERIFY_OP_ON(m_toParts[i_part], !=, 0, "bad to part");
              const CellTopologyData *ctopo = eMesh.get_cell_topology(*m_toParts[i_part]);
              if (getToTopology() != ctopo)
                {
                  if (ctopo == 0 && m_toParts[i_part]->subsets().size() == 0)
                    {
                      if (getPrimaryEntityRank() != m_toParts[i_part]->primary_entity_rank())
                        continue;
                    }
                  else
                    continue;
                }
              if (!ctopo)
                {
                  stk::mesh::Part& root_part = eMesh.get_fem_meta_data()->get_topology_root_part(stk::mesh::get_topology(getToTopology(), eMesh.get_fem_meta_data()->spatial_dimension()));
                  add_parts.push_back( &root_part);
                }
              eMesh.get_bulk_data()->change_entity_parts( newElement, add_parts, remove_parts );
              if (add_parts.size() != 1)
                add_parts.resize(1);

              found = true;
            }
        }
      if (!found)
        {
          std::cout << "stacktrace= " << eMesh.demangled_stacktrace() << std::endl;
          throw std::runtime_error("URP::change_entity_parts couldn't find part");
        }
    }

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

    int UniformRefinerPatternBase::
    getPermutation(PerceptMesh& eMesh, int num_verts, stk::mesh::Entity element, shards::CellTopology& cell_topo, unsigned rank_of_subcell, unsigned ordinal_of_subcell)
    {
      if (rank_of_subcell == 0 || rank_of_subcell == 3) return 0;

      //! We choose to define the "global baseline" as an imaginary face that has its nodes sorted on their identifiers.
      //! The main part of this is to find the minimum node index.
      //! Once sorted, the node id's go into the vector_sdcell_global_baseline array.

      static std::vector<unsigned> vector_sdcell_global_baseline(4);
      static std::vector<unsigned> subCell_from_element(4);

      const percept::MyPairIterRelation elem_nodes (eMesh, element, stk::topology::NODE_RANK);

      const unsigned * inodes = cell_topo.getCellTopologyData()->subcell[rank_of_subcell][ordinal_of_subcell].node;
      int num_subcell_verts = cell_topo.getCellTopologyData()->subcell[rank_of_subcell][ordinal_of_subcell].topology->vertex_count;

      // tmp
      //vector_sdcell_global_baseline.resize(num_subcell_verts);
      //subCell_from_element.resize(num_subcell_verts);
      // tmp end

      stk::mesh::EntityId minNodeId = 0;
      for (int iv = 0; iv < num_subcell_verts; iv++)
        {
          stk::mesh::EntityId nid = eMesh.identifier(elem_nodes[inodes[iv]].entity());
          if (iv == 0)
            minNodeId = nid;
          else
            minNodeId = std::min(minNodeId, nid);
          subCell_from_element[iv] = nid;
        }

      int perm = -1;

      /// for tri or quad faces we search for the min node, then look at its two neighbors along edges
      ///   - if the first edge is the next node in line (in terms of its id) we use the ordering as 0,1,2,3
      ///   else it is flipped and we reverse the ordering
      ///   - once the ordering is determined, the actual permutation can be deduced from shards
      ///

      {
        // quad or tri
        if (0 && num_subcell_verts==3)
          {
            std::cout << "tmp b4 element 1= " << element << " cell_topo= " << cell_topo.getName()
                      << " rank_of_subcell= " << rank_of_subcell << std::endl;
            std::cout << "tmp b4 vector_sdcell_global_baseline= " << vector_sdcell_global_baseline << std::endl;
            std::cout << "tmp b4 subCell_from_element = " << subCell_from_element << std::endl;
          }

        //! extract the minimal node index
        //set_sdcell_global_baseline_iter = set_sdcell_global_baseline.begin();
        stk::mesh::EntityId i0 = minNodeId;

        //! find the rotation to get to the minimal node
        int j0 = -1;
        for (int iv = 0; iv < num_subcell_verts; iv++)
          {
            if (i0 == subCell_from_element[iv])
              {
                j0 = iv;
                break;
              }
          }

        if (j0 < 0) throw std::logic_error("j0 < 0 ");

        int j1 = (j0 + 1) % num_subcell_verts;
        int j2 = (j0 + (num_subcell_verts-1)) % num_subcell_verts;  // adds 3 for quads, or 2 for tris to pickup the neigh node

        //! see if we need to reverse the order to make it match up; save the newly oriented nodes in vector_sdcell_global_baseline
        if (subCell_from_element[j1] < subCell_from_element[j2])
          {
            for (int iv = 0; iv < num_subcell_verts; iv++)
              {
                vector_sdcell_global_baseline[iv] = subCell_from_element[(j0 + iv) % num_subcell_verts];
              }
          }
        else
          {
            for (int iv = 0; iv < num_subcell_verts; iv++)
              {
                vector_sdcell_global_baseline[(num_subcell_verts - iv) % num_subcell_verts] =
                  subCell_from_element[(j0 + iv) % num_subcell_verts];
              }
          }

        //! now we have a set of nodes in the right order, use Shards to get the actual permutation
        perm = shards::findPermutation(cell_topo.getCellTopologyData()->subcell[rank_of_subcell][ordinal_of_subcell].topology,
                                       &vector_sdcell_global_baseline[0], &subCell_from_element[0]);

        //std::cout << "tmp perm = " << perm << std::endl;

        if ( perm < 0)
          {
            std::cout << "tmp aft element 1= " << element << " cell_topo= " << cell_topo.getName() << " rank_of_subcell= " << rank_of_subcell << std::endl;
            std::cout << "tmp aft vector_sdcell_global_baseline= " << vector_sdcell_global_baseline << std::endl;
            std::cout << "tmp aft subCell_from_element = " << subCell_from_element << std::endl;
            throw std::logic_error("getPermutation: perm < 0");
          }

        if (0 && num_subcell_verts==3)
          {
            const unsigned *perm_array = cell_topo.getCellTopologyData()->subcell[rank_of_subcell][ordinal_of_subcell].topology->permutation[perm].node;
            for (int iv = 0; iv < num_subcell_verts; iv++)
              {
                std::cout << "tmp perm_array[" << iv << "]=  " << perm_array[iv] << std::endl;
              }
          }

      }

      if (perm < 0)
        {
          std::cout << "tmp element 1= " << element << " cell_topo= " << cell_topo.getName() << " rank_of_subcell= " << rank_of_subcell << std::endl;
          std::cout << "tmp vector_sdcell_global_baseline= " << vector_sdcell_global_baseline << std::endl;
          std::cout << "tmp subCell_from_element = " << subCell_from_element << std::endl;
          throw std::logic_error("getPermutation 2: perm < 0");
        }

      return perm;
    }


    /// sets the needed number of nodes on each sub-entity to 1 - this is just a helper - in general, edges and faces have 1 new node
    /// for linear elements, and multiple new nodes in the case of quadratic elements
    void UniformRefinerPatternBase::
    setToOne(std::vector<NeededEntityType>& needed_entities)
    {
      for (unsigned i = 0; i < needed_entities.size(); i++)
        {
          needed_entities[i].second = 1u;
        }
    }
    double * UniformRefinerPatternBase::
    midPoint(const double *p1, const double *p2, int spatialDim, double *x)
    {
      x[0] = 0.5*(p1[0]+p2[0]);
      x[1] = 0.5*(p1[1]+p2[1]);
      if (spatialDim == 3)
        x[2] = 0.5*(p1[2]+p2[2]);
      return x;
    }

    double * UniformRefinerPatternBase::
    getCentroid( double* pts[], int len, int spatialDim, double *x)
    {
      double dlen = double(len);
      for (int jsp = 0; jsp < spatialDim; jsp++)
        {
          x[jsp] = 0.0;
        }
      for (int ipt = 0; ipt < len; ipt++)
        {
          for (int jsp = 0; jsp < spatialDim; jsp++)
            {
              x[jsp] += pts[ipt][jsp] / dlen;
            }
        }
      return x;
    }

    //static
    int UniformRefinerPatternBase::
    getTopoDim(shards::CellTopology& cell_topo)
    {
      int topoDim = cell_topo.getDimension();
      unsigned cell_topo_key = cell_topo.getKey();

      switch (cell_topo_key)
        {
        case base_s_shell_line_2_key:
        case base_s_shell_line_3_key:
        case base_s_beam_2_key:
        case base_s_beam_3_key:
          topoDim = 1;
          break;

        case base_s_shell_tri_3_key:
        case base_s_shell_tri_6_key:
        case base_s_shell_quad_4_key:
        case base_s_shell_quad_9_key:
        case base_s_shell_quad_8_key:
          topoDim = 2;
          break;
        }

      return topoDim;
    }

    /// optionally overridden (must be overridden if sidesets are to work properly) to provide info on which sub pattern
    /// should be used to refine side sets (and edge sets)
    void UniformRefinerPatternBase::setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
    {
      /// default is only this pattern
      bp.resize(1); // = std::vector<UniformRefinerPatternBase *>(1u, 0);
      bp[0] = this;
    }

    void UniformRefinerPatternBase::setSubPatternsForSetNeededParts( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
    {
      setSubPatterns(bp, eMesh);
    }

    /// for i/o to work properly, supply string replacements such as for hex-->tet breaking, you would supply "quad"-->"tri" etc. string maps
    StringStringMap UniformRefinerPatternBase::fixSurfaceAndEdgeSetNamesMap()
    {
      // provide a null implementation
      StringStringMap map;
      return map;
    }

    size_t UniformRefinerPatternBase::estimateNumberOfNewElements(percept::PerceptMesh& eMesh, stk::mesh::EntityRank rank, NodeRegistry& nodeRegistry, size_t num_elem_not_ghost)
    {
      EXCEPTWATCH;

      vector<NeededEntityType> needed_entity_ranks;
      this->fillNeededEntities(needed_entity_ranks);

      size_t num_elem_marked = 0;

      const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( rank );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          const unsigned num_elements_in_bucket = bucket.size();
          const CellTopologyData * const cell_topo_data = eMesh.get_cell_topology(bucket);
          if (cell_topo_data == 0 || 
              cell_topo_data->name != getFromTopoPartName())
            {
              continue;
            }
          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
            {
              stk::mesh::Entity element = bucket[iElement];

              if (eMesh.isGhostElement(element))
                continue;
              if (eMesh.hasFamilyTree(element) && eMesh.isParentElement(element, true))
                continue;

              bool foundMarkedEdge = false;

              //const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);
              for (unsigned ineed_ent=0; ineed_ent < needed_entity_ranks.size(); ineed_ent++)
                {
                  unsigned numSubDimNeededEntities = 0;

                  if (this->edgeMarkIsEnough() && needed_entity_ranks[ineed_ent].first != eMesh.edge_rank())
                    continue;

                  // special case of face in 3d or edge in 2d
                  if (needed_entity_ranks[ineed_ent].first == eMesh.entity_rank(element))
                    {
                      numSubDimNeededEntities = 1;
                    }
                  else if (needed_entity_ranks[ineed_ent].first == eMesh.edge_rank())
                    {
                      numSubDimNeededEntities = cell_topo_data->edge_count;
                    }
                  else if (needed_entity_ranks[ineed_ent].first == eMesh.face_rank())
                    {
                      numSubDimNeededEntities = cell_topo_data->side_count;
                    }
                  else if (needed_entity_ranks[ineed_ent].first == stk::topology::ELEMENT_RANK)
                    {
                      numSubDimNeededEntities = 1;
                    }

                  for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                    {
                      NodeIdsOnSubDimEntityType* nodeIds_onSE_ptr = nodeRegistry.getNewNodesOnSubDimEntity(element, needed_entity_ranks[ineed_ent].first, iSubDimOrd);
                      if (nodeIds_onSE_ptr == 0)
                        {
                          continue;
                        }
                      NodeIdsOnSubDimEntityType& nodeIds_onSE = *nodeIds_onSE_ptr;

                      bool doMark = true;
                      if (needed_entity_ranks[ineed_ent].third.size())
                        {
                          VERIFY_OP_ON(needed_entity_ranks[ineed_ent].third.size(), ==, numSubDimNeededEntities, "bad size");
                          if (!needed_entity_ranks[ineed_ent].third[iSubDimOrd])
                            doMark = false;
                        }

                      if (doMark && nodeIds_onSE.size() != 0)
                        {
                          foundMarkedEdge = true;
                          break;
                        }
                    }
                }
              if (foundMarkedEdge)
                ++num_elem_marked;
            }
        }

      VERIFY_OP_ON(this->getNumNewElemPerElem(), >, 0, "bad getNumNewElemPerElem");
      unsigned num_elem_needed = num_elem_marked * this->getNumNewElemPerElem();
      return num_elem_needed;
    }

    void UniformRefinerPatternBase::mergeOrAddParts(UniformRefinerPatternBase *bp_from, UniformRefinerPatternBase *bp_to, bool merge)
    {
      stk::mesh::PartVector& fromParts = bp_from->getFromParts();
      stk::mesh::PartVector& toParts = bp_from->getToParts();
      if (fromParts.size() != toParts.size()) {
        printParts(bp_from);
        printParts(bp_to);
      }
      VERIFY_OP_ON(fromParts.size(), ==, toParts.size(), "from parts/to parts size mismatch");
      for (unsigned ii=0; ii < fromParts.size(); ii++)
        {
          if (!merge || std::find(bp_to->getFromParts().begin(), bp_to->getFromParts().end(), fromParts[ii]) == bp_to->getFromParts().end())
            {
              bp_to->getFromParts().push_back(fromParts[ii]);
              bp_to->getToParts().push_back(toParts[ii]);
            }
        }
    }


    void UniformRefinerPatternBase::printParts(UniformRefinerPatternBase *bp, bool printAllParts)
    {
      std::cout << "\n===================================\n=====================================\ntmp printParts " << PerceptMesh::demangle(typeid(*bp).name()) << std::endl;
      if (bp->getFromParts().size() != bp->getToParts().size())
        {
          for (unsigned ii=0; ii < bp->getFromParts().size(); ii++)
            {
              std::cout << "tmp printParts ii, fromParts= " << ii << " " << (bp->getFromParts()[ii] ? bp->getFromParts()[ii]->name() : "Null From Part!!!") << std::endl;
            }
          for (unsigned ii=0; ii < bp->getToParts().size(); ii++)
            {
              std::cout << "tmp printParts ii, toParts= " << ii << " " << (bp->getToParts()[ii] ? bp->getToParts()[ii]->name() : "Null To Part!!!") << std::endl;
            }
        }
      else
        {
          for (unsigned ii=0; ii < bp->getFromParts().size(); ii++)
            {
              std::cout << "tmp printParts from/to parts[" << ii << "] = "
                        << (bp->getFromParts()[ii] ? bp->getFromParts()[ii]->name() : "Null From Part!!!")
                        << " ==> " << (bp->getToParts()[ii] ? bp->getToParts()[ii]->name() : "Null To Part!!!")
                        << std::endl;
            }
        }

      if (printAllParts && bp->getFromParts().size())
        {
          std::cout << std::endl;
          const stk::mesh::MetaData& metaData = stk::mesh::MetaData::get(*bp->getFromParts()[0]);
          const std::vector< stk::mesh::Part * > & parts = metaData.get_parts();
          for (unsigned ii=0; ii < parts.size(); ++ii)
            {
              std::cout << "tmp part[" << ii << "] = " << parts[ii]->name() << std::endl;
            }
        }
      std::cout << "\n\n";
    }


    void UniformRefinerPatternBase::create_side_element(PerceptMesh& eMesh, bool use_declare_element_side, stk::mesh::Entity *nodes, unsigned nodes_size, stk::mesh::Entity& newElement)
    {
      stk::topology topo = stk::mesh::get_topology( shards::CellTopology (getToTopology()));

      stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(0);
      std::vector<stk::mesh::Entity> new_nodes(nodes_size);

      if (use_declare_element_side)
        {
          stk::mesh::Entity element_found = stk::mesh::Entity();
          unsigned side_ord_found = 0;
          perm = eMesh.find_side_element_permutation_and_appropriate_element(nodes, nodes_size, topo, new_nodes.data(), element_found, side_ord_found);
          if (perm == stk::mesh::INVALID_PERMUTATION)
            {
              std::cout << eMesh.rank() << " perm= " << perm << std::endl;
              perm = eMesh.find_side_element_permutation_and_appropriate_element(nodes, nodes_size, topo, new_nodes.data(), element_found, side_ord_found, true);
            }
          VERIFY_OP_ON(perm, !=, stk::mesh::INVALID_PERMUTATION, "bad permutation");

          VERIFY_OP_ON(eMesh.is_valid(element_found), ==, true, "newElement");

          if (1)
            {
              newElement = eMesh.get_bulk_data()->declare_element_side(element_found, side_ord_found, stk::mesh::ConstPartVector{});
              VERIFY_OP_ON(eMesh.is_valid(newElement), ==, true, "bad newElement");
            }
        }
      else
        {
          for (unsigned inode=0; inode < nodes_size; ++inode)
            {
              eMesh.get_bulk_data()->declare_relation(newElement, nodes[inode], inode);
            }
        }

    }

  }


