// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <adapt/ElementRefinePredicate.hpp>

  namespace percept {

    /// Return DO_REFINE, DO_UNREFINE, DO_NOTHING
    int ElementRefinePredicate::operator()(const stk::mesh::Entity entity)
    {
      int mark = 0;

      if (getMarkNone())
        return DO_NOTHING;

      stk::mesh::EntityRank rank = m_eMesh.entity_rank(entity);
      if (rank != m_eMesh.element_rank())
        return DO_NOTHING;

      bool isParent = m_eMesh.hasFamilyTree(entity) && m_eMesh.isParentElement(entity, false);

      RefineFieldType::value_type *fdata = 0;
      if (m_field && m_field->entity_rank() == m_eMesh.entity_rank(entity))
        {
          //fdata = m_eMesh.field_data( *static_cast<const TransitionElementType *>(m_field) , entity );
          fdata = stk::mesh::field_data( *dynamic_cast<const RefineFieldType *>(m_field) , entity );
        }
      bool selected = (m_eb_selector==0 || (*m_eb_selector)(m_eMesh.bucket(entity)));

      bool ref_field_criterion = (fdata  && fdata[0] == 1);

      bool unref_field_criterion = (fdata && fdata[0] < 0);

      if (getRefineStage() == -10)
        {
          if (isParent)
            {
              ref_field_criterion = (fdata  && fdata[0] == 2);
              unref_field_criterion = false;
            }
          else
            {
              if (unref_field_criterion)
                {
                  // can only unrefine elements with parents
                  if (!m_eMesh.hasFamilyTree(entity))
                    unref_field_criterion = false;
                  else
                    {
                      stk::mesh::Entity parent = m_eMesh.getParent(entity, true);
                      if (!m_eMesh.is_valid(parent))
                        {
                          unref_field_criterion = false;
                        }
                    }
                }
            }
        }
      else
        {
          if (isParent) return DO_NOTHING;
          if (unref_field_criterion)
            {
              // can only unrefine elements with parents
              if (!m_eMesh.hasFamilyTree(entity))
                unref_field_criterion = false;
              else
                {
                  stk::mesh::Entity parent = m_eMesh.getParent(entity, true);
                  if (!m_eMesh.is_valid(parent))
                    {
                      unref_field_criterion = false;
                    }
                }
            }
        }

      if (selected && ref_field_criterion)
        {
          mark |= DO_REFINE;
        }
      if (selected && unref_field_criterion)
        {
          mark |= DO_UNREFINE;
        }

      return mark;
    }


  }

