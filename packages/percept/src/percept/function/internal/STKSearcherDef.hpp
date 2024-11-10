// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_STKSearcherDef_hpp
#define percept_STKSearcherDef_hpp

#include <percept/Percept.hpp>
#include <percept/util/Loops.hpp>
#include <percept/FieldTypes.hpp>

#include <Intrepid2_FunctionSpaceTools.hpp>
#include "STKSearcher.hpp"

  namespace percept
  {

    STKSearcher::STKSearcher(stk::mesh::BulkData *bulk) : m_bulk(bulk), m_boxes()
    {
    }

    STKSearcher::~STKSearcher() {}

    void
    STKSearcher::setupSearch()
    {
      const stk::mesh::MetaData& metaData = m_bulk->mesh_meta_data();
      stk::mesh::BulkData& bulkData = *m_bulk;
      CoordinatesFieldType *coords_field = metaData.get_field<double>(stk::topology::NODE_RANK, "coordinates");
      PerceptMesh meshUtil(&metaData, &bulkData);

      BBB buildBoundingBoxes(m_boxes, coords_field);
      elementOpLoop(*meshUtil.get_bulk_data(), buildBoundingBoxes, coords_field);
    }

    void
    STKSearcher::tearDownSearch()
    {
      m_boxes.clear();
    }

    /**
     *  Dimensions of input_phy_points = ([P]=1, [D])
     *  Dimensions of found_parametric_coordinates = ([P]=1, [D])
     */

    const stk::mesh::Entity
    STKSearcher::findElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates,
                                         unsigned& found_it, const stk::mesh::Entity hint_element )
    {
      //return 0;
      const stk::mesh::MetaData& metaData = m_bulk->mesh_meta_data();
      stk::mesh::BulkData& bulkData = *m_bulk;

      //CoordinatesFieldType *coords_field = metaData.get_field<CoordinatesFieldType::value_type >(stk::topology::NODE_RANK, "coordinates");

      PerceptMesh meshUtil(&metaData, &bulkData);

      double pts[3] = {};
      for (unsigned iDim = 0; iDim < metaData.spatial_dimension(); iDim++)
        {
          pts[iDim] = input_phy_points(0, iDim);
        }
      BPoint pointBoundingBox;
      pointBoundingBox.second.set_id(123);  // FIXME for multiple points
      pointBoundingBox.first[0] = pts[0];
      pointBoundingBox.first[1] = pts[1];
      pointBoundingBox.first[2] = pts[2];
      std::vector<BPoint> points(1, pointBoundingBox);

      if (EXTRA_PRINT)
        {
          bool box_p = stk::search::intersects(m_boxes[0].first, pointBoundingBox.first);
          bool p_box = stk::search::intersects(pointBoundingBox.first, m_boxes[0].first);

          std::cout << "STKSearcher::findElement: m_boxes[0]=  (" << m_boxes[0].first << "," << m_boxes[0].second << ")"  << std::endl;
          std::cout << "STKSearcher::findElement: pointBoundingBox=  (" << pointBoundingBox.first << "," << pointBoundingBox.second << ")" << std::endl;
          std::cout << "STKSearcher::findElement: box_p=  " << box_p << std::endl;
          std::cout << "STKSearcher::findElement: p_box=  " << p_box << std::endl;
        }

      if (EXTRA_PRINT) std::cout << "STKSearcher::findElement: nboxes=  " << m_boxes.size()  << std::endl;

      IdentProcRelation relation;
      stk::search::coarse_search(points, m_boxes, stk::search::KDTREE, bulkData.parallel(), relation);

      if (EXTRA_PRINT) std::cout << "STKSearcher::findElement: found  " << relation.size() << " containing bboxes"  << std::endl;

      if (relation.size())
        {
          IsInElement isIn(input_phy_points, found_parametric_coordinates);

          for (unsigned i = 0; i < relation.size(); i++)
            {
              if (EXTRA_PRINT)
                std::cout << "relation[ " << i << "]= {" << relation[i].first << "} --> { " << relation[i].second << "}" << std::endl;
              stk::mesh::Entity element = bulkData.get_entity(stk::topology::ELEMENT_RANK, relation[i].second.id());
              //bool loop_break = ... intentionally ignoring return value
              isIn(element, bulkData);
              if (EXTRA_PRINT) std::cout << "STKSearcher::findElement: found it= " << isIn.m_found_it << std::endl;
              if (isIn.m_found_it)
                {
                  found_it = 1;
                  return isIn.m_foundElement;
                }
              else
                {
                  found_it = 0;
                  return stk::mesh::Entity();
                }
            }
        }

      return stk::mesh::Entity();
    }

  }

#endif
