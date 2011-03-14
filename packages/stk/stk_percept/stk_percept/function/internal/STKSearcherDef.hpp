#ifndef stk_percept_STKSearcherDef_hpp
#define stk_percept_STKSearcherDef_hpp

#include <stk_percept/Percept.hpp>

#include <Intrepid_FunctionSpaceTools.hpp>
#include "STKSearcher.hpp"

namespace stk
{
  namespace percept
  {

    template<int SpatialDim>
    STKSearcher<SpatialDim>::STKSearcher(FieldFunction *ff) : m_fieldFunction(ff), m_boxes()
    {
    }

    template<int SpatialDim>
    void
    STKSearcher<SpatialDim>::setupSearch()
    {
      mesh::MetaData& metaData = MetaData::get(*(m_fieldFunction->getField()));
      mesh::BulkData& bulkData = *m_fieldFunction->getBulkData();
      VectorFieldType *coords_field = metaData.get_field<VectorFieldType >("coordinates");
      PerceptMesh meshUtil(&metaData, &bulkData);

      BBB buildBoundingBoxes(m_boxes, coords_field);
      meshUtil.elementOpLoop(buildBoundingBoxes, coords_field);
    }

    template<int SpatialDim>
    void
    STKSearcher<SpatialDim>::tearDownSearch()
    {
      m_boxes.clear();
    }

    /**
     *  Dimensions of input_phy_points = ([P]=1, [D])
     *  Dimensions of found_parametric_coordinates = ([P]=1, [D])
     */

    template<int SpatialDim>
    const stk::mesh::Entity *
    STKSearcher<SpatialDim>::findElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates,
                                         unsigned& found_it, const mesh::Entity *hint_element )
    {
      //return 0;
      mesh::MetaData& metaData = MetaData::get(*(m_fieldFunction->getField()));
      mesh::BulkData& bulkData = *m_fieldFunction->getBulkData();

      //VectorFieldType *coords_field = metaData.get_field<VectorFieldType >("coordinates");

      PerceptMesh meshUtil(&metaData, &bulkData);

      double pts[SpatialDim];
      for (unsigned iDim = 0; iDim < SpatialDim; iDim++)
        {
          pts[iDim] = input_phy_points(0, iDim);
        }
      BPoint pointBoundingBox;
      pointBoundingBox.key.ident = 123;  // FIXME for multiple points
      pointBoundingBox.set_center(pts);
      std::vector<BPoint> points(1, pointBoundingBox);

      stk::search::FactoryOrder order;
      order.m_communicator = bulkData.parallel();
      order.m_algorithm = stk::search::FactoryOrder::BIHTREE;

      if (0 || EXTRA_PRINT)
        {
          bool box_p = m_boxes[0].intersect(pointBoundingBox);
          bool p_box = pointBoundingBox.intersect(m_boxes[0]);

          std::cout << "STKSearcher::findElement: m_boxes[0]=  " << m_boxes[0] << std::endl;
          std::cout << "STKSearcher::findElement: pointBoundingBox=  " << pointBoundingBox << std::endl;
          std::cout << "STKSearcher::findElement: box_p=  " << box_p << std::endl;
          std::cout << "STKSearcher::findElement: p_box=  " << p_box << std::endl;
        }

      if (0 || EXTRA_PRINT) std::cout << "STKSearcher::findElement: nboxes=  " << m_boxes.size()  << std::endl;

      IdentProcRelation relation;
      stk::search::coarse_search(relation,  m_boxes, points, order);
      //stk::search::coarse_search(relation,   points, m_boxes, order);

      if (0 || EXTRA_PRINT) std::cout << "STKSearcher::findElement: found  " << relation.size() << " containing bboxes"  << std::endl;

      if (relation.size())
        {
          IsInElement isIn(input_phy_points, found_parametric_coordinates);

          for (unsigned i = 0; i < relation.size(); i++)
            {
              if (0 || EXTRA_PRINT)
                std::cout << "relation[ " << i << "]= {" << relation[i].first << "} --> { " << relation[i].second << "}" << std::endl;
              mesh::Entity *element = bulkData.get_entity(mesh::Element, relation[i].second.ident);
              //bool loop_break = ... intentionally ignoring return value
              isIn(*element, m_fieldFunction->getField(), bulkData);
              if (0 || EXTRA_PRINT) std::cout << "STKSearcher::findElement: found it= " << isIn.m_found_it << std::endl;
              if (isIn.m_found_it)
                {
                  found_it = 1;
                  return isIn.m_foundElement;
                }
              else
                {
                  found_it = 0;
                  return 0;
                }
            }
        }

      return 0;
    }

  }
}
#endif
