// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>

#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_FunctionSpaceTools.hpp>

#include <percept/FieldTypes.hpp>

#include <percept/function/FieldFunction.hpp>
#include <percept/function/internal/SimpleSearcher.hpp>
#include <percept/util/Loops.hpp>

  namespace percept
  {

    SimpleSearcher::SimpleSearcher(stk::mesh::BulkData *bulk) : m_bulk(bulk) {}

    /**
     *  Dimensions of input_phy_points = ([P]=1, [D])
     *  Dimensions of found_parametric_coordinates = ([P]=1, [D])
     */

    const stk::mesh::Entity SimpleSearcher::findElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates,
                                                         unsigned& found_it, const stk::mesh::Entity hint_element )
    {
      VERIFY_OP(input_phy_points.rank(), ==, found_parametric_coordinates.rank(), "SimpleSearcher::findElement bad dims");
      //VERIFY_OP(input_phy_points.rank(), ==, 2, "SimpleSearcher::findElement bad rank");

      const stk::mesh::MetaData& metaData = m_bulk->mesh_meta_data();
      stk::mesh::BulkData& bulkData = *m_bulk;

      // FIXME consider caching the coords_field
      CoordinatesFieldType *coords_field = metaData.get_field<double>(stk::topology::NODE_RANK, "coordinates");

      PerceptMesh meshUtil(&metaData, &bulkData);

      IsInElement isIn(input_phy_points, found_parametric_coordinates);

      // check first using the hint
      if (bulkData.is_valid(hint_element))
        {
          IsInElement isIn_hint(input_phy_points, found_parametric_coordinates);
          isIn_hint(hint_element,  bulkData);

          //if (EXTRA_PRINT) std::cout << "SimpleSearcher::findElement: hint found it= " << isIn_hint.m_found_it << std::endl;
          if (isIn_hint.m_found_it)
            {
              found_it = 1;
              return isIn_hint.m_foundElement;
            }
        }

      elementOpLoop(*meshUtil.get_bulk_data(), isIn, coords_field);
      //if (EXTRA_PRINT) std::cout << "SimpleSearcher::findElement: found it= " << isIn.m_found_it << std::endl;

      if (isIn.m_found_it)
        {
          found_it = 1;
          return isIn.m_foundElement;
        }
      else
        {
          found_it = 0;
        }
      return stk::mesh::Entity();
    }
  }
