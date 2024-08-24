// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>

#include <percept/function/internal/IsInElement.hpp>


#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_FunctionSpaceTools.hpp>

#include <percept/norm/IntrepidManager.hpp>
#include <percept/FieldTypes.hpp>

  namespace percept
  {

    IsInElement::IsInElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates) :
      m_found_it(false), m_input_phy_points(input_phy_points), m_found_parametric_coordinates(found_parametric_coordinates),
      m_foundElement()
    {}

    void IsInElement::init_elementOp()
    {
      m_found_it=false;
    }
    void IsInElement::fini_elementOp()
    {
    }

    bool IsInElement::operator()(const stk::mesh::Entity element, const stk::mesh::BulkData& bulkData)
    {

      unsigned found_it;
      isInElement(m_input_phy_points, m_found_parametric_coordinates, found_it, element, bulkData);
      //if (EXTRA_PRINT) std::cout << "IsInElement::operator() found_it = " << found_it << std::endl;
      // break out of enclosing element loop
      if (found_it)
        {
          m_found_it = true;
          m_foundElement = element;
          return true;
        }
      else
        return false;
    }
      
    bool IsInElement::operator()(const stk::mesh::Entity element, stk::mesh::FieldBase* field, const stk::mesh::BulkData& bulkData)
    {
      return (*this)(element, bulkData);
    }

    /**
     *  Dimensions of input_phy_points = ([P]=1, [D])
     *  Dimensions of found_parametric_coordinates = ([P]=1, [D])
     */
    void IsInElement::isInElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates, unsigned& found_it, const stk::mesh::Entity element,
                                  const stk::mesh::BulkData& bulkData)
    {
      IntrepidManager::isInElement(input_phy_points, found_parametric_coordinates, found_it, element, bulkData);
    }
  }
