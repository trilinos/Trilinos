
#include <stk_percept/Percept.hpp>


#include <Intrepid_CellTools.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>

#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/internal/SimpleSearcher.hpp>

namespace stk
{
  namespace percept
  {

    SimpleSearcher::SimpleSearcher(FieldFunction *ff) : m_fieldFunction(ff) {}

    /**
     *  Dimensions of input_phy_points = ([P]=1, [D]) 
     *  Dimensions of found_parametric_coordinates = ([P]=1, [D])
     */

    const stk::mesh::Entity *SimpleSearcher::findElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates, 
                                                         unsigned& found_it, const mesh::Entity *hint_element )
    {
      VERIFY_OP(input_phy_points.rank(), ==, found_parametric_coordinates.rank(), "SimpleSearcher::findElement bad dims");
      VERIFY_OP(input_phy_points.rank(), ==, 2, "SimpleSearcher::findElement bad rank");

      mesh::MetaData& metaData = m_fieldFunction->getField()->mesh_meta_data();
      mesh::BulkData& bulkData = *m_fieldFunction->getBulkData();
        
      // FIXME consider caching the coords_field 
      VectorFieldType *coords_field = metaData.get_field<VectorFieldType >("coordinates");

      PerceptMesh meshUtil(&metaData, &bulkData);

      IsInElement isIn(input_phy_points, found_parametric_coordinates);

      // check first using the hint
      if (hint_element)
        {
          IsInElement isIn_hint(input_phy_points, found_parametric_coordinates);
          isIn_hint(*hint_element, m_fieldFunction->getField(), bulkData);

          //if (EXTRA_PRINT) std::cout << "SimpleSearcher::findElement: hint found it= " << isIn_hint.m_found_it << std::endl;
          if (isIn_hint.m_found_it)
            {
              found_it = 1;
              return isIn_hint.m_foundElement;
            }
        }

      meshUtil.elementOpLoop(isIn, coords_field);
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
      return 0;
    }


  }
}

