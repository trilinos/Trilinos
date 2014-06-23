
#include <stk_percept/Percept.hpp>


#include <Intrepid_CellTools.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>

#include <stk_percept/function/FieldFunction.hpp>
#include <stk_percept/function/internal/SimpleSearcher.hpp>

namespace stk_classic
{
  namespace percept
  {

    SimpleSearcher::SimpleSearcher(stk_classic::mesh::BulkData *bulk) : m_bulk(bulk) {}

    /**
     *  Dimensions of input_phy_points = ([P]=1, [D])
     *  Dimensions of found_parametric_coordinates = ([P]=1, [D])
     */

    const stk_classic::mesh::Entity *SimpleSearcher::findElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates,
                                                         unsigned& found_it, const mesh::Entity *hint_element )
    {
      VERIFY_OP(input_phy_points.rank(), ==, found_parametric_coordinates.rank(), "SimpleSearcher::findElement bad dims");
      VERIFY_OP(input_phy_points.rank(), ==, 2, "SimpleSearcher::findElement bad rank");

      mesh::fem::FEMMetaData& metaData = stk_classic::mesh::fem::FEMMetaData::get( *m_bulk );
      mesh::BulkData& bulkData = *m_bulk;

      // FIXME consider caching the coords_field
      VectorFieldType *coords_field = metaData.get_field<VectorFieldType >("coordinates");

      PerceptMesh meshUtil(&metaData, &bulkData);

      IsInElement isIn(input_phy_points, found_parametric_coordinates);

      // check first using the hint
      if (hint_element)
        {
          IsInElement isIn_hint(input_phy_points, found_parametric_coordinates);
          isIn_hint(*hint_element,  bulkData);

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

