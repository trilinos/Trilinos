#ifndef stk_encr_Searcher_hpp
#define stk_encr_Searcher_hpp

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <stk_percept/function/Function.hpp>

namespace stk
{
  namespace percept
  {
    
    class Searcher
    {
    public:

      /** Find the element containing this physical point and return if found (also set the found_it flag to 1, else 0).
       *  If hint_element is non-null, use it to check first if it contains the point to potentially avoid a more costly search.
       * 
       *  Dimensions of input_phy_points = ([P]=1, [D]) 
       *  Dimensions of found_parametric_coordinates = ([P]=1, [D])
       */
      virtual const stk::mesh::Entity *findElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates, 
                                                   unsigned& found_it, const mesh::Entity *hint_element )=0;
      virtual void setupSearch() {}
      virtual void tearDownSearch() {}
      virtual ~Searcher() {}
    };
  }
}
#endif
