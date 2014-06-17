#ifndef stk_encr_IsInElement_hpp
#define stk_encr_IsInElement_hpp

#include <cmath>
#include <math.h>
#include <vector>
#include <utility>


#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>

#include <Shards_CellTopology.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/function/Function.hpp>
#include <stk_percept/function/ElementOp.hpp>

namespace stk_classic
{
  namespace percept
  {

    class IsInElement : public ElementOp
    {

    public:
      bool m_found_it;
      MDArray& m_input_phy_points;
      MDArray& m_found_parametric_coordinates;
      const stk_classic::mesh::Entity *m_foundElement;
    public:
      IsInElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates);

      bool operator()(const stk_classic::mesh::Entity& element, const mesh::BulkData& bulkData);
      bool operator()(const stk_classic::mesh::Entity& element, stk_classic::mesh::FieldBase* field, const mesh::BulkData& bulkData) 
      {
        return (*this)(element, bulkData);
      }
      void init_elementOp();
      void fini_elementOp();

    private:


      /**
       *  Dimensions of input_phy_points = ([P]=1, [D]) 
       *  Dimensions of found_parametric_coordinates = ([P]=1, [D])
       */
      void isInElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates, unsigned& found_it, const mesh::Entity& element,
                       const mesh::BulkData& bulkData);

    };

  }
}
#endif
