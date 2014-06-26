#ifndef stk_encr_SimpleSearcher_hpp
#define stk_encr_SimpleSearcher_hpp

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

#include <stk_percept/function/internal/Searcher.hpp>
#include <stk_percept/function/internal/IsInElement.hpp>

namespace stk_classic
{
  namespace percept
  {
    typedef mesh::Field<double>                     ScalarFieldType ;
    typedef mesh::Field<double, mesh::Cartesian>    VectorFieldType ;

    class FieldFunction;

    class SimpleSearcher : public Searcher
    {
      stk_classic::mesh::BulkData *m_bulk;
    public:
      SimpleSearcher(stk_classic::mesh::BulkData *bulk);

      virtual ~SimpleSearcher() {}

      /**
       *  Dimensions of input_phy_points = ([P]=1, [D]) 
       *  Dimensions of found_parametric_coordinates = ([P]=1, [D])
       */

      virtual const stk_classic::mesh::Entity *findElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates, 
                                                   unsigned& found_it, const mesh::Entity *hint_element );
    };

  }
}

#endif
