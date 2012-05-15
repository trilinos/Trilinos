#ifndef stk_percept_STKSearcher_hpp
#define stk_percept_STKSearcher_hpp

#include <stk_percept/Percept.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/diag/IdentProc.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>


#include <Intrepid_FunctionSpaceTools.hpp>

#include <Shards_CellTopology.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_percept/function/FieldFunction.hpp>

#include <stk_percept/function/internal/Searcher.hpp>
#include <stk_percept/function/internal/BuildBoundingBoxes.hpp>

#define EXTRA_PRINT 0

namespace stk
{
  namespace percept
  {
    typedef mesh::Field<double>                     ScalarFieldType ;
    typedef mesh::Field<double, mesh::Cartesian>    VectorFieldType ;

    template<int SpatialDim>
    class STKSearcher : public Searcher
    {
      typedef BuildBoundingBoxes<SpatialDim> BBB;
      typedef typename BBB ::AABoundingBox BBox;
      typedef typename BBB ::BoundingPoint BPoint;
      typedef typename BBB ::IdentProc IdProc;
      typedef std::vector<std::pair<IdProc, IdProc> > IdentProcRelation;

      stk::mesh::BulkData *m_bulk;
      //std::vector< template<> BuildBoundingBoxes<SpatialDim>::BoundingBox > m_boxes;
      std::vector< BBox > m_boxes;

    public:

      STKSearcher(stk::mesh::BulkData *bulk);

      virtual ~STKSearcher() {}

      void setupSearch();

      void tearDownSearch();

      /**
       *  Dimensions of input_phy_points = ([P]=1, [D]) 
       *  Dimensions of found_parametric_coordinates = ([P]=1, [D])
       */

      virtual const stk::mesh::Entity *findElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates, 
                                                   unsigned& found_it, const mesh::Entity *hint_element );

    private:


    };

  }
}

#include "STKSearcherDef.hpp"

#endif
