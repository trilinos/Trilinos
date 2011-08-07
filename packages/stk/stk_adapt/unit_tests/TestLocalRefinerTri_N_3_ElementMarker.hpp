#ifndef stk_adapt_TestLocalRefinerTri_N_3_ElementMarker_hpp
#define stk_adapt_TestLocalRefinerTri_N_3_ElementMarker_hpp

#include <stk_adapt/ElementMarker.hpp>

namespace stk {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * A test implementation as a use case for ElementMarker
     */
    class TestLocalRefinerTri_N_3_ElementMarker : public ElementMarker
    {
    public:
      TestLocalRefinerTri_N_3_ElementMarker(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0);

      //virtual ElementUnrefineCollection  buildUnrefList();

    protected:

      /// Client supplies these methods - given an element, which edge, and the nodes on the edge, return instruction on what to do to the edge,
      ///    0 (nothing), 1 (refine)

      virtual int mark(const stk::mesh::Entity& element);

      ///    -1 (unrefine), 0 (nothing)
      virtual int markUnrefine(const stk::mesh::Entity& element);

    };

    // This is a very specialized test that is used in unit testing only (see unit_localRefiner/break_tri_to_tri_N_3_ElementMarker in UnitTestLocalRefiner.cpp)

    TestLocalRefinerTri_N_3_ElementMarker::TestLocalRefinerTri_N_3_ElementMarker(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk::mesh::FieldBase *proc_rank_field) : 
      ElementMarker(eMesh, bp, proc_rank_field)
    {
    }


    int TestLocalRefinerTri_N_3_ElementMarker::
    mark(const stk::mesh::Entity& element)
    {
      // vertical line position
      const double vx = 0.21;

      // horizontal line position
      const double vy = 1.21;

      const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

      for (unsigned inode=0; inode < elem_nodes.size(); inode++)
        {
          stk::mesh::Entity *node = elem_nodes[inode].entity();
          stk::mesh::Entity *node1 = elem_nodes[(inode + 1) % 3].entity();
          double *coord0 = stk::mesh::field_data( *m_eMesh.getCoordinatesField(), *node );
          double *coord1 = stk::mesh::field_data( *m_eMesh.getCoordinatesField(), *node1 );

          // choose to refine or not 
          if (
              ( std::fabs(coord0[0]-coord1[0]) > 1.e-3 &&
                ( (coord0[0] < vx && vx < coord1[0]) || (coord1[0] < vx && vx < coord0[0]) )
                )
              ||
              ( std::fabs(coord0[1]-coord1[1]) > 1.e-3 &&
                ( (coord0[1] < vy && vy < coord1[1]) || (coord1[1] < vy && vy < coord0[1]) )
                )
              )
            {
              return 1;
            }
        }

      return 0;
    }

    int TestLocalRefinerTri_N_3_ElementMarker::
    markUnrefine(const stk::mesh::Entity& element)
    {
      const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

      bool found = true;
      for (unsigned inode=0; inode < elem_nodes.size(); inode++)
        {
          stk::mesh::Entity *node = elem_nodes[inode].entity();
          double *coord = stk::mesh::field_data( *m_eMesh.getCoordinatesField(), *node );
          //if (coord[0] > 2.1 || coord[1] > 2.1)
          if (coord[0] > 1.0001 || coord[1] > 1.0001)
            {
              found = false;
              break;
            }
        }

      if (found)
        {
          return -1;
        }
      else
        return 0;
    }

  }
}
#endif
