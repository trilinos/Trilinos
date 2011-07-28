#ifndef stk_adapt_TestLocalRefinerTri_N_3_EdgeMarker_hpp
#define stk_adapt_TestLocalRefinerTri_N_3_EdgeMarker_hpp

#include <stk_adapt/EdgeMarker.hpp>

namespace stk {
  namespace adapt {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * A test implementation as a use case for EdgeMarker
     */
    class TestLocalRefinerTri_N_3_EdgeMarker : public EdgeMarker
    {
    public:
      TestLocalRefinerTri_N_3_EdgeMarker(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0);

      //virtual ElementUnrefineCollection  buildUnrefList();

    protected:

      /// Client supplies these methods - given an element, which edge, and the nodes on the edge, return instruction on what to do to the edge,
      ///    0 (nothing), 1 (refine)

      virtual int mark(const stk::mesh::Entity& element, unsigned which_edge, stk::mesh::Entity & node0, stk::mesh::Entity & node1,
                           double *coord0, double *coord1, std::vector<int>& existing_edge_marks) ;

      ///    -1 (unrefine), 0 (nothing)
      virtual int markUnrefine(const stk::mesh::Entity& element, unsigned which_edge, stk::mesh::Entity & node0, stk::mesh::Entity & node1,
                                double *coord0, double *coord1);

    };

    // This is a very specialized test that is used in unit testing only (see unit_localRefiner/break_tri_to_tri_N_3_EdgeMarker in UnitTestLocalRefiner.cpp)

    TestLocalRefinerTri_N_3_EdgeMarker::TestLocalRefinerTri_N_3_EdgeMarker(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk::mesh::FieldBase *proc_rank_field) : 
      EdgeMarker(eMesh, bp, proc_rank_field)
    {
    }


    int TestLocalRefinerTri_N_3_EdgeMarker::
    mark(const stk::mesh::Entity& element, unsigned which_edge, stk::mesh::Entity & node0, stk::mesh::Entity & node1,
             double *coord0, double *coord1, std::vector<int>& existing_edge_marks)
    {
      // vertical line position
      const double vx = 0.21;

      // horizontal line position
      const double vy = 1.21;

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


      return 0;
    }

    int TestLocalRefinerTri_N_3_EdgeMarker::
    markUnrefine(const stk::mesh::Entity& element, unsigned which_edge, stk::mesh::Entity & node0, stk::mesh::Entity & node1,
                  double *coord0, double *coord1)
    {
      // choose to unrefine or not
      if (coord0[0] > 1.0001 || coord0[1] > 1.0001 || coord1[0] > 1.0001 || coord1[1] > 1.0001)
        {
          return 0;
        }
      else
        return -1;
    }

#if 0
    int TestLocalRefinerTri_N_3_EdgeMarker::
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

      return 0;
    }
#endif

  }
}
#endif
